"""Testi za Fazo 4 — post-procesiranje.

Testi ne potrebujejo FEMM zagona; uporabljajo sintetične NoLoadResult /
TorqueResult dataclasse.
"""

import numpy as np
import pytest

from src.inputs import MachineInputs, MaterialParams
from src.losses import IronLossModel
from src.analytical import MotorDesignGenes, analyze
from src.femm_sim import NoLoadResult, TorqueResult
from src.post import (
    femm_iron_losses,
    uind_spectrum,
    torque_thd,
    compare,
    summary_dataframe,
)


@pytest.fixture(scope="module")
def baseline_design():
    machine = MachineInputs()
    material = MaterialParams()
    loss = IronLossModel.fit_default()
    g = MotorDesignGenes(D_r=0.150, l_to_D=0.8, B_delta=0.9, B_ds=1.6,
                         B_sy=1.3, J_cu_s=6.0, J_cu_r=4.0, N_r=40, q=2.0)
    return analyze(g, machine, material, loss), machine, material, loss


def _synthetic_no_load() -> NoLoadResult:
    theta = np.linspace(0, 120, 13)  # 13 pozicij
    flux = 0.05 * np.cos(np.deg2rad(theta) * 3.0)  # 3 elektr. period v 120° mech
    # U_ind = -dPsi/dt; dPsi cca sinus
    omega = 2 * np.pi * 7000 / 60
    dt = np.diff(np.deg2rad(theta)) / omega
    uind = np.diff(flux) / dt
    return NoLoadResult(
        theta_deg_mech=theta,
        flux_A=flux,
        U_ind_inst=uind,
        U_ind_amp=float(np.max(np.abs(uind))),
        U_ind_rms=float(np.sqrt(np.mean(uind ** 2))),
        B_max_tooth=1.5,
        B_max_yoke=1.05,
    )


def _synthetic_torque(M_amp=68.0, ripple_amp=2.0) -> TorqueResult:
    theta = np.linspace(0, 120, 13)
    # M = M_amp · sin(3·θ_mech) + ripple · sin(6·θ_mech)
    rad = np.deg2rad(theta)
    M = M_amp * np.sin(3 * rad) + ripple_amp * np.sin(6 * rad)
    sig = M[:-1]
    spec = np.fft.rfft(sig)
    amps = np.abs(spec) / len(sig)
    if len(amps) > 1:
        amps[1:-1] *= 2.0
    return TorqueResult(
        theta_deg_mech=theta, torque=M,
        M_avg=float(np.mean(M)),
        M_max=float(np.max(np.abs(M))),
        M_fundamental=float(amps[1]) if len(amps) > 1 else 0.0,
        M_second=float(amps[2]) if len(amps) > 2 else 0.0,
        M_ripple_pp=float(M.max() - M.min()),
        M_ripple_percent=0.0,
        fft_harmonics=np.arange(len(amps)),
        fft_amplitudes=amps,
    )


class TestFemmIronLosses:
    def test_returns_positive_losses(self, baseline_design):
        design, machine, material, loss = baseline_design
        nl = _synthetic_no_load()
        fl = femm_iron_losses(design, nl, machine, material, loss)
        assert fl.P_fe_tooth > 0
        assert fl.P_fe_yoke > 0
        assert fl.P_fe_total == pytest.approx(fl.P_fe_tooth + fl.P_fe_yoke)
        assert 0.5 < fl.eta_femm < 1.0

    def test_efficiency_in_expected_range(self, baseline_design):
        design, machine, material, loss = baseline_design
        nl = _synthetic_no_load()
        fl = femm_iron_losses(design, nl, machine, material, loss)
        # Za 50 kW stroj naj eta v območju [0.85, 0.99]
        assert 0.85 < fl.eta_femm < 0.99


class TestUindSpectrum:
    def test_returns_correct_shape(self):
        nl = _synthetic_no_load()
        us = uind_spectrum(nl)
        assert us.harmonics.shape == us.amplitudes.shape
        assert us.fundamental >= 0
        assert us.thd_percent >= 0

    def test_pure_sinusoid_has_low_thd(self):
        theta = np.linspace(0, 120, 25)
        omega = 2 * np.pi * 7000 / 60
        rad = np.deg2rad(theta) * 3.0   # 3 polov para = osnovno
        flux = 0.05 * np.cos(rad)
        dt = np.diff(np.deg2rad(theta)) / omega
        uind = np.diff(flux) / dt
        nl = NoLoadResult(
            theta_deg_mech=theta, flux_A=flux,
            U_ind_inst=uind,
            U_ind_amp=float(np.max(np.abs(uind))),
            U_ind_rms=float(np.sqrt(np.mean(uind ** 2))),
        )
        us = uind_spectrum(nl)
        assert us.thd_percent < 50.0  # numerični diff doda nekaj harmonikov, a < 50%


class TestTorqueThd:
    def test_pure_sinusoid_thd_zero(self):
        tq = _synthetic_torque(M_amp=68.0, ripple_amp=0.0)
        assert torque_thd(tq) == pytest.approx(0.0, abs=1.0)

    def test_with_ripple_thd_nonzero(self):
        tq = _synthetic_torque(M_amp=68.0, ripple_amp=10.0)
        thd = torque_thd(tq)
        # Z amplitudo 10 Nm 2. harmonske in 68 Nm 1. harm: THD ≈ 10/68 ≈ 15 %
        assert 10 < thd < 20


class TestCompare:
    def test_returns_8_rows(self, baseline_design):
        design, machine, material, loss = baseline_design
        nl = _synthetic_no_load()
        tq = _synthetic_torque()
        fl = femm_iron_losses(design, nl, machine, material, loss)
        rows = compare(design, nl, tq, fl, machine)
        # 8 vrstic: U_ind, B_zob, B_jarem, Navor M_c (cilj), P_Fe×3, η.
        # (Nekdanja ločena vrstica "Navor (napoved)" je odstranjena skupaj z
        #  navorno korekcijo toka — navor primerjamo le proti nazivnemu M_c.)
        assert len(rows) == 8

    def test_diff_percent_finite(self, baseline_design):
        design, machine, material, loss = baseline_design
        nl = _synthetic_no_load()
        tq = _synthetic_torque()
        fl = femm_iron_losses(design, nl, machine, material, loss)
        rows = compare(design, nl, tq, fl, machine)
        for r in rows:
            assert np.isfinite(r.diff_percent)


class TestSummaryDataframe:
    def test_one_row_per_design(self, baseline_design):
        design, machine, material, loss = baseline_design
        nl = _synthetic_no_load()
        tq = _synthetic_torque()
        fl = femm_iron_losses(design, nl, machine, material, loss)
        df = summary_dataframe(
            designs=[design, design],
            no_loads=[nl, nl],
            torques=[tq, tq],
            losses=[fl, fl],
            labels=["mid", "max_eta"],
            machine=machine,
        )
        assert len(df) == 2

    def test_required_columns_present(self, baseline_design):
        design, machine, material, loss = baseline_design
        nl = _synthetic_no_load()
        tq = _synthetic_torque()
        fl = femm_iron_losses(design, nl, machine, material, loss)
        df = summary_dataframe([design], [nl], [tq], [fl], ["mid"], machine)
        required = [
            "design_id", "label", "D_r_mm", "L_r_mm", "D_se_mm",
            "q", "Q_s", "N_s", "N_r", "I_n_A", "I_m_A",
            "P_Fe_analit_W", "P_Cu_W", "eta_analit",
            "P_Fe_femm_W", "eta_femm",
            "U_ind_femm_V", "M_femm_Nm", "M_2harm_Nm",
            "B_max_zob_T", "B_max_jarem_T", "thd_torque_percent",
            "V_active_cm3",
        ]
        for col in required:
            assert col in df.columns, f"Manjka stolpec: {col}"
