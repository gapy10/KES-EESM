"""Smoke testi za Fazo 1 (analitični izračun)."""

import math
import pytest

from src.inputs import MachineInputs, MaterialParams, q_to_qs_and_kw1
from src.losses import IronLossModel, _B_TABLE, _F_TABLE, _P_TABLE
from src.analytical import MotorDesignGenes, analyze


# -----------------------------------------------------------------------------
# inputs.py
# -----------------------------------------------------------------------------
class TestInputs:
    def test_nominal_machine_params(self):
        m = MachineInputs()
        assert m.m == 3
        assert m.p == 3
        assert m.U_f == 200.0
        assert m.P_c == 50_000.0
        assert m.n_c == 7000.0

    def test_freq_c_is_350Hz(self):
        m = MachineInputs()
        assert m.freq_c == pytest.approx(350.0)

    def test_torque_at_corner(self):
        m = MachineInputs()
        # M = P / omega_mech;  omega = 2π·7000/60
        omega = 2 * math.pi * 7000 / 60
        expected = 50_000 / omega
        assert m.torque_c == pytest.approx(expected)
        assert m.torque_c == pytest.approx(68.21, rel=1e-3)

    def test_sigma_cu_temperature(self):
        mat = MaterialParams()
        # Pri T=20 °C velja sigma = sigma_0
        assert mat.sigma_cu(20.0) == pytest.approx(34e6, rel=1e-6)
        # Pri T=80 °C: sigma = sigma_0 / (1 + α·60) ≈ 27.7 MS/m
        s_80 = mat.sigma_cu(80.0)
        assert 27e6 < s_80 < 28e6

    @pytest.mark.parametrize("q,Qs", [(1.0, 18), (1.5, 27), (2.0, 36), (2.5, 45), (3.0, 54)])
    def test_q_to_Qs(self, q, Qs):
        Qs_calc, kw1 = q_to_qs_and_kw1(q)
        assert Qs_calc == Qs
        assert 0.94 < kw1 < 0.97   # tipično za 6-polni 3-fazni stroj

    def test_q_invalid_raises(self):
        with pytest.raises(ValueError):
            q_to_qs_and_kw1(1.7)


# -----------------------------------------------------------------------------
# losses.py
# -----------------------------------------------------------------------------
class TestLosses:
    def test_fit_default_returns_model(self):
        m = IronLossModel.fit_default()
        assert m.coeffs.shape[0] > 0
        # RMSE pod 2 W/kg na celem območju 0–70 W/kg.
        assert m.rmse_train < 2.0

    @pytest.mark.parametrize("B,f,P_table,tol_abs", [
        (1.0,  50, 1.22, 0.3),
        (1.5,  50, 2.94, 0.3),
        (1.0, 100, 3.07, 0.5),
        (1.5, 200, 19.6, 1.5),
    ])
    def test_table_points_within_tolerance(self, B, f, P_table, tol_abs):
        m = IronLossModel.fit_default()
        p = m.loss_density(B, f)
        assert p == pytest.approx(P_table, abs=tol_abs)

    def test_loss_density_nonnegative(self):
        m = IronLossModel.fit_default()
        for B in [0, 0.1, 0.5, 1.0, 1.5]:
            for f in [50, 100, 200, 350, 400]:
                p = m.loss_density(B, f)
                assert p >= 0.0, f"Negativna izguba pri B={B}, f={f}: {p}"

    def test_total_loss_uses_k_form(self):
        m = IronLossModel.fit_default()
        mass = 5.0
        p1 = m.total_loss(B=1.0, f=350, mass_kg=mass, k_form=1.0)
        p2 = m.total_loss(B=1.0, f=350, mass_kg=mass, k_form=2.0)
        assert p2 == pytest.approx(2.0 * p1)


# -----------------------------------------------------------------------------
# analytical.py
# -----------------------------------------------------------------------------
class TestAnalytical:
    @pytest.fixture(scope="class")
    def baseline(self):
        machine = MachineInputs()
        material = MaterialParams()
        loss = IronLossModel.fit_default()
        genes = MotorDesignGenes(
            D_r=0.150, l_to_D=0.8, B_delta=0.9,
            B_ds=1.6, B_sy=1.3, J_cu_s=6.0, J_cu_r=4.0,
            N_r=40, q=2.0,
        )
        return analyze(genes, machine, material, loss)

    def test_geometry_positive(self, baseline):
        d = baseline
        assert d.L_r > 0
        assert d.tau_p > 0
        assert d.delta > 0
        assert d.b_ds > 0 and d.h_ds > 0 and d.h_ys > 0

    def test_Qs_consistent_with_q(self, baseline):
        # q = 2, p = 3, m = 3  =>  Q_s = 2·p·m·q = 36
        assert baseline.Q_s == 36

    def test_Zq_is_integer_and_positive(self, baseline):
        assert isinstance(baseline.Z_q, int)
        assert baseline.Z_q > 0

    def test_efficiency_reasonable(self, baseline):
        # Za 50 kW stroj pričakujemo η med 0.85 in 0.98.
        assert 0.85 < baseline.eta < 0.98

    def test_current_densities_match_genes(self, baseline):
        # J_actual mora biti enako J_gene (po izračunu preseka).
        assert baseline.J_cu_s_actual == pytest.approx(baseline.genes.J_cu_s, rel=1e-3)
        assert baseline.J_cu_r_actual == pytest.approx(baseline.genes.J_cu_r, rel=1e-3)

    def test_feasible_baseline(self, baseline):
        assert baseline.feasible, f"Pričakovano feasible, kršitve: {baseline.infeasibility_reasons}"

    def test_infeasible_when_l_to_D_too_small(self):
        machine = MachineInputs()
        material = MaterialParams()
        loss = IronLossModel.fit_default()
        genes = MotorDesignGenes(
            D_r=0.300, l_to_D=0.2,    # preploščat
            B_delta=0.9, B_ds=1.6, B_sy=1.3,
            J_cu_s=6.0, J_cu_r=4.0, N_r=40, q=2.0,
        )
        d = analyze(genes, machine, material, loss)
        assert not d.feasible
        assert any("preplo" in r.lower() or "l/d" in r.lower() for r in d.infeasibility_reasons)
