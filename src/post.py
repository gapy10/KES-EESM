"""
Post-procesiranje rezultatov FEMM simulacij (Faza 4).

Implementira:
    - izračun železovih izgub iz FEMM B_max + nazivne frekvence,
    - izračun izkoristka iz FEMM rezultatov,
    - FFT inducirane napetosti, THD,
    - THD navorne karakteristike (1. + višje harmonske),
    - primerjavo analitičnih in FEMM vrednosti (tabela %-odstopanj),
    - povzetek vseh 5 Pareto rešitev v eno DataFrame / CSV.

Modul ne kliče FEMM — sprejme `NoLoadResult` in `TorqueResult` (že izračunane).
Glavna referenca: Pyrhönen §11 (železove izgube), §6.4 (THD).
"""

from __future__ import annotations

from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Optional, Sequence

import numpy as np
import pandas as pd

from .analytical import MotorDesign
from .inputs import MachineInputs, MaterialParams
from .losses import IronLossModel
from .femm_sim import NoLoadResult, TorqueResult


# =============================================================================
# Železove izgube iz FEMM
# =============================================================================

@dataclass
class FemmLossBreakdown:
    """Železove + bakrene izgube preračunane iz FEMM rezultatov."""
    P_fe_tooth: float
    P_fe_yoke: float
    P_fe_total: float
    P_cu_total: float       # iz analitičnega (FEMM ne računa AC bakrenih izgub)
    P_loss_total: float
    eta_femm: float
    P_in_femm: float        # P_c + P_loss


def femm_iron_losses(
    design: MotorDesign,
    no_load: NoLoadResult,
    machine: MachineInputs,
    material: MaterialParams,
    loss_model: IronLossModel,
) -> FemmLossBreakdown:
    """Izračuna železove izgube iz B_max iz FEMM in P_Cu iz analitičnega.

    P_Fe = k_form · m · p_Fe(B_max, f_c)

    Bakrene izgube ostanejo iz analitičnega izračuna (FEMM ne računa AC
    odvisnih izgub v navitju brez full transient solver).
    """
    f_c = machine.freq_c
    P_fe_tooth = loss_model.total_loss(
        B=no_load.B_max_tooth, f=f_c,
        mass_kg=design.m_fe_tooth_s,
        k_form=material.k_fe_tooth,
    )
    P_fe_yoke = loss_model.total_loss(
        B=no_load.B_max_yoke, f=f_c,
        mass_kg=design.m_fe_yoke_s,
        k_form=material.k_fe_yoke,
    )
    P_fe_total = P_fe_tooth + P_fe_yoke
    P_cu_total = design.P_cu_total
    P_loss_total = P_fe_total + P_cu_total
    P_in = machine.P_c + P_loss_total
    eta_femm = machine.P_c / P_in if P_in > 0 else 0.0
    return FemmLossBreakdown(
        P_fe_tooth=P_fe_tooth,
        P_fe_yoke=P_fe_yoke,
        P_fe_total=P_fe_total,
        P_cu_total=P_cu_total,
        P_loss_total=P_loss_total,
        eta_femm=eta_femm,
        P_in_femm=P_in,
    )


# =============================================================================
# FFT inducirane napetosti
# =============================================================================

@dataclass
class UindSpectrum:
    """Spekter inducirane napetosti."""
    harmonics: np.ndarray         # indeksi harmonikov 0..N/2
    amplitudes: np.ndarray        # amplitude [V]
    fundamental: float            # amplituda osnovne harmonske [V]
    thd_percent: float            # THD = √(Σ_{n≥2} A_n^2) / A_1 [%]


def uind_spectrum(no_load: NoLoadResult) -> UindSpectrum:
    """Izračuna FFT spekter U_ind(t).

    Vhodni signal je U_ind_inst (np.diff(Ψ_A) / dt). Predpostavljamo, da je
    enakomerno vzorčen čez en električni period (= 360°/p mech).
    """
    sig = np.asarray(no_load.U_ind_inst)
    if sig.size == 0:
        return UindSpectrum(np.array([0]), np.array([0.0]), 0.0, 0.0)

    spec = np.fft.rfft(sig)
    amps = np.abs(spec) / len(sig)
    if len(amps) > 1:
        amps[1:-1] *= 2.0
    harmonics = np.arange(len(amps))
    fundamental = float(amps[1]) if len(amps) > 1 else 0.0
    # THD = sqrt(sum A_n^2 for n>=2) / A_1
    higher = amps[2:] if len(amps) > 2 else np.array([0.0])
    thd = float(np.sqrt(np.sum(higher ** 2)) / fundamental * 100.0) if fundamental > 1e-6 else 0.0
    return UindSpectrum(harmonics, amps, fundamental, thd)


# =============================================================================
# THD navora
# =============================================================================

def torque_thd(tq: TorqueResult) -> float:
    """THD navorne karakteristike [%] = √(Σ_{n≥2} A_n^2) / A_1.

    Vir: Pyrhönen §6.7 (definicija THD za izmenične signale).
    """
    amps = np.asarray(tq.fft_amplitudes)
    if amps.size < 2:
        return 0.0
    A_1 = amps[1]
    if A_1 < 1e-6:
        return 0.0
    higher = amps[2:]
    return float(np.sqrt(np.sum(higher ** 2)) / A_1 * 100.0)


# =============================================================================
# Primerjava analitično vs FEMM
# =============================================================================

@dataclass
class AnalyticalVsFemm:
    """Primerjava ključnih količin (analitično, FEMM, % odstopanje)."""
    quantity: str
    unit: str
    analytical: float
    femm: float
    diff_percent: float

    def as_row(self) -> dict:
        return asdict(self)


def compare(
    design: MotorDesign,
    no_load: NoLoadResult,
    tq: TorqueResult,
    fem_loss: FemmLossBreakdown,
    machine: MachineInputs,
) -> list[AnalyticalVsFemm]:
    """Vrne seznam primerjav za poročilo."""
    def diff(a, f):
        return ((f - a) / a * 100.0) if abs(a) > 1e-12 else 0.0

    rows = [
        AnalyticalVsFemm(
            "U_ind (RMS)", "V",
            machine.U_f, no_load.U_ind_rms,
            diff(machine.U_f, no_load.U_ind_rms),
        ),
        AnalyticalVsFemm(
            "B_max v zobu", "T",
            design.genes.B_ds, no_load.B_max_tooth,
            diff(design.genes.B_ds, no_load.B_max_tooth),
        ),
        AnalyticalVsFemm(
            "B_max v jarmu", "T",
            design.genes.B_sy, no_load.B_max_yoke,
            diff(design.genes.B_sy, no_load.B_max_yoke),
        ),
        AnalyticalVsFemm(
            "Nazivni navor", "Nm",
            machine.torque_c, tq.M_fundamental,
            diff(machine.torque_c, tq.M_fundamental),
        ),
        AnalyticalVsFemm(
            "P_Fe (zobje)", "W",
            design.P_fe_tooth, fem_loss.P_fe_tooth,
            diff(design.P_fe_tooth, fem_loss.P_fe_tooth),
        ),
        AnalyticalVsFemm(
            "P_Fe (jarem)", "W",
            design.P_fe_yoke, fem_loss.P_fe_yoke,
            diff(design.P_fe_yoke, fem_loss.P_fe_yoke),
        ),
        AnalyticalVsFemm(
            "P_Fe (skupno)", "W",
            design.P_fe_total, fem_loss.P_fe_total,
            diff(design.P_fe_total, fem_loss.P_fe_total),
        ),
        AnalyticalVsFemm(
            "Izkoristek η", "%",
            design.eta * 100, fem_loss.eta_femm * 100,
            diff(design.eta, fem_loss.eta_femm),
        ),
    ]
    return rows


def comparison_table_text(rows: Sequence[AnalyticalVsFemm], title: str = "") -> str:
    """Lepa tekstovna tabela primerjave."""
    lines = []
    if title:
        lines.append(title)
        lines.append("=" * len(title))
    hdr = f"{'Količina':<18} {'Enota':<6} {'Analit.':>10} {'FEMM':>10} {'Δ%':>8}"
    lines.append(hdr)
    lines.append("-" * len(hdr))
    for r in rows:
        lines.append(
            f"{r.quantity:<18} {r.unit:<6} {r.analytical:>10.3f} "
            f"{r.femm:>10.3f} {r.diff_percent:>+7.1f}"
        )
    return "\n".join(lines)


# =============================================================================
# Povzetna tabela za vseh 5 rešitev
# =============================================================================

def summary_dataframe(
    designs: Sequence[MotorDesign],
    no_loads: Sequence[NoLoadResult],
    torques: Sequence[TorqueResult],
    losses: Sequence[FemmLossBreakdown],
    labels: Sequence[str],
    machine: MachineInputs,
    design_ids: Optional[Sequence[str]] = None,
    pareto_idx: Optional[Sequence[int]] = None,
    replaced: Optional[Sequence[bool]] = None,
) -> pd.DataFrame:
    """Sestavi povzetno tabelo za 5 izbranih rešitev.

    Stolpci natanko po FSD razdelek 6.3.1 (`results.csv` shema), z dodatnimi
    FEMM stolpci. Neobvezno: `design_ids` (npr. "D05_min_V") in `pareto_idx`
    omogočita zapis izvirne pozicije na Pareto fronti — uporabno, ko je
    bila rešitev zamenjana s sosedom (`replaced=True`).
    """
    rows = []
    for k, (d, nl, tq, fl, lab) in enumerate(zip(designs, no_loads, torques, losses, labels), start=1):
        did = design_ids[k - 1] if design_ids is not None else f"D{k:02d}_{lab}"
        row = {
            "design_id": did,
            "label": lab,
            "D_r_mm": d.genes.D_r * 1e3,
            "L_r_mm": d.L_r * 1e3,
            "D_se_mm": d.D_se * 1e3,
            "q": d.genes.q,
            "Q_s": d.Q_s,
            "N_s": d.N_s,
            "N_r": d.genes.N_r,
            "I_n_A": d.I_n,
            "I_m_A": d.I_m,
            # analitično:
            "P_Fe_analit_W": d.P_fe_total,
            "P_Cu_W": d.P_cu_total,
            "eta_analit": d.eta,
            "U_ind_analit_V": machine.U_f,
            "M_analit_Nm": machine.torque_c,
            # FEMM:
            "P_Fe_femm_W": fl.P_fe_total,
            "eta_femm": fl.eta_femm,
            "U_ind_femm_V": nl.U_ind_rms,
            "M_femm_Nm": tq.M_fundamental,
            "M_2harm_Nm": tq.M_second,
            "M_ripple_pp_Nm": tq.M_ripple_pp,
            "B_max_zob_T": nl.B_max_tooth,
            "B_max_jarem_T": nl.B_max_yoke,
            "thd_torque_percent": torque_thd(tq),
            "V_active_cm3": d.V_active * 1e6,
        }
        if pareto_idx is not None:
            row["pareto_idx"] = int(pareto_idx[k - 1])
        if replaced is not None:
            row["replaced"] = bool(replaced[k - 1])
        rows.append(row)
    return pd.DataFrame(rows)


# =============================================================================
# Pretvorba iz JSON nazaj v dataklasi (za --reprocess način)
# =============================================================================

def load_femm_results(json_path: str | Path) -> list[dict]:
    """Naloži femm_pareto_results.json za ponovno analizo brez FEMM zagona."""
    import json
    return json.loads(Path(json_path).read_text(encoding="utf-8"))
