"""
FEMM simulacijski moduli: prosti tek (točka L) in navor (točka M).

Implementira:
    L: simulacija prostega teka — I_stator = 0, vzbujanje = I_m, vrtenje
       rotorja za en polov par; izračun Ψ_A(θ) → U_ind pri n_c.
    M: simulacija navora — I_A = √2·I_n, I_B = -I_A/2, I_C = -I_A/2,
       vzbujanje = I_m, vrtenje za en polov par; mo_blockintegral(22).
       FFT navorne karakteristike.

Vse funkcije sprejmejo `MotorDesign` in delajo s pyFEMM. Predpostavljajo,
da je FEMM že odprt in da je model naložen.

Vir:
    Meeker, D., FEMM 4.2 — `mo_blockintegral(22)` = Steady-state weighted
    stress tensor torque about (0, 0).
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import femm
import numpy as np

from .analytical import MotorDesign


# =============================================================================
# Rezultati simulacij
# =============================================================================

@dataclass
class NoLoadResult:
    """Rezultat simulacije prostega teka."""
    theta_deg_mech: np.ndarray       # mehanski kot rotorja [°]
    flux_A: np.ndarray               # Ψ_A(θ) [Wb]
    U_ind_inst: np.ndarray           # trenutna U_ind (interpolirana) [V]
    U_ind_amp: float                 # amplituda U_ind [V]
    U_ind_rms: float                 # RMS vrednost U_ind [V]
    B_max_tooth: float = 0.0         # max B v statorskem zobu [T]
    B_max_yoke: float = 0.0          # max B v statorskem jarmu [T]


@dataclass
class TorqueResult:
    """Rezultat simulacije navora.

    Opomba o interpretaciji:
        Med rotacijo rotorja za en električen period in fiksnimi statorskimi
        tokovi je navor sinusoidalen z osnovno frekvenco = 1 cikel/period.
        Zato je `M_avg` (povprečje) ≈ 0, `M_fundamental` pa amplituda 1.
        harmonske, ki ustreza nazivnemu navoru (= maksimum sinusoide).
        `M_max` je dejanski maksimum vzorčenega signala (lahko manjši od
        `M_fundamental`, če vzorčevanje ne zadane vrha).
    """
    theta_deg_mech: np.ndarray       # mehanski kot rotorja [°]
    torque: np.ndarray               # M(θ) [Nm]
    M_avg: float                     # povprečni navor (DC) [Nm]
    M_max: float                     # max(|navor|) [Nm]
    M_fundamental: float             # amplituda 1. harmonske FFT [Nm] (= nazivni navor)
    M_second: float                  # amplituda 2. harmonske FFT [Nm] (reluktančna)
    M_ripple_pp: float               # peak-to-peak [Nm]
    M_ripple_percent: float          # % valovitosti glede na M_fundamental
    fft_harmonics: np.ndarray        # indeksi harmonikov [-]
    fft_amplitudes: np.ndarray       # amplitude FFT [Nm]


# =============================================================================
# Vrtenje rotorja
# =============================================================================

def _rotate_rotor(deg: float) -> None:
    """Zavrti skupino 1 (rotor) za podan kot v mehanskih stopinjah."""
    if abs(deg) < 1e-9:
        return
    femm.mi_selectgroup(1)
    femm.mi_moverotate(0, 0, deg)
    femm.mi_clearselected()


# =============================================================================
# Prosti tek (točka L)
# =============================================================================

def no_load(
    design: MotorDesign,
    *,
    rpm: Optional[float] = None,
    step_deg: float = 5.0,
    save_fem_path: Optional[str | Path] = None,
) -> NoLoadResult:
    """Simulira prosti tek: vzbujanje = I_m, statorski tok = 0, vrtenje za 360°/p.

    Args:
        design: stroj iz `analyze()`.
        rpm: mehanska hitrost za pretvorbo Ψ(θ) → U_ind(t). Privzeto n_c.
        step_deg: kotni inkrement vrtenja [mech °].
        save_fem_path: če podan, shrani končni `.fem`.

    Returns:
        NoLoadResult z Ψ_A(θ), U_ind(θ) in amplitudo/RMS.

    Predpostavlja, da je FEMM že odprt in model naložen.
    """
    p = 3  # 6-polni stroj
    if rpm is None:
        rpm = 7000.0  # privzeta vogalna hitrost; je v MachineInputs

    # Postavi vzbujanje in nič stator:
    femm.mi_setcurrent("DC", design.I_m)
    femm.mi_setcurrent("A", 0.0)
    femm.mi_setcurrent("B", 0.0)
    femm.mi_setcurrent("C", 0.0)

    # Vrtenje za 360°/p mech (en električen period) v korakih step_deg.
    span = 360.0 / p
    n_steps = int(round(span / step_deg)) + 1
    angles = np.linspace(0.0, span, n_steps)

    psi_A = np.zeros(n_steps)

    # Prva analiza pri kotu 0 (brez rotacije):
    femm.mi_analyze()
    femm.mi_loadsolution()
    _, _, psi = femm.mo_getcircuitproperties("A")
    psi_A[0] = psi

    # Vsak korak: zavrti za step_deg, ponovno analiziraj.
    for k in range(1, n_steps):
        delta = angles[k] - angles[k - 1]
        _rotate_rotor(delta)
        femm.mi_analyze()
        femm.mi_loadsolution()
        _, _, psi = femm.mo_getcircuitproperties("A")
        psi_A[k] = psi

    # Izračun U_ind(θ) = dΨ/dt.
    # dt = dθ_mech / ω_mech ; ω_mech = 2π·rpm/60
    omega_mech = 2.0 * math.pi * rpm / 60.0
    dtheta_rad = np.deg2rad(np.diff(angles))
    dt = dtheta_rad / omega_mech
    U_ind = np.diff(psi_A) / dt   # n_steps-1 vrednosti

    U_ind_amp = float(np.max(np.abs(U_ind)))
    U_ind_rms = float(np.sqrt(np.mean(U_ind ** 2)))

    # Izračun max B v statorskem zobu in jarmu (po vzoru izgube_fem.txt):
    B_tooth, B_yoke = sample_b_field(design)

    # Vrni rotor v začetno pozicijo (negativna kumulativna rotacija):
    _rotate_rotor(-span)

    if save_fem_path:
        femm.mi_saveas(str(save_fem_path))

    return NoLoadResult(
        theta_deg_mech=angles,
        flux_A=psi_A,
        U_ind_inst=U_ind,
        U_ind_amp=U_ind_amp,
        U_ind_rms=U_ind_rms,
        B_max_tooth=B_tooth,
        B_max_yoke=B_yoke,
    )


# =============================================================================
# Navor (točka M)
# =============================================================================

def torque(
    design: MotorDesign,
    *,
    step_deg: float = 5.0,
    save_fem_path: Optional[str | Path] = None,
) -> TorqueResult:
    """Simulira nazivni navor: vzbujanje + statorski tok, vrtenje za en polov par.

    Statorski tok (po vzoru primer1):
        I_A = √2 · I_n  (amplituda)
        I_B = -I_A / 2
        I_C = -I_A / 2

    Vrne časovni potek navora in FFT.
    """
    p = 3
    I_amp = math.sqrt(2.0) * design.I_n
    femm.mi_setcurrent("DC", design.I_m)
    femm.mi_setcurrent("A", I_amp)
    femm.mi_setcurrent("B", -I_amp / 2.0)
    femm.mi_setcurrent("C", -I_amp / 2.0)

    span = 360.0 / p
    n_steps = int(round(span / step_deg)) + 1
    angles = np.linspace(0.0, span, n_steps)
    torques = np.zeros(n_steps)

    femm.mi_analyze()
    femm.mi_loadsolution()
    femm.mo_groupselectblock(1)
    torques[0] = femm.mo_blockintegral(22)
    femm.mo_clearblock()

    for k in range(1, n_steps):
        delta = angles[k] - angles[k - 1]
        _rotate_rotor(delta)
        femm.mi_analyze()
        femm.mi_loadsolution()
        femm.mo_groupselectblock(1)
        torques[k] = femm.mo_blockintegral(22)
        femm.mo_clearblock()

    M_avg = float(np.mean(torques))
    M_max = float(np.max(np.abs(torques)))
    M_ripple_pp = float(torques.max() - torques.min())

    # FFT: zavržemo zadnjo točko (= prva po enem električnem periodu)
    sig = torques[:-1]
    spec = np.fft.rfft(sig)
    amps = np.abs(spec) / len(sig)
    if len(amps) > 1:
        amps[1:-1] *= 2.0
    harmonics = np.arange(len(amps))
    M_fundamental = float(amps[1]) if len(amps) > 1 else 0.0
    M_second = float(amps[2]) if len(amps) > 2 else 0.0
    M_ripple_percent = (
        float(M_ripple_pp / M_fundamental * 100.0) if M_fundamental > 1e-6 else 0.0
    )

    _rotate_rotor(-span)
    if save_fem_path:
        femm.mi_saveas(str(save_fem_path))

    return TorqueResult(
        theta_deg_mech=angles,
        torque=torques,
        M_avg=M_avg,
        M_max=M_max,
        M_fundamental=M_fundamental,
        M_second=M_second,
        M_ripple_pp=M_ripple_pp,
        M_ripple_percent=M_ripple_percent,
        fft_harmonics=harmonics,
        fft_amplitudes=amps,
    )


# =============================================================================
# Pomožno: max B v statorju
# =============================================================================

def sample_b_field(design: MotorDesign) -> tuple[float, float]:
    """Po vzoru izgube_fem.txt — vzorči B v eni točki na zobu in jarmu okrog
    celotnega oboda statorja, vrne maksimuma.

    Predpostavlja, da je FEMM v post-procesnem načinu (mo_*).
    """
    Q_s = design.Q_s
    D_r_mm = design.genes.D_r * 1e3
    delta_mm = design.delta * 1e3
    h_ds_mm = design.h_ds * 1e3
    h_ys_mm = design.h_ys * 1e3

    Rr = D_r_mm / 2.0
    Rsn = Rr + delta_mm
    # Točka na sredini zoba:
    r_tooth = Rsn + h_ds_mm / 2.0
    # Točka v jarmu:
    r_yoke = Rsn + h_ds_mm + h_ys_mm / 2.0

    B_max_tooth = 0.0
    B_max_yoke = 0.0
    n_samples = Q_s

    # Vzorčimo SREDINO STATORSKIH ZOB (med dvema sosednjima utoroma):
    # utori so centrirani na kotih 90° + k·360°/Q_s; zobje so na vmesnih kotih.
    for k in range(n_samples):
        ang = math.pi / 2.0 + (k + 0.5) * 2.0 * math.pi / Q_s
        cx = math.cos(ang)
        sy = math.sin(ang)
        bx, by = femm.mo_getpointvalues(r_tooth * cx, r_tooth * sy)[1:3]
        B = math.hypot(bx, by)
        if B > B_max_tooth:
            B_max_tooth = B
        bx, by = femm.mo_getpointvalues(r_yoke * cx, r_yoke * sy)[1:3]
        B = math.hypot(bx, by)
        if B > B_max_yoke:
            B_max_yoke = B
    return B_max_tooth, B_max_yoke
