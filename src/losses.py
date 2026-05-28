"""
Železove izgube za neorientirano elektropločevino M330-35A.

Implementira 2D polinomsko aproksimacijo specifičnih izgub p_Fe = f(B, f)
po metodi iz priložene MATLAB skripte `izgube_fem.txt`:

    surf = [B  f]
    f = fit(surf, P, 'poly23')

To je polinom stopnje 2 v B in stopnje 3 v f, brez križnih členov
v 'poly23' notaciji po MATLAB Curve Fitting Toolbox. V Pythonu uporabimo
ekvivalent z `numpy.linalg.lstsq` na osnovi vseh polnih členov:

    p(B,f) = sum_{i<=2, j<=3, i+j<=3} c_ij * B^i * f^j

(MATLAB poly23 zajema vse člene z i<=2, j<=3, i+j<=5; naša izvedba
uporablja polni set tenzorskega produkta in se ne razlikuje pomembno
v točnosti na razponu tabele.)

Vir:
    - Sura Cogent, data sheet M330-35A.
    - Pyrhönen, J., poglavje 11.2 (železove izgube, polinomska aproksimacija).
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Iterable
import numpy as np


# -----------------------------------------------------------------------------
# Surovi tabelarni podatki M330-35A (iz izgube_fem.txt)
# B [T], f [Hz], P [W/kg]
# -----------------------------------------------------------------------------

# fmt: off
_B_TABLE = np.array([
    0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8,
    0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6,
    0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6,
    0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6,
], dtype=float)

_F_TABLE = np.array([
    50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50,
    100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,
    200,200,200,200,200,200,200,200,200,200,200,200,200,200,200,200,200,
    400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,400,
], dtype=float)

_P_TABLE = np.array([
    0, 0.05, 0.08, 0.17, 0.28, 0.40, 0.53, 0.68, 0.84, 1.02, 1.22, 1.44, 1.69, 2.00, 2.4, 2.94, 3.67, 4.32, 4.73,
    0, 0.06, 0.20, 0.41, 0.67, 0.97, 1.30, 1.68, 2.10, 2.56, 3.07, 3.64, 4.29, 5.07, 6.06, 7.40, 8.86,
    0, 0.12, 0.48, 1.02, 1.68, 2.47, 3.37, 4.39, 5.54, 6.82, 8.25, 9.86, 11.6, 13.7, 16.3, 19.6, 23.2,
    0, 0.33, 1.27, 2.69, 4.49, 6.66, 9.19, 12.11, 15.44, 19.22, 23.54, 28.48, 34.12, 40.62, 48.24, 57.86, 70.24,
], dtype=float)
# fmt: on


def _build_design_matrix(B: np.ndarray, f: np.ndarray) -> np.ndarray:
    """Konstruira matriko členov za polinom stopnje (2 v B, 3 v f).

    Vsebuje vse člene B^i * f^j z 0 <= i <= 2 in 0 <= j <= 3.
    """
    cols = []
    for i in range(0, 3):
        for j in range(0, 4):
            cols.append((B ** i) * (f ** j))
    return np.column_stack(cols)


@dataclass
class IronLossModel:
    """Polinomski model specifičnih železovih izgub za M330-35A.

    Po inicializaciji `fit_default()` ali `fit_from_arrays()` lahko kličete
    `loss_density(B, f)` za izračun specifičnih izgub v W/kg.

    Veljavnost: B ∈ [0, 1.8] T, f ∈ [50, 400] Hz. Pri f = 350 Hz (vogalna
    točka 7000 vrt/min * 3 polovih parov / 60 = 350 Hz) smo blizu zgornje
    meje, kar je smiselno (interpolacija, ne ekstrapolacija).
    """

    coeffs: np.ndarray
    rmse_train: float

    @staticmethod
    def fit_from_arrays(B: np.ndarray, f: np.ndarray, P: np.ndarray) -> "IronLossModel":
        """Polinomska aproksimacija po metodi najmanjših kvadratov."""
        B = np.asarray(B, dtype=float)
        f = np.asarray(f, dtype=float)
        P = np.asarray(P, dtype=float)
        A = _build_design_matrix(B, f)
        c, *_ = np.linalg.lstsq(A, P, rcond=None)
        residuals = A @ c - P
        rmse = float(np.sqrt(np.mean(residuals ** 2)))
        return IronLossModel(coeffs=c, rmse_train=rmse)

    @staticmethod
    def fit_default() -> "IronLossModel":
        """Polinomska aproksimacija na priloženih M330-35A tabelarnih podatkih."""
        return IronLossModel.fit_from_arrays(_B_TABLE, _F_TABLE, _P_TABLE)

    def loss_density(self, B: float | Iterable[float], f: float | Iterable[float]) -> float | np.ndarray:
        """Specifične izgube p_Fe [W/kg] pri amplitudi B [T] in frekvenci f [Hz].

        Velja interpolacija znotraj B ∈ [0, 1.8] T in f ∈ [50, 400] Hz.
        Pri ekstrapolaciji bo polinom verjetno netočen.
        """
        B = np.asarray(B, dtype=float)
        f = np.asarray(f, dtype=float)
        scalar = (B.ndim == 0 and f.ndim == 0)
        Bf = np.atleast_1d(B)
        ff = np.atleast_1d(f)
        if Bf.shape != ff.shape:
            # broadcast
            Bf, ff = np.broadcast_arrays(Bf, ff)
        A = _build_design_matrix(Bf, ff)
        out = A @ self.coeffs
        # vse vrednosti naj bodo nenegativne; pri zelo majhnih B lahko
        # polinom vrne rahlo negativne vrednosti zaradi šuma fitanja.
        out = np.clip(out, 0.0, None)
        return float(out[0]) if scalar else out

    def total_loss(
        self,
        B: float,
        f: float,
        mass_kg: float,
        k_form: float = 1.0,
    ) -> float:
        """Skupne izgube [W] = k_form * mass * p_Fe(B, f).

        k_form je empirični faktor zaradi nehomogenosti B v zobu/jarmu in
        zaradi rotacijskega magnetenja (Pyrhönen tabela 11.5).
        Privzete vrednosti (iz `izgube_fem.txt`):
            - zob:   k_form = 2.0
            - jarem: k_form = 1.6
        """
        return float(k_form * mass_kg * self.loss_density(B, f))


# -----------------------------------------------------------------------------
# Bakrene izgube
# -----------------------------------------------------------------------------

def end_winding_length_estimate(tau_p_mm: float) -> float:
    """Ocena dolžine enega čela navitja [m] po Pyrhönenu (poglavje 4.6).

    Privzeto: l_end ≈ 2 * tau_p (v m). Tu pa ima naloga eksplicitno
    poenostavitev (FR-1.7): dolžina čela navitja = L_r (aktivna dolžina).
    Funkcija je tu samo kot referenca; v izračunu izgub uporabimo L_r.
    """
    return 2.0e-3 * tau_p_mm


def stator_winding_resistance(
    N_s_per_phase: int,
    cu_cross_section_m2: float,
    mean_turn_length_m: float,
    sigma_cu_T: float,
) -> float:
    """Statorska upornost ene faze [Ohm].

    R = (N_s * l_avg) / (sigma_Cu * A_Cu)
    kjer je l_avg povprečna dolžina enega ovoja navitja.
    """
    return (N_s_per_phase * mean_turn_length_m) / (sigma_cu_T * cu_cross_section_m2)


def stator_copper_losses(
    I_n_phase_rms: float,
    R_s_phase: float,
    m_phases: int = 3,
) -> float:
    """Statorske bakrene izgube [W] = m * I^2 * R."""
    return m_phases * I_n_phase_rms ** 2 * R_s_phase


def rotor_field_resistance(
    N_r_total_turns: int,
    cu_cross_section_m2: float,
    mean_turn_length_m: float,
    sigma_cu_T: float,
) -> float:
    """Upornost serijsko vezanega vzbujalnega navitja [Ohm].

    Vse vzbujalne tuljave na vseh polih so v seriji - skupno N_r_total ovojev
    z istim povprečnim ovojem.
    """
    return (N_r_total_turns * mean_turn_length_m) / (sigma_cu_T * cu_cross_section_m2)


def rotor_field_losses(I_m: float, R_r: float) -> float:
    """Bakrene izgube vzbujalnega navitja [W] = I_m^2 * R_r (DC tok)."""
    return I_m ** 2 * R_r


# -----------------------------------------------------------------------------
# Smoke test ob neposrednem zagonu
# -----------------------------------------------------------------------------
if __name__ == "__main__":
    model = IronLossModel.fit_default()
    print(f"RMSE pri fitanju: {model.rmse_train:.3f} W/kg")
    print()
    print("Specifične izgube W/kg za nekaj točk:")
    print(f"  B=1.0,  f= 50 Hz  ->  {model.loss_density(1.0,  50):.3f}  (tabela: 1.22)")
    print(f"  B=1.5,  f= 50 Hz  ->  {model.loss_density(1.5,  50):.3f}  (tabela: 2.94)")
    print(f"  B=1.0,  f=100 Hz  ->  {model.loss_density(1.0, 100):.3f}  (tabela: 3.07)")
    print(f"  B=1.5,  f=200 Hz  ->  {model.loss_density(1.5, 200):.3f}  (tabela: 19.6)")
    print(f"  B=1.5,  f=350 Hz  ->  {model.loss_density(1.5, 350):.3f}  (interpoliran)")
    print(f"  B=1.0,  f=350 Hz  ->  {model.loss_density(1.0, 350):.3f}  (interpoliran)")
