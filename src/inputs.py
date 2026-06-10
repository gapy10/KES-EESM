"""
Vhodni podatki za načrtovanje 6-polnega sinhronskega motorja z vzbujalnim navitjem.

Vsebuje:
    - razred MachineInputs: nazivni podatki iz seminarske naloge,
    - razred MaterialParams: lastnosti M330-35A in bakra,
    - razred DesignBounds: meje optimizacijskih spremenljivk za GA,
    - korelacijo q -> Qs, kw1 iz EMETOR (https://www.emetor.com/).

Vir:
    Pyrhönen, J., Jokinen, T., Hrabovcová, V.,
    "Design of Rotating Electrical Machines", 2nd ed., Wiley, 2014.
    Poglavje 2 in 6: nazivni izračun ter faktorji navitja.
"""

from __future__ import annotations

from dataclasses import dataclass, field, asdict
from typing import Optional

import math
from pathlib import Path

try:
    import yaml  # PyYAML — neobvezno, samo za branje konfiguracije
    _HAS_YAML = True
except ImportError:
    _HAS_YAML = False


# -----------------------------------------------------------------------------
# Nazivni podatki stroja (Sinhronc_opis.md, razdelek "Vhodni podatki")
# -----------------------------------------------------------------------------
@dataclass
class MachineInputs:
    # Topologija
    m: int = 3                  # število faz [-]
    p: int = 3                  # število polovih parov [-] (6-polni stroj)

    # Električni nazivni podatki
    U_f: float = 200.0          # nazivna fazna napetost, RMS [V]
    P_c: float = 50_000.0       # nazivna moč v vogalni točki [W] (analitično)
    P_c_min_FEMM: float = 50_000.0  # MINIMALNA P (= M·ω) ki jo mora design
                                    # dejansko doseči v FEMM (specifikacija
                                    # naloge: 50 kW pri 7000 rpm)
    n_c: float = 7000.0         # mehanska hitrost v vogalni točki [vrt/min]
    n_max: float = 14_000.0     # najvišja obratovalna hitrost [vrt/min]

    # Predpostavke obratovanja
    cos_phi: float = 0.95       # faktor delavnosti (predpostavka)
    eta_init: float = 0.95      # začetna ocena izkoristka za dimenzioniranje

    # Material
    material: str = "M330-35A"  # neorientirana elektropločevina

    @property
    def number_of_poles(self) -> int:
        """Število polov 2p."""
        return 2 * self.p

    @property
    def omega_mech(self) -> float:
        """Mehanska kotna hitrost v vogalni točki [rad/s]."""
        return 2.0 * math.pi * self.n_c / 60.0

    @property
    def omega_el(self) -> float:
        """Električna kotna hitrost v vogalni točki [rad/s]."""
        return self.omega_mech * self.p

    @property
    def freq_c(self) -> float:
        """Električna frekvenca v vogalni točki [Hz].
        f = n * p / 60  =  7000/60 * 3 = 350 Hz."""
        return self.n_c * self.p / 60.0

    @property
    def freq_max(self) -> float:
        """Električna frekvenca pri n_max [Hz] (= 700 Hz pri 14000 rpm)."""
        return self.n_max * self.p / 60.0

    @property
    def torque_c(self) -> float:
        """Nazivni navor v vogalni točki [Nm] = P_c / omega_mech."""
        return self.P_c / self.omega_mech


# -----------------------------------------------------------------------------
# Lastnosti materialov
# -----------------------------------------------------------------------------
@dataclass
class MaterialParams:
    # Železo M330-35A (Sura/Cogent)
    rho_fe: float = 7650.0          # gostota železa [kg/m^3]
    k_fe: float = 0.97              # polnilni faktor lameliranja [-]
    # Empirična ojačitvena faktorja izgub (iz izgube_fem.txt + Pyrhönen
    # str. 553-556) — pokrivata vpliv nehomogenosti B v zobu (sredina ima
    # manj, robovi več) in rotacijskega magnetenja v jarmu. Te faktorji
    # se uporabljajo TAKO v analitiki KOT v FEMM post-procesiranju, zato
    # ne ustvarjajo razlike med njima.
    k_fe_tooth: float = 2.0
    k_fe_yoke: float = 1.6

    # FEMM KALIBRACIJA — napove dejansko FEMM obnašanje na podlagi
    # opazovanj iz 5+ designov:
    #
    # 1) Faktor saturacije B v železu: analitika cilja B_ds=1.7, FEMM pa
    #    izmeri ~1.46 (razmerje 0.86), ker rotor-stator nelinearnost ne
    #    podpre projektirane gostote. Posledica: P_Fe v izgube_fem
    #    polinomu, ki je ~B², bo realno ~0.86² = 0.74-krat manjši.
    k_B_femm_factor: float = 0.86   # B_FEMM / B_target za predikcijo izgub

    # FEMM-verificiran navor na enoto statorskega toka glede na analitično
    # napoved. Idealni tok I_n_ideal = P_c/(η·m·cosφ·U_f) v FEMM proizvede le
    # ~0.96-kratnik nazivnega navora (analitika rahlo preceni navor na amper:
    # nasičenje, harmoniki, faktor navitja). Zato nazivni tok rahlo korigiramo
    # I_n = I_n_ideal / k_torque_femm (≈ +3.7 %), da stroj v FEMM DEJANSKO doda
    # nazivnih 50 kW (M_c = 68.2 Nm) — umerjeno tako, da FEMM navor doseže
    # nazivnega. To je kalibracija na POŠTENIH 50 kW, NE napihovanje (prejšnja
    # ~40 % korekcija, ki je silila navor 20–35 % nad nazivnega, je odstranjena).
    # U_ind ostane nespremenjen (prosti tek), J_cu_s ostane = genu (presek raste
    # z I), J_cu_r nedotaknjen — oba znotraj mej.
    k_torque_femm: float = 0.964

    # Baker
    sigma_cu_20: float = 34e6       # specifična prevodnost pri 20 °C [S/m]
    alpha_cu: float = 0.00381       # temperaturni koeficient [1/K]
    T_operating: float = 80.0       # delovna temperatura navitij [°C]
    T_reference: float = 20.0       # referenčna temperatura [°C]
    K_Cu_stator: float = 0.38       # faktor zapolnitve statorskega utora [-]
    K_Cu_rotor: float = 0.60        # faktor zapolnitve rotorskega utora [-] (assumed)

    # Magnetne omejitve (mehke; uporabljene kot začetne predpostavke)
    B_delta_init: float = 0.9       # ciljna B v zr. reži [T]
    B_tooth_max: float = 1.7        # max B v zobu [T] (saturacija ~ 1.8 T)
    B_yoke_max: float = 1.4         # max B v jarmu [T]

    # Tokovne omejitve (Sinhronc_opis.md razdelek A)
    J_cu_s_max: float = 10.0        # A/mm^2
    J_cu_r_max: float = 5.0         # A/mm^2

    def sigma_cu(self, T_celsius: Optional[float] = None) -> float:
        """Prevodnost bakra pri temperaturi T [S/m].
        sigma(T) = sigma(20) / (1 + alpha * (T - 20))
        """
        T = self.T_operating if T_celsius is None else T_celsius
        return self.sigma_cu_20 / (1.0 + self.alpha_cu * (T - self.T_reference))

    def resistivity_cu(self, T_celsius: Optional[float] = None) -> float:
        """Specifična upornost bakra [Ohm * m]."""
        return 1.0 / self.sigma_cu(T_celsius)


# -----------------------------------------------------------------------------
# EMETOR korelacija q -> (Qs, kw1)
# https://www.emetor.com/edit/windings/  (3-fazno, 2p = 6, single-layer / double-layer)
# Vrednosti za 6-polni 3-fazni stroj, izbranih je 5 dovoljenih q vrednosti.
# -----------------------------------------------------------------------------
EMETOR_Q_TABLE: dict[float, dict] = {
    1.0: {"Qs": 18,  "kw1": 0.9598},   # 6 polov × 3 faz × 1
    1.5: {"Qs": 27,  "kw1": 0.9452},   # frakcijsko navitje
    2.0: {"Qs": 36,  "kw1": 0.9659},   # standardno integralno
    2.5: {"Qs": 45,  "kw1": 0.9531},   # frakcijsko
    3.0: {"Qs": 54,  "kw1": 0.9598},   # integralno z dvojnim utorom
}
# Opomba: kw1 vrednosti so izračunane po standardnih obrazcih:
#   kw1 = kd1 * kp1 (faktor razdeljenosti × faktor koraka).
# Za 6-polni 3-fazni stroj se kw1 giblje med 0.94 in 0.97; nad 1.0 ni možno.


def q_to_qs_and_kw1(q: float) -> tuple[int, float]:
    """Vrne (Qs, kw1) za dovoljen q ∈ {1, 1.5, 2, 2.5, 3}.

    Args:
        q: število utorov na pol in fazo.

    Returns:
        (Qs, kw1) - število statorskih utorov in faktor navitja osnovne harmonske.

    Raises:
        ValueError če q ni v dovoljeni množici.
    """
    if q not in EMETOR_Q_TABLE:
        allowed = sorted(EMETOR_Q_TABLE.keys())
        raise ValueError(f"q = {q} ni dovoljen. Dovoljene vrednosti: {allowed}")
    entry = EMETOR_Q_TABLE[q]
    return entry["Qs"], entry["kw1"]


# -----------------------------------------------------------------------------
# Meje GA optimizacijskih spremenljivk
# -----------------------------------------------------------------------------
@dataclass
class DesignBounds:
    """Spodnje in zgornje meje za 9 optimizacijskih spremenljivk GA."""

    # Geometrija — meje izbrane za "dolge in vitke" stroje:
    #   * spodnja D_r = 0.145 m zagotavlja, da FEMM uspešno mreži najmanjše rešitve
    #   * zgornja D_r = 0.220 m omeji preveč debele rotorje
    #   * spodnja l/D = 0.9 sili dolžino (prej je GA padla na 0.4 → preploščati stroji)
    D_r_min: float = 0.145          # premer rotorja [m]
    D_r_max: float = 0.220
    l_to_D_min: float = 0.9         # l/D razmerje — minimum sili "kolonadno" geometrijo
    l_to_D_max: float = 2.5

    # Magnetne obremenitve
    B_delta_min: float = 0.6        # [T]
    B_delta_max: float = 1.05
    B_ds_min: float = 1.2           # B v statorskem zobu [T]
    B_ds_max: float = 1.7
    B_sy_min: float = 1.0           # B v statorskem jarmu [T]
    B_sy_max: float = 1.5

    # Tokovne gostote
    J_cu_s_min: float = 3.0         # [A/mm^2]
    J_cu_s_max: float = 10.0
    J_cu_r_min: float = 2.0
    J_cu_r_max: float = 5.0

    # Vzbujanje
    N_r_min: int = 10               # ovojev na rotorski pol
    N_r_max: int = 80

    # Zračna reža (10. optimizacijska spremenljivka)
    delta_min: float = 0.7e-3       # [m]
    delta_max: float = 1.5e-3       # [m]

    # q (diskretno)
    q_choices: list[float] = field(default_factory=lambda: [1.0, 1.5, 2.0, 2.5, 3.0])


# -----------------------------------------------------------------------------
# Pomožne funkcije za nalaganje YAML konfiguracije
# -----------------------------------------------------------------------------
def load_from_yaml(path: str | Path) -> tuple[MachineInputs, MaterialParams, DesignBounds]:
    """Naloži vse tri dataclassove iz YAML datoteke.

    Datoteka lahko vsebuje delne nadgradnje; manjkajoča polja ostanejo privzeta.
    """
    path = Path(path)
    if not path.exists():
        return MachineInputs(), MaterialParams(), DesignBounds()
    if not _HAS_YAML:
        raise ImportError(
            "PyYAML ni nameščen. Zaženi: pip install pyyaml  (ali uporabi privzete vrednosti)."
        )

    with path.open("r", encoding="utf-8") as f:
        cfg = yaml.safe_load(f) or {}

    machine = MachineInputs(**cfg.get("machine", {}))
    material = MaterialParams(**cfg.get("material", {}))

    bounds_cfg = cfg.get("bounds", {})
    bounds = DesignBounds(**bounds_cfg) if bounds_cfg else DesignBounds()

    return machine, material, bounds


def to_dict(*objs) -> dict:
    """Pripravi združen slovar za logiranje."""
    out: dict = {}
    for obj in objs:
        if hasattr(obj, "__dataclass_fields__"):
            out[type(obj).__name__] = asdict(obj)
    return out


if __name__ == "__main__":
    m = MachineInputs()
    mat = MaterialParams()
    print(f"P_c = {m.P_c/1000:.1f} kW")
    print(f"n_c = {m.n_c} rpm  ->  omega_mech = {m.omega_mech:.2f} rad/s")
    print(f"f_c = {m.freq_c:.1f} Hz  (vogalna točka)")
    print(f"f_max = {m.freq_max:.1f} Hz  (pri n_max)")
    print(f"Nazivni navor M_c = {m.torque_c:.2f} Nm")
    print(f"sigma_Cu pri {mat.T_operating} C = {mat.sigma_cu()/1e6:.2f} MS/m")
    for q in [1.0, 1.5, 2.0, 2.5, 3.0]:
        Qs, kw1 = q_to_qs_and_kw1(q)
        print(f"  q = {q}  ->  Qs = {Qs:2d},  kw1 = {kw1:.4f}")
