"""
Analitični izračun 6-polnega sinhronskega motorja z vzbujalnim navitjem.

Implementira točke A–F seminarske naloge:
    A: izračun celotne geometrije iz vhodnih parametrov,
    B: B_δ(I_m, θ_m),
    C: polinomska aproksimacija izgub (modul `losses.py`),
    D: P_Fe v statorskih zobeh in jarmu,
    E: P_Cu v statorju in rotorju,
    F: izkoristek η.

Glavni vstop:
    analyze(genes, machine, material, loss_model) -> MotorDesign

Vsi izračuni so v SI enotah (m, A, T, V). Pretvorba v mm je potrebna le
ob klicih FEMM (modul `femm_model.py`).

Glavna referenca:
    Pyrhönen, J., Jokinen, T., Hrabovcová, V.,
    "Design of Rotating Electrical Machines", 2nd ed., Wiley, 2014.
    Citati v obliki "Pyrhönen §X.Y" se nanašajo na ta učbenik.

Sekundarna referenca:
    primer1/python/KES_geometrija.py — letošnja referenčna implementacija
    (6-polni stroj kolega). Empirične izbire (npr. 4e-7 magnetna konstanta
    v Carterjevi formuli) so prevzete od tam zaradi konsistence z lansko
    nalogo, vendar so dodatno obrazložene s sklicem na knjigo.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Optional

from .inputs import MachineInputs, MaterialParams, q_to_qs_and_kw1
from .losses import (
    IronLossModel,
    stator_winding_resistance,
    stator_copper_losses,
    rotor_field_resistance,
    rotor_field_losses,
)


MU_0 = 4.0 * math.pi * 1e-7  # magnetna konstanta vakuuma [H/m]


# -----------------------------------------------------------------------------
# Vhodni parametri (geni za GA)
# -----------------------------------------------------------------------------
@dataclass
class MotorDesignGenes:
    """Devet optimizacijskih spremenljivk, ki definirajo en stroj."""

    D_r: float          # premer rotorja [m]
    l_to_D: float       # razmerje L_r / D_r [-]
    B_delta: float      # ciljna amplituda B v zr. reži [T]
    B_ds: float         # ciljna B v statorskem zobu [T]
    B_sy: float         # ciljna B v statorskem jarmu [T]
    J_cu_s: float       # tokovna gostota statorskega navitja [A/mm^2]
    J_cu_r: float       # tokovna gostota rotorskega navitja [A/mm^2]
    N_r: int            # ovojev na rotorski pol [-]
    q: float            # utori/(pol·fazo) ∈ {1, 1.5, 2, 2.5, 3} [-]

    def as_tuple(self) -> tuple:
        return (self.D_r, self.l_to_D, self.B_delta, self.B_ds, self.B_sy,
                self.J_cu_s, self.J_cu_r, self.N_r, self.q)


# -----------------------------------------------------------------------------
# Polni rezultat analitičnega izračuna
# -----------------------------------------------------------------------------
@dataclass
class MotorDesign:
    """Analitično izračunan stroj. Vse dimenzije v SI (m), vse mase v kg."""

    # --- vhodni geni ---
    genes: MotorDesignGenes

    # --- osnovna geometrija ---
    L_r: float = 0.0            # aktivna dolžina paketa [m]
    tau_p: float = 0.0          # polov korak v zračni reži [m]
    delta: float = 0.0          # minimalna zračna reža (sredina pola) [m]
    K_c: float = 1.0            # Carterjev koeficient (skupni) [-]

    # --- statorske dimenzije ---
    D_si: float = 0.0           # notranji premer statorja = D_r + 2·delta [m]
    D_se: float = 0.0           # zunanji premer statorja [m]
    Q_s: int = 0                # število statorskih utorov [-]
    kw1: float = 0.0            # faktor navitja osnovne harmonske [-]
    N_s: int = 0                # ovojev/fazo (serijsko) [-]
    Z_q: int = 0                # ovojev v enem utoru [-]
    tau_u: float = 0.0          # statorski utorni korak v zr. reži [m]
    b_ds: float = 0.0           # širina statorskega zoba [m]
    b_s_avg: float = 0.0        # povprečna širina statorskega utora [m]
    h_ds: float = 0.0           # višina statorskega zoba (utora) [m]
    h_ys: float = 0.0           # višina statorskega jarma [m]
    b_1s: float = 0.0           # širina statorske utorne odprtine [m]

    # --- rotorske dimenzije ---
    b_dr: float = 0.0           # širina rotorskega zoba (osnova) [m]
    h_yr: float = 0.0           # višina rotorskega jarma [m]
    b_1r: float = 0.0           # širina rotorske utorne odprtine [m]

    # --- elektromagnetne količine ---
    Phi_pole: float = 0.0       # max fluks na pol [Wb]
    I_n: float = 0.0            # nazivni statorski fazni tok, RMS [A]
    I_m: float = 0.0            # vzbujalni tok rotorja, DC [A]
    F_m: float = 0.0            # magnetna napetost (ampere-ovojev) [A]

    # --- preseki bakra ---
    S_cu_s: float = 0.0         # presek enega statorskega vodnika [m^2]
    S_cu_r: float = 0.0         # presek enega rotorskega vodnika [m^2]
    S_us: float = 0.0           # površina enega statorskega utora [m^2]
    S_ur: float = 0.0           # površina enega rotorskega utora [m^2]

    # --- upornosti pri T_op ---
    R_s_phase: float = 0.0      # upornost ene faze statorja [Ohm]
    R_r_total: float = 0.0      # upornost vseh serijsko vezanih rotorskih tuljav [Ohm]

    # --- mase ---
    m_fe_tooth_s: float = 0.0   # masa vseh statorskih zob [kg]
    m_fe_yoke_s: float = 0.0    # masa statorskega jarma [kg]

    # --- izgube in izkoristek ---
    P_fe_tooth: float = 0.0     # izgube v statorskih zobeh [W]
    P_fe_yoke: float = 0.0      # izgube v statorskem jarmu [W]
    P_fe_total: float = 0.0     # skupne izgube v železu [W]
    P_cu_stator: float = 0.0    # izgube v statorskem navitju [W]
    P_cu_rotor: float = 0.0     # izgube v rotorskem (vzbujalnem) navitju [W]
    P_cu_total: float = 0.0     # skupne bakrene izgube [W]
    P_loss_total: float = 0.0   # P_Fe + P_Cu [W]
    eta: float = 0.0            # izkoristek [-]

    # --- aktivni volumen (cilj GA optimizacije) ---
    V_active: float = 0.0       # π * (D_se/2)^2 * L_r [m^3]

    # --- gostote toka po izračunu (preverba) ---
    J_cu_s_actual: float = 0.0  # A/mm^2
    J_cu_r_actual: float = 0.0  # A/mm^2

    # --- izvedljivost in razlog za morebitno zavrnitev ---
    feasible: bool = True
    infeasibility_reasons: list[str] = field(default_factory=list)


# =============================================================================
# Glavni izračun
# =============================================================================

def analyze(
    genes: MotorDesignGenes,
    machine: MachineInputs,
    material: MaterialParams,
    loss_model: IronLossModel,
    *,
    slot_opening_stator_mm: float = 4.0,
    slot_opening_rotor_frac: float = 0.5,
    saturation_factor: float = 1.10,
    end_winding_length_eq_L: bool = True,
) -> MotorDesign:
    """Polni analitični izračun enega stroja.

    Args:
        genes: optimizacijske spremenljivke (vhod).
        machine: nazivni podatki stroja.
        material: lastnosti M330-35A in bakra.
        loss_model: aproksimacija specifičnih železovih izgub.
        slot_opening_stator_mm: širina statorske utorne odprtine (b_1s) [mm].
            Privzeto 4 mm (kot v primer1).
        slot_opening_rotor_frac: delež širine rotorskega vmesnega utora
            (b_1r) glede na τ_p; privzeto 0.5 (kot v primer1, b_1r ≈ τ_p - cevelj).
        saturation_factor: k_sat za Carterjev preračun (Pyrhönen §3.5).
        end_winding_length_eq_L: če True, dolžina čela navitij = L_r (zahteva
            naloge, FR-1.7).

    Returns:
        polni MotorDesign s `feasible` zastavico in seznamom kršitev.
    """
    d = MotorDesign(genes=genes)
    reasons: list[str] = []

    # -- 1) Osnovna geometrija aktivnega dela --------------------------------
    D_r = genes.D_r                                # [m]
    L_r = genes.l_to_D * D_r                       # [m]
    d.L_r = L_r
    d.tau_p = math.pi * D_r / (2.0 * machine.p)    # Pyrhönen eq. 2.1

    # -- 2) Q_s in kw1 iz EMETOR -----------------------------------------------
    try:
        Q_s, kw1 = q_to_qs_and_kw1(genes.q)
    except ValueError as e:
        reasons.append(str(e))
        d.feasible = False
        d.infeasibility_reasons = reasons
        return d
    d.Q_s = Q_s
    d.kw1 = kw1

    # -- 3) Magnetni in tokovni preračun (Pyrhönen §6.1, eq. 6.4–6.6) --------
    # Začetna ocena Carter koeficienta (iterativno popravljeno spodaj).
    K_c = 1.17
    d.K_c = K_c

    # Tangencialna obremenitev (specific torque shear):
    #   σ_F = M / (π * D_r^2 * L_r / 2)
    V_r = math.pi * (D_r / 2.0) ** 2 * L_r
    sigma_F = machine.torque_c / (2.0 * V_r)        # [N/m^2]

    # Tokovna obloga (linear current loading), peak vrednost:
    #   A = √2 * σ_F / (B_δ * cosφ)
    # (Pyrhönen eq. 6.5 — RMS variant; tu uporabimo peak po primer1.)
    A_strom = math.sqrt(2.0) * sigma_F / (genes.B_delta * machine.cos_phi)  # [A/m]

    # Začetna zr. reža iz Carter zveze (primer1 vrstice 39–40):
    delta = (1.0 / K_c) * 4e-7 * d.tau_p * A_strom / genes.B_delta          # [m]
    # Pojasnilo konstante 4e-7: Pyrhönen §3.5 podaja
    #   δ ≥ (γ · D · A · k_sat) / (B_δ · ...)
    # kjer γ ≈ 4e-7 m·m/A za zračne reže v sinhronskih strojih (tipično).
    # Empirična izbira; popravljena z faktorjem nasičenja spodaj.

    delta *= saturation_factor                                              # k_sat
    d.delta = delta

    # -- 4) Carterjev koeficient (popravek za odprtine utorov) ---------------
    # Statorska utorna odprtina:
    b_1s = slot_opening_stator_mm * 1e-3            # [m]
    d.b_1s = b_1s
    K_bs = _carter_K_b(b_1s, delta)
    bes = K_bs * b_1s
    tau_u_air = math.pi * (D_r + 2.0 * delta) / Q_s     # utorni korak v zr. reži, statorska stran
    K_cs = tau_u_air / (tau_u_air - bes) if tau_u_air > bes else 1.5
    d.tau_u = tau_u_air

    # Rotorska utorna odprtina (med poloma):
    b_1r = slot_opening_rotor_frac * d.tau_p        # [m] — približek
    d.b_1r = b_1r
    K_br = _carter_K_b(b_1r, delta)
    ber = K_br * b_1r
    K_cr = d.tau_p / (d.tau_p - ber) if d.tau_p > ber else 1.5

    K_c = K_cs * K_cr
    d.K_c = K_c

    # Ponovni izračun zr. reže z novim K_c (eno iteracijo zadošča):
    delta = (1.0 / K_c) * 4e-7 * d.tau_p * A_strom / genes.B_delta * saturation_factor
    delta = max(delta, 0.3e-3)  # tehnološka spodnja meja 0.3 mm
    d.delta = delta

    # -- 5) Statorsko navitje: N_s, Z_q ---------------------------------------
    # E1f (RMS) ≈ √2 · π · f · k_w1 · N_s · α_p · B_δ · τ_p · L_r
    # Rešeno za N_s, ob predpostavki E ≈ U_f:
    alpha_p = 2.0 / math.pi                          # razmerje pole-arc/pole-pitch (sinusoidno polje)
    N_s_float = (math.sqrt(2.0) * machine.U_f) / (
        machine.omega_el * alpha_p * genes.B_delta * d.tau_p * L_r * kw1
    )
    # Zaokrožimo Z_q (ovojev/utor) na celo število, nato preračunamo N_s:
    Z_q_float = 2.0 * machine.m * N_s_float / Q_s
    Z_q = max(1, int(math.ceil(Z_q_float)))
    N_s = Z_q * Q_s // (2 * machine.m)
    d.Z_q = Z_q
    d.N_s = N_s

    # Korigirana B_δ glede na zaokrožen N_s (manjši odstopki):
    B_delta_actual = (math.sqrt(2.0) * machine.U_f) / (
        machine.omega_el * alpha_p * N_s * d.tau_p * L_r * kw1
    )

    # -- 6) Statorski tok in presek vodnikov ----------------------------------
    # I_s_phase ≈ P / (m · η · cosφ · U_f)
    I_n = machine.P_c / (machine.eta_init * machine.m * machine.cos_phi * machine.U_f)
    d.I_n = I_n

    S_cu_s = (I_n / genes.J_cu_s) * 1e-6                # presek ene žice [m^2]
    d.S_cu_s = S_cu_s
    d.J_cu_s_actual = (I_n / (S_cu_s * 1e6))

    # -- 7) Geometrija statorskega zoba in jarma ------------------------------
    # B v zobu homogeno = B_δ * τ_u / (b_ds · k_fe)   (Pyrhönen §7.4)
    b_ds = (B_delta_actual * d.tau_u) / (genes.B_ds * material.k_fe)
    d.b_ds = b_ds

    # Površina utora potrebna za navitje:
    S_us = Z_q * S_cu_s / material.K_Cu_stator
    d.S_us = S_us

    # Pravokotna aproksimacija:
    b_s_top = d.tau_u - b_ds                # širina utora pri zračni reži (~ τ_u - b_ds)
    h_ds = S_us / max(b_s_top, 1e-4)
    d.h_ds = h_ds
    d.b_s_avg = b_s_top

    # Statorski jarem (vrne polovica fluksa skozi vsako polovico jarma):
    Phi_pole = alpha_p * B_delta_actual * d.tau_p * L_r
    d.Phi_pole = Phi_pole
    h_ys = Phi_pole / (2.0 * material.k_fe * L_r * genes.B_sy)
    d.h_ys = h_ys

    # Notranji in zunanji premer statorja:
    D_si = D_r + 2.0 * delta
    D_se = D_si + 2.0 * h_ds + 2.0 * h_ys
    d.D_si = D_si
    d.D_se = D_se
    d.V_active = math.pi * (D_se / 2.0) ** 2 * L_r

    # -- 8) Rotorska geometrija ----------------------------------------------
    # V primer1 je b_dr ≈ 45 mm, h_yr ≈ 22 mm za D_r = 190 mm (delež ~ 0.24, 0.12).
    # Tukaj pa uporabimo magnetni pretok in dovoljen B_dr za rotorski zob:
    B_dr = genes.B_ds                                   # uporabimo isto kot stator (predpostavka)
    b_dr = (B_delta_actual * d.tau_p * alpha_p) / (B_dr * material.k_fe)
    d.b_dr = b_dr

    # Rotorski jarem — pretok skozi polovico polovega pasu:
    h_yr = Phi_pole / (2.0 * material.k_fe * L_r * genes.B_sy)
    d.h_yr = h_yr

    # -- 9) Vzbujanje (točka B): I_m iz F_m = N_r * I_m ----------------------
    # F_m mora pokriti zr. režo: F_m = H_δ * δ_eff, kjer H_δ = B_δ / μ_0
    # in δ_eff = K_c * δ * k_sat.
    delta_eff = K_c * delta * saturation_factor
    F_m = (B_delta_actual / MU_0) * delta_eff           # [A] (ampere-ovojev na pol)
    d.F_m = F_m
    I_m = F_m / max(genes.N_r, 1)                       # tok pri N_r ovojih/pol
    d.I_m = I_m

    # Presek rotorskega vodnika:
    S_cu_r = (I_m / genes.J_cu_r) * 1e-6                # [m^2]
    d.S_cu_r = S_cu_r
    d.J_cu_r_actual = (I_m / (S_cu_r * 1e6))

    # Površina ene rotorske tuljave (en pol): 2 stranici dveh blokov.
    # Privzeto: 2 bloka po polu (kot v primer1, vrstici 272–290).
    # Potrebni navojni volumen na pol: N_r * S_cu_r / K_Cu_rotor
    S_ur = (genes.N_r * S_cu_r) / material.K_Cu_rotor
    d.S_ur = S_ur

    # -- 10) Upornosti in bakrene izgube --------------------------------------
    # Povprečna dolžina ovoja:
    #   l_avg = 2·L_r + 2·tau_u + 2·b_s_avg                              (stator)
    #   l_avg_rotor = 2·L_r + 2·b_dr + 2·širina_navitja                  (rotor)
    if end_winding_length_eq_L:
        l_end = L_r  # zahteva naloge, FR-1.7
    else:
        l_end = 2.0 * d.tau_p   # standardna ocena

    l_avg_stator = 2.0 * L_r + 2.0 * l_end
    R_s = stator_winding_resistance(
        N_s_per_phase=N_s,
        cu_cross_section_m2=S_cu_s,
        mean_turn_length_m=l_avg_stator,
        sigma_cu_T=material.sigma_cu(),
    )
    d.R_s_phase = R_s
    d.P_cu_stator = stator_copper_losses(I_n, R_s, m_phases=machine.m)

    # Rotorske tuljave: N_r ovojev na pol × 2p polov = 2·p·N_r serijsko vezanih ovojev.
    N_r_total = 2 * machine.p * genes.N_r
    l_avg_rotor = 2.0 * L_r + 2.0 * l_end   # poenostavitev po FR-1.7
    R_r = rotor_field_resistance(
        N_r_total_turns=N_r_total,
        cu_cross_section_m2=S_cu_r,
        mean_turn_length_m=l_avg_rotor,
        sigma_cu_T=material.sigma_cu(),
    )
    d.R_r_total = R_r
    d.P_cu_rotor = rotor_field_losses(I_m, R_r)
    d.P_cu_total = d.P_cu_stator + d.P_cu_rotor

    # -- 11) Železove izgube --------------------------------------------------
    # Mase: zobje + jarem.
    # Volumen vseh statorskih zob: Q_s * b_ds * h_ds * L_r * k_fe
    V_tooth = Q_s * b_ds * h_ds * L_r * material.k_fe
    V_yoke = math.pi * L_r * ((D_se / 2.0) ** 2 - (D_si / 2.0 + h_ds) ** 2) * material.k_fe
    V_yoke = max(V_yoke, 0.0)
    m_tooth = V_tooth * material.rho_fe
    m_yoke = V_yoke * material.rho_fe
    d.m_fe_tooth_s = m_tooth
    d.m_fe_yoke_s = m_yoke

    # Frekvenca v vogalni točki:
    f_c = machine.freq_c

    d.P_fe_tooth = loss_model.total_loss(
        B=genes.B_ds, f=f_c, mass_kg=m_tooth, k_form=material.k_fe_tooth
    )
    d.P_fe_yoke = loss_model.total_loss(
        B=genes.B_sy, f=f_c, mass_kg=m_yoke, k_form=material.k_fe_yoke
    )
    d.P_fe_total = d.P_fe_tooth + d.P_fe_yoke

    # -- 12) Izkoristek -------------------------------------------------------
    d.P_loss_total = d.P_fe_total + d.P_cu_total
    d.eta = machine.P_c / (machine.P_c + d.P_loss_total)

    # -- 13) Preverjanje izvedljivosti ---------------------------------------
    if d.J_cu_s_actual > material.J_cu_s_max + 1e-3:
        reasons.append(f"J_cu_s = {d.J_cu_s_actual:.1f} > {material.J_cu_s_max} A/mm^2")
    if d.J_cu_r_actual > material.J_cu_r_max + 1e-3:
        reasons.append(f"J_cu_r = {d.J_cu_r_actual:.1f} > {material.J_cu_r_max} A/mm^2")
    if b_ds <= 0 or h_ds <= 0 or h_ys <= 0:
        reasons.append("Nepozitivna geometrija (b_ds, h_ds ali h_ys ≤ 0)")
    if b_ds + b_s_top > d.tau_u * 1.01:
        reasons.append("Zob + utor presega utorni korak")
    if genes.l_to_D < 0.4:
        reasons.append(f"l/D = {genes.l_to_D:.2f} < 0.4 (preploščat motor)")
    if d.delta < 0.3e-3:
        reasons.append(f"delta = {d.delta*1e3:.2f} mm < 0.3 mm (tehnološko)")
    if d.eta < 0.5:
        reasons.append(f"eta = {d.eta:.3f} nerealno nizek")

    d.feasible = len(reasons) == 0
    d.infeasibility_reasons = reasons
    return d


def _carter_K_b(b_1: float, delta: float) -> float:
    """Carterjev pomožni faktor (Pyrhönen eq. 3.62)."""
    if delta <= 0:
        return 0.0
    r = b_1 / delta
    return r / (5.0 + r)


# =============================================================================
# Pomožni izpis za hitri pregled
# =============================================================================

def pretty_print(d: MotorDesign) -> str:
    g = d.genes
    lines = [
        f"=== MotorDesign (feasible = {d.feasible}) ===",
        f"Vhod (gen):",
        f"  D_r = {g.D_r*1e3:.1f} mm,  l/D = {g.l_to_D:.2f}  ->  L_r = {d.L_r*1e3:.1f} mm",
        f"  B_δ = {g.B_delta:.2f} T,  B_ds = {g.B_ds:.2f} T,  B_sy = {g.B_sy:.2f} T",
        f"  J_Cu_s = {g.J_cu_s:.1f} A/mm²,  J_Cu_r = {g.J_cu_r:.1f} A/mm²",
        f"  N_r = {g.N_r},  q = {g.q}  ->  Q_s = {d.Q_s},  k_w1 = {d.kw1:.4f}",
        f"",
        f"Geometrija:",
        f"  τ_p     = {d.tau_p*1e3:.2f} mm",
        f"  τ_u     = {d.tau_u*1e3:.2f} mm",
        f"  δ       = {d.delta*1e3:.2f} mm   (K_c = {d.K_c:.3f})",
        f"  D_si    = {d.D_si*1e3:.1f} mm",
        f"  D_se    = {d.D_se*1e3:.1f} mm",
        f"  b_ds    = {d.b_ds*1e3:.2f} mm,   h_ds  = {d.h_ds*1e3:.2f} mm",
        f"  h_ys    = {d.h_ys*1e3:.2f} mm",
        f"  b_dr    = {d.b_dr*1e3:.2f} mm,   h_yr  = {d.h_yr*1e3:.2f} mm",
        f"",
        f"Električne količine:",
        f"  N_s = {d.N_s} ovojev/fazo,  Z_q = {d.Z_q} ovojev/utor",
        f"  I_n = {d.I_n:.2f} A   (statorski fazni RMS)",
        f"  I_m = {d.I_m:.2f} A   (vzbujalni DC),   F_m = {d.F_m:.0f} A",
        f"  J_Cu_s_actual = {d.J_cu_s_actual:.2f} A/mm²",
        f"  J_Cu_r_actual = {d.J_cu_r_actual:.2f} A/mm²",
        f"  R_s = {d.R_s_phase*1e3:.2f} mΩ/fazo",
        f"  R_r_total = {d.R_r_total*1e3:.1f} mΩ",
        f"",
        f"Izgube:",
        f"  P_Fe_zob   = {d.P_fe_tooth:7.1f} W   (masa {d.m_fe_tooth_s:.2f} kg)",
        f"  P_Fe_jarem = {d.P_fe_yoke:7.1f} W   (masa {d.m_fe_yoke_s:.2f} kg)",
        f"  P_Fe       = {d.P_fe_total:7.1f} W",
        f"  P_Cu_s     = {d.P_cu_stator:7.1f} W",
        f"  P_Cu_r     = {d.P_cu_rotor:7.1f} W",
        f"  P_Cu       = {d.P_cu_total:7.1f} W",
        f"  P_loss     = {d.P_loss_total:7.1f} W",
        f"",
        f"Volumen aktivnega dela: {d.V_active*1e6:.1f} cm³",
        f"Izkoristek η = {d.eta*100:.2f} %",
    ]
    if not d.feasible:
        lines.append("")
        lines.append("Kršitve izvedljivosti:")
        for r in d.infeasibility_reasons:
            lines.append(f"  - {r}")
    return "\n".join(lines)


# =============================================================================
# Smoke test
# =============================================================================
if __name__ == "__main__":
    machine = MachineInputs()
    material = MaterialParams()
    loss = IronLossModel.fit_default()

    # Začetna ocena geometrije (po vzoru primer1 / Pyrhönen tabel 6.2):
    initial = MotorDesignGenes(
        D_r=0.150,           # 150 mm
        l_to_D=0.8,          # L = 120 mm
        B_delta=0.90,
        B_ds=1.60,
        B_sy=1.30,
        J_cu_s=6.0,
        J_cu_r=4.0,
        N_r=40,
        q=2.0,
    )

    design = analyze(initial, machine, material, loss)
    print(pretty_print(design))
