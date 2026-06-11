"""
Parametrični izris geometrije sinhronskega motorja v FEMM 4.2 preko pyFEMM.

Implementira točki I in J seminarske naloge:
    I: parametričen numeričen FEMM model iz analitičnih dimenzij.
    J: definicija materialov, navitja in robnih pogojev (priprava na mreženje).

Ključne lastnosti modela:
    * 2D planarno, enote v mm (FEMM standard), globina = L_r.
    * Stator: razširjena glava zoba ob zr. reži (cevelj), zaokrožena radij.
    * Rotor: sinusno aproksimirana zr. reža (zr(α) = zr_min / cos(α·κ)),
      ozka v sredini pola, široka ob robovih → sinusno inducirana napetost.
    * Vzbujalno navitje: 2 ločena bloka 'Copper' na rotorski pol (pozitiven,
      negativen smer ovojev), serijsko vezana v krogotok 'DC'.
    * Statorsko navitje: razporeditev generirana iz Q_s, p in m za izbiran q;
      podporo za q ∈ {1, 2, 3} (integralna navitja); frakcijska navitja
      (q ∈ {1.5, 2.5}) uporabijo isto razporeditev kot najbližja celostna,
      kar je za 2D simulacijo zadostna aproksimacija.
    * Robni pogoj A = 0 na zunanjem statorskem obroču.

Vir:
    Meeker, D., FEMM 4.2 User Manual.
    Meeker, D., Octave-FEMM Documentation (pyFEMM uporablja iste ALC ukaze).
"""

from __future__ import annotations

import math
import sys
import time
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence, Optional

import femm

from .analytical import MotorDesign


# =============================================================================
# Pomožne funkcije
# =============================================================================

def kill_stale_femm(verbose: bool = True) -> None:
    """Pobije morebitne osirotele femm.exe procese pred odpiranjem FEMM.

    pyfemm `openfemm()` se preko COM/ActiveX prilepi na *že odprto* FEMM
    instanco. Če je ta po prekinjenem/sesutem prejšnjem zagonu obtičala
    (odprt modalni dialog, sesutje sredi operacije), vsi `mi_*` ukazi
    spodletijo in zgradnja stroja pade pri vseh poskusih. Zato pred zagonom
    počistimo. Velja le na Windows; drugje je tih No-Op.
    """
    if not sys.platform.startswith("win"):
        return
    try:
        res = subprocess.run(
            ["taskkill", "/F", "/IM", "femm.exe"],
            capture_output=True, text=True,
        )
    except FileNotFoundError:
        return
    # taskkill vrne 0 ob uspehu, sicer (ni procesa) pustimo pri miru.
    if res.returncode == 0:
        if verbose:
            print("[info] Pobiti obstoječi femm.exe procesi pred zagonom.")
        # Kratka pavza, da OS sprosti COM registracijo in datotečne ročice.
        time.sleep(1.0)

def _mirror_rotate_point(x: float, y: float, deg: float) -> tuple[float, float]:
    """Zrcali x okrog y-osi, nato zavrti za deg stopinj (primer1 vzorec)."""
    rad = deg * math.pi / 180.0
    xm = -x
    xr = xm * math.cos(rad) - y * math.sin(rad)
    yr = xm * math.sin(rad) + y * math.cos(rad)
    return xr, yr


def generate_stator_winding_layout(Q_s: int, p: int, m: int = 3) -> list[str]:
    """Generira razporeditev faz po Q_s statorskih utorih za 3-fazni stroj.

    Vrne seznam dolžine Q_s s simboli {'A','a','B','b','C','c'}; velike
    črke pomenijo pozitivno smer, male negativno smer ovojev.

    Algoritem: standardna 60° pasovnica. Vsak utor dobi električni kot
    α_k = k · 360° · p / Q_s (mod 360°). Glede na sektor (0..6 × 60°) se
    dodeli faza:
        0..60°    -> A
        60..120°  -> -C ('c')
        120..180° -> B
        180..240° -> -A ('a')
        240..300° -> C
        300..360° -> -B ('b')
    """
    if m != 3:
        raise ValueError("Trenutno podprt le 3-fazni primer.")
    sectors = ["A", "c", "B", "a", "C", "b"]
    layout: list[str] = []
    for k in range(Q_s):
        alpha = (k * 360.0 * p / Q_s) % 360.0
        idx = int(alpha // 60.0) % 6
        layout.append(sectors[idx])
    return layout


# =============================================================================
# Parametri risanja, ki niso del MotorDesign
# =============================================================================

@dataclass
class FemmDrawOptions:
    """Konstrukcijske dimenzije, ki niso optimizacijske spremenljivke.

    Vse rotorske dimenzije so brez-dimenzionalne (delež τ_p oz. R_r), da se
    avtomatsko skalirajo na različne premere stroja. Primer1 hard-kodira
    tau_cevlja_r=90 mm in offset=16 mm za D_r=190 mm, kar daje deleže ~0.90
    (pole arc / τ_p) in ~0.17 (offset / R_r).
    """

    slot_opening_stator_mm: float = 4.0       # b_1s — odprtina statorskega utora [mm]
    slot_lining_mm: float = 1.0               # l_reza_utora — debelina izolacije [mm]
    slot_wedge_mm: float = 2.0                # zagozda — zalivanje glave utora [mm]

    # Rotor (brez-dimenzionalno):
    rotor_pole_arc_frac: float = 0.92         # širina pola / τ_p (povišano z 0.85
                                              # — pol čevelj sega bližje sosednjima
                                              # poloma)
    rotor_pole_height_frac: float = 0.17      # globina pola od površine / R_r
    sinusoidal_curve_factor: float = 85.0 / 90.0   # κ koeficient (povišano z 80/90
                                              # — sinusna kapica se na robovih
                                              # globlje potopi v rotor)
    smartmesh_level: int = 0                  # 0 = privzeti FEMM smart mesh
    precision: float = 1e-8                   # natančnost FEMM solverja (FEMM default; ne sprejme < 1e-8)
    mineangle: float = 30                     # min kot mreženja (30° = privzeto FEMM)


# =============================================================================
# Glavna funkcija
# =============================================================================

def build_motor(
    design: MotorDesign,
    out_fem_path: str | Path,
    bh_curve_path: str | Path = "data/BH_M330-35A.txt",
    options: Optional[FemmDrawOptions] = None,
    *,
    initial_angle_deg: float = 0.0,
    open_femm: bool = True,
    close_femm: bool = False,
    finalize_geometry: bool = True,
    kill_existing: bool = True,
) -> Path:
    """Zgradi FEMM model enega stroja in ga shrani v `.fem` datoteko.

    Args:
        design: rezultat `analyze()` z geometrijo v SI (m).
        out_fem_path: izhodna pot.
        bh_curve_path: pot do `BH_M330-35A.txt`.
        options: konstrukcijski parametri (b_1s, izolacija, ...).
        initial_angle_deg: zasuk rotorja po zaključenem izrisu (-35° v primer1
            ustreza nastavitvi kolesnega kota na 0°).
        open_femm: True → odpri FEMM (potrebno, če še ni odprt).
        close_femm: True → po shranjevanju zapri FEMM.
        kill_existing: True (in open_femm) → pred odprtjem pobije morebitne
            osirotele femm.exe procese, da se pyfemm ne prilepi na zataknjeno
            instanco iz prejšnjega zagona.
        finalize_geometry: True → izvede `koncaj_geometrijo` (materiali,
            krogotoki, robni pogoji). False uporabi se za testiranje izrisa.
    """
    options = options or FemmDrawOptions()
    out_fem_path = Path(out_fem_path)
    out_fem_path.parent.mkdir(parents=True, exist_ok=True)

    # ---- 1) Otvoritev FEMM in problemska definicija ------------------------
    if open_femm:
        if kill_existing:
            kill_stale_femm()
        femm.openfemm()
    femm.newdocument(0)  # 0 = magnetic problem
    # globina problema = L_r (mm), tipično 30° min angle, brez precision lim
    L_r_mm = design.L_r * 1e3
    femm.mi_probdef(
        0,                # frekvenca [Hz]; 0 = statika
        "millimeters",
        "planar",
        options.precision,
        L_r_mm,           # depth (globina paketa)
        options.mineangle,
        0,                # ac solver
    )
    femm.smartmesh(options.smartmesh_level)

    # ---- 2) Statorska geometrija --------------------------------------------
    narisi_stator(design, options)

    # ---- 3) Rotorska geometrija ---------------------------------------------
    narisi_rotor(design, options)

    # ---- 4) Materiali, krogotoki, robni pogoji ------------------------------
    if finalize_geometry:
        koncaj_geometrijo(design, options, bh_curve_path, initial_angle_deg)

    # ---- 5) Shrani --------------------------------------------------------
    femm.mi_zoomnatural()
    femm.mi_saveas(str(out_fem_path))
    if close_femm:
        femm.closefemm()
    return out_fem_path


# =============================================================================
# 2.1 STATOR
# =============================================================================

def narisi_stator(design: MotorDesign, opt: FemmDrawOptions) -> None:
    """Nariše statorski paket: zunanji obod, utore z razširjenim zobom."""
    # Vse v mm.
    D_r = design.genes.D_r * 1e3
    delta = design.delta * 1e3
    Rr = D_r / 2.0
    Rsn = Rr + delta                                # notranji radij statorja
    h_ds = design.h_ds * 1e3
    h_ys = design.h_ys * 1e3
    b_ds = design.b_ds * 1e3
    Q_s = design.Q_s
    b_1s = opt.slot_opening_stator_mm
    l_reza = opt.slot_lining_mm
    zagozda = opt.slot_wedge_mm

    # Zunanji radij statorja (skupna debelina utora + jarem):
    Rsz = Rsn + h_ys + h_ds + l_reza + zagozda

    # Korak utora ob zr. reži:
    tau_u = math.pi * (D_r + 2.0 * delta) / Q_s

    # --- 2.1.1 Zgradimo en utor pri "vrhu" (kot π/2 = 12 ura) -------------
    # Točke gradimo po vzoru primer1 KES_geometrija.py.
    # Stator nodes:
    #   0: točka na notranjem oboku, kjer se konča stranica zoba.
    #   1: zgornji rob odprtine utora.
    #   2: rob čez l_reza (izolacijo).
    #   3: spodnji rob zagozde / vstop v utor.
    nodes: list[tuple[float, float]] = []
    nodes.append((Rsn * math.cos(85 / 180 * math.pi), Rsn * math.sin(85 / 180 * math.pi)))
    op = (b_1s / 2.0, math.sqrt(max(Rsn ** 2 - (b_1s / 2.0) ** 2, 1e-9)))
    nodes.append(op)
    op2 = (op[0], op[1] + l_reza)
    nodes.append(op2)
    op3 = (op2[0], op2[1] + zagozda)

    # --- 2.1.2 Iterativno določimo širino utora znotraj zoba ------------
    # Spodnja širina utora se popravlja, dokler razdalja med utornima
    # robovoma pri preslikavi za -360/Qs ni enaka b_ds.
    utor = [op3[0], op3[1]]
    spodnja_bds = 1e9
    iters = 0
    while spodnja_bds > b_ds and iters < 5000:
        utor[0] += 0.001
        kopija = _mirror_rotate_point(utor[0], utor[1], -360.0 / Q_s)
        spodnja_bds = math.hypot(utor[0] - kopija[0], utor[1] - kopija[1])
        iters += 1

    # Zgornja širina utora (na koncu utora):
    utor_z = [utor[0], utor[1] + h_ds]
    zgornja_bds = 1e9
    iters = 0
    while zgornja_bds > b_ds and iters < 5000:
        utor_z[0] += 0.001
        kopija = _mirror_rotate_point(utor_z[0], utor_z[1], -360.0 / Q_s)
        zgornja_bds = math.hypot(utor_z[0] - kopija[0], utor_z[1] - kopija[1])
        iters += 1

    nodes.append(tuple(utor))
    nodes.append(tuple(utor_z))
    nodes.append((0.0, utor_z[1]))
    nodes.append((0.0, Rsn + l_reza))
    nodes.append((0.0, utor[1]))

    # Vrini vozlišča v FEMM:
    for x, y in nodes:
        femm.mi_addnode(x, y)

    # Loki in segmenti za en utor (vzorec primer1 vrstice 138–144):
    arc_len = 4.0 / 360.0 * 2 * math.pi * Rsn
    femm.mi_drawarc(nodes[0][0], nodes[0][1], nodes[1][0], nodes[1][1], 8, arc_len)
    femm.mi_addsegment(nodes[1][0], nodes[1][1], nodes[2][0], nodes[2][1])
    femm.mi_addsegment(nodes[3][0], nodes[3][1], nodes[2][0], nodes[2][1])
    femm.mi_addsegment(nodes[3][0], nodes[3][1], nodes[4][0], nodes[4][1])
    femm.mi_addsegment(nodes[5][0], nodes[5][1], nodes[4][0], nodes[4][1])
    femm.mi_addsegment(nodes[3][0], nodes[3][1], nodes[7][0], nodes[7][1])

    # Zaobli glavo zoba (cevelj):
    femm.mi_createradius(nodes[4][0], nodes[4][1], 1.0)

    # --- 2.1.3 Zrcali + razmnoži za vse Q_s utore ------------------------
    femm.mi_selectcircle(0, 0, Rsz, 1)
    femm.mi_mirror2(0, 0, Rsz * math.cos(math.pi / 2), Rsz * math.sin(math.pi / 2), 1)
    femm.mi_selectcircle(0, 0, Rsz, 3)
    femm.mi_mirror2(0, 0, Rsz * math.cos(math.pi / 2), Rsz * math.sin(math.pi / 2), 3)

    femm.mi_selectcircle(0, 0, Rsz, 1)
    femm.mi_copyrotate2(0, 0, 360.0 / Q_s, Q_s - 1, 1)
    femm.mi_selectcircle(0, 0, Rsz, 3)
    femm.mi_copyrotate2(0, 0, -360.0 / Q_s, Q_s - 1, 3)

    # --- 2.1.4 Zunanji statorski obroč -----------------------------------
    femm.mi_addnode(-Rsz, 0)
    femm.mi_addnode(Rsz, 0)
    femm.mi_drawarc(Rsz, 0, -Rsz, 0, 180, math.pi * Rsz)
    femm.mi_drawarc(-Rsz, 0, Rsz, 0, 180, math.pi * Rsz)


# =============================================================================
# 2.2 ROTOR
# =============================================================================

def narisi_rotor(design: MotorDesign, opt: FemmDrawOptions) -> None:
    """Nariše rotor: sinusno zaokrožen zob + jarem + utore za vzbujalno navitje."""
    D_r = design.genes.D_r * 1e3
    Rr = D_r / 2.0
    delta = design.delta * 1e3
    Rsn = Rr + delta
    p = 3                                        # 6-polni → 3 polovih para
    tau_p_mm = math.pi * D_r / (2.0 * p)         # polov korak na obodu rotorja, [mm]
    obseg_rotorja = math.pi * D_r
    bdr_mm = design.b_dr * 1e3                   # širina rotorskega zoba [mm]
    h_yr_mm = design.h_yr * 1e3                  # višina rotorskega jarma [mm]

    # Pole shoe arc dolžina v mm (delež τ_p):
    tau_cevlja_r = opt.rotor_pole_arc_frac * tau_p_mm
    # Globina pola (od površine do osnove pola) v mm:
    rotor_yoke_inner_offset = opt.rotor_pole_height_frac * Rr

    # Skupina 1 = rotor (rotira med simulacijo).
    # Sinusno zaokrožen rotorski zob (primer1 vrstice 173–193):
    tocke_a: list[tuple[float, float]] = []
    cevelj_rad = (tau_cevlja_r / 2.0 / obseg_rotorja) * 2.0 * math.pi
    krivulja_rad = opt.sinusoidal_curve_factor * math.pi / 2.0
    rad_a = cevelj_rad / (tau_cevlja_r / 2.0)
    krivulja_a_rad = krivulja_rad / (tau_cevlja_r / 2.0)

    # Generiraj točke iz sredine pola (a=0) proti robu pola (a = tau_cevlja_r/2):
    n_pts = max(int(tau_cevlja_r / 2.0) + 1, 6)
    last_R = Rsn  # placeholder; bo prepisan v zanki
    for a in range(n_pts):
        zr_a = delta / max(math.cos(a * krivulja_a_rad), 1e-6)
        R_zoba = Rsn - zr_a
        last_R = R_zoba
        x = R_zoba * math.cos(math.pi / 2.0 - a * rad_a)
        y = R_zoba * math.sin(math.pi / 2.0 - a * rad_a)
        tocke_a.append((x, y))

    # Dodatna točka navzdol (proti rotorskemu jarmu):
    x_last = tocke_a[-1][0]
    y_link = Rr - rotor_yoke_inner_offset
    link_pt = (x_last, y_link)
    tocke_a.append(link_pt)

    for x, y in tocke_a:
        femm.mi_addnode(x, y)
    for i in range(len(tocke_a) - 1):
        femm.mi_addsegment(tocke_a[i][0], tocke_a[i][1],
                           tocke_a[i + 1][0], tocke_a[i + 1][1])

    # Rotorski zobni stebriček (osnova zoba; trapez s polovično širino bdr_mm/2):
    base_x = bdr_mm / 2.0
    base_y_inner = bdr_mm * math.cos(math.pi / 6.0)  # vzorec primer1
    rotor_nodes = [
        (link_pt[0], link_pt[1]),
        (base_x, link_pt[1]),
        (base_x, base_y_inner),
    ]
    for x, y in rotor_nodes:
        femm.mi_addnode(x, y)
    femm.mi_addsegment(rotor_nodes[0][0], rotor_nodes[0][1],
                       rotor_nodes[1][0], rotor_nodes[1][1])
    femm.mi_addsegment(rotor_nodes[2][0], rotor_nodes[2][1],
                       rotor_nodes[1][0], rotor_nodes[1][1])
    femm.mi_addsegment(rotor_nodes[2][0], rotor_nodes[2][1],
                       rotor_nodes[0][0], rotor_nodes[0][1])

    # Zrcali en pol → pol postane simetričen okrog y-osi.
    femm.mi_selectcircle(0, 0, Rr + 1, 4)
    femm.mi_mirror2(0, 0, 0, Rr, 1)

    # Razmnoži za vseh 6 polov:
    femm.mi_selectcircle(0, 0, Rr + 1, 4)
    femm.mi_copyrotate2(0, 0, 60.0, 5, 1)

    # Notranji rotorski jarem (krog povezav med poli):
    R_yoke_outer = math.hypot(base_x, base_y_inner)
    R_yoke_inner = R_yoke_outer - h_yr_mm
    if R_yoke_inner < 5.0:
        R_yoke_inner = 5.0  # minimalna luknja za os
    # 6 vozlišč v sredini jarma:
    for i in range(6):
        ang = i * math.pi / 3.0 + math.pi / 6.0
        femm.mi_addnode(R_yoke_inner * math.sin(ang),
                        R_yoke_inner * math.cos(ang))


# =============================================================================
# 2.3 KONČNI KORAKI (materiali, krogotoki, robni pogoji)
# =============================================================================

def koncaj_geometrijo(
    design: MotorDesign,
    opt: FemmDrawOptions,
    bh_curve_path: str | Path,
    initial_angle_deg: float = 0.0,
) -> None:
    """Doda robni pogoj, krogotoke, materiale, navitja in nastavi začetni kot."""
    D_r = design.genes.D_r * 1e3
    delta = design.delta * 1e3
    Rr = D_r / 2.0
    Rsn = Rr + delta
    h_ds = design.h_ds * 1e3
    h_ys = design.h_ys * 1e3
    Rsz = Rsn + h_ys + h_ds + opt.slot_lining_mm + opt.slot_wedge_mm
    Q_s = design.Q_s

    # ---- Robni pogoj A = 0 ---------------------------------------------------
    femm.mi_addboundprop("A=0", 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)

    # ---- Krogotoki -----------------------------------------------------------
    femm.mi_addcircprop("A", 0, 1)   # serijsko
    femm.mi_addcircprop("B", 0, 1)
    femm.mi_addcircprop("C", 0, 1)
    femm.mi_addcircprop("DC", 0, 1)  # vzbujanje rotorja

    # ---- Materiali -----------------------------------------------------------
    # M330-35A z naloženo B-H krivuljo.
    femm.mi_addmaterial("M330-35A", 1, 1, 0, 0, 0, 0)
    with open(bh_curve_path, "r", encoding="utf-8") as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            parts = s.split()
            if len(parts) >= 2:
                B, H = float(parts[0]), float(parts[1])
                femm.mi_addbhpoint("M330-35A", B, H)
    femm.mi_getmaterial("Air")
    femm.mi_addmaterial("Copper", 1, 1, 0)

    # ---- Block labels --------------------------------------------------------
    # Statorski jarem (med h_ds in h_ys):
    R_yoke = Rsn + h_ds + opt.slot_lining_mm + opt.slot_wedge_mm + h_ys / 2.0
    femm.mi_addblocklabel(0, R_yoke)
    femm.mi_selectcircle(0, R_yoke, 5, 2)
    femm.mi_setblockprop("M330-35A", 1, 0, "<None>", 0, 0, 1)
    femm.mi_clearselected()

    # Rotor (središče):
    femm.mi_addblocklabel(0, 0)
    femm.mi_selectcircle(0, 0, 5, 2)
    femm.mi_setblockprop("M330-35A", 1, 0, "<None>", 0, 1, 1)
    femm.mi_clearselected()

    # Zrak: postavi label v MEDPOLNI ZALIV (kot 60° = med polom @90° in polom @30°),
    # kjer je zrak širok in lokacija ni občutljiva na zelo majhne δ.
    # Vsi zračni regioni (zr. reža + 6 medpolnih zalivov + 36 utornih cevljev) so
    # topološko povezani v en region, zato en label zadostuje.
    R_air = (Rr + Rsn) / 2.0
    ang_air = math.pi / 3.0       # 60° — sredina medpolnega prostora
    femm.mi_addblocklabel(R_air * math.cos(ang_air), R_air * math.sin(ang_air))
    femm.mi_selectcircle(R_air * math.cos(ang_air), R_air * math.sin(ang_air), 2, 2)
    femm.mi_setblockprop("Air", 1, 0, "<None>", 0, 0, 1)
    femm.mi_clearselected()

    # ---- Vzbujalno navitje rotorja (6 polov × 2 bloka) ----------------------
    # Postavi par blokov med dvema sosednjima polovima (kot v primer1):
    # radij na sredini med rotorskim jarmom (notranji rob pola) in spojem
    # pola s cevljem (zgornji rob).
    bdr_mm = design.b_dr * 1e3
    h_yr_mm = design.h_yr * 1e3
    R_pole_inner = math.hypot(bdr_mm / 2.0, bdr_mm * math.cos(math.pi / 6.0))
    R_pole_top = Rr - opt.rotor_pole_height_frac * Rr
    R_rotor_winding = (R_pole_inner + R_pole_top) / 2.0
    R_rotor_winding = max(R_rotor_winding, 10.0)

    for i in range(6):
        ang_pos = 3 * math.pi / 8.0 + i * 60.0 * math.pi / 180.0
        ang_neg = 5 * math.pi / 8.0 + i * 60.0 * math.pi / 180.0
        x_pos = R_rotor_winding * math.cos(ang_pos)
        y_pos = R_rotor_winding * math.sin(ang_pos)
        x_neg = R_rotor_winding * math.cos(ang_neg)
        y_neg = R_rotor_winding * math.sin(ang_neg)
        # Dva bloka: + / − ovojev, izmenično glede na pol:
        ovoji = design.genes.N_r if (i % 2 == 0) else -design.genes.N_r
        femm.mi_addblocklabel(x_pos, y_pos)
        femm.mi_selectcircle(x_pos, y_pos, 3, 2)
        femm.mi_setblockprop("Copper", 1, 0, "DC", 0, 1, ovoji)
        femm.mi_clearselected()

        femm.mi_addblocklabel(x_neg, y_neg)
        femm.mi_selectcircle(x_neg, y_neg, 3, 2)
        femm.mi_setblockprop("Copper", 1, 0, "DC", 0, 1, -ovoji)
        femm.mi_clearselected()

    # ---- Statorsko navitje --------------------------------------------------
    layout = generate_stator_winding_layout(Q_s, p=3, m=3)
    R_stator_winding = Rsn + opt.slot_lining_mm + opt.slot_wedge_mm + h_ds / 2.0
    Z_q = design.Z_q
    for i in range(Q_s):
        ang = math.pi / 2.0 + i * 2.0 * math.pi / Q_s
        x = R_stator_winding * math.cos(ang)
        y = R_stator_winding * math.sin(ang)
        sym = layout[i]
        circuit = sym.upper()
        ovoji = Z_q if sym.isupper() else -Z_q
        femm.mi_addblocklabel(x, y)
        femm.mi_selectcircle(x, y, 2, 2)
        femm.mi_setblockprop("Copper", 1, 0, circuit, 0, 0, ovoji)
        femm.mi_clearselected()

    # ---- Robni pogoj A = 0 na zunanjem statorskem oboku ---------------------
    femm.mi_selectarcsegment(0, Rsz)
    femm.mi_selectarcsegment(0, -Rsz)
    femm.mi_setarcsegmentprop(1, "A=0", 0, 0)
    femm.mi_clearselected()

    # ---- Vrtljiva skupina (rotor) -------------------------------------------
    femm.mi_selectcircle(0, 0, Rr + delta * 0.5, 1)
    femm.mi_setsegmentprop("", 0, 1, 0, 1)

    # ---- Začetni kolesni kot ------------------------------------------------
    if initial_angle_deg != 0.0:
        femm.mi_selectgroup(1)
        femm.mi_moverotate(0, 0, initial_angle_deg)
        femm.mi_clearselected()


# =============================================================================
# Smoke entrypoint
# =============================================================================
if __name__ == "__main__":
    import json
    import sys
    from .inputs import MachineInputs, MaterialParams
    from .losses import IronLossModel
    from .analytical import analyze, MotorDesignGenes, pretty_print

    # Naloži zadnjo `mid` rešitev iz selected5.json ali zaženi default:
    sel_path = Path("outputs/selected5.json")
    if sel_path.exists():
        data = json.loads(sel_path.read_text(encoding="utf-8"))
        # Vzemi 'mid' rešitev
        mid = next((d for d in data if d["label"] == "mid"), data[0])
        g = MotorDesignGenes(
            D_r=mid["genes"]["D_r_mm"] * 1e-3,
            l_to_D=mid["genes"]["l_to_D"],
            B_delta=mid["genes"]["B_delta"],
            B_ds=mid["genes"]["B_ds"],
            B_sy=mid["genes"]["B_sy"],
            J_cu_s=mid["genes"]["J_cu_s"],
            J_cu_r=mid["genes"]["J_cu_r"],
            N_r=int(mid["genes"]["N_r"]),
            q=mid["genes"]["q"],
        )
    else:
        g = MotorDesignGenes(D_r=0.150, l_to_D=0.8, B_delta=0.9,
                              B_ds=1.6, B_sy=1.3, J_cu_s=6.0,
                              J_cu_r=4.0, N_r=40, q=2.0)

    m = MachineInputs(); mat = MaterialParams(); lm = IronLossModel.fit_default()
    design = analyze(g, m, mat, lm)
    print(pretty_print(design))

    out = build_motor(design, "outputs/fem/smoke_test.fem")
    print(f"\nFEMM model shranjen: {out}")
