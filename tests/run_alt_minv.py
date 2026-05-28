"""Zažene FEMM simulacijo za alternativno min_V rešitev iz Pareto fronte.

Originalna D05 (D_r=135.6 mm) pade na geometrijski edge case. Izberemo
najmanjšo rešitev z D_r >= 145 mm kot zamenjavo.
"""

import json
import time
from pathlib import Path
import sys

import numpy as np
import femm

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from src.inputs import MachineInputs, MaterialParams, DesignBounds, q_to_qs_and_kw1
from src.losses import IronLossModel
from src.analytical import MotorDesignGenes, analyze, pretty_print
from src.femm_model import build_motor
from src.femm_sim import no_load, torque
from src.post import (
    femm_iron_losses, uind_spectrum, torque_thd, compare,
    comparison_table_text,
)
from src.plot import plot_uind, plot_torque, plot_fft, plot_flux, plot_comparison_bar


def main():
    machine = MachineInputs()
    material = MaterialParams()
    bounds = DesignBounds()
    loss = IronLossModel.fit_default()

    # Naloži Pareto fronto in poišči najboljši kandidat:
    data = np.load("outputs/pareto.npz")
    X = data["X"]
    V = data["V_active"] * 1e6   # cm3
    D_r_mm = X[:, 0] * 1e3

    # Najmanjši V z D_r >= 145 mm
    mask = D_r_mm >= 145.0
    candidates = np.where(mask)[0]
    k = candidates[np.argmin(V[mask])]
    x = X[k]
    q_idx = int(round(np.clip(x[8], 0, len(bounds.q_choices) - 1)))

    genes = MotorDesignGenes(
        D_r=float(x[0]),
        l_to_D=float(x[1]),
        B_delta=float(x[2]),
        B_ds=float(x[3]),
        B_sy=float(x[4]),
        J_cu_s=float(x[5]),
        J_cu_r=float(x[6]),
        N_r=int(round(x[7])),
        q=float(bounds.q_choices[q_idx]),
    )

    design = analyze(genes, machine, material, loss)
    print(pretty_print(design))
    print()

    label = "alt_min_V"
    design_id = f"D05_{label}"
    fem_dir = Path("outputs/fem")
    fem_dir.mkdir(exist_ok=True, parents=True)
    fig_dir = Path(f"outputs/figures/{design_id}")
    fig_dir.mkdir(exist_ok=True, parents=True)

    fem_built = fem_dir / f"{design_id}_excited.fem"
    fem_noload_path = fem_dir / f"{design_id}_noload.fem"
    fem_torque_path = fem_dir / f"{design_id}_torque.fem"

    # ---- Build ----
    t0 = time.time()
    build_motor(design, fem_built, open_femm=True, close_femm=False)
    print(f"[1/3] Build: {time.time()-t0:.1f}s")

    # ---- No-load (step=5 za hitrost) ----
    t0 = time.time()
    nl = no_load(design, rpm=machine.n_c, step_deg=5.0, save_fem_path=fem_noload_path)
    us = uind_spectrum(nl)
    print(f"[2/3] No-load: {time.time()-t0:.1f}s  U_ind={nl.U_ind_amp:.1f}V ({nl.U_ind_rms:.1f}Vrms) "
          f"THD={us.thd_percent:.1f}%  B_zob={nl.B_max_tooth:.2f}T B_jarem={nl.B_max_yoke:.2f}T")
    plot_uind(nl.theta_deg_mech[:-1], nl.U_ind_inst, fig_dir / "uind.png")
    plot_flux(nl.theta_deg_mech, nl.flux_A, fig_dir / "flux.png")
    plot_fft(us.harmonics[:15], us.amplitudes[:15], fig_dir / "uind_fft.png",
             title=f"Spekter U_ind ({label}); THD = {us.thd_percent:.1f}%",
             ylabel="Amplituda [V]")

    # ---- Torque ----
    t0 = time.time()
    tq = torque(design, step_deg=5.0, save_fem_path=fem_torque_path)
    thd = torque_thd(tq)
    print(f"[3/3] Torque: {time.time()-t0:.1f}s  M_1.harm={tq.M_fundamental:.2f}Nm "
          f"(cilj {machine.torque_c:.2f}Nm)  M_2.harm={tq.M_second:.2f}Nm  THD={thd:.1f}%")
    plot_torque(tq.theta_deg_mech, tq.torque, fig_dir / "navor.png",
                title=f"Navor M(θ) - {label}")
    plot_fft(tq.fft_harmonics[:15], tq.fft_amplitudes[:15], fig_dir / "navor_fft.png",
             title=f"Spekter navora ({label}); THD = {thd:.1f}%",
             ylabel="Amplituda [Nm]")

    # ---- Post ----
    fl = femm_iron_losses(design, nl, machine, material, loss)
    rows = compare(design, nl, tq, fl, machine)
    table_txt = comparison_table_text(rows, title=f"Primerjava — {design_id}")
    print()
    print(table_txt)
    (fig_dir / "comparison.txt").write_text(table_txt, encoding="utf-8")
    quantities = [r.quantity for r in rows if r.unit not in ("V", "Nm")]
    analyt = [r.analytical for r in rows if r.unit not in ("V", "Nm")]
    femmvals = [r.femm for r in rows if r.unit not in ("V", "Nm")]
    if quantities:
        plot_comparison_bar(quantities, analyt, femmvals, fig_dir / "comparison_bar.png",
                            title=f"Primerjava — {label}")

    femm.closefemm()

    # JSON izpis:
    out = {
        "label": label,
        "design_id": design_id,
        "design": {
            "D_r_mm": design.genes.D_r * 1e3,
            "L_r_mm": design.L_r * 1e3,
            "q": design.genes.q,
            "N_r": design.genes.N_r,
            "I_m_A": design.I_m,
            "I_n_A": design.I_n,
            "eta_analit": design.eta,
            "V_active_cm3": design.V_active * 1e6,
        },
        "no_load": {
            "U_ind_amp_V": nl.U_ind_amp,
            "U_ind_rms_V": nl.U_ind_rms,
            "U_ind_thd_percent": us.thd_percent,
            "B_max_tooth_T": nl.B_max_tooth,
            "B_max_yoke_T": nl.B_max_yoke,
        },
        "torque": {
            "M_1harm_Nm": tq.M_fundamental,
            "M_2harm_Nm": tq.M_second,
            "M_max_Nm": tq.M_max,
            "M_ripple_pp_Nm": tq.M_ripple_pp,
            "thd_percent": thd,
        },
        "femm_losses": {
            "P_Fe_total_W": fl.P_fe_total,
            "P_Cu_total_W": fl.P_cu_total,
            "eta_femm": fl.eta_femm,
        },
    }
    out_path = Path(f"outputs/{design_id}_result.json")
    out_path.write_text(json.dumps(out, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\nRezultat: {out_path}")


if __name__ == "__main__":
    main()
