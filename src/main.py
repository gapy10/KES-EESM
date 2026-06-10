"""
Orkestrator celotnega cevnega toka KES v7.

Tok delovanja (privzeto):
    1. Naloži vhodne podatke (CLI / YAML).
    2. Zaženi NSGA-II GA.
    3. Iz Pareto fronte izberi 5 reprezentativnih rešitev.
    4. (Faza 3) Za vsako: parametrično izrisi v FEMM in simuliraj.
    5. Izvozi rezultate v CSV in slike.

Zagon:
    python -m src.main                          # privzeti parametri
    python -m src.main --config inputs.yaml     # iz YAML
    python -m src.main --skip-femm              # samo GA + Pareto
    python -m src.main --pop 100 --gen 50

V tej fazi razvoja (po Fazi 2) je --skip-femm vedno aktiven, ker
Faza 3 (femm_model + femm_sim) še ni implementirana.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from .inputs import MachineInputs, MaterialParams, DesignBounds, load_from_yaml
from .losses import IronLossModel
from .optimization import run_nsga2, designs_from_result, save_pareto
from .pareto import select_five, summary_table
from .plot import plot_pareto


def parse_args(argv=None):
    p = argparse.ArgumentParser(prog="src.main", description="KES v7 — analitičen izračun + GA + Pareto")
    p.add_argument("--config", type=str, default=None, help="YAML konfiguracija (neobvezno)")
    p.add_argument("--pop", type=int, default=100, help="Velikost populacije GA")
    p.add_argument("--gen", type=int, default=50, help="Število generacij GA")
    p.add_argument("--seed", type=int, default=42, help="Seme za ponovljivost")
    p.add_argument("--out", type=str, default="outputs", help="Izhodni imenik")
    p.add_argument("--verbose", action="store_true", help="Verbose izpis GA")
    p.add_argument("--skip-femm", action="store_true", default=True, help="(privzeto) preskoči Fazo 3")
    return p.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)
    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    # ---- 1) Vhodni podatki ---------------------------------------------------
    if args.config:
        machine, material, bounds = load_from_yaml(args.config)
        print(f"[main] Konfiguracija naložena iz: {args.config}")
    else:
        machine, material, bounds = MachineInputs(), MaterialParams(), DesignBounds()
        print("[main] Uporabljam privzete vhodne parametre.")
    print(f"[main] P_c = {machine.P_c/1000:.1f} kW, n_c = {machine.n_c} rpm, "
          f"f_c = {machine.freq_c:.0f} Hz, M_c = {machine.torque_c:.1f} Nm")

    # ---- 2) Polinom železovih izgub ------------------------------------------
    loss = IronLossModel.fit_default()
    print(f"[main] Polinom železovih izgub: RMSE = {loss.rmse_train:.3f} W/kg")

    # ---- 3) NSGA-II ----------------------------------------------------------
    print(f"[main] NSGA-II: pop = {args.pop}, gen = {args.gen}, seed = {args.seed}")
    result, problem = run_nsga2(
        machine, material, bounds, loss,
        pop_size=args.pop, n_gen=args.gen, seed=args.seed, verbose=args.verbose,
    )

    if result.X is None or len(result.X) == 0:
        print("[main] !!! Pareto fronta prazna - vse rešitve so nedopustne.", file=sys.stderr)
        return 1

    designs = designs_from_result(result, problem)
    print(f"[main] Pareto fronta: {len(designs)} nedominiranih rešitev.")
    eta_arr = np.array([d.eta for d in designs])
    V_arr = np.array([d.V_active for d in designs])
    print(f"[main]   η razpon: {eta_arr.min()*100:.2f} % .. {eta_arr.max()*100:.2f} %")
    print(f"[main]   V razpon: {V_arr.min()*1e6:.0f} cm³ .. {V_arr.max()*1e6:.0f} cm³")

    # ---- 4) Shrani Pareto fronto --------------------------------------------
    pareto_npz = out_dir / "pareto.npz"
    save_pareto(designs, result.F, result.X, pareto_npz)
    print(f"[main] Pareto fronta shranjena: {pareto_npz}")

    # ---- 5) Izbor 5 reprezentantov ------------------------------------------
    selected, indices, labels = select_five(designs)
    print(f"[main] Izbranih {len(selected)} reprezentativnih rešitev:")
    print()
    print(summary_table(selected, labels))
    print()

    # ---- 6) Pareto graf ------------------------------------------------------
    png_path = plot_pareto(designs, indices, labels, out_path=out_dir / "pareto.png")
    print(f"[main] Pareto graf shranjen: {png_path}")

    # ---- 7) Selected5 JSON ---------------------------------------------------
    sel_path = out_dir / "selected5.json"
    sel_data = []
    for d, lab in zip(selected, labels):
        sel_data.append({
            "label": lab,
            "genes": {
                "D_r_mm": d.genes.D_r * 1e3,
                "l_to_D": d.genes.l_to_D,
                "B_delta": d.genes.B_delta,
                "B_ds": d.genes.B_ds,
                "B_sy": d.genes.B_sy,
                "J_cu_s": d.genes.J_cu_s,
                "J_cu_r": d.genes.J_cu_r,
                "N_r": d.genes.N_r,
                "q": d.genes.q,
                "delta_mm": d.genes.delta * 1e3,
            },
            "geometry_mm": {
                "L_r": d.L_r * 1e3,
                "D_si": d.D_si * 1e3,
                "D_se": d.D_se * 1e3,
                "delta": d.delta * 1e3,
                "tau_p": d.tau_p * 1e3,
                "b_ds": d.b_ds * 1e3,
                "h_ds": d.h_ds * 1e3,
                "h_ys": d.h_ys * 1e3,
                "b_dr": d.b_dr * 1e3,
                "h_yr": d.h_yr * 1e3,
            },
            "magnetic": {
                "K_cs": d.K_cs,     # Carterjev koeficient — statorska stran
                "K_cr": d.K_cr,     # Carterjev koeficient — rotorska stran
                "K_c": d.K_c,       # skupni Carterjev koeficient = K_cs · K_cr
            },
            "electrical": {
                "Q_s": d.Q_s,
                "N_s": d.N_s,
                "Z_q": d.Z_q,
                "I_n_A": d.I_n,
                "I_m_A": d.I_m,
                "F_m_A": d.F_m,
            },
            "losses_W": {
                "P_Fe_tooth": d.P_fe_tooth,
                "P_Fe_yoke": d.P_fe_yoke,
                "P_Fe_total": d.P_fe_total,
                "P_Cu_stator": d.P_cu_stator,
                "P_Cu_rotor": d.P_cu_rotor,
                "P_Cu_total": d.P_cu_total,
            },
            "eta": d.eta,
            "V_active_cm3": d.V_active * 1e6,
        })
    sel_path.write_text(json.dumps(sel_data, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"[main] selected5.json: {sel_path}")

    # ---- 8) Polni results.csv (vse Pareto rešitve) --------------------------
    rows = []
    for k, d in enumerate(designs):
        rows.append({
            "design_id": f"P{k:03d}",
            "D_r_mm": d.genes.D_r * 1e3,
            "L_r_mm": d.L_r * 1e3,
            "D_se_mm": d.D_se * 1e3,
            "q": d.genes.q,
            "Q_s": d.Q_s,
            "N_s": d.N_s,
            "N_r": d.genes.N_r,
            "I_n_A": d.I_n,
            "I_m_A": d.I_m,
            "P_Fe_W": d.P_fe_total,
            "P_Cu_W": d.P_cu_total,
            "eta_analit": d.eta,
            "V_active_cm3": d.V_active * 1e6,
            "selected_as": labels[indices.index(k)] if k in indices else "",
        })
    csv_path = out_dir / "results.csv"
    pd.DataFrame(rows).to_csv(csv_path, index=False)
    print(f"[main] results.csv: {csv_path}")

    if not args.skip_femm:
        print("[main] Faza 3 (FEMM) ni implementirana - --skip-femm je vedno aktiven.")

    return 0


if __name__ == "__main__":
    sys.exit(main())
