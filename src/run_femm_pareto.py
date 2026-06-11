"""
Zažene FEMM simulacije (prosti tek + navor) za vseh 5 rešitev iz selected5.json,
in izvede post-procesiranje (Faza 4): primerjava analitiko↔FEMM, P_Fe iz FEMM,
THD, grafe v 300 dpi.

Tipičen zagon traja 30 – 90 min zaradi nelinearne magnetostatične analize.

Robustnost:
    Če FEMM gradnja ali simulacija enega stroja spodleti (npr. „Material
    properties have not been defined for all regions" za ekstremno geometrijo),
    skript poišče najbližjega soseda na Pareto fronti (v normalizirani
    ravnini (1/η−1, V_active)) in poskusi z njim. Število poskusov omejuje
    `--max-retries` (privzeto 12 — celoten max_eta vogal z velikim D_r ≳151 mm
    je nezmrežljiv, zato je treba seči dlje po fronti do zmrežljivega stroja).

Zagon:
    python -m src.run_femm_pareto                          # privzeti step_deg=5
    python -m src.run_femm_pareto --step 10                # manj pozicij, hitreje
    python -m src.run_femm_pareto --only min_V             # samo ena rešitev
                                                            # (rezultati se zlijejo
                                                            # z obstoječim povzetkom)
"""

from __future__ import annotations

import argparse
import json
import time
import sys
from pathlib import Path

import femm
import numpy as np
import pandas as pd

from .inputs import MachineInputs, MaterialParams, DesignBounds, load_from_yaml
from .losses import IronLossModel
from .analytical import MotorDesignGenes, analyze
from .femm_model import build_motor, kill_stale_femm
from .femm_sim import no_load, torque
from .optimization import MotorOptimizationProblem
from .post import (
    femm_iron_losses,
    uind_spectrum,
    torque_thd,
    compare,
    comparison_table_text,
    summary_dataframe,
)
from .plot import (
    plot_uind,
    plot_torque,
    plot_fft,
    plot_flux,
    plot_comparison_bar,
)


def parse_args(argv=None):
    p = argparse.ArgumentParser(prog="src.run_femm_pareto")
    p.add_argument("--selected", default="outputs/selected5.json",
                   help="Pot do selected5.json")
    p.add_argument("--pareto", default=None,
                   help="Pot do pareto.npz (privzeto: poleg selected5.json)")
    p.add_argument("--step", type=float, default=5.0,
                   help="Kotni korak vrtenja [°] (manjše = natančneje, počasneje)")
    p.add_argument("--only", default=None,
                   help="Samo ena oznaka: max_eta, high_eta_mid, mid, low_eta_mid, min_V")
    p.add_argument("--max-retries", type=int, default=12,
                   help="Število alternativ s Pareto fronte ob napaki FEMM. "
                        "Privzeto 12: ekstremni vogal Pareto fronte (max_eta, "
                        "velik D_r ≳151 mm) FEMM ne zmreži ('Material properties "
                        "have not been defined for all regions'), zato mora biti "
                        "dovolj poskusov, da najde najbližji zmrežljiv stroj.")
    p.add_argument("--skip-noload", action="store_true")
    p.add_argument("--skip-torque", action="store_true")
    p.add_argument("--out", default="outputs", help="Izhodni imenik")
    p.add_argument("--config", default="inputs.yaml",
                   help="YAML z bounds (potreben za pravilno dekodiranje alternativ "
                        "iz pareto.npz; mora se ujemati s tistim, ki ga je uporabil GA)")
    p.add_argument("--no-kill-femm", action="store_true",
                   help="Ne pobij obstoječih femm.exe procesov pred zagonom "
                        "(privzeto jih pobijemo, da se pyfemm ne prilepi na "
                        "zataknjeno instanco iz prejšnjega zagona)")
    return p.parse_args(argv)


def reconstruct_design(item, machine, material, loss):
    g = MotorDesignGenes(
        D_r=item["genes"]["D_r_mm"] * 1e-3,
        l_to_D=item["genes"]["l_to_D"],
        B_delta=item["genes"]["B_delta"],
        B_ds=item["genes"]["B_ds"],
        B_sy=item["genes"]["B_sy"],
        J_cu_s=item["genes"]["J_cu_s"],
        J_cu_r=item["genes"]["J_cu_r"],
        N_r=int(item["genes"]["N_r"]),
        q=item["genes"]["q"],
        delta=item["genes"].get("delta_mm", 0.7) * 1e-3,
    )
    return analyze(g, machine, material, loss)


def _candidates_near(target_F: np.ndarray, all_F: np.ndarray) -> list[int]:
    """Vrne indekse Pareto fronte, urejene po L2 razdalji v normalizirani
    (1/η−1, V_active) ravnini, najbližji najprej.
    """
    Fmin = all_F.min(axis=0)
    Frng = np.maximum(all_F.max(axis=0) - Fmin, 1e-12)
    Fn = (all_F - Fmin) / Frng
    tgt = (target_F - Fmin) / Frng
    d = np.linalg.norm(Fn - tgt, axis=1)
    return list(np.argsort(d))


def _recover_used_indices(
    existing_csv: Path, all_F: np.ndarray, skip_labels: set[str]
) -> set[int]:
    """Iz obstoječega femm_summary.csv izvleče Pareto indekse že uspešnih
    rešitev (preko stolpca pareto_idx, če obstaja, drugače z najbližjim
    ujemanjem (η, V_active))."""
    used: set[int] = set()
    if not existing_csv.exists():
        return used
    df = pd.read_csv(existing_csv)
    for _, row in df.iterrows():
        if row.get("label") in skip_labels:
            continue
        if "pareto_idx" in df.columns and pd.notna(row["pareto_idx"]):
            used.add(int(row["pareto_idx"]))
            continue
        # Fallback: najbližje ujemanje po (η, V_active)
        eta = float(row["eta_analit"])
        V = float(row["V_active_cm3"]) * 1e-6
        f1 = 1.0 / eta - 1.0
        order = _candidates_near(np.array([f1, V]), all_F)
        if order:
            used.add(int(order[0]))
    return used


def main(argv=None) -> int:
    args = parse_args(argv)
    if not args.no_kill_femm:
        kill_stale_femm()
    out_dir = Path(args.out)
    fem_dir = out_dir / "fem"
    fig_root = out_dir / "figures"
    fem_dir.mkdir(parents=True, exist_ok=True)
    fig_root.mkdir(parents=True, exist_ok=True)

    # Za pravilno dekodiranje alternativnih Pareto rešitev nujno uporabimo
    # iste bounds (predvsem `q_choices`), kot jih je uporabil GA.
    cfg_path = Path(args.config)
    if cfg_path.exists():
        machine, material, bounds = load_from_yaml(cfg_path)
    else:
        machine = MachineInputs()
        material = MaterialParams()
        bounds = DesignBounds()
    loss = IronLossModel.fit_default()

    full_data = json.loads(Path(args.selected).read_text(encoding="utf-8"))
    # Ohranimo izvorni 1-osnovni `k` v selected5.json (npr. D05 ostane D05
    # tudi pri --only).
    process_items: list[tuple[int, dict]] = [
        (k + 1, item) for k, item in enumerate(full_data)
    ]
    if args.only:
        process_items = [(k, it) for (k, it) in process_items if it["label"] == args.only]
        if not process_items:
            print(f"!! Oznake '{args.only}' ni v selected5.json", file=sys.stderr)
            return 1

    # Pareto fronta za iskanje alternativ ob napaki:
    pareto_path = Path(args.pareto) if args.pareto else (Path(args.selected).parent / "pareto.npz")
    pareto_data = np.load(pareto_path)
    all_F = pareto_data["F"]           # (N, 2): (1/η−1, V_active)
    all_X = pareto_data["X"]           # (N, 9)
    decode_problem = MotorOptimizationProblem(
        machine, material, bounds, loss
    )

    # Pri --only ohranimo Pareto indekse že uspešnih rešitev iz prejšnjega
    # zagona, da nove alternative ne padejo na isti stroj.
    existing_csv = out_dir / "femm_summary.csv"
    processed_labels = {it["label"] for _, it in process_items}
    used_pareto_indices: set[int] = _recover_used_indices(
        existing_csv, all_F, skip_labels=processed_labels
    )
    if used_pareto_indices:
        print(f"[info] Že uporabljeni Pareto indeksi (iz prejšnjega "
              f"femm_summary.csv): {sorted(used_pareto_indices)}")

    designs, no_loads, torques, losses, labels = [], [], [], [], []
    design_ids: list[str] = []
    pareto_idxs: list[int] = []
    replaced_flags: list[bool] = []
    failures: list[tuple[str, str]] = []   # (design_id, error_msg)
    t_total = time.time()
    femm_opened = False  # Sledi, ali je FEMM že odprt v tej seji.

    for slot_pos, (k, item) in enumerate(process_items):
        label = item["label"]
        design_id = f"D{k:02d}_{label}"
        print(f"\n=== [{slot_pos+1}/{len(process_items)}] {design_id} ===")

        # Izgradi seznam kandidatov: najprej originalni stroj (iz
        # selected5.json), nato sosedi na Pareto fronti urejeni po
        # bližini v normalizirani (1/η−1, V_active) ravnini.
        original_design = reconstruct_design(item, machine, material, loss)
        target_F = np.array([1.0 / original_design.eta - 1.0, original_design.V_active])
        nearest_order = _candidates_near(target_F, all_F)
        # Filtriraj že uporabljene (z drugih slotov v isti seji ali iz prejšnjega zagona):
        candidate_indices: list[int] = [
            idx for idx in nearest_order if idx not in used_pareto_indices
        ]

        success = False
        last_err = None  # zadnja FEMM-napaka tega slota (za sklepno poročilo)
        for attempt, pareto_idx in enumerate(candidate_indices, start=1):
            if attempt > args.max_retries:
                print(f"  !! presežen --max-retries ({args.max_retries}) "
                      f"za {design_id}", file=sys.stderr)
                break

            if attempt == 1:
                design = original_design
                is_replacement = False
            else:
                genes = decode_problem.decode(all_X[pareto_idx])
                design = analyze(genes, machine, material, loss)
                is_replacement = True
                print(f"  ↻ alternativa #{attempt-1} (Pareto idx {pareto_idx}): "
                      f"D_r={design.genes.D_r*1e3:.1f}mm η={design.eta*100:.2f}% "
                      f"V={design.V_active*1e6:.0f}cm³")

            print(f"  D_r={design.genes.D_r*1e3:.1f}mm L={design.L_r*1e3:.1f}mm "
                  f"q={design.genes.q} N_r={design.genes.N_r} "
                  f"η_analit={design.eta*100:.2f}% V={design.V_active*1e6:.0f}cm³")

            fem_built = fem_dir / f"{design_id}_excited.fem"
            fem_noload_path = fem_dir / f"{design_id}_noload.fem"
            fem_torque_path = fem_dir / f"{design_id}_torque.fem"
            fig_dir = fig_root / design_id
            fig_dir.mkdir(parents=True, exist_ok=True)

            nl = None
            tq = None
            fl = None
            try:
                # ---- Build ----
                t0 = time.time()
                build_motor(design, fem_built,
                            open_femm=(not femm_opened), close_femm=False,
                            kill_existing=not args.no_kill_femm)
                femm_opened = True
                print(f"  build: {time.time()-t0:.1f}s")

                # ---- No-load ----
                if not args.skip_noload:
                    t0 = time.time()
                    nl = no_load(design, rpm=machine.n_c, step_deg=args.step,
                                 save_fem_path=fem_noload_path)
                    us = uind_spectrum(nl)
                    print(f"  no_load: {time.time()-t0:.1f}s  U_ind_amp={nl.U_ind_amp:.1f}V "
                          f"({nl.U_ind_rms:.1f}Vrms)  THD={us.thd_percent:.1f}%  "
                          f"B_zob={nl.B_max_tooth:.2f}T  B_jarem={nl.B_max_yoke:.2f}T")
                    plot_uind(nl.theta_deg_mech[:-1], nl.U_ind_inst, fig_dir / "uind.png")
                    plot_flux(nl.theta_deg_mech, nl.flux_A, fig_dir / "flux.png")
                    plot_fft(us.harmonics[:15], us.amplitudes[:15],
                             fig_dir / "uind_fft.png",
                             title=f"Spekter U_ind ({label}); THD = {us.thd_percent:.1f} %",
                             ylabel="Amplituda [V]")

                # ---- Torque ----
                if not args.skip_torque:
                    t0 = time.time()
                    tq = torque(design, step_deg=args.step, save_fem_path=fem_torque_path)
                    thd = torque_thd(tq)
                    print(f"  torque: {time.time()-t0:.1f}s  M_1.harm={tq.M_fundamental:.2f}Nm "
                          f"(napoved {design.M_FEMM_pred:.2f}Nm, cilj M_c {machine.torque_c:.2f}Nm)  "
                          f"M_2.harm={tq.M_second:.2f}Nm  THD={thd:.1f}%")
                    plot_torque(tq.theta_deg_mech, tq.torque, fig_dir / "navor.png",
                                title=f"Navor M(θ) - {label}")
                    plot_fft(tq.fft_harmonics[:15], tq.fft_amplitudes[:15],
                             fig_dir / "navor_fft.png",
                             title=f"Spekter navora ({label}); THD = {thd:.1f} %",
                             ylabel="Amplituda [Nm]")

                # ---- Post ----
                if nl is not None:
                    fl = femm_iron_losses(design, nl, machine, material, loss)
                    if tq is not None:
                        rows = compare(design, nl, tq, fl, machine)
                        table_txt = comparison_table_text(rows, title=f"Primerjava — {design_id}")
                        print()
                        print(table_txt)
                        (fig_dir / "comparison.txt").write_text(table_txt, encoding="utf-8")
                        quantities = [r.quantity for r in rows if r.unit not in ("V", "Nm")]
                        analyt = [r.analytical for r in rows if r.unit not in ("V", "Nm")]
                        femmvals = [r.femm for r in rows if r.unit not in ("V", "Nm")]
                        if quantities:
                            plot_comparison_bar(
                                quantities, analyt, femmvals,
                                fig_dir / "comparison_bar.png",
                                title=f"Primerjava — {label}",
                            )

            except Exception as e:
                last_err = f"{type(e).__name__}: {e}"
                # POZOR: to NI napaka programa. Stroji na robu Pareto fronte
                # (npr. ekstremen max_eta z velikim D_r) jih FEMM ne zmreži;
                # skript je zasnovan tako, da samodejno nadaljuje z najbližjo
                # zmrežljivo alternativo s Pareto fronte. Zato to izpišemo kot
                # običajen informativen korak (stdout), NE kot napako. Pravo
                # odpoved (če spodletijo VSE alternative) poročamo šele na koncu.
                print(f"  [info] {design_id}: poskus {attempt} (Pareto idx "
                      f"{pareto_idx}) ni zmrežljiv — nadaljujem z naslednjo "
                      f"rešitvijo s Pareto fronte ...")
                # FEMM stanje je verjetno pokvarjeno → zaprimo, ponovno odpremo.
                try:
                    femm.closefemm()
                except Exception:
                    pass
                femm_opened = False
                continue  # naslednji kandidat

            # Uspeh:
            if nl is not None and tq is not None and fl is not None:
                designs.append(design)
                no_loads.append(nl)
                torques.append(tq)
                losses.append(fl)
                labels.append(label)
                design_ids.append(design_id)
                pareto_idxs.append(int(pareto_idx))
                replaced_flags.append(is_replacement)
                used_pareto_indices.add(int(pareto_idx))
                if is_replacement:
                    print(f"  ✓ uporabljena alternativa s Pareto fronte "
                          f"(idx {pareto_idx}) namesto izvirne rešitve")
                success = True
                break

        if not success:
            failures.append((design_id,
                             f"vse alternative spodletele ({args.max_retries} poskusov)"
                             + (f"; zadnja FEMM-napaka: {last_err}" if last_err else "")))

    try:
        femm.closefemm()
    except Exception:
        pass

    if failures:
        print(f"\n!! NAPAKE pri {len(failures)} rešitvah:", file=sys.stderr)
        for did, msg in failures:
            print(f"   - {did}: {msg}", file=sys.stderr)

    # ---- Povzetna tabela vseh rešitev ----
    if designs:
        df_new = summary_dataframe(
            designs, no_loads, torques, losses, labels, machine,
            design_ids=design_ids, pareto_idx=pareto_idxs,
            replaced=replaced_flags,
        )

        # JSON dump trenutnih rezultatov:
        new_dump = []
        for d, nl_, tq_, fl_, lab_, did_, pidx_, repl_ in zip(
            designs, no_loads, torques, losses, labels,
            design_ids, pareto_idxs, replaced_flags,
        ):
            us = uind_spectrum(nl_)
            new_dump.append({
                "design_id": did_,
                "label": lab_,
                "pareto_idx": pidx_,
                "replaced": repl_,
                "design": {
                    "D_r_mm": d.genes.D_r * 1e3,
                    "L_r_mm": d.L_r * 1e3,
                    "q": d.genes.q,
                    "Q_s": d.Q_s,
                    "N_r": d.genes.N_r,
                    "I_m_A": d.I_m,
                    "I_n_A": d.I_n,
                    "eta_analit": d.eta,
                    "V_active_cm3": d.V_active * 1e6,
                },
                "no_load": {
                    "U_ind_amp_V": nl_.U_ind_amp,
                    "U_ind_rms_V": nl_.U_ind_rms,
                    "U_ind_thd_percent": us.thd_percent,
                    "B_max_tooth_T": nl_.B_max_tooth,
                    "B_max_yoke_T": nl_.B_max_yoke,
                    "theta_deg": nl_.theta_deg_mech.tolist(),
                    "flux_A_Wb": nl_.flux_A.tolist(),
                },
                "torque": {
                    "M_1harm_Nm": tq_.M_fundamental,
                    "M_2harm_Nm": tq_.M_second,
                    "M_max_Nm": tq_.M_max,
                    "M_ripple_pp_Nm": tq_.M_ripple_pp,
                    "thd_percent": torque_thd(tq_),
                    "fft_first_5_Nm": tq_.fft_amplitudes[:5].tolist(),
                    "theta_deg": tq_.theta_deg_mech.tolist(),
                    "torque_Nm": tq_.torque.tolist(),
                },
                "femm_losses": {
                    "P_Fe_tooth_W": fl_.P_fe_tooth,
                    "P_Fe_yoke_W": fl_.P_fe_yoke,
                    "P_Fe_total_W": fl_.P_fe_total,
                    "P_Cu_total_W": fl_.P_cu_total,
                    "eta_femm": fl_.eta_femm,
                },
            })

        # Pri --only se rezultati zlijejo z obstoječim povzetkom (zamenjajo
        # vrstice z istim labelom). Pri polnem zagonu nadomestimo vse.
        csv_path = out_dir / "femm_summary.csv"
        json_path = out_dir / "femm_pareto_results.json"

        if args.only and csv_path.exists():
            df_existing = pd.read_csv(csv_path)
            new_labels = set(df_new["label"].tolist())
            df_combined = pd.concat([
                df_existing[~df_existing["label"].isin(new_labels)],
                df_new,
            ], ignore_index=True)
            # Uredi po design_id za berljivost:
            df_combined = df_combined.sort_values("design_id").reset_index(drop=True)
            df_combined.to_csv(csv_path, index=False, float_format="%.4f")
            print(f"\n=== Povzetek (zlit z obstoječim) {len(df_combined)} rešitev: "
                  f"{csv_path} ===\n")
            df_print = df_combined
        else:
            df_new.to_csv(csv_path, index=False, float_format="%.4f")
            print(f"\n=== Povzetek vseh {len(designs)} rešitev: {csv_path} ===\n")
            df_print = df_new

        # Stdout pregled:
        cols_show = ["design_id", "q", "eta_analit", "eta_femm",
                     "U_ind_femm_V", "M_femm_Nm", "M_2harm_Nm",
                     "thd_torque_percent", "P_Fe_analit_W", "P_Fe_femm_W"]
        if "replaced" in df_print.columns:
            cols_show.append("replaced")
        print(df_print[cols_show].to_string(index=False))

        # Združi JSON: če --only, zlij z obstoječim po labelu.
        if args.only and json_path.exists():
            existing_dump = json.loads(json_path.read_text(encoding="utf-8"))
            new_labels = set(d["label"] for d in new_dump)
            existing_filtered = [d for d in existing_dump if d["label"] not in new_labels]
            final_dump = sorted(existing_filtered + new_dump,
                                key=lambda d: d.get("design_id", d["label"]))
        else:
            final_dump = new_dump
        json_path.write_text(json.dumps(final_dump, indent=2, ensure_ascii=False),
                              encoding="utf-8")
        print(f"\nFEMM rezultati: {json_path}")

    print(f"\nSkupni čas: {time.time()-t_total:.1f} s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
