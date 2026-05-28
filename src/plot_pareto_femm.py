"""
Nariše Pareto fronto z označenimi 5 končnimi (FEMM-validiranimi) rešitvami.
Zamenjane rešitve (`replaced=True` v femm_pareto_results.json) so prikazane
v drugi barvi, s puščico do izvirne (padle) izbire.

Zagon:
    python -m src.plot_pareto_femm
    python -m src.plot_pareto_femm --out outputs/pareto_femm.png
"""

from __future__ import annotations

import argparse
import json
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .inputs import MachineInputs, load_from_yaml
from .plot import (
    plot_pareto,
    plot_pareto_loss,
    plot_eta_comparison,
    plot_V_comparison,
)


@dataclass
class _LiteDesign:
    """Lahka ovojnica okrog (η, V_active) za potrebe plot_pareto."""
    eta: float
    V_active: float


def _front_designs(pareto_npz_path: Path) -> list[_LiteDesign]:
    data = np.load(pareto_npz_path)
    return [_LiteDesign(eta=float(e), V_active=float(v))
            for e, v in zip(data["eta"], data["V_active"])]


def _find_original_idx(
    label: str,
    selected: list[dict],
    all_F: np.ndarray,
) -> int:
    """Vrne Pareto indeks izvirne rešitve za label (po izboru iz selected5.json)."""
    item = next(d for d in selected if d["label"] == label)
    eta = item["eta"]
    V = item["V_active_cm3"] * 1e-6
    return _nearest_pareto_idx(eta, V, all_F)


def _nearest_pareto_idx(eta: float, V_active: float, all_F: np.ndarray) -> int:
    """Najbližji indeks na Pareto fronti po (1/η−1, V_active) v normalizirani ravnini."""
    f1 = 1.0 / eta - 1.0
    Fmin = all_F.min(axis=0)
    Frng = np.maximum(all_F.max(axis=0) - Fmin, 1e-12)
    Fn = (all_F - Fmin) / Frng
    tgt = (np.array([f1, V_active]) - Fmin) / Frng
    return int(np.argmin(np.linalg.norm(Fn - tgt, axis=1)))


def main(argv=None) -> int:
    p = argparse.ArgumentParser(prog="src.plot_pareto_femm")
    p.add_argument("--pareto", default="outputs/pareto.npz")
    p.add_argument("--selected", default="outputs/selected5.json")
    p.add_argument("--femm", default="outputs/femm_pareto_results.json")
    p.add_argument("--out", default="outputs/pareto_femm.png",
                   help="Pareto graf z η na y-osi")
    p.add_argument("--out-loss", default="outputs/pareto_femm_loss.png",
                   help="Pareto graf s P_loss na y-osi (minimizacija)")
    p.add_argument("--out-eta", default="outputs/eta_comparison.png",
                   help="Primerjalni graf η: analitično vs FEMM")
    p.add_argument("--out-V", default="outputs/V_comparison.png",
                   help="Primerjalni graf volumna aktivnega dela")
    p.add_argument("--config", default="inputs.yaml",
                   help="YAML s P_c (za pretvorbo η → P_loss)")
    p.add_argument("--title", default="Pareto fronta z 5 FEMM-validiranimi rešitvami")
    args = p.parse_args(argv)

    pareto_path = Path(args.pareto)
    sel_path = Path(args.selected)
    femm_path = Path(args.femm)
    if not (pareto_path.exists() and sel_path.exists() and femm_path.exists()):
        print("!! Manjkajoče vhodne datoteke. Preverite outputs/.", file=sys.stderr)
        return 1

    pareto_data = np.load(pareto_path)
    all_F = pareto_data["F"]
    designs = _front_designs(pareto_path)
    selected = json.loads(sel_path.read_text(encoding="utf-8"))
    femm_results = json.loads(femm_path.read_text(encoding="utf-8"))

    # Uredi po standardnem zaporedju oznak:
    label_order = ["max_eta", "high_eta_mid", "mid", "low_eta_mid", "min_V"]
    femm_by_label = {d["label"]: d for d in femm_results}

    selected_indices: list[int] = []
    replaced_flags: list[bool] = []
    original_indices: list[int] = []
    used_labels: list[str] = []
    for lab in label_order:
        if lab not in femm_by_label:
            continue
        fr = femm_by_label[lab]
        # pareto_idx je bil dodan šele v novejši verziji koda; če manjka,
        # ga rekonstruiramo iz analitičnih (η, V) iz JSON-a.
        if "pareto_idx" in fr:
            sel_idx = int(fr["pareto_idx"])
        else:
            eta_femm_design = fr["design"]["eta_analit"]
            V_femm_design = fr["design"]["V_active_cm3"] * 1e-6
            sel_idx = _nearest_pareto_idx(eta_femm_design, V_femm_design, all_F)
        selected_indices.append(sel_idx)
        replaced_flags.append(bool(fr.get("replaced", False)))
        original_indices.append(_find_original_idx(lab, selected, all_F))
        used_labels.append(lab)

    # Klasični Pareto graf (η na y-osi). Označimo samo vseh 5 kandidatov,
    # ki so uspešno prešli FEMM (vključno z morebitno zamenjavo); izvirne
    # padle rešitve ne risemo več, da graf ostane berljiv.
    out = plot_pareto(
        designs,
        selected_indices=selected_indices,
        selected_labels=used_labels,
        out_path=args.out,
        title=args.title,
    )
    print(f"Pareto (η) graf shranjen: {out}")

    # Pareto z izgubami na y-osi (minimizacijska perspektiva):
    cfg_path = Path(args.config)
    if cfg_path.exists():
        machine, _, _ = load_from_yaml(cfg_path)
    else:
        machine = MachineInputs()
    out_loss = plot_pareto_loss(
        designs,
        P_out_W=machine.P_c,
        selected_indices=selected_indices,
        selected_labels=used_labels,
        out_path=args.out_loss,
        title="Pareto fronta z 5 FEMM-validiranimi rešitvami (izgube)",
    )
    print(f"Pareto (izgube) graf shranjen: {out_loss}")

    # Primerjalni graf η: analitično vs FEMM:
    eta_analyt = [femm_by_label[lab]["design"]["eta_analit"] for lab in used_labels]
    eta_femm = [femm_by_label[lab]["femm_losses"]["eta_femm"] for lab in used_labels]
    out_eta = plot_eta_comparison(used_labels, eta_analyt, eta_femm,
                                  out_path=args.out_eta)
    print(f"η primerjalni graf shranjen: {out_eta}")

    # Primerjalni graf V:
    V_vals = [femm_by_label[lab]["design"]["V_active_cm3"] for lab in used_labels]
    out_V = plot_V_comparison(used_labels, V_vals, out_path=args.out_V)
    print(f"V primerjalni graf shranjen: {out_V}")

    print()
    for lab, sidx, oidx, repl in zip(used_labels, selected_indices,
                                     original_indices, replaced_flags):
        if repl:
            print(f"  [{lab}] zamenjan: izvirni Pareto idx {oidx} → "
                  f"uporabljen idx {sidx}")
        else:
            print(f"  [{lab}] Pareto idx {sidx}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
