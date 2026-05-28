"""
Izbor reprezentativnih rešitev iz Pareto fronte.

Implementira točko K seminarske naloge: iz Pareto fronte določi 5 rešitev,
ki bodo simulirane v FEMM (Faza 3).

Strategija izbora:
    1. točka A: maksimalen η (najnižja vrednost f1 = 1/η - 1),
    2. točka E: minimalen V_active (najnižja vrednost f2),
    3. točke B, C, D: enakomerno razporejene po loku fronte
       (po normalizirani razdalji vzdolž krivulje v ravnini f1-f2).

Vir:
    Deb, K., "Multi-Objective Optimization using Evolutionary Algorithms",
    Wiley, 2001, poglavje 8 (visualization & decision making).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np

from .analytical import MotorDesign


# =============================================================================
# Helpers
# =============================================================================

def _sort_by_eta(designs: Sequence[MotorDesign]) -> list[int]:
    """Vrne indekse rešitev, urejene po naraščajočem η (od najmanjšega)."""
    return sorted(range(len(designs)), key=lambda i: designs[i].eta)


def _normalize_front(F: np.ndarray) -> np.ndarray:
    """Normalizira (N, 2) ciljno matriko v [0, 1]^2."""
    F = np.asarray(F, dtype=float)
    rng = F.max(axis=0) - F.min(axis=0)
    rng = np.where(rng < 1e-12, 1.0, rng)
    return (F - F.min(axis=0)) / rng


# =============================================================================
# Selekcija 5 rešitev
# =============================================================================

def select_five(
    designs: Sequence[MotorDesign],
) -> tuple[list[MotorDesign], list[int], list[str]]:
    """Iz Pareto fronte izbere 5 reprezentativnih rešitev.

    Args:
        designs: seznam vseh nedominiranih MotorDesign rešitev.

    Returns:
        (selected_designs, indices, labels) — labels so opisi:
        ['min_V', 'low_eta_mid', 'mid', 'high_eta_mid', 'max_eta'].
    """
    n = len(designs)
    if n == 0:
        return [], [], []
    if n <= 5:
        # Premalo rešitev — vrnemo vse + ponovimo zadnjo.
        labels = ["pt_" + str(i + 1) for i in range(n)]
        return list(designs), list(range(n)), labels

    # Sestavi ciljno matriko F (f1 = 1/η - 1, f2 = V_active).
    F = np.array([[1.0 / d.eta - 1.0, d.V_active] for d in designs])
    Fn = _normalize_front(F)

    # Uredimo po naraščajočem f1 (= padajočem η) -> najprej najmanj učinkoviti.
    order = np.argsort(F[:, 0])

    # Izračunaj kumulativno razdaljo vzdolž zaporedja v normalizirani ravnini:
    pts = Fn[order]
    diffs = np.diff(pts, axis=0)
    dist = np.zeros(len(order))
    dist[1:] = np.cumsum(np.linalg.norm(diffs, axis=1))
    total = dist[-1] if dist[-1] > 0 else 1.0

    # Ciljne pozicije (0, 0.25, 0.5, 0.75, 1.0) * total
    targets = [0.0, 0.25, 0.50, 0.75, 1.0]
    chosen_in_order: list[int] = []
    for t in targets:
        idx_in_order = int(np.argmin(np.abs(dist - t * total)))
        chosen_in_order.append(idx_in_order)

    # Odpravimo duplikate (lahko se zgodi, da sta 2 cilja padla na isti idx):
    seen: set[int] = set()
    unique_in_order: list[int] = []
    for idx in chosen_in_order:
        if idx in seen:
            # poišči najbližji še nezaseden
            for delta in range(1, len(order)):
                for cand in (idx - delta, idx + delta):
                    if 0 <= cand < len(order) and cand not in seen:
                        idx = cand
                        break
                else:
                    continue
                break
        seen.add(idx)
        unique_in_order.append(idx)

    # Mapiraj nazaj v izvirne indekse:
    indices = [int(order[i]) for i in unique_in_order]
    selected = [designs[i] for i in indices]
    # order[] je po f1 naraščajoče = η padajoče: prvi element je max_eta,
    # zadnji je min_V. (Pareto: večji η = večji V.)
    labels = ["max_eta", "high_eta_mid", "mid", "low_eta_mid", "min_V"]
    return selected, indices, labels


# =============================================================================
# Pomožni izpis
# =============================================================================

def summary_table(
    selected: Sequence[MotorDesign],
    labels: Sequence[str],
) -> str:
    """Lepo formatirana tabela za izbranih 5 rešitev."""
    lines = []
    header = (
        f"{'#':<3} {'oznaka':<14} {'D_r':>7} {'L_r':>7} {'D_se':>7} "
        f"{'q':>4} {'N_r':>4} {'η':>7} {'V':>9} {'P_Fe':>7} {'P_Cu':>7} {'I_m':>7}"
    )
    units = (
        f"{'':<3} {'':<14} {'[mm]':>7} {'[mm]':>7} {'[mm]':>7} "
        f"{'':>4} {'':>4} {'[%]':>7} {'[cm³]':>9} {'[W]':>7} {'[W]':>7} {'[A]':>7}"
    )
    lines.append(header)
    lines.append(units)
    lines.append("-" * len(header))
    for i, (d, lab) in enumerate(zip(selected, labels), start=1):
        lines.append(
            f"{i:<3} {lab:<14} "
            f"{d.genes.D_r*1e3:7.1f} {d.L_r*1e3:7.1f} {d.D_se*1e3:7.1f} "
            f"{d.genes.q:>4} {d.genes.N_r:>4} "
            f"{d.eta*100:7.3f} {d.V_active*1e6:9.1f} "
            f"{d.P_fe_total:7.0f} {d.P_cu_total:7.0f} {d.I_m:7.1f}"
        )
    return "\n".join(lines)


# =============================================================================
# Smoke entrypoint
# =============================================================================
if __name__ == "__main__":
    from .inputs import MachineInputs, MaterialParams, DesignBounds
    from .losses import IronLossModel
    from .optimization import run_nsga2, designs_from_result

    machine = MachineInputs()
    material = MaterialParams()
    bounds = DesignBounds()
    loss = IronLossModel.fit_default()

    result, prob = run_nsga2(machine, material, bounds, loss,
                              pop_size=60, n_gen=30, seed=42, verbose=False)
    designs = designs_from_result(result, prob)
    selected, indices, labels = select_five(designs)

    print(f"\nPareto fronta: {len(designs)} rešitev")
    print(f"Izbranih: {len(selected)}\n")
    print(summary_table(selected, labels))
