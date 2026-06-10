"""
Večkriterijska optimizacija sinhronskega motorja z NSGA-II.

Implementira točke G in H seminarske naloge:
    G: analitičen preračun vključen v GA,
    H: dvokriterijska minimizacija (1/η in V_active) -> Pareto fronta.

Optimizacijski problem (9 spremenljivk):
    x[0] = D_r        ∈ [D_r_min, D_r_max]            (zvezno)
    x[1] = l_to_D     ∈ [l_to_D_min, l_to_D_max]      (zvezno)
    x[2] = B_delta    ∈ [B_delta_min, B_delta_max]    (zvezno)
    x[3] = B_ds       ∈ [B_ds_min, B_ds_max]          (zvezno)
    x[4] = B_sy       ∈ [B_sy_min, B_sy_max]          (zvezno)
    x[5] = J_cu_s     ∈ [J_cu_s_min, J_cu_s_max]      (zvezno)
    x[6] = J_cu_r     ∈ [J_cu_r_min, J_cu_r_max]      (zvezno)
    x[7] = N_r        ∈ [N_r_min, N_r_max]            (zaokroženo na int)
    x[8] = q_idx      ∈ [0, 4]                        (zaokroženo, ima 5 košev za q)
    x[9] = delta      ∈ [delta_min, delta_max]        (zvezno, zračna reža [m])

Cilji (minimizacija):
    f1 = 1/η - 1         (max izkoristek)
    f2 = V_active        (min volumen aktivnega dela)

Trde omejitve (g_i ≤ 0):
    g1 = J_cu_s_actual - 10
    g2 = J_cu_r_actual - 5
    g3 = feasibility flag iz analitične funkcije (≥1 če je infeasible)

Opomba: U_ind ≤ 0.95·U_grid (FR-3.8) je tavtologija v analitičnem modelu,
saj se N_s računa tako, da E = U_f. Ta omejitev se preverja šele v Fazi 3
(FEMM verifikacija) z dejansko FEMM induktivnostjo in zaokroženim N_s.
Spodnja meja l/D = 0.4 je vgrajena v `DesignBounds`, zato je ne potrebujemo
kot eksplicitno omejitev.

Algoritem: NSGA-II (Deb et al., 2002).

Vir:
    Deb, K., Pratap, A., Agarwal, S., Meyarivan, T.,
    "A Fast and Elitist Multiobjective Genetic Algorithm: NSGA-II",
    IEEE Trans. Evol. Comput., 6(2):182-197, 2002.
"""

from __future__ import annotations

import math
from dataclasses import asdict
from pathlib import Path
from typing import Optional

import numpy as np

from pymoo.algorithms.moo.nsga2 import NSGA2
from pymoo.core.problem import Problem
from pymoo.operators.crossover.sbx import SBX
from pymoo.operators.mutation.pm import PM
from pymoo.operators.sampling.lhs import LHS
from pymoo.optimize import minimize
from pymoo.termination import get_termination

from .inputs import MachineInputs, MaterialParams, DesignBounds
from .losses import IronLossModel
from .analytical import MotorDesignGenes, MotorDesign, analyze


# =============================================================================
# Problem definition
# =============================================================================

class MotorOptimizationProblem(Problem):
    """NSGA-II problem ovojnica za analitični stroj.

    9 spremenljivk, 2 cilja, 3 omejitve (J_cu_s, J_cu_r in penalty
    za ostale kršitve iz `MotorDesign.feasible`).
    """

    N_VAR = 10
    N_OBJ = 2
    N_CONSTR = 3    # J_cu_s, J_cu_r, feasibility-penalty
                    # (Inflacijski FEMM-navor filter g4 odstranjen: model
                    #  targetira poštene 50 kW, FEMM navor je rezultat, ne omejitev.)

    def __init__(
        self,
        machine: MachineInputs,
        material: MaterialParams,
        bounds: DesignBounds,
        loss_model: IronLossModel,
        U_ind_limit_factor: float = 0.95,
        **analyze_kwargs,
    ):
        self.machine = machine
        self.material = material
        self.design_bounds = bounds
        self.loss_model = loss_model
        self.U_ind_limit_factor = U_ind_limit_factor
        self.analyze_kwargs = analyze_kwargs

        # Meje za pymoo: vse spremenljivke zvezne; N_r in q_idx interpretirani diskretno.
        xl = np.array([
            bounds.D_r_min,
            bounds.l_to_D_min,
            bounds.B_delta_min,
            bounds.B_ds_min,
            bounds.B_sy_min,
            bounds.J_cu_s_min,
            bounds.J_cu_r_min,
            float(bounds.N_r_min),
            0.0,
            bounds.delta_min,
        ])
        xu = np.array([
            bounds.D_r_max,
            bounds.l_to_D_max,
            bounds.B_delta_max,
            bounds.B_ds_max,
            bounds.B_sy_max,
            bounds.J_cu_s_max,
            bounds.J_cu_r_max,
            float(bounds.N_r_max),
            float(len(bounds.q_choices) - 1),  # 0..4 za 5 vrednosti
            bounds.delta_max,
        ])
        super().__init__(n_var=self.N_VAR, n_obj=self.N_OBJ, n_ieq_constr=self.N_CONSTR,
                         xl=xl, xu=xu)

    # ------------------------------------------------------------------ helpers
    def decode(self, x: np.ndarray) -> MotorDesignGenes:
        """Dekodira en vektor x v MotorDesignGenes."""
        q_idx = int(round(np.clip(x[8], 0, len(self.design_bounds.q_choices) - 1)))
        return MotorDesignGenes(
            D_r=float(x[0]),
            l_to_D=float(x[1]),
            B_delta=float(x[2]),
            B_ds=float(x[3]),
            B_sy=float(x[4]),
            J_cu_s=float(x[5]),
            J_cu_r=float(x[6]),
            N_r=int(round(x[7])),
            q=float(self.design_bounds.q_choices[q_idx]),
            delta=float(x[9]),
        )

    def evaluate_one(self, x: np.ndarray) -> tuple[MotorDesign, np.ndarray, np.ndarray]:
        """Vrne (design, f, g) za en sample."""
        genes = self.decode(x)
        d = analyze(genes, self.machine, self.material, self.loss_model,
                    **self.analyze_kwargs)

        # Cilji (minimizacija):
        if d.eta > 0:
            f1 = 1.0 / d.eta - 1.0
        else:
            f1 = 10.0  # zelo kazenska vrednost
        f2 = d.V_active

        # Trde omejitve g ≤ 0 (NSGA-II konvencija):
        g1 = d.J_cu_s_actual - self.material.J_cu_s_max
        g2 = d.J_cu_r_actual - self.material.J_cu_r_max
        # Ostale geometrijske kršitve iz `analyze()` — vsaka po 1.0 v CV:
        g3 = -1.0 if d.feasible else float(len(d.infeasibility_reasons))
        # OPOMBA: prejšnja omejitev g4 (M_FEMM_pred >= P_c_min_FEMM/ω) je bila
        # odstranjena. Bila je vezana na inflacijo P_c na 70 kW; pri poštenih
        # 50 kW vsak analitični design po konstrukciji targetira M_c = 68.21 Nm,
        # dejanski FEMM navor pa je rezultat verifikacije, ne projektni filter.

        f = np.array([f1, f2])
        g = np.array([g1, g2, g3])
        return d, f, g

    # ------------------------------------------------------------------ pymoo API
    def _evaluate(self, X, out, *args, **kwargs):
        F = np.zeros((len(X), self.N_OBJ))
        G = np.zeros((len(X), self.N_CONSTR))
        designs: list[MotorDesign] = []
        for i, x in enumerate(X):
            d, f, g = self.evaluate_one(x)
            F[i] = f
            G[i] = g
            designs.append(d)
        out["F"] = F
        out["G"] = G
        # `out["designs"]` lahko prenesemo za poznejšo uporabo (pymoo
        # vmes ne uporablja, vendar ga lahko pohranimo prek callback-a).


# =============================================================================
# Run NSGA-II
# =============================================================================

def run_nsga2(
    machine: MachineInputs,
    material: MaterialParams,
    bounds: DesignBounds,
    loss_model: IronLossModel,
    pop_size: int = 100,
    n_gen: int = 50,
    seed: int = 42,
    verbose: bool = False,
):
    """Zažene NSGA-II in vrne rezultat.

    Returns:
        (result, problem) — `result.X` je `(N, 9)` matrika nedominiranih
        odločitvenih vektorjev, `result.F` je `(N, 2)` ciljnih vrednosti.
    """
    problem = MotorOptimizationProblem(machine, material, bounds, loss_model)

    algorithm = NSGA2(
        pop_size=pop_size,
        sampling=LHS(),
        crossover=SBX(prob=0.9, eta=15),
        mutation=PM(prob=1.0 / problem.N_VAR, eta=20),
        eliminate_duplicates=True,
    )

    termination = get_termination("n_gen", n_gen)

    result = minimize(
        problem,
        algorithm,
        termination,
        seed=seed,
        save_history=False,
        verbose=verbose,
    )
    return result, problem


# =============================================================================
# Save Pareto front to .npz + JSON
# =============================================================================

def designs_from_result(result, problem: MotorOptimizationProblem) -> list[MotorDesign]:
    """Iz rezultata NSGA-II rekonstruira MotorDesign objekte za vsako rešitev."""
    designs: list[MotorDesign] = []
    if result.X is None or len(result.X) == 0:
        return designs
    for x in np.atleast_2d(result.X):
        d, _f, _g = problem.evaluate_one(x)
        designs.append(d)
    return designs


def save_pareto(
    designs: list[MotorDesign],
    F: np.ndarray,
    X: np.ndarray,
    out_path: str | Path,
) -> None:
    """Shrani Pareto fronto v .npz."""
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    eta = np.array([d.eta for d in designs])
    V = np.array([d.V_active for d in designs])

    np.savez(
        out_path,
        F=np.asarray(F),
        X=np.asarray(X),
        eta=eta,
        V_active=V,
    )


# =============================================================================
# Smoke entrypoint
# =============================================================================
if __name__ == "__main__":
    machine = MachineInputs()
    material = MaterialParams()
    bounds = DesignBounds()
    loss = IronLossModel.fit_default()

    result, prob = run_nsga2(
        machine, material, bounds, loss,
        pop_size=40, n_gen=20, seed=42, verbose=True,
    )
    if result.X is None:
        print("Pareto fronta je prazna - vse rešitve so nedopustne.")
    else:
        print(f"\nPareto fronta: {len(result.X)} rešitev")
        print(f"η razpon: {1.0 / (1.0 + result.F[:, 0].max()):.4f} ... "
              f"{1.0 / (1.0 + result.F[:, 0].min()):.4f}")
        print(f"V razpon: {result.F[:, 1].min()*1e6:.0f} ... "
              f"{result.F[:, 1].max()*1e6:.0f} cm³")
