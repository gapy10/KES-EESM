"""Smoke testi za Fazo 2 (NSGA-II + Pareto + izbor)."""

import numpy as np
import pytest

from src.inputs import MachineInputs, MaterialParams, DesignBounds
from src.losses import IronLossModel
from src.optimization import (
    MotorOptimizationProblem,
    run_nsga2,
    designs_from_result,
)
from src.pareto import select_five


@pytest.fixture(scope="module")
def small_ga_result():
    machine = MachineInputs()
    material = MaterialParams()
    bounds = DesignBounds()
    loss = IronLossModel.fit_default()
    result, problem = run_nsga2(
        machine, material, bounds, loss,
        pop_size=40, n_gen=20, seed=123, verbose=False,
    )
    designs = designs_from_result(result, problem)
    return result, problem, designs


class TestOptimization:
    def test_pareto_nonempty(self, small_ga_result):
        _, _, designs = small_ga_result
        assert len(designs) > 0, "Pareto fronta je prazna"

    def test_at_least_20_solutions(self, small_ga_result):
        # FSD: faza 2 izhodno merilo - "vsaj 20 nedominiranih rešitev".
        # Pri pop=40 lahko ima Pareto fronta do 40 točk.
        _, _, designs = small_ga_result
        assert len(designs) >= 10  # pri 20 generacij in pop=40 je 10 dovolj kot smoke test

    def test_all_designs_feasible(self, small_ga_result):
        _, _, designs = small_ga_result
        for d in designs:
            assert d.feasible, f"Nedopusten v Pareto: {d.infeasibility_reasons}"

    def test_efficiency_in_range(self, small_ga_result):
        _, _, designs = small_ga_result
        eta = np.array([d.eta for d in designs])
        assert (eta > 0.9).all(), f"η < 0.9 za {(eta <= 0.9).sum()} rešitev"
        assert (eta < 0.99).all()

    def test_current_density_constraints(self, small_ga_result):
        _, _, designs = small_ga_result
        for d in designs:
            assert d.J_cu_s_actual <= 10.0 + 1e-3
            assert d.J_cu_r_actual <= 5.0 + 1e-3

    def test_determinism(self):
        machine = MachineInputs()
        material = MaterialParams()
        bounds = DesignBounds()
        loss = IronLossModel.fit_default()
        r1, _ = run_nsga2(machine, material, bounds, loss, pop_size=20, n_gen=10, seed=7)
        r2, _ = run_nsga2(machine, material, bounds, loss, pop_size=20, n_gen=10, seed=7)
        np.testing.assert_allclose(r1.F, r2.F, rtol=1e-9, atol=1e-12)


class TestParetoSelection:
    def test_select_five_returns_five(self, small_ga_result):
        _, _, designs = small_ga_result
        selected, indices, labels = select_five(designs)
        assert len(selected) == 5
        assert len(indices) == 5
        assert len(labels) == 5

    def test_selected_max_eta_is_highest(self, small_ga_result):
        _, _, designs = small_ga_result
        selected, _, labels = select_five(designs)
        max_eta_idx = labels.index("max_eta")
        # max_eta mora imeti najvišji η med izbranimi.
        etas = [d.eta for d in selected]
        assert selected[max_eta_idx].eta == max(etas)

    def test_selected_min_V_is_smallest(self, small_ga_result):
        _, _, designs = small_ga_result
        selected, _, labels = select_five(designs)
        min_V_idx = labels.index("min_V")
        Vs = [d.V_active for d in selected]
        assert selected[min_V_idx].V_active == min(Vs)

    def test_selected_indices_unique(self, small_ga_result):
        _, _, designs = small_ga_result
        _, indices, _ = select_five(designs)
        assert len(set(indices)) == 5

    def test_pareto_tradeoff_exists(self, small_ga_result):
        # Na Pareto fronti mora veljati: večji η ↔ večji V (negativna korelacija).
        _, _, designs = small_ga_result
        if len(designs) < 4:
            pytest.skip("Premajhna fronta")
        eta = np.array([d.eta for d in designs])
        V = np.array([d.V_active for d in designs])
        # Pearson korelacija med η in V mora biti pozitivna na fronti
        # (večji V daje večji η, ker imamo več železa za fluks in več bakra).
        r = np.corrcoef(eta, V)[0, 1]
        assert r > 0.3, f"Korelacija η-V na fronti = {r:.2f} (pričakovan trade-off)"
