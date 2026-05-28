"""Strukturni testi za Fazo 3.

FEMM klicev tu NE testiramo — to počne `tests/smoke_sim.py`, ki traja minute.
Tukaj preverjamo le komponente, ki ne potrebujejo odprte FEMM seje.
"""

import pytest
import math
import numpy as np

from src.femm_model import (
    generate_stator_winding_layout,
    _mirror_rotate_point,
    FemmDrawOptions,
)
from src.femm_sim import NoLoadResult, TorqueResult


class TestWindingLayout:
    @pytest.mark.parametrize("q,p,Qs", [(1.0, 3, 18), (2.0, 3, 36), (3.0, 3, 54)])
    def test_layout_length(self, q, p, Qs):
        layout = generate_stator_winding_layout(Qs, p=p, m=3)
        assert len(layout) == Qs

    def test_layout_contains_all_phases(self):
        layout = generate_stator_winding_layout(36, p=3, m=3)
        symbols = set(layout)
        # Vsebovati mora vse 3 faze v obeh smereh.
        assert symbols == {"A", "a", "B", "b", "C", "c"}

    def test_layout_balanced(self):
        # Vsaka faza (vključno z obema smerema) mora imeti enako število pojavov.
        layout = generate_stator_winding_layout(36, p=3, m=3)
        for sym in ["A", "a", "B", "b", "C", "c"]:
            assert layout.count(sym) == 6  # 36 / 6 = 6 utorov na simbol

    def test_layout_q2_starts_with_A(self):
        # Za q=2, Q_s=36, p=3: prvi utor (kot 0°) je v sektorju A.
        layout = generate_stator_winding_layout(36, p=3, m=3)
        assert layout[0] == "A"

    def test_rejects_non_3phase(self):
        with pytest.raises(ValueError):
            generate_stator_winding_layout(20, p=2, m=2)


class TestMirrorRotate:
    def test_mirror_zero_rotation(self):
        # Zrcalita x→-x; brez rotacije.
        x1, y1 = _mirror_rotate_point(5.0, 10.0, 0.0)
        assert x1 == pytest.approx(-5.0)
        assert y1 == pytest.approx(10.0)

    def test_mirror_180_rotation(self):
        # Po zrcaljenju x→-x in zasuku za 180°: (-5,10) → (5,-10).
        x1, y1 = _mirror_rotate_point(5.0, 10.0, 180.0)
        assert x1 == pytest.approx(5.0)
        assert y1 == pytest.approx(-10.0, abs=1e-9)


class TestFemmDrawOptions:
    def test_defaults_reasonable(self):
        opt = FemmDrawOptions()
        assert opt.slot_opening_stator_mm > 0
        assert opt.slot_lining_mm > 0
        assert opt.slot_wedge_mm > 0
        assert 0 < opt.rotor_pole_arc_frac < 1.0
        assert 0 < opt.rotor_pole_height_frac < 0.5
        assert opt.precision > 0

    def test_override(self):
        opt = FemmDrawOptions(slot_opening_stator_mm=3.0, precision=1e-6)
        assert opt.slot_opening_stator_mm == 3.0
        assert opt.precision == 1e-6


class TestResultDataclasses:
    """Preveri strukturo `NoLoadResult` in `TorqueResult` brez FEMM klicev."""

    def test_no_load_result_fields(self):
        n = 5
        res = NoLoadResult(
            theta_deg_mech=np.linspace(0, 120, n),
            flux_A=np.zeros(n),
            U_ind_inst=np.zeros(n - 1),
            U_ind_amp=0.0,
            U_ind_rms=0.0,
        )
        assert res.theta_deg_mech.shape == (n,)
        assert res.flux_A.shape == (n,)
        # U_ind je en korak manjši (np.diff)
        assert res.U_ind_inst.shape == (n - 1,)

    def test_torque_result_fields(self):
        n = 5
        res = TorqueResult(
            theta_deg_mech=np.linspace(0, 120, n),
            torque=np.zeros(n),
            M_avg=0.0,
            M_max=0.0,
            M_fundamental=0.0,
            M_second=0.0,
            M_ripple_pp=0.0,
            M_ripple_percent=0.0,
            fft_harmonics=np.arange(3),
            fft_amplitudes=np.zeros(3),
        )
        assert res.theta_deg_mech.shape == (n,)
        assert res.torque.shape == (n,)
        assert res.fft_amplitudes.shape == res.fft_harmonics.shape
