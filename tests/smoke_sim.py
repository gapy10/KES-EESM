"""Smoke test za femm_sim — zaženi prosti tek + navor na enem stroju.

To je dolgotrajen test (~5 min), ki uporablja FEMM 4.2.
Zaženi: python tests/smoke_sim.py
"""

import sys
import time
import json
from pathlib import Path

import numpy as np
import femm

# Pot k repu projekta:
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from src.inputs import MachineInputs, MaterialParams
from src.losses import IronLossModel
from src.analytical import MotorDesignGenes, analyze, pretty_print
from src.femm_model import build_motor
from src.femm_sim import no_load, torque


def main():
    machine = MachineInputs()
    material = MaterialParams()
    loss = IronLossModel.fit_default()

    # Naloži 'mid' rešitev:
    data = json.loads(Path("outputs/selected5.json").read_text(encoding="utf-8"))
    mid = next(d for d in data if d["label"] == "mid")
    genes = MotorDesignGenes(
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
    d = analyze(genes, machine, material, loss)
    print(pretty_print(d))
    print()

    fem_excited = "outputs/fem/mid_excited.fem"
    fem_noload = "outputs/fem/mid_noload.fem"
    fem_torque = "outputs/fem/mid_torque.fem"

    # ---- Build model ----
    t0 = time.time()
    build_motor(d, fem_excited, finalize_geometry=True, open_femm=True, close_femm=False)
    print(f"[1/3] Model zgrajen ({time.time()-t0:.1f}s) -> {fem_excited}")

    # ---- Prosti tek ----
    # step_deg=30 -> 5 pozicij (znotraj enega el. periodu), kar omogoči
    # smiselno dPsi/dt = U_ind izračun. Polni zagon je v src.run_femm_pareto.
    t0 = time.time()
    nl = no_load(d, rpm=machine.n_c, step_deg=30.0, save_fem_path=fem_noload)
    print(f"[2/3] Prosti tek ({time.time()-t0:.1f}s)")
    print(f"      U_ind amp = {nl.U_ind_amp:.1f} V   RMS = {nl.U_ind_rms:.1f} V")
    print(f"      U_grid = √2·U_f = {(2**0.5)*machine.U_f:.1f} V (amp), 0.95×={(2**0.5)*machine.U_f*0.95:.1f}")
    print(f"      B_max zob = {nl.B_max_tooth:.3f} T")
    print(f"      B_max jarem = {nl.B_max_yoke:.3f} T")

    # ---- Navor ----
    t0 = time.time()
    tq = torque(d, step_deg=30.0, save_fem_path=fem_torque)
    print(f"[3/3] Navor ({time.time()-t0:.1f}s)")
    print(f"      M_1.harm  = {tq.M_fundamental:.2f} Nm  (analitično: {machine.torque_c:.2f} Nm) "
          f"odstopanje {(tq.M_fundamental-machine.torque_c)/machine.torque_c*100:+.1f} %")
    print(f"      M_2.harm  = {tq.M_second:.2f} Nm")
    print(f"      M_max     = {tq.M_max:.2f} Nm   M_avg = {tq.M_avg:.2f} Nm")
    print(f"      FFT [{', '.join(f'{a:.2f}' for a in tq.fft_amplitudes[:5])}] Nm")

    femm.closefemm()

    # Shrani rezultate v json:
    out = {
        "no_load": {
            "U_ind_amp": float(nl.U_ind_amp),
            "U_ind_rms": float(nl.U_ind_rms),
            "B_max_tooth": float(nl.B_max_tooth),
            "B_max_yoke": float(nl.B_max_yoke),
            "theta_deg_mech": nl.theta_deg_mech.tolist(),
            "flux_A": nl.flux_A.tolist(),
        },
        "torque": {
            "M_avg": float(tq.M_avg),
            "M_max": float(tq.M_max),
            "M_fundamental": float(tq.M_fundamental),
            "M_second": float(tq.M_second),
            "M_ripple_pp": float(tq.M_ripple_pp),
            "M_ripple_percent": float(tq.M_ripple_percent),
            "fft_amplitudes": tq.fft_amplitudes[:10].tolist(),
            "theta_deg_mech": tq.theta_deg_mech.tolist(),
            "torque": tq.torque.tolist(),
        },
    }
    Path("outputs/mid_sim.json").write_text(
        json.dumps(out, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    print(f"\nRezultati: outputs/mid_sim.json")


if __name__ == "__main__":
    main()
