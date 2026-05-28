"""
Risanje grafov za Fazo 2 (Pareto fronta) in Fazo 4 (U_ind, navor, FFT, primerjave).

Privzeti backend: matplotlib `Agg` (brez okna), tako da se shrani PNG.
Vse slike so generirane v 300 dpi za vključitev v poročilo (FSD razdelek 3.4
izhodno merilo).
"""

from __future__ import annotations

from pathlib import Path
from typing import Sequence, Optional

import matplotlib

matplotlib.use("Agg")  # neinteraktivni backend, primeren za batch zagone
import matplotlib.pyplot as plt
import numpy as np

from .analytical import MotorDesign

DPI = 300


def plot_pareto(
    designs: Sequence[MotorDesign],
    selected_indices: Sequence[int] = (),
    selected_labels: Sequence[str] = (),
    out_path: str | Path = "outputs/pareto.png",
    title: str = "Pareto fronta: η vs. V_active",
    replaced_flags: Optional[Sequence[bool]] = None,
    original_indices: Optional[Sequence[int]] = None,
) -> Path:
    """Nariše Pareto fronto v ravnini (V_active [cm³], η [%]) z označenimi 5 rešitvami.

    Args:
        replaced_flags: vzporeden seznam k `selected_indices`; True pomeni,
            da je bila rešitev po napaki FEMM zamenjana z najbližjim sosedom.
            Zamenjana rešitev se izriše z drugačno barvo zvezde.
        original_indices: vzporeden seznam k `selected_indices`. Za vsak
            slot, ki je bil zamenjan (`replaced_flags[i] == True`), ta indeks
            kaže na izvirno (padlo) rešitev na Pareto fronti. Izriše se kot
            prazen kvadratek + puščica do nove izbire.
    """
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    V = np.array([d.V_active * 1e6 for d in designs])  # cm³
    eta = np.array([d.eta * 100 for d in designs])     # %

    fig, ax = plt.subplots(figsize=(7, 5), dpi=DPI)
    ax.scatter(V, eta, s=22, c="#1f77b4", alpha=0.6, label="Pareto fronta")

    repl_handle = None
    sel_handle = None
    orig_handle = None
    for i, (idx, lab) in enumerate(zip(selected_indices, selected_labels)):
        is_repl = bool(replaced_flags[i]) if replaced_flags is not None else False
        color = "#ff7f0e" if is_repl else "#d62728"
        h = ax.scatter(V[idx], eta[idx], s=140, marker="*", c=color,
                       edgecolors="black", linewidths=0.6, zorder=5)
        if is_repl and repl_handle is None:
            repl_handle = h
        elif not is_repl and sel_handle is None:
            sel_handle = h
        # Po potrebi povezava do izvirne (padle) rešitve:
        if is_repl and original_indices is not None:
            orig_idx = int(original_indices[i])
            oh = ax.scatter(V[orig_idx], eta[orig_idx], s=110, marker="s",
                            facecolors="none", edgecolors="#555",
                            linewidths=1.0, zorder=4)
            if orig_handle is None:
                orig_handle = oh
            ax.annotate(
                "", xy=(V[idx], eta[idx]), xytext=(V[orig_idx], eta[orig_idx]),
                arrowprops=dict(arrowstyle="->", color="#555",
                                lw=0.8, shrinkA=4, shrinkB=6),
            )
        ax.annotate(
            lab,
            (V[idx], eta[idx]),
            xytext=(8, 4), textcoords="offset points",
            fontsize=8, color="#444",
        )

    ax.set_xlabel("Volumen aktivnega dela V [cm³]")
    ax.set_ylabel("Izkoristek η [%]")
    ax.set_title(title)
    ax.grid(True, alpha=0.3)

    handles = [ax.collections[0]]
    labels_legend = ["Pareto fronta"]
    if sel_handle is not None:
        handles.append(sel_handle); labels_legend.append("Izbrana (FEMM OK)")
    if repl_handle is not None:
        handles.append(repl_handle); labels_legend.append("Zamenjana (sosed)")
    if orig_handle is not None:
        handles.append(orig_handle); labels_legend.append("Izvirna (FEMM padla)")
    ax.legend(handles, labels_legend, loc="lower right", fontsize=8)
    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)
    return out_path


def plot_pareto_loss(
    designs: Sequence[MotorDesign],
    P_out_W: float,
    selected_indices: Sequence[int] = (),
    selected_labels: Sequence[str] = (),
    out_path: str | Path = "outputs/pareto_loss.png",
    title: str = "Pareto fronta: P_skupne izgube vs. V_active",
    replaced_flags: Optional[Sequence[bool]] = None,
    original_indices: Optional[Sequence[int]] = None,
) -> Path:
    """Pareto fronta v ravnini (V_active, P_loss).

    Y-os prikazuje skupne izgube P_loss = P_out * (1/η − 1), kar je
    "minimizacijska" različica η (smer optimuma: levo-dol).
    """
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    V = np.array([d.V_active * 1e6 for d in designs])              # cm³
    P_loss = np.array([P_out_W * (1.0 / d.eta - 1.0) for d in designs])  # W

    fig, ax = plt.subplots(figsize=(7, 5), dpi=DPI)
    ax.scatter(V, P_loss, s=22, c="#1f77b4", alpha=0.6)

    repl_handle = None
    sel_handle = None
    orig_handle = None
    for i, (idx, lab) in enumerate(zip(selected_indices, selected_labels)):
        is_repl = bool(replaced_flags[i]) if replaced_flags is not None else False
        color = "#ff7f0e" if is_repl else "#d62728"
        h = ax.scatter(V[idx], P_loss[idx], s=140, marker="*", c=color,
                       edgecolors="black", linewidths=0.6, zorder=5)
        if is_repl and repl_handle is None:
            repl_handle = h
        elif not is_repl and sel_handle is None:
            sel_handle = h
        if is_repl and original_indices is not None:
            oi = int(original_indices[i])
            oh = ax.scatter(V[oi], P_loss[oi], s=110, marker="s",
                            facecolors="none", edgecolors="#555",
                            linewidths=1.0, zorder=4)
            if orig_handle is None:
                orig_handle = oh
            ax.annotate(
                "", xy=(V[idx], P_loss[idx]), xytext=(V[oi], P_loss[oi]),
                arrowprops=dict(arrowstyle="->", color="#555",
                                lw=0.8, shrinkA=4, shrinkB=6),
            )
        ax.annotate(lab, (V[idx], P_loss[idx]),
                    xytext=(8, 4), textcoords="offset points",
                    fontsize=8, color="#444")

    ax.set_xlabel("Volumen aktivnega dela V [cm³]")
    ax.set_ylabel("Skupne izgube P_loss [W]")
    ax.set_title(title)
    ax.grid(True, alpha=0.3)

    handles = [ax.collections[0]]
    labels_legend = ["Pareto fronta"]
    if sel_handle is not None:
        handles.append(sel_handle); labels_legend.append("Izbrana (FEMM OK)")
    if repl_handle is not None:
        handles.append(repl_handle); labels_legend.append("Zamenjana (sosed)")
    if orig_handle is not None:
        handles.append(orig_handle); labels_legend.append("Izvirna (FEMM padla)")
    ax.legend(handles, labels_legend, loc="upper right", fontsize=8)
    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)
    return out_path


def plot_eta_comparison(
    labels: Sequence[str],
    eta_analyt: Sequence[float],
    eta_femm: Sequence[float],
    out_path: str | Path = "outputs/eta_comparison.png",
    title: str = "Primerjava izkoristka — analitično vs. FEMM",
) -> Path:
    """Stolpčni graf η_analitično vs. η_FEMM za 5 izbranih rešitev."""
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    n = len(labels)
    x = np.arange(n)
    w = 0.38
    a = np.asarray(eta_analyt) * 100.0
    f = np.asarray(eta_femm) * 100.0

    fig, ax = plt.subplots(figsize=(8, 4.5), dpi=DPI)
    b1 = ax.bar(x - w / 2, a, width=w, label="Analitično", color="#1f77b4")
    b2 = ax.bar(x + w / 2, f, width=w, label="FEMM", color="#ff7f0e")

    for rect, val in zip(b1, a):
        ax.text(rect.get_x() + rect.get_width() / 2, val + 0.02,
                f"{val:.2f}", ha="center", va="bottom", fontsize=7)
    for rect, val in zip(b2, f):
        ax.text(rect.get_x() + rect.get_width() / 2, val + 0.02,
                f"{val:.2f}", ha="center", va="bottom", fontsize=7)

    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=15, ha="right")
    ax.set_ylabel("Izkoristek η [%]")
    ax.set_title(title)
    ax.grid(True, axis="y", alpha=0.3)
    # Omeji prikaz na koristen razpon (vse vrednosti so visoke):
    lo = min(a.min(), f.min()) - 0.5
    hi = max(a.max(), f.max()) + 0.5
    ax.set_ylim(lo, hi)
    ax.legend(loc="lower right", fontsize=9)
    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)
    return out_path


def plot_V_comparison(
    labels: Sequence[str],
    V_cm3: Sequence[float],
    out_path: str | Path = "outputs/V_comparison.png",
    title: str = "Primerjava volumna aktivnega dela",
) -> Path:
    """Stolpčni graf V_active [cm³] za 5 izbranih rešitev."""
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    n = len(labels)
    x = np.arange(n)
    V = np.asarray(V_cm3)

    fig, ax = plt.subplots(figsize=(8, 4.5), dpi=DPI)
    bars = ax.bar(x, V, color="#2ca02c", width=0.6)
    for rect, val in zip(bars, V):
        ax.text(rect.get_x() + rect.get_width() / 2, val + V.max() * 0.01,
                f"{val:.0f}", ha="center", va="bottom", fontsize=8)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=15, ha="right")
    ax.set_ylabel("Volumen aktivnega dela V [cm³]")
    ax.set_title(title)
    ax.grid(True, axis="y", alpha=0.3)
    ax.set_ylim(0, V.max() * 1.10)
    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)
    return out_path


def plot_uind(
    theta_deg: np.ndarray,
    uind: np.ndarray,
    out_path: str | Path,
    title: str = "Inducirana napetost U_ind(θ)",
) -> Path:
    """Graf inducirane napetosti pri prostem teku."""
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(7, 4), dpi=DPI)
    ax.plot(theta_deg, uind, color="#1f77b4")
    ax.set_xlabel("Kot rotorja θ [°]")
    ax.set_ylabel("U_ind [V]")
    ax.set_title(title)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)
    return out_path


def plot_torque(
    theta_deg: np.ndarray,
    torque: np.ndarray,
    out_path: str | Path,
    title: str = "Navor M(θ)",
) -> Path:
    """Graf navorne karakteristike v odvisnosti od kolesnega kota."""
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(7, 4), dpi=DPI)
    ax.plot(theta_deg, torque, color="#2ca02c")
    ax.set_xlabel("Kolesni kot θ [°]")
    ax.set_ylabel("M [Nm]")
    ax.set_title(title)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)
    return out_path


def plot_fft(
    harmonics: np.ndarray,
    amplitudes: np.ndarray,
    out_path: str | Path,
    title: str = "Spekter navora (FFT)",
    ylabel: str = "Amplituda [Nm]",
) -> Path:
    """Palični graf harmonikov."""
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(7, 4), dpi=DPI)
    ax.bar(harmonics, amplitudes, color="#ff7f0e", width=0.5)
    ax.set_xlabel("Harmonik [-]")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(True, alpha=0.3, axis="y")
    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)
    return out_path


def plot_flux(
    theta_deg: np.ndarray,
    flux: np.ndarray,
    out_path: str | Path,
    title: str = "Verižni magnetni pretok Ψ_A(θ)",
) -> Path:
    """Graf verižnega magnetnega pretoka (za diagnostiko prostega teka)."""
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(7, 4), dpi=DPI)
    ax.plot(theta_deg, flux * 1e3, color="#9467bd", marker="o", markersize=4)
    ax.set_xlabel("Kot rotorja θ [°]")
    ax.set_ylabel("Ψ_A [mWb]")
    ax.set_title(title)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)
    return out_path


def plot_comparison_bar(
    quantities: Sequence[str],
    analytical: Sequence[float],
    femm: Sequence[float],
    out_path: str | Path,
    title: str = "Primerjava analitično vs FEMM",
    unit: str = "",
) -> Path:
    """Stolpičast graf primerjave analitičnih in FEMM vrednosti."""
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    n = len(quantities)
    x = np.arange(n)
    width = 0.35
    fig, ax = plt.subplots(figsize=(8, 4.5), dpi=DPI)
    b1 = ax.bar(x - width / 2, analytical, width, label="Analitično", color="#1f77b4")
    b2 = ax.bar(x + width / 2, femm, width, label="FEMM", color="#d62728")
    ax.set_xticks(x)
    ax.set_xticklabels(quantities, rotation=20, ha="right")
    ax.set_ylabel(unit)
    ax.set_title(title)
    ax.legend()
    ax.grid(True, alpha=0.3, axis="y")
    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)
    return out_path
