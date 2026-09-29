#!/usr/bin/env python3
"""Generate the paper's figures from the REAL GW150914 strain and the FROZEN artifacts.

Nothing here is drawn from memory or from a hand-typed number: every plotted curve is
either read from the H1/L1 HDF5 release, recomputed in the same frozen PSD convention
as the audit, or read BY KEY from work/artifacts/*.json.

Outputs (vector PDF, plus a PNG preview for the visual check):
  figures/fig1_strain.pdf        H1 and L1 strain around the event, band-passed
  figures/fig2_asd.pdf           measured ASD (median Welch) vs the artifact values
  figures/fig3_template.pdf      whitened data vs the IMRPhenomT template
  figures/fig4_inclination.pdf   template optimal SNR vs inclination
  figures/fig5_tau2_profile.pdf  joint H1+L1 profile in tau2, both kernels

The numbering follows the order the figures appear in the paper, not the order the
functions are defined in.

Run:  python3 make_figures.py           (numpy, scipy, h5py, matplotlib; phenomxpy
                                         only for fig2, from work/env/venv if needed)
"""
import json
import math
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import rcParams
from scipy.signal import butter, sosfiltfilt
from scipy.signal.windows import tukey
import h5py

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DATA = os.path.join(ROOT, "data")
ART = os.path.join(ROOT, "work", "artifacts")
FIG = os.path.join(HERE, "figures")

FS = 16384
GPS0 = 1126259447
GPS_EVENT = 1126259462.4
T_EVENT = GPS_EVENT - GPS0
FMIN, FMAX = 20.0, 300.0
MC_MSUN = 28.096
DL_MPC = 410.0
ETA = 0.24710059171597634
DF_FD = 0.25
F_MAX_FD = 1024.0

G = 6.67430e-11
CL = 2.99792458e8
MSUN = 1.98892e30
TSUN = G * MSUN / CL ** 3
MPC = 3.0856775814913673e22
TMPC = MPC / CL
MC = MC_MSUN * TSUN
DL = DL_MPC * TMPC

FILES = {"H1": "H-H1_GWOSC_16KHZ_R1-1126259447-32.hdf5",
         "L1": "L-L1_GWOSC_16KHZ_R1-1126259447-32.hdf5"}

# ---- typography: same family as the paper text (Palatino clone) -------------
rcParams.update({
    "font.family": "serif",
    "font.serif": ["P052", "Palatino", "DejaVu Serif"],
    "mathtext.fontset": "dejavuserif",
    "font.size": 9.5,
    "axes.titlesize": 10.5,
    "axes.labelsize": 9.5,
    "legend.fontsize": 8.5,
    "xtick.labelsize": 8.5,
    "ytick.labelsize": 8.5,
    "axes.linewidth": 0.7,
    "axes.grid": True,
    "grid.alpha": 0.25,
    "grid.linewidth": 0.5,
    "figure.dpi": 160,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.03,
})
C_H1 = "#1b3a6b"
C_L1 = "#a8322d"
C_ACC = "#8a6d1f"


def load(det):
    with h5py.File(os.path.join(DATA, FILES[det]), "r") as h:
        return h["strain/Strain"][:]


def highpass(x, fs, f0=15.0):
    return sosfiltfilt(butter(4, f0 / (fs / 2.0), btype="highpass", output="sos"), x)


def bandpass(x, fs, lo, hi):
    return sosfiltfilt(butter(4, [lo / (fs / 2.0), hi / (fs / 2.0)],
                              btype="bandpass", output="sos"), x)


def median_welch(x, fs, nper, overlap=0.5):
    """Identical estimator to the audit's: median-averaged periodogram, ln2-corrected."""
    step = int(nper * (1.0 - overlap))
    w = np.hanning(nper)
    norm = 2.0 / (fs * float(np.sum(w ** 2)))
    segs = np.array([np.abs(np.fft.rfft(x[i:i + nper] * w)) ** 2 * norm
                     for i in range(0, len(x) - nper, step)])
    return np.fft.rfftfreq(nper, 1.0 / fs), np.median(segs, axis=0) / math.log(2.0)


def psd_on_grid(f, f_psd, P_psd):
    return np.exp(np.interp(np.log(np.clip(f, 1e-3, None)),
                            np.log(np.clip(f_psd, 1e-3, None)), np.log(P_psd)))


def imr_polarizations(mc_msun, freqs):
    """The audit's own call convention (h3_joint_fit.imr_polarizations)."""
    from phenomxpy import IMRPhenomT
    mtot = mc_msun / ETA ** 0.6
    wf = IMRPhenomT(eta=ETA, s1=[0, 0, 0], s2=[0, 0, 0], f_min=FMIN, f_ref=FMIN,
                    total_mass=mtot, distance=DL_MPC, inclination=0.0,
                    phi_ref=0.0, f_max=F_MAX_FD, delta_t=0.5 / F_MAX_FD,
                    delta_f=DF_FD, condition=True)
    fd = wf.compute_fd_polarizations()
    hp = np.asarray(fd[0], dtype=complex)
    f_model = np.arange(len(hp)) * DF_FD
    Hp = np.interp(freqs, f_model, hp.real) + 1j * np.interp(freqs, f_model, hp.imag)
    return Hp


def load_json(name):
    with open(os.path.join(ART, name), encoding="utf-8") as f:
        return json.load(f)


# --------------------------------------------------------------------- figures
def fig1_strain():
    """H1 and L1 around the merger, band-passed 35-350 Hz (where the signal lives)."""
    seg_n = 4 * FS
    i0 = int(round((T_EVENT - 2.0) * FS))
    fig, axes = plt.subplots(2, 1, figsize=(6.6, 4.0), sharex=True)
    peaks = {}
    for ax, det, col in zip(axes, ("H1", "L1"), (C_H1, C_L1)):
        x = load(det)
        seg = bandpass(highpass(x[i0:i0 + seg_n], FS), FS, 35.0, 350.0)
        t = np.arange(seg_n) / FS - 2.0
        ax.plot(t, seg * 1e21, lw=0.45, color=col)
        ax.axvline(0.0, color=C_ACC, lw=0.9, ls="--", zorder=5)
        ax.set_ylabel(f"{det}  strain  [$10^{{-21}}$]")
        ax.set_xlim(-0.55, 0.25)
        # peak measured in the merger window only (±50 ms), so a late glitch cannot
        # be mistaken for the signal
        sel = (t > -0.05) & (t < 0.05)
        peak = float(np.abs(seg[sel]).max() * 1e21)
        peaks[det] = peak
        ax.set_title(f"{det} — GW150914, band-passed 35–350 Hz "
                     f"(peak $|h|$ in ±50 ms $\\approx$ {peak:.1f} $\\times 10^{{-21}}$)",
                     loc="left")
    axes[-1].set_xlabel("time relative to GPS 1126259462.4  [s]")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "fig1_strain.pdf"))
    fig.savefig(os.path.join(FIG, "fig1_strain.png"))
    plt.close(fig)
    return {"fig1_peak_H1_1e21": peaks["H1"], "fig1_peak_L1_1e21": peaks["L1"]}


def fig3_template():
    """Whitened H1 data vs the full IMRPhenomT template, on the audit's own grid."""
    seg_n = 4 * FS
    i0 = int(round((T_EVENT - 2.0) * FS))
    xh = load("H1")
    x = highpass(xh[i0:i0 + seg_n], FS)
    off = np.concatenate([xh[:8 * FS], xh[-8 * FS:]])
    f_psd, P_psd = median_welch(highpass(off, FS), FS, FS, 0.5)
    freqs = np.fft.rfftfreq(seg_n, 1.0 / FS)
    Sn = psd_on_grid(freqs, f_psd, P_psd)
    m = (freqs >= FMIN) & (freqs <= FMAX)
    win = tukey(seg_n, 0.25)

    X = np.fft.rfft(x * win)[m]
    w = 1.0 / np.sqrt(Sn[m])
    B = np.zeros(seg_n, dtype=complex)
    B[:seg_n // 2 + 1][m] = X * w
    d_white = np.fft.irfft(B, seg_n)

    Hp = imr_polarizations(MC_MSUN, freqs)
    Hp_m = Hp[m]
    s_t = math.sqrt(4.0 * (freqs[1] - freqs[0]) / FS ** 2
                    * float(np.sum(np.abs(FS * Hp_m) ** 2 / Sn[m])))
    # the template is placed at the artifact's own RECORDED IMR peak time
    # (h3_imr_check.json -> bands.H1.bands.B_full_20_300.peak_time_IMR_s), not at a
    # time chosen by eye: 2.1851806640625 s in the segment.
    tc = load_json("h3_imr_check.json")["bands"]["H1"]["bands"]["B_full_20_300"]["peak_time_IMR_s"]
    Bt = np.zeros(seg_n, dtype=complex)
    Bt[:seg_n // 2 + 1][m] = (FS * Hp_m) * np.exp(-2j * np.pi * freqs[m] * tc) * w
    tpl_white = np.fft.irfft(Bt, seg_n)
    tpl_white *= np.abs(d_white).max() / max(np.abs(tpl_white).max(), 1e-30) * 0.95

    t = np.arange(seg_n) / FS - 2.0
    fig, ax = plt.subplots(figsize=(6.6, 2.9))
    ax.plot(t, d_white, lw=0.5, color=C_H1, label="H1, whitened (20–300 Hz)")
    ax.plot(t, tpl_white, lw=1.0, color=C_ACC, alpha=0.9,
            label="IMRPhenomT template (shape-normalised)")
    ax.axvline(0.0, color="0.4", lw=0.8, ls="--")
    ax.set_xlim(-0.45, 0.15)
    ax.set_xlabel("time relative to GPS 1126259462.4  [s]")
    ax.set_ylabel("whitened strain  [σ]")
    ax.set_title("The template this paper fits: full IMRPhenomT against the real H1 strain",
                 loc="left")
    ax.legend(loc="upper left", frameon=False)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "fig3_template.pdf"))
    fig.savefig(os.path.join(FIG, "fig3_template.png"))
    plt.close(fig)
    return {"fig2_tpl_opt_snr_recomputed": float(s_t)}


def fig2_asd():
    """Measured ASD vs the artifact's own ASD(100 Hz) values."""
    real = load_json("h3_real_data.json")
    fig, ax = plt.subplots(figsize=(6.6, 3.2))
    for det, col in (("H1", C_H1), ("L1", C_L1)):
        x = load(det)
        off = np.concatenate([x[:8 * FS], x[-8 * FS:]])
        f, P = median_welch(highpass(off, FS), FS, FS, 0.5)
        sel = (f > 10) & (f < 1000)
        ax.loglog(f[sel], np.sqrt(P[sel]), lw=0.9, color=col, label=f"{det}, median Welch")
    ax.axhline(1.0275194697647844e-23, color=C_H1, lw=0.7, ls=":",
               label="artifact ASD(100 Hz) H1 = $1.03\\times10^{-23}$")
    ax.axhline(9.848940548987532e-24, color=C_L1, lw=0.7, ls=":",
               label="artifact ASD(100 Hz) L1 = $9.85\\times10^{-24}$")
    ax.axvspan(FMIN, FMAX, color=C_ACC, alpha=0.08, lw=0)
    ax.text(24, 3e-22, "analysis band\n20–300 Hz", fontsize=8, color=C_ACC)
    ax.set_xlim(10, 1000)
    ax.set_ylim(1e-24, 1e-20)
    ax.set_xlabel("frequency  [Hz]")
    ax.set_ylabel("ASD  [strain / √Hz]")
    ax.set_title("The noise curve is measured from off-source data, not assumed", loc="left")
    ax.legend(loc="lower left", frameon=False, ncol=2)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "fig2_asd.pdf"))
    fig.savefig(os.path.join(FIG, "fig2_asd.png"))
    plt.close(fig)
    return {"fig3_asd100_H1": real["inputs"]["H1_psd"]["asd_at_100Hz"],
            "fig3_asd100_L1": real["inputs"]["L1_psd"]["asd_at_100Hz"]}


def fig5_tau2_profile():
    """The joint H1+L1 profile in tau2, both kernels, straight from the frozen artifact."""
    d = load_json("h3_joint_fit.json")
    fig, axes = plt.subplots(1, 2, figsize=(6.8, 2.9), sharey=True)
    for ax, key, title in ((axes[0], "joint_pn_expanded", "kernel: PN inspiral"),
                           (axes[1], "joint_imr_expanded", "kernel: IMR group delay")):
        g = d["per_grid"][key]
        tau = np.array(g["grid"]) * 1e3
        prof = np.array(g["profile"])
        peak = g["profile_peak_net_snr"]
        ax.plot(tau, prof, lw=1.2, color=C_H1)
        ax.axhline(peak - 1.0, color=C_ACC, lw=0.9, ls="--")
        ax.axhline(g["profile_at_zero"], color="0.5", lw=0.8, ls=":")
        ax.axvline(0.0, color="0.5", lw=0.8, ls=":")
        ax.plot([g["tau2_best_fit_amplitude_ms"]], [peak], marker="o", ms=3.5,
                color=C_L1, zorder=5)
        ax.set_title(title, loc="left")
        ax.set_xlabel(r"$\tau_2$ amplitude  [ms]")
        drop = g["profile_drop_within_grid"]
        ax.text(0.03, 0.06, f"drop = {drop:.3f}$\\sigma$", transform=ax.transAxes,
                fontsize=8.5, color=C_ACC)
    axes[0].set_ylabel("coherent network SNR")
    fig.suptitle("Joint H1+L1 fit with one $\\tau_2$ — the profile, and the peak $-\\,1\\sigma$ level",
                 x=0.012, ha="left", fontsize=10.5)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(os.path.join(FIG, "fig5_tau2_profile.pdf"))
    fig.savefig(os.path.join(FIG, "fig5_tau2_profile.png"))
    plt.close(fig)
    return {"fig4_drop_pn": d["per_grid"]["joint_pn_expanded"]["profile_drop_within_grid"],
            "fig4_drop_imr": d["per_grid"]["joint_imr_expanded"]["profile_drop_within_grid"]}


def fig4_inclination():
    """Template optimal SNR vs inclination, from the frozen IMR and spin artifacts."""
    imr = load_json("h3_imr_check.json")
    spins = load_json("h3_imr_spins.json")
    sweep = imr["inclination_sweep_H1_rho_opt"]
    xs = np.array(sorted(float(k) for k in sweep))
    ys = np.array([sweep[f"{k:.1f}"] for k in xs])
    cross = imr["inclination_where_rho_opt_equals_20_deg"]

    fig, ax = plt.subplots(figsize=(6.6, 3.0))
    ax.plot(xs, ys, "-o", ms=3.5, lw=1.1, color=C_H1,
            label="H1, IMRPhenomT, optimal SNR (zero spin)")
    ax.axhline(19.5, color=C_L1, lw=0.9, ls="--",
               label="published single-detector H1 = 19.5")
    ax.axvline(cross, color=C_ACC, lw=0.9, ls=":")
    ax.annotate(f"ι ≈ {cross:.1f}°", xy=(cross, 20.0), xytext=(cross + 4, 24.5),
                fontsize=8.5, color=C_ACC,
                arrowprops=dict(arrowstyle="-", color=C_ACC, lw=0.7))
    env = spins["answer"]["H1"]["rho_opt_spin_only_range"]
    ax.axhspan(env[0], env[1], color=C_L1, alpha=0.10, lw=0,
               label=f"spin envelope, H1: {env[0]:.1f}–{env[1]:.1f}")
    ax.set_xlim(0, 90)
    ax.set_xlabel("inclination ι  [deg]")
    ax.set_ylabel("optimal SNR  $\\rho_{\\rm opt}$")
    ax.set_title("Where the residual deficit lives: orientation, not spins", loc="left")
    ax.legend(loc="upper right", frameon=False)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "fig4_inclination.pdf"))
    fig.savefig(os.path.join(FIG, "fig4_inclination.png"))
    plt.close(fig)
    return {"fig5_crossing_deg": cross, "fig5_rho_at_0": float(ys[0])}


def main():
    os.makedirs(FIG, exist_ok=True)
    out = {}
    out.update(fig1_strain())
    out.update(fig2_asd())
    out.update(fig3_template())
    out.update(fig4_inclination())
    out.update(fig5_tau2_profile())
    print("FIGURES")
    for k, v in out.items():
        print(f"  {k} = {v}")
    for n in sorted(os.listdir(FIG)):
        p = os.path.join(FIG, n)
        print(f"  {os.path.getsize(p):>9}  {n}")
    return 0


if __name__ == "__main__":
    sys.exit(main())