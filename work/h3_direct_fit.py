#!/usr/bin/env python3
"""H3 issue 2: fit the delay DIRECTLY to the real GW150914 strain, and decompose
the 37.3-vs-7.27 SNR gap the owner flagged.

Two things the previous turns did NOT do:

 (A) DIRECT FIT.  h3_real_data.py only PROJECTED the delay direction off the
     nuisance span -- it never asked the data "what is the best-fit quadratic
     delay, and how well is it determined?".  Here we scan tau2 over a grid,
     build the delay-perturbed template, maximise the matched-filter SNR over
     (t_c, phi_c) at each tau2, and read off the profile likelihood.  This is a
     genuine fit to the measured strain.

 (B) SNR GAP.  The template's own optimal SNR is 37.3 (H1) but the filter
     extracts 7.27 from the data.  That factor 5.1 is unexplained unless
     decomposed.  Candidate causes, each measured here:
       B1  PSD calibration error        -- test: injection with target SNR 20
       B2  amplitude scale the data prefers (best-fit scale)
       B3  band truncation (inspiral only, no merger/ringdown) -- measure how
           much of the template's sigma sits above the ISCO frequency
       B4  template imperfection (leading-order phase) -- not separable from B3
           with this data, and said so rather than claimed.

CEILING.  Without the true IMR waveform we cannot assign the gap to B3 vs B4.
What we CAN do is show that B1 is excluded and give B2 exactly.
"""
import json
import math
import os
import sys

import numpy as np
import h5py
from scipy.signal import butter, sosfiltfilt
from scipy.signal.windows import tukey

FS = 16384
GPS0 = 1126259447
GPS_EVENT = 1126259462.4
T_EVENT = GPS_EVENT - GPS0
FMIN, FMAX = 20.0, 300.0
MC_MSUN = 28.096
DL_MPC = 410.0
SEED = 20260927
N_CAL = 64

G = 6.67430e-11
CL = 2.99792458e8
MSUN = 1.98892e30
TSUN = G * MSUN / CL ** 3
MPC = 3.0856775814913673e22
TMPC = MPC / CL
MC = MC_MSUN * TSUN
DL = DL_MPC * TMPC

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(os.path.dirname(HERE), "data")
FILES = {"H1": "H-H1_GWOSC_16KHZ_R1-1126259447-32.hdf5",
         "L1": "L-L1_GWOSC_16KHZ_R1-1126259447-32.hdf5"}


def load(det):
    with h5py.File(os.path.join(DATA, FILES[det]), "r") as h:
        return h["strain/Strain"][:], int(h["meta/GPSstart"][()])


def highpass(x, fs, f0=15.0):
    return sosfiltfilt(butter(4, f0 / (fs / 2.0), btype="highpass", output="sos"), x)


def median_welch(x, fs, nper, overlap=0.5):
    step = int(nper * (1.0 - overlap))
    w = np.hanning(nper)
    norm = 2.0 / (fs * float(np.sum(w ** 2)))
    segs = np.array([np.abs(np.fft.rfft(x[i:i + nper] * w)) ** 2 * norm
                     for i in range(0, len(x) - nper, step)])
    return np.fft.rfftfreq(nper, 1.0 / fs), np.median(segs, axis=0) / math.log(2.0)


def psd_on_grid(f, f_psd, P_psd):
    return np.exp(np.interp(np.log(np.clip(f, 1e-3, None)),
                            np.log(np.clip(f_psd, 1e-3, None)), np.log(P_psd)))


def amp_cont(f, Mc, DL_):
    return (1.0 / DL_) * math.sqrt(5.0 / 24.0) * math.pi ** (-2.0 / 3.0) \
        * Mc ** (5.0 / 6.0) * f ** (-7.0 / 6.0)


def t_of_f(f, Mc):
    return -(5.0 / 256.0) * Mc ** (-5.0 / 3.0) * (np.pi * f) ** (-8.0 / 3.0)


def f_isco(Mc):
    """ISCO frequency of the total mass inferred from the chirp mass; for
    GW150914 the remnant forms near 250 Hz, and the inspiral-only template
    should not be trusted above it.

    In geometric units (G = c = 1) the total mass in SECONDS is
    M_s = (G M_kg)/c^3, and f_isco = 1/(6 sqrt(6) pi M_s).
    """
    q = 29.0 / 36.0
    M_s = Mc * (1 + q) ** 1.2 / q ** 0.6          # total mass in seconds
    return 1.0 / (6.0 * math.sqrt(6.0) * math.pi * M_s)


class Filter:
    def __init__(self, seg_n, fs, fmin, fmax, Sn, tau2=None, t_ref=0.0):
        self.n = seg_n
        self.fs = fs
        self.f = np.fft.rfftfreq(seg_n, 1.0 / fs)
        self.df = self.f[1] - self.f[0]
        self.m = (self.f >= fmin) & (self.f <= fmax)
        self.fm = self.f[self.m]
        self.Sn = Sn[self.m]
        self.C = 4.0 * self.df / fs ** 2
        psi = (3.0 / 128.0) * (np.pi * MC * self.fm) ** (-5.0 / 3.0) - np.pi / 4.0
        if tau2 is not None:
            t = t_of_f(self.fm, MC) - t_ref
            psi = psi + 2.0 * np.pi * self.fm * (tau2 * t ** 2)
        self.H0 = fs * amp_cont(self.fm, MC, DL) * np.exp(1j * psi)
        self.sigma = math.sqrt(self.C * float(np.sum(np.abs(self.H0) ** 2 / self.Sn)))
        self.win = tukey(seg_n, 0.25)

    def rho(self, data, tc_grid):
        x = data * self.win
        base = self.C * np.fft.rfft(x)[self.m] * np.conj(self.H0) / self.Sn
        out = np.empty(len(tc_grid))
        step = 8192
        for i in range(0, len(tc_grid), step):
            chunk = tc_grid[i:i + step]
            ph = np.exp(2j * np.pi * np.outer(chunk, self.fm))
            out[i:i + step] = np.abs(ph @ base) / self.sigma
        return out

    def rho_max(self, data, coarse=1.0 / 1024.0):
        tc_c = np.arange(0.0, self.n / self.fs, coarse)
        rc = self.rho(data, tc_c)
        k = int(np.argmax(rc))
        lo = max(0.0, tc_c[k] - 2 * coarse)
        hi = min(self.n / self.fs, tc_c[k] + 2 * coarse)
        tc_f = np.arange(lo, hi, 1.0 / self.fs)
        rf = self.rho(data, tc_f)
        j = int(np.argmax(rf))
        return (float(rf[j]), float(tc_f[j])) if rf[j] >= rc[k] else (float(rc[k]), float(tc_c[k]))

    def template_time_series(self, tc):
        H = np.zeros(len(self.f), dtype=complex)
        H[self.m] = self.H0 * np.exp(-2j * np.pi * self.fm * tc)
        return np.fft.irfft(H, self.n)


def main():
    out = {"inputs": {}, "gap_decomposition": {}, "direct_fit": {}, "controls": {},
           "selfchecks": {}, "thresholds": {}}
    seg_n = 4 * FS
    seg_t0 = T_EVENT - 2.0
    i0 = int(round(seg_t0 * FS))

    data = {}
    for det in ("H1", "L1"):
        x, gps0 = load(det)
        data[det] = x

    # ---------------- PSD from off-source, same estimator as h3_real_data.py
    psds = {}
    n_off = 8 * FS
    for det in ("H1", "L1"):
        off = np.concatenate([data[det][:n_off], data[det][-n_off:]])
        psds[det] = median_welch(highpass(off, FS), FS, FS, 0.5)

    # ---------------- per-detector: filter, peak, gap decomposition
    rng = np.random.default_rng(SEED)
    res = {}
    for det in ("H1", "L1"):
        f_psd, P_psd = psds[det]
        Sn = psd_on_grid(np.fft.rfftfreq(seg_n, 1.0 / FS), f_psd, P_psd)
        flt = Filter(seg_n, FS, FMIN, FMAX, Sn)
        seg = highpass(data[det][i0:i0 + seg_n], FS)
        rho_pk, t_pk = flt.rho_max(seg)

        # B1: PSD calibration -- injection with target SNR 20 must come back ~20
        def gen(rng_):
            W = np.fft.rfft(rng_.normal(0.0, 1.0, seg_n)) * np.sqrt(Sn * FS / 2.0)
            return highpass(np.fft.irfft(W, seg_n), FS)
        inj = flt.template_time_series(2.0)
        recs = []
        for _ in range(N_CAL):
            nz = gen(rng)
            r = flt.rho(nz + inj * (20.0 / flt.sigma), np.arange(0.0, 4.0, 2.0 / FS))
            recs.append(float(r[np.argmax(r)]))
        recs = np.array(recs)

        # B2: amplitude scale the data prefers
        scale = rho_pk / flt.sigma

        # B3: band truncation -- how much of sigma lies above the ISCO frequency
        fisco = f_isco(MC)
        m_hi = flt.m & (flt.f > fisco)
        sig_above = math.sqrt(flt.C * float(np.sum(np.abs(flt.H0[m_hi[flt.m]]) ** 2
                                                  / flt.Sn[m_hi[flt.m]])))
        frac_above = (sig_above / flt.sigma) ** 2

        res[det] = {
            "peak_snr": rho_pk,
            "peak_time_in_segment_s": t_pk,
            "coalescence_gps": GPS0 + seg_t0 + t_pk,
            "template_optimal_snr": flt.sigma,
            "amplitude_scale_preferred_by_data": scale,
            "implied_distance_Mpc_if_scale_correct": DL_MPC / scale,
            "B1_injection_recovered_mean": float(recs.mean()),
            "B1_injection_recovered_std": float(recs.std()),
            "B3_fisco_Hz": fisco,
            "B3_fraction_of_sigma2_above_isco": frac_above,
            "B3_sigma_if_truncated_at_isco": math.sqrt(
                flt.C * float(np.sum(np.abs(flt.H0[~m_hi[flt.m]]) ** 2
                                     / flt.Sn[~m_hi[flt.m]]))),
        }
    res["H1_minus_L1_ms"] = 1e3 * (res["H1"]["coalescence_gps"] - res["L1"]["coalescence_gps"])
    out["gap_decomposition"] = res

    # ---------------- DIRECT FIT of tau2 to the real strain
    # profile likelihood: for each tau2, best matched-filter SNR over (tc, phic).
    # The nuisance (Mc) is held at its published value; the degeneracy with Mc is
    # handled separately by the projection in h3_real_data.py.
    t_ref = float(t_of_f(np.array([FMIN]), MC)[0])
    tau2_paper = 1.2e-3 / 0.8444822952069388 ** 2
    grid = np.concatenate([
        np.linspace(-8.0, 8.0, 33) * tau2_paper,
    ])
    direct = {}
    for det in ("H1", "L1"):
        f_psd, P_psd = psds[det]
        Sn = psd_on_grid(np.fft.rfftfreq(seg_n, 1.0 / FS), f_psd, P_psd)
        seg = highpass(data[det][i0:i0 + seg_n], FS)
        tc_coarse = np.arange(0.0, seg_n / FS, 2.0 / FS)
        prof = []
        for t2 in grid:
            flt = Filter(seg_n, FS, FMIN, FMAX, Sn, tau2=t2, t_ref=t_ref)
            prof.append(float(np.max(flt.rho(seg, tc_coarse))))
        prof = np.array(prof)
        kbest = int(np.argmax(prof))
        # 1-sigma: where the profile drops by 1 below its max (for a peak SNR
        # statistic the likelihood-ratio drop of 1 corresponds to 1 sigma).
        # If the profile never drops by 1 anywhere on the grid, the data do not
        # constrain tau2 AT ALL on that range -- and that is a result, not a
        # failure, so it is recorded as such rather than returned as None.
        target = prof[kbest] - 1.0
        left = grid[:kbest][prof[:kbest] < target]
        right = grid[kbest:][prof[kbest:] < target]
        lo = float(left[-1]) if len(left) else None
        hi = float(right[0]) if len(right) else None
        span = float(grid[-1] - grid[0])
        direct[det] = {
            "tau2_paper_per_s2": tau2_paper,
            "tau2_best_fit_per_s2": float(grid[kbest]),
            "tau2_best_fit_amplitude_ms": float(grid[kbest] * 0.8444822952069388 ** 2 * 1e3),
            "profile_peak_snr": float(prof[kbest]),
            "profile_at_zero": float(prof[np.argmin(np.abs(grid))]),
            "profile_at_paper": float(prof[np.argmin(np.abs(grid - tau2_paper))]),
            "tau2_at_1sigma_low": lo,
            "tau2_at_1sigma_high": hi,
            "unconstrained_on_grid": bool(lo is None and hi is None),
            "grid_span_tau2_per_s2": span,
            "grid_span_amplitude_ms": float(span * 0.8444822952069388 ** 2 * 1e3),
            "profile_drop_within_grid": float(prof[kbest] - min(prof)),
            "one_sigma_amplitude_ms": None if (lo is None or hi is None) else
                float((hi - lo) / 2 * 0.8444822952069388 ** 2 * 1e3),
            "grid": grid.tolist(),
            "profile": prof.tolist(),
        }
    out["direct_fit"] = direct

    # ---------------- self-checks
    c = {}
    c["S1_B1_calibration_within_5pct"] = all(
        abs(res[d]["B1_injection_recovered_mean"] - 20.0) < 1.0 for d in ("H1", "L1"))
    c["S2_two_detectors_agree_within_5ms"] = abs(res["H1_minus_L1_ms"]) < 5.0
    c["S3_amplitude_scale_below_1"] = all(
        res[d]["amplitude_scale_preferred_by_data"] < 1.0 for d in ("H1", "L1"))
    c["S4_direct_fit_profile_finite"] = all(
        np.all(np.isfinite(direct[d]["profile"])) for d in ("H1", "L1"))
    out["selfchecks"] = {k: bool(v) for k, v in c.items()}

    T = dict(c)
    # The decisive DIRECT-FIT statements. Because the delay direction is
    # (nearly) degenerate, the expected outcome is that the data do NOT
    # constrain tau2: the profile stays flat and never drops by 1 sigma over a
    # grid spanning +-9.6 ms of amplitude. That is the result; it is checked
    # directly rather than dressed up as an upper limit.
    T["T1_profile_flat_over_wide_grid"] = all(
        direct[d]["unconstrained_on_grid"] for d in ("H1", "L1"))
    T["T2_profile_drop_below_1sigma"] = all(
        direct[d]["profile_drop_within_grid"] < 1.0 for d in ("H1", "L1"))
    T["T3_zero_and_paper_indistinguishable"] = all(
        abs(direct[d]["profile_at_zero"] - direct[d]["profile_at_paper"]) < 0.05
        for d in ("H1", "L1"))
    out["thresholds"] = {k: bool(v) for k, v in T.items()}

    os.makedirs(os.path.join(HERE, "artifacts"), exist_ok=True)
    with open(os.path.join(HERE, "artifacts", "h3_direct_fit.json"), "w") as fh:
        json.dump(out, fh, indent=2, sort_keys=True)
    print(json.dumps(out, indent=2, sort_keys=True))
    return 0 if all(T.values()) else 1


if __name__ == "__main__":
    sys.exit(main())
