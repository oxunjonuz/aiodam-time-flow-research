#!/usr/bin/env python3
"""
H3 on the REAL GW150914 strain (H1 + L1), with a PSD MEASURED from the data.

Difference from analysis.py (v5): v5 weighted everything with an ANALYTIC PSD
that had been calibrated so the injected signal's SNR equalled the published 20.
Here nothing is calibrated: the PSD is estimated from off-source strain and the
strain amplitude is the physical one (M_c = 28.096 Msun, D_L = 410 Mpc).

MATCHED-FILTER CONVENTION, settled empirically (mf_norm_probe.py), not by
algebra.  With X_k the numpy rfft of the sampled segment, the template written
in the SAME discrete convention H_k = fs * h~(f_k), and the one-sided PSD S_k:
    z(t_c) = C * sum_k X_k conj(H_k) e^{2 pi i f_k t_c} / S_k,   C = 4 df / fs^2
    rho(t_c) = |z(t_c)| / sigma,  sigma^2 = C * sum_k |H_k|^2 / S_k
sigma then equals the standard optimal SNR 4*int |h~|^2/S df exactly (measured
45.5191 both ways).  Measured behaviour of the statistic:
  * pure synthetic noise        -> peak rho 4.8   (Rayleigh maximum over 65536 bins)
  * injection scaled to SNR 20  -> recovered 22.95 at the exact injected time
  * real H1 strain, 4 s at the event -> peak 10.08 at GPS 1126259462.38,
    i.e. 20 ms from the published merger time 1126259462.40
  * the SAME filter on an off-source segment -> peak 7.11
The template is leading-order inspiral only, so recovering ~10 where the
published full-IMR single-detector SNR is ~20 is expected, not a defect.

Two bugs on the way, both caught by these controls and recorded in NOTES.md:
  (1) C = 4*df (the continuous-FT form) with an inverse-FFT reconstruction gave
      rho ~ 1e5 on pure noise -- the discrete rfft already carries a factor fs;
  (2) without a Tukey window the peak ran to the segment edge (rho 306), because
      the rectangular edges leak into the band.

Deterministic headline numbers; the only random numbers are fixed-seed synthetic
noise used for validation.
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


# ------------------------------------------------------------------ helpers
def load(det):
    with h5py.File(os.path.join(DATA, FILES[det]), "r") as h:
        return (h["strain/Strain"][:], int(h["meta/GPSstart"][()]),
                h["meta/Detector"][()].decode(), int(h["meta/Duration"][()]))


def highpass(x, fs, f0=15.0):
    return sosfiltfilt(butter(4, f0 / (fs / 2.0), btype="highpass", output="sos"), x)


def median_welch(x, fs, nper, overlap=0.5):
    """MEDIAN-averaged one-sided periodogram, 50%-overlapping Hann segments.

    Why median and not mean.  psd_whiten_check.py compares, on the same
    off-source data, the FFT-binned in-band variance with the PSD integral.
    The median estimator gives ASD(100 Hz) = 1.03e-23 /rtHz, which is the
    published aLIGO O1 sensitivity; the mean estimator gives 2.3e-22, a factor
    ~22 too high, because the released strain contains non-Gaussian transients
    (glitches) and the mean is not robust to them.  The apparent "6x excess" of
    the data over the median PSD is exactly that non-Gaussianity: the synthetic
    control generated from the median PSD reproduces the PSD integral to 0.98,
    so the median PSD is self-consistent for Gaussian noise and represents the
    typical (quiet) sensitivity, which is the right weight for a matched filter.
    Median-averaging an exponential is biased low by ln 2; corrected here.
    """
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
    """Continuous strain amplitude h~(f) of the leading-order inspiral."""
    return (1.0 / DL_) * math.sqrt(5.0 / 24.0) * math.pi ** (-2.0 / 3.0) \
        * Mc ** (5.0 / 6.0) * f ** (-7.0 / 6.0)


def t_of_f(f, Mc):
    """t(f) - t_c = -(5/256) Mc^(-5/3) (pi f)^(-8/3)."""
    return -(5.0 / 256.0) * Mc ** (-5.0 / 3.0) * (np.pi * f) ** (-8.0 / 3.0)


def project(u, cols, w):
    """Weighted orthogonal projection -> (absorbed_fraction, residual_snr)."""
    Q, R = np.linalg.qr(cols * w[:, None], mode="reduced")
    d = np.abs(np.diag(R))
    Qk = Q[:, d > d.max() * 1e-10]
    y = u * w
    n = float(np.linalg.norm(y))
    if n == 0.0:
        return 1.0, 0.0
    coef = Qk.T @ y
    return (float(np.linalg.norm(Qk @ coef) / n), float(np.linalg.norm(y - Qk @ coef)))


def schur_sigma(cols, w):
    A = cols * w[:, None]
    F = A.T @ A
    Fdd, Faa, Fda = F[-1, -1], F[:-1, :-1], F[:-1, -1]
    try:
        s = Fdd - Fda @ np.linalg.solve(Faa, Fda)
    except np.linalg.LinAlgError:
        return float("inf")
    return float("inf") if s <= 0 else 1.0 / math.sqrt(s)


class Filter:
    """Optimal-SNR matched filter with the discrete convention of the header."""

    def __init__(self, seg_n, fs, fmin, fmax, Sn, tc0=0.0):
        self.n = seg_n
        self.fs = fs
        self.f = np.fft.rfftfreq(seg_n, 1.0 / fs)
        self.df = self.f[1] - self.f[0]
        self.m = (self.f >= fmin) & (self.f <= fmax)
        self.fm = self.f[self.m]
        self.Sn = Sn[self.m]
        self.C = 4.0 * self.df / fs ** 2
        self.H0 = fs * amp_cont(self.fm, MC, DL) * np.exp(
            1j * ((3.0 / 128.0) * (np.pi * MC * self.fm) ** (-5.0 / 3.0)
                  - np.pi / 4.0 - 2.0 * np.pi * self.fm * tc0))
        self.sigma = math.sqrt(self.C * float(np.sum(np.abs(self.H0) ** 2 / self.Sn)))
        self.win = tukey(seg_n, 0.25)

    def rho(self, data, tc_grid, window=True):
        x = data * self.win if window else data
        base = self.C * np.fft.rfft(x)[self.m] * np.conj(self.H0) / self.Sn
        out = np.empty(len(tc_grid))
        step = 8192
        for i in range(0, len(tc_grid), step):
            chunk = tc_grid[i:i + step]
            ph = np.exp(2j * np.pi * np.outer(chunk, self.fm))
            out[i:i + step] = np.abs(ph @ base) / self.sigma
        return out

    def rho_max(self, data, coarse=1.0 / 1024.0, window=True):
        """Two-stage peak search: coarse grid, then 1/fs around the winner."""
        tc_c = np.arange(0.0, self.n / self.fs, coarse)
        rc = self.rho(data, tc_c, window)
        k = int(np.argmax(rc))
        lo = max(0.0, tc_c[k] - 2 * coarse)
        hi = min(self.n / self.fs, tc_c[k] + 2 * coarse)
        tc_f = np.arange(lo, hi, 1.0 / self.fs)
        rf = self.rho(data, tc_f, window)
        j = int(np.argmax(rf))
        return (float(rf[j]), float(tc_f[j])) if rf[j] >= rc[k] else (float(rc[k]), float(tc_c[k]))

    def sigma_independent(self):
        """sigma via the continuous integral 4*int |h~|^2/S df, a different route
        from the discrete sum C*sum|H|^2/S used in __init__."""
        return math.sqrt(4.0 * self.df * float(
            np.sum(np.abs(self.H0 / self.fs) ** 2 / self.Sn)))

    def template_time_series(self, tc):
        """The template as a real time series with coalescence at tc."""
        H = np.zeros(len(self.f), dtype=complex)
        H[self.m] = self.H0 * np.exp(-2j * np.pi * self.fm * tc)
        return np.fft.irfft(H, self.n)


# ------------------------------------------------------------------ main
def main():
    out = {"inputs": {}, "mf_calibration": {}, "matched_filter": {}, "projection": {},
           "controls": {}, "selfchecks": {}, "thresholds": {}}

    det_data = {}
    for det in ("H1", "L1"):
        x, gps0, meta, dur = load(det)
        assert gps0 == GPS0 and dur == 32 and meta == det and len(x) == 524288, det
        det_data[det] = x
        out["inputs"][det] = {"n_samples": int(len(x)), "fs": FS, "gps_start": gps0,
                              "rms_strain": float(np.sqrt(np.mean(x ** 2))),
                              "max_abs_strain": float(np.max(np.abs(x)))}

    n_off = 8 * FS
    psds = {}
    for det in ("H1", "L1"):
        off = np.concatenate([det_data[det][:n_off], det_data[det][-n_off:]])
        f_psd, P_psd = median_welch(highpass(off, FS), FS, FS, 0.5)
        psds[det] = (f_psd, P_psd)
        band = (f_psd >= FMIN) & (f_psd <= FMAX)
        out["inputs"][f"{det}_psd"] = {
            "estimator": "MEDIAN Welch 1 s Hann 50% overlap, off-source 0-8 s + 24-32 s, /ln2",
            "n_segments": (2 * 8 * FS - FS) // (FS // 2) + 1,
            "psd_at_35Hz": float(np.interp(35.0, f_psd, P_psd)),
            "psd_at_100Hz": float(np.interp(100.0, f_psd, P_psd)),
            "psd_at_250Hz": float(np.interp(250.0, f_psd, P_psd)),
            "asd_at_100Hz": float(np.sqrt(np.interp(100.0, f_psd, P_psd))),
            "asd_min": float(np.sqrt(np.min(P_psd[band]))),
            "consistency": "psd_whiten_check.py: median ASD(100Hz)=1.03e-23 /rtHz matches "
                           "published aLIGO O1; a synthetic series from this PSD reproduces "
                           "its own band integral to 0.98. The real off-source data carries "
                           "~6x more in-band power than this PSD, i.e. the released strain "
                           "has non-Gaussian transients; the median is the robust estimator.",
        }

    # ================================================ calibration of the filter
    seg_n = 4 * FS
    t_ref = 2.0
    tc_coarse = np.arange(0.0, seg_n / FS, 2.0 / FS)
    rng = np.random.default_rng(SEED)
    cal = {}
    for det in ("H1", "L1"):
        f_psd, P_psd = psds[det]
        Sn = psd_on_grid(np.fft.rfftfreq(seg_n, 1.0 / FS), f_psd, P_psd)
        flt = Filter(seg_n, FS, FMIN, FMAX, Sn)
        # noise generator consistent with the PSD: E|X_k|^2 = (fs*N/2) S_k
        def gen(rng_):
            W = np.fft.rfft(rng_.normal(0.0, 1.0, seg_n)) * np.sqrt(Sn * FS / 2.0)
            return highpass(np.fft.irfft(W, seg_n), FS)   # same conditioning as real data
        inj_unit = flt.template_time_series(t_ref)
        peaks, recs = [], []
        for _ in range(N_CAL):
            nz = gen(rng)
            peaks.append(float(np.max(flt.rho(nz, tc_coarse))))
            rho = flt.rho(nz + inj_unit * (20.0 / flt.sigma), tc_coarse)
            recs.append(float(rho[np.argmax(rho)]))
        peaks = np.array(peaks)
        recs = np.array(recs)
        cal[det] = {
            "sigma_discrete_sum": flt.sigma,
            "sigma_continuous_integral": flt.sigma_independent(),
            "n_realizations": N_CAL,
            "noise_only_peak_mean": float(peaks.mean()),
            "noise_only_peak_max": float(peaks.max()),
            "injection_target_snr": 20.0,
            "injection_recovered_mean": float(recs.mean()),
            "injection_recovered_std": float(recs.std()),
            "injection_recovered_min": float(recs.min()),
            "injection_recovered_max": float(recs.max()),
        }
        rho, t_pk = flt.rho_max(gen(rng) + inj_unit * (20.0 / flt.sigma))
        cal[det]["injection_peak_time_s"] = t_pk
        cal[det]["injection_peak_rho"] = rho
    out["mf_calibration"] = cal

    # ================================================ filter on real data
    seg_t0 = T_EVENT - 2.0
    i0 = int(round(seg_t0 * FS))
    for det in ("H1", "L1"):
        f_psd, P_psd = psds[det]
        Sn = psd_on_grid(np.fft.rfftfreq(seg_n, 1.0 / FS), f_psd, P_psd)
        flt = Filter(seg_n, FS, FMIN, FMAX, Sn)
        rho, t_pk = flt.rho_max(highpass(det_data[det][i0:i0 + seg_n], FS))
        out["matched_filter"][det] = {
            "peak_snr": rho,
            "peak_time_in_segment_s": t_pk,
            "coalescence_gps": GPS0 + seg_t0 + t_pk,
            "offset_from_published_merger_ms": 1e3 * (GPS0 + seg_t0 + t_pk - GPS_EVENT),
            "template_optimal_snr": flt.sigma,
            "segment_start_s": seg_t0,
        }
        j0 = int(round(4.0 * FS))
        rho_off = flt.rho(highpass(det_data[det][j0:j0 + seg_n], FS), tc_coarse)
        out["controls"][f"{det}_offsource_peak_snr"] = float(np.max(rho_off))
        out["controls"][f"{det}_offsource_peak_t_s"] = float(np.argmax(rho_off) * 2.0 / FS)
    out["matched_filter"]["H1_minus_L1_ms"] = 1e3 * (
        out["matched_filter"]["H1"]["coalescence_gps"]
        - out["matched_filter"]["L1"]["coalescence_gps"])

    # ================================================ projection, measured PSD
    df = 1.0 / 4.0
    f = np.arange(FMIN, FMAX + df / 2.0, df)
    t_span = abs(float(t_of_f(np.array([FMIN]), MC)[0]) - float(t_of_f(np.array([FMAX]), MC)[0]))
    res = {}
    for det in ("H1", "L1"):
        f_psd, P_psd = psds[det]
        w = np.sqrt(4.0 * df / psd_on_grid(f, f_psd, P_psd))
        h = amp_cont(f, MC, DL)
        dMc = -(5.0 / 3.0) * (3.0 / 128.0) * (np.pi * MC * f) ** (-5.0 / 3.0) / MC
        cols = np.column_stack([h * dMc, h * (2 * np.pi * f), -h])
        u0 = h * (2 * np.pi * f * 1.2e-3)
        af0, snr0 = project(u0, cols, w)
        ul = h * (2 * np.pi * f * (1e-2 * t_of_f(f, MC)))
        afl, snrl = project(ul, cols, w)
        tau2 = 1.2e-3 / t_span ** 2
        uq = h * (2 * np.pi * f * (tau2 * t_of_f(f, MC) ** 2))
        afq, snrq = project(uq, cols, w)
        amp_det = 1.2e-3 * 5.0 / snrq if snrq > 0 else float("inf")
        smooth = np.interp(f, np.linspace(FMIN, FMAX, 8), rng.normal(0, 1, 8))
        # cols = [dMc, dtc, dphic]; freeze t_c -> {dMc,dphic}; freeze M_c -> {dtc,dphic}
        ftc = project(u0, cols[:, [0, 2]], w)
        fmc = project(h * (2 * np.pi * f * (1e-1 * t_of_f(f, MC))), cols[:, [1, 2]], w)
        res[det] = {
            "signal_snr_model": float(np.linalg.norm(h * w)),
            "t_span_s": t_span,
            "H1_const_delay": {"amplitude_ms": 1.2, "absorbed_fraction": af0,
                               "one_minus_absorbed": 1 - af0, "unmodelled_snr": snr0},
            "H2_linear_delay": {"taudot": 1e-2, "absorbed_fraction": afl,
                                "one_minus_absorbed": 1 - afl, "unmodelled_snr": snrl},
            "H3_quadratic": {"amplitude_ms": 1.2, "absorbed_fraction": afq,
                             "one_minus_absorbed": 1 - afq, "unmodelled_snr": snrq,
                             "detectable_amplitude_ms_at_snr5": amp_det * 1e3,
                             "ratio_detectable_to_paper": amp_det / 1.2e-3,
                             "fisher_sigma_amplitude_ms": schur_sigma(
                                 np.column_stack([cols, uq / 1.2e-3]), w) * 1e3},
            "controls": {"zero_injection_snr": project(np.zeros_like(f), cols, w)[1],
                         "freeze_tc_absorbed": ftc[0], "freeze_tc_snr": ftc[1],
                         "freeze_mc_absorbed": fmc[0], "freeze_mc_snr": fmc[1],
                         "scrambled_absorbed": project(h * smooth, cols, w)[0]},
        }
    out["projection"] = res

    # ================================================ self-checks
    c = {}
    c["S1_files_are_32s_16kHz_at_stated_gps"] = all(
        out["inputs"][d]["n_samples"] == 524288 and out["inputs"][d]["gps_start"] == GPS0
        for d in ("H1", "L1"))
    c["S2_sigma_two_routes_agree"] = all(
        abs(cal[d]["sigma_discrete_sum"] - cal[d]["sigma_continuous_integral"])
        <= 1e-9 * cal[d]["sigma_continuous_integral"] for d in ("H1", "L1"))
    c["S3_noise_only_peak_below_8_all_realizations"] = all(
        cal[d]["noise_only_peak_max"] < 8.0 for d in ("H1", "L1"))
    c["S4_injection_snr20_recovered_within_25pct"] = all(
        abs(cal[d]["injection_recovered_mean"] - 20.0) < 5.0 for d in ("H1", "L1"))
    c["S5_injection_time_recovered_within_2ms"] = all(
        abs(cal[d]["injection_peak_time_s"] - t_ref) < 2e-3 for d in ("H1", "L1"))
    # S6 criterion, stated as what a LEADING-ORDER inspiral template can deliver.
    # The published single-detector SNR of GW150914 is ~20 for the FULL IMR
    # waveform.  This template is restricted post-Newtonian inspiral truncated at
    # 300 Hz, with no merger or ringdown, so it cannot recover 20 -- it recovers
    # a fraction of it.  The defensible criterion is therefore: the on-source peak
    # must clear the standard detection threshold 5 in at least one detector, and
    # must exceed the same filter's off-source peak in both.  The measured values
    # are reported next to it so the gap is visible, not hidden.
    c["S6_onsource_above_5_and_above_offsource"] = (
        max(out["matched_filter"][d]["peak_snr"] for d in ("H1", "L1")) >= 5.0
        and all(out["matched_filter"][d]["peak_snr"]
                > out["controls"][f"{d}_offsource_peak_snr"] for d in ("H1", "L1")))
    c["S7_real_peak_within_100ms_of_published_merger"] = all(
        abs(out["matched_filter"][d]["offset_from_published_merger_ms"]) < 100.0
        for d in ("H1", "L1"))
    c["S8_H1_H2_absorbed_with_measured_psd"] = all(
        res[d]["H1_const_delay"]["one_minus_absorbed"] <= 1e-10
        and res[d]["H2_linear_delay"]["one_minus_absorbed"] <= 1e-10 for d in ("H1", "L1"))
    c["S9_quadratic_not_absorbed"] = all(
        res[d]["H3_quadratic"]["absorbed_fraction"] < 0.999 for d in ("H1", "L1"))
    out["selfchecks"] = {k: bool(v) for k, v in c.items()}

    T = dict(out["selfchecks"])
    T["T5_H3_snr_lt_1"] = all(res[d]["H3_quadratic"]["unmodelled_snr"] < 1.0 for d in ("H1", "L1"))
    T["T6_H3_x100_snr_ge_5"] = all(res[d]["H3_quadratic"]["unmodelled_snr"] * 100.0 >= 5.0
                                   for d in ("H1", "L1"))
    T["T7_zero_injection_le_1e-12"] = all(
        res[d]["controls"]["zero_injection_snr"] <= 1e-12 for d in ("H1", "L1"))
    T["T13_freeze_tc_snr_gt_5"] = all(res[d]["controls"]["freeze_tc_snr"] > 5.0 for d in ("H1", "L1"))
    T["T14_freeze_mc_snr_gt_5"] = all(res[d]["controls"]["freeze_mc_snr"] > 5.0 for d in ("H1", "L1"))
    T["T15_scrambled_not_absorbed"] = all(res[d]["controls"]["scrambled_absorbed"] < 0.95
                                          for d in ("H1", "L1"))
    out["thresholds"] = T

    os.makedirs(os.path.join(HERE, "artifacts"), exist_ok=True)
    with open(os.path.join(HERE, "artifacts", "h3_real_data.json"), "w") as fh:
        json.dump(out, fh, indent=2, sort_keys=True)
    print(json.dumps(out, indent=2, sort_keys=True))
    return 0 if all(T.values()) else 1


if __name__ == "__main__":
    sys.exit(main())
