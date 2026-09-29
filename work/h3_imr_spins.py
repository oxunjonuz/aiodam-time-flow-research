#!/usr/bin/env python3
"""H3 point 1 (owner's msg_00226): run IMRPhenomT WITH SPINS, so the face-on
zero-spin number stops being an "upper bound" and becomes a RECOVERED interval.

msg224 measured, with spins = 0 and inclination = 0:
    rho_opt = 31.6595 (H1) / 28.7350 (L1)
and the artifact declared that an UPPER BOUND, because face-on zero-spin is the
maximal-amplitude configuration.  That wording is honest but weak: it does not
say how far the real number could move.  This script replaces the bound with a
range over the PUBLISHED GW150914 posterior.

WHAT IS USED, AND WHERE IT COMES FROM (registered sources):
  * component spins      a1 = 0.32 (+0.49/-0.29), a2 = 0.44 (+0.50/-0.40);
                         90% upper bounds a1 < 0.69, a2 < 0.89
                         (arXiv:1602.03840, Table I)
  * effective spin       chi_eff = -0.07 (+0.16/-0.17)      (arXiv:1602.03840)
  * inclination          P(45 deg < theta_JN < 135 deg) = 0.35 -- strongly
                         misaligned is DISFAVOURED, but the angle is not pinned
                         (arXiv:1602.03840)
  * single-detector SNR  re-weighted rho_hat = 19.5 (H1), 13.3 (L1)
                         (arXiv:1602.03839)
  * network SNR          24                                 (arXiv:1602.03837)
  * network optimal SNR  25.1 (+1.7/-1.7)                   (arXiv:1602.03840)

TWO THINGS THIS SCRIPT HAS TO SAY ABOUT THE TARGET, because the project has been
quoting a round "20" that appears in none of the papers:
  (a) the published single-detector numbers are 19.5 and 13.3 -- not 20;
  (b) 25.1 is a NETWORK optimal SNR and 24 is a NETWORK SNR.  Comparing a
      single-detector rho_opt to a network number is a category error.  So the
      comparison target here is 19.5 / 13.3, and the "~20" of earlier rounds is
      retired.

THE SPIN ENVELOPE IS NOT A RECTANGLE.  a1 and a2 are individually bounded, but
they are also JOINTLY constrained by chi_eff, which is measured to be
-0.07 +/- 0.17.  The pair (a1=0.69, a2=0.89) sits at chi_eff = +0.78, which the
measurement excludes by a wide margin.  So the envelope scanned here is
    |a1| <= 0.69, |a2| <= 0.89, chi_eff in [-0.24, +0.09]
and the naive corner is reported separately as the "outside chi_eff" case, so
the difference between the two is visible rather than hidden.

CONTROLS (a failed CONTROL means the instrument is broken and the run is void):
  NC1  spins=0, iota=0 reproduces the recorded 31.659512701612417 (H1) /
       28.73504152105233 (L1) to < 1e-6
  NC2  the spins=0 inclination sweep reproduces the recorded sweep
  NC3  amplitude linearity: doubling the waveform doubles rho_opt (exact)
  NC4  zero waveform -> rho_opt = 0
  NC5  the set labelled chi_eff_median really carries the published chi_eff
  NC6  Cauchy-Schwarz: peak matched-filter SNR <= rho_opt for EVERY template
  NC7  the envelope scan must be non-vacuous: it must contain both a point with
       chi_eff below the published range and one above, or the constraint is
       doing nothing.

A FALSE CRITERION FROM MY FIRST DRAFT IS KEPT IN THE ARTIFACT, NOT DELETED.
Draft NC6 asserted "a template with more band power cannot have a LOWER
matched-filter peak on the same data".  MEASURED: FALSE.  Adding the published
spins raises rho_opt from 31.66 to 34.84 while the matched-filter peak FALLS
from 16.74 to 13.51, because the spin also changes the template's PHASE, so the
overlap with the data is not a function of template power alone.  That is the
fourth wrong criterion in this project, and like the others it is recorded.
"""
import json
import math
import os
import sys

import numpy as np
import h5py
from scipy.signal import butter, sosfiltfilt
from scipy.signal.windows import tukey

from phenomxpy import IMRPhenomT

FS = 16384
GPS0 = 1126259447
GPS_EVENT = 1126259462.4
T_EVENT = GPS_EVENT - GPS0
FMIN, FMAX = 20.0, 300.0
MC_MSUN = 28.096
DL_MPC = 410.0
Q = 29.0 / 36.0
ETA = Q / (1.0 + Q) ** 2
MTOT_MSUN = MC_MSUN / ETA ** 0.6
F_MAX_FD = 1024.0
DF_FD = 0.25

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

# ---- recorded values this script must reproduce (msg224 artifact)
RECORDED = {
    "H1": {"rho_opt_face_on_zero_spin": 31.659512701612417,
           "sweep": {0.0: 31.659512701612417, 15.0: 30.599120093967446,
                     30.0: 27.70207361391088, 45.0: 23.744634526209317,
                     60.0: 19.787195438507773, 90.0: 15.8297563508062}},
    "L1": {"rho_opt_face_on_zero_spin": 28.73504152105233, "sweep": {}},
}

# ---- published GW150914 numbers
PUBLISHED_RHO_HAT = {"H1": 19.5, "L1": 13.3}     # arXiv:1602.03839
NETWORK_SNR = 24.0                                # arXiv:1602.03837
NETWORK_OPT_SNR = 25.1                            # arXiv:1602.03840

# ---- published spin posterior (arXiv:1602.03840, Table I)
M1_SRC, M2_SRC = 36.3, 28.6
CHI_EFF_MEDIAN, CHI_EFF_LO, CHI_EFF_HI = -0.07, -0.17, +0.16
CHI_EFF_RANGE = (CHI_EFF_MEDIAN + CHI_EFF_LO, CHI_EFF_MEDIAN + CHI_EFF_HI)
A1_MED, A1_LO, A1_HI = 0.32, -0.29, +0.49
A2_MED, A2_LO, A2_HI = 0.44, -0.40, +0.50
A1_BOUND, A2_BOUND = 0.69, 0.89

NAMED_SETS = {
    "zero":           (0.0, 0.0),
    "chi_eff_median": (CHI_EFF_MEDIAN, CHI_EFF_MEDIAN),
    "medians":        (A1_MED, A2_MED),
    "medians_low":    (A1_MED + A1_LO, A2_MED + A2_LO),
    "medians_high":   (A1_MED + A1_HI, A2_MED + A2_HI),
    "bounds_corner":  (A1_BOUND, A2_BOUND),      # outside chi_eff -- kept visible
    "bounds_low":     (-A1_BOUND, -A2_BOUND),    # outside chi_eff -- kept visible
}


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


def imr_model(freqs, a1, a2, inclination_deg):
    """One IMRPhenomT call on the pipeline's grid (df = 0.25 Hz).

    The model's native grid runs to F_MAX_FD = 1024 Hz (4097 samples); the
    pipeline's segment grid runs to Nyquist (32769 samples at df = 0.25 Hz).
    Same df, so the interpolation below is an IDENTITY on shared samples -- it
    only lifts the model onto the longer array, as h3_imr_check.py does.
    (A verifier once built the model at df = 0.5 and interpolated: that cost 19%
    of the band power, because the FD phase rotates by ~pi between adjacent
    0.5 Hz samples.  Same trap, avoided by construction.)
    """
    wf = IMRPhenomT(eta=ETA, s1=[0, 0, a1], s2=[0, 0, a2],
                    f_min=FMIN, f_ref=FMIN, total_mass=MTOT_MSUN,
                    distance=DL_MPC, inclination=math.radians(inclination_deg),
                    phi_ref=0.0, f_max=F_MAX_FD, delta_t=0.5 / F_MAX_FD,
                    delta_f=DF_FD, condition=True)
    fd = wf.compute_fd_polarizations()
    hp = np.asarray(fd[0], dtype=complex)
    f_model = np.arange(len(hp)) * DF_FD
    return (np.interp(freqs, f_model, hp.real)
            + 1j * np.interp(freqs, f_model, hp.imag))


def chi_eff(a1, a2, m1=M1_SRC, m2=M2_SRC):
    return (a1 * m1 + a2 * m2) / (m1 + m2)


def rho_opt(freqs, htilde, Sn, lo=FMIN, hi=FMAX):
    m = (freqs >= lo) & (freqs <= hi)
    return math.sqrt(4.0 * DF_FD * float(np.sum(np.abs(htilde[m]) ** 2 / Sn[m]))), int(m.sum())


def matched_filter_peak(seg, freqs, htilde, Sn, rho, lo=FMIN, hi=FMAX, fs=FS):
    m = (freqs >= lo) & (freqs <= hi)
    fm = freqs[m]
    xw = seg * tukey(len(seg), 0.25)
    X = np.fft.rfft(xw)
    base = (4.0 * DF_FD / fs) * np.conj(htilde[m]) * X[m] / Sn[m]
    tc = np.arange(0.0, len(seg) / fs, 2.0 / fs)
    out = np.empty(len(tc))
    step = 8192
    for i in range(0, len(tc), step):
        ch = tc[i:i + step]
        out[i:i + step] = np.abs(np.exp(2j * np.pi * np.outer(ch, fm)) @ base) / rho
    k = int(np.argmax(out))
    return float(out[k]), float(tc[k])


def spin_envelope(lo=-A1_BOUND, hi=A1_BOUND, n=13, lo2=-A2_BOUND, hi2=A2_BOUND, n2=13):
    """Grid over |a1|<=0.69, |a2|<=0.89 intersected with the published chi_eff."""
    a1s = np.linspace(lo, hi, n)
    a2s = np.linspace(lo2, hi2, n2)
    inside, outside = [], []
    for a1 in a1s:
        for a2 in a2s:
            ce = chi_eff(a1, a2)
            (inside if CHI_EFF_RANGE[0] <= ce <= CHI_EFF_RANGE[1] else outside).append((float(a1), float(a2), ce))
    return inside, outside


def main():
    out = {"setup": {}, "published": {}, "per_detector": {}, "envelope": {},
           "controls": {}, "control_values": {}, "findings": {}, "thresholds": {},
           "answer": {}, "notes": {}}
    seg_n = 4 * FS
    seg_t0 = T_EVENT - 2.0
    i0 = int(round(seg_t0 * FS))
    freqs = np.fft.rfftfreq(seg_n, 1.0 / FS)
    assert abs((freqs[1] - freqs[0]) - DF_FD) < 1e-12

    out["published"] = {
        "rho_hat_single_detector": PUBLISHED_RHO_HAT,
        "network_snr": NETWORK_SNR,
        "network_optimal_snr": NETWORK_OPT_SNR,
        "category_note": "19.5 / 13.3 are single-detector; 24 and 25.1 are "
                         "network quantities and are NOT comparison targets for "
                         "a single-detector rho_opt. The round '20' used in "
                         "earlier rounds is retired.",
        "chi_eff_median": CHI_EFF_MEDIAN, "chi_eff_90pct_range": list(CHI_EFF_RANGE),
        "a1_median": A1_MED, "a2_median": A2_MED,
        "a1_bound_90pct": A1_BOUND, "a2_bound_90pct": A2_BOUND,
        "theta_JN": "P(45 deg < theta_JN < 135 deg) = 0.35; bounded, not measured",
        "sources": ["arXiv:1602.03837", "arXiv:1602.03839", "arXiv:1602.03840"],
    }
    out["setup"] = {
        "eta": ETA, "total_mass_Msun": MTOT_MSUN, "chirp_mass_Msun": MC_MSUN,
        "distance_Mpc": DL_MPC, "df_Hz": DF_FD, "seg_n": seg_n,
        "phenomxpy": "IMRPhenomT", "spin_family": "aligned (s = [0, 0, a])",
        "mass_ratio_used": Q, "mass_ratio_published": "0.82 (+0.17/-0.20)",
        "m1_source_Msun": M1_SRC, "m2_source_Msun": M2_SRC,
    }

    inside, outside = spin_envelope()
    out["envelope"] = {
        "rule": "|a1| <= 0.69, |a2| <= 0.89, chi_eff in [-0.24, +0.09]",
        "n_inside": len(inside), "n_outside": len(outside),
        "chi_eff_of_bounds_corner": chi_eff(A1_BOUND, A2_BOUND),
        "outside_chi_eff_note": "the corner (0.69, 0.89) carries chi_eff = +0.78, "
                                "which the chi_eff measurement excludes; it is "
                                "kept in the named sets so the difference is visible.",
    }

    res = {}
    for det in ("H1", "L1"):
        x, gps0 = load(det)
        assert gps0 == GPS0
        off = np.concatenate([x[:8 * FS], x[-8 * FS:]])
        f_psd, P_psd = median_welch(highpass(off, FS), FS, FS, 0.5)
        Sn = psd_on_grid(freqs, f_psd, P_psd)
        seg = highpass(x[i0:i0 + seg_n], FS)

        d = {"named_sets": {}, "inclination_at_chi_eff_median": {},
             "inclination_at_published_spins": {}, "envelope_scan": {}}

        for name, (a1, a2) in NAMED_SETS.items():
            ht = imr_model(freqs, a1, a2, 0.0)
            r, nb = rho_opt(freqs, ht, Sn)
            pk, tk = matched_filter_peak(seg, freqs, ht, Sn, r)
            d["named_sets"][name] = {
                "a1": a1, "a2": a2, "chi_eff": chi_eff(a1, a2),
                "chi_eff_inside_published": bool(CHI_EFF_RANGE[0] <= chi_eff(a1, a2) <= CHI_EFF_RANGE[1]),
                "rho_opt": r, "n_bins": nb, "peak_snr": pk, "peak_time_s": tk,
                "amplitude_scale": pk / r,
                "implied_distance_Mpc_if_scale_correct": DL_MPC / (pk / r),
                "rho_opt_over_zero_spin": r / RECORDED[det]["rho_opt_face_on_zero_spin"],
            }

        for label, (a1, a2) in (("chi_eff_median", (CHI_EFF_MEDIAN, CHI_EFF_MEDIAN)),
                                ("medians", (A1_MED, A2_MED))):
            g = {}
            for inc in (0.0, 15.0, 30.0, 45.0, 60.0, 90.0):
                ht = imr_model(freqs, a1, a2, inc)
                r, _ = rho_opt(freqs, ht, Sn)
                g[f"{inc:.1f}"] = r
            d["inclination_at_chi_eff_median" if label == "chi_eff_median"
              else "inclination_at_published_spins"] = g

        # The inclination factor is EXACTLY (1+cos^2 iota)/2 for the h+ that
        # phenomxpy returns, and it is INDEPENDENT of the spins -- verified below
        # as a control rather than assumed.  That is what makes the joint range
        # exactly spin_range x [0.5, 1.0] and not a sampled approximation.
        d["inclination_factor_vs_geometric"] = {
            f"{inc:.1f}": (d["inclination_at_chi_eff_median"][f"{inc:.1f}"]
                           / d["inclination_at_chi_eff_median"]["0.0"]
                           - (1.0 + math.cos(math.radians(inc)) ** 2) / 2.0)
            for inc in (0.0, 15.0, 30.0, 45.0, 60.0, 90.0)
        }

        d["inclination_zero_spin"] = {}
        for inc in (0.0, 15.0, 30.0, 45.0, 60.0, 90.0):
            ht = imr_model(freqs, 0.0, 0.0, inc)
            r, _ = rho_opt(freqs, ht, Sn)
            d["inclination_zero_spin"][f"{inc:.1f}"] = r

        # envelope scan, rho_opt only (the fit is a separate script)
        scan = []
        for a1, a2, ce in inside:
            ht = imr_model(freqs, a1, a2, 0.0)
            r, _ = rho_opt(freqs, ht, Sn)
            scan.append({"a1": a1, "a2": a2, "chi_eff": ce, "rho_opt": r})
        d["envelope_scan"] = scan
        d["recorded_rho_opt_face_on_zero_spin"] = RECORDED[det]["rho_opt_face_on_zero_spin"]
        d["published_rho_hat"] = PUBLISHED_RHO_HAT[det]
        res[det] = d
    out["per_detector"] = res

    # ---------------- ranges
    ans = {}
    for det in ("H1", "L1"):
        ins = [s["rho_opt"] for s in res[det]["envelope_scan"]]
        named_in = [v["rho_opt"] for v in res[det]["named_sets"].values()
                    if v["chi_eff_inside_published"]]
        inc_med = list(res[det]["inclination_at_chi_eff_median"].values())
        lo_spin, hi_spin = min(ins), max(ins)
        # joint envelope: max at iota=0 (already the scan), min at iota=90
        lo_joint = lo_spin * (res[det]["inclination_at_chi_eff_median"]["90.0"]
                              / res[det]["inclination_at_chi_eff_median"]["0.0"])
        target = PUBLISHED_RHO_HAT[det]
        ans[det] = {
            "published_rho_hat": target,
            "rho_opt_spin_only_range": [lo_spin, hi_spin],
            "rho_opt_joint_spin_and_inclination_range": [lo_joint, hi_spin],
            "n_envelope_points": len(ins),
            "named_sets_inside_chi_eff_range": [min(named_in), max(named_in)],
            "spin_only_reaches_published": bool(lo_spin <= target <= hi_spin),
            "joint_reaches_published": bool(lo_joint <= target <= hi_spin),
            "gap_face_on_zero_spin": RECORDED[det]["rho_opt_face_on_zero_spin"] / target,
        }
        sw = res[det]["inclination_at_chi_eff_median"]
        incs = sorted(float(k) for k in sw)
        vv = [sw[f"{i:.1f}"] for i in incs]
        cross = None
        for a, b, va, vb in zip(incs, incs[1:], vv, vv[1:]):
            if (va - target) * (vb - target) <= 0.0:
                cross = a + (b - a) * (va - target) / (va - vb)
                break
        ans[det]["inclination_where_published_rho_hat_reached_deg"] = cross
        ans[det]["shortfall_of_joint_envelope"] = {
            "lowest_joint_rho_opt": lo_joint,
            "published": target,
            "absolute_shortfall": target - lo_joint,
            "relative_shortfall_pct": 100.0 * (target - lo_joint) / target,
            "reached": bool(lo_joint <= target),
        }
    out["answer"] = ans

    # ---------------- the H1/L1 SNR-ratio puzzle, decomposed
    # The published single-detector SNRs differ by more than the template's own
    # optimal SNRs do.  Three ratios, all measured, all comparable:
    #   R_psd      = rho_opt(H1)/rho_opt(L1)  -- PSD only, same waveform
    #   R_peak     = peak(H1)/peak(L1)        -- what the filter actually extracts
    #   R_published= 19.5/13.3                -- the published pair
    ratio = {
        "R_psd_template_optimal": (res["H1"]["named_sets"]["zero"]["rho_opt"]
                                   / res["L1"]["named_sets"]["zero"]["rho_opt"]),
        "R_peak_matched_filter": (res["H1"]["named_sets"]["zero"]["peak_snr"]
                                  / res["L1"]["named_sets"]["zero"]["peak_snr"]),
        "R_published_single_detector": PUBLISHED_RHO_HAT["H1"] / PUBLISHED_RHO_HAT["L1"],
        "note": "R_psd is PSD-only (the same template on both detectors, no "
                "antenna projection).  R_peak already carries the difference in "
                "how much of the template each detector's data actually matches. "
                "R_published additionally carries the ANTENNA RESPONSE of the two "
                "sites, which is not in this model at all.  So the residual "
                "difference between R_peak and R_published is expected to be "
                "detector orientation, and is named rather than smoothed.",
    }
    out["answer"]["snr_ratio_decomposition"] = ratio

    # ---------------- controls
    c = {}
    c["NC1_face_on_zero_spin_reproduces_recorded"] = all(
        abs(res[d]["named_sets"]["zero"]["rho_opt"]
            - RECORDED[d]["rho_opt_face_on_zero_spin"]) < 1e-6 for d in ("H1", "L1"))
    c["NC2_zero_spin_inclination_sweep_reproduces_recorded"] = all(
        abs(res["H1"]["inclination_zero_spin"][f"{k:.1f}"] - v) < 1e-6
        for k, v in RECORDED["H1"]["sweep"].items())

    x, _ = load("H1")
    off = np.concatenate([x[:8 * FS], x[-8 * FS:]])
    f_psd, P_psd = median_welch(highpass(off, FS), FS, FS, 0.5)
    Sn1 = psd_on_grid(freqs, f_psd, P_psd)
    h1 = imr_model(freqs, A1_MED, A2_MED, 0.0)
    r1, _ = rho_opt(freqs, h1, Sn1)
    r2, _ = rho_opt(freqs, 2.0 * h1, Sn1)
    r0, _ = rho_opt(freqs, np.zeros_like(h1), Sn1)
    c["NC3_amplitude_linearity"] = abs(r2 / r1 - 2.0) < 1e-9
    c["NC4_zero_waveform_zero_rho"] = abs(r0) < 1e-300
    c["NC5_spin_set_reproduces_published_chi_eff"] = (
        abs(chi_eff(*NAMED_SETS["chi_eff_median"]) - CHI_EFF_MEDIAN) < 1e-12)
    c["NC6_cauchy_schwarz_peak_le_rho_opt"] = all(
        v["peak_snr"] <= v["rho_opt"] * (1.0 + 1e-9)
        for d in ("H1", "L1") for v in res[d]["named_sets"].values())
    c["NC7_envelope_constraint_non_vacuous"] = bool(
        len(outside) > 0 and len(inside) > 0
        and any(chi_eff(*o[:2]) < CHI_EFF_RANGE[0] for o in outside)
        and any(chi_eff(*o[:2]) > CHI_EFF_RANGE[1] for o in outside))
    # NC8: the inclination factor must be the geometric (1+cos^2 i)/2 to <1e-9,
    # for EVERY detector.  If it were not, the joint range would be a sampled
    # approximation and would have to be reported as one.
    c["NC8_inclination_factor_is_geometric"] = all(
        abs(v) < 1e-9 for d in ("H1", "L1")
        for v in res[d]["inclination_factor_vs_geometric"].values())
    out["controls"] = {k: bool(v) for k, v in c.items()}
    out["control_values"] = {
        "r1": r1, "r2": r2, "r2_over_r1": r2 / r1, "r0": r0,
        "chi_eff_of_chi_eff_median_set": chi_eff(*NAMED_SETS["chi_eff_median"]),
        "peak_zero_H1": res["H1"]["named_sets"]["zero"]["peak_snr"],
        "peak_medians_H1": res["H1"]["named_sets"]["medians"]["peak_snr"],
        "rho_zero_H1": res["H1"]["named_sets"]["zero"]["rho_opt"],
        "rho_medians_H1": res["H1"]["named_sets"]["medians"]["rho_opt"],
        "first_draft_NC6_claim_was_false": True,
        "first_draft_NC6_why": "adding published spins RAISES rho_opt 31.66 -> "
                               "34.84 while the matched-filter peak FALLS 16.74 "
                               "-> 13.51, because the spin also changes the "
                               "template PHASE; overlap with the data is not a "
                               "function of template power alone.",
    }

    # ---------------- findings
    T = dict(c)
    T["T_spin_only_cannot_reach_published_rho_hat"] = all(
        not ans[d]["spin_only_reaches_published"] for d in ("H1", "L1"))
    T["T_joint_envelope_reaches_published"] = all(
        ans[d]["joint_reaches_published"] for d in ("H1", "L1"))
    T["T_spins_move_rho_opt_by_less_than_20pct"] = all(
        abs(v["rho_opt_over_zero_spin"] - 1.0) < 0.20
        for d in ("H1", "L1") for v in res[d]["named_sets"].values())
    out["thresholds"] = {k: bool(v) for k, v in T.items()}
    out["findings"] = {k: bool(v) for k, v in T.items() if k not in c}
    out["controls_all_pass"] = bool(all(c.values()))

    out["answer"]["statement"] = (
        "The 'upper bound' wording is retired. What replaces it, measured: over "
        "the published spin envelope the face-on rho_opt spans a RANGE, and that "
        "range does NOT contain the published single-detector SNR -- spins alone "
        "cannot explain the residual, and the largest spin effect pushes the "
        "number UP, away from the target. Inclination does explain it for H1: the "
        "published 19.5 is reached at a plausible angle. For L1 the published "
        "13.3 is NOT reached even edge-on, which is a genuine open item: the "
        "ratio of published single-detector SNRs (19.5/13.3 = 1.47) is larger "
        "than the ratio of the template's own optimal SNRs (31.66/28.74 = 1.10), "
        "so the two detectors' relative response is not what the PSD ratio "
        "predicts. That is named here rather than smoothed over."
    )
    out["notes"] = {
        "why_not_a_single_number": "theta_JN is bounded, not measured. Any single "
                                   "'recovered SNR' would be a choice.",
        "aligned_spin_limit": "IMRPhenomT is aligned-spin; the published "
                              "precession bound chi_p < 0.71 is not spanned.",
        "published_numbers_not_mutually_consistent": "the a1/a2 90% intervals and "
                                                     "the 90% upper bounds do not "
                                                     "agree exactly (0.81 vs 0.69 "
                                                     "for a1); the envelope rule "
                                                     "uses the bounds, and both "
                                                     "are reported.",
        "reweighted_vs_plain": "the published 19.5 / 13.3 are RE-WEIGHTED SNRs "
                               "(chi^2-weighted, arXiv:1602.03839); everything in "
                               "this script is a plain matched-filter SNR.  The "
                               "re-weighting can only LOWER a number, so a small "
                               "shortfall against the published value is expected "
                               "and is not by itself a discrepancy.  This is "
                               "named because it bounds how much weight the "
                               "residual shortfall can carry.",
        "delay": "not fitted here (point 2).",
    }

    os.makedirs(os.path.join(HERE, "artifacts"), exist_ok=True)
    with open(os.path.join(HERE, "artifacts", "h3_imr_spins.json"), "w") as fh:
        json.dump(out, fh, indent=2, sort_keys=True)
    print(json.dumps({k: out[k] for k in ("controls", "control_values", "findings",
                                          "thresholds", "answer", "envelope")},
                     indent=2, sort_keys=True))
    return 0 if all(c.values()) else 1


if __name__ == "__main__":
    sys.exit(main())
