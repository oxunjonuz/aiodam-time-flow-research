#!/usr/bin/env python3
"""Step (a) continued: at what detector sensitivity does the quadratic delay
become visible?  Scaling estimate for ET / CE.

METHOD, AND WHY IT IS A SCALING ESTIMATE AND NOT A SIMULATION.
A real ET-D / CE sensitivity curve is not available inside this container (the
ET site and the LALSuite data files are unreachable; only a documented headline
"~factor of 10 better than existing detectors, band from a few Hz to >2 kHz" is
citable, source src_e8a31a47672c). So this file does NOT pretend to use ET-D.
It does something narrower and defensible: it takes the SAME aLIGO PSD SHAPE
used in analysis.py and asks how the unabsorbed SNR of the quadratic delay
responds to two independent knobs:

  * the overall noise scale, ASD -> ASD / S   (S = 1 aLIGO, S = 10 ET headline,
    S = 20 a deliberately optimistic CE-like value);
  * the low-frequency edge fmin (20 Hz aLIGO, 10 Hz and 5 Hz for a
    third-generation band).

Both knobs are exact in the statistic: the projection residual is linear in the
data and the weighting is 1/S_n, so scaling S_n by a constant scales every
unmodelled SNR by exactly S. Only the fmin dependence is a genuine re-computation.

The structural result -- that the const and linear delays are absorbed exactly,
at any S and any band -- is independent of all of this, because it is an
algebraic identity (see pn_orders.py).

CEILING ON THE CLAIM. The numbers below are a sensitivity study of the
STATISTIC, not a forecast for ET or CE. A real forecast needs the real curve and
a full IMR waveform, neither of which is available here.
"""
import json
import math
import os
import sys

import numpy as np

G = 6.67430e-11
C = 2.99792458e8
MSUN = 1.98892e30
TSUN = G * MSUN / C ** 3
MPC = 3.0856775814913673e22
TMPC = MPC / C

M1_S, M2_S = 36.0, 29.0
M1, M2 = M1_S * TSUN, M2_S * TSUN
MTOT = M1 + M2
ETA = M1 * M2 / MTOT ** 2
MCHIRP = MTOT * ETA ** 0.6
DL = 410.0 * TMPC
HERE = os.path.dirname(os.path.abspath(__file__))

FMAX, NF = 300.0, 4000
TAU0 = 1.2e-3          # the paper's claimed delay scale
SNR_TARGET = 20.0


def psi_chirp(f, Mc):
    return (3.0 / 128.0) * (np.pi * Mc * f) ** (-5.0 / 3.0)


def t_of_f(f, Mc):
    return -(5.0 / 256.0) * Mc ** (-5.0 / 3.0) * (np.pi * f) ** (-8.0 / 3.0)


def strain_amp(f, Mc, DL_):
    return (1.0 / DL_) * math.sqrt(5.0 / 24.0) * math.pi ** (-2.0 / 3.0) \
        * Mc ** (5.0 / 6.0) * f ** (-7.0 / 6.0)


def psd_shape(f):
    x = f / 215.0
    return 1e-49 * (x ** -4.14 - 5.0 * x ** -2.0
                    + 111.0 * (1.0 - x**2 + 0.5 * x**4) / (1.0 + 0.5 * x**2))


def project(u, cols, w):
    Q, R = np.linalg.qr(cols * w[:, None], mode="reduced")
    d = np.abs(np.diag(R))
    Qk = Q[:, d > d.max() * 1e-10]
    y = u * w
    n = float(np.linalg.norm(y))
    if n == 0.0:
        return 1.0, 0.0
    coef = Qk.T @ y
    return (float(np.linalg.norm(Qk @ coef) / n), float(np.linalg.norm(y - Qk @ coef)))


TAU2_FIXED = 1.2e-3 / 0.8444822952069388 ** 2   # tau2 fixed at the 20 Hz band value


def study(fmin, nscale, convention="fixed_peak_amplitude"):
    """Return the projection results for one (fmin, noise scale, convention) cell.

    CONVENTION MATTERS, AND THAT IS ITSELF A RESULT. The paper does not state the
    functional form of tau(t), so "the quadratic term" is not defined until a
    normalisation is chosen:

      fixed_peak_amplitude -- tau2 = 1.2 ms / t_span^2, i.e. the quadratic
          reaches 1.2 ms at the band edge.  A wider band (lower fmin) means a
          longer t_span, so the SAME 1.2 ms is spread over more time and tau2
          SHRINKS as t_span^-2.
      fixed_tau2 -- hold the coefficient tau2 fixed at its 20 Hz value and just
          widen the band.

    Under the first convention the residual DROPS as the band opens (the span
    effect beats the x^-13/2 growth); under the second it RISES. The sign of the
    ET/CE band benefit therefore depends on an unspecified modelling choice,
    which is exactly the gap the owner pointed at. Both are reported.
    """
    f = np.linspace(fmin, FMAX, NF)
    df = f[1] - f[0]
    h = strain_amp(f, MCHIRP, DL)
    k_cal = (float(np.linalg.norm(h * np.sqrt(4.0 * df / psd_shape(f)))) / SNR_TARGET) ** 2
    Sn = psd_shape(f) * k_cal / nscale ** 2          # better detector -> lower Sn
    w = np.sqrt(4.0 * df / Sn)

    dMc = -(5.0 / 3.0) * psi_chirp(f, MCHIRP) / MCHIRP
    cols = np.column_stack([h * dMc, h * 2 * np.pi * f, -h])
    t = t_of_f(f, MCHIRP)
    t_ref = float(t[0])
    t_span = abs(float(t[0]) - float(t[-1]))
    if convention == "fixed_peak_amplitude":
        tau2 = TAU0 / t_span ** 2
    elif convention == "fixed_tau2":
        tau2 = TAU2_FIXED
    else:
        raise ValueError(convention)

    out = {"fmin_Hz": fmin, "noise_scale": nscale, "convention": convention,
           "signal_snr": float(np.linalg.norm(h * w)), "t_span_s": t_span,
           "tau2": tau2}
    for name, tau in (("const", TAU0 * np.ones_like(t)),
                      ("linear", 1e-2 * (t - t_ref)),
                      ("quadratic", tau2 * (t - t_ref) ** 2)):
        u = h * (2 * np.pi * f * tau)
        af, snr = project(u, cols, w)
        out[name] = {"absorbed_fraction": af, "unmodelled_snr": snr,
                     "detectable_amplitude_ms_at_snr5":
                         (TAU0 * 5.0 / snr * 1e3) if snr > 0 else None}
    return out


def main():
    out = {"cells": {}, "notes": {}}
    for conv in ("fixed_peak_amplitude", "fixed_tau2"):
        out["cells"][conv] = {}
        for fmin in (5.0, 10.0, 20.0):
            for nscale in (1.0, 10.0, 20.0):
                key = f"fmin{fmin:g}_S{nscale:g}"
                out["cells"][conv][key] = study(fmin, nscale, conv)

    # the linearity check: scaling the noise by S must scale the unmodelled SNR
    # by exactly S at fixed fmin -- verified, not assumed
    lin = {}
    base = out["cells"]["fixed_peak_amplitude"]["fmin20_S1"]["quadratic"]["unmodelled_snr"]
    for nscale in (10.0, 20.0):
        got = out["cells"]["fixed_peak_amplitude"][f"fmin20_S{nscale:g}"]["quadratic"][
            "unmodelled_snr"]
        lin[f"S{nscale:g}_ratio"] = got / base
    out["noise_scaling_ratios"] = lin

    T = {}
    allc = [v for conv in out["cells"].values() for v in conv.values()]
    T["all_const_absorbed"] = all(v["const"]["absorbed_fraction"] > 1 - 1e-10 for v in allc)
    T["all_linear_absorbed"] = all(v["linear"]["absorbed_fraction"] > 1 - 1e-10 for v in allc)
    T["noise_scaling_is_exact"] = all(
        abs(r - s) < 1e-6 for s, r in ((10.0, lin["S10_ratio"]), (20.0, lin["S20_ratio"])))

    qf = {k: v["quadratic"]["unmodelled_snr"]
          for k, v in out["cells"]["fixed_peak_amplitude"].items()}
    qt = {k: v["quadratic"]["unmodelled_snr"]
          for k, v in out["cells"]["fixed_tau2"].items()}
    T["quadratic_grows_with_sensitivity"] = qf["fmin20_S20"] > qf["fmin20_S1"]
    T["band_effect_flips_sign_with_convention"] = (
        (qf["fmin5_S1"] < qf["fmin20_S1"]) and (qt["fmin5_S1"] > qt["fmin20_S1"]))
    out["thresholds"] = {k: bool(v) for k, v in T.items()}

    out["notes"]["interpretation"] = (
        "The const and linear delays are absorbed exactly (to 1e-10) in EVERY cell: "
        "any band, any sensitivity. That part is an algebraic identity, not a number. "
        "The quadratic residual scales EXACTLY linearly with detector sensitivity "
        f"(ratio {lin['S10_ratio']:.9f} for S=10, {lin['S20_ratio']:.9f} for S=20), so a "
        "10x-better detector multiplies it by 10. Its response to opening the band "
        "DOWNWARD flips sign with the normalisation convention: with the amplitude "
        f"pinned at 1.2 ms at the band edge it DROPS ({qf['fmin20_S1']:.4f} at 20 Hz vs "
        f"{qf['fmin5_S1']:.4f} at 5 Hz, because t_span^2 grows), while with tau2 held "
        f"fixed it RISES ({qt['fmin20_S1']:.4f} vs {qt['fmin5_S1']:.4f}). Since the paper "
        "never states tau(t)'s functional form, the sign of the ET/CE band benefit is "
        "not determined by the paper. That is the honest ceiling on this step."
    )

    os.makedirs(os.path.join(HERE, "artifacts"), exist_ok=True)
    with open(os.path.join(HERE, "artifacts", "et_ce_scaling.json"), "w") as fh:
        json.dump(out, fh, indent=2, sort_keys=True)
    print(json.dumps(out, indent=2, sort_keys=True))
    return 0 if all(T.values()) else 1


if __name__ == "__main__":
    sys.exit(main())
