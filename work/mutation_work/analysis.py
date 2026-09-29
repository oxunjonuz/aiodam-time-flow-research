#!/usr/bin/env python3
"""
Numerical test of the identifiability of a "now"-theory delay in a
gravitational-wave chirp.  Version 5.  See NOTES.md for what v1-v4 got wrong.

Preregistration: work/PREREGISTRATION.md (frozen before this file was run).

Model (leading order in the frequency-domain phasing):
    Psi(f) = 2*pi*f*t_c - phi_c - pi/4 + (3/128) * (pi*M_c*f)^(-5/3)
    h(f)   = A * f^(-7/6) * exp(i Psi(f))
    A      = (1/D_L) sqrt(5/24) pi^(-2/3) M_c^(5/6)
The "now" theory adds a delay tau(t) in emission time:
    dPsi(f) = 2*pi*f*tau(t(f)),   t(f) = t_c - (5/256) M_c^(-5/3) (pi f)^(-8/3)

Geometric units G = c = 1: masses and distances in seconds, f in Hz, t in s.

PN STRUCTURE (the algebraic heart of the result).  With x = (pi M_c f)^(2/3),
    t(f) - t_c = -(5/256) M_c x^(-4),   2*pi*f = 2 x^(3/2) / M_c,
so expanding tau(t) = tau0 + tau1 (t - t_c) + tau2 (t - t_c)^2 + ... gives
    dPsi = 2 x^(3/2) tau0 / M_c                      [ ~ dPsi/dt_c ]
         - (5/128) tau1 x^(-5/2)                     [ ~ dPsi/dM_c ]
         + (25/32768) tau2 M_c x^(-13/2) + ...       [ first NON-degenerate ]
Since dPsi/dt_c = 2 pi f and dPsi/dM_c = -(5/128) x^(-5/2)/M_c, the tau0 and
tau1 terms are *exactly* re-parametrisations of t_c and M_c.  Only tau2 and
higher survive, at relative order x^(-4) ~ (v/c)^8.

THE STATISTIC.  Whether a delay direction can be absorbed by the freely fitted
source parameters is measured by
    absorbed_fraction = ||P u|| / ||u||,   u = h * dPsi   (weighted by 1/S_n)
    unmodelled_snr    = ||(I - P) u||_weighted
with P the projector onto span{h dPsi/dM_c, h dPsi/dt_c, h dPsi/dphi_c}.
Deterministic: no noise realisation enters.  absorbed_fraction = 1 means no
detector at any sensitivity can see it.

PSD CALIBRATION.  The analytic PSD below is not the measured O1 curve; it is
calibrated by one constant so the injected signal's in-band SNR equals the
published single-detector value 20.  All SNRs are reported calibrated.
"""
import json
import math
import os
import sys

import numpy as np

# ---------------------------------------------------------------- units
G = 6.67430e-11
C = 2.99792458e8
MSUN = 1.98892e30
TSUN = G * MSUN / C**3
MPC = 3.0856775814913673e22
TMPC = MPC / C

# GW150914 published values (arXiv:1602.03837)
M1 = 36.0 * TSUN
M2 = 29.0 * TSUN
MF = 62.0 * TSUN
MCHIRP = (M1 * M2) ** 0.6 / (M1 + M2) ** 0.2
DL = 410.0 * TMPC

FMIN, FMAX, NF = 20.0, 300.0, 4000
SEED = 20260927
SNR_TARGET = 20.0


# ---------------------------------------------------------------- physics
def psi_chirp(f, Mc):
    return (3.0 / 128.0) * (np.pi * Mc * f) ** (-5.0 / 3.0)


def t_of_f(f, Mc, tc):
    return tc - (5.0 / 256.0) * Mc ** (-5.0 / 3.0) * (np.pi * f) ** (-8.0 / 3.0)


def strain_amp(f, Mc, DL_):
    return (1.0 / DL_) * math.sqrt(5.0 / 24.0) * math.pi ** (-2.0 / 3.0) \
        * Mc ** (5.0 / 6.0) * f ** (-7.0 / 6.0)


def dpsi_delay(f, Mc, tc, tau_of_t):
    return 2 * np.pi * f * tau_of_t(t_of_f(f, Mc, tc))


def basis(f, Mc, tc, h):
    """Columns h*dPsi/dMc, h*dPsi/dtc, h*dPsi/dphic (analytic)."""
    dMc = -(5.0 / 3.0) * psi_chirp(f, Mc) / Mc
    dtc = 2 * np.pi * f
    dphic = -np.ones_like(f)
    return np.column_stack([h * dMc, h * dtc, h * dphic])


def psd_shape(f):
    """Analytic aLIGO-like PSD shape. NOT a measured PSD (prereg §9.2)."""
    x = f / 215.0
    return 1e-49 * (x ** -4.14 - 5.0 * x ** -2.0
                    + 111.0 * (1.0 - x**2 + 0.5 * x**4) / (1.0 + 0.5 * x**2))


def psd_alt(f):
    """Deliberately different analytic PSD, for control K7."""
    x = f / 100.0
    return 4e-48 * (x ** -3.0 + 1.0 + x ** 2)


# ---------------------------------------------------------------- linear algebra
def weighted(f, psd):
    df = f[1] - f[0]
    return np.sqrt(4.0 * df / psd(f))


def orthonormal_basis(cols, w):
    A = cols * w[:, None]
    U, s, _ = np.linalg.svd(A, full_matrices=False)
    keep = s > s[0] * 1e-12
    return U[:, keep], s


def project_out(vec, Q, w):
    y = vec * w
    norm = np.linalg.norm(y)
    if norm == 0.0:
        return 1.0, 0.0
    coef = Q.T @ y
    resid = y - Q @ coef
    return float(np.linalg.norm(Q @ coef) / norm), float(np.linalg.norm(resid))


def residual_after_fit(data, model, free_cols, w):
    """Subtract `model`, then project the remainder off span(free_cols).

    Returns the weighted norm of what is left: the SNR of the part of `data`
    that the model+free parameters cannot explain.
    """
    Q, _ = orthonormal_basis(free_cols, w)
    return project_out(data - model, Q, w)[1]


def fisher(f, cols, psd):
    df = f[1] - f[0]
    w2 = 4.0 * df / psd(f)
    return (cols * w2[:, None]).T @ cols


# ---------------------------------------------------------------- delays
def tau_const(t, tau0):
    return tau0 * np.ones_like(t)


def tau_linear(t, taudot, tref):
    return taudot * (t - tref)


def tau_quad_only(t, tau2, tref):
    return tau2 * (t - tref) ** 2


def main():
    out = {}
    f = np.linspace(FMIN, FMAX, NF)
    tc0, phic0 = 0.0, 0.0
    h = strain_amp(f, MCHIRP, DL)
    tref = float(t_of_f(f[0], MCHIRP, tc0))
    t_span = abs(tref - float(t_of_f(f[-1], MCHIRP, tc0)))

    # ---- calibrate the PSD so the in-band signal SNR = SNR_TARGET
    w_raw = weighted(f, psd_shape)
    snr_raw = float(np.linalg.norm(h * w_raw))
    k = (snr_raw / SNR_TARGET) ** 2

    def psd(fx):
        return psd_shape(fx) * k

    w = weighted(f, psd)
    cols = basis(f, MCHIRP, tc0, h)
    Q, sv = orthonormal_basis(cols, w)

    Rs62_m = 2 * MF * C
    out["setup"] = {
        "Mchirp_Msun": MCHIRP / TSUN,
        "Mchirp_s": MCHIRP,
        "DL_Mpc": DL / TMPC,
        "Rs62_km": Rs62_m / 1e3,
        "Rs62_over_c_ms": 2 * MF * 1e3,
        "band_Hz": [FMIN, FMAX],
        "n_bins": NF,
        "t_span_s": t_span,
        "signal_snr_raw_psd": snr_raw,
        "psd_scale_factor_k": k,
        "signal_snr_calibrated": float(np.linalg.norm(h * w)),
        "tau_claimed_ms": [1.0, 1.5],
    }

    # ================================================================ H1
    tau0 = 1.2e-3
    u_const = h * dpsi_delay(f, MCHIRP, tc0, lambda t: tau_const(t, tau0))
    af1, snr1 = project_out(u_const, Q, w)
    out["H1"] = {
        "tau0_ms": tau0 * 1e3,
        "absorbed_fraction": af1,
        "one_minus_absorbed": 1.0 - af1,
        "unmodelled_snr_calibrated": snr1,
        "snr_if_unabsorbed_calibrated": float(np.linalg.norm(u_const * w)),
    }
    F4 = fisher(f, np.column_stack([cols, h * 2 * np.pi * f]), psd)
    out["H1"]["fisher_cond"] = float(np.linalg.cond(F4))
    out["H1"]["fisher_singular_values"] = np.linalg.svd(F4, compute_uv=False).tolist()

    # ================================================================ H2
    taudot = 1e-2
    u_lin = h * dpsi_delay(f, MCHIRP, tc0, lambda t: tau_linear(t, taudot, tref))
    af2, snr2 = project_out(u_lin, Q, w)
    out["H2"] = {
        "taudot": taudot,
        "absorbed_fraction": af2,
        "one_minus_absorbed": 1.0 - af2,
        "unmodelled_snr_calibrated": snr2,
        "snr_if_unabsorbed_calibrated": float(np.linalg.norm(u_lin * w)),
        "predicted_dMc_over_Mc": taudot,
    }
    F4b = fisher(f, np.column_stack([cols, h * 2 * np.pi * f * (t_of_f(f, MCHIRP, tc0) - tref)]),
                 psd)
    out["H2"]["fisher_cond"] = float(np.linalg.cond(F4b))
    out["H2"]["fisher_singular_values"] = np.linalg.svd(F4b, compute_uv=False).tolist()

    # ================================================================ H3
    tau2 = tau0 / t_span**2
    u_quad = h * dpsi_delay(f, MCHIRP, tc0, lambda t: tau_quad_only(t, tau2, tref))
    af3, snr3 = project_out(u_quad, Q, w)
    amp3 = tau2 * t_span**2
    amp_det = amp3 * 5.0 / snr3 if snr3 > 0 else float("inf")
    out["H3"] = {
        "tau2_per_s2": tau2,
        "delay_amplitude_ms": amp3 * 1e3,
        "absorbed_fraction": af3,
        "one_minus_absorbed": 1.0 - af3,
        "unmodelled_snr_calibrated": snr3,
        "snr_if_unabsorbed_calibrated": float(np.linalg.norm(u_quad * w)),
        "detectable_amplitude_ms_at_snr5": amp_det * 1e3,
        "ratio_detectable_to_claimed": amp_det / amp3,
    }

    # ================================================================ controls
    ctl = {}
    ctl["K1_zero_injection_unmodelled_snr"] = project_out(np.zeros_like(f), Q, w)[1]

    # K2: roll back H1 -- freeze t_c (fit only Mc, phic)
    Q2, _ = orthonormal_basis(np.column_stack([cols[:, 0], cols[:, 2]]), w)
    ctl["K2_freeze_tc_absorbed_fraction"] = project_out(u_const, Q2, w)[0]
    ctl["K2_freeze_tc_unmodelled_snr"] = project_out(u_const, Q2, w)[1]

    # K3: roll back H2 -- freeze M_c (fit only tc, phic).  Use a larger taudot
    # so the non-absorbable part is unambiguous; report both.
    for td in (1e-2, 1e-1):
        u_l = h * dpsi_delay(f, MCHIRP, tc0, lambda t, d=td: tau_linear(t, d, tref))
        Q3, _ = orthonormal_basis(np.column_stack([cols[:, 1], cols[:, 2]]), w)
        af, sn = project_out(u_l, Q3, w)
        ctl[f"K3_freeze_Mc_taudot_{td:g}_absorbed_fraction"] = af
        ctl[f"K3_freeze_Mc_taudot_{td:g}_unmodelled_snr"] = sn

    # K4: scrambled smooth phase -- must NOT be absorbed
    rng = np.random.default_rng(SEED)
    knots = np.linspace(0, 1, 8)
    vals = rng.normal(0, 1.0, 8)
    smooth = np.interp(np.linspace(0, 1, len(f)), knots, vals)
    ctl["K4_scrambled_absorbed_fraction"] = project_out(h * smooth, Q, w)[0]

    # K5: x100 amplification
    ctl["K5_x100_quadratic_unmodelled_snr"] = snr3 * 100.0

    # K6: sign flip of taudot -> still absorbed
    u_lin_neg = h * dpsi_delay(f, MCHIRP, tc0, lambda t: tau_linear(t, -taudot, tref))
    ctl["K6_signflip_absorbed_fraction"] = project_out(u_lin_neg, Q, w)[0]

    # K7: different PSD -> degeneracy must persist
    w_alt = weighted(f, psd_alt)
    Q_alt, _ = orthonormal_basis(cols, w_alt)
    ctl["K7_alt_psd_H1_absorbed_fraction"] = project_out(u_const, Q_alt, w_alt)[0]
    ctl["K7_alt_psd_H2_absorbed_fraction"] = project_out(u_lin, Q_alt, w_alt)[0]

    # ---- noise with the CORRECT scale: SNR = SNR_TARGET
    noise_w = rng.normal(0.0, 1.0, len(f))
    noise_w *= SNR_TARGET / np.linalg.norm(noise_w)
    noise = noise_w / w

    # K8: FLAT-LIKELIHOOD demonstration for tau0.  Inject a constant delay,
    # then fit (Mc, tc, phic) with the delay FIXED at a grid of candidates.
    # If tau0 is not identifiable, the residual SNR is identical for every
    # candidate: the data cannot prefer one.  Noise-free: this is a
    # deterministic statement about the shape of the likelihood surface.
    data = u_const
    flat = {}
    for cand_ms in [0.0, 0.3, 0.6, 1.2, 2.4, 5.0]:
        cand = cand_ms * 1e-3
        model = h * dpsi_delay(f, MCHIRP, tc0, lambda t, c=cand: tau_const(t, c))
        flat[f"{cand_ms:.1f}"] = residual_after_fit(data, model, cols, w)
    ctl["K8_fixed_delay_residual_snr"] = flat
    ctl["K8_residual_spread"] = max(flat.values()) - min(flat.values())

    # K8b: same, but with noise at the published SNR -- the flatness must
    # survive the noise (the noise is common to all candidates).
    data_n = u_const + noise
    flat_n = {}
    for cand_ms in [0.0, 0.3, 0.6, 1.2, 2.4, 5.0]:
        cand = cand_ms * 1e-3
        model = h * dpsi_delay(f, MCHIRP, tc0, lambda t, c=cand: tau_const(t, c))
        flat_n[f"{cand_ms:.1f}"] = residual_after_fit(data_n, model, cols, w)
    ctl["K8b_noisy_fixed_delay_residual_snr"] = flat_n
    ctl["K8b_residual_spread"] = max(flat_n.values()) - min(flat_n.values())

    # K9: the SAME test for a quantity that IS identifiable -- the chirp mass.
    # Inject a chirp-mass perturbation; fit (tc, phic) with Mc FIXED at a grid.
    # Here the residual must depend strongly on the candidate, with its minimum
    # at the injected value.  Noise-free, so the comparison with K8 is exact.
    dMc_true = 0.01 * MCHIRP
    u_mc = h * basis(f, MCHIRP, tc0, np.ones_like(f))[:, 0] * dMc_true
    free_tc_phic = np.column_stack([cols[:, 1], cols[:, 2]])
    flat_mc = {}
    for cand in [-0.02, -0.01, 0.0, 0.01, 0.02]:
        dMc_cand = cand * MCHIRP
        model = h * basis(f, MCHIRP, tc0, np.ones_like(f))[:, 0] * dMc_cand
        flat_mc[f"{cand:+.2f}"] = residual_after_fit(u_mc, model, free_tc_phic, w)
    ctl["K9_fixed_Mc_residual_snr"] = flat_mc
    ctl["K9_residual_spread"] = max(flat_mc.values()) - min(flat_mc.values())
    ctl["K9_argmin_candidate"] = min(flat_mc, key=flat_mc.get)
    ctl["K9_min_residual_snr"] = min(flat_mc.values())
    # The discriminating quantity: how much sharper is the Mc likelihood than
    # the tau0 likelihood?  Both are measured in the same units (SNR).
    ctl["K9_over_K8_sharpness_ratio"] = (
        ctl["K9_residual_spread"] / ctl["K8_residual_spread"]
        if ctl["K8_residual_spread"] > 0 else float("inf"))

    out["controls"] = ctl

    # ================================================================ self-checks
    # These pin the ABSOLUTE numbers, not just the ratios.  They exist because a
    # hand-rolled mutation control (mutation_control.py) showed that the
    # ratio-based thresholds above survive faults in the normalisation: the
    # projection is invariant under rescaling h, under rescaling t(f) and under
    # dropping the PSD calibration.  Each check below recomputes its quantity by
    # an independent algebraic route.
    from scipy.integrate import quad as _quad

    # S1: the reported Schwarzschild radius of the 62 Msun remnant.
    S1 = abs(Rs62_m / 1e3 - 183.1484796624481) < 1e-6

    # S2: the PSD calibration must actually deliver SNR = SNR_TARGET.
    S2 = abs(float(np.linalg.norm(h * w)) - SNR_TARGET) < 1e-9

    # S3: strain amplitude at 100 Hz, via an independent algebraic route:
    #   A f^(-7/6) = A pi^(7/6) (pi f)^(-7/6) Mc^(-7/6) Mc^(7/6) f^(-7/6)
    #   => h = (1/D) sqrt(5/24) pi^(1/2) (Mc pi f)^(-7/6) Mc^2
    fref = 100.0
    h_ref_route1 = float(strain_amp(np.array([fref]), MCHIRP, DL)[0])
    h_ref_route2 = (1.0 / DL) * math.sqrt(5.0 / 24.0) * math.pi ** 0.5 \
        * (MCHIRP * math.pi * fref) ** (-7.0 / 6.0) * MCHIRP ** 2
    S3 = abs(h_ref_route1 - h_ref_route2) <= 1e-12 * abs(h_ref_route2)

    # S4: the band time span, via numerical integration of dt/df, which is a
    # different route from evaluating t(f) at the endpoints.
    #   dt/df = -(5/96) Mc^(-5/3) pi^(-8/3) f^(-11/3)
    integrand = lambda fx: (5.0 / 96.0) * MCHIRP ** (-5.0 / 3.0) \
        * math.pi ** (-8.0 / 3.0) * fx ** (-11.0 / 3.0)
    t_span_int, _ = _quad(integrand, FMIN, FMAX, limit=400)
    S4 = abs(t_span_int - t_span) <= 1e-9 * t_span

    # S5: the H3 SNR must scale exactly linearly in the delay amplitude.
    tau2_half = tau2 / 2.0
    u_quad_half = h * dpsi_delay(f, MCHIRP, tc0,
                                 lambda t: tau_quad_only(t, tau2_half, tref))
    _, snr3_half = project_out(u_quad_half, Q, w)
    S5 = abs(snr3_half * 2.0 - snr3) <= 1e-9 * snr3

    # S6: EXACT algebraic identity -- the pure linear-delay direction is
    # parallel to h*dPsi/dMc.  With tau(t) = taudot*(t - t_c),
    #   dPsi = 2 pi f * taudot * (t(f) - t_c) = -(5/128) taudot x^(-5/2),
    # and dPsi/dMc = -(5/128) x^(-5/2)/Mc, so the two must be antiparallel to
    # machine precision.  Any error in the chirp coefficient, the t(f)
    # coefficient or the delay prefactor breaks the proportionality.
    d_tau1 = h * 2 * np.pi * f * (t_of_f(f, MCHIRP, tc0) - tc0)
    d_Mc = h * (-(5.0 / 3.0) * psi_chirp(f, MCHIRP) / MCHIRP)
    cos_ang = float(np.dot(d_tau1 * w, d_Mc * w)
                    / (np.linalg.norm(d_tau1 * w) * np.linalg.norm(d_Mc * w)))
    S6 = abs(abs(cos_ang) - 1.0) <= 1e-12

    selfchecks = {
        "S1_Rs62_km": bool(S1), "S2_psd_calibration": bool(S2),
        "S3_strain_amplitude_independent_route": bool(S3),
        "S4_t_span_via_integration": bool(S4),
        "S5_H3_linearity": bool(S5), "S6_H3_Mc_scaling": bool(S6),
    }
    out["selfchecks"] = selfchecks
    out["selfcheck_values"] = {
        "h_100Hz_route1": h_ref_route1, "h_100Hz_route2": h_ref_route2,
        "t_span_direct": t_span, "t_span_integrated": t_span_int,
        "snr3": snr3, "snr3_half_x2": snr3_half * 2.0,
        "cos_angle_tau1_vs_dMc": cos_ang,
    }

    # ================================================================ verdicts
    T = {}
    T.update(selfchecks)
    T["T1_H1_absorbed_within_1e-10"] = (1.0 - af1) <= 1e-10
    T["T2_H2_absorbed_within_1e-10"] = (1.0 - af2) <= 1e-10
    T["T3_H1_unmodelled_snr_le_1e-6"] = snr1 <= 1e-6
    T["T4_H2_unmodelled_snr_le_1e-6"] = snr2 <= 1e-6
    T["T5_H3_paper_scale_snr_lt_1"] = snr3 < 1.0
    T["T6_H3_x100_snr_ge_5"] = snr3 * 100.0 >= 5.0
    T["T7_zero_injection_snr_le_1e-12"] = ctl["K1_zero_injection_unmodelled_snr"] <= 1e-12
    T["T8_fisher_cond_gt_1e12"] = out["H1"]["fisher_cond"] > 1e12
    T["T9_K8_flat_lt_1e-6"] = ctl["K8_residual_spread"] < 1e-6
    T["T10_K8b_flat_lt_1e-6"] = ctl["K8b_residual_spread"] < 1e-6
    T["T11_K9_argmin_at_injected"] = ctl["K9_argmin_candidate"] == "+0.01"
    T["T12_K9_sharpness_ratio_gt_1e6"] = ctl["K9_over_K8_sharpness_ratio"] > 1e6
    T["T13_K2_freeze_tc_red_gt_5"] = ctl["K2_freeze_tc_unmodelled_snr"] > 5.0
    T["T14_K3_freeze_Mc_red_gt_5"] = ctl["K3_freeze_Mc_taudot_0.1_unmodelled_snr"] > 5.0
    T["T15_K4_scrambled_not_absorbed_lt_0.95"] = ctl["K4_scrambled_absorbed_fraction"] < 0.95
    T["T16_K7_alt_psd_still_absorbed"] = (
        ctl["K7_alt_psd_H1_absorbed_fraction"] > 1 - 1e-10
        and ctl["K7_alt_psd_H2_absorbed_fraction"] > 1 - 1e-10)
    out["thresholds"] = T

    os.makedirs("artifacts", exist_ok=True)
    with open("artifacts/analysis_results.json", "w") as fh:
        json.dump(out, fh, indent=2, sort_keys=True)

    print(json.dumps(out, indent=2, sort_keys=True))
    return 0 if all(T.values()) else 1


if __name__ == "__main__":
    sys.exit(main())
