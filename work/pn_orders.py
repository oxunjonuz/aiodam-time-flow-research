#!/usr/bin/env python3
"""Higher-PN orders, and whether the Lean criterion is genuinely instantiated.

Owner's questions answered here:
 (1) the Lean criterion is stated for an abstract linear model A theta + B tau --
     is it instantiated by the concrete physics vectors?
 (2) H1/H2 were proved at leading order -- does the degeneracy survive to 3.5PN?

TWO EXACT IDENTITIES, valid at EVERY PN order for any smooth chirp phase:

 (I)  dPsi/dt_c = 2 pi f                       (t_c enters linearly)
 (II) f dPsi_chirp/df = M_c dPsi_chirp/dM_c

 (II) is Euler's theorem in one variable, because Psi_chirp depends on M_c and f
 ONLY through v = (pi M_c eta^(-3/5) f)^(1/3). No PN coefficient enters, so no
 truncation order can break it. It is checked SYMBOLICALLY (sympy) on the
 truncated series with a NEGATIVE CONTROL that adds a spurious explicit M_c
 dependence and must make it fail.

Consequence. The stationary-phase time is t(f) - t_c = (1/2pi) dPsi_chirp/df, so
a linear delay tau(t) = taudot (t - t_ref) gives
    dPsi_lin = 2 pi f taudot (t(f) - t_ref)
             = taudot [ M_c dPsi/dM_c + (t_c - t_ref) dPsi/dt_c ],
EXACTLY a combination of the nuisance directions at every PN order.

METHOD NOTE. Each PN order is tested in its OWN basis and with its OWN t(f):
the previous version mixed a 3.5PN delay with a leading-order basis and reported
a spurious 4.8e-9 leak. Fixed here. The delay directions are referenced to
t_ref = t(f_min) as in analysis.py, so the quadratic term vanishes at the band
edge; the unabsorbed SNR is invariant under that choice (a shift of reference
only adds absorbed constant+linear pieces), and that invariance is checked.

Also recorded in NOTES.md: a first version put the 3PN log term -(6848/21)ln(4v)
OUTSIDE the v^6 factor, which made the 3.5PN phase at 20 Hz wrong by ~27x and the
reported band time span 32.5 s instead of 0.84 s.
"""
import json
import math
import os
import sys

import numpy as np
import sympy as sp

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

FMIN, FMAX, NF = 20.0, 300.0, 4000
HERE = os.path.dirname(os.path.abspath(__file__))
GAMMA_E = 0.5772156649015329


def alpha(k, eta):
    if k == 0:
        return 1.0
    if k == 1:
        return 0.0
    if k == 2:
        return 3715.0 / 756.0 + 55.0 * eta / 9.0
    if k == 3:
        return -16.0 * math.pi
    if k == 4:
        return 4.0 * (15293365.0 / 508032.0 + 27145.0 * eta / 504.0
                      + 3085.0 * eta ** 2 / 72.0)
    if k == 5:
        return math.pi * (38645.0 / 756.0 - 65.0 * eta / 9.0)
    if k == 6:
        return (11583231236531.0 / 4694215680.0 - 640.0 * math.pi ** 2 / 3.0
                - 6848.0 * GAMMA_E / 21.0) \
            + eta * (-15737765635.0 / 3048192.0 + 2255.0 * math.pi ** 2 / 12.0) \
            + eta ** 2 * (76055.0 / 1728.0) - eta ** 3 * (127825.0 / 1296.0)
    if k == 7:
        return math.pi * (77096675.0 / 254016.0 + 378515.0 * eta / 1512.0
                          - 74045.0 * eta ** 2 / 756.0)
    return 0.0


def S_of_v(v, eta, pnmax):
    s = np.zeros_like(v)
    for k in range(0, pnmax + 1):
        a = alpha(k, eta)
        if k == 5:
            s = s + a * v ** k * (1.0 + 3.0 * np.log(v))
        elif k == 6:
            s = s + (a - (6848.0 / 21.0) * np.log(4.0 * v)) * v ** k
        else:
            s = s + a * v ** k
    return s


def v_of_f(f, Mc, eta):
    return (np.pi * Mc * eta ** (-0.6) * f) ** (1.0 / 3.0)


def psi_chirp(f, Mc, eta, pnmax=7):
    v = v_of_f(f, Mc, eta)
    return (3.0 / (128.0 * eta * v ** 5)) * S_of_v(v, eta, pnmax)


def dpsi_df(f, Mc, eta, pnmax=7):
    """Analytic dPsi_chirp/df, chain rule through v, difference only in v."""
    v = v_of_f(f, Mc, eta)
    h = 1e-7

    def S(vv):
        return S_of_v(vv, eta, pnmax)

    Sp = (S(v * (1 + h)) - S(v * (1 - h))) / (2 * h * v)
    dpsi_dv = (3.0 / (128.0 * eta)) * (-5.0 * v ** -6 * S(v) + v ** -5 * Sp)
    return dpsi_dv * (v / (3.0 * f))


def t_of_f(f, Mc, eta, pnmax=7):
    return dpsi_df(f, Mc, eta, pnmax) / (2.0 * np.pi)


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


# ------------------------------------------------------------------ symbolic
def symbolic_identity():
    """Prove f dPsi/df = M_c dPsi/dM_c as an identity in v, and show it FAILS
    when a spurious explicit M_c dependence is added (negative control)."""
    out = {}
    Mc_s, f_s, eta_s = sp.symbols("Mc f eta", positive=True)
    v_s = (sp.pi * Mc_s * eta_s ** sp.Rational(-3, 5) * f_s) ** sp.Rational(1, 3)

    def make_S(vv, pnmax, spurious=False):
        S = sp.Integer(0)
        for k in range(0, pnmax + 1):
            a = alpha(k, 1)  # numeric eta; the check does not depend on its value
            if k == 5:
                S += a * vv ** k * (1 + 3 * sp.log(vv))
            elif k == 6:
                S += (a - sp.Rational(6848, 21) * sp.log(4 * vv)) * vv ** k
            else:
                S += a * vv ** k
        if spurious:
            # A term that depends on M_c and f NOT only through v. Note
            # log(M_c * f) would NOT work: M_c*f is a function of v alone, so it
            # would satisfy the identity too and the control would be vacuous.
            # log(M_c / f) is genuinely not a function of v, so it must break it.
            S += sp.log(Mc_s / f_s)
        return S

    for pnmax in (0, 2, 7):
        Psi = sp.Rational(3, 128) * v_s ** -5 * make_S(v_s, pnmax)
        A = sp.simplify(Mc_s * sp.diff(Psi, Mc_s))
        B = sp.simplify(f_s * sp.diff(Psi, f_s))
        out[f"pn{pnmax}_identity_holds"] = bool(sp.simplify(A - B) == 0)

        Psi_bad = sp.Rational(3, 128) * v_s ** -5 * make_S(v_s, pnmax, spurious=True)
        A_bad = sp.simplify(Mc_s * sp.diff(Psi_bad, Mc_s))
        B_bad = sp.simplify(f_s * sp.diff(Psi_bad, f_s))
        out[f"pn{pnmax}_negative_control_breaks"] = bool(sp.simplify(A_bad - B_bad) != 0)
    return out


def main():
    out = {"symbolic": {}, "identities": {}, "projection": {}, "reference_invariance": {}}
    f = np.linspace(FMIN, FMAX, NF)
    h_strain = strain_amp(f, MCHIRP, DL)
    df = f[1] - f[0]
    k_cal = (float(np.linalg.norm(h_strain * np.sqrt(4.0 * df / psd_shape(f)))) / 20.0) ** 2
    w = np.sqrt(4.0 * df / (psd_shape(f) * k_cal))

    out["symbolic"] = symbolic_identity()

    # ------------------------------------------------ identity (I): LO t(f) closed form
    t_lo = t_of_f(f, MCHIRP, ETA, 0)
    t_closed = -(5.0 / 256.0) * MCHIRP ** (-5.0 / 3.0) * (np.pi * f) ** (-8.0 / 3.0)
    idents = {"I_LO_t_of_f_matches_closed_form": bool(
        np.max(np.abs(t_lo - t_closed)) <= 1e-6 * np.max(np.abs(t_closed)))}

    # ------------------------------------------------ identity (II) numeric, per order
    for pnmax in (0, 7):
        ddf = dpsi_df(f, MCHIRP, ETA, pnmax)
        eps = 1e-7
        dMc = (psi_chirp(f, MCHIRP * (1 + eps), ETA, pnmax)
               - psi_chirp(f, MCHIRP * (1 - eps), ETA, pnmax)) / (2 * eps * MCHIRP)
        rel = float(np.max(np.abs(f * ddf - MCHIRP * dMc)) / np.max(np.abs(MCHIRP * dMc)))
        out[f"II_numeric_rel_pn{pnmax}"] = rel
        idents[f"II_numeric_pn{pnmax}"] = bool(rel < 1e-5)
    out["identities"] = idents

    # ------------------------------------------------ per-order projection
    tau0 = 1.2e-3
    taudot = 1e-2
    for pnmax in (0, 7):
        t_pn = t_of_f(f, MCHIRP, ETA, pnmax)
        t_span = abs(float(t_pn[0]) - float(t_pn[-1]))
        t_ref = float(t_pn[0])
        eps = 1e-6
        dMc = (psi_chirp(f, MCHIRP * (1 + eps), ETA, pnmax)
               - psi_chirp(f, MCHIRP * (1 - eps), ETA, pnmax)) / (2 * eps * MCHIRP)
        cols = [h_strain * dMc, h_strain * 2 * np.pi * f, -h_strain]
        if pnmax == 7:
            d_eta = (psi_chirp(f, MCHIRP, ETA * (1 + eps), pnmax)
                     - psi_chirp(f, MCHIRP, ETA * (1 - eps), pnmax)) / (2 * eps * ETA)
            cols = [h_strain * dMc, h_strain * d_eta, h_strain * 2 * np.pi * f, -h_strain]
        A = np.column_stack(cols)
        tau2 = tau0 / t_span ** 2
        dirs = {
            "const": h_strain * (2 * np.pi * f * tau0),
            "linear": h_strain * (2 * np.pi * f * taudot * (t_pn - t_ref)),
            "quadratic": h_strain * (2 * np.pi * f * tau2 * (t_pn - t_ref) ** 2),
        }
        label = f"pn{pnmax}" + ("_Mc_eta_tc_phic" if pnmax == 7 else "_Mc_tc_phic")
        out["projection"][label] = {"t_span_s": t_span}
        for dname, u in dirs.items():
            af, snr = project(u, A, w)
            out["projection"][label][dname] = {
                "absorbed_fraction": af, "one_minus_absorbed": 1.0 - af,
                "unmodelled_snr": snr, "total_snr": float(np.linalg.norm(u * w))}

    # ------------------------------------------------ reference invariance of the residual
    # shifting the quadratic reference only adds absorbed constant+linear pieces,
    # so the UNMODELLED SNR must be invariant.
    t_pn = t_of_f(f, MCHIRP, ETA)
    t_span = abs(float(t_pn[0]) - float(t_pn[-1]))
    tau2 = tau0 / t_span ** 2
    eps = 1e-6
    dMc = (psi_chirp(f, MCHIRP * (1 + eps), ETA) - psi_chirp(f, MCHIRP * (1 - eps), ETA)) \
        / (2 * eps * MCHIRP)
    d_eta = (psi_chirp(f, MCHIRP, ETA * (1 + eps)) - psi_chirp(f, MCHIRP, ETA * (1 - eps))) \
        / (2 * eps * ETA)
    A = np.column_stack([h_strain * dMc, h_strain * d_eta, h_strain * 2 * np.pi * f, -h_strain])
    for ref_name, ref in (("band_start", float(t_pn[0])), ("coalescence", 0.0),
                          ("band_end", float(t_pn[-1]))):
        u = h_strain * (2 * np.pi * f * tau2 * (t_pn - ref) ** 2)
        af, snr = project(u, A, w)
        out["reference_invariance"][ref_name] = {
            "absorbed_fraction": af, "unmodelled_snr": snr}

    snrs = [v["unmodelled_snr"] for v in out["reference_invariance"].values()]
    ref_spread = (max(snrs) - min(snrs)) / max(snrs)
    # The ABSORBED FRACTION, unlike the residual, is NOT invariant under the
    # choice of reference.  This is a finding about the metric, not about the
    # physics: a nearly-degenerate direction can be 99.9997% absorbed at one
    # reference and 95.0% at another while its residual is bit-identical. So the
    # absorbed fraction must never be quoted as the evidence for H3 -- the
    # residual SNR must.
    afs = [v["absorbed_fraction"] for v in out["reference_invariance"].values()]
    af_spread = max(afs) - min(afs)

    T = {}
    T.update(out["symbolic"])
    T.update(idents)
    T["A_const_absorbed_3.5PN"] = \
        out["projection"]["pn7_Mc_eta_tc_phic"]["const"]["one_minus_absorbed"] <= 1e-10
    T["B_linear_absorbed_3.5PN"] = \
        out["projection"]["pn7_Mc_eta_tc_phic"]["linear"]["one_minus_absorbed"] <= 1e-10
    T["C_const_absorbed_LO"] = \
        out["projection"]["pn0_Mc_tc_phic"]["const"]["one_minus_absorbed"] <= 1e-10
    T["D_linear_absorbed_LO"] = \
        out["projection"]["pn0_Mc_tc_phic"]["linear"]["one_minus_absorbed"] <= 1e-10
    T["E_quadratic_residual_lt_1_3.5PN"] = \
        out["projection"]["pn7_Mc_eta_tc_phic"]["quadratic"]["unmodelled_snr"] < 1.0
    T["F_quadratic_residual_lt_1_LO"] = \
        out["projection"]["pn0_Mc_tc_phic"]["quadratic"]["unmodelled_snr"] < 1.0
    T["G_residual_invariant_under_reference"] = bool(ref_spread < 1e-6)
    T["H_absorbed_fraction_NOT_reference_invariant"] = bool(af_spread > 1e-3)
    T["I_t_span_plausible"] = bool(0.5 < t_span < 2.0)

    out["thresholds"] = {kk: bool(vv) for kk, vv in T.items()}
    out["reference_spread"] = ref_spread
    out["setup"] = {"Mchirp_Msun": MCHIRP / TSUN, "eta": ETA, "t_span_3.5PN_s": t_span,
                    "tau0_ms": tau0 * 1e3, "taudot": taudot}

    os.makedirs(os.path.join(HERE, "artifacts"), exist_ok=True)
    with open(os.path.join(HERE, "artifacts", "pn_orders.json"), "w") as fh:
        json.dump(out, fh, indent=2, sort_keys=True)
    print(json.dumps(out, indent=2, sort_keys=True))
    return 0 if all(T.values()) else 1


if __name__ == "__main__":
    sys.exit(main())
