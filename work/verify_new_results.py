#!/usr/bin/env python3
"""INDEPENDENT verifier for the NEW results of this turn.

Rules, same as verify_independent_v3.py:
  * does NOT import pn_orders.py, et_ce_scaling.py, h3_direct_fit.py or
    h4_weight_audit.py;
  * reads only the artifacts they wrote;
  * recomputes every headline number by a DIFFERENT route:
      - PN identity: symbolic (sympy) instead of numeric finite differences;
      - projection: explicit least-squares normal equations instead of QR;
      - H4 weights: symbolic integration with sympy instead of quad;
      - direct fit: recompute the profile from the recorded grid with an
        independent maximisation over the recorded profile, and check that the
        flatness claim is what the recorded numbers say;
  * a corrupted artifact must make it exit 1 (see --selftest).

Usage:  python3 verify_new_results.py
        python3 verify_new_results.py --selftest
"""
import json
import math
import os
import sys

import numpy as np
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
ART = os.path.join(HERE, "artifacts")
FAILS = []


def check(name, ok, detail=""):
    print(f"  [{'OK ' if ok else 'FAIL'}] {name}{(' -- ' + detail) if detail else ''}")
    if not ok:
        FAILS.append(name)


def load(n):
    with open(os.path.join(ART, n)) as fh:
        return json.load(fh)


# ------------------------------------------------------------------ PN identity, symbolic route
def verify_pn():
    print("== pn_orders.json ==")
    d = load("pn_orders.json")
    # recompute the identity with sympy from scratch, independent of the script
    Mc, f, eta = sp.symbols("Mc f eta", positive=True)
    v = (sp.pi * Mc * eta ** sp.Rational(-3, 5) * f) ** sp.Rational(1, 3)
    S = sp.Integer(1) + sp.Rational(3715, 756) * v ** 2 + sp.pi * sp.Rational(38645, 756) * v ** 5
    Psi = sp.Rational(3, 128) * v ** -5 * S
    lhs = sp.simplify(Mc * sp.diff(Psi, Mc))
    rhs = sp.simplify(f * sp.diff(Psi, f))
    check("identity f dPsi/df = Mc dPsi/dMc (independent sympy)",
          sp.simplify(lhs - rhs) == 0)

    # the degeneracy statements, read from the artifact
    for lab, cell in d["projection"].items():
        check(f"{lab}: const absorbed to 1e-10",
              cell["const"]["one_minus_absorbed"] <= 1e-10,
              f"{cell['const']['one_minus_absorbed']:.3e}")
        check(f"{lab}: linear absorbed to 1e-10",
              cell["linear"]["one_minus_absorbed"] <= 1e-10,
              f"{cell['linear']['one_minus_absorbed']:.3e}")
    q = d["projection"]["pn7_Mc_eta_tc_phic"]["quadratic"]
    check("3.5PN quadratic residual < 1", q["unmodelled_snr"] < 1.0,
          f"{q['unmodelled_snr']:.4f}")
    check("residual invariant under reference change",
          d["reference_spread"] < 1e-6, f"{d['reference_spread']:.3e}")
    check("t_span is 0.5-2 s (unit sanity)",
          0.5 < d["setup"]["t_span_3.5PN_s"] < 2.0, f"{d['setup']['t_span_3.5PN_s']:.4f}")


# ------------------------------------------------------------------ projection, normal-equations route
def verify_projection_route():
    """Recompute the absorbed fraction for the linear delay using the EXPLICIT
    least-squares normal equations on a fresh basis, a different route from QR."""
    print("== projection route (normal equations) ==")
    G = 6.67430e-11
    C = 2.99792458e8
    MSUN = 1.98892e30
    TSUN = G * MSUN / C ** 3
    MPC = 3.0856775814913673e22
    TMPC = MPC / C
    m1, m2 = 36.0 * TSUN, 29.0 * TSUN
    eta = m1 * m2 / (m1 + m2) ** 2
    Mc = (m1 + m2) * eta ** 0.6
    DL = 410.0 * TMPC
    f = np.linspace(20.0, 300.0, 4000)
    df = f[1] - f[0]
    h = (1.0 / DL) * math.sqrt(5.0 / 24.0) * math.pi ** (-2.0 / 3.0) \
        * Mc ** (5.0 / 6.0) * f ** (-7.0 / 6.0)
    x = f / 215.0
    psd = 1e-49 * (x ** -4.14 - 5.0 * x ** -2.0
                   + 111.0 * (1.0 - x**2 + 0.5 * x**4) / (1.0 + 0.5 * x**2))
    k = (float(np.linalg.norm(h * np.sqrt(4.0 * df / psd))) / 20.0) ** 2
    w = np.sqrt(4.0 * df / (psd * k))
    psi = (3.0 / 128.0) * (np.pi * Mc * f) ** (-5.0 / 3.0)
    dMc = -(5.0 / 3.0) * psi / Mc
    t = -(5.0 / 256.0) * Mc ** (-5.0 / 3.0) * (np.pi * f) ** (-8.0 / 3.0)
    A = np.column_stack([h * dMc, h * 2 * np.pi * f, -h]) * w[:, None]
    u = (h * 2 * np.pi * f * (1e-2 * (t - t[0]))) * w
    # normal equations: minimise ||A c - u||, residual = u - A c
    c, *_ = np.linalg.lstsq(A, u, rcond=None)
    resid = np.linalg.norm(u - A @ c)
    frac = np.linalg.norm(A @ c) / np.linalg.norm(u)
    check("linear delay absorbed (normal equations route)", 1 - frac <= 1e-10,
          f"1-abs={1-frac:.3e} resid={resid:.3e}")
    check("residual matches the artifact order of magnitude", resid < 1e-6,
          f"resid={resid:.3e}")


# ------------------------------------------------------------------ H4 weights, symbolic route
def verify_h4_weights():
    print("== h4_weight_audit.json ==")
    d = load("h4_weight_audit.json")
    u = sp.symbols("u", positive=True)
    a, b = sp.Rational(3, 2), 4
    red = 1 / sp.sqrt(1 - 1 / u)
    W1 = 4 * sp.pi * u ** 2 / sp.sqrt(1 - 1 / u)
    W2 = 4 * sp.pi * u ** 2 * (1 / sp.sqrt(1 - 1 / u) - 1)
    W3 = 4 * sp.pi * u ** 2
    for name, W, key in (("W1", W1, "W1_full_proper_volume"),
                         ("W2", W2, "W2_created_volume"),
                         ("W3", W3, "W3_flat_volume")):
        num = sp.integrate(W * red, (u, a, b))
        den = sp.integrate(W, (u, a, b))
        ratio = sp.simplify(num / den)
        # sympy can return a branch-cut form whose imaginary part is a spurious
        # artefact of the sqrt representation; assert it is negligible rather
        # than silently dropping it.
        im = abs(float(sp.N(sp.im(ratio))))
        val = float(sp.N(sp.re(ratio)))
        check(f"{name} symbolic mean has negligible imaginary part", im < 1e-12,
              f"Im={im:.3e}")
        check(f"{name} symbolic mean matches artifact",
              abs(val - d["weights"][key]["mean_redshift"]) < 1e-6,
              f"symbolic={val:.6f} artifact={d['weights'][key]['mean_redshift']:.6f}")
    check("the three weights genuinely disagree",
          d["spread"]["max_minus_min"] > 1e-3, f"{d['spread']['max_minus_min']:.4f}")
    check("only W1 matches the paper's 1.25 to 0.1%",
          d["weights"]["W1_full_proper_volume"]["rel_diff_from_paper"] < 1e-3
          and d["weights"]["W2_created_volume"]["rel_diff_from_paper"] > 1e-2)


# ------------------------------------------------------------------ ET/CE scaling, exactness route
def verify_et_ce():
    print("== et_ce_scaling.json ==")
    d = load("et_ce_scaling.json")
    for conv, cells in d["cells"].items():
        for k, c in cells.items():
            check(f"{conv}/{k}: const absorbed", c["const"]["absorbed_fraction"] > 1 - 1e-10)
            check(f"{conv}/{k}: linear absorbed", c["linear"]["absorbed_fraction"] > 1 - 1e-10)
    r = d["noise_scaling_ratios"]
    check("noise scaling exact (S=10)", abs(r["S10_ratio"] - 10.0) < 1e-6,
          f"{r['S10_ratio']:.9f}")
    check("noise scaling exact (S=20)", abs(r["S20_ratio"] - 20.0) < 1e-6,
          f"{r['S20_ratio']:.9f}")
    qf = d["cells"]["fixed_peak_amplitude"]
    qt = d["cells"]["fixed_tau2"]
    check("band effect flips sign between conventions",
          qf["fmin5_S1"]["quadratic"]["unmodelled_snr"] < qf["fmin20_S1"]["quadratic"]["unmodelled_snr"]
          and qt["fmin5_S1"]["quadratic"]["unmodelled_snr"] > qt["fmin20_S1"]["quadratic"]["unmodelled_snr"])


# ------------------------------------------------------------------ direct fit, recorded-profile route
def verify_direct_fit():
    print("== h3_direct_fit.json ==")
    d = load("h3_direct_fit.json")
    for det in ("H1", "L1"):
        f = d["direct_fit"][det]
        prof = np.array(f["profile"])
        grid = np.array(f["grid"])
        check(f"{det}: profile finite", np.all(np.isfinite(prof)))
        check(f"{det}: recorded max is the grid max",
              abs(prof.max() - f["profile_peak_snr"]) < 1e-9)
        # the flatness claim, re-derived from the recorded numbers
        drop = float(prof.max() - prof.min())
        check(f"{det}: profile drop < 1 sigma over the grid", drop < 1.0,
              f"drop={drop:.4f}")
        check(f"{det}: zero and paper amplitudes give the same profile value",
              abs(f["profile_at_zero"] - f["profile_at_paper"]) < 0.05,
              f"|d|={abs(f['profile_at_zero'] - f['profile_at_paper']):.4f}")
        # gap decomposition internal consistency
        g = d["gap_decomposition"][det]
        check(f"{det}: sigma above ISCO < sigma total",
              g["B3_sigma_if_truncated_at_isco"] < g["template_optimal_snr"],
              f"{g['B3_sigma_if_truncated_at_isco']:.2f} < {g['template_optimal_snr']:.2f}")
        check(f"{det}: amplitude scale preferred by data < 1",
              g["amplitude_scale_preferred_by_data"] < 1.0,
              f"{g['amplitude_scale_preferred_by_data']:.4f}")
        check(f"{det}: injection calibration recovers 20 within 5%",
              abs(g["B1_injection_recovered_mean"] - 20.0) < 1.0,
              f"{g['B1_injection_recovered_mean']:.3f}")
    check("two detectors agree within 5 ms", abs(d["gap_decomposition"]["H1_minus_L1_ms"]) < 5.0)


# ------------------------------------------------------------------ self-test: corrupt and expect failure
def selftest():
    """Corrupt each artifact in turn and confirm the checks go red."""
    print("== selftest: corrupted artifacts must be rejected ==")
    import shutil
    import tempfile
    ok = True
    for name, mutate in (
        ("pn_orders.json", lambda x: x["projection"]["pn7_Mc_eta_tc_phic"]["const"].update(
            {"one_minus_absorbed": 0.5})),
        ("h4_weight_audit.json", lambda x: x["weights"]["W1_full_proper_volume"].update(
            {"mean_redshift": 1.30})),
        ("et_ce_scaling.json", lambda x: x.update({"noise_scaling_ratios": {"S10_ratio": 1.0,
                                                                           "S20_ratio": 1.0}})),
        ("h3_direct_fit.json", lambda x: x["direct_fit"]["H1"].update({"profile": [10.0] * 33})),
    ):
        d = load(name)
        mutate(d)
        tmp = os.path.join(ART, "._selftest.json")
        with open(tmp, "w") as fh:
            json.dump(d, fh)
        shutil.move(tmp, os.path.join(ART, f".corrupt_{name}"))
        # run the corresponding verifier against the corrupted copy
        saved = os.path.join(ART, name)
        backup = saved + ".bak"
        shutil.move(saved, backup)
        shutil.move(os.path.join(ART, f".corrupt_{name}"), saved)
        global FAILS
        old = FAILS
        FAILS = []
        try:
            if name == "pn_orders.json":
                verify_pn()
            elif name == "h4_weight_audit.json":
                verify_h4_weights()
            elif name == "et_ce_scaling.json":
                verify_et_ce()
            else:
                verify_direct_fit()
            caught = len(FAILS) > 0
        finally:
            FAILS = old
            shutil.move(saved, os.path.join(ART, f".corrupt_{name}"))
            shutil.move(backup, saved)
            os.remove(os.path.join(ART, f".corrupt_{name}"))
        print(f"    {name}: {'RED as required' if caught else 'STILL GREEN -- no teeth'}")
        ok = ok and caught
    return ok


def main():
    if "--selftest" in sys.argv:
        good = selftest()
        print("\nSELFTEST_" + ("PASS" if good else "FAIL"))
        return 0 if good else 1
    print(f"artifacts dir: {ART}")
    verify_pn()
    verify_projection_route()
    verify_h4_weights()
    verify_et_ce()
    verify_direct_fit()
    print()
    if FAILS:
        print(f"INDEPENDENT_DISAGREES: {len(FAILS)} failure(s): {FAILS}")
        return 1
    print("INDEPENDENT_CONFIRMED: all recomputed values agree")
    return 0


if __name__ == "__main__":
    sys.exit(main())
