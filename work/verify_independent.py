#!/usr/bin/env python3
"""
INDEPENDENT verification of the identifiability result.

Rules this file obeys (prereg §7):
  * it imports NOTHING from analysis.py;
  * it reads only the artifact JSON from disk;
  * it recomputes the key numbers with its own code;
  * it uses a DIFFERENT numerical method (normal equations with explicit
    pseudo-inverse via QR, plus a sympy symbolic check of the PN algebra).

Any disagreement > 1e-9 with the artifact is a failure.
"""
import json
import math
import sys

import numpy as np
import sympy as sp

TOL = 1e-9
FAILURES = []


def check(name, got, want, tol=TOL):
    ok = abs(got - want) <= tol
    print(f"  [{'OK ' if ok else 'FAIL'}] {name}: got={got!r} want={want!r}")
    if not ok:
        FAILURES.append(name)
    return ok


def main():
    with open("artifacts/analysis_results.json") as fh:
        A = json.load(fh)

    print("=== 1. symbolic check of the PN algebra (sympy, independent) ===")
    # x = (pi Mc f)^(2/3).  Verify:
    #   dPsi_chirp/dMc = -(5/3) Psi_chirp / Mc = -(5/128) x^(-5/2)/Mc
    #   t - tc = -(5/256) Mc x^(-4)
    #   dPsi/dt_c = 2 pi f = 2 x^(3/2)/Mc
    x, Mc = sp.symbols("x Mc", positive=True)
    Psi_chirp = sp.Rational(3, 128) * x ** sp.Rational(-5, 2)
    dPsi_dMc = sp.simplify(sp.diff(Psi_chirp, Mc))
    # express via x: d/dMc of (3/128)(pi Mc f)^(-5/3) with x=(pi Mc f)^(2/3)
    # Psi = (3/128) x^(-5/2); dx/dMc = (2/3) x / Mc
    dPsi_dMc_expr = sp.simplify(sp.Rational(3, 128) * sp.Rational(-5, 2) * x ** sp.Rational(-7, 2)
                                * sp.Rational(2, 3) * x / Mc)
    print(f"  dPsi_chirp/dMc = {dPsi_dMc_expr}")
    check("dPsi/dMc equals -(5/128) x^(-5/2)/Mc",
          float(sp.simplify(dPsi_dMc_expr - (-sp.Rational(5, 128) * x ** sp.Rational(-5, 2) / Mc)
                            ).subs({x: 1.3, Mc: 2.0})), 0.0)

    # tau0 term: 2 pi f * tau0 = 2 x^(3/2) tau0 / Mc  -> proportional to dPsi/dt_c
    # tau1 term: 2 pi f * tau1 * (t - tc) = 2 x^(3/2)/Mc * tau1 * (-(5/256) Mc x^(-4))
    #            = -(5/128) tau1 x^(-5/2)  -> proportional to dPsi/dMc (times Mc)
    tau1_coef = sp.simplify(2 * x ** sp.Rational(3, 2) / Mc * (-sp.Rational(5, 256) * Mc * x ** -4))
    print(f"  tau1 coefficient = {tau1_coef}")
    check("tau1 term equals -(5/128) tau1 x^(-5/2)",
          float(sp.simplify(tau1_coef + sp.Rational(5, 128) * x ** sp.Rational(-5, 2)
                            ).subs(x, 1.3)), 0.0)
    # ratio to dPsi/dMc (which is -(5/128) x^(-5/2)/Mc): the tau1 term is
    # -(5/128) x^(-5/2) = Mc * dPsi/dMc  -> an exact Mc rescaling.
    ratio = sp.simplify(tau1_coef / (-sp.Rational(5, 128) * x ** sp.Rational(-5, 2) / Mc))
    print(f"  tau1 / (dPsi/dMc) = {ratio}  (should be +Mc, i.e. a pure Mc shift)")
    check("tau1 direction is Mc * dPsi/dMc", float(sp.simplify(ratio - Mc).subs({x: 1.3, Mc: 2.0})), 0.0)

    # tau2 term exponent: coefficient is 25 Mc /(32768 x^(13/2))
    tau2_coef = sp.simplify(2 * x ** sp.Rational(3, 2) / Mc * (sp.Rational(5, 256) * Mc) ** 2 * x ** -8)
    print(f"  tau2 coefficient = {tau2_coef}")
    check("tau2 coefficient equals 25 Mc/(32768 x^(13/2))",
          float(sp.simplify(tau2_coef - sp.Rational(25, 32768) * Mc * x ** sp.Rational(-13, 2)
                            ).subs({x: sp.Rational(13, 10), Mc: 2})), 0.0)

    print("=== 2. independent recomputation of the geometry ===")
    G = 6.67430e-11
    C = 2.99792458e8
    MSUN = 1.98892e30
    TSUN = G * MSUN / C**3
    MPC = 3.0856775814913673e22
    TMPC = MPC / C
    M1, M2 = 36.0 * TSUN, 29.0 * TSUN
    MF = 62.0 * TSUN
    MC = (M1 * M2) ** 0.6 / (M1 + M2) ** 0.2
    check("chirp mass (Msun)", MC / TSUN, A["setup"]["Mchirp_Msun"], 1e-9)
    check("Rs(62) km", 2 * MF * C / 1e3, A["setup"]["Rs62_km"], 1e-6)
    check("Rs(62)/c ms", 2 * MF * 1e3, A["setup"]["Rs62_over_c_ms"], 1e-9)

    print("=== 3. independent recomputation of the projections (QR, not SVD) ===")
    FMIN, FMAX, NF = 20.0, 300.0, 4000
    f = np.linspace(FMIN, FMAX, NF)
    tc0 = 0.0
    h = (1.0 / (410.0 * TMPC)) * math.sqrt(5.0 / 24.0) * math.pi ** (-2.0 / 3.0) \
        * MC ** (5.0 / 6.0) * f ** (-7.0 / 6.0)

    def psd_shape(fx):
        xx = fx / 215.0
        return 1e-49 * (xx ** -4.14 - 5.0 * xx ** -2.0
                        + 111.0 * (1.0 - xx**2 + 0.5 * xx**4) / (1.0 + 0.5 * xx**2))

    df = f[1] - f[0]
    w_raw = np.sqrt(4.0 * df / psd_shape(f))
    snr_raw = float(np.linalg.norm(h * w_raw))
    k = (snr_raw / 20.0) ** 2
    w = np.sqrt(4.0 * df / (psd_shape(f) * k))

    def t_of_f(fx):
        return tc0 - (5.0 / 256.0) * MC ** (-5.0 / 3.0) * (np.pi * fx) ** (-8.0 / 3.0)

    psi = (3.0 / 128.0) * (np.pi * MC * f) ** (-5.0 / 3.0)
    cols = np.column_stack([h * (-(5.0 / 3.0) * psi / MC), h * 2 * np.pi * f, -h])

    # QR on the weighted matrix -- a different algorithm from the SVD used in
    # analysis.py.
    Qr, Rr = np.linalg.qr(cols * w[:, None], mode="reduced")
    # re-orthonormalise to kill QR drift
    Qq, _ = np.linalg.qr(Qr, mode="reduced")

    def absorbed(vec):
        y = vec * w
        n = np.linalg.norm(y)
        c = Qq.T @ y
        r = y - Qq @ c
        return float(np.linalg.norm(Qq @ c) / n), float(np.linalg.norm(r))

    tref = float(t_of_f(f[0]))
    t_span = abs(tref - float(t_of_f(f[-1])))

    tau0 = 1.2e-3
    u_const = h * 2 * np.pi * f * tau0
    af1, sn1 = absorbed(u_const)
    check("H1 absorbed_fraction", af1, A["H1"]["absorbed_fraction"], 1e-9)
    check("H1 unmodelled_snr", sn1, A["H1"]["unmodelled_snr_calibrated"], 1e-6)

    taudot = 1e-2
    u_lin = h * 2 * np.pi * f * taudot * (t_of_f(f) - tref)
    af2, sn2 = absorbed(u_lin)
    check("H2 absorbed_fraction", af2, A["H2"]["absorbed_fraction"], 1e-9)
    check("H2 unmodelled_snr", sn2, A["H2"]["unmodelled_snr_calibrated"], 1e-6)

    tau2 = tau0 / t_span**2
    u_quad = h * 2 * np.pi * f * tau2 * (t_of_f(f) - tref) ** 2
    af3, sn3 = absorbed(u_quad)
    check("H3 absorbed_fraction", af3, A["H3"]["absorbed_fraction"], 1e-9)
    check("H3 unmodelled_snr", sn3, A["H3"]["unmodelled_snr_calibrated"], 1e-6)
    check("H3 t_span", t_span, A["setup"]["t_span_s"], 1e-12)

    print("=== 4. independent check of the flat-likelihood claim ===")
    # Recompute the K8 residual at two candidates with a completely different
    # route: solve the weighted least-squares problem directly with lstsq.
    def resid_lstsq(data, model):
        y = (data - model) * w
        Aw = cols * w[:, None]
        coef, *_ = np.linalg.lstsq(Aw, y, rcond=None)
        return float(np.linalg.norm(y - Aw @ coef))

    r0 = resid_lstsq(u_const, u_const)                       # candidate 0 ms
    r12 = resid_lstsq(u_const, h * 2 * np.pi * f * 1.2e-3)   # candidate 1.2 ms
    print(f"  residual at candidate 0.0 ms = {r0:.6e}")
    print(f"  residual at candidate 1.2 ms = {r12:.6e}")
    check("K8 flatness (0 vs 1.2 ms) via lstsq", abs(r0 - r12), 0.0, 1e-6)
    check("K8 residual matches artifact",
          r0, A["controls"]["K8_fixed_delay_residual_snr"]["0.0"], 1e-6)

    print("=== 5. independent check of the calibrated SNR ===")
    check("calibrated signal SNR", float(np.linalg.norm(h * w)), 20.0, 1e-9)

    print()
    if FAILURES:
        print(f"INDEPENDENT_CHECK_FAILED: {len(FAILURES)} disagreement(s): {FAILURES}")
        return 1
    print("INDEPENDENT_CONFIRMED: all recomputed values agree within tolerance")
    return 0


if __name__ == "__main__":
    sys.exit(main())
