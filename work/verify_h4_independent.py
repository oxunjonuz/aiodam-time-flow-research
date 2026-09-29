#!/usr/bin/env python3
"""
Independent verification of the CORRECTED H4 redshift sub-claim.

This file imports NOTHING from h4_paper_arithmetic.py.  It reads only the
artifact JSON from disk and recomputes the claim by a different route:

  * the root of 1/sqrt(1 - 1/r) = 1.25 is solved SYMBOLICALLY with sympy
    (the audited script locates it analytically by hand);
  * the reachability question is answered by solving the equation, not by
    scanning a grid;
  * a negative control solves the same equation on [1.5, 2.5] and must find
    no root there.

Disagreement with the artifact beyond tolerance = FAIL.
"""
import json
import sys

import sympy as sp

TOL = 1e-9


def main():
    with open("artifacts/h4_paper_arithmetic.json") as fh:
        art = json.load(fh)

    r = sp.symbols("r", positive=True)
    # 1 + z = 1/sqrt(1 - 1/r) = 1.25  ->  solve exactly
    sols = sp.solve(sp.Eq(1 / sp.sqrt(1 - 1 / r), sp.Rational(125, 100)), r)
    roots = [s for s in sols if s.is_real and s > 1]
    ok = True

    def check(name, got, want, tol):
        nonlocal ok
        good = abs(float(got) - float(want)) <= tol
        ok = ok and good
        print(f"  [{'OK ' if good else 'BAD'}] {name}: got={got!r} want={want!r}")

    print("=== 1. symbolic root of 1+z(r) = 1.25 ===")
    print("    sympy solutions:", sols)
    if len(roots) != 1:
        print("    [BAD] expected exactly one real root > 1")
        ok = False
    else:
        r_root = float(roots[0])
        check("root r/Rs (sympy)", r_root,
              art["r_over_Rs_needed_for_factor_1.25"], TOL)
        check("1+z at root", float(1 / sp.sqrt(1 - 1 / roots[0])), 1.25, TOL)
        lo, hi = art["paper_uses_r_range"]
        reachable = lo <= r_root <= hi
        print(f"    root {r_root:.6f} inside stated range [{lo}, {hi}] -> {reachable}")
        if reachable != art["factor_125_root_inside_stated_range"]:
            print("    [BAD] reachability disagrees with artifact")
            ok = False
        if not reachable:
            print("    [BAD] root must be inside the stated range")
            ok = False

    print("=== 2. negative control: no root on [1.5, 2.5] ===")
    nc_lo, nc_hi = 1.5, 2.5
    # 1+z is strictly decreasing in r, so its range on [a,b] is [1+z(b), 1+z(a)].
    f_hi = float(1 / sp.sqrt(1 - 1 / sp.Rational(15, 10)))   # at r=1.5
    f_lo = float(1 / sp.sqrt(1 - 1 / sp.Rational(25, 10)))   # at r=2.5
    contains = f_lo <= 1.25 <= f_hi
    print(f"    1+z on [{nc_lo}, {nc_hi}] spans [{f_lo:.6f}, {f_hi:.6f}]")
    print(f"    contains 1.25 -> {contains}")
    if contains != art["negative_control"]["contains_1.25"]:
        print("    [BAD] negative control disagrees with artifact")
        ok = False
    if contains:
        print("    [BAD] negative control must NOT contain the root")
        ok = False

    print("=== 3. main H4 result must be unchanged ===")
    check("geometric factor ~2", art["paper_choice"]["geometric_factor"],
          2.0083004031215705, 1e-9)
    check("t_volume_ms", art["paper_choice"]["t_volume_ms"],
          1.2269060002076415, 1e-9)

    print()
    if ok:
        print("H4_INDEPENDENT_CONFIRMED: corrected sub-claim agrees by an "
              "independent symbolic route")
        return 0
    print("H4_INDEPENDENT_DISAGREES")
    return 1


if __name__ == "__main__":
    sys.exit(main())
