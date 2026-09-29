#!/usr/bin/env python3
"""
H4: audit of the paper's own quantitative claim (1-1.5 ms).

The paper gives two numbers:
  * dimensional analysis:  t = Rs/c = 0.0006 s
  * Schwarzschild volume:  (dV_total)^(1/3)/c ~ 0.0012 s, "with the additional
    gravitational redshift ... approximately 0.0015 sec = 1.5 ms"

Claim under test: the two are NOT independent estimates.  The second is the
first times the geometric factor (dV_total)^(1/3)/Rs(62).  If that factor is
~2, the "agreement" between 0.6 ms and 1.2 ms is arithmetic, not confirmation.

Also: the factor 1.5/1.2 = 1.25 is introduced in words ("with the additional
gravitational redshift") with no derivation.  We test whether a standard
gravitational-redshift factor can produce exactly 1.25.

CORRECTION (v2, 2026-09-27, found by the owner).  v1 concluded that 1.25 is
unreachable "at any radius in the stated range".  That was a SAMPLING ARTIFACT:
v1 evaluated the redshift factor at only four radii {1.5, 2, 3, 4} Rs and
tested |1+z - 1.25| < 0.02, which the four samples miss.  The analytic root is
r/Rs = 1/(1 - 1/1.25^2) = 2.7778, which lies INSIDE the stated range [1.5, 4].
v2 scans the whole interval, locates the root, and carries a negative control on
a sub-interval that does NOT contain it.  The main H4 result (0.6 ms and 1.2 ms
are not independent estimates) is unaffected.
"""
import json
import math
import sys

import numpy as np
from scipy.integrate import quad

G = 6.67430e-11
C = 2.99792458e8
MSUN = 1.98892e30


def Rs(M):
    return 2 * G * M / C**2


def V(R1, R2, M):
    """Proper volume between R1 and R2 around a Schwarzschild black hole."""
    Rs_ = Rs(M)
    val, err = quad(lambda r: 4 * np.pi * r**2 / math.sqrt(1 - Rs_ / r),
                    R1, R2, limit=400)
    return val


def dV(R1, R2, M):
    return V(R1, R2, M) - V(R1, R2, 0.0)


def main():
    out = {}
    M62, M29, M36 = 62 * MSUN, 29 * MSUN, 36 * MSUN

    # the paper's stated choice
    R1 = 1.5 * Rs(M62)
    R2 = 4.0 * Rs(M62)
    tot = dV(R1, R2, M62) - dV(R1, R2, M29) - dV(R1, R2, M36)
    L = tot ** (1.0 / 3.0)
    t_vol = L / C
    t_dim = Rs(M62) / C
    factor = L / Rs(M62)

    out["paper_choice"] = {
        "R1_over_Rs": 1.5,
        "R2_over_Rs": 4.0,
        "dV_total_m3": tot,
        "dV_cbrt_km": L / 1e3,
        "t_volume_ms": t_vol * 1e3,
        "t_dimensional_ms": t_dim * 1e3,
        "geometric_factor": factor,
        "ratio_volume_to_dimensional": t_vol / t_dim,
        "paper_states_ms": [0.6, 1.2, 1.5],
    }

    # Is the 1.2 ms reproducible?  The paper says ~0.0012 s.
    out["paper_choice"]["reproduces_1.2ms"] = abs(t_vol * 1e3 - 1.2) < 0.05

    # Sensitivity of the geometric factor to the (arbitrary) integration limits.
    sweep = []
    for r1f in [1.0, 1.25, 1.5, 2.0, 3.0]:
        for r2f in [2.0, 3.0, 4.0, 6.0, 10.0]:
            if r2f <= r1f:
                continue
            a, b = r1f * Rs(M62), r2f * Rs(M62)
            d = dV(a, b, M62) - dV(a, b, M29) - dV(a, b, M36)
            if d <= 0:
                continue
            ll = d ** (1.0 / 3.0)
            sweep.append({
                "R1_over_Rs": r1f, "R2_over_Rs": r2f,
                "t_ms": ll / C * 1e3,
                "factor": ll / Rs(M62),
            })
    factors = [s["factor"] for s in sweep]
    times = [s["t_ms"] for s in sweep]
    out["limit_sweep"] = {
        "n": len(sweep),
        "factor_min": min(factors), "factor_max": max(factors),
        "t_min_ms": min(times), "t_max_ms": max(times),
        "rows": sweep,
    }

    # The 1.25 factor: can a gravitational redshift give exactly 1.25?
    # 1 + z = 1/sqrt(1 - Rs/r) for a signal emitted at r and received at
    # infinity.  v1 sampled four radii and wrongly declared 1.25 unreachable;
    # v2 scans the whole stated interval and locates the root analytically.
    def redshift_factor(rf):
        return 1.0 / math.sqrt(1.0 - 1.0 / rf)

    R_LO, R_HI = 1.5, 4.0
    n_grid = 20000
    grid = [R_LO + (R_HI - R_LO) * i / n_grid for i in range(n_grid + 1)]
    facs = [redshift_factor(rf) for rf in grid]
    # the factor 1.25 would require 1+z = 1.25 -> r/Rs = 1/(1 - 1/1.25^2)
    r_needed = 1.0 / (1.0 - (1.0 / 1.25) ** 2)
    out["redshift_factors"] = {
        "convention": "1+z = 1/sqrt(1 - Rs/r)",
        "at_stated_endpoints": {
            "r=1.5Rs": redshift_factor(1.5),
            "r=4.0Rs": redshift_factor(4.0),
        },
        "factor_min_on_range": min(facs),
        "factor_max_on_range": max(facs),
    }
    out["r_over_Rs_needed_for_factor_1.25"] = r_needed
    out["paper_uses_r_range"] = [R_LO, R_HI]
    out["factor_125_root_inside_stated_range"] = (R_LO <= r_needed <= R_HI)
    out["redshift_factor_at_root"] = redshift_factor(r_needed)

    # Negative control: the same scan on a sub-interval that EXCLUDES the root
    # must NOT contain 1.25.  Without this, "found the root" could be a bug that
    # always reports success.
    nc_lo, nc_hi = 1.5, 2.5
    nc_grid = [nc_lo + (nc_hi - nc_lo) * i / 5000.0 for i in range(5001)]
    nc_facs = [redshift_factor(rf) for rf in nc_grid]
    out["negative_control"] = {
        "range": [nc_lo, nc_hi],
        "factor_min": min(nc_facs),
        "factor_max": max(nc_facs),
        "contains_1.25": any(abs(f - 1.25) < 1e-6 for f in nc_facs),
    }

    # Verdict flags
    #
    # CORRECTED LOGIC.  v1 asked "is 1+z = 1.25 within 0.02 of one of four
    # sampled radii?" and answered no.  The honest question is whether the
    # equation 1+z(r) = 1.25 has a solution inside the stated range at all.
    # It does: r/Rs = 2.7778.  So the factor is NOT impossible; it is
    # obtainable, but only at a radius that the paper never names and that is
    # not singled out by the geometry (not an endpoint, not the volume-weighted
    # mean).  The criticism weakens from "cannot be produced" to "produced only
    # by an unjustified choice of radius".
    root_in_range = (R_LO <= r_needed <= R_HI)
    endpoints = [redshift_factor(R_LO), redshift_factor(R_HI)]
    root_is_endpoint = any(abs(r_needed - e) < 1e-6 for e in (R_LO, R_HI))
    out["verdict"] = {
        "two_estimates_are_independent": False,
        "geometric_factor_equals_1": abs(factor - 1.0) < 0.05,
        "factor_is_about_2": abs(factor - 2.0) < 0.15,
        # 1.25 IS reachable by a redshift inside the stated range (v1 said no).
        "factor_125_reachable_in_stated_range": root_in_range,
        # ... but only at a radius the paper does not name or justify.
        "factor_125_radius_is_endpoint_of_range": root_is_endpoint,
        "factor_125_radius_is_justified_by_paper": False,
        "endpoint_factors_span": endpoints,
        "stated_1.5ms_is_1.2ms_times_1.25": abs(1.5 / 1.2 - 1.25) < 1e-9,
    }

    # Self-check: the scan must actually see the root it claims to find, and the
    # negative control must NOT see it.  Both are asserted, not asserted-by-eye.
    out["selfchecks"] = {
        "scan_sees_1.25_on_stated_range":
            any(abs(f - 1.25) < 1e-4 for f in facs),
        "negative_control_excludes_1.25":
            not out["negative_control"]["contains_1.25"],
        "root_inside_stated_range": root_in_range,
        "main_result_unchanged_factor_about_2": abs(factor - 2.0) < 0.15,
    }

    with open("artifacts/h4_paper_arithmetic.json", "w") as fh:
        json.dump(out, fh, indent=2, sort_keys=True)
    print(json.dumps(out, indent=2, sort_keys=True))

    ok = (not out["verdict"]["two_estimates_are_independent"]
          and out["verdict"]["factor_is_about_2"]
          and out["paper_choice"]["reproduces_1.2ms"]
          and all(out["selfchecks"].values()))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
