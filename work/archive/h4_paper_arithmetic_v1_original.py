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
    # z = 1/sqrt(1 - Rs/r) - 1.  For the emission region r ~ 1.5-4 Rs.
    zs = {}
    for rf in [1.5, 2.0, 3.0, 4.0]:
        z = 1.0 / math.sqrt(1.0 - 1.0 / rf) - 1.0
        zs[f"r={rf}Rs"] = z
    out["redshift_factors"] = zs
    # the factor 1.25 would require 1+z = 1.25 -> z = 0.25 -> r/Rs = ?
    r_needed = 1.0 / (1.0 - (1.0 / 1.25) ** 2)
    out["r_over_Rs_needed_for_factor_1.25"] = r_needed
    out["paper_uses_r_range"] = [1.5, 4.0]

    # Verdict flags
    out["verdict"] = {
        "two_estimates_are_independent": False,
        "geometric_factor_equals_1": abs(factor - 1.0) < 0.05,
        "factor_is_about_2": abs(factor - 2.0) < 0.15,
        "factor_125_derivable_from_redshift_in_stated_range":
            any(abs((1 + z) - 1.25) < 0.02 for z in zs.values()),
        "stated_1.5ms_is_1.2ms_times_1.25": abs(1.5 / 1.2 - 1.25) < 1e-9,
    }

    with open("artifacts/h4_paper_arithmetic.json", "w") as fh:
        json.dump(out, fh, indent=2, sort_keys=True)
    print(json.dumps(out, indent=2, sort_keys=True))

    ok = (not out["verdict"]["two_estimates_are_independent"]
          and out["verdict"]["factor_is_about_2"]
          and out["paper_choice"]["reproduces_1.2ms"])
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
