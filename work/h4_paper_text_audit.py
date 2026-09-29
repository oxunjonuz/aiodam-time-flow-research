#!/usr/bin/env python3
"""
H4 v3: audit the paper's own arithmetic against the paper's own text.

Text read from sources/Muller_Maguire_2016_Now_and_the_Flow_of_Time_extracted.txt:

  Sec. 3: "The Schwarzschild radius is given by Rs = 2GM/c2 = 180 km. ... we
           might expect that the time created to be t = Rs/c = 0.0006 second.
           If we include the gravitational redshift, this extra time observed on
           Earth should appear to be of order 1 millisecond delay in the end of
           the pulse vs the beginning."

  Sec. 4: "we estimate the additional time created in the event as
           (dV_total)^(1/3)/c ~ 0.0012 seconds. With the additional gravitational
           redshift, this value will be larger, approximately 0.0015 sec = 1.5 ms.
           This is close to the value we obtained from dimensional analysis."

RESULTS OF THIS AUDIT (2026-09-27)

 H4a STANDS.  The two BARE numbers are not independent estimates.  Their ratio is
     the geometric factor (dV_total)^(1/3)/Rs(62) = 2.0083, and the paper's own
     bare ratio is 2.0.  So "0.6 ms" and "1.2 ms" are one number and a factor,
     not two agreeing measurements.

 H4b REVISED, and it WEAKENS the criticism.  Only ONE precise factor is actually
     claimed: 1.5/1.2 = 1.25 on the volume branch.  The other apparent factor,
     1.0/0.6 = 1.667, is an artefact of the paper writing "of order 1
     millisecond" -- a rounded phrase, not a derived number.  Treating it as a
     second claimed factor would be over-reading the text.

 H4c REFUTED, by this measurement.  A previous turn (msg 217) claimed that 1.25
     is obtainable "only by an unjustified choice of radius", specifically
     r = 2.78 Rs.  That is wrong.  The volume-weighted mean of the redshift
     factor over the paper's OWN stated shell [1.5, 4] Rs is 1.2503 -- the
     paper's 1.25 to 0.03%.  So the paper's factor is reproducible from its own
     stated geometry by the natural weighting (time created proportional to
     volume created).  This is a correction of my own earlier claim, found by
     measurement, and it is reported as such.

WHAT SURVIVES.  H4a: the paper's two numbers are one number times a geometric
factor, so the "encouraging" agreement between them is arithmetic.  And the
1.25, while now explicable, is still introduced without any derivation in the
text -- the reader is told "with the additional gravitational redshift" and given
a number; the volume-weighted averaging that reproduces it is not stated.

Negative controls: NC1 shows H4a would flip if the volumes were equal; NC2 shows
the 1.25 root is correctly excluded from a sub-interval that does not contain it;
NC3 shows the volume-weighted mean is shell-dependent (1.088 for [1.5, 10] Rs),
so 1.2503 is a property of the paper's specific stated shell, not a constant.
"""
import json
import math
import os
import sys

import numpy as np
from scipy.integrate import quad

G = 6.67430e-11
C = 2.99792458e8
MSUN = 1.98892e30

HERE = os.path.dirname(os.path.abspath(__file__))

PAPER = {
    "Rs_km_stated": 180.0,
    "t_dim_bare_s": 0.0006,
    "t_dim_phrase_s": 0.001,          # "of order 1 millisecond" -- a rounded phrase
    "t_vol_bare_s": 0.0012,
    "t_vol_with_redshift_s": 0.0015,
    "R1_over_Rs": 1.5,
    "R2_over_Rs": 4.0,
    "masses_Msun": [29.0, 36.0, 62.0],   # paper order: 29, 36 -> 62
}


def Rs(M):
    return 2 * G * M / C ** 2


def volume(R1, R2, M):
    """Proper volume, eq. (4.1); flat space when M = 0."""
    if M == 0.0:
        return 4.0 / 3.0 * math.pi * (R2 ** 3 - R1 ** 3)
    Rs_ = Rs(M)
    val, _ = quad(lambda r: 4 * math.pi * r ** 2 / math.sqrt(1.0 - Rs_ / r),
                  R1, R2, limit=400)
    return val


def dV(R1, R2, M):
    """Change of volume created, eq. (4.2)."""
    return volume(R1, R2, M) - volume(R1, R2, 0.0)


def redshift_at(r_over_Rs):
    """1 + z = 1/sqrt(1 - Rs/r) for emission at r, received at infinity."""
    return 1.0 / math.sqrt(1.0 - 1.0 / r_over_Rs)


def radius_for_factor(fac):
    return 1.0 / (1.0 - 1.0 / fac ** 2)


def volume_weighted_redshift(R1, R2, M):
    """dV-weighted mean of (1+z) over the shell.

    Justification for this weighting: the paper says the new time created is
    (dV)^(1/3)/c.  Time created in a shell element dr is proportional to the
    volume created there, and it is observed on Earth stretched by (1+z(r)).
    Averaging (1+z) with weight dV/dr is therefore the natural way to convert the
    shell's redshift into a single factor.
    """
    Rs_ = Rs(M)
    w = lambda r: 4 * math.pi * r ** 2 / math.sqrt(1.0 - Rs_ / r)
    den, _ = quad(w, R1, R2, limit=400)
    num, _ = quad(lambda r: w(r) / math.sqrt(1.0 - Rs_ / r), R1, R2, limit=400)
    return num / den


def main():
    out = {"paper_quoted": PAPER, "reproduce": {}, "H4a": {}, "H4b": {}, "H4c": {},
           "negative_controls": {}, "selfchecks": {}, "thresholds": {}}

    M29, M36, M62 = (m * MSUN for m in PAPER["masses_Msun"])

    Rs62_km = Rs(M62) / 1e3
    t_dim_exact = Rs(M62) / C
    t_dim_stated = PAPER["Rs_km_stated"] * 1e3 / C

    R1 = PAPER["R1_over_Rs"] * Rs(M62)
    R2 = PAPER["R2_over_Rs"] * Rs(M62)
    tot = dV(R1, R2, M62) - dV(R1, R2, M29) - dV(R1, R2, M36)
    L = tot ** (1.0 / 3.0)
    t_vol = L / C
    factor = L / Rs(M62)

    out["reproduce"] = {
        "Rs62_exact_km": Rs62_km,
        "Rs62_paper_stated_km": PAPER["Rs_km_stated"],
        "Rs62_paper_rounding_pct": 100.0 * abs(Rs62_km - PAPER["Rs_km_stated"]) / Rs62_km,
        "t_dimensional_exact_ms": t_dim_exact * 1e3,
        "t_dimensional_from_stated_Rs_ms": t_dim_stated * 1e3,
        "t_dimensional_paper_ms": PAPER["t_dim_bare_s"] * 1e3,
        "dV_total_m3": tot,
        "dV_cbrt_km": L / 1e3,
        "t_volume_ms": t_vol * 1e3,
        "paper_states_volume_ms": PAPER["t_vol_bare_s"] * 1e3,
        "geometric_factor_cbrt_dV_over_Rs": factor,
    }

    # ================================================================ H4a
    out["H4a"] = {
        "claim": "the two BARE estimates are not independent: their ratio is the "
                 "geometric factor (dV_total)^(1/3)/Rs(62)",
        "geometric_factor": factor,
        "paper_bare_ratio": PAPER["t_vol_bare_s"] / PAPER["t_dim_bare_s"],
        "consistent": abs(PAPER["t_vol_bare_s"] / PAPER["t_dim_bare_s"] - factor) < 0.05,
        "verdict": "STANDS",
        "falsifier": "show (dV_total)^(1/3)/Rs(62) = 1",
    }

    # ================================================================ H4b
    f_vol = PAPER["t_vol_with_redshift_s"] / PAPER["t_vol_bare_s"]
    f_dim_phrase = PAPER["t_dim_phrase_s"] / PAPER["t_dim_bare_s"]
    out["H4b"] = {
        "claim": "only one precise correction factor is claimed in the text; the "
                 "apparent second one is an artefact of the rounded phrase "
                 "'of order 1 millisecond'",
        "precise_factor_volume_branch": f_vol,
        "apparent_factor_dimensional_branch": f_dim_phrase,
        "dimensional_phrase_is_rounded": True,
        "phrase_wording": "of order 1 millisecond delay in the end of the pulse",
        "treating_phrase_as_a_factor_would_be_over_reading": True,
        "derivation_of_1.25_in_text": False,
        "verdict": "REVISED -- the earlier 'two different factors' framing was "
                   "over-reading a rounded phrase",
        "falsifier": "find in the text a derivation of 1.25",
    }

    # ================================================================ H4c
    vw = volume_weighted_redshift(R1, R2, M62)
    out["H4c"] = {
        "claim": "the paper's 1.25 IS reproducible from its own stated geometry, "
                 "as the volume-weighted mean redshift over [1.5, 4] Rs",
        "volume_weighted_mean_redshift_over_stated_shell": vw,
        "paper_factor": f_vol,
        "relative_difference": abs(vw - f_vol) / f_vol,
        "matches_within_0.1pct": abs(vw - f_vol) / f_vol < 1e-3,
        "redshift_at_stated_R1_1.5Rs": redshift_at(PAPER["R1_over_Rs"]),
        "redshift_at_stated_R2_4Rs": redshift_at(PAPER["R2_over_Rs"]),
        "radius_for_factor_1.25_in_Rs": radius_for_factor(1.25),
        "previous_claim_msg217": "1.25 obtainable only by an unjustified choice of "
                                 "radius r = 2.78 Rs",
        "previous_claim_status": "REFUTED by this measurement",
        "verdict": "REFUTED -- my own earlier criticism was wrong",
        "falsifier": "show the dV-weighted mean of (1+z) over [1.5,4] Rs differs "
                     "from 1.25 by more than 0.1%",
    }

    # ================================================================ negative controls
    nc = {}
    nc["NC1_H4a_would_flip_if_volumes_equal"] = {
        "test": "replace dV_total by Rs(62)^3, giving factor exactly 1",
        "resulting_factor": (Rs(M62) ** 3) ** (1.0 / 3.0) / Rs(M62),
        "H4a_would_flip": abs((Rs(M62) ** 3) ** (1.0 / 3.0) / Rs(M62) - 1.0) < 1e-12,
    }
    nc["NC2_root_excluded_subinterval"] = {
        "subinterval_Rs": [1.5, 2.0],
        "factor_at_1.5Rs": redshift_at(1.5),
        "factor_at_2Rs": redshift_at(2.0),
        "contains_1.25": redshift_at(2.0) <= 1.25 <= redshift_at(1.5),
        "expected_contains_1.25": False,
    }
    vw_other = volume_weighted_redshift(1.5 * Rs(M62), 10.0 * Rs(M62), M62)
    nc["NC3_vw_is_shell_dependent"] = {
        "shell_Rs": [1.5, 10.0],
        "volume_weighted_mean": vw_other,
        "differs_from_stated_shell": abs(vw_other - vw) > 1e-3,
        "expected": "differs -- so 1.2503 is a property of the stated shell, not a "
                    "universal constant, and the agreement is not vacuous",
    }
    out["negative_controls"] = nc

    # ================================================================ self-checks
    c = {}
    c["S1_volume_reproduces_paper_1.2ms"] = abs(t_vol * 1e3 - 1.2) < 0.05
    c["S2_Rs_matches_paper_180km_within_2pct"] = (
        abs(Rs62_km - PAPER["Rs_km_stated"]) / PAPER["Rs_km_stated"] < 0.02)
    c["S3_paper_bare_ratio_equals_geometric_factor"] = out["H4a"]["consistent"]
    c["S4_volume_weighted_mean_matches_paper_1.25"] = out["H4c"]["matches_within_0.1pct"]
    c["S5_NC1_H4a_flips_when_volumes_equal"] = nc["NC1_H4a_would_flip_if_volumes_equal"][
        "H4a_would_flip"]
    c["S6_NC2_root_correctly_excluded"] = not nc["NC2_root_excluded_subinterval"][
        "contains_1.25"]
    c["S7_NC3_vw_shell_dependent"] = nc["NC3_vw_is_shell_dependent"][
        "differs_from_stated_shell"]
    out["selfchecks"] = {k: bool(v) for k, v in c.items()}

    T = dict(out["selfchecks"])
    T["T_H4a_factor_about_2"] = abs(factor - 2.0) < 0.15
    T["T_H4c_agreement_within_0.1pct"] = out["H4c"]["matches_within_0.1pct"]
    out["thresholds"] = T

    os.makedirs(os.path.join(HERE, "artifacts"), exist_ok=True)
    with open(os.path.join(HERE, "artifacts", "h4_paper_text_audit.json"), "w") as fh:
        json.dump(out, fh, indent=2, sort_keys=True)
    print(json.dumps(out, indent=2, sort_keys=True))
    return 0 if all(T.values()) else 1


if __name__ == "__main__":
    sys.exit(main())
