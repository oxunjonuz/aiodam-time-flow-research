#!/usr/bin/env python3
"""H4 issue 3: the 1.25 factor depends on an UNJUSTIFIED CHOICE OF WEIGHT.

Owner's point: h4_paper_text_audit.py averages the redshift over the FULL proper
volume, but its own docstring justifies the weight as the CREATED volume. Those
are different quantities. The full-proper-volume weighting gives 1.2503 (the
paper's 1.25 to 0.03%); the created-volume weighting gives a different number.

This script computes BOTH, plus the flat-volume weighting, on the paper's own
stated shell [1.5, 4] Rs(62), and reports the spread. It does not decide which
is physically right -- it shows that the agreement with the paper's text depends
on the choice, which is exactly what the owner said. The honest conclusion is
therefore weaker than "1.25 is reproduced": it is "1.25 is reproduced by ONE
natural weighting, not by all".

Weights compared (all on [1.5, 4] Rs(62), redshift 1 + z = 1/sqrt(1 - Rs/r)):
  W1 full proper volume element   dV_proper/dr  = 4 pi r^2 / sqrt(1 - Rs/r)
  W2 created volume element       dV_created/dr = 4 pi r^2 (1/sqrt(1-Rs/r) - 1)
  W3 coordinate (flat) volume     dV_flat/dr    = 4 pi r^2

W3 is the null control: with no curvature weight at all, the mean redshift is
still > 1, so the agreement is not automatic -- but it is also not 1.25.
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

R1_OVER_RS, R2_OVER_RS = 1.5, 4.0
PAPER_FACTOR = 0.0015 / 0.0012        # 1.25, the only precise factor in the text


def Rs(M):
    return 2 * G * M / C ** 2


def red(r, Rs_):
    """1 + z at coordinate radius r (metres), for emission received at infinity."""
    return 1.0 / math.sqrt(1.0 - Rs_ / r)


def weighted_mean(weight, a, b, Rs_):
    den, _ = quad(weight, a, b, limit=400)
    num, _ = quad(lambda r: weight(r) * red(r, Rs_), a, b, limit=400)
    return num / den


def main():
    out = {}
    M62 = 62.0 * MSUN
    Rs62 = Rs(M62)
    a, b = R1_OVER_RS * Rs62, R2_OVER_RS * Rs62

    W1 = lambda r: 4 * math.pi * r ** 2 / math.sqrt(1 - Rs62 / r)
    W2 = lambda r: 4 * math.pi * r ** 2 * (1.0 / math.sqrt(1 - Rs62 / r) - 1.0)
    W3 = lambda r: 4 * math.pi * r ** 2

    m1, m2, m3 = (weighted_mean(W, a, b, Rs62) for W in (W1, W2, W3))

    out["paper_factor"] = PAPER_FACTOR
    out["weights"] = {
        "W1_full_proper_volume": {
            "mean_redshift": m1,
            "rel_diff_from_paper": abs(m1 - PAPER_FACTOR) / PAPER_FACTOR},
        "W2_created_volume": {
            "mean_redshift": m2,
            "rel_diff_from_paper": abs(m2 - PAPER_FACTOR) / PAPER_FACTOR},
        "W3_flat_volume": {
            "mean_redshift": m3,
            "rel_diff_from_paper": abs(m3 - PAPER_FACTOR) / PAPER_FACTOR},
    }
    out["spread"] = {
        "max_minus_min": max(m1, m2, m3) - min(m1, m2, m3),
        "relative_to_paper": (max(m1, m2, m3) - min(m1, m2, m3)) / PAPER_FACTOR,
    }
    out["conclusion"] = (
        "The paper's 1.25 is reproduced by the FULL-PROPER-VOLUME weighting "
        f"({m1:.5f}) but NOT by the CREATED-volume weighting ({m2:.5f}) that "
        "h4_paper_text_audit.py's own docstring names as the justification, nor by "
        f"the flat-volume weighting ({m3:.5f}). The agreement with the paper's text "
        "is therefore a property of one particular (natural but unstated) choice of "
        "weight. Claim weakened: '1.25 is reproducible by a natural weighting', not "
        "'1.25 is derived'."
    )

    c = {}
    c["NC1_weights_disagree"] = abs(m1 - m2) > 1e-3
    c["NC2_means_exceed_1"] = m1 > 1.0 and m2 > 1.0 and m3 > 1.0
    c["NC3_flat_weight_is_smallest"] = m3 < m1 and m3 < m2
    out["negative_controls"] = {k: bool(v) for k, v in c.items()}

    T = dict(c)
    T["T_W1_matches_paper_within_0.1pct"] = \
        out["weights"]["W1_full_proper_volume"]["rel_diff_from_paper"] < 1e-3
    T["T_W2_does_not_match_paper"] = \
        out["weights"]["W2_created_volume"]["rel_diff_from_paper"] > 1e-2
    out["thresholds"] = {k: bool(v) for k, v in T.items()}

    os.makedirs(os.path.join(HERE, "artifacts"), exist_ok=True)
    with open(os.path.join(HERE, "artifacts", "h4_weight_audit.json"), "w") as fh:
        json.dump(out, fh, indent=2, sort_keys=True)
    print(json.dumps(out, indent=2, sort_keys=True))
    return 0 if all(T.values()) else 1


if __name__ == "__main__":
    sys.exit(main())
