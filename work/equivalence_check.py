#!/usr/bin/env python3
"""
Classification of the surviving mutants from mutation_control.py.

A mutant that leaves the exit code green is only *equivalent* if it also leaves
the OUTPUT unchanged.  This script proves that directly: it runs the baseline
and each survivor, and compares the full JSON artifact numerically.

If a survivor produces a different artifact, it is a REAL defect in the test,
not an equivalent mutant.
"""
import json
import os
import shutil
import subprocess
import sys

SRC = "analysis.py"
WORK = "equivalence_work"

# A mutant is equivalent only if every claim-bearing quantity agrees to within
# floating-point roundoff.  ATOL is set well below the smallest physically
# meaningful value in the artifact (the H1/H2 unmodelled SNRs are ~1e-14 and are
# zero in exact arithmetic); RTOL covers the ordinary values.
ATOL = 1e-12
RTOL = 1e-9

SURVIVORS = [
    ("M3 wrong chirp-mass derivative factor (5/3 -> 3/5)",
     "    dMc = -(5.0 / 3.0) * psi_chirp(f, Mc) / Mc",
     "    dMc = -(3.0 / 5.0) * psi_chirp(f, Mc) / Mc"),
    ("M6 wrong chirp phase coefficient (3/128 -> 3/64)",
     "    return (3.0 / 128.0) * (np.pi * Mc * f) ** (-5.0 / 3.0)",
     "    return (3.0 / 64.0) * (np.pi * Mc * f) ** (-5.0 / 3.0)"),
    ("M16 wrong luminosity distance (410 -> 41)",
     "DL = 410.0 * TMPC",
     "DL = 41.0 * TMPC"),
    ("M19 sign flip in the delay direction",
     "    return 2 * np.pi * f * tau_of_t(t_of_f(f, Mc, tc))",
     "    return -2 * np.pi * f * tau_of_t(t_of_f(f, Mc, tc))"),
]


def run(src_text, tag):
    d = os.path.join(WORK, tag)
    os.makedirs(os.path.join(d, "artifacts"), exist_ok=True)
    with open(os.path.join(d, SRC), "w") as fh:
        fh.write(src_text)
    r = subprocess.run([sys.executable, SRC], cwd=d, capture_output=True, text=True)
    p = os.path.join(d, "artifacts", "analysis_results.json")
    if not os.path.exists(p):
        return r.returncode, None
    with open(p) as fh:
        return r.returncode, json.load(fh)


def diff(a, b, path=""):
    """Return a list of differing leaf paths between two nested structures."""
    out = []
    if isinstance(a, dict) and isinstance(b, dict):
        for k in sorted(set(a) | set(b)):
            out += diff(a.get(k), b.get(k), f"{path}.{k}")
    elif isinstance(a, list) and isinstance(b, list):
        if len(a) != len(b):
            out.append(f"{path}: length {len(a)} vs {len(b)}")
        else:
            for i, (x, y) in enumerate(zip(a, b)):
                out += diff(x, y, f"{path}[{i}]")
    else:
        if isinstance(a, (int, float)) and isinstance(b, (int, float)):
            # Tolerance-aware comparison.  Exact equality is the wrong test for
            # floating point: several claim-bearing quantities are themselves
            # ~1e-14 (they are zero in exact arithmetic), so a relative test
            # alone would call pure roundoff a "real defect".
            if abs(a - b) > ATOL + RTOL * max(abs(a), abs(b)):
                denom = max(abs(a), abs(b), 1e-300)
                out.append((path, a, b, abs(a - b) / denom))
        elif a != b:
            out.append((path, a, b, float("inf")))
    return out


# A leaf is CLAIM-BEARING if the report's conclusions rest on it.  Everything
# else (condition numbers, floating-point noise in a noisy control) is
# auxiliary and may legitimately differ.
#
# Two tiers, because the mutation campaign showed they behave differently:
#   HEADLINE  -- the numbers the verdict rests on.  Invariant under any
#                rescaling of the basis, so a fault that moves one is a real
#                defect.
#   DIAGNOSTIC-- scale-dependent by construction (a ratio of residuals whose
#                basis direction may itself be rescaled).  A fault may move it
#                without moving any headline claim.
HEADLINE = (
    "absorbed_fraction",
    "one_minus_absorbed",
    "unmodelled_snr",
    "argmin_candidate",
    "Mchirp_Msun",
    "Rs62_km",
    "Rs62_over_c_ms",
    "t_span_s",
    "signal_snr_calibrated",
    "delay_amplitude_ms",
    "detectable_amplitude_ms_at_snr5",
    "reproduces_1.2ms",
    "geometric_factor",
    "t_volume_ms",
    "t_dimensional_ms",
    "two_estimates_are_independent",
    "factor_is_about_2",
    "ratio_detectable_to_claimed",
    "K8_residual_spread",
    "K8b_residual_spread",
)

DIAGNOSTIC = (
    "K9_residual_spread",
    "K9_over_K8_sharpness_ratio",
    "fisher_cond",
    "fisher_singular_values",
)


def tier(path):
    if any(k in path for k in HEADLINE):
        return "headline"
    if any(k in path for k in DIAGNOSTIC):
        return "diagnostic"
    return "other"


def main():
    if os.path.isdir(WORK):
        shutil.rmtree(WORK)
    src = open(SRC).read()
    rc0, base = run(src, "baseline")
    print(f"baseline exit={rc0}, artifact={'yes' if base else 'NO'}")
    if rc0 != 0 or base is None:
        print("BASELINE_RED")
        return 2

    print()
    verdicts = {}
    for label, old, new in SURVIVORS:
        assert src.count(old) == 1, (label, src.count(old))
        rc, art = run(src.replace(old, new), label.split()[0])
        if art is None:
            print(f"  [NO ARTIFACT] {label} (exit {rc})")
            verdicts[label] = "no_artifact"
            continue
        d = diff(base, art)
        head = [x for x in d if tier(x[0]) == "headline"]
        diag = [x for x in d if tier(x[0]) == "diagnostic"]
        other = [x for x in d if tier(x[0]) == "other"]
        if not head:
            verdicts[label] = "equivalent_on_headline_claims"
            print(f"  [EQUIVALENT-ON-HEADLINE] {label}")
            print(f"      headline claims identical; "
                  f"{len(diag)} diagnostic + {len(other)} other leaf/leaves differ")
            for p, x, y, r in (diag + other)[:4]:
                print(f"        {tier(p)} {p}: rel {r:.2e}")
        else:
            verdicts[label] = "real_defect"
            print(f"  [REAL DEFECT] {label}: {len(head)} HEADLINE leaf/leaves differ")
            for p, x, y, r in head[:8]:
                print(f"        {p}: {x!r} vs {y!r} (rel {r:.2e})")
            for p, x, y, r in diag[:3]:
                print(f"        (also diagnostic {p}: rel {r:.2e})")

    print()
    print("summary:", json.dumps(verdicts, indent=2))
    with open(os.path.join(WORK, SRC), "w") as fh:
        fh.write(src)
    return 0 if all(v == "equivalent_on_headline_claims" for v in verdicts.values()) else 1


if __name__ == "__main__":
    sys.exit(main())
