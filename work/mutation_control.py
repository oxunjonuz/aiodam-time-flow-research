#!/usr/bin/env python3
"""
Hand-rolled mutation control for analysis.py.

The automated campaign (mut_dcd2ff7529de) returned a null result: its coverage
tracer reported "0 executed lines", so it classified every mutant as
"unexecuted" and never ran one. That is an instrument failure, not evidence.
This script does the job directly: apply one textual mutation at a time to a
copy, run the real analysis, and require the exit code to go RED.

A mutation that leaves the exit code GREEN is a hole in the test.
"""
import os
import shutil
import subprocess
import sys

SRC = "analysis.py"
WORK = "mutation_work"

# (label, old_text, new_text) -- each must break a specific claim
MUTATIONS = [
    ("M1 drop the t_c column from the fit basis",
     "    return np.column_stack([h * dMc, h * dtc, h * dphic])",
     "    return np.column_stack([h * dMc, h * dphic])"),
    ("M2 drop the M_c column from the fit basis",
     "    return np.column_stack([h * dMc, h * dtc, h * dphic])",
     "    return np.column_stack([h * dtc, h * dphic])"),
    ("M3 wrong chirp-mass derivative factor (5/3 -> 3/5)",
     "    dMc = -(5.0 / 3.0) * psi_chirp(f, Mc) / Mc",
     "    dMc = -(3.0 / 5.0) * psi_chirp(f, Mc) / Mc"),
    ("M4 wrong t(f) coefficient (5/256 -> 5/128)",
     "    return tc - (5.0 / 256.0) * Mc ** (-5.0 / 3.0) * (np.pi * f) ** (-8.0 / 3.0)",
     "    return tc - (5.0 / 128.0) * Mc ** (-5.0 / 3.0) * (np.pi * f) ** (-8.0 / 3.0)"),
    ("M5 wrong strain amplitude exponent (-7/6 -> -7/5)",
     "        * Mc ** (5.0 / 6.0) * f ** (-7.0 / 6.0)",
     "        * Mc ** (5.0 / 6.0) * f ** (-7.0 / 5.0)"),
    ("M6 wrong chirp phase coefficient (3/128 -> 3/64)",
     "    return (3.0 / 128.0) * (np.pi * Mc * f) ** (-5.0 / 3.0)",
     "    return (3.0 / 64.0) * (np.pi * Mc * f) ** (-5.0 / 3.0)"),
    ("M7 delay phase missing the 2*pi factor",
     "    return 2 * np.pi * f * tau_of_t(t_of_f(f, Mc, tc))",
     "    return f * tau_of_t(t_of_f(f, Mc, tc))"),
    ("M8 PSD not calibrated (drop the scale factor)",
     "    def psd(fx):\n        return psd_shape(fx) * k",
     "    def psd(fx):\n        return psd_shape(fx)"),
    ("M9 projection keeps only 2 of 3 basis directions",
     "    keep = s > s[0] * 1e-12",
     "    keep = s > s[0] * 1e-1"),
    ("M10 absorbed fraction reported as the residual fraction",
     "    return float(np.linalg.norm(Q @ coef) / norm), float(np.linalg.norm(resid))",
     "    return float(np.linalg.norm(resid) / norm), float(np.linalg.norm(resid))"),
    ("M11 H3 uses the linear instead of the quadratic delay",
     "    u_quad = h * dpsi_delay(f, MCHIRP, tc0, lambda t: tau_quad_only(t, tau2, tref))",
     "    u_quad = h * dpsi_delay(f, MCHIRP, tc0, lambda t: tau_linear(t, tau2, tref))"),
    ("M12 K8 grid collapses to a single candidate",
     "    for cand_ms in [0.0, 0.3, 0.6, 1.2, 2.4, 5.0]:",
     "    for cand_ms in [1.2, 1.2, 1.2, 1.2, 1.2, 1.2]:"),
    ("M13 K9 grid collapses to a single candidate",
     "    for cand in [-0.02, -0.01, 0.0, 0.01, 0.02]:",
     "    for cand in [0.01, 0.01, 0.01, 0.01, 0.01]:"),
    ("M14 wrong Rs (drop the factor 2)",
     "    Rs62_m = 2 * MF * C",
     "    Rs62_m = MF * C"),
    ("M15 wrong chirp mass formula (3/5 -> 2/5 exponent)",
     "MCHIRP = (M1 * M2) ** 0.6 / (M1 + M2) ** 0.2",
     "MCHIRP = (M1 * M2) ** 0.4 / (M1 + M2) ** 0.2"),
    ("M16 wrong luminosity distance (410 -> 41)",
     "DL = 410.0 * TMPC",
     "DL = 41.0 * TMPC"),
    ("M17 wrong frequency band (20-300 -> 20-3000)",
     "FMIN, FMAX, NF = 20.0, 300.0, 4000",
     "FMIN, FMAX, NF = 20.0, 3000.0, 4000"),
    ("M18 wrong solar mass in seconds (drop G)",
     "TSUN = G * MSUN / C**3",
     "TSUN = MSUN / C**3"),
    ("M19 sign flip in the delay direction",
     "    return 2 * np.pi * f * tau_of_t(t_of_f(f, Mc, tc))",
     "    return -2 * np.pi * f * tau_of_t(t_of_f(f, Mc, tc))"),
    ("M20 threshold T5 inverted (lt 1 -> gt 1)",
     '    T["T5_H3_paper_scale_snr_lt_1"] = snr3 < 1.0',
     '    T["T5_H3_paper_scale_snr_lt_1"] = snr3 > 1.0'),
]


def main():
    if os.path.isdir(WORK):
        shutil.rmtree(WORK)
    os.makedirs(os.path.join(WORK, "artifacts"))
    src = open(SRC).read()

    # baseline must be green
    with open(os.path.join(WORK, SRC), "w") as fh:
        fh.write(src)
    r = subprocess.run([sys.executable, SRC], cwd=WORK,
                       capture_output=True, text=True)
    print(f"baseline exit = {r.returncode}  (must be 0)")
    if r.returncode != 0:
        print("BASELINE_RED -- cannot score mutations")
        return 2

    killed, survived, invalid = [], [], []
    for label, old, new in MUTATIONS:
        if src.count(old) != 1:
            invalid.append((label, f"anchor count = {src.count(old)}"))
            print(f"  [INVALID] {label}: anchor count = {src.count(old)}")
            continue
        with open(os.path.join(WORK, SRC), "w") as fh:
            fh.write(src.replace(old, new))
        r = subprocess.run([sys.executable, SRC], cwd=WORK,
                           capture_output=True, text=True)
        if r.returncode != 0:
            killed.append(label)
            print(f"  [KILLED ] {label}")
        else:
            survived.append(label)
            print(f"  [SURVIVED] {label}   <-- hole in the test")

    print()
    print(f"killed   = {len(killed)}/{len(MUTATIONS)}")
    print(f"survived = {len(survived)}")
    print(f"invalid  = {len(invalid)}")
    for s in survived:
        print(f"  SURVIVOR: {s}")
    with open(os.path.join(WORK, SRC), "w") as fh:
        fh.write(src)
    return 0 if not survived and not invalid else 1


if __name__ == "__main__":
    sys.exit(main())
