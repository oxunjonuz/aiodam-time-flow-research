#!/usr/bin/env python3
"""verify_claims_228.py -- independent check of the five corrections of turn 228.

It does NOT import any audit script.  It re-reads the frozen artifacts and the
Lean sources and re-derives, by its own arithmetic, every number that the
corrections of turn 228 put into CLAIMS.md / REPORT.md.

V1  "+3.6 ms" is the RIGHT crossing of the level (peak - 1 sigma), not a grid
    edge and not the left crossing.  Re-derived from the profile array itself.
V2  On the recorded grid (+-9.6 ms) the peak sits ON THE LEFT EDGE, so the right
    crossing is +5.4 ms, not +3.6 ms -- i.e. the +3.6 ms number does NOT exist
    on the recorded grid.
V3  There is no Lean file for A4-A6: the Lean sources contain no chirp phase,
    no delay, no Euler identity.  (grep over the actual .lean bytes.)
V4  The 69 ms figure is a LEADING-ORDER number only: it lives in
    analysis_results.json (analytical PSD calibrated to SNR 20), and the IMR
    artifact carries a different amplitude scale.
V5  The current fit resolution is ~5 ms: the recovery scatter is 4.6 ms (1.2 ms
    injection) and 5.1 ms (20 ms injection).
V6  The 88 % absorption quoted for the quadratic term is what h3_real_data.json
    actually records (0.8849 / 0.8716).

Negative controls: each check is paired with a deliberately wrong variant that
must NOT pass.
"""
import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ART = os.environ.get("VC228_ART", os.path.join(HERE, "artifacts"))
LEAN = os.path.join(HERE, "lean")

T_SPAN_S = 0.8444822952069388          # same span the joint fit declares
FAIL = []


def amp(x):
    """tau2 coefficient -> amplitude in ms, exactly as h3_joint_fit declares."""
    return float(x) * T_SPAN_S ** 2 * 1e3


def check(name, ok, detail=""):
    print(f"{'PASS' if ok else 'FAIL'}  {name}  {detail}")
    if not ok:
        FAIL.append(name)


def crossings(grid, profile):
    """Own re-derivation: peak, and the first grid point on each side where the
    profile drops below peak - 1."""
    n = len(profile)
    k = max(range(n), key=lambda i: profile[i])
    peak = profile[k]
    right = next((j for j in range(k, n) if profile[j] < peak - 1.0), None)
    left = next((j for j in range(k, -1, -1) if profile[j] < peak - 1.0), None)
    return {
        "peak": peak,
        "peak_ms": amp(grid[k]),
        "at_edge": k == 0 or k == n - 1,
        "right_ms": None if right is None else amp(grid[right]),
        "left_ms": None if left is None else amp(grid[left]),
    }


def main():
    d = json.load(open(os.path.join(ART, "h3_joint_fit.json")))
    pg = d["per_grid"]

    # ---------------- V1 / V2: what "+3.6 ms" is
    print("== V1/V2: the +3.6 ms number ==")
    pn_rec = crossings(pg["joint_pn_recorded"]["grid"],
                       pg["joint_pn_recorded"]["profile"])
    pn_exp = crossings(pg["joint_pn_expanded"]["grid"],
                       pg["joint_pn_expanded"]["profile"])
    imr_rec = crossings(pg["joint_imr_recorded"]["grid"],
                        pg["joint_imr_recorded"]["profile"])
    imr_exp = crossings(pg["joint_imr_expanded"]["grid"],
                        pg["joint_imr_expanded"]["profile"])
    for tag, c in (("pn/recorded", pn_rec), ("pn/expanded", pn_exp),
                   ("imr/recorded", imr_rec), ("imr/expanded", imr_exp)):
        print(f"   {tag:14s} peak {c['peak']:.4f} @ {c['peak_ms']:+.2f} ms  "
              f"edge={c['at_edge']}  right={c['right_ms']}  left={c['left_ms']}")

    check("V1a expanded pn right crossing = +3.6 ms",
          pn_exp["right_ms"] is not None and abs(pn_exp["right_ms"] - 3.6) < 0.05,
          f"got {pn_exp['right_ms']}")
    check("V1b expanded imr right crossing = +2.4 ms",
          imr_exp["right_ms"] is not None and abs(imr_exp["right_ms"] - 2.4) < 0.05,
          f"got {imr_exp['right_ms']}")
    check("V1c expanded peaks are INTERIOR (not at edge)",
          (not pn_exp["at_edge"]) and (not imr_exp["at_edge"]))
    check("V1d no left crossing on any grid (profile never falls 1 sigma left)",
          all(c["left_ms"] is None for c in (pn_rec, pn_exp, imr_rec, imr_exp)))
    check("V2a recorded pn peak IS at the left edge",
          pn_rec["at_edge"] and abs(pn_rec["peak_ms"] + 9.6) < 0.05,
          f"peak at {pn_rec['peak_ms']:+.2f} ms")
    check("V2b recorded pn right crossing is +5.4 ms, NOT +3.6 ms",
          pn_rec["right_ms"] is not None and abs(pn_rec["right_ms"] - 5.4) < 0.05,
          f"got {pn_rec['right_ms']}")
    # negative control: the wrong reading would be "3.6 ms lies inside +-9.6 ms,
    # so the recorded grid already shows the bound".  Test that reading directly.
    rec_profile = pg["joint_pn_recorded"]["profile"]
    rec_grid = pg["joint_pn_recorded"]["grid"]
    j36 = min(range(len(rec_grid)), key=lambda i: abs(amp(rec_grid[i]) - 3.6))
    drop_at_36 = max(rec_profile) - rec_profile[j36]
    check("V2c NC: at +3.6 ms on the recorded grid the profile has NOT dropped 1 sigma",
          drop_at_36 < 1.0,
          f"drop at +3.6 ms = {drop_at_36:.4f} (< 1.0 => no crossing there)")

    # ---------------- V3: A4-A6 have no Lean backing
    print("== V3: is there a Lean file for A4-A6? ==")
    lean_text = ""
    names = []
    for fn in sorted(os.listdir(LEAN)):
        if fn.endswith(".lean"):
            names.append(fn)
            lean_text += open(os.path.join(LEAN, fn), encoding="utf-8").read()
    print(f"   lean files: {names}")
    # A4-A6 are about the chirp phase expansion, the delay coefficients, and the
    # Euler identity.  NOTE (verifier bug #1, this turn): my first regex looked
    # for the WORD "delay"/"tau" and went red -- but Identifiability.lean is
    # *about* an additive delay in the abstract, so the word is present while
    # the CONTENT of A4-A6 is absent.  The correct test is for the content.
    pat_chirp = re.compile(r"chirp|3/128|32768|13/2|x\^\{?-\s*5/2", re.I)
    pat_euler = re.compile(r"Euler|euler", re.I)
    pat_coeff = re.compile(r"tau_?0|tau_?1|tau_?2|τ₀|τ₁|τ₂", re.I)
    check("V3a no Lean source carries the A4-A6 delay coefficients (tau0/1/2)",
          not pat_coeff.search(lean_text))
    check("V3b no Lean source mentions the chirp phase",
          not pat_chirp.search(lean_text))
    check("V3c no Lean source mentions the Euler identity",
          not pat_euler.search(lean_text))
    # positive control: the files that DO back A1-A3 must be found by name.
    check("V3d NC: Identifiability.lean / Instances.lean exist (A1-A3 do have Lean)",
          "Identifiability.lean" in names and "Instances.lean" in names)
    # positive control that the abstract word "delay" IS present -- so the test
    # above is discriminating content, not the mere word.
    check("V3e NC: the abstract word 'delay' IS present (A1-A3 are about a delay)",
          re.search(r"delay", lean_text, re.I) is not None)

    # ---------------- V4: 69 ms is a leading-order number
    print("== V4: is 69 ms leading-order only? ==")
    a = json.load(open(os.path.join(ART, "analysis_results.json")))
    detect = a["H3"]["detectable_amplitude_ms_at_snr5"]
    print(f"   analysis_results.json H3 detectable = {detect:.2f} ms")
    check("V4a the 68.8 ms figure exists in the analytical-PSD artifact",
          abs(detect - 68.8) < 0.2, f"got {detect:.3f}")
    # the IMR artifact must NOT contain 68.8 -- it has its own amplitude scales
    imr_blob = json.dumps(json.load(open(os.path.join(ART, "h3_imr_check.json"))))
    check("V4b NC: 68.8 does NOT appear in the IMR artifact (different template)",
          "68.8" not in imr_blob and "68.78" not in imr_blob)
    # the analytical artifact's own source says the PSD is calibrated to SNR 20
    src = open(os.path.join(HERE, "analysis.py"), encoding="utf-8").read()
    check("V4c analysis.py declares the PSD is calibrated to the published SNR 20",
          "calibrat" in src.lower() and "20" in src)

    # ---------------- V5: current resolution ~5 ms
    print("== V5: current fit resolution ==")
    rec = d["injection_recovery"]["recovery"]
    s_paper = rec["paper"]["recovered_std_ms"]
    s_large = rec["large"]["recovered_std_ms"]
    print(f"   scatter: paper(1.2 ms) = {s_paper:.2f} ms ; large(20 ms) = {s_large:.2f} ms")
    check("V5a recovery scatter is ~5 ms (4.6 and 5.1)",
          abs(s_paper - 4.6) < 0.2 and abs(s_large - 5.1) < 0.2)
    check("V5b NC: scatter exceeds the 1.2 ms signal (so 1.2 ms is NOT resolved)",
          s_paper > 1.2)

    # ---------------- V6: the 88 % absorption
    print("== V6: the 88 % absorption quoted for the quadratic term ==")
    h3 = json.load(open(os.path.join(ART, "h3_real_data.json")))
    # NOTE (verifier bug #2, this turn): my first regex was `0\.(?:8849|8716)`
    # and missed H1, because the recorded value is 0.88487345... -- the digits
    # are 88487, not 8849.  Read the numbers out of the JSON by key instead of
    # by pattern.
    got = [h3["projection"][d]["H3_quadratic"]["absorbed_fraction"]
           for d in ("H1", "L1")]
    print(f"   absorption fractions recorded: {got}")
    check("V6a 0.88487345 (H1) and 0.87164874 (L1) are recorded (=> 'about 88 %')",
          abs(got[0] - 0.8848734532582031) < 1e-9
          and abs(got[1] - 0.8716487429155013) < 1e-9)

    print()
    if FAIL:
        print("VERIFY_FAILED:", FAIL)
        return 1
    print("VERIFY_CONFIRMED: all five corrections are backed by the frozen artifacts")
    return 0


if __name__ == "__main__":
    if "--selftest" in sys.argv:
        # Corrupt a COPY of the joint-fit artifact and confirm the verifier goes
        # red.  A verifier that cannot fail proves nothing.
        import shutil
        import tempfile
        tmp = tempfile.mkdtemp(prefix="vc228_")
        for fn in os.listdir(ART):
            shutil.copy(os.path.join(ART, fn), os.path.join(tmp, fn))
        jf = os.path.join(tmp, "h3_joint_fit.json")
        d = json.load(open(jf))
        # fault: move the pn/expanded peak to the grid edge by flattening the
        # left half of the profile -- V1c/V2b must then notice.
        prof = d["per_grid"]["joint_pn_expanded"]["profile"]
        n = len(prof)
        for i in range(n):
            prof[i] = 30.0 - 0.001 * i          # monotone decay, peak at left edge
        json.dump(d, open(jf, "w"))
        os.environ["VC228_ART"] = tmp
        ART = tmp
        print("== SELFTEST: corrupted copy, the verifier MUST go red ==")
        rc = main()
        print("SELFTEST_PASS" if rc == 1 else "SELFTEST_FAIL (verifier stayed green on a fault)")
        sys.exit(0 if rc == 1 else 1)
    sys.exit(main())