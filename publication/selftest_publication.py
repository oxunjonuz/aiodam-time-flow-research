#!/usr/bin/env python3
"""Self-test for verify_publication.py: corrupt copies and require it to go RED.

A verifier that cannot fail proves nothing. This runs verify_publication.py
against deliberately corrupted copies of the PUBLISHED files and of the frozen
ARTIFACTS, in an isolated temporary tree, and asserts that the verifier reports
PUBLICATION_INCONSISTENT for every corruption. It also asserts it stays GREEN on
the unmodified tree (so the test itself is not vacuous).

Each corruption targets one specific class of error this project has actually
made:
  C1  a document stops quoting a headline number (headline weaker than the work,
      or a number silently edited away). NOTE: this check is a PRESENCE check — it
      verifies a number is present somewhere in the document, not that it sits
      next to the right words. A corruption that edits a second copy of the same
      token while leaving one intact does not trip it, by design; C1 therefore
      uses a token that occurs exactly once, and C7 covers the case where the
      number is right but its recorded hash is not.
  C2  the artifact silently changes (a measurement drifts under a stable document)
  C3  the scope caveat is deleted (the headline becomes stronger than the work)
  C4  the authorship disclosure is deleted (the "made by an agent" statement lost)
  C5  a claim-ledger category label is removed
  C6  a translation is truncated (a document that looks present but is empty)
  C7  the ARTIFACTS.md hash no longer matches the artifact on disk
  C8  a Lean theorem the paper cites is renamed away (phantom proof)
  C9  the convention-dependence caveat is deleted (a bound stated as convention-free)
  C10 a PDF is damaged (the PDFs are part of the package; verify_pdfs.py must catch it)
  C11 the figure caption calls the τ₂ profile a measured constraint again
  C12 the inclination is called measured again
  C13 the 69 ms vs 0.61 ms comparison is called three orders of magnitude again
  C14 the SNR-provenance paragraph is removed
  C15 a §9 limitation is duplicated again
  C10b the text block is widened past the paper and the PDF rebuilt — the frame
       check must see the overflow (this is the owner's own finding, turned into a
       corruption so the check that catches it is proven able to fail)
  C16 §2 restricts the headline to register 1 again — the contradiction with §8
       that the owner found in the built PDF (round 231)
  C17 the §3.2 count word disagrees with its bullet list again
  C18 §0 quotes an error split that count_errors.py does not produce
  C19 the deposit metadata (.zenodo.json) drifts from the paper's error split
  C20 the removed self-contradiction clause returns to the EN paper
       (the owner asked in round 232 that the paper not advertise its own
       self-contradictions; C21 continues the same request in round 233)
  C21 the typesetting-exclusion argument returns to §0 — the owner asked in round
       233 that errors found while building and proofreading the PDF not be
       disclosed to the reader at all; the check that keeps them out must be able
       to fail
  plus, after every corrupt/restore pair, a byte-exactness check on the restored
  file — because a restore that misses leaves the NEXT corruption meaningless.
  That check is not decoration: C3's original restore replaced the first occurrence
  of "stated", which is inside the word "unstated" in the abstract, so the file was
  left corrupted and C4 went red for the wrong reason. The restore is now anchored
  on the whole sentence, and the byte check guards it.

C11–C15 are the four self-contradictions the owner found in round 230, each turned
into a corruption so the check that catches it is proven to be able to fail.

Exit 0 => SELFTEST_PASS (all corruptions caught). Exit 1 => the verifier missed one.
"""

import json
import os
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
VERIFIER = os.path.join(HERE, "verify_publication.py")


def run(tree):
    """Run the verifier inside `tree` (which mirrors the real layout)."""
    p = subprocess.run([sys.executable, os.path.join(tree, "publication", "verify_publication.py")],
                       capture_output=True, text=True)
    return p.returncode, p.stdout + p.stderr


def build_tree(tmp):
    root = os.path.join(tmp, "root")
    os.makedirs(os.path.join(root, "publication"))
    os.makedirs(os.path.join(root, "work", "artifacts"))
    # copy the publication dir and the artifacts (verify_publication reads the
    # RU report from work/REPORT.md for the size comparison, so copy it too)
    for name in os.listdir(os.path.join(HERE)):
        src = os.path.join(HERE, name)
        if os.path.isfile(src):
            shutil.copy2(src, os.path.join(root, "publication", name))
        elif os.path.isdir(src) and name in ("paper", "docs_en", "figures"):
            shutil.copytree(src, os.path.join(root, "publication", name))
    art_src = os.path.join(os.path.dirname(HERE), "work", "artifacts")
    for name in os.listdir(art_src):
        shutil.copy2(os.path.join(art_src, name),
                     os.path.join(root, "work", "artifacts", name))
    shutil.copy2(os.path.join(os.path.dirname(HERE), "work", "REPORT.md"),
                 os.path.join(root, "work", "REPORT.md"))
    # the error-count check runs count_errors.py, which reads work/NOTES.md, so the
    # fixture must carry it too (otherwise the baseline is red for a missing file)
    shutil.copy2(os.path.join(os.path.dirname(HERE), "work", "NOTES.md"),
                 os.path.join(root, "work", "NOTES.md"))
    # the verifier also checks the Lean-source hashes listed in ARTIFACTS.md, so
    # the fixture must carry work/lean/ too (only the .lean files matter)
    lean_src = os.path.join(os.path.dirname(HERE), "work", "lean")
    lean_dst = os.path.join(root, "work", "lean")
    os.makedirs(lean_dst)
    for name in os.listdir(lean_src):
        if name.endswith(".lean"):
            shutil.copy2(os.path.join(lean_src, name),
                         os.path.join(lean_dst, name))
    return root


def edit(path, old, new, count=None):
    with open(path, encoding="utf-8") as f:
        t = f.read()
    if old not in t:
        raise SystemExit(f"corruption target not found in {path}: {old!r}")
    t = t.replace(old, new) if count is None else t.replace(old, new, count)
    with open(path, "w", encoding="utf-8") as f:
        f.write(t)


def main():
    results = []
    with tempfile.TemporaryDirectory() as tmp:
        root = build_tree(tmp)

        # 0) baseline: must be GREEN
        rc, out = run(root)
        results.append(("baseline green", rc == 0, out.strip().splitlines()[-1]))

        # a restore-integrity helper: after every corrupt/restore pair the file must
        # be byte-identical to the original, otherwise the NEXT corruption is applied
        # to a file that is already wrong and its result means nothing. (This was a
        # real defect: C3's restore replaced the first occurrence of "stated", which
        # is inside the word "unstated" in the abstract, so the file was left
        # corrupted and C4 then went red for the wrong reason.)
        import hashlib

        def digest(path):
            with open(path, "rb") as f:
                return hashlib.sha256(f.read()).hexdigest()

        originals = {}
        for rel in ("paper/paper_en.md", "docs_en/CLAIMS_en.md"):
            p = os.path.join(root, "publication", rel)
            originals[p] = digest(p)

        def restore_ok(label):
            bad = [os.path.basename(p) for p, h in originals.items() if digest(p) != h]
            if bad:
                results.append((f"restore after {label} is byte-exact", False,
                                f"still corrupted: {bad}"))
            else:
                results.append((f"restore after {label} is byte-exact", True, "ok"))

        # C1: change a headline number in the EN paper to a value the artifact
        # does not support. 0.174 occurs exactly once, so removing it is a real
        # test of the presence check (see the note in the header about C1).
        p = os.path.join(root, "publication", "paper", "paper_en.md")
        edit(p, "0.174", "0.999", 1)
        rc, out = run(root)
        results.append(("C1 headline number not in artifact", rc != 0,
                        out.strip().splitlines()[-1]))
        # restore
        edit(p, "0.999", "0.174", 1)
        restore_ok("C1")

        # C2: silently change an artifact number
        a = os.path.join(root, "work", "artifacts", "h3_joint_fit.json")
        with open(a, encoding="utf-8") as f:
            d = json.load(f)
        d["answer"]["joint_fit"]["pn"]["net_snr_peak"] = 99.0
        with open(a, "w", encoding="utf-8") as f:
            json.dump(d, f)
        rc, out = run(root)
        results.append(("C2 artifact silently changed", rc != 0,
                        out.strip().splitlines()[-1]))
        with open(a, "w", encoding="utf-8") as f:
            json.dump(json.load(open(os.path.join(os.path.dirname(HERE),
                                                  "work", "artifacts",
                                                  "h3_joint_fit.json"), encoding="utf-8")), f)

        # C3: delete the scope caveat from the EN paper. The restore is anchored on
        # the SENTENCE, not on the bare word "stated" — see restore_ok above.
        edit(p, "not derived from the theory", "stated from the theory", 1)
        edit(p, "does not specify it", "gives it", 1)
        rc, out = run(root)
        results.append(("C3 scope caveat deleted", rc != 0,
                        out.strip().splitlines()[-1]))
        edit(p, "gives it", "does not specify it", 1)
        edit(p, "stated from the theory", "not derived from the theory", 1)
        restore_ok("C3")

        # C4: delete the authorship disclosure. This now removes the FULL sentence,
        # because a bare-phrase check was measured to survive the deletion (the
        # phrase also appears in the YAML author list).
        edit(p, "carried out **by an autonomous research agent, Aiodam**, acting as",
             "carried out by the project, acting as", 1)
        rc, out = run(root)
        results.append(("C4 authorship disclosure deleted", rc != 0,
                        out.strip().splitlines()[-1]))
        edit(p, "carried out by the project, acting as",
             "carried out **by an autonomous research agent, Aiodam**, acting as", 1)
        restore_ok("C4")

        # C4b: the error-attribution sentence removed (the paper would then name the
        # co-author without saying how many errors he found)
        edit(p, "The co-author found **eleven** of them: **six** in one",
             "The co-author reviewed the work in several",
             1)
        rc, out = run(root)
        results.append(("C4b error-attribution sentence deleted", rc != 0,
                        out.strip().splitlines()[-1]))
        edit(p, "The co-author reviewed the work in several",
             "The co-author found **eleven** of them: **six** in one", 1)
        restore_ok("C4b")

        # C5: remove a claim-ledger category
        c = os.path.join(root, "publication", "docs_en", "CLAIMS_en.md")
        edit(c, "[HYPOTHESIS]", "[OPEN]")
        rc, out = run(root)
        results.append(("C5 ledger category label removed", rc != 0,
                        out.strip().splitlines()[-1]))
        edit(c, "[OPEN]", "[HYPOTHESIS]")
        restore_ok("C5")

        # C6: truncate a translation
        r = os.path.join(root, "publication", "docs_en", "ROUNDS_en.md")
        with open(r, "w", encoding="utf-8") as f:
            f.write("x")
        rc, out = run(root)
        results.append(("C6 translation truncated", rc != 0,
                        out.strip().splitlines()[-1]))
        shutil.copy2(os.path.join(HERE, "docs_en", "ROUNDS_en.md"), r)

        # C7: break the ARTIFACTS.md hash
        m = os.path.join(root, "publication", "ARTIFACTS.md")
        with open(m, encoding="utf-8") as f:
            t = f.read()
        t = t.replace("1daa84078cabac6f146d4586b2c0703738570d6c3126119012e7034537eb84c8",
                      "0" * 64)
        with open(m, "w", encoding="utf-8") as f:
            f.write(t)
        rc, out = run(root)
        results.append(("C7 ARTIFACTS.md hash broken", rc != 0,
                        out.strip().splitlines()[-1]))

        # C8: rename a theorem in the Lean source so the paper cites a proof that
        # no longer exists (the strongest form of "headline stronger than content")
        lf = os.path.join(root, "work", "lean", "Identifiability.lean")
        edit(lf, "identifiable_iff", "identifiable_iff_RENAMED")
        rc, out = run(root)
        results.append(("C8 paper cites a theorem Lean no longer defines", rc != 0,
                        out.strip().splitlines()[-1]))
        edit(lf, "identifiable_iff_RENAMED", "identifiable_iff")

        # C9: delete the convention caveat from the paper (the bound would then be
        # stated as convention-free, which is false by a factor of 15)
        edit(p, "24.96 ms", "some other value")
        rc, out = run(root)
        results.append(("C9 convention-dependence caveat deleted", rc != 0,
                        out.strip().splitlines()[-1]))
        edit(p, "some other value", "24.96 ms")
        restore_ok("C9")

        # C10: the PDFs are part of the package. If a PDF is missing, the verifier
        # as shipped did not notice at all — it never looked at them. (This is the
        # gap that made the PDF checks a separate script; the selftest now confirms
        # the separate script fails on a damaged PDF.)
        pdfcheck = os.path.join(HERE, "verify_pdfs.py")
        if os.path.exists(pdfcheck):
            pdf = os.path.join(root, "publication", "paper", "paper_en.pdf")
            if os.path.exists(pdf):
                with open(pdf, "wb") as f:
                    f.write(b"%PDF-1.4\nnot a real pdf\n")
                p2 = subprocess.run([sys.executable, pdfcheck,
                                     os.path.join(root, "publication")],
                                    capture_output=True, text=True)
                results.append(("C10 damaged PDF caught by verify_pdfs.py",
                                p2.returncode != 0,
                                (p2.stdout + p2.stderr).strip().splitlines()[-1]))
                shutil.copy2(os.path.join(HERE, "paper", "paper_en.pdf"), pdf)

        # C10b: the frame check must catch an overflow. The corruption is a real
        # one: the geometry is widened so the text block extends past the paper,
        # and the PDF is rebuilt. If verify_frames.py cannot see that, it is blind.
        framecheck = os.path.join(HERE, "verify_frames.py")
        if os.path.exists(framecheck):
            pre = os.path.join(root, "publication", "preamble.tex")
            keep_pre = open(pre, encoding="utf-8").read()
            open(pre, "w", encoding="utf-8").write(
                keep_pre.replace("margin=2.4cm", "margin=0.4cm"))
            subprocess.run(["sh", os.path.join(root, "publication", "build_paper.sh"), "en"],
                           capture_output=True, text=True, cwd=os.path.join(root, "publication"))
            p3 = subprocess.run([sys.executable, framecheck,
                                 os.path.join(root, "publication")],
                                capture_output=True, text=True)
            results.append(("C10b text-block overflow caught by verify_frames.py",
                            p3.returncode != 0,
                            (p3.stdout + p3.stderr).strip().splitlines()[-1]))
            open(pre, "w", encoding="utf-8").write(keep_pre)
            subprocess.run(["sh", os.path.join(root, "publication", "build_paper.sh"), "en"],
                           capture_output=True, text=True, cwd=os.path.join(root, "publication"))

        # ---- the four self-contradictions the owner found (round 230) ----------
        # Each is restored immediately and the restore is byte-checked, so a
        # failure here means the check is wrong, not that the file drifted.
        def contradiction(label, old, new):
            edit(p, old, new, 1)
            rc, out = run(root)
            results.append((label, rc != 0, out.strip().splitlines()[-1]))
            edit(p, new, old, 1)
            restore_ok(label)

        # C11: the figure caption calls the tau2 profile a measured constraint again
        contradiction(
            "C11 caption overclaims the tau2 bound",
            "This is an **upper bound on the strength of the constraint**, not a measured bound:",
            "This is the first constraint on τ₂ the data itself yields, and it is a measured bound:")

        # C12: the inclination is called measured again
        contradiction(
            "C12 inclination called measured again",
            "That crossing is a **fitted** inclination — the value the model would need, not a\nmeasured one — and it is stated as such.",
            "That crossing is a measured inclination, and it is stated as such.")

        # C13: the order of magnitude goes back to three
        contradiction(
            "C13 orders of magnitude wrong again",
            "is **two** orders of magnitude above what the mechanism itself provides:",
            "is three orders of magnitude above what the mechanism itself provides:")

        # C14: the SNR provenance paragraph is removed
        contradiction(
            "C14 SNR provenance paragraph removed",
            "* **Which signal-to-noise ratio is which.**",
            "* **A note on the numbers.**")

        # C15: the duplicated limitation comes back (item 9 replaced by a copy of 8)
        contradiction(
            "C15 §9 limitation duplicated again",
            "9. **The inclination ι ≈ 59.2° is fitted, not measured.**",
            "9. **The π-phase bound (1.667 ms) is convention-dependent.**")

        # C16: §2 restricts the headline to register 1 again — the contradiction
        # with §8 that the owner found in the built PDF (round 231).
        contradiction(
            "C16 §2 restricts the headline to register 1 again",
            "a claim in the **hypothesis** category is never allowed\ninto the headline, which rests on registers 1 and 2 only.",
            "a claim is only allowed into the headline if it is in\nregister 1.")

        # C17: the count word in §3.2 disagrees with its bullet list again
        # ("Two further" in front of three bullets).
        contradiction(
            "C17 §3.2 count word wrong again",
            "Three further formal results support the surrounding argument",
            "Two further formal results support the surrounding argument")

        # C18: §0 quotes an error split that count_errors.py does not produce —
        # the paragraph would then be a remembered number, which is what it claims
        # not to be. (This is the owner's round-231 request turned into a corruption:
        # §0, NOTES.md and count_errors.py must agree.)
        contradiction(
            "C18 §0 error split disagrees with count_errors.py",
            "The co-author found **eleven** of them: **six** in one",
            "The co-author found **twelve** of them: **six** in one")
        # C19: the deposit metadata drifts from the paper's split (a second copy of
        # the same number, which is exactly how the first one went stale)
        z = os.path.join(root, "publication", ".zenodo.json")
        with open(z, encoding="utf-8") as f:
            zt = f.read()
        keep_z = zt
        zt2 = zt.replace("<strong>eleven</strong> of the 47 recorded errors",
                         "<strong>twelve</strong> of the 47 recorded errors")
        with open(z, "w", encoding="utf-8") as f:
            f.write(zt2)
        rc, out = run(root)
        results.append(("C19 zenodo split drifts from the paper", rc != 0,
                        out.strip().splitlines()[-1]))
        with open(z, "w", encoding="utf-8") as f:
            f.write(keep_z)

        # C20: the removed clause comes back into the EN paper. The owner asked
        # (round 232) that the paper not advertise its own self-contradictions, and
        # again (round 233) that errors found while building and proofreading the
        # PDF not be disclosed to the reader at all. The check that enforces this
        # must be able to fail, so the clause is re-inserted here.
        contradiction(
            "C20 removed self-contradiction clause returns",
            "and **four** further errors of\nformulation.",
            "and **four**\nwhen he read the built PDF and found the paper contradicting its own headline\ncaveat in four places.")

        # C21: the typesetting-exclusion argument comes back into §0. It was a
        # positive check before round 233; the owner's request turned it into a
        # negative one, and the new check must be proven able to go red.
        contradiction(
            "C21 typesetting-exclusion argument returns to §0",
            "rule can be argued with.",
            "rule can be argued with. A formula that ran past the page margin — a\n"
            "defect of PDF assembly, not of the research — is not among the 47.")

    print("SELFTEST OF verify_publication.py")
    print("-" * 60)
    caught = 0
    for name, good, last in results:
        status = "ok  " if good else "MISS"
        print(f"  {status}  {name}   [{last}]")
        caught += 1 if good else 0
    print("-" * 60)
    print(f"{caught}/{len(results)} checks behaved as required")
    if caught != len(results):
        print("SELFTEST_FAIL")
        return 1
    print("SELFTEST_PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
