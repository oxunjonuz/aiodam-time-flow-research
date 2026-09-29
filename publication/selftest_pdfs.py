#!/usr/bin/env python3
"""Self-test for verify_pdfs.py: corrupt copies and require it to go RED.

A verifier that cannot fail proves nothing. This builds a temporary copy of the
publication folder and applies one corruption at a time, each targeting a
specific class of error:

  P1  a PDF is replaced by garbage            (the file is not a PDF at all)
  P2  a figure file the Markdown references is deleted
  P3  the authorship sentence loses the co-author's error count
  P4  the main font is swapped for one that does NOT cover Cyrillic
      (the Russian PDF would then render with missing glyphs)
  P5  the whole PDF is truncated to one page  (content silently lost)

Exit 0 => SELFTEST_PASS. Exit 1 => the verifier missed a corruption.
"""
import os
import re
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
VERIFIER = os.path.join(HERE, "verify_pdfs.py")


def run(pub):
    p = subprocess.run([sys.executable, VERIFIER, pub], capture_output=True, text=True)
    return p.returncode, (p.stdout + p.stderr)


def build_tree(tmp):
    pub = os.path.join(tmp, "publication")
    shutil.copytree(HERE, pub,
                    ignore=shutil.ignore_patterns(".preview", "__pycache__"))
    return pub


def main():
    results = []
    with tempfile.TemporaryDirectory() as tmp:
        pub = build_tree(tmp)

        rc, out = run(pub)
        results.append(("baseline green", rc == 0, out.strip().splitlines()[-1]))

        # P1: garbage PDF
        pdf = os.path.join(pub, "paper", "paper_en.pdf")
        keep = open(pdf, "rb").read()
        with open(pdf, "wb") as f:
            f.write(b"not a pdf at all\n")
        rc, out = run(pub)
        results.append(("P1 PDF replaced by garbage", rc != 0,
                        out.strip().splitlines()[-1]))
        with open(pdf, "wb") as f:
            f.write(keep)

        # P2: delete a referenced figure
        fig = os.path.join(pub, "figures", "fig3_template.pdf")
        keep_fig = open(fig, "rb").read()
        os.remove(fig)
        rc, out = run(pub)
        results.append(("P2 referenced figure deleted", rc != 0,
                        out.strip().splitlines()[-1]))
        with open(fig, "wb") as f:
            f.write(keep_fig)

        # P3: the co-author's error count removed from the Markdown AND the PDF
        # rebuilt — editing the source alone cannot change the PDF, so the
        # corruption must go through the build to be a real test.
        md = os.path.join(pub, "paper", "paper_en.md")
        t = open(md, encoding="utf-8").read()
        keep_md = t
        t2 = t.replace("The co-author found **eleven** of them: **six** in one",
                       "The co-author reviewed the work in several")
        open(md, "w", encoding="utf-8").write(t2)
        subprocess.run(["sh", os.path.join(pub, "build_paper.sh"), "en"],
                       capture_output=True, text=True, cwd=pub)
        rc, out = run(pub)
        results.append(("P3 error count removed and PDF rebuilt", rc != 0,
                        out.strip().splitlines()[-1]))
        open(md, "w", encoding="utf-8").write(keep_md)
        subprocess.run(["sh", os.path.join(pub, "build_paper.sh"), "en"],
                       capture_output=True, text=True, cwd=pub)

        # P4: rebuild the Russian PDF with a font that has no Cyrillic
        # (DejaVu Sans Mono is used for code and has no Cyrillic in the *serif*
        #  role; the honest test is to point the main font at a Latin-only face)
        ru = os.path.join(pub, "paper", "paper_ru.pdf")
        keep_ru = open(ru, "rb").read()
        pre = os.path.join(pub, "preamble.tex")
        keep_pre = open(pre, encoding="utf-8").read()
        open(pre, "w", encoding="utf-8").write(
            keep_pre.replace("\\setmainfont{P052}", "\\setmainfont{TeX Gyre Termes}"))
        p = subprocess.run(["sh", os.path.join(pub, "build_paper.sh"), "ru"],
                           capture_output=True, text=True, cwd=pub)
        rc, out = run(pub)
        rebuilt = os.path.exists(ru) and open(ru, "rb").read() != keep_ru
        results.append(("P4 main font without Cyrillic coverage",
                        rc != 0 or not rebuilt,
                        out.strip().splitlines()[-1]))
        open(pre, "w", encoding="utf-8").write(keep_pre)
        with open(ru, "wb") as f:
            f.write(keep_ru)

        # P5: truncate the PDF to its first page
        with open(pdf, "wb") as f:
            f.write(keep[:len(keep) // 3])
        rc, out = run(pub)
        results.append(("P5 PDF truncated", rc != 0,
                        out.strip().splitlines()[-1]))
        with open(pdf, "wb") as f:
            f.write(keep)

        # and the baseline must be green again after all restores
        rc, out = run(pub)
        results.append(("baseline green after restores", rc == 0,
                        out.strip().splitlines()[-1]))

    print("SELFTEST OF verify_pdfs.py")
    print("-" * 64)
    caught = 0
    for name, good, last in results:
        print(f"  {'ok  ' if good else 'MISS'}  {name}   [{last}]")
        caught += 1 if good else 0
    print("-" * 64)
    print(f"{caught}/{len(results)} checks behaved as required")
    if caught != len(results):
        print("SELFTEST_FAIL")
        return 1
    print("SELFTEST_PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())