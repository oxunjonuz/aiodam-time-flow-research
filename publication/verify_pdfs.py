#!/usr/bin/env python3
"""Independent check of the two built PDFs.

Does NOT read the Markdown sources for its verdict: it reads the PDFs themselves.
Checks:
  * both PDFs parse, are non-trivial, and are A4;
  * every font used is EMBEDDED (an unembedded font renders differently on the
    reader's machine — the whole point of choosing a font);
  * the glyphs the main font lacks are actually PRESENT in the extracted text
    (a missing glyph would silently drop out of the PDF);
  * every figure file referenced by the Markdown is present on disk AND its
    caption text reaches the PDF;
  * the authorship paragraph names both people and the number of errors the
    co-author found, and the count matches count_errors.py.

Exit 0 => PDFS_OK.
"""
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
# an optional first argument overrides the package root, so the self-test can run
# this script against a deliberately damaged copy in a temporary tree
ROOT = os.path.abspath(sys.argv[1]) if len(sys.argv) > 1 else HERE
PAPER = os.path.join(ROOT, "paper")
FIG = os.path.join(ROOT, "figures")

ok, bad = [], []


def check(name, cond, detail=""):
    (ok if cond else bad).append(name if cond else f"{name}: {detail}")


def run(cmd):
    p = subprocess.run(cmd, capture_output=True, text=True)
    return p.returncode, p.stdout + p.stderr


# the co-author's error count, taken from the counter rather than typed here.
# count_errors.py reads work/NOTES.md, which the temporary tree does not carry, so
# if it is absent the counts are taken from the real package (the PDFs are what is
# under test here, not the notes).
counter = os.path.join(HERE, "count_errors.py")
if not os.path.exists(os.path.join(os.path.dirname(ROOT), "work", "NOTES.md")):
    counter = os.path.join(HERE, "count_errors.py")
rc, cnt = run([sys.executable, counter])
m = re.search(r"found by the co-author:\s+(\d+)", cnt)
co_author = int(m.group(1)) if m else None
m2 = re.search(r"found by the agent:\s+(\d+)", cnt)
agent = int(m2.group(1)) if m2 else None
check("count_errors.py runs", rc == 0, cnt[-300:])
check("count_errors.py reports a co-author count", co_author is not None, cnt[-300:])
check("count_errors.py reports an agent count", agent is not None, cnt[-300:])

for lang in ("en", "ru"):
    pdf = os.path.join(PAPER, f"paper_{lang}.pdf")
    check(f"{lang}: pdf exists", os.path.exists(pdf), "missing")
    if not os.path.exists(pdf):
        continue
    n = os.path.getsize(pdf)
    check(f"{lang}: pdf non-trivial", n > 60_000, f"{n} bytes")

    rc, info = run(["pdfinfo", pdf])
    check(f"{lang}: pdfinfo ok", rc == 0, info[:200])
    m = re.search(r"Pages:\s+(\d+)", info)
    pages = int(m.group(1)) if m else 0
    check(f"{lang}: has pages", pages >= 6, f"pages={pages}")
    check(f"{lang}: A4 page size", "595" in info and "841" in info, info.splitlines()[-3:])

    rc, fonts = run(["pdffonts", pdf])
    check(f"{lang}: pdffonts ok", rc == 0, fonts[:200])
    lines = [l for l in fonts.splitlines()[2:] if l.strip()]
    check(f"{lang}: fonts are used", len(lines) > 0, "no font rows")
    not_embedded = [l.split()[0] for l in lines if len(l.split()) > 3
                    and l.split()[-4] != "yes"]
    check(f"{lang}: every font embedded", not not_embedded, f"not embedded: {not_embedded}")
    fams = {l.split()[0].split("+")[-1].split("-")[0] for l in lines}
    check(f"{lang}: main font is P052", "P052" in fams, f"families={sorted(fams)}")

    rc, txt = run(["pdftotext", "-layout", pdf, "-"])
    check(f"{lang}: pdftotext ok", rc == 0, txt[:200])

    def norm(s):
        """Same normalisation on both sides: strip punctuation, collapse space."""
        return " ".join(re.sub(r"[^0-9A-Za-z\u0400-\u04FF ]", " ", s).split()).lower()

    flat = norm(txt)

    # the glyphs the main font lacks must survive into the extracted text AS
    # THEMSELVES — a fallback that renders but loses the character would fail here
    for ch, name in (("\u207b", "superscript minus"), ("\u2080", "subscript zero"),
                     ("\u2082", "subscript two"), ("\u208a", "subscript plus"),
                     ("\u209b", "subscript s"), ("\u27f9", "Longrightarrow"),
                     ("\u2293", "sqcap")):
        check(f"{lang}: glyph {name} present in PDF text", ch in txt, f"{ch!r} missing")
    check(f"{lang}: tau present", "\u03c4" in txt, "no tau in text")
    # the combining dot must be present AND attached to a tau (that is the point)
    check(f"{lang}: combining dot present", "\u0307" in txt, "U+0307 missing")
    check(f"{lang}: combining dot follows a tau",
          "\u03c4\u0307" in txt, "no tau+U+0307 pair in the text layer")

    # authorship: both names, the co-author's count, and the agent's count
    check(f"{lang}: names Aiodam", "Aiodam" in txt)
    check(f"{lang}: names Ubaydullayev", "Ubaydullayev" in txt)
    if lang == "en":
        _w = {11: "eleven", 12: "twelve", 13: "thirteen", 14: "fourteen",
              15: "fifteen", 10: "ten", 9: "nine", 8: "eight", 7: "seven",
              6: "six", 5: "five", 4: "four", 3: "three", 2: "two", 1: "one"}
        check(f"en: states the co-author found {co_author} errors",
              re.search(rf"found {_w[co_author]}", flat) is not None,
              f"no 'found {_w[co_author]}' phrase")
        check(f"en: states the agent found {agent}",
              re.search(rf"remaining {agent}", flat) is not None,
              f"agent count {agent} not stated")
        check(f"en: total {agent + co_author} stated",
              f"{agent + co_author} errors" in flat,
              f"no '{agent + co_author} errors'")
    else:
        _w_ru = {11: "одиннадцать", 12: "двенадцать", 13: "тринадцать",
                 14: "четырнадцать", 15: "пятнадцать", 10: "десять",
                 9: "девять", 8: "восемь", 7: "семь", 6: "шесть", 5: "пять",
                 4: "четыре", 3: "три", 2: "два", 1: "один"}
        check(f"ru: states the co-author found {co_author} errors",
              _w_ru[co_author] in flat,
              f"no '{_w_ru[co_author]}' phrase")
        check(f"ru: states the agent found {agent}",
              re.search(rf"остальные {agent}", flat) is not None,
              f"agent count {agent} not stated")
        check(f"ru: total {agent + co_author} stated",
              f"{agent + co_author} ошибок" in flat,
              f"no '{agent + co_author} ошибок'")

    # figures: every referenced file exists, and its caption reaches the PDF
    md = open(os.path.join(PAPER, f"paper_{lang}.md"), encoding="utf-8").read()
    refs = re.findall(r"!\[[^\]]*\]\(([^)]+)\)", md)
    check(f"{lang}: paper references figures", len(refs) >= 3, f"refs={refs}")
    for r in refs:
        p = os.path.normpath(os.path.join(PAPER, r))
        check(f"{lang}: figure file exists {r}", os.path.exists(p), "missing")
    caps = re.findall(r"!\[([^\]]*)\]\([^)]+\)", md)
    for c in caps:
        key = norm(c)[:60]
        if key:
            check(f"{lang}: caption text in PDF ({key[:30]}…)",
                  key in flat, f"missing {key!r}")
    # a figure that failed to include would leave a literal '![' in the text
    check(f"{lang}: no unrendered image markup", "![" not in txt, "raw markdown left")

print("PDF CHECK")
print("-" * 64)
for o in ok:
    print("  ok    ", o)
for b in bad:
    print("  FAIL  ", b)
print("-" * 64)
print(f"{len(ok)} ok, {len(bad)} failed")
if bad:
    print("VERDICT: PDFS_BAD")
    sys.exit(1)
print("VERDICT: PDFS_OK")