#!/usr/bin/env python3
"""Generate MANIFEST.md: sha256 of every published file, from the bytes on disk.

Run from the publication folder. Writes MANIFEST.md in place.

ORDER MATTERS: MANIFEST.md is itself checked against the bytes by
verify_publication.py, and the manifest covers verify_publication.py. So the
manifest must be regenerated AFTER any edit to a published file, and
verify_publication.py must be run after that. The sequence is:

    make_figures.py -> build_paper.sh -> make_manifest.py -> verify_publication.py

(the figures and the PDFs are published files too, so they come before the
manifest, not after it).
"""
import hashlib
import os

HERE = os.path.dirname(os.path.abspath(__file__))

# files that are part of the publication package (this script writes the list into
# MANIFEST.md, so it must not include MANIFEST.md itself)
INCLUDE_DIRS = ["", "paper", "docs_en", "figures"]
SKIP = {"MANIFEST.md", "make_manifest.py"}


def sha(path):
    with open(path, "rb") as f:
        return hashlib.sha256(f.read()).hexdigest(), os.path.getsize(path)


rows = []
for d in INCLUDE_DIRS:
    folder = os.path.join(HERE, d) if d else HERE
    for name in sorted(os.listdir(folder)):
        p = os.path.join(folder, name)
        if not os.path.isfile(p):
            continue
        rel = os.path.relpath(p, HERE).replace(os.sep, "/")
        if rel in SKIP:
            continue
        h, n = sha(p)
        rows.append((rel, n, h))

total = sum(n for _, n, _ in rows)

lines = []
lines.append("# MANIFEST — published files and their SHA-256")
lines.append("")
lines.append("Generated from the bytes on disk by `make_manifest.py`. "
             f"{len(rows)} files, {total} bytes total.")
lines.append("")
lines.append("*Сгенерировано из байтов на диске скриптом `make_manifest.py`.*")
lines.append("")
lines.append("| file | bytes | SHA-256 |")
lines.append("|---|---:|---|")
for rel, n, h in rows:
    lines.append(f"| `{rel}` | {n} | `{h}` |")
lines.append("")
lines.append("## What is NOT in this manifest")
lines.append("")
lines.append("The research tree itself (`work/`, `data/`, `sources/`) lives one level "
             "up and is not duplicated here. Its measurement artifacts are indexed, "
             "with their own hashes, in `ARTIFACTS.md`; the Lean formalization "
             "sources are indexed there too. The two PDFs in `paper/` are built from "
             "the two Markdown files beside them by `pandoc` + `xelatex` "
             "(`build_paper.sh` + `preamble.tex`), and the figures in `figures/` are "
             "drawn from the real GW150914 strain and the frozen artifacts by "
             "`make_figures.py` — so the PDFs and the figures are derived, not "
             "independent, artifacts.")
lines.append("")
lines.append("The PNG files beside each figure are low-resolution previews written for "
             "the visual check; the PDFs are what the paper includes.")
lines.append("")

with open(os.path.join(HERE, "MANIFEST.md"), "w", encoding="utf-8") as f:
    f.write("\n".join(lines))

print(f"wrote MANIFEST.md: {len(rows)} files, {total} bytes")
for rel, n, h in rows:
    print(f"  {h[:12]}…  {n:>8}  {rel}")
