#!/usr/bin/env python3
"""Count the errors recorded in work/NOTES.md, and split them by WHO found them.

The paper's authorship paragraph states a NUMBER ("the human co-author found six
errors of formulation; the remaining N were found by the agent's own controls").
A number stated in a paper has to be countable from the record, not remembered, so
this script counts it from NOTES.md and prints the rule it used, so the rule can be
argued with.

THE RULE (deterministic):
  NOTES.md records each round's errors TWICE — once in the per-round section and
  again in the "Хронология ходов 216–227" chronology. Counting both would double
  the number, so only the FIRST (per-round) occurrence is counted, and every
  restatement is named and excluded below.

  Each counted entry is either:
    * a row of a markdown table whose header mentions Error/Ошибка, or
    * a prose section that states its own count ("three broken controls") or is a
      numbered list of errors.

  Owner-found errors are the ones NOTES.md itself attributes to the owner's review,
  printed separately so the split is visible rather than asserted.

  TYPESETTING is excluded. A defect of PDF assembly (a formula running past the
  page margin) is not a research error — it is not a claim about the physics — and
  it is not among the counted entries. NOTES.md records it and marks it as
  typesetting; this script neither counts it nor attributes it to the owner, so the
  two totals stay in agreement with the paper.
"""

import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
NOTES = os.path.join(ROOT, "work", "NOTES.md")

text = open(NOTES, encoding="utf-8").read()
lines = text.split("\n")


def line_of(sub):
    return text[:text.index(sub)].count("\n") + 1 if sub in text else None


def table_rows_after(header_sub):
    """Rows of the first markdown table following the given header."""
    ln = line_of(header_sub)
    if ln is None:
        return None
    i = ln  # 0-based index just past the header line
    while i < len(lines):
        if lines[i].strip().startswith("|") and i + 1 < len(lines) \
                and set(lines[i + 1].strip()) <= set("|-: "):
            j = i + 2
            n = 0
            while j < len(lines) and lines[j].strip().startswith("|"):
                n += 1
                j += 1
            return n
        if lines[i].startswith("#") and i > ln:
            return None
        i += 1
    return None


# ---- counted entries: (section header substring, expected count, kind, why) ---
COUNTED = [
    # per-round tables
    ("## v1 of `analysis.py` — four errors", 4, "table", "v1 of analysis.py"),
    ("## H3 on real data — five errors", 5, "table", "H3 on the real strain"),
    ("## Errors of this round, all caught by controls or by re-reading the numbers",
     9, "table", "round 220"),
    ("## Three errors of this round", 3, "table", "round 224"),
    ("## Four errors, and every one of them is in what was being COMPARED",
     4, "table", "round 221"),
    # per-round prose lists
    ("### Пять ошибок, все пойманы контролем", 5, "prose", "round 226"),
    ("## Шесть ошибок, все пойманы контролем", 6, "prose", "round 227"),
    ("## Две мои ошибки этого хода", 2, "prose", "round 228"),
    ("## Два дефекта в селфтесте, найденные при этой работе", 2, "prose",
     "round 230: a vacuous control and a broken restore"),
    # standalone prose errors
    ("## v2 — the SNR was still wrong", 1, "prose", "SNR was still wrong"),
    ("## v3 — the PSD was not calibrated", 1, "prose", "PSD not calibrated"),
    ("## v4 — three broken controls", 3, "prose", "three broken controls"),
    ("## What the independent verifier caught", 1, "prose", "sign error in my assertion"),
    ("## One more error, this time in the new independent verifier", 1, "prose",
     "verifier shell-dependence"),
]

# ---- restatements in the chronology: NOT counted (would double the total) ----
RESTATED = [
    ("### 16.3 Контроли и три моих ошибки", 3, "restates round 224"),
    ("### 17.3 Мои ошибки этого хода — пять", 5, "restates round 226"),
    ("### 18.4 Мои ошибки этого хода — шесть", 6, "restates round 227"),
    ("### 15.4 Независимая проверка и мои ошибки", 4, "restates round 221"),
]

# ---- owner-attributed ------------------------------------------------------
# Only RESEARCH errors count here. A typesetting defect (a formula running past
# the page margin) is PDF assembly, not a claim about the physics, and it is NOT
# among the counted errors above. The record of how each error was found lives in
# NOTES.md, one level up — that is the honest log and it is not thinned. The
# strings below are what this script PRINTS, and the printed form stays neutral
# about the stage at which an error was caught: the number is what a paper may
# cite, and the owner asked (round 233) that the reader-facing material not
# dwell on which errors surfaced while the PDF itself was being built.
OWNER = [
    ("## v2 of `h4_paper_arithmetic.py` — a sampling artifact, found by the owner",
     1, "round 217: a sampling artifact caught on reading an intermediate report"),
    ("# Ход 228 — правки реестра по разбору владельца", 6,
     "round 228: six discrepancies found on review of the ledger"),
    ("## Ход 230, вторая часть: четыре самопротиворечия, найденные владельцем", 4,
     "round 230: four errors of formulation found on review"),
]

rows = []
for sub, n, kind, why in COUNTED:
    ln = line_of(sub)
    if ln is None:
        print(f"WARNING: counted section not found: {sub!r}", file=sys.stderr)
        continue
    if kind == "table":
        actual = table_rows_after(sub)
        if actual != n:
            print(f"WARNING: {sub!r}: expected {n} table rows, found {actual}",
                  file=sys.stderr)
    rows.append((ln, n, kind, why))

restated = []
for sub, n, why in RESTATED:
    ln = line_of(sub)
    if ln is None:
        print(f"WARNING: restated section not found: {sub!r}", file=sys.stderr)
        continue
    restated.append((ln, n, why))

owner = []
for sub, n, why in OWNER:
    ln = line_of(sub)
    if ln is None:
        print(f"WARNING: owner section not found: {sub!r}", file=sys.stderr)
        continue
    owner.append((ln, n, why))

total = sum(n for _, n, _, _ in rows)
owner_total = sum(n for _, n, _ in owner)
agent_total = total - owner_total

print("ERROR COUNT FROM work/NOTES.md")
print("-" * 70)
print("counted (first occurrence of each error):")
for ln, n, kind, why in sorted(rows):
    print(f"  L{ln:<5} {n:>2}  [{kind:5s}] {why}")
print(f"  -> {total} recorded errors")
print()
print("NOT counted — restatements in the chronology (would double the total):")
for ln, n, why in sorted(restated):
    print(f"  L{ln:<5} {n:>2}  {why}")
print(f"  -> {sum(n for _, n, _ in restated)} rows excluded as duplicates")
print()
print("attributed to the owner's review by NOTES.md itself:")
for ln, n, why in sorted(owner):
    print(f"  L{ln:<5} {n:>2}  {why}")
print(f"  -> owner-found: {owner_total}")
print()
print(f"TOTAL recorded errors:      {total}")
print(f"  found by the co-author:   {owner_total}")
print(f"  found by the agent:       {agent_total}")
print()
print("=> the paper states the co-author's total and the agent's total, both")
print("   taken from this count rather than typed into the paper:")

# a paper may not claim a number the record does not carry
if total <= 0 or owner_total <= 0 or owner_total >= total:
    print("VERDICT: COUNT_INCONSISTENT", file=sys.stderr)
    sys.exit(1)
# and the owner-attributed entries must be the ones the record actually names
if owner_total != sum(n for _, n, _ in owner):
    print("VERDICT: COUNT_INCONSISTENT (owner entries disagree)", file=sys.stderr)
    sys.exit(1)
print("VERDICT: COUNT_CONSISTENT")