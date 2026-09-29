#!/bin/sh
# Build the two paper PDFs from the two Markdown sources.
#
#   sh build_paper.sh            build both
#   sh build_paper.sh en         build only paper_en
#
# The Markdown sources are NOT modified: the build works on a copy in a temp
# directory, and the only substitution applied to that copy is
#     τ̇  ->  $\dot{\tau}$
# because U+0307 (combining dot above) is absent from the main font and cannot be
# redefined generically by newunicodechar. Everything else is handled by the
# preamble's explicit glyph map, so the PDF text stays identical to the source.
#
# Requires: pandoc, xelatex, P052 + TeX Gyre Heros + DejaVu Sans Mono + Latin Modern Math.
set -e

HERE=$(cd "$(dirname "$0")" && pwd)
PAPER="$HERE/paper"
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT

build() {
  lang="$1"
  src="$PAPER/paper_${lang}.md"
  out="$PAPER/paper_${lang}.pdf"
  [ -f "$src" ] || { echo "missing $src" >&2; exit 1; }
  rm -f "$out"   # so a stale PDF can never be mistaken for a fresh build

  # the only textual substitution, applied to the COPY
  python3 - "$src" "$TMP/paper_${lang}.md" <<'PY'
import sys
src, dst = sys.argv[1], sys.argv[2]
t = open(src, encoding="utf-8").read()
n = t.count("\u03c4\u0307")
t = t.replace("\u03c4\u0307", r"\tauDot{}")
open(dst, "w", encoding="utf-8").write(t)
print(f"  {src}: {n} occurrence(s) of tau+combining-dot mapped to \\tauDot")
PY

  echo "  pandoc -> xelatex  ($lang)"
  pandoc "$TMP/paper_${lang}.md" \
    -o "$out" \
    --pdf-engine=xelatex \
    --pdf-engine-opt=-interaction=nonstopmode \
    --include-in-header="$HERE/preamble.tex" \
    --toc --toc-depth=2 \
    -V documentclass=article \
    -V fontsize=11pt \
    -V microtypeoptions=protrusion=false \
    -f markdown+pipe_tables+raw_tex \
    --resource-path="$PAPER:$HERE" \
    2>&1 | grep -v -E "^(Overfull|Underfull|\[[0-9]+\]|LaTeX Warning: Reference)" || true

  [ -s "$out" ] || { echo "FAILED: $out not produced" >&2; exit 1; }
  echo "  wrote $out ($(wc -c < "$out") bytes)"
}

case "${1:-both}" in
  en)   build en ;;
  ru)   build ru ;;
  both) build en; build ru ;;
  *)    echo "usage: sh build_paper.sh [en|ru|both]" >&2; exit 2 ;;
esac

echo "done."