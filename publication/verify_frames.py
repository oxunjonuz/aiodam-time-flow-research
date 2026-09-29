#!/usr/bin/env python3
"""Check that no INK leaves the text block, on every page of both PDFs.

The owner found a formula running past the right margin. This answers that
question the way a reader sees it: it rasterises each page and looks for dark
pixels outside the text block.

Why ink and not word boxes. `pdftotext -bbox` reports a glyph's ADVANCE box, not
its ink: for a math glyph like the minus in `10⁻¹⁴` it reports x0 = 64.8 pt while
the ink actually starts at ~70 pt, inside the block. A check built on those boxes
therefore reports a false overflow. Measured, not assumed — the crop of that line
shows the mark inside the margin. Ink is the ground truth here, so ink is what is
measured.

Geometry is taken from the preamble (a4paper, margin 2.4 cm, top 2.6 cm,
bottom 2.4 cm). The header/footer sit in the margin by design, so the band above
the top margin and below the bottom margin is excluded from the horizontal test;
only the left and right edges are tested there, and the top/bottom tests use the
physical page edge.

Exit 0 => FRAME_OK. Exit 1 => ink leaves the frame.
"""
import os
import subprocess
import sys
import tempfile

from PIL import Image

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(sys.argv[1]) if len(sys.argv) > 1 else HERE
PAPER = os.path.join(ROOT, "paper")

PT_PER_CM = 72.0 / 2.54
MARGIN = 2.4 * PT_PER_CM
TOP = 2.6 * PT_PER_CM
BOTTOM = 2.4 * PT_PER_CM
DPI = 200
# a glyph's antialiased edge can reach a fraction of a pixel past the box; allow
# one pixel (0.36 pt at 200 dpi) and no more
TOL_PX = 1


def dark_mask(img):
    """True where the page has ink (luma below a threshold)."""
    g = img.convert("L")
    return g.point(lambda v: 255 if v < 200 else 0)


def main():
    ok, bad = [], []
    for lang in ("en", "ru"):
        pdf = os.path.join(PAPER, f"paper_{lang}.pdf")
        if not os.path.exists(pdf):
            bad.append(f"{lang}: missing {pdf}")
            continue
        with tempfile.TemporaryDirectory() as tmp:
            p = subprocess.run(["pdftoppm", "-png", "-r", str(DPI), pdf,
                                os.path.join(tmp, "pg")],
                               capture_output=True, text=True)
            if p.returncode != 0:
                bad.append(f"{lang}: pdftoppm failed: {p.stderr[:200]}")
                continue
            pages = sorted(f for f in os.listdir(tmp) if f.endswith(".png"))
            ok.append(f"{lang}: {len(pages)} pages rasterised at {DPI} dpi")
            scale = DPI / 72.0
            worst = 0
            offenders = []
            for name in pages:
                idx = int(name.split("-")[-1].split(".")[0])
                im = Image.open(os.path.join(tmp, name))
                W, H = im.size
                m = dark_mask(im)
                px = m.load()
                left = int(MARGIN * scale) - TOL_PX
                right = int((W / scale - MARGIN) * scale) + TOL_PX
                top = int(TOP * scale) - TOL_PX
                bottom = int((H / scale - BOTTOM) * scale) + TOL_PX
                # scan the whole page for ink; a hit outside the block is reported
                for y in range(0, H, 1):
                    row = None
                    for x in (list(range(0, max(left, 1))) +
                              list(range(min(right + 1, W), W))):
                        if px[x, y]:
                            row = x
                            break
                    if row is not None:
                        over = max(left - row, row - right)
                        worst = max(worst, over)
                        offenders.append((idx, "horizontal", row, y, over))
                # vertical: ink above the top margin or below the bottom margin is
                # fine for the running head/foot, but ink past the PHYSICAL edge is
                # not. The page edge is checked by the rasteriser itself (the image
                # is exactly the page), so only a full-bleed check is meaningful
                # here: report ink touching the outermost 2 px.
                for y in list(range(0, 2)) + list(range(H - 2, H)):
                    for x in range(W):
                        if px[x, y]:
                            offenders.append((idx, "page-edge", x, y, 99))
                            break
            if offenders:
                bad.append(f"{lang}: {len(offenders)} ink pixel(s) outside the block")
                seen = set()
                for idx, kind, x, y, over in offenders:
                    if (idx, kind) in seen:
                        continue
                    seen.add((idx, kind))
                    bad.append(f"    page {idx}: {kind} at x={x} y={y} "
                               f"(over by {over}px)")
                    if len(seen) > 6:
                        break
            else:
                ok.append(f"{lang}: no ink outside the text block "
                          f"(worst {worst}px ≤ {TOL_PX}px tolerance)")

    print("FRAME CHECK (ink, not advance boxes)")
    print("-" * 64)
    for o in ok:
        print("  ok    ", o)
    for b in bad:
        print("  FAIL  ", b)
    print("-" * 64)
    print(f"{len(ok)} ok, {len(bad)} failed")
    if bad:
        print("VERDICT: FRAME_BAD")
        return 1
    print("VERDICT: FRAME_OK")
    return 0


if __name__ == "__main__":
    sys.exit(main())