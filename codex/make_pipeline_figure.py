"""Superseded 2x2 redraw of the MATCH pipeline figure.

The author requested minimal text edits on the original four-panel layout.
The current source of `paper/pipline.png` is `codex/edit_pipeline_text.py`;
this file is retained only as the superseded redraw alternative.
"""

from pathlib import Path
from PIL import Image, ImageDraw, ImageFont


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "paper" / "pipline.png"
PREVIEW = ROOT / "codex" / "validation" / "match_pipeline_20261020.png"

W, H = 1800, 1300
MARGIN = 30
GAP = 22
PANEL_W = (W - 2 * MARGIN - GAP) // 2
PANEL_H = (H - 2 * MARGIN - GAP) // 2

FONT_PATH = "/System/Library/Fonts/Supplemental/Arial Unicode.ttf"
FONT_BOLD = "/System/Library/Fonts/Supplemental/Arial Bold.ttf"

INK = (31, 42, 55)
MUTED = (82, 96, 109)
LINE = (176, 189, 199)
BLUE = (42, 92, 138)
TEAL = (33, 121, 132)
ORANGE = (190, 112, 42)
GREEN = (67, 122, 78)
PANEL_FILL = (247, 249, 251)
BOX_FILL = (255, 255, 255)
BLUE_FILL = (232, 241, 248)
TEAL_FILL = (229, 243, 244)
ORANGE_FILL = (251, 240, 226)
GREEN_FILL = (235, 244, 234)


def font(size, bold=False):
    return ImageFont.truetype(FONT_BOLD if bold else FONT_PATH, size)


def wrap(draw, text, fnt, max_width):
    words = text.split()
    lines = []
    current = ""
    for word in words:
        candidate = word if not current else current + " " + word
        if draw.textlength(candidate, font=fnt) <= max_width:
            current = candidate
        else:
            if current:
                lines.append(current)
            current = word
    if current:
        lines.append(current)
    return lines


def text_block(draw, box, text, size=27, bold=False, color=INK, align="center",
               line_gap=5):
    x, y, w, h = box
    fnt = font(size, bold=bold)
    lines = []
    for paragraph in text.split("\n"):
        lines.extend(wrap(draw, paragraph, fnt, w - 24))
    line_h = size + line_gap
    total_h = len(lines) * line_h - line_gap
    cy = y + (h - total_h) / 2
    for line in lines:
        if align == "center":
            tx = x + w / 2
            anchor = "mm"
        else:
            tx = x + 12
            anchor = "lm"
        draw.text((tx, cy + line_h / 2), line, font=fnt, fill=color, anchor=anchor)
        cy += line_h


def box(draw, box_xywh, text, fill=BOX_FILL, outline=LINE, size=27, bold=False,
        color=INK, width=2):
    x, y, w, h = box_xywh
    draw.rounded_rectangle(
        [x, y, x + w, y + h], radius=11, fill=fill, outline=outline, width=width
    )
    text_block(draw, box_xywh, text, size=size, bold=bold, color=color)


def arrow(draw, x1, y1, x2, y2, color=MUTED, width=4, head=10):
    draw.line([(x1, y1), (x2, y2)], fill=color, width=width)
    if abs(y2 - y1) >= abs(x2 - x1):
        direction = 1 if y2 >= y1 else -1
        pts = [(x2, y2), (x2 - head, y2 - direction * head * 1.6),
               (x2 + head, y2 - direction * head * 1.6)]
    else:
        direction = 1 if x2 >= x1 else -1
        pts = [(x2, y2), (x2 - direction * head * 1.6, y2 - head),
               (x2 - direction * head * 1.6, y2 + head)]
    draw.polygon(pts, fill=color)


def panel(draw, x, y, w, h, tag, title):
    draw.rounded_rectangle(
        [x, y, x + w, y + h], radius=16, fill=PANEL_FILL, outline=(210, 218, 224),
        width=2
    )
    draw.text((x + 22, y + 16), tag, font=font(34, bold=True), fill=TEAL)
    draw.text((x + 22, y + 56), title, font=font(31, bold=True), fill=INK)


def footer(draw, x, y, w, text, color=MUTED, size=22):
    draw.rounded_rectangle(
        [x, y, x + w, y + 64], radius=8, fill=(255, 255, 255),
        outline=(220, 226, 231), width=1
    )
    text_block(draw, (x, y, w, 64), text, size=size, color=color)


def main():
    img = Image.new("RGB", (W, H), "white")
    draw = ImageDraw.Draw(img)

    xs = [MARGIN, MARGIN + PANEL_W + GAP]
    ys = [MARGIN, MARGIN + PANEL_H + GAP]
    pad = 24
    iw = PANEL_W - 2 * pad

    # (a) Input and assignment space
    x, y0 = xs[0], ys[0]
    panel(draw, x, y0, PANEL_W, PANEL_H, "(a)", "Input and assignment space")
    yy = y0 + 110
    box(draw, (x + pad, yy, iw, 64), "Embedded net + repeat cell", BLUE_FILL, BLUE, size=30)
    arrow(draw, x + PANEL_W / 2, yy + 64, x + PANEL_W / 2, yy + 82)
    yy += 82
    box(draw, (x + pad, yy, iw, 64), "Quotient graph Q = Q/T", BLUE_FILL, BLUE, size=30)
    arrow(draw, x + PANEL_W / 2, yy + 64, x + PANEL_W / 2, yy + 82)
    yy += 82
    box(draw, (x + pad, yy, iw, 64), "Degree-typed node pools C_d", TEAL_FILL, TEAL, size=30)
    arrow(draw, x + PANEL_W / 2, yy + 64, x + PANEL_W / 2, yy + 82)
    yy += 82
    box(draw, (x + pad, yy, iw, 78), "Node assignments X: one node label per vertex",
        TEAL_FILL, TEAL, size=30)
    yy += 96
    box(draw, (x + pad, yy, iw, 78),
        "Fixed edge template 1B_4H_Ch\n(reconstruction input; not in X)",
        ORANGE_FILL, ORANGE, size=26)
    footer(draw, x + pad, y0 + PANEL_H - 76, iw,
           "Node labels are the enumerated variable", color=TEAL, size=22)

    # (b) Symmetry and candidate identity
    x, y0 = xs[1], ys[0]
    panel(draw, x, y0, PANEL_W, PANEL_H, "(b)", "Symmetry and candidate identity")
    yy = y0 + 110
    box(draw, (x + pad, yy, iw, 78), "Proper crystallographic operations S_emb^+",
        BLUE_FILL, BLUE, size=29)
    yy += 98
    arrow(draw, x + PANEL_W / 2, yy, x + PANEL_W / 2, yy + 20)
    yy += 20
    box(draw, (x + pad, yy, iw, 64), "Γ = Φ(S_emb^+) ≤ Sym(V)", TEAL_FILL, TEAL, size=30)
    yy += 84
    arrow(draw, x + PANEL_W / 2, yy, x + PANEL_W / 2, yy + 20)
    yy += 20
    box(draw, (x + pad, yy, iw, 70), "Orbits Γ·c; canonical representative κ(c)",
        TEAL_FILL, TEAL, size=29)
    yy += 90
    arrow(draw, x + PANEL_W / 2, yy, x + PANEL_W / 2, yy + 20)
    yy += 20
    box(draw, (x + pad, yy, iw, 88),
        "N_Γ = |X/Γ|\nρ(c) = |Stab(c)| / |Γ|", GREEN_FILL, GREEN, size=30)
    footer(draw, x + pad, y0 + PANEL_H - 76, iw,
           "Reflections remain distinct unless a proper operation maps them",
           color=ORANGE, size=22)

    # (c) Subgroup-lattice enumeration
    x, y0 = xs[0], ys[1]
    panel(draw, x, y0, PANEL_W, PANEL_H, "(c)", "Subgroup-lattice enumeration")
    # lattice icon on the left
    cx = x + 190
    top = y0 + 126
    draw.ellipse([cx - 30, top - 4, cx + 30, top + 56], fill=TEAL_FILL, outline=TEAL, width=2)
    draw.text((cx, top + 26), "Γ", font=font(27, bold=True), fill=TEAL, anchor="mm")
    hy = top + 100
    for dx, lab in [(-82, "H1"), (82, "H2")]:
        draw.line([(cx, top + 56), (cx + dx, hy)], fill=LINE, width=3)
        draw.ellipse([cx + dx - 30, hy - 4, cx + dx + 30, hy + 56],
                     fill=BLUE_FILL, outline=BLUE, width=2)
        draw.text((cx + dx, hy + 26), lab, font=font(25, bold=True), fill=BLUE, anchor="mm")
    ey = hy + 92
    for dx in (-82, 82):
        draw.line([(cx + dx, hy + 56), (cx, ey)], fill=LINE, width=3)
    draw.ellipse([cx - 32, ey - 4, cx + 32, ey + 58], fill=ORANGE_FILL, outline=ORANGE, width=2)
    draw.text((cx, ey + 27), "{e}", font=font(25, bold=True), fill=ORANGE, anchor="mm")
    draw.text((cx, ey + 80), "subgroup lattice", font=font(21), fill=MUTED, anchor="mm")

    # enumeration steps on the right
    rx = x + 370
    rw = PANEL_W - 370 - pad
    yy = y0 + 120
    box(draw, (rx, yy, rw, 70), "Fixed-point domains Fix(H)", TEAL_FILL, TEAL, size=27)
    yy += 92
    box(draw, (rx, yy, rw, 70), "DFS + prefix rejection", BLUE_FILL, BLUE, size=27)
    yy += 92
    box(draw, (rx, yy, rw, 86), "Optional exact-rho quotas B_s for Ω_s",
        ORANGE_FILL, ORANGE, size=25)
    yy += 110
    draw.text((rx, yy), "Exact-rho strata (vba):", font=font(23, bold=True), fill=INK)
    draw.text((rx, yy + 30), "1, 1/2, 1/4, 1/8, 1/16", font=font(23), fill=INK)
    draw.text((rx, yy + 62), "retained across strata, not filtered",
              font=font(21), fill=MUTED)

    # (d) Oriented reconstruction
    x, y0 = xs[1], ys[1]
    panel(draw, x, y0, PANEL_W, PANEL_H, "(d)", "Oriented reconstruction")
    yy = y0 + 110
    box(draw, (x + pad, yy, iw, 64), "Selected c* → node templates", TEAL_FILL, TEAL, size=30)
    yy += 84
    arrow(draw, x + PANEL_W / 2, yy, x + PANEL_W / 2, yy + 20)
    yy += 20
    box(draw, (x + pad, yy, iw, 64), "Kabsch proper rotation; det R = 1", BLUE_FILL, BLUE, size=29)
    yy += 84
    arrow(draw, x + PANEL_W / 2, yy, x + PANEL_W / 2, yy + 20)
    yy += 20
    box(draw, (x + pad, yy, iw, 80),
        "Linker placement; bridge repair; X fusion < 0.8 Å", ORANGE_FILL, ORANGE, size=29)
    yy += 100
    arrow(draw, x + PANEL_W / 2, yy, x + PANEL_W / 2, yy + 20)
    yy += 20
    box(draw, (x + pad, yy, iw, 80),
        "Remove X; merge sites < 0.2 Å; export CIF", GREEN_FILL, GREEN, size=27)
    footer(draw, x + pad, y0 + PANEL_H - 76, iw,
           "Pre-relaxation periodic models for downstream screening", color=GREEN, size=24)

    img.save(OUT, dpi=(300, 300))
    img.save(PREVIEW, dpi=(300, 300))
    print(OUT)
    print(PREVIEW)


if __name__ == "__main__":
    main()
