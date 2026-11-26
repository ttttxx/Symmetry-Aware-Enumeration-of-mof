"""Superseded raster-only text edits for the MATCH pipeline figure.

The current source of truth is `paper/流水线.drawio`, edited by
`codex/edit_pipeline_drawio.py` and exported with draw.io. This raster script
is kept only to document the earlier PNG-only pass.
"""

from pathlib import Path
from PIL import Image, ImageDraw, ImageFont


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "paper" / "pipline.png"
OUT = ROOT / "paper" / "pipline.png"
PREVIEW = ROOT / "codex" / "validation" / "match_pipeline_textedit_20261020.png"

ARIAL = "/System/Library/Fonts/Supplemental/Arial.ttf"
ARIAL_BOLD = "/System/Library/Fonts/Supplemental/Arial Bold.ttf"
MENLO = "/System/Library/Fonts/Menlo.ttc"


def erase(draw, box, fill):
    draw.rectangle(box, fill=fill)


def main():
    img = Image.open(SRC).convert("RGB")
    draw = ImageDraw.Draw(img)

    # Panel (d) header: remove unsupported "seconds" and "AIMD" claims.
    erase(draw, (1265, 5, 1645, 50), (255, 255, 255))
    header_font = ImageFont.truetype(ARIAL, 17)
    draw.text((1452, 20), "Zero-redundancy library of pre-relaxation",
              font=header_font, fill=(85, 85, 85), anchor="mm")
    draw.text((1452, 40), "models for downstream evaluation.",
              font=header_font, fill=(85, 85, 85), anchor="mm")

    # Panel (c) code box: 0.5 A -> 0.2 A.
    erase(draw, (976, 466, 1008, 484), (245, 245, 245))
    code_font = ImageFont.truetype(MENLO, 12)
    draw.text((978, 467), "0.2)", font=code_font, fill=(34, 34, 34))

    # Panel (a): the paper's Gamma is the proper-operation permutation image,
    # not the full graph automorphism group.
    erase(draw, (20, 462, 255, 512), (255, 255, 255))
    a_font = ImageFont.truetype(ARIAL, 12)
    draw.text((136, 472), "Proper-Operation Permutation",
              font=a_font, fill=(31, 42, 55), anchor="mm")
    draw.text((136, 492), "Group Γ = Φ(S_emb^+)",
              font=a_font, fill=(31, 42, 55), anchor="mm")

    # Panel (b): replace the unexplained correlation plot with a text block
    # describing the task output; the search steps are already in the text above.
    erase(draw, (270, 252, 530, 270), (255, 255, 255))
    erase(draw, (270, 246, 510, 263), (255, 255, 255))
    erase(draw, (285, 270, 710, 528), (255, 255, 255))
    draw.rounded_rectangle(
        [285, 262, 700, 430], radius=8, fill=(247, 249, 251),
        outline=(201, 211, 219), width=1
    )
    draw.text((300, 278), "task output", font=ImageFont.truetype(ARIAL_BOLD, 12),
              fill=(31, 42, 55), anchor="lm")
    out_font = ImageFont.truetype(ARIAL, 11.5)
    out_lines = [
        "canonical leaf κ(c)",
        "Ω: one representative per orbit",
        "merge identical P_H; global dedup",
        "optional per-stratum quota B_s",
    ]
    for i, line in enumerate(out_lines):
        draw.text((302, 302 + 30 * i), line, font=out_font,
                  fill=(31, 42, 55), anchor="lm")

    label_font = ImageFont.truetype(ARIAL_BOLD, 14)
    label_color = (26, 26, 26)
    label_bg = (250, 250, 250)

    # Top-left: isotropic limit -> exact rho = 1.
    erase(draw, (1256, 56, 1408, 73), label_bg)
    draw.text((1332, 64.5), "Exact rho = 1",
              font=label_font, fill=label_color, anchor="mm")

    # Top-right: tier 2 -> exact-rho range.
    erase(draw, (1504, 56, 1655, 73), label_bg)
    draw.text((1579.5, 64.5), "Exact rho in (1/|Γ|, 1)",
              font=label_font, fill=label_color, anchor="mm")

    # Bottom-right: tier 3 -> trivial-stabilizer exact rho.
    erase(draw, (1499, 295, 1672, 313), label_bg)
    draw.text((1585.5, 304.5), "Exact rho = 1/|Γ|",
              font=label_font, fill=label_color, anchor="mm")

    # Bottom-right note: no high-entropy filtering; representatives are retained.
    erase(draw, (1503, 480, 1655, 510), label_bg)
    note_font = ImageFont.truetype(ARIAL, 12)
    draw.text((1579, 489), "Trivial stabilizer.",
              font=note_font, fill=(51, 51, 51), anchor="mm")
    draw.text((1579, 503), "Representatives retained.",
              font=note_font, fill=(51, 51, 51), anchor="mm")

    img.save(OUT, dpi=(300, 300))
    img.save(PREVIEW, dpi=(300, 300))
    print(OUT)
    print(PREVIEW)


if __name__ == "__main__":
    main()
