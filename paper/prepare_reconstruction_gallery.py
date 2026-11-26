from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
from PIL import Image


ROOT = Path(__file__).resolve().parent
DESKTOP = Path("/Users/jing/Desktop")
OUTDIR = ROOT / "reconstruction_gallery"
CONTACT_SHEET = ROOT / "reconstruction_gallery_contact_sheet.png"


@dataclass(frozen=True)
class PanelSpec:
    name: str
    filename: str


SPECS = [
    PanelSpec("ffc", "ffc.png"),
    PanelSpec("fjh", "fjh.png"),
    PanelSpec("vba", "vba.png"),
]

CANVAS_SIZE = (1400, 860)
CONTENT_MARGIN = 56
MASK_THRESHOLD = 12
ARROW_MASK_WIDTH = 280
ARROW_MASK_HEIGHT = 340
BOUNDARY_BRIGHTNESS = 232
BOUNDARY_SATURATION = 26
BOUNDARY_COLOR = np.array([110, 110, 110], dtype=np.uint8)
FINAL_CORNER_MASK_WIDTH = 320
FINAL_CORNER_MASK_HEIGHT = 210


def find_structure_bbox(image: Image.Image) -> tuple[int, int, int, int]:
    arr = np.asarray(image.convert("RGBA"))
    rgb = arr[..., :3]
    alpha = arr[..., 3]
    active = (rgb.max(axis=2) > MASK_THRESHOLD) & (alpha > 0)

    height, width = active.shape
    active[height - ARROW_MASK_HEIGHT :, :ARROW_MASK_WIDTH] = False

    ys, xs = np.nonzero(active)
    if len(xs) == 0 or len(ys) == 0:
        return (0, 0, width, height)

    x0 = max(int(xs.min()) - 16, 0)
    x1 = min(int(xs.max()) + 17, width)
    y0 = max(int(ys.min()) - 16, 0)
    y1 = min(int(ys.max()) + 17, height)
    return (x0, y0, x1, y1)


def recolor_unit_cell(arr: np.ndarray) -> np.ndarray:
    rgb = arr[..., :3]
    brightness = rgb.max(axis=2)
    saturation = rgb.max(axis=2) - rgb.min(axis=2)
    mask = (brightness >= BOUNDARY_BRIGHTNESS) & (saturation <= BOUNDARY_SATURATION)
    rgb[mask] = BOUNDARY_COLOR
    return arr


def normalize_panel(image_path: Path, out_path: Path) -> None:
    original = Image.open(image_path).convert("RGBA")
    cropped = original.crop(find_structure_bbox(original))

    arr = np.asarray(cropped).copy()
    arr = recolor_unit_cell(arr)
    cropped = Image.fromarray(arr, mode="RGBA")

    canvas_w, canvas_h = CANVAS_SIZE
    inner_w = canvas_w - 2 * CONTENT_MARGIN
    inner_h = canvas_h - 2 * CONTENT_MARGIN

    scale = min(inner_w / cropped.width, inner_h / cropped.height)
    resized = cropped.resize(
        (max(1, int(round(cropped.width * scale))), max(1, int(round(cropped.height * scale)))),
        Image.Resampling.LANCZOS,
    )

    canvas = Image.new("RGBA", CANVAS_SIZE, (255, 255, 255, 255))
    paste_x = (canvas_w - resized.width) // 2
    paste_y = (canvas_h - resized.height) // 2
    canvas.alpha_composite(resized, (paste_x, paste_y))

    # Clear any residual axis labels that survived anti-aliased rendering.
    canvas_arr = np.asarray(canvas).copy()
    canvas_arr[canvas_h - FINAL_CORNER_MASK_HEIGHT :, :FINAL_CORNER_MASK_WIDTH, :3] = 255
    canvas_arr[canvas_h - FINAL_CORNER_MASK_HEIGHT :, :FINAL_CORNER_MASK_WIDTH, 3] = 255
    canvas = Image.fromarray(canvas_arr, mode="RGBA")
    canvas.save(out_path)


def make_contact_sheet(panel_paths: list[Path]) -> None:
    columns = 3
    gap = 30
    thumb_w = 620
    thumb_h = int(round(thumb_w * CANVAS_SIZE[1] / CANVAS_SIZE[0]))
    rows = (len(panel_paths) + columns - 1) // columns

    sheet_w = columns * thumb_w + (columns + 1) * gap
    sheet_h = rows * thumb_h + (rows + 1) * gap
    sheet = Image.new("RGB", (sheet_w, sheet_h), (255, 255, 255))

    for idx, panel_path in enumerate(panel_paths):
        row = idx // columns
        col = idx % columns
        x = gap + col * (thumb_w + gap)
        y = gap + row * (thumb_h + gap)
        panel = Image.open(panel_path).convert("RGB")
        panel = panel.resize((thumb_w, thumb_h), Image.Resampling.LANCZOS)
        sheet.paste(panel, (x, y))

    sheet.save(CONTACT_SHEET)


def main() -> None:
    OUTDIR.mkdir(exist_ok=True)
    outputs: list[Path] = []

    for spec in SPECS:
        input_path = DESKTOP / spec.filename
        output_path = OUTDIR / f"{spec.name}_panel.png"
        normalize_panel(input_path, output_path)
        outputs.append(output_path)
        print(f"wrote {output_path}")

    make_contact_sheet(outputs)
    print(f"wrote {CONTACT_SHEET}")


if __name__ == "__main__":
    main()
