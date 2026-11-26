"""Superseded drawio edit generated before the author's own final drawing.

The current source of truth is the author-edited `paper/流水线.drawio` and
`paper/流水线.png` (copied to `paper/pipline.png`). Do not rerun this script:
it starts from the backup and would overwrite the author's final layout/text.
"""

from pathlib import Path
import xml.etree.ElementTree as ET


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "codex" / "validation" / "流水线_original_20261020.drawio"
OUT = ROOT / "paper" / "流水线.drawio"


def geom(cell):
    return cell.find("mxGeometry")


def set_xy(cell, x=None, y=None):
    g = geom(cell)
    if g is None:
        return
    if x is not None:
        g.set("x", f"{x:g}")
    if y is not None:
        g.set("y", f"{y:g}")


def set_size(cell, w=None, h=None):
    g = geom(cell)
    if g is None:
        return
    if w is not None:
        g.set("width", f"{w:g}")
    if h is not None:
        g.set("height", f"{h:g}")


def add_text(root, cell_id, value, x, y, w, h, font_size, bold=False):
    style = (
        "text;html=1;strokeColor=none;fillColor=none;align=left;"
        "verticalAlign=top;whiteSpace=wrap;rounded=0;"
        f"fontSize={font_size};"
    )
    cell = ET.SubElement(root, "mxCell", {
        "id": cell_id,
        "parent": "1",
        "style": style,
        "value": value,
        "vertex": "1",
    })
    g = ET.SubElement(cell, "mxGeometry", {
        "height": f"{h:g}",
        "width": f"{w:g}",
        "x": f"{x:g}",
        "y": f"{y:g}",
        "as": "geometry",
    })
    return cell


def add_box(root, cell_id, value, x, y, w, h, fill, stroke, font_size=10.5):
    style = (
        "rounded=1;whiteSpace=wrap;html=1;"
        f"fillColor={fill};strokeColor={stroke};fontSize={font_size};"
        "align=center;verticalAlign=middle;arcSize=8;"
    )
    cell = ET.SubElement(root, "mxCell", {
        "id": cell_id,
        "parent": "1",
        "style": style,
        "value": value,
        "vertex": "1",
    })
    ET.SubElement(cell, "mxGeometry", {
        "height": f"{h:g}",
        "width": f"{w:g}",
        "x": f"{x:g}",
        "y": f"{y:g}",
        "as": "geometry",
    })
    return cell


def main():
    tree = ET.parse(SRC)
    root = tree.getroot()
    model = root.find("diagram").find("mxGraphModel")
    cells = model.find("root")

    by_id = {c.get("id"): c for c in cells.findall("mxCell")}

    # Narrow panel (b) and its label.
    for cell_id in ("RH0wHt1T8WRH78aJM92P-415", "RH0wHt1T8WRH78aJM92P-416"):
        set_size(by_id[cell_id], w=390)

    # Move the lattice nodes left and down so the text block can sit below.
    lattice = {
        "RH0wHt1T8WRH78aJM92P-420": (-1685, 320),  # Gamma
        "RH0wHt1T8WRH78aJM92P-421": (-1760, 420),  # H1
        "RH0wHt1T8WRH78aJM92P-422": (-1680, 420),  # Hk
        "RH0wHt1T8WRH78aJM92P-423": (-1600, 420),  # Hn
        "RH0wHt1T8WRH78aJM92P-424": (-1685, 510),  # {e}
    }
    for cell_id, (x, y) in lattice.items():
        set_xy(by_id[cell_id], x, y)

    # Remove old panel-b text and the correlation plot image.
    for cell in list(cells):
        cell_id = cell.get("id") or ""
        if cell_id.startswith("mod-b-") or cell_id == "JLeGss2_ZZTua6th3Dmy-149":
            cells.remove(cell)

    # Title above the lattice.
    add_text(
        cells, "mod-b-title",
        '<font style="font-size: 15px;"><b>Subgroup Lattice L(&#915;)</b></font>',
        -1820, 270, 340, 28, 15,
    )

    # Three subgroup-task rows, then the merge/dedup output and remaining details.
    add_box(
        cells, "mod-b-task1",
        "&#915;&nbsp;&nbsp;&#8594;&nbsp;&nbsp;DFS on Fix(&#915;)&nbsp;&nbsp;&#8594;&nbsp;&nbsp;&#954;(c)",
        -1820, 540, 340, 24, "#EAF2F8", "#4A7BA7", 10.5,
    )
    add_box(
        cells, "mod-b-task2",
        "H<sub>i</sub>&nbsp;&nbsp;&#8594;&nbsp;&nbsp;DFS on Fix(H<sub>i</sub>)"
        "&nbsp;&nbsp;&#8594;&nbsp;&nbsp;&#954;(c)",
        -1820, 572, 340, 24, "#EAF2F8", "#4A7BA7", 10.5,
    )
    add_box(
        cells, "mod-b-task3",
        "{e}&nbsp;&nbsp;&#8594;&nbsp;&nbsp;DFS on Fix({e}) = X"
        "&nbsp;&nbsp;&#8594;&nbsp;&nbsp;&#954;(c)",
        -1820, 604, 340, 24, "#EAF2F8", "#4A7BA7", 10.5,
    )
    add_box(
        cells, "mod-b-output",
        "merge identical P<sub>H</sub> &#8594; canonicalize &#8594; "
        "global dedup &#8594; &#937;",
        -1820, 644, 340, 28, "#FBF0E2", "#BE702A", 9.5,
    )
    add_text(
        cells, "mod-b-details",
        '<font style="font-size: 10.5px;">'
        'all subgroups; min-heap by |P<sub>H</sub>|<br>'
        'degree-count |C<sub>d</sub>|; canonicity<br>'
        '{e}: per-stratum quota B<sub>s</sub>'
        '</font>',
        -1820, 686, 340, 70, 10.5,
    )

    # Shift panels (c) and (d) left by 70, preserving their internal layout.
    for cell in cells.findall("mxCell"):
        if cell.get("parent") != "1":
            continue
        g = geom(cell)
        if g is None:
            continue
        try:
            x = float(g.get("x", "0"))
        except ValueError:
            continue
        if -1370 <= x <= -400:
            g.set("x", f"{x - 70:g}")

    # Text corrections that were previously applied only to the raster export.
    by_id["RH0wHt1T8WRH78aJM92P-414"].set(
        "value",
        '<font style="font-size: 16px;"><b style="">Discrete Algebraic Space</b><br>'
        'Proper-Operation Permutation<br>'
        'Group \\(\\Gamma = \\Phi(\\mathsf{S}_{\\mathrm{emb}}^{+})\\)</font>'
    )
    by_id["JLeGss2_ZZTua6th3Dmy-147"].set(
        "value",
        by_id["JLeGss2_ZZTua6th3Dmy-147"].get("value").replace(
            "struct.merge_sites(tol=0.5)", "struct.merge_sites(tol=0.2)"
        )
    )
    by_id["RH0wHt1T8WRH78aJM92P-472"].set(
        "value",
        '<font style="font-size: 17px;">Zero-redundancy library of pre-relaxation<br>'
        'models for downstream evaluation.</font>'
    )
    by_id["RH0wHt1T8WRH78aJM92P-474"].set(
        "value",
        '<span style="font-size: 14px;"><font style="color: rgb(29, 29, 29);">'
        'Exact rho = 1</font></span>'
    )
    by_id["RH0wHt1T8WRH78aJM92P-478"].set("value", "Exact rho in (1/|Γ|,1)")
    by_id["JLeGss2_ZZTua6th3Dmy-25"].set(
        "value",
        '<font style="color: rgb(29, 29, 29);">Exact rho = 1/|Γ|</font>'
    )
    by_id["JLeGss2_ZZTua6th3Dmy-29"].set(
        "value", "Trivial stabilizer.<br>Representatives retained."
    )

    tree.write(OUT, encoding="UTF-8", xml_declaration=True)
    print(OUT)


if __name__ == "__main__":
    main()
