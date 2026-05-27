"""
Streamlit frontend for MOF Structure Enumerator.
Run with:  streamlit run app.py
"""
from __future__ import annotations
import contextlib
import json
import os
import sys
import traceback
from io import StringIO
from pathlib import Path

import streamlit as st

# ensure imports resolve relative to this file's directory
PROJECT_ROOT = Path(__file__).parent.resolve()
os.chdir(PROJECT_ROOT)
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

# ── topology auto-analysis ────────────────────────────────────────────────
from main import calculate_topology_requirements, load_graph

# ── 3-D visualisation ─────────────────────────────────────────────────────────
try:
    import py3Dmol  # noqa: F401
    HAS_PY3DMOL = True
except ImportError:
    HAS_PY3DMOL = False

try:
    from stmol import showmol  # noqa: F401
    HAS_STMOL = True
except ImportError:
    HAS_STMOL = False


# ── helpers ───────────────────────────────────────────────────────────────────

def _topology_list() -> list[str]:
    return sorted(d.name for d in (PROJECT_ROOT / "data").iterdir() if d.is_dir())


def _node_list() -> list[str]:
    return sorted(f.stem for f in (PROJECT_ROOT / "nodes").glob("*.cif"))


def _parse_degree(name: str) -> int | None:
    """Extract coordination number from a node name like '4c_Cu_1_Ch' -> 4."""
    import re
    m = re.match(r'(\d+)c_', name)
    return int(m.group(1)) if m else None


def _nodes_by_degree() -> dict[int, list[str]]:
    """Group available node CIFs by their coordination degree."""
    groups: dict[int, list[str]] = {}
    for name in _node_list():
        d = _parse_degree(name)
        if d is not None:
            groups.setdefault(d, []).append(name)
    return groups


def _edge_list() -> list[str]:
    return sorted(f.stem for f in (PROJECT_ROOT / "edges").glob("*.cif"))


def _show_cif_3d(cif_path: str, style: str = "Ball and Stick", height: int = 420) -> None:
    """Render a CIF file as an interactive 3-D viewer inside Streamlit."""
    with open(cif_path) as fh:
        cif_str = fh.read()

    # json.dumps escapes safely
    cif_json = json.dumps(cif_str)

    # 匹配 MOF 对应的 3Dmol 渲染样式
    if style == "Stick":
        style_js = "stick: { radius: 0.15, colorscheme: 'Jmol' }"
    elif style == "Sphere":
        style_js = "sphere: { radius: 1.0, colorscheme: 'Jmol' }"
    elif style == "Line":
        style_js = "line: { colorscheme: 'Jmol' }"
    else: # 默认 Ball and Stick
        style_js = "stick: { radius: 0.12, colorscheme: 'Jmol' }, sphere: { radius: 0.28, colorscheme: 'Jmol' }"

    html = f"""<!DOCTYPE html>
<html>
<head>
  <script src="https://cdnjs.cloudflare.com/ajax/libs/3Dmol/2.1.0/3Dmol-min.js"></script>
  <style>
    html, body {{ margin: 0; padding: 0; width: 100%; height: 100%; overflow: hidden;
                  background: #ffffff; }}
    #viewer    {{ width: 100%; height: {height}px; position: relative; }}
  </style>
</head>
<body>
  <div id="viewer"></div>
  <script>
    var viewer = $3Dmol.createViewer(
      document.getElementById('viewer'),
      {{ backgroundColor: '#ffffff' }}
    );
    viewer.addModel({cif_json}, 'cif');
    viewer.setStyle({{}}, {{
      {style_js}
    }});
    viewer.zoomTo();   // centre + fit before render
    viewer.render();
  </script>
</body>
</html>"""

    st.components.v1.html(html, height=height + 2, scrolling=False)

# ── page config ───────────────────────────────────────────────────────────────

st.set_page_config(
    page_title="MOF Structure Enumerator",
    page_icon="🔬",
    layout="wide",
    initial_sidebar_state="expanded",
)

st.title("MOF Structure Enumerator")
st.caption(
    "Enumerate symmetry-unique MOF colorings, verify chemistry, and visualise top-ranked CIF structures."
)

# ── sidebar inputs ────────────────────────────────────────────────────────────

with st.sidebar:
    st.header("Parameters")

    topologies = _topology_list()
    default_topo_idx = topologies.index("ffc") if "ffc" in topologies else 0
    topology = st.selectbox("Topology (DATA_PREFIX)", topologies, index=default_topo_idx)

    # ── 拓扑需求自动分析 + 动态配体选择 ──────────────────────────────────
    user_selected_ligands: dict[int, list[str]] = {}
    nodes_by_deg = _nodes_by_degree()

    if topology:
        graph_path = os.path.join("data", topology, f"{topology}_graph.json")
        if os.path.exists(graph_path):
            with open(graph_path) as f:
                gdata = json.load(f)
            topo_reqs = calculate_topology_requirements(gdata)

            req_text = "Current topology requires: " + ", ".join(
                f"{cnt} × {deg}-c node{'s' if cnt > 1 else ''}"
                for deg, cnt in sorted(topo_reqs.items())
            )
            st.caption(req_text)

            for degree, required_count in sorted(topo_reqs.items()):
                available = nodes_by_deg.get(degree, [])
                if not available:
                    st.warning(f"No {degree}-c ligands available in library!")
                    continue
                # 智能默认：优先选中 Cu 配体
                default = [n for n in available if "Cu" in n]
                if not default:
                    default = available[:1]

                selected = st.multiselect(
                    f"**{degree}-c slots** (need {required_count})",
                    options=available,
                    default=default,
                    key=f"ligand_sel_{degree}",
                )
                if selected:
                    user_selected_ligands[degree] = selected

    # ── 边缘 (Edge) 选择 ──────────────────────────────────────────────────
    all_edges = _edge_list()
    default_edge_idx = all_edges.index("1B_4H_Ch") if "1B_4H_Ch" in all_edges else 0
    forced_edge = st.selectbox(
        "Forced Edge (FORCED_EDGE)",
        all_edges,
        index=default_edge_idx,
        help="The linker CIF used for edge assembly.",
    )

    st.divider()
    max_cifs = st.slider("Top-K CIF structures to generate", min_value=1, max_value=10, value=5)
    max_results = st.number_input(
        "Max enumeration results",
        min_value=100,
        max_value=100_000,
        value=10_000,
        step=1_000,
    )

    st.divider()
    run_btn = st.button("Generate Structures", type="primary", use_container_width=True)

# ── session state ─────────────────────────────────────────────────────────────

if "run_results" not in st.session_state:
    st.session_state.run_results = None

# ── run enumeration ───────────────────────────────────────────────────────────

if run_btn:
    if not user_selected_ligands:
        st.error("Please select at least one ligand for each degree slot.")
        st.stop()

    # Validate every topology requirement has a ligand selection
    log_buf = StringIO()
    progress_placeholder = st.empty()

    try:
        with progress_placeholder.container():
            st.info("Running enumeration — this may take a few minutes...")

        from main import run_enumeration

        with contextlib.redirect_stdout(log_buf):
            cif_paths, json_path, output_data = run_enumeration(
                data_prefix=topology,
                user_selected_ligands=user_selected_ligands,
                forced_edge=forced_edge,
                max_cifs=int(max_cifs),
                max_results=int(max_results),
            )

        st.session_state.run_results = {
            "cif_paths": cif_paths,
            "json_path": json_path,
            "output_data": output_data,
            "log": log_buf.getvalue(),
        }
        progress_placeholder.success(
            f"Done! Generated {len(cif_paths)} CIF file(s). "
            f"({len(output_data.get('structures', []))} valid structures total)"
        )

    except FileNotFoundError as exc:
        progress_placeholder.error(str(exc))
        st.stop()
    except Exception:
        progress_placeholder.error("Enumeration failed — see traceback below.")
        st.code(traceback.format_exc(), language="python")
        if log_buf.getvalue():
            with st.expander("Partial run log"):
                st.text(log_buf.getvalue())
        st.stop()

    with st.expander("Run log", expanded=False):
        st.text(st.session_state.run_results["log"])

# ── display results ───────────────────────────────────────────────────────────

if st.session_state.run_results:
    res = st.session_state.run_results
    output_data = res["output_data"]
    meta = output_data.get("meta", {})
    structures = output_data.get("structures", [])
    cif_paths = res["cif_paths"]

    # ── summary metrics ──────────────────────────────────────────────────────
    st.divider()
    c1, c2, c3, c4 = st.columns(4)
    c1.metric("Topology", meta.get("topology", "-"))
    c2.metric("Group Size |G|", meta.get("group_size", "-"))
    c3.metric("Valid Structures", len(structures))
    c4.metric("CIFs Generated", len(cif_paths))

    verify = meta.get("verification", {})
    if verify:
        with st.expander("Verification stats"):
            st.json(verify)

    # ── 3-D viewers ──────────────────────────────────────────────────────────
    st.subheader("3-D Structure Viewer — Top Ranked CIFs")

    if not cif_paths:
        st.warning("No CIF files were generated for this run.")
    else:
        # 新增：3D 渲染样式切换菜单
        render_style = st.selectbox(
            "Select 3D Rendering Style",
            ["Ball and Stick", "Stick", "Sphere", "Line"],
            index=0,
            help="Sphere is useful for visualizing pore sizes, while Line/Stick is better for topology."
        )
        st.write("")  # 增加一点留白

        # Layout: up to 3 columns per row
        n_cols = min(3, len(cif_paths))
        for row_start in range(0, len(cif_paths), n_cols):
            batch = cif_paths[row_start: row_start + n_cols]
            cols = st.columns(len(batch))
            for col, path in zip(cols, batch):
                with col:
                    fname = Path(path).name
                    rank = fname.split("_rank")[-1].split("_")[0] if "_rank" in fname else "?"
                    sym = fname.split("_sym")[-1].replace(".cif", "") if "_sym" in fname else "?"
                    st.markdown(f"**Rank {rank}** &nbsp; | &nbsp; Stabilizer = {sym}")

                    # 传入选中的样式
                    _show_cif_3d(path, style=render_style, height=380)

                    # per-structure download button
                    with open(path, "rb") as fh:
                        st.download_button(
                            label=f"Download {fname}",
                            data=fh,
                            file_name=fname,
                            mime="chemical/x-cif",
                            key=f"dl_cif_{path}",
                        )
    # ── structure table ───────────────────────────────────────────────────────
    st.subheader("All Valid Structures")

    table_rows = [
        {
            "id": s.get("id", i),
            "Coloring": " | ".join(s.get("coloring", [])),
            "Stabilizer Size": s.get("stabilizer_size"),
            "Orbit Size": s.get("orbit_size"),
        }
        for i, s in enumerate(structures)
    ]
    st.dataframe(table_rows, use_container_width=True, height=300)

    # ── JSON display + download ───────────────────────────────────────────────
    st.subheader("Full JSON Output")

    with st.expander("Preview JSON (first 20 structures)", expanded=False):
        preview = dict(output_data)
        preview["structures"] = structures[:20]
        st.json(preview)

    json_bytes = json.dumps(output_data, indent=2).encode()
    json_filename = Path(res["json_path"]).name

    st.download_button(
        label="Download Complete JSON",
        data=json_bytes,
        file_name=json_filename,
        mime="application/json",
        type="primary",
    )

    # ── run log ───────────────────────────────────────────────────────────────
    with st.expander("Run log", expanded=False):
        st.text(res["log"])
