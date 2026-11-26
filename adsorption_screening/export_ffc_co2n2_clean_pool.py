#!/usr/bin/env python3
"""Export a large strict-clean ffc CIF pool for CO2/N2 active learning."""

from __future__ import annotations

import argparse
import contextlib
import csv
import io
import json
import logging
import shutil
import sys
import time
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait
from collections import Counter
from pathlib import Path
from typing import Iterable

sys.path.insert(0, ".")

from pymatgen.core import Structure

from assembler import StructureAssembler
from adsorption_screening.export_ffc_co2n2_candidates import (
    FORCED_EDGE,
    MAX_ATOMS,
    NodeInfo,
    _check_passed,
    _copy_node_snapshots,
    _enumerate_ffc,
    _manifest_row,
    _node_info_map,
    _node_rows,
    _sorted_bucket,
    _write_csv,
    screen_nodes,
)


DEFAULT_BATCH_ID = "batch_002_clean_pool_5000"
DEFAULT_TARGET_CLEAN = 5000
DEFAULT_OUT_ROOT = Path("results/co2n2_ffc_screening")

CLEAN_BIASED_NODE_CHOICES = {
    3: ["3c_B_Ch", "3c_pyrrole_Ch", "3c_triazine_Ch"],
    4: [
        "4c_B_Ch",
        "4c_Cd_1_Ch",
        "4c_Cu_1_Ch",
        "4c_In_1_Ch",
        "4c_Ni_1_Ch",
        "4c_V_1_Ch",
        "4c_bicyclooctane_Ch",
        "4c_porphyrin_Ch",
        "4c_pyrene_Ch",
        "4c_tetrathiafulvalene_Ch",
    ],
}

FAST_CLEAN_NODE_CHOICES = {
    3: CLEAN_BIASED_NODE_CHOICES[3],
    4: [
        "4c_B_Ch",
        "4c_Cd_1_Ch",
        "4c_Cu_1_Ch",
        "4c_In_1_Ch",
        "4c_Ni_1_Ch",
        "4c_V_1_Ch",
        "4c_bicyclooctane_Ch",
        "4c_tetrathiafulvalene_Ch",
    ],
}

FULL_NODE_CHOICES = {
    3: ["3c_B_Ch", "3c_Cu_imadazoline_Ch", "3c_pyrrole_Ch", "3c_triazine_Ch"],
    4: CLEAN_BIASED_NODE_CHOICES[4],
}

NODE_SETS = {
    "clean-biased": CLEAN_BIASED_NODE_CHOICES,
    "fast-clean": FAST_CLEAN_NODE_CHOICES,
    "full": FULL_NODE_CHOICES,
}


def _read_csv(path: Path) -> list[dict[str, str]]:
    if not path.exists() or path.stat().st_size == 0:
        return []
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _append_csv(path: Path, row: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    exists = path.exists() and path.stat().st_size > 0
    with path.open("a", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(row.keys()))
        if not exists:
            writer.writeheader()
        writer.writerow(row)


def _record_key_from_row(row: dict[str, str]) -> tuple[int, tuple[str, ...], int]:
    return (
        int(row["source_enum_id"]),
        tuple(json.loads(row["coloring_json"])),
        int(row["stabilizer_size"]),
    )


def _record_key(record: dict) -> tuple[int, tuple[str, ...], int]:
    return (
        int(record["source_enum_id"]),
        tuple(record["coloring"]),
        int(record["stabilizer_size"]),
    )


def _candidate_order(records: Iterable[dict], mode: str) -> list[dict]:
    high = [record for record in records if record["sym_class"] == "High"]
    med = [record for record in records if record["sym_class"] == "Med"]
    low = [record for record in records if record["sym_class"] == "Low"]

    if mode == "symmetry_first":
        return _sorted_bucket(high) + _sorted_bucket(med) + _sorted_bucket(low)
    if mode == "balanced":
        buckets = [_sorted_bucket(high), _sorted_bucket(med), _sorted_bucket(low)]
        ordered: list[dict] = []
        while any(buckets):
            for bucket in buckets:
                if bucket:
                    ordered.append(bucket.pop(0))
        return ordered
    raise ValueError(f"Unsupported order mode: {mode}")


def _selected_node_infos(node_names: set[str]) -> list[NodeInfo]:
    return [info for info in screen_nodes() if info.name in node_names]


def _validation_warning(metadata: dict | None) -> str:
    if not metadata:
        return "assembly"
    validation = metadata.get("validation", {})
    if not isinstance(validation, dict):
        return "validation"
    failed = []
    for name in ("clash", "topology", "symmetry"):
        if _check_passed(validation, name) is False:
            failed.append(name)
    return ";".join(failed) or "validation"


def _write_summary(out_dir: Path, manifest: list[dict], rejected: list[dict], *,
                   node_choices: dict[int, list[str]], total_records: int, attempts: int,
                   elapsed: float, complete: bool) -> None:
    by_bucket = Counter(row["selection_bucket"] for row in manifest)
    by_class = Counter(row["sym_class"] for row in manifest)
    by_reject = Counter(row["reason"] for row in rejected)
    summary = {
        "batch_id": out_dir.name,
        "topology": "ffc",
        "edge": FORCED_EDGE,
        "complete": complete,
        "strict_clean_count": len(manifest),
        "attempts": attempts,
        "clean_rate": len(manifest) / attempts if attempts else 0.0,
        "elapsed_seconds": elapsed,
        "total_enumerated_records": total_records,
        "max_atoms": MAX_ATOMS,
        "node_choices": node_choices,
        "bucket_counts": dict(sorted(by_bucket.items())),
        "sym_class_counts": dict(sorted(by_class.items())),
        "rejected_count": len(rejected),
        "rejected_reason_counts": dict(sorted(by_reject.items())),
    }
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n")

    lines = [
        f"# {out_dir.name} ffc clean CO2/N2 active-learning pool",
        "",
        f"- Strict-clean CIFs: `{len(manifest)}`",
        f"- Attempts: `{attempts}`",
        f"- Clean rate: `{summary['clean_rate']:.3f}`",
        f"- Complete: `{complete}`",
        f"- Node choices: `{node_choices}`",
        f"- Bucket counts: `{dict(sorted(by_bucket.items()))}`",
        f"- Symmetry classes: `{dict(sorted(by_class.items()))}`",
        f"- Rejected candidates: `{len(rejected)}`",
        "",
        "This batch is a structural candidate pool for downstream CO2/N2 active learning. "
        "It does not run GCMC and does not make adsorption-performance claims.",
    ]
    (out_dir / "summary.md").write_text("\n".join(lines) + "\n")


def _assemble_one(task: dict) -> dict:
    record = task["record"]
    graph_data = task["graph_data"]
    tmp_path = Path(task["tmp_path"])
    structure_map = {str(i): name for i, name in enumerate(record["coloring"])}
    logging.getLogger().setLevel(logging.ERROR)
    assembler = StructureAssembler(output_dir=str(tmp_path.parent))

    metadata = None
    try:
        with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            metadata = assembler.assemble_and_save(
                str(tmp_path),
                structure_map,
                graph_data,
                forced_edge=FORCED_EDGE,
                return_metadata=True,
            )
        structure = Structure.from_file(tmp_path)
        if len(structure) > MAX_ATOMS:
            return {"record": record, "tmp_path": str(tmp_path), "metadata": metadata, "reason": f"n_atoms>{MAX_ATOMS}"}
        if not metadata.get("validation_passed", False):
            return {"record": record, "tmp_path": str(tmp_path), "metadata": metadata, "reason": _validation_warning(metadata)}
        return {"record": record, "tmp_path": str(tmp_path), "metadata": metadata, "reason": ""}
    except Exception as exc:
        return {
            "record": record,
            "tmp_path": str(tmp_path),
            "metadata": metadata,
            "reason": f"exception:{type(exc).__name__}:{exc}",
        }


def _reject_row(record: dict, rejected_index: int, reason: str, gamma: int) -> dict:
    return {
        "rejected_id": f"ffc_co2n2_rejected_{rejected_index:06d}",
        "source_enum_id": record["source_enum_id"],
        "composition": record["composition"],
        "rho": f"{record['rho']:.6f}",
        "sym_class": record["sym_class"],
        "selection_bucket": record["selection_bucket"],
        "stabilizer_size": record["stabilizer_size"],
        "gamma_size": gamma,
        "coloring_json": json.dumps(record["coloring"], sort_keys=False),
        "reason": reason,
    }


def _handle_result(result: dict, *, manifest: list[dict], rejected: list[dict],
                   manifest_path: Path, rejected_path: Path, cif_dir: Path,
                   node_map: dict[str, NodeInfo], gamma: int, clean_index: int,
                   rejected_index: int, keep_rejected_cifs: bool) -> tuple[int, int]:
    record = result["record"]
    tmp_path = Path(result["tmp_path"])
    reason = result["reason"]
    metadata = result.get("metadata")

    if reason:
        rejected_index += 1
        if result.get("tmp_path") and tmp_path.is_file() and not keep_rejected_cifs:
            tmp_path.unlink(missing_ok=True)
        row = _reject_row(record, rejected_index, reason, gamma)
        rejected.append(row)
        _append_csv(rejected_path, row)
        return clean_index, rejected_index

    clean_index += 1
    candidate_id = f"ffc_co2n2_cleanpool_{clean_index:05d}"
    clean_path = cif_dir / f"{candidate_id}_{record['selection_bucket']}_rho{record['rho']:.4f}_{record['sym_class']}.cif"
    tmp_path.replace(clean_path)
    structure = Structure.from_file(clean_path)
    row = _manifest_row(record, candidate_id, clean_path, structure, node_map, metadata)
    row["gamma_size"] = gamma
    manifest.append(row)
    _append_csv(manifest_path, row)
    return clean_index, rejected_index


def export_clean_pool(batch_id: str = DEFAULT_BATCH_ID,
                      target_clean: int = DEFAULT_TARGET_CLEAN,
                      out_root: Path = DEFAULT_OUT_ROOT,
                      node_set: str = "clean-biased",
                      order: str = "balanced",
                      max_attempts: int | None = None,
                      resume: bool = False,
                      keep_rejected_cifs: bool = False,
                      progress_every: int = 100,
                      workers: int = 1) -> list[dict]:
    node_choices = NODE_SETS[node_set]
    records, graph_data, gamma = _enumerate_ffc(node_choices)
    ordered = _candidate_order(records, order)
    node_names = {name for names in node_choices.values() for name in names}
    node_infos = _selected_node_infos(node_names)
    node_map = _node_info_map(node_infos)

    out_dir = out_root / batch_id
    cif_dir = out_dir / "cifs"
    rejected_dir = out_dir / "rejected_cifs"
    manifest_path = out_dir / "manifest.csv"
    rejected_path = out_dir / "rejected_manifest.csv"
    out_dir.mkdir(parents=True, exist_ok=True)
    cif_dir.mkdir(parents=True, exist_ok=True)
    if keep_rejected_cifs:
        rejected_dir.mkdir(parents=True, exist_ok=True)

    if not resume and (manifest_path.exists() or rejected_path.exists()):
        raise FileExistsError(f"{out_dir} already contains output; pass --resume or choose a new --batch-id")

    _copy_node_snapshots(node_infos, out_dir / "nodes_snapshot")
    _write_csv(out_dir / "nodes_screening.csv", _node_rows(node_infos))

    manifest = _read_csv(manifest_path)
    rejected = _read_csv(rejected_path)
    seen = {_record_key_from_row(row) for row in manifest + rejected}
    attempts = len(manifest) + len(rejected)
    start = time.time()
    clean_index = len(manifest)
    rejected_index = len(rejected)

    logging.getLogger().setLevel(logging.ERROR)
    if workers > 1:
        record_iter = iter(ordered)
        futures = {}
        submitting = True
        with ProcessPoolExecutor(max_workers=workers) as pool:
            while futures or submitting:
                while submitting and len(futures) < workers and len(manifest) < target_clean:
                    if max_attempts is not None and attempts >= max_attempts:
                        submitting = False
                        break
                    for base_record in record_iter:
                        if _record_key(base_record) not in seen:
                            break
                    else:
                        submitting = False
                        break

                    attempts += 1
                    bucket = f"clean_pool_{base_record['sym_class'].lower()}"
                    record = {**base_record, "selection_bucket": bucket}
                    tmp_path = (
                        rejected_dir / f"attempt_{attempts:06d}.cif"
                        if keep_rejected_cifs
                        else cif_dir / f"_attempt_{attempts:06d}.cif"
                    )
                    future = pool.submit(
                        _assemble_one,
                        {"record": record, "graph_data": graph_data, "tmp_path": str(tmp_path)},
                    )
                    futures[future] = record
                    seen.add(_record_key(record))

                if not futures:
                    break

                done, _pending = wait(futures, return_when=FIRST_COMPLETED)
                for future in done:
                    record = futures.pop(future)
                    previous_clean = len(manifest)
                    try:
                        result = future.result()
                    except Exception as exc:
                        result = {
                            "record": record,
                            "tmp_path": "",
                            "metadata": None,
                            "reason": f"exception:{type(exc).__name__}:{exc}",
                        }
                    clean_index, rejected_index = _handle_result(
                        result,
                        manifest=manifest,
                        rejected=rejected,
                        manifest_path=manifest_path,
                        rejected_path=rejected_path,
                        cif_dir=cif_dir,
                        node_map=node_map,
                        gamma=gamma,
                        clean_index=clean_index,
                        rejected_index=rejected_index,
                        keep_rejected_cifs=keep_rejected_cifs,
                    )
                    if (
                        progress_every
                        and len(manifest) > previous_clean
                        and len(manifest) % progress_every == 0
                    ):
                        print(
                            f"[clean-pool] clean={len(manifest)} attempts={attempts} "
                            f"rate={len(manifest) / attempts:.3f}",
                            flush=True,
                        )
                if len(manifest) >= target_clean:
                    submitting = False
        elapsed = time.time() - start
        complete = len(manifest) >= target_clean
        (out_dir / "manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n")
        _write_summary(
            out_dir,
            manifest,
            rejected,
            node_choices=node_choices,
            total_records=len(records),
            attempts=attempts,
            elapsed=elapsed,
            complete=complete,
        )
        return manifest

    assembler = StructureAssembler(output_dir=str(cif_dir))
    for record in ordered:
        if len(manifest) >= target_clean:
            break
        if max_attempts is not None and attempts >= max_attempts:
            break
        if _record_key(record) in seen:
            continue

        attempts += 1
        bucket = f"clean_pool_{record['sym_class'].lower()}"
        record = {**record, "selection_bucket": bucket}
        structure_map = {str(i): name for i, name in enumerate(record["coloring"])}
        tmp_path = rejected_dir / f"attempt_{attempts:06d}.cif" if keep_rejected_cifs else cif_dir / f"_attempt_{attempts:06d}.cif"
        metadata = None

        try:
            with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
                metadata = assembler.assemble_and_save(
                    str(tmp_path),
                    structure_map,
                    graph_data,
                    forced_edge=FORCED_EDGE,
                    return_metadata=True,
                )
            structure = Structure.from_file(tmp_path)
            if len(structure) > MAX_ATOMS:
                reason = f"n_atoms>{MAX_ATOMS}"
            elif not metadata.get("validation_passed", False):
                reason = _validation_warning(metadata)
            else:
                clean_index += 1
                candidate_id = f"ffc_co2n2_cleanpool_{clean_index:05d}"
                clean_path = cif_dir / f"{candidate_id}_{bucket}_rho{record['rho']:.4f}_{record['sym_class']}.cif"
                tmp_path.replace(clean_path)
                row = _manifest_row(record, candidate_id, clean_path, structure, node_map, metadata)
                row["gamma_size"] = gamma
                manifest.append(row)
                _append_csv(manifest_path, row)
                seen.add(_record_key(record))
                if progress_every and len(manifest) % progress_every == 0:
                    print(
                        f"[clean-pool] clean={len(manifest)} attempts={attempts} "
                        f"rate={len(manifest) / attempts:.3f}",
                        flush=True,
                    )
                continue
        except Exception as exc:
            reason = f"exception:{type(exc).__name__}:{exc}"

        rejected_index += 1
        if tmp_path.exists() and not keep_rejected_cifs:
            tmp_path.unlink(missing_ok=True)
        rejected_row = _reject_row(record, rejected_index, reason, gamma)
        rejected.append(rejected_row)
        _append_csv(rejected_path, rejected_row)
        seen.add(_record_key(record))

    elapsed = time.time() - start
    complete = len(manifest) >= target_clean
    (out_dir / "manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n")
    _write_summary(
        out_dir,
        manifest,
        rejected,
        node_choices=node_choices,
        total_records=len(records),
        attempts=attempts,
        elapsed=elapsed,
        complete=complete,
    )
    return manifest


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--batch-id", default=DEFAULT_BATCH_ID)
    parser.add_argument("--target-clean", type=int, default=DEFAULT_TARGET_CLEAN)
    parser.add_argument("--out-root", default=str(DEFAULT_OUT_ROOT))
    parser.add_argument("--node-set", choices=sorted(NODE_SETS), default="clean-biased")
    parser.add_argument("--order", choices=["balanced", "symmetry_first"], default="balanced")
    parser.add_argument("--max-attempts", type=int)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--keep-rejected-cifs", action="store_true")
    parser.add_argument("--progress-every", type=int, default=100)
    parser.add_argument("--workers", type=int, default=1)
    args = parser.parse_args(argv)

    rows = export_clean_pool(
        batch_id=args.batch_id,
        target_clean=args.target_clean,
        out_root=Path(args.out_root),
        node_set=args.node_set,
        order=args.order,
        max_attempts=args.max_attempts,
        resume=args.resume,
        keep_rejected_cifs=args.keep_rejected_cifs,
        progress_every=args.progress_every,
        workers=args.workers,
    )
    print(f"[clean-pool] prepared {len(rows)} strict-clean CIFs")
    print(f"[clean-pool] output: {Path(args.out_root) / args.batch_id}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
