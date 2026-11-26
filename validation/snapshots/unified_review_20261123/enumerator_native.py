"""Python scheduler + C++ search with compact full-representative Trie.

Build explicitly with ``python tools/build_enumeration_core.py``. ``auto`` uses
the shared library when available; ``python`` forces the tested Python engine;
``cpp`` requires the native engine. No compilation occurs while importing.
"""

import ctypes as ct
from pathlib import Path
import sys
import time

from enumerator_storage import HashStorageMOFEnumerator

_I = ct.c_int32
_U = ct.c_uint64
_U32 = ct.c_uint32
_PI, _PU = ct.POINTER(_I), ct.POINTER(_U)


def _load_core():
    suffix = ".dylib" if sys.platform == "darwin" else ".so"
    path = Path(__file__).resolve().parent / "native" / ("libmof_enumeration" + suffix)
    lib = ct.CDLL(str(path))
    lib.mof_core_version.restype = _I
    if lib.mof_core_version() != 4:
        raise OSError("Native core ABI mismatch; rebuild the optional library")
    lib.mof_core_error.restype = ct.c_char_p
    lib.mof_core_create.argtypes = [_I, _I, _PI]
    lib.mof_core_create.restype = ct.c_void_p
    lib.mof_core_destroy.argtypes = [ct.c_void_p]
    lib.mof_core_target_order.argtypes = [ct.c_void_p, _I]
    lib.mof_core_target_order.restype = _I
    lib.mof_core_count.argtypes = [ct.c_void_p]
    lib.mof_core_count.restype = _U
    lib.mof_core_storage_mode.argtypes = [ct.c_void_p, _I]
    lib.mof_core_storage_mode.restype = _I
    lib.mof_core_storage_stats.argtypes = [ct.c_void_p, _PU]
    lib.mof_core_storage_stats.restype = None
    lib.mof_core_trie_copy_nodes.argtypes = [ct.c_void_p, _U, ct.POINTER(_U32)]
    lib.mof_core_trie_copy_nodes.restype = _I
    lib.mof_core_trie_lookup.argtypes = [ct.c_void_p, ct.POINTER(ct.c_ubyte)]
    lib.mof_core_trie_lookup.restype = ct.c_int64
    lib.mof_core_prefix_mode.argtypes = [ct.c_void_p, _I, _U]
    lib.mof_core_prefix_mode.restype = _I
    lib.mof_core_prefix_stats.argtypes = [ct.c_void_p, _PU]
    lib.mof_core_prefix_stats.restype = None
    lib.mof_core_run.argtypes = [ct.c_void_p, _I, _PI, _PI, _PI, _PI, _I,
                                 _PU, _PU, _PI, _PI, _PI, _PU, ct.c_int64,
                                 ct.c_double, _I, _PU]
    lib.mof_core_run.restype = _I
    lib.mof_core_copy.argtypes = [ct.c_void_p, _U, _U, ct.POINTER(ct.c_ubyte), _PI]
    lib.mof_core_copy.restype = _I
    return lib


def native_available():
    try:
        _load_core()
        return True
    except OSError:
        return False


def _ints(values):
    values = list(values)
    return (_I * len(values))(*values)


def _words(values, width):
    values = [((value >> (64 * word)) & ((1 << 64)-1))
              for value in values for word in range(width)]
    return (_U * len(values))(*values)


class HybridMOFEnumerator(HashStorageMOFEnumerator):
    def __init__(self, *args, backend="auto", native_kernel="dfs", native_store="trie", **kwargs):
        if backend not in ("auto", "python", "cpp"):
            raise ValueError("backend must be auto, python or cpp")
        if native_kernel not in ("dfs", "counter"):
            raise ValueError("native_kernel must be dfs or counter")
        if native_store not in ("trie", "hash"):
            raise ValueError("native_store must be trie or hash")
        super().__init__(*args, **kwargs)
        self.backend = backend
        self.native_kernel = native_kernel
        self.native_store = native_store
        self.active_store = "python_hash"
        self.result_index_stats = {}
        self._native_trie_nodes = None
        if self.prefix_trie and native_kernel != "dfs":
            raise ValueError("Completed-prefix Trie is supported by the DFS kernel only")
        self.active_backend = "python"
        self.backend_reason = "Python explicitly selected"
        self._native_context = None
        self.native_task_reports = []

    def enumerate_subgroups(self, atom_types):
        self.active_backend = "python"
        self.active_store = "python_hash"
        self.result_index_stats = {}
        self._native_trie_nodes = None
        self.native_task_reports = []
        self.backend_reason = "Python explicitly selected" if self.backend == "python" else "Native library unavailable"
        supported = 0 < self.n <= 512 and len(self._code_to_label) <= 256
        lib = None
        if self.backend != "python" and supported:
            try:
                lib = _load_core()
            except OSError as exc:
                self.backend_reason = str(exc)
        elif not supported:
            self.backend_reason = "Native engine supports 1..512 sites and at most 256 labels"
        if self.backend == "cpp" and lib is None:
            raise RuntimeError(self.backend_reason + "; build with python tools/build_enumeration_core.py")
        if lib is None:
            if self.native_kernel != "dfs":
                raise ValueError("The comparison counter kernel requires backend='cpp'")
            return self._run_schedule(atom_types)
        inverse_data = _ints(value for inverse in self.perms_inv for value in inverse)
        context = lib.mof_core_create(self.n, len(self.perms), inverse_data)
        if not context:
            raise RuntimeError(lib.mof_core_error().decode())
        self._native_context = context
        self._core = lib
        lib.mof_core_storage_mode(context, int(self.native_store == "trie"))
        lib.mof_core_prefix_mode(context, int(self.prefix_trie), self.prefix_trie_max_nodes)
        self.active_backend = "cpp"
        self.active_store = "cpp_" + self.native_store
        self.backend_reason = "Optional native library loaded"
        try:
            return self._run_schedule(atom_types)
        finally:
            try:
                if self.prefix_trie:
                    self._update_prefix_stats()
                    self.prefix_trie_stats["retained_nodes_before_release"] = self.prefix_trie_stats["live_nodes"]
                stats = (_U * 7)()
                lib.mof_core_storage_stats(context, stats)
                self.result_index_stats = dict(zip(
                    ("trie_enabled", "nodes", "capacity", "node_bytes", "capacity_bytes",
                     "hash_keys", "hash_buckets"), map(int, stats)))
                # Retain the flat index for the existing post-run lookup API.
                # Edges refer to classes by stable insertion ID; no Python
                # full-key hash mirror is built in Trie mode.
                if self.native_store == "trie" and lib.mof_core_count(context) == len(self.classes):
                    self._native_trie_nodes = (_U32 * (4 * stats[1]))()
                    if lib.mof_core_trie_copy_nodes(context, stats[1], self._native_trie_nodes):
                        raise RuntimeError(lib.mof_core_error().decode())
            finally:
                lib.mof_core_destroy(context)
                self._native_context = None
                if self.prefix_trie:
                    self.prefix_trie_stats["live_nodes"] = 0

    def _run_schedule(self, atom_types):
        return super().enumerate_subgroups(atom_types)

    def _trie_lookup(self, canonical):
        if self.active_backend != "cpp" or self.native_store != "trie":
            return super()._trie_lookup(canonical)
        key = self._storage_key(canonical)
        if key is None or len(key) != self.n:
            return None
        if self._native_context is not None:
            raw = (ct.c_ubyte * self.n).from_buffer_copy(bytes(key))
            result = self._core.mof_core_trie_lookup(self._native_context, raw)
            return int(result) if result >= 0 else None
        nodes = self._native_trie_nodes
        if nodes is None:
            return None
        nil, parent = (1 << 32) - 1, 0
        while True:
            depth, child = nodes[parent * 4 + 1], nodes[parent * 4 + 2]
            while child != nil:
                sample, end = nodes[child * 4], nodes[child * 4 + 1]
                coloring = self.classes[sample]["coloring"]
                first = coloring[depth]
                if not self._labels_encoded:
                    first = self._label_to_code[first]
                if first >= key[depth]:
                    break
                child = nodes[child * 4 + 3]
            if child == nil:
                return None
            for position in range(depth, end):
                value = coloring[position]
                if not self._labels_encoded:
                    value = self._label_to_code[value]
                if value != key[position]:
                    return None
            if end == self.n:
                return int(sample)
            parent = child

    def _update_prefix_stats(self):
        stats = (_U * 10)()
        self._core.mof_core_prefix_stats(self._native_context, stats)
        self.prefix_trie_stats = dict(zip(
            ("queries", "hits", "early_hits", "closed_prefixes", "nodes_created",
             "live_nodes", "peak_nodes", "budget_skips", "snapshots", "forests"),
            map(int, stats)))

    def _search_compiled_task(self, orbits, atom_types, coloring,
                              active_perm_indices, current_degree_counts,
                              identity_start, identity_budget, deadline):
        if self._native_context is None:
            return super()._search_compiled_task(
                orbits, atom_types, coloring, active_perm_indices, current_degree_counts,
                identity_start, identity_budget, deadline)
        if self.topology_requirements and self.node_degrees:
            counts = dict(current_degree_counts or {degree: 0 for degree in self.topology_requirements})
            for orbit in orbits:
                degree = self.node_degrees[orbit[0]]
                counts[degree] += len(orbit)
                if counts[degree] > self.topology_requirements.get(degree, counts[degree]):
                    return
        target_order = getattr(self, "_active_rho_order", None) or 0
        if self._core.mof_core_target_order(self._native_context, target_order):
            raise RuntimeError(self._core.mof_core_error().decode())
        domains = [self._choices(orbit, atom_types) for orbit in orbits]
        nodes, offsets, choices, choice_offsets = [], [0], [], [0]
        for orbit, domain in zip(orbits, domains):
            nodes.extend(orbit)
            offsets.append(len(nodes))
            choices.extend(domain)
            choice_offsets.append(len(choices))
        kinds, others, masks, event_offsets = [], [], [], [0]
        for events in self._exclusion_events:
            for (kind, other), mask in events:
                kinds.append(0 if kind == "constant" else 1)
                others.append(other)
                masks.append(mask)
            event_offsets.append(len(kinds))
        initial = self._exclusion_stack[0]
        width = max(1, (initial.bit_length()+63)//64)
        limits = [self.max_total_results]
        if identity_start is not None and identity_budget is not None:
            limits.append(identity_start + identity_budget)
        before = len(self.classes)
        finite_limits = [value for value in limits if value is not None]
        if finite_limits and min(finite_limits) <= before:
            return
        limit = min(min(finite_limits), (1 << 63)-1) if finite_limits else -1
        if self._core.mof_core_count(self._native_context) != before:
            raise RuntimeError("Python/native result indexes diverged")
        arguments = (
            self._native_context, len(orbits), _ints(nodes), _ints(offsets),
            _ints(choices), _ints(choice_offsets), width, _words([initial], width),
            _words(self._decided_masks, width), _ints(event_offsets), _ints(kinds),
            _ints(others), _words(masks, width))
        stats = (_U * 7)()
        remaining = max(0., deadline-time.monotonic()) if deadline is not None else -1.
        started = time.perf_counter()
        status = self._core.mof_core_run(*arguments, limit, remaining,
                                        0 if self.native_kernel == "dfs" else 1, stats)
        kernel_seconds = time.perf_counter()-started
        if status:
            raise RuntimeError(self._core.mof_core_error().decode())
        count = self._core.mof_core_count(self._native_context)-before
        key_buffer = (ct.c_ubyte * (count*self.n))()
        stabilizers = (_I * count)()
        if self._core.mof_core_copy(self._native_context, before, count, key_buffer, stabilizers):
            raise RuntimeError(self._core.mof_core_error().decode())
        raw = bytes(key_buffer)
        for i, stabilizer in enumerate(stabilizers):
            key = raw[i*self.n:(i+1)*self.n]
            cid = before+i
            if self.native_store == "hash":
                self._canonical_index[key if self._use_bytes_key else tuple(key)] = cid
            self.classes[cid] = dict(id=cid, coloring=list(key), stabilizer_size=stabilizer,
                                     orbit_size=len(self.perms)//stabilizer)
        for name, value in zip(("leaf_visits", "duplicate_leaves", "pruned_prefixes"), stats):
            self.task_dedup_stats[name] += int(value)
        self.task_dedup_stats["search_nodes"] = self.task_dedup_stats.get("search_nodes", 0) + int(stats[6])
        if self.prefix_trie:
            self._update_prefix_stats()
        self.native_task_reports.append(dict(kernel_s=kernel_seconds, added=count,
            candidates=int(stats[3]), timed_out=bool(stats[4]), output_capped=bool(stats[5]),
            search_nodes=int(stats[6])))
