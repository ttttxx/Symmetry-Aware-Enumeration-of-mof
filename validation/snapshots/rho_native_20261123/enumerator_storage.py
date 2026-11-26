"""Full-key storage alternatives for the experimental subgroup DFS.

Only the completed-canonical membership index changes; scheduling, prefix
pruning, result records, sorting and label decoding remain inherited. These
indexes have no subtree-coverage meaning and do not replace DFS pruning.
"""

from enumerator_casm_inspired import CasmInspiredMOFEnumerator


class HashStorageMOFEnumerator(CasmInspiredMOFEnumerator):
    def __init__(self, *args, storage_key="bytes", **kwargs):
        if storage_key not in ("tuple", "bytes"):
            raise ValueError("storage_key must be 'tuple' or 'bytes'")
        super().__init__(*args, **kwargs)
        self.storage_key = storage_key
        self._canonical_index = {}
        self._use_bytes_key = False

    def enumerate_subgroups(self, atom_types):
        self._canonical_index = {}
        # Integer codes 0..255 fit in one byte; larger label sets retain exact
        # tuple keys. No digest or truncated key is used, so no collision-based
        # loss of representatives is possible.
        self._use_bytes_key = self.storage_key == "bytes" and len(self._code_to_label) <= 256
        return super().enumerate_subgroups(atom_types)

    def _storage_key(self, canonical):
        if hasattr(self, "_label_to_code") and not self._labels_encoded:
            try:
                canonical = tuple(self._label_to_code[label] for label in canonical)
            except KeyError:
                return None
        return bytes(canonical) if self._use_bytes_key else tuple(canonical)

    def _trie_lookup(self, canonical):
        key = self._storage_key(canonical)
        return self._canonical_index.get(key) if key is not None else None

    def _trie_insert(self, canonical, class_id, stabilizer_size):
        key = self._storage_key(canonical)
        if key is None:
            raise ValueError("Canonical coloring contains an unknown label")
        self._canonical_index[key] = class_id
