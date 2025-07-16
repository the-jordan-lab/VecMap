# flake8: noqa: E203
"""Numba-accelerated mapper functions for VecMap.

This module demonstrates a prototype of VecMap 2.0 where the seed index
construction and candidate scanning loops are compiled with Numba to
remove Python interpreter overhead. It is not fully optimized but
serves as a starting point for GPU/SIMD exploration.
"""

from typing import Iterable, List, Tuple

import numpy as np
from numba import njit, types
from numba.typed import Dict, List as NumbaList


@njit
def _build_seed_index_nb(ref: str, seed_len: int) -> Dict:
    index = Dict.empty(
        key_type=types.unicode_type,
        value_type=types.ListType(types.int64),
    )
    for i in range(len(ref) - seed_len + 1):
        seed = ref[i : i + seed_len]
        if seed in index:
            index[seed].append(i)
        else:
            lst = NumbaList.empty_list(types.int64)
            lst.append(i)
            index[seed] = lst
    return index


@njit
def _get_candidate_starts_nb(
    index: Dict,
    read: str,
    read_len: int,
    seed_offsets: Iterable[int],
    seed_len: int,
    ref_len: int,
) -> NumbaList:
    starts = NumbaList.empty_list(types.int64)
    for k in range(len(seed_offsets)):
        offset = seed_offsets[k]
        seed = read[offset : offset + seed_len]
        if seed in index:
            hits = index[seed]
            for j in range(len(hits)):
                start = hits[j] - offset
                if 0 <= start and start + read_len <= ref_len:
                    starts.append(start)
    return starts


def build_seed_index(ref: str, seed_len: int) -> dict:
    """Build seed index using numba-accelerated loop."""
    return _build_seed_index_nb(ref, seed_len)


def vecmap_numba(
    ref: str,
    reads: Iterable[Tuple[str, str]],
    read_len: int,
    seed_len: int = 20,
    seed_offsets: Tuple[int, ...] = (0, 20, 40, 60, 80),
) -> List[Tuple[int, int, str]]:
    """Vectorized mapping with numba-accelerated candidate search."""
    index = build_seed_index(ref, seed_len)
    ref_arr = np.array(list(ref))
    mappings: List[Tuple[int, int, str]] = []
    for read, read_id in reads:
        cand_starts_nb = _get_candidate_starts_nb(
            index,
            read,
            read_len,
            np.array(seed_offsets, dtype=np.int64),
            seed_len,
            len(ref),
        )
        candidate_list = sorted(set(cand_starts_nb))
        if not candidate_list:
            mappings.append((-1, -1, read_id))
            continue
        read_arr = np.array(list(read))
        starts = np.array(candidate_list)
        substrs = ref_arr[starts[:, None] + np.arange(read_len)]
        mismatches_arr = (substrs != read_arr).sum(axis=1)
        best_idx = mismatches_arr.argmin()
        mappings.append(
            (
                int(starts[best_idx]),
                int(mismatches_arr[best_idx]),
                read_id,
            )
        )
    return mappings
