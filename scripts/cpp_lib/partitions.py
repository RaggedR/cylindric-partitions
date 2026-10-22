"""
partitions.py — Partition utilities for cylindric partition computations.

Pure Python, no Sage dependency. These are hot-path functions used in
the innermost enumeration loops, so they stay as plain Python for speed.
"""


def get_part(nu, j):
    """0-indexed j-th part of partition nu; 0 if j >= len(nu)."""
    return nu[j] if j < len(nu) else 0


def gen_partitions(max_part, max_size):
    """Generate all partitions with largest part <= max_part, |lambda| <= max_size."""
    def _gen(prev_val, remaining):
        if remaining <= 0 or prev_val <= 0:
            yield ()
            return
        yield ()
        for p in range(min(prev_val, remaining), 0, -1):
            for rest in _gen(p, remaining - p):
                yield (p,) + rest
    return _gen(max_part, max_size)


def gen_compatible(nu_prev, shift, max_part, max_size):
    """
    Generate partitions nu with nu_prev[j] >= nu[j+shift] for all j >= 0.

    Equivalently: nu[k] <= nu_prev[k-shift] for k >= shift, and
    nu[k] <= max_part for k < shift (partition condition only).
    Also: parts <= max_part, total size <= max_size.
    """
    max_len = len(nu_prev) + shift
    if max_len == 0:
        max_len = min(max_part, max_size) if max_part > 0 and max_size > 0 else 0
        if shift == 0:
            yield ()
            return
        else:
            max_len = shift

    def _gen(pos, prev_val, remaining):
        if remaining <= 0 or prev_val <= 0:
            yield ()
            return
        if pos >= max_len:
            yield ()
            return
        yield ()
        if pos >= shift:
            ub = min(prev_val, get_part(nu_prev, pos - shift), remaining)
        else:
            ub = min(prev_val, max_part, remaining)
        for p in range(ub, 0, -1):
            for rest in _gen(pos + 1, p, remaining - p):
                yield (p,) + rest

    init_ub = min(max_part, max_size)
    if shift > 0:
        init_ub = min(init_ub, max_part)
    else:
        init_ub = min(init_ub, get_part(nu_prev, 0) if nu_prev else 0)
        if init_ub <= 0:
            yield ()
            return

    for parts in _gen(0, init_ub, max_size):
        yield parts


def conjugate(lam):
    """Conjugate of a partition (transpose the Young diagram).

    Input: partition as tuple in decreasing order, e.g. (5, 3, 2).
    Output: conjugate partition as tuple in decreasing order.
    """
    if not lam:
        return ()
    cols = []
    for i in range(lam[0]):
        col_len = sum(1 for part in lam if part > i)
        cols.append(col_len)
    return tuple(cols)


def partitions_of(n, max_part=None):
    """Generate all partitions of n as tuples in decreasing order."""
    if max_part is None:
        max_part = n
    if n == 0:
        yield ()
        return
    for first in range(min(n, max_part), 0, -1):
        for rest in partitions_of(n - first, first):
            yield (first,) + rest
