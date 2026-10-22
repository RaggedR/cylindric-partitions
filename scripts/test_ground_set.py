#!/usr/bin/env python3
"""
test_ground_set.py — Enumerate the ground set G_{n,c} explicitly.

For the REGULAR case (ordinary partitions, no cylindric constraint):
  h_m = q^m, so H_m has one element of weight m.
  The ground set reduces to subsets S ⊆ {1,...,n}, with weight σ(S) + (n-|S|).
  Q_n = 0 (complete cancellation).

For the CYLINDRIC case:
  h_{m,c}(q) is a nonneg polynomial (computed from GK_c).
  The ground set has richer structure.
  Q_{n,c}(q) ≥ 0 (the conjecture).

We enumerate elements, group by weight, sum signs, and verify the result
matches Q_{n,c}(q).
"""

from collections import defaultdict
from itertools import combinations


# ---------------------------------------------------------------------------
# Regular partitions: h_m = q^m
# ---------------------------------------------------------------------------

def ground_set_regular(n):
    """
    Enumerate the ground set for regular partitions.
    h_m = q^m, so H_m = {single element of weight m}.

    Elements: (S, m) where S ⊆ {1,...,n}, m = n - |S|.
    Sign: (-1)^|S|
    Weight: σ(S) + m = σ(S) + n - |S|
    """
    elements = []
    for j in range(n + 1):
        for S in combinations(range(1, n + 1), j):
            S = frozenset(S)
            m = n - j
            sigma_S = sum(S)
            weight = sigma_S + m
            sign = (-1) ** j
            elements.append((S, m, sign, weight))
    return elements


def test_regular(n):
    """Test complete cancellation for regular partitions."""
    print(f"\n{'='*60}")
    print(f"REGULAR PARTITIONS: n = {n}")
    print(f"{'='*60}")

    elts = ground_set_regular(n)
    print(f"Ground set size: {len(elts)} = 2^{n}")

    # Group by weight, sum signs
    by_weight = defaultdict(int)
    by_weight_detail = defaultdict(list)
    for S, m, sign, w in elts:
        by_weight[w] += sign
        by_weight_detail[w].append((S, m, sign))

    print(f"\nSigned sum by weight (should all be 0 for n ≥ 1):")
    max_w = max(by_weight.keys()) if by_weight else 0
    all_zero = True
    for w in range(max_w + 1):
        val = by_weight.get(w, 0)
        if val != 0:
            all_zero = False
        if n <= 4 or val != 0:
            print(f"  q^{w}: {by_weight.get(w, 0):+d}")

    if all_zero and n >= 1:
        print(f"\n  ✓ Complete cancellation (Q_{n} = 0)")
    elif n == 0:
        print(f"\n  Q_0 = 1 (as expected)")
    else:
        print(f"\n  ✗ INCOMPLETE cancellation!")

    # Show the pairing for small n
    if n <= 3:
        print(f"\nGround set elements:")
        for S, m, sign, w in sorted(elts, key=lambda x: (x[3], -x[2])):
            s_str = '{' + ','.join(str(x) for x in sorted(S)) + '}'
            print(f"  S={s_str:<12} τ_wt={m}  sign={sign:+d}  "
                  f"total_wt={w}")

    # Show the "flip element 1" involution
    if n <= 4:
        print(f"\nFlip-1 involution pairing:")
        for S, m, sign, w in sorted(elts, key=lambda x: (x[3], len(x[0]))):
            if 1 not in S:
                S2 = frozenset(S | {1})
                m2 = n - len(S2)
                w2 = sum(S2) + m2
                s1 = '{' + ','.join(str(x) for x in sorted(S)) + '}'
                s2 = '{' + ','.join(str(x) for x in sorted(S2)) + '}'
                print(f"  ({s1:<8}, wt={w}) ↔ ({s2:<8}, wt={w2})"
                      f"  signs {sign:+d} and {(-1)**len(S2):+d}")


# ---------------------------------------------------------------------------
# Cylindric partitions: use compute_gk
# ---------------------------------------------------------------------------

def compute_h_m_all(c, max_m, max_w):
    """
    Compute h_{m,c}(q) = (q)_m · f_{m,c}(q) for m = 0, ..., max_m.

    Returns dict {m: list of int coefficients [a_0, ..., a_{max_w}]}.
    """
    from cpp_lib import compute_gk, extract_z_coefficient, compute_h_m

    gk = compute_gk(c, max_n=max_m, max_w=max_w)

    result = {}
    for m in range(max_m + 1):
        g_m = extract_z_coefficient(gk, m)
        h = compute_h_m(c, m, max_w, g_m_dict=g_m)
        result[m] = h

    return result


def ground_set_cylindric(c, n, max_w):
    """
    Enumerate the ground set G_{n,c} for cylindric partitions.

    Elements: (S, τ) where S ⊆ {1,...,n}, τ ∈ H_{n-|S|,c}.
    H_m has elements indexed by weight w with multiplicity h_{m,c}[w].
    Sign: (-1)^|S|
    Weight: σ(S) + wt(τ)

    Returns elements as (S, tau_weight, multiplicity, sign, total_weight).
    """
    h = compute_h_m_all(c, max_m=n, max_w=max_w)

    elements = []
    for j in range(n + 1):
        m = n - j
        h_m = h[m]
        for S in combinations(range(1, n + 1), j):
            S = frozenset(S)
            sigma_S = sum(S)
            sign = (-1) ** j
            for tau_w in range(max_w + 1):
                mult = int(h_m[tau_w])
                if mult == 0:
                    continue
                total_w = sigma_S + tau_w
                if total_w <= max_w:
                    elements.append((S, tau_w, mult, sign, total_w))

    return elements, h


def test_cylindric(c, n, max_w):
    """Test the ground set for cylindric partitions."""
    d = sum(c)
    print(f"\n{'='*60}")
    print(f"CYLINDRIC: c = {c}, d = {d}, n = {n}")
    print(f"{'='*60}")

    elts, h = ground_set_cylindric(c, n, max_w)

    # Show h_m polynomials
    print(f"\nh_m polynomials:")
    for m in range(n + 1):
        h_m = h[m]
        terms = [(w, int(h_m[w])) for w in range(max_w + 1) if h_m[w] != 0]
        poly_str = " + ".join(f"{c}q^{w}" if c != 1 else f"q^{w}"
                              for w, c in terms) if terms else "0"
        print(f"  h_{m} = {poly_str}")

    # Group by total weight, sum (sign × multiplicity)
    by_weight = defaultdict(int)
    for S, tau_w, mult, sign, total_w in elts:
        by_weight[total_w] += sign * mult

    # Also compute Q_{n,c} via cpp_lib for comparison
    from cpp_lib import compute_gk as _compute_gk, compute_Q
    gk = _compute_gk(c, max_n=n + max_w, max_w=max_w + n*(n+1))
    Q = compute_Q(c, n, max_w, gk_dict=gk)

    print(f"\nSigned sum vs Q_{n},{c}(q):")
    max_wt = max(max(by_weight.keys(), default=0),
                 max(w for w in range(len(Q)) if Q[w] != 0) if any(Q[w] != 0 for w in range(len(Q))) else 0)

    all_match = True
    for w in range(min(max_wt + 1, max_w + 1)):
        ground = by_weight.get(w, 0)
        q_val = int(Q[w]) if w < len(Q) else 0
        marker = "  ✓" if ground == q_val else "  ✗ MISMATCH"
        if ground != 0 or q_val != 0:
            print(f"  q^{w}: ground_set={ground:+d}  Q={q_val:+d}{marker}")
        if ground != q_val:
            all_match = False

    if all_match:
        print(f"\n  ✓ Ground set signed sum matches Q_{n},{c}(q)")
    else:
        print(f"\n  ✗ MISMATCH detected!")

    # Show ground set elements for small cases
    total_pos = sum(mult for _, _, mult, sign, _ in elts if sign > 0)
    total_neg = sum(mult for _, _, mult, sign, _ in elts if sign < 0)
    print(f"\nGround set: {total_pos} positive, {total_neg} negative "
          f"(total {total_pos + total_neg})")

    # Show elements for small n
    if n <= 2 and d <= 5:
        print(f"\nGround set elements (S, τ_wt × mult, sign, total_wt):")
        for S, tau_w, mult, sign, total_w in sorted(
                elts, key=lambda x: (x[4], -x[3])):
            s_str = '{' + ','.join(str(x) for x in sorted(S)) + '}'
            print(f"  S={s_str:<8} τ_wt={tau_w} ×{mult}  "
                  f"sign={sign:+d}  total_wt={total_w}")


if __name__ == '__main__':
    # Regular partitions
    for n in range(5):
        test_regular(n)

    # Cylindric: d=4, c=(2,2,0), n=1 (simplest non-trivial case)
    test_cylindric((2, 2, 0), n=1, max_w=20)

    # Cylindric: d=4, c=(2,2,0), n=2
    test_cylindric((2, 2, 0), n=2, max_w=30)

    # Cylindric: d=1, c=(1,0,0), n=1
    test_cylindric((1, 0, 0), n=1, max_w=10)

    # Cylindric: d=2, c=(1,1,0), n=1 (Rogers-Ramanujan)
    test_cylindric((1, 1, 0), n=1, max_w=10)
