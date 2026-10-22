#!/usr/bin/env python3
"""
fixed_point_structure.py — Study the structure of U-toggle fixed points.

After the U-toggle involution, the remaining fixed points are
(S, ∅, ν) triples that encode the q-binomial alternating sum:
    Q_n = Σ_j (-1)^j q^{j(j+1)/2} [n,j]_q · g_{n-j}(q)

This script analyzes:
1. The j-stratification: fixed points grouped by |S| = j
2. Within each j-layer: which CPPs ν contribute?
3. Can a second involution on S handle the remaining cancellation?
4. What is the "fixed point of fixed point" — the survivors that give Q_n?
"""

from collections import defaultdict
from itertools import combinations
from cpp_lib import enumerate_cpps_with_max, compute_g_m, parse_args, get_part
import sys


def analyze_fixed_points(c, n, max_weight):
    """
    Analyze the structure of U-toggle fixed points (S, ∅, ν).
    """
    d = sum(c)
    print(f"Fixed point structure: c={c}, d={d}, n={n}")
    print(f"max_weight={max_weight}")
    print("=" * 70)

    # Compute g_m for m = 0,...,n
    print("\n--- g_m(q) polynomials ---")
    g = {}
    for m in range(n + 1):
        g[m] = compute_g_m(c, m, max_weight)
        terms = []
        for w in sorted(g[m].keys())[:10]:
            terms.append(f"{g[m][w]}q^{w}" if w > 0 else str(g[m][w]))
        if len(g[m]) > 10:
            terms.append("...")
        print(f"  g_{m}(q) = {' + '.join(terms)}")
        print(f"    ({len(g[m])} terms, g_{m}(1) = {sum(g[m].values())})")

    # Enumerate fixed points: (S, ∅, ν) for each j = |S|
    print(f"\n--- Fixed points by j = |S| ---")

    all_fixed = []
    for j in range(n + 1):
        m = n - j
        sign = (-1) ** j

        for S in combinations(range(1, n + 1), j):
            S = frozenset(S)
            sig_S = sum(S)

            if m == 0:
                cpps = [((), (), (), 0)]
            else:
                cpps = enumerate_cpps_with_max(c, m, max_weight)

            for nu1, nu2, nu3, cpp_wt in cpps:
                total_wt = sig_S + cpp_wt
                all_fixed.append({
                    'S': S,
                    'nu': (nu1, nu2, nu3),
                    'sign': sign,
                    'q_weight': total_wt,
                    'j': j,
                    'm': m,
                    'cpp_weight': cpp_wt,
                })

    # Group by j
    by_j = defaultdict(list)
    for f in all_fixed:
        by_j[f['j']].append(f)

    total_pos = sum(1 for f in all_fixed if f['sign'] > 0)
    total_neg = sum(1 for f in all_fixed if f['sign'] < 0)
    print(f"Total fixed points: {len(all_fixed)} "
          f"({total_pos} pos, {total_neg} neg)")

    for j in range(n + 1):
        items = by_j[j]
        m = n - j
        sign = (-1) ** j
        sign_str = "+" if sign > 0 else "−"
        # Count distinct S choices
        n_S = len(set(f['S'] for f in items))
        n_cpps = len(items) // n_S if n_S > 0 else 0
        print(f"\n  j={j}: m={m}, sign={sign_str}, "
              f"{n_S} choices of S × {n_cpps} CPPs = {len(items)} triples")

        # Show S contributions: Σ_{S,|S|=j} q^{σ(S)}
        S_contrib = defaultdict(int)
        for S in combinations(range(1, n + 1), j):
            S_contrib[sum(S)] += 1
        S_terms = []
        for w in sorted(S_contrib.keys()):
            S_terms.append(f"{S_contrib[w]}q^{w}" if w > 0 else str(S_contrib[w]))
        print(f"    Σ q^{{σ(S)}} = {' + '.join(S_terms)}")
        print(f"    = q^{{{j*(j+1)//2}}} [{n},{j}]_q")

    # Compute Q_n from fixed points
    Q_from_fixed = defaultdict(int)
    for f in all_fixed:
        Q_from_fixed[f['q_weight']] += f['sign']

    # Show Q_n
    print(f"\n--- Q_{n} from fixed points ---")
    terms = []
    for w in sorted(Q_from_fixed.keys()):
        c_w = Q_from_fixed[w]
        if c_w != 0:
            terms.append(f"{c_w}q^{w}" if w > 0 else str(c_w))
    print(f"Q_{n} = {' + '.join(terms[:20])}")

    # Analyze the j-layered cancellation
    print(f"\n--- Layer-by-layer cancellation ---")
    safe_max_w = n * n * d // 3 + n * d + 5  # generous bound

    for w in range(safe_max_w + 1):
        layers = []
        total = 0
        for j in range(n + 1):
            m = n - j
            sign = (-1) ** j
            layer_contrib = 0
            for S in combinations(range(1, n + 1), j):
                s = sum(S)
                if w - s >= 0 and (w - s) in g[m]:
                    layer_contrib += g[m][w - s]
            layer_contrib *= sign
            total += layer_contrib
            if layer_contrib != 0:
                layers.append(f"j={j}:{layer_contrib:+d}")

        if total != 0 or layers:
            if w <= safe_max_w:
                print(f"  q^{w}: {' '.join(layers)} → {total}")

    # Show the "survivor" structure
    print(f"\n--- Survivor analysis ---")
    print(f"The survivors (net positive) at each weight are CPPs from j=0 layer")
    print(f"that survive cancellation from higher j layers.")
    print()

    for w in range(safe_max_w + 1):
        q_coeff = Q_from_fixed.get(w, 0)
        if q_coeff <= 0:
            continue

        j0_count = g[n].get(w, 0)
        print(f"  q^{w}: Q_n coeff = {q_coeff}")
        print(f"    j=0 layer: {j0_count} CPPs with max={n}, weight={w}")

        for j in range(1, n + 1):
            m = n - j
            sign = (-1) ** j
            cancel = 0
            for S in combinations(range(1, n + 1), j):
                s = sum(S)
                if w - s >= 0 and (w - s) in g[m]:
                    cancel += g[m][w - s]
            if cancel > 0:
                print(f"    j={j} layer: -{cancel} (sign {'+' if sign>0 else '-'})")


def main():
    n, c, max_weight = parse_args(sys.argv,
                                  defaults={'n': 2, 'c': (2, 1, 1), 'max_weight': 30})
    analyze_fixed_points(c, n, max_weight)


if __name__ == '__main__':
    main()
