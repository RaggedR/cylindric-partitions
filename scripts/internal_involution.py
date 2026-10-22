#!/usr/bin/env python3
"""
internal_involution.py — Construct an explicit sign-reversing involution
on (S, U, ν) triples that preserves the QSym monomial (frequency vector).

Since 99.9% of cancellation is internal to QSym classes, finding this
involution solves almost all of the cancellation problem. The remaining
~0.1% requires cross-QSym cancellation at low q-weights.

STRATEGY: Within a fixed frequency vector α, the triples differ only
in their (S, U) decomposition and the specific CPP ν. The natural
involution pairs:
    (S, U, ν)  ↔  (S, U △ {k}, ν')
where △ denotes symmetric difference and ν' is obtained by
adding/removing parts from ν to compensate for the weight change.
"""

from collections import defaultdict
from cpp_lib import (enumerate_cpps_with_max, frequency_vector,
                     enumerate_all_triples, parse_args)
import sys


def try_U_toggle_involution(triples, c, n):
    """
    Try to construct involution by toggling elements of U.

    The idea: for a triple (S, U, ν) with max_entry m and k ∈ {1,...,m},
    toggle k in U (add if absent, remove if present). This flips the sign.
    To preserve q-weight: toggling k in U changes weight by ±k.
    We need to compensate by changing ν by ∓k. Specifically:

    If k ∉ U: add k to U (weight +k from σ(U)), need to decrease |ν| by k.
              Remove k parts equal to 1 from ν? No, that changes freq.
              We need to change ν without changing the frequency vector!

    The issue: to keep freq fixed, the CPP parts can't change in value.
    Only the number of parts at each value can change. But adding k to U
    already adds one count at value k, so ν must lose one count at value k.

    So: toggle k in U ↔ add/remove a part of value k from ν.
    - (S, U, ν) with k ∉ U, ν has at least one part = k
      ↔ (S, U∪{k}, ν' where ν' has one fewer part = k)
    - Both have the same freq: U gains +1 at k, ν loses -1 at k → net 0.
    - Weight: σ(U) +k, |ν| -k → net 0. ✓
    - Sign: (-1)^{|U|+1} = -(-1)^{|U|}. ✓ (sign flips)

    But we need ν' to be a valid CPP with max = m!
    """
    print(f"\n{'='*70}")
    print(f"U-TOGGLE INVOLUTION")
    print(f"{'='*70}")

    # Group by freq vector
    by_freq = defaultdict(list)
    for i, t in enumerate(triples):
        by_freq[t['freq']].append(i)

    paired = [None] * len(triples)
    n_paired = 0
    n_fixed = 0

    for freq, indices in by_freq.items():
        # For each freq class, try to pair triples
        for i in indices:
            if paired[i] is not None:
                continue

            t = triples[i]
            m = t['max_entry']
            nu1, nu2, nu3 = t['nu']

            found_partner = False

            # Try toggling each k ∈ {1,...,m} in U
            for k in range(1, m + 1):
                if k in t['U']:
                    # k ∈ U: remove k from U, add a part k to ν
                    new_U = t['U'] - {k}
                    target_cpp_wt = sum(nu1) + sum(nu2) + sum(nu3) + k

                    for j in indices:
                        if j == i or paired[j] is not None:
                            continue
                        t2 = triples[j]
                        if t2['S'] != t['S']:
                            continue
                        if t2['U'] != new_U:
                            continue
                        if t2['max_entry'] != m:
                            continue
                        cpp_wt2 = sum(t2['nu'][0]) + sum(t2['nu'][1]) + sum(t2['nu'][2])
                        if cpp_wt2 != target_cpp_wt:
                            continue
                        parts_k_orig = sum(1 for p in nu1+nu2+nu3 if p == k)
                        parts_k_new = sum(1 for p in t2['nu'][0]+t2['nu'][1]+t2['nu'][2] if p == k)
                        if parts_k_new == parts_k_orig + 1:
                            if t['q_weight'] == t2['q_weight']:
                                paired[i] = j
                                paired[j] = i
                                n_paired += 1
                                found_partner = True
                                break

                    if found_partner:
                        break

                else:
                    # k ∉ U: add k to U, remove a part k from ν
                    parts_k = sum(1 for p in nu1+nu2+nu3 if p == k)
                    if parts_k == 0:
                        continue

                    new_U = t['U'] | {k}
                    target_cpp_wt = sum(nu1) + sum(nu2) + sum(nu3) - k

                    for j in indices:
                        if j == i or paired[j] is not None:
                            continue
                        t2 = triples[j]
                        if t2['S'] != t['S']:
                            continue
                        if t2['U'] != new_U:
                            continue
                        if t2['max_entry'] != m:
                            continue
                        cpp_wt2 = sum(t2['nu'][0]) + sum(t2['nu'][1]) + sum(t2['nu'][2])
                        if cpp_wt2 != target_cpp_wt:
                            continue
                        parts_k_new = sum(1 for p in t2['nu'][0]+t2['nu'][1]+t2['nu'][2] if p == k)
                        if parts_k_new == parts_k - 1:
                            if t['q_weight'] == t2['q_weight']:
                                paired[i] = j
                                paired[j] = i
                                n_paired += 1
                                found_partner = True
                                break

                    if found_partner:
                        break

            if not found_partner:
                n_fixed += 1

    # Verify
    errors = 0
    for i in range(len(triples)):
        j = paired[i]
        if j is not None and j > i:
            ti, tj = triples[i], triples[j]
            if ti['q_weight'] != tj['q_weight']:
                errors += 1
            if ti['sign'] == tj['sign']:
                errors += 1
            if ti['freq'] != tj['freq']:
                errors += 1

    print(f"Paired: {n_paired} pairs ({2*n_paired} triples)")
    print(f"Fixed points: {n_fixed} triples")
    if errors == 0:
        print(f"All pairs verified: same weight, opposite sign, same freq ✓")
    else:
        print(f"ERRORS: {errors}")

    # Analyze fixed points
    fixed_by_weight = defaultdict(int)
    fixed_items = []
    for i in range(len(triples)):
        if paired[i] is None:
            t = triples[i]
            fixed_by_weight[t['q_weight']] += t['sign']
            fixed_items.append(t)

    # Fixed points should give Q_n
    print(f"\nFixed point contribution (should give Q_n):")
    terms = []
    for w in sorted(fixed_by_weight.keys()):
        c_w = fixed_by_weight[w]
        if c_w != 0:
            terms.append(f"{c_w}q^{w}" if w > 0 else str(c_w))
    print(f"  {' + '.join(terms[:20]) if terms else '0'}")

    # Compare with actual Q_n
    actual_Q = defaultdict(int)
    for t in triples:
        actual_Q[t['q_weight']] += t['sign']

    ok = True
    for w in set(list(fixed_by_weight.keys()) + list(actual_Q.keys())):
        if fixed_by_weight.get(w, 0) != actual_Q.get(w, 0):
            print(f"  MISMATCH at q^{w}: fixed={fixed_by_weight.get(w,0)}, "
                  f"actual={actual_Q.get(w,0)}")
            ok = False
    if ok:
        print(f"  Matches Q_n ✓")

    # Show fixed points at low weights
    print(f"\nFixed point details (first 30):")
    for t in sorted(fixed_items, key=lambda t: t['q_weight'])[:30]:
        sign_str = "+" if t['sign'] == +1 else "−"
        S_str = "{" + ",".join(str(s) for s in sorted(t['S'])) + "}" \
                if t['S'] else "∅"
        U_str = "{" + ",".join(str(u) for u in sorted(t['U'])) + "}" \
                if t['U'] else "∅"
        print(f"  {sign_str} q^{t['q_weight']} z^{t['freq']} "
              f"S={S_str} U={U_str} ν={t['nu']}")

    return paired, n_paired, n_fixed


def main():
    n, c, max_weight = parse_args(sys.argv,
                                  defaults={'n': 1, 'c': (2, 1, 1), 'max_weight': 30})

    d = sum(c)
    print(f"Internal involution search: c={c}, d={d}, n={n}")
    print(f"max_weight={max_weight}")

    triples = enumerate_all_triples(c, n, max_weight)
    print(f"Total triples: {len(triples)}")

    paired, n_paired, n_fixed = try_U_toggle_involution(triples, c, n)


if __name__ == '__main__':
    main()
