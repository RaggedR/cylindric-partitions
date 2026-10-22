#!/usr/bin/env python3
"""
qsym_Q_polynomial.py — Compute the "QSym Q-polynomial" Q_n^{QSym}
and test L-positivity.

Q_n^{QSym} = Σ_{(S,U,ν)} (-1)^{|S|+|U|} M_{freq(S,U,ν)}

where the sum is over all triples in the ground set for Q_{n,c}.

KEY QUESTION: Is Q_n^{QSym} L-positive (nonneg in Gessel's fundamental basis)?
If yes, then L-positivity of Q_n^{QSym} implies Q_n(q) ≥ 0.

RESULT: NO — Q_n^{QSym} is NOT L-positive (nor even M-positive for n >= 2).
This is a different object from GK_tilde (which IS L-positive).
"""

from collections import defaultdict
from itertools import combinations
from cpp_lib import (enumerate_cpps_with_max, frequency_vector,
                     compositions_of, refines, coarsens,
                     q_specialization, M_to_L, parse_args)
import sys


def compute_qsym_Q(c, n, max_weight=60):
    """
    Compute Q_n^{QSym} as a dict: composition -> coefficient.
    Also compute Q_n(q) for comparison.
    """
    M_coeffs = defaultdict(int)  # composition -> coefficient
    q_coeffs = defaultdict(int)  # q-weight -> coefficient

    for j in range(n + 1):
        m = n - j
        for S in combinations(range(1, n + 1), j):
            S = frozenset(S)
            sig_S = sum(S)
            if m == 0:
                cpps = [((), (), (), 0)]
            else:
                cpps = enumerate_cpps_with_max(c, m, max_weight)
            for nu1, nu2, nu3, cpp_wt in cpps:
                for u_size in range(m + 1):
                    for U in combinations(range(1, m + 1), u_size):
                        U = frozenset(U)
                        total_wt = sig_S + sum(U) + cpp_wt
                        sign = (-1) ** (len(S) + len(U))
                        freq = frequency_vector(S, U, nu1, nu2, nu3)
                        M_coeffs[freq] += sign
                        q_coeffs[total_wt] += sign

    # Remove zero entries
    M_coeffs = {k: v for k, v in M_coeffs.items() if v != 0}
    q_coeffs = {k: v for k, v in q_coeffs.items() if v != 0}

    return M_coeffs, q_coeffs


def main():
    n, c, max_weight = parse_args(sys.argv,
                                  defaults={'n': 1, 'c': (2, 1, 1), 'max_weight': 30})

    d = sum(c)
    print(f"Computing Q_{n}^{{QSym}} for c={c}, d={d}")
    print(f"max_weight={max_weight}")
    print("=" * 70)

    M_coeffs, q_coeffs = compute_qsym_Q(c, n, max_weight)

    # Show Q_n(q)
    max_w = max(q_coeffs.keys()) if q_coeffs else 0
    terms = []
    for w in range(max_w + 1):
        c_w = q_coeffs.get(w, 0)
        if c_w != 0:
            terms.append(f"{c_w}q^{w}" if w > 0 else str(c_w))
    print(f"\nQ_{n}(q) = {' + '.join(terms) if terms else '0'}")

    # Determine safe range: Q_n should be a polynomial of degree ~ n^2 * d/3
    expected_deg = n * n * d // 3 + n * d  # generous upper bound
    safe_max_w = min(max_weight - n, expected_deg + 5)
    print(f"Expected max degree of Q_{n}: ~{n*n*d//3}")
    print(f"Safe range for analysis: q^0 to q^{safe_max_w}")

    # Filter M-coefficients to safe range
    safe_M = {}
    for alpha, coeff in M_coeffs.items():
        qw = q_specialization(alpha)
        if qw <= safe_max_w:
            safe_M[alpha] = coeff

    print(f"\n{'='*70}")
    print(f"M-BASIS EXPANSION (safe range)")
    print(f"{'='*70}")
    print(f"Q_{n}^{{QSym}} = Σ a_α M_α  where:")

    # Group by q-weight for display
    by_qw = defaultdict(list)
    for alpha, coeff in safe_M.items():
        qw = q_specialization(alpha)
        by_qw[qw].append((alpha, coeff))

    n_neg_M = 0
    n_pos_M = 0
    for qw in sorted(by_qw.keys()):
        items = by_qw[qw]
        items.sort(key=lambda x: x[0])
        for alpha, coeff in items:
            sign = "+" if coeff > 0 else ""
            print(f"  a_{alpha} = {sign}{coeff}")
            if coeff > 0:
                n_pos_M += 1
            elif coeff < 0:
                n_neg_M += 1

    print(f"\nM-basis: {n_pos_M} positive, {n_neg_M} negative")
    if n_neg_M == 0:
        print(">>> Q_n^{QSym} IS M-POSITIVE! <<<")
    else:
        print("Q_n^{QSym} is NOT M-positive")

    # Compute L-basis expansion
    print(f"\n{'='*70}")
    print(f"L-BASIS EXPANSION (safe range)")
    print(f"{'='*70}")

    L_coeffs = M_to_L(safe_M)
    print(f"Q_{n}^{{QSym}} = Σ b_α L_α  where:")

    # Group by q-weight
    by_qw_L = defaultdict(list)
    for alpha, coeff in L_coeffs.items():
        qw = q_specialization(alpha)
        by_qw_L[qw].append((alpha, coeff))

    n_neg_L = 0
    n_pos_L = 0
    for qw in sorted(by_qw_L.keys()):
        items = by_qw_L[qw]
        items.sort(key=lambda x: x[0])
        for alpha, coeff in items:
            sign = "+" if coeff > 0 else ""
            print(f"  b_{alpha} = {sign}{coeff}")
            if coeff > 0:
                n_pos_L += 1
            elif coeff < 0:
                n_neg_L += 1

    print(f"\nL-basis: {n_pos_L} positive, {n_neg_L} negative")
    if n_neg_L == 0:
        print(">>> Q_n^{QSym} IS L-POSITIVE! <<<")
    else:
        print("Q_n^{QSym} is NOT L-positive")
        # Show negatives
        print("\nNegative L-coefficients:")
        for alpha, coeff in sorted(L_coeffs.items()):
            if coeff < 0:
                print(f"  b_{alpha} = {coeff}")

    # Verify q-specialization
    print(f"\n{'='*70}")
    print(f"VERIFICATION: q-specialization")
    print(f"{'='*70}")

    q_from_M = defaultdict(int)
    for alpha, coeff in safe_M.items():
        qw = q_specialization(alpha)
        q_from_M[qw] += coeff

    ok = True
    for w in range(safe_max_w + 1):
        from_M = q_from_M.get(w, 0)
        from_Q = q_coeffs.get(w, 0)
        if from_M != from_Q:
            print(f"  MISMATCH at q^{w}: M gives {from_M}, Q gives {from_Q}")
            ok = False
    if ok:
        print("  All q-specialization checks pass ✓")


if __name__ == '__main__':
    main()
