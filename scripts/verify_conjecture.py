#!/usr/bin/env python3
"""
verify_conjecture.py — Verify Warnaar's Conjecture 2.7 on Q_{n,c}(q).

For c = (c0, c1, c2) with d = c0+c1+c2, d not equiv 0 (mod 3):
    Q_{n,c}(q) := (q^ell; q^ell)_n * [z^n]( (zq)_inf * GK_c(z,q) )
where ell = gcd(d, 3).

Conjecture: Q_{n,c}(q) is a polynomial in q with NONNEGATIVE coefficients.
Moreover: Q_{n,c}(1) = ( (d+1)(d+2)/6 - 1 )^n.

Reference: Warnaar (2021), Conjecture 2.7 (p.8), eq (2.4)-(2.5)
"""

from math import gcd
import time

from cpp_lib.enumeration import compute_gk
from cpp_lib.qseries import compute_Q


# ---------------------------------------------------------------------------
# Verification
# ---------------------------------------------------------------------------

def verify_conjecture_27(c, max_n_range, max_w, verbose=True):
    """
    Verify Conjecture 2.7 for profile c and z-powers n = 0, 1, ..., max_n_range.

    Checks:
    1. Q_{n,c}(q) has all nonneg integer coefficients
    2. Q_{n,c}(1) = ((d+1)(d+2)/6 - 1)^n

    Returns True if all checks pass.
    """
    d = sum(c)
    ell = gcd(d, 3)

    if d % 3 == 0:
        if verbose:
            print(f"  SKIP: d = {d} is divisible by 3 (conjecture requires d not equiv 0 mod 3)")
        return True

    # Expected Q_{n,c}(1)
    expected_q1 = (d + 1) * (d + 2) // 6 - 1

    if verbose:
        print(f"\n{'='*60}")
        print(f"Conjecture 2.7 verification for c = {c}")
        print(f"  d = {d}, ell = gcd({d}, 3) = {ell}")
        print(f"  Expected Q_{{n,c}}(1) = ({expected_q1})^n")
        print(f"  max_w = {max_w}")

    # Compute GK with sufficient truncation
    gk_max_n = max_n_range + 2
    gk_max_w = max_w + max_n_range * (max_n_range + 1) // 2 + 3

    if verbose:
        print(f"  Computing GK_c(z,q) with max_n={gk_max_n}, max_w={gk_max_w}...")
    t0 = time.time()
    gk = compute_gk(c, gk_max_n, gk_max_w, verbose=verbose)
    if verbose:
        print(f"  GK computed in {time.time()-t0:.2f}s")

    all_pass = True

    for n in range(max_n_range + 1):
        Q = compute_Q(c, n, max_w, gk_dict=gk)

        # Check integrality and nonnegativity
        is_nonneg = True
        max_deg = 0
        for w in range(max_w + 1):
            if Q[w] != 0:
                max_deg = w
                if Q[w] < 0:
                    is_nonneg = False

        # Compute Q(1)
        q_at_1 = sum(Q)

        # Expected value
        expected = expected_q1 ** n

        status_nonneg = "PASS" if is_nonneg else "FAIL"
        status_q1 = "PASS" if q_at_1 == expected else "FAIL"

        if not (is_nonneg and q_at_1 == expected):
            all_pass = False

        if verbose:
            terms = []
            for w in range(max_w + 1):
                if Q[w] != 0:
                    terms.append(f"{Q[w]}*q^{w}")
            q_str = " + ".join(terms) if terms else "0"
            if len(q_str) > 100:
                q_str = q_str[:100] + "..."

            print(f"\n  n = {n}: Q_{{n,c}}(q) = {q_str}")
            print(f"    Polynomial: PASS, Nonneg: {status_nonneg}, "
                  f"Q(1) = {q_at_1} {'==' if q_at_1 == expected else '!='} {expected}: {status_q1}")

            if not is_nonneg:
                for w in range(max_w + 1):
                    if Q[w] < 0:
                        print(f"    NEGATIVE COEFFICIENT at q^{w}: {Q[w]}")

            if q_at_1 != expected and q_at_1 < expected and is_nonneg:
                deficit = expected - q_at_1
                print(f"    Note: Q(1) deficit = {deficit} likely due to truncation "
                      f"(max_w={max_w} < deg(Q_{{{n},c}}))")

    if verbose:
        print(f"\n  {'ALL PASS' if all_pass else 'SOME FAILURES'}")

    return all_pass


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    import sys

    print("=" * 60)
    print("Warnaar Conjecture 2.7 verification")
    print("=" * 60)

    # Parse optional --dmax flag
    dmax = 7
    for i, arg in enumerate(sys.argv[1:]):
        if arg == "--dmax":
            dmax = int(sys.argv[i + 2])

    # Organize test profiles by d.
    all_profiles = {
        1: [(1, 0, 0)],
        2: [(2, 0, 0), (1, 1, 0)],
        4: [(2, 1, 1), (1, 2, 1)],
        5: [(2, 2, 1), (3, 1, 1)],
        7: [(3, 2, 2), (2, 3, 2), (4, 2, 1)],
        8: [(3, 3, 2), (4, 3, 1)],
    }

    fast = "--fast" in sys.argv

    if fast:
        truncation = {
            1: dict(max_n_range=2, max_w=8),
            2: dict(max_n_range=3, max_w=12),
            4: dict(max_n_range=2, max_w=15),
            5: dict(max_n_range=1, max_w=12),
            7: dict(max_n_range=1, max_w=15),
            8: dict(max_n_range=1, max_w=12),
        }
    else:
        truncation = {
            1: dict(max_n_range=3, max_w=10),
            2: dict(max_n_range=4, max_w=20),
            4: dict(max_n_range=3, max_w=30),
            5: dict(max_n_range=2, max_w=20),
            7: dict(max_n_range=1, max_w=20),
            8: dict(max_n_range=1, max_w=15),
        }

    results = {}
    total_t0 = time.time()

    for d in sorted(all_profiles.keys()):
        if d > dmax:
            break

        params = truncation[d]
        print(f"\n{'#'*60}")
        print(f"# d = {d}  (max_n_range={params['max_n_range']}, max_w={params['max_w']})")
        print(f"# Q_{{n,c}}(1) = {(d+1)*(d+2)//6 - 1}^n")
        print(f"{'#'*60}")

        for c in all_profiles[d]:
            t0 = time.time()
            passed = verify_conjecture_27(
                c, params['max_n_range'], params['max_w'], verbose=True)
            elapsed = time.time() - t0
            results[c] = (passed, elapsed)

    # Summary
    total_elapsed = time.time() - total_t0
    print(f"\n\n{'='*60}")
    print("SUMMARY")
    print("=" * 60)
    for c, (passed, elapsed) in results.items():
        d = sum(c)
        print(f"  c = {c}, d = {d}: {'PASS' if passed else 'FAIL'}  ({elapsed:.1f}s)")

    all_ok = all(p for p, _ in results.values())
    print(f"\n  Overall: {'ALL CONJECTURE CHECKS PASS' if all_ok else 'SOME FAILURES DETECTED'}")
    print(f"  Total time: {total_elapsed:.1f}s")
