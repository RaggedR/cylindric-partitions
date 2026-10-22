#!/usr/bin/env python3
"""
compute_gk.py — Compute GK_c(z,q) for rank-3 cylindric partitions.

Thin wrapper around cpp_lib providing validation and the __main__ test suite.
All core computation lives in cpp_lib.enumeration and cpp_lib.qseries.

Reference: Warnaar (2021), eqs (3.1), (3.3), (3.4a), (3.11)
           Langer thesis (2021), Chapter 1a
"""

from collections import defaultdict

from cpp_lib.partitions import get_part
from cpp_lib.profiles import profile_to_mu
from cpp_lib.enumeration import (compute_gk, extract_z_coefficient,
                                  gk_z1_coefficients)
from cpp_lib.qseries import borodin_product


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------

def validate_gk(c, max_n, max_w, verbose=True):
    """Compute GK_c(z,q) and validate against Borodin product at z=1."""
    d = sum(c)
    mu = profile_to_mu(c)

    if verbose:
        print(f"\n{'='*60}")
        print(f"Profile c = {c}, d = {d}, mu(c) = {mu}")

    gk = compute_gk(c, max_n, max_w, verbose=verbose)

    # z=1 specialization from enumeration
    enum_z1 = gk_z1_coefficients(gk, max_w)

    # Borodin product
    product = borodin_product(c, max_w + 5)

    safe_w = max_n

    all_match = True
    for w in range(min(safe_w + 1, max_w + 1)):
        if enum_z1[w] != product[w]:
            if verbose:
                print(f"  MISMATCH at q^{w}: enum={enum_z1[w]}, product={product[w]}")
            all_match = False

    if verbose:
        if all_match:
            print(f"  PASS: GK_c(1,q) matches Borodin product through q^{safe_w}")
        else:
            print(f"  FAIL: mismatch found in safe range q^0..q^{safe_w}")

    # Print first few GK terms
    if verbose and len(gk) <= 200:
        print(f"\n  GK_{c}(z,q) terms:")
        for (n, w) in sorted(gk.keys())[:50]:
            if gk[(n, w)] != 0:
                print(f"    z^{n} q^{w} : {gk[(n, w)]}")
        if len(gk) > 50:
            print(f"    ... ({len(gk)} total terms)")

    return gk


def check_cyclic_symmetry(c, max_n, max_w, verbose=True):
    """Check GK_{(c0,c1,c2)} = GK_{(c2,c0,c1)} (Warnaar eq 3.12)."""
    c0, c1, c2 = c
    c_rot = (c2, c0, c1)

    gk1 = compute_gk(c, max_n, max_w)
    gk2 = compute_gk(c_rot, max_n, max_w)

    match = (gk1 == gk2)
    if verbose:
        print(f"  Cyclic: GK_{c} = GK_{c_rot}: {'PASS' if match else 'FAIL'}")
    return match


def check_d1_identity(max_n, max_w, verbose=True):
    """
    Check GK_{(1,0,0)}(z,q) = 1/(zq;q)_inf (Warnaar eq 1.9 / d=1 case).
    """
    if verbose:
        print(f"\n{'='*60}")
        print(f"d=1 identity: GK_{{(1,0,0)}}(z,q) = 1/(zq;q)_inf")

    gk = compute_gk((1, 0, 0), max_n, max_w, verbose=verbose)

    # Compute 1/(zq;q)_inf = prod_{k>=1} 1/(1-zq^k)
    # [z^m q^w] = number of partitions of w into exactly m parts (each >= 1)
    expected = {}
    for m in range(max_n + 1):
        dp = [0] * (max_w + 1)
        dp[0] = 1
        for part in range(1, m + 1):
            for k in range(part, max_w - m + 1):
                dp[k] += dp[k - part]
        for w_reduced in range(max_w - m + 1):
            w = w_reduced + m
            if w <= max_w and dp[w_reduced] > 0:
                expected[(m, w)] = expected.get((m, w), 0) + dp[w_reduced]

    safe_w = max_n
    all_match = True
    for (n, w), cnt in expected.items():
        if n <= max_n and w <= safe_w:
            gk_val = gk.get((n, w), 0)
            if gk_val != cnt:
                if verbose:
                    print(f"  MISMATCH: z^{n}q^{w}: gk={gk_val}, expected={cnt}")
                all_match = False

    for (n, w), cnt in gk.items():
        if n <= max_n and w <= safe_w:
            exp_val = expected.get((n, w), 0)
            if cnt != exp_val:
                if verbose:
                    print(f"  MISMATCH: z^{n}q^{w}: gk={cnt}, expected={exp_val}")
                all_match = False

    if verbose:
        print(f"  {'PASS' if all_match else 'FAIL'}: d=1 identity through q^{safe_w}")

    return all_match


def check_reversal_symmetry(c, max_n, max_w, verbose=True):
    """Check GK_{(c0,c1,c2)} = GK_{(c2,c1,c0)} (reversal symmetry, experimental)."""
    c_rev = tuple(reversed(c))
    if c_rev == c:
        if verbose:
            print(f"  Reversal: c = {c} is palindromic, trivially holds")
        return True

    c0, c1, c2 = c
    orbit = {(c0, c1, c2), (c2, c0, c1), (c1, c2, c0)}
    if c_rev in orbit:
        if verbose:
            print(f"  Reversal: c_rev = {c_rev} is in cyclic orbit of {c}, trivially holds")
        return True

    gk1 = compute_gk(c, max_n, max_w)
    gk2 = compute_gk(c_rev, max_n, max_w)

    safe_w = max_n
    all_match = True
    all_keys = set(gk1.keys()) | set(gk2.keys())
    for (n, w) in all_keys:
        if w <= safe_w and n <= max_n:
            v1 = gk1.get((n, w), 0)
            v2 = gk2.get((n, w), 0)
            if v1 != v2:
                if verbose and all_match:
                    print(f"  MISMATCH at z^{n}q^{w}: GK_{c}={v1}, GK_{c_rev}={v2}")
                all_match = False

    if verbose:
        status = 'PASS' if all_match else 'FAIL (not expected to hold in general)'
        print(f"  Reversal: GK_{c} vs GK_{c_rev}: {status}")

    return all_match


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    print("=" * 60)
    print("GK_c(z,q) computation for rank-3 cylindric partitions")
    print("=" * 60)

    # d=1 tests
    for c in [(1, 0, 0), (0, 1, 0), (0, 0, 1)]:
        validate_gk(c, max_n=10, max_w=10)

    # d=1 closed-form identity
    check_d1_identity(max_n=10, max_w=10)

    # d=2 tests
    for c in [(2, 0, 0), (1, 1, 0), (1, 0, 1), (0, 1, 1), (0, 2, 0), (0, 0, 2)]:
        validate_gk(c, max_n=8, max_w=8)

    # Cyclic symmetry
    print(f"\n{'='*60}")
    print("Cyclic symmetry checks")
    check_cyclic_symmetry((1, 1, 0), max_n=6, max_w=8)
    check_cyclic_symmetry((2, 1, 0), max_n=5, max_w=7)
    check_cyclic_symmetry((2, 1, 1), max_n=4, max_w=6)

    # Reversal symmetry (experimental — not universally valid)
    print(f"\n{'='*60}")
    print("Reversal symmetry checks (experimental)")
    check_reversal_symmetry((1, 1, 0), max_n=6, max_w=8)
    check_reversal_symmetry((2, 1, 0), max_n=5, max_w=7)
    check_reversal_symmetry((2, 1, 1), max_n=4, max_w=6)
    check_reversal_symmetry((2, 2, 1), max_n=3, max_w=5)

    # d=4 test
    validate_gk((2, 1, 1), max_n=6, max_w=8)

    # d=5 tests
    for c in [(2, 2, 1), (2, 1, 2), (1, 2, 2), (3, 1, 1)]:
        validate_gk(c, max_n=4, max_w=6)

    print("\n" + "=" * 60)
    print("All tests complete.")
