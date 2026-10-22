#!/usr/bin/env python3
"""
hook_multiplicities.py — Compute cylindric hook multiplicities for each profile.

For a cylindric shape with profile c, the Borodin product formula gives:
    GK_c(1,q) = prod_{s in hat{lambda}} 1/(1 - q^{h(s)}) = prod_{k>=1} (1-q^k)^{-n_k}
where n_k = number of cells with hook length k.

This script:
1. Computes GK_c(1,q) via Borodin's product formula (high precision)
2. Extracts n_k by iterative division (power series arithmetic)
3. Verifies against the analytical formula from theta/Pochhammer residues
4. Displays the hook multiplicity pattern for each profile

The key finding: n_k is eventually periodic with period T = d+3.
"""

import sys
sys.path.insert(0, '.')

from cpp_lib.profiles import profile_to_mu, all_profiles, representative_profiles
from cpp_lib.qseries import borodin_product

from sage.all import ZZ, PowerSeriesRing


def extract_hook_multiplicities(borodin_coeffs, max_k=None):
    """
    Given GK_c(1,q) as a list of coefficients, extract n_k such that
    GK_c(1,q) = prod_{k>=1} (1-q^k)^{-n_k}.

    Algorithm: iteratively divide by (1-q^k)^{-n_k} for k=1,2,...
    At step k, n_k = coefficient of q^k in the remaining series.
    """
    N = len(borodin_coeffs)
    if max_k is None:
        max_k = N // 2  # safe upper bound

    PS = PowerSeriesRing(ZZ, 'q', default_prec=N)
    q = PS.gen()

    # Build the power series from coefficients
    f = PS(borodin_coeffs)

    multiplicities = {}
    for k in range(1, max_k + 1):
        # n_k is the coefficient of q^k in the current series
        n_k = int(f[k])
        multiplicities[k] = n_k

        if n_k != 0:
            # Multiply by (1-q^k)^{n_k} to cancel this factor
            factor = (1 - q**k)**n_k
            f = f * factor
            # After multiplication, f should have [q^k] = 0

    return multiplicities


def analytical_hook_multiplicities(c, max_k):
    """
    Compute hook multiplicities analytically from the Borodin formula structure.

    For rank 3, GK_c(1,q) = Num/Den where:
        Den = (q;q)_inf^3  (exponent -3 for each k)
        Num = (q^T;q^T)_inf^2 * prod_{i<j} theta(q^{a_ij}; q^T)

    The exponent of (1-q^k) in GK_c(1,q) is -n_k = (numerator exp) - 3.

    For k large enough, n_k depends only on k mod T.
    """
    r = len(c)
    d = sum(c)
    T = d + r  # period
    mu = profile_to_mu(c)

    # Compute the a_{ij} parameters
    a_values = []
    for i in range(r):
        for j in range(i + 1, r):
            a = mu[i] - mu[j] + j - i
            a_values.append(a)

    # Residues from theta functions: each theta(q^a; q^T) contributes
    # (q^a; q^T)_inf * (q^{T-a}; q^T)_inf
    # So residues are: {a mod T, (T-a) mod T} for each a
    theta_residues = []
    for a in a_values:
        theta_residues.append(a % T)
        theta_residues.append((T - a) % T)

    # Count residue multiplicities
    residue_count = {}
    for res in theta_residues:
        residue_count[res] = residue_count.get(res, 0) + 1

    # Add (q^T;q^T)_inf^2 contribution (residue 0, multiplicity 2)
    residue_count[0] = residue_count.get(0, 0) + 2

    # For each k, n_k = 3 - (numerator exponent at residue k mod T)
    # But only for k >= a (the Pochhammer starts at q^a, not q^0)
    # For the eventual periodic part:
    periodic_n = {}
    for res in range(T):
        nu = residue_count.get(res, 0)
        periodic_n[res] = r - nu  # r=3 for rank 3

    # Now compute actual n_k accounting for initial terms
    # The theta factor (q^a; q^T)_inf starts at q^a, so for k < a it doesn't contribute
    multiplicities = {}
    for k in range(1, max_k + 1):
        res = k % T
        # Start with base: 3 from denominator (q;q)_inf^3
        n_k = r

        # Subtract numerator contributions
        # (q^T;q^T)_inf^2: contributes when T | k, i.e., k = T, 2T, ...
        if k % T == 0 and k > 0:
            n_k -= 2

        # Theta functions: each (q^a; q^T)_inf contributes when k = a, a+T, a+2T, ...
        # i.e., k >= a and k ≡ a mod T
        for a in a_values:
            a_mod = a % T
            if a_mod == 0:
                a_eff = T  # (q^T; q^T) — but this is already in the (q^T;q^T)^2 factor
                # Actually, if a ≡ 0 mod T, then (q^a; q^T) = (q^T; q^T) shifted
                # Need to be more careful here
                pass
            # (q^a; q^T)_inf: contributes at k = a, a+T, a+2T, ...
            if k >= a and (k - a) % T == 0:
                n_k -= 1
            # (q^{T-a}; q^T)_inf: contributes at k = T-a, 2T-a, ...
            ta = T - a
            if ta > 0 and k >= ta and (k - ta) % T == 0:
                n_k -= 1

        multiplicities[k] = n_k

    return multiplicities, periodic_n, a_values, T


def main():
    print("=" * 70)
    print("Hook multiplicities for rank-3 cylindric shapes")
    print("=" * 70)

    max_k = 30  # compute multiplicities up to hook length 30
    N = max_k * 3  # precision for Borodin formula

    # Test profiles organized by d
    test_profiles = {
        1: [(1, 0, 0)],
        2: [(2, 0, 0), (1, 1, 0)],
        4: [(2, 1, 1), (1, 2, 1), (4, 0, 0), (2, 2, 0)],
        5: [(2, 2, 1), (3, 1, 1)],
    }

    for d in sorted(test_profiles.keys()):
        for c in test_profiles[d]:
            print(f"\n{'─' * 70}")
            print(f"Profile c = {c}, d = {d}, T = d+3 = {d+3}")
            print(f"{'─' * 70}")

            mu = profile_to_mu(c)
            print(f"  mu = {mu}")

            # Compute a_ij values
            r = 3
            a_vals = []
            for i in range(r):
                for j in range(i + 1, r):
                    a = mu[i] - mu[j] + j - i
                    a_vals.append(a)
            print(f"  a values: a_01={a_vals[0]}, a_02={a_vals[1]}, a_12={a_vals[2]}")
            T = d + r
            print(f"  Theta residues mod {T}: {sorted([a % T for a in a_vals] + [(T - a) % T for a in a_vals])}")

            # Method 1: Extract from Borodin product (numerical)
            bp = borodin_product(c, N)
            n_numerical = extract_hook_multiplicities(bp, max_k)

            # Method 2: Analytical formula
            n_analytical, periodic_n, _, _ = analytical_hook_multiplicities(c, max_k)

            # Display
            print(f"\n  Hook multiplicities n_k (k=1..{max_k}):")
            print(f"  {'k':>3} {'n_k(num)':>8} {'n_k(ana)':>8} {'match':>6}")
            print(f"  {'─'*30}")

            all_match = True
            for k in range(1, max_k + 1):
                nn = n_numerical.get(k, '?')
                na = n_analytical.get(k, '?')
                match = "✓" if nn == na else "✗"
                if nn != na:
                    all_match = False
                # Only print first 2 periods + a few extra
                if k <= 2 * T + 3:
                    print(f"  {k:>3} {nn:>8} {na:>8} {match:>6}")

            if max_k > 2 * T + 3:
                print(f"  ... (remaining {max_k - 2*T - 3} values all match)" if all_match else "  ... (some mismatches!)")

            # Show periodic pattern
            print(f"\n  Eventual periodic pattern (period T={T}):")
            print(f"  residue:      {' '.join(f'{r:>3}' for r in range(T))}")
            vals = [periodic_n.get(r, '?') for r in range(T)]
            print(f"  n_k:          {' '.join(f'{v:>3}' for v in vals)}")
            print(f"  Sum per period: {sum(v for v in vals if isinstance(v, int))}")

            # Check: is it constant?
            unique_vals = set(v for v in vals if isinstance(v, int))
            if len(unique_vals) == 1:
                print(f"  ** CONSTANT: n_k = {unique_vals.pop()} for all k (in periodic part)")
            else:
                print(f"  ** NON-CONSTANT: n_k takes values {sorted(unique_vals)}")


def analyze_module_structure():
    """
    Analyze the Sym-module structure of prod_k (1-z_k)^{-n_k}.

    Key observations:
    1. n_k is eventually periodic with period T = d+3
    2. Within each residue class mod T, the generating function is symmetric
    3. The "stabilizer" of the multiplicity sequence determines the module

    The GF decomposes as:
       prod_k (1-z_k)^{-n_k} = prod_{r=0}^{T-1} sigma[n_r * Z_r]
    where Z_r = z_r + z_{r+T} + z_{r+2T} + ... is the r-th "column" alphabet.
    """
    print("\n" + "=" * 70)
    print("Module structure analysis")
    print("=" * 70)

    for d in [1, 2, 3, 4, 5]:
        profiles = all_profiles(d)
        T = d + 3
        print(f"\n{'━' * 70}")
        print(f"d = {d}, T = {T}, all {len(profiles)} profiles")
        print(f"{'━' * 70}")

        # Collect periodic patterns for all profiles
        patterns = {}
        for c in profiles:
            _, periodic_n, a_values, _ = analytical_hook_multiplicities(c, T)
            pattern = tuple(periodic_n.get(r, 0) for r in range(T))
            patterns[c] = pattern

        # Group by sorted pattern (same multiset)
        from collections import defaultdict
        multiset_groups = defaultdict(list)
        for c, pat in patterns.items():
            multiset_groups[tuple(sorted(pat))].append((c, pat))

        for ms, group in sorted(multiset_groups.items()):
            print(f"\n  Multiplicity multiset: {ms}")
            print(f"    Sum = {sum(ms)}, expected 3d+1 = {3*d+1}")
            for c, pat in group:
                mu = profile_to_mu(c)
                print(f"    c={c}, mu={mu}: pattern={pat}")

            # Analyze stabilizer structure
            # The stabilizer of the pattern under S_T (permuting residue classes)
            # is the Young subgroup S_lambda where lambda groups equal values
            from collections import Counter
            val_counts = Counter(ms)
            stab_parts = sorted(val_counts.values(), reverse=True)
            print(f"    Stabilizer: S_{stab_parts} in S_{T}")
            print(f"    Orbit size (# rearrangements): {T}")
            # Multinomial: T! / prod(stab_parts!)
            from math import factorial
            orbit_size = factorial(T)
            for p in stab_parts:
                orbit_size //= factorial(p)
            print(f"    Actual orbit size: {orbit_size}")

        # Check which patterns are cyclic rotations of each other
        print(f"\n  Cyclic orbit analysis:")
        seen = set()
        for c in profiles:
            if c in seen:
                continue
            pat = patterns[c]
            # Find all cyclic rotations of c that give the same pattern
            cyclic_group = []
            for c2 in profiles:
                if patterns[c2] == pat:
                    cyclic_group.append(c2)
                    seen.add(c2)
            if len(cyclic_group) > 1:
                print(f"    Same pattern {pat}: {cyclic_group}")


def analyze_weighted_alphabet():
    """
    For each profile, express the weighted alphabet sum_k n_k z_k
    in terms of the periodic structure.
    """
    print("\n" + "=" * 70)
    print("Weighted alphabet decomposition")
    print("=" * 70)

    for d in [1, 2, 4, 5]:
        T = d + 3
        reps = representative_profiles().get(d, [])
        if not reps:
            continue

        print(f"\n{'─' * 70}")
        print(f"d = {d}, T = {T}")

        for c in reps:
            _, periodic_n, a_values, _ = analytical_hook_multiplicities(c, T)
            pattern = [periodic_n.get(r, 0) for r in range(T)]

            print(f"\n  Profile c = {c}")
            print(f"  Periodic multiplicities: {pattern}")

            # Decompose: group residue classes by multiplicity
            from collections import defaultdict
            groups = defaultdict(list)
            for r in range(T):
                groups[pattern[r]].append(r)

            print(f"  Alphabet decomposition:")
            for mult, residues in sorted(groups.items()):
                res_str = ", ".join(str(r) for r in residues)
                print(f"    Multiplicity {mult}: residues {{{res_str}}}")
                # In the eventual periodic part:
                # These residues contribute n_r * (z_r + z_{r+T} + z_{r+2T} + ...)
                # = mult * sum of variables in these columns

            # The generating function is:
            # prod_{mult m} prod_{r with n_r=m} sigma_1[m * Z_r]
            # = prod_{mult m} sigma_1[m * (sum of Z_r for r with n_r=m)]
            # = prod_{mult m} Sym^{(m)} evaluated at the combined alphabet
            print(f"  Structure: ", end="")
            parts = []
            for mult, residues in sorted(groups.items()):
                n_cols = len(residues)
                parts.append(f"Sym^({mult}) in {n_cols} columns")
            print(" ⊗ ".join(parts))


if __name__ == '__main__':
    main()
    analyze_module_structure()
    analyze_weighted_alphabet()
