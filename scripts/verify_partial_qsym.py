#!/usr/bin/env python3
"""
verify_partial_qsym.py — Verify partial quasi-symmetry of g_tilde_n.

Claim: g_tilde_n (CPPs with max entry = n) lives in z_n · QSym(z_1,...,z_{n-1})[z_n].
That is, for each fixed power z_n^m, the coefficient polynomial in z_1,...,z_{n-1}
is quasi-symmetric.

Why: The relabelling argument for values {1,...,n-1} works with n held fixed.
Interlacing among non-n entries depends only on relative order, and any
non-n entry v < n, so inequalities involving n-entries are preserved.

Also checks L-positivity of each z_n^m slice.

Usage:
    python3 verify_partial_qsym.py              # all d <= 5, n = 1..3
    python3 verify_partial_qsym.py --fast       # all d <= 3, n = 1..2
    python3 verify_partial_qsym.py --dmax 4     # all d <= 4, n = 1..3
"""

import sys
from collections import defaultdict
from cpp_lib import (enumerate_cpps_with_max, compositions_of,
                     coarsens, M_to_L)


def M_to_L_bounded(M_coeffs, max_len):
    """
    Convert M-basis to L-basis in QSym_r (r = max_len variables).

    In finitely many variables, compositions with length > r have
    M_alpha = L_alpha = 0. We restrict the Mobius inversion to
    compositions with length <= max_len.
    """
    # Filter M_coeffs to those with length <= max_len
    filtered = {alpha: c for alpha, c in M_coeffs.items()
                if len(alpha) <= max_len}

    if max_len is None or max_len >= max(len(a) for a in M_coeffs):
        return M_to_L(filtered)

    # Do bounded Mobius inversion manually
    L_coeffs = defaultdict(int)
    all_totals = set(sum(c) for c in filtered.keys())

    for n in sorted(all_totals):
        if n == 0:
            continue
        comps_n = [c for c in compositions_of(n) if len(c) <= max_len]

        for alpha in comps_n:
            coeff = 0
            for beta in comps_n:
                if beta not in filtered:
                    continue
                if coarsens(beta, alpha):
                    sign = (-1) ** (len(beta) - len(alpha))
                    coeff += sign * filtered[beta]
            if coeff != 0:
                L_coeffs[alpha] = coeff

    return dict(L_coeffs)


def verify_partial_qsym_profile(c, max_entry, max_weight, verbose=False):
    """
    Verify partial quasi-symmetry for a single profile c at max_entry n.

    For each z_n^m slice (fixed number of parts = n), group CPPs by their
    reduced composition alpha = (a_1,...,a_k) of non-n value frequencies.
    The quasi-symmetry claim is: count depends only on alpha, not on which
    values in {1,...,n-1} carry those frequencies.

    Returns: (qsym_ok, l_pos_ok, n_slices, details)
    """
    n = max_entry
    r = n - 1  # number of non-n variables

    if n == 0:
        return True, True, 0, "trivial (n=0)"

    cpps = enumerate_cpps_with_max(c, n, max_weight)

    # Group by (m, placement, reduced_comp)
    # placement = tuple of actual non-n values (v_1, ..., v_k) with v_i < n
    by_m_placement_comp = defaultdict(int)

    for nu1, nu2, nu3, wt in cpps:
        freq = defaultdict(int)
        for part in nu1 + nu2 + nu3:
            if part > 0:
                freq[part] += 1
        m = freq.get(n, 0)
        non_max = sorted(v for v in freq if v < n and freq[v] > 0)
        placement = tuple(non_max)
        reduced_comp = tuple(freq[v] for v in non_max)
        by_m_placement_comp[(m, placement, reduced_comp)] += 1

    # Verify QSym: for fixed (m, reduced_comp), count is placement-independent
    by_m_alpha = defaultdict(dict)
    for (m, placement, alpha), count in by_m_placement_comp.items():
        by_m_alpha[(m, alpha)][placement] = count

    qsym_failures = []
    for (m, alpha), placement_counts in by_m_alpha.items():
        counts = list(placement_counts.values())
        if len(set(counts)) > 1:
            qsym_failures.append((m, alpha, placement_counts))

    qsym_ok = len(qsym_failures) == 0

    # Check L-positivity of each z_n^m slice (only meaningful if r >= 1)
    l_pos_ok = True
    l_neg_details = []

    # Collect all m values
    all_m = sorted(set(m for (m, _, _) in by_m_placement_comp))

    for m in all_m:
        if not qsym_ok:
            break

        # Get M-coefficients: for each alpha, pick any placement count
        M_coeffs = {}
        for (m2, alpha), placement_counts in by_m_alpha.items():
            if m2 != m:
                continue
            if not alpha:  # empty composition (only n-entries)
                continue
            if len(alpha) > r:  # can't have more distinct values than {1,...,r}
                continue
            # All placement counts should be equal (verified above)
            M_coeffs[alpha] = next(iter(placement_counts.values()))

        if not M_coeffs:
            continue

        if r < 2:
            # With 0 or 1 variables, L-positivity is trivially equivalent
            # to M-positivity (which holds since counts are nonneg)
            continue

        # Convert to L-basis in QSym_r
        L_coeffs = M_to_L_bounded(M_coeffs, max_len=r)
        for alpha, coeff in L_coeffs.items():
            if coeff < 0:
                l_pos_ok = False
                l_neg_details.append((m, alpha, coeff))

    n_slices = len(all_m)

    details = ""
    if not qsym_ok:
        details = f"{len(qsym_failures)} QSym failures"
        if verbose:
            for m, alpha, pc in qsym_failures:
                details += f"\n  m={m}, alpha={alpha}: {pc}"
    if not l_pos_ok:
        details += f"; {len(l_neg_details)} L-negative terms"
        if verbose:
            for m, alpha, coeff in l_neg_details[:10]:
                details += f"\n  m={m}, L_{alpha} = {coeff}"
            if len(l_neg_details) > 10:
                details += f"\n  ... and {len(l_neg_details) - 10} more"

    return qsym_ok, l_pos_ok, n_slices, details


def generate_profiles(d):
    """Generate all rank-3 profiles (c0, c1, c2) with c0+c1+c2 = d."""
    profiles = []
    for c0 in range(d + 1):
        for c1 in range(d - c0 + 1):
            c2 = d - c0 - c1
            profiles.append((c0, c1, c2))
    return profiles


def main():
    fast = '--fast' in sys.argv
    dmax = 5
    nmax = 3
    max_weight = 18

    for arg in sys.argv[1:]:
        if arg.startswith('--dmax='):
            dmax = int(arg.split('=')[1])
        elif arg.startswith('--nmax='):
            nmax = int(arg.split('=')[1])
        elif arg.startswith('--maxw='):
            max_weight = int(arg.split('=')[1])

    if fast:
        dmax = min(dmax, 3)
        nmax = min(nmax, 2)
        max_weight = min(max_weight, 12)

    print("=" * 70)
    print("Partial quasi-symmetry verification for g_tilde_n")
    print("=" * 70)
    print(f"Parameters: dmax={dmax}, nmax={nmax}, max_weight={max_weight}")
    print()
    print("Claim: For CPPs with max entry = n, each z_n^m slice")
    print("is quasi-symmetric in z_1,...,z_{n-1}.")
    print()

    total_profiles = 0
    total_qsym_ok = 0
    total_l_pos_ok = 0
    total_tests = 0
    results = []

    for d in range(1, dmax + 1):
        profiles = generate_profiles(d)
        for c in profiles:
            total_profiles += 1
            profile_qsym = True
            profile_lpos = True
            profile_slices = 0

            for n in range(1, nmax + 1):
                qsym_ok, l_pos_ok, n_slices, details = \
                    verify_partial_qsym_profile(c, n, max_weight, verbose=True)
                total_tests += 1
                profile_slices += n_slices
                if not qsym_ok:
                    profile_qsym = False
                if not l_pos_ok:
                    profile_lpos = False

                status_q = "OK" if qsym_ok else "FAIL"
                status_l = "OK" if l_pos_ok else "FAIL"
                label = f"c={c}, n={n}"
                print(f"  {label:20s}  QSym: {status_q:4s}  L-pos: {status_l:4s}  "
                      f"slices: {n_slices}", end="")
                if details:
                    print(f"  [{details}]")
                else:
                    print()

            if profile_qsym:
                total_qsym_ok += 1
            if profile_lpos:
                total_l_pos_ok += 1

            results.append((c, d, profile_qsym, profile_lpos, profile_slices))

    print()
    print("=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"Profiles tested:     {total_profiles}")
    print(f"QSym verified:       {total_qsym_ok}/{total_profiles}")
    print(f"L-positive:          {total_l_pos_ok}/{total_profiles}")
    print(f"Total (profile, n):  {total_tests}")
    print()

    print(f"{'d':>3s}  {'Profiles':>8s}  {'QSym OK':>7s}  {'L-pos OK':>8s}")
    print("-" * 32)
    for d in range(1, dmax + 1):
        d_results = [r for r in results if r[1] == d]
        n_prof = len(d_results)
        n_qsym = sum(1 for r in d_results if r[2])
        n_lpos = sum(1 for r in d_results if r[3])
        print(f"{d:>3d}  {n_prof:>8d}  {n_qsym:>7d}  {n_lpos:>8d}")

    total_qsym_fail = total_profiles - total_qsym_ok
    total_lpos_fail = total_profiles - total_l_pos_ok
    print("-" * 32)
    print(f"Tot  {total_profiles:>8d}  {total_qsym_ok:>7d}  {total_l_pos_ok:>8d}")

    if total_qsym_fail == 0 and total_lpos_fail == 0:
        print("\nAll checks passed: g_tilde_n is partially quasi-symmetric")
        print("and each z_n^m slice is L-positive in z_1,...,z_{n-1}.")
    elif total_qsym_fail == 0:
        print(f"\nPartial quasi-symmetry verified for all profiles.")
        print(f"L-positivity: {total_lpos_fail} profile(s) have negative L-coefficients.")
    else:
        print(f"\nFAILURES: {total_qsym_fail} QSym, {total_lpos_fail} L-pos")


if __name__ == '__main__':
    main()
