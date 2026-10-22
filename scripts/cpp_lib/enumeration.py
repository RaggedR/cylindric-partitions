"""
enumeration.py — Core CPP enumeration engine for rank-3 cylindric partitions.

Pure Python with int arithmetic only. This is the performance-critical hot path:
the triple-loop over (nu1, nu2, nu3) can enumerate millions of partitions.
No Sage types should enter these loops.

Exports:
    compute_gk(c, max_n, max_w, verbose) -> dict{(z,q): count}
    enumerate_cpps_with_max(c, max_entry, max_weight) -> list of triples
    extract_z_coefficient(gk_dict, n) -> dict{q_exp: coeff}
    gk_z1_coefficients(gk_dict, max_w) -> list
    compute_g_m(c, m, max_weight) -> dict{weight: count}
    compute_gk_both_stats(c, max_n, max_w, verbose) -> (gk_dict, gk_omega_dict)
"""

from collections import defaultdict
import time

from .partitions import get_part, gen_partitions, gen_compatible


def compute_gk(c, max_n, max_w, verbose=False):
    """
    Compute GK_c(z,q) for rank 3 by enumerating multipartitions.

    Parameters:
        c: (c0, c1, c2) -- profile
        max_n: max value of max(pi) to enumerate (z-degree truncation)
        max_w: max total weight |pi| to enumerate (q-degree truncation)
        verbose: print progress info

    Returns: dict { (n, w): count } where GK_c(z,q) = sum count * z^n * q^w
    """
    c0, c1, c2 = c

    result = defaultdict(int)
    total_count = 0

    t0 = time.time()

    for nu1 in gen_partitions(max_n, max_w):
        s1 = sum(nu1)
        rem1 = max_w - s1

        for nu2 in gen_compatible(nu1, c1, max_n, rem1):
            s2 = sum(nu2)
            rem2 = rem1 - s2

            for nu3 in gen_compatible(nu2, c2, max_n, rem2):
                s3 = sum(nu3)
                # Cyclic closure: nu3[j] >= nu1[j+c0] for all j
                ok = True
                for j in range(len(nu3)):
                    if nu3[j] < get_part(nu1, j + c0):
                        ok = False
                        break
                if ok and len(nu1) > len(nu3) + c0:
                    ok = False
                if not ok:
                    continue

                total = s1 + s2 + s3
                mx = max(get_part(nu1, 0), get_part(nu2, 0), get_part(nu3, 0))
                result[(mx, total)] += 1
                total_count += 1

    if verbose:
        elapsed = time.time() - t0
        print(f"  Enumerated {total_count} cylindric partitions in {elapsed:.2f}s")

    return dict(result)


def enumerate_cpps_with_max(c, max_entry, max_weight):
    """
    Enumerate all rank-3 CPPs with profile c and EXACT max entry = max_entry.

    Returns: list of (nu1, nu2, nu3, total_weight) tuples.
    """
    c0, c1, c2 = c
    results = []

    for nu1 in gen_partitions(max_entry, max_weight):
        s1 = sum(nu1)
        rem1 = max_weight - s1
        for nu2 in gen_compatible(nu1, c1, max_entry, rem1):
            s2 = sum(nu2)
            rem2 = rem1 - s2
            for nu3 in gen_compatible(nu2, c2, max_entry, rem2):
                s3 = sum(nu3)
                ok = True
                for j in range(len(nu3)):
                    if nu3[j] < get_part(nu1, j + c0):
                        ok = False
                        break
                if ok and len(nu1) > len(nu3) + c0:
                    ok = False
                if not ok:
                    continue
                mx = max(get_part(nu1, 0), get_part(nu2, 0),
                         get_part(nu3, 0))
                if mx == max_entry:
                    total = s1 + s2 + s3
                    results.append((nu1, nu2, nu3, total))

    return results


def extract_z_coefficient(gk_dict, n):
    """Extract [z^n] GK_c(z,q) as a dict {q_exp: coeff}."""
    result = {}
    for (zn, qw), coeff in gk_dict.items():
        if zn == n and coeff != 0:
            result[qw] = result.get(qw, 0) + coeff
    return result


def gk_z1_coefficients(gk_dict, max_w):
    """Compute GK_c(1,q) = sum_n [z^n] GK_c(z,q), returning list of max_w+1 coefficients."""
    result = [0] * (max_w + 1)
    for (zn, qw), coeff in gk_dict.items():
        if qw <= max_w:
            result[qw] += coeff
    return result


def compute_g_m(c, m, max_weight):
    """
    Compute g_m(q) = [z^m] GK_c(z,q) restricted to CPPs with exact max = m.

    Returns: dict {weight: count}.
    """
    if m == 0:
        return {0: 1}
    poly = defaultdict(int)
    cpps = enumerate_cpps_with_max(c, m, max_weight)
    for nu1, nu2, nu3, wt in cpps:
        poly[wt] += 1
    return dict(poly)


def compute_gk_both_stats(c, max_n, max_w, verbose=False):
    """
    Enumerate cylindric partitions once, recording both z-statistics:
      GK_c(z,q): z = max first part (standard)
      GK_c^omega(z,q): z = max length (omega-conjugate statistic)

    Returns (gk_dict, gk_omega_dict).
    """
    c0, c1, c2 = c
    gk = defaultdict(int)
    gk_omega = defaultdict(int)
    total_count = 0

    t0 = time.time()

    for nu1 in gen_partitions(max_n, max_w):
        s1 = sum(nu1)
        rem1 = max_w - s1

        for nu2 in gen_compatible(nu1, c1, max_n, rem1):
            s2 = sum(nu2)
            rem2 = rem1 - s2

            for nu3 in gen_compatible(nu2, c2, max_n, rem2):
                s3 = sum(nu3)
                ok = True
                for j in range(len(nu3)):
                    if nu3[j] < get_part(nu1, j + c0):
                        ok = False
                        break
                if ok and len(nu1) > len(nu3) + c0:
                    ok = False
                if not ok:
                    continue

                total = s1 + s2 + s3
                mx_part = max(get_part(nu1, 0), get_part(nu2, 0), get_part(nu3, 0))
                mx_len = max(len(nu1), len(nu2), len(nu3))

                gk[(mx_part, total)] += 1
                gk_omega[(mx_len, total)] += 1
                total_count += 1

    if verbose:
        elapsed = time.time() - t0
        print(f"  Enumerated {total_count} CPPs in {elapsed:.2f}s")

    return dict(gk), dict(gk_omega)
