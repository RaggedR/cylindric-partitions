"""
qsym.py — Quasi-symmetric function utilities.

Provides M-to-L basis conversion (Mobius inversion on the refinement poset),
composition enumeration, and refinement testing.

Pure Python (no Sage dependency) — the hot paths are combinatorial tests
on small compositions where Sage overhead isn't worth it.
"""

from collections import defaultdict


def compositions_of(n, max_len=None):
    """Generate all compositions of n (ordered sequences summing to n)."""
    if n == 0:
        yield ()
        return
    if max_len is None:
        max_len = n
    if max_len <= 0:
        return
    for first in range(1, n + 1):
        if max_len == 1:
            if first == n:
                yield (first,)
        else:
            for rest in compositions_of(n - first, max_len - 1):
                yield (first,) + rest


def refines(beta, alpha):
    """
    Check if composition beta refines alpha.

    beta refines alpha if alpha can be obtained by merging consecutive
    parts of beta. E.g., (1,2,1) refines (3,1) since 1+2=3, 1=1.
    """
    if sum(beta) != sum(alpha):
        return False
    j = 0
    for a in alpha:
        s = 0
        while s < a and j < len(beta):
            s += beta[j]
            j += 1
        if s != a:
            return False
    return j == len(beta)


def coarsens(alpha, beta):
    """Check if alpha coarsens beta (i.e., beta refines alpha)."""
    return refines(beta, alpha)


def q_specialization(alpha):
    """
    q-weight of composition alpha under z_k -> q^k.

    Returns sum_{i} (i+1) * alpha_i = sum_k k * f_k where f_k = alpha_{k-1}.
    """
    return sum((i + 1) * alpha[i] for i in range(len(alpha)))


def M_to_L(M_coeffs):
    """
    Convert M-basis coefficients to L-basis coefficients via Mobius
    inversion on the refinement poset of compositions.

    L_alpha = sum_{beta refines alpha} M_beta, so inverting:
    L_coeff(alpha) = sum_{beta coarsens alpha} (-1)^{l(beta)-l(alpha)} M_coeff(beta)
    """
    L_coeffs = defaultdict(int)
    all_totals = set(sum(c) for c in M_coeffs.keys())

    for n in sorted(all_totals):
        if n == 0:
            continue
        comps_n = list(compositions_of(n))

        for alpha in comps_n:
            coeff = 0
            for beta in comps_n:
                if beta not in M_coeffs:
                    continue
                if coarsens(beta, alpha):
                    sign = (-1) ** (len(beta) - len(alpha))
                    coeff += sign * M_coeffs[beta]
            if coeff != 0:
                L_coeffs[alpha] = coeff

    return dict(L_coeffs)


def M_to_L_bounded(M_coeffs, max_len):
    """
    Convert M-basis to L-basis in QSym_r (r = max_len variables).

    In finitely many variables, compositions with length > r have
    M_alpha = L_alpha = 0. We restrict the Mobius inversion to
    compositions with length <= max_len.
    """
    filtered = {alpha: c for alpha, c in M_coeffs.items()
                if len(alpha) <= max_len}

    if not filtered:
        return {}

    if max_len >= max(len(a) for a in M_coeffs):
        return M_to_L(filtered)

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
