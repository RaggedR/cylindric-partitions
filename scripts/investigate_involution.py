#!/usr/bin/env python3
"""
investigate_involution.py — The omega involution and Andrews-Gordon connection.

Key insight: The omega involution on symmetric functions sends e_k <-> h_k.
In q-series:
  E(z) = (zq;q)_inf  (elementary symmetric function GF)
  H(z) = 1/(zq;q)_inf  (complete homogeneous GF)

The E-H factorization: GK_c(z,q) = H(z) * F_c(z,q), where
  F_c(z,q) = E(z) * GK_c(z,q) = (zq;q)_inf * GK_c(z,q)

Reference: Warnaar (2021), Langer thesis (2021)
"""

from collections import defaultdict
import time

from cpp_lib.enumeration import (compute_gk, extract_z_coefficient,
                                  compute_gk_both_stats)
from cpp_lib.partitions import get_part
from cpp_lib.qseries import compute_Q, compute_h_m, compute_qpoch_m, qbinomial


# ---------------------------------------------------------------------------
# h_m analysis
# ---------------------------------------------------------------------------

def analyze_h_m(c, max_m, max_w, verbose=True):
    """
    Compute h_m for m = 0, ..., max_m and analyze:
    1. Is h_m a polynomial? (do coefficients become 0?)
    2. Are all coefficients nonneg?
    3. What is the degree?
    """
    gk = compute_gk(c, max_m + 2, max_w + 10, verbose=verbose)

    if verbose:
        d = sum(c)
        print(f"\n{'='*60}")
        print(f"h_m = (q)_m * g_m  for c = {c}, d = {d}")

    results = {}

    for m in range(max_m + 1):
        g_m = extract_z_coefficient(gk, m)
        h = compute_h_m(c, m, max_w, g_m_dict=g_m)

        # Find degree
        deg = 0
        for w in range(max_w, -1, -1):
            if h[w] != 0:
                deg = w
                break

        trailing_zeros = max_w - deg
        nonneg = all(h[w] >= 0 for w in range(max_w + 1))
        terms = [(w, h[w]) for w in range(max_w + 1) if h[w] != 0]

        results[m] = dict(h=h, deg=deg, trailing_zeros=trailing_zeros,
                          nonneg=nonneg, terms=terms)

        if verbose:
            term_str = " + ".join(f"{c}*q^{w}" for w, c in terms[:20])
            if len(terms) > 20:
                term_str += "..."
            poly_note = f"  (poly? trailing zeros={trailing_zeros})"
            nonneg_note = " nonneg" if nonneg else " HAS NEGATIVES"
            print(f"  m={m}: h_{m} = {term_str}")
            print(f"    deg={deg},{poly_note},{nonneg_note}")

    return results


# ---------------------------------------------------------------------------
# Rogers-Ramanujan verification for d=2
# ---------------------------------------------------------------------------

def verify_rogers_ramanujan(max_n, max_w, verbose=True):
    """
    Verify the Rogers-Ramanujan connection for d=2:
      c=(1,1,0): Q_m = q^{m^2}    (first Rogers-Ramanujan)
      c=(2,0,0): Q_m = q^{m(m+1)} (second Rogers-Ramanujan)
    """
    if verbose:
        print(f"\n{'='*60}")
        print("Rogers-Ramanujan verification for d=2")

    for c, quad_form, name in [
        ((1, 1, 0), lambda m: m * m, "q^{m^2}  (first RR)"),
        ((2, 0, 0), lambda m: m * (m + 1), "q^{m(m+1)}  (second RR)"),
    ]:
        if verbose:
            print(f"\n  c = {c}: expecting Q_m = {name}")

        gk = compute_gk(c, max_n + 2, max_w + max_n + 5, verbose=False)

        all_Q_pass = True
        for m in range(max_n + 1):
            Q_m = compute_Q(c, m, max_w, gk_dict=gk)
            expected_exp = quad_form(m)

            is_mono = True
            for w in range(max_w + 1):
                expected_coeff = 1 if w == expected_exp else 0
                if int(Q_m[w]) != expected_coeff:
                    is_mono = False
                    break

            g_m = extract_z_coefficient(gk, m)
            h = compute_h_m(c, m, max_w, g_m_dict=g_m)
            h_terms = [(w, h[w]) for w in range(max_w + 1) if h[w] != 0]
            h_nonneg = all(h[w] >= 0 for w in range(max_w + 1))
            h_str = " + ".join(f"{cv}*q^{w}" for w, cv in h_terms[:8])
            if len(h_terms) > 8:
                h_str += "..."

            if verbose:
                q_status = "PASS" if is_mono else "FAIL"
                h_status = "nonneg" if h_nonneg else "HAS NEGATIVES"
                print(f"    m={m}: Q_{m} = q^{expected_exp}  {q_status}")
                print(f"           h_{m} = {h_str}  ({h_status})")

            if not is_mono:
                all_Q_pass = False

        if verbose:
            status = "CONFIRMED" if all_Q_pass else "FAILED"
            print(f"  {status}: Q_m is a monomial ({name})")


# ---------------------------------------------------------------------------
# q-binomial decomposition
# ---------------------------------------------------------------------------

def qbinomial_decomposition(c, n, max_w, gk_dict=None, verbose=True):
    """
    Decompose Q_{n,c}(q) using the q-binomial formula:
    Q_n = sum_{j=0}^n (-1)^j q^{binom(j+1,2)} [n choose j]_q h_{n-j}
    """
    d = sum(c)
    if gk_dict is None:
        gk_max_n = n + 2
        gk_max_w = max_w + n * (n + 1) // 2 + 10
        gk_dict = compute_gk(c, gk_max_n, gk_max_w, verbose=False)

    if verbose:
        print(f"\n{'='*60}")
        print(f"q-binomial decomposition of Q_{{{n},{c}}}(q)")

    # Compute h_m for m = 0, ..., n
    h = {}
    for m in range(n + 1):
        g_m = extract_z_coefficient(gk_dict, m)
        h[m] = compute_h_m(c, m, max_w, g_m_dict=g_m)

    # Compute each term in the sum
    Q_reconstructed = [0] * (max_w + 1)

    for j in range(n + 1):
        sign = (-1) ** j
        shift = j * (j + 1) // 2
        qbin = qbinomial(n, j, max_w=max_w)
        hval = h[n - j]

        term = [0] * (max_w + 1)
        for w_q in range(max_w + 1):
            if qbin[w_q] == 0:
                continue
            for w_h in range(max_w + 1):
                if hval[w_h] == 0:
                    continue
                w_total = w_q + w_h + shift
                if w_total <= max_w:
                    term[w_total] += qbin[w_q] * hval[w_h]

        for w in range(max_w + 1):
            Q_reconstructed[w] += sign * term[w]

        if verbose:
            h_terms = [(w, hval[w]) for w in range(max_w + 1) if hval[w] != 0]
            h_str = " + ".join(f"{c}*q^{w}" for w, c in h_terms[:10])
            if len(h_terms) > 10:
                h_str += "..."

            term_terms = [(w, sign * term[w]) for w in range(max_w + 1) if term[w] != 0]
            t_str = " + ".join(f"{c}*q^{w}" for w, c in term_terms[:10])
            if len(term_terms) > 10:
                t_str += "..."

            sign_str = "+" if sign > 0 else "-"
            print(f"\n  j={j}: {sign_str} q^{shift} * [{n} choose {j}]_q * h_{n-j}")
            print(f"    h_{n-j} = {h_str}")
            print(f"    contribution = {t_str}")

    # Verify against compute_Q
    Q_actual = compute_Q(c, n, max_w, gk_dict=gk_dict)
    match = all(Q_reconstructed[w] == int(Q_actual[w])
                for w in range(max_w + 1))

    if verbose:
        q_str = " + ".join(f"{Q_reconstructed[w]}*q^{w}"
                           for w in range(max_w + 1) if Q_reconstructed[w] != 0)
        print(f"\n  Q_{n} (reconstructed) = {q_str}")
        print(f"  Matches compute_Q: {'PASS' if match else 'FAIL'}")

    return Q_reconstructed, match


# ---------------------------------------------------------------------------
# GK vs GK^omega comparison
# ---------------------------------------------------------------------------

def compare_gk_omega(c, max_n, max_w, verbose=True):
    """Compare GK_c(z,q) (z = max first part) with GK_c^omega(z,q) (z = max length)."""
    if verbose:
        print(f"\n{'='*60}")
        print(f"Comparing GK vs GK^omega for c = {c}, d = {sum(c)}")

    gk, gk_omega = compute_gk_both_stats(c, max_n, max_w, verbose=verbose)

    z1_gk = defaultdict(int)
    z1_omega = defaultdict(int)
    for (n, w), cnt in gk.items():
        z1_gk[w] += cnt
    for (n, w), cnt in gk_omega.items():
        z1_omega[w] += cnt

    z1_match = all(z1_gk[w] == z1_omega[w] for w in range(max_w + 1))
    same = (gk == gk_omega)

    if verbose:
        print(f"  GK(1,q) = GK^omega(1,q): {'PASS' if z1_match else 'FAIL (BUG!)'}")
        print(f"  GK(z,q) = GK^omega(z,q): {'YES (self-conjugate)' if same else 'NO (different z-refinement)'}")

        if not same:
            max_z = min(max_n, 5)
            for n in range(max_z + 1):
                gk_n = {w: cnt for (zn, w), cnt in gk.items() if zn == n}
                om_n = {w: cnt for (zn, w), cnt in gk_omega.items() if zn == n}
                gk_sum = sum(gk_n.values())
                om_sum = sum(om_n.values())
                print(f"    [z^{n}]: GK has {gk_sum} terms, GK^omega has {om_sum} terms")

    return same, z1_match


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    import sys
    fast = "--fast" in sys.argv

    print("=" * 60)
    print("Omega involution investigation")
    print("=" * 60)

    # 1. Rogers-Ramanujan verification
    print("\n" + "#" * 60)
    print("# 1. Rogers-Ramanujan verification (d=2)")
    print("#" * 60)
    verify_rogers_ramanujan(
        max_n=4 if not fast else 3,
        max_w=20 if not fast else 12
    )

    # 2. h_m polynomiality and nonnegativity
    print("\n" + "#" * 60)
    print("# 2. Corrected functions h_m = (q)_m * g_m")
    print("#" * 60)

    for c in [(1, 0, 0), (2, 1, 1), (2, 2, 1)]:
        mw = 15 if fast else 25
        mm = 2 if fast else 3
        if c == (1, 0, 0):
            mm = 3
            mw = 12
        analyze_h_m(c, max_m=mm, max_w=mw)

    # 3. q-binomial decomposition
    print("\n" + "#" * 60)
    print("# 3. q-binomial decomposition")
    print("#" * 60)
    for c in [(2, 1, 1), (2, 2, 1)]:
        qbinomial_decomposition(c, n=1, max_w=15)
    if not fast:
        qbinomial_decomposition((2, 1, 1), n=2, max_w=20)

    # 4. GK vs GK^omega
    print("\n" + "#" * 60)
    print("# 4. GK vs GK^omega (max-first-part vs max-length)")
    print("#" * 60)
    for c in [(1, 0, 0), (1, 1, 0), (2, 1, 1)]:
        mn = 4 if c == (1, 0, 0) else 3
        mw = 8 if fast else 10
        compare_gk_omega(c, max_n=mn, max_w=mw)

    print("\n" + "=" * 60)
    print("Investigation complete.")
