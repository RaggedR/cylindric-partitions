#!/usr/bin/env python3
"""
analyze_Q.py — Coefficient pattern analysis for Q_{n,c}(q).

For each computed Q_{n,c}(q):
- Test unimodality and log-concavity
- Check palindromic symmetry
- Decompose the alternating sum structure
- Compare across different profiles with the same d
- Kostka-Foulkes polynomial comparison

Reference: Warnaar (2021), Conjecture 2.7
"""

from collections import defaultdict
from math import gcd
import time

from cpp_lib.enumeration import compute_gk, extract_z_coefficient
from cpp_lib.profiles import profile_to_mu
from cpp_lib.qseries import compute_Q, compute_qpoch_m
from cpp_lib.partitions import partitions_of


# ---------------------------------------------------------------------------
# Coefficient analysis tools
# ---------------------------------------------------------------------------

def poly_to_list(Q, max_w):
    """Extract nonzero coefficients as list of (degree, coeff) pairs."""
    return [(w, int(Q[w])) for w in range(max_w + 1) if Q[w] != 0]


def is_unimodal(coeffs):
    """Test if a sequence of coefficients is unimodal."""
    if len(coeffs) <= 2:
        return True
    increasing = True
    for i in range(1, len(coeffs)):
        if coeffs[i] < coeffs[i - 1]:
            increasing = False
        elif coeffs[i] > coeffs[i - 1] and not increasing:
            return False
    return True


def is_log_concave(coeffs):
    """Test if a_i^2 >= a_{i-1} * a_{i+1} for all interior indices."""
    for i in range(1, len(coeffs) - 1):
        if coeffs[i] ** 2 < coeffs[i - 1] * coeffs[i + 1]:
            return False
    return True


def is_palindromic(coeffs):
    """Test if coefficient sequence is palindromic."""
    n = len(coeffs)
    for i in range(n // 2):
        if coeffs[i] != coeffs[n - 1 - i]:
            return False
    return True


def degree_and_valuation(Q, max_w):
    """Return (valuation, degree) of polynomial Q."""
    val = None
    deg = None
    for w in range(max_w + 1):
        if Q[w] != 0:
            if val is None:
                val = w
            deg = w
    return (val if val is not None else 0, deg if deg is not None else 0)


# ---------------------------------------------------------------------------
# Detailed analysis
# ---------------------------------------------------------------------------

def analyze_profile(c, max_n_range, max_w, verbose=True):
    """Compute and analyze Q_{n,c}(q) for a given profile."""
    d = sum(c)
    ell = gcd(d, 3)
    expected_q1 = (d + 1) * (d + 2) // 6 - 1

    if verbose:
        print(f"\n{'='*60}")
        print(f"Analysis of Q_{{n,c}}(q) for c = {c}, d = {d}")
        print(f"  ell = {ell}, Q(1) = {expected_q1}^n")

    # Compute GK
    gk_max_n = max_n_range + 2
    gk_max_w = max_w + max_n_range * (max_n_range + 1) // 2 + 3
    gk = compute_gk(c, gk_max_n, gk_max_w, verbose=verbose)

    results = {}

    for n in range(max_n_range + 1):
        Q = compute_Q(c, n, max_w, gk_dict=gk)

        terms = poly_to_list(Q, max_w)
        val, deg = degree_and_valuation(Q, max_w)

        if terms:
            coeff_seq = [int(Q[w]) for w in range(val, deg + 1)]
        else:
            coeff_seq = [0]

        total = sum(int(Q[w]) for w in range(max_w + 1))
        num_terms = len([c for c in coeff_seq if c != 0])

        unimodal = is_unimodal(coeff_seq) if coeff_seq else True
        log_conc = is_log_concave(coeff_seq) if len(coeff_seq) >= 3 else True
        palindrome = is_palindromic(coeff_seq) if coeff_seq else True

        is_mono = (num_terms <= 1)
        gap_free = all(c != 0 for c in coeff_seq) if coeff_seq else True

        res = dict(
            n=n, terms=terms, val=val, deg=deg,
            coeff_seq=coeff_seq, total=total, num_terms=num_terms,
            is_monomial=is_mono, unimodal=unimodal,
            log_concave=log_conc, palindromic=palindrome,
            gap_free=gap_free,
        )
        results[n] = res

        if verbose:
            print(f"\n  n = {n}:")
            if is_mono:
                if num_terms == 0:
                    print(f"    Q = 0")
                else:
                    print(f"    Q = {terms[0][1]}*q^{terms[0][0]}  (monomial)")
            else:
                q_str = " + ".join(f"{c}*q^{w}" for w, c in terms)
                if len(q_str) > 120:
                    q_str = q_str[:120] + "..."
                print(f"    Q = {q_str}")
                print(f"    val = {val}, deg = {deg}, #terms = {num_terms}")
                print(f"    Q(1) = {total}")
                print(f"    Unimodal: {unimodal}, Log-concave: {log_conc}, "
                      f"Palindromic: {palindrome}, Gap-free: {gap_free}")
                if len(coeff_seq) <= 40:
                    print(f"    Coefficients: {coeff_seq}")

    return results


def compare_profiles(profiles_list, max_n_range, max_w, verbose=True):
    """Compare Q_{n,c}(q) across profiles with the same d."""
    all_results = {}
    for c in profiles_list:
        all_results[c] = analyze_profile(c, max_n_range, max_w, verbose=verbose)

    if verbose and len(profiles_list) > 1:
        print(f"\n{'='*60}")
        print("Cross-profile comparison")
        for n in range(max_n_range + 1):
            print(f"\n  n = {n}:")
            for c in profiles_list:
                res = all_results[c][n]
                if res['is_monomial']:
                    desc = f"monomial q^{res['val']}" if res['num_terms'] > 0 else "0"
                else:
                    desc = (f"val={res['val']}, deg={res['deg']}, "
                            f"#terms={res['num_terms']}, Q(1)={res['total']}")
                print(f"    c={c}: {desc}")

    return all_results


# ---------------------------------------------------------------------------
# Soft positivity counterexample
# ---------------------------------------------------------------------------

def check_soft_positivity(max_n, max_w, verbose=True):
    """
    Test whether (q;q)_n * [z^n]((zq)_inf * f(z,q)) >= 0 for arbitrary f.

    Counterexample: f(z,q) = 1 + zq.
    """
    if verbose:
        print(f"\n{'='*60}")
        print("Soft positivity test: f(z,q) = 1 + zq")
        print("  (zq)_inf * (1+zq) = (zq^2;q)_inf")

    for n in range(max_n + 1):
        tri_shift = n * (n + 3) // 2
        sign = (-1) ** n

        result_exp = tri_shift
        result_coeff = sign

        if verbose:
            print(f"  n={n}: (q;q)_n * [z^n]((zq^2)_inf) = {result_coeff}*q^{result_exp}"
                  f"  {'>=0' if result_coeff >= 0 else '< 0 *** NEGATIVE ***'}")

    if verbose:
        print("  Conclusion: Positivity does NOT hold for arbitrary f(z,q) with nonneg coeffs.")
        print("  The cylindric partition structure of GK_c(z,q) is essential.")


# ---------------------------------------------------------------------------
# Kostka-Foulkes polynomials (charge statistic on SSYT)
# ---------------------------------------------------------------------------

def kostka_foulkes(lam, mu, max_deg=50):
    """
    Compute the Kostka-Foulkes polynomial K_{lam,mu}(q) using the charge
    statistic on semistandard Young tableaux.

    Returns: list of coefficients [a_0, a_1, ..., a_{max_deg}].
    """
    if sum(lam) != sum(mu):
        return [0] * (max_deg + 1)

    n_entries = len(mu)
    content = list(mu)
    rows = list(lam)
    n_rows = len(rows)

    result = [0] * (max_deg + 1)

    def _fill(tableau, row, col, remaining):
        if row >= n_rows:
            ch = _charge(tableau, lam, mu)
            if ch <= max_deg:
                result[ch] += 1
            return

        if col >= rows[row]:
            _fill(tableau, row + 1, 0, remaining)
            return

        min_entry = tableau[row][col - 1] if col > 0 else 1
        max_entry = n_entries

        for e in range(min_entry, max_entry + 1):
            if remaining[e - 1] <= 0:
                continue
            if row > 0 and col < rows[row - 1]:
                if e <= tableau[row - 1][col]:
                    continue
            remaining[e - 1] -= 1
            tableau[row][col] = e
            _fill(tableau, row, col + 1, remaining)
            remaining[e - 1] += 1

    tableau = [[0] * rows[i] for i in range(n_rows)]
    _fill(tableau, 0, 0, list(content))

    return result


def _charge(tableau, lam, mu):
    """Compute the charge of a semistandard Young tableau."""
    word = []
    for row in tableau:
        for entry in reversed(row):
            word.append(entry)

    if not word:
        return 0

    remaining = list(word)
    total_charge = 0

    while remaining:
        subword = []
        used = [False] * len(remaining)
        target = 1
        max_target = max(remaining)

        while target <= max_target:
            for i in range(len(remaining)):
                if not used[i] and remaining[i] == target:
                    subword.append((i, target))
                    used[i] = True
                    target += 1
                    break
            else:
                target += 1

        if not subword:
            break

        ch = 0
        running = 0
        for idx in range(1, len(subword)):
            if subword[idx][0] < subword[idx - 1][0]:
                running += 1
            ch += running

        total_charge += ch
        remaining = [remaining[i] for i in range(len(remaining)) if not used[i]]

    return total_charge


# ---------------------------------------------------------------------------
# Franklin-type involution analysis
# ---------------------------------------------------------------------------

def analyze_alternating_sum(c, n, max_w, gk_dict, verbose=True):
    """Decompose [z^n]((zq)_inf * GK_c(z,q)) to understand the sign structure."""
    d = sum(c)
    ell = gcd(d, 3)

    g = {}
    for m in range(n + 1):
        g[m] = extract_z_coefficient(gk_dict, m)

    pos_coeffs = defaultdict(int)
    neg_coeffs = defaultdict(int)

    for j in range(n + 1):
        tri = j * (j + 1) // 2
        sign = (-1) ** j

        # 1/(q)_j coefficients via DP
        inv_qj = _compute_inv_qj_list(j, max_w)

        gm = g.get(n - j, {})

        for w_gm, c_gm in gm.items():
            for w_inv, c_inv in enumerate(inv_qj):
                w_total = w_gm + tri + w_inv
                if w_total > max_w:
                    break
                contribution = c_gm * c_inv
                if sign > 0:
                    pos_coeffs[w_total] += contribution
                else:
                    neg_coeffs[w_total] += contribution

    if verbose:
        print(f"\n  Alternating sum decomposition for n={n}, c={c}:")
        print(f"    f_n(q) = P(q) - N(q) where:")

        max_show = min(max_w, 25)
        p_terms = [(w, pos_coeffs[w]) for w in range(max_show + 1) if pos_coeffs[w] != 0]
        n_terms = [(w, neg_coeffs[w]) for w in range(max_show + 1) if neg_coeffs[w] != 0]

        p_str = " + ".join(f"{c}*q^{w}" for w, c in p_terms[:15])
        n_str = " + ".join(f"{c}*q^{w}" for w, c in n_terms[:15])
        print(f"    P(q) = {p_str}{'...' if len(p_terms) > 15 else ''}")
        print(f"    N(q) = {n_str}{'...' if len(n_terms) > 15 else ''}")

        p_total = sum(pos_coeffs[w] for w in range(max_w + 1))
        n_total = sum(neg_coeffs[w] for w in range(max_w + 1))
        if p_total > 0:
            ratio = float(n_total / p_total)
            print(f"    P(1) = {p_total}, N(1) = {n_total}, cancellation ratio = {ratio:.4f}")

    return pos_coeffs, neg_coeffs


def _compute_inv_qj_list(j, max_w):
    """Compute 1/(q)_j as a list of int coefficients (partitions into parts <= j)."""
    if j == 0:
        result = [0] * (max_w + 1)
        result[0] = 1
        return result
    result = [0] * (max_w + 1)
    result[0] = 1
    for part in range(1, j + 1):
        for k in range(part, max_w + 1):
            result[k] += result[k - part]
    return result


def cancellation_ratio_table(profiles, max_n_range, max_w, verbose=True):
    """Compute cancellation ratio N(1)/P(1) across profiles and n values."""
    if verbose:
        print(f"\n{'='*60}")
        print("Cancellation ratio N(1)/P(1) table")
        header = "  c \\ n       " + "".join(f"  n={n:2d}  " for n in range(1, max_n_range + 1))
        print(header)
        print("  " + "-" * len(header))

    for c in profiles:
        gk_max_n = max_n_range + 2
        gk_max_w = max_w + max_n_range * (max_n_range + 1) // 2 + 3
        gk = compute_gk(c, gk_max_n, gk_max_w)

        ratios = []
        for n in range(1, max_n_range + 1):
            pos, neg = analyze_alternating_sum(c, n, max_w, gk, verbose=False)
            p_total = sum(pos[w] for w in range(max_w + 1))
            n_total = sum(neg[w] for w in range(max_w + 1))
            ratio = float(n_total / p_total) if p_total > 0 else 0.0
            ratios.append(ratio)

        if verbose:
            ratio_strs = [f"  {r:.4f}" for r in ratios]
            print(f"  {str(c):14s}" + "".join(ratio_strs))


# ---------------------------------------------------------------------------
# Kostka-Foulkes comparison
# ---------------------------------------------------------------------------

def compare_with_kostka(c, n, max_w, verbose=True):
    """Compare Q_{n,c}(q) with Kostka-Foulkes polynomials."""
    d = sum(c)
    if d % 3 == 0:
        return

    gk_max_n = n + 2
    gk_max_w = max_w + n * (n + 1) // 2 + 3
    gk = compute_gk(c, gk_max_n, gk_max_w)
    Q = compute_Q(c, n, max_w, gk_dict=gk)

    Q_terms = [(w, int(Q[w])) for w in range(max_w + 1) if Q[w] != 0]

    if verbose:
        print(f"\n  Kostka-Foulkes comparison for c={c}, n={n}:")
        q_str = " + ".join(f"{coeff}*q^{w}" for w, coeff in Q_terms)
        if len(q_str) > 100:
            q_str = q_str[:100] + "..."
        print(f"    Q_{{{n},{c}}}(q) = {q_str}")

    if not Q_terms:
        if verbose:
            print(f"    Q = 0, nothing to compare")
        return

    val = Q_terms[0][0]
    deg = Q_terms[-1][0]

    if verbose:
        total_size = deg
        if total_size <= 12:
            print(f"    Computing K_{{lam,(1^{total_size})}}(q) for small lam:")
            mu_col = tuple([1] * total_size) if total_size > 0 else (0,)
            for lam in partitions_of(total_size):
                if len(lam) > 5:
                    continue
                kf = kostka_foulkes(lam, mu_col, max_deg=max_w)
                kf_terms = [(w, kf[w]) for w in range(max_w + 1) if kf[w] != 0]
                if kf_terms:
                    kf_str = " + ".join(f"{c}*q^{w}" for w, c in kf_terms[:10])
                    print(f"      K_{{{lam}, {mu_col}}}(q) = {kf_str}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    import sys

    fast = "--fast" in sys.argv

    print("=" * 60)
    print("Q_{n,c}(q) coefficient analysis")
    print("=" * 60)

    # Soft positivity counterexample
    print("\n" + "#" * 60)
    print("# Soft positivity counterexample")
    print("#" * 60)
    check_soft_positivity(max_n=5, max_w=30)

    # d = 2
    print("\n" + "#" * 60)
    print("# d = 2: Q_{n,c}(q) is a monomial")
    print("#" * 60)
    compare_profiles(
        [(2, 0, 0), (1, 1, 0)],
        max_n_range=4 if not fast else 3,
        max_w=20 if not fast else 12
    )

    # d = 4
    print("\n" + "#" * 60)
    print("# d = 4: first nontrivial positivity")
    print("#" * 60)
    results_d4 = compare_profiles(
        [(2, 1, 1), (1, 2, 1)],
        max_n_range=2 if fast else 3,
        max_w=15 if fast else 25
    )

    # Alternating sum decomposition for d=4
    print("\n" + "#" * 60)
    print("# Alternating sum decomposition (d=4)")
    print("#" * 60)
    c = (2, 1, 1)
    mw = 15 if fast else 25
    gk = compute_gk(c, 5, mw + 10, verbose=True)
    for n in range(1, 3 if fast else 4):
        analyze_alternating_sum(c, n, mw, gk, verbose=True)

    # d = 5
    print("\n" + "#" * 60)
    print("# d = 5")
    print("#" * 60)
    results_d5 = compare_profiles(
        [(2, 2, 1), (3, 1, 1)],
        max_n_range=1 if fast else 2,
        max_w=12 if fast else 20
    )

    # d = 7
    print("\n" + "#" * 60)
    print("# d = 7")
    print("#" * 60)
    results_d7 = compare_profiles(
        [(3, 2, 2), (4, 2, 1)],
        max_n_range=1,
        max_w=12 if fast else 18
    )

    # Alternating sum for d=5
    if not fast:
        print("\n" + "#" * 60)
        print("# Alternating sum decomposition (d=5)")
        print("#" * 60)
        for c in [(2, 2, 1), (3, 1, 1)]:
            gk = compute_gk(c, 4, 25, verbose=True)
            for n in range(1, 3):
                analyze_alternating_sum(c, n, 20, gk, verbose=True)

    # Cancellation ratio table
    print("\n" + "#" * 60)
    print("# Cancellation ratio N(1)/P(1)")
    print("#" * 60)
    ratio_profiles = [(2, 1, 1), (2, 2, 1), (3, 1, 1)]
    if not fast:
        ratio_profiles.append((3, 2, 2))
    cancellation_ratio_table(
        ratio_profiles,
        max_n_range=2 if fast else 3,
        max_w=12 if fast else 20
    )

    # Kostka-Foulkes comparison for d=4
    print("\n" + "#" * 60)
    print("# Kostka-Foulkes polynomial comparison")
    print("#" * 60)
    compare_with_kostka((2, 1, 1), n=1, max_w=15)
    compare_with_kostka((2, 2, 1), n=1, max_w=15)
    if not fast:
        compare_with_kostka((2, 1, 1), n=2, max_w=20)

    print("\n" + "=" * 60)
    print("Analysis complete.")
