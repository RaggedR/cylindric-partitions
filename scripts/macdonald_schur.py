"""macdonald_schur.py — Check Macdonald Schur positivity via Sage.

Key question: In the Macdonald-Borodin identity, do the Schur coefficients
(after expanding Macdonald P and Q functions in the Schur basis) have
coefficients in Z_≥0[q,t]?

Tests:
1. P_mu → Schur: Kostka-Macdonald K_{lam,mu}(q,t) — known Z_≥0[q,t] (Haiman)
2. Q_mu → Schur: K'_{lam,mu}(q,t) — check
3. Products P_lam · Q_mu → Schur — check
4. Sum of Q_mu (RPP generating function level by level) — check
"""

from sage.all import (QQ, PolynomialRing, Partitions, Partition,
                       SymmetricFunctions)


def setup():
    R = PolynomialRing(QQ, ['q', 't'])
    q, t = R.gens()
    F = R.fraction_field()
    Sym = SymmetricFunctions(F)
    s = Sym.s()
    h = Sym.h()
    Mac = Sym.macdonald(q=q, t=t)
    P = Mac.P()
    Q = Mac.Q()
    Ht = Mac.Ht()
    return R, q, t, F, Sym, s, h, P, Q, Ht


def check(coeff, R):
    """Returns (is_polynomial, is_nonneg, display_string)."""
    try:
        poly = R(coeff)
        coeffs = poly.coefficients()
        nonneg = all(c >= 0 for c in coeffs)
        return True, nonneg, poly
    except (TypeError, ValueError):
        return False, False, coeff


def test_P_schur():
    """Macdonald P in Schur basis = Kostka-Macdonald polynomials."""
    R, q, t, F, Sym, s, h, P, Q, Ht = setup()
    print("TEST 1: P_mu in Schur basis (Kostka-Macdonald)")
    print("─" * 50)
    ok = True
    for n in range(1, 5):
        for mu in Partitions(n):
            for lam, c in s(P(mu)):
                ip, nn, p = check(c, R)
                if not ip or not nn:
                    print(f"  FAIL: P_{mu} -> s_{lam}: {p}")
                    ok = False
    print(f"  K_{{lam,mu}}(q,t) ∈ Z_≥0[q,t] for |mu| ≤ 4: {'PASS' if ok else 'FAIL'}")


def test_Q_schur():
    """Macdonald Q in Schur basis."""
    R, q, t, F, Sym, s, h, P, Q, Ht = setup()
    print("\nTEST 2: Q_mu in Schur basis")
    print("─" * 50)
    for n in range(1, 4):
        for mu in Partitions(n):
            print(f"  Q_{mu} =")
            for lam, c in s(Q(mu)):
                ip, nn, p = check(c, R)
                tag = "✓" if (ip and nn) else ("poly" if ip else "RAT")
                print(f"    {tag} s_{lam}: {p}")


def test_Ht_schur():
    """Modified Macdonald Ht_mu in Schur basis (= modified Kostka)."""
    R, q, t, F, Sym, s, h, P, Q, Ht = setup()
    print("\nTEST 3: Ht_mu (modified Macdonald) in Schur basis")
    print("─" * 50)
    for n in range(1, 5):
        for mu in Partitions(n):
            print(f"  Ht_{mu} =")
            for lam, c in s(Ht(mu)):
                ip, nn, p = check(c, R)
                tag = "✓" if (ip and nn) else ("poly" if ip else "RAT")
                print(f"    {tag} s_{lam}: {p}")


def test_product_PQ():
    """Product P_lam · Q_mu in Schur basis."""
    R, q, t, F, Sym, s, h, P, Q, Ht = setup()
    print("\nTEST 4: Products P_lam · Q_mu in Schur basis")
    print("─" * 50)
    pairs = [
        ([1], [1]), ([2], [1]), ([1, 1], [1]),
        ([2], [2]), ([2, 1], [1]),
        ([2], [1, 1]), ([1, 1], [1, 1]),
    ]
    for l, m in pairs:
        lam, mu = Partition(l), Partition(m)
        product = s(P(lam) * Q(mu))
        print(f"\n  P_{lam} · Q_{mu} =")
        all_ok = True
        for kappa, c in product:
            ip, nn, p = check(c, R)
            tag = "✓" if (ip and nn) else ("NEG" if ip else "RAT")
            if not ip or not nn:
                all_ok = False
            print(f"    {tag} s_{kappa}: {p}")
        if all_ok:
            print(f"    → poly, nonneg ✓")


def test_sum_Q():
    """Sum of Q_mu for fixed |mu| = RPP generating function per level."""
    R, q, t, F, Sym, s, h, P, Q, Ht = setup()
    print("\nTEST 5: sum_{|mu|=n} Q_mu (RPP Macdonald GF, level n)")
    print("─" * 50)
    for n in range(1, 5):
        total = s.zero()
        for mu in Partitions(n):
            total += s(Q(mu))
        print(f"\n  Level {n}:")
        for lam, c in total:
            ip, nn, p = check(c, R)
            tag = "✓" if (ip and nn) else ("NEG" if ip else "RAT")
            print(f"    {tag} s_{lam}: {p}")


def test_skew_Q():
    """Skew Macdonald Q_{lam/mu} in Schur basis.

    The Borodin identity involves skew Macdonald functions.
    Check if their Schur expansions are q,t-positive.
    """
    R, q, t, F, Sym, s, h, P, Q, Ht = setup()
    print("\nTEST 6: Skew Q_{lam/mu} in Schur basis")
    print("─" * 50)
    # Q_lam.skew_by(Q_mu) = sum_nu c_{mu,nu}^{lam} Q_nu
    # In Schur basis, this should give q,t-polynomials.
    skew_pairs = [
        ([2], [1]), ([2, 1], [1]), ([2, 1], [2]),
        ([3], [1]), ([3, 1], [2]), ([2, 2], [1, 1]),
        ([3, 2, 1], [2, 1]),
    ]
    for l, m in skew_pairs:
        lam, mu = Partition(l), Partition(m)
        skewed = Q(lam).skew_by(Q(mu))
        schur_exp = s(skewed)
        print(f"\n  Q_{lam}/{mu} =")
        all_ok = True
        for kappa, c in schur_exp:
            ip, nn, p = check(c, R)
            tag = "✓" if (ip and nn) else ("NEG" if ip else "RAT")
            if not ip or not nn:
                all_ok = False
            print(f"    {tag} s_{kappa}: {p}")
        if all_ok:
            print(f"    → poly, nonneg ✓")
        else:
            print(f"    → FAILS")


def test_product_skew_Q():
    """Product of two skew Q functions in Schur basis.

    The Borodin identity for CPPs gives GK as a sum of PRODUCTS of skew
    Macdonald functions. Check if the product remains Schur positive.
    """
    R, q, t, F, Sym, s, h, P, Q, Ht = setup()
    print("\nTEST 7: Products of skew Q functions in Schur basis")
    print("These appear directly in the Macdonald-Borodin identity")
    print("─" * 50)

    # Product Q_{lam1/mu1} * Q_{lam2/mu2}
    test_cases = [
        (([2], [1]), ([2], [1])),         # Q_{(2)/(1)} * Q_{(2)/(1)}
        (([2, 1], [1]), ([2], [1])),       # Q_{(2,1)/(1)} * Q_{(2)/(1)}
        (([3], [1]), ([2, 1], [1])),       # Q_{(3)/(1)} * Q_{(2,1)/(1)}
        (([2, 1], [1]), ([2, 1], [1])),    # Q_{(2,1)/(1)} * Q_{(2,1)/(1)}
    ]

    for (l1, m1), (l2, m2) in test_cases:
        lam1, mu1, lam2, mu2 = Partition(l1), Partition(m1), Partition(l2), Partition(m2)
        q1 = Q(lam1).skew_by(Q(mu1))
        q2 = Q(lam2).skew_by(Q(mu2))
        product = s(q1 * q2)
        print(f"\n  Q_{lam1}/{mu1} · Q_{lam2}/{mu2} =")
        all_ok = True
        for kappa, c in product:
            ip, nn, p = check(c, R)
            tag = "✓" if (ip and nn) else ("NEG" if ip else "RAT")
            if not ip or not nn:
                all_ok = False
            print(f"    {tag} s_{kappa}: {p}")
        if all_ok:
            print(f"    → poly, nonneg ✓")
        else:
            print(f"    → FAILS positivity")


if __name__ == '__main__':
    test_P_schur()
    test_Q_schur()
    test_Ht_schur()
    test_product_PQ()
    test_sum_Q()
    test_skew_Q()
    test_product_skew_Q()
