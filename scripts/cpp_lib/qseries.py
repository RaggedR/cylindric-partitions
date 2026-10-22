"""
qseries.py — q-series algebra using SageMath.

Replaces ~200 lines of hand-rolled polynomial arithmetic (Fraction-based
truncated power series, manual q-Pochhammer, etc.) with clean Sage
polynomial ring operations.

Design principle: enumeration outputs dict{int: int}, then gk_dict_to_poly()
converts to a Sage polynomial. All algebra after that point uses Sage.

Exports:
    qpoch(m) -> (q)_m as Sage polynomial
    qbinomial(n, k) -> [n choose k]_q as Sage polynomial
    compute_qpoch_m(m, max_w) -> list of int coefficients (pure Python fallback)
    gk_dict_to_poly(coeffs) -> Sage polynomial
    compute_h_m(c, m, max_w, g_m_dict) -> list of int coefficients
    compute_Q(c, n, max_w, gk_dict, verbose) -> list of int coefficients
    borodin_product(c, N) -> list of int coefficients
"""

from sage.all import ZZ, PolynomialRing, PowerSeriesRing, prod
from sage.combinat.q_analogues import q_binomial as sage_q_binomial

from .enumeration import compute_gk, extract_z_coefficient, compute_g_m, enumerate_cpps_with_max
from .profiles import profile_to_mu

# Global Sage rings
_R = PolynomialRing(ZZ, 'q')
q = _R.gen()

_PS = PowerSeriesRing(ZZ, 'q')


def qpoch(m):
    """(q)_m = prod_{k=1}^m (1 - q^k) as a Sage polynomial."""
    return prod(1 - q**k for k in range(1, m + 1))


def qbinomial(n, k, max_w=None):
    """
    Gaussian binomial [n choose k]_q as a Sage polynomial.

    If max_w is given, returns a list of int coefficients of length max_w+1.
    Otherwise returns a Sage polynomial.
    """
    result = sage_q_binomial(n, k)
    if max_w is not None:
        poly = _R(result)
        coeffs = [0] * (max_w + 1)
        for exp, c in poly.dict().items():
            if exp <= max_w:
                coeffs[exp] = int(c)
        return coeffs
    return _R(result)


def compute_qpoch_m(m, max_w):
    """
    Compute (q)_m = prod_{k=1}^m (1 - q^k) as a list of int coefficients.

    Pure Python fallback for contexts that don't need full Sage polynomials.
    """
    result = [0] * (max_w + 1)
    result[0] = 1
    for k in range(1, m + 1):
        if k > max_w:
            break
        new = list(result)
        for i in range(max_w + 1 - k):
            new[i + k] -= result[i]
        result = new
    return result


def gk_dict_to_poly(coeffs):
    """Convert a {q_power: count} dict to a Sage polynomial in q."""
    if not coeffs:
        return _R.zero()
    return _R(coeffs)


def compute_h_m(c, m, max_w, g_m_dict=None):
    """
    Compute h_m(q) = (q)_m * g_m(q) as a list of int coefficients.

    h_m are the "corrected" GF coefficients: nonneg polynomials for all
    tested profiles, and Q_n = sum (-1)^j q^{j(j+1)/2} [n,j]_q h_{n-j}.

    Returns: list of int coefficients [a_0, a_1, ..., a_{max_w}].
    """
    if g_m_dict is None:
        g_m_dict = compute_g_m(c, m, max_w)

    qm = compute_qpoch_m(m, max_w)

    h_m = [0] * (max_w + 1)
    for w_g, c_g in g_m_dict.items():
        if w_g > max_w:
            continue
        for w_q in range(max_w + 1 - w_g):
            if qm[w_q] != 0:
                h_m[w_g + w_q] += c_g * qm[w_q]

    return h_m


def compute_Q(c, n, max_w, gk_dict=None, verbose=False):
    """
    Compute Q_{n,c}(q) for given profile c and z-power n.

    Q_{n,c}(q) = (q^ell; q^ell)_n * [z^n]((zq;q)_inf * GK_c(z,q))

    where ell = gcd(d, 3), d = c0+c1+c2.

    Uses Sage power series for exact polynomial arithmetic.

    Returns: list of int coefficients [a_0, a_1, ..., a_{max_w}]
    """
    from math import gcd
    d = sum(c)
    ell = gcd(d, 3)

    if gk_dict is None:
        gk_dict = compute_gk(c, max_n=n + max_w, max_w=max_w + n * (n + 1), verbose=verbose)

    # Extract g_m(q) for m = 0, ..., n
    g = {}
    for m in range(n + 1):
        g[m] = extract_z_coefficient(gk_dict, m)

    # Compute [z^n]((zq;q)_inf * GK_c(z,q))
    # = sum_{j=0}^n (-1)^j * q^{j(j+1)/2} / (q)_j * g_{n-j}(q)
    #
    # Use power series truncated to O(q^{max_w+1}).
    prec = max_w + 1
    PS = PowerSeriesRing(ZZ, 'q', default_prec=prec)
    q_ps = PS.gen()

    fn = PS.zero()
    for j in range(n + 1):
        tri = j * (j + 1) // 2
        sign = (-1) ** j

        # 1/(q)_j as power series = generating function for partitions
        # into parts <= j
        inv_qj = _compute_inv_qj_series(j, PS, prec)

        # g_{n-j}(q) as power series
        gm = g.get(n - j, {})
        gm_ps = PS(gm) if gm else PS.zero()

        fn += sign * q_ps**tri * inv_qj * gm_ps

    # Multiply by (q^ell; q^ell)_n = prod_{k=1}^n (1 - q^{k*ell})
    qell_n = PS.one()
    for k in range(1, n + 1):
        exp = k * ell
        if exp < prec:
            qell_n *= (1 - q_ps**exp)

    Q_ps = fn * qell_n

    # Extract coefficients as integers
    Q = [0] * (max_w + 1)
    for w in range(max_w + 1):
        Q[w] = int(Q_ps[w])

    return Q


def _compute_inv_qj_series(j, PS, prec):
    """
    Compute 1/(q)_j as a power series truncated to O(q^prec).

    The coefficient of q^k is the number of partitions of k into parts <= j.
    """
    if j == 0:
        return PS.one()

    # Build as a product of geometric series
    result = PS.one()
    q_ps = PS.gen()
    for k in range(1, j + 1):
        if k < prec:
            result *= ~(1 - q_ps**k)  # Sage inversion in power series ring

    return result


def borodin_product(c, N):
    """
    GK_c(1,q) via Borodin's product formula (Warnaar eq 3.15).

    Returns list of N int coefficients.
    """
    r = len(c)
    d = sum(c)
    mu = profile_to_mu(c)
    qdr = d + r

    PS = PowerSeriesRing(ZZ, 'q', default_prec=N)
    q_ps = PS.gen()

    def theta_ps(a, m):
        """Jacobi theta: (q^a; q^m)_inf * (q^{m-a}; q^m)_inf as power series."""
        t1 = PS.one()
        t2 = PS.one()
        k = 0
        while a + k * m < N:
            t1 *= (1 - q_ps**(a + k * m))
            k += 1
        k = 0
        while (m - a) + k * m < N:
            t2 *= (1 - q_ps**((m - a) + k * m))
            k += 1
        return t1 * t2

    def qpoch_inf_ps(a, m):
        """(q^a; q^m)_inf as power series."""
        result = PS.one()
        k = 0
        while a + k * m < N:
            result *= (1 - q_ps**(a + k * m))
            k += 1
        return result

    # Numerator: (q^{d+r})_inf^{r-1}
    qdr_inf = qpoch_inf_ps(qdr, qdr)
    num = qdr_inf ** (r - 1)

    # Theta functions
    for i in range(r):
        for j in range(i + 1, r):
            a = mu[i] - mu[j] + j - i
            num *= theta_ps(a, qdr)

    # Denominator: (q)_inf^r
    q_inf = qpoch_inf_ps(1, 1)
    denom = q_inf ** r

    result_ps = num * ~denom  # Sage power series inversion

    return [int(result_ps[k]) for k in range(N)]
