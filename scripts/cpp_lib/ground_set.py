"""
ground_set.py — (S, U, nu) triple enumeration for Q_{n,c}.

The ground set G_{n,c} consists of triples (S, U, nu) where:
    S subset {1,...,n}, |S| = j
    U subset {1,...,n-j}
    nu is a CPP with profile c and max entry = n-j

Each triple has sign (-1)^{|S|+|U|} and q-weight sigma(S) + sigma(U) + |nu|.

Pure Python with int arithmetic.
"""

from collections import defaultdict
from itertools import combinations

from .enumeration import enumerate_cpps_with_max


def frequency_vector(S, U, nu1, nu2, nu3):
    """
    Compute the frequency vector of a (S, U, nu) triple.

    Counts how many objects appear at each integer value k >= 1,
    across elements of S, elements of U, and parts of (nu1, nu2, nu3).

    Returns: tuple (f_1, f_2, ..., f_max) where f_k = count at value k.
             Returns () for the empty triple.
    """
    freq = defaultdict(int)
    for s in S:
        freq[s] += 1
    for u in U:
        freq[u] += 1
    for part in nu1 + nu2 + nu3:
        freq[part] += 1
    if not freq:
        return ()
    max_val = max(freq.keys())
    return tuple(freq.get(k, 0) for k in range(1, max_val + 1))


def enumerate_all_triples(c, n, max_weight=60):
    """
    Enumerate all (S, U, nu) triples in the ground set for Q_{n,c}.

    Returns: list of dicts with keys S, U, nu, sign, q_weight, freq, max_entry.
    """
    triples = []
    for j in range(n + 1):
        m = n - j
        for S in combinations(range(1, n + 1), j):
            S = frozenset(S)
            if m == 0:
                cpps = [((), (), (), 0)]
            else:
                cpps = enumerate_cpps_with_max(c, m, max_weight)
            for nu1, nu2, nu3, cpp_wt in cpps:
                for u_size in range(m + 1):
                    for U in combinations(range(1, m + 1), u_size):
                        U = frozenset(U)
                        total_wt = sum(S) + sum(U) + cpp_wt
                        sign = (-1) ** (len(S) + len(U))
                        freq = frequency_vector(S, U, nu1, nu2, nu3)
                        triples.append({
                            'S': S,
                            'U': U,
                            'nu': (nu1, nu2, nu3),
                            'sign': sign,
                            'q_weight': total_wt,
                            'freq': freq,
                            'max_entry': m,
                        })
    return triples
