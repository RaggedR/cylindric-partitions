#!/usr/bin/env python3
"""Test Q_{n,c}(q) for d ≡ 0 mod 3 profiles (d=3, d=6).
These are excluded from Warnaar's Conjecture 2.7.
We want to see: does Q actually go negative?"""

from math import gcd
import time

from cpp_lib.enumeration import compute_gk
from cpp_lib.qseries import compute_Q

# Profiles with d ≡ 0 mod 3
test_profiles = {
    3: [(1, 1, 1), (2, 1, 0), (3, 0, 0)],
    6: [(2, 2, 2), (3, 2, 1), (4, 1, 1)],
}

for d, profiles in sorted(test_profiles.items()):
    ell = gcd(d, 3)
    print(f"\n{'='*60}")
    print(f"d = {d}, ell = gcd({d},3) = {ell}")
    print(f"{'='*60}")

    for c in profiles:
        print(f"\n  Profile c = {c}")
        max_n_range = 3
        max_w = 20

        gk_max_n = max_n_range + 2
        gk_max_w = max_w + max_n_range * (max_n_range + 1) // 2 + 3

        t0 = time.time()
        gk = compute_gk(c, gk_max_n, gk_max_w, verbose=False)
        print(f"    GK computed in {time.time()-t0:.2f}s")

        for n in range(max_n_range + 1):
            Q = compute_Q(c, n, max_w, gk_dict=gk)

            # Check for negatives
            neg_coeffs = [(w, Q[w]) for w in range(max_w+1) if Q[w] < 0]
            q_at_1 = sum(Q)

            terms = []
            for w in range(max_w+1):
                if Q[w] != 0:
                    terms.append(f"{Q[w]}*q^{w}")
            q_str = " + ".join(terms) if terms else "0"
            if len(q_str) > 120:
                q_str = q_str[:120] + "..."

            status = "NEGATIVE!" if neg_coeffs else "nonneg"
            print(f"    n={n}: Q = {q_str}")
            print(f"         [{status}] Q(1) = {q_at_1}")
            if neg_coeffs:
                print(f"         Negative at: {neg_coeffs}")
