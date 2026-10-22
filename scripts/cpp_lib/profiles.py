"""
profiles.py — Profile utilities for rank-3 cylindric partitions.

A profile c = (c0, c1, c2) with d = c0+c1+c2 determines the interlacing
conditions for the multipartition representation. Pure Python, no Sage.
"""


def profile_to_mu(c):
    """mu(c) = (c1+c2+...+c_{r-1}, c2+...+c_{r-1}, ..., c_{r-1}, 0)."""
    r = len(c)
    return tuple(sum(c[j] for j in range(i + 1, r)) for i in range(r))


def all_profiles(d):
    """Generate all rank-3 profiles (c0, c1, c2) with c0+c1+c2 = d."""
    profiles = []
    for c0 in range(d + 1):
        for c1 in range(d - c0 + 1):
            c2 = d - c0 - c1
            profiles.append((c0, c1, c2))
    return profiles


def representative_profiles():
    """
    One profile per cyclic orbit, organized by d.

    Cyclic symmetry: GK_{(c0,c1,c2)} = GK_{(c2,c0,c1)} = GK_{(c1,c2,c0)}.
    We pick one representative per orbit.
    """
    return {
        1: [(1, 0, 0)],
        2: [(2, 0, 0), (1, 1, 0)],
        4: [(2, 1, 1), (1, 2, 1)],
        5: [(2, 2, 1), (3, 1, 1)],
        7: [(3, 2, 2), (2, 3, 2), (4, 2, 1)],
        8: [(3, 3, 2), (4, 3, 1)],
    }
