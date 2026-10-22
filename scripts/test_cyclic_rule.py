#!/usr/bin/env python3
"""
test_cyclic_rule.py — Test the cyclic descent injection rule for max-1 CPPs.

RULE: Given a max-1 CPP triple (l1,l2,l3) with slack s = (s1,s2,s3),
      increment component k where:
      k = argmin{s_{(k-2) mod 3} : k in {1,2,3}, s_k > s_{(k-2) mod 3}}

This is the "cyclic descent with min-predecessor" rule.
"""

from collections import defaultdict, Counter
import sys


def enumerate_triples(c, max_w):
    """Enumerate max-1 CPP triples (l1,l2,l3) by weight."""
    c0, c1, c2 = c
    by_weight = defaultdict(list)
    for w in range(1, max_w + 1):
        for l1 in range(w + 1):
            for l2 in range(w - l1 + 1):
                l3 = w - l1 - l2
                if l2 > l1 + c1 or l3 > l2 + c2 or l1 > l3 + c0:
                    continue
                by_weight[w].append((l1, l2, l3))
    return by_weight


def slack(c, triple):
    """Compute slack vector."""
    c0, c1, c2 = c
    l1, l2, l3 = triple
    return (l3 + c0 - l1, l1 + c1 - l2, l2 + c2 - l3)


def cyclic_rule_minprev(s):
    """
    Cyclic descent rule with min-predecessor tie-breaking.
    Find k with s_k > s_{(k-2)%3}, pick k with smallest s_{(k-2)%3}, then smallest k.
    """
    candidates = []
    for k in range(1, 4):
        prev_idx = (k - 2) % 3
        if s[k - 1] > s[prev_idx]:
            candidates.append((s[prev_idx], k))
    if not candidates:
        return None
    candidates.sort()
    return candidates[0][1]


def test_injection(c, max_w):
    """Test injection property only (no Q comparison). Returns (ok, fps_by_weight)."""
    d = sum(c)
    by_weight = enumerate_triples(c, max_w)

    rule = {}
    for s1 in range(d + 1):
        for s2 in range(d - s1 + 1):
            s3 = d - s1 - s2
            s = (s1, s2, s3)
            comp = cyclic_rule_minprev(s)
            if comp is not None:
                rule[s] = comp

    all_ok = True
    fps_by_weight = {1: by_weight.get(1, [])}
    errors = []

    for w in range(2, max_w + 1):
        prev = by_weight.get(w - 1, [])
        curr_set = set(by_weight.get(w, []))

        if not prev:
            fps_by_weight[w] = list(curr_set)
            continue

        images = []
        for p in prev:
            s = slack(c, p)
            comp = rule.get(s)
            if comp is None:
                errors.append(f"  w={w}: No rule for slack {s}")
                all_ok = False
                continue

            t = list(p)
            t[comp - 1] += 1
            t = tuple(t)

            c0, c1, c2 = c
            if t[1] > t[0] + c1 or t[2] > t[1] + c2 or t[0] > t[2] + c0:
                errors.append(f"  w={w}: Invalid {p}→{t}")
                all_ok = False
                continue

            images.append(t)

        if len(set(images)) != len(images):
            dupes = {t: cnt for t, cnt in Counter(images).items() if cnt > 1}
            errors.append(f"  w={w}: DUPLICATES: {dupes}")
            all_ok = False

        matched = set(images)
        fps_by_weight[w] = [t for t in by_weight.get(w, []) if t not in matched]

    return all_ok, fps_by_weight, errors


def main():
    fast = "--fast" in sys.argv
    max_w = 15 if fast else 30

    # All rank-3 profiles with d <= 7
    profiles = []
    for d in range(2, 8 if not fast else 6):
        for c0 in range(d + 1):
            for c1 in range(d - c0 + 1):
                c2 = d - c0 - c1
                if c0 + c1 + c2 != d:
                    continue
                profiles.append((c0, c1, c2))

    print(f"{'='*70}")
    print(f"Testing cyclic descent injection rule for max-1 CPPs")
    print(f"Rule: k = argmin{{s_{{k-2 mod 3}} : s_k > s_{{k-2 mod 3}}}}")
    print(f"Testing {len(profiles)} profiles, max_w={max_w}")
    print(f"{'='*70}")

    pass_count = 0
    fail_count = 0

    for c in profiles:
        d = sum(c)
        ok, fps, errors = test_injection(c, max_w)

        # Count stable fps
        by_weight = enumerate_triples(c, max_w)
        stable_w = max_w
        count_at_stable = len(by_weight.get(stable_w, []))
        count_prev = len(by_weight.get(stable_w - 1, []))
        stable_fps = len(fps.get(stable_w, []))

        if ok:
            pass_count += 1
            print(f"  c={c}, d={d}: PASS "
                  f"(stable: {count_prev}→{count_at_stable}, fps={stable_fps})")
        else:
            fail_count += 1
            print(f"  c={c}, d={d}: FAIL")
            for e in errors[:5]:
                print(f"    {e}")

    print(f"\n  Results: {pass_count} PASS, {fail_count} FAIL out of {len(profiles)}")

    # Now compare with Q_1 for key profiles
    if pass_count > 0:
        print(f"\n{'='*70}")
        print(f"Comparing injection fps with Q_1 coefficients")
        print(f"{'='*70}")

        from cpp_lib import compute_gk, compute_Q

        key_profiles = [(2, 1, 1), (2, 2, 1), (1, 1, 0)]
        if not fast:
            key_profiles += [(3, 2, 2), (4, 2, 1)]

        for c in key_profiles:
            d = sum(c)
            ok, fps, _ = test_injection(c, max_w)
            if not ok:
                continue

            gk = compute_gk(c, 6, max_w + 20, verbose=False)
            Q = compute_Q(c, 1, max_w, gk_dict=gk)

            print(f"\n  c={c}, d={d}:")
            print(f"  {'w':>3}  {'#triple':>8}  {'#fps':>6}  {'Q[w]':>6}  {'fps-Q':>6}")
            for w in range(1, min(max_w + 1, 20)):
                n_triple = len(enumerate_triples(c, w).get(w, []))
                n_fps = len(fps.get(w, []))
                q_w = int(Q[w]) if w < len(Q) else 0
                diff = n_fps - q_w
                marker = "" if diff == 0 else f" (diff={diff})"
                if n_fps > 0 or q_w > 0:
                    print(f"  {w:>3}  {n_triple:>8}  {n_fps:>6}  {q_w:>6}  {diff:>6}{marker}")

    # Show rule for c=(2,1,1)
    print(f"\n{'='*70}")
    print(f"Complete rule for c=(2,1,1), d=4:")
    print(f"{'='*70}")
    d = 4
    for s1 in range(d + 1):
        for s2 in range(d - s1 + 1):
            s3 = d - s1 - s2
            s = (s1, s2, s3)
            comp = cyclic_rule_minprev(s)
            if comp:
                # Show the slack effect
                new_s = list(s)
                new_s[comp - 1] -= 1
                new_s[(comp - 2) % 3] += 1
                print(f"  s={s} → incr l_{comp} → s'={tuple(new_s)}")


if __name__ == "__main__":
    main()
