#!/usr/bin/env python3
"""
analyze_qsym_cancellation.py — Study how (S, U, ν) triples cancel
in the QSym refinement.

KEY QUESTION: Can we find a sign-reversing involution on (S, U, ν)
triples that preserves the QSym monomial (frequency vector)?

The answer turns out to be NO — some frequency classes have net
negative coefficient, so cancellation MUST cross between different
QSym classes. The q-specialization z_k → q^k merges different
frequency vectors into the same q-power, enabling cross-cancellation.

This script analyzes:
1. Which frequency vectors cancel internally (same QSym class)
2. Which require cross-cancellation (different QSym → same q-power)
3. The structure of the cross-cancellation patterns
"""

from collections import defaultdict
from cpp_lib import (enumerate_all_triples, frequency_vector,
                     q_specialization, parse_args)
import sys


def analyze_cross_cancellation(triples, n, c):
    """
    Main analysis: group by q-weight, then by frequency vector.
    Identify internal vs cross-cancellation.
    """
    d = sum(c)
    print(f"Profile c = {c}, d = {d}, n = {n}")
    print(f"Total triples: {len(triples)}")

    # Group by q-weight
    by_qweight = defaultdict(list)
    for t in triples:
        by_qweight[t['q_weight']].append(t)

    # Q-polynomial
    Q_coeffs = {}
    for w, items in by_qweight.items():
        net = sum(t['sign'] for t in items)
        Q_coeffs[w] = net

    max_w = max(Q_coeffs.keys()) if Q_coeffs else 0
    terms = []
    for w in range(max_w + 1):
        c_w = Q_coeffs.get(w, 0)
        if c_w != 0:
            terms.append(f"{c_w}q^{w}" if w > 0 else str(c_w))
    print(f"Q_{n} = {' + '.join(terms) if terms else '0'}")

    # Now analyze cancellation structure at each q-weight
    print(f"\n{'='*70}")
    print(f"CANCELLATION STRUCTURE BY q-WEIGHT")
    print(f"{'='*70}")

    total_internal = 0
    total_cross = 0

    for w in sorted(by_qweight.keys()):
        items = by_qweight[w]
        if not items:
            continue

        total_pos = sum(1 for t in items if t['sign'] == +1)
        total_neg = sum(1 for t in items if t['sign'] == -1)
        net = total_pos - total_neg

        # Group by frequency vector within this q-weight
        by_freq = defaultdict(list)
        for t in items:
            by_freq[t['freq']].append(t)

        # Internal cancellation: positive and negative within same freq
        internal_cancel = 0
        cross_excess = defaultdict(int)  # freq -> excess after internal cancel

        for freq, fitems in by_freq.items():
            fp = sum(1 for t in fitems if t['sign'] == +1)
            fn = sum(1 for t in fitems if t['sign'] == -1)
            internal_cancel += min(fp, fn)
            cross_excess[freq] = fp - fn

        # Cross-cancellation: excess positive from one freq cancels excess
        # negative from another
        pos_excess = sum(v for v in cross_excess.values() if v > 0)
        neg_excess = sum(-v for v in cross_excess.values() if v < 0)
        cross_cancel = min(pos_excess - net if net > 0 else pos_excess,
                          neg_excess)

        total_internal += internal_cancel
        total_cross += cross_cancel

        # Only show weights with cross-cancellation
        if cross_cancel > 0 or w <= 8:
            print(f"\nq^{w}: {total_pos}+ {total_neg}− net={net}")
            # Show freq breakdown
            for freq in sorted(by_freq.keys()):
                fitems = by_freq[freq]
                fp = sum(1 for t in fitems if t['sign'] == +1)
                fn = sum(1 for t in fitems if t['sign'] == -1)
                fnet = fp - fn
                freq_str = f"z^{freq}"
                ic = min(fp, fn)
                mark = "" if fnet == 0 else f"  ← excess {fnet:+d}"
                if cross_cancel > 0 and fnet != 0:
                    mark += " (CROSS)"
                print(f"  {freq_str:20s}  {fp}+ {fn}−  "
                      f"internal={ic} net={fnet}{mark}")

            if cross_cancel > 0:
                print(f"  → Internal: {internal_cancel}, "
                      f"Cross: {cross_cancel}")

                # Show which freqs pair in cross-cancellation
                pos_freqs = {f: v for f, v in cross_excess.items() if v > 0}
                neg_freqs = {f: v for f, v in cross_excess.items() if v < 0}
                print(f"  → Pos excess: {dict(pos_freqs)}")
                print(f"  → Neg excess: {dict(neg_freqs)}")

                # Check: do the crossing freqs map to same q-weight?
                for f in list(pos_freqs) + list(neg_freqs):
                    qw = sum((k+1) * f[k] for k in range(len(f)))
                    assert qw == w, f"Freq {f} has q-weight {qw} != {w}"

    print(f"\n{'='*70}")
    print(f"SUMMARY")
    print(f"{'='*70}")
    print(f"Total internal cancellations: {total_internal}")
    print(f"Total cross-cancellations: {total_cross}")
    print(f"Cross fraction: {total_cross}/{total_internal + total_cross}"
          f" = {total_cross/(total_internal + total_cross):.1%}"
          if total_internal + total_cross > 0 else "")


def analyze_freq_structure(triples, n, c):
    """
    Analyze the frequency vectors that need cross-cancellation.
    Which compositions map to the same q-weight?
    """
    print(f"\n{'='*70}")
    print(f"FREQUENCY VECTOR → q-WEIGHT MAP")
    print(f"{'='*70}")

    # Group triples by frequency vector
    by_freq = defaultdict(list)
    for t in triples:
        by_freq[t['freq']].append(t)

    # For each freq, compute its q-weight and net coefficient
    freq_data = {}
    for freq, items in by_freq.items():
        qw = sum((k+1) * freq[k] for k in range(len(freq)))
        net = sum(t['sign'] for t in items)
        freq_data[freq] = {'q_weight': qw, 'net': net, 'count': len(items)}

    # Group freqs by q-weight
    qw_to_freqs = defaultdict(list)
    for freq, data in freq_data.items():
        qw_to_freqs[data['q_weight']].append((freq, data))

    # Show partitions of the same q-weight
    for qw in sorted(qw_to_freqs.keys()):
        freqs = qw_to_freqs[qw]
        if len(freqs) <= 1:
            continue
        # Check if any cross-cancellation needed
        nets = [d['net'] for _, d in freqs]
        if all(n == 0 for n in nets):
            continue

        print(f"\nq^{qw}: {len(freqs)} frequency vectors")
        total_net = sum(d['net'] for _, d in freqs)
        for freq, data in sorted(freqs, key=lambda x: x[0]):
            print(f"  {str(freq):20s}  net={data['net']:+3d}  "
                  f"({data['count']} triples)")
        print(f"  Total net: {total_net}")

    # COMPOSITION PAIRS: which compositions share the same q-weight?
    print(f"\n{'='*70}")
    print(f"COMPOSITION PAIRS (same q-weight, different QSym monomial)")
    print(f"{'='*70}")
    print(f"These pairs of compositions α, β have Σ i·α_i = Σ i·β_i")
    print(f"and their net (S,U,ν) contributions must cross-cancel.")
    print()

    for qw in sorted(qw_to_freqs.keys()):
        freqs = qw_to_freqs[qw]
        if len(freqs) <= 1:
            continue
        nets = [d['net'] for _, d in freqs]
        if all(n == 0 for n in nets):
            continue

        pos_freqs = [(f, d) for f, d in freqs if d['net'] > 0]
        neg_freqs = [(f, d) for f, d in freqs if d['net'] < 0]

        if pos_freqs and neg_freqs:
            for pf, pd in pos_freqs:
                for nf, nd in neg_freqs:
                    print(f"  q^{qw}: z^{pf} (net={pd['net']:+d}) ↔ "
                          f"z^{nf} (net={nd['net']:+d})")
                    changes = []
                    max_len = max(len(pf), len(nf))
                    for k in range(max_len):
                        pk = pf[k] if k < len(pf) else 0
                        nk = nf[k] if k < len(nf) else 0
                        if pk != nk:
                            changes.append(f"  val {k+1}: {pk} → {nk}")
                    if changes:
                        print("    " + ", ".join(changes))


def analyze_S_U_patterns(triples, n, c):
    """
    Analyze the (S, U) patterns that appear at each frequency vector.
    """
    print(f"\n{'='*70}")
    print(f"(S, U) PATTERNS AT KEY FREQUENCY VECTORS")
    print(f"{'='*70}")

    # Group by frequency
    by_freq = defaultdict(list)
    for t in triples:
        by_freq[t['freq']].append(t)

    # Show patterns at interesting frequency vectors
    for freq in sorted(by_freq.keys()):
        items = by_freq[freq]
        net = sum(t['sign'] for t in items)
        if net == 0:
            continue  # Only show non-cancelling ones

        print(f"\nz^{freq}  (net = {net:+d})")
        for t in sorted(items, key=lambda t: (-t['sign'], len(t['S']),
                                               len(t['U']))):
            sign_str = "+" if t['sign'] == +1 else "−"
            S_str = "{" + ",".join(str(s) for s in sorted(t['S'])) + "}" \
                    if t['S'] else "∅"
            U_str = "{" + ",".join(str(u) for u in sorted(t['U'])) + "}" \
                    if t['U'] else "∅"
            nu_str = str(t['nu'])
            print(f"  {sign_str} S={S_str:8s} U={U_str:8s} "
                  f"ν={nu_str}")


def main():
    n, c, max_weight = parse_args(sys.argv,
                                  defaults={'n': 2, 'c': (1, 1, 0), 'max_weight': 40})

    d = sum(c)
    print(f"Analyzing cross-cancellation: c={c}, d={d}, n={n}")
    print(f"{'='*70}")

    triples = enumerate_all_triples(c, n, max_weight)
    analyze_cross_cancellation(triples, n, c)
    analyze_freq_structure(triples, n, c)
    analyze_S_U_patterns(triples, n, c)


if __name__ == '__main__':
    main()
