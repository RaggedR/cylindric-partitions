#!/usr/bin/env python3
"""Littlewood-Richardson puzzle generator.

Given three partitions encoded as binary strings on the boundary of a triangle,
generates all valid LR puzzles using the Zinn-Justin integrable tiling model.

Encoding (following ZJ / Knutson-Tao):
  - Binary string of length n with k zeros and (n-k) ones.
  - '1' = step right, '0' = step up in a lattice path from (0,0) to (n-k, k).
  - The partition is the area below the path in the k x (n-k) box.
  - The cumulative x-coordinate at the y-th up step gives the width at height y.
  - Partition parts are in English convention: largest (top) first.

Boundary assignment for c^lambda_{mu,nu}:
  - Bottom edge (lambda): read left to right
  - Left edge (nu): read bottom to top (green particles exit left)
  - Right edge (mu^T): read top to bottom (red particles exit right)
"""


# ---------------------------------------------------------------------------
# Partition <-> binary string conversion
# ---------------------------------------------------------------------------

def binary_to_partition(s):
    """Convert binary string to partition (English convention: largest part first).

    The width at height y = cumulative x-coordinate at the y-th up step.
    Reverse for English convention.
    """
    widths = []  # widths from bottom to top
    x = 0
    for c in s:
        if c == '1':
            x += 1
        else:  # '0' = up step
            widths.append(x)
    # Reverse: English convention (top row = largest first)
    widths.reverse()
    while widths and widths[-1] == 0:
        widths.pop()
    return tuple(widths)


def partition_to_binary(partition, n, k):
    """Convert partition to binary string of length n with k zeros.

    Reverse partition to get widths bottom-to-top (non-decreasing),
    then emit incremental right steps at each height.
    """
    parts = list(partition) + [0] * k
    parts = parts[:k]
    # Widths bottom to top: reverse of English convention
    widths = parts[::-1]  # widths[0] = bottom (smallest), widths[k-1] = top (largest)
    s = []
    prev = 0
    for y in range(k):
        rights = widths[y] - prev
        s.extend(['1'] * rights)
        s.append('0')
        prev = widths[y]
    remaining = (n - k) - prev
    s.extend(['1'] * remaining)
    return ''.join(s)


def conjugate_partition(partition):
    """Compute the conjugate (transpose) of a partition."""
    if not partition:
        return ()
    cols = []
    for j in range(partition[0]):
        col_height = sum(1 for p in partition if p > j)
        cols.append(col_height)
    return tuple(cols)


# ---------------------------------------------------------------------------
# Row evolution (the core algorithm)
# ---------------------------------------------------------------------------

def evolve_row(row):
    """Evolve a row one step up the triangle.

    Returns list of (next_row, left_exit, right_exit) tuples.
    """
    m = len(row)
    if m == 0:
        return []
    if m == 1:
        val = row[0]
        if val == 1:
            return [([], False, True)]
        elif val == -1:
            return [([], True, False)]
        else:
            return [([], False, False)]

    n_slots = m - 1

    # Compute natural targets
    targets = [[] for _ in range(n_slots)]
    left_forced = False
    right_forced = False

    for j in range(m):
        val = row[j]
        if val == 1:  # red -> slot j
            if j < n_slots:
                targets[j].append((1, j))
            else:
                right_forced = True
        elif val == -1:  # green -> slot j-1
            if j > 0:
                targets[j - 1].append((-1, j))
            else:
                left_forced = True

    # Find collision slots
    collisions = []
    for slot in range(n_slots):
        has_red = any(t == 1 for t, _ in targets[slot])
        has_green = any(t == -1 for t, _ in targets[slot])
        if has_red and has_green:
            collisions.append(slot)

    n_collisions = len(collisions)
    results = []

    for choices in range(1 << n_collisions):
        next_row = [0] * n_slots
        extra_left = False
        extra_right = False
        occupied = set()

        # Place non-colliding particles
        for slot in range(n_slots):
            if slot not in collisions:
                for t, _ in targets[slot]:
                    next_row[slot] = t
                    occupied.add(slot)

        # Resolve collisions left to right
        for i, coll_slot in enumerate(collisions):
            green_displaced = not (choices & (1 << i))

            if green_displaced:
                next_row[coll_slot] = 1
                occupied.add(coll_slot)
                placed = False
                for s in range(coll_slot - 1, -1, -1):
                    if s not in occupied:
                        next_row[s] = -1
                        occupied.add(s)
                        placed = True
                        break
                if not placed:
                    extra_left = True
            else:
                next_row[coll_slot] = -1
                occupied.add(coll_slot)
                placed = False
                for s in range(coll_slot + 1, n_slots):
                    if s not in occupied:
                        next_row[s] = 1
                        occupied.add(s)
                        placed = True
                        break
                if not placed:
                    extra_right = True

        results.append((next_row, left_forced or extra_left, right_forced or extra_right))

    return results


# ---------------------------------------------------------------------------
# Puzzle generation
# ---------------------------------------------------------------------------

def generate_puzzles(bottom, target_left, target_right):
    """Generate all LR puzzles for given boundary conditions.

    Args:
        bottom: binary string for lambda (bottom edge, left to right)
        target_left: binary string for nu (left edge, bottom to top)
        target_right: binary string for mu^T (right edge, top to bottom)
    """
    n = len(bottom)
    assert len(target_left) == n
    assert len(target_right) == n

    initial_row = [1 if c == '1' else -1 for c in bottom]

    # Left: '0' at step r means green exits
    target_left_exits = [c == '0' for c in target_left]
    # Right: read top to bottom, reverse for step indexing (bottom to top)
    target_right_reversed = list(target_right)[::-1]
    target_right_exits = [c == '1' for c in target_right_reversed]

    puzzles = []

    def backtrack(row, step, rows_so_far, left_exits, right_exits):
        if len(row) == 0:
            if left_exits == target_left_exits and right_exits == target_right_exits:
                puzzles.append({
                    'rows': rows_so_far[:],
                    'left_exits': left_exits[:],
                    'right_exits': right_exits[:],
                })
            return

        for i in range(step):
            if left_exits[i] != target_left_exits[i]:
                return
            if right_exits[i] != target_right_exits[i]:
                return

        for next_row, left_exit, right_exit in evolve_row(row):
            backtrack(next_row, step + 1, rows_so_far + [next_row],
                      left_exits + [left_exit], right_exits + [right_exit])

    backtrack(initial_row, 0, [initial_row], [], [])
    return puzzles


def count_lr(bottom):
    """Count all puzzles for a given bottom edge, grouped by boundary.

    Returns dict mapping (left_partition, right_partition) -> count.
    """
    n = len(bottom)
    initial_row = [1 if c == '1' else -1 for c in bottom]
    counts = {}

    def backtrack(row, left_exits, right_exits):
        if len(row) == 0:
            left_str = ''.join('0' if e else '1' for e in left_exits)
            right_list = ['1' if e else '0' for e in right_exits][::-1]
            right_str = ''.join(right_list)
            nu = binary_to_partition(left_str)
            mu_t = binary_to_partition(right_str)
            mu = conjugate_partition(mu_t)
            key = (mu, nu)
            counts[key] = counts.get(key, 0) + 1
            return

        for next_row, left_exit, right_exit in evolve_row(row):
            backtrack(next_row, left_exits + [left_exit],
                      right_exits + [right_exit])

    backtrack(initial_row, [], [])
    return counts


# ---------------------------------------------------------------------------
# High-level interface
# ---------------------------------------------------------------------------

def lr_coefficient(lam, mu, nu, n=None, k=None):
    """Compute c^lambda_{mu, nu} directly from partitions.

    If n, k not given, chooses the smallest n, k that fit all three partitions.
    """
    if k is None:
        k = max(len(lam), len(mu), len(nu))
    if n is None:
        max_width = max(
            (lam[0] if lam else 0),
            (mu[0] if mu else 0),
            (nu[0] if nu else 0),
        )
        n = k + max_width

    bottom = partition_to_binary(lam, n, k)
    mu_t = conjugate_partition(mu)
    left = partition_to_binary(nu, n, k)
    right = partition_to_binary(mu_t, n, k)

    puzzles = generate_puzzles(bottom, left, right)
    return len(puzzles), puzzles


# ---------------------------------------------------------------------------
# TikZ output
# ---------------------------------------------------------------------------

import math

def puzzle_to_tikz(puzzle, n, lam=None, mu=None, nu=None):
    """Generate TikZ code for a puzzle.

    Draws colored dots on a triangular lattice: red for 1, green for -1.
    """
    rows = puzzle['rows']

    h = math.sqrt(3) / 2
    scale = 1.2

    lines = []
    lines.append(r"\begin{tikzpicture}[scale=%.1f]" % scale)
    lines.append(r"  \tikzstyle{reddot}=[circle, fill=red!80, inner sep=2.5pt]")
    lines.append(r"  \tikzstyle{greendot}=[circle, fill=green!60!black, inner sep=2.5pt]")

    def pos(row_idx, col_idx):
        x = col_idx + row_idx * 0.5
        y = row_idx * h
        return (x, y)

    # Draw dots: red for 1, green for -1, skip 0
    for r in range(len(rows)):
        for j, val in enumerate(rows[r]):
            x, y = pos(r, j)
            if val == 1:
                lines.append(r"  \node[reddot] at (%.3f, %.3f) {};" % (x, y))
            elif val == -1:
                lines.append(r"  \node[greendot] at (%.3f, %.3f) {};" % (x, y))

    # Partition labels
    if lam or mu or nu:
        mid_x = (n - 1) / 2.0
        top_x, top_y = pos(n - 1, 0)
        lines.append(r"  \node[below] at (%.3f, -0.3) {$\lambda = %s$};" % (mid_x, _partition_str(lam)))
        lines.append(r"  \node[left] at (-0.5, %.3f) {$\nu = %s$};" % (top_y / 2, _partition_str(nu)))
        right_mid_x = (n - 1) * 0.75 + 0.8
        lines.append(r"  \node[right] at (%.3f, %.3f) {$\mu = %s$};" % (
            right_mid_x, top_y / 2, _partition_str(mu)))

    lines.append(r"\end{tikzpicture}")
    return '\n'.join(lines)


def _partition_str(p):
    """Format a partition for LaTeX."""
    if not p:
        return r"\varnothing"
    return "(%s)" % ", ".join(str(x) for x in p)


def puzzles_to_tikz_document(puzzles, lam, mu, nu):
    """Generate a complete LaTeX document showing all puzzles."""
    lines = []
    lines.append(r"\documentclass[border=10pt]{standalone}")
    lines.append(r"\usepackage{tikz}")
    lines.append(r"\begin{document}")

    n = len(puzzles[0]['rows'][0]) if puzzles else 0
    count = len(puzzles)

    # Title
    lines.append(r"\begin{tabular}{c}")
    lines.append(r"$c^{%s}_{%s,\,%s} = %d$ \\[1em]" % (
        _partition_str(lam), _partition_str(mu), _partition_str(nu), count))

    for i, puzzle in enumerate(puzzles):
        lines.append(r"\textbf{Puzzle %d} \\[0.5em]" % (i + 1))
        # Only show boundary labels on first puzzle
        if i == 0:
            lines.append(puzzle_to_tikz(puzzle, n, lam, mu, nu))
        else:
            lines.append(puzzle_to_tikz(puzzle, n))
        if i < len(puzzles) - 1:
            lines.append(r"\\[2em]")
        lines.append("")

    lines.append(r"\end{tabular}")
    lines.append(r"\end{document}")
    return '\n'.join(lines)


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

def test_conversions():
    """Verify partition <-> binary string conversions."""
    cases = [
        # (partition, n, k, expected_string)
        ((2,), 3, 1, "110"),
        ((1,), 3, 1, "101"),
        ((), 3, 1, "011"),
        ((2, 1), 5, 2, "10101"),   # NOT "11010"!
        ((3, 1), 6, 2, "101101"),  # paper's +-++-+
        ((3, 2), 5, 2, "11010"),
        ((1, 1), 4, 2, "1001"),
        ((2,), 4, 2, "0110"),
    ]
    for part, n, k, expected in cases:
        got = binary_to_partition(expected)
        assert got == part, f"binary_to_partition('{expected}') = {got}, expected {part}"
        got_str = partition_to_binary(part, n, k)
        assert got_str == expected, f"partition_to_binary({part}, {n}, {k}) = '{got_str}', expected '{expected}'"

    # Test conjugate
    assert conjugate_partition((2, 1)) == (2, 1)
    assert conjugate_partition((3, 1)) == (2, 1, 1)
    assert conjugate_partition((2, 2)) == (2, 2)
    assert conjugate_partition((3,)) == (1, 1, 1)
    assert conjugate_partition(()) == ()

    print("All conversion tests passed!")


def test_evolve():
    """Test row evolution with Robin's example."""
    row = [0, 0, 1, 1, -1, 0]
    results = evolve_row(row)
    next_rows = [r for r, _, _ in results]
    print(f"Evolving {row}:")
    for r, le, re in results:
        print(f"  -> {r}  left_exit={le}  right_exit={re}")
    assert [0, -1, 1, 1, 0] in next_rows
    assert [0, 0, 1, -1, 1] in next_rows
    print("Evolution test passed!\n")


def test_basic():
    """Test c^{(2)}_{(1),(1)} = 1."""
    count, puzzles = lr_coefficient((2,), (1,), (1,))
    print(f"c^(2)_((1),(1)) = {count}")
    assert count == 1, f"Expected 1, got {count}"
    print("Basic test passed!\n")


def test_21():
    """Test all LR coefficients for lambda=(2,1)."""
    # Use n=5, k=2 (fits in 2x3 box)
    bottom = partition_to_binary((2, 1), 5, 2)
    print(f"lambda=(2,1), bottom='{bottom}'")
    counts = count_lr(bottom)
    print(f"All nonzero c^(2,1)_(mu, nu):")
    for (mu, nu), c in sorted(counts.items(), key=lambda x: str(x)):
        if c > 0:
            print(f"  c^(2,1)_({mu}, {nu}) = {c}")

    # Known values
    expected = {
        ((2,), (1,)): 1,
        ((1, 1), (1,)): 1,
        ((1,), (2,)): 1,
        ((1,), (1, 1)): 1,
        ((2, 1), ()): 1,
        ((), (2, 1)): 1,
    }
    for (mu, nu), exp_c in expected.items():
        got = counts.get((mu, nu), 0)
        status = "OK" if got == exp_c else "FAIL"
        print(f"  [{status}] c^(2,1)_({mu}, {nu}) = {got} (expected {exp_c})")
    print()


def test_c321():
    """Test c^{(3,2,1)}_{(2,1),(2,1)} = 2."""
    count, puzzles = lr_coefficient((3, 2, 1), (2, 1), (2, 1))
    assert count == 2, f"Expected 2, got {count}"
    print(f"c^(3,2,1)_((2,1),(2,1)) = {count}  (expected 2)")

    # Verify puzzles differ
    assert puzzles[0]['rows'] != puzzles[1]['rows'], "Two puzzles should be distinct"
    print("c^(3,2,1) test passed!\n")


def test_tikz():
    """Test TikZ generation compiles (structural check)."""
    lam, mu, nu = (3, 2, 1), (2, 1), (2, 1)
    count, puzzles = lr_coefficient(lam, mu, nu)

    # Single puzzle TikZ
    tikz = puzzle_to_tikz(puzzles[0], 6, lam, mu, nu)
    assert r"\begin{tikzpicture}" in tikz
    assert r"\end{tikzpicture}" in tikz
    assert r"\node[reddot]" in tikz
    assert r"\node[greendot]" in tikz
    assert r"\lambda = (3, 2, 1)" in tikz

    # Full document
    doc = puzzles_to_tikz_document(puzzles, lam, mu, nu)
    assert r"\documentclass" in doc
    assert r"c^{(3, 2, 1)}_{(2, 1)" in doc
    assert doc.count(r"\begin{tikzpicture}") == 2

    print("TikZ generation test passed!\n")


if __name__ == "__main__":
    test_conversions()
    test_evolve()
    test_basic()
    test_21()
    test_c321()
    test_tikz()
