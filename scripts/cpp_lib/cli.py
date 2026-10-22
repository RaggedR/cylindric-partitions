"""
cli.py — CLI argument parsing for investigation scripts.
"""


def parse_args(argv, defaults=None):
    """
    Parse common CLI arguments: n, c, max_weight.

    Usage: script [n] [c0,c1,c2] [max_weight]
    Returns: (n, c, max_weight)
    """
    if defaults is None:
        defaults = {'n': 1, 'c': (2, 1, 1), 'max_weight': 30}

    n = defaults['n']
    c = defaults['c']
    max_weight = defaults['max_weight']

    if len(argv) > 1:
        n = int(argv[1])
    if len(argv) > 2:
        parts = argv[2].split(',')
        c = tuple(int(x) for x in parts)
    if len(argv) > 3:
        max_weight = int(argv[3])

    return n, c, max_weight
