import itertools
n=8
def monos(D): return [m for d in range(D+1) for m in itertools.combinations_with_replacement(range(n),d)]
def key(m): return tuple(sorted(m))
