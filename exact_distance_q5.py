"""Exact CSS distances for the q=5,s=4 example in AEAQEC_revised.tex.

This dependency-free script fixes explicit MDS constituent generators.  It checks
the prescribed intersection dimensions and enumerates all words of increasing
weight in the two relative CSS spaces
    U^perp \\ (V cap U^perp) and V^perp \\ (U cap V^perp).
Because the first word occurs at the reported weight, the output is an exact
distance certificate for these fixed constituents (not for every MDS choice).
"""
from itertools import combinations, product

P = 5
N = 4
S = 4
X = [1, 2, 3, 4]
FOURIER = [
    [1, 1, 1, 1],
    [1, 2, 4, 3],
    [1, 4, 1, 4],
    [1, 3, 4, 2],
]
IDENTITY = [[int(i == j) for j in range(S)] for i in range(S)]


def inv(a):
    return pow(a % P, -1, P)


def rank(M):
    M = [[x % P for x in row] for row in M]
    if not M:
        return 0
    rows, cols = len(M), len(M[0])
    pivot_row = 0
    for col in range(cols):
        pivot = next((r for r in range(pivot_row, rows) if M[r][col]), None)
        if pivot is None:
            continue
        M[pivot_row], M[pivot] = M[pivot], M[pivot_row]
        z = inv(M[pivot_row][col])
        M[pivot_row] = [(x * z) % P for x in M[pivot_row]]
        for r in range(rows):
            if r != pivot_row and M[r][col]:
                z = M[r][col]
                M[r] = [(u - z * v) % P for u, v in zip(M[r], M[pivot_row])]
        pivot_row += 1
    return pivot_row


def nullspace(M):
    """A row basis for the nullspace of M over F_5."""
    M = [[x % P for x in row] for row in M]
    cols = len(M[0])
    pivots = []
    pivot_row = 0
    for col in range(cols):
        pivot = next((r for r in range(pivot_row, len(M)) if M[r][col]), None)
        if pivot is None:
            continue
        M[pivot_row], M[pivot] = M[pivot], M[pivot_row]
        z = inv(M[pivot_row][col])
        M[pivot_row] = [(x * z) % P for x in M[pivot_row]]
        for r in range(len(M)):
            if r != pivot_row and M[r][col]:
                z = M[r][col]
                M[r] = [(u - z * v) % P for u, v in zip(M[r], M[pivot_row])]
        pivots.append(col)
        pivot_row += 1
    free = [c for c in range(cols) if c not in pivots]
    basis = []
    for f in free:
        v = [0] * cols
        v[f] = 1
        for r, c in enumerate(pivots):
            v[c] = (-M[r][f]) % P
        basis.append(v)
    return basis


def grs(dim, multipliers):
    return [[multipliers[j] * pow(X[j], r, P) % P for j in range(N)]
            for r in range(dim)]


def in_span(v, basis):
    return rank(basis) == rank(basis + [v])


def intersection_dimension(A, B):
    return len(A) + len(B) - rank(A + B)


def mp_generator(constituents, matrix):
    """Generator rows for [C_1,...,C_s] matrix, with contiguous blocks."""
    out = []
    for i, basis in enumerate(constituents):
        for c in basis:
            word = [0] * (S * N)
            for pos in range(N):
                for block in range(S):
                    word[block * N + pos] = (
                        word[block * N + pos] + c[pos] * matrix[i][block]
                    ) % P
            out.append(word)
    return out


def first_relative_word(G, excluded, max_weight):
    r"""Return the first word in span(G)\span(excluded), by Hamming weight."""
    length = len(G[0])
    for weight in range(1, max_weight + 1):
        for support in combinations(range(length), weight):
            for values in product(range(1, P), repeat=weight):
                word = [0] * length
                for pos, value in zip(support, values):
                    word[pos] = value
                if in_span(word, G) and not in_span(word, excluded):
                    return word
    return None


def exact_pair(matrix, k, ell, target_t, max_weight=3):
    # C_i^perp=E_i are [4,4-k_i] RS codes.  D_i are GRS codes with
    # multiplier (1,1,1,2), giving the required t_i values.
    C, D = [], []
    for ki, li, ti in zip(k, ell, target_t):
        E = grs(N - ki, [1] * N)
        Di = grs(li, [1, 1, 1, 2])
        assert intersection_dimension(E, Di) == ti
        C.append(nullspace(E))
        D.append(Di)
    U = mp_generator(C, matrix)
    V = mp_generator(D, matrix)
    U_perp, V_perp = nullspace(U), nullspace(V)
    z = first_relative_word(U_perp, V, max_weight)
    x = first_relative_word(V_perp, U, max_weight)
    assert z is not None and x is not None
    return len([x for x in z if x]), len([x for x in x if x]), z, x


def main():
    # The first row of the manuscript: [[16,10,3/3;4]]_5 versus its
    # identity-map baseline [[16,10,2/2;4]]_5.
    k0, ell0, t0 = [1, 2, 1, 1], [2, 1, 1, 1], [1, 0, 0, 0]
    twisted = exact_pair(FOURIER, k0, ell0, t0)
    untwisted = exact_pair(IDENTITY, k0, ell0, t0)
    assert twisted[:2] == (3, 3)
    assert untwisted[:2] == (2, 2)
    print("baseline twisted exact (d_Z,d_X) =", twisted[:2], "witnesses:", twisted[2], twisted[3])
    print("baseline identity exact (d_Z,d_X) =", untwisted[:2], "witnesses:", untwisted[2], untwisted[3])

    # A Pareto point found by finite order/constituent search.  Zero-dimensional
    # constituents are allowed by the manuscript's conventions and produce c=0
    # while preserving K=10 and the same exact distances.
    k1, ell1, t1 = [0, 2, 1, 0], [2, 0, 0, 1], [2, 0, 0, 1]
    zero_ebit = exact_pair(FOURIER, k1, ell1, t1)
    assert zero_ebit[:2] == (3, 3)
    print("zero-ebit twisted exact (d_Z,d_X) =", zero_ebit[:2], "witnesses:", zero_ebit[2], zero_ebit[3])
    print("PASS: q=5,s=4 fixed-constituent and zero-ebit exact-distance checks")


if __name__ == "__main__":
    main()
