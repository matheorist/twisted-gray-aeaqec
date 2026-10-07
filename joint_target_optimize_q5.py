"""Target-aware order/profile optimization for the q=5, n=4 catalogue.

The script runs the fixed-target dynamic program over a finite grid of
(delta_Z, delta_X) pairs and reports the largest certified delta_Z attainable
for each delta_X at K >= 10 and c <= 6.  It is dependency-free.
"""

from collections import defaultdict


N = 4
C_MAX = 6
K_TARGET = 10


def constituent_options():
    """All MDS dimension/intersection triples for n=4, with k,l<n."""
    out = []
    for k in range(N):
        for ell in range(N):
            for t in range(max(0, ell - k), min(ell, N - k) + 1):
                out.append((k, ell, t))
    return out


def separate_cost(delta, s=4):
    """Unavoidable separate constituent-dimension cost T_s(delta)."""
    return sum((delta + r - 1) // r - 1 for r in range(1, s + 1))


def optimize(s, rho, options, delta_z, delta_x):
    """Return K_max over exact ebit counts c <= C_MAX for one record."""
    allowed = []
    for j in range(s):
        mp_pos = rho[j]
        z_min = (delta_z + (s - mp_pos - 1)) // (s - mp_pos) - 1
        x_min = (delta_x + mp_pos) // (mp_pos + 1) - 1
        allowed.append([
            item for item in options
            if item[0] >= z_min and item[1] >= x_min
        ])

    dp = {0: 0}
    for slot in allowed:
        nxt = {}
        for used, cost in dp.items():
            for k, ell, t in slot:
                new_c = used + ell - t
                if new_c <= C_MAX:
                    new_cost = cost + k + t
                    nxt[new_c] = min(nxt.get(new_c, 10**9), new_cost)
        dp = nxt
    return max((s * N - cost for c, cost in dp.items() if c <= C_MAX),
               default=-1)


def main():
    options = constituent_options()
    catalogue = [(2, (1, 2)), (2, (2, 1)), (4, (1, 4, 3, 2))]
    observed = {}
    for delta_x in range(1, N + 1):
        best_delta_z = 0
        for delta_z in range(1, N + 1):
            best_k = max(
                optimize(s, tuple(x - 1 for x in rho), options,
                         delta_z, delta_x)
                for s, rho in catalogue
            )
            if best_k >= K_TARGET:
                best_delta_z = delta_z
        observed[delta_x] = best_delta_z

    expected = {1: 4, 2: 4, 3: 3, 4: 2}
    assert observed == expected, (observed, expected)
    feasible = set()
    for delta_x in range(1, N + 1):
        for delta_z in range(1, N + 1):
            best_k = max(
                optimize(s, tuple(x - 1 for x in rho), options,
                         delta_z, delta_x)
                for s, rho in catalogue
            )
            if best_k >= K_TARGET:
                feasible.add((delta_z, delta_x))

    # Feasibility is downward closed in both certified distance targets.
    for delta_z, delta_x in feasible:
        for lower_z in range(1, delta_z + 1):
            for lower_x in range(1, delta_x + 1):
                assert (lower_z, lower_x) in feasible
    boundary = {
        pair for pair in feasible
        if not any(
            (z >= pair[0] and x >= pair[1] and (z, x) != pair)
            for z, x in feasible
        )
    }
    assert boundary == {(4, 2), (3, 3), (2, 4)}, boundary
    print("feasible target pairs on 1..4 grid:", len(feasible))
    print("target antichain boundary:", sorted(boundary))

    # Matrix-independent prefilter from Corollary cor:targetprefilter.
    threshold = 4 * N + C_MAX - K_TARGET
    wide_grid = range(1, 2 * N + 1)
    pruned = [
        (delta_z, delta_x)
        for delta_z in wide_grid
        for delta_x in wide_grid
        if separate_cost(delta_z) + separate_cost(delta_x) > threshold
    ]
    assert (7, 2) in pruned
    assert (4, 4) not in pruned
    print("max delta_Z at K>=10,c<=6:", observed)
    print("prefilter threshold:", threshold)
    print("prefilter pruned (1..8)^2 target pairs:", len(pruned))
    print("PASS: matrix-independent target prefilter")
    print("PASS: target-aware q=5 order/profile optimization")


if __name__ == "__main__":
    main()

