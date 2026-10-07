"""Joint order/profile optimization for the q=5, n=4 catalogue.

The script searches the available type-I orders for s=2 and s=4 and uses a
component-wise dynamic program to maximize K at each exact ebit count under
the certified targets d_Z,d_X >= 3 across the complete admissible budget
range.  It is dependency-free and illustrates
the algorithm stated in Section 4 of the manuscript.
"""

from collections import defaultdict


Q = 5
N = 4
DELTA_Z = 3
DELTA_X = 3
C_MAX = 12


def constituent_options():
    """All MDS dimension/intersection triples for n=4, with k,l<n."""
    out = []
    for k in range(N):
        for ell in range(N):
            for t in range(max(0, ell - k), min(ell, N - k) + 1):
                out.append((k, ell, t))
    return out


def optimize(s, rho, options):
    """Return the best K for each exact c for one (s,rho) record."""
    # rho is zero-based.  Constituent slot j receives both forced thresholds
    # from the MP position rho(j), because rho is an involution.
    allowed = []
    for j in range(s):
        mp_pos = rho[j]
        z_min = (DELTA_Z + (s - mp_pos - 1)) // (s - mp_pos) - 1
        x_min = (DELTA_X + mp_pos) // (mp_pos + 1) - 1
        allowed.append([
            item for item in options
            if item[0] >= z_min and item[1] >= x_min
        ])

    # W_j(c) is the minimum sum of k_i+t_i after j slots.
    dp = {0: 0}
    for slot in allowed:
        nxt = {}
        for used, cost in dp.items():
            for k, ell, t in slot:
                c_i = ell - t
                new_c = used + c_i
                if new_c <= C_MAX:
                    new_cost = cost + k + t
                    nxt[new_c] = min(nxt.get(new_c, 10**9), new_cost)
        dp = nxt

    return {c: s * N - cost for c, cost in sorted(dp.items())}


def main():
    options = constituent_options()
    # Catalogue records use one-line permutation notation.  For s=4 the
    # Fourier matrix has rho=(1,4,3,2); s=2 records cover both involutions.
    catalogue = [
        (2, (1, 2)),
        (2, (2, 1)),
        (4, (1, 4, 3, 2)),
    ]
    best = defaultdict(lambda: (-1, None))
    for s, rho_one_line in catalogue:
        rho = tuple(x - 1 for x in rho_one_line)
        curve = optimize(s, rho, options)
        for c, K in curve.items():
            if K > best[c][0]:
                best[c] = (K, (s, rho_one_line))

    expected = {c: (10 if c <= 6 else 16 - c) for c in range(13)}
    assert {c: best[c][0] for c in range(13)} == expected
    for c in range(13):
        K, record = best[c]
        print(f"c={c}: K_max={K}, selected (s,rho)={record}")
    print("PASS: joint q=5 order/profile optimization")


if __name__ == "__main__":
    main()
