"""Joint order/profile optimization for the q=5, n=4 catalogue."""

from collections import defaultdict

Q = 5
N = 4
DELTA_Z = 3
DELTA_X = 3
C_MAX = 6


def constituent_options():
    out = []
    for k in range(N):
        for ell in range(N):
            for t in range(max(0, ell - k), min(ell, N - k) + 1):
                out.append((k, ell, t))
    return out


def optimize(s, rho, options):
    allowed = []
    for i in range(s):
        z_min = (DELTA_Z + (s - i - 1)) // (s - i) - 1
        x_min = (DELTA_X + i) // (i + 1) - 1
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
                    nxt[new_c] = min(nxt.get(new_c, 10**9), cost + k + t)
        dp = nxt
    return {c: s * N - cost for c, cost in sorted(dp.items())}


def main():
    options = constituent_options()
    catalogue = [(2, (1, 2)), (2, (2, 1)), (4, (1, 4, 3, 2))]
    best = defaultdict(lambda: (-1, None))
    for s, rho_one_line in catalogue:
        rho = tuple(x - 1 for x in rho_one_line)
        for c, K in optimize(s, rho, options).items():
            if K > best[c][0]:
                best[c] = (K, (s, rho_one_line))

    expected = {c: 10 for c in range(7)}
    assert {c: best[c][0] for c in range(7)} == expected
    for c in range(7):
        print(f"c={c}: K_max={best[c][0]}, selected (s,rho)={best[c][1]}")
    print("PASS: joint q=5 order/profile optimization")


if __name__ == "__main__":
    main()
