"""Exact inverse ebit-budget check for the q=5 target catalogue.

For each target pair, this script computes the smallest exact ebit count at
which the catalogue reaches K >= 9. It also checks the matrix-independent
separate-cost lower bound and the catalogue gate deficit.
"""

N = 4
K_TARGET = 9
C_MAX = 12


def constituent_options():
    out = []
    for k in range(N):
        for ell in range(N):
            for t in range(max(0, ell - k), min(ell, N - k) + 1):
                out.append((k, ell, t))
    return out


def separate_cost(delta, s=4):
    return sum((delta + r - 1) // r - 1 for r in range(1, s + 1))


def frontier(s, rho, options, delta_z, delta_x):
    allowed = []
    for j in range(s):
        mp_pos = rho[j]
        z_min = (delta_z + (s - mp_pos - 1)) // (s - mp_pos) - 1
        x_min = (delta_x + mp_pos) // (mp_pos + 1) - 1
        allowed.append([item for item in options
                        if item[0] >= z_min and item[1] >= x_min])
    dp = {0: 0}
    for slot in allowed:
        nxt = {}
        for used, cost in dp.items():
            for k, ell, t in slot:
                c = used + ell - t
                if c <= C_MAX:
                    nxt[c] = min(nxt.get(c, 10**9), cost + k + t)
        dp = nxt
    return {c: s * N - cost for c, cost in dp.items()}


def catalogue_ceiling(delta_z, delta_x, options, catalogue):
    ceilings = []
    for s, rho_one_line in catalogue:
        rho = tuple(x - 1 for x in rho_one_line)
        curve = frontier(s, rho, options, delta_z, delta_x)
        if curve:
            ceilings.append(max(curve.values()))
    return max(ceilings, default=None)


def minimum_budget(delta_z, delta_x, options, catalogue):
    records = []
    for s, rho_one_line in catalogue:
        rho = tuple(x - 1 for x in rho_one_line)
        curve = frontier(s, rho, options, delta_z, delta_x)
        feasible = [c for c, k in curve.items() if k >= K_TARGET]
        if feasible:
            records.append((min(feasible), (s, rho_one_line)))
    return min(records) if records else (None, None)


def target_boundary(K_target, budget, options, catalogue):
    feasible = set()
    for delta_x in range(1, N + 1):
        for delta_z in range(1, N + 1):
            for s, rho_one_line in catalogue:
                rho = tuple(x - 1 for x in rho_one_line)
                curve = frontier(s, rho, options, delta_z, delta_x)
                if any(k >= K_target and c <= budget
                       for c, k in curve.items()):
                    feasible.add((delta_z, delta_x))
                    break
    return {
        pair for pair in feasible
        if not any(
            (z >= pair[0] and x >= pair[1] and (z, x) != pair)
            for z, x in feasible
        )
    }


def main():
    options = constituent_options()
    catalogue = [(2, (1, 2)), (2, (2, 1)), (4, (1, 4, 3, 2))]
    table = {}
    for delta_x in range(1, N + 1):
        for delta_z in range(1, N + 1):
            table[(delta_z, delta_x)] = minimum_budget(
                delta_z, delta_x, options, catalogue)[0]
    expected = {
        (1, 1): 0, (2, 1): 0, (3, 1): 0, (4, 1): 0,
        (1, 2): 0, (2, 2): 0, (3, 2): 0, (4, 2): 0,
        (1, 3): 0, (2, 3): 0, (3, 3): 0, (4, 3): 1,
        (1, 4): 0, (2, 4): 0, (3, 4): 1, (4, 4): None,
    }
    assert table == expected, (table, expected)
    for delta_x in range(1, N + 1):
        row = [table[(delta_z, delta_x)] for delta_z in range(1, N + 1)]
        print(f"delta_X={delta_x}: c_min(delta_Z=1..4)={row}")
    boundary0 = target_boundary(9, 0, options, catalogue)
    boundary1 = target_boundary(9, 1, options, catalogue)
    boundary12 = target_boundary(9, 12, options, catalogue)
    assert boundary0 == {(4, 2), (3, 3), (2, 4)}
    assert boundary1 == {(4, 3), (3, 4)}
    assert boundary12 == boundary1
    print("K=9 boundary at cmax=0:", sorted(boundary0))
    print("K=9 boundary at cmax=1:", sorted(boundary1))
    print("K=9 boundary at cmax=12:", sorted(boundary12))
    lower = max(0, K_TARGET - 4 * N +
                separate_cost(4) + separate_cost(3))
    assert lower == 1
    c43 = table[(4, 3)]
    c44 = table[(4, 4)]
    assert c43 == lower and c44 is None
    gate_44 = catalogue_ceiling(4, 4, options, catalogue)
    assert gate_44 == 8 and K_TARGET - gate_44 == 1
    print("matrix-independent lower bound for (4,3):", lower)
    print("catalogue minimum for (4,3):", c43)
    print("catalogue minimum for (4,4):", c44)
    print("(4,4) catalogue ceiling:", gate_44)
    print("(4,4) gate deficit:", K_TARGET - gate_44)
    print("PASS: exact inverse ebit-budget q=5 catalogue")


if __name__ == "__main__":
    main()


