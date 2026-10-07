"""Exhaustive small-profile check for the integer envelope proposition."""

from itertools import product


def envelope(n, alpha, beta):
    eta = [max(0, a + b - n) for a, b in zip(alpha, beta)]
    m = [min(a, b) for a, b in zip(alpha, beta)]
    M = [max(a, b) for a, b in zip(alpha, beta)]
    Q, E, H = sum(eta), sum(m), sum(M)
    Ta, Tb = sum(alpha), sum(beta)
    Cmax = len(alpha) * (n - 1)
    expected = {}
    for c in range(Q, Cmax + 1):
        if c <= E:
            expected[c] = len(alpha) * n - Ta - Tb + c
        elif c <= H:
            expected[c] = len(alpha) * n - H
        else:
            expected[c] = len(alpha) * n - c
    return eta, m, M, Q, E, H, Cmax, expected


def greedy(n, alpha, beta, c):
    eta, m, M, Q, E, H, Cmax, expected = envelope(n, alpha, beta)
    if not (Q <= c <= Cmax):
        return None
    y = eta[:]
    remaining = c - Q
    for lower, upper in ((eta, m), (m, M), (M, [n - 1] * len(alpha))):
        for i in range(len(alpha)):
            capacity = upper[i] - lower[i]
            take = min(remaining, capacity)
            y[i] += take
            remaining -= take
    assert remaining == 0
    return y


def check_profile(n, alpha, beta):
    eta, m, M, Q, E, H, Cmax, expected = envelope(n, alpha, beta)
    best = {}
    ranges = [range(e, n) for e in eta]
    for y in product(*ranges):
        c = sum(y)
        value = sum(max(a, yi) + max(b, yi) - yi
                    for a, b, yi in zip(alpha, beta, y))
        best[c] = min(best.get(c, 10**9), value)
    for c, lower_sum in expected.items():
        assert best[c] == n * len(alpha) - lower_sum
        y = greedy(n, alpha, beta, c)
        assert sum(y) == c
        for a, b, yi in zip(alpha, beta, y):
            k, l, t = max(a, yi), max(b, yi), max(b, yi) - yi
            assert 0 <= k < n and 0 <= l < n
            assert max(0, l - k) <= t <= min(l, n - k)


def main():
    profiles = 0
    for n in range(1, 5):
        for s in range(1, 4):
            for alpha in product(range(n), repeat=s):
                for beta in product(range(n), repeat=s):
                    check_profile(n, alpha, beta)
                    profiles += 1
    print(f"PASS: integer envelope and greedy checks for {profiles} profiles")


if __name__ == "__main__":
    main()
