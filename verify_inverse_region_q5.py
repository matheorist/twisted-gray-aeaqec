"""Verify the inverse distance--ebit--dimension table in the manuscript."""

from math import ceil


def profile(delta_z, delta_x, n=4, s=4, rho=(1, 4, 3, 2)):
    a = [0] * (s + 1)
    b = [0] * (s + 1)
    for i in range(1, s + 1):
        a[rho[i - 1]] = ceil(delta_z / (s - i + 1)) - 1
        b[rho[s - i]] = ceil(delta_x / (s - i + 1)) - 1
    tz, tx = sum(a[1:]), sum(b[1:])
    h = sum(max(a[i], b[i]) for i in range(1, s + 1))
    e = tz + tx - h
    q = sum(max(0, a[i] + b[i] - n) for i in range(1, s + 1))
    return a, b, tz, tx, h, e, q


def main():
    n = s = 4
    k0, cmax = 10, 6
    expected = {1: 4, 2: 4, 3: 3, 4: 2}
    observed = {}
    for delta_x in range(1, n + 1):
        best = 0
        for delta_z in range(1, n + 1):
            a, b, tz, tx, h, e, q = profile(delta_z, delta_x, n, s)
            if max(a[1:] + b[1:]) >= n or q > cmax:
                continue
            budget_ceiling = s * n - tz - tx + min(cmax, e)
            if budget_ceiling >= k0:
                best = delta_z
        observed[delta_x] = best
    assert observed == expected, (observed, expected)

    # The symmetric (3,3) target has a zero-ebit ceiling point.
    _, _, tz, tx, h, e, q = profile(3, 3, n, s)
    assert (tz, tx, h, e, q) == (3, 3, 6, 0, 0)
    assert s * n - tz - tx + min(cmax, e) == k0
    print("PASS: inverse q=5,s=4 target table", observed)


if __name__ == "__main__":
    main()
