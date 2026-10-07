"""Exhaustive check of the q=5, s=4 certified fixed-ebit frontier.

The script enumerates all MDS dimension/intersection profiles with n=4,
the Fourier-induced order rho=(1,4,3,2), and certified targets dZ,dX >= 3.
It compares the largest K at each exact ebit count with the closed formula
K_max(c)=10 for 0<=c<=6 and K_max(c)=16-c for 6<=c<=12.
"""

from itertools import product

n = 4
s = 4
rho = (0, 3, 2, 1)  # one-based (1,4,3,2)

options = []
for k in range(n):
    for ell in range(n):
        for t in range(max(0, ell - k), min(ell, n - k) + 1):
            options.append((k, ell, t))

best = {}
for profile in product(options, repeat=s):
    ks = tuple(x[0] for x in profile)
    ls = tuple(x[1] for x in profile)
    ts = tuple(x[2] for x in profile)
    mz = min((s - i) * (ks[rho[i]] + 1) for i in range(s))
    mx = min((s - i) * (ls[rho[s - 1 - i]] + 1) for i in range(s))
    if mz < 3 or mx < 3:
        continue
    c = sum(ls[i] - ts[i] for i in range(s))
    K = s * n - sum(ks[i] + ts[i] for i in range(s))
    best[c] = max(best.get(c, -1), K)

expected = {c: (10 if c <= 6 else 16 - c) for c in range(13)}
assert best == expected, (best, expected)
print("PASS: q=5,s=4 exhaustive certified frontier", best)
