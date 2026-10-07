"""Verify a five-point exact-distance segment of the q=5,s=4 frontier.

The fixed Fourier/GRS family realizes the certified plateau point K=10 with
exact distances (d_Z,d_X)=(3,3) at ebit counts c=0,1,2,3,4.
"""

from exact_distance_q5 import FOURIER, exact_pair


CASES = {
    0: ((0, 2, 1, 0), (2, 0, 0, 1), (2, 0, 0, 1)),
    1: ((0, 2, 1, 1), (2, 0, 0, 1), (2, 0, 0, 0)),
    2: ((0, 2, 1, 1), (2, 0, 1, 1), (2, 0, 0, 0)),
    3: ((0, 2, 1, 1), (2, 1, 1, 1), (2, 0, 0, 0)),
    4: ((1, 2, 1, 1), (2, 1, 1, 1), (1, 0, 0, 0)),
}


def main():
    for c, (k, ell, t) in CASES.items():
        observed_c = sum(ell_i - t_i for ell_i, t_i in zip(ell, t))
        observed_k = 16 - sum(k_i + t_i for k_i, t_i in zip(k, t))
        assert observed_c == c
        assert observed_k == 10
        distances = exact_pair(FOURIER, k, ell, t)[:2]
        assert distances == (3, 3), (c, distances)
        print(f"c={c}: K=10, exact (d_Z,d_X)={distances}")
    print("PASS: five-point q=5,s=4 exact witness coverage")


if __name__ == "__main__":
    main()
