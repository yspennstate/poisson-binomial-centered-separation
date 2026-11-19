"""Section 4 (six summands): check, in exact rational arithmetic at random centered feasible pairs with n = 6,
    D_3 - D_4 = L_6 / 8   and   (D_3 + D_4) / 2 = (5 Delta - e_3(d) + e_5(d) - Q + H) / 16,
and the value b_6(g) = (30 g - 20 g^3 + 6 g^5) / 16. Standard library only; stops with an error if a check fails."""
from fractions import Fraction as F
from itertools import combinations
import random


def esym(v, k):
    tot = F(0)
    for I in combinations(range(len(v)), k):
        p = F(1)
        for i in I:
            p *= v[i]
        tot += p
    return tot


def pmf(p):
    dist = [F(1)]
    for pi in p:
        new = [F(0)] * (len(dist) + 1)
        for k, w in enumerate(dist):
            new[k] += w * (1 - pi)
            new[k + 1] += w * pi
        dist = new
    return dist


def tails(p, q, k):
    PX, PY = pmf(p), pmf(q)
    return sum(PX[k:]) - sum(PY[k:])


def forms(s, d):
    n = len(s)
    Delta = sum(d)
    Q = sum((d[i] - Delta / 2) * s[i] ** 2 for i in range(n))
    H = sum(d[i] * esym([s[j] for j in range(n) if j != i], 4) for i in range(n))
    for I in combinations(range(n), 3):
        pr = F(1)
        for i in I:
            pr *= d[i]
        H += pr * esym([s[j] for j in range(n) if j not in I], 2)
    sp = [s[i] + d[i] for i in range(n)]
    sm = [s[i] - d[i] for i in range(n)]
    O = lambda j: (esym(sp, j) - esym(sm, j)) / 2
    L6 = sum(s[i] * d[i] for i in range(n)) + O(4) - 5 * O(6)
    return Q, H, L6


def main(trials=200, seed=3):
    rng = random.Random(seed)
    for _ in range(trials):
        while True:
            s = [F(rng.randint(-40, 40), 100) for _ in range(5)]
            s.append(-sum(s))
            if all(abs(x) < 1 for x in s):
                break
        d = [F(rng.randint(0, int((1 - abs(x)) * 100)), 100) for x in s]
        p = [(1 + s[i] + d[i]) / 2 for i in range(6)]
        q = [(1 + s[i] - d[i]) / 2 for i in range(6)]
        D3, D4 = tails(p, q, 3), tails(p, q, 4)
        Q, H, L6 = forms(s, d)
        assert D3 - D4 == L6 / 8
        assert (D3 + D4) / 2 == (5 * sum(d) - esym(d, 3) + esym(d, 5) - Q + H) / 16
    for g in [F(k, 17) for k in range(18)]:
        p, q = [(1 + g) / 2] * 6, [(1 - g) / 2] * 6
        assert tails(p, q, 3) == (30 * g - 20 * g ** 3 + 6 * g ** 5) / 16
    print(f"six-summand identities hold at {trials} random centered pairs; b_6 formula holds at 18 values of g")


if __name__ == "__main__":
    main()
