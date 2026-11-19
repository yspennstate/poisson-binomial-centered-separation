"""Section 4 (six summands): exact checks of the three examples, and a floating-point search for the inequality
Q - H <= |L_6| over centered feasible pairs with n = 6.
(1) SymPy: the formulas for Q - H and L_6 on the family s = (beta, -beta, 0, 0, 0, 0), and on the family
    s = (lam, -lam/5, ..., -lam/5), d = (mu, 1 - lam/5, ..., 1 - lam/5), with |L_6| / (Q - H) -> 2 as lam -> 0.
(2) Exact rational arithmetic: the pair s = (-3/4, 1/2, 2/5, 2/5, -19/20, 2/5), d = (0, 0, 0, 0, 1/20, 0) has
    largest s_i^2 / S_2 = 361/878, Q - H = 1/20000 and |L_6| = 167/2500.
(3) The search draws 25,000 points in five equal classes. Classes one to three: s uniform with sum zero and
    |s_i| < 1, d_i = (1 - |s_i|) u_i with u_i uniform (one), with u_i in {0, 1, uniform} (two), and with s then
    scaled by a factor in [0, 0.1] (three). Class four: one dominant coordinate s_1 = x, the others negative,
    d_1 in [(1-x)/2, 1-x] and the other gaps small. Class five: s with no s_i^2 above S_2/2 that has a vertex
    d_i in {0, 1 - |s_i|} with Q - H > 0, and d drawn near such a vertex or uniformly in the box.
Evidence, not a proof. Needs SymPy; stops with an error if a check fails."""
from fractions import Fraction as F
from itertools import combinations, product
import random
import sympy as sp


def esym(v, k):
    return sp.Add(*[sp.Mul(*[v[i] for i in I]) for I in combinations(range(len(v)), k)])


def esym_f(v, k):
    tot = 0
    for I in combinations(range(len(v)), k):
        p = 1
        for i in I:
            p *= v[i]
        tot += p
    return tot


def forms(s, d, e):
    n = len(s)
    Delta = sum(d)
    Q = sum((d[i] - Delta / 2) * s[i] ** 2 for i in range(n))
    H = sum(d[i] * e([s[j] for j in range(n) if j != i], 4) for i in range(n))
    for I in combinations(range(n), 3):
        pr = 1
        for i in I:
            pr = pr * d[i]
        H = H + pr * e([s[j] for j in range(n) if j not in I], 2)
    spl = [s[i] + d[i] for i in range(n)]
    smi = [s[i] - d[i] for i in range(n)]
    O = lambda j: (e(spl, j) - e(smi, j)) / 2
    L6 = sum(s[i] * d[i] for i in range(n)) + O(4) - 5 * O(6)
    return Q, H, L6


def families():
    b, p, q = sp.symbols('beta p q')
    r = list(sp.symbols('tau3:7'))
    Q, H, L6 = forms([b, -b, 0, 0, 0, 0], [p, q, *r], esym)
    assert sp.expand(Q - H - b ** 2 * (esym(r, 3) - esym(r, 1))) == 0
    assert sp.expand(L6 - b * (p - q) * (1 - esym(r, 2) + 5 * esym(r, 4))) == 0
    lam, mu = sp.symbols('lam mu', positive=True)
    c = 1 - lam / 5
    Q, H, L6 = forms([lam] + [-lam / 5] * 5, [mu] + [c] * 5, esym)
    A = Q - H
    assert sp.expand(A - lam ** 2 / 625 * (500 * (1 - mu) - lam * (1000 - 300 * mu) + lam ** 2 * (365 - 35 * mu) - 37 * lam ** 3)) == 0
    assert sp.expand(L6 - 8 * lam ** 2 * (5 - lam) ** 2 / 625 * (5 * (1 - mu) - lam * (5 - 2 * mu) + 2 * lam ** 2)) == 0
    assert sp.simplify(sp.limit(L6 / A, lam, 0) - 2) == 0


def example():
    s = [F(-3, 4), F(1, 2), F(2, 5), F(2, 5), F(-19, 20), F(2, 5)]
    d = [F(0)] * 4 + [F(1, 20), F(0)]
    assert sum(s) == 0 and all(0 <= d[i] <= 1 - abs(s[i]) for i in range(6))
    S2 = sum(x * x for x in s)
    assert max(x * x for x in s) / S2 == F(361, 878)
    Q, H, L6 = forms(s, d, esym_f)
    assert Q - H == F(1, 20000) and abs(L6) == F(167, 2500)


def flat_with_positive_vertex(rng):
    while True:
        s = [rng.uniform(-1, 1) for _ in range(5)]
        s.append(-sum(s))
        if max(abs(x) for x in s) >= 1:
            continue
        S2 = sum(x * x for x in s)
        if max(x * x for x in s) > S2 / 2:
            continue
        caps = [1 - abs(x) for x in s]
        pos = []
        for bits in product([0, 1], repeat=6):
            Q, H, _ = forms(s, [c if bt else 0.0 for c, bt in zip(caps, bits)], esym_f)
            if Q - H > 1e-12:
                pos.append(bits)
        if pos:
            return s, caps, pos


def search(n=25000, seed=7):
    rng = random.Random(seed)
    worst, positive, ratio = float("inf"), 0, float("inf")
    flat, positive_flat = 0, 0
    cur = None
    for t in range(n):
        mode = t % 5
        if mode == 4:
            if t % 250 == 4:
                cur = flat_with_positive_vertex(rng)
            s, caps, pos = cur
            if rng.random() < 0.5:
                bits = rng.choice(pos)
                d = [c * min(1.0, max(0.0, bt + rng.gauss(0, 0.25))) for c, bt in zip(caps, bits)]
            else:
                d = [c * rng.random() for c in caps]
        elif mode == 3:
            x = rng.uniform(0.01, 0.99)
            w = [rng.uniform(0.5, 1.5) for _ in range(5)]
            s = [x] + [-x * wi / sum(w) for wi in w]
            d = [(1 - x) * rng.uniform(0.5, 1.0)] + [(1 - abs(y)) * rng.uniform(0, 0.2) for y in s[1:]]
        else:
            while True:
                s = [rng.uniform(-0.9, 0.9) for _ in range(5)]
                s.append(-sum(s))
                if abs(s[-1]) < 1:
                    break
            d = []
            for x in s:
                u = rng.random()
                if mode == 1:
                    u = rng.choice([0.0, 1.0, u])
                d.append((1 - abs(x)) * u)
            if mode == 2:
                f = rng.uniform(0, 0.1)
                s = [x * f for x in s]
        assert abs(sum(s)) < 1e-9, "a drawn point is not centered"
        Q, H, L6 = forms(s, d, esym_f)
        A = Q - H
        scale = max(abs(A), abs(L6), 1e-300)
        worst = min(worst, (abs(L6) - A) / scale)
        S2 = sum(x * x for x in s)
        is_flat = max(x * x for x in s) <= S2 / 2
        flat += is_flat
        if A > 0:
            positive += 1
            positive_flat += is_flat
            ratio = min(ratio, abs(L6) / A)
    return worst, positive, ratio, flat, positive_flat


def main():
    families()
    example()
    worst, positive, ratio, flat, positive_flat = search()
    print(f"family formulas and the rational example hold; {flat} of 25000 points have no s_i^2 above S_2/2; "
          f"Q - H > 0 at {positive} points, {positive_flat} of them with no s_i^2 above S_2/2; "
          f"smallest (|L6| - (Q - H)) / max(|Q - H|, |L6|) = {worst:.3g}; smallest |L6| / (Q - H) where Q - H > 0 = {ratio:.4f}")
    assert worst > -1e-9


if __name__ == "__main__":
    main()
