"""Section 4 (six summands): exact checks behind the gap inequality
    5 Delta - e_3(d) + e_5(d) >= 30 g - 20 g^3 + 6 g^5,  g = Delta / 6,  d in [0,1]^6.
(1) With d = (d_i, d_j, eta), the identity
    e_3(d) - e_5(d) = d_i d_j (e_1(eta) - e_3(eta)) + (d_i + d_j)(e_2(eta) - e_4(eta)) + e_3(eta)
holds as a polynomial identity (SymPy).
(2) At d = (g, ..., g) the two sides agree: e_3 = 20 g^3 and e_5 = 6 g^5.
(3) e_1(eta) - e_3(eta) >= 0 is checked exactly at random rational points of [0,1]^4 (it is proved in the paper).
Needs SymPy; stops with an error if a check fails."""
from fractions import Fraction as F
from itertools import combinations
import random
import sympy as sp


def esym(v, k):
    return sp.Add(*[sp.Mul(*[v[i] for i in I]) for I in combinations(range(len(v)), k)])


def main():
    di, dj, g = sp.symbols('d_i d_j g')
    eta = list(sp.symbols('r1:5'))
    d = [di, dj, *eta]
    lhs = esym(d, 3) - esym(d, 5)
    rhs = di * dj * (esym(eta, 1) - esym(eta, 3)) + (di + dj) * (esym(eta, 2) - esym(eta, 4)) + esym(eta, 3)
    assert sp.expand(lhs - rhs) == 0
    eq = [g] * 6
    assert sp.expand(esym(eq, 3) - 20 * g ** 3) == 0 and sp.expand(esym(eq, 5) - 6 * g ** 5) == 0
    rng = random.Random(5)
    for _ in range(2000):
        r = [F(rng.randint(0, 50), 50) for _ in range(4)]
        e1 = sum(r)
        e3 = sum(r[a] * r[b] * r[c] for a, b, c in combinations(range(4), 3))
        assert e1 >= e3
    print("gap identity, equal-gap values and e_3(eta) <= e_1(eta) at 2000 rational points: all hold")


if __name__ == "__main__":
    main()
