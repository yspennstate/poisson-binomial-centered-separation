"""Check the identities of the paper in its own notation (SymPy), and the equal-gap formulas in exact rationals."""
import sys, itertools
from fractions import Fraction as Fr
from math import comb
import sympy as sp
sys.stdout.reconfigure(encoding="utf-8")
ok = True
def chk(name, expr):
    global ok
    e = sp.expand(expr)
    print(name, "OK" if e == 0 else f"FAIL {e}")
    ok &= (e == 0)

# majority M = (3e1 - e3)/8 on {-1,1}^4
for z in itertools.product([-1, 1], repeat=4):
    pos = sum(1 for v in z if v > 0)
    M = 1 if pos > 2 else (0 if pos == 2 else -1)
    e1 = sum(z); e3 = sum(z[i]*z[j]*z[k] for i, j, k in itertools.combinations(range(4), 3))
    assert Fr(3*e1 - e3, 8) == M
print("majority OK")

s = sp.symbols("s1:5"); d = sp.symbols("d1:5")
sub = {s[3]: -(s[0]+s[1]+s[2])}
p = [(1 + s[i] + d[i])/2 for i in range(4)]; q = [(1 + s[i] - d[i])/2 for i in range(4)]
t_ = sp.Symbol("t")
def pgf(pr): return sp.expand(sp.prod([(1-x) + x*t_ for x in pr]))
GX, GY = pgf(p), pgf(q)
PX = [GX.coeff(t_, k) for k in range(5)]; PY = [GY.coeff(t_, k) for k in range(5)]
D2 = (PX[2]+PX[3]+PX[4]) - (PY[2]+PY[3]+PY[4]); D3 = (PX[3]+PX[4]) - (PY[3]+PY[4])
Delta = sum(d); S2 = sum(x**2 for x in s)
e3d = sum(d[i]*d[j]*d[k] for i, j, k in itertools.combinations(range(4), 3))
Q = sum((d[i] - Delta/2)*s[i]**2 for i in range(4))
L = sum(s[i]*d[i] for i in range(4)) + 3*sum(d[i]*sp.prod([s[j] for j in range(4) if j != i]) + s[i]*sp.prod([d[j] for j in range(4) if j != i]) for i in range(4))
chk("avg4", ((D2+D3)/2 - (3*Delta - e3d - Q)/8).subs(sub))
chk("diff4", ((D2-D3) - L/4).subs(sub))
chk("gap", Delta**3/16 - e3d - sp.Rational(1, 48)*sum((d[i]-d[j])**2*(5*Delta - 4*(d[i]+d[j])) for i, j in itertools.combinations(range(4), 2)))
# b4
g = sp.Symbol("g")
b4 = sum(comb(4, k)*(((1+g)/2)**k*((1-g)/2)**(4-k) - ((1-g)/2)**k*((1+g)/2)**(4-k)) for k in range(2, 5))
chk("b4", b4 - (sp.Rational(3, 2)*g - g**3/2))
chk("b4Delta", b4.subs(g, sp.Symbol("D")/4) - (sp.Rational(3, 8)*sp.Symbol("D") - sp.Symbol("D")**3/128))

# coordinates of section 3: s=(x,-x y1,-x y2,-x y3), d=(t,u1,u2,u3), y1+y2+y3=1
x, t = sp.symbols("x t"); y = sp.symbols("y1:4"); u = sp.symbols("u1:4")
ysub = {y[2]: 1 - y[0] - y[1]}
sc = [x, -x*y[0], -x*y[1], -x*y[2]]; dc = [t, u[0], u[1], u[2]]
rep = dict(zip(s, sc)); rep.update(dict(zip(d, dc)))
Qc = sum((dc[i] - sum(dc)/2)*sc[i]**2 for i in range(4))
rho = sum(v**2 for v in y); kap = (1 - rho)/2; th = [(1 + rho)/2 - v**2 for v in y]
chk("Q coords", (Qc - x**2*(kap*t - sum(th[i]*u[i] for i in range(3)))).subs(ysub))
Lc = L.subs(rep)
C = 1 - 3*x**2*y[0]*y[1]*y[2] - 3*sum(y[i]*sp.prod([u[j] for j in range(3) if j != i]) for i in range(3))
G = sum(u[i]*(3*x**2*sp.prod([y[j] for j in range(3) if j != i]) - y[i]) for i in range(3)) + 3*u[0]*u[1]*u[2]
chk("affine", (Lc/x - (t*C + G)).subs(ysub))
sig = sum(u)
v = [u[i]/sig for i in range(3)]
U = sum(v[i]*(1 - y[i]) for i in range(3))
F = sum(y[i]*sp.prod([v[j] for j in range(3) if j != i]) for i in range(3)) - v[0]*v[1]*v[2]
V = sum(v[i]*(1 - y[i])*sp.prod([y[j] for j in range(3) if j != i]) for i in range(3))
mid = sum(u[i]*(1 - y[i]) for i in range(3)) + 3*x**2*sum(u[i]*(1 - y[i])*sp.prod([y[j] for j in range(3) if j != i]) for i in range(3)) - 3*(sig*sum(y[i]*sp.prod([u[j] for j in range(3) if j != i]) for i in range(3)) - u[0]*u[1]*u[2])
chk("regroup1", sp.simplify((sig*C + G - mid).subs(ysub)))
chk("regroup2", sp.simplify((mid - sig*(U + 3*x**2*V - 3*sig**2*F)).subs(ysub)))
# vertex expansion at y=(-1/3,1/3,1), v homogenized
v1, v2, v3 = sp.symbols("v1:4")
yv = [sp.Rational(-1, 3), sp.Rational(1, 3), 1]; vv = [v1, v2, v3]; S = v1 + v2 + v3
Uv = sum(vv[i]*(1 - yv[i]) for i in range(3))*S**2
Fv = sum(yv[i]*sp.prod([vv[j] for j in range(3) if j != i]) for i in range(3))*S - v1*v2*v3
rhs = 4*v1**3 + v2*(v1 - v2/2)**2 + sp.Rational(7, 4)*v2**3 + 5*v1**2*v3 + 7*v2**2*v3 + v1*v3**2 + 5*v2*v3**2 + 12*v1*v2*v3
chk("vertex", 3*(Uv - 3*Fv) - rhs)
# sharpness example
tt = sp.Symbol("tt"); xx = 1 - tt
ex = {s[0]: xx, s[1]: -xx/3, s[2]: -xx/3, s[3]: -xx/3, d[0]: tt, d[1]: 0, d[2]: 0, d[3]: 0}
chk("sharp Q", Q.subs(ex) - tt*xx**2/3)
chk("sharp L", L.subs(ex) - (xx*tt - tt*xx**3/9))

# equal-gap theorem: jumps formula, exact remainder, bound, exact rationals
def binpmf(n, pp, k): return comb(n, k)*pp**k*(1-pp)**(n-k) if 0 <= k <= n else Fr(0)
def hfun(k, m): return -1 if k < m else (0 if k == m else 1)
for m in range(2, 7):
    for gg in (Fr(1, 7), Fr(1, 2), Fr(5, 6)):
        a = 1 - gg
        phi = []
        for w in range(2*m + 1):
            e1 = sum(binpmf(2*m - w, gg, k)*hfun(w + k, m) for k in range(2*m - w + 1))
            e2 = sum(binpmf(w, a, k)*hfun(k, m) for k in range(w + 1))
            phi.append((e1 - e2)/2)
        assert phi[m] == 1 - a**m
        A = [phi[m - j] - phi[m - j - 1] for j in range(m)]
        assert A[0] == Fr(m - 1, 2)*gg*a**m
        for j in range(1, m):
            assert A[j] == Fr(1, 2)*comb(m + j, j)*gg**j*a**m*(a + Fr(m, j + 1)*gg), (m, gg, j)
        for w in range(2*m + 1):
            assert phi[w] == phi[2*m - w]
            assert phi[w] == phi[0] + sum(A[j] for j in range(m) if abs(w - m) <= j)
        bound = Fr(m + 1, 6)*gg**(m - 1)*a**m
        assert all(A[j]/(2*j + 1) >= bound for j in range(m))
print("equal-gap formulas OK")
print("ALL OK" if ok else "SOME FAIL")
