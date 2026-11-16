#!/usr/bin/env python3
"""Independent sparse-polynomial verification of the centered four-bit proof."""
from fractions import Fraction as F
from itertools import combinations
import json,time
N=8;Z=(0,)*N

def const(v):return {Z:F(v)} if v else {}
def var(i):
    t=list(Z);t[i]=1;return {tuple(t):F(1)}
def add(*xs):
    out={}
    for x in xs:
        for k,v in x.items():out[k]=out.get(k,F(0))+v
    return {k:v for k,v in out.items() if v}
def scale(x,v):return {k:a*v for k,a in x.items() if a*v}
def sub(x,y):return add(x,scale(y,-1))
def mul(x,y):
    out={}
    for a,av in x.items():
        for b,bv in y.items():
            k=tuple(i+j for i,j in zip(a,b));out[k]=out.get(k,F(0))+av*bv
    return {k:v for k,v in out.items() if v}
def prod(xs):
    out=const(1)
    for x in xs:out=mul(out,x)
    return out
def power(x,n):return prod([x]*n)
def es(xs,n):return add(*(prod(ts) for ts in combinations(xs,n)))
def pmf(ps):
    out=[const(1)]
    for p in ps:
        q=sub(const(1),p);new=[{} for _ in range(len(out)+1)]
        for k,v in enumerate(out):new[k]=add(new[k],mul(v,q));new[k+1]=add(new[k+1],mul(v,p))
        out=new
    return out

def QL(s,d):
    delta=add(*d)
    Q=add(*(mul(sub(di,scale(delta,F(1,2))),power(si,2)) for si,di in zip(s,d)))
    L=add(*(mul(si,di) for si,di in zip(s,d)),scale(add(*(add(mul(d[i],prod(s[:i]+s[i+1:])),mul(s[i],prod(d[:i]+d[i+1:]))) for i in range(4))),3))
    return Q,L

def main():
    start=time.monotonic();checks=[]
    s=[var(i) for i in range(3)];s.append(scale(add(*s),-1));d=[var(i) for i in range(3,7)]
    p=[scale(add(const(1),si,di),F(1,2)) for si,di in zip(s,d)]
    q=[scale(add(const(1),si,scale(di,-1)),F(1,2)) for si,di in zip(s,d)]
    mp,mq=pmf(p),pmf(q);D=[sub(add(*mp[k:]),add(*mq[k:])) for k in range(1,5)]
    Q,L=QL(s,d);delta=add(*d)
    assert sub(scale(add(D[1],D[2]),F(1,2)),scale(add(scale(delta,3),scale(es(d,3),-1),scale(Q,-1)),F(1,8)))=={}
    assert sub(sub(D[1],D[2]),scale(L,F(1,4)))=={}
    checks.append('central average and signed contrast: zero polynomial residuals')
    gap=scale(add(*(mul(power(sub(d[i],d[j]),2),sub(scale(delta,5),scale(add(d[i],d[j]),4))) for i,j in combinations(range(4),2))),F(1,48))
    assert sub(sub(scale(power(delta,3),F(1,16)),es(d,3)),gap)=={}
    checks.append('gap Maclaurin sum-of-squares: zero polynomial residual')
    # Simplex lemma at a determining vertex of its r-polytope, homogenized.
    a,b,c=var(0),var(1),var(2);one=add(a,b,c)
    A=scale(add(scale(a,4),scale(b,2)),F(1,3))
    Fhom=sub(scale(add(scale(mul(b,c),-1),mul(a,c),scale(mul(a,b),3)),F(1,3)),{})
    target=scale(add(mul(A,power(one,2)),scale(mul(Fhom,one),-3),scale(prod([a,b,c]),3)),3)
    certificate=add(scale(power(a,3),4),mul(b,power(sub(a,scale(b,F(1,2))),2)),scale(power(b,3),F(7,4)),scale(prod([a,a,c]),5),scale(prod([b,b,c]),7),prod([a,c,c]),scale(prod([b,c,c]),5),scale(prod([a,b,c]),12))
    assert sub(target,certificate)=={}
    checks.append('six-vertex simplex lemma: exact nonnegative cubic certificate')
    # General transformed Q and L; r3=1-r1-r2.
    x,r1,r2,t,e1,e2,e3=[var(i) for i in range(7)];r=[r1,r2,sub(const(1),add(r1,r2))];e=[e1,e2,e3];h=add(*e)
    s=[x]+[scale(mul(x,ri),-1) for ri in r];d=[t]+e
    Q,L=QL(s,d);norm=add(*(power(ri,2) for ri in r));cc=scale(sub(const(1),norm),F(1,2))
    bb=[sub(scale(add(const(1),norm),F(1,2)),power(ri,2)) for ri in r]
    assert sub(Q,mul(power(x,2),sub(mul(cc,t),add(*(mul(bi,ei) for bi,ei in zip(bb,e))))))=={}
    slope=add(const(1),scale(mul(power(x,2),prod(r)),-3),scale(add(*(mul(r[i],prod(e[:i]+e[i+1:])) for i in range(3))),-3))
    free=add(*(mul(e[i],add(scale(r[i],-1),scale(mul(power(x,2),prod(r[:i]+r[i+1:])),3))) for i in range(3)),scale(prod(e),3))
    assert sub(L,mul(x,add(mul(t,slope),free)))=={}
    A=add(*(mul(ei,sub(const(1),ri)) for ei,ri in zip(e,r)))
    B=add(*(prod([e[i],sub(const(1),r[i]),*r[:i],*r[i+1:]]) for i in range(3)))
    F3=sub(mul(h,add(*(mul(r[i],prod(e[:i]+e[i+1:])) for i in range(3)))),prod(e))
    assert sub(add(mul(h,slope),free),add(A,scale(mul(power(x,2),B),3),scale(F3,-3)))=={}
    checks.append('dominant midpoint/gap substitutions and affine decomposition: zero polynomial residuals')
    out={'checks':checks,'nonzero_residual_monomials':0,'seconds':time.monotonic()-start}
    print(json.dumps(out,indent=2));open('certificate_checks.json','w').write(json.dumps(out,indent=2)+'\n')
if __name__=='__main__':main()
