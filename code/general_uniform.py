#!/usr/bin/env python3
"""Exact forcing-kernel identities and stability checks; standard library only."""
from fractions import Fraction as F
from math import comb
from itertools import combinations_with_replacement
from random import Random
import json, time

def pmf(ps):
    v=[F(1)]
    for p in ps:
        v=[(v[k] if k<len(v) else 0)*(1-p)+(v[k-1] if k else 0)*p for k in range(len(v)+1)]
    return v

def binpmf(n,p):
    return [F(comb(n,k))*p**k*(1-p)**(n-k) for k in range(n+1)]

def tails(ps,qs):
    p,q=pmf(ps),pmf(qs)
    return [sum(p[k:])-sum(q[k:]) for k in range(1,len(p))]

def benchmark(n,g):
    p,q=binpmf(n,(1+g)/2),binpmf(n,(1-g)/2)
    return sum(p[n//2:])-sum(q[n//2:])

def phi(n,g,w):
    m=n//2
    h=lambda k: (k>m)-(k<m)
    return (sum(h(w+k)*v for k,v in enumerate(binpmf(n-w,g)))-sum(h(k)*v for k,v in enumerate(binpmf(w,1-g))))/2

def jumps(n,g):
    m=n//2;a=1-g
    return [F(m-1,2)*g*a**m]+[F(comb(m+j,j),2)*g**j*a**m*(a+F(m,j+1)*g) for j in range(1,m)]

def stability(n,g):
    if n==2:return F(0)
    m=n//2
    return F(m+1,24)*g**(m-1)*(1-g)**(m-2)

def check(rs,g):
    n=len(rs);m=n//2;a=1-g
    assert sum(rs)==m and n%2==0
    ss=[a*(2*r-1) for r in rs]
    ds=tails([g+a*r for r in rs],[a*r for r in rs])
    avg=(ds[m-1]+ds[m])/2;base=benchmark(n,g)
    wp=pmf(rs);bp=binpmf(n,F(1,2))
    delta=[sum(wp[m-j:m+j+1])-sum(bp[m-j:m+j+1]) for j in range(m)]
    aa=jumps(n,g)
    assert all(v>=0 for v in delta)
    assert all(phi(n,g,m-j)-phi(n,g,m-j-1)==aa[j] for j in range(m))
    assert all(phi(n,g,w)==phi(n,g,n-w) for w in range(n+1))
    assert avg==sum(phi(n,g,w)*wp[w] for w in range(n+1))
    assert avg-base==sum(x*y for x,y in zip(aa,delta))
    assert sum((2*j+1)*v for j,v in enumerate(delta))==sum((r-F(1,2))**2 for r in rs)
    assert avg-base>=stability(n,g)*sum(s*s for s in ss)
    if n>=4 and 0<g<1:
        sharp=min(aa[j]/(2*j+1) for j in range(m))/(4*a*a)
        assert sharp>=stability(n,g)
        assert (avg==base)==all(r==F(1,2) for r in rs)
    return {'n':n,'g':str(g),'r':[str(v) for v in rs],'average':str(avg),'benchmark':str(base),'remainder':str(avg-base),'S2':str(sum(s*s for s in ss)),'c_simple':str(stability(n,g)), 'c_min':str(min(aa[j]/(2*j+1) for j in range(m))/(4*a*a)) if a else None}

def main():
    start=time.monotonic();results=[];rng=Random(4911116)
    for n in [2,4,6,8,10,12,16,20]:
        count=0
        for v in combinations_with_replacement(range(5),n) if n<=12 else []:
            if sum(v)!=2*n:continue
            for g in [F(0),F(1,10),F(1,2),F(9,10),F(1)]:
                check([F(t,4) for t in v],g);count+=1
        for trial in range(40):
            # Random balancing moves give centered, generally unpaired rational vectors.
            r=[F(1,2)]*n
            for step in range(2*n):
                i,j=rng.sample(range(n),2)
                lo=max(-r[i],r[j]-1);hi=min(1-r[i],r[j])
                stepval=lo+(hi-lo)*F(rng.randrange(9),8)
                r[i]+=stepval;r[j]-=stepval
            check(r,F(rng.randrange(1,16),16));count+=1
        results.append({'n':n,'exact_cases':count})
    # Small g and asymmetric baseline configurations; endpoint degeneracy is explicit.
    examples=[]
    for n in [4,6,8,10,12]:
        r=[F(0),F(1,3),F(2,3),F(1)]+[F(1,2)]*(n-4)
        for g in [F(1,1000),F(1,2),F(999,1000)]:examples.append(check(r,g))
    output={'seed':4911116,'checks':results,'examples':examples,'seconds':time.monotonic()-start}
    print(json.dumps(output,indent=2));open('uniform_checks.json','w').write(json.dumps(output,indent=2)+'\n')
if __name__=='__main__':main()
