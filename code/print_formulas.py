#!/usr/bin/env python3
"""Print every displayed low-dimensional coefficient as an exact rational."""
from fractions import Fraction as F
from math import comb
import json
out=[]
for n in [2,4,6,8,10,12]:
    m=n//2
    multiplier=F(n*comb(n-1,m-1),2**(n-1))
    coefficients={str(2*j+1):str(multiplier*F((-1)**j*comb(m-1,j),2*j+1)) for j in range(m)}
    # Jump ratios after extracting (1-g)^m/2; values encode g^j * (a + m*g/(j+1)).
    jump_ratios=[{'g_power':1,'coefficient':str(m-1),'factor':'1'}]+[{'g_power':j,'coefficient':str(F(comb(m+j,j),2*j+1)),'factor':f'(1-g)+{F(m,j+1)}*g'} for j in range(1,m)]
    out.append({'n':n,'b_polynomial_power_to_coefficient':coefficients,'c_simple_prefactor':str(F(m+1,24)) if m>=2 else '0','c_simple_g_power':m-1 if m>=2 else None,'c_simple_1_minus_g_power':m-2 if m>=2 else None,'jump_ratio_terms':jump_ratios})
print(json.dumps(out,indent=2))
