import sympy as sp, mpmath as mp
from parse_rur import load
mp.mp.dps=200
R=load(); P=R[1][5][1]
fc=[mp.mpf(c) for c in reversed(P[0][1])]; qc=list(reversed(P[1][1]))
nums=[(list(reversed(a[0][1])),a[1]) for a in P[2]]
roots=mp.polyroots(fc,maxsteps=500,extraprec=2000)
from cf8 import hank_sym, ys
import numpy as np
H=[sp.lambdify(ys,h,'mpmath') for h in hank_sym(list(ys))]
pv=lambda cs,x: mp.polyval([mp.mpf(c) for c in cs],x)
for s in [1,-1]:
    res=[]
    for r in roots[:3]:
        y=[s*pv(n,r)/(c*pv(qc,r)) for n,c in nums]
        res.append(max(abs(h(*y)) for h in H))
    print("sign",s,[mp.nstr(v,3) for v in res])
print("real roots:",[mp.nstr(r,12) for r in roots if abs(mp.im(r))<mp.mpf(10)**-100])
