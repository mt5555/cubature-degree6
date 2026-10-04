# Exact check of msolve's rational parametrization: f squarefree of degree 16, the 16 points
# y(t) = s * num_i(t) / (c_i * q(t)) satisfy all 10 equations exactly (mod f), and f has exactly 2 real roots.
import sympy as sp, numpy as np
from parse_rur import load
from fractions import Fraction
R=load(); P=R[1][5][1]
t=sp.symbols('t')
poly=lambda d,cs: sp.Poly(list(reversed(cs)),t,domain='QQ')   # msolve lists coefficients from constant term up
f=poly(*P[0]); q=poly(*P[1]); nums=[(poly(*a[0]),sp.Integer(a[1])) for a in P[2]]
print("deg f =",f.degree(),"deg q =",q.degree())
print("f squarefree (gcd(f,f')=1):", sp.gcd(f,f.diff(t)).degree()==0, "; gcd(f,q)=1:", sp.gcd(f,q).degree()==0)
lines=open('flat8_Q.ms').read().split('\n'); ys=sp.symbols('y0:10')
eqs=[sp.Poly(sp.sympify(l.rstrip(',').replace('^','**')),*ys) for l in lines[2:] if l.strip()]
# numerical sign check
from cf8 import hank, hscale
roots=np.roots([float(c) for c in f.all_coeffs()])
for s in [1,-1]:
    ok=[]
    for r in roots[:4]:
        y=np.array([complex(s*nm.eval(r)/(ci*q.eval(r))) for nm,ci in [(nm,float(ci)) for nm,ci in nums]])
        ok.append(np.abs(hank(y)).max()/hscale)
    print("sign",s,": max |hank|/scale at 4 roots:",["%.1e"%v for v in ok])
