# Exact verification of the rational parametrization: E_k(y(t)) * q(t)^2 == 0 mod f(t) for all 10 equations,
# with y_i = -n_i(t)/(c_i q(t)); f squarefree of degree 16 with exactly 2 real roots (Sturm).
import sympy as sp, time
from parse_rur import load
t0=time.time(); t=sp.symbols('t')
R=load(); P=R[1][5][1]
poly=lambda d,cs: sp.Poly(list(reversed(cs)),t,domain='ZZ')
f=poly(*P[0]); q=poly(*P[1]); nums=[(poly(*a[0]),int(a[1])) for a in P[2]]
lines=open('flat8_Q.ms').read().split('\n'); ys=sp.symbols('y0:10')
eqs=[sp.Poly(sp.sympify(l.rstrip(',').replace('^','**')),*ys) for l in lines[2:] if l.strip()]
# homogenize by q: y_i = -n_i/(c_i q)  =>  y_i * q = -n_i / c_i  =: Y_i (polynomial in t over Q)
Y=[sp.Poly(-nm.as_expr()/c, t, domain='QQ') for nm,c in nums]
qQ=sp.Poly(q.as_expr(),t,domain='QQ'); fQ=sp.Poly(f.as_expr(),t,domain='QQ')
allzero=True
for k,E in enumerate(eqs):
    acc=sp.Poly(0,t,domain='QQ')
    for mon,co in E.terms():
        deg=sum(mon); term=sp.Poly(co,t,domain='QQ')*qQ**(2-deg)
        for i,e in enumerate(mon):
            if e: term=term*Y[i]**e
        acc=acc+term
    r=acc.rem(fQ); allzero&= r.is_zero
    print(f"  equation {k}: remainder mod f is zero: {r.is_zero}",flush=True)
print("all 10 equations vanish on the 16 parametrized points:",allzero)
print("f squarefree:",sp.gcd(f,f.diff(t)).degree()==0,"; gcd(f,q)=1:",sp.gcd(f,q).degree()==0)
print("number of real roots of f (exact, Sturm):",sp.Poly(f.as_expr(),t).count_roots(), f"  [{time.time()-t0:.0f}s]")
