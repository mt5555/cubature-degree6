# Lagrange system for critical points of l on the smooth part of V: grad l = sum lam_k grad h_k, h = 0
import sympy as sp, pickle
d=pickle.load(open('data6.pkl','rb'))
ys=sp.symbols('y0:8'); lam=sp.symbols('L0:6')
H=[sp.sympify(e) for e in d['H']]
ell=[3,-5,7,2,-4,6,-1,8]
eqs=[ell[i]-sum(lam[k]*sp.diff(H[k],ys[i]) for k in range(6)) for i in range(8)]+H
lines=[]
for e in eqs:
    P=sp.Poly(sp.expand(e),*ys,*lam); g=sp.gcd_list([sp.Integer(c) for c in P.coeffs()])
    lines.append(str((P.as_expr()/g)).replace('**','^'))
open('lagrange_Q.ms','w').write(','.join(map(str,ys+lam))+'\n0\n'+',\n'.join(lines)+'\n')
print("wrote lagrange_Q.ms")
