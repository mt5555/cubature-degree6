# critical points of l on the boundary curve V ∩ {det G_g = 0}:  grad l = sum lam_k grad h_k + mu grad det_g
import sympy as sp, pickle, sys
d=pickle.load(open('data6.pkl','rb')); g=sys.argv[1]
ys=sp.symbols('y0:8'); lam=sp.symbols('L0:6'); mu=sp.symbols('M')
H=[sp.sympify(e) for e in d['H']]; Dg=sp.sympify(d['DET'][g])
ell=[3,-5,7,2,-4,6,-1,8]
eqs=[ell[i]-sum(lam[k]*sp.diff(H[k],ys[i]) for k in range(6))-mu*sp.diff(Dg,ys[i]) for i in range(8)]+H+[Dg]
lines=[]
for e in eqs:
    P=sp.Poly(sp.expand(e),*ys,*lam,mu); c=sp.gcd_list([sp.Integer(x) for x in P.coeffs()])
    lines.append(str(P.as_expr()/c).replace('**','^'))
open(f'bd_{g}.ms','w').write(','.join(map(str,list(ys)+list(lam)+[mu]))+'\n0\n'+',\n'.join(lines)+'\n')
