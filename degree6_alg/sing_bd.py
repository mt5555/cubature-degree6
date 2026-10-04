# singular points of the boundary curve V ∩ {det G_g = 0}:  h = 0, det_g = 0, grad det_g = sum c_k grad h_k
import sympy as sp, pickle, sys
d=pickle.load(open('data6.pkl','rb')); g=sys.argv[1]
ys=sp.symbols('y0:8'); c=sp.symbols('C0:6')
H=[sp.sympify(e) for e in d['H']]; Dg=sp.sympify(d['DET'][g])
eqs=[sp.diff(Dg,ys[i])-sum(c[k]*sp.diff(H[k],ys[i]) for k in range(6)) for i in range(8)]+H+[Dg]
lines=[]
for e in eqs:
    P=sp.Poly(sp.expand(e),*ys,*c); cc=sp.gcd_list([sp.Integer(x) for x in P.coeffs()])
    lines.append(str(P.as_expr()/cc).replace('**','^'))
open(f'sbd_{g}.ms','w').write(','.join(map(str,list(ys)+list(c)))+'\n0\n'+',\n'.join(lines)+'\n')
