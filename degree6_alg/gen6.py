# exact degree-6 data in y-coordinates (from ../cf2.py): 6 flatness quadrics h, 3 localizing determinants det G_g
import sys, sympy as sp
sys.path.insert(0,'..')
from cf2 import hank_sym, G_sym, ys
def intpoly(e):
    num,den=sp.fraction(sp.together(sp.expand(e))); P=sp.Poly(sp.expand(num),*ys)
    g=sp.gcd_list([sp.Integer(c) for c in P.coeffs()]); return sp.Poly(P.as_expr()/g,*ys)
H=[intpoly(h) for h in hank_sym(list(ys))]
Gs=G_sym(list(ys)); DET={g:intpoly(Gs[g].det()) for g in Gs}
def s(P): return str(P.as_expr()).replace('**','^')
if __name__=="__main__":
    import pickle
    pickle.dump({'H':[str(P.as_expr()) for P in H],'DET':{g:str(P.as_expr()) for g,P in DET.items()}},open('data6.pkl','wb'))
    print("h degrees",[P.total_degree() for P in H]," det degrees",{g:P.total_degree() for g,P in DET.items()})
    ring="ring r = 0,(y0,y1,y2,y3,y4,y5,y6,y7),dp;\n"
    with open('V.sing','w') as f:
        f.write('LIB "modstd.lib";\n'+ring+"ideal I =\n"+",\n".join(s(P) for P in H)+";\n")
        f.write('ideal G = modStd(I,1);\n"dim V (Krull):"; dim(G);\n"degree of V:"; mult(G);\n"Hilbert polynomial check (hilb):"; hilb(G,2);\nquit;\n')
