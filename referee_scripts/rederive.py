# Independent exact re-derivation: 10x10 localizing matrices + flatness, compare to sosdata.pkl
import sympy as sp, pickle, itertools
from fractions import Fraction as Fr
from math import factorial
Q=sp.Rational
def mom(a,b): return Q(factorial(a)*factorial(b), factorial(a+b+2))
y=sp.symbols('y0:8'); t=sp.symbols('t0:8')
def L(a,b):
    assert a>=0 and b>=0 and a+b<=7
    return y[b] if a+b==7 else mom(a,b)
P3=[(a,d-a) for d in range(4) for a in range(d,-1,-1)]        # 10 monomials, P2 first (6) then degree 3 (4)
D4=[(4-j,j) for j in range(5)]
gdef={'x':{(1,0):1},'y':{(0,1):1},'h':{(0,0):1,(1,0):-1,(0,1):-1}}
def Lg(g,a,b): return sum(c*L(a+da,b+db) for (da,db),c in gdef[g].items())
M3g={g:sp.Matrix(10,10,lambda i,j: Lg(g,P3[i][0]+P3[j][0],P3[i][1]+P3[j][1])) for g in gdef}
# Schur complement on the P2 block
Gmine={}
for g in gdef:
    M=M3g[g]; A=M[:6,:6]; F=M[:6,6:]; D=M[6:,6:]
    assert A.free_symbols==set() and F.free_symbols==set()
    # A must be PD (exact): check leading principal minors
    assert all(A[:k,:k].det()>0 for k in range(1,7)), g
    Gmine[g]=(D-F.T*A.inv()*F).applyfunc(sp.expand)
# flatness: M3 (Lebesgue, fixed), B=M[P3,D4], C=B^T M3^{-1} B, Hankel differences
M3=sp.Matrix(10,10,lambda i,j: mom(P3[i][0]+P3[j][0],P3[i][1]+P3[j][1]))
assert all(M3[:k,:k].det()>0 for k in range(1,11))
B=sp.Matrix(10,5,lambda i,j: L(P3[i][0]+D4[j][0],P3[i][1]+D4[j][1]))
C=(B.T*M3.inv()*B).applyfunc(sp.expand)
# the (i,j) entry of C is the moment of x^{8-(i+j)} y^{i+j}: group by s=i+j
groups={}
for i in range(5):
    for j in range(i,5): groups.setdefault(i+j,[]).append((i,j))
hmine=[]
for s,lst in sorted(groups.items()):
    for q in lst[1:]: hmine.append(sp.expand(C[lst[0]]-C[q]))
assert len(hmine)==6 and len(groups)==9
# compare with cf2's symbolic objects (the write-up's G_g = D_g - K_g and hank_sym)
import cf2
Gcf=cf2.G_sym(list(y)); hcf=cf2.hank_sym(list(y))
print("(a) Schur 4x4 == cf2 G_sym:",all(sp.simplify(Gmine[g]-Gcf[g])==sp.zeros(4,4) for g in gdef))
print("(b) hankel eqs == cf2 hank_sym (as polynomials, same order):",[sp.expand(a-b)==0 for a,b in zip(hmine,hcf)])
# now compare with sosdata.pkl: GT (after y=c+hw t and congruence by W) and EQ (up to scalar)
dat=pickle.load(open('sosdata.pkl','rb')); GT=dat['GT']; EQ=dat['EQ']; c=dat['c']; hw=dat['hw']
# W is not stored in sosdata.pkl: recompute as the exact dyadic Fr(float(Wm)) from cf2 the way sosdata does
W={g:sp.Matrix(4,4,lambda i,j: Q(Fr(float(cf2.Wm[g][i,j])))) for g in cf2.Wm}
sub={y[k]:Q(c[k])+Q(hw[k])*t[k] for k in range(8)}
okGT=True
for g in gdef:
    Gt=(W[g]*Gmine[g].subs(sub)*W[g].T).applyfunc(sp.expand)
    B0=sp.Matrix(4,4,lambda i,j: Q(GT[g][0][i][j])); Bk=[sp.Matrix(4,4,lambda i,j: Q(GT[g][1][k][i][j])) for k in range(8)]
    rec=B0+sum((t[k]*Bk[k] for k in range(8)),sp.zeros(4,4))
    okGT&= (Gt-rec).applyfunc(sp.expand)==sp.zeros(4,4)
    print(f"   W_{g} det = {W[g].det()} (nonsingular: {W[g].det()!=0}), lower-triangular: {all(W[g][i,j]==0 for i in range(4) for j in range(i+1,4))}")
print("(c) GT in sosdata.pkl == W * Schur(y=c+hw t) * W^T exactly:",okGT)
def todict(expr):
    P=sp.Poly(sp.expand(expr),*t,domain='QQ'); d={}
    for mon,co in P.terms():
        d[tuple(i for i in range(8) for _ in range(mon[i]))]=Fr(int(co.p),int(co.q))
    return d
mine=[todict(h.subs(sub)) for h in hmine]
matched=[]
for e in EQ:
    found=None
    for idx,md in enumerate(mine):
        if set(md)!=set(e): continue
        m0=next(iter(e)); r=e[m0]/md[m0]
        if all(e[m]==r*md[m] for m in e): found=(idx,r); break
    matched.append(found)
print("(d) each EQ_k equals (nonzero scalar) * own hankel eq:",matched)
print("    all matched, distinct, scalars nonzero:",all(m is not None and m[1]!=0 for m in matched) and len({m[0] for m in matched if m})==6)
print("    EQ coefficient types all Fraction:",all(type(v) is Fr for e in EQ for v in e.values()), " c,hw Fraction:",all(type(v) is Fr for v in c+hw), " hw nonzero:",all(v!=0 for v in hw))
