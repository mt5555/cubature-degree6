# Independent exact re-derivation of the data in sosdata.pkl (adapted from the referee's rederive.py).
# Builds the 10x10 localizing matrices M3(g L) and the flatness condition symbolically in the degree-7
# moments y, then checks EXACTLY that sosdata.pkl's GT and EQ are these objects after y = c + hw*t,
# congruence by W (GT), and nonzero rational scaling (EQ).
import sympy as sp, pickle
from fractions import Fraction as Fr
from math import factorial
Q=sp.Rational
def mom(a,b): return Q(factorial(a)*factorial(b), factorial(a+b+2))   # Lebesgue moments on T
y=sp.symbols('y0:8'); t=sp.symbols('t0:8')
def L(a,b):
    assert a>=0 and b>=0 and a+b<=7
    return y[b] if a+b==7 else mom(a,b)
P3=[(a,d-a) for d in range(4) for a in range(d,-1,-1)]      # P2 first (6), then degree 3 (4)
D4=[(4-j,j) for j in range(5)]
gdef={'x':{(1,0):1},'y':{(0,1):1},'h':{(0,0):1,(1,0):-1,(0,1):-1}}
def Lg(g,a,b): return sum(c*L(a+da,b+db) for (da,db),c in gdef[g].items())
G={}
for g in gdef:
    M=sp.Matrix(10,10,lambda i,j: Lg(g,P3[i][0]+P3[j][0],P3[i][1]+P3[j][1]))
    A=M[:6,:6]; F=M[:6,6:]; D=M[6:,6:]
    assert A.free_symbols==set() and F.free_symbols==set()
    assert all(A[:k,:k].det()>0 for k in range(1,7)), g          # A_g PD, so the Schur step is valid
    G[g]=(D-F.T*A.inv()*F).applyfunc(sp.expand)
M3=sp.Matrix(10,10,lambda i,j: mom(P3[i][0]+P3[j][0],P3[i][1]+P3[j][1]))
assert all(M3[:k,:k].det()>0 for k in range(1,11))                # M3 PD
Bm=sp.Matrix(10,5,lambda i,j: L(P3[i][0]+D4[j][0],P3[i][1]+D4[j][1]))
C=(Bm.T*M3.inv()*Bm).applyfunc(sp.expand)
groups={}
for i in range(5):
    for j in range(i,5): groups.setdefault(i+j,[]).append((i,j))
h=[sp.expand(C[lst[0]]-C[q]) for s,lst in sorted(groups.items()) for q in lst[1:]]
assert len(h)==6
dat=pickle.load(open('sosdata.pkl','rb')); GT=dat['GT']; EQ=dat['EQ']; c=dat['c']; hw=dat['hw']; W=dat['W']
assert all(v!=0 for v in hw)
sub={y[k]:Q(c[k])+Q(hw[k])*t[k] for k in range(8)}
okGT=True
for g in gdef:
    Wg=sp.Matrix(4,4,lambda i,j: Q(W[g][i][j]))
    Gt=(Wg*G[g].subs(sub)*Wg.T).applyfunc(sp.expand)
    rec=sp.Matrix(4,4,lambda i,j: Q(GT[g][0][i][j]))+sum((t[k]*sp.Matrix(4,4,lambda i,j: Q(GT[g][1][k][i][j])) for k in range(8)),sp.zeros(4,4))
    okGT&=(Gt-rec).applyfunc(sp.expand)==sp.zeros(4,4)
print("GT == W * SchurComplement(M3(gL))(y=c+hw t) * W^T exactly, all g:",okGT)
def todict(expr):
    P=sp.Poly(sp.expand(expr),*t,domain='QQ')
    return {tuple(i for i in range(8) for _ in range(mon[i])):Fr(int(co.p),int(co.q)) for mon,co in P.terms()}
mine=[todict(hk.subs(sub)) for hk in h]
ok=True
for e,md in zip(EQ,mine):
    m0=next(iter(e)); r=e[m0]/md[m0]
    ok&= set(md)==set(e) and r!=0 and all(e[m]==r*md[m] for m in e)
print("EQ_k == (nonzero rational) * k-th Hankel condition on B^T M3^-1 B, all k:",ok)
print("DATA VERIFIED:",okGT and ok)
