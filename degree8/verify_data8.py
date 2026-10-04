# Independent exact re-derivation for degree 8 (direct 15x15 moment matrices, no block formulas):
# (a) M4 PD; the 10 Hankel conditions on C = B^T M4^{-1} B equal the msolve input equations up to nonzero scalars;
# (b) Schur complements of the 15x15 localizing matrices on the P3 block equal cf8's G_g = D_g - K_g, with A_g PD.
import sympy as sp
from math import factorial
Q=sp.Rational
def mom(a,b): return Q(factorial(a)*factorial(b),factorial(a+b+2))
y=sp.symbols('y0:10')
def L(a,b):
    assert a+b<=9; return y[b] if a+b==9 else mom(a,b)
P4=[(a,d-a) for d in range(5) for a in range(d,-1,-1)]; D5=[(5-j,j) for j in range(6)]
M4=sp.Matrix(15,15,lambda i,j: mom(P4[i][0]+P4[j][0],P4[i][1]+P4[j][1]))
assert all(M4[:k,:k].det()>0 for k in range(1,16)); print("M4 positive definite: True")
B=sp.Matrix(15,6,lambda i,j: L(P4[i][0]+D5[j][0],P4[i][1]+D5[j][1]))
C=(B.T*M4.inv()*B).applyfunc(sp.expand)
groups={}
for i in range(6):
    for j in range(i,6): groups.setdefault(i+j,[]).append((i,j))
h=[sp.expand(C[l[0]]-C[q]) for s,l in sorted(groups.items()) for q in l[1:]]
assert len(h)==10
lines=open('flat8_Q.ms').read().split('\n')
E=[sp.Poly(sp.sympify(l.rstrip(',').replace('^','**')),*y) for l in lines[2:] if l.strip()]
ok=True
for hk,Ek in zip(h,E):
    Hp=sp.Poly(hk,*y); m0=Ek.monoms()[0]; r=Ek.coeff_monomial(m0)/Hp.coeff_monomial(m0)
    ok&= r!=0 and (Ek-Hp*r).is_zero
print("msolve input equations == nonzero rational multiples of the Hankel conditions:",ok)
import cf8
gdef={'x':{(1,0):1},'y':{(0,1):1},'h':{(0,0):1,(1,0):-1,(0,1):-1}}
okG=True
for g,gd in gdef.items():
    M=sp.Matrix(15,15,lambda i,j: sum(c*L(P4[i][0]+P4[j][0]+da,P4[i][1]+P4[j][1]+db) for (da,db),c in gd.items()))
    A=M[:10,:10]; F=M[:10,10:]; D=M[10:,10:]
    assert A.free_symbols==set() and F.free_symbols==set()
    assert all(A[:k,:k].det()>0 for k in range(1,11))
    okG&= ((D-F.T*A.inv()*F)-cf8.G_sym(list(y))[g]).applyfunc(sp.expand)==sp.zeros(5,5)
print("Schur complements of M4(gL) on the P3 block == cf8 G_g (A_g PD), all g:",okG)
