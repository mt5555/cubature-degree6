# Standalone exact verification of the infeasibility certificate (no SDP solving, no float arithmetic).
# Inputs: sosdata.pkl (exact data GT, EQ), rigbox.pkl (box certificates Z), sos_num.pkl (certificate).
# Step 1: re-verify the 16 box identities  s*t_k + c0 = sum_g <Z_g, Gt_g(t)>,  Z_g exactly PD  =>  |t_k| <= B_k on K.
# Step 2: P(t) := sigma0'(t) + sum_g <S_g, Gt_g(t) (x) v1 v1^T> + sum lam_k(t) h_k(t) + gamma.
#   On K ∩ V:  sigma0'>=0, <S_g,.> >= 0, h=0  =>  P(t*) >= gamma.  If sup_box |P| < gamma, K ∩ V is empty.
import pickle, time
from fractions import Fraction as Fr
from exact import ldl_pd
from polyutil import monos, key
t0=time.time()
dat=pickle.load(open('sosdata.pkl','rb')); GT=dat['GT']; EQ=dat['EQ']
box=pickle.load(open('rigbox.pkl','rb')); num=pickle.load(open('sos_num.pkl','rb'))
gl=['x','y','h']
def ip(A,B): return sum(A[i][j]*B[i][j] for i in range(4) for j in range(4))
# ---- Step 1: box certificates
assert all(type(v) is Fr for g in gl for M in [GT[g][0]]+GT[g][1] for r in M for v in r)
c0={}
for (k,s),Z in box['Z'].items():
    assert all(type(v) is Fr for g in gl for r in Z[g] for v in r)
    assert all(Z[g][i][j]==Z[g][j][i] for g in gl for i in range(4) for j in range(4))
    assert all(sum(ip(Z[g],GT[g][1][j]) for g in gl)==(s if j==k else 0) for j in range(8)), (k,s)
    assert all(ldl_pd(Z[g])[0] for g in gl), (k,s)
    c0[(k,s)]=sum(ip(Z[g],GT[g][0]) for g in gl)
B=[max(c0[(k,1)],c0[(k,-1)]) for k in range(8)]
assert len(c0)==16 and all(b>0 for b in B)
print(f"Step 1: 16 box certificates verified exactly; max |t_k| bound = {max(float(b) for b in B):.6f}  [{time.time()-t0:.1f}s]")
# ---- Step 2: main certificate
v2=monos(2); v1=monos(1)
delta=Fr(1,10**8)
q=lambda x: Fr(float(x))     # exact conversion of each stored float to a rational constant
S0=[[ (q(num['S0'][a,b])+q(num['S0'][b,a]))/2 + (delta if a==b else 0) for b in range(45)] for a in range(45)]
Sg={g:[[ (q(num['Sg'][g][a,b])+q(num['Sg'][g][b,a]))/2 for b in range(36)] for a in range(36)] for g in gl}
lam=[q(x) for x in num['lam']]; gam=q(num['gam'])
ok0,p0=ldl_pd(S0); assert ok0
print(f"Step 2: S0+delta*I exactly PD (min pivot {p0:.2e})",end="; ")
for g in gl:
    ok,p=ldl_pd(Sg[g]); assert ok; print(f"S_{g} PD (min pivot {p:.2e})",end="; ")
print()
P={}
def add(m,v):
    m=key(m); P[m]=P.get(m,Fr(0))+v
for a in range(45):
    for b in range(45):
        if S0[a][b]: add(v2[a]+v2[b],S0[a][b])
rows=[(i,mm) for i in range(4) for mm in v1]          # Kronecker order: (matrix index i, monomial mm)
for g in gl:
    B0,Bk=GT[g]
    for A,(i,ma) in enumerate(rows):
        for Bb,(j,mb) in enumerate(rows):
            s=Sg[g][A][Bb]
            if not s: continue
            add(ma+mb, s*B0[i][j])
            for k in range(8):
                if Bk[k][i][j]: add(ma+mb+(k,), s*Bk[k][i][j])
col=0
for e in EQ:
    for beta in v2:
        if lam[col]:
            for mm,co in e.items(): add(mm+beta, lam[col]*co)
        col+=1
add((),gam)
assert all(type(v) is Fr for v in P.values())
bound=Fr(0)
for m,v in P.items():
    w=abs(v)
    for i in m: w*=B[i]
    bound+=w
print(f"gamma = {float(gam):.6e};  sup over certified box of |P(t)| <= {float(bound):.6e}  ({len(P)} monomials)")
print("CERTIFICATE VALID (bound < gamma):", bound<gam, f"  [{time.time()-t0:.1f}s]")
