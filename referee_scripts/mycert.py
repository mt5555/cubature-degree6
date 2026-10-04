# Independent exact re-check of the final certificate (own polynomial representation, own PD test)
import pickle, itertools, numpy as np
from fractions import Fraction as Fr
dat=pickle.load(open('sosdata.pkl','rb')); GT=dat['GT']; EQ=dat['EQ']
num=pickle.load(open('sos_num.pkl','rb')); rb=pickle.load(open('rigbox.pkl','rb')); B=rb['B']
print("rigbox bounds types/values:",{k:(type(v).__name__,float(v)) for k,v in rb['bounds'].items()})
print("B:",[float(b) for b in B], "all Fraction and >0:",all(type(b) is Fr and b>0 for b in B))
n=8
def monos(D): return [m for d in range(D+1) for m in itertools.combinations_with_replacement(range(n),d)]
v2=monos(2); v1=monos(1); assert len(v2)==45 and len(v1)==9
def ex(m):  # exponent vector
    e=[0]*n
    for i in m: e[i]+=1
    return tuple(e)
def addp(P,e,v): P[e]=P.get(e,Fr(0))+v
def mul(e1,e2): return tuple(a+b for a,b in zip(e1,e2))
delta=Fr(1,10**8)
S0=[[Fr(float(num['S0'][a,b]))/2+Fr(float(num['S0'][b,a]))/2+(delta if a==b else 0) for b in range(45)] for a in range(45)]
Sg={g:[[Fr(float(num['Sg'][g][a,b]))/2+Fr(float(num['Sg'][g][b,a]))/2 for b in range(36)] for a in range(36)] for g in 'xyh'}
lam=[Fr(float(x)) for x in num['lam']]; gam=Fr(float(num['gam']))
def cholesky_pd(M):
    # exact rational Cholesky-type test: returns True iff M (symmetric) is PD, via LDL^T pivots
    nn=len(M); A=[row[:] for row in M]; Lm=[[Fr(0)]*nn for _ in range(nn)]; Dv=[]
    for j in range(nn):
        d=A[j][j]-sum(Lm[j][k]**2*Dv[k] for k in range(j))
        if d<=0: return False,float(d)
        Dv.append(d); Lm[j][j]=Fr(1)
        for i in range(j+1,nn):
            Lm[i][j]=(A[i][j]-sum(Lm[i][k]*Lm[j][k]*Dv[k] for k in range(j)))/d
    return True,float(min(Dv))
print("S0+delta I symmetric:",all(S0[a][b]==S0[b][a] for a in range(45) for b in range(45)))
print("S0 PD (own LDL):",cholesky_pd(S0)," float min eig:",np.linalg.eigvalsh(np.array(S0,dtype=float)).min())
for g in 'xyh': print(f"S_{g} PD (own LDL):",cholesky_pd(Sg[g])," float min eig:",np.linalg.eigvalsh(np.array(Sg[g],dtype=float)).min())
P={}
for a in range(45):
    for b in range(45):
        if S0[a][b]: addp(P,mul(ex(v2[a]),ex(v2[b])),S0[a][b])
rows=[(i,m) for i in range(4) for m in v1]
for g in 'xyh':
    B0,Bk=GT[g]
    for A,(i,ma) in enumerate(rows):
        for Bb,(j,mb) in enumerate(rows):
            s=Sg[g][A][Bb]
            if s==0: continue
            e=mul(ex(ma),ex(mb))
            if B0[i][j]: addp(P,e,s*B0[i][j])
            for k in range(8):
                if Bk[k][i][j]: addp(P,mul(e,ex((k,))),s*Bk[k][i][j])
col=0
for e in EQ:
    for beta in v2:
        if lam[col]:
            for mm,co in e.items(): addp(P,mul(ex(mm),ex(beta)),lam[col]*co)
        col+=1
assert col==len(lam)
addp(P,ex(()),gam)
P={e:v for e,v in P.items() if v!=0}
bound=sum(abs(v)*np.prod([B[i]**e[i] for i in range(n)],dtype=object) for e,v in P.items())
print("monomials:",len(P)," max degree:",max(sum(e) for e in P))
print("gamma = %.6e   bound = %.6e   bound<gamma: %s"%(float(gam),float(bound),bound<gam))
# how much of the bound is the delta shift?
print("delta contribution: %.3e"%float(sum(delta*np.prod([B[i]**(2*ex(m)[i]) for i in range(n)],dtype=object) for m in v2)))
# sanity: evaluate P at random points in box (float) -- should be ~<=bound
rng=np.random.default_rng(1)
mx=0
for _ in range(2000):
    tt=rng.uniform(-1,1,8)*np.array([float(b) for b in B])
    val=sum(float(v)*np.prod(tt**np.array(e)) for e,v in P.items()); mx=max(mx,abs(val))
print("max |P| at 2000 random box points: %.3e"%mx)
