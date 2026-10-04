# Facially reduced SOS: S0 = U^T S' U, U rational basis of span{n_p}^perp
import numpy as np, cvxpy as cp, scipy.sparse as sps, pickle, warnings
from fractions import Fraction as Fr
warnings.filterwarnings('ignore')
from sos import *          # A0, Ag, Ah, nl, v2, I4, NC
dat=pickle.load(open('sosdata.pkl','rb')); hw=dat['hw']
low=[a for a,m in enumerate(v2) if len(m)<2]
groups={}
for a,m in enumerate(v2):
    if len(m)==2: groups.setdefault(m[0]+m[1],[]).append(a)
nval=lambda a: Fr(1)/(hw[v2[a][0]]*hw[v2[a][1]])
Urows=[]
for a in low:
    r=[Fr(0)]*45; r[a]=Fr(1); Urows.append(r)
for p,mem in sorted(groups.items()):
    a1=mem[0]
    for al in mem[1:]:
        r=[Fr(0)]*45; sc=max(abs(nval(al)),abs(nval(a1))); r[a1]=nval(al)/sc; r[al]=-nval(a1)/sc; Urows.append(r)
Uq=Urows; U=np.array([[float(x) for x in r] for r in Uq]); nU=len(Uq)
print("reduced Gram size:",nU)
# map S' -> coefficients: A0 vec(U^T S' U)
kron=np.kron(U.T,U.T)   # vec_F(U^T S U) = (U^T kron U^T) vec_F(S)
A0r=sps.csr_matrix(A0@kron)
def solve2(stage,gam_fix=None):
    S1=cp.Variable((nU,nU),PSD=True); Sg={g:cp.Variable((Ag[g][1],Ag[g][1]),PSD=True) for g in Ag}
    lam=cp.Variable(nl); gam=cp.Variable(); tau=cp.Variable()
    e0=np.zeros(NC); e0[I4[()]]=1
    expr=A0r@cp.vec(S1,order='F')+sum(Ag[g][0]@cp.vec(Sg[g],order='F') for g in Ag)+Ah@lam+gam*e0
    cons=[expr==0, cp.trace(S1)+sum(cp.trace(Sg[g]) for g in Sg)<=1]
    if stage=='a': obj=cp.Maximize(gam)
    else:
        cons+=[gam==gam_fix, S1>>tau*np.eye(nU)]+[Sg[g]>>tau*np.eye(Ag[g][1]) for g in Sg]; obj=cp.Maximize(tau)
    p=cp.Problem(obj,cons); p.solve(solver='CLARABEL')
    return p,S1.value,{g:Sg[g].value for g in Sg},lam.value,gam.value,tau.value
if __name__=="__main__":
    p,S1,Sg,lam,gam,_=solve2('a'); print("stage a:",p.status,"gamma*=%.6e"%gam)
    p,S1,Sg,lam,g2,tau=solve2('b',gam/2); print("stage b:",p.status,"tau=%.4e"%tau)
    print("min eigs: S' %.3e"%np.linalg.eigvalsh(S1).min(),{g:"%.3e"%np.linalg.eigvalsh(Sg[g]).min() for g in Sg})
    M=sps.hstack([A0r]+[Ag[g][0] for g in Ag]+[Ah]).toarray()
    print("rank of full coefficient map:",np.linalg.matrix_rank(M,tol=1e-9),"of",NC)
    Mb=sps.hstack([A0r,Ah]).toarray(); print("rank of [A0r, Ah] (vars solved exactly):",np.linalg.matrix_rank(Mb,tol=1e-9))
    pickle.dump({'S1':S1,'Sg':Sg,'lam':lam,'gam':g2,'Uq':Uq},open('sos2_num.pkl','wb'))
