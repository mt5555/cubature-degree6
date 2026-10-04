# Find: sigma0 + sum_g <S_g, Gt_g(t) (x) v1 v1^T> + sum_k lam_k(t) h_k(t) = -gamma   (degree 4)
import numpy as np, cvxpy as cp, scipy.sparse as sps, pickle, sys, warnings
warnings.filterwarnings('ignore')
from polyutil import monos, key, n
dat=pickle.load(open('sosdata.pkl','rb')); GT=dat['GT']; EQ=dat['EQ']
v2=monos(2); v1=monos(1); V4=monos(4); I4={m:i for i,m in enumerate(V4)}; NC=len(V4)
def gram_map(basis):
    m=len(basis); r=[];cc=[];v=[]
    for a in range(m):
        for b in range(m):
            r.append(I4[key(basis[a]+basis[b])]); cc.append(a+b*m); v.append(1.0)
    return sps.csr_matrix((v,(r,cc)),shape=(NC,m*m))
def loc_map(B0,Bk):
    rows=[(i,mm) for i in range(4) for mm in v1]; m=len(rows); r=[];cc=[];v=[]
    for A,(i,ma) in enumerate(rows):
        for Bb,(j,mb) in enumerate(rows):
            col=A+Bb*m
            if B0[i][j]!=0: r.append(I4[key(ma+mb)]); cc.append(col); v.append(float(B0[i][j]))
            for k in range(8):
                if Bk[k][i][j]!=0: r.append(I4[key(ma+mb+(k,))]); cc.append(col); v.append(float(Bk[k][i][j]))
    return sps.csr_matrix((v,(r,cc)),shape=(NC,m*m)), m
def eq_map():
    r=[];cc=[];v=[];col=0
    for e in EQ:
        for beta in v2:
            for mm,co in e.items(): r.append(I4[key(mm+beta)]); cc.append(col); v.append(float(co))
            col+=1
    return sps.csr_matrix((v,(r,cc)),shape=(NC,col)), col
A0=gram_map(v2); Ag={}; 
for g in GT: Ag[g]=loc_map(*GT[g])
Ah,nl=eq_map()
def solve(stage, gam_fix=None, solver='CLARABEL'):
    S0=cp.Variable((45,45),PSD=True); Sg={g:cp.Variable((Ag[g][1],Ag[g][1]),PSD=True) for g in Ag}
    lam=cp.Variable(nl); gam=cp.Variable(); tau=cp.Variable()
    e0=np.zeros(NC); e0[I4[()]]=1
    expr=A0@cp.vec(S0,order='F')+sum(Ag[g][0]@cp.vec(Sg[g],order='F') for g in Ag)+Ah@lam+gam*e0
    cons=[expr==0, cp.trace(S0)+sum(cp.trace(Sg[g]) for g in Sg)<=1]
    if stage=='a': obj=cp.Maximize(gam)
    else:
        cons+=[gam==gam_fix, S0>>tau*np.eye(45)]+[Sg[g]>>tau*np.eye(Ag[g][1]) for g in Sg]; obj=cp.Maximize(tau)
    p=cp.Problem(obj,cons); p.solve(solver=solver,verbose=False)
    return p, S0.value, {g:Sg[g].value for g in Sg}, lam.value, gam.value, tau.value
if __name__=="__main__":
    for solver in ['CLARABEL','SCS']:
        try:
            p,S0,Sg,lam,gam,_=solve('a',solver=solver); print(solver,"stage a:",p.status,"gamma* = %.6e"%gam); break
        except Exception as ex: print(solver,"failed",ex)
    p,S0,Sg,lam,gam,tau=solve('b',gam_fix=gam/2,solver=solver)
    print("stage b:",p.status,"gamma=%.4e tau=%.4e"%(gam,tau))
    print("min eigs: S0 %.3e"%np.linalg.eigvalsh(S0).min(), {g:"%.3e"%np.linalg.eigvalsh(Sg[g]).min() for g in Sg})
    pickle.dump({'S0':S0,'Sg':Sg,'lam':lam,'gam':gam},open('sos_num.pkl','wb'))
