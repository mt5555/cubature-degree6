# P1: rigorous bounding box  |t_k| <= B_k on K  via exact LMI duality certificates
import numpy as np, cvxpy as cp, pickle, warnings
from fractions import Fraction as Fr
from exact import ldl_pd, solve_exact
warnings.filterwarnings('ignore')
GT=pickle.load(open('sosdata.pkl','rb'))['GT']; gl=['x','y','h']
B0f={g:np.array(GT[g][0],dtype=float) for g in gl}; Bkf={g:[np.array(GT[g][1][j],dtype=float) for j in range(8)] for g in gl}
def ip(A,B): return sum(A[i][j]*B[i][j] for i in range(4) for j in range(4))
Gam=[[sum(ip(GT[g][1][i],GT[g][1][j]) for g in gl) for j in range(8)] for i in range(8)]
bounds={}; Zstore={}
for k in range(8):
    for s in [1,-1]:
        Z={g:cp.Variable((4,4),PSD=True) for g in gl}; eps=cp.Variable()
        cons=[sum(cp.trace(Z[g]@Bkf[g][j]) for g in gl)==(s if j==k else 0) for j in range(8)]
        cons+=[sum(cp.trace(Z[g]@B0f[g]) for g in gl)<=1.03]+[Z[g]>>eps*np.eye(4) for g in gl]
        cp.Problem(cp.Maximize(eps),cons).solve(solver='CLARABEL')
        Zq={g:[[Fr(float((Z[g].value[i,j]+Z[g].value[j,i])/2)) for j in range(4)] for i in range(4)] for g in gl}
        r=[(Fr(s) if j==k else Fr(0))-sum(ip(Zq[g],GT[g][1][j]) for g in gl) for j in range(8)]
        al=solve_exact(Gam,r)
        Zc={g:[[Zq[g][i][j]+sum(al[l]*GT[g][1][l][i][j] for l in range(8)) for j in range(4)] for i in range(4)] for g in gl}
        assert all(sum(ip(Zc[g],GT[g][1][j]) for g in gl)==(s if j==k else 0) for j in range(8))
        pd=[ldl_pd(Zc[g]) for g in gl]; assert all(p[0] for p in pd), pd
        c0=sum(ip(Zc[g],GT[g][0]) for g in gl)    # s*t_k + c0 = sum <Z,Gt(t)> >= 0 on K
        bounds[(k,s)]=c0; Zstore[(k,s)]=Zc
        print(f"t_{k} {'>=' if s==1 else '<='} {'-' if s==1 else '+'}{float(c0):.6f}   (eps={eps.value:.2e}, min pivots {[f'{p[1]:.1e}' for p in pd]})",flush=True)
B=[max(bounds[(k,1)],bounds[(k,-1)]) for k in range(8)]
pickle.dump({'bounds':bounds,'B':B,'Z':Zstore},open('rigbox.pkl','wb'))
print("rigorous |t_k| bounds:",[f"{float(b):.6f}" for b in B])
