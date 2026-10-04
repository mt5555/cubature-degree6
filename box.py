# bounding box of the convex set K = {y : G_g(y) >= 0} (whitened LMIs)
import numpy as np, cvxpy as cp
from cf2 import *
def Gw_affine():
    # whitened G_g(y) = A0 + sum_k y_k A_k (exact affine), returned per g
    out={}
    base=G(np.zeros(8))
    for g in base:
        A0=Wm[g]@base[g]@Wm[g].T
        Ak=[]
        for k in range(8):
            e=np.zeros(8); e[k]=1; Ak.append(Wm[g]@(G(e)[g]-base[g])@Wm[g].T)
        out[g]=(A0,Ak)
    return out
if __name__=="__main__":
    aff=Gw_affine(); u=cp.Variable(8); y=cp.multiply(yLf,1+u)
    cons=[aff[g][0] + sum(y[k]*aff[g][1][k] for k in range(8)) >> 0 for g in aff]
    lo=[];hi=[]
    for k in range(8):
        for sgn,lst in [(1,lo),(-1,hi)]:
            p=cp.Problem(cp.Minimize(sgn*u[k]),cons); p.solve(solver='CLARABEL')
            lst.append(yLf[k]*(1+u.value[k])); print(k,sgn,p.status)
    lo=np.array(lo); hi=np.array(hi)
    np.save('box.npy',np.vstack([lo,hi]))
    for k in range(8): print(k,"lo %.6e hi %.6e  Leb %.6e  (rel width %.3f)"%(lo[k],hi[k],yLf[k],(hi[k]-lo[k])/yLf[k]))
