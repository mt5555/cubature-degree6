# re-verify min residual for 10 inside nodes with SLSQP + linear inequality constraints
import numpy as np, sys
from scipy.optimize import minimize
from common import resid
N=10; rng=np.random.default_rng(int(sys.argv[1])); out=[]
r0=np.load('res_10_inside_3.npy')
def start(t):
    if t < 40:   # from previous trf results
        z=r0[t,3:]; s=z[:30].reshape(N,3)**2; lam=s/s.sum(1,keepdims=True); return np.concatenate([lam[:,0],lam[:,1],z[30:]**2])
    lam=rng.dirichlet([1,1,1],N); return np.concatenate([lam[:,0],lam[:,1],rng.uniform(0.02,0.08,N)])
A=np.zeros((3*N+N,3*N)); b=np.zeros(4*N)
for i in range(N):
    A[i,i]=1; A[N+i,N+i]=1; A[2*N+i,i]=-1; A[2*N+i,N+i]=-1; b[2*N+i]=1; A[3*N+i,2*N+i]=1
cons=[{'type':'ineq','fun':lambda u: A@u+b,'jac':lambda u: A}]
f=lambda u: 0.5*np.sum(resid(u[:N],u[N:2*N],u[2*N:])**2)
for t in range(int(sys.argv[2])):
    r=minimize(f,start(t),method='SLSQP',constraints=cons,options={'ftol':1e-20,'maxiter':3000})
    u=r.x; viol=max(0,-(A@u+b).min())
    out.append((np.sqrt(2*f(u)),viol,t))
out.sort()
print("seed",sys.argv[1],"best:",["%.5f(viol %.0e,start %d)"%o for o in out[:6]])
print("distinct values:",np.unique(np.round([o[0] for o in out],4))[:8])
