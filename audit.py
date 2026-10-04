import numpy as np, pickle, scipy.linalg as sl
from fractions import Fraction as Fr
import cf, cf2
from sosdata import monos, key
rng=np.random.default_rng(0)
# (a) Schur 4x4 form vs direct 10x10 localizing matrices: generalized eigenvalues vs Lebesgue must match
y7L=cf2.yLf; ML=cf.locs(y7L)
worst=0
for trial in range(20):
    y=y7L*(1+0.03*rng.standard_normal(8))
    direct=[np.sort(sl.eigh(Mg,MLg,eigvals_only=True)) for Mg,MLg in zip(cf.locs(y),ML)]
    schur=cf2.geigs(y)
    for d,g in zip(direct,['x','y','h']):
        # direct has 6 eigenvalues equal to 1 (Pi_2 block) plus the 4 Schur ones
        s=np.sort(np.concatenate([schur[g],np.ones(6)])); worst=max(worst,np.abs(s-d).max())
print("(a) max |gen.eig(10x10) - gen.eig(Schur 4x4 + six 1's)| over 20 random y: %.2e"%worst)
# also hank vs cf.hank (independent implementations) at random y
print("    hank: |cf2.hank - cf.hank-based| :", end=" ")
mx=0
for trial in range(20):
    y=y7L*(1+0.03*rng.standard_normal(8))
    C=cf.C_of(y); hk=np.array([C[p]-C[q] for p,q in cf.pairs])
    mx=max(mx,np.abs(np.sort(np.abs(hk))-np.sort(np.abs(cf2.hank(y)))).max()/np.abs(hk).max())
print("%.2e (relative, as sets)"%mx)
# (b) exact data vs float formulation at the extremal outside rule
dat=pickle.load(open('sosdata.pkl','rb')); GT=dat['GT']; EQ=dat['EQ']; c=dat['c']; hw=dat['hw']
u=np.load('extremal.npy'); N=10; x,yy,w=u[:N],u[N:2*N],u[2*N:3*N]
y7=np.array([np.sum(w*x**(7-k)*yy**k) for k in range(8)]); t=np.array([(y7[k]-float(c[k]))/float(hw[k]) for k in range(8)])
Gt={g:np.array(GT[g][0],dtype=float)+sum(t[k]*np.array(GT[g][1][k],dtype=float) for k in range(8)) for g in GT}
Wn={g:cf2.Wm[g]@cf2.G(y7)[g]@cf2.Wm[g].T for g in GT}
print("(b) max |Gtilde_exact(t) - W G(y) W^T| at extremal rule: %.2e"%max(np.abs(Gt[g]-Wn[g]).max() for g in GT))
print("    exact h~_k at extremal rule:",["%.1e"%sum(float(co)*np.prod([t[i] for i in m]) for m,co in e.items()) for e in EQ])
# (c) negative control: evaluate certificate pieces at the extremal (outside) rule
num=pickle.load(open('sos_num.pkl','rb')); v2=monos(2); v1=monos(1)
vv2=np.array([np.prod([t[i] for i in m]) for m in v2]); vv1=np.array([np.prod([t[i] for i in m]) for m in v1])
sig=vv2@(num['S0']+1e-8*np.eye(45))@vv2
loc={g:np.sum(num['Sg'][g]*np.kron(Gt[g],np.outer(vv1,vv1))) for g in GT}
print("(c) at extremal outside rule: sigma0=%.3e, loc terms="%sig,{g:"%.3e"%v for g,v in loc.items()},
      " sum+gamma=%.2e (should be ~0)"%(sig+sum(loc.values())+num['gam']))
print("    min eig of Gtilde_h there: %.3f -> the h-term is what goes negative, as expected"%np.linalg.eigvalsh(Gt['h']).min())
