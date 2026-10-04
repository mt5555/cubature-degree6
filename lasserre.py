# Moment (Lasserre) relaxation in t in [-1,1]^8, y = c + h*t
import numpy as np, sympy as sp, cvxpy as cp, scipy.sparse as sps, itertools, sys, time
from cf2 import *
from box import Gw_affine
lo,hi=np.load('box.npy'); c=(lo+hi)/2; hw=(hi-lo)/2
n=8
def monos(D): return [m for d in range(D+1) for m in itertools.combinations_with_replacement(range(n),d)]
def key(m): return tuple(sorted(m))
# whitened LMIs in t: Gw(t) = B0 + sum t_k Bk
aff=Gw_affine(); LMI={}
for g,(A0,Ak) in aff.items():
    B0=A0+sum(c[k]*Ak[k] for k in range(8)); Bk=[hw[k]*Ak[k] for k in range(8)]
    LMI[g]=(B0,Bk)
# hank in t (polynomial dicts), each scaled to max-coef 1
ts=sp.symbols('t0:8')
hs=hank_sym([sp.Float(c[k],30)+sp.Float(hw[k],30)*ts[k] for k in range(8)])
EQ=[]
for h in hs:
    P=sp.Poly(sp.expand(h),*ts); d={}
    for mon,co in P.terms():
        m=tuple(i for i in range(n) for _ in range(mon[i])); d[m]=float(co)
    s=max(abs(v) for v in d.values()); EQ.append({m:v/s for m,v in d.items()})
def build(dlev, s_margin=0.0, box=True, verbose=True):
    D=2*dlev; M=monos(D); idx={m:i for i,m in enumerate(M)}; N=len(M)
    z=cp.Variable(N); cons=[z[idx[()]]==1]
    def mat(rows, cols, entry):   # entry(a,b)-> list of (monomial, coef)
        r=[];cidx=[];v=[]
        nr=len(rows)
        for a in range(nr):
            for b in range(nr):
                for m,co in entry(rows[a],rows[b]):
                    r.append(a*nr+b); cidx.append(idx[key(m)]); v.append(co)
        A=sps.csr_matrix((v,(r,cidx)),shape=(nr*nr,N))
        X=cp.reshape(A@z,(nr,nr),order='C'); return X
    vd=monos(dlev); vd1=monos(dlev-1)
    cons.append(mat(vd,None,lambda a,b:[(a+b,1.0)])>>0)
    for g,(B0,Bk) in LMI.items():
        rows=[(i,m) for i in range(4) for m in vd1]
        def ent(ra,rb,B0=B0,Bk=Bk):
            (i,ma),(j,mb)=ra,rb; out=[(ma+mb,B0[i,j]-(s_margin if i==j else 0))]
            out+= [(ma+mb+(k,),Bk[k][i,j]) for k in range(8)]; return out
        cons.append(mat(rows,None,ent)>>0)
    if box:
        for k in range(n):
            cons.append(mat(vd1,None,lambda a,b,k=k:[(a+b,1.0),(a+b+(k,k),-1.0)])>>0)
    rowsE=[];colsE=[];valsE=[];r=0
    for e in EQ:
        for beta in monos(D-2):
            for m,co in e.items():
                rowsE.append(r); colsE.append(idx[key(m+beta)]); valsE.append(co)
            r+=1
    AE=sps.csr_matrix((valsE,(rowsE,colsE)),shape=(r,N)).toarray()
    # drop dependent rows (orthonormal basis of row space)
    U,sv,Vt=np.linalg.svd(AE,full_matrices=False); rk=int(np.sum(sv>1e-10*sv[0]))
    if verbose: print(f"  eq rows {r}, rank {rk}, smallest kept sv {sv[rk-1]:.2e}",flush=True)
    cons.append(Vt[:rk]@z==0)
    if verbose: print(f"level {dlev}: {N} moments, moment matrix {len(vd)}, loc blocks {4*len(vd1)}, {r} eq rows",flush=True)
    return cp.Problem(cp.Minimize(0),cons), z, idx
if __name__=="__main__":
    dlev=int(sys.argv[1]); smarg=[float(a) for a in sys.argv[2:]] or [0.0]
    for s in smarg:
        t0=time.time(); p,z,idx=build(dlev,s)
        try: p.solve(solver='CLARABEL',max_iter=500)
        except Exception as ex: print("solver error",ex)
        print(f"  margin s={s:+.3f}: status {p.status}  ({time.time()-t0:.1f}s)",flush=True)
