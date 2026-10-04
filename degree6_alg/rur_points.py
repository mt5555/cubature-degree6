# Real points of an msolve parametrization (-P 1 output): take msolve's isolating interval for the
# parametrizing variable t, refine t by Newton on f at high precision, evaluate all coordinates.
import sys, mpmath as mp, numpy as np
from parse_rur import load
sys.path.insert(0,'..')
from cf2 import yLf, geigs
mp.mp.dps=300
def points(fn):
    R=load(fn); names=R[1][3]; P=R[1][5][1]
    # the parametrizing variable t is the last listed variable (msolve may add a variable A for genericity)
    fc=[mp.mpf(c) for c in reversed(P[0][1])]; dfc=[c*(len(fc)-1-i) for i,c in enumerate(fc[:-1])]
    qc=[mp.mpf(c) for c in reversed(P[1][1])]; nums=[([mp.mpf(c) for c in reversed(a[0][1])],a[1]) for a in P[2]]
    sols=R[2][1] if len(R)>2 else []
    ts=[]
    if sols and len(sols[0])==len(names):          # msolve gave intervals for every variable incl. t
        for so in sols:
            lo,hi=[mp.mpf(z.numerator)/z.denominator for z in so[-1]]; t=(lo+hi)/2
            for _ in range(200): t=t-mp.polyval(fc,t)/mp.polyval(dfc,t)
            assert lo<=t<=hi, "Newton left the isolating interval"
            ts.append(t)
    else:                                            # t is an added variable: real roots of f directly
        ts=[mp.re(r) for r in mp.polyroots(fc,maxsteps=500,extraprec=5000) if abs(mp.im(r))<mp.mpf(10)**-100]
        assert len(ts)==len(sols), (len(ts),len(sols))
    pts=[]
    for t in ts:
        vals={names[i]:-mp.polyval(n,t)/(c*mp.polyval(qc,t)) for i,(n,c) in enumerate(nums)}
        vals[names[-1]]=t; pts.append(vals)
    return pts
if __name__=="__main__":
    for fn in sys.argv[1:]:
        pts=points(fn); print(f"{fn}: {len(pts)} real points")
        for v in pts:
            y=np.array([float(v[f'y{i}']) for i in range(8)]); dev=np.abs(y/yLf-1).max()
            me=min(e.min() for e in geigs(y).values()) if dev<1 else float('nan')
            print(f"   max|y/yL-1| = {dev:.3e}   min whitened eig = {me:9.3f}   in K: {bool(dev<0.07 and me>=0)}")
