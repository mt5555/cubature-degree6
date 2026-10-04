import pickle
from fractions import Fraction as Fr
dat=pickle.load(open('sosdata.pkl','rb')); GT=dat['GT']
rb=pickle.load(open('rigbox_dump.pkl','rb')); Z=rb['Z']; bounds=rb['bounds']
def ldl(M):
    n=len(M); A=[r[:] for r in M]; L=[[Fr(0)]*n for _ in range(n)]; D=[]
    for j in range(n):
        d=A[j][j]-sum(L[j][k]**2*D[k] for k in range(j))
        if d<=0: return False
        D.append(d)
        for i in range(j+1,n): L[i][j]=(A[i][j]-sum(L[i][k]*L[j][k]*D[k] for k in range(j)))/d
    return True
ok=True
for (k,s),Zg in Z.items():
    # identity: sum_g <Z_g, B0_g + sum_j t_j Bk_g[j]> == s*t_k + c0  (coefficientwise, exact)
    lin=[sum(sum(Zg[g][i][j]*GT[g][1][m][i][j] for i in range(4) for j in range(4)) for g in 'xyh') for m in range(8)]
    c0=sum(sum(Zg[g][i][j]*GT[g][0][i][j] for i in range(4) for j in range(4)) for g in 'xyh')
    sym=all(Zg[g][i][j]==Zg[g][j][i] for g in 'xyh' for i in range(4) for j in range(4))
    pd=all(ldl(Zg[g]) for g in 'xyh')
    good= all(lin[m]==(s if m==k else 0) for m in range(8)) and c0==bounds[(k,s)] and sym and pd and all(type(v) is Fr for g in 'xyh' for r in Zg[g] for v in r)
    ok&=good
    print(k,s,"identity exact:",all(lin[m]==(s if m==k else 0) for m in range(8)),"c0 match:",c0==bounds[(k,s)],"Z sym:",sym,"Z PD:",pd,"c0=%.9f"%float(c0))
print("ALL BOX CERTIFICATES VERIFIED INDEPENDENTLY:",ok)
print("B (rerun) =",[float(b) for b in rb['B']])
