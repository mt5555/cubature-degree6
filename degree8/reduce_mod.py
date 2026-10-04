import sys, re
src, p, dst = sys.argv[1], int(sys.argv[2]), sys.argv[3]
L=open(src).read().split('\n')
body='\n'.join(L[2:])
# reduce every integer coefficient (numbers not preceded by 'y' or '^') mod p
body=re.sub(r'(?<![y\^\d])(\d+)', lambda m: str(int(m.group(1)) % p), body)
open(dst,'w').write(L[0]+'\n'+str(p)+'\n'+body)
