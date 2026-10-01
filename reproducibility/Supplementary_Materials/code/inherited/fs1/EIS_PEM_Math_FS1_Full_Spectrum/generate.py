"""Generate exact analytic test pairs; this does NOT run a spectrum estimator."""
from pathlib import Path
from fractions import Fraction as F
from math import comb, factorial
import json
from exact import *
ROOT=Path(__file__).resolve().parent

def generate_pair(r,h,c,label,N=2000000):
    n=2*r;nodes=[c+(l-r)*h for l in range(n+1)];ts=[1+x/2 for x in nodes]
    ref=[(nodes[l],F(comb(n,l),2**(n-1))) for l in range(n+1) if l%2]
    alt=[(nodes[l],F(comb(n,l),2**(n-1))) for l in range(n+1) if not l%2]
    C=(F(factorial(n),2**(2*n-1))*h**n)**2
    A=lambda v:sum(t*t*v/(1+t*t*v) for t in ts)
    L=F(n)/max(ts)**2;U=F(n)/min(ts)**2
    for _ in range(90):
        mid=(L+U)/2
        if A(mid)<n:L=mid
        else:U=mid
    B=lambda v:C*v**n/product(1+t*t*v for t in ts)
    sqlo=B((L+U)/2);squp=C*U**n/product(1+t*t*L for t in ts)
    lo=sqrt_bracket(sqlo,220)[0];eps=sqrt_bracket(squp,220)[1]
    flo=sqrt_bracket(L,150)[0];fhi=sqrt_bracket(U,150)[1]
    theta=[F(1,10),F(1,100)];freq=[1,2,3,4,5,6];y=[]
    for om in freq:
        y.extend([theta[0]+sum(w/(1+om*om*(1+x/2)**2) for x,w in ref),
                  -om*theta[1]+sum(w*om*(1+x/2)/(1+om*om*(1+x/2)**2) for x,w in ref)])
    tg=[dict(s=s,normalizer=factorial(s),sign=(-1)**r) for s in [2,3,4]]
    gaps={}
    for s in [2,3,4]:
        fun=lambda x:F((-1)**r,factorial(s))*h**(s+2)/(h*h+(x-c)**2)
        gaps[str(s)]=sum(w*fun(x) for x,w in alt)-sum(w*fun(x) for x,w in ref)
    st=dict(sigma=F(1,1000),alpha=F(1,20),n=N)
    inp=dict(format='FS1_positive_pair_v1',r=r,h=h,c=c,reference=ref,alternative=alt,theta=theta,
             frequencies=freq,y=y,epsilon=eps,targets=tg,statistics=st)
    kl=N*eps**2/(2*st['sigma']**2)
    conf={s:(1-2*st['alpha']-F(1,4))*g for s,g in gaps.items()} if kl<=F(1,8) else {}
    cert=dict(input_hash=digest(inp),v_lo=L,v_hi=U,norm_square_lo=sqlo,norm_square_hi=squp,
              omega_lo=flo,omega_hi=fhi,norm_lo=lo,norm_hi=eps,target_gaps=gaps,
              kl_upper=kl,kl_at_most_eighth=kl<=F(1,8),honest_expected_length_lower=conf)
    dump(ROOT/'inputs'/f'{label}.json',inp);dump(ROOT/'certificates'/f'{label}.json',cert)

def reproduction(m):
    freq=list(range(1,m+1));D=[F(1)];factors=[[F(1),F(0),F(w*w)] for w in freq]
    for f in factors:D=mul(D,f)
    D1=val(D,F(1));nums=[]
    for i,w in enumerate(freq):
        p=[F(1)]
        for j,f in enumerate(factors):
            if j!=i:p=mul(p,f)
        nums.extend([scale(p,1/D1),scale(mul(p,[F(0),F(w)]),1/D1)])
    A=[[p[i] if i<len(p) else F(0) for p in nums] for i in range(2*m)]
    lam=[]
    for j in range(2*(m-2)+1):
        rhs=mul([F(0),F(0),F(1)],power([F(-2),F(2)],j));rhs+=[F(0)]*(2*m-len(rhs))
        lam.append(solve(A,rhs))
    dump(ROOT/'inputs'/f'model_m{m}.json',dict(format='FS1_reproduction_v1',frequencies=freq,**{'lambda':lam}))

if __name__=='__main__':
    for r in range(1,5):
        for j in [3,4,5,6,7]:generate_pair(r,F(1,2**j),F(0),f'r{r}_h{j}')
    for c,tag in [(F(-1,4),'minus'),(F(1,4),'plus')]:generate_pair(3,F(1,16),c,f'r3_shift_{tag}')
    for m in [3,4,5,6]:reproduction(m)
    print('Generated 22 positive-pair certificates and 4 exact operator maps.')
