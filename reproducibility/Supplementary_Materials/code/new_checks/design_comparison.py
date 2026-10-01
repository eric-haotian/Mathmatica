"""Same-hypothesis-pair, same-total-effort information comparison.
Balanced finite designs are compared with a pair-informed frequency oracle.
This does not simulate an estimator or implement an adaptive sampling strategy.
"""
from pathlib import Path
from fractions import Fraction as F
from math import factorial,comb
import json
ROOT=Path(__file__).resolve().parent

def pack(x):
    if isinstance(x,F): return str(x)
    if isinstance(x,dict): return {k:pack(v) for k,v in x.items()}
    if isinstance(x,(tuple,list)): return [pack(v) for v in x]
    return x

def numerator(n,h): return (F(factorial(n),2**(n-1))*(h/2)**n)**2

def qprod(ts,v):
    ans=F(1)
    for t in ts: ans*=1+t*t*v
    return ans

def direct_gap(atoms0,atoms1,w):
    def obs(atoms):
        a=b=F(0)
        for x,m in atoms:
            t=1+x/2; den=1+w*w*t*t
            a+=m/den; b-=m*w*t/den
        return a,b
    p=obs(atoms0); q=obs(atoms1)
    return sum((a-b)**2 for a,b in zip(p,q))

def create(r,h):
    n=2*r;xs=[F(l-r)*h for l in range(n+1)];ts=[1+x/2 for x in xs]
    mu0=[(xs[l],F(comb(n,l),2**(n-1))) for l in range(n+1) if l%2]
    mu1=[(xs[l],F(comb(n,l),2**(n-1))) for l in range(n+1) if l%2==0]
    L=F(n)/max(ts)**2;U=F(n)/min(ts)**2
    phi=lambda v:sum(t*t*v/(1+t*t*v) for t in ts)
    for _ in range(110):
        mid=(L+U)/2
        if phi(mid)<n:L=mid
        else:U=mid
    C=numerator(n,h);mid=(L+U)/2
    lo=C*mid**n/qprod(ts,mid);hi=C*U**n/qprod(ts,L)
    result={'format':'same_effort_information_v1','r':r,'h':h,'center':F(0),
            'reference':mu0,'alternative':mu1,'v_lo':L,'v_hi':U,
            'max_gap_square_lower':lo,'max_gap_square_upper':hi,'designs':{}}
    for name,ws in [('baseline',list(map(F,range(1,7)))),('high_band',[F(16*i) for i in range(1,7)])]:
        vals=[direct_gap(mu0,mu1,w) for w in ws];A=sum(vals)/len(vals)
        cn=F(factorial(n)**2,2**(4*n-2))
        limfinite=cn*sum(w**(2*n)/(1+w*w)**(n+1) for w in ws)/len(ws)
        limoracle=cn*F(n**n,(1+n)**(n+1))
        result['designs'][name]={'frequencies':ws,'point_gap_squares':vals,'mean_gap_square':A,
            'oracle_efficiency_lower':A/hi,'oracle_efficiency_upper':A/lo,
            'limit_efficiency':limfinite/limoracle}
    key=f'information_r{r}_h{h.denominator}'
    (ROOT/'certificates'/f'{key}.json').write_text(json.dumps(pack(result),indent=2))
    return result

if __name__=='__main__':
    rows=[]
    for r in range(1,5):
        for h in [F(1,8),F(1,16),F(1,32)]:
            a=create(r,h);rows.append(a)
            print(r,str(h),*[float(a['designs'][k]['oracle_efficiency_lower']) for k in ['baseline','high_band']])
    (ROOT/'results/design_comparison.json').write_text(json.dumps(pack(rows),indent=2))
