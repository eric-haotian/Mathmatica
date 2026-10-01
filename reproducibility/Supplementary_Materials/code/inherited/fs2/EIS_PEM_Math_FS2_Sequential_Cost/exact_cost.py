"""Exact rational certificates for stopped precision-cost lower bounds.

This module does not simulate experiments.  The statistical conclusion uses
Theorem 2 of Proof_FS2; the functions below check its numerical premises.
"""
from __future__ import annotations
from fractions import Fraction as F
from math import comb, factorial
import hashlib, json
from pathlib import Path

def require(ok: bool, message: str) -> None:
    if not ok:
        raise ValueError(message)

def rational(x) -> F:
    require(isinstance(x, (str, int, F)) and not isinstance(x, bool), 'exact rational input required')
    return F(x)

def canonical(x):
    if isinstance(x, F): return str(x)
    if isinstance(x, dict): return {k: canonical(v) for k,v in x.items()}
    if isinstance(x, (list,tuple)): return [canonical(v) for v in x]
    return x

def save(path, obj):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(canonical(obj),indent=2,sort_keys=True)+'\n')

def digest(obj):
    return hashlib.sha256(json.dumps(canonical(obj),sort_keys=True,separators=(',',':')).encode()).hexdigest()

def trim(p):
    p=list(p)
    while len(p)>1 and p[-1]==0:p.pop()
    return p

def add(p,q):
    return trim([(p[i] if i<len(p) else F(0))+(q[i] if i<len(q) else F(0)) for i in range(max(len(p),len(q)))])

def mul(p,q):
    out=[F(0)]*(len(p)+len(q)-1)
    for i,a in enumerate(p):
        for j,b in enumerate(q):out[i+j]+=a*b
    return trim(out)

def evaluate(p,x):
    ans=F(0)
    for a in reversed(p):ans=ans*x+a
    return ans

def _log_unit(x: F, terms: int):
    require(1<=x<=2 and terms>=1, 'log series domain')
    z=(x-1)/(x+1);z2=z*z;power=z;low=F(0)
    for j in range(terms):
        low+=2*power/(2*j+1);power*=z2
    tail=2*power/((2*terms+1)*(1-z2))
    return low,low+tail

def log_interval(x, terms=64):
    """Enclose log(x) by atanh series; only integer/rational arithmetic."""
    x=rational(x);require(x>0,'log argument must be positive')
    if x<1:
        lo,hi=log_interval(1/x,terms);return -hi,-lo
    n=0
    while x>=2:x/=2;n+=1
    lo,hi=_log_unit(x,terms);l2,u2=_log_unit(F(2),terms)
    return lo+n*l2,hi+n*u2

def binary_information(alpha):
    alpha=rational(alpha);require(0<alpha<F(1,4),'confidence error must be in (0,1/4)')
    l,u=log_interval((1-alpha)/alpha)
    return (1-2*alpha)*l,(1-2*alpha)*u

def pair_spec(r,h,c=0):
    require(type(r) is int and 1<=r<=8, 'r integer in [1,8]')
    h=rational(h);c=rational(c);n=2*r
    nodes=[c+(j-r)*h for j in range(n+1)]
    require(h>0 and min(nodes)>=-1 and max(nodes)<=1,'positive scale with support on [-1,1]')
    weights=[F(comb(n,j),2**(n-1)) for j in range(n+1)]
    return nodes,weights

def full_norm(r,h,c=0,bisections=104):
    """Prove the resolvent numerator identity and enclose the all-frequency norm."""
    nodes,weights=pair_spec(r,h,c);h=rational(h);n=2*r;ts=[1+x/2 for x in nodes]
    # Reconstruct response numerator in z = i omega, without using the claimed identity.
    num=[F(0)]
    for j,w in enumerate(weights):
        p=[F(1)]
        for k,t in enumerate(ts):
            if k!=j:p=mul(p,[F(1),t])
        num=add(num,[(-1)**j*w*a for a in p])
    coefficient=F(factorial(n),2**(n-1))*(h/2)**n
    require(num==[F(0)]*n+[coefficient], 'all-frequency numerator identity')
    for k in range(n):
        require(sum((-1)**j*weights[j]*nodes[j]**k for j in range(n+1))==0,'moment cancellation')
    L=F(n)/(max(ts)**2);U=F(n)/(min(ts)**2)
    phi=lambda v:sum(t*t*v/(1+t*t*v) for t in ts)
    for _ in range(bisections):
        mid=(L+U)/2
        if phi(mid)<n:L=mid
        else:U=mid
    require(phi(L)<=n<=phi(U),'maximum bracket')
    Q=[F(1)]
    for t in ts:Q=mul(Q,[F(1),t*t])
    mid=(L+U)/2;A=coefficient**2
    lo=A*mid**n/evaluate(Q,mid)
    hi=A*U**n/evaluate(Q,L)
    require(0<lo<=hi and hi/lo<1+F(1,2**85),'squared norm precision')
    return dict(v_lo=L,v_hi=U,norm_sq_lo=lo,norm_sq_hi=hi)

def target_gap(r,h,c,s):
    nodes,weights=pair_spec(r,h,c);h=rational(h);c=rational(c)
    require(type(s) is int and 1<=s<=6,'integer target order')
    # f(x)=(-1)^r h^(s+2)/(s! [h^2+(x-c)^2]).
    f=lambda x:(-1)**r*h**(s+2)/(factorial(s)*(h*h+(x-c)**2))
    require(h<=1,'target normalization uses h <= 1')
    for k in range(s+1):require(F(factorial(k),factorial(s))*h**(s-k)<=1,'C^s derivative bound')
    gap=sum((-1)**j*weights[j]*f(x) for j,x in enumerate(nodes))
    require(gap>0,'positive functional separation')
    return gap

def cost_certificate(inp):
    require(set(inp)=={'format','r','h','c','target_order','sigma','alpha','requested_width','noise_model','stopping_requirement'},'input schema')
    require(inp['format']=='FS2_stopped_cost_input_v1','input format')
    require(inp['noise_model']=='independent Gaussian channels; variance sigma^2/effort; cost effort','noise/precision law')
    require(inp['stopping_requirement']=='finite on every finite positive spectrum; terminal width bounded; uniform unconditional coverage','stopping semantics')
    r=inp['r'];h=rational(inp['h']);c=rational(inp['c']);s=inp['target_order']
    sig=rational(inp['sigma']);alpha=rational(inp['alpha']);w=rational(inp['requested_width'])
    require(sig>0 and w>0,'positive sigma and width')
    peak=full_norm(r,h,c);gap=target_gap(r,h,c,s)
    require(w<gap,'terminal width must be strictly less than the target gap')
    kl_lo,kl_hi=binary_information(alpha)
    bound_lo=2*sig*sig*kl_lo/peak['norm_sq_hi']
    bound_hi=2*sig*sig*kl_hi/peak['norm_sq_lo']
    return dict(format='FS2_stopped_cost_certificate_v1',input_hash=digest(inp),
                target_gap=gap,**peak,kl_lo=kl_lo,kl_hi=kl_hi,
                expected_effort_lower=bound_lo,exact_threshold_upper=bound_hi)

def verify_cost(inp,cert):
    ans=cost_certificate(inp)
    require(canonical(ans)==cert,'stopped-cost certificate mismatch')
    return ans

def ceil_fraction(x):
    x=rational(x);return -(-x.numerator//x.denominator)

if __name__=='__main__':
    import argparse
    ap=argparse.ArgumentParser();ap.add_argument('input');ap.add_argument('certificate');ap.add_argument('--output')
    a=ap.parse_args();inp=json.loads(Path(a.input).read_text());cert=json.loads(Path(a.certificate).read_text())
    ans=verify_cost(inp,cert)
    if a.output:save(a.output,ans)
    else:print('VERIFIED: fixed-confidence stopped experimental cost')
