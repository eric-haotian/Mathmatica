"""Independent verifier of all-frequency positive-pair and adaptive-KL certificates.

The generator is not imported. Identities are checked as rational polynomials,
not by sampling a frequency grid. Numerical certificates verify special finite
instances; they do not prove the general statistical theorem or priority.
"""
from pathlib import Path
from fractions import Fraction as F
from math import comb, factorial
import argparse,json
from exact import frac,need,trim,add,scale,mul,power,val,product,dump,digest

def read_atoms(a):
    out=[(frac(x),frac(w)) for x,w in a]
    need(all(-1<=x<=1 and w>0 for x,w in out),'positive atoms on J')
    need(len(set(x for x,w in out))==len(out),'duplicate atoms')
    return out

def raw(atoms,theta,frequencies):
    ans=[]
    for om in map(frac,frequencies):
        need(om>0,'positive frequency')
        re=theta[0];nim=-om*theta[1]
        for x,w in atoms:
            t=1+x/2;den=1+om*om*t*t
            re+=w/den;nim+=w*om*t/den
        ans.extend([re,nim])
    return ans

def full_spectrum_pair(inp,cert):
    need(set(inp)=={'format','r','h','c','reference','alternative','theta','frequencies','y','epsilon','targets','statistics'},'input fields')
    need(inp['format']=='FS1_positive_pair_v1','format')
    need(cert['input_hash']==digest(inp),'input binding')
    r=inp['r'];need(type(r) is int and 1<=r<=20,'r')
    n=2*r;h=frac(inp['h']);c=frac(inp['c']);eps=frac(inp['epsilon'])
    need(h>0 and eps>0,'positive scales')
    ref=read_atoms(inp['reference']);alt=read_atoms(inp['alternative'])
    need(len(ref)==r and len(alt)==r+1,'atomic counts')
    need(sum(w for x,w in ref)==sum(w for x,w in alt)==1,'unit physical mass')
    theta=list(map(frac,inp['theta']));need(len(theta)==2 and min(theta)>=0,'nuisance positivity')
    nodes=[c+(l-r)*h for l in range(n+1)];ts=[1+x/2 for x in nodes]
    expected_ref=[(nodes[l],F(comb(n,l),2**(n-1))) for l in range(n+1) if l%2]
    expected_alt=[(nodes[l],F(comb(n,l),2**(n-1))) for l in range(n+1) if not l%2]
    need(ref==expected_ref and alt==expected_alt,'binomial pair structure')
    for k in range(n):
        need(sum(w*x**k for x,w in alt)==sum(w*x**k for x,w in ref),'matched low moments')
    # Independent reconstruction of the rational-function numerator in z=i*omega.
    signed=[F((-1)**l*comb(n,l),2**(n-1)) for l in range(n+1)]
    numerator=[F(0)]
    for l,w in enumerate(signed):
        p=[F(1)]
        for j,t in enumerate(ts):
            if j!=l:p=mul(p,[F(1),t])
        numerator=add(numerator,scale(p,w))
    expected=[F(0)]*n+[F(factorial(n),2**(n-1))*(h/2)**n]
    need(trim(numerator)==trim(expected),'full frequency numerator identity')
    Q=[F(1)]
    for t in ts:Q=mul(Q,[F(1),t*t])
    C=expected[-1]**2
    L=frac(cert['v_lo']);U=frac(cert['v_hi']);need(0<L<U,'frequency-square bracket')
    A=lambda v:sum(t*t*v/(1+t*t*v) for t in ts)
    need(A(L)<=n<=A(U),'unique stationary point bracket')
    need(U-L<=max(F(1),L)/2**70,'bracket precision')
    vlo=(L+U)/2
    sqlo=C*vlo**n/val(Q,vlo)
    squp=C*U**n/val(Q,L)
    need(frac(cert['norm_square_lo'])==sqlo and frac(cert['norm_square_hi'])==squp,'norm bound recomputation')
    need(0<sqlo<=squp<=eps*eps,'all-frequency noise radius')
    need(squp/sqlo<=1+F(1,2**60),'norm relative precision')
    flo=frac(cert['omega_lo']);fhi=frac(cert['omega_hi'])
    need(0<flo<fhi and flo*flo<=L and fhi*fhi>=U,'frequency bracket')
    nlo=frac(cert['norm_lo']);nhi=frac(cert['norm_hi'])
    need(0<nlo<=nhi and nlo*nlo<=sqlo and nhi*nhi>=squp,'norm interval')
    need(nhi==eps,'epsilon is the certified norm upper bound')
    y=raw(ref,theta,inp['frequencies']);need(y==list(map(frac,inp['y'])),'raw input values')
    ya=raw(alt,theta,inp['frequencies'])
    need(all(abs(a-b)<=eps for a,b in zip(y,ya)),'unchanged raw observation box')
    gaps={}
    for target in inp['targets']:
        need(set(target)=={'s','normalizer','sign'},'target fields')
        s=target['s'];need(type(s) is int and s in [2,3,4],'tested target orders')
        S=frac(target['normalizer']);sg=target['sign']
        need(S==factorial(s) and sg==(-1)**r,'uniform target normalization/sign')
        # Partial-fraction derivative bound: |f^(j)| <= j!*h^(s-j)/S.
        need(h<=1 and all(F(factorial(j),1)*h**(s-j)<=S for j in range(s+1)),'physical C^s normalization')
        f=lambda x:sg*h**(s+2)/(S*(h*h+(x-c)**2))
        gap=sum(w*f(x) for x,w in alt)-sum(w*f(x) for x,w in ref)
        need(gap>0,'positive target separation')
        need(gap==frac(cert['target_gaps'][str(s)]),'target gap recomputation')
        gaps[str(s)]=gap
    st=inp['statistics'];need(set(st)=={'sigma','alpha','n'},'statistical model fields')
    sig=frac(st['sigma']);alpha=frac(st['alpha']);N=st['n']
    need(sig>0 and 0<alpha<F(1,4) and type(N) is int and N>0,'statistical parameters')
    kl=F(N)*eps*eps/(2*sig*sig)
    need(kl==frac(cert['kl_upper']),'adaptive KL bound')
    need(cert['kl_at_most_eighth']==(kl<=F(1,8)),'KL flag')
    lower={}
    if kl<=F(1,8):
        for s,g in gaps.items():lower[s]=(1-2*alpha-F(1,4))*g
    need({k:frac(v) for k,v in cert['honest_expected_length_lower'].items()}==lower,'confidence length consequence')
    return dict(r=r,h=h,c=c,epsilon=eps,omega_lo=flo,omega_hi=fhi,norm_lo=nlo,norm_hi=nhi,
                target_gaps=gaps,kl_upper=kl,honest_expected_length_lower=lower,
                raw_peak_ratio=max(abs(a-b) for a,b in zip(y,ya))/eps,
                normalized_norm_lo=nlo/h**n,normalized_norm_hi=nhi/h**n,
                identity_degree=n,relative_norm_bracket=squp/sqlo-1)

def model_reproduction(obj):
    need(obj['format']=='FS1_reproduction_v1','model format')
    freq=list(map(frac,obj['frequencies']));m=len(freq)
    need(len(set(freq))==m and min(freq)>0,'distinct frequencies')
    factors=[[F(1),F(0),w*w] for w in freq];D=[F(1)]
    for fac in factors:D=mul(D,fac)
    D1=val(D,F(1));nums=[]
    for i,w in enumerate(freq):
        p=[F(1)]
        for j,fac in enumerate(factors):
            if i!=j:p=mul(p,fac)
        nums.extend([scale(p,1/D1),scale(mul(p,[F(0),w]),1/D1)])
    for degree,row0 in enumerate(obj['lambda']):
        row=list(map(frac,row0));need(len(row)==2*m,'map dimensions')
        p=[F(0)]
        for a,q in zip(row,nums):p=add(p,scale(q,a))
        need(p==mul([F(0),F(0),F(1)],power([F(-2),F(2)],degree)),'reproduction identity')
        need(sum(row[::2])==0 and sum(row[2*i+1]*freq[i] for i in range(m))==0,'nuisance annihilation')
    need(len(obj['lambda'])==2*(m-2)+1,'degree 2K coverage')
    return dict(frequencies=m,max_degree=len(obj['lambda'])-1,polynomial_identities=len(obj['lambda']))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('input');ap.add_argument('certificate',nargs='?');ap.add_argument('--output');ap.add_argument('--model',action='store_true');a=ap.parse_args()
    inp=json.loads(Path(a.input).read_text())
    if a.model:ans=model_reproduction(inp)
    else:
        need(a.certificate is not None,'certificate path');ans=full_spectrum_pair(inp,json.loads(Path(a.certificate).read_text()))
    if a.output:dump(a.output,ans)
    else:print('VERIFIED')
