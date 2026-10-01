"""Certificate verifier for exact positive RC ambiguity, standard library only.
Inputs define the task. A certificate cannot change the observation set or calibration mode.
"""
from fractions import Fraction as F
from math import factorial
from pathlib import Path
import json,sys
from exact import frac,need,add,mul,power,scale,val,product,digest,dump,trim

def product_poly(seq):
    out=[F(1)]
    for q in seq:out=mul(out,q)
    return out

def atoms(raw):
    need(isinstance(raw,list) and len(raw)>0,'nonempty atomic spectrum')
    out=[]
    for row in raw:
        need(isinstance(row,list) and len(row)==2,'atom shape')
        x,a=map(frac,row)
        need(-1<x<1 and a>0,'positive interior atom')
        out.append((x,a))
    need(len(set(x for x,a in out))==len(out),'distinct atom sites')
    return out

def channels(mu,theta,ws):
    ans=[]
    for w in ws:
        re,ni=theta[0],-w*theta[1]
        for x,a in mu:
            t=F(1)+x/2; z=w*t; dd=F(1)+z*z
            re+=a/dd;ni+=a*z/dd
        ans.extend([re,ni])
    return ans

def as_t(pol):
    out=[F(0)]
    for j,c in enumerate(pol):out=add(out,scale(power([F(-2),F(2)],j),c))
    return out

def verify(data,cert):
    need(data.get('schema')=='calibration_threshold_v1','schema')
    need(cert.get('input_sha256')==digest(data),'input binding')
    K=data['K'];need(isinstance(K,int) and not isinstance(K,bool) and K>=1,'K')
    mode=data['mode'];need(mode in ('known','unknown'),'calibration mode')
    known=mode=='known';ws=list(map(frac,data['frequencies']));wn=frac(data['additional_frequency'])
    need(len(ws)==(K if known else K+1),'claimed below-threshold count')
    need(min(ws+[wn])>0 and len(set(ws+[wn]))==len(ws)+1,'positive distinct frequencies')
    m0=atoms(data['reference']);m1=atoms(data['alternative'])
    need(len(m0)==K and len(m1)==K+1,'atom counts')
    need(sum(a for x,a in m0)==1,'reference is a probability measure')
    t0=list(map(frac,data['reference_calibration']));t1=list(map(frac,data['alternative_calibration']))
    need(len(t0)==len(t1)==2 and min(t0+t1)>0,'positive calibration')
    if known:need(t0==t1,'known calibration cannot change')
    need(frac(data['error_radius'])==0,'noiseless task')
    yy=list(map(frac,data['raw_data']));need(len(yy)==2*len(ws),'observation count')
    need(channels(m0,t0,ws)==yy==channels(m1,t1,ws),'original predictions exactly equal')
    W=list(map(frac,cert['annihilator']));f=list(map(frac,data['target_coefficients']));S=frac(cert['normalizer'])
    expected=product_poly([power([-x,F(1)],2) for x,a in m0])
    need(trim(W)==expected,'annihilator product')
    need(S>0 and trim(scale(f,S))==expected,'target normalization identity')
    norms=[sum(abs(f[j])*F(factorial(j),factorial(j-q)) for j in range(q,len(f))) for q in range(3)]
    need(max(norms)<=1,'physical C2 target norm')
    v0=sum(a*val(f,x) for x,a in m0);v1=sum(a*val(f,x) for x,a in m1)
    need(v0==0 and v1>0,'strict positive functional separation')
    need(frac(cert['target_gap'])==v1-v0,'exact target gap')
    # Independent check of weighted moment cancellation, not needed merely to check raw equality.
    D1=product(1+w*w for w in ws)
    def rho(x):
        t=1+x/2;d=product(1+w*w*t*t for w in ws)
        return d/D1 if known else d/(D1*t*t)
    for j in range(2*K):
        a0=sum(a*x**j/rho(x) for x,a in m0);a1=sum(a*x**j/rho(x) for x,a in m1)
        need(a0==a1,'low weighted moments agree')
    ext=ws+[wn]
    z0=channels(m0,t0,ext);z1=channels(m1,t1,ext)
    extra=[z1[-2]-z0[-2],z1[-1]-z0[-1]]
    need(extra==list(map(frac,cert['additional_row_difference'])) and any(extra),'extra row separates')
    # Whole-interval exposing dual identity. No sampled positivity decision.
    lam=list(map(frac,cert['extended_design_annihilator_dual']));need(len(lam)==2*len(ext),'dual length')
    lhs=[F(0)]
    for i,w in enumerate(ext):
        other=product_poly([[F(1),F(0),v*v] for j,v in enumerate(ext) if i!=j])
        lhs=add(lhs,scale(other,lam[2*i]));lhs=add(lhs,scale(mul([F(0),F(1)],other),w*lam[2*i+1]))
    rhs=scale(as_t(W),product(1+w*w for w in ext))
    if not known:
        rhs=mul([F(0),F(0),F(1)],rhs)
        need(sum(lam[::2])==0 and sum(w*lam[2*i+1] for i,w in enumerate(ext))==0,'dual nuisance annihilation')
    need(trim(lhs)==trim(rhs),'continuous reproducing identity')
    if known:
        zz0=channels(m0,[F(0),F(0)],ext);zz1=channels(m1,[F(0),F(0)],ext)
    else:zz0=z0;zz1=z1
    expose0=sum(a*b for a,b in zip(lam,zz0));expose1=sum(a*b for a,b in zip(lam,zz1))
    need(expose0==0 and expose1>0,'additional-design exposed face')
    return {'status':'PASS','K':K,'calibration_mode':mode,'frequencies':ws,
            'raw_difference':F(0),'target_gap':v1-v0,'reference_value':v0,
            'alternative_value':v1,'reference_mass':F(1),'alternative_mass':sum(a for x,a in m1),
            'additional_row_difference':extra,'extended_dual_excess':expose1,
            'reference_calibration':t0,'alternative_calibration':t1}

if __name__=='__main__':
    if len(sys.argv)!=3:raise SystemExit('Usage: python verify.py INPUT CERTIFICATE')
    data=json.loads(Path(sys.argv[1]).read_text());cert=json.loads(Path(sys.argv[2]).read_text())
    from exact import pack
    print(json.dumps(pack(verify(data,cert)),sort_keys=True))
