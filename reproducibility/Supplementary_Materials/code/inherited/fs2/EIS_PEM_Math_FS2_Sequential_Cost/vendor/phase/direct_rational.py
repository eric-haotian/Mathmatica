"""v0.4 general-rational specialization of the inherited endpoint engine.

L: inherited root-free localized polynomial; G: global geometric polynomial;
R: direct rational target, hence zero target approximation error.
All arms use the same exact LP engine, anchor, nodes policy, dual mass shift,
and exact continuous nonnegativity checker. This is not a reproduction of
any external author's complete algorithm. The shared primitives are inherited
from the archived v0.1 software without modification.
"""
from pathlib import Path
import sys, time, json, argparse, platform
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'software'))
from core import *

def target_proxy(inp, method):
    mo,y,de,tau,M=data_context(inp)
    spec=inp['target']
    need(method=='direct' and set(spec)=={'kind','numerator','denominator'},'rational schema')
    need(spec['kind']=='rational','rational target kind')
    num=trim(spec['numerator']);den=trim(spec['denominator'])
    need(positive(den,maxdepth=48) is not None,'positive target denominator')
    return {'numerator':num,'denominator':den,'E':F(0),
            'degree':None,'moment_scalars':None,'localization':None}

def obj(pr,x):return val(pr['numerator'],x)/val(pr['denominator'],x)
def residual(mo,g,pr,sg):
    su=[F(0)]
    for t,p in zip(g,mo['nums']):su=add(su,scale(p,t))
    return sub(mul(su,pr['denominator']),scale(mul(mo['den'],pr['numerator']),sg))

def check(inp,cert):
    """Exact replay; no floating optimizer, polynomial roots or hidden source."""
    need(cert['format']=='V04_rational_endpoint_v1','certificate format')
    need(cert['input_hash']==digest(inp),'input binding')
    mo,y,de,tau,M=data_context(inp);pr=target_proxy(inp,cert['method'])
    need(cert['proxy_hash']==digest(pr),'proxy binding')
    sg=cert['sign'];need(sg in (-1,1),'sign')
    nodes=[frac(x) for x in cert['nodes']]
    need(nodes==sorted(set(nodes)) and nodes[0]==-1 and nodes[-1]==1,'node domain')
    atoms=[[frac(x),frac(w)] for x,w in cert['atoms']];theta=list(map(frac,cert['theta']))
    anc=[[frac(x),frac(w)] for x,w in cert['anchor_atoms']];ath=list(map(frac,cert['anchor_theta']))
    for aset,th,strict in [(atoms,theta,False),(anc,ath,True)]:
        z=raw(mo,aset,th)
        margins=[d-abs(a-b) for a,b,d in zip(z,y,de)]
        need(min(margins)>0 if strict else min(margins)>=0,'original observation box')
        need(all(x in nodes for x,w in aset),'support included')
    g=list(map(frac,cert['gamma']));need(len(g)==mo['m'],'dual dimension')
    need(all(sum(g[i]*mo['B'][i][j] for i in range(mo['m']))>=0 for j in range(mo['k'])),'nuisance cone')
    need(all(dot(g,obs(mo,x))>=sg*obj(pr,x) for x in nodes),'nodal dual')
    # All denominators are also checked on the full continuous interval.
    need(positive(mo['den'],maxdepth=48) is not None,'model denominator')
    need(positive(pr['denominator'],maxdepth=48) is not None,'target denominator')
    shift=frac(cert['shift']);need(shift==tau/(32*M),'common correction budget')
    gc=[a+shift*b for a,b in zip(g,mo['lam'][0])]
    pos=positive(residual(mo,gc,pr,sg),maxdepth=cert['maxdepth'])
    P=sg*sum(w*obj(pr,x) for x,w in atoms)
    D0=support(g,y,de);D=support(gc,y,de);E=pr['E']
    need(D0==P,'finite LP optimality equality')
    need(P<=D<=P+shift*M,'continuous correction cost')
    gap=D-P+2*E;need(0<=gap<=tau/2,'full endpoint tolerance')
    return {'lower':P-E,'upper':D+E,'gap':gap,'gap_over_tau':gap/tau,
            'finite_gap':D0-P,'continuous_correction':D-P,'approximation_pair':2*E,
            'degree':pr['degree'],'moment_scalars':pr['moment_scalars'],
            'nodes':len(nodes),'atoms':len(atoms),'positive':pos,'E':E,
            'max_primal_bits':max((max(v.numerator.bit_length(),v.denominator.bit_length()) for x,w in atoms for v in (x,w)),default=0),
            'max_dual_bits':max(max(v.numerator.bit_length(),v.denominator.bit_length()) for v in g)}

def proposals(mo,g,pr,sg,bits):
    import numpy as np
    from numpy.polynomial import polynomial as NP
    pol=residual(mo,g,pr,sg);den=mul(mo['den'],pr['denominator'])
    dp=sub(mul(deriv(pol),den),mul(pol,deriv(den)))
    ff=np.array([float(v) for v in dp]);mag=max(abs(ff),default=0)
    if not mag:return []
    out=[]
    for z in NP.polyroots(ff/mag):
        if abs(z.imag)<1e-6 and -1<z.real<1:
            x=F(round(float(z.real)*(1<<bits)),1<<bits)
            r=val(pol,x)/val(den,x)
            if r<0:out.append((r,x))
    return sorted(out)

def endpoint(inp,method,sg,anchor,maxrounds=28):
    from generate import dense_nodes, lp_rows
    from exact_lp import solve_lp
    mo,y,de,tau,M=data_context(inp);t0=time.perf_counter();pr=target_proxy(inp,method)
    formation=time.perf_counter()-t0
    X=dense_nodes(3)|{x for x,w in anchor['atoms']};previous=None;trace=[]
    for it in range(maxrounds):
        X|=dense_nodes(3+it//6);nodes=sorted(X);A=lp_rows(mo,nodes)
        AA=A+[[-v for v in row] for row in A]
        bb=[b+d for b,d in zip(y,de)]+[d-b for b,d in zip(y,de)]
        c=[sg*obj(pr,x)*val(mo['den'],x) for x in nodes]+[F(0)]*mo['k']
        warm=None if previous is None else [previous[0].get(x,F(0))/val(mo['den'],x) for x in nodes]+previous[1]
        rr,st=solve_lp(AA,bb,c,complete=False,warm=warm)
        rec={'iteration':it,'nodes':len(nodes),'lp':st};trace.append(rec)
        if rr is None:continue
        sol,u=rr;g=[u[i]-u[i+mo['m']] for i in range(mo['m'])]
        atoms=[[x,w*val(mo['den'],x)] for x,w in zip(nodes,sol[:len(nodes)]) if w]
        theta=sol[len(nodes):];previous=({x:w for x,w in atoms},theta)
        cert={'format':'V04_rational_endpoint_v1','method':method,'input_hash':digest(inp),
              'proxy_hash':digest(pr),'sign':sg,'nodes':nodes,'atoms':atoms,'theta':theta,
              'anchor_atoms':anchor['atoms'],'anchor_theta':anchor['theta'],'gamma':g,
              'shift':tau/(32*M),'maxdepth':36+it}
        try:
            frozen=jsonable(cert);res=check(inp,frozen);rec['accepted']=True
            return frozen,res,trace,{'formation_seconds':formation,'endpoint_seconds':time.perf_counter()-t0}
        except ValueError as ex:
            rec['accepted']=False;rec['reason']=str(ex)
            if str(ex).startswith('negative polynomial at '):X.add(F(str(ex).split(' at ')[1]))
        try:cand=proposals(mo,g,pr,sg,40+it//4)
        except Exception as ex:cand=[];rec['proposal_error']=repr(ex)
        X.update(x for v,x in cand[:16]);rec['proposed_roots']=len(cand[:16])
    return None,None,trace,{'formation_seconds':formation,'endpoint_seconds':time.perf_counter()-t0,'status':'UNRESOLVED'}


def main():
    import argparse
    from generate import find_anchor
    parser=argparse.ArgumentParser()
    parser.add_argument('key');parser.add_argument('--sign',type=int,choices=[-1,1])
    args=parser.parse_args()
    inp=load(HERE/'inputs'/(args.key+'.json'))
    mo,y,de,tau,M=data_context(inp)
    t0=time.perf_counter();anchor=find_anchor(mo,y,de)
    original_hash=digest(inp)
    scale_factor=F(863656)
    scaled=dict(inp);scaled['target']=dict(inp['target'])
    scaled['target']['numerator']=scale(trim(inp['target']['numerator']),scale_factor)
    scaled['tau']=tau*scale_factor
    out={'input_hash':original_hash,'target_and_tolerance_scale':scale_factor,
         'anchor':anchor,'endpoints':[]}
    for sg in ([args.sign] if args.sign else [1,-1]):
        cert,res,trace,timing=endpoint(scaled,'direct',sg,anchor)
        item={'sign':sg,'trace':trace,'timing':timing,'status':'UNRESOLVED'}
        if cert is not None:
            cert['input_hash']=original_hash
            cert['proxy_hash']=digest(target_proxy(inp,'direct'))
            cert['gamma']=[frac(x)/scale_factor for x in cert['gamma']]
            cert['shift']=frac(cert['shift'])/scale_factor
            frozen=jsonable(cert)
            result=check(inp,frozen)
            save(HERE/'certificates'/f'{args.key}_s{sg}.json',frozen)
            need(check(inp,load(HERE/'certificates'/f'{args.key}_s{sg}.json'))==result,'disk replay')
            item.update(status='CERTIFIED',result=result)
            print(args.key,sg,'nodes',result['nodes'],'gap/tau',float(result['gap_over_tau']),flush=True)
        out['endpoints'].append(item)
        out['elapsed_seconds']=time.perf_counter()-t0
        save(HERE/'results'/f'{args.key}{"_sign"+str(args.sign) if args.sign else ""}.json',out)

if __name__=='__main__':main()
