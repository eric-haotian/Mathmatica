"""Truth-free R15 proposal engine + certified endpoint stopping.
Fast LP/roots only propose. Every LP has an exact Bland-simplex fallback. --complete enables
unbounded fair refinement; this guarantees termination under the theorem assumptions,
not an efficient running time. Default caps return UNRESOLVED, never success.
"""
from core import *
from exact_lp import solve_lp
from time import perf_counter

def dense_nodes(level):return {F(-1)+F(2*j,2**level) for j in range(2**level+1)}
def columns(mo,nodes):
    return [[val(p,x) for p in mo['nums']] for x in nodes]+list(map(list,zip(*mo['B']))) if mo['k'] else [[val(p,x) for p in mo['nums']] for x in nodes]
def lp_rows(mo,nodes):return list(map(list,zip(*columns(mo,nodes))))

def find_anchor(mo,y,de,complete=False,maxlevel=10,force_exact=False):
    extra={F(0)}
    if len(mo['lam'])>=3:
        m=[dot(l,y) for l in mo['lam'][:3]]
        if m[0]>0:extra.add(max(F(-1),min(F(1),F(round((m[1]/m[0])*(1<<40)),1<<40))))
        if m[1]!=0:extra.add(max(F(-1),min(F(1),F(round((m[2]/m[1])*(1<<40)),1<<40))))
    level=2;history=[]
    while complete or level<=maxlevel:
        nodes=sorted(dense_nodes(level)|extra);A=lp_rows(mo,nodes);n=len(nodes)+mo['k']
        AA=[row+[F(1)] for row in A]+[[-v for v in row]+[F(1)] for row in A]+[[F(0)]*n+[F(1)]]
        bb=[b+d for b,d in zip(y,de)]+[d-b for b,d in zip(y,de)]+[min(de)]
        c=[F(0)]*n+[F(1)]
        rr,st=solve_lp(AA,bb,c,complete=complete,force_exact=force_exact)
        history.append(dict(level=level,nodes=len(nodes),lp=st))
        if rr is not None:
            x,u=rr
            if x[-1]>0:
                atoms=[[z,w*val(mo['den'],z)] for z,w in zip(nodes,x[:len(nodes)]) if w]
                theta=x[len(nodes):len(nodes)+mo['k']]
                zz=raw(mo,atoms,theta);margin=min(d-abs(a-b) for a,b,d in zip(zz,y,de))
                need(margin>0,'anchor verified strict slack')
                return dict(atoms=atoms,theta=theta,margin=margin,trace=history)
        level+=1
    raise RuntimeError('UNRESOLVED: no certified strict anchor within configured cap')

def candidates(mo,g,p,sg,bits=40):
    import numpy as np
    from numpy.polynomial import polynomial as NP
    pol=dual_numerator(mo,g,p,sg)
    derivative=sub(mul(deriv(pol),mo['den']),mul(pol,deriv(mo['den'])))
    ff=np.array([float(x) for x in derivative]);z=max(abs(ff),default=0)
    if not z:return []
    rr=NP.polyroots(ff/z);out=[]
    for r in rr:
        if abs(r.imag)<1e-6 and -1<r.real<1:
            x=F(round(float(r.real)*(1<<bits)),1<<bits)
            value=val(pol,x)/val(mo['den'],x)
            if value<0:out.append((value,x))
    return sorted(out)

def endpoint(inp,sign,anchor,complete=False,maxrounds=32,force_exact=False):
    mo,y,de,tau,M=data_context(inp);pr=proxy(mo,y,de,tau,M,inp['target']);p=pr['poly']
    shift=tau/(32*M);X=dense_nodes(3)|{x for x,w in anchor['atoms']};hist=[];it=0;previous=None
    while complete or it<maxrounds:
        # A fair uniform refinement subsequence makes the mathematical fallback complete.
        lev=3+it//6;X|=dense_nodes(lev);nodes=sorted(X);A=lp_rows(mo,nodes)
        AA=A+[[-v for v in row] for row in A]
        bb=[b+d for b,d in zip(y,de)]+[d-b for b,d in zip(y,de)]
        c=[sign*val(p,x)*val(mo['den'],x) for x in nodes]+[F(0)]*mo['k']
        warm=None if previous is None else [previous[0].get(x,F(0))/val(mo['den'],x) for x in nodes]+previous[1]
        rr,st=solve_lp(AA,bb,c,complete=complete,force_exact=force_exact,warm=warm)
        if rr is None:
            hist.append(dict(iteration=it,nodes=len(nodes),lp=st));it+=1
            continue
        sol,u=rr;g=[u[i]-u[i+mo['m']] for i in range(mo['m'])]
        atoms=[[x,w*val(mo['den'],x)] for x,w in zip(nodes,sol[:len(nodes)]) if w];theta=sol[len(nodes):];previous=({x:w for x,w in atoms},theta)
        corr=[a+shift*b for a,b in zip(g,mo['lam'][0])]
        maxdepth=36+it
        cert=dict(format='R15_endpoint_v1',input_hash=digest(inp),proxy_hash=digest(pr),sign=sign,
                  atoms=atoms,theta=theta,anchor_atoms=anchor['atoms'],anchor_theta=anchor['theta'],
                  nodes=nodes,gamma=g,shift=shift,maxdepth=maxdepth)
        status=''
        try:
            pp=dual_numerator(mo,corr,p,sign);proof=positive(pp,maxdepth=maxdepth)
            frozen=jsonable(cert);res=verify(inp,frozen)
            hist.append(dict(iteration=it,nodes=len(nodes),lp=st,accepted=True,positivity=proof))
            return frozen,res,hist
        except ValueError as ex:
            status=str(ex)
            if status.startswith('negative polynomial at '):X.add(F(status.split(' at ')[1]))
        try: cand=candidates(mo,g,p,sign,bits=40+it//4)
        except Exception: cand=[]  # uniform refinement remains the complete fallback
        for _,x in cand[:16]:X.add(x)
        hist.append(dict(iteration=it,nodes=len(nodes),lp=st,accepted=False,reason=status,added_roots=len(cand[:16])))
        it+=1
    raise RuntimeError('UNRESOLVED: endpoint did not meet certified tolerance within cap; last='+str(hist[-1] if hist else None))

def run(input_path,complete=False,maxrounds=32,force_exact=False):
    start=perf_counter();inp=load(input_path);mo,y,de,tau,M=data_context(inp)
    anchor=find_anchor(mo,y,de,complete=complete,force_exact=force_exact)
    result={'input':Path(input_path).name,'input_hash':digest(inp),'anchor':anchor,'endpoints':[]}
    key=Path(input_path).stem
    for sg in [1,-1]:
        cert,res,trace=endpoint(inp,sg,anchor,complete,maxrounds,force_exact)
        path=ROOT/'certificates'/f'{key}_s{sg}.json';save(path,cert)
        checked=verify(inp,load(path));need(checked==res,'readback result equality')
        result['endpoints'].append(dict(certificate=path.name,verified=res,trace=trace))
        print(key,sg,'degree',res['proxy_degree'],'vars',res['moment_scalars'],'nodes',res['finite_nodes'],
              'gap/tau',float(res['gap_over_tau']),flush=True)
    rp,rm=[x['verified'] for x in result['endpoints']]
    result['information_range_outer']=[-rm['upper'],rp['upper']]
    result['information_width_upper']=rp['upper']+rm['upper']
    result['information_width_lower']=max(F(0),rp['lower']+rm['lower'])
    result['numerical_excess_budget']=rp['gap']+rm['gap']
    need(result['numerical_excess_budget']<=tau,'full range excess')
    result['generation_seconds']=perf_counter()-start
    save(ROOT/'results'/f'{key}.json',result)
    return result
if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('input');p.add_argument('--complete',action='store_true');p.add_argument('--force-exact',action='store_true');p.add_argument('--max-rounds',type=int,default=32)
    a=p.parse_args();run(a.input,a.complete,a.max_rounds,a.force_exact)
