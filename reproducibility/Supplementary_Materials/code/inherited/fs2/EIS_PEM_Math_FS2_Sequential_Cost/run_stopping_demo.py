"""Seeded Gaussian sample-mean prototype with exact endpoint certificates.

Efforts, error spending, physical target and rounding are bound to each stage.
The terminal range test is exact. Proposal caps may yield UNRESOLVED; the demo
is not the unbounded fair-refinement implementation of the completeness proof.
Samples are simulated batch means, not millions of collected observations.
"""
from pathlib import Path
from fractions import Fraction as F
import sys, json, math, time, argparse
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'vendor'/'phase'));sys.path.insert(0,str(ROOT/'vendor'/'software'))
import core as C
import direct_rational as R
from generate import find_anchor
from exact_cost import log_interval, ceil_fraction, save, digest, canonical

def dense_sequence():
    level=1
    while True:
        for k in range(1,2**level,2):yield F(1)+F(k,2**level)
        level+=1

def physical_raw(atoms,theta,ws):
    y=[]
    for w in ws:
        re=theta[0];ni=-w*theta[1]
        for x,mass in atoms:
            t=1+x/2;den=1+w*w*t*t
            re+=mass/den;ni+=mass*w*t/den
        y.extend([re,ni])
    return y

def weighted_target(ws):
    # Same physical target f=rho_(1,2,3,4)^(-1)/(2-x) at EVERY stage.
    # rho_current/rho_baseline = product of newly added factors / their values at t=1.
    num=[F(1)]
    for w in ws[4:]:
        fac=[1+w*w,w*w,w*w/4]
        num=C.mul(num,[a/(1+w*w) for a in fac])
    return dict(kind='rational',numerator=num,denominator=[F(2),F(-1)])

def run(case,seed,maxstages=7):
    import numpy as np
    rng=np.random.default_rng(seed)
    sigma=F(1,1000);alpha=F(1,20);width=F(1,20);b0=F(1,4096)
    if case=='two_atoms':atoms=[(F(-1,8),F(1,2)),(F(1,8),F(1,2))]
    elif case=='three_atoms':atoms=[(F(-1,4),F(3,16)),(F(0),F(5,8)),(F(1,4),F(3,16))]
    else:raise ValueError('case')
    theta=[F(1,10),F(1,100)];ws=list(map(F,[1,2,3,4]));extras=dense_sequence()
    out={'format':'FS2_sequential_demo_v1','case':case,'seed':seed,'sigma':sigma,'alpha':alpha,
         'requested_width':width,'b0':b0,'terminal_status':'UNRESOLVED','stages':[],
         'simulation':'independent Gaussian sample means from numpy PCG64; no physical experiment',
         'cost_counts':'one complex unit-precision observation; two real channels',
         'class_constant_tuning':'b0 chosen for mechanism demonstration, not from a worst-class sharp constant'}
    save(ROOT/'provenance'/f'{case}_{seed}.json',{'physical_atoms':atoms,'theta':theta})
    start=time.perf_counter()
    for j in range(1,maxstages+1):
        if j>1:ws.append(next(extras))
        b=b0/2**(j-1);risk=alpha/2**j;m=len(ws)
        _,loghi=log_interval(F(4*m)/risk)
        n=ceil_fraction(2*sigma*sigma*loghi/(b*b))
        means=physical_raw(atoms,theta,ws)
        draws=rng.standard_normal(2*m)
        # Round observed values to the b/16 grid. The mathematical error cap is b.
        step=b/16
        yf=[float(z)+float(sigma)/math.sqrt(n)*q for z,q in zip(means,draws)]
        yr=[F(round(z/float(step)))*step for z in yf]
        inp={'format':'R15_raw_data_v1','model':{'kind':'rc','parameters':ws},
             'y':yr,'delta':[3*b]*(2*m),'target':weighted_target(ws),'tau':width/4}
        inp=canonical(inp);key=f'{case}_{seed}_stage{j}'
        save(ROOT/'sequential'/'inputs'/f'{key}.json',inp)
        stage={'stage':j,'frequencies':list(ws),'b':b,'risk_spend':risk,'repetitions_per_frequency':n,
               'new_complex_observations':m*n,'input':key+'.json','input_hash':digest(inp),
               'rounded_mean_step':step,'reported_mean_binary64':[z.hex() for z in yf],
               'status':'UNRESOLVED','endpoints':[]}
        mo,y,de,tau,M=C.data_context(inp)
        try:
            anchor=find_anchor(mo,y,de,maxlevel=7)
            for sign in [1,-1]:
                cert,res,tr,ti=R.endpoint(inp,'direct',sign,anchor,maxrounds=18)
                if cert is None:raise RuntimeError('endpoint proposal cap')
                fn=f'{key}_s{sign}.json';save(ROOT/'sequential'/'certificates'/fn,cert)
                result=R.check(inp,json.loads((ROOT/'sequential'/'certificates'/fn).read_text()))
                stage['endpoints'].append({'certificate':fn,'result':result,'trace':tr,'timing':ti})
            plus,minus=[e['result'] for e in stage['endpoints']]
            interval=[-minus['upper'],plus['upper']];length=interval[1]-interval[0]
            stage.update(status='CERTIFIED_STOP' if length<=width else 'CERTIFIED_CONTINUE',
                         interval=interval,certified_length=length,
                         numerical_excess_budget=plus['gap']+minus['gap'])
        except (ValueError,RuntimeError) as ex:
            stage['reason']=str(ex)
        out['stages'].append(stage);out['elapsed_seconds']=time.perf_counter()-start
        if stage['status']=='CERTIFIED_STOP':out['terminal_status']='CERTIFIED_STOP'
        save(ROOT/'results'/f'sequential_{case}_{seed}.json',out)
        print(key,stage['status'],'n',n,'width',float(stage.get('certified_length',-1)),flush=True)
        if out['terminal_status']=='CERTIFIED_STOP':break
    out['total_complex_observations']=sum(s['new_complex_observations'] for s in out['stages'])
    out['elapsed_seconds']=time.perf_counter()-start
    save(ROOT/'results'/f'sequential_{case}_{seed}.json',out)
    return out

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('case',choices=['two_atoms','three_atoms']);ap.add_argument('--seed',type=int,default=20260916);ap.add_argument('--maxstages',type=int,default=7)
    a=ap.parse_args();run(a.case,a.seed,a.maxstages)
