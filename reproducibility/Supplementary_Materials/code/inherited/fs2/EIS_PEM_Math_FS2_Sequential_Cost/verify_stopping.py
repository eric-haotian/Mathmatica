"""Standard-library verification of staged data, stopping, and endpoint bounds.

Does not load the simulator, random generator, synthetic truth, or an optimizer.
It checks certificate arithmetic and alpha spending, not a Monte Carlo coverage rate.
"""
from pathlib import Path
from fractions import Fraction as F
import json, sys
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT))
sys.path.insert(0,str(ROOT/'vendor'/'phase'));sys.path.insert(0,str(ROOT/'vendor'/'software'))
import core as C
import direct_rational as R
from exact_cost import require, rational, log_interval, digest, canonical, save

def frequency_plan(stage):
    ws=list(map(F,[1,2,3,4]));level=1
    while len(ws)<stage+3:
        for i in range(1,2**level,2):
            if len(ws)>=stage+3:break
            ws.append(F(1)+F(i,2**level))
        level+=1
    return ws

def verify_stage(trace,stage,inp,certs):
    j=stage['stage'];require(type(j) is int and j>=1,'positive stage index')
    b0=rational(trace['b0']);b=b0/2**(j-1);alpha=rational(trace['alpha']);sigma=rational(trace['sigma'])
    w=rational(trace['requested_width']);require(0<alpha<F(1,4) and sigma>0 and b>0 and w>0,'experiment parameters')
    ws=frequency_plan(j);m=len(ws);risk=alpha/2**j
    require(list(map(rational,stage['frequencies']))==ws,'fixed dense frequency schedule')
    require(rational(stage['b'])==b and rational(stage['risk_spend'])==risk,'stage radius/risk binding')
    require(stage['input_hash']==digest(inp),'stage-input binding')
    require(inp['model']=={'kind':'rc','parameters':list(map(str,ws))},'model frequencies')
    require(list(map(rational,inp['delta']))==[3*b]*(2*m),'unaltered rounded-data box')
    require(rational(inp['tau'])==w/4,'numerical tolerance allocation')
    step=b/16;require(rational(stage['rounded_mean_step'])==step,'rounding scale')
    yf=[F.from_float(float.fromhex(z)) for z in stage['reported_mean_binary64']]
    yr=list(map(rational,inp['y']));require(len(yf)==len(yr)==2*m,'data dimensions')
    require(all(abs(a-z)<=b for a,z in zip(yr,yf)),'rounding error cap')
    _,lh=log_interval(F(4*m)/risk)
    n=stage['repetitions_per_frequency'];require(type(n) is int and n>0,'positive integer repetitions')
    require(F(n)>=2*sigma*sigma*lh/(b*b),'confidence-spending effort bound')
    require(stage['new_complex_observations']==m*n,'stage acquisition cost')
    mo=C.from_spec(inp['model']);base=C.get_model('rc',tuple(map(F,[1,2,3,4])))
    p=list(map(rational,inp['target']['numerator']));q=list(map(rational,inp['target']['denominator']))
    # Binding to one fixed physical target; not merely one weighted-coordinate formula.
    lhs=C.mul(C.mul(C.mul(p,mo['rhod']),base['rhon']),[F(2),F(-1)])
    rhs=C.mul(C.mul(q,mo['rhon']),base['rhod'])
    require(lhs==rhs,'physical target unchanged across designs')
    require(min(C.bernstein(q,-1,1))>0,'strict target denominator')
    results=[R.check(inp,c) for c in certs]
    require([c['sign'] for c in certs]==[1,-1],'two endpoint signs')
    interval=[-results[1]['upper'],results[0]['upper']];length=interval[1]-interval[0]
    require(canonical(interval)==stage['interval'],'outer interval arithmetic')
    require(str(length)==stage['certified_length'],'interval length')
    extra=sum(v['gap'] for v in results)
    require(extra<=w/4 and str(extra)==stage['numerical_excess_budget'],'complete numerical allowance')
    expected='CERTIFIED_STOP' if length<=w else 'CERTIFIED_CONTINUE'
    require(stage['status']==expected,'terminal width decision')
    for res,old in zip(results,stage['endpoints']):require(canonical(res)==old['result'],'endpoint replay equality')
    return {'stage':j,'status':expected,'length':length,'numerical_excess':extra,'cost':m*n,'risk':risk,'nodes':[z['nodes'] for z in results]}

def verify_trace(path):
    trace=json.loads(Path(path).read_text());require(trace['format']=='FS2_sequential_demo_v1','trace format')
    out=[]
    for j,stage in enumerate(trace['stages'],1):
        require(stage['stage']==j,'uninterrupted stage indices')
        require(stage['status'] in ['CERTIFIED_STOP','CERTIFIED_CONTINUE'],'unresolved not certified')
        inp=json.loads((ROOT/'sequential'/'inputs'/stage['input']).read_text())
        certs=[json.loads((ROOT/'sequential'/'certificates'/e['certificate']).read_text()) for e in stage['endpoints']]
        out.append(verify_stage(trace,stage,inp,certs))
        if stage['status']=='CERTIFIED_STOP':require(j==len(trace['stages']),'no stages after stopping')
    require(out and out[-1]['status']=='CERTIFIED_STOP' and trace['terminal_status']=='CERTIFIED_STOP','successful stopping trace')
    require(sum(o['cost'] for o in out)==trace['total_complex_observations'],'total acquisition count')
    require(sum(o['risk'] for o in out)<=rational(trace['alpha']),'spent risk no more than global alpha')
    return {'case':trace['case'],'stages':out,'total_cost':sum(o['cost'] for o in out),
            'spent_alpha':sum(o['risk'] for o in out),'final_interval':trace['stages'][-1]['interval']}

if __name__=='__main__':
    import argparse
    ap=argparse.ArgumentParser();ap.add_argument('trace');ap.add_argument('--output');a=ap.parse_args()
    ans=verify_trace(a.trace)
    if a.output:save(a.output,ans)
    else:print('VERIFIED: sequential decisions and continuous endpoint certificates')
