"""Named exact checks and mutation tests. These do not prove coverage by simulation."""
from pathlib import Path
from fractions import Fraction as F
from copy import deepcopy
import json,hashlib
from exact_cost import *
from verify_stopping import verify_stage,frequency_plan
ROOT=Path(__file__).resolve().parent
checks=[]
def check(name,ok):
    require(ok,name);checks.append({'name':name,'status':'PASS'})
def reject(name,fun):
    try:fun()
    except (ValueError,ZeroDivisionError,KeyError,IndexError):checks.append({'name':name,'status':'REJECTED_AS_REQUIRED'});return
    raise AssertionError('accepted mutation: '+name)

l,u=log_interval(F(19));check('log19 bracket width',u-l<F(1,10**50))
check('log19 elementary decimal enclosure',F(2944,1000)<l<u<F(2945,1000))
l2,u2=log_interval(F(2));l4,u4=log_interval(F(4));check('log4 equals twice log2 interval overlap',max(l4,2*l2)<=min(u4,2*u2))
ln,un=log_interval(F(1,19));check('log reciprocal negates interval',ln==-u and un==-l)
lo=[]
for a in [F(1,10),F(1,20),F(1,100)]:lo.append(binary_information(a)[0])
check('confidence information increases at 90 95 99 percent',lo[0]<lo[1]<lo[2])
for alpha in [F(1,10),F(1,20),F(1,100)]:
    q=alpha/2;x=4*q
    partial=sum(F(j*j)*x**(j-1) for j in range(1,25))
    closed=(1+x)/(1-x)**3
    check(f'cost geometric-square upper identity alpha={alpha}',partial<closed<=12)
    risk=sum(alpha/F(2**j) for j in range(1,25))
    check(f'exact risk tail alpha={alpha}',risk+alpha/F(2**24)==alpha)
    check(f'uniform geometric cost bound alpha={alpha}',1/(1-q)**2<=F(64,49))

inp=json.loads((ROOT/'inputs'/'r3_h8_alpha20.json').read_text());cert=json.loads((ROOT/'certificates'/'r3_h8_alpha20.json').read_text())
ans=verify_cost(inp,cert)
check('main physical gap exactly 9/5120',ans['target_gap']==F(9,5120))
check('main expected cost lower exceeds 48614350',ans['expected_effort_lower']>48614350)
changed=deepcopy(inp);changed['sigma']='1/500';rr=cost_certificate(changed)
check('noise doubling costs factor four',rr['expected_effort_lower']==4*ans['expected_effort_lower'])
# The pair bound is constant for any requested width strictly below its fixed gap.
changed=deepcopy(inp);changed['requested_width']='1/2000'
check('fixed pair width threshold',cost_certificate(changed)['expected_effort_lower']==ans['expected_effort_lower'])
for field,value,label in [('requested_width','9/5120','width equals separation'),('sigma','-1/1000','negative noise'),('alpha','1/2','unsupported confidence'),('noise_model','relative noise','unproved noise law'),('stopping_requirement','coverage only on sparse truth','narrowed coverage class')]:
    bad=deepcopy(inp);bad[field]=value
    reject(label,lambda bad=bad:cost_certificate(bad))
bad=deepcopy(cert);bad['expected_effort_lower']=str(2*rational(bad['expected_effort_lower']))
reject('inflated experimental lower bound',lambda:verify_cost(inp,bad))
bad=deepcopy(cert);bad['norm_sq_hi']='1/1000000000000000'
reject('false whole-frequency bound',lambda:verify_cost(inp,bad))
bad=deepcopy(inp);bad['h']='1/2'
reject('nodes outside compact physical domain',lambda:cost_certificate(bad))

trace=json.loads((ROOT/'results'/'sequential_two_atoms_20260916.json').read_text());stage=trace['stages'][0]
sinp=json.loads((ROOT/'sequential'/'inputs'/stage['input']).read_text());cs=[json.loads((ROOT/'sequential'/'certificates'/z['certificate']).read_text()) for z in stage['endpoints']]
check('first valid stage must continue',verify_stage(trace,stage,sinp,cs)['status']=='CERTIFIED_CONTINUE')
for field,value,label in [('status','CERTIFIED_STOP','premature stopping'),('repetitions_per_frequency',1,'unpaid statistical precision'),('risk_spend','1/10','incorrect risk spending'),('b','1/16','incorrect stage radius')]:
    bad=deepcopy(stage);bad[field]=value
    reject(label,lambda bad=bad:verify_stage(trace,bad,sinp,cs))
bad_inp=deepcopy(sinp);bad_inp['target']['numerator']=['2'];bad_stage=deepcopy(stage);bad_stage['input_hash']=digest(bad_inp)
reject('changed physical target after data',lambda:verify_stage(trace,bad_stage,bad_inp,cs))
bad_inp=deepcopy(sinp);bad_inp['delta']=[str(2*rational(z)) for z in sinp['delta']];bad_stage=deepcopy(stage);bad_stage['input_hash']=digest(bad_inp)
reject('unreported enlargement of raw boxes',lambda:verify_stage(trace,bad_stage,bad_inp,cs))
bad_cs=deepcopy(cs);bad_cs[0]['gamma'][0]=str(rational(bad_cs[0]['gamma'][0])+1)
reject('corrupted continuous dual',lambda:verify_stage(trace,stage,sinp,bad_cs))

old=json.loads((ROOT/'checks'/'pre_repair_input_hashes.json').read_text())
new={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (ROOT/'sequential'/'inputs').glob('*.json')}
check('metadata alias repair leaves numerical inputs byte identical',old==new)
check('third-stage model appends dense locations',frequency_plan(3)==list(map(F,[1,2,3,4]))+[F(3,2),F(5,4)])

save(ROOT/'results'/'mechanism_checks.json',{'passed':True,'checks':checks,'count':len(checks),
                                        'rejection_checks':sum(x['status']=='REJECTED_AS_REQUIRED' for x in checks)})
print(len(checks),'named checks,',sum(x['status']=='REJECTED_AS_REQUIRED' for x in checks),'rejections')
