"""Named exact checks and adversarial certificate rejection tests."""
from pathlib import Path
from copy import deepcopy
from fractions import Fraction as F
from math import comb,factorial
import json
from exact import *
from verify import full_spectrum_pair,model_reproduction
ROOT=Path(__file__).resolve().parent
checks=[]
def ck(name,condition):
    need(condition,name);checks.append(dict(name=name,result='PASS'))

def rejected(name,change):
    ii=deepcopy(base_i);cc=deepcopy(base_c);change(ii,cc)
    # Re-bind mutated inputs to prevent relying on hash mismatch alone.
    cc['input_hash']=digest(ii)
    try:full_spectrum_pair(ii,cc)
    except (ValueError,TypeError,KeyError,ZeroDivisionError,OverflowError):checks.append(dict(name=name,result='REJECTED'));return
    raise AssertionError('tamper accepted: '+name)

base_i=json.loads((ROOT/'inputs/r3_h3.json').read_text());base_c=json.loads((ROOT/'certificates/r3_h3.json').read_text())
base=full_spectrum_pair(base_i,base_c)
ck('two_million_adaptive_queries_KL_at_most_1_over_8',frac(base['kl_upper'])<=F(1,8))
ck('physical_C2_gap_9_over_5120',base['target_gaps']['2']==F(9,5120))
ck('95_percent_expected_length_lower_117_over_102400',base['honest_expected_length_lower']['2']==F(117,102400))
for r in range(1,5):
    n=2*r
    c=F(0);h=F(1,32)
    v=[F((-1)**l*comb(n,l),2**(n-1)) for l in range(n+1)]
    for k in [0,n-1]:ck(f'r{r}_moment_{k}',sum(v[l]*(l-r)**k for l in range(n+1))==0)
    ck(f'r{r}_first_nonzero_moment',sum(v[l]*(l-r)**n for l in range(n+1))==F(factorial(n),2**(n-1)))
    # Unbounded-band Fourier control: at omega*h=pi the magnitudes sum to 2.
    ck(f'r{r}_Fourier_countermodel_amplitude_two',abs(sum(v[l]*(-1)**l for l in range(n+1)))==2)
    ck(f'r{r}_resolvent_target_sign',(-1)**r*sum(v[l]/(1+(l-r)**2) for l in range(n+1))>0)
    # Geometric scaling of the exact rational target separation.
    aa=json.loads((ROOT/'certificates'/f'r{r}_h4.json').read_text())
    bb=json.loads((ROOT/'certificates'/f'r{r}_h5.json').read_text())
    for s in [2,3,4]:ck(f'r{r}_target_s{s}_exact_halving',frac(bb['target_gaps'][str(s)])/frac(aa['target_gaps'][str(s)])==F(1,2**s))
rejected('negative_mass',lambda i,c:i['reference'][0].__setitem__(1,'-1'))
rejected('out_of_domain_atom',lambda i,c:i['alternative'][0].__setitem__(0,'-2'))
rejected('nuisance_negative',lambda i,c:i['theta'].__setitem__(0,'-1'))
rejected('radius_halved',lambda i,c:i.__setitem__('epsilon',str(frac(i['epsilon'])/2)))
rejected('reference_weight_changed',lambda i,c:i['reference'][0].__setitem__(1,'1/9'))
rejected('frequency_square_reversed',lambda i,c:c.__setitem__('v_hi',c['v_lo']))
rejected('frequency_bracket_elsewhere',lambda i,c:c.__setitem__('v_lo','100'))
rejected('norm_upper_zero',lambda i,c:c.__setitem__('norm_square_hi','0'))
rejected('norm_interval_falsified',lambda i,c:c.__setitem__('norm_hi','1'))
rejected('wrong_raw_y',lambda i,c:i['y'].__setitem__(0,'0'))
rejected('target_norm_weakened',lambda i,c:i['targets'][0].__setitem__('normalizer','1'))
rejected('target_sign_reversed',lambda i,c:i['targets'][0].__setitem__('sign',1))
rejected('gap_inflated',lambda i,c:c['target_gaps'].__setitem__('2','1'))
rejected('KL_artificially_zero',lambda i,c:c.__setitem__('kl_upper','0'))
rejected('confidence_bound_inflated',lambda i,c:c['honest_expected_length_lower'].__setitem__('2','1'))
rejected('unlisted_noise_model',lambda i,c:i['statistics'].__setitem__('heteroscedastic',True))
rejected('floating_point_radius',lambda i,c:i.__setitem__('epsilon',1e-8))
rejected('sigma_zero',lambda i,c:i['statistics'].__setitem__('sigma','0'))
ck('at_least_50_checks',len(checks)>=50)
dump(ROOT/'results/mechanism_checks.json',dict(status='PASS',checks=checks,count=len(checks),rejections=sum(a['result']=='REJECTED' for a in checks)))
print('PASS',len(checks),'checks;',sum(a['result']=='REJECTED' for a in checks),'rejections')
