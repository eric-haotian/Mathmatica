"""Build new stopping-cost certificates from stated binomial positive pairs."""
from pathlib import Path
from fractions import Fraction as F
import sys, json, hashlib
from exact_cost import *
ROOT=Path(__file__).resolve().parent

def main():
    rows=[]
    # Six representative pairs; a case is not a stochastic experiment.
    cases=[(1,F(1,8)),(2,F(1,8)),(3,F(1,8)),(3,F(1,16)),(4,F(1,8)),(4,F(1,16))]
    for r,h in cases:
        for alpha in [F(1,10),F(1,20),F(1,100)]:
            gap=target_gap(r,h,0,2)
            width=F(1,1000) if r==3 and h==F(1,8) else gap/2
            obj=dict(format='FS2_stopped_cost_input_v1',r=r,h=h,c=F(0),target_order=2,
                     sigma=F(1,1000),alpha=alpha,requested_width=width,
                     noise_model='independent Gaussian channels; variance sigma^2/effort; cost effort',
                     stopping_requirement='finite on every finite positive spectrum; terminal width bounded; uniform unconditional coverage')
            key=f'r{r}_h{h.denominator}_alpha{alpha.denominator}'
            cert=cost_certificate(obj)
            save(ROOT/'inputs'/f'{key}.json',obj);save(ROOT/'certificates'/f'{key}.json',cert)
            rows.append(dict(case=key,r=r,h=h,alpha=alpha,requested_width=width,target_gap=gap,
                             expected_effort_lower=cert['expected_effort_lower'],
                             threshold_upper=cert['exact_threshold_upper'],
                             displayed_floor=int(cert['expected_effort_lower'])))
    save(ROOT/'results'/'cost_summary.json',rows)
    for row in rows:
        if row['r']==3 and row['h']==F(1,8):
            print(row['case'],'gap',float(row['target_gap']),'cost lower',float(row['expected_effort_lower']))

if __name__=='__main__':main()
