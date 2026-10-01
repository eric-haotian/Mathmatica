"""Summarize verified rational results; displays are directed, not rounded inward."""
from pathlib import Path
from fractions import Fraction as F
from decimal import Decimal,localcontext,ROUND_FLOOR,ROUND_CEILING
import json,csv
ROOT=Path(__file__).resolve().parent

def dec(x,places,up=False):
    x=F(x)
    with localcontext() as c:
        c.prec=150
        d=Decimal(x.numerator)/Decimal(x.denominator)
        return str(d.quantize(Decimal(10)**(-places),rounding=ROUND_CEILING if up else ROUND_FLOOR))

def sci(x,digits=7,up=False):
    x=F(x)
    with localcontext() as c:
        c.prec=150
        d=Decimal(x.numerator)/Decimal(x.denominator)
        e=d.adjusted()
        d=d.quantize(Decimal(10)**(e-digits+1),rounding=ROUND_CEILING if up else ROUND_FLOOR)
        return f'{d:.{digits-1}E}'

r=json.loads((ROOT/'results/replay.json').read_text())
pairs=[x for x in r['results'] if x['type']=='positive_pair']
rows=[]
for item in pairs:
    z=item['result']
    rows.append(dict(case=item['case'],r=z['r'],h=z['h'],c=z['c'],
      omega_lo=dec(z['omega_lo'],10),omega_hi=dec(z['omega_hi'],10,True),
      response_norm_lo=sci(z['norm_lo'],9),response_norm_hi=sci(z['norm_hi'],9,True),
      C2_gap=z['target_gaps']['2'],C3_gap=z['target_gaps']['3'],C4_gap=z['target_gaps']['4'],
      KL_upper=sci(z['kl_upper'],9,True),eligible=bool(z['honest_expected_length_lower']),
      honest_C2_expected_length_lower=z['honest_expected_length_lower'].get('2','')))
with (ROOT/'results/summary.csv').open('w') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
example=next(x['result'] for x in pairs if x['case']=='r3_h3')
summary=dict(pair_count=len(rows),operator_count=4,target_pairs=3*len(rows),
  KL_eligible_pairs=sum(x['eligible'] for x in rows),KL_ineligible_pairs=sum(not x['eligible'] for x in rows),
  maximum_raw_component_residual_ratio_upper=dec(max(F(x['result']['raw_peak_ratio']) for x in pairs),12,True),
  maximum_relative_norm_interval_width_upper=sci(max(F(x['result']['relative_norm_bracket']) for x in pairs),9,True),
  example_r3_h3=example,rows=rows,
  limitations='No estimator or Monte Carlo coverage experiment; exact model-pair certificates and analytic adaptive-KL consequences only.')
(ROOT/'results/summary.json').write_text(json.dumps(summary,indent=2)+'\n')
# TeX body used directly in the evidence document, so table values trace to exact replay.
lines=[]
for x in rows:
    def tx(y):
        a,b=y.split('E');return a+r'\times10^{'+str(int(b))+'}'
    z=next(i['result'] for i in pairs if i['case']==x['case'])
    olo,ohi=dec(z['omega_lo'],6),dec(z['omega_hi'],6,True)
    line=f"{x['r']} & ${x['h']}$ & ${x['c']}$ & $[{olo},{ohi}]$ & ${tx(x['response_norm_hi'])}$ & {'yes' if x['eligible'] else 'no'} "
    lines.append(line+chr(92)*2)
(ROOT/'documents/number_table.tex').write_text('\n'.join(lines)+'\n')
print({k:v for k,v in summary.items() if k not in ['rows','example_r3_h3']})
print('example',next(x for x in rows if x['case']=='r3_h3'))
