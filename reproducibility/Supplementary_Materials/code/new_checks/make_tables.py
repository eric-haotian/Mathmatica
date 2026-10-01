"""Regenerate main-text tables from exact certificates; all bounds round outward."""
from pathlib import Path
from fractions import Fraction as F
import json
ROOT=Path(__file__).resolve().parents[1]
def fixed(n:int,digits:int=6)->str:
    sign='-' if n<0 else ''; n=abs(n); scale=10**digits
    return f'{sign}{n//scale}.{n%scale:0{digits}d}'
def interval(a,b,digits=6):
    a,b=F(a),F(b); scale=10**digits
    lo=(a.numerator*scale)//a.denominator
    hi=-((-b.numerator*scale)//b.denominator)
    return '['+fixed(lo,digits)+', '+fixed(hi,digits)+']'
designs=json.loads((ROOT/'new_checks/results/design_comparison.json').read_text())
lines=[r'\begin{table}[tb]',r'\centering',r'\caption{Same-pair Gaussian information efficiency $e_\Omega$ at equal total effort. Intervals are outward rounded. The oracle knows both hypotheses and is not a spectrum-recovery algorithm.}\label{tab:efficiency}',r'\begin{tabular}{cccc}\toprule',r'$r$ & $h$ & $\Omega_0=(1,\ldots,6)$ & $16\Omega_0$ \\\midrule']
for r in [2,3,4]:
    for den in [8,32]:
        row=next(x for x in designs if x['r']==r and x['h']==f'1/{den}')
        vals=[]
        for key in ['baseline','high_band']:
            d=row['designs'][key]; vals.append(interval(d['oracle_efficiency_lower'],d['oracle_efficiency_upper']))
        lines.append(f"{r} & $1/{den}$ & ${vals[0]}$ & ${vals[1]}$ \\")
lines += [r'\bottomrule\end{tabular}',r'\end{table}']
# Each row needs two literal TeX backslashes.
lines=[line+'\\' if line.endswith(' \\') and not line.endswith(' \\\\') else line for line in lines]
(ROOT/'manuscript/efficiency_table.tex').write_text('\n'.join(lines)+'\n')
rp=json.loads((ROOT/'new_checks/results/replay.json').read_text())
rows={x['result']['input']:x['result'] for x in rp['results'] if x['kind']=='guard'}
lines=[r'\begin{table}[tb]',r'\centering',r'\caption{Sufficient guard-grid nodes on unchanged archived stage inputs. Linear counts are calculated, not executed. All six quadratic grids produced exactly verified strict anchors.}\label{tab:guard}',r'\begin{tabular}{crrrr}\toprule',r' & \multicolumn{2}{c}{Two-atom data} & \multicolumn{2}{c}{Three-atom data} \\',r'Stage & Linear & Quadratic & Linear & Quadratic \\\midrule']
for stage in [1,2,3]:
    a=rows[f'two_atoms_20260916_stage{stage}.json']; b=rows[f'three_atoms_20260917_stage{stage}.json']
    lines.append(f"{stage} & {a['old_nodes']:,} & {a['nodes']:,} & {b['old_nodes']:,} & {b['nodes']:,} \\")
lines += [r'\bottomrule\end{tabular}',r'\end{table}']
lines=[line+'\\' if line.endswith(' \\') and not line.endswith(' \\\\') else line for line in lines]
(ROOT/'manuscript/guard_table.tex').write_text('\n'.join(lines)+'\n')
print('wrote two tables')
