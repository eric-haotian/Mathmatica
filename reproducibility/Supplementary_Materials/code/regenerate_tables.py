"""Rebuild manuscript tables from frozen FS2 rational outputs; never round a bound inward."""
from pathlib import Path
from fractions import Fraction
import json

ROOT=Path(__file__).resolve().parent
DATA=ROOT/'inherited/fs2/EIS_PEM_Math_FS2_Sequential_Cost/results'
OUT=ROOT/'manuscript'
def upward(value: str, digits: int=9) -> str:
    q=Fraction(value); s=10**digits
    n=-((-q.numerator*s)//q.denominator)
    return f'{n//s}.{n%s:0{digits}d}'

def main() -> None:
    costs=json.loads((DATA/'cost_summary.json').read_text())
    costs=sorted([c for c in costs if c['r']==3 and c['h']=='1/8'],key=lambda c:-Fraction(c['alpha']))
    lines=[r'\begin{table}[tb]',r'\centering',r'\caption{Expected effort lower bounds for the pair in \eqref{eq:example}, $w=\sigma=10^{-3}$. Integer displays are rounded downward. The interval must cover all positive spectra, including the four-atom alternative.}\label{tab:cost}',r'\begin{tabular}{cr}\toprule',r'Coverage $1-\alpha$ & Certified lower bound for $\E_0 C_T$\\\midrule']
    for c in costs:
        q=Fraction(c['expected_effort_lower']); floor=q.numerator//q.denominator
        assert floor==c['displayed_floor']
        lines.append(f"{float(1-Fraction(c['alpha']))*100:.0f}\\% & $>{floor:,}$"+r'\\')
    lines += [r'\bottomrule\end{tabular}',r'\end{table}']
    (OUT/'cost_table.tex').write_text('\n'.join(lines)+'\n')
    lines=[r'\begin{table}[tb]',r'\centering',r'\caption{Two seeded stopping paths at requested width $w=0.05$. Width and excess-budget displays are rounded upward. The common acquisition counts at stages 1, 2, 3 are 868, 5823, 32409 cumulatively.}\label{tab:stopping}',r'\begin{l|r}']
    # ordinary table, not an image or a chart
    lines[-1]=r'\begin{tabular}{lcrrl}\toprule'
    lines += [r'Reference & Stage & Outer width & Numerical excess & Decision\\\midrule']
    for fn,name in [('sequential_two_atoms_20260916.json','Two atoms'),('sequential_three_atoms_20260917.json','Three atoms')]:
        path=json.loads((DATA/fn).read_text())
        for row in path['stages']:
            width=Fraction(row['certified_length'])
            stop=width<=Fraction(path['requested_width'])
            assert ('STOP' in row['status'])==stop
            lines.append(f"{name if row['stage']==1 else ''} & {row['stage']} & {upward(row['certified_length'])} & {upward(row['numerical_excess_budget'])} & {'stop' if stop else 'continue'}"+r'\\')
    lines += [r'\bottomrule\end{tabular}',r'\end{table}']
    (OUT/'stopping_table.tex').write_text('\n'.join(lines)+'\n')
    print('PASS: regenerated 3 cost rows and 6 stopping rows from exact fractions')
if __name__=='__main__': main()
