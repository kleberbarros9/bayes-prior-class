#!/usr/bin/env python3
"""Read-only full-workflow monitor. Closing it does not stop the analysis."""
import json,os,time
from pathlib import Path
root=Path(__file__).resolve().parents[1]
def read(p):
 try:return json.loads(p.read_text())
 except (OSError,ValueError):return {}
def bar(done,total):
 n=int(25*done/max(1,total));return '['+'█'*n+'░'*(25-n)+f'] {done}/{total}'
s=read(root/'results/workflow_status.json');print('FULL MANUSCRIPT REPRODUCTION');print('Status:',s.get('status','not started'),'| Section:',s.get('section','—'))
for label,path in [('Simulation','simulation'),('Meta-analysis','meta-analysis'),('Survival-model','survival-model')]:print(label+':',read(root/f'results/{path}/complete.json').get('status','pending'))
for name in ['Quadratic-Weibull','QROG-Weibull','Weibull-PH','Exponential-PH','Lognormal-AFT','Loglogistic-AFT','Piecewise-RW1','Mspline-PH','P-spline-PH']:
 d=root/('results/pspline_lung' if name=='P-spline-PH' else 'results/full_lung')/name
 attempts=sorted(d.glob('attempt_*'));a=attempts[-1] if attempts else None
 chains=sum(read(a/f'chain_{c:02d}.json').get('status')=='complete' for c in range(1,5)) if a else 0
 stages=sum(read(d/(step+'_complete.json')).get('status')=='complete' for step in ['tables','curves','ppc','loo'])
 print(f'{name:20s} chains {bar(chains,4)} | postprocessing {bar(stages,4)}')
print('Bars count completed checkpoints, not elapsed computational cost. No ETA is inferred across unequal stages.')
print('Log: tail -n 40 -F results/full_run.log')
