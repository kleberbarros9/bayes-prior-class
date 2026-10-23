#!/usr/bin/env python3
"""Paired RQMC refinement plus direct held-out posterior integration, resumable."""
import os
for n in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'): os.environ[n]='1'
import argparse,hashlib,importlib.util,json,math
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import special,stats
ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('integrator',ROOT/'python/refine_integrator.py');mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
PAIRS=[(1,0,'QROG minus Weibull'),(3,0,'Quadratic minus Weibull'),(1,3,'QROG minus Quadratic')]
def atomic(p,obj):
 t=p.with_suffix(p.suffix+'.tmp');t.write_text(json.dumps(obj,indent=2,allow_nan=False));os.replace(t,p)
def aggregate(records,out):
 df=pd.DataFrame(records);df.to_csv(out/'replicates.csv',index=False);rows=[]
 for (power,rep),g in df.groupby(['power','replicate']):
  b=g.set_index('q')
  for q,ref,label in PAIRS:
   rows.append(dict(power=power,replicate=rep,comparison=label,delta_ELPD=b.loc[q,'ELPD']-b.loc[ref,'ELPD'],improvement_WAIC=b.loc[ref,'WAIC']-b.loc[q,'WAIC'],improvement_DIC=b.loc[ref,'DIC']-b.loc[q,'DIC']))
 pairs=pd.DataFrame(rows);pairs.to_csv(out/'paired_replicates.csv',index=False);summary=[]
 for (power,label),g in pairs.groupby(['power','comparison']):
  r=dict(power=int(power),comparison=label,replicates=len(g))
  for col in ['delta_ELPD','improvement_WAIC','improvement_DIC']:
   r[col]=g[col].mean();r['MCSE_'+col]=g[col].std(ddof=1)/math.sqrt(len(g)) if len(g)>1 else None
  summary.append(r)
 pd.DataFrame(summary).to_csv(out/'paired_summary.csv',index=False)
 return df,pairs

def direct(m,test_model,unit,power,seed):
 """Integrate the TRAIN posterior using geometry fitted without the held-out record."""
 mode,cov=m.geometry();L=np.linalg.cholesky(cov)
 u=mod.qmc.Sobol(m.dim+1,scramble=True,seed=seed).random_base2(power)
 p=mode+(special.ndtri(u[:,:m.dim])@L.T)/np.sqrt(stats.chi2.ppf(u[:,-1],7)[:,None]/7)
 prop=stats.multivariate_t.logpdf(p,loc=mode,shape=cov,df=7)
 total=np.empty(len(p));held=np.empty(len(p))
 for j in range(0,len(p),512):
  total[j:j+512]=m.ll(p[j:j+512]).sum(axis=1);held[j:j+512]=test_model.ll(p[j:j+512])[:,unit]
 rows=[]
 for q in (0,1,3):
  lw=total+m.prior(p,q)-prop;lw=np.where(np.isfinite(lw),lw,-np.inf);norm=special.logsumexp(lw)
  rows.append(dict(unit=unit+1,q=q,ELPD=float(special.logsumexp(lw+held)-norm),ESS=float(np.exp(2*norm-special.logsumexp(2*lw)))))
 return rows

def main():
 a=argparse.ArgumentParser();a.add_argument('--powers',default='18,20');a.add_argument('--reps',type=int,default=12);a.add_argument('--direct-power',type=int,default=18);a.add_argument('--direct-reps',type=int,default=8);a.add_argument('--direct-units',type=int,default=6);args=a.parse_args()
 source=ROOT/'prepared/lung_python.csv';data=pd.read_csv(source);m=mod.Model(data)
 meta=json.loads((ROOT/'prepared/lung_python.json').read_text())
 if meta.get('frailty'):raise RuntimeError('Frailty not supported')
 out=ROOT/'results/lung_refinement_v1/precision';out.mkdir(parents=True,exist_ok=True)
 key=hashlib.sha256(source.read_bytes()+(ROOT/'python/refine_integrator.py').read_bytes()+Path(__file__).read_bytes()+json.dumps(vars(args),sort_keys=True).encode()).hexdigest()
 records=[];powers=sorted(map(int,args.powers.split(',')))
 for power in powers:
  for rep in range(args.reps):
   cp=out/f'p{power}_r{rep:02d}.json';pp=out/f'p{power}_r{rep:02d}_pointwise.csv.gz'
   if cp.exists():
    o=json.loads(cp.read_text());assert o['key']==key,'Checkpoint mismatch';rows=o['rows']
   else:
    # Same scramble seeds across resolutions: nested points, paired refinement.
    rows,points=mod.run(m,power,1,20300000+rep)
    for r in rows:r['replicate']=rep
    pd.DataFrame(points).to_csv(pp,index=False,compression='gzip');atomic(cp,dict(key=key,rows=rows))
   records+=rows;aggregate(records,out);print('PRECISION',power,rep+1,'/',args.reps,flush=True)
 df,pairs=aggregate(records,out)
 refined=[]
 for label,g in pairs.groupby('comparison'):
  lo=g[g.power==powers[-2]].set_index('replicate') if len(powers)>1 else None;hi=g[g.power==powers[-1]].set_index('replicate')
  for col in ['delta_ELPD','improvement_WAIC','improvement_DIC']:
   v=hi[col];se=v.std(ddof=1)/np.sqrt(len(v));shift=(v-lo[col]) if lo is not None else v*float('nan');ds=shift.std(ddof=1)/np.sqrt(len(v))
   refined.append(dict(comparison=label,metric=col,estimate=v.mean(),MCSE=se,change_from_previous=shift.mean(),MCSE_change=ds,resolution_stable=bool(lo is not None and abs(shift.mean())<=3*max(ds,1e-15)),sign_resolved=bool(abs(v.mean())>3*se)))
 pd.DataFrame(refined).to_csv(out/'resolution_audit.csv',index=False)
 # Check largest pointwise discrepancies between saved MCMC/PSIS and integrated scores.
 pts=pd.concat([pd.read_csv(out/f'p{powers[-1]}_r{r:02d}_pointwise.csv.gz') for r in range(args.reps)])
 avg=pts.groupby(['q','unit']).ELPD.mean().unstack('q');rank=np.zeros(m.n)
 for q,_,name in [(0,0,'Weibull-PH'),(1,0,'QROG-Weibull'),(3,0,'Quadratic-Weibull')]:
  original=pd.read_csv(ROOT/f'results/full_lung/{name}/pointwise_final.csv').ELPD.to_numpy()
  rank=np.maximum(rank,np.abs(original-avg[q].to_numpy()))
 selected=np.argsort(rank)[-args.direct_units:][::-1] if args.direct_units else []
 direct_rows=[]
 for unit in selected:
  train=mod.Model(data.drop(data.index[unit]).reset_index(drop=True))
  for rep in range(args.direct_reps):
   cp=out/f'direct_u{unit+1:03d}_r{rep:02d}.json'
   if cp.exists():
    o=json.loads(cp.read_text());assert o['key']==key;rows=o['rows']
   else:
    rows=direct(train,m,int(unit),args.direct_power,20400000+100*int(unit)+rep);atomic(cp,dict(key=key,rows=rows))
   for r in rows:direct_rows.append(dict(**r,replicate=rep))
   print('DIRECT LOO',unit+1,rep+1,'/',args.direct_reps,flush=True)
 if direct_rows:
  dd=pd.DataFrame(direct_rows);dd.to_csv(out/'direct_loo_replicates.csv',index=False);rr=[]
  for (unit,q),g in dd.groupby(['unit','q']):
   vals=pts[(pts.unit==unit)&(pts.q==q)].ELPD
   rr.append(dict(unit=int(unit),q=int(q),direct_ELPD=g.ELPD.mean(),direct_MCSE=g.ELPD.std(ddof=1)/np.sqrt(len(g)),full_posterior_ELPD=vals.mean(),full_posterior_MCSE=vals.std(ddof=1)/np.sqrt(len(vals)),difference=g.ELPD.mean()-vals.mean()))
  pd.DataFrame(rr).to_csv(out/'direct_loo_check.csv',index=False)
 atomic(out/'complete.json',dict(status='complete',key=key,settings=vars(args),scope='Numerical checks; unresolved resolution changes remain flagged, no forced winner.'))
if __name__=='__main__':main()
