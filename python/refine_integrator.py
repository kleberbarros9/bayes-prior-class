#!/usr/bin/env python3
"""Independent common-point RQMC verification of six unclustered Weibull priors.
Not a replacement for PSIS diagnostics or validation of the seven other families.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
import argparse,json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import special,stats,optimize
from scipy.stats import qmc
ROOT=Path(__file__).resolve().parents[1];R=-np.log(.05)/np.log(2);XI=np.log(3.2)/np.log(8)
NAMES=['Weibull-PH','QROG-Weibull','Power-Weibull','Quadratic-Weibull','Exponential-q-Weibull','Loglinear-Weibull']
def calibration(q):
 r=3.2
 if q==1:k=(R-r)/(r*r/(1+r)-R/2);q1=1+k/2
 elif q==2:k=np.log(R)/np.log(r);q1=1.
 elif q==3:k=(R-r)/(r*r-R);q1=1+k
 elif q==4:k=optimize.brentq(lambda k:np.expm1(k*r)/np.expm1(k)-R,1e-8,2);q1=np.expm1(k)/k
 elif q==5:k=(R-r)/(r*np.log1p(r)-R*np.log(2));q1=1+k*np.log(2)
 else:return 0,0
 return k,np.log(2)/q1
def prior_alpha(a,q):
 if q==0:return stats.norm.logpdf(a,0,1.2642168)
 k,l=calibration(q);lu=XI*a;u=np.exp(lu)
 with np.errstate(over='ignore',invalid='ignore'):
  if q==1:v=u+k*u*u/(1+u);ld=np.log(1+k-k/(1+u)**2)
  if q==2:return np.log(l*k*XI)+k*lu-l*np.exp(k*lu)
  if q==3:v=u+k*u*u;ld=np.log1p(2*k*u)
  if q==4:v=np.expm1(np.minimum(k*u,700))/k;ld=k*u
  if q==5:v=u*(1+k*np.log1p(u));ld=np.log(1+k*(np.log1p(u)+u/(1+u)))
  return np.log(l*XI)+lu+ld-l*v
class Model:
 def __init__(self,d):
  self.X=d.filter(regex='^x_[0-9]+$').to_numpy();self.lt=np.log(d.model_time.to_numpy());self.event=d.event.to_numpy();self.ratio=float(d.time_ratio.iloc[0]);self.n=len(d);self.dim=self.X.shape[1]+2
 def ll(self,p):
  p=np.atleast_2d(p);a=p[:,0,None];ls=p[:,1,None];nu=np.exp(ls);eta=p[:,2:]@self.X.T
  with np.errstate(over='ignore',invalid='ignore'):
   return self.event*(a+ls+(nu-1)*self.lt+eta+np.log(self.ratio))-np.exp(a+nu*self.lt+eta)
 def prior(self,p,q):
  p=np.atleast_2d(p);return prior_alpha(p[:,0],q)+stats.norm.logpdf(p[:,1],0,.5)+stats.norm.logpdf(p[:,2:],0,.7).sum(axis=1)
 def geometry(self):
  def f(p):return -float(self.ll(p).sum()+self.prior(p,0)[0])
  fit=optimize.minimize(f,np.zeros(self.dim),method='BFGS',options={'gtol':1e-6});mode=fit.x
  h=2e-4;delta=np.eye(self.dim)*h
  H=np.array([[(f(mode+delta[i]+delta[j])-f(mode+delta[i]-delta[j])-f(mode-delta[i]+delta[j])+f(mode-delta[i]-delta[j]))/(4*h*h) for j in range(self.dim)] for i in range(self.dim)])
  if np.linalg.eigvalsh(H).min()<=0:raise RuntimeError('Nonpositive Hessian; independent integration not certified')
  return mode,np.linalg.inv(H)
def run(m,power,reps,seed):
 mode,cov=m.geometry();L=np.linalg.cholesky(cov);records=[];points=[]
 for rep in range(reps):
  u=qmc.Sobol(m.dim+1,scramble=True,seed=seed+rep).random_base2(power);z=special.ndtri(u[:,:m.dim]);c=stats.chi2.ppf(u[:,-1],7)
  p=mode+(z@L.T)/np.sqrt(c[:,None]/7);logprop=stats.multivariate_t.logpdf(p,loc=mode,shape=cov,df=7)
  total=np.empty(len(p))
  for j in range(0,len(p),512):total[j:j+512]=m.ll(p[j:j+512]).sum(axis=1)
  scores={}
  for q in (0,1,3):
   lw=total+m.prior(p,q)-logprop;lw=np.where(np.isfinite(lw),lw,-np.inf);norm=special.logsumexp(lw);wn=np.exp(lw-norm)
   minus=np.full(m.n,-np.inf);plus=np.full(m.n,-np.inf);mu=np.zeros(m.n);second=np.zeros(m.n)
   for j in range(0,len(p),512):
    ll=m.ll(p[j:j+512]);w=wn[j:j+512]
    minus=np.logaddexp(minus,special.logsumexp(lw[j:j+512,None]-ll,axis=0));plus=np.logaddexp(plus,special.logsumexp(lw[j:j+512,None]+ll,axis=0))
    mu+=w@ll;second+=w@(ll*ll)
   loo=norm-minus;scores[q]=loo;waic=-2*np.sum(plus-norm-(second-mu**2));pm=wn@p;Dbar=-2*float(wn@total);Dhat=-2*float(m.ll(pm).sum())
   records.append(dict(power=power,replicate=rep,q=q,model=NAMES[q],ELPD=float(loo.sum()),LOOIC=float(-2*loo.sum()),WAIC=float(waic),DIC=float(2*Dbar-Dhat),ESS=float(np.exp(2*norm-special.logsumexp(2*lw))),mean_alpha=float(pm[0]),mean_log_shape=float(pm[1])))
  for r in records[-3:]:r['delta_vs_Weibull']=float((scores[r['q']]-scores[0]).sum())
  for q,loo in scores.items():
   for i,v in enumerate(loo):points.append(dict(power=power,replicate=rep,q=q,unit=i+1,ELPD=float(v)))
  print('RQMC',power,rep+1,'/',reps,flush=True)
 return records,points
def main():
 a=argparse.ArgumentParser();a.add_argument('dataset');a.add_argument('--powers',default='14,16');a.add_argument('--replicates',type=int,default=8);args=a.parse_args()
 meta=json.loads((ROOT/'prepared'/f'{args.dataset}_python.json').read_text())
 if meta['frailty']:raise SystemExit('Frailty dataset: use shared-prior R verification and direct LOO; this independent integrator does not cover frailty')
 m=Model(pd.read_csv(ROOT/'prepared'/f'{args.dataset}_python.csv'));out=ROOT/'results'/'verification'/args.dataset;out.mkdir(parents=True,exist_ok=True)
 records=[];points=[]
 for power in map(int,args.powers.split(',')):
  r,p=run(m,power,args.replicates,20261022);records+=r;points+=p
  pd.DataFrame(records).to_csv(out/'replicates.csv',index=False)
 pd.DataFrame(points).to_csv(out/'pointwise.csv',index=False)
 df=pd.DataFrame(records);summary=[]
 for (power,q),g in df.groupby(['power','q']):
  summary.append(dict(power=int(power),model=NAMES[q],ELPD=g.ELPD.mean(),MCSE_ELPD=g.ELPD.std(ddof=1)/np.sqrt(len(g)),delta_vs_Weibull=g.delta_vs_Weibull.mean(),MCSE_delta=g.delta_vs_Weibull.std(ddof=1)/np.sqrt(len(g)),WAIC=g.WAIC.mean(),DIC=g.DIC.mean(),minimum_ESS=g.ESS.min()))
 pd.DataFrame(summary).to_csv(out/'summary.csv',index=False);print(pd.DataFrame(summary).to_string(index=False))
if __name__=='__main__':main()
