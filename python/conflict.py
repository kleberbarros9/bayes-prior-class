"""Independent revision calculations. Run from any directory with Python 3 + numpy/scipy/pandas.

No original Stan chains are required or reconstructed. Survival inference uses
randomized quasi-Monte Carlo importance integration, with common proposal draws
for paired prior comparisons, not MCMC or the original PSIS calculation.
"""
from pathlib import Path
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
import json, re, hashlib, platform
import numpy as np
import pandas as pd
import scipy
from scipy import stats, special, optimize, integrate
from scipy.stats import qmc

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/simulation/conflict'; OUT.mkdir(parents=True,exist_ok=True)
R=-np.log(.05)/np.log(2)

def calibration(family, r=3.2):
    if family=='QROG':
        if not r <= R < 2*r*r/(1+r): raise ValueError('Inadmissible finite QROG pair')
        k=(R-r)/(r*r/(1+r)-R/2); q1=1+k/2
    elif family=='Quadratic':
        if not r <= R < r*r: raise ValueError('Inadmissible finite quadratic pair')
        k=(R-r)/(r*r-R); q1=1+k
    elif family=='Exponential':
        if r>R: raise ValueError('Inadmissible exponential pair')
        k=optimize.brentq(lambda k:np.expm1(k*r)/np.expm1(k)-R,1e-9,4)
        q1=np.expm1(k)/k
    else: raise ValueError(family)
    return np.log(2)/q1,k

def generated_log_z(z,family,r=3.2,median=1.,upper=8.):
    lam,k=calibration(family,r)
    xi=np.log(r)/np.log(upper/median)
    logu=xi*(np.asarray(z)-np.log(median));u=np.exp(np.clip(logu,-700,300))
    if family=='QROG':
        q=u+k*u*u/(1+u);qp=1+k-k/(1+u)**2; lqp=np.log(qp)
    elif family=='Quadratic': q=u+k*u*u;lqp=np.log1p(2*k*u)
    else:
        ku=k*u;q=np.expm1(np.minimum(ku,700))/k;lqp=ku
    return np.log(lam*xi)+logu+lqp-lam*q

def conventional_z(z,median,upper,family):
    if family=='Lognormal':return stats.norm.logpdf(z,np.log(median),np.log(upper/median)/stats.norm.ppf(.95))
    if family=='Gamma':
        a=optimize.brentq(lambda a:stats.gamma.ppf(.95,a)/stats.gamma.ppf(.5,a)-upper/median,.01,100)
        rate=stats.gamma.ppf(.5,a)/median
        return a*np.log(rate)-special.gammaln(a)+a*z-rate*np.exp(z)
    raise ValueError(family)

def qsummary(z,logw):
    w=np.exp(logw-special.logsumexp(logw)); x=np.exp(z)
    c=np.cumsum(w)-.5*w;qs=np.interp([.025,.5,.975],c,x)
    return w,qs

def conflict():
    z=np.linspace(-18,6,24577);theta=np.exp(z);dz=z[1]-z[0];trap=np.ones(len(z));trap[[0,-1]]=.5
    priors={f:generated_log_z(z,f,3.2,.2,.8) for f in ['QROG','Quadratic','Exponential']}
    priors.update({f:conventional_z(z,.2,.8,f) for f in ['Gamma','Lognormal']})
    rng=np.random.default_rng(20261024);rows=[];pairs=[];replicates=[]
    for prob in [.95,.99,.999]:
        lam,k=calibration('QROG');v=-np.log1p(-prob)/lam
        u=(v-1+np.sqrt((1-v)**2+4*(1+k)*v))/(2*(1+k))
        truth=float(.2*u**(np.log(4)/np.log(3.2)))
        for E in [5.,20.]:
            y=rng.poisson(E*truth,size=2000);metrics={}
            for fam,lp in priors.items():
                cache={}
                for yy in np.unique(y):
                    lw=lp+yy*z-E*theta+np.log(trap*dz);w,qs=qsummary(z,lw)
                    mean=float(w@theta); tail=float(w@(theta>.8));pred=float(w@np.exp(-5*theta))
                    # Held-out future count expected log score under true Poisson law.
                    fut=np.arange(int(stats.poisson.ppf(1-1e-10,5*truth))+1)
                    logpred=np.array([special.logsumexp(lw+kk*np.log(5*theta)-5*theta-special.gammaln(kk+1))-special.logsumexp(lw) for kk in fut])
                    score=float(stats.poisson.pmf(fut,5*truth)@logpred)
                    cache[int(yy)]=[mean-truth,(mean-truth)**2,float(qs[0]<=truth<=qs[2]),qs[2]-qs[0],tail,pred,score]
                arr=np.array([cache[int(v)] for v in y]);metrics[fam]=arr
                rep=pd.DataFrame(arr,columns=['error','squared_error','covered','interval_length','posterior_tail_gt_08','predictive_survival_5','expected_future_log_score'])
                rep['replicate']=np.arange(1,len(y)+1);rep['y']=y;rep['family']=fam
                rep['reference_prior_quantile']=prob;rep['truth']=truth;rep['exposure']=E
                replicates.append(rep)
                vals=arr.mean(axis=0);ses=arr.std(axis=0,ddof=1)/np.sqrt(len(y))
                label=dict(family=fam,reference_prior_quantile=prob,truth=truth,exposure=E,M=len(y),prior_tail_at_truth=float(integrate.trapezoid(np.exp(lp)[theta>truth],theta[theta>truth])/1) if False else float(np.sum(np.exp(lp[theta>truth])*dz)))
                row=dict(**label)
                for k,value,se in zip(['bias','MSE','coverage','interval_length','posterior_tail_gt_08','predictive_survival_5','expected_future_log_score'],vals,ses):row[k]=value;row[k+'_MCSE']=se
                row['RMSE']=np.sqrt(vals[1]);row['RMSE_MCSE']=ses[1]/(2*np.sqrt(vals[1]));rows.append(row)
            for fam,arr in metrics.items():
                delta=arr-metrics['Gamma'];pairs.append(dict(family=fam,reference_prior_quantile=prob,exposure=E,delta_MSE_vs_Gamma=delta[:,1].mean(),MCSE_delta_MSE=delta[:,1].std(ddof=1)/np.sqrt(len(y)),delta_coverage_vs_Gamma=delta[:,2].mean(),MCSE_delta_coverage=delta[:,2].std(ddof=1)/np.sqrt(len(y))))
    pd.DataFrame(rows).to_csv(OUT/'upper_tail_conflict.csv',index=False);pd.DataFrame(pairs).to_csv(OUT/'upper_tail_conflict_paired.csv',index=False)
    pd.concat(replicates,ignore_index=True).to_csv(OUT/'upper_tail_conflict_replicates.csv.gz',index=False,compression='gzip')
    return dict(seed=20261024,M=2000,truth_definition='QROG common calibrated reference quantiles 0.95, 0.99, 0.999; same generated Poisson observations across five methods.')


if __name__ == "__main__":
    print(json.dumps(conflict(), indent=2))
