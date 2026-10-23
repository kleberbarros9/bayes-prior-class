"""Record deviations from previous runs; Monte Carlo trajectories need not match."""
from pathlib import Path
import pandas as pd,numpy as np
root=Path(__file__).resolve().parents[1]
a=pd.read_csv(root/'results/survival-model/model_comparison.csv');b=pd.read_csv(root/'reference_outputs/lung_model_comparison.csv')
m=a.merge(b,on='model',suffixes=('_new','_reference'))
for metric in ['ELPD','WAIC','DIC']:m['change_'+metric]=m[metric+'_new']-m[metric+'_reference']
m['combined_ELPD_MCSE']=np.hypot(m.MCSE_ELPD_new,m.MCSE_ELPD_reference)
m.to_csv(root/'results/survival-model/reference_comparison.csv',index=False)
print(m[['model','change_ELPD','combined_ELPD_MCSE']].to_string(index=False))
