"""Lightweight numerical identities for the independent calculation."""
import importlib.util,sys
from pathlib import Path
import numpy as np
from scipy import integrate
root=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('conflict',root/'python/conflict.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
z=np.linspace(-60,12,65537)
for family in ['QROG','Quadratic','Exponential']:
 with np.errstate(over='ignore',invalid='ignore'):
  density=np.exp(m.generated_log_z(z,family,3.2,.2,.8))
 assert abs(integrate.trapezoid(density,z)-1)<1e-7
 lam,k=m.calibration(family);assert lam>0 and k>0
print('Python generated-prior normalization and calibration checks passed.')
# Exercise the independent integration API on a validation grid, using source CSV data.
spec=importlib.util.spec_from_file_location('precision',root/'python/refine_integrator.py');precision=importlib.util.module_from_spec(spec);spec.loader.exec_module(precision)
import pandas as pd
d=pd.read_csv(root/'data/lung.csv').dropna(subset=['time','status','age','sex','ph.ecog']).copy()
x=pd.DataFrame({'model_time':d.time/365.25,'event':(d.status==2).astype(int),'time_ratio':.25,
 'x_1':(d.age-d.age.mean())/d.age.std(ddof=1),'x_2':(d.sex==2).astype(int),
 'x_3':(d['ph.ecog']-d['ph.ecog'].mean())/d['ph.ecog'].std(ddof=1)})
model=precision.Model(x);rows,points=precision.run(model,9,1,2201)
assert len(rows)==3 and all(np.isfinite(row['ELPD']) for row in rows)
print('Independent integration API passed on a validation grid; production precision is not certified by this check.')
