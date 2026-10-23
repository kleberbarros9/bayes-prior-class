#!/usr/bin/env python3
"""Run the final, full manuscript analyses. No dataset or transformation search."""
import argparse, concurrent.futures, contextlib, fcntl, hashlib, json, os, shutil, socket, subprocess, sys, time
from pathlib import Path
ROOT=Path(__file__).resolve().parent
BASE=['Quadratic-Weibull','QROG-Weibull','Weibull-PH','Exponential-PH','Lognormal-AFT','Loglogistic-AFT','Piecewise-RW1','Mspline-PH']
ALL=BASE+['P-spline-PH']
def read(path):
    try:return json.loads(path.read_text())
    except (OSError,ValueError):return {}
def atomic(path,data):
    path.parent.mkdir(parents=True,exist_ok=True);tmp=path.with_suffix(path.suffix+'.tmp');tmp.write_text(json.dumps(data,indent=2));os.replace(tmp,path)
def complete(path):return read(path).get('status')=='complete'
def run(command,label,allowed=(0,)):
    log=ROOT/'results/logs'/f'{label}.log';log.parent.mkdir(parents=True,exist_ok=True)
    print('START',label,flush=True)
    with log.open('a') as stream:
        stream.write('\nSTART '+time.strftime('%Y-%m-%d %H:%M:%S %z')+'\n');stream.flush()
        result=subprocess.run(command,cwd=ROOT,stdout=stream,stderr=subprocess.STDOUT)
    if result.returncode not in allowed:raise RuntimeError(f'{label} exited {result.returncode}; inspect {log}')
    print('END',label,result.returncode,flush=True);return result.returncode
def r(script,args=(),label=None,allowed=(0,)):
    return run([os.environ.get('RSCRIPT','Rscript'),'R/'+script,*map(str,args)],label or Path(script).stem,allowed)
def guard(section):
    files=sorted([*ROOT.glob('R/*.R'),*ROOT.glob('models/*.stan'),*ROOT.glob('data/*'),*ROOT.glob('python/*.py'),ROOT/'run.py'])
    key=hashlib.sha256(b''.join(str(p.relative_to(ROOT)).encode()+p.read_bytes() for p in files if p.is_file())).hexdigest()
    path=ROOT/'results'/section/'configuration.json';old=read(path)
    if old and old.get('source_sha256')!=key:raise RuntimeError(f'Source/data changed since {section} started. Preserve results and use a clean checkout/output tree.')
    atomic(path,dict(source_sha256=key,profile='full',created=old.get('created',time.time())))
def chain_paths(name):return ROOT/('results/pspline_lung' if name=='P-spline-PH' else 'results/full_lung')/name
def fit(name,workers):
    out=chain_paths(name);prefix='pspline' if name=='P-spline-PH' else 'full_lung'
    if complete(out/'fit_status.json'):return
    for attempt in range(1,(2 if prefix=='pspline' else 3)+1):
        with concurrent.futures.ThreadPoolExecutor(max_workers=min(workers,4)) as pool:
            jobs=[]
            for ch in range(1,5):
                cp=out/f'attempt_{attempt:02d}/chain_{ch:02d}.json'
                if complete(cp) and cp.with_suffix('.rds').exists():continue
                jobs.append(pool.submit(r,prefix+'_chain.R',[name,attempt,ch],f'{name}_attempt{attempt}_chain{ch}'))
            for job in jobs:job.result()
        rc=r(prefix+'_assemble.R',[name,attempt],f'{name}_assemble{attempt}',(0,2))
        if rc==0:return
    raise RuntimeError(f'{name} failed the declared diagnostics after all attempts; results remain saved and no final comparison is certified.')
def post(name):
    prefix='pspline' if name=='P-spline-PH' else 'full_lung';out=chain_paths(name)
    if prefix=='pspline' and not complete(out/'quadrature_complete.json'):r('pspline_verify.R',label='pspline_quadrature_and_inverse')
    if not complete(out/'loo_complete.json'):r(prefix+'_loo.R',[name],name+'_loo')
    for step in ['tables','curves','ppc']:
        if not complete(out/(step+'_complete.json')):r(prefix+'_post.R',[name,step],name+'_'+step)
def survival(workers):
    r('prepare_survival.R')
    for name in BASE:fit(name,workers);post(name)
    r('pspline_prepare.R');fit('P-spline-PH',workers);post('P-spline-PH')
    r('compare_survival.R')
    # Generate figures after all curves exist so every overlay includes all nine models.
    for name in ALL:r('pspline_figures.R',[name],name+'_figures')
    run([sys.executable,'python/refine_precision.py'],'weibull_independent_precision')
    run([sys.executable,'python/compare_reference.py'],'compare_reference_outputs')
def validate():
    r('validate.R')
    run([sys.executable,'tests/test_python.py'],'python_numerical_validation')
    for script in ROOT.glob('python/*.py'):compile(script.read_text(),str(script),'exec')
    print('Validation passed: R syntax, Stan syntax, data rules and numerical identities. No posterior experiment was substituted for a full run.')
def recover_stale_chain_locks():
    for parent in ['full_lung','pspline_lung']:
        for lock in (ROOT/'results'/parent).glob('*/attempt_*/chain_*.rds.lock'):
            owner=read(lock/'owner.json')
            if owner.get('host')!=socket.gethostname() or not isinstance(owner.get('pid'),int):
                raise RuntimeError(f'Cannot establish whether {lock} is stale. Check its recorded process before removing this lock directory only.')
            try:os.kill(owner['pid'],0)
            except ProcessLookupError:shutil.rmtree(lock)
            else:raise RuntimeError(f"Chain process {owner['pid']} is still active. Wait for it to finish before restarting the runner.")
def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--section',choices=['all','simulation','meta-analysis','survival-model'],default='all');p.add_argument('--workers',type=int,default=4);p.add_argument('--validate',action='store_true',help='Validate code, data and numerical identities without running the scientific experiments');a=p.parse_args()
    if a.workers<1:p.error('--workers must be positive')
    os.environ['AJS_PROJECT_ROOT']=str(ROOT);os.environ['AJS_WORKERS']=str(a.workers)
    for name in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:os.environ[name]='1'
    if not shutil.which(os.environ.get('RSCRIPT','Rscript')):p.error('Rscript is unavailable; run scripts/setup_linux.sh first.')
    if not os.environ.get('R_LIBS_USER'):os.environ['R_LIBS_USER']=str(ROOT/'.Rlibrary')
    out=ROOT/'results';out.mkdir(exist_ok=True)
    with (out/'RUN.lock').open('a+') as lock:
        try:fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError:p.error('Another workflow holds results/RUN.lock. Keep one runner per checkout.')
        if a.validate:validate();return
        recover_stale_chain_locks()
        sections=['simulation','meta-analysis','survival-model'] if a.section=='all' else [a.section]
        state=out/'workflow_status.json';atomic(out/'RUNNING.pid',{'pid':os.getpid()})
        for section in sections:
            guard(section);marker=out/section/'complete.json'
            if complete(marker):print('REUSE completed section:',section,flush=True);continue
            atomic(state,dict(status='running',section=section,pid=os.getpid(),started=time.time()))
            try:
                if section=='simulation':
                    r('run_simulation.R');run([sys.executable,'python/conflict.py'],'directed_conflict')
                elif section=='meta-analysis':r('run_meta_analysis.R')
                else:survival(a.workers)
                atomic(marker,dict(status='complete',finished=time.time()))
            except BaseException as error:
                atomic(state,dict(status='failed',section=section,error=str(error),finished=time.time()));raise
        atomic(state,dict(status='complete',sections=sections,finished=time.time()))
        print('Requested full analyses complete. Run scripts/export_results.sh for compact and full archives.',flush=True)
if __name__=='__main__':main()
