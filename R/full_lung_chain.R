#!/usr/bin/env Rscript
source(file.path(dirname(sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1])),'full_lung_common.R'))
a<-commandArgs(TRUE);name<-a[1];attempt<-as.integer(a[2]);chain<-as.integer(a[3]);z<-full_load(name);cfg<-full_cfg(attempt)
path<-file.path(z$path,sprintf('attempt_%02d',attempt));dir.create(path,showWarnings=FALSE)
key<-digest::digest(list(z$key,cfg,chain),algo='sha256');file<-file.path(path,sprintf('chain_%02d.rds',chain));old<-full_read(file,key)
if(!is.null(old)&&identical(old$status,'complete')){full_json(list(status='complete',model=name,attempt=attempt,chain=chain,reused=TRUE),file.path(path,sprintf('chain_%02d.json',chain)));cat('REUSE',name,attempt,chain,'\n');quit(status=0)}
lock<-paste0(file,'.lock');if(!dir.create(lock,showWarnings=FALSE))stop('Chain locked: ',lock)
analysis_json(list(pid=Sys.getpid(),host=Sys.info()[['nodename']]),file.path(lock,'owner.json'))
tryCatch({
 full_json(list(status='sampling',model=name,attempt=attempt,chain=chain,started=as.character(Sys.time())),file.path(path,sprintf('chain_%02d.json',chain)))
 sd<-surv_stan_data(z$d,z$m,cfg$gh_points);surv_validate_inputs(sd)
 init<-surv_initial_values(sd,chain);init<-init[vapply(init,length,integer(1))>0L]
 st<-rstan::sampling(z$compiled$compiled,data=sd,chains=1L,cores=1L,iter=cfg$iterations,warmup=cfg$warmup,chain_id=chain,
 seed=cfg$seed+10000L*z$b$block+100L*attempt+chain,init=list(init),refresh=500L,
 control=list(adapt_delta=cfg$adapt_delta,max_treedepth=cfg$max_treedepth))
 arr<-rstan::extract(st,pars=surv_parameter_names(sd),permuted=FALSE)
 ll<-rstan::extract(st,pars='log_lik',permuted=FALSE);sam<-rstan::get_sampler_params(st,inc_warmup=FALSE)[[1]]
 if(dim(arr)[1]!=cfg$iterations-cfg$warmup||dim(ll)[3]!=sd$U||any(!is.finite(arr))||any(!is.finite(ll)))stop('Invalid chain dimensions or values')
 analysis_atomic(list(key=key,input_key=z$key,status='complete',model=name,attempt=attempt,chain=chain,cfg=cfg,parameters=arr,log_lik=ll,sampler=sam,session=capture.output(sessionInfo())),file)
 full_json(list(status='complete',model=name,attempt=attempt,chain=chain,retained=dim(arr)[1],divergences=sum(sam[,'divergent__']),finished=as.character(Sys.time())),file.path(path,sprintf('chain_%02d.json',chain)))
 unlink(lock,recursive=TRUE);cat('CHAIN COMPLETE',name,attempt,chain,'\n')
},error=function(e){unlink(lock,recursive=TRUE);full_json(list(status='error',error=conditionMessage(e),model=name,attempt=attempt,chain=chain),file.path(path,sprintf('chain_%02d.json',chain)));stop(e)})
