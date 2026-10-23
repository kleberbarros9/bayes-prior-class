#!/usr/bin/env Rscript
source(file.path(dirname(sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1])),'pspline_common.R'))
a<-commandArgs(TRUE);z<-full_load(a[1]);attempt<-as.integer(a[2]);cfg<-full_cfg(attempt);out<-file.path(z$path,'fit.rds')
old<-full_read(out,z$key);if(!is.null(old)&&identical(old$status,'complete')){cat('REUSE FIT\n');quit(status=0)}
objects<-lapply(1:4,function(ch){p<-file.path(z$path,sprintf('attempt_%02d/chain_%02d.rds',attempt,ch));x<-full_read(p,digest::digest(list(z$key,cfg,ch),algo='sha256'));if(is.null(x)||x$status!='complete')stop('Missing completed chain ',ch);x})
n<-dim(objects[[1]]$parameters)[1];pars<-dimnames(objects[[1]]$parameters)[[3]];U<-length(z$d$starts)
ar<-array(NA_real_,c(n,4,length(pars)),dimnames=list(NULL,paste0('chain',1:4),pars));ll<-array(NA_real_,c(n,4,U))
for(ch in 1:4){ar[,ch,]<-objects[[ch]]$parameters[,1,];ll[,ch,]<-objects[[ch]]$log_lik[,1,]}
su<-as.data.frame(posterior::summarise_draws(posterior::as_draws_array(ar)))
# Exclude deterministic simplex[1] if present, otherwise all model parameters.
active<-is.finite(su$sd)&su$sd>0;dg<-su[active,,drop=FALSE]
if(!nrow(dg))stop('No nonconstant parameters')
sampler<-lapply(objects,`[[`,'sampler')
bfmi<-vapply(sampler,function(s)mean(diff(s[,'energy__'])^2)/var(s[,'energy__']),0.0)
diag<-list(rank_Rhat=max(dg$rhat),max_Rhat=max(dg$rhat),ESS_bulk=min(dg$ess_bulk),ESS_tail=min(dg$ess_tail),min_ESS=min(dg$ess_bulk),divergences=sum(vapply(sampler,function(s)sum(s[,'divergent__']),0.0)),depth_hits=sum(vapply(sampler,function(s)sum(s[,'treedepth__']>=cfg$max_treedepth),0.0)),min_BFMI=min(bfmi))
diag$valid<-all(is.finite(unlist(diag)))&&diag$rank_Rhat<=1.01&&diag$ESS_bulk>=cfg$min_ess&&diag$ESS_tail>=cfg$min_ess&&diag$divergences==0&&diag$depth_hits==0&&diag$min_BFMI>=.3
write.csv(su,file.path(z$path,sprintf('diagnostics_attempt_%02d.csv',attempt)),row.names=FALSE)
full_json(list(status=if(diag$valid)'diagnostics_valid'else'diagnostics_failed',model=a[1],attempt=attempt,diagnostics=diag,key=z$key),file.path(z$path,'fit_status.json'))
if(!diag$valid){cat('DIAGNOSTICS FAILED',a[1],attempt,'\n');quit(status=2)}
pp<-matrix(ar,ncol=length(pars),dimnames=list(NULL,pars))
f<-list(status='posterior_saved',diagnostics=diag,chains=ar,parameter_draws=pp,log_lik=ll,iterations=cfg$iterations,history=list(diag))
# Save combined posterior before any metrics: a scoring failure never refits chains.
analysis_atomic(list(key=z$key,status='posterior_saved',model=z$m,cfg=cfg,attempt=attempt,fit=f),out)
f<-analysis_metrics(f,z$d,z$m,cfg);if(!identical(f$status,'complete'))stop('Metric failure: ',f$error)
analysis_atomic(list(key=z$key,status='complete',model=z$m,cfg=cfg,attempt=attempt,fit=f),out)
row<-list(model=a[1],ELPD=sum(f$pointwise_elpd),LOOIC=-2*sum(f$pointwise_elpd),WAIC=unname(f$waic$estimates['waic','Estimate']),DIC=f$DIC,MCSE_ELPD=f$MCSE_total,max_k=max(f$pareto_k),high_k=sum(f$pareto_k>.7),diagnostics=diag,attempt=attempt,iterations=cfg$iterations,chains=4,key=z$key)
full_json(row,file.path(z$path,'metrics.json'));full_json(list(status='complete',key=z$key,attempt=attempt,model=a[1],diagnostics=diag),file.path(z$path,'fit_status.json'))
write.csv(data.frame(unit=seq_along(f$pointwise_elpd),ELPD=f$pointwise_elpd,pareto_k=f$pareto_k),file.path(z$path,'pointwise.csv'),row.names=FALSE)
cat('FIT COMPLETE',a[1],'ELPD',row$ELPD,'\n')
