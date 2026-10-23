#!/usr/bin/env Rscript
source(file.path(dirname(sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1])),'full_lung_common.R'))
a<-commandArgs(TRUE);z<-full_load(a[1]);obj<-full_read(file.path(z$path,'fit.rds'),z$key);if(is.null(obj)||obj$status!='complete')stop('No valid full fit')
f<-obj$fit;units<-which(f$pareto_k>.7);if(length(units)>20L)stop('More than 20 high-k patients: targeted review needed; no final ranking certified.');out<-file.path(z$path,'direct_loo');dir.create(out,showWarnings=FALSE)
subset_data<-function(d,keep){d$X<-d$X[keep,,drop=FALSE];for(n in c('time','event','model_time'))d[[n]]<-d[[n]][keep];d$analyzed<-d$analyzed[keep,,drop=FALSE];d$cluster<-seq_along(keep);d$starts<-seq_along(keep);d$ends<-seq_along(keep);d}
rows<-list()
for(u in units){
 key<-digest::digest(list(z$key,obj$attempt,u,'direct_loo_v1'),algo='sha256');path<-file.path(out,paste0('unit_',u,'.rds'));old<-full_read(path,key)
 if(!is.null(old)&&old$status=='complete'){rows[[length(rows)+1]]<-old$summary;next}
 train<-subset_data(z$d,setdiff(seq_along(z$d$event),u));test<-subset_data(z$d,u);cfg<-analysis_cfg(final=TRUE);cfg$cores<-1L
 ATTEMPT_DIR<-file.path(out,paste0('unit_',u,'_attempts'))
 save<-function(fit)analysis_atomic(list(key=key,status='posterior_saved',fit=fit),path)
 fit<-if(!is.null(old)&&isTRUE(old$fit$diagnostics$valid))old$fit else surv_fit(z$compiled$compiled,train,z$m,cfg,seed=cfg$seed+90000000L+10000L*z$b$block+u,save_posterior=save)
 if(!isTRUE(fit$diagnostics$valid))stop('Direct LOO diagnostics failed at unit ',u)
 score_ll<-vapply(seq_len(nrow(fit$parameter_draws)),function(j)sum(analysis_loglik(fit$parameter_draws[j,],z$m,test,cfg$gh_points)),0.0)
 score<-analysis_lse(score_ll)-log(length(score_ll));dd<-dim(fit$log_lik);mc<-analysis_mcse(matrix(exp(score_ll-score)-1,nrow=dd[1],ncol=dd[2]));row<-data.frame(unit=u,ELPD=score,MCSE=tail(mc$MCSE,1));analysis_atomic(list(key=key,status='complete',summary=row,fit=fit),path);rows[[length(rows)+1]]<-row
}
point<-f$pointwise_elpd;mc<-f$MCSE_total
if(length(rows)){dr<-do.call(rbind,rows);point[dr$unit]<-dr$ELPD;mc<-mc+sum(dr$MCSE);write.csv(dr,file.path(out,'scores.csv'),row.names=FALSE)}
write.csv(data.frame(unit=seq_along(point),ELPD=point,original_pareto_k=f$pareto_k),file.path(z$path,'pointwise_final.csv'),row.names=FALSE)
full_json(list(status='complete',key=z$key,refits=length(units),ELPD=sum(point),LOOIC=-2*sum(point),MCSE_ELPD_conservative=mc),file.path(z$path,'loo_complete.json'))
cat('LOO COMPLETE',a[1],'direct refits',length(units),'\n')
