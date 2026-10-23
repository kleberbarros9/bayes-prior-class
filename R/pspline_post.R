#!/usr/bin/env Rscript
source(file.path(dirname(sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1])),'pspline_common.R'))
a<-commandArgs(TRUE);name<-a[1];step<-a[2];z<-full_load(name);obj<-full_read(file.path(z$path,'fit.rds'),z$key)
if(is.null(obj)||obj$status!='complete'||!isTRUE(obj$fit$diagnostics$valid))stop('No valid full fit')
f<-obj$fit;p<-f$parameter_draws;cfg<-obj$cfg;d<-z$d;model<-z$m;n<-length(d$event);years<-full_years(d);ref_years<-d$reference_scale/365.25
key<-digest::digest(list(z$key,obj$attempt,cfg,'post_v1'),algo='sha256')
rows_survival<-function(v,t) {
 pr<-survival_profile(v,model$family);eta<-as.numeric(d$X%*%v[paste0('beta[',seq_len(ncol(d$X)),']')]);tm<-full_model_time(d,t)
 if(model$family==3L)return(pnorm((log(tm)+v['alpha']+eta)/pr$shape,lower.tail=FALSE))
 if(model$family==4L)return(exp(-survival_softplus(pr$shape*(log(tm)+v['alpha']+eta))))
 exp(-exp(v['alpha']+eta)*survival_baseline(tm,pr,model$family)$H)
}
curve_draw<-function(v,grid,basis=NULL) {
 pr<-survival_profile(v,model$family);eta<-as.numeric(d$X%*%v[paste0('beta[',seq_len(ncol(d$X)),']')]);tm<-full_model_time(d,grid)
 if(model$family==3L){ans<-rowMeans(matrix(pnorm(outer(log(pmax(tm,1e-300)),v['alpha']+eta,'+')/pr$shape,lower.tail=FALSE),nrow=length(grid)));ans[tm==0]<-1;return(ans)}
 if(model$family==4L){ans<-rowMeans(exp(-survival_softplus(pr$shape*outer(log(pmax(tm,1e-300)),v['alpha']+eta,'+'))));ans[tm==0]<-1;return(ans)}
 H<-survival_baseline(tm,pr,model$family,basis)$H;rowMeans(exp(-outer(H,exp(v['alpha']+eta))))
}
make_censor<-function() {
 age<-if('age'%in%names(d$analyzed))as.numeric(d$analyzed$age) else d$X[,1]
 sex<-if('sex'%in%names(d$analyzed))as.character(d$analyzed$sex)else rep('all',n)
 ecog<-if('ph.ecog'%in%names(d$analyzed))as.character(d$analyzed$ph.ecog)else rep('all',n)
 br<-unique(quantile(age,seq(0,1,length.out=5),names=FALSE));ag<-if(length(br)>1)as.character(cut(age,br,include.lowest=TRUE))else rep('all',n)
 full<-interaction(sex,ag,ecog,drop=TRUE);mid<-interaction(sex,ag,drop=TRUE)
 labels<-character(n);laws<-list();audit<-list()
 for(i in seq_len(n)) {
  candidates<-list(which(full==full[i]),which(mid==mid[i]),which(sex==sex[i]),seq_len(n))
  k<-which(vapply(candidates,length,integer(1))>=10L)[1];if(is.na(k))k<-4L
  ids<-candidates[[k]];id<-paste(ids,collapse=',');labels[i]<-id
  if(is.null(laws[[id]])){
   rev<-survival::survfit(survival::Surv(years[ids],1L-d$event[ids])~1)
   prob<-pmax(0,c(1,head(rev$surv,-1))-rev$surv);prob<-c(prob,max(0,tail(rev$surv,1)))
   laws[[id]]<-list(times=c(rev$time,max(years)),prob=prob/sum(prob))
   audit[[length(audit)+1]]<-data.frame(law=length(laws),source=c('sex_agequartile_ecog','sex_agequartile','sex','pooled')[k],training_n=length(ids),training_censored=sum(1L-d$event[ids]))
  }
 }
 for(id in names(laws))laws[[id]]$indices<-which(labels==id)
 list(laws=laws,audit=do.call(rbind,audit),note='Reverse-KM conditional approximation; sparse strata (<10 records) fall back from sex x age quartile x ECOG to sex x age quartile, sex, then pooled. Remaining mass at maximum observed follow-up; this assumption is not validated by event-model diagnostics.')
}
inverse_events<-function(target,pr,sg=NULL) survival_inverse(target,pr,model$family,sg)
if(step=='curves') {
 ix<-unique(round(seq(1,nrow(p),length.out=min(cfg$curve_draws,nrow(p)))))
 grid<-sort(unique(c(seq(0,max(years),length.out=cfg$curve_points),c(.5,1,2)[c(.5,1,2)<=max(years)])))
 basis<-if(model$family==6L)survival_splines(full_model_time(d,grid))else if(model$family==8L)ps_quad(full_model_time(d,grid),.PS,32L)else NULL
 dir<-file.path(z$path,'curve_chunks');dir.create(dir,showWarnings=FALSE);parts<-split(seq_along(ix),ceiling(seq_along(ix)/250L))
 for(k in seq_along(parts)){
  path<-file.path(dir,sprintf('chunk_%04d.rds',k));old<-full_read(path,key)
  if(is.null(old)){ii<-parts[[k]];mat<-t(vapply(ix[ii],function(j)curve_draw(p[j,],grid,basis),numeric(length(grid))));if(any(!is.finite(mat))||any(mat<0|mat>1)||any(apply(mat,1,function(x)any(diff(x)>1e-8))))stop('Invalid survival curves');analysis_atomic(list(key=key,index=ii,selected=ix[ii],curves=mat),path)}
  full_json(list(status='curves',done=max(parts[[k]]),total=length(ix),model=name),file.path(z$path,'post_status.json'))
 }
 mat<-do.call(rbind,lapply(seq_along(parts),function(k)full_read(file.path(dir,sprintf('chunk_%04d.rds',k)),key)$curves))
 qu<-apply(mat,2,quantile,c(.025,.5,.975));tab<-data.frame(model=name,time_years=grid,mean=colMeans(mat),lower=qu[1,],median=qu[2,],upper=qu[3,])
 write.csv(tab,file.path(z$path,'survival_curve.csv'),row.names=FALSE);analysis_atomic(list(key=key,selected=ix,grid=grid,curves=mat),file.path(z$path,'curves.rds'))
 full_json(list(status='complete',key=key,draws=length(ix)),file.path(z$path,'curves_complete.json'))
}
if(step=='ppc') {
 censor<-make_censor();write.csv(censor$audit,file.path(z$path,'censor_audit.csv'),row.names=FALSE)
 analysis_atomic(list(key=key,censor=censor),file.path(z$path,'censor_model.rds'))
 dir<-file.path(z$path,'predictive_chunks');dir.create(dir,showWarnings=FALSE);parts<-split(seq_len(cfg$predictive_draws),ceiling(seq_len(cfg$predictive_draws)/cfg$predictive_chunk))
 sg<-NULL;if(model$family==8L){tt<-seq(0,.PS$upper,length.out=8193);sg<-list(time=tt,B=ps_basis(tt,.PS))};if(model$family==6L){tt<-sort(unique(c(seq(0,20,length.out=8193),.5,1,2)));sg<-list(time=tt,basis=survival_splines(tt))}
 for(k in seq_along(parts)) {
  path<-file.path(dir,sprintf('chunk_%04d.rds',k));old<-full_read(path,key)
  if(is.null(old)){
   set.seed(cfg$seed+20000000L+z$b$block*1000L+k);positions<-parts[[k]];B<-length(positions);selected<-sample.int(nrow(p),B,replace=TRUE)
   latent<-observed<-matrix(NA_real_,B,n);event<-matrix(0L,B,n);future<-numeric(B)
   for(j in seq_len(B)){
    v<-p[selected[j],];pr<-survival_profile(v,model$family);eta<-as.numeric(d$X%*%v[paste0('beta[',seq_len(ncol(d$X)),']')])
    if(model$family==3L)te<-exp(-v['alpha']-eta+rnorm(n,0,pr$shape)) else if(model$family==4L)te<-exp(-v['alpha']-eta+qlogis(runif(n))/pr$shape) else te<-inverse_events(rexp(n)*exp(-v['alpha']-eta),pr,sg)
    te<-te*ref_years;ct<-numeric(n)
    for(law in censor$laws)ct[law$indices]<-sample(law$times,length(law$indices),replace=TRUE,prob=law$prob)
    if(any(!is.finite(te))||any(te<0))stop('Invalid replicated event times')
    latent[j,]<-te;observed[j,]<-pmin(te,ct);event[j,]<-as.integer(te<=ct);future[j]<-te[sample.int(n,1L)]
   }
   stats<-data.frame(replication=positions,parameter_draw=selected,deaths=rowSums(event),followup_person_years=rowSums(observed),future_event_time_years=future)
   analysis_atomic(list(key=key,summary=stats,latent_event_time=latent,observed_time=observed,event=event),path)
  }
  full_json(list(status='ppc',done=max(parts[[k]]),total=cfg$predictive_draws,model=name),file.path(z$path,'post_status.json'))
 }
 pp<-do.call(rbind,lapply(seq_along(parts),function(k)full_read(file.path(dir,sprintf('chunk_%04d.rds',k)),key)$summary))
 con<-gzfile(file.path(z$path,'posterior_predictive.csv.gz'),'wt');write.csv(pp,con,row.names=FALSE);close(con)
 psummary<-function(x,observed,quantity){q<-quantile(x,c(.025,.5,.975));pv<-mean(x>=observed);data.frame(quantity=quantity,observed=observed,mean=mean(x),lower=q[1],median=q[2],upper=q[3],upper_tail_probability=pv,simulation_MCSE_conditional=sqrt(pv*(1-pv)/length(x)))}
 tab<-rbind(psummary(pp$deaths,sum(d$event),'death_count'),psummary(pp$followup_person_years,sum(years),'followup_person_years'))
 write.csv(tab,file.path(z$path,'posterior_predictive_summary.csv'),row.names=FALSE)
 qu<-quantile(pp$future_event_time_years,c(.025,.5,.975));future<-data.frame(model=name,lower=qu[1],median=qu[2],upper=qu[3],Pr_T_gt_0.5=mean(pp$future_event_time_years>.5),Pr_T_gt_1=mean(pp$future_event_time_years>1),Pr_T_gt_2=mean(pp$future_event_time_years>2),Pr_T_gt_5=mean(pp$future_event_time_years>5),Pr_T_gt_10=mean(pp$future_event_time_years>10),observed_max_followup_years=max(years))
 write.csv(future,file.path(z$path,'future_patient_summary.csv'),row.names=FALSE)
 full_json(list(status='complete',key=key,replications=nrow(pp)),file.path(z$path,'ppc_complete.json'))
}
if(step=='tables') {
 ar<-f$chains;su<-as.data.frame(posterior::summarise_draws(posterior::as_draws_array(ar),'mean','sd','rhat','ess_bulk','ess_tail','mcse_mean'));qu<-t(apply(p,2,quantile,c(.025,.5,.975)))
 tab<-data.frame(parameter=colnames(p),mean=colMeans(p),SD=apply(p,2,sd),lower=qu[,1],median=qu[,2],upper=qu[,3],MCSE_mean=su$mcse_mean[match(colnames(p),su$variable)],row.names=NULL)
 write.csv(tab,file.path(z$path,'parameters.csv'),row.names=FALSE)
 derived<-list(baseline_cumulative_hazard_at_reference=if(model$family==3L)-pnorm(-p[,'alpha']/exp(p[,'log_shape[1]']),log.p=TRUE)else if(model$family==4L)survival_softplus(exp(p[,'log_shape[1]'])*p[,'alpha'])else exp(p[,'alpha']))
 if('log_shape[1]'%in%colnames(p))derived$shape_or_logtime_SD<-exp(p[,'log_shape[1]'])
 for(j in seq_len(ncol(d$X))){label<-colnames(d$X)[j];scale<-1;if(label%in%names(d$standardization))scale<-unname(d$standardization[[label]]['scale']);mult<-if(label=='age')10/scale else 1/scale;sign<-if(model$family%in%c(3L,4L))-1 else 1;derived[[paste0(if(model$family%in%c(3L,4L))'time_ratio_'else'HR_',label,if(label=='age')'_per_10_years'else'_per_original_unit')]]<-exp(p[,paste0('beta[',j,']')]*mult*sign)}
 tt<-do.call(rbind,lapply(names(derived),function(nm){v<-derived[[nm]];q<-quantile(v,c(.025,.5,.975));data.frame(parameter=nm,mean=mean(v),SD=sd(v),lower=q[1],median=q[2],upper=q[3],MCSE_mean=posterior::mcse_mean(matrix(v,nrow=dim(ar)[1],ncol=4)))}))
 write.csv(tt,file.path(z$path,'derived_parameters.csv'),row.names=FALSE)
 write.csv(d$analyzed,file.path(FULL_ROOT,'analyzed_records.csv'),row.names=FALSE)
 full_json(list(status='complete',key=key),file.path(z$path,'tables_complete.json'))
}
cat('POST COMPLETE',name,step,'\n')
