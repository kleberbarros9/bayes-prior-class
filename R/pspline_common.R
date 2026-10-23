source(file.path(dirname(sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1])),'full_lung_common.R'))
source(file.path(ROOT,'R','pspline_math.R'))
.PS_OUT<-file.path(ROOT,'results','pspline_lung');.PS<-readRDS(file.path(.PS_OUT,'spec.rds'));.ps_cache<-new.env(parent=emptyenv())
.original_load<-full_load
ps_old_load<-function(name){saved<-FULL_ROOT;FULL_ROOT<<-file.path(ROOT,'results','full_lung');on.exit(FULL_ROOT<<-saved);.original_load(name)}
FULL_ROOT<-.PS_OUT
full_names<-c('Quadratic-Weibull','QROG-Weibull','Weibull-PH','Exponential-PH','Lognormal-AFT','Loglogistic-AFT','Piecewise-RW1','Mspline-PH','P-spline-PH')
full_cfg<-function(attempt=1L)list(version='pspline_ph_v1',mode='exec',attempts=2L,chains=4L,cores=1L,iterations=16000L*2L^(attempt-1L),warmup=8000L*2L^(attempt-1L),min_ess=4000L,adapt_delta=if(attempt==1L).995 else .999,max_treedepth=13L,seed=20261026L,gh_points=45L,curve_draws=8000L,curve_points=161L,predictive_draws=40000L,predictive_chunk=200L,quadrature=16L)
full_load<-function(name){
 if(name!='P-spline-PH')return(ps_old_load(name))
 z<-ps_old_load('Weibull-PH');if(!identical(.PS$data_hash,digest::digest(z$d,algo='sha256')))stop('Prepared lung data changed')
 co<-readRDS(file.path(.PS_OUT,'compiled.rds'));h<-digest::digest(paste(readLines(file.path(ROOT,'models','pspline_ph.stan')),collapse='\n'),algo='sha256')
 if(!identical(h,co$source_hash))stop('P-spline source changed')
 z$m$family<-8L;z$m$q<-0L;z$m$model<-'P-spline-PH';z$b$block<-999L;z$compiled<-co
 z$key<-digest::digest(list('pspline_ph_v1',.PS,h,z$d$input_hash,z$d$X),algo='sha256');z$path<-file.path(.PS_OUT,'P-spline-PH');dir.create(z$path,recursive=TRUE,showWarnings=FALSE);z
}
.old_stan_data<-surv_stan_data
surv_stan_data<-function(dataset,model,gh_points=45L){
 if(model$family!=8L)return(.old_stan_data(dataset,model,gh_points))
 d<-dataset;q<-ps_quad(c(d$model_time,1),.PS,16L)
 list(N=length(d$event),U=length(d$event),P=ncol(d$X),K=.PS$K,R=ncol(.PS$F),J=q$J,L=q$L,C=ncol(q$piece_ids),X=d$X,event=as.integer(d$event),B=ps_basis(d$model_time,.PS),Bq=q$B,log_weights=q$log_weights,piece_ids=q$piece_ids,piece_n=q$piece_n,F=.PS$F,trend=.PS$trend,log_time_ratio=log(d$time_ratio),family=8L,is_pspline=1L,centered=0L)
}
.old_validate<-surv_validate_inputs
surv_validate_inputs<-function(sd){if(is.null(sd$is_pspline))return(.old_validate(sd));stopifnot(sd$N>0,sd$P==3,all(sd$piece_n>=1),all(is.finite(sd$Bq)),all(sd$event%in%0:1));invisible(TRUE)}
.old_init<-surv_initial_values
surv_initial_values<-function(sd,chain_id=1L){if(is.null(sd$is_pspline))return(.old_init(sd,chain_id));list(alpha=.01*(chain_id-2.5),beta=rep(.01*(chain_id-2.5),sd$P),slope=.01*(chain_id-2.5),log_tau=log(.3),z_smooth=rep(.01*(chain_id-2.5),sd$R))}
.old_names<-surv_parameter_names
surv_parameter_names<-function(sd){if(is.null(sd$is_pspline))return(.old_names(sd));c('alpha','beta','slope','log_tau','z_smooth','gamma')}
.old_rawll<-analysis_loglik
analysis_loglik<-function(parameters,model,dataset,gh_points=45L){if(model$family!=8L)return(.old_rawll(parameters,model,dataset,gh_points));as.numeric(ps_ll_matrix(parameters,dataset,.PS,64L))}
.old_profile<-survival_profile
survival_profile<-function(parameters,family){if(family!=8L)return(.old_profile(parameters,family));ps_profile(parameters,.PS)}
.old_baseline<-survival_baseline
survival_baseline<-function(time,profile,family,basis=NULL){if(family!=8L)return(.old_baseline(time,profile,family,basis));ps_baseline(time,profile,.PS,basis)}
.old_inverse<-survival_inverse
survival_inverse<-function(target,profile,family,spline_grid=NULL){if(family!=8L)return(.old_inverse(target,profile,family,spline_grid));ps_inverse(target,profile,.PS,spline_grid)}
.old_fit<-surv_fit
surv_fit<-function(compiled,dataset,model,cfg,seed,save_posterior=NULL){
 if(model$family!=8L)return(.old_fit(compiled,dataset,model,cfg,seed,save_posterior))
 sd<-surv_stan_data(dataset,model,cfg$gh_points);history<-list()
 for(attempt in 1:2){
  st<-rstan::sampling(compiled,data=sd,chains=4L,cores=cfg$cores,iter=cfg$iterations*2L^(attempt-1L),warmup=cfg$iterations*2L^(attempt-1L)%/%2L,seed=seed+attempt,init=function(chain_id=1L)surv_initial_values(sd,chain_id),refresh=500L,control=list(adapt_delta=if(attempt==1L).995 else .999,max_treedepth=13L,metric='dense_e'))
  dg<-surv_diagnostics(st,cfg,surv_parameter_names(sd));history[[attempt]]<-dg
  if(isTRUE(dg$valid))break
 }
 if(!isTRUE(dg$valid))return(list(status='diagnostics_failed',diagnostics=dg,history=history))
 ar<-rstan::extract(st,pars=surv_parameter_names(sd),permuted=FALSE);pp<-matrix(ar,ncol=dim(ar)[3],dimnames=list(NULL,dimnames(ar)[[3]]))
 ll<-ps_ll_matrix(pp,dataset,.PS,64L);ll<-array(ll,dim=c(dim(ar)[1],4L,length(dataset$event)))
 ans<-list(status='posterior_saved',diagnostics=dg,history=history,chains=ar,parameter_draws=pp,log_lik=ll,iterations=cfg$iterations*2L^(attempt-1L))
 if(!is.null(save_posterior))save_posterior(ans)
 surv_compute_metrics(ans,cfg)
}
