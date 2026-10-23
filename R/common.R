ROOT <- Sys.getenv("AJS_PROJECT_ROOT")
if (!nzchar(ROOT)) stop("Run through run.py or the supplied scripts.")
ROOT <- normalizePath(ROOT, mustWork=TRUE)
source(file.path(ROOT,"R","engine.R")); source(file.path(ROOT,"R","likelihood.R"))
.surv_root <- ROOT
analysis_atomic <- function(x,path) {
 dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE)
 tmp<-tempfile("partial_",tmpdir=dirname(path));on.exit(unlink(tmp),add=TRUE)
 saveRDS(x,tmp,compress=FALSE);if(!file.rename(tmp,path))stop("Atomic rename failed: ",path)
}
analysis_json <- function(x,path) {
 dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE)
 tmp<-tempfile("partial_",tmpdir=dirname(path));on.exit(unlink(tmp),add=TRUE)
 jsonlite::write_json(x,tmp,auto_unbox=TRUE,pretty=TRUE,digits=16,null="null")
 if(!file.rename(tmp,path))stop("JSON rename failed")
}
analysis_lse <- survival_log_sum
analysis_stan <- function(aft="matched") paste(readLines(file.path(ROOT,"models","survival.stan")),collapse="\n")
analysis_cfg <- function(final=TRUE,test=FALSE) full_cfg()
analysis_loglik <- function(parameters,model,dataset,gh_points=45L) survival_loglik(parameters,model,dataset)
analysis_prior_alpha <- function(alpha,q) {
  if(q==0)return(dnorm(alpha,0,1.2642168,log=TRUE))
  cal<-surv_q_calibration(q); k<-cal$kappa;lambda<-cal$lambda;xi<-cal$baseline_shape;lu<-xi*alpha;u<-exp(lu)
  if(q==1){v<-u+k*u*u/(1+u);ld<-log(1+k-k/(1+u)^2)}
  if(q==2)return(log(lambda*k*xi)+k*lu-lambda*exp(k*lu))
  if(q==3){v<-u+k*u*u;ld<-log1p(2*k*u)}
  if(q==4){v<-expm1(pmin(k*u,700))/k;ld<-k*u}
  if(q==5){v<-u*(1+k*log1p(u));ld<-log(1+k*(log1p(u)+u/(1+u)))}
  log(lambda*xi)+lu+ld-lambda*v
}
analysis_mcse <- function(series) {
  # series: iteration x chain. Covariances across residents already included.
  n<-nrow(series);out<-lapply(c(50L,100L,250L,500L,1000L,2000L),function(b){
    nb<-n%/%b;if(nb<4L)return(NULL)
    means<-unlist(lapply(seq_len(ncol(series)),function(c)colMeans(matrix(series[seq_len(nb*b),c],nrow=b))))
    data.frame(batch=b,batches=length(means),MCSE=sd(means)/sqrt(length(means)))
  });out<-do.call(rbind,out);out
}
analysis_metrics <- function(fit,data,model,cfg) {
  if(!identical(fit$status,'complete') || is.null(fit$loo) || is.null(fit$waic) || is.null(fit$pointwise_elpd)) fit<-surv_compute_metrics(fit,cfg)
  if(fit$status!='complete')return(fit)
  d<-dim(fit$log_lik);ll<-matrix(fit$log_lik,ncol=d[3]);draws<-fit$parameter_draws
  check<-unique(round(seq(1,nrow(draws),length.out=5)))
  err<-max(vapply(check,function(j)max(abs(analysis_loglik(draws[j,],model,data,cfg$gh_points)-ll[j,])),0.0))
  if(!is.finite(err)||err>1e-6)stop('Independent R likelihood check failed: ',err)
  Dbar<- -2*mean(rowSums(ll));Dhat<- -2*sum(analysis_loglik(colMeans(draws),model,data,cfg$gh_points))
  fit$DIC<-2*Dbar-Dhat;fit$p_DIC<-Dbar-Dhat;fit$likelihood_check<-err
  lm<-apply(-ll,2,analysis_lse)-log(nrow(ll));psi<-rowSums(1-exp(sweep(-ll,2,lm,'-')))
  fit$analysis_ELPD<- -sum(lm);fit$MCSE_batches<-analysis_mcse(matrix(psi,nrow=d[1],ncol=d[2]))
  fit$MCSE_total<-if(nrow(fit$MCSE_batches))tail(fit$MCSE_batches$MCSE,1) else NA_real_
  fit$PSIS_analysis_difference<-sum(fit$pointwise_elpd)-fit$analysis_ELPD
  if(isTRUE(data$frailty)) {
    high<-analysis_loglik(colMeans(draws),model,data,45L);low<-analysis_loglik(colMeans(draws),model,data,25L)
    fit$frailty_GH_meanparameter_change<-max(abs(high-low))
  }
  fit
}
