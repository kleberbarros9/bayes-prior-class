source(file.path(Sys.getenv("AJS_PROJECT_ROOT"),"R","common.R"))
FULL_ROOT<-file.path(ROOT,"results","full_lung");dir.create(FULL_ROOT,recursive=TRUE,showWarnings=FALSE)
surv_models <- function() data.frame(
 model=c("Quadratic-Weibull","QROG-Weibull","Weibull-PH","Exponential-PH","Lognormal-AFT","Loglogistic-AFT","Piecewise-RW1","Mspline-PH"),
 family=c(2L,2L,2L,1L,3L,4L,5L,6L),q=c(3L,1L,rep(0L,6)),
 proposed=c(TRUE,TRUE,rep(FALSE,6)), stringsAsFactors=FALSE)
full_names<-surv_models()$model
full_cfg<-function(attempt=1L) list(version="lung_full_repro_v2",mode="exec",attempts=3L,chains=4L,cores=1L,
 iterations=16000L*2L^(attempt-1L),warmup=8000L*2L^(attempt-1L),min_ess=4000L,
 adapt_delta=.999,max_treedepth=13L,seed=20261023L,gh_points=45L,curve_draws=8000L,
 curve_points=161L,predictive_draws=40000L,predictive_chunk=200L)
full_load<-function(name) {
 m<-surv_models();idx<-match(name,m$model);if(is.na(idx))stop("Unknown final model: ",name)
 m<-m[idx,,drop=FALSE];d<-readRDS(file.path(ROOT,"prepared","lung.rds"))
 co<-readRDS(file.path(ROOT,"models","compiled.rds"))
 if(!identical(co$code_hash,digest::digest(analysis_stan(),algo="sha256")))stop("Stan source changed; use a new output directory")
 key<-digest::digest(list("lung_full_repro_v2",d$input_hash,d$X,d$event,d$model_time,d$time_ratio,m,co$code_hash),algo="sha256")
 path<-file.path(FULL_ROOT,name);dir.create(path,recursive=TRUE,showWarnings=FALSE)
 b<-data.frame(block=c(51L,49L,55L,54L,56L,57L,58L,59L)[idx],dataset="lung",model=name,model_index=idx)
 list(b=b,d=d,m=m,compiled=co,key=key,path=path)
}
full_json<-analysis_json
full_read<-function(p,key=NULL) {if(!file.exists(p))return(NULL);x<-readRDS(p);if(!is.null(key)&&!identical(x$key,key))stop("Checkpoint inputs changed; preserved: ",p);x}
full_years<-function(d)d$model_time*d$reference_scale/365.25
full_model_time<-function(d,t)t/(d$reference_scale/365.25)
