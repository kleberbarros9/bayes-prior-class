root<-Sys.getenv('AJS_PROJECT_ROOT');setwd(root)
for(p in list.files('R',pattern='[.]R$',full.names=TRUE))invisible(parse(p))
packages<-c('rstan','loo','posterior','ggplot2','gridExtra','splines2','digest','jsonlite','survival')
missing<-packages[!vapply(packages,requireNamespace,TRUE,quietly=TRUE)]
if(length(missing))stop('Missing packages: ',paste(missing,collapse=', '),'. Run scripts/setup_linux.sh.')
source('R/common.R');source('R/pspline_math.R');source('R/meta_analysis.R');.univ_root<-root
lung<-read.csv('data/lung.csv');oncology<-univ_load_data()
stopifnot(nrow(lung)==228L,nrow(oncology)==8L)
keep<-complete.cases(lung[,c('age','sex','ph.ecog','time','status')])&lung$time>0
stopifnot(sum(keep)==227L,sum(lung$status[keep]==2L)==164L)
for(q in c(0L,1L,3L)) {
 mass<-integrate(function(a)exp(analysis_prior_alpha(a,q)),-80,12,rel.tol=1e-8)$value
 stopifnot(abs(mass-1)<1e-7)
}
for(f in list.files('models',pattern='[.]stan$',full.names=TRUE)) {
 compiled<-rstan::stanc(file=f,model_name=gsub('[^a-zA-Z0-9]','_',basename(f)))
 if(!isTRUE(compiled$status))stop('Stan syntax failed: ',f)
}
spec<-list(K=7L,upper=3,knots=c(.75,1.5,2.25));times<-c(0,.1,1,2.9,3,8)
G<-matrix(rep(0,7),nrow=1);Q<-ps_quad(times,spec,16L)
stopifnot(max(abs(as.numeric(ps_integrals(G,Q))-times))<1e-11)
pr<-list(shape=1,weights=rep(0,7),normalization=1)
stopifnot(max(abs(ps_inverse(times,pr,spec)-times))<1e-7)
stopifnot(all(is.finite(ps_quad(0,spec)$log_weights)))
cat('R syntax, both Stan models, input snapshots, prior normalization and P-spline quadrature/inversion passed.\n')
write.csv(data.frame(package=packages,version=vapply(packages,function(p)as.character(packageVersion(p)),'')), 'installed_R_packages.csv',row.names=FALSE)
# Compare the normalized zero-spline baseline with a unit exponential model.
time<-c(.1,.5,1,2,4);base<-ps_baseline(time,pr,spec)
stopifnot(max(abs(base$H-time))<1e-11,max(abs(base$log_h))<1e-11)
# Independently check AFT location-to-log-cumulative-hazard Jacobians.
for(a in c(-1,.2,1))for(shape in c(.7,1.3)) {
 f_ln<-function(v)log(-pnorm(-v/shape,log.p=TRUE))
 f_ll<-function(v)log(log1p(exp(shape*v)))
 eps<-1e-5
 for(fam in c('LN','LL')) {
  fn<-if(fam=='LN')f_ln else f_ll
  numerical<-(fn(a+eps)-fn(a-eps))/(2*eps)
  if(fam=='LN')analytic<-exp(dnorm(a/shape,log=TRUE)-log(shape)-pnorm(-a/shape,log.p=TRUE)-f_ln(a))else analytic<-shape*plogis(shape*a)/log1p(exp(shape*a))
  stopifnot(abs(numerical-analytic)<1e-7)
 }
}
cat('Common-target AFT Jacobian and zero-spline hazard identities passed.\n')
