source(file.path(dirname(sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1])),'pspline_common.R'))
z<-full_load('P-spline-PH');o<-full_read(file.path(z$path,'fit.rds'),z$key);if(is.null(o)||o$status!='complete')stop('No accepted P-spline fit')
marker<-file.path(z$path,'quadrature_complete.json');if(file.exists(marker)){cat('REUSE QUADRATURE\n');quit(status=0)}
p<-o$fit$parameter_draws;native<-matrix(o$fit$log_lik,ncol=length(z$d$event))
ll64<-ps_ll_matrix(p,z$d,.PS,64L)
ix<-unique(c(round(seq(1,nrow(p),length.out=min(512,nrow(p)))),order(abs(p[,'log_tau']),decreasing=TRUE)[1:min(64,nrow(p))]))
ll32<-ps_ll_matrix(p[ix,,drop=FALSE],z$d,.PS,32L)
point<-max(abs(ll64-native));total<-max(abs(rowSums(ll64-native)));higher<-max(abs(ll32-ll64[ix,,drop=FALSE]))
audit<-list(status='checked',draws_checked=nrow(p),sampling_order=16L,scoring_order=64L,max_pointwise_difference=point,max_total_loglik_difference=total,max_32_vs_64=higher,valid=is.finite(point)&&point<1e-8&&total<1e-7&&higher<1e-9,DIC_parameterization='Posterior mean of alpha, beta and physical spline coefficients gamma')
full_json(audit,file.path(z$path,'quadrature_audit.json'));if(!isTRUE(audit$valid))stop('Quadrature accuracy failed; original fit preserved, no final comparison authorized')
# Retain native chain RDS; replace only this new model's combined scoring object.
o$fit$log_lik<-array(ll64,dim=dim(o$fit$log_lik));o$fit$status<-'posterior_saved';o$fit$loo<-NULL;o$fit$waic<-NULL;o$fit$pointwise_elpd<-NULL
o$fit<-analysis_metrics(o$fit,z$d,z$m,o$cfg);if(o$fit$status!='complete')stop('High-order scoring failed')
o$quadrature<-audit;analysis_atomic(o,file.path(z$path,'fit.rds'))
f<-o$fit;row<-list(model='P-spline-PH',ELPD=sum(f$pointwise_elpd),LOOIC=-2*sum(f$pointwise_elpd),WAIC=unname(f$waic$estimates['waic','Estimate']),DIC=f$DIC,MCSE_ELPD=f$MCSE_total,max_k=max(f$pareto_k),high_k=sum(f$pareto_k>.7),diagnostics=f$diagnostics,attempt=o$attempt,iterations=o$cfg$iterations,chains=4,quadrature=audit,key=z$key)
full_json(row,file.path(z$path,'metrics.json'));write.csv(data.frame(unit=seq_along(f$pointwise_elpd),ELPD=f$pointwise_elpd,pareto_k=f$pareto_k),file.path(z$path,'pointwise.csv'),row.names=FALSE)
# Check inverse sampling at fixed probabilities and posterior tail draws.
ri<-unique(round(seq(1,nrow(p),length.out=min(100,nrow(p)))));err<-0
sg<-list(time=seq(0,.PS$upper,length.out=8193));sg$B<-ps_basis(sg$time,.PS)
for(j in ri){pr<-ps_profile(p[j,],.PS);target<- -log(c(.999,.9,.5,.1,.01,.001));t<-ps_inverse(target,pr,.PS,sg);H<-ps_baseline(t,pr,.PS,ps_quad(t,.PS,64L))$H;err<-max(err,max(abs(H-target)/pmax(1,target)))}
full_json(list(status='complete',draws_checked=length(ri),max_scaled_H_error=err,tolerance=1e-6),file.path(z$path,'inverse_audit.json'));if(err>1e-6)stop('Inverse audit failed')
full_json(list(status='complete',key=z$key,quadrature=audit),marker);cat('QUADRATURE AND INVERSE COMPLETE\n')
