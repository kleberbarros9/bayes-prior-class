# Installation check only; these draws are never used in a manuscript result.
code<-'parameters {real z;} model {z ~ normal(0,1);}'
model<-rstan::stan_model(model_code=code,model_name='installation_check',auto_write=FALSE)
fit<-rstan::sampling(model,chains=1,cores=1,iter=40,warmup=20,seed=1701,refresh=0,init=list(list(z=0)))
stopifnot(all(is.finite(as.matrix(fit,pars='z'))))
cat('C++ compilation and Stan execution passed (installation check only).\n')
