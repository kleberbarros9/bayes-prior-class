source(file.path(Sys.getenv("AJS_PROJECT_ROOT"),"R","full_lung_common.R"))
dir.create(file.path(ROOT,"prepared"),showWarnings=FALSE)
original<-read.csv(file.path(ROOT,"data","lung.csv"));d<-original
d$time_value<-d$time;d$event<-as.integer(d$status==2L);d$sex<-factor(d$sex,levels=c(1,2))
keep<-complete.cases(d[,c("time_value","event","age","sex","ph.ecog")])&d$time_value>0
excluded<-which(!keep);d<-droplevels(d[keep,,drop=FALSE]);rownames(d)<-NULL
standardization<-lapply(c("age","ph.ecog"),function(n)c(center=mean(d[[n]]),scale=sd(d[[n]])))
names(standardization)<-c("age","ph.ecog")
for(n in names(standardization))d[[n]]<-(d[[n]]-standardization[[n]][1])/standardization[[n]][2]
X<-model.matrix(~age+sex+ph.ecog,d)[,-1,drop=FALSE];N<-nrow(d)
calendar<-list(original_unit="days",analysis_unit="quarters",divisor=91.3125,original_max=max(d$time),transformed_max=max(d$time)/91.3125)
x<-list(name="lung",original=original,analyzed=d,excluded_rows=excluded,formula="~ age + sex + ph.ecog",X=X,
 time=d$time/91.3125,event=d$event,cluster=seq_len(N),starts=seq_len(N),ends=seq_len(N),frailty=FALSE,
 standardization=standardization,time_scale=91.3125,calendar=calendar,model_time=d$time/365.25,
 reference_scale=365.25,time_ratio=.25,input_hash=digest::digest(original,algo="sha256"))
stopifnot(N==227L,sum(x$event)==164L,length(excluded)==1L)
f<-file.path(ROOT,"prepared","lung.rds")
if(file.exists(f)){if(!identical(digest::digest(readRDS(f)),digest::digest(x)))stop("Prepared data changed; preserve previous results in a separate directory")}else analysis_atomic(x,f)
write.csv(data.frame(model_time=x$model_time,event=x$event,time_ratio=.25,x_1=X[,1],x_2=X[,2],x_3=X[,3]),file.path(ROOT,"prepared","lung_python.csv"),row.names=FALSE)
analysis_json(list(frailty=FALSE,columns=colnames(X),time_ratio=.25),file.path(ROOT,"prepared","lung_python.json"))
write.csv(d,file.path(FULL_ROOT,"analyzed_records.csv"),row.names=FALSE)
analysis_json(list(observations=N,deaths=sum(x$event),excluded_rows=excluded,calendar=calendar,standardization=standardization,reference_original_days=365.25,configuration=full_cfg(),session=capture.output(sessionInfo())),file.path(FULL_ROOT,"data_and_configuration.json"))
code<-analysis_stan();h<-digest::digest(code,algo="sha256");f<-file.path(ROOT,"models","compiled.rds")
if(file.exists(f)){if(!identical(readRDS(f)$code_hash,h))stop("Compiled source changed; remove only compiled.rds after preserving incompatible results")}else {
 rstan::rstan_options(auto_write=FALSE);co<-rstan::stan_model(model_code=code,model_name="lung_full_repro_v2",auto_write=FALSE)
 analysis_atomic(list(compiled=co,code_hash=h,aft_prior="matched"),f)
}
for(name in full_names){z<-full_load(name);surv_validate_inputs(surv_stan_data(z$d,z$m,45L))}
cat("Survival preparation complete: 227 patients, eight base models plus P-spline-PH.\n")
