root<-normalizePath(commandArgs(TRUE)[1],mustWork=TRUE)
lib<-file.path(root,'.Rlibrary');dir.create(lib,showWarnings=FALSE);.libPaths(c(lib,.libPaths()))
options(repos=c(CRAN='https://cloud.r-project.org'),timeout=1200)
if(getRversion()<'4.3.0')stop('R >= 4.3 is required; R 4.5.0 is the documented container version.')
if(!requireNamespace('renv',quietly=TRUE))install.packages('renv',lib=lib)
# Root package versions match the archived full-analysis session where available.
pkgs<-c('RcppParallel@5.1.7','rstan@2.32.7','loo@2.8.0','posterior@1.6.1','ggplot2@3.5.1','gridExtra@2.3','splines2@0.5.4','digest@0.6.37','jsonlite@1.9.1','survival@3.8-3')
lock<-file.path(root,'renv.lock')
if(file.exists(lock))renv::restore(project=root,lockfile=lock,library=lib,prompt=FALSE)else {
 renv::install(pkgs,library=lib,prompt=FALSE)
 renv::snapshot(project=root,library=lib,lockfile=lock,type='all',prompt=FALSE)
}
write.csv(as.data.frame(installed.packages(lib.loc=.libPaths()))[,c('Package','Version')],file.path(root,'installed_R_packages.csv'),row.names=FALSE)
