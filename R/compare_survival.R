source(file.path(Sys.getenv('AJS_PROJECT_ROOT'),'R','pspline_common.R'))
rows<-list();points<-list()
for(name in full_names) {
 z<-full_load(name);o<-full_read(file.path(z$path,'fit.rds'),z$key)
 if(is.null(o)||o$status!='complete'||!isTRUE(o$fit$diagnostics$valid))stop('Missing valid model: ',name)
 final<-jsonlite::read_json(file.path(z$path,'loo_complete.json'),simplifyVector=TRUE)
 if(final$status!='complete')stop('LOO incomplete: ',name)
 f<-o$fit;v<-read.csv(file.path(z$path,'pointwise_final.csv'))$ELPD;points[[name]]<-v;dg<-f$diagnostics
 rows[[name]]<-data.frame(model=name,ELPD=sum(v),LOOIC=-2*sum(v),WAIC=unname(f$waic$estimates['waic','Estimate']),DIC=f$DIC,MCSE_ELPD=final$MCSE_ELPD_conservative,Rhat=dg$rank_Rhat,ESS_bulk=dg$ESS_bulk,ESS_tail=dg$ESS_tail,divergences=dg$divergences,depth_hits=dg$depth_hits,BFMI=dg$min_BFMI)
}
out<-file.path(ROOT,'results','survival-model');dir.create(out,recursive=TRUE,showWarnings=FALSE)
tab<-do.call(rbind,rows);tab<-tab[order(-tab$ELPD),];write.csv(tab,file.path(out,'model_comparison.csv'),row.names=FALSE)
paired<-list()
for(name in setdiff(full_names,'Weibull-PH')) {
 d<-points[[name]]-points[['Weibull-PH']];a<-tab[tab$model==name,];b<-tab[tab$model=='Weibull-PH',]
 paired[[name]]<-data.frame(model=name,reference='Weibull-PH',delta_ELPD=sum(d),paired_patient_SE=sqrt(length(d)*var(d)),independent_run_MCSE=sqrt(a$MCSE_ELPD^2+b$MCSE_ELPD^2))
 write.csv(data.frame(unit=seq_along(d),delta_ELPD=d),file.path(full_load(name)$path,'delta_vs_Weibull.csv'),row.names=FALSE)
}
write.csv(do.call(rbind,paired),file.path(out,'paired_differences.csv'),row.names=FALSE)
tex<-c('\\begin{table}[H]','\\centering\\small','\\caption{Full Bayesian survival comparison for lung.}','\\begin{tabular}{lrrrr}','\\hline','Model & ELPD & LOOIC & WAIC & DIC \\\\','\\hline')
for(i in seq_len(nrow(tab)))tex<-c(tex,sprintf('%s & %.3f & %.3f & %.3f & %.3f \\\\',tab$model[i],tab$ELPD[i],tab$LOOIC[i],tab$WAIC[i],tab$DIC[i]))
writeLines(c(tex,'\\hline','\\end{tabular}','\\end{table}'),file.path(out,'model_comparison.tex'))
writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
print(tab,row.names=FALSE)

writeLines(c("FULL LUNG SURVIVAL ANALYSIS", "Nine fixed final models; no dataset search is executed.",
"The application remains exploratory after the original selection process.",
"All scores use censored observed-data likelihoods in a common time unit.",
"ELPD and LOOIC encode the same criterion; WAIC is a related predictive estimate; DIC is parameterization dependent.",
"Inspect paired patient SE, Monte Carlo SE and the independent precision calculation before interpreting tiny ranks.",
capture.output(print(tab,row.names=FALSE))),file.path(out,'REPORT.txt'))
