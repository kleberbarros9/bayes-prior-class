root<-Sys.getenv('AJS_PROJECT_ROOT');setwd(root)
source(file.path(root,'R','meta_analysis.R'));.univ_root<-root
result<-run_univariate_rover('exec',prediction='pred',mc_draws=40000L)
writeLines(capture.output(sessionInfo()),file.path(root,'results','meta_analysis','oncology','sessionInfo.txt'))
