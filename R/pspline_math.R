ps_basis<-function(t,spec){if(!length(t))return(matrix(numeric(),0,spec$K));as.matrix(splines::bs(pmax(0,pmin(spec$upper,t)),knots=spec$knots,Boundary.knots=c(0,spec$upper),degree=3,intercept=TRUE))}
ps_gauss<-function(n){j<-seq_len(n-1L);A<-matrix(0,n,n);A[cbind(j,j+1L)]<-j/sqrt(4*j*j-1);A<-A+t(A);e<-eigen(A,symmetric=TRUE);i<-order(e$values);list(x=e$values[i],w=2*e$vectors[1,i]^2)}
ps_quad<-function(times,spec,order=16L){
 # Unique full/partial knot intervals, reused across patients and reference.
 g<-ps_gauss(order);cuts<-c(0,spec$knots,spec$upper);intervals<-list();lookup<-new.env(parent=emptyenv());pieces<-vector('list',length(times));tail<-pmax(0,times-spec$upper)
 for(i in seq_along(times)){
  t<-min(spec$upper,max(0,times[i]));ids<-integer()
  for(s in seq_len(length(cuts)-1L)){
   a<-cuts[s];b<-min(t,cuts[s+1L]);if(b<=a)next
   key<-paste(sprintf('%.17g',a),sprintf('%.17g',b),sep=':')
   if(!exists(key,envir=lookup,inherits=FALSE)){intervals[[length(intervals)+1L]]<-c(a,b);assign(key,length(intervals),envir=lookup)}
   ids<-c(ids,get(key,envir=lookup,inherits=FALSE))
  };pieces[[i]]<-ids
 }
 nodes<-unlist(lapply(intervals,function(v)(v[1]+v[2])/2+(v[2]-v[1])*g$x/2),use.names=FALSE)
 weights<-unlist(lapply(intervals,function(v)(v[2]-v[1])*g$w/2),use.names=FALSE)
 ids<-matrix(1L,length(times),length(cuts)-1L);for(i in seq_along(pieces))if(length(pieces[[i]]))ids[i,seq_along(pieces[[i]])]<-pieces[[i]]
 list(B=ps_basis(nodes,spec),log_weights=log(as.numeric(weights)),piece_ids=ids,piece_n=vapply(pieces,length,integer(1)),L=as.integer(order),J=length(intervals),tail=tail,at_end=ps_basis(spec$upper,spec))
}
ps_rows_lse<-function(x){if(is.null(dim(x)))return(x);m<-apply(x,1,max);m+log(rowSums(exp(x-m)))}
ps_gamma<-function(p,spec){
 if(is.null(dim(p)))p<-matrix(p,nrow=1,dimnames=list(NULL,names(p)))
 gn<-paste0('gamma[',seq_len(spec$K),']');if(all(gn%in%colnames(p)))return(p[,gn,drop=FALSE])
 zz<-p[,paste0('z_smooth[',seq_len(ncol(spec$F)),']'),drop=FALSE]
 sweep(zz%*%t(spec$F),1,exp(p[,'log_tau']),'*')+outer(p[,'slope'],spec$trend)
}
ps_integrals<-function(G,q){
 if(q$J==0L)return(matrix(0,nrow(G),nrow(q$piece_ids)))
 lq<-G%*%t(q$B);lq<-sweep(lq,2,q$log_weights,'+');nr<-nrow(G)
 ints<-matrix(0,nr,q$J)
 for(j in seq_len(q$J)){ix<-seq.int((j-1L)*q$L+1L,j*q$L);ints[,j]<-ps_rows_lse(lq[,ix,drop=FALSE])}
 result<-matrix(0,nr,nrow(q$piece_ids))
 for(i in seq_len(ncol(result))){if(q$piece_n[i])result[,i]<-exp(ps_rows_lse(ints[,q$piece_ids[i,seq_len(q$piece_n[i])],drop=FALSE]))}
 if(any(q$tail>0))result<-result+outer(exp(as.numeric(G%*%t(q$at_end))),q$tail)
 result
}
ps_ll_matrix<-function(p,d,spec,order=16L,block=128L){
 if(is.null(dim(p)))p<-matrix(p,nrow=1,dimnames=list(NULL,names(p)))
 out<-matrix(NA_real_,nrow(p),length(d$event));q<-ps_quad(c(d$model_time,1),spec,order);B<-ps_basis(d$model_time,spec);N<-length(d$event)
 for(first in seq.int(1L,nrow(p),by=block)){
  ix<-first:min(nrow(p),first+block-1L);v<-p[ix,,drop=FALSE];G<-ps_gamma(v,spec);I<-ps_integrals(G,q);ln<-log(I[,N+1L]);eta<-v[,paste0('beta[',seq_len(ncol(d$X)),']'),drop=FALSE]%*%t(d$X)
  level<-sweep(eta,1,v[,'alpha']-ln,'+');lh<-G%*%t(B)+level+log(d$time_ratio)
  out[ix,]<-sweep(lh,2,d$event,'*')-I[,seq_len(N),drop=FALSE]*exp(level)
 };out
}
ps_profile<-function(v,spec){G<-ps_gamma(v,spec);q<-ps_quad(1,spec,32L);norm<-ps_integrals(G,q)[1,1];list(shape=1,weights=as.numeric(G),normalization=norm)}
ps_baseline<-function(t,pr,spec,basis=NULL){q<-if(is.null(basis))ps_quad(t,spec,32L)else basis;I<-as.numeric(ps_integrals(matrix(pr$weights,nrow=1),q));list(H=I/pr$normalization,log_h=as.numeric(ps_basis(t,spec)%*%pr$weights)-log(pr$normalization))}
ps_inverse<-function(target,pr,spec,sg=NULL){
 if(is.null(sg)){grid<-seq(0,spec$upper,length.out=8193);sg<-list(time=grid,B=ps_basis(grid,spec))}
 grid<-sg$time;h<-exp(as.numeric(sg$B%*%pr$weights))/pr$normalization
 H<-c(0,cumsum(diff(grid)*(head(h,-1)+tail(h,-1))/2));exact_end<-ps_baseline(spec$upper,pr,spec)$H;H<-H*(exact_end/tail(H,1));keep<-!duplicated(H)
 ans<-approx(H[keep],grid[keep],xout=pmin(target,exact_end),rule=2,ties='ordered')$y
 beyond<-target>exact_end;ans[beyond]<-spec$upper+(target[beyond]-exact_end)/tail(h,1)
 check<-ps_baseline(ans,pr,spec,ps_quad(ans,spec,16L))$H
 bad<-!is.finite(check)|abs(check-target)>1e-6*pmax(1,target)
 if(any(bad))for(i in which(bad)){
  if(target[i]>exact_end){ans[i]<-spec$upper+(target[i]-exact_end)/tail(h,1)}else ans[i]<-uniroot(function(t)ps_baseline(t,pr,spec)$H-target[i],c(0,spec$upper),tol=1e-10)$root
 }
 if(any(!is.finite(ans))||any(ans<0))stop('Invalid P-spline inverse');ans
}
