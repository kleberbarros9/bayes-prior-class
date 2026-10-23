functions {
  real generated_alpha_lpdf(real alpha, int q, real kappa, real lambda, real baseline_shape) {
    real log_u = baseline_shape * alpha;
    real u = exp(log_u);
    real value;
    real derivative = 1;
    if (q == 1) {
      value = u * (1 + kappa * inv_logit(log_u));
      derivative = 1 + kappa - kappa / square(1 + u);
    } else if (q == 2) {
      value = exp(kappa * log_u);
      return log(lambda) + log(kappa) + kappa * log_u + log(baseline_shape) - lambda * value;
    } else if (q == 3) {
      value = u + kappa * square(u);
      derivative = 1 + 2 * kappa * u;
    } else if (q == 4) {
      if (kappa * u > 700) return negative_infinity();
      value = expm1(kappa * u) / kappa;
      return log(lambda) + kappa * u + log(baseline_shape) + log_u - lambda * value;
    } else {
      value = u * (1 + kappa * log1p(u));
      derivative = 1 + kappa * (log1p(u) + u / (1 + u));
    }
    return log(lambda) + log(derivative) + log(baseline_shape) + log_u - lambda * value;
  }
  real record_loglik(int family, real t, int event, real alpha, real shape,
      real eta, row_vector M, row_vector I, row_vector exposure,
      int bin, vector hazard_weights, real normalization) {
    real log_t = log(t);
    if (family == 1 || family == 2) {
      real log_h = alpha + log(shape) + (shape - 1) * log_t + eta;
      real H = exp(alpha + shape * log_t + eta);
      return (event == 1 ? log_h : 0) - H;
    } else if (family == 3) {
      real mu = -alpha - eta;
      if (event == 1) return normal_lpdf(log_t | mu, shape) - log_t;
      return normal_lccdf(log_t | mu, shape);
    } else if (family == 4) {
      real z = shape * (log_t + alpha + eta);
      if (event == 1) return log(shape) - log_t + z - 2 * log1p_exp(z);
      return -log1p_exp(z);
    } else if (family == 5) {
      real log_h = alpha + log(hazard_weights[bin]) + eta - log(normalization);
      real H = exp(alpha + eta) * dot_product(exposure, hazard_weights) / normalization;
      return (event == 1 ? log_h : 0) - H;
    } else if (family == 6) {
      real h = dot_product(M, hazard_weights) / normalization;
      real H = exp(alpha + eta) * dot_product(I, hazard_weights) / normalization;
      return (event == 1 ? alpha + eta + log(h) : 0) - H;
    }
    reject("Unsupported survival family");
    return negative_infinity();
  }
}
data {
 int<lower=1> N; int<lower=1> U; int<lower=1> P; int<lower=0,upper=1> S;
 int<lower=1> K; int<lower=1,upper=6> family; int<lower=0,upper=5> q;
 vector<lower=0>[N] time; array[N] int<lower=0,upper=1> event; matrix[N,P] X;
 matrix[N,K] M; matrix[N,K] I; matrix[N,K] exposure;
 vector[K] reference; array[N] int<lower=1,upper=K> bin;
 real<lower=0> kappa; real<lower=0> lambda; real<lower=0> baseline_shape;
 real<lower=0> time_ratio;
}
parameters {
 real alpha; vector[P] beta; vector[S] log_shape;
 vector[K-1] increments; vector<lower=0>[K>1 ? 1:0] smooth_sd;
}
transformed parameters {
 real shape=1; vector[N] eta=X*beta;
 vector[K] log_profile=rep_vector(0,K); vector[K] hazard_weights=rep_vector(1,K);
 real normalization=1;
 if(S==1) shape=exp(log_shape[1]);
 if(K>1) {
  for(k in 2:K) log_profile[k]=log_profile[k-1]+smooth_sd[1]*increments[k-1];
  if(family==6) hazard_weights=softmax(log_profile); else hazard_weights=exp(log_profile);
  normalization=dot_product(reference,hazard_weights);
 }
}
model {
 if(q>0) target+=generated_alpha_lpdf(alpha | q,kappa,lambda,baseline_shape);
 else if(family==3) {
  real logS=normal_lcdf(-alpha/shape | 0,1); real logH=log(-logS);
  target+=normal_lpdf(logH | 0,1.2642168)+normal_lpdf(alpha/shape | 0,1)-log(shape)-logS-logH;
 } else if(family==4) {
  real z=shape*alpha; real logH=log(log1p_exp(z));
  target+=normal_lpdf(logH | 0,1.2642168)+log(shape)+log_inv_logit(z)-logH;
 } else alpha~normal(0,1.2642168);
 beta~normal(0,.7); log_shape~normal(0,.5); increments~std_normal(); smooth_sd~normal(0,.5);
 for(i in 1:N) target+=record_loglik(family,time[i]*time_ratio,event[i],alpha,shape,eta[i],M[i],I[i],exposure[i],bin[i],hazard_weights,normalization)+event[i]*log(time_ratio);
}
generated quantities {
 vector[N] log_lik;
 for(i in 1:N) log_lik[i]=record_loglik(family,time[i]*time_ratio,event[i],alpha,shape,eta[i],M[i],I[i],exposure[i],bin[i],hazard_weights,normalization)+event[i]*log(time_ratio);
}
