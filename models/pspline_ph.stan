data {
  int<lower=1> N;
  int<lower=1> P;
  int<lower=3> K;
  int<lower=1> R;
  int<lower=1> J;
  int<lower=2> L;
  int<lower=1> C;
  matrix[N,P] X;
  array[N] int<lower=0,upper=1> event;
  matrix[N,K] B;
  matrix[J*L,K] Bq;
  vector[J*L] log_weights;
  array[N+1,C] int<lower=1,upper=J> piece_ids;
  array[N+1] int<lower=1,upper=C> piece_n;
  matrix[K,R] F;
  vector[K] trend;
  real log_time_ratio;
}
parameters {
  real alpha;
  vector[P] beta;
  real slope;
  real log_tau;
  vector[R] z_smooth;
}
transformed parameters {
  vector[K] gamma = slope * trend + exp(log_tau) * F * z_smooth;
  vector[J] log_integrals;
  vector[N] eta = X * beta;
  vector[N] log_H;
  real log_norm;
  {
    vector[J*L] lq = Bq * gamma + log_weights;
    for (j in 1:J) log_integrals[j] = log_sum_exp(segment(lq,(j-1)*L+1,L));
    log_norm = log_sum_exp(log_integrals[piece_ids[N+1,1:piece_n[N+1]]]);
    for (i in 1:N) log_H[i] = alpha + eta[i] - log_norm
      + log_sum_exp(log_integrals[piece_ids[i,1:piece_n[i]]]);
  }
}
model {
  alpha ~ normal(0,1.2642168);
  beta ~ normal(0,.7);
  slope ~ normal(0,.7);
  log_tau ~ normal(log(.3),.5);
  z_smooth ~ std_normal();
  for (i in 1:N)
    target += event[i] * (alpha + eta[i] + B[i] * gamma - log_norm + log_time_ratio)
      - exp(log_H[i]);
}
generated quantities {
  vector[N] log_lik;
  for (i in 1:N) log_lik[i] = event[i] * (alpha + eta[i] + B[i] * gamma - log_norm + log_time_ratio)
    - exp(log_H[i]);
}
