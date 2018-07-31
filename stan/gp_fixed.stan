
// Fit the hyperparameters of a Gaussian process with an exponentiated
// quadratic kernel and analytically predict out-of-sample observations
functions{
  matrix distance(real[] x1,
                  real[] x2,
                  int N1,
                  int N2,
                  real length_scale) {
    matrix[N1, N2] dist;
    for (n1 in 1:N1) {
      for (n2 in 1:N2) {
        dist[n1,n2] = (x1[n1]- x2[n2])/length_scale/length_scale;
      }
    }       
    return dist;
  }
  matrix dx1_cov_exp_quad(real[] x1d,
                          real[] x2,
                          int N1d,
                          int N2,
                          real alpha,
                          real length_scale) {
    matrix[N1d, N2] dist;
    matrix[N1d, N2] Sigma_reg;
    matrix[N1d, N2] Sigma;
    dist = distance(x1d,x2, N1d, N2, length_scale);
    Sigma_reg = cov_exp_quad(x1d, x2, alpha, length_scale);
    Sigma = -1. * dist .* Sigma_reg;
    return Sigma;
  }
  matrix dx1dx2_cov_exp_quad(real[] x1d,
                             real[] x2d,
                             int N1d,
                             int N2d,
                             real alpha,
                             real length_scale) {
    matrix[N1d, N2d] dist;
    matrix[N1d, N2d] Sigma_reg;
    matrix[N1d, N2d] Sigma;
    dist = distance(x1d,x2d, N1d, N2d, length_scale);
    Sigma_reg = cov_exp_quad(x1d, x2d, alpha, length_scale);
    Sigma = -1*Sigma_reg .* dist .* dist + Sigma_reg ./ length_scale ./ length_scale;
    return Sigma;
  }
  matrix cov_exp_quad_full(real[] x,
                           real[] xd,
                           int N,
                           int Nd,
                           real alpha,
                           real length_scale) {
    matrix[N + Nd, N + Nd] Sigma;
    Sigma[1:N,1:N] = cov_exp_quad(x,alpha,length_scale);
    Sigma[N+1:N+Nd, 1:N] = dx1_cov_exp_quad(xd,x, Nd,N,alpha,length_scale);
    Sigma[1:N, N+1:N+Nd] = -1.0*dx1_cov_exp_quad(x,xd,N,Nd,alpha,length_scale);
    Sigma[N+1:N+Nd, N+1:N+Nd] = dx1dx2_cov_exp_quad(xd,xd,Nd,Nd,alpha, length_scale);
    return Sigma;
  }
  // Gives latent f given all real + fake data. x contains the true values, xd contains the derivative values.
  vector latent_f(real rho, real alpha, real delta, vector eta, real[] x, real[] xd, int n, int nd) {
    matrix[n+nd, n+nd] K = cov_exp_quad_full(x, xd, n, nd, alpha, rho);
    matrix[n+nd, n+nd] L_K;
    for (j in 1:n+nd)
      K[j, j] = K[j, j] + delta;
    L_K = cholesky_decompose(K);
    return L_K * eta;
  }
}
data {
  int<lower=0> N1;
  real x1[N1];
  int<lower=0> N1d;
  real x1d[N1d];
  int<lower=0> N1df;
  real x1df[N1df];
  vector[N1] y1;
  vector[N1df] ydf;
  int<lower=0> N2;
  real x2[N2];
  real nu;
  real rho;
  real alpha;
  real sigma;
  real<lower=0> sigma_d_f;
}
transformed data {
  real delta = 1e-9;
  real x[N1+N2];
  real xd[N1d+N1df];
  int ones[N1d];
  for(i in 1:N1d) ones[i] = 1;
  x[1:N1] = x1;
  x[N1+1:N1+N2] = x2;
  xd[1:N1d] = x1d;
  xd[1+N1d:N1d+N1df] = x1df;
}
parameters {
  vector[N1+N2+N1d+N1df] eta;
}
transformed parameters {
  vector[N1+N2+N1d+N1df] f;
  f = latent_f(rho, alpha, delta, eta, x, xd, N1+N2, N1d+N1df);
}
model {
  eta ~ normal(0,1);
  y1 ~ normal(f[1:N1], sigma);
  ydf ~ normal(f[1+N1+N2+N1d:N1+N2+N1d+N1df], sigma_d_f);
  ones ~  bernoulli(inv_logit(nu * f[N1+N2+1:N1+N2+N1d]));
}
generated quantities {
  vector[N2] f2;
  vector[N2] y2;
  f2 = f[N1+1:N1+N2];
  for (n2 in 1:N2)
    y2[n2] = normal_rng(f2[n2], sigma);
}
