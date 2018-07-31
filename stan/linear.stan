// Fit the hyperparameters of a Gaussian process with an exponentiated
// quadratic kernel and analytically predict out-of-sample observations

data {
  int<lower=1> N1;
  real x1[N1];
  vector[N1] y1;
  int<lower=1> N2;
  real x2[N2];
}
parameters {
  real b;
  real k;
  real sigma;
}
transformed parameters {
  vector[N1] f1;
  for (n1 in 1:N1)
    f1[n1] = b+k*x1[n1];
}
model {
  y1 ~ normal(f1,sigma);
}
generated quantities {
  vector[N2] f2;
  vector[N2] y2;
  for (n2 in 1:N2)
    f2[n2] = b+k*x2[n2];
  for (n2 in 1:N2)
    y2[n2] = normal_rng(f2[n2], sigma);
}
