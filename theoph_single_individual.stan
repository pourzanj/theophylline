functions {
  real[] dz_dt(real t,       // time
               real[] z,     // system state {prey, predator}
               real[] theta, // parameters
               real[] x_r,   // unused data
               int[] x_i) {
    real A = z[1];
    real c = z[2];

    real Ka = theta[1];
    real Cl = theta[2];
    real V = theta[3];

    real dA_dt = -Ka*A;
    real dc_dt = (1/V)*(Ka*A - (Cl/V)*(c*V));

    return { dA_dt, dc_dt };
  }
}
data {
  int<lower = 0> Nt;          // number of measurement times
  real ts[Nt-1];                // measurement times > 0
  real y_init[2];            // initial measured populations
  real<lower = 0> y[Nt-1];   // measured populations
}
parameters {
  real<lower=0> b;          // bioavailablility  
  real<lower = 0> theta[3];   // { Ka,Cl,V }
  real<lower = 0> c_init;  // initial concentration in blood
  real<lower = 0> sigma;   // measurement error
}
transformed parameters {
  real<lower = 0> z_init[2] = {y_init[1]*b, c_init};
  real z[Nt-1, 2]
    = integrate_ode_rk45(dz_dt, z_init, 0, ts, theta,
                         rep_array(0.0, 0), rep_array(0, 0),
                         1e-5, 1e-3, 5e2);
}
model {
  
  b ~ normal(0, 0.3);
  c_init ~ student_t(3, 0, 0.01);
  
  theta[1] ~ normal(0, 10); // Ka
  theta[2] ~ normal(0, 10); // Cl
  theta[3] ~ normal(5.1, 0.2); // V usually 4.7 to 5.5 liters in the human body
  
  sigma ~ normal(0, 0.5);
  
  y_init[2] ~ lognormal(log(z_init[2]), sigma);
  y ~ lognormal(log(z[, 2]), sigma);
  
}
generated quantities {
  real y_init_rep;
  real y_rep[Nt-1];
  
  y_init_rep = lognormal_rng(log(z_init[2]), sigma);
  for (n in 1:(Nt-1)) {
      y_rep[n] = lognormal_rng(log(z[n, 2]), sigma);
  }
}