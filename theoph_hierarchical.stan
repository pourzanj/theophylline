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
  int<lower=0> N; //number of subjects
  int<lower = 0> Nt;          // number of measurement times
  real ts[N,Nt-1];                // measurement times > 0
  real y_init[N,2];            // initial measured populations
  real<lower = 0> y[N,Nt-1];   // measured populations
}
parameters {
  real<lower=0> b[N];          // bioavailablility  
  real<lower = 0> theta[N,3];   // { Ka,Cl,V }
  real<lower = 0> c_init[N];  // initial concentration in blood
  real<lower = 0> sigma[N];   // measurement error
}
transformed parameters {
  real<lower = 0> z_init[N,2];
  real z[N,Nt-1, 2];
  
  for(n in 1:N) {
    z_init[n,] = {y_init[n,1]*b[n], c_init[n] + 1e-6};
    z[n,,] = integrate_ode_rk45(dz_dt, z_init[n,], 0, ts[n,], theta[n,],
                         rep_array(0.0, 0), rep_array(0, 0),
                         1e-8, 1e-8, 5e2);
  }

  for(i in 1:N) {
    for(j in 1:(Nt-1)) {
      for(k in 1:2) {
        z[i,j,k] = z[i,j,k] + 1e-6;
      }
    }
  }
}
model {
  
  b ~ normal(0, 0.3);
  c_init ~ student_t(3, 0, 0.01);
  
  theta[,1] ~ normal(0, 10); // Ka
  theta[,2] ~ normal(0, 10); // Cl
  theta[,3] ~ normal(5.1, 0.2); // V usually 4.7 to 5.5 liters in the human body
  
  sigma ~ normal(0, 0.5);
  
  for(n in 1:N) {
    y_init[n,2] ~ lognormal(log(z_init[n,2]), sigma[n]);
    y[n,] ~ lognormal(log(z[n,,2]), sigma[n]);
  }
  
}
generated quantities {
  real y_init_rep[N];
  real y_rep[N,Nt-1];
  
  for(n in 1:N) {
    y_init_rep[n] = lognormal_rng(log(z_init[n,2]), sigma[n]);
    for (nt in 1:(Nt-1)) {
        y_rep[n,nt] = lognormal_rng(log(z[n,nt, 2]), sigma[n]);
    }
  }
}