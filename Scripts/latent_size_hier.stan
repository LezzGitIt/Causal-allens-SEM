// Hierarchical latent-size SEM across species and studies (STANDARDIZED factor).
//
// The latent Size factor has variance fixed to 1 (there is no free sigma_size to
// collapse toward zero — that collapse, with loadings running to infinity, is the
// funnel that wrecked the unstandardized version on weakly-coupled species).
//
// Per species s, on data z-scored WITHIN species:
//   Size_i = kappa[s]*T_i + sqrt(1 - kappa[s]^2) * eta_i,  eta_i ~ N(0,1)   // Var(Size)=1; kappa = cor(Size,T) = Bergmann
//   y_k,i  = a_k[s]*Size_i + d_k[s]*T_i + e_k,i,  e ~ N(0, psi_k[s])        // d_mass = 0 (mass anchors, no direct effect)
// with y = (mass, wing, app2). Marginalizing Size, each individual is
//   y_i ~ MVN(c_s * T_i, Sigma_s),
//   c_s    = ( a_mass*kappa,  a_wing*kappa + d_wing,  a_app2*kappa + d_app2 ),
//   Sigma_s = (1 - kappa^2) * a a' + diag(psi),   a = (a_mass, a_wing, a_app2).
// The likelihood uses per-species SUFFICIENT STATISTICS (n, sum T^2, sum T*y, sum y y').
//
// a_* are STANDARDIZED loadings (correlations of each indicator with Size), so weak
// allometric coupling is just a_k ~ 0 — bounded and harmless, no runaway. Sign is
// fixed by a_mass > 0. Estimand: d_wing[s] per species, mu_dw[study] per study.
//
// Pooling: wing loading a_wing GLOBAL; 2nd-appendage loading, Bergmann corr kappa,
// and direct effects pooled WITHIN study (keeps per-year vs per-degC effects apart).

data {
  int<lower=1> S;
  int<lower=1> G;
  array[S] int<lower=1, upper=G> study;
  array[S] int<lower=1> n;
  array[S] real<lower=0> Sxx;
  array[S] vector[3] Syx;               // order: mass, wing, app2
  array[S] matrix[3, 3] Syy;
}

parameters {
  // standardized wing loading — pooled GLOBALLY (non-centred)
  real mu_aw;  real<lower=0> tau_aw;  vector[S] z_aw;
  // standardized 2nd-appendage loading — pooled BY STUDY
  vector[G] mu_a2; real<lower=0> tau_a2; vector[S] z_a2;
  // direct effects — pooled BY STUDY (non-centred)
  vector[G] mu_dw; real<lower=0> tau_dw; vector[S] z_dw;   // wing (ESTIMAND)
  vector[G] mu_d2; real<lower=0> tau_d2; vector[S] z_d2;   // 2nd appendage (nuisance)
  // mass loading (sign anchor) and Bergmann correlation, both bounded
  vector<lower=0, upper=1>[S] a_mass;
  vector<lower=-1, upper=1>[S] kappa;
  vector[G] mu_kappa; real<lower=0> tau_kappa;
  // residual variances
  vector<lower=0>[S] psi_mass;
  vector<lower=0>[S] psi_wing;
  vector<lower=0>[S] psi_app2;
}

transformed parameters {
  vector[S] a_wing = mu_aw         + tau_aw * z_aw;
  vector[S] a_app2 = mu_a2[study]  + tau_a2 * z_a2;
  vector[S] d_wing = mu_dw[study]  + tau_dw * z_dw;   // per-species direct wing effect
  vector[S] d_app2 = mu_d2[study]  + tau_d2 * z_d2;
}

model {
  mu_aw ~ normal(0.4, 0.3);   tau_aw ~ normal(0, 0.3);
  mu_a2 ~ normal(0.4, 0.3);   tau_a2 ~ normal(0, 0.3);
  mu_dw ~ normal(0, 0.5);     tau_dw ~ normal(0, 0.3);
  mu_d2 ~ normal(0, 0.5);     tau_d2 ~ normal(0, 0.3);
  z_aw ~ std_normal(); z_a2 ~ std_normal(); z_dw ~ std_normal(); z_d2 ~ std_normal();
  a_mass ~ normal(0.6, 0.25);                       // truncated to (0, 1)
  kappa  ~ normal(mu_kappa[study], tau_kappa);      // truncated to (-1, 1)
  mu_kappa ~ normal(0, 0.5); tau_kappa ~ normal(0, 0.3);
  psi_mass ~ normal(0, 0.7); psi_wing ~ normal(0, 0.7); psi_app2 ~ normal(0, 0.7);

  for (s in 1:S) {
    real k = kappa[s];
    vector[3] a = [a_mass[s], a_wing[s], a_app2[s]]';
    vector[3] c = [a_mass[s] * k,
                   a_wing[s] * k + d_wing[s],
                   a_app2[s] * k + d_app2[s]]';
    matrix[3, 3] Sig = (1 - square(k)) * (a * a')
                       + diag_matrix([psi_mass[s], psi_wing[s], psi_app2[s]]');
    matrix[3, 3] Sinv = inverse_spd(Sig);
    real logdet = log_determinant(Sig);
    real quad = trace(Sinv * Syy[s])
                - 2 * dot_product(c, Sinv * Syx[s])
                + Sxx[s] * quad_form(Sinv, c);
    target += -0.5 * (n[s] * (3 * log(2 * pi()) + logdet) + quad);
  }
}
