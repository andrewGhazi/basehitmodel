// This is adapted from reduce_redundancy_rs_nonsquare_prior_updates.stan but for one protein at a 
// time. That involves on the model side:
// 1. Making the protein baseline the intercept
// 
// And on the practical side:
// 1. Removing the reduce sum functionality.
//
// The other consequence is that the dispersion will be estimated on a per-protein basis.

data {
  int<lower = 1> n_strain;
  int<lower = 1> n_bc;
  int<lower = 1> n_ps;
  int<lower = 1> N; //total number of count observations
  
  int<lower=1, upper = n_strain> strain_i[N];
  
  int n_eta; 
  int<lower = 1, upper = n_strain> eta_ps_i[n_eta];
  vector[n_eta] eta_presum; // precomputed sum of depth & log(pre_count)
  int<lower = 1> eta_i[N]; // eta index onto count_obs
  
  int<lower = 0> count_obs[N];
}

parameters {
  real<lower = 0, upper = 1> theta[n_strain]; // the fraction of zeros by strain
  real protein_baseline;
  
  vector[n_ps] prot_strain;
  real<lower = 0> phi; 
}

model {
  int grainsize = 1;
  real bern_terms[n_strain,2];
  vector[n_eta] etas;
  etas = eta_presum + protein_baseline + prot_strain[eta_ps_i];
  
  protein_baseline ~ normal(-2,2);
  theta ~ beta(3,1);
  phi ~ gamma(2,2);
  prot_strain ~ normal(0,.15);
  
  // It would help to vectorize the contribution from the zero counts, but I 
  // can't figure out how to vectorize the log_sum_exp operation
  for (s in 1:n_strain){
    bern_terms[s,1] = bernoulli_lpmf(1 | theta[s]);
    bern_terms[s,2] = bernoulli_lpmf(0 | theta[s]);
  }
  
  for (i in 1:N){
    if (count_obs[i] == 0){
      target += log_sum_exp(bern_terms[strain_i[i],1], 
                            bern_terms[strain_i[i],2] + 
                              neg_binomial_2_log_lpmf(0 | etas[eta_i[i]], phi));
    } else {
      target += bern_terms[strain_i[i], 2] + 
                neg_binomial_2_log_lpmf(count_obs[i] | etas[eta_i[i]], phi);
    }
  }
}

