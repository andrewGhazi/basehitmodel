// This is adapted from reduce_redundancy_rs_nonsquare_prior_updates.stan but for one protein at a
// time. That involves on the model side:
// 1. Making the protein baseline the intercept
//
// And on the practical side:
// 1. Removing the reduce sum functionality.
//
// The other consequence is that the dispersion will be estimated on a per-protein basis.
//
// Some of the speedup tricks are discussed in this thread on the Stan discourse forum:
// https://discourse.mc-stan.org/t/speeding-up-a-hierarchical-negative-binomial-regression-with-a-huge-number-of-groups/19837
data {
  int<lower = 1> n_strain;
  int<lower = 1> n_bc;
  int<lower = 1> N; //total number of count observations

  array[N] int<lower=1, upper = n_strain> strain_i; // index array giving the strain for each count observation

  int n_eta; // number of unique eta components
  array[n_eta] int<lower = 1, upper = n_strain> eta_ps_i; // index array of etas onto protein:strain interactions. (It only goes up to n_strain because this model is fit to one protein at a time.)
  vector[n_eta] eta_presum; // precomputed sum of depth & log(pre_count)
  array[N] int<lower = 1> eta_i; // eta index onto count_obs

  real<lower=0> ixn_prior_width; // width of normal prior on interaction scores (default: 0.15)
  array[N] int<lower = 0> count_obs; // output count observations
}

parameters {
  array[n_strain] real<lower = 0, upper = 1> theta; // the fraction of zeros by strain
  real protein_baseline; // baseline output for the protein

  vector[n_strain] prot_strain; // vector of protein:strain ixn scores (length n_strain because the model is fit to only one protein)
  real<lower = 0> phi; // dispersion of counts in the negative binomial component
}

model {
  array[n_strain,2] real bern_terms;
  vector[n_eta] etas;
  etas = eta_presum + protein_baseline + prot_strain[eta_ps_i];

  protein_baseline ~ normal(-2,2);
  theta ~ beta(3,1);
  phi ~ gamma(2,2);
  prot_strain ~ normal(0, ixn_prior_width);

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

