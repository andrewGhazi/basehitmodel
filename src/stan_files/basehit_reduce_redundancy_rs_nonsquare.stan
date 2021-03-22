// This is basehit_reduce_redundancy_rs.stan but for nonsquare data: some interactions are cut out because they don't 
// have >= 3 nonzero counts.
functions {
  real partial_sum_lpmf(int[] slice_counts,
                         int start, int end,
                         int[] strain_i,
                         real[,] bern_terms,
                         int[] eta_i,
                         vector etas,
                         real phi) {

    real target_component = 0;
    int slice_i = 0;
    
    // Turns out it's currently impossible to vectorize mixture models:
    // https://mc-stan.org/docs/2_25/stan-users-guide/vectorizing-mixtures.html
    for (i in start:end){
      slice_i += 1;
      if (slice_counts[slice_i] == 0){
        target_component += log_sum_exp(bern_terms[strain_i[i],1], 
                                        bern_terms[strain_i[i],2] + 
                                          neg_binomial_2_log_lpmf(0 | etas[eta_i[i]],
                                                                      phi));
      } else{
        target_component += bern_terms[strain_i[i], 2] + neg_binomial_2_log_lpmf(slice_counts[slice_i] | etas[eta_i[i]], phi);
      }
      
    }

    return target_component;
  }
}

data {
  int<lower = 1> n_prot;
  int<lower = 1> n_strain;
  int<lower = 1> n_bc;
  int<lower = 1> n_ps;
  int<lower = 1> N; //total number of count observations
  
  int<lower=1, upper = n_strain> strain_i[N];
  
  int n_eta; // number of unique combinations of pre_count, protein, and protein:strain
  int<lower = 1, upper = n_prot> eta_p_i[n_eta];
  int<lower = 1, upper = n_prot*n_strain> eta_ps_i[n_eta];
  vector[n_eta] eta_presum; // precomputed sum of depth & log(pre_count)
  int<lower = 1> eta_i[N]; // eta index onto count_obs
  
  int<lower = 0> count_obs[N];
}

parameters {
  real<lower = 0, upper = 1> theta[n_strain]; // the fraction of zeros by strain
  
  real baseline_binding; // this should ideally be estimated from the counts in the streptavidin beads
  vector[n_prot] protein;
  vector[n_ps] prot_strain;
  real<lower = 0> phi; 
}

model {
  int grainsize = 1;
  real bern_terms[n_strain,2];
  vector[n_eta] etas;
  etas = eta_presum + baseline_binding + protein[eta_p_i] + prot_strain[eta_ps_i];
  
  baseline_binding ~ normal(-2,2);
  theta ~ beta(1,1);
  phi ~ gamma(2,2);
  protein ~ normal(0,.33);
  prot_strain ~ normal(0,.2);
  
  // It would help to vectorize the contribution from the zero counts, but I 
  // can't figure out how to vectorize the log_sum_exp operation
  for (s in 1:n_strain){
    bern_terms[s,1] = bernoulli_lpmf(1 | theta[s]);
    bern_terms[s,2] = bernoulli_lpmf(0 | theta[s]);
  }
  
  target += reduce_sum(partial_sum_lpmf, count_obs, grainsize,
                       strain_i, bern_terms, eta_i, etas, phi);
                       
}

