data{
  int<lower=1> S; // Number of subjects
  int<lower=0> ns[S]; // counts of stimulated cells
  int<lower=0> nu[S]; // counts of unstimulated cells
  int<lower=0> Ns[S]; // counts of total cells in stimulated samples
  int<lower=0> Nu[S]; //counts of total cells in unstimulated samples
  real<lower=0, upper=1> a;
  real<lower=0,upper=1> b;
}

parameters{
  real<lower=0,upper=1> proportion; //estimated proportion of responders
  positive_ordered[2] mu;
  real<lower=1,upper=100000> phi_s;
  real<lower=1,upper=100000> phi_u;
}

transformed parameters{
  vector[2] log_ps[S];
  #beta_s = phi_s*(1-mu_s)
  #alpha_s = mu_s*phi_s
  #beta_u = phi_u*(1-mu_u)
  #alpha_u = mu_u*phi_u

  for(s in 1:S){
    log_ps[s,1] = bernoulli_lpmf(0|proportion) + beta_binomial_lpmf(ns[s]|Ns[s],mu[1]*phi_u,phi_u*(1-mu[1])) + beta_binomial_lpmf(nu[s]|Nu[s],mu[1]*phi_u,phi_u*(1-mu[1]));

      log_ps[s,2] = bernoulli_lpmf(1|proportion)  + beta_binomial_lpmf(ns[s]|Ns[s],mu[2]*phi_s,phi_s*(1-mu[2])) + beta_binomial_lpmf(nu[s]|Nu[s],mu[1]*phi_u,phi_u*(1-mu[1]));  
  }
}

model{
  mu ~ beta(1,1);
  proportion ~ beta(1,1);
  phi_u ~ exponential(a);
  phi_s ~ exponential(b);
  // nu ~ beta_binomial(Nu,mu[1]*phi_u,phi_u*(1-mu[1]));
	for(s in 1:S){
		target += log_sum_exp(log_ps[s]);
	}
}
generated quantities{
  real z[S];
  real<lower=0,upper=1> ps[S]; // proportion  of positive cells in stimulated sample
  real<lower=0,upper=1> pu[S]; //proportion of positive cells in unstimulated sample
	for(s in 1:S){
		vector[2] prob;
		prob = softmax(log_ps[s]);
		z[s] = bernoulli_rng(prob[2]);
		if(z[s]==1){
		    ps[s] = beta_rng(mu[2]*phi_s,phi_s*(1-mu[2]));
		}else{
		  	ps[s] = beta_rng(mu[1]*phi_u,phi_u*(1-mu[1]));
		}
		pu[s] = beta_rng(mu[1]*phi_u,phi_u*(1-mu[1]));
	}
}
 