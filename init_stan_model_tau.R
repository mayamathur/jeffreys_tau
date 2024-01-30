
# 2024-01-30 - PRIOR ON TAU ONLY ----------------------

model.text <- "

functions{

	real jeffreys_prior(real mu, real tau, int k, real[] sei){

		real kmm;
		real kms;
		real kss;
		matrix[2,2] fishinfo;
    real sigma;
		// will just be set to 1
		int n;

    real Si;


		// this will be the TOTALS for all observations
		matrix[2,2] fishinfototal;
		fishinfototal[1,1] = 0;
  	fishinfototal[1,2] = 0;
  	fishinfototal[2,1] = 0;
  	fishinfototal[2,2] = 0;


		// build a Fisher info matrix for EACH observation
		for (i in 1:k) {

		  Si = sqrt( tau^2 + sei[i]^2 );
      
      // only difference from joint Jeffreys prior is that here, kmm=1
      kmm = 1;
      kms = 0;
      kss = tau^2 * Si^(-4); 

  		fishinfo[1,1] = -kmm;
      fishinfo[1,2] = -kms;
      fishinfo[2,1] = -kms;
      fishinfo[2,2] = -kss;

  		// add the new fisher info to the total one
  		fishinfototal = fishinfototal + fishinfo;
		}

		return sqrt(determinant(fishinfototal));
	}
}

data{
	int<lower=0> k;
  real sei[k];
	real y[k];
}

parameters{
  real mu;
	real<lower=0> tau;
}


model{
  // this is to remove prior, as a sanity check:
  // target += 0;
  //see 'foundational ideas' here: https://vasishth.github.io/bayescogsci/book/sec-firststan.html
	target += log( jeffreys_prior(mu, tau, k, sei) );
	for(i in 1:k) {
      y[i] ~ normal( mu, sqrt(tau^2 + sei[i]^2) );
	}
}

// this chunk doesn't actually affect the model that's being fit to the data;
//  it's just re-calculating the prior, lkl, and post to return to user
// Stan docs: 'Nothing in the generated quantities block affects the sampled parameter values. The block is executed only after a sample has been generated'

generated quantities{
  real log_lik = 0;
  real log_prior = log( jeffreys_prior(mu, tau, k, sei) );
  real log_post;

  // versions that are evaluated at a SPECIFIC (mu=2, tau=2) so that we can compare
  //  to R functions for MAP, MLE, etc.
  real log_lik_sanity = 0;
  real log_prior_sanity = log( jeffreys_prior(2, 2, k, sei) );

  for ( i in 1:k ){
    log_lik += normal_lpdf( y[i] | mu, sqrt(tau^2 + sei[i]^2) );
    log_lik_sanity += normal_lpdf( y[i] | 2, sqrt(2^2 + sei[i]^2) );
  }
  log_post = log_prior + log_lik;
}
"



# necessary to prevent ReadRDS errors in which cores try to work with other cores' intermediate results
# https://groups.google.com/g/stan-users/c/8snqQTTfWVs?pli=1
options(mc.cores = parallel::detectCores())

# "isystem" arg is just a placeholder to avoid Stan's not understanding special characters
#  in getwd(), even though we don't actually use the dir at all
# note: removing the isystem arg does NOT fix the very sporadic "error reading from connection" on cluster
cat( paste("\n init_stan_model_tau: about to call stan_model_tau") )
stan.model.tau <- stan_model(model_code = model.text,
                         isystem = "~/Desktop")


cat( paste("\n init_stan_model_tau: done calling stan_model_tau") )