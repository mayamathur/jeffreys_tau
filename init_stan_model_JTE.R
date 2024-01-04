
# version that assumes only nonaffirmative results:
# https://github.com/mayamathur/phacking/blob/e335a054aedc448d9a1bc0f6e82300a0fb876583/Sherlock%20code/init_stan_model_SAPH.R

# # to test if this compiles in interactive session:
# # helper code
# path = "/home/groups/manishad/SAPH"
# setwd(path)
# source("init_stan_model_applied_SAPH.R")


# 2022-5-3: SIMPLIFIED JEFFREYS PRIOR WRT TAU ----------------------


model.text <- "

functions{

	real jeffreys_prior(real mu, real tau, int k, real[] sei, real[] tcrit, real[] affirm){

		real mustarL;
		real mustarU;
		real alphaL;
		real alphaU;
		real kmm;
		real kms;
		real kss;
		matrix[2,2] fishinfo;
    real sigma;
		real LL;
		real UU;
		// will just be set to 1
		int n;

    real Si;
    real cz;
    real dnor;
    real pnor;
    real r;

		// this will be the TOTALS for all observations
		matrix[2,2] fishinfototal;
		fishinfototal[1,1] = 0;
  	fishinfototal[1,2] = 0;
  	fishinfototal[2,1] = 0;
  	fishinfototal[2,2] = 0;


		// build a Fisher info matrix for EACH observation
		for (i in 1:k) {

		  Si = sqrt( tau^2 + sei[i]^2 );
      cz = (sei[i] * tcrit[i] - mu) / Si;
      dnor = exp( normal_lpdf(cz | 0, 1 ) );
      pnor = exp( normal_lcdf(cz | 0, 1 ) );
      r = dnor/pnor;
      
      kmm = Si^(-2)*(cz*r + r^2 - 1);
      kms = tau*Si^(-3)*r*( cz^2 + cz*r + 1 );
      kss = ( tau^2 * Si^(-4) ) * ( cz^3*r + cz^2*r^2 + cz*r - 2 ); 

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
  real tcrit[k];
  real affirm[k];
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
	target += log( jeffreys_prior(mu, tau, k, sei, tcrit, affirm) );
	for(i in 1:k) {
      if ( affirm[i] == 0 ) {
        y[i] ~ normal( mu, sqrt(tau^2 + sei[i]^2) ) T[ , tcrit[i] * sei[i] ];
      } else if ( affirm[i] == 1 ) {
        y[i] ~ normal( mu, sqrt(tau^2 + sei[i]^2) ) T[ tcrit[i] * sei[i] , ];
      }
	}
}

// this chunk doesn't actually affect the model that's being fit to the data;
//  it's just re-calculating the prior, lkl, and post to return to user
// Stan docs: 'Nothing in the generated quantities block affects the sampled parameter values. The block is executed only after a sample has been generated'

generated quantities{
  real log_lik = 0;
  real log_prior = log( jeffreys_prior(mu, tau, k, sei, tcrit, affirm) );
  real log_post;
  // this is just an intermediate quantity for log_lik
  // will be equal to UU or LL above, depending on affirm status
  real critScaled;

  // versions that are evaluated at a SPECIFIC (mu=2, tau=2) so that we can compare
  //  to R functions for MAP, MLE, etc.
  real log_lik_sanity = 0;
  real log_prior_sanity = log( jeffreys_prior(2, 2, k, sei, tcrit, affirm) );

  for ( i in 1:k ){
    log_lik += normal_lpdf( y[i] | mu, sqrt(tau^2 + sei[i]^2) );
    log_lik_sanity += normal_lpdf( y[i] | 2, sqrt(2^2 + sei[i]^2) );

    critScaled = tcrit[i] * sei[i];

    // https://mc-stan.org/docs/2_20/reference-manual/sampling-statements-section.html
    // see 'Truncation with upper bounds in Stan' section
    // nonaffirm case:
    if ( y[i] <= critScaled ) {
    // from sanity checks in doParallel, I know this matches joint_nll_2
      log_lik += -1 * normal_lcdf(critScaled | mu, sqrt(tau^2 + sei[i]^2) );
      log_lik_sanity += -1 * normal_lcdf(critScaled | 2, sqrt(2^2 + sei[i]^2) );

    // affirm case:
    } else if ( y[i] > critScaled ) {
      log_lik += -1 * log( 1 - normal_cdf( critScaled, mu, sqrt(tau^2 + sei[i]^2) ) );
      log_lik_sanity += -1 * log( 1 - normal_cdf( critScaled, 2, sqrt(2^2 + sei[i]^2) ) );
    }
  }
  log_post = log_prior + log_lik;
}
"


# UNTRUNCATED PRIOR ------------------------------

# # see prior_untrunc in “2022-4-20 New Jeffreys prior using numerical derivatives wrt tau.R” 
# model.text <- "
# 
# functions{
# 
# 	real jeffreys_prior(real mu, real tau, int k, real[] sei, real[] tcrit, real[] affirm){
# 
# 		real mustarL;
# 		real mustarU;
# 		real alphaL;
# 		real alphaU;
# 		real kmm;
# 		real kms;
# 		real kss;
# 		matrix[2,2] fishinfo;
#     real sigma;
# 		real LL;
# 		real UU;
# 		// will just be set to 1
# 		int n;
# 
#   real e1;
# 		real e2;
#     real e3;
#     real e4;
#     real e5;
#     real e6;
#     real e7;
#     real e8;
#     real e9;
#     real e10;
#     real e11;
#     real e12;
#     real e13;
#     real e14;
#     real e15;
#     real expectation1;
#     real expectation2;
#     real fishdet;
# 
# 		// this will be the TOTALS for all observations
# 		matrix[2,2] fishinfototal;
# 		fishinfototal[1,1] = 0;
#   	fishinfototal[1,2] = 0;
#   	fishinfototal[2,1] = 0;
#   	fishinfototal[2,2] = 0;
# 
# 
# 		// build a Fisher info matrix for EACH observation
#   for (i in 1:k) {
#     
#     // for this observation
#  		// this will be the TOTALS for all observations
#  		fishinfo[1,1] = 0;
#    	fishinfo[1,2] = 0;
#    	fishinfo[2,1] = 0;
#    	fishinfo[2,2] = 0;
#     
#     kmm = -(1/(sei[i]^2 + tau^2));
#     
#     // this is just 0
#     kms = -(2 * (tau * (0)/(sei[i]^2 + tau^2)^2));
#     
#     // from body of R's get_D22_num:
#     e1 = tau^2;
#     e3 = sei[i]^2 + e1;
#     // expectation of (.yi - mu)^2 is just marginal variance when not truncated
#     //e5 = (.yi - mu)^2
#     e5 = tau^2 + sei[i]^2;
#     kss = ((e5 * (1 - 2 * (e1/e3)) - 2 * (e1 * (e5/e3 - 1)))/e3 - 
#         1)/e3;
# 
# 
#     fishinfo[1,1] = -kmm;
#     fishinfo[1,2] = -kms;
#     fishinfo[2,1] = -kms;
#     fishinfo[2,2] = -kss;
#     
#     // add the new fisher info to the total one
#     fishinfototal = fishinfototal + fishinfo;
#   }
#   
# 	
# 		// can't use 0 here because we'll be taking log of prior
# 		fishdet = determinant(fishinfototal);
# 		//if ( fishdet < 0 ) fishdet = 0.01;
# 		return sqrt(fishdet);
# 	}
# }
# 
# data{
# 	int<lower=0> k;
#   real sei[k];
#   real tcrit[k];
#   real affirm[k];
# 	real y[k];
# }
# 
# parameters{
#   real mu;
# 	real<lower=0> tau;
# }
# 
# 
# model{
#   // this is to remove prior, as a sanity check:
#   // target += 0;
#   //see 'foundational ideas' here: https://vasishth.github.io/bayescogsci/book/sec-firststan.html
# 	target += log( jeffreys_prior(mu, tau, k, sei, tcrit, affirm) );
# 	for(i in 1:k) {
#       if ( affirm[i] == 0 ) {
#         y[i] ~ normal( mu, sqrt(tau^2 + sei[i]^2) ) T[ , tcrit[i] * sei[i] ];
#       } else if ( affirm[i] == 1 ) {
#         y[i] ~ normal( mu, sqrt(tau^2 + sei[i]^2) ) T[ tcrit[i] * sei[i] , ];
#       }
# 	}
# }
# 
# // this chunk doesn't actually affect the model that's being fit to the data;
# //  it's just re-calculating the prior, lkl, and post to return to user
# // Stan docs: 'Nothing in the generated quantities block affects the sampled parameter values. The block is executed only after a sample has been generated'
# 
# generated quantities{
#   real log_lik = 0;
#   real log_prior = log( jeffreys_prior(mu, tau, k, sei, tcrit, affirm) );
#   real log_post;
#   // this is just an intermediate quantity for log_lik
#   // will be equal to UU or LL above, depending on affirm status
#   real critScaled;
# 
#   // versions that are evaluated at a SPECIFIC (mu=2, tau=2) so that we can compare
#   //  to R functions for MAP, MLE, etc.
#   real log_lik_sanity = 0;
#   real log_prior_sanity = log( jeffreys_prior(2, 2, k, sei, tcrit, affirm) );
# 
#   for ( i in 1:k ){
#     log_lik += normal_lpdf( y[i] | mu, sqrt(tau^2 + sei[i]^2) );
#     log_lik_sanity += normal_lpdf( y[i] | 2, sqrt(2^2 + sei[i]^2) );
# 
#     critScaled = tcrit[i] * sei[i];
# 
#     // https://mc-stan.org/docs/2_20/reference-manual/sampling-statements-section.html
#     // see 'Truncation with upper bounds in Stan' section
#     // nonaffirm case:
#     if ( y[i] <= critScaled ) {
#     // from sanity checks in doParallel, I know this matches joint_nll_2
#       log_lik += -1 * normal_lcdf(critScaled | mu, sqrt(tau^2 + sei[i]^2) );
#       log_lik_sanity += -1 * normal_lcdf(critScaled | 2, sqrt(2^2 + sei[i]^2) );
# 
#     // affirm case:
#     } else if ( y[i] > critScaled ) {
#       log_lik += -1 * log( 1 - normal_cdf( critScaled, mu, sqrt(tau^2 + sei[i]^2) ) );
#       log_lik_sanity += -1 * log( 1 - normal_cdf( critScaled, 2, sqrt(2^2 + sei[i]^2) ) );
#     }
#   }
#   log_post = log_prior + log_lik;
# }
# "

# NEW JEFFREYS PRIOR W/ CHANGE OF VARS FROM R'S DERIV -----------------------------------
# See 2022-4-20 Check Mathematica using R

# model.text <- "
# 
# functions{
# 
# 	real jeffreys_prior(real mu, real tau, int k, real[] sei, real[] tcrit, real[] affirm){
# 
# 		real mustarL;
# 		real mustarU;
# 		real alphaL;
# 		real alphaU;
# 		real kmm;
# 		real kms;
# 		real kss;
# 		matrix[2,2] fishinfo;
#     real sigma;
# 		real LL;
# 		real UU;
# 		// will just be set to 1
# 		int n;
# 
#   real e1;
# 		real e2;
#     real e3;
#     real e4;
#     real e5;
#     real e6;
#     real e7;
#     real e8;
#     real e9;
#     real e10;
#     real e11;
#     real e12;
#     real e13;
#     real e14;
#     real e15;
#     real expectation1;
#     real expectation2;
# 
# 		// this will be the TOTALS for all observations
# 		matrix[2,2] fishinfototal;
# 		fishinfototal[1,1] = 0;
#   	fishinfototal[1,2] = 0;
#   	fishinfototal[2,1] = 0;
#   	fishinfototal[2,2] = 0;
# 
# 
# 		// build a Fisher info matrix for EACH observation
# 		for (i in 1:k) {
# 
# 		  // from body of R's get_D11_num:
#        e2 = sei[i]^2 + tau^2;
#         e3 = sqrt(e2);
#         e5 = sei[i] * tcrit[i] - mu;
#         e6 = e5/e3;
#         // R version:
#         //e7 = dnorm(e6, 0, 1);
#         e7 = exp( normal_lpdf(e6 | 0, 1) );
#         // R version:
#         //e8 = pnorm(e6);
#         e8 = exp( normal_lcdf(e6 | 0, 1 ) );
# 		    kmm = -(1/e2 - (e5/(e2 * e8) + e7 * e3/(e8 * e3)^2) * e7/e3);
# 
# 		  // from body of R's get_D12_num:
# 		    e2 = sei[i]^2 + tau^2;
#         e3 = sqrt(e2);
#         e5 = sei[i] * tcrit[i] - mu;
#         // e6 is scaled critical value:
#         e6 = e5/e3;
#         //e7 = pnorm(e6);
#         e7 = exp( normal_lcdf(e6 | 0, 1 ) );
#         e8 = e2^2;
#         //e9 = dnorm(e6, 0, 1);
#         e9 = exp( normal_lpdf(e6 | 0, 1) );
# 
#         // my own expectation of .yi - .mu:
#         expectation1 = -sqrt(sei[i]^2 + tau^2) * e9/e7;
#         kms = -(tau * (((e7/e3 - e5 * e9/e2)/(e7 * e3)^2 - e5^2/(e8 *
#             e7 * e3)) * e9 + 2 * ((expectation1)/e8)));
# 
# 
#   		// from body of R's get_D22_num:
#   		  e1 = tau^2;
#         e3 = sei[i]^2 + e1;
#         e5 = sei[i] * tcrit[i] - mu;
#         e6 = sqrt(e3);
#         // e7 is scaled crit value:
#         e7 = e5/e6;
#         //e8 = pnorm(e7)
#         e8 = exp( normal_lcdf(e7 | 0, 1 ) );
#         //e9 = dnorm(e7, 0, 1)
#         e9 = exp( normal_lpdf(e7 | 0, 1 ) );
#         e10 = e5 * e9;
#         e11 = e8 * e6;
#         e13 = e10/e11;
#         // *replace this one with its expectation:
#         //e15 = (.yi - .mu)^2/e3
#         // expectation of (.yi - .mu)^2:
#         expectation2 = (sei[i]^2 + tau^2)*(1 - e7 * e9/e8);
#         e15 = expectation2/e3;
# 
#         kss = (e13 + e15 - (e1 * (e5 * ((e8/e6 - e10/e3)/e11^2 -
#             e5^2/(e3^2 * e8 * e6)) * e9 + 2 * ((e13 + 2 * e15 -
#             1)/e3)) + 1))/e3;
# 
# 
#   		fishinfo[1,1] = -kmm;
#       fishinfo[1,2] = -kms;
#       fishinfo[2,1] = -kms;
#       fishinfo[2,2] = -kss;
# 
#   		// add the new fisher info to the total one
#   		fishinfototal = fishinfototal + fishinfo;
# 		}
# 		
# 		
# 		//TEMP ONLY: NUMERICAL ISSUES
# 		if ( sqrt(determinant(fishinfototal)) < 0 ) {
# 		return 0;
# 		} else {
# 		  return( sqrt(determinant(fishinfototal)) );
# 		}
# 		//return sqrt(determinant(fishinfototal));
# 	}
# }
# 
# data{
# 	int<lower=0> k;
#   real sei[k];
#   real tcrit[k];
#   real affirm[k];
# 	real y[k];
# }
# 
# parameters{
#   real mu;
# 	real<lower=0> tau;
# }
# 
# 
# model{
#   // this is to remove prior, as a sanity check:
#   // target += 0;
#   //see 'foundational ideas' here: https://vasishth.github.io/bayescogsci/book/sec-firststan.html
# 	target += log( jeffreys_prior(mu, tau, k, sei, tcrit, affirm) );
# 	for(i in 1:k) {
#       if ( affirm[i] == 0 ) {
#         y[i] ~ normal( mu, sqrt(tau^2 + sei[i]^2) ) T[ , tcrit[i] * sei[i] ];
#       } else if ( affirm[i] == 1 ) {
#         y[i] ~ normal( mu, sqrt(tau^2 + sei[i]^2) ) T[ tcrit[i] * sei[i] , ];
#       }
# 	}
# }
# 
# // this chunk doesn't actually affect the model that's being fit to the data;
# //  it's just re-calculating the prior, lkl, and post to return to user
# // Stan docs: 'Nothing in the generated quantities block affects the sampled parameter values. The block is executed only after a sample has been generated'
# 
# generated quantities{
#   real log_lik = 0;
#   real log_prior = log( jeffreys_prior(mu, tau, k, sei, tcrit, affirm) );
#   real log_post;
#   // this is just an intermediate quantity for log_lik
#   // will be equal to UU or LL above, depending on affirm status
#   real critScaled;
# 
#   // versions that are evaluated at a SPECIFIC (mu=2, tau=2) so that we can compare
#   //  to R functions for MAP, MLE, etc.
#   real log_lik_sanity = 0;
#   real log_prior_sanity = log( jeffreys_prior(2, 2, k, sei, tcrit, affirm) );
# 
#   for ( i in 1:k ){
#     log_lik += normal_lpdf( y[i] | mu, sqrt(tau^2 + sei[i]^2) );
#     log_lik_sanity += normal_lpdf( y[i] | 2, sqrt(2^2 + sei[i]^2) );
# 
#     critScaled = tcrit[i] * sei[i];
# 
#     // https://mc-stan.org/docs/2_20/reference-manual/sampling-statements-section.html
#     // see 'Truncation with upper bounds in Stan' section
#     // nonaffirm case:
#     if ( y[i] <= critScaled ) {
#     // from sanity checks in doParallel, I know this matches joint_nll_2
#       log_lik += -1 * normal_lcdf(critScaled | mu, sqrt(tau^2 + sei[i]^2) );
#       log_lik_sanity += -1 * normal_lcdf(critScaled | 2, sqrt(2^2 + sei[i]^2) );
# 
#     // affirm case:
#     } else if ( y[i] > critScaled ) {
#       log_lik += -1 * log( 1 - normal_cdf( critScaled, mu, sqrt(tau^2 + sei[i]^2) ) );
#       log_lik_sanity += -1 * log( 1 - normal_cdf( critScaled, 2, sqrt(2^2 + sei[i]^2) ) );
#     }
#   }
#   log_post = log_prior + log_lik;
# }
# "



# NEW JEFFREYS PRIOR W/ CHANGE OF VARIABLES - FROM MANUAL THEORY ON IPAD  -----------------------------------
# 2022-4-19
# # LL and UU: cutpoints on RAW scale, not Z-scores
# # tau: SD, not variance
# # READY TO TRY THIS :)
# model.text <- "
# 
# functions{
# 
# 	real jeffreys_prior(real mu, real tau, int k, real[] sei, real[] tcrit, real[] affirm){
# 
# 	  // these will be overwritten for EACH observation
# 		real mustarL;
# 		real mustarU;
# 		real alphaL;
# 		real alphaU;
# 		real kmm;
# 		real kms;
# 		real kss;
# 		matrix[2,2] fishinfo;
#     real sigma;
# 		real LL;
# 		real UU;
# 		// will just be set to 1
# 		int n;
# 		
# 		real dS_dtau;
# 		real d2S_dtau2;
# 		real expectation2;
# 
# 
# 		// this will be the TOTALS for all observations
# 		matrix[2,2] fishinfototal;
# 		fishinfototal[1,1] = 0;
#   	fishinfototal[1,2] = 0;
#   	fishinfototal[2,1] = 0;
#   	fishinfototal[2,2] = 0;
# 
# 
# 		// build a Fisher info matrix for EACH observation
# 		for (i in 1:k) {
# 
# 		  // MARGINAL SD for this one observation
# 		  sigma = sqrt(tau^2 + sei[i]^2);
# 
# 		  // depending on whether study is affirmative, set truncation limits
# 		  // for THIS study, given its SE
# 		  if ( affirm[i] == 0 ) {
# 		  		UU = tcrit[i] * sei[i];
# 		  		// standardized truncation limits
# 		  		mustarL = -999;
#   		    mustarU = (UU - mu) / sigma;
# 		  } else if ( affirm[i] == 1 ) {
# 		      LL = tcrit[i] * sei[i];
# 		      // standardized truncation limits
# 		  		mustarL = (LL - mu) / sigma;
#   		    mustarU = 999;
# 		  }
# 		  
#   		// beginning of stuff that is not modified at all from TNE,
#   		//  *except* for the change-of-variables terms in kms and kss applied in a
#   		//   final code block
#   		// note that normal_lpdf, etc., parameterize in terms of SD, not var
#   		//  the (0,1) below are *not* start values for MCMC
#   		alphaU = exp( normal_lpdf(mustarU | 0, 1) -
#    	                log_diff_exp( normal_lcdf(mustarU | 0, 1),
#    	                normal_lcdf(mustarL | 0, 1) ) );
#    	        
#    	  // for change of vars from sigma to tau:        
#    	  // first derivative (blue)
#   		dS_dtau = tau/sigma;
#   		// second derivative (purple)
#   		d2S_dtau2 = sigma^(-1) * ( 1 - tau^2 * sigma^(-2) );
# 
#   		// second derivatives for Fisher info
#   		// wrt tau (not sigma) via change-of-variables terms
#   		
#   		//2022-4-19: updated :)
#   		kmm = -1/sigma^2 + 1/sigma^2 * (alphaU^2 + alphaU*mustarU);
#   		
#   		kms = ( 2*sigma^(-3)*alphaU + sigma^(-3)*( alphaU*mustarU^2 + alphaU^2*mustarU - alphaU ) ) * tau;
#   		
#   		expectation2 = sigma^2*(1 - alphaU*mustarU);
#       kss = ( -sigma^(-1) + sigma^(-3)*expectation2 + sigma^(-1)*mustarU*alphaU ) * d2S_dtau2^2 +
#             ( sigma^(-2) - 3*sigma^(-4)*expectation2 + sigma^(-2)*( mustarU*alphaU*(mustarU^2 - 2) + (mustarU*alphaU)^2 ) ) * dS_dtau^2;
#   			  				
#   		fishinfo[1,1] = -kmm;
#       fishinfo[1,2] = -kms;
#       fishinfo[2,1] = -kms;
#       fishinfo[2,2] = -kss;
# 
#   		// MM: add the new fisher info to the total one
#   		fishinfototal = fishinfototal + fishinfo;
# 		}
# 		return sqrt(determinant(fishinfototal));
# 	}
# }
# 
# data{
# 	int<lower=0> k;
#   real sei[k];
#   real tcrit[k];
#   real affirm[k];
# 	real y[k];
# }
# 
# parameters{
#   real mu;
# 	real<lower=0> tau;
# }
# 
# 
# model{
#   // this is to remove prior, as a sanity check:
#   // target += 0;
#   //see 'foundational ideas' here: https://vasishth.github.io/bayescogsci/book/sec-firststan.html
# 	target += log( jeffreys_prior(mu, tau, k, sei, tcrit, affirm) );
# 	for(i in 1:k) {
#       if ( affirm[i] == 0 ) {
#         y[i] ~ normal( mu, sqrt(tau^2 + sei[i]^2) ) T[ , tcrit[i] * sei[i] ];
#       } else if ( affirm[i] == 1 ) {
#         y[i] ~ normal( mu, sqrt(tau^2 + sei[i]^2) ) T[ tcrit[i] * sei[i] , ];
#       }
# 	}
# }
# 
# // this chunk doesn't actually affect the model that's being fit to the data;
# //  it's just re-calculating the prior, lkl, and post to return to user
# // Stan docs: 'Nothing in the generated quantities block affects the sampled parameter values. The block is executed only after a sample has been generated'
# 
# generated quantities{
#   real log_lik = 0;
#   real log_prior = log( jeffreys_prior(mu, tau, k, sei, tcrit, affirm) );
#   real log_post;
#   // this is just an intermediate quantity for log_lik
#   // will be equal to UU or LL above, depending on affirm status
#   real critScaled;
# 
#   // versions that are evaluated at a SPECIFIC (mu=2, tau=2) so that we can compare
#   //  to R functions for MAP, MLE, etc.
#   real log_lik_sanity = 0;
#   real log_prior_sanity = log( jeffreys_prior(2, 2, k, sei, tcrit, affirm) );
# 
#   for ( i in 1:k ){
#     log_lik += normal_lpdf( y[i] | mu, sqrt(tau^2 + sei[i]^2) );
#     log_lik_sanity += normal_lpdf( y[i] | 2, sqrt(2^2 + sei[i]^2) );
# 
#     critScaled = tcrit[i] * sei[i];
# 
#     // https://mc-stan.org/docs/2_20/reference-manual/sampling-statements-section.html
#     // see 'Truncation with upper bounds in Stan' section
#     // nonaffirm case:
#     if ( y[i] <= critScaled ) {
#     // from sanity checks in doParallel, I know this matches joint_nll_2
#       log_lik += -1 * normal_lcdf(critScaled | mu, sqrt(tau^2 + sei[i]^2) );
#       log_lik_sanity += -1 * normal_lcdf(critScaled | 2, sqrt(2^2 + sei[i]^2) );
# 
#     // affirm case:
#     } else if ( y[i] > critScaled ) {
#       log_lik += -1 * log( 1 - normal_cdf( critScaled, mu, sqrt(tau^2 + sei[i]^2) ) );
#       log_lik_sanity += -1 * log( 1 - normal_cdf( critScaled, 2, sqrt(2^2 + sei[i]^2) ) );
#     }
#   }
#   log_post = log_prior + log_lik;
# }
# "



# JUST 1/TAU^2  -----------------------------------
# # great when no heterogeneity; otherwise very bad
# model.text <- "
# 
# functions{
# 
# 	real simple_prior(real mu, real tau, int k, real[] sei, real[] tcrit, real[] affirm){
# 		return 1 / ( tau^2 );
# 	}
# }
# 
# data{
# 	int<lower=0> k;
#   real sei[k];
#   real tcrit[k];
#   real affirm[k];
# 	real y[k];
# }
# 
# parameters{
#   real mu;
# 	real<lower=0> tau;
# }
# 
# 
# model{
#   // this is to remove prior, as a sanity check:
#   // target += 0;
#   //see 'foundational ideas' here: https://vasishth.github.io/bayescogsci/book/sec-firststan.html
# 	target += log( simple_prior(mu, tau, k, sei, tcrit, affirm) );
# 	for(i in 1:k) {
#       if ( affirm[i] == 0 ) {
#         y[i] ~ normal( mu, sqrt(tau^2 + sei[i]^2) ) T[ , tcrit[i] * sei[i] ];
#       } else if ( affirm[i] == 1 ) {
#         y[i] ~ normal( mu, sqrt(tau^2 + sei[i]^2) ) T[ tcrit[i] * sei[i] , ];
#       }
# 	}
# }
# 
# // this chunk doesn't actually affect the model that's being fit to the data;
# //  it's just re-calculating the prior, lkl, and post to return to user
# // Stan docs: 'Nothing in the generated quantities block affects the sampled parameter values. The block is executed only after a sample has been generated'
# 
# generated quantities{
#   real log_lik = 0;
#   real log_prior = log( simple_prior(mu, tau, k, sei, tcrit, affirm) );
#   real log_post;
#   // this is just an intermediate quantity for log_lik
#   // will be equal to UU or LL above, depending on affirm status
#   real critScaled;
# 
#    // versions that are evaluated at a SPECIFIC (mu=2, tau=2) so that we can compare
#    //  to R functions for MAP, MLE, etc.
#    real log_lik_sanity = 0;
#    real log_prior_sanity = log( simple_prior(2, 2, k, sei, tcrit, affirm) );
# 
#   for ( i in 1:k ){
#     log_lik += normal_lpdf( y[i] | mu, sqrt(tau^2 + sei[i]^2) );
#     log_lik_sanity += normal_lpdf( y[i] | 2, sqrt(2^2 + sei[i]^2) );
# 
#     critScaled = tcrit[i] * sei[i];
# 
#     // https://mc-stan.org/docs/2_20/reference-manual/sampling-statements-section.html
#     // see 'Truncation with upper bounds in Stan' section
#     // nonaffirm case:
#     if ( y[i] <= critScaled ) {
#     // from sanity checks in doParallel, I know this matches joint_nll_2
#       log_lik += -1 * normal_lcdf(critScaled | mu, sqrt(tau^2 + sei[i]^2) );
#       log_lik_sanity += -1 * normal_lcdf(critScaled | 2, sqrt(2^2 + sei[i]^2) );
# 
#     // affirm case:
#     } else if ( y[i] > critScaled ) {
#       log_lik += -1 * log( 1 - normal_cdf( critScaled, mu, sqrt(tau^2 + sei[i]^2) ) );
#       log_lik_sanity += -1 * log( 1 - normal_cdf( critScaled, 2, sqrt(2^2 + sei[i]^2) ) );
#     }
#   }
#   log_post = log_prior + log_lik;
# }
# "






# SILLIMAN EQ 2 (SIMPLIFIED FOR CASE WITH ALL SIGMA EQUAL)  -----------------------------------
# model.text <- "
#
# functions{
#
# 	real simple_prior(real mu, real tau, int k, real[] sei, real[] tcrit, real[] affirm){
# 		return 1 / ( sei[1]^2 + tau^2 );
# 	}
# }
#
# data{
# 	int<lower=0> k;
#   real sei[k];
#   real tcrit[k];
#   real affirm[k];
# 	real y[k];
# }
#
# parameters{
#   real mu;
# 	real<lower=0> tau;
# }
#
#
# model{
#   // this is to remove prior, as a sanity check:
#   // target += 0;
#   //see 'foundational ideas' here: https://vasishth.github.io/bayescogsci/book/sec-firststan.html
# 	target += log( simple_prior(mu, tau, k, sei, tcrit, affirm) );
# 	for(i in 1:k) {
#       if ( affirm[i] == 0 ) {
#         y[i] ~ normal( mu, sqrt(tau^2 + sei[i]^2) ) T[ , tcrit[i] * sei[i] ];
#       } else if ( affirm[i] == 1 ) {
#         y[i] ~ normal( mu, sqrt(tau^2 + sei[i]^2) ) T[ tcrit[i] * sei[i] , ];
#       }
# 	}
# }
#
# // this chunk doesn't actually affect the model that's being fit to the data;
# //  it's just re-calculating the prior, lkl, and post to return to user
# // Stan docs: 'Nothing in the generated quantities block affects the sampled parameter values. The block is executed only after a sample has been generated'
#
# generated quantities{
#   real log_lik = 0;
#   real log_prior = log( simple_prior(mu, tau, k, sei, tcrit, affirm) );
#   real log_post;
#   // this is just an intermediate quantity for log_lik
#   // will be equal to UU or LL above, depending on affirm status
#   real critScaled;
#
#    // versions that are evaluated at a SPECIFIC (mu=2, tau=2) so that we can compare
#    //  to R functions for MAP, MLE, etc.
#    real log_lik_sanity = 0;
#    real log_prior_sanity = log( simple_prior(2, 2, k, sei, tcrit, affirm) );
#
#   for ( i in 1:k ){
#     log_lik += normal_lpdf( y[i] | mu, sqrt(tau^2 + sei[i]^2) );
#     log_lik_sanity += normal_lpdf( y[i] | 2, sqrt(2^2 + sei[i]^2) );
#
#     critScaled = tcrit[i] * sei[i];
#
#     // https://mc-stan.org/docs/2_20/reference-manual/sampling-statements-section.html
#     // see 'Truncation with upper bounds in Stan' section
#     // nonaffirm case:
#     if ( y[i] <= critScaled ) {
#     // from sanity checks in doParallel, I know this matches joint_nll_2
#       log_lik += -1 * normal_lcdf(critScaled | mu, sqrt(tau^2 + sei[i]^2) );
#       log_lik_sanity += -1 * normal_lcdf(critScaled | 2, sqrt(2^2 + sei[i]^2) );
#
#     // affirm case:
#     } else if ( y[i] > critScaled ) {
#       log_lik += -1 * log( 1 - normal_cdf( critScaled, mu, sqrt(tau^2 + sei[i]^2) ) );
#       log_lik_sanity += -1 * log( 1 - normal_cdf( critScaled, 2, sqrt(2^2 + sei[i]^2) ) );
#     }
#   }
#   log_post = log_prior + log_lik;
# }
# "



# SIMPLE PRIOR - SILLIMAN  -----------------------------------
# LL and UU: cutpoints on RAW scale, not Z-scores
# tau: SD, not variance
# to do:
#  - can get rid of mustarL and LL; they are just placeholders now
# this model.text is from "2022-2-22 stan for Sherlock"
# model.text <- "
# 
# functions{
# 
# 	real simple_prior(real mu, real tau, int k, real[] sei, real[] tcrit, real[] affirm){
# 
#     // initialize the sum
#     real termA;
#     termA = 0;
# 
#     // Silliman prior in Eq (2), but they parameterize in terms of tau^2 - IS THAT OK?
#     for ( i in 1:k ) {
#     //Silliman version
#       termA = termA + ( 1 / ( sei[i]^2 + tau^2 )^2 );
#       //termA = termA + ( 1 / ( sei[i]^2 + tau^2 ) );
#     }
# 
# 		return (0.5*termA)^0.5;
# 	}
# }
# 
# data{
# 	int<lower=0> k;
#   real sei[k];
#   real tcrit[k];
#   real affirm[k];
# 	real y[k];
# }
# 
# parameters{
#   real mu;
# 	real<lower=0> tau;
# }
# 
# 
# model{
#   // this is to remove prior, as a sanity check:
#   // target += 0;
#   //see 'foundational ideas' here: https://vasishth.github.io/bayescogsci/book/sec-firststan.html
# 	target += log( simple_prior(mu, tau, k, sei, tcrit, affirm) );
# 	for(i in 1:k) {
#       if ( affirm[i] == 0 ) {
#         y[i] ~ normal( mu, sqrt(tau^2 + sei[i]^2) ) T[ , tcrit[i] * sei[i] ];
#       } else if ( affirm[i] == 1 ) {
#         y[i] ~ normal( mu, sqrt(tau^2 + sei[i]^2) ) T[ tcrit[i] * sei[i] , ];
#       }
# 	}
# }
# 
# // this chunk doesn't actually affect the model that's being fit to the data;
# //  it's just re-calculating the prior, lkl, and post to return to user
# // Stan docs: 'Nothing in the generated quantities block affects the sampled parameter values. The block is executed only after a sample has been generated'
# 
# generated quantities{
#   real log_lik = 0;
#   real log_prior = log( simple_prior(mu, tau, k, sei, tcrit, affirm) );
#   real log_post;
#   // this is just an intermediate quantity for log_lik
#   // will be equal to UU or LL above, depending on affirm status
#   real critScaled;
# 
#    // versions that are evaluated at a SPECIFIC (mu=2, tau=2) so that we can compare
#    //  to R functions for MAP, MLE, etc.
#    real log_lik_sanity = 0;
#    real log_prior_sanity = log( simple_prior(2, 2, k, sei, tcrit, affirm) );
# 
#   for ( i in 1:k ){
#     log_lik += normal_lpdf( y[i] | mu, sqrt(tau^2 + sei[i]^2) );
#     log_lik_sanity += normal_lpdf( y[i] | 2, sqrt(2^2 + sei[i]^2) );
# 
#     critScaled = tcrit[i] * sei[i];
# 
#     // https://mc-stan.org/docs/2_20/reference-manual/sampling-statements-section.html
#     // see 'Truncation with upper bounds in Stan' section
#     // nonaffirm case:
#     if ( y[i] <= critScaled ) {
#     // from sanity checks in doParallel, I know this matches joint_nll_2
#       log_lik += -1 * normal_lcdf(critScaled | mu, sqrt(tau^2 + sei[i]^2) );
#       log_lik_sanity += -1 * normal_lcdf(critScaled | 2, sqrt(2^2 + sei[i]^2) );
# 
#     // affirm case:
#     } else if ( y[i] > critScaled ) {
#       log_lik += -1 * log( 1 - normal_cdf( critScaled, mu, sqrt(tau^2 + sei[i]^2) ) );
#       log_lik_sanity += -1 * log( 1 - normal_cdf( critScaled, 2, sqrt(2^2 + sei[i]^2) ) );
#     }
#   }
#   log_post = log_prior + log_lik;
# }
# "

# MAIN VERSION (JEFFREYS PRIOR W/O CHANGE OF VARIABLES)  -----------------------------------
# LL and UU: cutpoints on RAW scale, not Z-scores
# tau: SD, not variance
# to do:
#  - can get rid of mustarL and LL; they are just placeholders now
# this model.text is from "2022-2-22 stan for Sherlock"
# model.text <- "
# 
# functions{
# 
# 	real jeffreys_prior(real mu, real tau, int k, real[] sei, real[] tcrit, real[] affirm){
# 
# 	  // these will be overwritten for EACH observation
# 		real mustarL;
# 		real mustarU;
# 		real alphaL;
# 		real alphaU;
# 		real kmm;
# 		real kms;
# 		real kss;
# 		matrix[2,2] fishinfo;
#     real sigma;
# 		real LL;
# 		real UU;
# 		// will just be set to 1
# 		int n;
# 
# 
# 		// this will be the TOTALS for all observations
# 		matrix[2,2] fishinfototal;
# 		fishinfototal[1,1] = 0;
#   	fishinfototal[1,2] = 0;
#   	fishinfototal[2,1] = 0;
#   	fishinfototal[2,2] = 0;
# 
# 
# 		// MM: build a Fisher info matrix for EACH observation
# 		for (i in 1:k) {
# 
# 		  // MARGINAL SD for this one observation
# 		  sigma = sqrt(tau^2 + sei[i]^2);
# 
# 		  // depending on whether study is affirmative, set truncation limits
# 		  // for THIS study, given its SE
# 		  if ( affirm[i] == 0 ) {
# 		  		UU = tcrit[i] * sei[i];
# 		  		// standardized truncation limits
# 		  		mustarL = -999;
#   		    mustarU = (UU - mu) / sigma;
# 		  } else if ( affirm[i] == 1 ) {
# 		      LL = tcrit[i] * sei[i];
# 		      // standardized truncation limits
# 		  		mustarL = (LL - mu) / sigma;
#   		    mustarU = 999;
# 		  }
# 
#   		// because EACH fisher info below has n=1 only
#   		n = 1;
# 
#   		// beginning of stuff that is not modified at all from TNE,
#   		//  *except* for the change-of-variables terms in kms and kss applied in a
#   		//   final code block
#   		// note that normal_lpdf, etc., parameterize in terms of SD, not var
#   		//  the (0,1) below are *not* start values for MCMC
#   		alphaL = exp( normal_lpdf(mustarL | 0, 1) -
#   	                log_diff_exp( normal_lcdf(mustarU | 0, 1),
#   	                normal_lcdf(mustarL | 0, 1) ) );
# 
#   		alphaU = exp( normal_lpdf(mustarU | 0, 1) -
#    	                log_diff_exp( normal_lcdf(mustarU | 0, 1),
#    	                normal_lcdf(mustarL | 0, 1) ) );
# 
#   		// second derivatives for Fisher info
#   		// wrt sigma, not tau
#   		kmm = -n/sigma^2 + n/sigma^2 * ((alphaU-alphaL)^2 + alphaU*mustarU- alphaL*mustarL);
#   		kms = -2*n/sigma^2 * (alphaL - alphaU) +
#   	   		  n/sigma^2 * (alphaL - alphaU + (alphaU*mustarU^2 - alphaL*mustarL^2) +
#   			  				(alphaL-alphaU) * (alphaL*mustarL - alphaU*mustarU));
#   		kss = n/sigma^2 - 3*n/sigma^2 * (1 + mustarL*alphaL - mustarU*alphaU) +
#   	   			n/sigma^2 * (mustarU*alphaU*(mustarU^2 - 2) - mustarL*alphaL*(mustarL^2 - 2) +
#   								(alphaU*mustarU - alphaL*mustarL)^2);
# 
# 
#   		// change of variables to get derivs wrt tau
#   		//kms = kms * tau / sqrt(tau^2 + sei[i]^2);
#   		//kss = kss * tau^2 / (tau^2 + sei[i]^2);
# 
#   		fishinfo[1,1] = -kmm;
#       fishinfo[1,2] = -kms;
#       fishinfo[2,1] = -kms;
#       fishinfo[2,2] = -kss;
# 
#   		// MM: add the new fisher info to the total one
#   		fishinfototal = fishinfototal + fishinfo;
# 		}
# 		return sqrt(determinant(fishinfototal));
# 	}
# }
# 
# data{
# 	int<lower=0> k;
#   real sei[k];
#   real tcrit[k];
#   real affirm[k];
# 	real y[k];
# }
# 
# parameters{
#   real mu;
# 	real<lower=0> tau;
# }
# 
# 
# model{
#   // this is to remove prior, as a sanity check:
#   // target += 0;
#   //see 'foundational ideas' here: https://vasishth.github.io/bayescogsci/book/sec-firststan.html
# 	target += log( jeffreys_prior(mu, tau, k, sei, tcrit, affirm) );
# 	for(i in 1:k) {
#       if ( affirm[i] == 0 ) {
#         y[i] ~ normal( mu, sqrt(tau^2 + sei[i]^2) ) T[ , tcrit[i] * sei[i] ];
#       } else if ( affirm[i] == 1 ) {
#         y[i] ~ normal( mu, sqrt(tau^2 + sei[i]^2) ) T[ tcrit[i] * sei[i] , ];
#       }
# 	}
# }
# 
# // this chunk doesn't actually affect the model that's being fit to the data;
# //  it's just re-calculating the prior, lkl, and post to return to user
# // Stan docs: 'Nothing in the generated quantities block affects the sampled parameter values. The block is executed only after a sample has been generated'
# 
# generated quantities{
#   real log_lik = 0;
#   real log_prior = log( jeffreys_prior(mu, tau, k, sei, tcrit, affirm) );
#   real log_post;
#   // this is just an intermediate quantity for log_lik
#   // will be equal to UU or LL above, depending on affirm status
#   real critScaled;
# 
#   // versions that are evaluated at a SPECIFIC (mu=2, tau=2) so that we can compare
#   //  to R functions for MAP, MLE, etc.
#   real log_lik_sanity = 0;
#   real log_prior_sanity = log( jeffreys_prior(2, 2, k, sei, tcrit, affirm) );
# 
#   for ( i in 1:k ){
#     log_lik += normal_lpdf( y[i] | mu, sqrt(tau^2 + sei[i]^2) );
#     log_lik_sanity += normal_lpdf( y[i] | 2, sqrt(2^2 + sei[i]^2) );
# 
#     critScaled = tcrit[i] * sei[i];
# 
#     // https://mc-stan.org/docs/2_20/reference-manual/sampling-statements-section.html
#     // see 'Truncation with upper bounds in Stan' section
#     // nonaffirm case:
#     if ( y[i] <= critScaled ) {
#     // from sanity checks in doParallel, I know this matches joint_nll_2
#       log_lik += -1 * normal_lcdf(critScaled | mu, sqrt(tau^2 + sei[i]^2) );
#       log_lik_sanity += -1 * normal_lcdf(critScaled | 2, sqrt(2^2 + sei[i]^2) );
# 
#     // affirm case:
#     } else if ( y[i] > critScaled ) {
#       log_lik += -1 * log( 1 - normal_cdf( critScaled, mu, sqrt(tau^2 + sei[i]^2) ) );
#       log_lik_sanity += -1 * log( 1 - normal_cdf( critScaled, 2, sqrt(2^2 + sei[i]^2) ) );
#     }
#   }
#   log_post = log_prior + log_lik;
# }
# "

# R implementation of affirm part that I know matches joint_nll2 and dtruncnorm:
# term1 = dnorm(x = dpa$yi,
#               mean = 2,
#               sd = sqrt(2^2 + dpa$vi),
#               log = TRUE )
# critScaled = dpa$tcrit * sqrt(dpa$vi)
# term2 = log( 1 - pnorm(q = critScaled,
#                        mean = 2,
#                        sd = sqrt(2^2 + dpa$vi) ) )


# necessary to prevent ReadRDS errors in which cores try to work with other cores' intermediate results
# https://groups.google.com/g/stan-users/c/8snqQTTfWVs?pli=1
options(mc.cores = parallel::detectCores())

# "isystem" arg is just a placeholder to avoid Stan's not understanding special characters
#  in getwd(), even though we don't actually use the dir at all
# note: removing the isystem arg does NOT fix the very sporadic "error reading from connection" on cluster
cat( paste("\n init_stan_model: about to call stan_model") )
stan.model <- stan_model(model_code = model.text,
                         isystem = "~/Desktop")


cat( paste("\n init_stan_model: done calling stan_model") )