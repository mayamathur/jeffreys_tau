
# NOTES ---------------------------------------------------------------

# "draw" always refers to within-study draw
# "study set" refers to all draws, observed and observed, for a given study


# ESTIMATION METHOD FNS ----------------------------------------------

# ~ Estimation Methods Structured for Use Inside run_method_safe --------

# Because these fns are run inside run_method_safe, the latter will handle editing rep.res
#  All of these fns should take get.CIs as an argument and return CIs as c(NA, NA) if not wanted
# Fns in this category need to return a dataframe with the below structure, although it's okay if they don't return all of these names since run_method_safe will handle that. Note that this is a LIST containing a dataframe called "stats", not just the dataframe; this allows easy extension in case you want to return other objects, like model objects.

# .yi: published point estimates
# .sei: their SEs
# .tcrit: critical values on t or z scale for each study; can just use qnorm(.975) by default
# .Mu.start: optimizer starting value for meta-analysis Mu
# .Tt.start: optimizer starting value for meta-analysis tau
# .stan.adapt_delta: passed to rstan
# .stan.maxtreedepth: same
#  we should later use ellipsis to allow passing arbitrary args to rstan
estimate_jeffreys = function(.yi,
                             .sei,
                             
                             .Mu.start,
                             .Tt.start,
                             .stan.adapt_delta = 0.8,
                             .stan.maxtreedepth = 10 ) {
  
  # stan.model (used later) is compiled OUTSIDE this fn in doParallel to avoid 
  #  issues with nodes competing with one another
  
  
  # prepare to capture warnings from Stan
  stan.warned = 0
  stan.warning = NA
  
  # set start values for sampler
  init.fcn = function(o){ list(mu = .Mu.start,
                               tau = .Tt.start ) }
  
  
  # like tryCatch, but captures warnings without stopping the function from
  #  returning its results
  withCallingHandlers({
    
    # necessary to prevent ReadRDS errors in which cores try to work with other cores' intermediate results
    # https://groups.google.com/g/stan-users/c/8snqQTTfWVs?pli=1
    options(mc.cores = parallel::detectCores())
    
    cat( paste("\n estimate_jeffreys flag 2: about to call sampling") )
    
    post = sampling(stan.model,
                    cores = 1,
                    refresh = 0,
                    data = list( k = length(.yi),
                                 sei = .sei,
                                 y = .yi ),
                    
                    #iter = p$stan.iter,   
                    control = list(max_treedepth = p$stan.maxtreedepth,
                                   adapt_delta = p$stan.adapt_delta),
                    
                    init = init.fcn)
    
    
  }, warning = function(condition){
    stan.warned <<- 1
    stan.warning <<- condition$message
  } )
  
  cat( paste("\n estimate_jeffreys flag 3: about to call postSumm") )
  postSumm = summary(post)$summary
  if (is.null(postSumm)) stop("In stan, postSumm is null")
  
  # pull out best iterate to pass to MAP optimization later
  ext = rstan::extract(post) # a vector of all post-WU iterates across all chains
  best.ind = which.max(ext$log_post)  # single iterate with best log-posterior should be very close to MAP
  
  
  # posterior means, posterior medians, modes, and max-LP iterate
  Mhat = c( postSumm["mu", "mean"],
            median( rstan::extract(post, "mu")[[1]] ),
            ext$mu[best.ind] )
  
  Shat = c( postSumm["tau", "mean"],
            median( rstan::extract(post, "tau")[[1]] ),
            ext$tau[best.ind] )
  
  # sanity check
  #expect_equal( Mhat[1], mean( rstan::extract(post, "mu")[[1]] ) )
  
  
  # SEs
  MhatSE = postSumm["mu", "se_mean"]
  ShatSE = postSumm["tau", "se_mean"]
  # how Stan estimates the SE: https://discourse.mc-stan.org/t/se-mean-in-print-stanfit/2869
  expect_equal( postSumm["mu", "sd"],
                sd( rstan::extract(post, "mu")[[1]] ) )
  expect_equal( MhatSE,
                postSumm["mu", "sd"] / sqrt( postSumm["mu", "n_eff"] ) )
  
  # CI limits
  S.CI = c( postSumm["tau", "2.5%"], postSumm["tau", "97.5%"] )
  M.CI = c( postSumm["mu", "2.5%"], postSumm["mu", "97.5%"] )
  # sanity check:
  myMhatCI = as.numeric( c( quantile( rstan::extract(post, "mu")[[1]], 0.025 ),
                            quantile( rstan::extract(post, "mu")[[1]], 0.975 ) ) )
  expect_equal(M.CI, myMhatCI)
  
  
  # the point estimates are length 2 (post means, then medians),
  #  but the inference is the same for each type of point estimate
  return( list( stats = data.frame( 
    
    Mhat = Mhat,
    Shat = Shat,
    
    MhatSE = MhatSE,
    ShatSE = ShatSE,
    
    # this will use same CI limits for all pt estimates
    MLo = M.CI[1],
    MHi = M.CI[2],
    
    SLo = S.CI[1],
    SHi = S.CI[2],
    
    stan.warned = stan.warned,
    stan.warning = stan.warning,
    MhatRhat = postSumm["mu", "Rhat"],
    ShatRhat = postSumm["tau", "Rhat"] ),
    
    post = post,
    postSumm = postSumm ) )
  
}


mw_est <- function(yi,  # the UNQUOTED name of the yi column in data, not the vector itself
                   vi,  # similar
                   data,
                   lab,
                   method = c("DL", "DL2"),
                   ci = .95,
                   Q.profile = FALSE,
                   pi = FALSE,
                   hksj = FALSE) {
  
  ## General Method-of-Moments Estimate for Tau-Squared (Eq. 6)
  .tau2MM <- function(yi, vi, ai) {
    yw <- sum(ai * yi) / sum(ai)
    
    tau2 <- (sum(ai * (yi - yw)^2) - (sum(ai * vi) - sum(ai^2 * vi) / sum(ai))) / (sum(ai) - sum(ai^2) / sum(ai))
    
    # Return only non-negative values
    max(0, tau2)
  }
  
  # Keep complete cases
  data <- data %>%
    select(c({{lab}}, {{ yi }}, {{ vi }})) %>%
    filter_all(all_vars(!is.na(.)))
  
  # Extract variables
  yi <- pull(data, {{ yi }})
  vi <- pull(data, {{ vi }})
  
  # Calculate DL tau2
  ai <- 1 / vi
  tau2 <- .tau2MM(yi = yi, vi = vi, ai = ai)
  
  # Calculate DL2 if requested
  if (method == "DL2") {
    ai <- 1 / (vi + tau2)
    tau2 <- .tau2MM(yi = yi, vi = vi, ai = ai)
  }
  
  ## Calculate inverse variance weights
  wi <- (1 / (vi + tau2)) # Eq. 2
  
  wy <- wi * yi
  
  mw <- sum(wy) / sum(wi) # Eq. 1
  
  # Get number of studies
  k <- length(yi)
  
  # Degrees of freedom
  df <- k - 1
  
  ### Heterogeneity
  
  ## Cochran's Q-statistic
  
  # Use DL inverse-variance weight (ai = 1 / vi)
  ai <- 1 / vi
  
  # Calculate Q-statistic
  Q <- sum(ai * ((yi - (sum(ai * yi) / sum(ai)))^2))
  
  # Calculate p-value
  Q.p <- pchisq(Q, df = df, lower.tail = FALSE)
  
  z.cdf <- qnorm(1 - ((1 - ci) / 2))
  
  # Method used to calculate heterogeneity and confidence intervals
  if (Q.profile) {
    
    ## Q-profile method to calculate heterogeneity and confidence 
    ## intervals (Eq. 14)
    
    ## Calculate I-squared and H
    
    # a. Use DL inverse-variance weight (ai = 1 / vi)
    ai <- 1 / vi
    
    # b. Calculate V
    V <- (k - 1) * sum(ai) / (sum(ai)^2 - sum(ai^2))
    
    # c. Calculate I-squared
    I.sq <- 100 * tau2 / (tau2 + V)
    
    # d. Calculate H-squared
    H <- (tau2 + V) / V
    
    # Calculate confidence intervals if k > 2
    if (k > 2) {
      
      ## Generalised Q-statistics bounds
      Q.lb <- qchisq((1 - ci) / 2, df, lower = TRUE)
      Q.ub <- qchisq((1 - ci) / 2, df, lower = FALSE)
      
      ## Lower confidence intervals
      # Initialise values
      tau2.tmp <- 0
      F.tau2 <- 1
      
      while (F.tau2 > 0) {
        
        # i. calculate weight for tau2
        W <- 1 / (vi + tau2.tmp)
        
        yW <- sum(yi * W) / sum(W)
        
        
        # ii. Calculate delta tau2
        F.tau2 <- sum(W * (yi - yW)^2) - Q.ub
        F.tau2 <- max(F.tau2, 0)
        
        delta.tau2 <- F.tau2 / sum((W^2) * (yi - yW)^2)
        
        # iii. Next iterative step
        tau2.tmp <- tau2.tmp + delta.tau2
      }
      
      tau2.tmp <- max(0, tau2.tmp)
      I.lower <- 100 * tau2.tmp / (tau2.tmp + V)
      H.lower <- (tau2.tmp + V) / V
      
      ## Upper confidence intervals
      # Initialise values
      tau2.tmp <- 0
      F.tau2 <- 1
      
      while (F.tau2 > 0) {
        
        # i. calculate weight for tau2
        W <- 1 / (vi + tau2.tmp)
        
        yW <- sum(yi * W) / sum(W)
        
        # ii. Calculate delta tau2
        F.tau2 <- sum(W * (yi - yW)^2) - Q.lb
        F.tau2 <- max(F.tau2, 0)
        
        delta.tau2 <- F.tau2 / sum((W^2) * (yi - yW)^2)
        
        # iii. Next iterative step
        tau2.tmp <- tau2.tmp + delta.tau2
      }
      
      tau2.tmp <- max(0, tau2.tmp)
      I.upper <- 100 * tau2.tmp / (tau2.tmp + V)
      H.upper <- (tau2.tmp + V) / V
    } else {
      I.lower <- NA
      I.upper <- NA
      H.lower <- NA
      H.upper <- NA
      
      pi.res <- NULL
    }
  } else {
    
    ## I-squared
    I.sq <- max(0, 100 * ((Q - df) / Q))
    
    ## H
    H <- sqrt(Q / df)
    
    # Calculate the SE of ln(H) aka B
    if (k > 2) {
      if (Q > k) {
        B <- 0.5 * ((log(Q) - log(df)) / (sqrt(2 * Q) - sqrt(2 * df - 1)))
      } else {
        B <- sqrt(1 / ((2 * df - 1) * (1 - (1 / (3 * (df - 1)^2)))))
      }
      
      # H confidence intervals
      H.lower <- max(1, exp(log(H) - z.cdf * B))
      H.upper <- exp(log(H) + z.cdf * B)
      
      # I-Squared
      I.L <- exp(0.5 * log(Q / df) - (z.cdf * B))
      I.U <- exp(0.5 * log(Q / df) + (z.cdf * B))
      
      I.lower <- max(0, 100 * ((I.L^2 - 1) / I.L^2))
      I.upper <- max(0, 100 * ((I.U^2 - 1) / I.U^2))
    } else {
      I.lower <- NA
      I.upper <- NA
      
      H.lower <- NA
      H.upper <- NA
      
      pi.res <- NULL
    }
  }
  
  # Save all heterogeneity statistics in a list
  het.test <- list(
    Q = Q, Q.p = Q.p,
    tau2 = tau2,
    I.sq = I.sq, I.lb = I.lower, I.ub = I.upper,
    H = H, H.lb = H.lower, H.ub = H.upper
  )
  
  # HKSJ adjustment
  if (hksj & k >= 3) {
    
    # Calculate HKSJ-adjusted standard error
    sew <- sqrt(sum(wi * (yi - mw)^2) / (df * sum(wi)))
    
    # Calculate t-value
    t <- mw / sew
    
    # Calculate p-value based on t-statistics and k - 1 df
    t.p <- 2 * pt(-abs(t), df = k - 1)
    
    # Output test for summary effect size
    mw.test <- list(t = t, p = t.p)
    
    # Confidence Intervals
    mw.lower <- mw - sew * qt(1 - ((1 - ci) / 2), df = df)
    mw.upper <- mw + sew * qt(1 - ((1 - ci) / 2), df = df)
    
    # Wald-type statistics
  } else {
    
    # Calculate standard error
    sew <- 1 / sqrt(sum(wi)) # Eq. 3
    
    # Calculate z-statistic and p-value
    z <- mw / sew
    z.p <- 2 * pnorm(-abs(z))
    
    # Confidence Intervals
    mw.lower <- mw - z.cdf * sew
    mw.upper <- mw + z.cdf * sew
    
    # Output test for summary effect size
    mw.test <- list(z = z, p = z.p)
  }
  
  # Prediction Interval (based on t-distribution)
  if (all(k > 2, pi)) {
    if (hksj) {
      pi.lower <- mw - sqrt(tau2 + sew^2) * qt(1 - ((1 - ci) / 2), df = df)
      pi.upper <- mw + sqrt(tau2 + sew^2) * qt(1 - ((1 - ci) / 2), df = df)
    } else {
      pi.lower <- mw - sqrt(tau2 + sew^2) * z.cdf
      pi.upper <- mw + sqrt(tau2 + sew^2) * z.cdf
    }
    
    # Output prediction interval
    pi.res <- list(pi.lb = pi.lower, pi.ub = pi.upper)
  } else {
    pi.res <- NULL
  }
  # Print data set
  cat(
    "\n", "Data used for meta-analysis",
    "\n",
    "\n",
    sep = ""
  )
  
  print(data)
  
  cat(
    "\n",
    sep = ""
  )
  
  # Output results as a list
  return(
    c(
      method = method,
      est = mw, se.est = sew,
      mw.test, ci.lb = mw.lower, ci.ub = mw.upper,
      het.test, pi.res
    )
  )
}

# example
if ( FALSE ){
  
  d = data.frame(yi = c(0.693147180559945, -0.510825623765991, 0.405465108108165, 
                        -0.693147180559945, -1.09861228866811, -2.19722457733622, -1.94591014905531, 
                        0, 0.405465108108164, -1.38629436111989),
                 vi = c(1.36666666666667, 0.483333333333333, 0.361111111111111, 0.675925925925926, 
                        2.52380952380952, 2.13526570048309, 2.17460317460317, 1.92, 0.738095238095238, 
                        1.16666666666667) )
  
  mw_est(yi = yi,
         vi = vi,

         data = d,
         method = "DL2",
         hksj = TRUE)
  
}







estimate_DL2 = function(.yi,
                             .sei,
                             
                             .Mu.start,
                             .Tt.start,
                             .stan.adapt_delta = 0.8,
                             .stan.maxtreedepth = 10 ) {
  
  ## General Method-of-Moments Estimate for Tau-Squared (Eq. 6)
  .tau2MM <- function(yi, vi, ai) {
    yw <- sum(ai * yi) / sum(ai)
    
    tau2 <- (sum(ai * (yi - yw)^2) - (sum(ai * vi) - sum(ai^2 * vi) / sum(ai))) / (sum(ai) - sum(ai^2) / sum(ai))
    
    # Return only non-negative values
    max(0, tau2)
  }
  
  
  # Calculate DL tau2
  vi= .sei^2
  ai <- 1 / vi
  tau2_init <- .tau2MM(yi = .yi, vi = vi, ai = ai)
  
  # Calculate DL2
  ai <- 1 / (vi + tau2_init)
  tau2 <- .tau2MM(yi = .yi, vi = vi, ai = ai)
  
  ## Calculate inverse variance weights
  wi <- (1 / (vi + tau2)) # Eq. 2
  
  wy <- wi * .yi
  
  mw <- sum(wy) / sum(wi) # Eq. 1
  
  # Get number of studies
  k <- length(.yi)
  
  # Degrees of freedom
  df <- k - 1
  
  ### Heterogeneity
  
  ## Cochran's Q-statistic
  
  # Use DL inverse-variance weight (ai = 1 / vi)
  ai <- 1 / vi
  
  # Calculate Q-statistic
  Q <- sum(ai * ((yi - (sum(ai * yi) / sum(ai)))^2))
  
  # Calculate p-value
  Q.p <- pchisq(Q, df = df, lower.tail = FALSE)
  
  z.cdf <- qnorm(1 - ((1 - ci) / 2))
  
  # Method used to calculate heterogeneity and confidence intervals
  if (Q.profile) {
    
    ## Q-profile method to calculate heterogeneity and confidence 
    ## intervals (Eq. 14)
    
    ## Calculate I-squared and H
    
    # a. Use DL inverse-variance weight (ai = 1 / vi)
    ai <- 1 / vi
    
    # b. Calculate V
    V <- (k - 1) * sum(ai) / (sum(ai)^2 - sum(ai^2))
    
    # c. Calculate I-squared
    I.sq <- 100 * tau2 / (tau2 + V)
    
    # d. Calculate H-squared
    H <- (tau2 + V) / V
    
    # Calculate confidence intervals if k > 2
    if (k > 2) {
      
      ## Generalised Q-statistics bounds
      Q.lb <- qchisq((1 - ci) / 2, df, lower = TRUE)
      Q.ub <- qchisq((1 - ci) / 2, df, lower = FALSE)
      
      ## Lower confidence intervals
      # Initialise values
      tau2.tmp <- 0
      F.tau2 <- 1
      
      while (F.tau2 > 0) {
        
        # i. calculate weight for tau2
        W <- 1 / (vi + tau2.tmp)
        
        yW <- sum(yi * W) / sum(W)
        
        
        # ii. Calculate delta tau2
        F.tau2 <- sum(W * (yi - yW)^2) - Q.ub
        F.tau2 <- max(F.tau2, 0)
        
        delta.tau2 <- F.tau2 / sum((W^2) * (yi - yW)^2)
        
        # iii. Next iterative step
        tau2.tmp <- tau2.tmp + delta.tau2
      }
      
      tau2.tmp <- max(0, tau2.tmp)
      I.lower <- 100 * tau2.tmp / (tau2.tmp + V)
      H.lower <- (tau2.tmp + V) / V
      
      ## Upper confidence intervals
      # Initialise values
      tau2.tmp <- 0
      F.tau2 <- 1
      
      while (F.tau2 > 0) {
        
        # i. calculate weight for tau2
        W <- 1 / (vi + tau2.tmp)
        
        yW <- sum(yi * W) / sum(W)
        
        # ii. Calculate delta tau2
        F.tau2 <- sum(W * (yi - yW)^2) - Q.lb
        F.tau2 <- max(F.tau2, 0)
        
        delta.tau2 <- F.tau2 / sum((W^2) * (yi - yW)^2)
        
        # iii. Next iterative step
        tau2.tmp <- tau2.tmp + delta.tau2
      }
      
      tau2.tmp <- max(0, tau2.tmp)
      I.upper <- 100 * tau2.tmp / (tau2.tmp + V)
      H.upper <- (tau2.tmp + V) / V
    } else {
      I.lower <- NA
      I.upper <- NA
      H.lower <- NA
      H.upper <- NA
      
      pi.res <- NULL
    }
  } else {
    
    ## I-squared
    I.sq <- max(0, 100 * ((Q - df) / Q))
    
    ## H
    H <- sqrt(Q / df)
    
    # Calculate the SE of ln(H) aka B
    if (k > 2) {
      if (Q > k) {
        B <- 0.5 * ((log(Q) - log(df)) / (sqrt(2 * Q) - sqrt(2 * df - 1)))
      } else {
        B <- sqrt(1 / ((2 * df - 1) * (1 - (1 / (3 * (df - 1)^2)))))
      }
      
      # H confidence intervals
      H.lower <- max(1, exp(log(H) - z.cdf * B))
      H.upper <- exp(log(H) + z.cdf * B)
      
      # I-Squared
      I.L <- exp(0.5 * log(Q / df) - (z.cdf * B))
      I.U <- exp(0.5 * log(Q / df) + (z.cdf * B))
      
      I.lower <- max(0, 100 * ((I.L^2 - 1) / I.L^2))
      I.upper <- max(0, 100 * ((I.U^2 - 1) / I.U^2))
    } else {
      I.lower <- NA
      I.upper <- NA
      
      H.lower <- NA
      H.upper <- NA
      
      pi.res <- NULL
    }
  }
  
  # Save all heterogeneity statistics in a list
  het.test <- list(
    Q = Q, Q.p = Q.p,
    tau2 = tau2,
    I.sq = I.sq, I.lb = I.lower, I.ub = I.upper,
    H = H, H.lb = H.lower, H.ub = H.upper
  )
  
  # HKSJ adjustment
  if (hksj & k >= 3) {
    
    # Calculate HKSJ-adjusted standard error
    sew <- sqrt(sum(wi * (yi - mw)^2) / (df * sum(wi)))
    
    # Calculate t-value
    t <- mw / sew
    
    # Calculate p-value based on t-statistics and k - 1 df
    t.p <- 2 * pt(-abs(t), df = k - 1)
    
    # Output test for summary effect size
    mw.test <- list(t = t, p = t.p)
    
    # Confidence Intervals
    mw.lower <- mw - sew * qt(1 - ((1 - ci) / 2), df = df)
    mw.upper <- mw + sew * qt(1 - ((1 - ci) / 2), df = df)
    
    # Wald-type statistics
  } else {
    
    # Calculate standard error
    sew <- 1 / sqrt(sum(wi)) # Eq. 3
    
    # Calculate z-statistic and p-value
    z <- mw / sew
    z.p <- 2 * pnorm(-abs(z))
    
    # Confidence Intervals
    mw.lower <- mw - z.cdf * sew
    mw.upper <- mw + z.cdf * sew
    
    # Output test for summary effect size
    mw.test <- list(z = z, p = z.p)
  }
  
  # Prediction Interval (based on t-distribution)
  if (all(k > 2, pi)) {
    if (hksj) {
      pi.lower <- mw - sqrt(tau2 + sew^2) * qt(1 - ((1 - ci) / 2), df = df)
      pi.upper <- mw + sqrt(tau2 + sew^2) * qt(1 - ((1 - ci) / 2), df = df)
    } else {
      pi.lower <- mw - sqrt(tau2 + sew^2) * z.cdf
      pi.upper <- mw + sqrt(tau2 + sew^2) * z.cdf
    }
    
    # Output prediction interval
    pi.res <- list(pi.lb = pi.lower, pi.ub = pi.upper)
  } else {
    pi.res <- NULL
  }
  # Print data set
  cat(
    "\n", "Data used for meta-analysis",
    "\n",
    "\n",
    sep = ""
  )
  
  print(data)
  
  cat(
    "\n",
    sep = ""
  )
  
  # Output results as a list
  return(
    c(
      method = method,
      est = mw, se.est = sew,
      mw.test, ci.lb = mw.lower, ci.ub = mw.upper,
      het.test, pi.res
    )
  )
  
  
  # # the point estimates are length 2 (post means, then medians),
  # #  but the inference is the same for each type of point estimate
  # return( list( stats = data.frame( 
  #   
  #   Mhat = Mhat,
  #   Shat = Shat,
  #   
  #   MhatSE = MhatSE,
  #   ShatSE = ShatSE,
  #   
  #   # this will use same CI limits for all pt estimates
  #   MLo = M.CI[1],
  #   MHi = M.CI[2],
  #   
  #   SLo = S.CI[1],
  #   SHi = S.CI[2],
  #   
  #   stan.warned = stan.warned,
  #   stan.warning = stan.warning,
  #   MhatRhat = postSumm["mu", "Rhat"],
  #   ShatRhat = postSumm["tau", "Rhat"] ),
  #   
  #   post = post,
  #   postSumm = postSumm ) )
  
}


# nicely report a metafor or robumeta object with optional suffix to denote which model
report_meta = function(.mod,
                       .mod.type = "rma",  # "rma" or "robu"
                       .suffix = "") {
  
  if ( !is.null(.mod) ) {
    
    
    if ( .mod.type == "rma" ) {
      tau.CI = tau_CI(.mod)
      .res = data.frame( .mod$b,
                         .mod$ci.lb,
                         .mod$ci.ub,
                         
                         sqrt(.mod$tau2),
                         tau.CI[1],
                         tau.CI[2] )
    } 
    
    
    if ( .mod.type == "robu" ) {
      
      .res = data.frame( .mod$b.r,
                         .mod$reg_table$CI.L,
                         .mod$reg_table$CI.U,
                         
                         sqrt(.mod$mod_info$tau.sq),
                         NA,
                         NA )
    } 
    
  } else {
    .res = data.frame( rep(NA, 6) )
  }
  
  
  names(.res) = paste( c("Mhat", "MLo", "MHi", "Shat", "SLo", "SHi"), .suffix, sep = "" )
  row.names(.res) = NULL
  
  return( list(stats = .res) )
}



# HELPERS FOR ABOVE ESTIMATION METHODS ----------------------------



# Fisher info when taking derivatives wrt mu and tau
# This fn ONLY handles nonaffirm results, but could easily be adapted to handle
#  affirms by changing tcrit.
# important: note that in this fn, critical value is on t/z scale, NOT raw scale
#  vs. in E_fisher_TNE, .b is on raw scale
E_fisher_RTMA = function( .sei, .Mu, .Tt, .tcrit = qnorm(0.975) ) {
  
  Efish.list = lapply( X = as.list(.sei),
                       FUN = function(.s) {
                         
                         # for this observation
                         sei = .s
                         mu = .Mu
                         tau = .Tt
                         tcrit = .tcrit  # currently assumed to be a scalar
                         if ( length(tcrit) > 1 ) tcrit = tcrit[1] #OBVIOUSLY NEEDS TO BE GENERALIZED
                         
                         fishinfo = matrix( NA, nrow = 2, ncol = 2 )
                         
                         # from body of R's get_D11_num:
                         e2 = sei^2 + tau^2
                         e3 = sqrt(e2)
                         e5 = sei * tcrit - mu
                         e6 = e5/e3
                         e7 = dnorm(e6, 0, 1)
                         # Stan version:
                         # e7 = exp( normal_lpdf(e6 | 0, 1) )
                         e8 = pnorm(e6)
                         #e8 = exp( normal_lcdf(e6 | 0, 1 ) )
                         kmm = -(1/e2 - (e5/(e2 * e8) + e7 * e3/(e8 * e3)^2) * e7/e3)
                         
                         # from body of R's get_D12_num:
                         e2 = sei^2 + tau^2
                         e3 = sqrt(e2)
                         e5 = sei * tcrit - mu
                         # e6 is scaled critical value:
                         e6 = e5/e3
                         e7 = pnorm(e6)
                         # e7 = exp( normal_lcdf(e6 | 0, 1 ) )
                         e8 = e2^2
                         e9 = dnorm(e6, 0, 1)
                         #e9 = exp( normal_lpdf(e6 | 0, 1) )
                         
                         # my own expectation of .yi - .mu:
                         expectation1 = -sqrt(sei^2 + tau^2) * e9/e7
                         kms = -(tau * (((e7/e3 - e5 * e9/e2)/(e7 * e3)^2 - e5^2/(e8 *
                                                                                    e7 * e3)) * e9 + 2 * ((expectation1)/e8)))
                         
                         
                         # from body of R's get_D22_num:
                         e1 = tau^2
                         e3 = sei^2 + e1
                         e5 = sei * tcrit - mu
                         e6 = sqrt(e3)
                         # e7 is scaled crit value:
                         e7 = e5/e6
                         e8 = pnorm(e7)
                         # e8 = exp( normal_lcdf(e7 | 0, 1 ) )
                         e9 = dnorm(e7, 0, 1)
                         # e9 = exp( normal_lpdf(e7 | 0, 1 ) )
                         e10 = e5 * e9
                         e11 = e8 * e6
                         e13 = e10/e11
                         # *replace this one with its expectation:
                         # e15 = (.yi - .mu)^2/e3
                         # expectation of (.yi - .mu)^2:
                         expectation2 = (sei^2 + tau^2)*(1 - e7 * e9/e8)
                         e15 = (expectation2)/e3
                         
                         kss = (e13 + e15 - (e1 * (e5 * ((e8/e6 - e10/e3)/e11^2 -
                                                           e5^2/(e3^2 * e8 * e6)) * e9 + 2 * ((e13 + 2 * e15 -
                                                                                                 1)/e3)) + 1))/e3
                         
                         
                         fishinfo[1,1] = -kmm
                         fishinfo[1,2] = -kms
                         fishinfo[2,1] = -kms
                         fishinfo[2,2] = -kss
                         
                         return(fishinfo)
                         
                         
                       })
  
  # add all the matrices entrywise
  # https://stackoverflow.com/questions/11641701/sum-a-list-of-matrices
  Efish.all = Reduce('+', Efish.list) 
  
  # cat("\nFirst observation Efish:")
  # print(Efish.list[[1]])
  
  return(Efish.all)
}


lprior = function(.sei, .Mu, .Tt, .tcrit) {
  Efish = E_fisher_RTMA( .sei = .sei, .Mu = .Mu, .Tt = .Tt, .tcrit = .tcrit )
  log( sqrt( det(Efish) ) )
}



# ~ Other Helpers ---------------

# taken from TNE 2022-2-26
get_optimx_dataframe = function( .yi,
                                 .sei,
                                 .tcrit,
                                 .usePrior,
                                 .par2is,
                                 .Mu.start,
                                 .par2.start ) {
  
  
  ox.methods <- c('Nelder-Mead', 'BFGS', 'CG', 'L-BFGS-B', 'nlm', 'nlminb', 'spg', 'ucminf',
                  'newuoa', 'bobyqa', 'nmkb', 'hjkb', 'Rcgmin', 'Rvmmin')
  
  l = optimx( par = c(.Mu.start, .par2.start),
              fn = function(..pars) as.numeric( nlpost_jeffreys_RTMA( .pars = ..pars,
                                                                      .par2is = .par2is,
                                                                      .yi = .yi,
                                                                      .sei = .sei,
                                                                      .tcrit = .tcrit,
                                                                      .usePrior = .usePrior ) ),
              method = ox.methods )
  
  l$opt.method = row.names(l)
  
  # transform second parameter so it's always Shat instead of Vhat
  if ( .par2is == "T2t" ) { l$p2 = sqrt(l$p2) }
  
  l2 = l %>% select(opt.method, p1, p2, convcode, value, kkt1, kkt2) 
  
  l2 = l2 %>% rename( Mhat = p1, Shat = p2, nll = value )
  
  w = pivot_wider(l2, 
                  names_from = "opt.method",
                  values_from = c("Mhat", "Shat", "convcode", "nll", "kkt1", "kkt2"),
                  names_glue = "optimx.{opt.method}.{.value}")
  
  
  if ( length( l$p1[ l$convcode == 0 ] ) > 0 ){
    
    # only keep the ones that had values for Mhat, Shat (not ones that didn't even give a value)
    l = l[ !is.na(l$p1) & !is.na(l$p2), ]
    
    
    #**optimizers that converged AND
    # had a small gradient (kkt1) AND
    # had a positive-definite Hessian (kkt2)
    lc = l[ l$convcode == 0 & l$kkt1 == TRUE & l$kkt2 == TRUE, ]
    
    # index of optimizer with the best nll
    lc.winner.ind = which.min(lc$value)
    
    
    # Mhat.winner is the Mhat of the optimizer with the best nll, OF converged ones
    # catch case in which no optimizers converged
    if ( length(lc.winner.ind > 0) ) {
      Mhat.winner = lc$p1[lc.winner.ind]
      Shat.winner = lc$p2[lc.winner.ind]
    } else {
      Mhat.winner = Shat.winner = NA
    }
    
    
    # **note that this is the criterion for agreement
    l$agree.Mhat = abs(l$p1 - Mhat.winner) < 0.01
    l$agree.Shat = abs(l$p2 - Shat.winner) < 0.01
    
    # sanity check: look at differences of non-agreers from Mhat.winner
    #l$p1[ l$agree.Mhat == FALSE ] - Mhat.winner
    
    
    # get lc again now that we have the agreement indicator
    lc = l[ l$convcode == 0 & l$kkt1 == TRUE & l$kkt2 == TRUE, ]
    w$optimx.Nconvergers = nrow(lc)
    w$optimx.convergers = paste( lc$opt.method, collapse = " ")
    
    w$optimx.Mhat.winner = Mhat.winner
    w$optimx.Pagree.Mhat.winner = sum(l$agree.Mhat)/nrow(l)
    # number and proportion of optimizers that converged that agreed with mode:
    w$optimx.Nagree.of.convergers.Mhat.winner = sum(lc$agree.Mhat)
    w$optimx.Pagree.of.convergers.Mhat.winner = sum(lc$agree.Mhat)/nrow(lc)
    w$optimx.Mhat.agreers = paste( l$opt.method[ l$agree.Mhat == TRUE ], collapse = " ")
    w$optimx.Mhat.convergers.agreers = paste( lc$opt.method[ lc$agree.Mhat == TRUE ], collapse = " ")
    
    w$optimx.Shat.winner = Shat.winner
    w$optimx.Pagree.Shat.winner = sum(l$agree.Shat)/nrow(l)
    w$optimx.Nagree.of.convergers.Shat.winner = sum(lc$agree.Shat)
    w$optimx.Pagree.of.convergers.Shat.winner = sum(lc$agree.Shat)/nrow(lc)
    w$optimx.Shat.agreers = paste( l$opt.method[ l$agree.Shat == TRUE ], collapse = " ")
    w$optimx.Shat.convergers.agreers = paste( lc$opt.method[ lc$agree.Shat == TRUE ], collapse = " ")
    
  } else {
    w$optimx.Nconvergers = NA
    w$optimx.convergers = NA
    w$optimx.Mhat.winner = NA
    w$optimx.Pagree.Mhat.winner = NA
    w$optimx.Nagree.of.convergers.Mhat.winner = NA
    w$optimx.Pagree.of.convergers.Mhat.winner = NA
    w$optimx.Mhat.agreers = NA
    w$optimx.Mhat.convergers.agreers = NA
    
    w$optimx.Shat.winner = NA
    w$optimx.Nagree.of.convergers.Shat.winner = NA
    w$optimx.Pagree.of.convergers.Shat.winner = NA
    w$optimx.Pagree.Shat.winner = NA
    w$optimx.Shat.agreers = NA
    w$optimx.Shat.convergers.agreers = NA
  }
  
  return(w)
} 


# ANALYSIS FNS ---------------------------------------------------------------


# Notes from TNE:
# In order to catch errors from individual estimation methods safely and informatively,
#  in general the estimation method fns are structured st they can be run within the
#  fn run_method_safe, which automatically runs the method within a tryCatch loop, 
#  records any error messages, and writes a results row to global var rep.res whether 
#  or not the estimation method threw an error.

# ~~ Wrapper Fn to Safely Run a Method -------

# See note at the beginning of this script
#  this fn automatically runs the method within a tryCatch loop, 
#  records any error messages, and writes a results row to global var rep.res whether 
#  or not the estimation method threw an error

# Important: this fn works if method.fn() returns multiple rows
# BUT in that case, it assumes that the CIs are shared for all rows of that method

# expects global vars: all.errors, rep.res
# directly edits res via superassignment
run_method_safe = function( method.label,
                            method.fn,
                            .rep.res ) {
  
  cat( paste("\n run_method_safe flag 1: about to try running method", method.label) )
  
  
  tryCatch({
    
    method.output = method.fn()
    new.rows = method.output$stats
    
    if ( !exists("new.rows") ) {
      cat("\n\n**** Object new.rows didn't exist for method", method.label)
      cat("\nHere is method.output:\n")
      print(method.output)
    }
    
    cat( paste("\n run_method_safe flag 2: done calling method.fn() for", method.label) )
    
    error = NA
    
  }, error = function(err) {
    # needs to be superassignment because inside the "error" fn
    error <<- err$message
    
    # only need one variable in the blank dataframe since bind_rows below
    #  will fill in the rest
    new.rows <<- data.frame( method = method.label )
    
  })
  
  new.rows = new.rows %>% add_column( method = method.label, .before = 1 )
  new.rows$overall.error = error
  
  # optimx.dataframe is itself a df, so needs to be handled differently
  # if ( !is.null(optimx.dataframe) ) new.row = bind_cols(new.row, optimx.dataframe)
  
  if ( nrow(.rep.res) == 0 ) .rep.res = new.rows else .rep.res = bind_rows(.rep.res, new.rows)
  return(.rep.res) 
  
}

# example of how to call it when method.fn takes args
# all.errors = c()
# if  exists("rep.res") ) r("rep.re("rep.re
# run_method_safe( method = "mle",
#                  method.fn = function() estimate_mles(x = x, get.CIs = TRUE ) )

# #### Sanity checks
# # fake method for estimating the moments, but it breaks if x<0 or x>5
# crappy_method = function(x) {
#   if ( x > 0 & x < 5 ) return( list(Mhat = x+1,
#                                     Vhat = x-1,
#                                     M.CI = c(NA, NA),
#                                     V.CI = c(NA, NA) ) )
#   if ( x <= 0 ) stop("Fake error A generated by method.fn!")
#   if ( x >= 5 ) stop("Fake error B generated by method.fn!")
# }
# 
# all.errors = c()
# if( exists("rep.res") ) rm(rep.res)
# 
# run_method_safe( "mle",
#                  method.fn = function() { crappy_method(-1) } )
# 
# # no error on this one
# run_method_safe( "mle",
#                  method.fn = function() { crappy_method(4) } )
# 
# 
# # this one will have a different error
# # no error on this one
# run_method_safe( "mle", method.fn = function() { crappy_method(40) } )
# 
# expect_equal( all.errors, c( "mle: Fake error A generated by method.fn!",
#                             "mle: Fake error B generated by method.fn!" ) )
# 
# expect_equal( rep.res,
#               data.frame( method = rep("mle", 3),
#                           Mhat = c(NA, 5, NA),
#                           Vhat = c(NA, 3, NA),
#                           MLo = rep(NA, 3),
#                           MHi = rep(NA, 3),
#                           VLo = rep(NA, 3),
#                           VHi = rep(NA, 3) ) )
# #### end sanity checks



# DATA SIMULATION ---------------------------------------------------------------


# - Mu: overall mean for meta-analysis (SMD if Ytype = "cont-SMD"; log-RR if Ytype = "bin-RR"; log-OR if Ytype = "log-OR")
# - t2a: across-study heterogeneity (NOT total heterogeneity)
# - true.dist: "norm" or "expo"
# - muN, minN: mean and lower limit of uniform dist from which to draw sample sizes
# - Ytype: "cont-SMD" (SMD effect size), "bin-RR" (RR effect size), "bin-OR" (OR effect size)
# - p0: P(Y=0 | X=0); only needed if Ytype is binary 
sim_meta = function(k.pub,
                    Mu,  
                    t2a,  
                    true.dist,
                    
                    # within-study parameters
                    N.expr,
                    Ytype,
                    p0) {
  
  
  # collect arguments
  .args = mget(names(formals()), sys.frame(sys.nframe()))
  # remove unnecessary args for sim_one_study_set
  .args = .args[ !names(.args) %in% c("k.pub")]
  
  #browser()

  for( i in 1:k.pub ) {
    
    newRow = do.call( sim_one_study, .args )
    
    # add study ID
    newRow = newRow %>% add_column( .before = 1,
                                      study = i )
    
    if ( i == 1 ) .dat = newRow else .dat = rbind( .dat, newRow )
  }
  
  # add more info to dataset
  .dat$mean_pY = meanNA(.dat$pY)
  .dat$mean_yi = meanNA(.dat$yi)
  
  return(.dat)
}

# test
if (FALSE) {
  # example: continuous Y
  d = sim_meta(k.pub = 10,
               
               Mu = 0.5,  
               t2a = 0.1^2,  
               true.dist = "norm",
               
               # within-study parameters
               # muN = 1000,
               # minN = 1000,
               N.expr = "round( runif(n = 1, min = 40, max = 400) )",
               Ytype = "cont-SMD",
               p0 = NA)
  
  # example: binary Y
  # evil scen 105
  d = sim_meta(  k.pub = 100,
                 t2a = 0.0001,
                 Mu = 0,
                 true.dist = "norm",
                 p0 = 0.05,
                 Ytype = "bin-OR",
                 N.expr = "40" )
  hist(d$sei)
  
  mean(d$pY)
  mean(d$pY0) # should match p0
}



# ~ Simulate a single study ----------------- 

# see sim_meta for args
sim_one_study = function( Mu,  # overall mean for meta-analysis
                          t2a,  # across-study heterogeneity
                          true.dist,
                          
                          # within-study sample size parameters
                          N.expr, 
                          
                          sd.w = 1,  # within-study SD(Y|X); only needed for cont outcome
                          
                          Ytype,  
                          p0 = NULL # P(Y | X=0); only needed for binary outcome
) {  
  
  # for testing
  if (FALSE){
    true.dist = "expo"
    Mu = 0.5
    t2a = 0
    muN = minN = 40
    Ytype = "bin-OR"
    sd.w = 1
    p0 = 0.01
  }
  
  # ~~ Mean for this study set -------------------------------------------------
  
  if( !true.dist %in% c("norm", "expo") ) stop("true.dist not recognized")
  
  if ( true.dist == "norm" ){
    mui = Mu + rnorm(mean = 0,
                     sd = sqrt(t2a),
                     n = 1)
  }
  
  if ( true.dist == "expo" ){
    # set the rate so the heterogeneity is correct
    mui = rexp( n = 1, rate = sqrt(1/t2a) )
    # now the mean is sqrt(t2a) rather than Mu
    # shift to have the correct mean (in expectation)
    mui = mui + ( Mu - sqrt(t2a))
  }
  
  # simulate total N for this study
  #if ( muN < minN ) stop("Should not have muN < minN")
  #N = round( runif( n = 1, min = minN, max = minN + 2*( muN - minN ) ) ) # draw from uniform centered on muN
  N = eval( parse( text = N.expr ) )
  muN = N  #@TEMP FOR USING N.EXPR INSTEAD OF MUN, MINN (needed for sanchecks)
  
  # ~~ Simulate individual subject data -------------------------------------------------
  
  # as in MRM helper code
  if ( Ytype == "cont-SMD" ) {
    # group assignments
    X = c( rep( 0, N/2 ), rep( 1, N/2 ) )
    
    ### Continuous Y ###
    # 2-group study of raw mean difference with means 0 and Mi in each group
    # and same SD
    Y = c( rnorm( n = N/2, mean = 0, sd = sd.w ),
           rnorm( n = N/2, mean = mui, sd = sd.w ) )
    
    # calculate ES for this study using metafor (see Viechtbauer "Conducting...", pg 10)
    ES = escalc( measure="SMD",   
                 n1i = N/2, 
                 n2i = N/2,
                 m1i = mean( Y[X==1] ),
                 m2i = mean( Y[X==0] ),
                 sd1i = sd( Y[X==1] ),
                 sd2i = sd( Y[X==0] ) ) 
    
    # only here to be consistent with binary Y case below
    pY = pY1 = pY0 = nY1 = nY1_theory = nY0 = nY0_theory = NA
  }
  
  
  # similar to Metasens helper code
  if ( Ytype %in% c("bin-RR", "bin-OR") ) {
    
    if ( is.null(p0) ) stop("Must specify p0 for binary Y")
    
    # group assignments
    X = c( rep( 0, N/2 ), rep( 1, N/2 ) )
    
    ### Binary Y; odds ratio ###
    if (Ytype == "bin-RR") {
      
      # check that args are ok
      if ( p0 * exp(mui) > 1 ) stop("Theoretical P(Y=1 | X=1) > 1. Adjust p0 or Mu.")
      if ( p0 * exp(mui) < 0 ) stop("Theoretical P(Y=1 | X=1) < 0. Adjust p0 or Mu.")
      
      linpred = log(p0) + mui*X  # mui is already on log scale
      # exp here because log-RR model
      Y = rbinom( size=1, n=N, prob=exp(linpred) ) 
      
      # sanity check to be returned
      nY0_theory = p0 * (muN/2)
      nY1_theory = p0 * exp(mui) * (muN/2)
      
      # sanity check
      if (FALSE){
        coef( glm(Y ~ X, family=binomial(link = "log")) )[["X"]]; mui
        mean(Y[X==0]); p0
        mean(Y[X==1]); p0 * exp(mui)
      }
    
    ### Binary Y; risk ratio ###
    } else if (Ytype == "bin-OR") {
      
      # no need to check that args are ok as above, since expit in [0,1]
      
      linpred = logit(p0) + mui*X 
      Y = rbinom( size=1, n=N, prob=expit(linpred) ) 
      
      # sanity check to be returned
      nY0_theory = p0 * (muN/2)
      nY1_theory = expit( logit(p0) + mui ) * (muN/2)

      # sanity check
      if (FALSE){
        coef( glm(Y ~ X, family=binomial(link = "logit")) )[["X"]]; mui
        mean(Y[X==0]); p0
        mean(Y[X==1]); expit( logit(p0) + mui )
      }
    }
    
    
    # calculate deaths (Y=1) and sample sizes in each group
    n1 = sum(X)  # number deaths among X=1
    n0 = length(X) - sum(X)
    y1 = sum(Y[X==1])  # number of deaths in Tx group
    y0 = sum(Y[X==0])  # number of deaths in control group
    
    # calculate log-RR for this study using metafor (see Viechtbauer "Conducting...", pg 10)
    if (Ytype == "bin-RR") {
      ES = escalc( measure="RR", ai=y1, bi=n1-y1, ci=y0, di=n0-y0 )  # returns on log scale
      
    } else if (Ytype == "bin-OR") {
      
      ES = escalc( measure="OR", ai=y1, bi=n1-y1, ci=y0, di=n0-y0 )  # returns on log scale
    }
    
    # sanity checks: summary stats
    pY = mean(Y)
    pY1 = mean(Y[X==1])
    pY0 = mean(Y[X==0])
    nY1 = y1
    nY0 = y0
  }
  
  yi = ES$yi
  vi = ES$vi
  sei = sqrt(vi)
  
  # ~~ One row for meta-analytic dataset  -------------------------------------------------
  d = data.frame( yi = yi, 
                  vi = vi, 
                  sei = sei, 
                  N = N,
                  
                  pY = pY,
                  pY1 = pY1,
                  pY0 = pY0,
                  
                  nY1 = nY1,
                  nY1_theory = nY1_theory,
                  
                  nY0 = nY0,
                  nY0_theory = nY0_theory)
  
  return(d)
  
}



# DATA WRANGLING ---------------------------------------------------------------

# corrObject: something returned by correct_dataset_phack
# looks for (or makes) global object, "res"
add_method_result_row = function(repRes = NA,
                                 corrObject,
                                 methName) {
  
  # newRow = bind_cols( corrObject$metaCorr,
  #                 corrObject$sanityChecks )
  #TEMP: DON'T KEEP THE SANITY CHECKS BECAUSE CORRECT_META_PHACK2 doesn't have it
  newRow = corrObject$metaCorr
  
  newRow = newRow %>% add_column(.before = 1,
                                 methName = methName )
  
  
  # "if" condition is hacky way to deal with repRes = NA case
  if ( is.null( nrow(repRes) ) ) repRes = newRow else repRes = bind_rows(repRes, newRow)
  return(repRes)
}



# quickly look at results when running doParallel locally
srr = function(rep.res) {
  
  if( "optimx.Mhat.winner" %in% names(rep.res) ) {
    cat("\n")
    print( rep.res %>% select(method, Mhat, MLo, MHi,
                              Shat,
                              optim.converged,
                              optimx.Mhat.winner,
                              optimx.Nconvergers,
                              optimx.Pagree.of.convergers.Mhat.winner) %>%
             mutate_if(is.numeric, function(x) round(x,2)) )
    cat("\n")
  } else {
    cat("\n")
    print( rep.res %>%
             mutate_if(is.numeric, function(x) round(x,2)) )
    cat("\n")
  }
}

# SMALL GENERIC HELPERS ---------------------

# take the logit of a probability, but truncate
#  to avoid infinities
truncLogit <- function(p) {
  p[p==0] = 0.001
  p[p==1] = 0.999
  log(p/(1-p))
}


logit <- function(p) {
  log(p/(1-p))
}


expit = function(x) {
  exp(x) / (1 + exp(x))
}

# calculate I^2 from t^2 and N
I2 = function(t2, N) {
  t2 / (t2 + 4/N)
}

# quick mean with NAs removed
meanNA = function(x){
  mean(x, na.rm = TRUE)
}


# check CI coverage
covers = function( truth, lo, hi ) {
  return( (lo <= truth) & (hi >= truth) )
}

# get names of dataframe containing a string
namesWith = function(pattern, dat){
  names(dat)[ grepl(pattern = pattern, x = names(dat) ) ]
}


# quick length(unique)
nuni = function(x) {
  length(unique(x))
}

# (re-)install package AND its dependencies
# useful for stupid rstan issues in which rstan itself it UTD but not its dependencies
# https://stackoverflow.com/questions/21010705/update-a-specific-r-package-and-its-dependencies
instPkgPlusDeps <- function(pkg, install = FALSE,
                            which = c("Depends", "Imports", "LinkingTo"),
                            inc.pkg = TRUE) {
  stopifnot(require("tools")) ## load tools
  ap <- available.packages() ## takes a minute on first use
  ## get dependencies for pkg recursively through all dependencies
  deps <- package_dependencies(pkg, db = ap, which = which, recursive = TRUE)
  ## the next line can generate warnings; I think these are harmless
  ## returns the Priority field. `NA` indicates not Base or Recommended
  pri <- sapply(deps[[1]], packageDescription, fields = "Priority")
  ## filter out Base & Recommended pkgs - we want the `NA` entries
  deps <- deps[[1]][is.na(pri)]
  ## install pkg too?
  if (inc.pkg) {
    deps = c(pkg, deps)
  }
  ## are we installing?
  if (install) {
    install.packages(deps)
  }
  deps ## return dependencies
}

# example
# instPkgPlusDeps("fields")


# CLUSTER FNS ---------------------------------------------------------------

# DO NOT CHANGE THE INDENTATION IN THE BELOW OR ELSE SLURM 
#  WILL SILENTLY IGNORE THE BATCH COMMANDS DUE TO EXTRA WHITESPACE!!
# DO NOT CHANGE THE INDENTATION IN THE BELOW OR ELSE SLURM 
#  WILL SILENTLY IGNORE THE BATCH COMMANDS DUE TO EXTRA WHITESPACE!!
sbatch_skeleton <- function() {
  return(
    "#!/bin/bash
#################
#set a job name  
#SBATCH --job-name=JOBNAME
#################  
#a file for job output, you can check job progress
#SBATCH --output=OUTFILE
#################
# a file for errors from the job
#SBATCH --error=ERRORFILE
#################
#time you think you need; default is one hour
#SBATCH --time=JOBTIME
#################
#quality of service; think of it as job priority
#SBATCH --qos=QUALITY
#################
#submit to both owners and normal partition
#SBATCH -p normal,owners,qsu
#################
#number of nodes you are requesting
#SBATCH --nodes=NODENUMBER
#################
#memory per node; default is 4000 MB
#SBATCH --mem=MEMPERNODE
#you could use --mem-per-cpu; they mean what we are calling cores
#################
#get emailed about job BEGIN, END, and FAIL
#SBATCH --mail-type=MAILTYPE
#################
#who to send email to; please change to your email
#SBATCH  --mail-user=USER_EMAIL
#################
#task to run per node; each node has 16 cores
#SBATCH --ntasks=TASKS_PER_NODE
#################
#SBATCH --cpus-per-task=CPUS_PER_TASK
#now run normal batch commands

ml load v8
ml load R/4.2.0
ml load jags
R -f PATH_TO_R_SCRIPT ARGS_TO_R_SCRIPT")
}



generateSbatch <- function(sbatch_params,
                           runfile_path = NA,
                           run_now = F) {
  
  #sbatch_params is a data frame with the following columns
  #jobname: string, specifies name associated with job in SLURM queue
  #outfile: string, specifies the name of the output file generated by job
  #errorfile: string, specifies the name of the error file generated by job
  #jobtime: string in hh:mm:ss format, max (maybe soft) is 48:00:00 
  #specifies the amoung of time job resources should be allocated
  #jobs still running after this amount of time will be aborted
  #quality: kind of like priority, normal works
  #node_number, integer: the number of nodes (computers w/16 cpus each) to allocate 
  #mem_per_node, integer: RAM, in MB, to allocate to each node
  #mailtype, string: ALL, BEGIN, END, FAIL: what types of events should you be notified about via email
  #user_email string: email address: email address to send notifications
  #tasks_per_node: integer, number of tasks, you should probably use 1
  #cpus_per_task: integer, 1-16, number of cpus to use, corresponds to number of available cores per task
  #path_to_r_script: path to r script on sherlock
  #args_to_r_script: arguments to pass to r script on command line
  #write_path: where to write the sbatch file
  #server_sbatch_path: where sbatch files will be stored on sherlock
  #runfile_path is a string containing a path at which to write an R script that can be used to run
  #the batch files generated by this function. 
  #if NA, no runfile will be written
  #run_now is a boolean specifying whether batch files should be run as they are generated
  
  sbatches <- list()
  if (!is.na(runfile_path)) {
    outfile_lines <- c(paste0("# Generated on ",  Sys.time()))
  }
  for (sbatch in 1:nrow(sbatch_params) ) {
    gen_batch <- sbatch_skeleton()
    #set job name
    if (is.null(sbatch_params$jobname[sbatch])) { 
      gen_batch <- gsub("JOBNAME", "unnamed", gen_batch) 
    } else { 
      gen_batch <- gsub("JOBNAME", sbatch_params$jobname[sbatch], gen_batch) 
    }
    #set outfile name
    if (is.null(sbatch_params$outfile[sbatch])) { 
      gen_batch <- gsub("OUTFILE", "unnamed", gen_batch) 
    } else { 
      gen_batch <- gsub("OUTFILE", sbatch_params$outfile[sbatch], gen_batch) 
    }
    #set errorfile name
    if (is.null(sbatch_params$errorfile[sbatch])) { 
      gen_batch <- gsub("ERRORFILE", "unnamed", gen_batch) 
    } else { 
      gen_batch <- gsub("ERRORFILE", sbatch_params$errorfile[sbatch], gen_batch) 
    }
    #set jobtime
    if (is.null(sbatch_params$jobtime[sbatch])) { 
      gen_batch <- gsub("JOBTIME", "unnamed", gen_batch) 
    } else { 
      gen_batch <- gsub("JOBTIME", sbatch_params$jobtime[sbatch], gen_batch) 
    }
    #set quality
    if (is.null(sbatch_params$quality[sbatch])) { 
      gen_batch <- gsub("QUALITY", "unnamed", gen_batch) 
    } else { 
      gen_batch <- gsub("QUALITY", sbatch_params$quality[sbatch], gen_batch) 
    }
    #set number of nodes
    if (is.null(sbatch_params$node_number[sbatch])) { 
      gen_batch <- gsub("NODENUMBER", "unnamed", gen_batch) 
    } else { 
      gen_batch <- gsub("NODENUMBER", sbatch_params$node_number[sbatch], gen_batch) 
    }
    #set memory per node
    if (is.null(sbatch_params$mem_per_node[sbatch])) { 
      gen_batch <- gsub("MEMPERNODE", "unnamed", gen_batch) 
    } else { 
      gen_batch <- gsub("MEMPERNODE", sbatch_params$mem_per_node[sbatch], gen_batch) 
    }
    #set requested mail message types
    if (is.null(sbatch_params$mailtype[sbatch])) { 
      gen_batch <- gsub("MAILTYPE", "unnamed", gen_batch) 
    } else { 
      gen_batch <- gsub("MAILTYPE", sbatch_params$mailtype[sbatch], gen_batch) 
    }
    #set email at which to receive messages
    if (is.null(sbatch_params$user_email[sbatch])) { 
      gen_batch <- gsub("USER_EMAIL", "unnamed", gen_batch) 
    } else { 
      gen_batch <- gsub("USER_EMAIL", sbatch_params$user_email[sbatch], gen_batch) 
    }
    #set tasks per node
    if (is.null(sbatch_params$tasks_per_node[sbatch])) { 
      gen_batch <- gsub("TASKS_PER_NODE", "unnamed", gen_batch) 
    } else { 
      gen_batch <- gsub("TASKS_PER_NODE", sbatch_params$tasks_per_node[sbatch], gen_batch) 
    }
    #set cpus per task
    if (is.null(sbatch_params$cpus_per_task[sbatch])) { 
      gen_batch <- gsub("CPUS_PER_TASK", "unnamed", gen_batch) 
    } else { 
      gen_batch <- gsub("CPUS_PER_TASK", sbatch_params$cpus_per_task[sbatch], gen_batch) 
    }
    #set path to r script
    if (is.null(sbatch_params$path_to_r_script[sbatch])) { 
      gen_batch <- gsub("PATH_TO_R_SCRIPT", "unnamed", gen_batch) 
    } else { 
      gen_batch <- gsub("PATH_TO_R_SCRIPT", sbatch_params$path_to_r_script[sbatch], gen_batch) 
    }
    #set args to r script
    if (is.null(sbatch_params$args_to_r_script[sbatch])) { 
      gen_batch <- gsub("ARGS_TO_R_SCRIPT", "unnamed", gen_batch) 
    } else { 
      gen_batch <- gsub("ARGS_TO_R_SCRIPT", sbatch_params$args_to_r_script[sbatch], gen_batch) 
    }
    
    #write batch file
    if (is.null(sbatch_params$write_path[sbatch])) { 
      cat(gen_batch, file = paste0("~/sbatch_generated_at_", gsub(" |:|-", "_", Sys.time()) ), append = F)
    } else { 
      cat(gen_batch, file = sbatch_params$write_path[sbatch], append = F)
    }
    
    if (!is.na(sbatch_params$server_sbatch_path[sbatch])) {
      outfile_lines <- c(outfile_lines, paste0("system(\"sbatch ", sbatch_params$server_sbatch_path[sbatch], "\")"))
    } 
    sbatches[[sbatch]] <- gen_batch
  }
  if (!is.na(runfile_path)) {
    cat(paste0(outfile_lines, collapse = "\n"), file = runfile_path)
  }
  if(run_now) { system(paste0("R -f ", runfile_path)) } 
  
  return(sbatches)
}


# looks at results files to identify sbatches that didn't write a file
# .max.sbatch.num: If not passed, defaults to largest number in actually run jobs.

sbatch_not_run = function(.results.singles.path,
                          .results.write.path,
                          .name.prefix,
                          .max.sbatch.num = NA ) {
  
  # get list of all files in folder
  all.files = list.files(.results.singles.path, full.names=TRUE)
  
  # we only want the ones whose name includes .name.prefix
  keepers = all.files[ grep( .name.prefix, all.files ) ]
  
  # extract job numbers
  sbatch.nums = as.numeric( unlist( lapply( strsplit( keepers, split = "_"), FUN = function(x) x[5] ) ) )
  
  # check for missed jobs before the max one
  if ( is.na(.max.sbatch.num) ) .max.sbatch.num = max(sbatch.nums)
  all.nums = 1 : .max.sbatch.num
  missed.nums = all.nums[ !all.nums %in% sbatch.nums ]
  
  # give info
  print( paste("The max job number is: ", max(sbatch.nums) ) )
  print( paste( "Number of jobs that weren't run: ",
                ifelse( length(missed.nums) > 0, length(missed.nums), "none" ) ) )
  
  if( length(missed.nums) > 0 ) {
    setwd(.results.write.path)
    write.csv(missed.nums, "missed_job_nums.csv")
  }
  
  return(missed.nums)
  
}

# FN: STITCH RESULTS FILES -------------------------------------

# given a folder path for results and a common beginning of file name of results files
#   written by separate workers in parallel, stitch results files together into a
#   single csv.

stitch_files = function(.results.singles.path, .results.stitched.write.path=.results.singles.path,
                        .name.prefix, .stitch.file.name="stitched_model_fit_results.csv") {
  
  # .results.singles.path = "/home/groups/manishad/MRM/sim_results/long"
  # .results.stitched.write.path = "/home/groups/manishad/MRM/sim_results/overall_stitched"
  # .name.prefix = "long_results"
  # .stitch.file.name="stitched.csv"
  
  # get list of all files in folder
  all.files = list.files(.results.singles.path, full.names=TRUE)
  
  # we only want the ones whose name includes .name.prefix
  keepers = all.files[ grep( .name.prefix, all.files ) ]
  
  # grab variable names from first file
  names = names( read.csv(keepers[1] )[-1] )
  
  # read in and rbind the keepers
  tables <- lapply( keepers, function(x) read.csv(x, header= TRUE) )
  s <- do.call(rbind, tables)
  
  names(s) = names( read.csv(keepers[1], header= TRUE) )
  
  if( is.na(s[1,1]) ) s = s[-1,]  # delete annoying NA row
  write.csv(s, paste(.results.stitched.write.path, .stitch.file.name, sep="/") )
  return(s)
}


# quickly look at results from job #1

res1 = function() {
  setwd("/home/groups/manishad/SAPH/long_results")
  rep.res = fread("long_results_job_1_.csv")
  srr()
  
  cat("\nErrors by method:" )
  print( rep.res %>% group_by(method) %>%
           summarise(prop.error = mean( overall.error != "" ) ) )
  
  #table(rep.res$overall.error)
  
  cat("\n\nDim:", dim(rep.res))
  cat("\n\nReps completed:", nrow(rep.res)/nuni(rep.res$method))
}




# TRASH -------------------------------------

# # 2022-3-19: I NO LONGER THINK IT MAKES SENSE TO USE THIS PLOT WITH VARIABLE SEI'S. 
# #  BUT SAVE IT B/C USEFUL WHEN SEI'S ARE THE SAME.
# #  Since the sei's differ, so do the cutpoints, and even if we plot the full normal dist, 
# #  its height doesn't align properly with the density when the density is truncated,
# #  giving the impression of a poor fit.
# # Args:
# #  - d: dataset with var names "yi", "vi", "affirm"
# #  - Mhat: RTMA estimate of underlying mean effect size
# #  - Shat: Same for heterogeneity
# #  - showAffirms: should it show all studies, even affirms?
# 
# # Returned plot:
# #  - black line = LOESS density of nonaffirms
# #  - red line = MLE from RTMA (parametric counterpart to the above)
# #  - blue line = LOESS density of all tstats (including affirms)
# 
# # IMPORTANT:
# # Note that truncation point shown in plots is only approximate because it’s set to 1.96 for all Zi.tilde, whereas actually each Zi.tilde has its own trunc point (see Zi_tilde_cdf).
# plot_trunc_densities_RTMA = function(d,
#                                      Mhat,
#                                      Shat,
#                                      showAffirms = FALSE) {
#   
#   # #TEST ONLY
#   # d = dp
#   # showAffirms = FALSE
#   # Mhat = 0.45  # FAKE for now
#   # Shat = 0.10
#   
#   # add Z-scores
#   # these are standardized ACROSS studies using the ESTIMATED mean and heterogeneity
#   # so they should look truncated N(0,1)
#   d$Zi.tilde = (d$yi - Mhat) / sqrt(Shat^2 + d$vi)
#   
#   # already has affirmative indicator
#   dn = d %>% filter(affirm == FALSE)
#   
#   
#   xmin = floor(min(dn$Zi.tilde))
#   xmax = ceiling(max(dn$Zi.tilde))
#   
#   p = ggplot(data = data.frame(x = c(xmin, 3)),
#              aes(x)) +
#     
#     geom_vline(xintercept = 0,
#                lwd = 1,
#                color = "gray") +
#     
#     # estimated density of estimates
#     geom_density( data = dn,
#                   aes(x = Zi.tilde),
#                   adjust = .3 ) +
#     
#     # estimated density from meta-analysis
#     # stat_function( fun = dnorm,
#     #                n = 101,
#     #                args = list( 
#     #                             mean = 0,
#     #                             sd = 1 ),
#     #                #aes(y = .25 * ..count..),  # doesn't work
#     #                lwd = 1.2,
#     #                color = "red") +
#     stat_function( fun = dtrunc,
#                    n = 101,
#                    args = list( spec = "norm",
#                                 mean = 0,
#                                 sd = 1,
#                                 b = qnorm(.975) ),
#                    #aes(y = .25 * ..count..),  # doesn't work
#                    lwd = 1.2,
#                    color = "red") +
#     #     
#     
#     ylab("") +
#     #scale_x_continuous( breaks = seq(xmin, 3, 0.5)) +
#     xlab("Across-study Z-score using (Mhat, Shat)") +
#     theme_minimal() +
#     scale_y_continuous(breaks = NULL) +
#     theme(text = element_text(size=16),
#           axis.text.x = element_text(size=16))
#   
#   
#   # also show density of all t-stats, not just the nonaffirms
#   if ( showAffirms == TRUE ) {
#     p = p + geom_density( data = d,
#                           aes(x = Zi.tilde),
#                           color = "blue",
#                           adjust = .3 )
#   }
#   
#   
#   return(p)
#   
# }
# # OLDER VERSION - MAYBE SAVE?
# # plot empirical data
# 
# # .obj: object returned by correct_meta_phack2
# # showAffirms: should it show all studies, even affirms?
# # black line = LOESS density of nonaffirms
# # red line = MLE from RTMA (parametric counterpart to the above)
# # blue line = LOESS density of all tstats (including affirms)
# plot_trunc_densities = function(.obj,
#                                 showAffirms = FALSE) {
#   
#   # already has affirmative indicator
#   d = .obj$data
#   dn = d[d$affirm == FALSE,]
#   
#   tstatMeanMLE = .obj$sanityChecks$tstatMeanMLE
#   tstatVarMLE = .obj$sanityChecks$tstatVarMLE
#   
#   xmin = floor(min(dn$tstat))
#   xmax = ceiling(max(dn$tstat))
#   
#   p = ggplot(data = data.frame(x = c(xmin, 3)),
#              aes(x)) +
#     
#     geom_vline(xintercept = 0,
#                lwd = 1,
#                color = "gray") +
#     
#     # geom_vline(xintercept = tstatMeanMLE,
#     #            lty = 2,
#     #            lwd = 1,
#     #            color = "red") +
#     
#     
#     # estimated density of estimates
#     geom_density( data = dn,
#                   aes(x = tstat),
#                   adjust = .3 ) + 
#     
#     
#     
#     # estimated density from meta-analysis
#     stat_function( fun = dtrunc,
#                    n = 101,
#                    args = list( spec = "norm",
#                                 mean = tstatMeanMLE,
#                                 sd = sqrt(tstatVarMLE),
#                                 b = .obj$crit),
#                    #aes(y = .25 * ..count..),  # doesn't work
#                    lwd = 1.2,
#                    color = "red") +
#     
#     
#     
#     ylab("") +
#     #scale_x_continuous( breaks = seq(xmin, 3, 0.5)) +
#     xlab("t-stat") +
#     theme_minimal() +
#     scale_y_continuous(breaks = NULL) +
#     theme(text = element_text(size=16),
#           axis.text.x = element_text(size=16))
#   
#   
#   # also show density of all t-stats, not just the nonaffirms
#   if ( showAffirms == TRUE ) {
#     p = p + geom_density( data = d,
#                           aes(x = tstat),
#                           color = "blue",
#                           adjust = .3 )
#   }
#   
#   
#   return(p)
#   
# }

# OLD VERSION (for sanity-checking the one below)

# 2022-3-12
# nonaffirms only
# Zi_tilde_cdf_OLD = function(x, .SE, .Shat) {
#   
#   # calculate cutpoint for EACH Zi.tilde
#   # **reasoning:
#   #  we observe a truncated sample st yi > 1.96*SE
#   # therefore Zi.tilde = (yi - mu) / ( sqrt(Shat^2 + SE^2) ) > (1.96*SE - mu) / ( sqrt(Shat^2 + SE^2) )
#   # and we know that Zi.tilde ~ N(0,1) prior to truncation
#   Zi.tilde.crit = ( qnorm(.975) * .SE ) / sqrt(.Shat^2 + .SE^2)
#   
#   ptruncnorm(q = x,
#              a = -Inf,
#              b = Zi.tilde.crit,
#              mean = 0,
#              sd = 1)
# }

# # 2022-3-12
# # fit diagnostics
# # get CDF of (non-iid) marginal Z-scores (Zi.tilde)
# #  given a fitted Shat
# # .affirm: VECTOR with same length as x for affirm status
# #  including the affirms is useful for 2PSM
# Zi_tilde_cdf = function(.Zi.tilde,
#                         .SE,
#                         .Shat,
#                         .affirm) {
#   
#   
#   #if ( length(.Zi.tilde) > 1 ) stop("Length of .Zi.tilde must be 1")
#   if ( length(.affirm) != length(.Zi.tilde) ) stop(".affirm must have same length as x")
#   
#   # calculate cutpoint for this Zi.tilde
#   # **reasoning:
#   #  we observe a truncated sample st yi > 1.96*SE
#   # therefore Zi.tilde = (yi - mu) / ( sqrt(Shat^2 + .SE^2) ) > (1.96*SE - mu) / ( sqrt(Shat^2 + .SE^2) )
#   # and we know that Zi.tilde ~ N(0,1) prior to truncation
#   Zi.tilde.crit = ( qnorm(.975) * .SE ) / sqrt(.Shat^2 + .SE^2)
#   
#   
#   dat = data.frame(Zi.Tilde = .Zi.tilde,
#                    Zi.Tilde.Crit = Zi.tilde.crit,
#                    Affirm = .affirm)
#   
#   # if ( .affirm == FALSE ) expect_equal( .Zi.tilde < Zi.tilde.crit, TRUE )
#   # if ( .affirm == FALSE ) expect_equal( .Zi.tilde < Zi.tilde.crit, TRUE )
#   
#   # if ( Affirm == FALSE ) {
#   #   return( ptruncnorm(q = Zi.Tilde,
#   #                      a = -Inf,
#   #                      b = Zi.Tilde.Crit,
#   #                      mean = 0,
#   #                      sd = 1) )
#   # } else if ( Affirm == TRUE ) {
#   #   return( ptruncnorm(q = Zi.Tilde,
#   #                      a = Zi.Tilde.Crit,
#   #                      b = Inf,
#   #                      mean = 0,
#   #                      sd = 1) )
#   # }
#   
#   dat$cdfi = NA
#   
#   if ( any(dat$Affirm == FALSE) ) {
#     dat$cdfi[ dat$Affirm == FALSE ] = ptruncnorm(q = dat$Zi.Tilde[ dat$Affirm == FALSE ],
#                                                  a = -Inf,
#                                                  b = dat$Zi.Tilde.Crit[ dat$Affirm == FALSE ],
#                                                  mean = 0,
#                                                  sd = 1)
#   }
#   
#   if ( any(dat$Affirm == TRUE) ) {
#     dat$cdfi[ dat$Affirm == TRUE ] = ptruncnorm(q = dat$Zi.Tilde[ dat$Affirm == TRUE ],
#                                                 a = dat$Zi.Tilde.Crit[ dat$Affirm == TRUE ],
#                                                 b = Inf,
#                                                 mean = 0,
#                                                 sd = 1)
#   }
#   
#   return(dat$cdfi)
#   
#   
#   # 
#   # dat = dat %>% rowwise() %>%
#   #   mutate( cdfi = function(Zi.Tilde, Affirm){
#   #     if ( Affirm == FALSE ) {
#   #       return( ptruncnorm(q = Zi.Tilde,
#   #                          a = -Inf,
#   #                          b = Zi.Tilde.Crit,
#   #                          mean = 0,
#   #                          sd = 1) )
#   #     } else if ( Affirm == TRUE ) {
#   #       return( ptruncnorm(q = Zi.Tilde,
#   #                          a = Zi.Tilde.Crit,
#   #                          b = Inf,
#   #                          mean = 0,
#   #                          sd = 1) )
#   #     }
#   #     
#   #     
#   #   } )
#   
# }


# ### Sanity checks ###
# 
# # sanity-check nonaffirmatives (taken from Kvarven's Belle meta-analysis)
# Zi.tilde = c(-1.19493855966001, -0.330782431293096, 0.18135714022493, -0.355495378590117)
# sei = c(0.183, 0.169, 0.269999994444444, 0.222)
# Shat = 0.23  # from 2PSM
# cdfi.mine = ptruncnorm(q = Zi.tilde,
#                            a = -Inf,
#                            b = ( qnorm(.975) * sei ) / sqrt(Shat^2 + sei^2),
#                            mean = 0,
#                            sd = 1)
# 
# 
# cdfi = Zi_tilde_cdf(.Zi.tilde = Zi.tilde,
#                     .SE = sei,
#                     .Shat = Shat,
#                     .affirm = rep(FALSE, length(Zi.tilde)))
# expect_equal( cdfi.mine, cdfi)
# 
# # sanity-check affirmatives (also from Belle)
# Zi.tilde = c(2.10452951219512, 2.28191480814403, 4.29142844422158, 2.33576635036496, 
#              2.707112988856, 2.86885247443619, 3.10843366971151, 2.37931028904956, 
#              2.6049383759669, 2.06417123549919, 2.031579)
# sei = c(0.287, 0.188000002659574, 0.175000002857143, 0.137, 0.23900000209205, 
#         0.304999998360656, 0.249000002008032, 0.261000001915709, 0.242999997942387, 
#         0.186999994652406, 0.19)
# Shat = 0.23
# cdfi.mine = ptruncnorm(q = Zi.tilde,
#                        a = ( qnorm(.975) * sei ) / sqrt(Shat^2 + sei^2),
#                        b = Inf,
#                        mean = 0,
#                        sd = 1)
# 
# 
# cdfi = Zi_tilde_cdf(.Zi.tilde = Zi.tilde,
#                     .SE = sei,
#                     .Shat = Shat,
#                     .affirm = rep(TRUE, length(Zi.tilde)))
# expect_equal( cdfi.mine, cdfi)



# # 2022-3-12
# # test for RTMA fit
# # yi, sei: can be for all published studies or for just nonaffirms or just affirms
# #  including the affirms is useful for 2PSM but not RTMA
# 
# my_ks_test_RTMA = function(yi,
#                            sei,
#                            Mhat,
#                            Shat) {
#   
#   if ( is.na(Mhat) | is.na(Shat) ) return(NA)
#   
#   affirm = (yi/sei) > qnorm(.975)
#   
#   # retain only nonaffirmatives
#   # nonaffirm = yi/sei < qnorm(.975)
#   # yi = yi[ nonaffirm == TRUE ]
#   # sei = sei[ nonaffirm == TRUE ]
#   
#   Zi.tilde = (yi - Mhat) / sqrt(Shat^2 + sei^2)
#   
#   res = ks.test(Zi.tilde, function(x) Zi_tilde_cdf(.Zi.tilde = Zi.tilde,
#                                                    .SE = sei,
#                                                    .Shat = Shat,
#                                                    .affirm = affirm) )
#   res$p.value
#   
# }




# # 2022-3-12
# # test for RTMA fit
# # yi, sei: can be for all published studies or for just nonaffirms; 
# #  fn will automatically retain only nonaffirms
# my_ks_test_RTMA = function(yi,
#                            sei,
#                            Mhat,
#                            Shat) {
#   
#   # retain only nonaffirmatives
#   nonaffirm = yi/sei < qnorm(.975)
#   yi = yi[ nonaffirm == TRUE ]
#   sei = sei[ nonaffirm == TRUE ]
#   
#   Zi.tilde = (yi - Mhat) / sqrt(Shat^2 + sei^2)
#   
#   res = ks.test(Zi.tilde, function(x) Zi_tilde_cdf(x,
#                                                    .SE = sei,
#                                                    .Shat = Shat) )
#   res$p.value
#   
# }
