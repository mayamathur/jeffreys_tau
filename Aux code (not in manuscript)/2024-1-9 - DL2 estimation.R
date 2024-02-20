
library(metafor)
library(tibble)

# try the example given here:
# https://sammancuso.com/2020/02/12/meta-analysis-with-the-2-step-dersimonian-and-laird-estimator/

d <- tibble(
  study = c("Baum", "Anand", "CLASP", "Bradley", "Lovet", "Teo"),
  n1i = c(156, 303, 565, 1570, 103, 4659),
  n2i = c(74, 303, 477, 1565, 105, 4650),
  ai = c(5, 5, 12, 69, 9, 313),
  ci = c(8, 17, 9, 94, 11, 352),
  yi = c(-1.2976, -1.2649, 0.1208, -0.3294, -0.2007, -0.1285),
  vi = c(0.3468, 0.2657, 0.1984, 0.0265, 0.2233, 0.0065)
)

# first, regular DL:
rma.uni(yi = yi,
        vi = vi, 
        data = d,
        method = "DL",
        knha = TRUE)


rma.uni(yi = yi,
        vi = vi, 
        data = d,
        method = "DLIT",
        knha = TRUE)


# c.f. function from that website:
mw_est <- function(yi,
                   vi,
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

mw_est(
  yi = yi,
  vi = vi,
  data = dat_clasp,
  lab = study,
  method = "DL2",
  pi = FALSE,
  Q.profile = TRUE,
  hksj = FALSE
)

# they don't agree. :(