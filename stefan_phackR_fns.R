
# This code is from a repo I forked (with no edits) from Stefan's on 2023-05-06.
#   my fork: https://github.com/mayamathur/phacking_compendium
#   Stefan: astefan1/phacking_compendium


# SUBGROUP SELECTIVE REPORTING  -------------------------------------------------

# https://github.com/mayamathur/phacking_compendium/blob/master/phackR/R/subgroupAnalysis.R

# ==============================================================================
# Subgroup Analyses
# ==============================================================================

#' Simulate data with subgroups
#' @description Outputs data frame with multiple binary variables from which subgroups can be extracted
#' @param nobs.group Vector giving number of observations per group
#' @param nsubvars Integer specifying number of variables for potential subgroups

.sim.subgroup <- function(nobs.group, nsubvars){
  
  dat <- .sim.data(nobs.group)
  
  # Observations per group and total observations
  if(length(nobs.group) == 1) nobs.group <- rep(nobs.group, 2)
  nobs <- sum(nobs.group)
  
  subvars <- matrix(NA, nrow = nobs, ncol = nsubvars)
  for(i in 1:nsubvars){
    subvars[,i] <- sample(c(0, 1), size = nobs, replace = TRUE)
  }
  
  res <- cbind(dat, subvars)
  
  return(res)
  
}

#' P-Hacking function for multiple subgroups analysis
#' @description Outputs a p-hacked p-value and a vector of all p-values that were computed in the process
#' @param df A matrix or data frame containing all relevant data
#' @param iv Integer specifying the location of the binary independent variable in the data frame
#' @param dv Integer specifying the location of the dependent variable in the data frame
#' @param subvars Vector specifying the location of the subgroup variables in the data frame
#' @param alternative Direction of the t-test ("two.sided", "less", "greater")
#' @param strategy String value: One out of "firstsig", "smallest", "smallest.sig"
#' @param alpha Significance level of the t-test
#' @importFrom dplyr group_by_at do
#' @importFrom stats t.test
#' @importFrom dplyr "%>%"
#' @importFrom rlang .data

.subgroupHack <- function(df, iv, dv, subvars, alternative = "two.sided", strategy = "firstsig", alpha = 0.05){
  
  # Prepare data frame
  ttest.df <- cbind(df[,iv], df[,dv])
  subvars.df <- cbind(df[, subvars])
  dfnew <- as.data.frame(cbind(ttest.df, subvars.df))
  
  # Compute p-values, R^2, Cohen's d
  
  # Not p-hacked
  mod.orig <- stats::t.test(ttest.df[,2] ~ ttest.df[,1], var.equal = TRUE, alternative = alternative)
  p.orig <- mod.orig$p.value
  r2.orig <- .compR2t(ttest.df[ttest.df[,1] == unique(ttest.df[,1])[1],2],
                      ttest.df[ttest.df[,1] == unique(ttest.df[,1])[2],2])
  d.orig <- .compCohensD(unname(mod.orig$statistic), nrow(ttest.df)/2)
  
  
  # p-hacked
  ps <- list()
  ds <- list()
  r2s <- list()
  
  for(i in 1:length(subvars)){
    
    tmp <- dplyr::group_by_at(dfnew, subvars[i]) %>%
      dplyr::do(as.data.frame(stats::t.test(.data$V2 ~ .data$V1, var.equal = TRUE, alternative = alternative)[c("p.value", "statistic")]))
    tmp2 <- dplyr::group_by_at(dfnew, subvars[i]) %>%
      dplyr::do(as.data.frame(table(.data$V1)))
    tmp3 <- dplyr::group_by_at(dfnew, subvars[i]) %>% do(as.data.frame(.compR2t(.data$V2[.data$V1 == unique(.data$V1)[1]], .data$V2[.data$V1 == unique(.data$V1)[2]])))
    
    ps[[i]] <- tmp[[2]]
    ds[[i]] <- c(tmp[[3]][1]*sqrt(sum(1/tmp2[[3]][1:2])), tmp[[3]][2]*sqrt(sum(1/tmp2[[3]][3:4])))
    r2s[[i]] <- tmp3[[2]]
    
  }
  
  ps <- c(p.orig, unlist(ps))
  r2s <- c(r2.orig, unlist(r2s))
  ds <- c(d.orig, unlist(ds))
  
  # Select final p-hacked p-value based on strategy
  p.final <- .selectpvalue(ps = ps, strategy = strategy, alpha = alpha)
  r2.final <- unique(r2s[ps == p.final])
  d.final <- unique(ds[ps == p.final])
  
  return(list(p.final = p.final,
              ps = ps,
              r2.final = r2.final,
              r2s = r2s,
              d.final = d.final,
              ds = ds))
  
}

#' Simulate p-hacking with multiple subgroups
#' Outputs a matrix containing the p-hacked p-values (\code{ps.hack}) and the original p-values (\code{ps.orig}) from all iterations
#' @param nobs.group Vector giving number of observations per group
#' @param nsubvars Integer specifying number of variables for potential subgroups
#' @param alternative Direction of the t-test ("two.sided", "less", "greater")
#' @param strategy String value: One out of "firstsig", "smallest", "smallest.sig"
#' @param alpha Significance level of the t-test
#' @param iter Number of simulation iterations
#' @param shinyEnv Is the function run in a Shiny session? TRUE/FALSE
#' @export

sim.subgroupHack <- function(nobs.group, nsubvars, alternative = "two.sided", strategy = "firstsig", alpha = 0.05, iter = 1000, shinyEnv = FALSE){
  
  # Simulate as many datasets as desired iterations
  dat <- list()
  for(i in 1:iter){
    dat[[i]] <- .sim.subgroup(nobs.group = nobs.group, nsubvars = nsubvars)
  }
  
  # Apply p-hacking procedure to each dataset
  .subgroupHackList <- function(x){
    .subgroupHack(df = x, iv = 1, dv = 2, subvars = c(3:(2+nsubvars)),
                  alternative = alternative, strategy = strategy, alpha = alpha)
  }
  
  if(!shinyEnv){
    res <- pbapply::pblapply(dat, .subgroupHackList)
  }
  
  if(shinyEnv){
    percentage <- 0
    withProgress(message = "Running simulation", value = 0, {
      res = lapply(dat, function(x){
        percentage <<- percentage + 1/length(dat)*100
        incProgress(1/length(dat), detail = paste0("Progress: ",round(percentage,2), "%"))
        .subgroupHack(df = x, iv = 1, dv = 2, subvars = c(3:(2+nsubvars)),
                      alternative = alternative, strategy = strategy, alpha = alpha)
      })
    })
  }
  
  ps.hack <- NULL
  ps.orig <- NULL
  r2s.hack <- NULL
  r2s.orig <- NULL
  ds.hack <- NULL
  ds.orig <- NULL
  
  for(i in 1:iter){
    ps.hack[i] <- res[[i]][["p.final"]]
    ps.orig[i] <- res[[i]][["ps"]][1]
    r2s.hack[i] <- res[[i]][["r2.final"]]
    r2s.orig[i] <- res[[i]][["r2s"]][1]
    ds.hack[i] <- res[[i]][["d.final"]]
    ds.orig[i] <- res[[i]][["ds"]][1]
  }
  
  res <- cbind(ps.hack, ps.orig, r2s.hack, r2s.orig, ds.hack, ds.orig)
  
  return(res)
  
}

# DV SELECTIVE REPORTING  -------------------------------------------------

# https://github.com/mayamathur/phacking_compendium/blob/master/phackR/R/selectiveReportingDV.R

# ==============================================================================
# Selective Reporting of the Dependent Variable
# ==============================================================================

#' Simulate dataset with multiple dependent variables
#' @description Outputs data frame with a grouping variable and multiple correlated dependent variables
#' @param nobs.group Vector giving number of observations per group
#' @param nvar Number of dependent variables in the data frame
#' @param r Desired correlation between the dependent variables (scalar)

.sim.multDV <- function(nobs.group, nvar, r){
  
  # Observations per group
  if(length(nobs.group) == 1) nobs.group <- rep(nobs.group, 2)
  
  # Generate group vector
  group <- rep(1:length(nobs.group), nobs.group)
  
  # Generate dependent variables
  dvs <- .sim.multcor(nobs = sum(nobs.group), nvar = nvar, r = r)
  
  # Generate data frame
  res <- cbind(group, dvs)
  
  return(res)
}

#' P-Hacking function for multiple dependent variables
#' @description Outputs a p-hacked p-value and a vector of all p-values that were computed in the process
#' @param df Data frame with one group variable and multiple dependent variables
#' @param dvs Vector defining the DV columns (will be checked in given order)
#' @param group Scalar defining grouping column
#' @param strategy String value: One out of "firstsig", "smallest", "smallest.sig"
#' @param alternative Direction of the t-test ("two.sided", "less", "greater")
#' @param alpha Significance level of the t-test
#' @importFrom stats t.test

.multDVhack <- function(df, dvs, group, strategy = "firstsig", alternative = "two.sided", alpha = 0.05){
  
  # Prepare data frame
  dvs <- as.matrix(df[, dvs], ncol = length(dvs))
  group <- df[, group]
  mod <- list()
  r2s <- NULL
  
  # Compute t-tests
  for(i in 1:ncol(dvs)){
    
    mod[[i]] <- stats::t.test(dvs[, i] ~ group,
                              var.equal = TRUE, alternative = alternative)
    r2s[i] <- .compR2t(dvs[group == unique(group)[1], i],
                       dvs[group == unique(group)[2], i])
  }
  
  ps <- unlist(simplify2array(mod)["p.value", ])
  ds <- .compCohensD(unlist(simplify2array(mod)["statistic", ]), length(df[, group])/2)
  
  # Select final p-hacked p-value based on strategy
  p.final <- .selectpvalue(ps = ps, strategy = strategy, alpha = alpha)
  r2.final <- unique(r2s[ps == p.final])
  d.final <- unique(ds[ps == p.final])
  
  return(list(p.final = p.final,
              ps = ps,
              r2.final = r2.final,
              r2s = r2s,
              d.final = d.final,
              ds = ds))
  
}

#' Simulate p-Hacking with multiple dependent variables
#' @description Outputs a matrix containing the p-hacked p-values (\code{ps.hack}) and the original p-values (\code{ps.orig}) from all iterations
#' @param nobs.group Vector giving number of observations per group
#' @param nvar Number of dependent variables (columns) in the data frame
#' @param r Desired correlation between the dependent variables (scalar)
#' @param strategy String value: One out of "firstsig", "smallest", "smallest.sig"
#' @param iter Number of simulation iterations
#' @param alternative Direction of the t-test ("two.sided", "less", "greater")
#' @param alpha Significance level of the t-test (default: 0.05)
#' @param shinyEnv Is the function run in a Shiny session? TRUE/FALSE
#' @export

sim.multDVhack <- function(nobs.group, nvar, r, strategy = "firstsig", iter = 1000, alternative = "two.sided", alpha = 0.05, shinyEnv = FALSE){
  
  # Simulate as many datasets as desired iterations
  dat <- list()
  for(i in 1:iter){
    dat[[i]] <- .sim.multDV(nobs.group = nobs.group, nvar = nvar, r = r)
  }
  
  # Apply p-hacking procedure to each dataset
  
  if(!shinyEnv){
    res <- pbapply::pblapply(dat, .multDVhack, dvs = c(2:(nvar+1)), group = 1,
                             strategy = strategy, alternative = alternative, alpha = alpha)
  }
  
  if(shinyEnv){
    percentage <- 0
    withProgress(message = "Running simulation", value = 0, {
      res = lapply(dat, function(x){
        percentage <<- percentage + 1/length(dat)*100
        incProgress(1/length(dat), detail = paste0("Progress: ",round(percentage,2), "%"))
        .multDVhack(df=x, dvs = c(2:(nvar+1)), group = 1,
                    strategy = strategy, alternative = alternative, alpha = alpha)
      })
    })
  }
  
  ps.hack <- NULL
  ps.orig <- NULL
  r2s.hack <- NULL
  r2s.orig <- NULL
  ds.hack <- NULL
  ds.orig <- NULL
  
  for(i in 1:iter){
    ps.hack[i] <- res[[i]][["p.final"]]
    ps.orig[i] <- res[[i]][["ps"]][1]
    r2s.hack[i] <- res[[i]][["r2.final"]]
    r2s.orig[i] <- res[[i]][["r2s"]][1]
    ds.hack[i] <- res[[i]][["d.final"]]
    ds.orig[i] <- res[[i]][["ds"]][1]
  }
  
  res <- cbind(ps.hack, ps.orig, r2s.hack, r2s.orig, ds.hack, ds.orig)
  
  return(res)
}

# OPTIONAL STOPPING  -------------------------------------------------

# https://github.com/mayamathur/phacking_compendium/blob/master/phackR/R/optionalStopping.R

# ==============================================================================
# Optional Stopping Based on Significance
# ==============================================================================

# Generic sampling function .sim.data() can be used

#' Optional Stopping based on existing dataset
#' @description Returns a p-hacked p-value and a non-p-hacked p-value based on the maximum sample size
#' @param df Data frame
#' @param group group Scalar defining grouping column
#' @param dv Scalar defining location of dependent variable in the data frame
#' @param n.min Minimum sample size
#' @param n.max Maximum sample size
#' @param step Step size of the optional stopping (default is 1)
#' @param peek Determines how often one peeks at the data. Overrides step argument if not NULL.
#' @param alternative Direction of the t-test ("two.sided", "less", "greater")
#' @param alpha Significance level of the t-test (default: 0.05)
#' @importFrom stats t.test
#' @importFrom utils tail

.optstop <- function(df, group, dv, n.min, n.max, step = 1, peek = NULL, alternative = "two.sided", alpha = 0.05){
  
  # Extract group variables
  g1 <- df[df[,group] == unique(df[,group])[1], dv]
  g2 <- df[df[,group] == unique(df[,group])[2], dv]
  
  # Sanity check: Enough data?
  stopifnot(length(g1) >= n.max && length(g2) >= n.max)
  
  # Determine places of peeks
  if(is.null(peek)){
    peeks <- seq(n.min, n.max, by=step)
    if(step > (n.max-n.min)) peeks <- c(n.min, n.max)
  } else {
    peeks <- round(seq(n.min, n.max, length.out = peek))
  }
  
  # Compute t-tests
  mod <- sapply(peeks, FUN = function(x) {stats::t.test(g1[1:x], g2[1:x], var.equal = TRUE, alternative = alternative)})
  ps <- simplify2array(mod["p.value",])
  r2s <- sapply(peeks, FUN = function(x) {.compR2t(g1[1:x], g2[1:x])})
  ds <- .compCohensD(simplify2array(mod["statistic",]), peeks)
  
  # Do the p-hacking
  if(any(ps < alpha) == FALSE){
    p.final <- utils::tail(ps, 1)
    r2.final <- utils::tail(r2s, 1)
    d.final <- utils::tail(ds, 1)
  } else if (any(ps < alpha) == TRUE) {
    p.final <- ps[which(ps < alpha)][1]
    r2.final <- unique(r2s[ps == p.final])
    d.final <- unique(ds[ps == p.final])
  }
  
  return(list(p.final = p.final,
              ps = ps,
              r2.final = r2.final,
              r2s = r2s,
              d.final = d.final,
              ds = ds))
}

#' Simulate p-hacking with incorrect rounding
#' @param n.min Minimum sample size
#' @param n.max Maximum sample size
#' @param step Step size of the optional stopping (default is 1)
#' @param peek Determines how often one peeks at the data. Overrides step argument if not NULL.
#' @param alternative Direction of the t-test ("two.sided", "less", "greater")
#' @param iter Number of iterations
#' @param alpha Significance level of the t-test (default: 0.05)
#' @param shinyEnv Is the function run in a Shiny session? TRUE/FALSE
#' @importFrom utils tail
#' @export
#'

sim.optstop <- function(n.min, n.max, step = 1, peek = NULL, alternative = "two.sided", iter = 1000, alpha = 0.05, shinyEnv = FALSE){
  
  # Simulate as many datasets as desired iterations
  dat <- list()
  for(i in 1:iter){
    dat[[i]] <- .sim.data(nobs.group = n.max)
  }
  
  # Apply p-hacking procedure to each dataset
  if(!shinyEnv){
    res <- pbapply::pblapply(dat, .optstop, group = 1, dv = 2,
                             n.min = n.min, n.max = n.max, step = step, peek = peek,
                             alternative = alternative, alpha = alpha)
  }
  
  if(shinyEnv){
    percentage <- 0
    withProgress(message = "Running simulation", value = 0, {
      res = lapply(dat, function(x){
        percentage <<- percentage + 1/length(dat)*100
        incProgress(1/length(dat), detail = paste0("Progress: ",round(percentage,2), "%"))
        .optstop(df=x, group = 1, dv = 2,
                 n.min = n.min, n.max = n.max, step = step,
                 alternative = alternative, alpha = alpha)
      })
    })
  }
  
  ps.hack <- NULL
  ps.orig <- NULL
  r2s.hack <- NULL
  r2s.orig <- NULL
  ds.hack <- NULL
  ds.orig <- NULL
  
  for(i in 1:iter){
    ps.hack[i] <- res[[i]][["p.final"]]
    ps.orig[i] <- utils::tail(res[[i]][["ps"]], 1)
    r2s.hack[i] <- res[[i]][["r2.final"]]
    r2s.orig[i] <- utils::tail(res[[i]][["r2s"]], 1)
    ds.hack[i] <- res[[i]][["d.final"]]
    ds.orig[i] <- utils::tail(res[[i]][["ds"]], 1)
  }
  
  res <- cbind(ps.hack, ps.orig, r2s.hack, r2s.orig, ds.hack, ds.orig)
  
  return(res)
  
}


# HELPERS.R  -------------------------------------------------

# https://github.com/mayamathur/phacking_compendium/blob/master/phackR/R/helpers.R


# ==============================================================================
# Helpers
# ==============================================================================

#' Simulate multivariate correlated data for continuous variables
#' @description Outputs a data frame with correlated variables of defined length
#' @param nobs Number of observations (rows) in the simulated data frame
#' @param nvar Number of variables (columns) in the data frame
#' @param r Desired correlation between the variables (integer)
#' @param mu Mean of the random data
#' @param sd Standard deviation of the random data
#' @param missing Proportion of missing values per variable (e.g., 0.2 = 20 percent)
#' @importFrom stats rnorm

.sim.multcor <- function(nobs, nvar, r, mu = 0, sd = 1, missing = 0){
  
  # set up correlation matrix
  R <- matrix(rep(r, nvar**2), nrow = nvar)
  diag(R) <- rep(1, nvar)
  
  # transposed Cholesky decomposition of correlation matrix
  U <- t(chol(R))
  
  # create random noise matrix
  random.normal <- matrix(stats::rnorm(nvar*nobs, mu, sd), nrow=nvar, ncol=nobs)
  
  # create raw data from matrix multiplication of U and random noise
  X <- as.data.frame(t(U %*% random.normal))
  
  # add missing values
  if(missing > 0){
    if(missing * nobs < 2){
      navalues <- as.data.frame(t(replicate(nvar, sample(1:nobs, missing*nobs))))
    } else {
      navalues <- as.data.frame(replicate(nvar, sample(1:nobs, missing*nobs)))
    }
    for(i in 1:nvar){
      X[unlist(navalues[,i]),i] <- NA
    }
  }
  
  return(X)
  
}

#' Generic sampling function
#' @description Outputs a data frame with two columns
#' @param nobs.group Number of observations per group. Either a scalar or a vector with two elements.
#' @importFrom stats rnorm

.sim.data <- function(nobs.group){
  
  if(length(nobs.group) == 1) nobs.group <- rep(nobs.group, 2)
  V1 <- stats::rnorm(nobs.group[1], 0, 1)
  V2 <- stats::rnorm(nobs.group[2], 0, 1)
  group <- c(rep(1, nobs.group[1]), rep(2, nobs.group[2]))
  
  res <- cbind(group, c(V1, V2))
  return(res)
  
}

#' Create data frames without outliers
#' @description Inputs data frame and two sets of outlier values, outputs list with three data frames
#' @param x Original vector of x values
#' @param y Original vector of y values
#' @param outsx Outlier values to be removed from x
#' @param outsy Outlier values to be removed from y


.extractoutlier <- function(x, y, outsx, outsy){
  
  # Remove x outliers from x and y
  if(length(outsx) > 0){
    x1 <- x[!x %in% outsx]
    y1 <- y[!x %in% outsx]
  } else {
    x1 <- x
    y1 <- y
  }
  xy1 <- unname(cbind(x1, y1))
  
  # Remove y outliers from x and y
  if(length(outsy) > 0){
    x2 <- x[!y %in% outsy]
    y2 <- y[!y %in% outsy]
  } else {
    x2 <- x
    y2 <- y
  }
  xy2 <- unname(cbind(x2, y2))
  
  # Remove x and y outliers from x and y
  if(length(outsx) > 0 && length(outsy) > 0){
    x3 <- x[!x %in% outsx & !y %in% outsy]
    y3 <- y[!x %in% outsx & !y %in% outsy]
  } else {
    x3 <- x
    y3 <- y
  }
  xy3 <- unname(cbind(x3, y3))
  
  # Combine results
  res <- unname(list(xy1, xy2, xy3))
  res <- unique(res)
  
  return(res)
  
}

#' Select a p-value from a vector of p-hacked p-values
#' @description Takes a vector of p-values and selects the smallest, first significant, or smallest significant p-value.
#' @param ps Vector of p values
#' @param strategy String value: One out of "firstsig", "smallest", "smallest.sig"
#' @param alpha Significance level (default: 0.05)

.selectpvalue <- function(ps, strategy, alpha){
  
  p.final <- NA
  p.orig <- ps[1]
  
  # Select smallest significant p-value
  if(strategy == "smallest.sig"){
    
    if(min(ps) < alpha){
      p.final <- min(ps)
    } else {
      p.final <- p.orig
    }
    
    # Select first significant p-value
  } else if (strategy == "firstsig") {
    
    if(min(ps) < alpha){
      p.final <- ps[which(ps < alpha)[1]]
    } else {
      p.final <- p.orig
    }
    
    # Select smallest p-value
  } else if (strategy == "smallest") {
    p.final <- min(ps)
  }
  
  return(p.final)
  
}

#' Compute R squared for the t-test
#' @param x values of group 1
#' @param y values of group 2

.compR2t <- function(x, y){
  grandmean <- mean(c(x, y))
  sst <- sum((c(x,y)-grandmean)^2)
  sse <- sum((x-mean(x))^2)+sum((y-mean(y))^2)
  return(1-(sse/sst))
}

#' Compute Cohen's d
#' @description Compute Cohen's d from t-value with equal sized groups of size n
#' @param t t-value
#' @param n sample size per group

.compCohensD <- function(t, n){
  t*sqrt(2/n)
}
