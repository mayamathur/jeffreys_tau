# FROM METASENS (FUNCTIONS_V2.R):

########################### FN: SIMULATE 1 STUDY, SCHLESSELMAN-STYLE ###########################
# p1 = P(X=1)
# p.int = P(Y=1 | X=0, U=0), i.e., intercept probability for logistic model

sim_one_study = function(.muB, .sigB, .Mt, .V, .p1, .p.int, .mu.N) {
  
  # simulate total N for each study
  min.N = 150
  N = round( runif( n = 1, min = min.N, max = min.N + 2*( .mu.N - min.N ) ) ) # draw from uniform centered on .mu.N
  
  # draw log-bias
  B = rnorm( n = 1, mean = .muB, sd = .sigB )
  
  # draw population true effect
  Mi = rnorm( n=1, mean=.Mt, sd=sqrt(.V) )
  
  # calculate population confounded effect
  Mi.c = Mi + B
  
  ###### Simulate Data For Individual Subjects ######
  # simulate X
  X = rbinom(size=1, n=N, prob=.p1)

  # assume RR_XU = RR_UY and call it gamma
  # compute it from definition of Tyler's bias factor
  # ~~~~~~~~~~~~ NOTE: THIS COULD BE NA IF B < 0; NEED TO BE CAREFUL CHOOSING PARAMETERS
  gam = exp(B) + sqrt( exp(B)^2 - exp(B) )

  # P( U=1 | X=1 ) - treatment group confounder prevalence
  p.u1x1 = 1
  
  # P( U=1 | X=0 ) - control group confounder prevalence
  # calculate using Schlesselman assumption (my notes, pg 65)
  p.u1x0 = ( exp(Mi) * (1 + (gam-1)*p.u1x1 ) - exp(Mi.c) ) / ( (gam-1)*exp(Mi.c) )
  
  # sanity check: these should match by Schlesselman
  # exp(Mi) * ( 1 + (gam-1)*p.u1x1 ) / (1 + (gam-1)*p.u1x0 )
  # exp(Mi.c)
  
  # generate U
  probs = rep( p.u1x0, length(X) )
  probs[X==1] = p.u1x1  # vector of P( U=1 | X=x )
  table(probs, X)  # check it
  U = rbinom(size=1, n=N, prob=probs)
  table(U, X) # check it
  
  # generate D - outcome itself
  linpred = log(p.int) + log(gam)*U + Mi*X  # Mi is already on log scale
  
  # exp here because log-RR model
  D = rbinom( size=1, n=N, prob=exp(linpred) )  
  
  # sanity check: see whether unconfounded RR generated here is correct
#   mod = glm(D ~ X + U, family=binomial)
#   exp(coef(mod)); exp(Mi)
#   # and confounded one
#   mod = glm(D ~ X, family=binomial)
#   exp(coef(mod)); exp(Mi.c)
  
  # sanity check -- look at outcome probs for each stratum of X and U
  #temp = data.frame( U, X, prob=expit(linpred) )
  #require(ggplot2)
  #ggplot(aes(y = prob, x = as.factor(X), color = as.factor(U) ), data = temp) + geom_point(size=3) + ylab("P(D=1)") + theme_bw()
  
  # calculate deaths and sample sizes in each group
  n1 = sum(X)  # number with X=1
  n0 = length(X) - sum(X)
  d1 = sum(D[X==1])  # number of deaths in Tx group
  d0 = sum(D[X==0])  # number of deaths in control group
  ###### end of individual subject data simulation ######

  # calculate confounded ES for this study using metafor (see Viechtbauer "Conducting...", pg 10)
  require(metafor)
  temp = escalc( measure="RR", ai=d1, bi=n1-d1, ci=d0, di=n0-d0 )  # returns on log scale
  yi.c = temp$yi
  vyi = temp$vi
  
  ##### adjusted RR
  # also calculate its adjusted RR (including the unmeasured confounder)
  # make 3D contingency table stratified by U
  library(car)
  X2 = recode( X, "0 = 'b.No'; 1 = 'a.Yes'")
  D2 = recode( D, "0 = 'b.No'; 1 = 'a.Yes'")
  U2 = recode( U, "0 = 'b.No'; 1 = 'a.Yes'")
  tab = table(D2, X2, U2)  # set up as on Modern Epi, pg 247

  # fix zero cells by replacing with a count of 0.5 observation
  tab[ tab==0 ] = 0.5
  # note to self: if we don't do this, then MH risk ratio has 0 in denom
  # if BOTH strata have zero cells, in which case we get an error when 
  # trying to meta-analyze

# check if variables are ALWAYS same value
#   if( length(levels(as.factor(D))) == 1 | length(levels(as.factor(X))) == 1 | length(levels(as.factor(U))) == 1 ){
#     cat("\nOne of the variables is always same value")
#     browser()
#   }

# # test only - compute stratum RRs by hand
# 
#   RR.U0 = ( tab[1,1,1] / (tab[1,1,1] + tab[1,2,1]) ) / ( tab[2,1,1] / (tab[2,1,1] + tab[2,2,1]) )
#   RR.U1 = ( tab[1,1,2] / (tab[1,1,2] + tab[1,2,2]) ) / ( tab[2,1,2] / (tab[2,1,2] + tab[2,2,2]) )

#   # expand the table
#   df = as.data.frame(tab)
#   df2 = df[ rep(row.names(df), df$Freq) , 1:3]

#   m = glm( D ~ X + U, data = df2, family = binomial )
#   OR = coef(m)[["X"]]

  # compute Mantel-Haenzsel risk ratio
  MH = MH_risk_ratio(tab)
  yi.lr = log( MH$RR )
  vyi.lr = MH$var.log.RR

  if( is.infinite(yi.lr) | is.na(yi.lr) ) browser()

  return( data.frame(Mi, yi.c, vyi, yi.lr, vyi.lr) )
  #return( data.frame(Mi, yi.c, vyi) )
}


# FROM MRM (HELPER_MRM):

########################### FN: SIMULATE 1 STUDY ###########################

# mu = true effect size as raw mean difference
# V = true variance of true effects
# muN = mean sample size in each study
# minN = minimum sample size 
# sd.w = SD within each group (2-group experiment)


# potentially with clustering
sim_one_study2 = function(b0, # intercept
                          bc, # effect of continuous moderator
                          bb, # effect of binary moderator
                          V, 
                          Vzeta, # used to calcuate within-cluster variance
                          zeta1,  # scalar cluster random intercept for this study's cluster
                          muN,
                          minN,
                          sd.w,
                          true.effect.dist = "normal"
) {
  
  # # @test for m=1 case
  # # TEST ONLY
  # b0 = 0.5 # intercept
  # bc = 0.5 # effect of continuous moderator
  # bb = 1 # effect of binary moderator
  # V = .5
  # Vzeta = 0.25
  # zeta1 = -0.2
  # muN = 100
  # minN = 50
  # sd.w = 1
  # true.effect.dist = "normal"
  
  if( !true.effect.dist %in% c("normal", "expo") ) stop("True effect dist not recognized")
  
  ##### Simulate Sample Size and Fixed Design Matrix for This Study #####
  # simulate total N for this study
  N = round( runif( n = 1, min = minN, max = minN + 2*( muN - minN ) ) ) # draw from uniform centered on muN
  
  # simulate study-level moderators (each a scalar)
  Zc = rnorm( n = 1, mean = 0, sd = 1)
  Zb = rbinom( n = 1, size = 1, prob = 0.5)
  
  # mean (i.e., linear predictor) conditional on the moderators and cluster membership
  mu = b0 + zeta1 + bc*Zc + bb*Zb
  # all that follows is that same as in NPPhat, except incorporating clustering as in SAPB
  
  ##### Draw a Single Population True Effect for This Study #####
  if ( true.effect.dist == "normal" ) {
    Mi = rnorm( n=1, mean=mu, sd=sqrt(V - Vzeta) )
  }
  if ( true.effect.dist == "expo" ) {
    # within-cluster variance = total - between
    Vwithin = V - Vzeta
    # set the rate so the heterogeneity is correct
    Mi = rexp( n = 1, rate = sqrt(1/Vwithin) )
    # now the mean is sqrt(V) rather than mu
    # shift to have the correct mean (in expectation)
    Mi = Mi + (mu - sqrt(Vwithin))
  }
  
  ###### Simulate Data For Individual Subjects ######
  # group assignments
  X = c( rep( 0, N/2 ), rep( 1, N/2 ) )
  
  # simulate continuous outcomes
  # 2-group study of raw mean difference with means 0 and Mi in each group
  # and same SD
  Y = c( rnorm( n = N/2, mean = 0, sd = sd.w ),
         rnorm( n = N/2, mean = Mi, sd = sd.w ) )
  
  # calculate ES for this study using metafor (see Viechtbauer "Conducting...", pg 10)
  require(metafor)
  ES = escalc( measure="SMD",   
               n1i = N/2, 
               n2i = N/2,
               m1i = mean( Y[X==1] ),
               m2i = mean( Y[X==0] ),
               sd1i = sd( Y[X==1] ),
               sd2i = sd( Y[X==0] ) ) 
  yi = ES$yi
  vyi = ES$vi
  
  return( data.frame( Mi, # study's true effect size; if within-cluster heterogeneity is zero, will be equal to mu
                      mu, # study's linear predictor conditional on the moderators and cluster membership
                      zeta1,
                      Zc,
                      Zb,
                      yi,
                      vyi ) )
}