
# Goal: Confirm understanding of why some small-N, binary-Y scenarios have bias even in individual studies.

# PRELIMINARIES  -------------------------------------------------


# helper fns
code.dir = here()
setwd(code.dir)
source("helper_JTE.R")


setwd(data.dir)
agg = fread("agg.csv")

# bad scens only: 
b = fread("agg_excluded_scens_biased_yi.csv")
nuni(b$scen.name)

# first row per iterate only
bf = b[ !duplicated(b$scen.name), ]
nrow(bf)

# initialize global variables that describe estimate and outcome names, etc.
# this must be after calling wrangle_agg_local
init_var_names(.agg = agg)


# EXPLORE BIAS -------------------------------------------------

t = bf %>% select( scen.name,
                   param.vars.manip2, 
                   all_of(starts_with("sancheck")) )

# scens with large k are more interesting since biased yi wouldn't just be chance
View(t %>% filter(k.pub == 100))


# ~ Isolate an interesting scen  -------------------------------------------------

### Scen with p0 = 0.05
# here, the sanchecks on nY0, nY1, pY0 are fine, but yi itself is biased
scen = 2109
View(t %>% filter(scen.name == scen))

# look at this scen's sim parameters
( p = as.data.frame( bf %>% filter( scen.name == scen ) %>%
  select( all_of(param.vars.manip2) ) ) )

### Scen with p0 > 0.05
scen = 2169
View(t %>% filter(scen.name == scen))

# look at this scen's sim parameters
( p = as.data.frame( bf %>% filter( scen.name == scen ) %>%
                       select( all_of(param.vars.manip2) ) ) )


# BINARY Y  -------------------------------------------------

# ~ Langan's N distributions  -------------------------------------------------

# this is reasonable, though less variation than seen in SAPB-E
d = sim_meta(  k.pub = p$k.pub,
               t2a = p$t2a,
               Mu = p$Mu,
               true.dist = p$true.dist,
               p0 = p$p0,
               Ytype = p$Ytype,
               N.expr = p$N.expr )

mean(d$yi) # biased upward, as in agg data (2.44-2.47 vs. truth 2.3) 


hist(d$yi)









