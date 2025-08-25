# Contains functions for simulations

# Libraries ====
library(simstudy)

# Functions ====

# generate data function
# multinomial distribution
snsp_data_gen = function(n, def) {
  data = genData(n, def)
  data$disease = data$test = 0
  data[data$td == 1 | data$td == 2, "disease"] = 1
  data[data$td == 1 | data$td == 3, "test"] = 1
  return(data)
}
# can also done in stages
# i.e. gen f(D), then f(T|D)
# but slower based on test, but clearer
# see tests2.R in simulated data/ folder

# create pvb problem, biased verification probability
# MAR
pvb_gen = function(complete_data, p1, p0, gen_na = TRUE) {  # p1: P(V=1|T=1), p0: P(V=1|T=0)
  data = complete_data
  # verification
  data$verify = NA
  t1 = which(data$test == 1); t0 = which(data$test == 0)
  data$verify[t1] = rbinom(nrow(data[t1]), 1, p1)
  data$verify[t0] = rbinom(nrow(data[t0]), 1, p0)
  # missing disease status
  if (gen_na == TRUE) {  # default
    data$disease[data$verify == 0] = NA
    return(data)
  } else if (gen_na == FALSE) {
    return(data)
  }  # use FALSE just to gen verification status, without specifying NA to disease for unverified
     # use this to check logit(V|T) significance
}  # note two calls to RNG here, make sure to set.seed()
# For MAR, may also extend to for P(V|T,X) --, this was done in deGroot 2011b, albeit less clearly in the study
# For the moment stick to P(V=T), while the inclusion of Xs is meant to provide better info about D
# for imputation & getting predicted Pr
pvb_mar_gen = pvb_gen

