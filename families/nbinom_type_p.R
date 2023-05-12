## Negative Binomial distribution parameterized by mean (mu), overdispersion scaling parameter (kappa), and shape (exponent) parameter p.
## This parameterization is referred to as NEGBIN type p by Greene (2008) [https://doi.org/10.1016/j.econlet.2007.10.015]
##
## x ~ nbinom_typep(mu, kappa), where E(x) = mu, Var(x) = mu (1 + kappa mu^(p-1)) = mu + kappa mu^p
## This should not be confused with the mu and shape parameterization of nbinom in R or the 'alternative' NB (neg_binomial_2_...) in stan
## Note using kappa instead of disp because using disp is a reserved term in `brms`

library(brms) 

run_examples <- TRUE

# helper functions for post-processing of the family
# cribbed from lognormal_natural.R

log_lik_nbinom_typep <- function(i, prep) {
  mu <- prep$dpars$mu[, i]
  if(NCOL(prep$dpars$kappa)==1) {
      kappa <- prep$dpars$kappa
  } else {
    kappa <- prep$dpars$kappa[, i]
  }   ## [, i] if kappa is modelled, without otherwise

  if(NCOL(prep$dpars$p)==1) {
      p <- prep$dpars$p
  } else {
    p <- prep$dpars$p[, i]
  }   ## [, i] if p is modelled, without otherwise
  
  y <- prep$data$Y[i]
  llik <- function(y, mu, kappa, p) {
    mu2mp <- mu^(2-p)
    s <- mu/(mu + mu2mp/kappa)
    lpmf <- lgamma(mu2mp/kappa + y) + (mu2mp/kappa) * log(s) + y * log(1-s) - lgamma(y + 1) - lgamma(mu2mp/kappa) - 2 * log(kappa)
    return(lpmf)
  }
  
  llik(y, mu, kappa, p)
  Vectorize(llik(y, mu, kappa, p) )
  
}

log_lik_nbinom_type1 <- function(i, prep) {
  prep$dpars$p <- 1
  return(log_lik_nbinom_typep(i, prep))
}

log_lik_nbinom_type2 <- function(i, prep) {
  prep$dpars$p <- 2
  return(log_lik_nbinom_typep(i, prep))
}

posterior_predict_nbinom_typep <- function(i, prep, ...) {
  mu <- prep$dpars$mu[, i]
  if(NCOL(prep$dpars$kappa)==1) {
    kappa <- prep$dpars$kappa
  } else {
    kappa <- prep$dpars$kappa[, i]
  }
  if(NCOL(prep$dpars$p)==1) {
      p <- prep$dpars$p
  } else {
    p <- prep$dpars$p[, i]
  }
  size = 1/(kappa*mu^(p-2))
  rnbinom(n = prep$ndraws, mu = mu, size = size)
}

posterior_predict_nbinom_type1 <- function(i, prep, ...) {
  prep$dpars$p <- 1
  return(posterior_predict_nbinom_typep(i, prep))
}

posterior_predict_nbinom_type2 <- function(i, prep, ...) {
  prep$dpars$p <- 2
  return(posterior_predict_nbinom_typep(i, prep))
}


posterior_epred_nbinom_typep <- function(prep) {
  mu <- prep$dpars$mu
  return(mu)
}

posterior_epred_nbinom_type1 <- function(prep) {
  mu <- prep$dpars$mu
  return(mu)
}

posterior_epred_nbinom_type2 <- function(prep) {
  mu <- prep$dpars$mu
  return(mu)
}

# define custom families
nbinom_typep <- function(link_mu = "identity", link_kappa = "identity", link_p = "identity")
    custom_family(
 name = "nbinom_typep", 
 dpars = c("mu", "kappa", "p"), #expectd value and variance scale (> 0)
 links = c(link_mu, link_kappa), 
 lb = c(0, 0, 0),
 ub = c(NA, NA, NA),
 type = "int", #category of response variable
 log_lik = log_lik_nbinom_typep,
 posterior_predict = posterior_predict_nbinom_typep,
 posterior_epred = posterior_epred_nbinom_typep
 )

nbinom_type1 <- function(link_mu = "identity", link_kappa = "identity")
    custom_family(
 name = "nbinom_type1", 
 dpars = c("mu", "kappa"), #expectd value and variance scale (> 0)
 links = c(link_mu, link_kappa), 
 lb = c(0, 0),
 ub = c(NA, NA),
 type = "int", #category of response variable
 log_lik = log_lik_nbinom_type1,
 posterior_predict = posterior_predict_nbinom_type1,
 posterior_epred = posterior_epred_nbinom_type1
 )

nbinom_type2 <- function(link_mu = "identity", link_kappa = "identity")
    custom_family(
 name = "nbinom_type2", 
 dpars = c("mu", "kappa"), #expectd value and variance scale (> 0)
 links = c(link_mu, link_kappa), 
 lb = c(0, 0),
 ub = c(NA, NA),
 type = "int", #category of response variable
 log_lik = log_lik_nbinom_type2,
 posterior_predict = posterior_predict_nbinom_type2,
 posterior_epred = posterior_epred_nbinom_type2
)

# additionally required Stan code
# Not clear if the arguments should be ints and reals (as in the stan documentation) instead of int and real
stan_nbinom_typep <- "
 real nbinom_typep_lpmf(int y, real mu, real kappa, real p) {
    // Formula based on eq 2-14 in Greene (2008) which uses the theta formulation for the pdf
    // NEGBIN 1 is when the generalized NEGBIN P has p = 1
    // Note: mu = lambda
    //       kappa = 1/theta
    // and g^(-1)(kappa) = theta = 1/kappa
    //     d g^(-1)/d kappa = -1/kappa^2
    // Thus, 
    real mu2mp = mu^(2 - p);
    real s = mu/(mu + mu2mp/kappa);
    real lpmf = lgamma(mu2mp/kappa + y)
                 + (mu2mp/kappa) * log(s)
                 + y * log(1-s)
                 - lgamma(y + 1)
                 - lgamma(mu2mp/kappa)
                 - 2 * log(kappa);
     return lpmf;
  }

  real nbinom_type1_lpmf(int y, real mu, real kappa) {
       return nbinom_typep_lpmf(y | mu, kappa, 1);
  }

real nbinom_type2_lpmf(int y, real mu, real kappa) {
       return nbinom_typep_lpmf(y | mu, kappa, 2);
}

  int nbinom_typep_rng(real mu, real kappa, real p) {
    real phi = 1/(kappa * mu^(p-2));
    return neg_binomial_2_rng(mu, phi);
  }

int nbinom_type1_rng(real mu, real kappa) {
    return nbinom_typep_rng(mu, kappa, 1);
}

int nbinom_type2_rng(real mu, real kappa) {
    return nbinom_typep_rng(mu, kappa, 2);
}
"

# example

if(run_examples) {
  
  library(dplyr)
  library(tidyr)

  ## ensure ncores =1 if there's <=2 or 4 if there's > 6
  ncores = min(max(parallel::detectCores()-2, 1), 4)
  nchains = ncores
  options(mc.cores = ncores)
  ## create data

  ### define variables

  x <- (0:50)*4
  group <- c(1)
  kappa <- c(2)

  ## Create data
  data1 <- crossing(group, rep = 1:10, x = x) %>%
    mutate(mu = 10*(1 + x), kappa = kappa[group]) %>%
    rowwise() %>%
    mutate(y = rnbinom(n = 1, mu = mu, size = mu/kappa),
           group = factor(group)) %>%
    select(-c(mu, kappa, rep))





  ## fit models
  stanvar <- stanvar(scode = stan_nbinom_typep, block = "functions")

  centered_intercept = FALSE
  
  if(centered_intercept) {
    formula <- formula(y ~ 1 + x)
  } else {
    formula <- formula(y ~ 0 + Intercept + x)
    }
  ## Std regression model
  get_prior(formula = formula, 
            family = nbinom_type1(),
            data = data1,
            stanvar = stanvar
            )

  bprior1 <- prior(gamma(0.01, 0.01), class = "b", lb = 0) +
    prior(gamma(0.01, 0.01), class = "kappa", lb = 0) +
    ## class intercept only works if you use `1 + ...`
    ## You cannot use it with `0 + Intercept + ...`
    if(centered_intercept) {
      prior(gamma(0.01, 0.01), class = "Intercept", lb = 0)
    } else {
      NULL #prior(gamma(0.01, 0.01), class = "b", coef = "Intercept", lb = 0)
    }

  start_time <- Sys.time()
  
  brm(formula = formula,
      family = nbinom_type1(),
      data = data1,
      stanvar = stanvar,
      prior = bprior1,
      control = list(adapt_delta = 0.95,
                     max_treedepth = 12),
      chains = nchains,
      cores = ncores,
      iter = 5000
      ) -> model_test1


  end_time <- Sys.time()
  total_time <- end_time - start_time
print(paste("Elapsed wall time", round(total_time, 2), "min"))

  summary(model_test1)
  model_test1 |> conditional_effects() |> plot(points = TRUE)
  model_test1 |> conditional_effects(method = "posterior_predict")



## Two groups
  group <- c(1,2)
  kappa <- c(2, 5)
  data2 <- crossing(group, rep = 1:10, x = x) %>%
    mutate(mu = 10*(1 + x), kappa = kappa[group]) %>%
    rowwise() %>%
    mutate(y = rnbinom(n = 1, mu = mu, size = mu/kappa),
           group = factor(group)) %>%
    select(-c(mu, kappa, rep))
  
  get_prior(bf(y ~ 0 + Intercept + x,
               kappa ~ 0 + group),
            family = nbinom_type1(),
            data = data,
            stanvar = stanvar
            )

  bprior2 <- prior(gamma(0.01, 0.01), class = "b", lb = 0) +
    prior(gamma(0.01, 0.01), dpar = "kappa", lb = 0)

  start_time <- Sys.time()
    
  brm(bf(y ~ 0 + Intercept + x,
         kappa ~ 0 + group),
      family = nbinom_type1(),
      data = data2,
      stanvar = stanvar,
      prior = bprior2,
      control = list(adapt_delta = 0.95,
                     max_treedepth = 14),
      chains = nchains,
      cores = ncores,
      ) -> model_test2

  end_time <- Sys.time()
  total_time2 <- end_time - start_time
  print(paste("Elapsed wall time", round(total_time2, 2), "min"))
  
  summary(model_test2)


  model_test2 |> conditional_effects() |> plot(points = TRUE)
  model_test2 |> conditional_effects(method = "posterior_predict")

#  loo(model_test1, model_test2)
}
