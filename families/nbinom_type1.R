## Negative Binomial distribution parameterized by mean (mu) and overdispersion parameter (theta).
## This parameterization is referred to as NEGBIN type I (Cameron and Trivedi, 1998) as cited by
## https://doi.org/10.1080/03610926.2018.1563164
##
## x ~ nbinom_type1(mu, theta), where E(x) = mu, Var(x) = (theta + 1) mu
## This should not be confused with the mu and shape parameterization of nbinom in R or the 'alternative' NB (neg_binomial_2_...) in stan
## Note using disp instead of theta because using theta is a reserved term in `brms`

library(brms) 
library(shinybrms)

run_examples <- TRUE

# helper functions for post-processing of the family
# cribbed from lognormal_natural.R
log_lik_nbinom_type1 <- function(i, prep) {
  mu <- prep$dpars$mu[, i]
  if(NCOL(prep$dpars$disp)==1) {
      disp <- prep$dpars$disp
  } else { disp <- prep$dpars$disp[, i] }   ## [, i] if disp is modelled, without otherwise
  y <- prep$data$Y[i]
  Vectorize(dnbinom)(x = y, mu = mu, size = mu/disp, log = TRUE)
}


posterior_predict_nbinom_type1 <- function(i, prep, ...) {
  mu <- prep$dpars$mu[, i]
  if(NCOL(prep$dpars$disp)==1){ disp <- prep$dpars$disp
  } else { disp <- prep$dpars$disp[, i] }   ## [, i] if disp is modelled, without otherwise
  size = mu/disp
  rnbinom(n = prep$ndraws, mu = mu, size = size)
}

posterior_epred_nbinom_type1 <- function(prep) {
  mu <- prep$dpars$mu
  return(mu)
}

# define custom family
nbinom_type1 <- function(link_mu = "identity", link_disp = "identity")
    custom_family(
 name = "nbinom_type1", 
 dpars = c("mu", "disp"), #expectd value and variance scale (> 0)
 links = c(link_mu, link_disp), 
 lb = c(0, 0),
 ub = c(NA, NA),
 type = "int", #category of response variable
 log_lik = log_lik_nbinom_type1,
 posterior_predict = posterior_predict_nbinom_type1,
 posterior_epred = posterior_epred_nbinom_type1
)

# additionally required Stan code
# Not clear if the arguments should be ints and reals (as in the stan documentation) instead of int and real

## Below doesn't seemto work well.  I believe there are some numerical issues
## Possibly due to using a real instead of an int for disp 
## return neg_binomial_2_lpmf(y | mu, mu/(disp)) + log(mu) - 2 * log(disp);
##
## Noting that Gamma(n + 1) = n!, it follows that
## (n + disp - 1)_choose_n = (n + disp - 1)!/(n! * (disp - 1)!)
## Gamma(n + disp)/(Gamma(n+1)*Gamma(disp)

stan_nbinom_type1 <- "
  real nbinom_type1_lpmf(int y, real mu, real disp) {
    return (lgamma(y + mu/disp) - (lgamma(y + 1) + lgamma(mu/disp)) + y * log(disp) - (y + mu/disp) * log(1 + disp));
  }

  int nbinom_type1_rng(real mu, real disp) {
    return neg_binomial_2_rng(mu, mu/(disp));
  }
"

# example

if(run_examples) {
  
  library(dplyr)
  library(tidyr)
  library(modelr)
  library(tidybayes)
  
  ## ensure ncores =1 if there's <=2 or 4 if there's > 6
  ncores = min(max(parallel::detectCores()-2, 1), 4)
  nchains = ncores
  options(mc.cores = ncores)

  centered_intercept = FALSE
  

  stanvar <- stanvar(scode = stan_nbinom_type1, block = "functions")

  ## create data
  ## Single population data
  ### define variables
  group <- c(1)
  x <- 0:50
  theta <- c(5)

  ## Create data
  data1 <- crossing(group, rep = 1:100, x = x) %>%
    mutate(mu = 10*(1 + x), theta = theta[group]) %>%
    rowwise() %>%
    mutate(y = rnbinom(n = 1, mu = mu, size = mu/theta),
           group = factor(group))
  #%>%
   # select(-c(mu, theta))

  data_tmp <- data1 %>% group_by(x) %>% summarize(var = var(y), mean = mean(y))
  ## fit models


  if(centered_intercept) {
    formula <- formula(y ~ 1 + x)
  } else {
    formula <- formula(y ~ 0 + Intercept + x)
  }
  
  bprior1 <- prior(uniform(0, 200), class = "b", lb = 0, ub=200) +
    ## prior(gamma(500, 100), class = "disp", lb = 3)
    prior(uniform(0.1, 200), class = "disp", lb = 0, ub=200)
  ## Setting lb = 0 is non-sensical given data is overdispersed
  ## class intercept only works if you use `1 + ...`
  ## You cannot use it with `0 + Intercept + ...`

  
  if(centered_intercept) {
      bprior1 <- bprior1 + prior(uniform(1, 200), class = "Intercept", lb = 1)
  }

  ## init: Initial values for the sampler. If ‘NULL’ (the default) or
  ## ‘"random"’, Stan will randomly generate initial values for
  ## parameters in a reasonable range. If ‘0’, all parameters are
  ## initialized to zero on the unconstrained space. 
  ##
  ## You can get the initial values for a fit using
  ## get_inits(fit)
  ##
  ##
  init <- 0
  #init <- rep(list(list(b = c(10, 10), disp = 5)), nchains)
  
  start_time <- Sys.time()
  
  brm(formula = formula,
      family = nbinom_type1(),
      data = data1,
      #data2 = list(alpha = 500, beta = 100), 
      stanvar = stanvar,
      prior = bprior1,
      control = list(adapt_delta = 0.8,
                     max_treedepth = 12),
      chains = nchains,
      cores = ncores,
      init = init, 
      iter = 1000
      ) -> model_test1


  
  end_time <- Sys.time()
  total_time <- end_time - start_time
  print(paste("Wall time", total_time))
  summary(model_test1)


  model <- model_test1
  #model <- update(model_test1, ...)

  plot(model)
  ## Generate epred and plot
  data_grid <- data1 %>%
    data_grid(x = seq_range(c(1, 51), n = 200))
  
  data_epred <- add_epred_draws(newdata = data_grid, object = model)
  plot_epred <- ggplot(data = data_epred, aes(x = x, y = .epred)) + 
    ## Combine Scatter Plots and Model vs Data Plots
    geom_point(data = data1,
               aes(x = x, y = y), color = "red") +
    stat_lineribbon(aes(y = .epred), .width = c(.95, .8, 0.5, 0.25), color = "#08519C") +
#    scale_fill_brewer(palette = "Greys") +
    labs(title = "Expected Values") + 
    ylim(0, max(data1$y)*1.05)

  plot_epred
  ## Generate predictions and plot
  data_pred <- add_predicted_draws(newdata = data_grid, object = model)
    
  
  plot_pred <- ggplot(data = data_pred, aes(x = x, y = .prediction)) + #, color = male)) +
    ## Combine Scatter Plots and Model vs Data Plots
    geom_jitter(data = data1, aes(x = x, y = y), color = "red") +
    stat_lineribbon(aes(y = .prediction), .width = c(.95, .8, .5, 0.25), color = "#08519C", alpha = 0.5) +
    scale_fill_brewer(palette = "Greys", direction =  -1) +
    labs(title = "Predicted Values & Data") + 
    ylim(0, max(data1$y)*1.1)

  plot_pred
  
  # End simple fit


  ## Fit Mixture data
  ### define variables
  group <- c(1:2)
  x <- 0:50
  theta <- c(2,5)

  ## Create data
  data2 <- crossing(group, rep = 1:100, x = x) %>%
    mutate(mu = 10*(1 + x), theta = theta[group]) %>%
    rowwise() %>%
    mutate(y = rnbinom(n = 1, mu = mu, size = mu/theta),
           group = factor(group)) %>%
    select(-c(mu, theta))

  start_time <- Sys.time()
  
  brm(formula = formula,
      family = nbinom_type1(),
      data = data2,
      stanvar = stanvar,
      prior = bprior1,
      control = list(adapt_delta = 0.8),
      chains = nchains,
      cores = ncores,
      iter = 1000
      ) -> model_test2


  end_time <- Sys.time()
  total_time <- end_time - start_time
  print(paste("Wall time", total_time))

  


  ## fit disp ~ group models
  start_time <- Sys.time()
  
  brm(formula = formula,
      family = nbinom_type1(),
      data = data2,
      stanvar = stanvar,
      prior = bprior1,
      control = list(adapt_delta = 0.8,
                     max_treedepth = 12),
      chains = nchains,
      cores = ncores,
      iter = 1000
      ) -> model_test2


  end_time <- Sys.time()
  total_time <- end_time - start_time
  print(paste("Wall time", total_time))
  
  ## Distributional regression model

  get_prior(bf(y ~ 0 + Intercept + x,
               disp ~ 0 + group),
            family = nbinom_type1(),
            data = data2,
            stanvar = stanvar
            )

  bprior3 <- prior(uniform(0, 200), class = "b", lb = 0, ub = 200) +
    prior(uniform(0, 200), dpar = "disp", lb = 0, ub = 200)

  brm(bf(y ~ 0 + Intercept + x,
         disp ~ 0 + group),
      family = nbinom_type1(),
      data = data2,
      stanvar = stanvar,
      prior = bprior3,
      control = list(adapt_delta = 0.8),
      chains = nchains,
      cores = ncores,
      ) -> model_test3

  summary(model_test2)
  summary(model_test3)

  model_test2 |> conditional_effects() |> plot(points = TRUE)
  model_test3 |> conditional_effects() |> plot(points = TRUE)

  model_test2 |> conditional_effects(method = "posterior_predict")
  model_test3 |> conditional_effects(method = "posterior_predict")


  ## Get 10,000s of warnings
  ## > Warning in formals(fun) : argument is not a function
  loo(model_test2, model_test3)
}
