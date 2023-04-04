## Negative Binomial distribution parameterized by mean (mu) and overdispersion parameter (theta).
## This parameterization is referred to as NEGBIN type I (Cameron and Trivedi, 1998) as cited by
## https://doi.org/10.1080/03610926.2018.1563164
##
## x ~ nbinom_type1(mu, theta), where E(x) = mu, Var(x) = (theta + 1) mu
## This should not be confused with the mu and shape parameterization of nbinom in R or the 'alternative' NB (neg_binomial_2_...) in stan
## Note using disp instead of theta because using theta gives the error
## > Error: Currently 'dirichlet' is the only valid prior for simplex parameters. See help(set_prior) for more details.
## when trying to fit the model.
library(brms) 

# helper functions for post-processing of the family
# cribbed from lognormal_natural.R
log_lik_nbinom_type1 <- function(i, prep) {
  mu <- prep$dpars$mu[, i]
  if(NCOL(prep$dpars$disp)==1){disp <- prep$dpars$disp}else
  {disp <- prep$dpars$disp[, i]}   ## [, i] if disp is modelled, without otherwise
  y <- prep$data$Y[i]
  size = mu/disp
  Vectorize(dnbinom)(x = y, mu = mu, size = size, log = TRUE)
}


posterior_predict_nbinom_type1 <- function(i, prep, ...) {
  mu <- prep$dpars$mu[, i]
  if(NCOL(prep$dpars$disp)==1){disp <- prep$dpars$disp}else
  {disp <- prep$dpars$disp[, i]}   ## [, i] if disp is modelled, without otherwise
  size = mu/disp
  rnbinom(n = 1, mu = mu, size = size)
}

posterior_epred_nbinom_type1 <- function(prep) {
  mu <- prep$dpars$mu
  return(mu)
}

# define custom family
nbinom_type1 <- function(link_mu = "identity", link_disp = "identity")
    custom_family(
 name = "nbinom_type1", 
 dpars = c("mu", "disp"), #expectd value and variance scale (>1)
 links = c(link_mu, link_disp), # only used if parameters are predicted
 lb = c(0, 0),
 ub = c(NA, NA),
 type = "int" #category of response variable
## log_lik = log_lik_nbionm_type1
## posterior_predict = posterior_predict_nbionm_type1,
## posterior_epred = posterior_epred_nbionm_type1
)

# additionally required Stan code
# Not clear if the arguments should be ints and reals (as in the stan documentation) instead of int and real
stan_nbinom_type1 <- "
  real nbinom_type1_lpmf(int y, real mu, real disp) {
    return neg_binomial_2_lpmf(y | mu, mu/(disp));
  }
  int nbinom_type1_rng(real mu, real disp) {
    return neg_binomial_2_rng(mu, mu/(disp));
  }
"

# example

if(FALSE) {
    
library(dplyr)
library(tidyr)

## create data

### define variables
group <- c(1:2)
x <- 0:50
theta <- c(2,5)

## Create data
data <- crossing(x = x, group) %>%
  mutate(mu = 10*(1 + x), theta = theta[group]) %>%
  rowwise() %>%
  mutate(y = rnbinom(n = 1, mu = mu, size = mu/theta), group = factor(group)) %>%
  select(-mu)





## fit models
nbinom_type1_vars <- stanvar(scode = stan_nbinom_type1, block = "functions")

## Std regression model
get_prior(y ~ 0 + Intercept + x,
    family = nbinom_type1(),
    data = data,
    stanvar = nbinom_type1_vars
    )

bprior1 <- prior(gamma(0.01, 0.01), class = "b", lb = 0) +
    prior(gamma(0.01, 0.01), class = "disp", lb = 0)



brm(y ~ 0 + Intercept + x,
    family = nbinom_type1(),
    data = data,
    stanvars = stanvars,
    prior = bprior1
    ) -> model_test1


get_prior(bf(y ~ 0 + Intercept + x,
    disp ~ 0 + group),
    family = nbinom_type1(),
    data = data,
    stanvars = stanvars
    )


## Distributional regression model

bprior2 <- prior(gamma(0.01, 0.01), class = "b", lb = 0)

brm(bf(y ~ 0 + Intercept + x,
    disp ~ 0 + group),
    family = nbinom_type1(),
    data = data,
    stanvars = stanvars,
    prior = bprior2
    ) -> model_test2

summary(model_test1)
summary(model_test2)

model_test1 |> conditional_effects() |> plot(points = TRUE)
model_test2 |> conditional_effects() |> plot(points = TRUE)

model_test1 |> conditional_effects(method = "posterior_predict")
model_test2 |> conditional_effects(method = "posterior_predict")

loo(model_test1, model_test2)
}
