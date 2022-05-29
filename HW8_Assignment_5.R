

library(aaltobda)
data("bioassay")

### 1.

## a)

# Hint! Compute with log-densities. Reasons are explained on page 261 of BDA3
# and lecture video 4.1. Remember that p1/p0 = exp(log(p1) − log(p0)). For your
# convenience we have provided functions that will evaluate the log-likelihood for given
# α and β (see bioassaylp() in the aaltobda package). Notice that you still need
# to add the prior yourself and remember the unnormalized log posterior is simply the
# sum of log-likelihood and log-prior. For evaluating the log of the Gaussian prior you
# can use the function dmvnorm from package aaltobda



density_ratio <- function(
    alpha_prop, 
    alpha_prev,
    beta_prop, 
    beta_prev, 
    x,
    y,
    n){
  
  unnormalized_log_posterior <- function(
    alpha,
    beta,
    mean,
    sigma,
    x,
    y,
    n){
    log_likelihood = bioassaylp(alpha = alpha, beta = beta, x = x, y = y, n = n)
    log_prior = dmvnorm(x = c(alpha, beta ), mean = mean, sigma = sigma, log = TRUE)
    return(log_likelihood + log_prior)
  }
  
  sigma_0 = matrix(data = c(4, 10, 10, 100), nrow = 2, ncol = 2)
  mu_0 = matrix(data = c(0 ,10), nrow = 2, ncol = 1)
  
  p1 = unnormalized_log_posterior(
    alpha = alpha_prop, 
    beta = beta_prop, 
    mean = mu_0, 
    sigma = sigma_0, 
    x = x, 
    y = y, 
    n = n)
  p0 = unnormalized_log_posterior(
    alpha = alpha_prev, 
    beta = beta_prev, 
    mean = mu_0, 
    sigma = sigma_0, 
    x = x,
    y = y, 
    n = n)
  ratio = exp(p1 - p0)
  return(ratio)
}

# test

density_ratio(alpha_prop = 1.89, alpha_prev = 0.374,
              beta_prop = 24.76, beta_prev = 20.04,
              x = bioassay$x, y = bioassay$y, n = bioassay$n)

density_ratio(alpha_prop = 0.374, alpha_prev = 1.89,
              beta_prop = 20.04, beta_prev = 24.76,
              x = bioassay$x, y = bioassay$y, n = bioassay$n)


## b) 

# Hint! Use a simple (normal) proposal distribution. Example proposals are α
# ∗ ∼
# N(αt−1, σ = 1) and β
# ∗ ∼ N(βt−1, σ = 5). There is no need to try to find optimal
# proposal but test some different values for the jump scale (σ). Remember to report
# the one you used. Efficient proposals are dicussed in BDA3 p. 295–297 (not part of
# the course). In real-life a pre-run could be made with an automatic adaptive control to
# adapt the proposal distribution.

# slide ch 11, p12
metropolis_bioassay <- function(
    alpha_0, 
    beta_0, 
    n_step, 
    n_warmup,
    bioassay){
  # data
  steps = matrix(rep(0, 2*(n_step+n_warmup)), ncol = 2)
  # c() : vector, column
  steps[1, ] <- c(alpha_0, beta_0)
  
  # Starting proposal 
  alpha_prev = alpha_0
  beta_prev = beta_0
  
  # t = 1, 2, ... n_step
  for (i in 1: (n_warmup+n_step) ){
    # (a) pick a proposal from the proposal distribution J 
    alpha_prop = rnorm(n = 1, mean = alpha_prev, sd = 1)
    beta_prop = rnorm(n = 1, mean = beta_prev , sd = 5)
    
    
    # (b) Calculate acceptance ratio
    acceptance_ratio = density_ratio(
      alpha_prop, 
      alpha_prev,
      beta_prop, 
      beta_prev, 
      bioassay$x, 
      bioassay$y,
      bioassay$n
    )
    # (c) set
    # rejection of a proposal step c is executed by generating a random number from U(0,1)
    # runif: The uniform distribution on the simplex
    if (runif(1) < acceptance_ratio){
      steps[i,] = c(alpha_prop, beta_prop)
      alpha_prev = alpha_prop
      beta_prev = beta_prop
    # rejection of a proposal
    } else {
      steps[i,] = c(alpha_prev, beta_prev)
    }
  }
  return(steps[n_warmup:(n_warmup+n_step),])
  #return(steps)
}
# Testing
n_step<-1000
n_warmup<-500
n_exp<-3
alpha_mat = matrix(rep(0, n_exp*(n_step)), ncol = n_exp)
beta_mat = matrix(rep(0, n_exp*(n_step)), ncol = n_exp)

steps = metropolis_bioassay(0.5, 5, n_step, n_warmup, bioassay)
alpha_mat[1:n_step,1]<-steps[1:n_step,1]
beta_mat[1:n_step,1]<-steps[1:n_step,2]

steps = metropolis_bioassay(1.5, 10, n_step, n_warmup, bioassay)
alpha_mat[1:n_step,2]<-steps[1:n_step,1]
beta_mat[1:n_step,2]<-steps[1:n_step,2]

steps = metropolis_bioassay(3, 20, n_step, n_warmup, bioassay)
alpha_mat[1:n_step,3]<-steps[1:n_step,1]
beta_mat[1:n_step,3]<-steps[1:n_step,2]

# Plot

matplot(alpha_mat, type = c("b"),pch=1, col = 1:n_exp)# plot
legend("topleft", legend = c(0.5,1.5,3), col=1:n_exp, pch=1) # optional legend

matplot(beta_mat, type = c("b"),pch=1, col = 1:n_exp)# plot
legend("topleft", legend = c(5,10,20), col=1:n_exp, pch=1) # optional legend
