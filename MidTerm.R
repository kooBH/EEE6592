## Midterm Replacement
## BDA assignment 4 

# requires Rtools : https://cran.rstudio.com/bin/windows/Rtools/
if (!require("remotes")) {
  install.packages("remotes")
  remotes::install_github("avehtari/BDA_course_Aalto", subdir = "rpackage", upgrade="never")
}
if (!require("pracma"))install.packages("pracma")

library(pracma)
library(aaltobda)

data("bioassay_posterior")

# b)

alpha = unlist(bioassay_posterior[1])
beta = unlist(bioassay_posterior[2])

mean_alpha = mean(alpha)
#mean_beta = sapply(beta, mean)
mean_beta = mean(beta)

var_alpha = var(alpha)
var_beta = var(beta)

MCSE_alpha = sqrt(var_alpha/4000)
MCSE_beta = sqrt(var_beta/4000)

quantile_alpha_05 = unlist(mcse_quantile(alpha, prob = 0.05 ))
quantile_alpha_95 = unlist(mcse_quantile(alpha, prob = 0.95 ))

quantile_beta_05 = unlist(mcse_quantile(beta, prob = 0.05 ))
quantile_beta_95 = unlist(mcse_quantile(beta, prob = 0.95 ))

print("+ alpha")
print(paste0("    mean     : ", mean_alpha ))
print(paste0("    MCSE     : ", MCSE_alpha ))
print(paste0(" 5% quantile : ", quantile_alpha_05 ))
print(paste0("95% quantile : ", quantile_alpha_95 ))

print("+ beta")
print(paste0("        mean : ", mean_beta))
print(paste0("        MCSE : ", MCSE_beta ))
print(paste0(" 5% quantile : ", quantile_beta_05 ))
print(paste0("95% quantile : ", quantile_beta_95 ))


# c) Implement a function for computing the log importance ratios
# when the importance sampling target distribution is the posterior dis -
# tribution, and the proposal distribution is the prior distribution from a)


log_importance_weights <- function(data, alpha , beta ) {
  x <-data$x
  y <-data$y
  n <-data$n
  
  return(bioassaylp(alpha = alpha, beta = beta, x = x, n = n, y = y))
}

# testing

alpha <- c(1.896, -3.6, 0.374, 0.964, -3.123, -1.581)
beta <- c(24.76, 20.04, 6.15, 18.65, 8.16, 17.4)
print(round(log_importance_weights(bioassay,alpha, beta), 2))

## [1] -8.95 -23.47 -6.02 -8.13 -16.61 -14.57

# d) normalized_importance_weights

# this function computes the normalized importance ratios
# from the unnormalized logratios in (c)

normalized_importance_weights <- function(data, alpha , beta){
  lw <- log_importance_weights(data,alpha,beta)
  
  # exponentiate the log ratio...
  elw <- exp(lw)
  
  # ... scale them such that they sum to one
  return(elw / sum(elw))
}

# testing
alpha <- c(1.896, -3.6, 0.374, 0.964, -3.123, -1.581)
beta <- c(24.76, 20.04, 6.15, 18.65, 8.16, 17.4)
print(round(normalized_importance_weights(bioassay, alpha = alpha, beta = beta), 2))

# e)  Sample 4000 draws of α and β from the prior distribution from a). 
# Compute and plot a histogram of the 4000 normalized importance ratios. 
# Use the functions you implemented in c) and d).

num_sample <- 4000
mean_a <- c(0,10)
cov_a  <- matrix(c(4,12,12,100),2,2)

data_e <- rmvnorm(num_sample,mean_a,cov_a)
alpha_e <- data_e[1:num_sample,1]
beta_e <- data_e[1:num_sample,2]

Normalized_Importance_Ratios <- normalized_importance_weights(bioassay,alpha = alpha_e, beta = beta_e)
hist(Normalized_Importance_Ratios)


# f)  Using the importance ratios, compute the importance sampling eective sample size Se and report it

S_eff <- function(data,alpha,beta) {
  # (10.4)
  return(1 / sum(normalized_importance_weights(data,alpha,beta)**2))
}

print(paste("S_eff 4000 : ",round(S_eff(bioassay, alpha = alpha_e, beta = beta_e),3)))

# test
alpha <- c(1.896, -3.6, 0.374, 0.964, -3.123, -1.581)
beta <- c(24.76, 20.04, 6.15, 18.65, 8.16, 17.4)
print(paste("S_eff test : ",round(S_eff(bioassay, alpha = alpha, beta = beta),3)))



# g) Explain in your own words what the importance sampling eective sample size represents.
# Also explain how the eective sample size is seen in the histogram of the weights that you plotted in e).



print(round(S_eff(bioassay, alpha = alpha_e, beta = beta_e),3))


# h) Implement a function for computing the posterior mean using importance sampling, and 
# compute the mean using your 4000 draws. Explain in your own words the computation for importance
# sampling. Below is an example how the function would work with the example values for alpha and beta above.
# Report the means for alpha and beta, and also the Monte Carlo standard errors (MCSEs) for the mean estimates. 
# Report the number of digits for the means based on the MCSEs


posterior_mean <- function(data, alpha, beta) {
  w = normalized_importance_weights(data, alpha, beta)
  sum_w <- sum(w)
  return(c(sum(alpha * w)/sum_w,sum(beta * w)/sum_w))
}




  # testing
  alpha <- c(1.896, -3.6, 0.374, 0.964, -3.123, -1.581)
  beta <- c(24.76, 20.04, 6.15, 18.65, 8.16, 17.4)
  posteriors_means <- posterior_mean(bioassay, alpha = alpha, beta = beta)
  print(round(posteriors_means,3))

  ## Report
  
posterior_squared_mean <- function(data, alpha, beta) {
  w <- normalized_importance_weights(data, alpha, beta)
  sum_w <- sum(w)
  return(c(sum(alpha**2 * w)/sum_w,sum(beta**2 * w)/sum_w))
  }
  
  num_sample <- 4000
  data_h <- rmvnorm(num_sample,mean_a,cov_a)
  alpha_h <- data_h[1:num_sample,1]
  beta_h <- data_h[1:num_sample,2]
  
  
  E <- posterior_mean(bioassay, alpha = alpha_h, beta = beta_h)
  squared_E <- posterior_squared_mean(bioassay, alpha = alpha_h, beta = beta_h)
  s_eff = S_eff(bioassay, alpha = alpha_h, beta = beta_h)
  
  var_alpha_h <- squared_E[1] - E[1]**2
  var_beta_h <- squared_E[2] - E[2]**2

  MCSE_alpha_h <- sqrt(var_alpha_h) / sqrt(s_eff)
  MCSE_beta_h <- sqrt(var_beta_h) / sqrt(s_eff)
  
  print(paste('alpha_mean: ', round(E[1],3), ' alpha mcse: ', MCSE_alpha_h))
  print(paste('beta_mean: ', round(E[2],3), ' beta mcse: ', MCSE_beta_h))

