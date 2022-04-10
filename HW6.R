## HW 6

# requires Rtools : https://cran.rstudio.com/bin/windows/Rtools/
if (!require("remotes")) {
  install.packages("remotes")
  remotes::install_github("avehtari/BDA_course_Aalto", subdir = "rpackage", upgrade="never")
}
if (!require("pracma"))install.packages("pracma")
  
library(pracma)
library(aaltobda)
data("algae")
# the data is now stored in the variable 'algae'

# test data to check correctness
algae_test <- c(0, 1, 1, 0, 0, 0)


# pi : prob 
# y  : observations in 'algea'
# Beta(2,10) prior
a <- 2
b <- 10 

data <- algae

## a) formulate followings, in the format Beta(,)

# (1) likelihood p(y|pi)
likelihood <- function(pi,data){
  n <- length(data)
  x <- sum(data)
  likelihood <- dbinom(x, size = n, prob = pi)
  return(likelihood)
}

# (2) prior p(pi)
prior <- function(a,b){
  seq <- linspace(0, 1, 1000)
  prior <- dbeta(seq,a,b)
  return(prior)
}



# (3) posterior p(pi|y)
posterior <-function(data){ 
  n <- length(data)
  y <- sum(data)
  seq <- linspace(0, 1, 1000)
  # Hw 3.2
  posterior = dBeta(seq,a + y, b + n - y)
  return(posterior)
}

## b)
# What can you say about the value of the unknown π according to the observations and your
# prior knowledge? Summarize your results with a point estimate (i.e. E(π|y)) and a 90%
# posterior interval. Note! Posterior intervals are also called credible intervals and are different
# from confidence intervals. Note! In your report, use the values from the data algae, not
# algae_test.
beta_point_est <- function(alpha,beta,data){
  y <- sum(data)
  n <- length(data)
  # HW3,3
  return ((alpha + y) / (alpha + beta + n))
}

beta_interval <- function(a,b,data,prob){
  y <- sum(data)
  n <- length(data)
  
  alpha <- a + y
  beta <- b + n - y
  
  seq <- seq(0, 1, by = 0.025)
  # beta quantile function and used t return quantile values of the function.
  seq <- qbeta(seq,alpha,beta)
  
  #return (seq)
  return (quantile(seq, c(0.05, 0.95)))
}

# check
beta_point_est(a,b,algae_test)
beta_interval(a,b,algae_test, 0.9)

# answer
beta_point_est(a,b,data)
beta_interval(a,b,data, 0.9)

 
## c)
# What is the probability that the proportion of monitoring sites with detectable algae levels π is
# smaller than π0 = 0.2 that is known from historical records?

beta_low <- function(alpha,beta,data,pi_0){
  y <- sum(data)
  n <- length(data)
  
  alpha <- a + y
  beta <- b + n - y
  
  # create cumulative distribution function of the beta distribution.
  return (pbeta(pi_0,alpha,beta))
}  
# check
beta_low(a,b,algae_test,0.2)

#answer
beta_low(a,b,data,0.2)


## d)
# What assumptions are required in order to use this kind of a model with this type of data? 
# (No  need to discuss exchangeability yet, as it is discussed in more detail in BDA Chapter 5 and
# Lecture 7)

# Assumption of binomial distribution. each event are indenpendent.

## e)
# Make prior sensitivity analysis by testing a couple of different reasonable priors and plot the
# different posteriors. Summarize the results by one or two sentences.

seq <- linspace(0, 1, 1000)

# 1) uniform distribution
post_uniform <- dbeta(seq,1+y,1+n-y)
plot(seq,post1, type="l",col="blue", xlab="pi", ylab="posterior")

# 2_) different beta
a2 = 10
b2 = 10
post_diff_beta <- dbeta(seq,a2+y,b2+n-y)

lines(seq,post_diff_beta,col="red")

post <- dbeta(seq,a+y,b+n-y)
lines(seq,post,col="green")

legend("topright",legend=c("B(2,10)","B(10,10)","uniform"),fill=c("green","red","blue"),border="white",box.lty=0,cex=1.5)

# With non-informative prior, posterior is costructed only by likelihood and for this data, posterior is similar to beta(2,10) prior
# But if beta parameter is set differently, center of posterior distribution is relatively more different that non-informative prior.
