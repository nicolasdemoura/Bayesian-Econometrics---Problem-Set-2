#####################################################################
# Topic: Bayesian Econometrics - Problem Set 2
# Instructor: Eduardo Mendes
# Course: Bayesian Econometrics
# Author: Nícolas de Moura
# Goal: Implement a Bayesian model to estimate the parameters of a non-linear model with Time-Varying Coefficients
#####################################################################
# Organize the working environment
#####################################################################
# Clear the workspace
rm(list=ls())

library(dplyr)
library(tidyr)
library(KFAS)
library(MASS)
library(mvtnorm)
library(coda)
library(ggplot2)
library(grid)       
library(progress)
library(stargazer)   

set.seed(20250420)

#####################################################################
# Simulate the data
#####################################################################

# Function to compute the Lambda matrix
get_Lambda <- function(taus, lambda) {
    Lambda <- matrix(0, nrow = length(taus), ncol = 3)
    for (i in 1:length(taus)) {
        Lambda[i, ] <- c(
            1,
            (1 - exp(-lambda * taus[i])) / (lambda * taus[i]),
            ((1 - exp(-lambda * taus[i])) / (lambda * taus[i])) - exp(-lambda * taus[i])
        )
    }
    return(Lambda)
}

simulate_data <- function(taus, n_obs, lambda, mu, A, H, Q) {
    Lambda <- get_Lambda(taus, lambda)

    beta <- matrix(0, nrow = n_obs, ncol = 3)
    y <- matrix(0, nrow = n_obs, ncol = length(taus))

    for (i in 2:n_obs) {
        eps <- mvrnorm(n = 1, mu = rep(0, length(taus)), Sigma = H)
        eta <- mvrnorm(n = 1, mu = rep(0, 3), Sigma = Q)
        beta[i, ] <- mu + A %*% (beta[i - 1, ] - mu) + eta
        y[i, ] <- Lambda %*% beta[i, ] + eps
    }

    y_df <- as.data.frame(y)
    colnames(y_df) <- paste0("y_", taus)
    beta_df <- as.data.frame(beta)

    return(list(y = y_df, beta = beta_df))
}

# Parameters
taus <- c(3, 6, 9, 12, 15, 18, 21, 24, 30, 36, 48, 60, 72, 84, 96, 108, 120, 180, 240, 360)
lambda_true <- 0.077
mu_true <- c(0.08, -1.5, 0)
A_true <- matrix(c(0.9, -0.1, -0.1, 
                   -0.1, 0.8, 0.05, 
                   -0.1, -0.1, 0.7), nrow = 3)
H_true <- diag(rep(0.05, length(taus)))
Q_true <- matrix(c(0.09, -0.01, 0.04, 
                   -0.01, 0.38, 0.01, 
                   0.04,  0.01, 0.80), nrow = 3)

# Simulate
ds <- simulate_data(taus, n_obs = 100, lambda_true, mu_true, A_true, H_true, Q_true)
y <- as.matrix(ds$y)
n_obs <- nrow(y)

taus        <- c(3,6,9,12,15,18,21,24,30,36,48,60,72,84,96,108,120,180,240,360)
n_obs       <- 100
sim         <- simulate_data(taus, n_obs, lambda_true, mu_true, A_true, H_true, Q_true)
y           <- as.matrix(sim$y)
beta_true   <- as.matrix(sim$beta)
T           <- n_obs; N <- ncol(y)

#################################################################
# Functions for sampling
#################################################################

# MH for A_unconstrained
# Adaptive proposal distribution for A
sample_A_unconstr <- function(A_curr, beta, mu, Q, sigma_prop = diag(0.01,9), accept_rate_target = 0.25) {
  # Sample from the proposal
  vecA <- as.numeric(A_curr)
  proposal <- as.numeric(rmvnorm(1, vecA, sigma_prop))
  A_prop <- matrix(proposal, 3, 3)

  # Rescale the eigenvalues to enforce stationarity
  e <- eigen(A_prop)
  vals <- e$values
  vals <- vals * min(0.999 / max(Mod(vals)), 1)
  A_prop <- Re(e$vectors %*% diag(vals) %*% solve(e$vectors))
  A_prop <- 0.8 * diag(3) + 0.2 * A_prop  # Ensure A is close to the identity matrix

  # Log posterior calculation
  logpost <- function(A) {
    res <- beta[-1,] - (mu + (beta[-nrow(beta),] - mu) %*% t(A))
    ll   <- sum(dmvnorm(res, sigma = Q, log = TRUE))
    lp   <- -sum(A^2) / (2 * sigma_A2)
    ll + lp
  }

  # Log acceptance ratio
  lr <- logpost(A_prop) - logpost(A_curr)
  
  # Update the proposal covariance adaptively
  acceptance_rate <- mean(lr > 0)  # Example calculation (you could track this over iterations)
  
  if (acceptance_rate < accept_rate_target) {
    sigma_prop <- 0.9 * sigma_prop  # Shrink the proposal covariance
  } else {
    sigma_prop <- 1.1 * sigma_prop  # Expand the proposal covariance
  }

  if (log(runif(1)) < lr) {
    return(list(A = A_prop, accept = TRUE, sigma_prop = sigma_prop))
  } else {
    return(list(A = A_curr, accept = FALSE, sigma_prop = sigma_prop))
  }
}

# MH for lambda
sample_lambda <- function(lambda_curr, y, beta, H, taus, sd_prop = 0.01, accept_rate_target = 0.30) {
  log_lambda_curr <- log(lambda_curr)
  log_lambda_prop <- rnorm(1, log_lambda_curr, sd_prop)
  lambda_prop <- exp(log_lambda_prop)

  loglike <- function(lam) {
    Z <- get_Lambda(taus, lam)
    sum(sapply(1:nrow(beta), function(t) {
      dmvnorm(y[t, ], mean = Z %*% beta[t, ], sigma = diag(H), log = TRUE)
    }))
  }

  log_accept_ratio <- (
    loglike(lambda_prop) - loglike(lambda_curr) +
    dgamma(lambda_prop, shape_l, rate_l, log = TRUE) - 
    dgamma(lambda_curr, shape_l, rate_l, log = TRUE) +
    log_lambda_curr - log_lambda_prop  # Jacobian correction
  )

  acceptance_rate <- mean(log_accept_ratio > 0)

  # Adapt step size based on acceptance rate
  if (acceptance_rate < accept_rate_target) {
    sd_prop <- 0.9 * sd_prop  # Shrink the step size
  } else {
    sd_prop <- 1.1 * sd_prop  # Expand the step size
  }

  if (log(runif(1)) < log_accept_ratio) {
    return(list(lambda = lambda_prop, accept = TRUE, sd_prop = sd_prop))
  } else {
    return(list(lambda = lambda_curr, accept = FALSE, sd_prop = sd_prop))
  }
}

# plot template function
gg_template_save <- function(p, filename) {
  p +
    theme_bw(base_size = 25) +
    theme(plot.margin = unit(c(5, 7, 2, 2), "mm"),
          legend.position = "bottom",
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 16),
          legend.key.size = unit(1, "cm"),
          legend.background = element_rect(color = "black", size = 0.5)) -> p2
  ggsave(filename, plot = p2, width = 10, height = 10)
}

# FFBS sampler for states
sample_beta <- function(y, A, mu, Q, H, lambda, taus) {
  Z <- get_Lambda(taus, lambda)
  mod <- SSModel(y ~ -1 + SSMcustom(Z = Z, T = A, Q = Q,
                                    P1 = diag(1e2, 3), a1 = mu), H = diag(H))
  simulateSSM(mod, type = "states", conditional = TRUE) %>%
    matrix(ncol = 3, byrow = TRUE)
}

#####################################################################
# Priors
#####################################################################
# Parameters: mu (3x1), A (3x3), Q (3x3 dense), H (NxN diagonal), lambda (scalar)
# Priors:
#   mu ~ N(m0, V0)
#   vec(A_unconstrained) ~ N(0, sigma_A^2 * I)
#   Q ~ Inverse-Wishart(nu_Q, S_Q)
#   H_diag_i ~ Inverse-Gamma(shape_h, rate_h)
#   lambda ~ Gamma(shape_l, rate_l)

m0       <- c(0.08, -1.5, 0)
V0       <- diag(100, 3)
sigma_A2 <- 0.5
nu_Q     <- 7;  S_Q <- diag(0.05, 3)
shape_h  <- 2;  rate_h <- 0.1
shape_l  <- 4;  rate_l <- 50

#-----------------------------#
# 3. MCMC sampler             #
#-----------------------------#
n_iter      <- 10000; burn <- 500
# storage
mu_chain     <- matrix(0, n_iter, 3)
A_chain      <- array(0, c(n_iter,3,3))
Q_chain      <- array(0, c(n_iter,3,3))
H_chain      <- matrix(0, n_iter, N)
lambda_chain<- numeric(n_iter)
beta_chain   <- array(0, c(n_iter, T, 3))
accept_A     <- accept_lam <- c(0,0)

# init
mu_curr     <- m0
A_curr      <- diag(0.5,3)
Q_curr      <- diag(0.1,3)
H_curr      <- rep(0.1, N)
lambda_curr<- 0.5

pb <- progress_bar$new(total = n_iter, format = "Sampling [:bar] :percent eta: :eta")
for(i in 1:n_iter) {
  pb$tick()
  # (a)
  beta_curr <- sample_beta(y, A_curr, mu_curr, Q_curr, H_curr, lambda_curr, taus)

  # (b)
  b0  <- beta_curr[-1,] - (beta_curr[-T,] - mu_curr) %*% t(A_curr)
  Vn  <- solve(solve(V0) + (T-1)*solve(Q_curr))
  mn  <- Vn %*% (solve(V0)%*%m0 + solve(Q_curr)%*%colSums(b0))
  mu_curr <- as.numeric(rmvnorm(1, mn, Vn))

  # (c)
  outA <- sample_A_unconstr(A_curr, beta_curr, mu_curr, Q_curr)
  A_curr <- outA$A; accept_A <- c(accept_A, outA$accept)

  # (d)
  E    <- beta_curr[-1,] - (mu_curr + (beta_curr[-T,] - mu_curr) %*% t(A_curr))
  S_post <- t(E) %*% E + S_Q; nu_post <- nu_Q + T - 1
  Q_curr <- solve(rWishart(1, nu_post, solve(S_post))[,,1])
  if (max(eigen(Q_curr)$values) > 100) {
    Q_curr <- Q_curr * (100 / max(eigen(Q_curr)$values))
  }

  # (e)
  res  <- y - beta_curr %*% t(get_Lambda(taus, lambda_curr))
  shape_p <- shape_h + T/2; rate_p <- rate_h + colSums(res^2)/2
  H_curr  <- 1 / rgamma(N, shape_p, rate_p)
  H_curr[H_curr > 100] <- 100

  # (f)
  outL <- sample_lambda(lambda_curr, y, beta_curr, H_curr, taus)
  lambda_curr <- outL$lambda; accept_lam <- c(accept_lam, outL$accept)

  # store
  mu_chain[i,]      <- mu_curr
  A_chain[i,,]      <- A_curr
  Q_chain[i,,]      <- Q_curr
  H_chain[i,]       <- H_curr
  lambda_chain[i]   <- lambda_curr
  beta_chain[i,,]   <- beta_curr

}

##################################################################
# Diagnostics
##################################################################

# prepare trace data
df_trace <- data.frame(
  iter = 1:n_iter,
  lambda = lambda_chain,
  a11    = A_chain[,1,1],
  h11    = H_chain[,1],
  q11    = Q_chain[,1,1]
  mu1 = mu_chain[,1],
)

# Remove burn-in
df_trace <- df_trace[-(1:burn),]

# trace plots
g1 <- ggplot(df_trace, aes(x=iter, y=lambda)) + geom_line(size=1.5) + labs(x="Iter", y=expression(lambda))
plot1 <- gg_template_save(g1, "figures/trace_lambda.png")

g2 <- ggplot(df_trace, aes(x=iter, y=a11)) + geom_line(size=1.5) + labs(x="Iter", y=expression(A[1,1]))
plot2 <- gg_template_save(g2, "figures/trace_A11.png")

g3 <- ggplot(df_trace, aes(x=iter, y=h11)) + geom_line(size=1.5) + labs(x="Iter", y=expression(H[1,1]))
plot3 <- gg_template_save(g3, "figures/trace_H11.png")

g4 <- ggplot(df_trace, aes(x=iter, y=q11)) + geom_line(size=1.5) + labs(x="Iter", y=expression(Q[1,1]))
plot4 <- gg_template_save(g4, "figures/trace_Q11.png")

g5 <- ggplot(df_trace, aes(x=iter, y=mu1)) + geom_line(size=1.5) + labs(x="Iter", y=expression(mu1[1]))
plot5 <- gg_template_save(g5, "figures/trace_mu1.png")

# Plot autocorrelations
acf_lambda <- acf(lambda_chain[-(1:burn)], plot = FALSE)
acf_A11    <- acf(A_chain[-(1:burn),1,1], plot = FALSE)
acf_H11    <- acf(H_chain[-(1:burn),1], plot = FALSE)
acf_Q11    <- acf(Q_chain[-(1:burn),1,1], plot = FALSE)
acf_mu1    <- acf(mu_chain[-(1:burn),1], plot = FALSE)

g1 <- ggplot(data.frame(lag = acf_lambda$lag, acf = acf_lambda$acf), aes(x=lag, y=acf)) +
  geom_bar(stat="identity") + labs(x="Lag", y=expression(lambda)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
plot5 <- gg_template_save(g1, "figures/acf_lambda.png")

g2 <- ggplot(data.frame(lag = acf_A11$lag, acf = acf_A11$acf), aes(x=lag, y=acf)) +
  geom_bar(stat="identity") + labs(x="Lag", y=expression(A[1,1])) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
plot6 <- gg_template_save(g2, "figures/acf_A11.png")

g3 <- ggplot(data.frame(lag = acf_H11$lag, acf = acf_H11$acf), aes(x=lag, y=acf)) +
  geom_bar(stat="identity") + labs(x="Lag", y=expression(H[1,1])) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
plot7 <- gg_template_save(g3, "figures/acf_H11.png")

g4 <- ggplot(data.frame(lag = acf_Q11$lag, acf = acf_Q11$acf), aes(x=lag, y=acf)) +
  geom_bar(stat="identity") + labs(x="Lag", y=expression(Q[1,1])) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
plot8 <- gg_template_save(g4, "figures/acf_Q11.png")

g5 <- ggplot(data.frame(lag = acf_mu1$lag, acf = acf_mu1$acf), aes(x=lag, y=acf)) +
  geom_bar(stat="identity") + labs(x="Lag", y=expression(mu[1])) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
plot9 <- gg_template_save(g5, "figures/acf_mu1.png")

# Create tex table for acceptance rates and effective sample size
effectiveSize <- function(x) {
  n <- length(x)
  acf_x <- acf(x, plot = FALSE)$acf[-1]
  var_x <- var(x)
  ess <- n / (1 + 2 * sum(acf_x) / var_x)
  return(ess)
}

acceptance_rates <- data.frame(
  Parameter = c("A", "lambda"),
  Acceptance_Rate = c(mean(accept_A[-(1:burn)]), mean(accept_lam[-(1:burn)])),
  ESS = c(effectiveSize(A_chain[-(1:burn),1,1]), effectiveSize(lambda_chain[-(1:burn)]))
)

stargazer(acceptance_rates, type = "latex", summary = FALSE, 
title = "Acceptance Rates and Effective Sample Size",
        digits = 3, out = "tables/acceptance_rates.tex")

# true vs estimated state
beta_mean <- apply(beta_chain[-(1:burn),,], 2:3, mean)
beta_lo   <- apply(beta_chain[-(1:burn),,], 2:3, quantile, .025)
beta_hi   <- apply(beta_chain[-(1:burn),,], 2:3, quantile, .975)

df_state <- data.frame(
  time = 1:T,
  true  = beta_true[,1],
  est   = beta_mean[,1],
  lo    = beta_lo[,1],
  hi    = beta_hi[,1]
)
p2 <- ggplot(df_state, aes(x=time)) +
  geom_line(aes(y=true, color="True"), size=1.5) +
  geom_line(aes(y=est,  color="Est"), size=1.5) +
  geom_ribbon(aes(ymin=lo, ymax=hi), alpha=0.2) +
  scale_color_manual(values=c("True"="#1E88E5","Est"="#D81B60")) +
  labs(x="Time", y="State 1", color="Series")
plot2 <- gg_template_save(p2, "figures/state_compare.png")

#####################################################################
# Posterior distributions
#####################################################################

for(r in 1:3) for(c in 1:3) {
  var <- A_chain[-(1:burn),r,c]
  dfA <- data.frame(val=var)
  p <- ggplot(dfA, aes(x=val)) + geom_histogram(aes(y=..density..), bins=50) +
       labs(x=paste0("A[",r,",",c,"]"), y="Density") + 
        geom_vline(xintercept=A_true[r,c], color="red", linetype="dashed", size=1)
  gg_template_save(p, paste0("figures/post_A",r,c,".png"))
}

#########################################################################
# Forecasting
#########################################################################

# Get A, Q, H, lambda from posterior means
A_mean     <- apply(A_chain[-(1:burn),,], 2:3, mean)
Q_mean     <- apply(Q_chain[-(1:burn),,], 2:3, mean)
H_mean     <- apply(H_chain[-(1:burn),], 2, mean)
lambda_mean <- mean(lambda_chain[-(1:burn)])
mu_mean     <- apply(mu_chain[-(1:burn),], 2, mean)

# Forecast horizon
h <- 1

# Create a new model with the estimated parameters
Z_forecast <- get_Lambda(taus, lambda_mean)

# Forecast the next h steps with 95% prediction intervals
model_forecast <- SSModel(y ~ -1 + SSMcustom(Z = Z_forecast, T = A_mean, Q = Q_mean,
                                             P1 = diag(0.01, 3), a1 = mu_mean), H = diag(H_mean))
pred <- predict(model_forecast, interval = "prediction", level = 0.95, nsim = 100)

pred <- as.data.frame(pred) 
reshaped_pred <- pred %>%
  pivot_longer(cols = everything(), 
               names_to = c("time", "measure"), 
               names_pattern = "y_(\\d+)\\.(fit|lwr|upr)", 
               values_to = "value") %>%
  group_by(time, measure) %>%
  summarize(mean_value = mean(value), .groups = "drop") %>%
  pivot_wider(names_from = "measure", values_from = "mean_value") %>%
  mutate(time = as.numeric(time)) %>%
  arrange(time)
reshaped_pred <- as.data.frame(reshaped_pred)

# Plot the forecasted values with 95% prediction intervals
gg <- ggplot(reshaped_pred, aes(x = time)) +
  geom_line(aes(y = fit), color = "blue", size = 1.5) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2) +
  labs(x = "Time", y = "Forecasted Values") +
  theme_bw(base_size = 25) +
  theme(plot.margin = unit(c(5, 7, 2, 2), "mm"),
        legend.position = "bottom",
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 16),
        legend.key.size = unit(1, "cm"),
        legend.background = element_rect(color = "black", size = 0.5))
gg_template_save(gg, "figures/forecast.png")

# Plot the forecasted values without
gg <- ggplot(reshaped_pred, aes(x = time)) +
  geom_line(aes(y = fit), color = "blue", size = 1.5) +
  labs(x = "Time", y = "Forecasted Values") +
  theme_bw(base_size = 25) +
  theme(plot.margin = unit(c(5, 7, 2, 2), "mm"),
        legend.position = "bottom",
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 16),
        legend.key.size = unit(1, "cm"),
        legend.background = element_rect(color = "black", size = 0.5))
gg_template_save(gg, "figures/forecast_no_intervals.png")
