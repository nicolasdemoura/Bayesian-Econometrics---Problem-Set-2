#####################################################################
# Topic: Bayesian Econometrics - Problem Set 2
# Instructor: Eduardo Mendes
# Course: Bayesian Econometrics
# Author: Nícolas de Moura
# Goal: Implement a Bayesian model to estimate the parameters of a non-linear model with Time-Varying Coefficients
#####################################################################
# Organize the working environment
#####################################################################

# Clean the working environment
rm(list = ls())
load.lib <- c("dplyr", "KFAS", "MASS", "progress", "mvtnorm", "ggplot2")
install.lib <- load.lib[!load.lib %in% installed.packages()]
for(lib in install.lib) install.packages(lib, dependencies = TRUE)
sapply(load.lib, require, character.only = TRUE)

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
lambda_true <- 1
mu_true <- c(20, -2, -1)
A_true <- matrix(c(0.4, 0.2, -0.1, 
                   -0.1, 0.4, 0.2, 
                   -0.1, -0.1, 0.2), nrow = 3)
H_true <- diag(taus)
Q_true <- diag(10, 3)

# Simulate
ds <- simulate_data(taus, n_obs = 100, lambda_true, mu_true, A_true, H_true, Q_true)
y <- as.matrix(ds$y)
n_obs <- nrow(y)

#####################################################################
# MH-within-Gibbs Sampler
#####################################################################

    # Define update function
    update_fn <- function(pars, model) {
        # Transition matrix A via scaled eigenvalue transformation
        A_unconstrained <- matrix(pars[1:9], 3, 3)
        eigs <- eigen(A_unconstrained)
        vals <- eigs$values
        vals <- 0.98 * vals / max(Mod(vals))    # rescale to ensure stationarity
        A_stable <- eigs$vectors %*% diag(vals) %*% solve(eigs$vectors)
        storage.mode(A_stable) <- "double"
        model$T[,,1] <- A_stable

        # Q = L %*% t(L), L lower triangular
        L <- matrix(0, 3, 3)
        L[1,1] <- pars[10]
        L[2,1] <- pars[11]; L[2,2] <- pars[12]
        L[3,1] <- pars[13]; L[3,2] <- pars[14]; L[3,3] <- pars[15]
        Q_est <- L %*% t(L)
        storage.mode(Q_est) <- "double"
        model$Q[,,1] <- Q_est

        # H diagonal with positive entries
        H_diag <- exp(pars[16:(15 + n_vars)])
        H_mat <- diag(H_diag)
        storage.mode(H_mat) <- "double"
        model$H[,,1] <- H_mat

        return(model)
    }

sample_AQH <- function(lambda, A, Q, H, y, taus) {

    # Observation matrix Z
    Z <- get_Lambda(taus, lambda)

    # Number of parameters to estimate:
    # - A: 9 elements (3x3)
    # - Q: 6 elements (Cholesky factor of 3x3)
    # - H: 20 elements (diagonal, log-transformed)
    n_params <- 9 + 6 + n_vars

    # Build model with placeholder matrices
    A_init <- A 
    Q_init <- Q 
    H_init <- H 

    model <- SSModel(y ~ -1 + SSMcustom(Z = Z, T = A_init, Q = Q_init, P1 = diag(100, 3), a1 = rep(0, 3)), H = H_init)

    # Initial values (some random starting values)
    init_vals <- c(rep(0.1, 9), rep(0.1, 6), rep(log(0.1), n_vars))

    # Estimate model
    fit <- fitSSM(model, inits = init_vals, updatefn = update_fn, method = "BFGS")
    kfs <- KFS(fit$model, smoothing = c("state"))


    return(list(A = fit$model$T[,,1], Q = fit$model$Q[,,1], H = fit$model$H[,,1], beta_hat = kfs$alphahat))
}

sample_lambda <- function(current_lambda, H, beta, taus, y) {
    # Propose from N(current_lambda, 0.0001) and take absolute value
    proposal <- abs(rnorm(1, mean = current_lambda, sd = 0.0001))

    loglike <- function(lambda) {
        Z <- get_Lambda(taus, lambda)
        ll <- 0
        for (t in 2:n_obs) {
            mu_t <- Z %*% beta[t, ]
            ll <- ll + dmvnorm(y[t, ], mean = mu_t, sigma = H, log = TRUE)
        }
        return(ll)
    }

    log_accept_ratio <- loglike(proposal) - loglike(current_lambda)
    if (log(runif(1)) < log_accept_ratio) {
        return(list(lambda=proposal, acceptance = 1))
    } else {
        return(list(lambda=current_lambda, acceptance = 0))
    }
}

mh_gibbs <- function(n_iter = 1000) {
    # Initialize parameter
    lambda_chain <- numeric(n_iter)
    A_chain <- array(0, dim = c(n_iter, 3, 3))
    Q_chain <- array(0, dim = c(n_iter, 3, 3))
    H_chain <- array(0, dim = c(n_iter, n_vars, n_vars))
    beta_chain <- matrix(0, nrow = n_iter, ncol = n_state*n_obs)
    acceptance <- matrix(0, nrow = n_iter, ncol = 1)

    # Initialize the first sample progress bar
    pb <- progress_bar$new(total = n_iter, format = "  Sampling [:bar] :percent eta: :eta")

    # Initialize the first sample
    lambda_curr <- 0.1
    A_curr <- diag(0.9, 3) 
    Q_curr <- diag(0.1, 3) 
    H_curr <- diag(0.1, n_vars) 

    # Store the first sample
    lambda_chain[1] <- lambda_curr
    A_chain[1, , ] <- A_curr
    Q_chain[1, , ] <- Q_curr
    H_chain[1, , ] <- H_curr
    beta_chain[1, ] <- rep(0, n_state*n_obs)
    acceptance[1] <- 0  

    for (i in 2:n_iter) {
        # Update progress bar
        pb$tick()

        # Sample theta = (A, Q, H) given lambda via MCMC block update 
        df <- sample_AQH(lambda_curr, A_curr, Q_curr, H_curr, y, taus)
        A_curr <- df$A
        Q_curr <- df$Q
        H_curr <- df$H
        beta_hat <- df$beta_hat

        # Sample lambda given theta via Metropolis-Hastings
        df <- sample_lambda(lambda_curr, H_curr, beta_hat, taus, y)
        lambda_curr <- df$lambda
        acceptance[i] <- df$acceptance

        # Store the samples
        lambda_chain[i] <- lambda_curr
        A_chain[i, , ] <- A_curr
        Q_chain[i, , ] <- Q_curr
        H_chain[i, , ] <- H_curr
        beta_chain[i, ] <- beta_hat
    }

    return(list(lambda_chain = lambda_chain, A_chain = A_chain, Q_chain = Q_chain, H_chain = H_chain, beta_chain = beta_chain, acceptance = acceptance))
}

n_state <- 3
n_vars <- length(taus)
n_obs <- nrow(y)
n_iter <- 100
n_burnin <- 2

# Run the MCMC sampler
mh_results <- mh_gibbs(n_iter = n_iter)

lambda_chain <- mh_results$lambda_chain[n_burnin:n_iter]
A_chain <- mh_results$A_chain[n_burnin:n_iter, , ]
Q_chain <- mh_results$Q_chain[n_burnin:n_iter, , ]
H_chain <- mh_results$H_chain[n_burnin:n_iter, , ]
beta_chain <- mh_results$beta_chain[n_burnin:n_iter, ]

#################################################################
# MCMC Diagnostics
#################################################################

# Ensure beta_chain and ds$beta have the correct dimensions
estimated_state <- beta_chain[nrow(beta_chain), 1:n_obs]
true_state <- ds$beta[1:n_obs, 1]

# Plot the last state estimate against the true state
ggplot() +
    geom_line(aes(x = 1:n_obs, y = estimated_state), color = "blue") +
    geom_line(aes(x = 1:n_obs, y = true_state), color = "red") +
    labs(title = "Estimated vs True State", x = "Time", y = "State") +
    theme_minimal() +
    scale_color_manual(values = c("blue", "red"), labels = c("Estimated State", "True State")) +
    theme(legend.title = element_blank())

# Plot the traceplots for the parameters lambda, A[1,1], Q[1,1], H[1,1] 

ggplot() +
    geom_line(aes(x = 1:length(lambda_chain), y = lambda_chain), color = "blue") +
    labs(title = "Traceplot of lambda", x = "Iteration", y = "lambda") +
    theme_minimal()

ggplot() +
    geom_line(aes(x = 1:length(A_chain[, 1, 1]), y = A_chain[, 1, 1]), color = "blue") +
    labs(title = "Traceplot of A[1,1]", x = "Iteration", y = "A[1,1]") +
    theme_minimal()

ggplot() +
    geom_line(aes(x = 1:length(Q_chain[, 1, 1]), y = Q_chain[, 1, 1]), color = "blue") +
    labs(title = "Traceplot of Q[1,1]", x = "Iteration", y = "Q[1,1]") +
    theme_minimal()

ggplot() +
    geom_line(aes(x = 1:length(H_chain[, 1, 1]), y = H_chain[, 1, 1]), color = "blue") +
    labs(title = "Traceplot of H[1,1]", x = "Iteration", y = "H[1,1]") +
    theme_minimal()

################################################################################
# Plotting the posterior distributions of A, Q, H
################################################################################


###########################################################################
# Predictive curve 
###########################################################################

# Estimate lambda, A, Q, H from the mean of the chains
lambda_est <- mean(lambda_chain)
A_est <- apply(A_chain, c(2, 3), mean)
Q_est <- apply(Q_chain, c(2, 3), mean)
H_est <- apply(H_chain, c(2, 3), mean)

# Get the mean beta from the last observation
beta_est <- apply(beta_chain, 2, mean)[(n_state * (n_obs - 1) + 1):(n_state * n_obs)]
beta_est <- matrix(beta_est, nrow = n_obs, ncol = n_state, byrow = TRUE)

# Nelson-Siegel Model
nelson_siegel <- function(tau, beta, lambda){
    # Tau : Time to maturity
    # beta[1] : Beta_0 - Level
    # beta[2] : Beta_1 - Slope
    # beta[3] : Beta_2 - Curvature
    # lambda : lambda - Factor of decay
    term1 <- (1 - exp(-tau / lambda)) / (tau / lambda)
    term2 <- term1 - exp(-tau / lambda)

    # Compute the Nelson-Siegel model
    y <- beta[1] + beta[2] * term1 + beta[3] * term2
    return(y)
}