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
load.lib <- c("dplyr", "KFAS", "MASS")
install.lib <- load.lib[!load.lib %in% installed.packages()]
for(lib in install.lib) install.packages(lib, dependencies = TRUE)
sapply(load.lib, require, character.only = TRUE)

set.seed(20250420)

#####################################################################
# Simulate the data
#####################################################################

simulate_data <- function(taus, N, lambda, mu, A, H, Q) {
    Lambda <- matrix(0, nrow = length(taus), ncol = 3)
    for (i in 1:length(taus)) {
        Lambda[i, ] <- c(
            1,
            (1 - exp(-lambda * taus[i])) / (lambda * taus[i]),
            ((1 - exp(-lambda * taus[i])) / (lambda * taus[i])) - exp(-lambda * taus[i])
        )
    }

    beta <- matrix(0, nrow = N, ncol = 3)
    y <- matrix(0, nrow = N, ncol = length(taus))

    for (i in 2:N) {
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
lambda_true <- 0.1
mu_true <- c(20, -2, -1)
A_true <- matrix(c(0.5, 0.1, 0.1, 
                   -0.1, 0.5, -0.1, 
                   -0.1, -0.1, 0.5), nrow = 3)
H_true <- diag(sqrt(taus) * 0.1)
Q_true <- diag(10, 3)

# Simulate
ds <- simulate_data(taus, 100, lambda_true, mu_true, A_true, H_true, Q_true)
y <- as.matrix(ds$y)
N <- nrow(y)

#####################################################################
# Model the data with the Kalman Filter
#####################################################################

n_state <- 3
n_obs <- length(taus)

# Observation matrix Z
Z <- matrix(0, nrow = n_obs, ncol = n_state)
for (i in 1:n_obs) {
    Z[i, ] <- c(
        1,
        (1 - exp(-lambda_true * taus[i])) / (lambda_true * taus[i]),
        ((1 - exp(-lambda_true * taus[i])) / (lambda_true * taus[i])) - exp(-lambda_true * taus[i])
    )
}

# Define zero matrices with appropriate dimensions
zero_state_obs <- matrix(0, nrow = n_obs, ncol = n_state)
zero_state_state <- matrix(0, nrow = n_state, ncol = n_state)

# Combine matrices to form the full Z matrix
Z_tilde <- rbind(cbind(Z, Z))

# Number of parameters to estimate:
# - A: 9 elements (3x3)
# - Q: 6 elements (Cholesky factor of 3x3)
# - H: 20 elements (diagonal, log-transformed)
n_params <- 9 + 6 + n_obs

# Build model with placeholder matrices
A_init <- diag(0.9, 3)
Q_init <- diag(0.1, 3)
H_init <- diag(0.1, n_obs)

A_tilde_init <- rbind(cbind(A_init, zero_state_state), cbind(zero_state_state, zero_state_state))
Q_tilde_init <- rbind(cbind(Q_init, zero_state_state), cbind(zero_state_state, zero_state_state))
H_tilde_init <- H_init

model <- SSModel(y ~ -1 + SSMcustom(Z = Z_tilde, T = A_tilde_init, Q = Q_tilde_init, P1 = diag(100, n_state*2), a1 = rep(0, n_state*2)), H = H_tilde_init)

# Define update function
update_fn <- function(pars, model) {
    # Transition matrix A via scaled eigenvalue transformation
    A_unconstrained <- matrix(pars[1:9], 3, 3)
    eigs <- eigen(A_unconstrained)
    vals <- eigs$values
    vals <- 0.98 * vals / max(Mod(vals))    # rescale to ensure stationarity
    A_stable <- eigs$vectors %*% diag(vals) %*% solve(eigs$vectors)
    A_tilde_stable <- rbind(cbind(A_stable, zero_state_state), cbind(zero_state_state, zero_state_state))
    storage.mode(A_tilde_stable) <- "double"
    model$T[,,1] <- A_tilde_stable

    # Q = L %*% t(L), L lower triangular
    L <- matrix(0, 3, 3)
    L[1,1] <- pars[10]
    L[2,1] <- pars[11]; L[2,2] <- pars[12]
    L[3,1] <- pars[13]; L[3,2] <- pars[14]; L[3,3] <- pars[15]
    Q_est <- L %*% t(L)
    storage.mode(Q_est) <- "double"
    Q_tilde_est <- rbind(cbind(Q_est, zero_state_state), cbind(zero_state_state, zero_state_state))
    model$Q[,,1] <- Q_tilde_est

    # H diagonal with positive entries
    H_diag <- exp(pars[16:(15 + n_obs)])
    H_mat <- diag(H_diag)
    storage.mode(H_mat) <- "double"
    H_tilde_mat <- H_mat
    model$H[,,1] <- H_tilde_mat

    return(model)
}

# Initial values (some random starting values)
init_vals <- c(rep(0.1, 9), rep(0.1, 6), rep(log(0.1), n_obs))

# Estimate model
fit <- fitSSM(model, inits = init_vals, updatefn = update_fn, method = "BFGS")

kfs <- KFS(fit$model, smoothing = c("state", "disturbance"))

# Plot the smoothed state estimates against the true values
ggplot2::ggplot() +
    geom_line(aes(x = 1:N, y = kfs$alphahat[, 1]+kfs$alphahat[, 4]), color = "blue") +
    geom_line(aes(x = 1:N, y = ds$beta[, 1]), color = "red") +
    labs(title = "Smoothed State Estimates vs True Values",
         x = "Time",
         y = "State Estimate") +
    theme_minimal()
