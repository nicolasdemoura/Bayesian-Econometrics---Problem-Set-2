estimate_kalman_em <- function(y, Lambda, mu, max_iter = 100, tol = 1e-5) {
    N <- nrow(y)
    n_states <- ncol(Lambda)  # 3
    n_obs <- ncol(y)

    # Initial guesses
    A <- diag(0.5, n_states)
    Q <- diag(1, n_states)
    H <- diag(0.1, n_obs)

    loglik_prev <- -Inf

    for (iter in 1:max_iter) {
        # Build state-space model with KFAS
        model <- SSModel(y ~ -1 +
                         SSMcustom(
                             Z = array(rep(Lambda, N), dim = c(n_obs, n_states, N)),
                             T = array(A, dim = c(n_states, n_states, N - 1)),
                             Q = array(Q, dim = c(n_states, n_states, N - 1)),
                             R = diag(n_states),
                             P1 = diag(1, n_states),
                             a1 = rep(0, n_states)
                         ), H = H)

        # Run Kalman smoother
        kfs <- KFS(model, smoothing = c("state", "disturbance"))

        beta_smoothed <- kfs$alphahat

        # E-step: calculate sufficient statistics
        beta_lag <- beta_smoothed[-N, ]
        beta_lead <- beta_smoothed[-1, ]
        beta_centered <- sweep(beta_lag, 2, mu)

        # M-step: update A
        A_new <- t(beta_lead - matrix(rep(mu, N - 1), ncol = n_states, byrow = TRUE)) %*% beta_centered
        A_new <- A_new %*% solve(t(beta_centered) %*% beta_centered)

        # Update Q
        eta_hat <- beta_lead - matrix(rep(mu, N - 1), ncol = n_states, byrow = TRUE) - (beta_centered %*% t(A_new))
        Q_new <- t(eta_hat) %*% eta_hat / (N - 1)

        # Update H
        eps_hat <- y - beta_smoothed %*% t(Lambda)
        H_new <- t(eps_hat) %*% eps_hat / N

        # Check for convergence
        loglik <- logLik(model)
        if (abs(loglik - loglik_prev) < tol) break
        loglik_prev <- loglik

        # Update parameters
        A <- A_new
        Q <- Q_new
        H <- H_new
    }

    return(list(A = A, Q = Q, H = H, beta_smoothed = beta_smoothed))
}

Lambda <- matrix(0, nrow = length(taus), ncol = 3)
for (i in 1:length(taus)) {
    Lambda[i, ] <- c(
        1,
        (1 - exp(-lambda_true * taus[i])) / (lambda_true * taus[i]),
        ((1 - exp(-lambda_true * taus[i])) / (lambda_true * taus[i])) - exp(-lambda_true * taus[i])
    )
}

result <- estimate_kalman_em(y = y, Lambda = Lambda, mu = mu_true)

# Extract estimates
A_est <- result$A
Q_est <- result$Q
H_est <- result$H
