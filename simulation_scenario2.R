library(SoftBart)
library(truncnorm)
library(MASS)

# ==============================================================================
# 1. CORE FUNCTIONS (DGP & UTILITIES)
# ==============================================================================

b2_func_true <- function(x_vals, z_val, w_val) {
  x1 <- x_vals[1]; x2 <- x_vals[2]; x3 <- x_vals[3]; x4 <- x_vals[4]; x5 <- x_vals[5]
  main_effect <- 1.0 * z_val + 0.7 * x1 - 0.5 * x2 + 0.3 * (z_val * x1 * exp(x2)) +
    0.5 * x3 - 0.2 * x4 + 0.1 * x5 + 0.4 * x1 * x2 * exp(x2) + 0.3 * x2 * x4 +
    0.2 * x1 * x3 * x4 + 0.35 * z_val * x2 * x3 + 0.25 * z_val * log(x1) * x4 +
    0.15 * z_val * x3 * x4 * x5 / 10 + 0.5 * sin(pi * x5 / 5) + 0.3 * sin(pi * x1 * x3) +
    0.2 * z_val * sin(pi * x2 * x4) + 0.5 * z_val * exp(x1) * log(x2)
  covariate_w_interaction <- 0.3 * x1 * w_val + 0.2 * x2 * w_val + 0.15 * x3 * w_val + 
    0.1 * (x5 - 5.5) * w_val + 0.25 * z_val * x1 * w_val
  return(main_effect + covariate_w_interaction)
}

scale_to_unit <- function(X) {
  X <- as.matrix(X)
  mins <- apply(X, 2, min); maxs <- apply(X, 2, max)
  ranges <- maxs - mins; ranges[ranges == 0] <- 1 
  return(sweep(sweep(X, 2, mins), 2, ranges, "/"))
}

create_dagar_precision <- function(adj_matrix, rho) {
  K <- nrow(adj_matrix); L <- diag(K); n_pi <- integer(K); B <- matrix(0, K, K)
  for (i in 1:K) {
    neighbors <- which(adj_matrix[i, ] == 1)
    dir_nb <- neighbors[neighbors < i]
    n_pi[i] <- length(dir_nb)
    if (n_pi[i] > 0) B[i, dir_nb] <- rho / (1 + (n_pi[i] - 1) * rho^2)
  }
  tau_i <- (1 + (n_pi - 1) * rho^2) / (1 - rho^2)
  L <- diag(K) - B
  return(t(L) %*% diag(tau_i) %*% L)
}

# ==============================================================================
# 2. THE REPLICATION WRAPPER
# ==============================================================================

run_replicate <- function(rep_id) {
  set.seed(2025 + rep_id)
  
  # --- Setup Grid & Data ---
  K <- 35; adj <- matrix(0, K, K)
  for(i in 1:K) {
    r <- ceiling(i/7); c <- i - (r-1)*7
    if(r > 1) adj[i, i-7] <- 1; if(r < 5) adj[i, i+7] <- 1
    if(c > 1) adj[i, i-1] <- 1; if(c < 7) adj[i, i+1] <- 1
  }
  
  W_true <- as.numeric(mvrnorm(1, rep(0, K), solve(create_dagar_precision(adj, 0.7))))
  cluster_id <- rep(1:K, each = 50); n <- length(cluster_id)
  X1 <- runif(n, 0, 1)
  X2 <- rbeta(n, 2, 6)
  X3 <- rbinom(n, 1, 0.8)
  X4 <- rbinom(n, 1, 0.7)
  X5 <- sample(1:10, n, replace = TRUE)
  X <- cbind(X1, X2, X3, X4, X5)
  b1_lin <- 0.5 * X1 - 0.3 * X2 + 0.2 * X3 - 0.1 * X4 + 0.05 * X2*(X5 - 5.5)+ X1*X2
  e_true <- pnorm(b1_lin)
  Z <- rbinom(n, 1, prob = e_true)
  
  # Survival Times
  logT_true <- sapply(1:n, function(i) b2_func_true(X[i,], Z[i], W_true[cluster_id[i]])) + rnorm(n, 0, sqrt(0.15))
  T_true <- exp(logT_true)
  
  # --- Introduce Censoring ---
  C_time <- rexp(n, rate = 0.1) # Adjust rate for desired censoring %
  y_obs_log <- pmin(logT_true, log(C_time))
  status <- as.integer(logT_true <= log(C_time))
  # --- Calculate and Print Censoring Percentage ---
  censoring_pct <- mean(status == 0) * 100
  cat("Replicate", rep_id, "- Censoring Rate:", round(censoring_pct, 2), "%\n")
  cens_idx <- which(status == 0)
  
  # Initial Scaling
  y_mean <- mean(y_obs_log); y_sd <- sd(y_obs_log)
  y_latent_log <- y_obs_log # Initialize latent vector with observed
  
  # --- Estimation ---
  X_s <- scale_to_unit(X)
  f_ps <- MakeForest(Hypers(X_s, qnorm(pmin(pmax(mean(Z), 0.1), 0.9)), sigma_hat=1), Opts(update_sigma=F), warn=F)
  for(i in 1:400) f_ps$do_gibbs(X_s, Z - 0.5, X_s, 1)
  e_hat <- pnorm(as.numeric(f_ps$do_predict(X_s)))
  
  X_mu <- scale_to_unit(cbind(X, e_hat)); X_tau <- scale_to_unit(X)
  f_mu <- MakeForest(Hypers(X_mu, (y_latent_log - y_mean)/y_sd, sigma_hat=1), Opts(), warn=F)
  f_tau <- MakeForest(Hypers(X_tau, rep(0, n), sigma_hat=0.5), Opts(), warn=F)
  
  W <- rep(0, K); sigma2 <- 1; Q <- create_dagar_precision(adj, 0.7)
  n_mcmc <- 1000; burn_in <- 500; thin <- 2
  n_save <- (n_mcmc - burn_in) / thin
  
  target_combos <- data.frame(county=c(4,4,10,10), t_star=c(5,15,5,15))
  crate_draws <- matrix(0, nrow=n_save, ncol=4)
  
  for (iter in 1:n_mcmc) {
    # Current scaled predictions
    mu_p_std <- as.numeric(f_mu$do_predict(X_mu))
    tau_p_std <- as.numeric(f_tau$do_predict(X_tau))
    
    # --- DATA AUGMENTATION STEP ---
    # Draw latent log-survival for censored observations
    if(length(cens_idx) > 0) {
      # Current mean on original log-scale
      cur_mu_log <- (mu_p_std * y_sd) + y_mean + (tau_p_std * y_sd) * Z + W[cluster_id]
      cur_sd_log <- sqrt(sigma2) * y_sd
      
      y_latent_log[cens_idx] <- rtruncnorm(length(cens_idx), 
                                           a = log(C_time[cens_idx]), 
                                           mean = cur_mu_log[cens_idx], 
                                           sd = cur_sd_log)
    }
    
    # Standardize updated latent outcome
    y_std <- (y_latent_log - y_mean) / y_sd
    
    # --- Spatial & Forest Updates ---
    for(i in 1:K) {
      idx <- which(cluster_id == i)
      prec_i <- length(idx)/sigma2 + Q[i,i]
      rem <- sum(y_std[idx] - mu_p_std[idx] - tau_p_std[idx]*Z[idx])/sigma2 - sum(Q[i,-i]*W[-i])
      W[i] <- rnorm(1, rem/prec_i, sqrt(1/prec_i))
    }
    
    f_mu$do_gibbs(X_mu, y_std - tau_p_std*Z - W[cluster_id], X_mu, 1)
    f_tau$do_gibbs(X_tau[Z==1,], (y_std - mu_p_std - W[cluster_id])[Z==1], X_tau, 1)
    sigma2 <- 1/rgamma(1, 1 + n/2, 1 + 0.5*sum((y_std - mu_p_std - tau_p_std*Z - W[cluster_id])^2))
    
    # --- Posterior Draws ---
    if (iter > burn_in && (iter-burn_in) %% thin == 0) {
      s_idx <- (iter-burn_in)/thin
      for(d in 1:4) {
        c_id <- target_combos$county[d]
        t_s <- target_combos$t_star[d]
        m_r <- (mu_p_std[cluster_id == c_id] * y_sd) + y_mean
        t_r <- (tau_p_std[cluster_id == c_id] * y_sd)
        w_r <- W[c_id] * y_sd
        s_r <- sqrt(sigma2) * y_sd
        s1 <- pnorm((log(t_s) - (m_r + t_r + w_r)) / s_r, lower.tail=F)
        s0 <- pnorm((log(t_s) - (m_r + w_r)) / s_r, lower.tail=F)
        crate_draws[s_idx, d] <- mean(s1 - s0)
      }
    }
  }
  
  # --- Summary Statistics ---
  rep_results <- data.frame()
  for(d in 1:4) {
    c_id <- target_combos$county[d]
    t_s <- target_combos$t_star[d]
    # True CRATE is counterfactual (independent of censoring)
    tr_s1 <- pnorm((log(t_s) - sapply(which(cluster_id==c_id), function(j) b2_func_true(X[j,], 1, W_true[c_id]))) / sqrt(0.15), lower.tail=F)
    tr_s0 <- pnorm((log(t_s) - sapply(which(cluster_id==c_id), function(j) b2_func_true(X[j,], 0, W_true[c_id]))) / sqrt(0.15), lower.tail=F)
    true_crate <- mean(tr_s1 - tr_s0)
    
    est_mean <- mean(crate_draws[,d]); ci <- quantile(crate_draws[,d], c(0.025, 0.975))
    rep_results <- rbind(rep_results, data.frame(County = c_id, Time = t_s, MAE = abs(est_mean - true_crate), Covered = (true_crate >= ci[1] & true_crate <= ci[2]), CI_Width = ci[2] - ci[1]))
  }
  return(rep_results)
}

# ==============================================================================
# 3. RUN & AGGREGATE
# ==============================================================================

all_reps <- lapply(1:5, function(r) {
  cat("Starting Replicate:", r, "\n")
  run_replicate(r)
})

final_df <- do.call(rbind, all_reps)
summary_stats <- aggregate(cbind(MAE, Covered, CI_Width) ~ County + Time, data = final_df, FUN = mean)
print(summary_stats)