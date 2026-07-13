library(SoftBart)
library(truncnorm)
library(MASS)

# ==============================================================================
# 1. CORE COVARIATE SCALING & DAGAR UTILITIES
# ==============================================================================

scale_to_unit <- function(X) {
  X <- as.matrix(X)
  mins <- apply(X, 2, min)
  maxs <- apply(X, 2, max)
  ranges <- maxs - mins
  ranges[ranges == 0] <- 1  
  X_scaled <- sweep(sweep(X, 2, mins), 2, ranges, "/")
  list(X_scaled = X_scaled, mins = mins, ranges = ranges)
}

scale_new <- function(x, params) {
  x <- as.matrix(x)
  if (ncol(x) != length(params$mins)) {
    stop("Dimension mismatch: x has ", ncol(x), " columns but params has ", length(params$mins))
  }
  sweep(sweep(x, 2, params$mins), 2, params$ranges, "/")
}

create_dagar_precision <- function(adj_matrix, ordering, rho) {
  K <- nrow(adj_matrix)
  if (rho < 0 || rho >= 1) stop("rho must be in [0, 1)")
  
  N_pi <- vector("list", K)
  n_pi <- integer(K)
  
  for (i in 1:K) {
    node <- ordering[i]
    neighbors <- unname(which(adj_matrix[node, ] == 1))
    if (length(neighbors) == 0) {
      N_pi[[node]] <- integer(0)
      n_pi[node] <- 0
    } else {
      filter_mask <- sapply(neighbors, function(j) which(ordering == j) < i)
      directed_neighbors <- neighbors[filter_mask]
      N_pi[[node]] <- directed_neighbors
      n_pi[node] <- length(directed_neighbors)
    }
  }
  
  B <- matrix(0, K, K)
  for (i in 1:K) {
    node <- ordering[i]
    if (n_pi[node] > 0) {
      B[node, N_pi[[node]]] <- rho / (1 + (n_pi[node] - 1) * rho^2)
    }
  }
  
  tau_i <- (1 + (n_pi - 1) * rho^2) / (1 - rho^2)
  F_mat <- diag(tau_i)
  L <- diag(K) - B
  Q <- t(L) %*% F_mat %*% L
  
  list(Q = Q, B = B, F_mat = F_mat, L = L, N_pi = N_pi, n_pi = n_pi, ordering = ordering)
}

# ==============================================================================
# 2. STAGE 1: PROPENSITY SCORE MODEL WITH ADAPTIVE DAGAR UPDATES
# ==============================================================================

fit_ps_dagar <- function(X, Z, cluster_id, adj_matrix, ordering = NULL,
                         n_mcmc = 1000, burn_in = 500, thin = 2, num_tree = 50,
                         rho_ps_init = 0.3, rho_ps_prop_sd = 0.05,
                         a_tau_ps = 1, b_tau_ps = 1, a_rho_ps = 2, b_rho_ps = 2) {
  
  n <- length(Z)
  unique_clusters <- sort(unique(cluster_id))
  K <- length(unique_clusters)
  cluster_indices <- split(1:n, cluster_id)
  n_i_vec <- sapply(cluster_indices, length)
  
  if (is.null(ordering)) ordering <- 1:K
  
  scaling_results <- scale_to_unit(X)
  X_scaled <- scaling_results$X_scaled
  X_scaled[is.nan(X_scaled)] <- 0
  X_scaled <- pmax(pmin(X_scaled, 10), -10)
  
  y_init <- qnorm(pmax(0.01, pmin(0.99, mean(Z) + 0.1 * (Z - mean(Z)))))
  hypers <- Hypers(X = X_scaled, Y = y_init, k = 2, num_tree = num_tree, sigma_hat = 1)
  opts <- Opts(update_sigma = FALSE, cache_trees = TRUE)
  forest <- MakeForest(hypers, opts, warn = FALSE)
  b1_current <- as.numeric(forest$do_predict(X_scaled))
  
  Z_star <- numeric(n)
  V <- rep(0, K)
  rho_ps <- rho_ps_init
  tau_v <- 1
  Q_ps <- create_dagar_precision(adj_matrix, ordering, rho_ps)$Q
  
  n_save <- floor((n_mcmc - burn_in) / thin)
  ps_samples <- matrix(0, nrow = n_save, ncol = n)
  accept_count_rho <- 0
  sample_idx <- 1
  
  for (iter in 1:n_mcmc) {
    V_expanded <- V[cluster_id]
    mean_vec <- b1_current + V_expanded
    Z_star <- rtruncnorm(n, a = ifelse(Z==1, 0, -Inf), b = ifelse(Z==1, Inf, 0), mean = mean_vec, sd = 1)
    
    b_vec_v <- sapply(cluster_indices, function(idx) sum(Z_star[idx] - b1_current[idx]))
    Prec_V <- diag(n_i_vec) + tau_v * Q_ps
    Cov_V <- solve(Prec_V)
    V <- as.numeric(mvrnorm(1, mu = Cov_V %*% b_vec_v, Sigma = Cov_V))
    
    tau_v <- rgamma(1, shape = a_tau_ps + K/2, rate = b_tau_ps + 0.5 * as.numeric(t(V) %*% Q_ps %*% V))
    
    # Metropolis-Hastings for Stage 1 rho_ps (Beta(2,2) prior)
    rho_ps_prop <- rho_ps + rnorm(1, 0, rho_ps_prop_sd)
    if (rho_ps_prop > 0 && rho_ps_prop < 1) {
      Q_prop <- create_dagar_precision(adj_matrix, ordering, rho_ps_prop)$Q
      log_alpha <- (0.5 * as.numeric(determinant(Q_prop, logarithm=T)$modulus) - 0.5 * tau_v * t(V) %*% Q_prop %*% V + (a_rho_ps-1)*log(rho_ps_prop) + (b_rho_ps-1)*log(1-rho_ps_prop)) -
        (0.5 * as.numeric(determinant(Q_ps, logarithm=T)$modulus) - 0.5 * tau_v * t(V) %*% Q_ps %*% V + (a_rho_ps-1)*log(rho_ps) + (b_rho_ps-1)*log(1-rho)  ) # Fixed typo here from rho to rho_ps
      if (log(runif(1)) < log_alpha) {
        rho_ps <- rho_ps_prop
        Q_ps <- Q_prop
        accept_count_rho <- accept_count_rho + 1
      }
    }
    
    # Adaptive tuning window for proposal variance
    if (iter <= burn_in && iter %% 50 == 0) {
      acc_rate <- accept_count_rho / 50
      rho_ps_prop_sd <- rho_ps_prop_sd * ifelse(acc_rate < 0.2, 0.8, ifelse(acc_rate > 0.45, 1.2, 1))
      accept_count_rho <- 0
    }
    
    forest$do_gibbs(X_scaled, Z_star - V[cluster_id], X_scaled, 1)
    b1_current <- as.numeric(forest$do_predict(X_scaled))
    
    if (iter > burn_in && ((iter - burn_in) %% thin == 0)) {
      ps_samples[sample_idx, ] <- pnorm(b1_current + V[cluster_id])
      sample_idx <- sample_idx + 1
    }
  }
  return(list(e_hat = colMeans(ps_samples)))
}

# ==============================================================================
# 3. STAGE 2: BCF OUTCOME MODEL WITH DATA AUGMENTATION & DAGAR
# ==============================================================================

fit_outcome_bcf_dagar <- function(time, status, X, Z, e_hat, cluster_id, adj_matrix, ordering = NULL,
                                  n_mcmc = 1000, burn_in = 500, thin = 2,
                                  num_tree_mu = 50, num_tree_tau = 25, k_mu = 2, k_tau = 3,
                                  rho_init = 0.5, rho_prop_sd = 0.05,
                                  a_tau = 1, b_tau = 1, a_sigma = 2, b_sigma = 1, a_rho = 2, b_rho = 2) {
  
  n <- length(time); unique_clusters <- sort(unique(cluster_id)); K <- length(unique_clusters)
  cluster_map <- setNames(1:K, as.character(unique_clusters)); cluster_idx <- cluster_map[as.character(cluster_id)]
  cluster_indices <- split(1:n, cluster_idx); n_i_vec <- sapply(cluster_indices, length)
  if (is.null(ordering)) ordering <- 1:K
  
  y_log <- log(ifelse(time <= 0, min(time[time > 0])/10, time))
  scaling_mu <- scale_to_unit(cbind(X, e_hat)); X_mu_scaled <- scaling_mu$X_scaled
  scaling_tau <- scale_to_unit(as.matrix(X)); X_tau_scaled <- scaling_tau$X_scaled
  
  tr_idx <- which(Z == 1)
  forest_mu <- MakeForest(Hypers(X_mu_scaled, y_log, k = k_mu, num_tree = num_tree_mu, sigma_hat = sd(y_log), normalize_Y = F), Opts(F, T))
  mu_curr <- as.numeric(forest_mu$do_predict(X_mu_scaled))
  forest_tau <- MakeForest(Hypers(X_tau_scaled[tr_idx,], rep(0, length(tr_idx)), k = k_tau, num_tree = num_tree_tau, sigma_hat = sd(y_log)/2, normalize_Y = F), Opts(F, T))
  tau_curr <- as.numeric(forest_tau$do_predict(X_tau_scaled))
  
  sigma2 <- var(y_log); tau_w <- 1; rho <- rho_init; W <- rep(0, K)
  Q <- create_dagar_precision(adj_matrix, ordering, rho)$Q
  tilde_y <- y_log; cens_idx <- which(status == 0)
  
  n_save <- floor((n_mcmc - burn_in) / thin)
  mu_samples <- matrix(0, n_save, n); tau_samples <- matrix(0, n_save, n); w_samples <- matrix(0, n_save, K)
  sigma2_samples <- rho_samples <- numeric(n_save); iter_indices <- integer(n_save)
  sample_idx <- 1; accept_count_rho <- 0
  
  for (iter in 1:n_mcmc) {
    # --- Data Augmentation Loop for Censored Survival Times ---
    mu_full <- mu_curr + tau_curr * Z + W[cluster_idx]
    if(length(cens_idx) > 0) {
      tilde_y[cens_idx] <- rtruncnorm(length(cens_idx), a = y_log[cens_idx], b = Inf, mean = mu_full[cens_idx], sd = sqrt(sigma2))
    }
    
    b_vec <- sapply(1:K, function(i) sum(tilde_y[cluster_indices[[i]]] - mu_curr[cluster_indices[[i]]] - tau_curr[cluster_indices[[i]]] * Z[cluster_indices[[i]]]) / sigma2)
    Prec_W <- diag(n_i_vec)/sigma2 + tau_w * Q; Cov_W <- solve(Prec_W)
    W <- as.numeric(mvrnorm(1, mu = Cov_W %*% b_vec, Sigma = Cov_W))
    
    tau_w <- rgamma(1, a_tau + K/2, b_tau + 0.5 * as.numeric(t(W) %*% Q %*% W))
    
    # Metropolis-Hastings for Stage 2 rho (Beta(2,2) prior)
    rho_prop <- rho + rnorm(1, 0, rho_prop_sd)
    if (rho_prop > 0 && rho_prop < 1) {
      Q_prop <- create_dagar_precision(adj_matrix, ordering, rho_prop)$Q
      log_alpha <- (0.5 * as.numeric(determinant(Q_prop, logarithm=T)$modulus) - 0.5 * tau_w * t(W) %*% Q_prop %*% W + (a_rho-1)*log(rho_prop) + (b_rho-1)*log(1-rho_prop)) -
        (0.5 * as.numeric(determinant(Q, logarithm=T)$modulus) - 0.5 * tau_w * t(W) %*% Q %*% W + (a_rho-1)*log(rho) + (b_rho-1)*log(1-rho))
      if (log(runif(1)) < log_alpha) { rho <- rho_prop; Q <- Q_prop; accept_count_rho <- accept_count_rho + 1 }
    }
    
    if (iter <= burn_in && iter %% 50 == 0) {
      acc_rate <- accept_count_rho / 50
      rho_prop_sd <- rho_prop_sd * ifelse(acc_rate < 0.2, 0.8, ifelse(acc_rate > 0.45, 1.2, 1))
      accept_count_rho <- 0
    }
    
    forest_mu$do_gibbs(X_mu_scaled, tilde_y - tau_curr * Z - W[cluster_idx], X_mu_scaled, 1)
    mu_curr <- as.numeric(forest_mu$do_predict(X_mu_scaled))
    forest_tau$do_gibbs(X_tau_scaled[tr_idx,], tilde_y[tr_idx] - mu_curr[tr_idx] - W[cluster_idx[tr_idx]], X_tau_scaled, 1)
    tau_curr <- as.numeric(forest_tau$do_predict(X_tau_scaled))
    
    sigma2 <- 1/rgamma(1, a_sigma + n/2, b_sigma + 0.5 * sum((tilde_y - mu_curr - tau_curr * Z - W[cluster_idx])^2))
    
    if (iter > burn_in && ((iter - burn_in) %% thin == 0)) {
      mu_samples[sample_idx,] <- mu_curr; tau_samples[sample_idx,] <- tau_curr; w_samples[sample_idx,] <- W
      sigma2_samples[sample_idx] <- sigma2; rho_samples[sample_idx] <- rho; iter_indices[sample_idx] <- iter; sample_idx <- sample_idx + 1
    }
  }
  return(list(forest_mu = forest_mu, forest_tau = forest_tau, mu_samples = mu_samples, tau_samples = tau_samples, 
              w_samples = w_samples, sigma2_samples = sigma2_samples, rho_samples = rho_samples, 
              scaling_mu = scaling_mu, scaling_tau = scaling_tau, cluster_map = cluster_map, iter_indices = iter_indices, n_save = n_save))
}

# ==============================================================================
# 4. DATA-GENERATING TRUE SURFACE WITH COVARIATE-SPATIAL INTERACTION
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

# ==============================================================================
# 5. CAUSAL EFFECT EXTRACTION ALGORITHM (CERM)
# ==============================================================================

estimate_cerm <- function(t_star, x, county_id, e_hat, fit) {
  n_draws <- fit$n_save
  CERM_draws <- numeric(n_draws)
  county_idx <- fit$cluster_map[as.character(county_id)]
  
  x_mu <- matrix(c(x, e_hat), nrow = 1)
  x_mu_scaled <- scale_new(x_mu, fit$scaling_mu)
  x_tau <- matrix(x, nrow = 1)
  x_tau_scaled <- scale_new(x_tau, fit$scaling_tau)
  
  for (m in 1:n_draws) {
    iter_m <- fit$iter_indices[m]
    mu_pred <- fit$forest_mu$predict_iteration(x_mu_scaled, iter = iter_m)
    tau_pred <- fit$forest_tau$predict_iteration(x_tau_scaled, iter = iter_m)
    
    W_m <- fit$w_samples[m, county_idx]
    sigma2_m <- fit$sigma2_samples[m]
    sigma_m <- sqrt(sigma2_m)
    
    mu_1 <- mu_pred + tau_pred + W_m   
    mu_0 <- mu_pred + W_m              
    
    log_t_star <- log(t_star)
    z_1 <- (log_t_star - mu_1) / sigma_m
    z_0 <- (log_t_star - mu_0) / sigma_m
    
    RM_1 <- exp(pmin(mu_1 + sigma2_m / 2, 700)) * pnorm(z_1 - sigma_m) + t_star * (1 - pnorm(z_1))
    RM_0 <- exp(pmin(mu_0 + sigma2_m / 2, 700)) * pnorm(z_0 - sigma_m) + t_star * (1 - pnorm(z_0))
    
    CERM_draws[m] <- pmin(RM_1, t_star) - pmin(RM_0, t_star)
  }
  return(CERM_draws)
}

# ==============================================================================
# 6. REPLICATION WRAPPER ENGINE
# ==============================================================================

run_replicate <- function(rep_id) {
  set.seed(2025 + rep_id)
  
  K <- 35; adj <- matrix(0, K, K)
  for(i in 1:K) {
    r <- ceiling(i/7); c <- i - (r-1)*7
    if(r > 1) adj[i, i-7] <- 1; if(r < 5) adj[i, i+7] <- 1
    if(c > 1) adj[i, i-1] <- 1; if(c < 7) adj[i, i+1] <- 1
  }
  ordering <- 1:K
  
  W_true <- as.numeric(mvrnorm(1, rep(0, K), solve(create_dagar_precision(adj, ordering, 0.7)$Q)))
  cluster_id <- rep(1:K, each = 50); n <- length(cluster_id)
  X1 <- runif(n, 0, 1)
  X2 <- rbeta(n, 2, 6)
  X3 <- rbinom(n, 1, 0.8)
  X4 <- rbinom(n, 1, 0.7)
  X5 <- sample(1:10, n, replace = TRUE)
  X <- cbind(X1, X2, X3, X4, X5)
  
  b1_lin <- 0.5 * X1 - 0.3 * X2 + 0.2 * X3 - 0.1 * X4 + 0.05 * X2*(X5 - 5.5) + X1*X2
  Z <- rbinom(n, 1, prob = pnorm(b1_lin))
  
  logT_true <- sapply(1:n, function(i) b2_func_true(X[i,], Z[i], W_true[cluster_id[i]])) + rnorm(n, 0, sqrt(0.15))
  C_time <- rexp(n, rate = 0.1) 
  y_obs_log <- pmin(logT_true, log(C_time))
  status <- as.integer(logT_true <= log(C_time))
  
  censoring_pct <- mean(status == 0) * 100
  cat("Replicate", rep_id, "- Censoring Rate:", round(censoring_pct, 2), "%\n")
  
  # --- Stage 1 Strategy: Keep all variables (No longer omitting X2) ---
  ps_fit <- fit_ps_dagar(X, Z, cluster_id, adj, ordering = ordering)
  e_hat <- ps_fit$e_hat
  
  # --- Stage 2 Outcome Strategy: Preserve all baseline variables ---
  out_fit <- fit_outcome_bcf_dagar(exp(y_obs_log), status, X, Z, e_hat, cluster_id, adj, ordering = ordering)
  
  target_combos <- data.frame(county=c(4,4,10,10), t_star=c(5,15,5,15))
  rep_results <- data.frame()
  
  for(d in 1:4) {
    c_id <- target_combos$county[d]
    t_s <- target_combos$t_star[d]
    
    # Mathematical Counterfactual Ground Truth (incorporating X * W interactions)
    tr_s1 <- pnorm((log(t_s) - sapply(which(cluster_id==c_id), function(j) b2_func_true(X[j,], 1, W_true[c_id]))) / sqrt(0.15), lower.tail=F)
    tr_s0 <- pnorm((log(t_s) - sapply(which(cluster_id==c_id), function(j) b2_func_true(X[j,], 0, W_true[c_id]))) / sqrt(0.15), lower.tail=F)
    true_crate <- mean(tr_s1 - tr_s0)
    
    # Aggregate draws across target cluster targets
    c_idx_pool <- which(cluster_id == c_id)
    draws_matrix <- matrix(0, nrow = out_fit$n_save, ncol = length(c_idx_pool))
    for(j in seq_along(c_idx_pool)) {
      subj_id <- c_idx_pool[j]
      draws_matrix[, j] <- estimate_cerm(t_s, X[subj_id,], c_id, e_hat[subj_id], out_fit)
    }
    acerm_draws <- rowMeans(draws_matrix)
    
    est_mean <- mean(acerm_draws)
    ci <- quantile(acerm_draws, c(0.025, 0.975))
    
    rep_results <- rbind(rep_results, data.frame(County = c_id, Time = t_s, MAE = abs(est_mean - true_crate), Covered = (true_crate >= ci[1] & true_crate <= ci[2]), CI_Width = ci[2] - ci[1]))
  }
  return(rep_results)
}

# ==============================================================================
# 7. RUN & AGGREGATE EXECUTION
# ==============================================================================

all_reps <- lapply(1:150, function(r) {
  cat("Starting Replicate:", r, "\n")
  run_replicate(r)
})

final_df <- do.call(rbind, all_reps)
summary_stats <- aggregate(cbind(MAE, Covered, CI_Width) ~ County + Time, data = final_df, FUN = mean)
print(summary_stats)
