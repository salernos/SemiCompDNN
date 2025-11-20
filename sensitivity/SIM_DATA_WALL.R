#===============================================================================
#
#  PROGRAM: SIM_DATA_WALL.R
#
#  AUTHOR:  Stephen Salerno, Zhilin Zhang
#
#  PURPOSE: Generate data for wall time sensitivity analysis
#
#  INPUT:   None
#
#  OUTPUT:  sim_linear_np_setting_{i}.RData
#
#  NOTES:   The simulated data are stored in twelve files, with each file
#           corresponding to one simulation setting.
#
#===============================================================================

#=== INITIALIZATION ============================================================

#--- NECESSARY PACKAGES --------------------------------------------------------

library(pacman)

p_load(MASS, tmvtnorm, here, tidyverse,TruncatedNormal, update = FALSE)

#=== HELPER FUNCTIONS ==========================================================

#--- DATA GENERATION -----------------------------------------------------------

simdat <- function(n, p,
                   theta = 0.5,
                   covariates = "normal",
                   risk = "linear",
                   cens = 0.25,
                   cens_admin = 5,
                   seed = NULL) {
  
  if (missing(n)) stop('Error: argument "n" is missing')
  if (missing(p)) stop('Error: argument "p" is missing')
  if (!is.null(seed)) set.seed(seed)
  
  if (theta == 0) {
    
    gamma <- rep(1, n)
  
  } else {
    
    gamma <- rgamma(n, shape = 1/theta, rate = 1/theta)
  }

  if (covariates == "uniform") {
    
    X <- matrix(runif(n * p, -1, 1), nrow = n, ncol = p)
  
  } else { 
    
    mu_X <- rep(0, p)
    
    Sigma_X <- matrix(0.2, nrow = p, ncol = p); diag(Sigma_X) <- 1
    
    X <- rtmvnorm(n, mean = mu_X, sigma = Sigma_X, 
                  lower = rep(-2, p), upper = rep(2, p))
  }
  
  colnames(X) <- paste0("X", seq_len(p))
  
  k_neg <- max(1, round(0.4 * p))
  
  beta  <- c(rep(-0.5, k_neg), rep(0.5, p - k_neg))
  
  if (risk == "none") {
    
    h1 <- h2 <- h3 <- 0
  
  } else if (risk == "linear") {
    
    h1 <- as.vector(X %*% beta)
    h2 <- as.vector(X %*% beta)
    h3 <- as.vector(X %*% beta)
  
  } else if (risk == "nonlinear") {
    
    stop("This script is for the linear case only. Set risk='linear'.")
  }
  
  if (cens == 0.25) {
    
    muc <- if (theta == 0.5) 0.5 else if (theta == 2) 0.05 else 0.5
    
  } else if (cens == 0.50) {
    
    muc <- if (theta == 0.5) 2 else if (theta == 2) 1 else 2
    
  } else {
    
    muc <- 0.5
  }
  
  Ti1 <- -log(runif(n)) / (gamma * exp(-h1)) / 2
  Ti2 <- -log(runif(n)) / (gamma * exp(-h2)) / 3
  Ti3 <- -log(runif(n)) / (gamma * exp(-h3)) / 2
  
  ind <- Ti1 < Ti2
  Ti2[ind] <- Ti1[ind] + Ti3[ind]
  
  s <- min(6, p)
  
  tail_idx <- (p - s + 1):p
  
  alt_coef <- rep(c(-0.5, 0.5), length.out = s)
  
  if (risk == "none") {
    
    hc <- 0 
  
  } else {
    
    hc <- as.vector(X[, tail_idx, drop = FALSE] %*% alt_coef)
  }
  
  Ci  <- pmin(-log(runif(n)) / (muc * exp(-hc)), cens_admin)
  Yi2 <- pmin(Ti2, Ci)
  Yi1 <- pmin(Ti1, Yi2)
  
  Di2 <- as.integer(Ti2 <= Ci)
  Di1 <- as.integer(Ti1 <= Ti2)
  
  ctime <- ifelse(Di1 == 0 & Di2 == 0, Yi1, pmax(Yi1, Yi2))
  cind  <- as.integer(Di1 == 0 & Di2 == 0)
  
  dat <- tibble(Yi1 = Yi1, Di1 = Di1, Yi2 = Yi2, Di2 = Di2) |>
    
    bind_cols(as_tibble(X)) |>
    
    mutate(h1 = h1, h2 = h2, h3 = h3, hc = hc, gamma = gamma,
      Ti1 = Ti1, Ti2 = Ti2, Ti3 = Ti3, ctime = ctime, cind = cind)
  
  return(dat)
}

#=== GENERATE DATA =============================================================

#--- DATA GENERATING MECHANISMS ------------------------------------------------

np_grid <- tibble::tibble(
  n = c(100, 100, 100, 100, 1000, 1000, 1000, 1000, 10000, 10000, 10000, 10000),
  p = c(1,   10,  100, 1000,1,    10,   100,  1000, 1,     10,    100,   1000)
)

nsims <- 1
theta_vals <- 0.5    
cens_vals  <- 0.25

settings <- tidyr::crossing(np_grid, theta = theta_vals, cens = cens_vals) |>
  mutate(risk = "linear")

#--- LOOP OVER SETTINGS --------------------------------------------------------

dir.create(here("Data"), showWarnings = FALSE)

for (i in seq_len(nrow(settings))) {
  
  n_i <- settings$n[i]; p_i <- settings$p[i]
  
  th  <- settings$theta[i]; cs <- settings$cens[i]
  
  cat("Setting", i, "— n =", n_i, "p =", p_i, "theta =", th, "cens =", cs, "\n")
  
  K <- p_i + 20
  
  dat_arr <- array(NA_real_, dim = c(n_i, K, nsims))
  
  for (sim in seq_len(nsims)) {
    
    if (sim %% 50 == 0) cat("  sim", sim, "\n")
    
    dat_sim <- simdat(n = n_i, p = p_i, theta = th, 
                      covariates = "uniform", risk = "linear",
                      cens = cs, seed = sim) |>
      
      mutate(setting = i, n = n_i, p = p_i, theta = th, risk = 0L, cens = cs)
    
    ord <- c(c("Yi1","Di1","Yi2","Di2"), paste0("X", seq_len(p_i)),
             c("h1","h2","h3","hc","gamma","Ti1","Ti2","Ti3","ctime","cind",
               "setting","n","p","theta","risk","cens"))
    
    dat_mat <- as.matrix(dat_sim[, ord])
    
    dat_arr[,,sim] <- dat_mat
  }
  
  save(dat_arr, 
       file = here("Data", sprintf("sim_linear_np_setting_%02d.RData", i)))
}

#=== END =======================================================================