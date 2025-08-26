#===============================================================================
#
#  PROGRAM: TUNE_01.R
#
#  AUTHORS: Stephen Salerno, Zhilin Zhang
#
#  PURPOSE: Code for hyperparameter tuning for Method 01: Lee et al. (2015)
#
#===============================================================================

#=== SETUP =====================================================================

mean_grid <- seq(0.1, 1, by = 0.1)

sd_grid <- seq(0.1, 1, by = 0.1)

results <- expand.grid(mean = mean_grid, sd = sd_grid)

results$bbs <- NA

myModel <- c("Markov", "Weibull")

form <- Formula(Yi1 + Di1 | Yi2 + Di2 ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 +

  X9 + X10 + X11 + X12 | X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 +

  X11 + X12 | X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + X11 + X12)

#=== TUNE OVER GRID OF PARAMETERS ==============================================

for (i in seq_len(nrow(results))) {
  
  mean_theta <- results$mean[i]
  
  sd_theta <- results$sd[i]

  theta.ab <- c(mean_theta, sd_theta)
  
  WB.ab1 <- c(0.5, 0.01)
  WB.ab2 <- c(0.5, 0.01)
  WB.ab3 <- c(0.5, 0.01)

  WB.cd1 <- c(0.5, 0.05)
  WB.cd2 <- c(0.5, 0.05)
  WB.cd3 <- c(0.5, 0.05)

  hyperParams <- list(theta = theta.ab,

    WB = list(WB.ab1 = WB.ab1, WB.ab2 = WB.ab2, WB.ab3 = WB.ab3,
              WB.cd1 = WB.cd1, WB.cd2 = WB.cd2, WB.cd3 = WB.cd3))

  numReps    <- 1000
  thin       <- 10
  burninPerc <- 0.5

  nGam_save <- 0
  storeV    <- rep(F, 3)

  mhProp_theta_var  <- 0.05
  mhProp_alphag_var <- c(0.01, 0.01, 0.01)

  mcmc.WB <- list(

    run = list(numReps = numReps, thin = thin, burninPerc = burninPerc),

    storage = list(nGam_save = nGam_save, storeV = storeV),

    tuning = list(mhProp_theta_var  = mhProp_theta_var,
                  mhProp_alphag_var = mhProp_alphag_var))

  dat_sim <- data.frame(dat[,, sim])

  colnames(dat_sim) <- c("Yi1", "Di1", "Yi2", "Di2", paste0("X", 1:12),

    "h1", "h2", "h3", "hc", "gamma", "Ti1", "Ti2", "Ti3", "ctime",

    "cind", "setting", "n", "theta", "risk", "cens")

  startValues <- initiate.startValues_HReg(

    form, dat_sim, model = myModel, nChain = 2)

  fit_01 <- BayesID_HReg(form, dat_sim, id = NULL, model = myModel,

    hyperParams, startValues, mcmc.WB)

  beta1 <- log(summary(fit_01)$coef[,1])
  beta2 <- log(summary(fit_01)$coef[,4])
  beta3 <- log(summary(fit_01)$coef[,7])

  h1 <- as.matrix(dat_sim[, paste0("X", 1:12)]) %*% beta1
  h2 <- as.matrix(dat_sim[, paste0("X", 1:12)]) %*% beta2
  h3 <- as.matrix(dat_sim[, paste0("X", 1:12)]) %*% beta3

  theta <- rep(summary(fit_01)$theta[1], n_sim)

  k1 <- exp(rep(summary(fit_01)$h0[1,1], n_sim))
  a1 <- exp(rep(summary(fit_01)$h0[2,1], n_sim))

  k2 <- exp(rep(summary(fit_01)$h0[1,4], n_sim))
  a2 <- exp(rep(summary(fit_01)$h0[2,4], n_sim))

  k3 <- exp(rep(summary(fit_01)$h0[1,7], n_sim))
  a3 <- exp(rep(summary(fit_01)$h0[2,7], n_sim))

  preds <- sapply(t_pred, function(t) {
    (1 +
       theta * k1 * t^a1 * exp(h1) +
       theta * k2 * t^a2 * exp(h2)
    )^(-1/theta)
  })

  res <- bbs1(test_dat, preds, t_pred)

  if (!inherits(sim_result, "try-error")) {
    results$bbs[i] <- res$ibs
  }
}

#=== DISPLAY BEST HYPERPARAMETER SETTINGS ======================================

best_idx <- which.min(results$bcindex)

best_params <- results[best_idx, ]

print(best_params)

#=== END =======================================================================
