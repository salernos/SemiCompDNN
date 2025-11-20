#===============================================================================
#
#  PROGRAM: TUNE_02.R
#
#  AUTHORS: Stephen Salerno, Zhilin Zhang
#
#  PURPOSE: Code for hyperparameter tuning for Method 02: Lee et al. (2017)
#
#===============================================================================

#=== SETUP =====================================================================

mean_grid <- seq(0.1, 1, by = 0.1)

sd_grid <- seq(0.1, 1, by = 0.1)

results <- expand.grid(mean = mean_grid, sd = sd_grid)

results$bbs <- NA

myModel <- "LN"

form <- Formula(LT | y1L + y1U | y2L + y2U ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 +
                  
  X8 + X9 + X10 + X11 + X12 | X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 +
                  
  X11 + X12 | X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + X11 + X12)

#=== TUNE OVER GRID OF PARAMETERS ==============================================

for (i in seq_len(nrow(results))) {
  
  mean_theta <- results$mean[i]
  
  sd_theta <- results$sd[i]

  theta.ab <- c(mean_theta, sd_theta)

  LN.ab1 <- c(0.3, 0.3)
  LN.ab2 <- c(0.3, 0.3)
  LN.ab3 <- c(0.3, 0.3)

  hyperParams <- list(theta = theta.ab,

    LN = list(LN.ab1 = LN.ab1, LN.ab2 = LN.ab2, LN.ab3 = LN.ab3))

  #--- MCMC SETTINGS -----------------------------------------------------------

  numReps    <- 300
  thin       <- 3
  burninPerc <- 0.5

  nGam_save   <- 10
  nY1_save    <- 10
  nY2_save    <- 10
  nY1.NA_save <- 10

  betag.prop.var <- c(0.01,0.01,0.01)
  mug.prop.var	 <- c(0.1,0.1,0.1)
  zetag.prop.var <- c(0.1,0.1,0.1)
  gamma.prop.var <- 0.01

  mcmcParams <- list(

    run = list(numReps = numReps, thin = thin, burninPerc = burninPerc),

    storage = list(nGam_save = nGam_save, nY1_save = nY1_save,

      nY2_save = nY2_save, nY1.NA_save = nY1.NA_save),

    tuning = list(betag.prop.var = betag.prop.var, mug.prop.var = mug.prop.var,

      zetag.prop.var = zetag.prop.var, gamma.prop.var = gamma.prop.var))

  dat_sim <- data.frame(dat[,, sim])

  colnames(dat_sim) <- c("Yi1", "Di1", "Yi2", "Di2", paste0("X", 1:12),

    "h1", "h2", "h3", "hc", "gamma", "Ti1", "Ti2", "Ti3", "ctime",

    "cind", "setting", "n", "theta", "risk", "cens")
  
  dat_sim$y1L <- dat_sim$y1U <- dat_sim[,1]

  dat_sim$y1U[which(dat_sim[,2] == 0)] <- Inf

  dat_sim$y2L <- dat_sim$y2U <- dat_sim[,3]

  dat_sim$y2U[which(dat_sim[,4] == 0)] <- Inf

  dat_sim$LT <- rep(0, dim(dat_sim)[1])

  startValues <- initiate.startValues_AFT(

    form, dat_sim, model = myModel, nChain = 2)

  fit_02 <- BayesID_AFT(form, dat_sim, model = myModel, hyperParams,

    startValues, mcmcParams)

  beta1 <- log(summary(fit_02)$coef[,1])
  beta2 <- log(summary(fit_02)$coef[,4])
  beta3 <- log(summary(fit_02)$coef[,7])

  h1 <- as.matrix(dat_sim[, paste0("X", 1:12)]) %*% beta1
  h2 <- as.matrix(dat_sim[, paste0("X", 1:12)]) %*% beta2
  h3 <- as.matrix(dat_sim[, paste0("X", 1:12)]) %*% beta3

  theta <- rep(summary(fit_02)$theta[1], n_sim)

  k1  <- rep(summary(fit_02)$h0[1,1], n_sim)
  a1 <- rep(summary(fit_02)$h0[2,1], n_sim)

  k2  <- rep(summary(fit_02)$h0[1,4], n_sim)
  a2 <- rep(summary(fit_02)$h0[2,4], n_sim)

  k3  <- rep(summary(fit_02)$h0[1,7], n_sim)
  a3 <- rep(summary(fit_02)$h0[2,7], n_sim)

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
