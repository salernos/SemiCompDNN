#===============================================================================
#
#  PROGRAM: SIM_TREE.R
#
#  AUTHORS: Stephen Salerno, Zhilin Zhang
#
#  PURPOSE: Additional sensitivity analysis results from tree-based method.
#
#  INPUT:   simdat{X}.RData
#
#  NOTES:   
#
#  The tree-based method for competing risks prediction is not directly 
#  comparable in our semi-competing risks setting in terms of parameter 
#  estimation. However we apply this this method in simulation to compare
#  in terms of our bivariate Brier score and bivariate concordance index.
#
#===============================================================================

#=== INITIALIZATION ============================================================

#--- NECESSARY PACKAGES --------------------------------------------------------

library(pacman)

p_load(randomForestSRC, progressr, doFuture, abind, here, survival, prodlim,
       
       riskRegression, pec, rtmvnorm, update = FALSE)

source("./Simulation/METRICS1.R")

args <- commandArgs(trailingOnly = TRUE)

k <- as.integer(args[1])  

i <- 8

nsims <- 500

load(here("Data", paste0("simdat", i, ".RData")))

n_sim <- dim(dat)[1]

cinds <- c()

ibs_true <- c()

cr_surv_from_trees <- function(pr1, pr2, t_pred) {
  
  t1 <- pr1$time.interest; H1 <- pr1$chf
  t2 <- pr2$time.interest; H2 <- pr2$chf
  
  t_grid <- sort(unique(c(t1, t2)))
  
  step_interp_mat <- function(times_from, mat_from, times_to) {
    t(apply(mat_from, 1, function(rowH){
      approx(x = times_from, y = rowH, xout = times_to,
             method = "constant", f = 0, rule = 2)$y
    }))
  }
  
  H1g <- step_interp_mat(t1, H1, t_grid)
  H2g <- step_interp_mat(t2, H2, t_grid)
  
  Htot <- H1g + H2g
  Sg   <- exp(-Htot)
  
  incr <- function(Hg) {
    cbind(Hg[,1,drop=FALSE], Hg[, -1, drop=FALSE] - Hg[, -ncol(Hg), drop=FALSE])
  }
  
  dH1 <- incr(H1g)
  dH2 <- incr(H2g)
  
  S_left <- cbind(1, Sg[, -ncol(Sg), drop=FALSE])
  
  CIF1g <- t(apply(S_left * dH1, 1, cumsum))
  CIF2g <- t(apply(S_left * dH2, 1, cumsum))
  
  step_interp_vecs <- function(times_from, mat_from, times_to) {
    t(apply(mat_from, 1, function(rowv){
      approx(x = times_from, y = rowv, xout = times_to,
             method = "constant", f = 0, rule = 2)$y
    }))
  }
  
  S_pred    <- step_interp_vecs(t_grid, Sg,    t_pred)
  CIF1_pred <- step_interp_vecs(t_grid, CIF1g, t_pred)
  CIF2_pred <- step_interp_vecs(t_grid, CIF2g, t_pred)
  
  colnames(S_pred)    <- paste0("S_t", t_pred)
  colnames(CIF1_pred) <- paste0("CIF1_t", t_pred)
  colnames(CIF2_pred) <- paste0("CIF2_t", t_pred)
  
  list(times = t_pred, S = S_pred, CIF1 = CIF1_pred, CIF2 = CIF2_pred)
}

for (sim in 1:nsims) {
  
  if (sim %% 1 == 0) cat("  Sim", sim, "of", nsims, "\n")
  
  dat_sim <- data.frame(dat[,, sim])
  
  colnames(dat_sim) <- c(
    "Yi1","Di1","Yi2","Di2",
    paste0("X", 1:12),
    "h1","h2","h3","hc","gamma",
    "Ti1","Ti2","Ti3","ctime",
    "cind","setting","n","theta","risk","cens"
  )
  
  # 80/20 train/test split 
  
  set.seed(1000 + sim)
  
  n <- nrow(dat_sim)
  
  train_idx <- sample.int(n, floor(0.8 * n))
  
  train <- dat_sim[train_idx, , drop = FALSE]
  test  <- dat_sim[-train_idx, , drop = FALSE]
  
  # Fit model on training split
  
  rsf_h1 <- rfsrc(
    Surv(Yi1, Di1) ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + X11 + X12,
    data = train, ntree = 500, splitrule = "logrank"
  )
  
  rsf_h2 <- rfsrc(
    Surv(Yi2, Di2) ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + X11 + X12,
    data = train, ntree = 500, splitrule = "logrank"
  )
  
  # Predict on testing split
  
  pr1 <- predict(rsf_h1, newdata = test)
  pr2 <- predict(rsf_h2, newdata = test)
  
  t_pred <- seq(0, max(c(test$Yi1, test$Yi2), na.rm = TRUE), length.out = 100)
  
  S_pred <- cr_surv_from_trees(pr1, pr2, t_pred)$S
  
  # Metrics on testing split
  
  res  <- bbs1(test, S_pred, t_pred)
  IBBS <- res$ibs
  
  test$event_all <- ifelse(test$Di1 == 1 | test$Di2 == 1, 1, 0)
  test$time_all  <- pmin(test$Yi1, test$Yi2, na.rm = TRUE)
  
  risk_score <- rowMeans(1 - S_pred) 
  cind <- survConcordance(Surv(time_all, event_all) ~ risk_score, data = test)
  
  cinds    <- c(cinds,    cind[[1]])
  ibs_true <- c(ibs_true, IBBS)
}

cat("Mean C-index:", mean(cinds), "\n")
cat("SD   C-index:", sd(cinds), "\n")
cat("Mean IBBS:   ", mean(ibs_true), "\n")
cat("SD   IBBS:   ", sd(ibs_true), "\n")

#=== END =======================================================================