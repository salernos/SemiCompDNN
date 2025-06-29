#===============================================================================
#
#  PROGRAM: METHOD_04.R
#
#  AUTHOR:  Stephen Salerno
#
#  PURPOSE: Generate simulation results for Method 04: Kats et al. (2022)
#
#  INPUT:   simdat{X}.RData
#
#  OUTPUT:  Sim_04_{X}.RData
#
#  NOTES:
#
#  The simulation results are stored in sixteen .RData files, with each file 
#  corresponding to one simulation setting. Each file is an n x 10 x nsims 
#  array, storing the estimated frailty variance, baseline hazard parameters, 
#  and log risk functions for each simulated dataset.
#
#===============================================================================

#=== INITIALIZATION ============================================================

#--- NECESSARY PACKAGES --------------------------------------------------------

library(pacman)

p_load(progressr, doFuture, abind, here, update = F)

source(here("semicompAFT-main", "semicompAFT.R"))

#--- SIMULATION SETTINGS -------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

k <- as.integer(args[1])  #- Setting Index

i <- 10

nsims <- 500

n <- c(1000, 5000)

theta <- c(0.5, 2.0)

risk <- c("Linear", "Nonlinear")

cens <- c(0.25, 0.50)

settings <- expand.grid(n = n, theta = theta, risk = risk, cens = cens)

cat("\nSetting", i, "of", nrow(settings), "---\n\n")

#=== RUN SIMULATION STUDY ======================================================

load(here("Data", paste0("simdat", i, ".RData")))

n_sim <- dim(dat)[1]

sim_array <- array(NA, dim = c(n_sim, 7, nsims))

for (sim in 1:nsims) {
  
  print(k)
  
  tryCatch({
    
    if (sim %% 1 == 0) cat("  Sim", sim, "of", nsims, "\n")
    
    dat_sim <- data.frame(dat[,, sim+k*5])
    
    colnames(dat_sim) <- c("Yi1", "Di1", "Yi2", "Di2", paste0("X", 1:12),
                           
      "h1", "h2", "h3", "hc", "gamma", "Ti1", "Ti2", "Ti3", "ctime",
                           
      "cind", "setting", "n", "theta", "risk", "cens")
    
    Z12 <- Z13 <- Z23 <- as.matrix(dat_sim[, paste0("X", 1:12)])
    
    V <- dat_sim[, "Yi1"]
    W <- dat_sim[, "Yi2"]
    
    delta1 <- dat_sim$Yi1 == dat_sim$Ti1
    delta2 <- dat_sim$Yi1 == dat_sim$Ti2
    delta3 <- dat_sim$Yi2 == dat_sim$Ti2
    
    fit_04 <- estimation_with_frailty(X01 = Z12, X02 = Z13, X12 = Z23, V = V,
                                      
      W = W, delta1 = delta1, delta2 = delta2, delta3 = delta3, B = 0,
                                      
      print = F)
    
    beta1 <- fit_04$ests[ 2:13, 1]
    beta2 <- fit_04$ests[14:25, 1]
    beta3 <- fit_04$ests[26:37, 1]
    
    h1 <- as.matrix(dat_sim[, paste0("X", 1:12)]) %*% beta1
    h2 <- as.matrix(dat_sim[, paste0("X", 1:12)]) %*% beta2
    h3 <- as.matrix(dat_sim[, paste0("X", 1:12)]) %*% beta3
    
    theta <- rep(fit_04$ests[1, 1], n_sim)
    
    lam01 <- fit_04$H0s[[1]]
    lam02 <- fit_04$H0s[[2]]
    lam03 <- fit_04$H0s[[3]]
    
    sim_array[,, sim] <- cbind(theta, lam01, lam02, lam03, h1, h2, h3)
    
  }, error = function(e) {
    
    cat("Error in simulation", i, ": ", e$message, "\n")
    
    sim_array[,, sim] <- matrix(nrow = n_sim, ncol = 7)
  })
}

save(sim_array, file = here("Results", paste0("Sim_04_", i,"_", k, ".RData")))

#=== END =======================================================================
