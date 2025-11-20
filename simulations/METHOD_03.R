#===============================================================================
#
#  PROGRAM: METHOD_03.R
#
#  AUTHORS: Stephen Salerno, Zhilin Zhang
#
#  PURPOSE: Generate simulation results for Method 03: Gorfine et al. (2020)
#
#  INPUT:   simdat{X}.RData
#
#  OUTPUT:  Sim_03_{X}.RData
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

source(here("frailty-LTRC-master", "Estimation Procedure.R"))

#--- SIMULATION SETTINGS -------------------------------------------------------

nsims <- 500

n <- c(1000, 5000)

theta <- c(0.5, 2.0)

risk <- c("Linear", "Nonlinear")

cens <- c(0.25, 0.50)

settings <- expand.grid(n, theta, risk, cens)

#=== RUN SIMULATION STUDY ======================================================

plan(multisession)

tic <- proc.time()

for (i in 1:nrow(settings)) {

  cat("\nSetting", i, "of", nrow(settings), "---\n\n")

  #--- LOAD DATA ---------------------------------------------------------------

  load(here("Data", paste0("simdat", i, ".RData")))

  n_sim <- dim(dat)[1]

  #--- RUN OVER NSIMS DATASETS -------------------------------------------------

  idx <- 1:nsims

  with_progress({

    p <- progressor(along = idx)

    sim_results <- foreach(sim = idx,

      .options.future = list(seed = T)) %dofuture% {

        p()

        dat_sim <- data.frame(dat[,, sim])

        colnames(dat_sim) <- c("Yi1", "Di1", "Yi2", "Di2", paste0("X", 1:12),

          "h1", "h2", "h3", "hc", "gamma", "Ti1", "Ti2", "Ti3", "ctime",

          "cind", "setting", "n", "theta", "risk", "cens")

        Z12 <- Z13 <- Z23 <- as.matrix(dat_sim[, paste0("X", 1:12)])

        V <- dat_sim[, "Yi1"]

        W <- dat_sim[, "Yi2"]

        delta1 <- dat_sim$Yi1 == dat_sim$Ti1

        delta2 <- dat_sim$Yi1 == dat_sim$Ti2

        delta3 <- dat_sim$Yi2 == dat_sim$Ti2

        fit_03 <- estimate_all(z12 = Z12, z13 = Z13, z23 = Z23,

          obs_times1 = V, obs_times2 = W,

          delta1 = delta1, delta2 = delta2, delta3 = delta3,

          save_address = "Results/Method_03/",

          name = paste0("simulation number_", seed))

        beta1 <- log(summary(fit_02)$coef[,1])
        beta2 <- log(summary(fit_02)$coef[,4])
        beta3 <- log(summary(fit_02)$coef[,7])

        h1 <- as.matrix(dat_sim[, paste0("X", 1:12)]) %*% beta1
        h2 <- as.matrix(dat_sim[, paste0("X", 1:12)]) %*% beta2
        h3 <- as.matrix(dat_sim[, paste0("X", 1:12)]) %*% beta3

        theta <- rep(summary(fit_02)$theta[1], n_sim)

        mu1  <- rep(summary(fit_02)$h0[1,1], n_sim)
        sig1 <- rep(summary(fit_02)$h0[2,1], n_sim)

        mu2  <- rep(summary(fit_02)$h0[1,4], n_sim)
        sig2 <- rep(summary(fit_02)$h0[2,4], n_sim)

        mu3  <- rep(summary(fit_02)$h0[1,7], n_sim)
        sig3 <- rep(summary(fit_02)$h0[2,7], n_sim)

        sim_result <- cbind(theta, mu1, sig1, mu2, sig2, mu3, sig3, h1, h2, h3)

        return(sim_result)
    }
  })

  sim_array <- array(unlist(sim_results), dim = c(n_sim, 10, nsims))

  save(sim_array, file = here("Results", paste0("Sim_02_", i, ".RData")))
}

toc <- proc.time(); toc - tic

#=== END =======================================================================
