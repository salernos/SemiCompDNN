#===============================================================================
#
#  PROGRAM: METHOD_01.R
#
#  AUTHORS: Stephen Salerno, Zhilin Zhang
#
#  PURPOSE: Generate simulation results for Method 01: Lee et al. (2015)
#
#  INPUT:   simdat{X}.RData
#
#  OUTPUT:  Sim_01_{X}.RData
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

p_load(SemiCompRisks, progressr, doFuture, abind, here, update = F)

#--- SIMULATION SETTINGS -------------------------------------------------------

nsims <- 500

n <- c(1000, 5000)

theta <- c(0.5, 2.0)

risk <- c("Linear", "Nonlinear")

cens <- c(0.25, 0.50)

settings <- expand.grid(n, theta, risk, cens)

#--- ANALYSIS SETTINGS ---------------------------------------------------------

myModel <- c("Markov", "Weibull")

form <- Formula(Yi1 + Di1 | Yi2 + Di2 ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 +

  X9 + X10 + X11 + X12 | X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 +

  X11 + X12 | X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + X11 + X12)

#--- HYPERPARAMETERS -----------------------------------------------------------

theta.ab <- c(0.7, 0.7)

WB.ab1 <- c(0.5, 0.01)
WB.ab2 <- c(0.5, 0.01)
WB.ab3 <- c(0.5, 0.01)

WB.cd1 <- c(0.5, 0.05)
WB.cd2 <- c(0.5, 0.05)
WB.cd3 <- c(0.5, 0.05)

hyperParams <- list(theta = theta.ab,

  WB = list(WB.ab1 = WB.ab1, WB.ab2 = WB.ab2, WB.ab3 = WB.ab3,
            WB.cd1 = WB.cd1, WB.cd2 = WB.cd2, WB.cd3 = WB.cd3))

#--- MCMC SETTINGS -------------------------------------------------------------

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

#=== RUN SIMULATION STUDY ======================================================

plan(multisession)

options(future.globals.maxSize = 600*1024^2)

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

        myPath <- paste0("Results/Method_01/Sim_", i, "_", sim, "/")

        startValues <- initiate.startValues_HReg(

          form, dat_sim, model = myModel, nChain = 2)

        fit_01 <- BayesID_HReg(form, dat_sim, id = NULL, model = myModel,

          hyperParams, startValues, mcmc.WB, path = myPath)

        beta1 <- log(summary(fit_01)$coef[,1])
        beta2 <- log(summary(fit_01)$coef[,4])
        beta3 <- log(summary(fit_01)$coef[,7])

        h1 <- as.matrix(dat_sim[, paste0("X", 1:12)]) %*% beta1
        h2 <- as.matrix(dat_sim[, paste0("X", 1:12)]) %*% beta2
        h3 <- as.matrix(dat_sim[, paste0("X", 1:12)]) %*% beta3

        theta <- rep(summary(fit_01)$theta[1], n_sim)

        kappa1 <- exp(rep(summary(fit_01)$h0[1,1], n_sim))
        alpha1 <- exp(rep(summary(fit_01)$h0[2,1], n_sim))

        kappa2 <- exp(rep(summary(fit_01)$h0[1,4], n_sim))
        alpha2 <- exp(rep(summary(fit_01)$h0[2,4], n_sim))

        kappa3 <- exp(rep(summary(fit_01)$h0[1,7], n_sim))
        alpha3 <- exp(rep(summary(fit_01)$h0[2,7], n_sim))

        sim_result <- cbind(

          theta, kappa1, alpha1, kappa2, alpha2, kappa3, alpha3, h1, h2, h3)

        return(sim_result)
    }
  })

  sim_array <- array(unlist(sim_results), dim = c(n_sim, 10, nsims))

  save(sim_array, file = here("Results", paste0("Sim_01_", i, ".RData")))
}

toc <- proc.time(); toc - tic

#=== END =======================================================================
