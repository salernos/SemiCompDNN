#===============================================================================
#
#  PROGRAM: METHOD_02.R
#
#  AUTHORS: Stephen Salerno, Zhilin Zhang
#
#  PURPOSE: Generate simulation results for Method 02: Lee et al. (2017)
#
#  INPUT:   simdat{X}.RData
#
#  OUTPUT:  Sim_02_{X}.RData
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

myModel <- "LN"

form <- Formula(LT | y1L + y1U | y2L + y2U ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 +

  X8 + X9 + X10 + X11 + X12 | X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 +

  X11 + X12 | X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + X11 + X12)

#--- HYPERPARAMETERS -----------------------------------------------------------

theta.ab <- c(0.5, 0.05)

LN.ab1 <- c(0.3, 0.3)
LN.ab2 <- c(0.3, 0.3)
LN.ab3 <- c(0.3, 0.3)

hyperParams <- list(theta=theta.ab,

  LN=list(LN.ab1=LN.ab1, LN.ab2=LN.ab2, LN.ab3=LN.ab3))

#--- MCMC SETTINGS -------------------------------------------------------------

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

        dat_sim$y1L <- dat_sim$y1U <- dat_sim[,1]

        dat_sim$y1U[which(dat_sim[,2] == 0)] <- Inf

        dat_sim$y2L <- dat_sim$y2U <- dat_sim[,3]

        dat_sim$y2U[which(dat_sim[,4] == 0)] <- Inf

        dat_sim$LT <- rep(0, dim(dat_sim)[1])

        myPath <- paste0("Results/Method_02/Sim_", i, "_", sim, "/")

        startValues <- initiate.startValues_AFT(

          form, dat_sim, model = myModel, nChain = 2)

        fit_02 <- BayesID_AFT(form, dat_sim, model = myModel, hyperParams,

          startValues, mcmcParams, path = myPath)

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
