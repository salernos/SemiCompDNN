#===============================================================================
#
#  PROGRAM: SIM_SPLINE.R
#
#  AUTHORS: Stephen Salerno, Zhilin Zhang
#
#  PURPOSE: Generate simulation results for the spline-based method.
#
#  INPUT:   simdat{X}.RData
#
#  OUTPUT:  Sim_06_{X}.RData
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

# The file functions.txt is available at:
# https://bmcmedresmethodol.biomedcentral.com/articles/10.1186/s12874-020-01203-8#Sec15
# We saved this as functions.R for convenience

source(here("functions.R")) 

#--- SIMULATION SETTINGS -------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

k <- as.integer(args[1])

i = 3

nsims = 500

#=== RUN SIMULATION STUDY ======================================================

load(here("Data", paste0("simdat", i, ".RData")))

n_sim <- dim(dat)[1]

sim_array <- array(NA, dim = c(n_sim, 13, nsims))

for (sim in k:k) {
  
  if (sim %% 1 == 0) cat("  Sim", sim, "of", nsims, "\n")
  
  dat_sim <- data.frame(dat[,, sim])
  
  colnames(dat_sim) <- c("Yi1", "Di1", "Yi2", "Di2", paste0("X", 1:12),
                         
    "h1", "h2", "h3", "hc", "gamma", "Ti1", "Ti2", "Ti3", "ctime",
                         
    "cind", "setting", "n", "theta", "risk", "cens")
  
  Y <- data.frame(
    y1     = dat_sim$Yi1,
    delta1 = dat_sim$Di1,
    y2     = dat_sim$Yi2,
    delta2 = dat_sim$Di2,
    L      = rep(0,n_sim)
  )
  
  Y$L <- fix_L(Y$L, Y$y1)
  
  lin.pred <- list(
    as.formula(paste("~", paste0("X", 1:12, collapse = " + "))), 
    as.formula(paste("~", paste0("X", 1:12, collapse = " + "))),  
    as.formula(paste("~", paste0("X", 1:12, collapse = " + ")))  
  )
  
  
  frailty <- TRUE 
  
  y1 <- Y$y1; delta1 <- Y$delta1
  y2 <- Y$y2; delta2 <- Y$delta2
  
  n.internalKnots.1 <- 1; Bspline.degree.1 <- 1
  n.internalKnots.2 <- 1; Bspline.degree.2 <- 1
  n.internalKnots.3 <- 1; Bspline.degree.3 <- 1
  
  num.Bspline.params.1 <- n.internalKnots.1 + (Bspline.degree.1 + 1)
  num.Bspline.params.2 <- n.internalKnots.2 + (Bspline.degree.2 + 1)
  num.Bspline.params.3 <- n.internalKnots.3 + (Bspline.degree.3 + 1)
  
  bdy.knots.b.1      <- c(0, max(y1))
  bdy.knots.b.2      <- c(0, max(y2))
  bdy.knots.b.3.y2my1<- c(0, max(y2 - y1))
  
  knot.loc.1 <- stats::quantile(y1[delta1 == 1], ppoints(n.internalKnots.1))
  b.1.event  <- bSpline(y1[delta1 == 1], knots = knot.loc.1, 
                        degree = Bspline.degree.1,
                        intercept = TRUE, Boundary.knots = bdy.knots.b.1)
  b.1        <- predict(b.1.event, y1)
  
  knot.loc.2 <- stats::quantile(y2[delta2 == 1], ppoints(n.internalKnots.2))
  b.2.event  <- bSpline(y2[delta2 == 1], knots = knot.loc.2, 
                        degree = Bspline.degree.2,
                        intercept = TRUE, Boundary.knots = bdy.knots.b.2)
  b.2        <- predict(b.2.event, y2)
  
  knot.loc.3 <- stats::quantile((y2 - y1)[delta1 == 1], 
                                ppoints(n.internalKnots.3))
  b.3.event  <- bSpline((y2 - y1)[delta1 == 1], knots = knot.loc.3, 
                        degree = Bspline.degree.3,
                        intercept = TRUE, Boundary.knots = bdy.knots.b.3.y2my1)
  b.3.y2my1  <- predict(b.3.event, y2 - y1)
  
  start1 <- get.Bspline.startVals(y1, delta1, lin.pred[[1]], 
                                  dat_sim, b.1, knot.loc.1, Bspline.degree.1)
  phi1.start <- start1[1:(length(knot.loc.1) + Bspline.degree.1 + 1)]
  beta1.start <- start1[(length(phi1.start)+1):length(start1)]
  
  start2 <- get.Bspline.startVals(y2, delta2, lin.pred[[2]], 
                                  dat_sim, b.2, knot.loc.2, Bspline.degree.2)
  phi2.start <- start2[1:(length(knot.loc.2) + Bspline.degree.2 + 1)]
  beta2.start <- start2[(length(phi2.start)+1):length(start2)]
  
  start3 <- get.Bspline.startVals((y2 - y1)[delta1 == 1], delta2[delta1 == 1],
                                  lin.pred[[3]], 
                                  dat_sim[delta1 == 1, , drop = FALSE],
                                  b.3.y2my1, knot.loc.3, Bspline.degree.3)
  phi3.start <- start3[1:(length(knot.loc.3) + Bspline.degree.3 + 1)]
  beta3.start <- start3[(length(phi3.start)+1):length(start3)]
  
  startVals <- c(phi1.start, phi2.start, phi3.start)
  if (frailty) startVals <- c(startVals, 0.5) 
  startVals <- c(startVals, beta1.start, beta2.start, beta3.start)
  
  ## Fit B-spline model
  
  fit.bs <- FreqID.LT.bSpline123.v3(
    Y=Y, lin.pred=lin.pred, data=dat_sim, startVals=startVals, frailty=frailty,
    b.1=b.1, b.2=b.2, b.3.y2my1=b.3.y2my1,
    bdy.knots.b.1=bdy.knots.b.1, bdy.knots.b.2=bdy.knots.b.2,
    bdy.knots.b.3.y2my1=bdy.knots.b.3.y2my1
  )
  
  beta1 <- fit.bs$estimate[11:22]
  beta2 <- fit.bs$estimate[23:34]
  beta3 <- fit.bs$estimate[35:46]
  
  h1 <- as.matrix(dat_sim[, paste0("X", 1:12)]) %*% beta1
  h2 <- as.matrix(dat_sim[, paste0("X", 1:12)]) %*% beta2
  h3 <- as.matrix(dat_sim[, paste0("X", 1:12)]) %*% beta3
  
  theta <- rep(exp(fit.bs$estimate[10]), n_sim)
  
  phi11 <- rep(fit.bs$estimate[1], n_sim)
  phi12 <- rep(fit.bs$estimate[2], n_sim)
  phi13 <- rep(fit.bs$estimate[3], n_sim)
  
  phi21 <- rep(fit.bs$estimate[4], n_sim)
  phi22 <- rep(fit.bs$estimate[5], n_sim)
  phi23 <- rep(fit.bs$estimate[6], n_sim)
  
  phi31 <- rep(fit.bs$estimate[7], n_sim)
  phi32 <- rep(fit.bs$estimate[8], n_sim)
  phi33 <- rep(fit.bs$estimate[9], n_sim)
  
  sim_array[,, 1] <- cbind(theta, h1, h2, h3, 
                           phi11, phi12, phi13, 
                           phi21, phi22, phi23,
                           phi31, phi32, phi33)
}

save(sim_array, file = here("Results", paste0("Sim_06_",i,"_",k, ".RData")))

#=== END =======================================================================
  