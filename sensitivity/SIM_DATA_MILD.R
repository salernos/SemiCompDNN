#===============================================================================
#
#  PROGRAM: SIM_DATA_MILD.R
#
#  AUTHOR:  Stephen Salerno, Zhilin Zhang
#
#  PURPOSE: Generate additional data for sensitivity analysis
#
#  INPUT:   None
#
#  OUTPUT:  simdat{X+16}.RData
#
#  NOTES:   The simulated data are stored in ten files, with each file
#           corresponding to one simulation setting. Each file is a size
#           `n` x `31` x `nsims` array, storing `nsims` generated datasets
#           of size `n` x `31`.
#
#===============================================================================

#=== INITIALIZATION ============================================================

#--- NECESSARY PACKAGES --------------------------------------------------------

library(pacman)

p_load(MASS, tmvtnorm, patchwork, ggpmisc, here, tidyverse, update = F)

#=== HELPER FUNCTIONS ==========================================================

#--- DATA GENERATION -----------------------------------------------------------

simdat <- function(n, theta = 0.5, covariates = "normal", risk = "linear",

  cens = 0.25, cens_admin = 5, seed = NULL) {

  if (is.null(n)) stop('Error: argument "n" is missing, with no default')

  if (!is.null(seed)) set.seed(seed)

  if (theta == 0) {

    gamma <- rep(1, n)

  } else if (theta > 0) {

    gamma <- rgamma(n, 1/theta, 1/theta)
  }

  if (covariates == "uniform") {

    X <- matrix(runif(n*12, -1, 1), nrow = n, ncol = 12)

  } else if (covariates == "normal") {

    mu_X <- rep(0, 12)

    Sigma_X <- matrix(0.2, nrow = 12, ncol = 12); diag(Sigma_X) <- rep(1, 12)

    X <- rtmvnorm(n, mu_X, Sigma_X, rep(-2, 12), rep(2, 12))
  }

  colnames(X) <- paste0("X", 1:12)

  if (risk == "none") {

    h1 <- h2 <- h3 <- 0

  } else if (risk == "linear") {

    h1 <- X[, 1:9] %*% c(rep(-0.5, 4), rep(0.5, 5))

    h2 <- X[, 1:9] %*% c(rep(-0.5, 4), rep(0.5, 5))

    h3 <- X[, 1:9] %*% c(rep(-0.5, 4), rep(0.5, 5))

    if (cens == 0.25) {

      if (theta == 0.5) {

        muc <- 0.5

      } else if (theta == 2) {

        muc <- 0.05
      }

    } else if (cens == 0.50) {

      if (theta == 0.5) {

        muc <- 2

      } else if (theta == 2) {

        muc <- 1
      }
    }

  } else if (risk == "nonlinear") {

    h1 <- 0.5*exp(X[,1] - X[,2]) - 0.5*log((X[,3] + X[,4])^2) +

      0.5*sin(X[,5]*X[,6]) - 0.5*(X[,7] - X[,8] + X[,9])^2

    h2 <- 0.5*exp(X[,1] - X[,2]) - 0.5*log((X[,3] + X[,4])^2) +

      0.5*sin(X[,5]*X[,6]) - 0.5*(X[,7] - X[,8] + X[,9])^2

    h3 <- 0.5*exp(X[,1] - X[,2]) - 0.5*log((X[,3] + X[,4])^2) +

      0.5*sin(X[,5]*X[,6]) - 0.5*(X[,7] - X[,8] + X[,9])^2

    if (cens == 0.25) {

      if (theta == 0.5) {

        muc <- 0.1

      } else if (theta == 2) {

        muc <- 0.025
      }

    } else if (cens == 0.50) {

      if (theta == 0.5) {

        muc <- 2

      } else if (theta == 2) {

        muc <- 0.5
      }
    }
    
  } else if (risk == "mildnonlinear") {
    
    lin_beta  <- c(rep(-0.5, 4), rep(0.5, 5))
    lin_part  <- drop(X[, 1:9, drop = FALSE] %*% lin_beta)
    
    Z <- scale(X[, 1:9, drop = FALSE], center = TRUE, scale = TRUE)
    Z <- as.matrix(Z)
    
    int_mat  <- cbind(Z[,1] * Z[,2], Z[,3] * Z[,4], Z[,5] * Z[,6])
    int_w    <- c(0.15, -0.10, 0.10)
    int_term <- drop(int_mat %*% int_w)

    quad_w    <- c(rep(0.05, 4), rep(-0.05, 5))
    quad_term <- drop((Z^2) %*% quad_w)

    h1 <- lin_part + 0.20 * int_term + 0.10 * quad_term
    h2 <- lin_part + 0.15 * int_term + 0.10 * quad_term
    h3 <- lin_part + 0.20 * int_term + 0.05 * quad_term

    if (cens == 0.25) {
      
      if (theta == 0.5) {
        
        muc <- 0.5
        
      } else if (theta == 2) {
        
        muc <- 0.05
      }
      
    } else if (cens == 0.50) {
      
      if (theta == 0.5) {
        
        muc <- 2
        
      } else if (theta == 2) {
        
        muc <- 1
      }
    }
     
  }

  Ti1 <- -log(runif(n)) / (gamma*exp(-h1)) / 2
  Ti2 <- -log(runif(n)) / (gamma*exp(-h2)) / 3
  Ti3 <- -log(runif(n)) / (gamma*exp(-h3)) / 2

  ind <- Ti1 < Ti2

  Ti2[ind] <- Ti1[ind] + Ti3[ind]

  if (risk == "none") {

    hc <- 0

  } else {

    hc <- X[, 7:12] %*% c(-0.5, 0.5, -0.5, 0.5, -0.5, 0.5)
  }

  Ci <- pmin(-log(runif(n)) / (muc*exp(-hc)), cens_admin)

  Yi2 <- pmin(Ti2, Ci)
  Yi1 <- pmin(Ti1, Yi2)

  Di2 <- as.numeric(Ti2 <= Ci)
  Di1 <- as.numeric(Ti1 <= Ti2)

  ctime <- ifelse(Di1 == 0 & Di2 == 0, Yi1, pmax(Yi1, Yi2))
  cind  <- ifelse(Di1 == 0 & Di2 == 0, 1, 0)

  dat <- data.frame(

    cbind(Yi1, Di1, Yi2, Di2, X, h1, h2, h3, hc, gamma,

      Ti1, Ti2, Ti3, ctime, cind))

  colnames(dat) <- c("Yi1", "Di1", "Yi2", "Di2", paste0("X", 1:12),

    "h1", "h2", "h3", "hc", "gamma", "Ti1", "Ti2", "Ti3", "ctime", "cind")

  return(dat)
}

#=== GENERATE DATA =============================================================

#--- DATA GENERATING MECHANISMS ------------------------------------------------

nsims <- 500

n <- c(1000, 5000)

theta <- 2.0

risk <- "mildnonlinear"

cens <- 0.25

settings <- expand.grid(n, theta, risk, cens)

#--- LOOP OVER SETTINGS --------------------------------------------------------

for (i in 1:nrow(settings)) {

  cat("Setting", i, "--\n")

  dat <- array(dim = c(settings[i, 1], 31, nsims))

  for (sim in 1:nsims) {

    if (sim %% 50 == 0) {

      cat("Simulated Dataset", sim, "\n")
    }

    dat_sim <- simdat(n = settings[i, 1], theta = settings[i, 2],

      risk = settings[i, 3], cens = settings[i, 4], seed = sim) |>

      mutate(setting = i, n = settings[i, 1], theta = settings[i, 2],

        risk = settings[i, 3], cens = settings[i, 4],

        risk = case_when(
          risk == "linear"         ~ 0,
          risk == "mildnonlinear"  ~ 1,
          risk == "nonlinear"      ~ 2
        )) |>

      as.matrix()

    dat[, , sim] <- dat_sim
  }

  save(dat, file = here("Datanew", paste0("simdat", i+16, ".RData")))
}

#--- EXAMPLE DATASETS ----------------------------------------------------------

ex_dat <- vector("list", nrow(settings))

for (i in 1:nrow(settings)) {

  load(here("Data", paste0("simdat", i, ".RData")))

  dat_i <- dat[,,1]

  colnames(dat_i) <- c("Yi1", "Di1", "Yi2", "Di2", paste0("X", 1:12),

    "h1", "h2", "h3", "hc", "gamma", "Ti1", "Ti2", "Ti3", "ctime", "cind",

    "setting", "n", "theta", "risk", "cens")

  ex_dat[[i]] <- dat_i |> as_tibble() |>

    rename("Sample Size" = n, "Frailty Variance" = theta, "Log Risk" = risk,

      "Censoring Rate" = cens)
}

ex_dat_all <- bind_rows(ex_dat)

ex_dat_summ <- ex_dat_all |>

  group_by(setting, `Sample Size`, `Frailty Variance`, `Log Risk`,

    `Censoring Rate`, Di1, Di2) |>

  count(name = "Count") |> group_by(setting) |>

  mutate(Percent = Count / sum(Count) * 100, name = "Yi1") |> ungroup() |>

  mutate(x = 1.5, y = rep(seq(0.5, 0.875, 0.125), 16),

    lab = paste0("Di1 =", Di1, ", Di2 =", Di2, " : ", Percent, "%"))

#--- PLOT EXAMPLE DATASETS -----------------------------------------------------

ex_plt <- ex_dat_all |>

  pivot_longer(c(Yi1, Yi2)) |>

  ggplot(aes(x = value, group = name, color = name, fill = name)) +

    theme_bw() +

    facet_grid(`Sample Size` + `Frailty Variance` ~

      `Log Risk` + `Censoring Rate`, labeller = label_both) +

    geom_density(aes(y = after_stat(scaled)), alpha = 0.25) +

    geom_text(aes(x = x, y = y, label = lab), data = ex_dat_summ,

      color = "black", hjust = 0) +

    scale_color_manual(values = c("#FFCB05", "#00274C")) +

    scale_fill_manual(values = c("#FFCB05", "#00274C")) +

    labs(x = "Simulated Observation Time", y = "Scaled Density",

      color = "Observation Time:", fill = 'Observation Time:') +

    theme(

      legend.position = "bottom",

      strip.background = element_rect(fill = "#00274C", color = "white"),

      strip.text = element_text(face = "bold", color = "white", size = 12))

ggsave(here("Results", "simulated_data_example.png"), ex_plt,

  width = 16, height = 16)

#=== END =======================================================================
