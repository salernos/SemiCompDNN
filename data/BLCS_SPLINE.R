#===============================================================================
#
#  PROGRAM: BLCS_SPLINE.R
#
#  AUTHORS: Stephen Salerno, Zhilin Zhang
#
#  PURPOSE: Real-data analysis of the BLSC data with the spline-based method
#
#  INPUT:   BLCS_Data_Updated_February2023.xlxs
#
#  NOTES: 
#
#  We cannot provide the real data here, however we provide the codes here for
#  review and reproducibility with similarly structured data. This script
#  implements the real data analysis for the spline-based method.
#
#===============================================================================

#=== INITIALIZATION ============================================================

#--- LOAD NECESSARY PACKAGES ---------------------------------------------------

library(pacman)

p_load(gtsummary, survival, survcomp, survminer, ggsci, patchwork,
       
  readxl, lubridate, here, tidyverse, SemiCompRisks, update = F)

#--- HELPER FUNCTIONS ----------------------------------------------------------

source(here("R", "Simulation", "METRICS1.R"))

bcindex <- function(dat, h1, h2, h3) {

  Y1 <- dat$Yi1 
  Y2 <- dat$Yi2 
  
  D1 <- dat$Di1  
  D2 <- dat$Di2  
  
  prank_h1 <- ecdf(h1)(h1)
  prank_h2 <- ecdf(h2)(h2)
  prank_h3 <- ecdf(h3)(h3)

  qval_h1 <- qnorm(prank_h1)
  qval_h2 <- qnorm(prank_h2)
  qval_h3 <- qnorm(prank_h3)

  risk <- (qval_h1 + qval_h2 + qval_h3) / 3
  
  numer <- 0
  denom <- 0
  n <- length(Y1)
  
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      if ((D1[i] == 1 | D2[i] == 1) & (D1[j] == 1 | D2[j] == 1)) {
        denom <- denom + as.numeric(Y1[i] < Y1[j] |Y2[i] < Y2[j])
        numer <- numer + as.numeric((Y1[i] < Y1[j] & risk[i] > risk[j]) | 
                                      (Y2[i] < Y2[j] & risk[i] > risk[j]))
      }
    }
  }
  
  bcidx <- numer / denom
  
  return(bcidx)
}

#=== DATA ======================================================================

#--- READ-IN DATA --------------------------------------------------------------

BLCS <- read_excel("./BLCS/BLCS_Data_Updated_February2023.xlsx")

#--- PRE-PROCESS DATA ----------------------------------------------------------

extrcdate <- ymd("2023-02-28")

BLCS_OUT1 <- BLCS |>
  
  #-- Correct Variable Types and Missingness
  
  mutate(
    
    across(where(is.character), ~na_if(., "")),
    
    across(contains("date"), ~as_date(.))) |>
  
  #-- Remove 58 Patients Without Diagnosis Date
  
  filter(!is.na(dxdate))

BLCS_OUT2 <- BLCS_OUT1 |>
  
  #-- Remove 56 Patients Diagnosed < 6 Months
  
  filter(difftime(extrcdate, dxdate, units = "days") > 180)

BLCS_OUT3 <- BLCS_OUT2 |>
  
  #-- Remove 212 Patients with SCLC and 6 Patients with Stage 0
  
  filter(!(cstage %in% c("0", "SCLC-ES", "SCLC-LS")))

BLCS_OUT4 <- BLCS_OUT3 |>
  
  mutate(
    
    #- Progression
    
    pdate = case_when(
      
      !is.na(progressiondate) ~ progressiondate,
      
      !is.na(relapsedate) ~ relapsedate,
      
      !is.na(progdate) & progdate > dxdate ~ progdate,
      
      T ~ NA_Date_),
    
    pevent = case_when(
      
      !is.na(progression) ~ progression,
      
      !is.na(relapse) ~ relapse,
      
      !is.na(prog) ~ prog,
      
      T ~ NA_real_),
    
    pevent = case_when(
      
      pevent %in% c(0, 2) ~ 0,
      
      pevent == 1 ~ 1,
      
      T ~ NA_real_),
    
    #- Censoring
    
    censordate = case_when(
      
      !(is.na(alivedate)) & alivedate <= dxdate &
        
        !(is.na(pdate) & is.na(deathdate)) ~ extrcdate,
      
      !(is.na(alivedate)) & !(is.na(deathdate)) &
        
        alivedate <= deathdate ~ extrcdate,
      
      T ~ pmin(alivedate, extrcdate, na.rm = T)),
    
    censorcheck = case_when(
      
      !(is.na(alivedate)) & alivedate <= dxdate &
        
        !(is.na(pdate) & is.na(deathdate)) ~ "Corrected Negative Survival",
      
      !(is.na(alivedate)) & !(is.na(deathdate)) &
        
        alivedate <= deathdate ~ "Corrected Censoring with Death",
      
      T ~ "Not Corrected"),
    
    #- Univariate Outcomes
    
    Ti1 = as.duration(dxdate %--% pdate)      / dyears(1),
    
    Ti2 = as.duration(dxdate %--% deathdate)  / dyears(1),
    
    Ci  = as.duration(dxdate %--% censordate) / dyears(1),
    
    ptime  = pmin(Ti1, Ci, na.rm = T),
    
    dtime  = pmin(Ti2, Ci, na.rm = T),
    
    pevent = if_else(Ti1 <= Ci, 1, 0, 0),
    
    devent = if_else(Ti2 <= Ci, 1, 0, 0),
    
    #- Semi-Competing Outcomes
    
    Yi1 = pmin(Ti1, Ti2, Ci, na.rm = T),
    
    Yi2 = pmin(Ti2, Ci, na.rm = T),
    
    YiS = pmin(Ti1, Ti2, Ci, na.rm = T),
    
    Di1 = if_else(Ti1 <= Yi2, 1, 0, 0),
    
    Di2 = if_else(Ti2 <= Ci,  1, 0, 0),
    
    DiS = if_else(pmin(Ti2, Ti2, na.rm = T) <= Ci, 1, 0, 0),
    
    #- Predictors
    
    sex	= factor(sex, 1:2, c("Male", "Female")),
    
    race = factor(race, 1:7,
                  
      c("White/Caucasian", "Native American/Alaska Native", "Asian",
        "Black/African American", "Native Hawaiian/Pacific Islander",
        "Multiracial", "Other")),
    
    ethnic = factor(ethnic, 0:1, c("Non-Hispanic", "Hispanic")),
    
    education = factor(education, 1:8,
                       
      c("Some grade school", "Some high school", "High school graduate",
        "Vocational/tech school after high school",
        "Some college or associate's degree",
        "College graduate", "Graduate or professional school", "Other")),
    
    bmi = wtkg / htm^2,
    
    smk	= factor(smk, 1:4,
                 
      c("Never smoker", "Former smoker", "Current smoker",
        "Smoker, status unknown")),
    
    trt = case_when(
      
      surg     == 1 ~ "Surgery",
      chemo    == 1 ~ "Chemotherapy",
      radio    == 1 ~ "Radiation",
      othertrt == 1 ~ "Other"),
    
    copd = factor(copd, 0:1, c("No", "Yes")),
    
    asthma = factor(asthma, 0:1, c("No", "Yes")),
    
    egfr = factor(egfr, 0:1, c("No", "Yes")) |>
      
      fct_na_value_to_level("Not Tested"),
    
    kras = factor(kras, 0:1, c("No", "Yes")) |>
      
      fct_na_value_to_level("Not Tested")) |>
  
  select(ptime, pevent, dtime, devent, Yi1, Yi2, YiS, Di1, Di2, DiS, agedx, sex,
         
    race, ethnic, education, bmi, smk, pkyrs, trt, ctype, cstage, copd, asthma,
         
    egfr, kras, everything())

#-- QA: PATIENTS WITH NEGATIVE SURVIVAL TIMES

BLCS_CHECK <- BLCS_OUT4 |>
  
  filter(Yi1 <= 0 | Yi2 <= 0) |>
  
  select(STIKRNUM, dxdate, relapsedate, progressiondate, progdate,
         
         nonprogressiondate, deathdate, alivedate, Yi1, Yi2)

#-- FINAL PRE-PROCESSED DATASET ------------------------------------------------

BLCS_CLEAN <- BLCS_OUT4 |>
  
  #-- Remove 25 patients with negative survival times
  
  filter(Yi1 > 0, Yi2 > 0)

#-- DATASET RE-FACTORED FOR MODEL FITTING --------------------------------------

BLCS_CLEAN2 <- BLCS_CLEAN |>
  
  mutate(
    
    stage_splt = if_else(cstage %in%
                           
      c("1", "1A", "1B", "2", "2A", "2B", "3", "3A"), 0, 1, 0),
    
    agedx = if_else(is.na(agedx), mean(agedx, na.rm = T), agedx),
    
    sex = sex |> fct_na_value_to_level("Unknown"),
    
    race = race |> fct_na_value_to_level("Unknown") |>
      
      fct_collapse(
        
        "White/Caucasian" = "White/Caucasian",
        
        "Other" = c(
          
          "Native American/Alaska Native", "Asian", "Black/African American",
          "Native Hawaiian/Pacific Islander", "Multiracial", "Other")),
    
    ethnic = ethnic |> fct_na_value_to_level("Unknown"),
    
    education = education |> fct_na_value_to_level("Unknown") |>
      
      fct_collapse("Other" = c("Other", "Unknown")),
    
    bmi = if_else(is.na(bmi), mean(bmi, na.rm = T), bmi),
    
    smk = smk |> fct_na_value_to_level("Unknown") |>
      
      fct_collapse(
        
        "Smoker, status unknown" = c("Smoker, status unknown", "Unknown")),
    
    pkyrs = if_else(is.na(pkyrs), mean(pkyrs, na.rm = T), pkyrs),
    
    surg = if_else(trt == "Surgery", 1, 0, 0)) |>
  
  select(Yi1, Yi2, YiS, Di1, Di2, DiS, Ti1, Ti2, Ci, stage_splt, cstage, trt,
         
         agedx:pkyrs, egfr:kras)

data <- BLCS_CLEAN |>
  
  mutate(
    
    stage_splt = if_else(cstage %in%
                           
      c("1", "1A", "1B", "2", "2A", "2B", "3", "3A"), 0, 1, 0),
    
    stage_splt = as.factor(stage_splt),
    
    agedx = if_else(is.na(agedx), mean(agedx, na.rm = T), agedx),
    
    sex = sex |> fct_na_value_to_level("Unknown"),
    
    race = race |> fct_na_value_to_level("Unknown") |>
      
      fct_collapse(
        
        "White/Caucasian" = "White/Caucasian",
        
        "Asian" = "Asian",
        
        "Black/African American" = "Black/African American",
        
        "Other" = c(
          
          "Native American/Alaska Native",
          "Native Hawaiian/Pacific Islander", "Multiracial", "Other")),
    
    ethnic = ethnic |> fct_na_value_to_level("Unknown"),
    
    education = education |> fct_na_value_to_level("Unknown") |>
      
      fct_collapse("Other" = c("Other", "Unknown")),
    
    bmi = if_else(is.na(bmi), mean(bmi, na.rm = T), bmi),
    
    smk = smk |> fct_na_value_to_level("Unknown") |>
      
      fct_collapse(
        
        "Smoker, status unknown" = c("Smoker, status unknown", "Unknown")),
    
    pkyrs = if_else(is.na(pkyrs), mean(pkyrs, na.rm = T), pkyrs),
    
    surg = if_else(trt == "Surgery", 1, 0, 0),
    
    ctype = as.factor(ctype) |> fct_na_value_to_level("Unknown"),
    
    trt = trt |> fct_na_value_to_level("Unknown"),
    
    copd = copd |> fct_na_value_to_level("Unknown"),
    
    asthma = asthma |> fct_na_value_to_level("Unknown")) |>
  
  select(Yi1, Yi2, Di1, Di2, stage_splt, ctype, trt, 
         
    agedx:pkyrs, egfr:kras, copd, asthma)

X <- model.matrix(~ agedx + sex + race + ethnic + copd + asthma + stage_splt +
                    trt + smk + pkyrs + egfr + kras,data = data)

data <- cbind(data, X)

form <- Formula(~ agedx + sexFemale + sexUnknown + raceOther + raceAsian + 
                  `raceBlack/African American` + raceUnknown + ethnicHispanic + 
                  ethnicUnknown + copdYes + copdUnknown + asthmaYes + 
                  asthmaUnknown + stage_splt1 + trtOther + trtRadiation + 
                  trtSurgery + trtUnknown + `smkFormer smoker` + 
                  `smkCurrent smoker` + `smkSmoker, status unknown` + pkyrs + 
                  egfrYes + `egfrNot Tested` + krasYes)

source(here("./R/functions.R")) 

K      <- 5

folds  <- createFolds(data$Di1, k = K)

t_pred <- seq(0, 5, length.out = 100)

bs_mat   <- matrix(NA, nrow = length(t_pred), ncol = K)
ibs_vec  <- numeric(K)

cind <- numeric(K)

for (k in 3:5) {
  
  test_idx  <- folds[[k]]
  train_dat <- data[-test_idx, ]
  test_dat  <- data[ test_idx, ]
  X_train <- X[-test_idx,-c(1,27)]
  X_test <- X[test_idx,-c(1,27)]
  
  Y <- data.frame(
    y1     = train_dat$Yi1,
    delta1 = train_dat$Di1,
    y2     = train_dat$Yi2,
    delta2 = train_dat$Di2,
    L      = rep(0,dim(train_dat)[1])
  )
  
  Y$L <- fix_L(Y$L, Y$y1)
  
  lin.pred <- list(
    form,  
    form,  
    form  
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

  start1 <- get.Bspline.startVals(y1, delta1, lin.pred[[1]], train_dat, b.1, 
                                  knot.loc.1, Bspline.degree.1)
  phi1.start <- start1[1:(length(knot.loc.1) + Bspline.degree.1 + 1)]
  beta1.start <- start1[(length(phi1.start)+1):length(start1)]
  
  start2 <- get.Bspline.startVals(y2, delta2, lin.pred[[2]], train_dat, b.2, 
                                  knot.loc.2, Bspline.degree.2)
  phi2.start <- start2[1:(length(knot.loc.2) + Bspline.degree.2 + 1)]
  beta2.start <- start2[(length(phi2.start)+1):length(start2)]
  
  start3 <- get.Bspline.startVals((y2 - y1)[delta1 == 1], delta2[delta1 == 1],
                                  lin.pred[[3]], 
                                  train_dat[delta1 == 1, , drop = FALSE],
                                  b.3.y2my1, knot.loc.3, Bspline.degree.3)
  phi3.start <- start3[1:(length(knot.loc.3) + Bspline.degree.3 + 1)]
  beta3.start <- start3[(length(phi3.start)+1):length(start3)]
  
  startVals <- c(phi1.start, phi2.start, phi3.start)
  if (frailty) startVals <- c(startVals, 0.5)  
  startVals <- c(startVals, beta1.start, beta2.start, beta3.start)
  
  fit.bs <- FreqID.LT.bSpline123.v3(
    Y=Y, lin.pred=lin.pred, data=train_dat, startVals=startVals, frailty=T,
    b.1=b.1, b.2=b.2, b.3.y2my1=b.3.y2my1,
    bdy.knots.b.1=bdy.knots.b.1, bdy.knots.b.2=bdy.knots.b.2,
    bdy.knots.b.3.y2my1=bdy.knots.b.3.y2my1
  )
  
  beta1 <- fit.bs$estimate[11:35]
  beta2 <- fit.bs$estimate[36:60]
  beta3 <- fit.bs$estimate[61:85]
  
  theta <- exp(fit.bs$estimate[10])
  
  phi11 <- fit.bs$estimate[1]
  phi12 <- fit.bs$estimate[2]
  phi13 <- fit.bs$estimate[3]
  
  phi21 <- fit.bs$estimate[4]
  phi22 <- fit.bs$estimate[5]
  phi23 <- fit.bs$estimate[6]
  
  phi31 <- fit.bs$estimate[7]
  phi32 <- fit.bs$estimate[8]
  phi33 <- fit.bs$estimate[9]
  
  h1 <- as.matrix(X_test) %*% beta1
  h2 <- as.matrix(X_test) %*% beta2
  h3 <- as.matrix(X_test) %*% beta3
  
  phi.1.hat <- unlist(list(phi11,phi12,phi13))
  
  phi.2.hat <- unlist(list(phi21,phi22,phi23))
  
  phi.3.hat <- unlist(list(phi31,phi32,phi33))
  
  preds <- sapply(t_pred, function(t) {
                         
    (1 + theta*exp(as.vector(matrix(phi.1.hat, nrow=1) %*% 
                               
      t(as.matrix(predict(b.1.event, t)))))*exp(h1) +
       
       theta*exp(as.vector(matrix(phi.2.hat, nrow=1) %*% 
                             
          t(as.matrix(predict(b.2.event, t)))))*exp(h2))^(-1/theta)
  })
  
  res        <- bbs1(test_dat, preds, t_pred)
  bs_mat[,k] <- res$bs
  ibs_vec[k] <- res$ibs
  
  test_dat$gamma <- rgamma(dim(test_dat)[1], 1/theta, 1/theta)
  
  cind[k] <- bcindex(test_dat,h1,h2,h3)
}

mean_bs <- rowMeans(bs_mat)
sd_bs   <- apply(bs_mat, 1, sd)
se_bs   <- sd_bs / sqrt(K)
t_q     <- qt(0.975, df = K-1)

ci_lower <- mean_bs - t_q * se_bs
ci_upper <- mean_bs + t_q * se_bs

mean(cind)
mean(cind) - t_q * sd(cind)
mean(cind) + t_q * sd(cind)

mean(ibs_vec)
mean(ibs_vec) - t_q * sd(ibs_vec)
mean(ibs_vec) + t_q * sd(ibs_vec)

plot_df <- data.frame(
  time  = t_pred,
  mean  = mean_bs,
  lower = ci_lower,
  upper = ci_upper
)

mean_ibs <- mean(ibs_vec)
se_ibs   <- sd(ibs_vec) / sqrt(K)
ci_ibs   <- mean_ibs + c(-1,1) * t_q * se_ibs

cat(sprintf("iBBS = %.4f (95%% CI: %.4f – %.4f)\n", 
            mean_ibs, ci_ibs[1], ci_ibs[2]))

#=== END =======================================================================