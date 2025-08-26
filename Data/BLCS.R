#===============================================================================
#
#  PROGRAM: BLCS.R
#
#  AUTHORS: Stephen Salerno, Zhilin Zhang
#
#  PURPOSE: Real-data analysis of the Boston Lung Cancer Study for Methods 0-4.
#
#  INPUT:   BLCS_Data_Updated_February2023.xlxs
#
#  OUTPUT:  
#
#  Descriptive Analysis
#
#    - Tab. 1: Distribution of Observed Events
#    - Fig. 1: Distribution of Observation Times
#    - Tab. 2: Descriptive Statistics for All Covariates
#    - Tab. 3: Descriptive Statistics for Refactored Model Covariates
#    - Tab. 4: Distribution of Observed Events by Stage
#
#  Modeling: 
#
#    - Method 00: Xu et al. (2010)
#    - Method 01: Lee et al. (2015)
#    - Method 02: Lee et al. (2017)
#    - Method 03: Gorfine et al. (2020)
#    - Method 04: Kats et al. (2022)
#
#  NOTES: 
#
#  We cannot provide the real data here, however we provide the codes here for
#  review and reproducibility with similarly structured data. This script
#  implements the real data analysis for the five methods we compare to, all 
#  of which are implemented in R.
#
#===============================================================================

#=== INITIALIZATION ============================================================

#--- LOAD NECESSARY PACKAGES ---------------------------------------------------

library(pacman)

p_load(

  #-- Packages for Analysis

  SemiCompRisks, caret,

  #-- Packages for Producing Tables and Figures

  gtsummary, survival, survcomp, survminer, ggsci, patchwork,

  #-- Packages for Reading and Pre-Processing Data

  readxl, lubridate, here, tidyverse,

  #-- Do Not Update Packages for Version Control

  update = F)

source(here("frailty-LTRC-master", "Estimation Procedure.R"))

source(here("semicompAFT-main", "semicompAFT.R"))

#=== DATA ======================================================================

#--- READ-IN DATA --------------------------------------------------------------

BLCS <- read_excel("./R/BLCS/BLCS_Data_Updated_February2023.xlsx")

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

      c("White/Caucasian", 
        "Native American/Alaska Native", 
        "Asian",
        "Black/African American", 
        "Native Hawaiian/Pacific Islander",
        "Multiracial", "Other")),

    ethnic = factor(ethnic, 0:1, c("Non-Hispanic", "Hispanic")),

    education = factor(education, 1:8,

      c("Some grade school", 
        "Some high school", 
        "High school graduate",
        "Vocational/tech school after high school",
        "Some college or associate's degree",
        "College graduate", 
        "Graduate or professional school", 
        "Other")),

    bmi = wtkg / htm^2,

    smk	= factor(smk, 1:4,

      c("Never smoker", 
        "Former smoker", 
        "Current smoker",
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

  dplyr::select(ptime, pevent, dtime, devent, Yi1, Yi2, YiS, Di1, Di2, DiS, 
                
    agedx, sex, race, ethnic, education, bmi, smk, pkyrs, trt, ctype, cstage, 
    
    copd, asthma, egfr, kras, everything())

#-- QA: PATIENTS WITH NEGATIVE SURVIVAL TIMES

BLCS_CHECK <- BLCS_OUT4 |>

  filter(Yi1 <= 0 | Yi2 <= 0) |>

  dplyr::select(STIKRNUM, dxdate, relapsedate, progressiondate, progdate,

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

          "Native American/Alaska Native", 
          "Asian", 
          "Black/African American",
          "Native Hawaiian/Pacific Islander", 
          "Multiracial", 
          "Other")),

    ethnic = ethnic |> fct_na_value_to_level("Unknown"),

    education = education |> fct_na_value_to_level("Unknown") |>

      fct_collapse("Other" = c("Other", "Unknown")),

    bmi = if_else(is.na(bmi), mean(bmi, na.rm = T), bmi),

    smk = smk |> fct_na_value_to_level("Unknown") |>

      fct_collapse(

        "Smoker, status unknown" = c("Smoker, status unknown", "Unknown")),

    pkyrs = if_else(is.na(pkyrs), mean(pkyrs, na.rm = T), pkyrs),

    surg = if_else(trt == "Surgery", 1, 0, 0)) |>

  dplyr::select(Yi1, Yi2, YiS, Di1, Di2, DiS, stage_splt, cstage, trt,

    agedx:pkyrs, egfr:kras)

#=== DESCRIPTIVE ANALYSIS ======================================================

#--- TABLE 1: OUTCOME DISTRIBUTION ---------------------------------------------

TAB1 <- BLCS_CLEAN |>

  mutate(

    Di1 = factor(Di1, 0:1,

      c("Progression Not Observed", "Progression Observed")),

    Di2 = factor(Di2, 0:1,

      c("Death not Observed", "Death Observed"))) |>

  group_by(Di1, Di2) |>

  count() |> ungroup() |>

  mutate(Percent = round(n / sum(n) * 100, 2)) |>

  gt::gt() |>

  gt::tab_caption("Distribution of Observed Events")

TAB1

#--- FIGURE 1: OUTCOME DISTRIBUTION --------------------------------------------

FIG1 <- BLCS_CLEAN |>

  mutate(

    cases = interaction(Di1, Di2) |>

      factor(labels = c("Both Events Censored", 
                        "Progression Observed",
                        "Death Observed", 
                        "Progression and Death Observed")),

    cens = case_when(cases == "Both Events Censored"           ~ 0.25,
                     cases == "Progression Observed"           ~ 1.0,
                     cases == "Death Observed"                 ~ 1.0,
                     cases == "Progression and Death Observed" ~ 1.0)) |>

  sample_frac(1L) |>

  ggplot(aes(x = Yi1, y = Yi2, color = cases, fill = cases)) +

  theme_bw() + coord_fixed() +

  geom_abline(slope = 1, intercept = 0, linetype = 2) +

  geom_point(size = 3, alpha = 0.5, stroke = 0.1, shape = 21) +

  scale_shape_manual(values = c(4,3,16,17)) +

  scale_color_manual(values = c("#989C97", "#CFC096", "#9A3324", "#00274C")) +

  scale_fill_manual(values = c("#989C97", "#CFC096", "#9A3324", "#00274C")) +

  labs(

    color = "Observed Data:",
    fill  = "Observed Data:",

    x = expression(paste(y[i1], "  (yrs.)")),
    y = expression(paste(y[i2], "  (yrs.)"))) +

  theme(

    legend.position = c(0.7, 0.2))

FIG1

#--- TABLE 2: DESCRIPTIVE CHARACTERISTICS --------------------------------------

tab2_labs <- list(

  agedx	    ~ "Age at Diagnosis (yrs.)",
  sex       ~ "Sex",
  race      ~ "Race",
  ethnic    ~ "Ethnicity",
  education ~ "Education",
  bmi       ~ "Body Mass Index",
  smk	      ~ "Smoking Status",
  pkyrs	    ~ "Pack-Years of Smoking",
  trt       ~ "First-Line Treatment",
  ctype     ~ "Histologic Type",
  cstage    ~ "Cancer Stage",
  copd	    ~ "COPD",
  asthma	  ~ "Asthma",
  egfr	    ~ "EGFR Mutation",
  kras	    ~ "KRAS Mutation")

# Define the categorization
BLCS_CLEAN <- BLCS_CLEAN %>%
  mutate(Stage_Group = case_when(
    cstage %in% c("1", "1A", "1B", "2", "2A", "2B", "3", "3A") ~ "Early Stage",
    cstage %in% c("3B", "4", "5", "6", "6A", "6B") ~ "Late Stage",
    TRUE ~ "Unknown"
  ))

TAB2 <- BLCS_CLEAN |>

  dplyr::select(agedx:kras) |>

  tbl_summary(label = tab2_labs) |>

  as_gt() |>

  gt::tab_caption("Descriptive Statistics")

TAB2 <- BLCS_CLEAN |>

  dplyr::select(Stage_Group, agedx:kras) |>

  tbl_summary(by = Stage_Group, label = tab2_labs) |>

  as_gt() |>

  gt::tab_caption("Descriptive Statistics")


TAB2

#--- TABLE 3: DESCRIPTIVE CHARACTERISTICS FOR MODEL DATA -----------------------

with(BLCS_CLEAN2, table(stage_splt, trt))

with(BLCS_CLEAN2, table(cstage, trt))

tab3_labs <- list(

  cstage    ~ "Stage at Diagnosis",
  agedx	    ~ "Age at Diagnosis (yrs.)",
  sex       ~ "Sex",
  race      ~ "Race",
  ethnic    ~ "Ethnicity",
  education ~ "Education",
  bmi       ~ "Body Mass Index",
  smk	      ~ "Smoking Status",
  pkyrs	    ~ "Pack-Years of Smoking",
  # trt       ~ "First-Line Treatment",
  egfr	    ~ "EGFR Mutation",
  kras	    ~ "KRAS Mutation")

TAB3 <- BLCS_CLEAN2 |>

  mutate(stage_splt = factor(stage_splt, 0:1, c("1-3A", "3B-4"))) |>

  dplyr::select(stage_splt, trt:kras) |>

  tbl_summary(stage_splt, label = tab3_labs) |>

  add_overall() |>

  as_gt() |>

  gt::tab_caption("Descriptive Statistics for Model Data, Stratified by Stage")

TAB3

TAB3_ALT <- BLCS_CLEAN2 |>

  filter(stage_splt == 0) |>

  mutate(

    cstage = factor(cstage) |>

      fct_collapse("Stage 1" = c("1", "1A", "1B"),

        "Stage 2" = c("2", "2A", "2B"),

        "Stage 3" = c("3", "3A")),

    trt = if_else(trt == "Surgery", "Surgery",

      "Chemotherapy, Radiation, or Other",

      "Chemotherapy, Radiation, or Other")) |>

  dplyr::select(cstage, trt:kras) |>

  tbl_summary(trt, label = tab3_labs) |>

  add_overall() |> bold_labels() |>

  add_p() |> bold_p() |>

  as_gt() |>

  gt::tab_caption("Descriptive Statistics for Model Data, Stratified by Stage")

TAB3_ALT

#--- TABLE 4: OUTCOME DISTRIBUTION BY STAGE ------------------------------------

TAB4 <- BLCS_CLEAN2 |>

  mutate(

    stage_splt = factor(stage_splt, 0:1,

      c("Stage 1-3A", "Stage 3B-4")),

    Di1 = factor(Di1, 0:1,

      c("Progression Not Observed", "Progression Observed")),

    Di2 = factor(Di2, 0:1,

      c("Death not Observed", "Death Observed"))) |>

  group_by(stage_splt, Di1, Di2) |>

  count() |> group_by(stage_splt) |>

  mutate(Percent = round(n / sum(n) * 100, 2)) |>

  gt::gt() |>

  gt::tab_caption("Distribution of Observed Events")

TAB4

#=== ANALYSIS ==================================================================

#--- ANALYTIC DATASET ----------------------------------------------------------

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

        "Other" = c(

          "Native American/Alaska Native", 
          "Asian", "Black/African American",
          "Native Hawaiian/Pacific Islander", 
          "Multiracial", 
          "Other")),

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

  dplyr::select(Yi1, Yi2, Di1, Di2, stage_splt, ctype, trt,

         agedx:pkyrs, egfr:kras, copd, asthma)

X <- model.matrix(~ agedx + sex + race + ethnic + copd + asthma + 
                    
  stage_splt + trt + smk + pkyrs + egfr + kras, data = data)

data <- cbind(data, X)

form <- Formula(
  
  #-- Outcomes
  
  Yi1 + Di1 | Yi2 + Di2 ~ 
    
  #-- Covariates
    
  #- Transition 1
                  
  agedx + sexFemale + sexUnknown + raceOther + raceUnknown + ethnicHispanic + 
    
    ethnicUnknown + copdYes + copdUnknown + asthmaYes + asthmaUnknown + 
    
    stage_splt1 + trtOther + trtRadiation + trtSurgery + trtUnknown + 
    
    `smkFormer smoker` + `smkCurrent smoker` + `smkSmoker, status unknown` + 
    
    pkyrs + egfrYes + `egfrNot Tested` + krasYes |
    
  #- Transition 2
                  
  agedx + sexFemale + sexUnknown + raceOther + raceUnknown + ethnicHispanic + 
    
    ethnicUnknown + copdYes + copdUnknown + asthmaYes + asthmaUnknown + 
    
    stage_splt1 + trtOther + trtRadiation + trtSurgery + trtUnknown + 
    
    `smkFormer smoker` + `smkCurrent smoker` + `smkSmoker, status unknown` + 
    
    pkyrs + egfrYes + `egfrNot Tested` + krasYes |
  
  #- Transition 3                
    
  agedx + sexFemale + sexUnknown + raceOther + raceUnknown + ethnicHispanic + 
    
    ethnicUnknown + copdYes + copdUnknown + asthmaYes + asthmaUnknown + 
    
    stage_splt1 + trtOther + trtRadiation + trtSurgery + trtUnknown + 
    
    `smkFormer smoker` + `smkCurrent smoker` + `smkSmoker, status unknown` + 
    
    pkyrs + egfrYes + `egfrNot Tested` + krasYes)

#--- METHOD 00 -----------------------------------------------------------------

fit0 <- nlm(logLike.weibull.SCR, 
            
  p = startVals * runif(length(startVals), 0.9, 1.1),

  y1 = y1, delta1 = delta1, y2 = y2, delta2 = delta2,

  Xmat1 = as.matrix(Xmat1), 
  
  Xmat2 = as.matrix(Xmat2), 
  
  Xmat3 = as.matrix(Xmat3),

  frailty = frailty,

  iterlim = 1000,  # NOTE: This can be changed, but it is not an issue

  hessian = TRUE,

  steptol = 1e-5)  # NOTE: Adding this argument to allow larger step sizes

fit0

set.seed(12345)

#- Prepare folds and prediction grid

K      <- 5
folds  <- createFolds(data$Di1, k = K)
t_pred <- seq(0, 5, length.out = 100)

# Storage

bs_mat  <- matrix(NA, nrow = length(t_pred), ncol = K)
ibs_vec <- numeric(K)

for (k in seq_along(folds)) {

  # Split train/test
  
  test_idx  <- folds[[k]]
  train_dat <- data[-test_idx, ]
  test_dat  <- data[ test_idx, ]
  X_train   <- X[-test_idx,-c(1,25)]
  X_test    <- X[test_idx,-c(1,25)]

  y1     <- as.vector(train_dat$Yi1)
  delta1 <- as.vector(train_dat$Di1)
  y2     <- as.vector(train_dat$Yi2)
  delta2 <- as.vector(train_dat$Di2)

  #- Covariates

  Xmat1 <- X_train
  Xmat2 <- X_train
  Xmat3 <- X_train

  #- Starting Values

  fit.survreg.1 <- survreg(as.formula(paste("Surv(y1, delta1) ",

    as.character(formula(form, lhs=0, rhs=1))[1],
    as.character(formula(form, lhs=0, rhs=1))[2])), 
    dist="weibull", data=train_dat)

  fit.survreg.2 <- survreg(as.formula(paste("Surv(y2, delta2) ",

    as.character(formula(form, lhs=0, rhs=2))[1],
    as.character(formula(form, lhs=0, rhs=2))[2])), 
    dist="weibull", data=train_dat)

  fit.survreg.3 <- survreg(as.formula(paste("Surv(y2, delta2) ",

    as.character(formula(form, lhs=0, rhs=3))[1],
    as.character(formula(form, lhs=0, rhs=3))[2])), 
    dist="weibull", data=train_dat)

  alpha1 <- 1 / fit.survreg.1$scale
  alpha2 <- 1 / fit.survreg.2$scale
  alpha3 <- 1 / fit.survreg.3$scale

  startVals <- c(-alpha1*coef(fit.survreg.1)[1], log(alpha1),
                 -alpha2*coef(fit.survreg.2)[1], log(alpha2),
                 -alpha3*coef(fit.survreg.3)[1], log(alpha3))

  frailty <- TRUE

  if(frailty == TRUE) startVals <- c(startVals, 0.5)

  startVals <- c(startVals,

    -coef(fit.survreg.1)[-1] * alpha1,
    -coef(fit.survreg.2)[-1] * alpha2,
    -coef(fit.survreg.3)[-1] * alpha3)

  #-  Fit Markov Model on Training Folds
  
  fit <- nlm(logLike.weibull.SCR,

    p = startVals * runif(length(startVals), 0.9, 1.1),

    y1 = y1, delta1 = delta1, y2 = y2, delta2 = delta2,

    Xmat1 = as.matrix(Xmat1),
    Xmat2 = as.matrix(Xmat2),
    Xmat3 = as.matrix(Xmat3),

    frailty = frailty,

    iterlim = 1000,

    hessian = TRUE,

    steptol = 1e-5)
  
  myLabels <- c("log(kappa1)", "log(alpha1)",
                "log(kappa2)", "log(alpha2)",
                "log(kappa3)", "log(alpha3)")

  if(frailty == TRUE) myLabels <- c(myLabels, "log(theta)")

  myLabels <- c(myLabels, colnames(Xmat1), colnames(Xmat2), colnames(Xmat3))

  nP <- c(ncol(Xmat1), ncol(Xmat2), ncol(Xmat3))

  value <- list(estimate = fit$estimate,
                Finv     = solve(fit$hessian),
                logLike  = -fit$minimum,
                myLabels = myLabels,
                frailty  = frailty,
                nP       = nP,
                Xmat     = list(Xmat1, Xmat2, Xmat3))
  
  value$class  <- c("Freq_HReg", "ID", "Ind", "WB", "Markov")
  
  class(value) <- "Freq_HReg"
  
  theta <- summary(value)$theta[1]
  
  k1 <- exp(summary(value)$h0[1,1]);  a1 <- exp(summary(value)$h0[2,1])
  k2 <- exp(summary(value)$h0[1,4]);  a2 <- exp(summary(value)$h0[2,4])
  k3 <- exp(summary(value)$h0[1,7]);  a3 <- exp(summary(value)$h0[2,7])

  p1 <- length(summary(value)$coef[,1])
  p2 <- length(summary(value)$coef[,1])
  p3 <- length(summary(value)$coef[,1])

  beta1 <- summary(value)$coef[,1]
  beta2 <- summary(value)$coef[,4]
  beta3 <- summary(value)$coef[,7]

  # Linear predictors
  
  h1 <- as.numeric(X_test %*% beta1)
  h2 <- as.numeric(X_test %*% beta2)
  h3 <- as.numeric(X_test %*% beta3)

  # Compute survival probability predictions at each t_pred

  preds <- sapply(t_pred, function(t) {
    
    (1 + theta * k1 * t^a1 * exp(h1) + theta * k2 * t^a2 * exp(h2))^(-1/theta)
  })

  # Compute the bbs() on the test fold
  
  res        <- bbs1(test_dat, preds, t_pred)
  bs_mat[,k] <- res$bs
  ibs_vec[k] <- res$ibs
}

# Point-wise mean & 95% CI of BBS

mean_bs <- rowMeans(bs_mat)
sd_bs   <- apply(bs_mat, 1, sd)
se_bs   <- sd_bs / sqrt(K)
t_q     <- qt(0.975, df = K-1)

ci_lower <- mean_bs - t_q * se_bs
ci_upper <- mean_bs + t_q * se_bs

# Assemble for plotting

plot_df <- data.frame(
  time  = t_pred,
  mean  = mean_bs,
  lower = ci_lower,
  upper = ci_upper
)

# Example plot

plot(plot_df$time, plot_df$mean, 
     
  type = "l", lwd = 2, xlab = "Time (yrs)", ylab = "BBS")

polygon(c(t_pred, rev(t_pred)), c(ci_lower, rev(ci_upper)),
        
  col = rgb(0,0,1,0.2), border = NA)

lines(plot_df$time, plot_df$mean, lwd = 2, col = "blue")

# Integrated BBS & 95% CI

mean_ibs <- mean(ibs_vec)
se_ibs   <- sd(ibs_vec) / sqrt(K)
ci_ibs   <- mean_ibs + c(-1,1) * t_q * se_ibs

cat(sprintf("iBBS = %.4f (95%% CI: %.4f – %.4f)\n",
            
  mean_ibs, ci_ibs[1], ci_ibs[2]))

#--- METHOD 01 -----------------------------------------------------------------

K      <- 5
folds  <- createFolds(data$Di1, k = K)
t_pred <- seq(0, 5, length.out = 100)

# Storage

bs_mat   <- matrix(NA, nrow = length(t_pred), ncol = K)
ibs_vec  <- numeric(K)

myModel  <- c("Markov", "Weibull")

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

for (k in seq_along(folds)) {

  # Split train/test
  
  test_idx  <- folds[[k]]
  train_dat <- data[-test_idx, ]
  test_dat  <- data[ test_idx, ]
  X_train   <- X[-test_idx,-c(1,25)]
  X_test    <- X[test_idx,-c(1,25)]

  y1     <- as.vector(train_dat$Yi1)
  delta1 <- as.vector(train_dat$Di1)
  y2     <- as.vector(train_dat$Yi2)
  delta2 <- as.vector(train_dat$Di2)

  startValues <- initiate.startValues_HReg(

    form, train_dat, model = myModel, nChain = 2)

  fit <- BayesID_HReg(form, train_dat, id = NULL, model = myModel,

    hyperParams, startValues, mcmc.WB)

  theta <- summary(fit)$theta[1]
  
  k1 <- exp(summary(fit)$h0[1,1]);  a1 <- exp(summary(fit)$h0[2,1])
  k2 <- exp(summary(fit)$h0[1,4]);  a2 <- exp(summary(fit)$h0[2,4])
  k3 <- exp(summary(fit)$h0[1,7]);  a3 <- exp(summary(fit)$h0[2,7])

  p1 <- length(summary(fit)$coef[,1])
  p2 <- length(summary(fit)$coef[,1])
  p3 <- length(summary(fit)$coef[,1])

  beta1 <- summary(fit)$coef[,1]
  beta2 <- summary(fit)$coef[,4]
  beta3 <- summary(fit)$coef[,7]

  # Linear predictors
  
  h1 <- as.numeric(X_test %*% beta1)
  h2 <- as.numeric(X_test %*% beta2)
  h3 <- as.numeric(X_test %*% beta3)

  # Compute survival probability predictions at each t_pred

  preds <- sapply(t_pred, function(t) {
    
    (1 + theta * k1 * t^a1 * exp(h1) + theta * k2 * t^a2 * exp(h2))^(-1/theta)
  })

  # Compute the bbs() on the test fold
  
  res        <- bbs1(test_dat, preds, t_pred)
  bs_mat[,k] <- res$bs
  ibs_vec[k] <- res$ibs
}

# Point-wise mean & 95% CI of BBS

mean_bs <- rowMeans(bs_mat)
sd_bs   <- apply(bs_mat, 1, sd)
se_bs   <- sd_bs / sqrt(K)
t_q     <- qt(0.975, df = K-1)

ci_lower <- mean_bs - t_q * se_bs
ci_upper <- mean_bs + t_q * se_bs

# Assemble for plotting

plot_df <- data.frame(
  time  = t_pred,
  mean  = mean_bs,
  lower = ci_lower,
  upper = ci_upper
)

# Example plot

plot(plot_df$time, plot_df$mean, 
     
  type = "l", lwd = 2, ylim = range(0, 0.8), xlab = "Time (yrs)", ylab = "BBS")

polygon(c(t_pred, rev(t_pred)), c(ci_lower, rev(ci_upper)),
        
  col = rgb(0,0,1,0.2), border = NA)

lines(plot_df$time, plot_df$mean, lwd = 2, col = "blue")

plot_df <- read.csv("/R/method0.csv")
bbs     <- read.csv("./R/method0_bbs.csv")

mean2 <- bbs$mean

ci_lower2 <- bbs$ci_lower
ci_upper2 <- bbs$ci_upper

polygon(c(t_pred, rev(t_pred)), c(ci_lower2, rev(ci_upper2)),
        
  col = rgb(1, 0, 0, 0.2), border = NA)

lines(t_pred, mean2, lwd = 2, col = "red")

# Integrated BBS & 95% CI

mean_ibs <- mean(ibs_vec)
se_ibs   <- sd(ibs_vec) / sqrt(K)
ci_ibs   <- mean_ibs + c(-1,1) * t_q * se_ibs

cat(sprintf("iBBS = %.4f (95%% CI: %.4f – %.4f)\n",

  mean_ibs, ci_ibs[1], ci_ibs[2]))

#--- METHOD 02 -----------------------------------------------------------------

K      <- 5
folds  <- createFolds(data$Di1, k = K)
t_pred <- seq(0, 5, length.out = 100)

# Storage

bs_mat  <- matrix(NA, nrow = length(t_pred), ncol = K)
ibs_vec <- numeric(K)

myModel <- "LN"

theta.ab <- c(0.5, 0.05)

LN.ab1 <- c(0.3, 0.3)
LN.ab2 <- c(0.3, 0.3)
LN.ab3 <- c(0.3, 0.3)

hyperParams <- list(theta = theta.ab, 
                    
  LN = list(LN.ab1 = LN.ab1, LN.ab2 = LN.ab2, LN.ab3 = LN.ab3))

# MCMC Settings

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

for (k in seq_along(folds)) {

  # Split train/test
  
  test_idx  <- folds[[k]]
  train_dat <- data[-test_idx, ]
  test_dat  <- data[ test_idx, ]
  X_train   <- X[-test_idx,-1]
  X_test    <- X[test_idx,-1]

  y1     <- as.vector(train_dat$Yi1)
  delta1 <- as.vector(train_dat$Di1)
  y2     <- as.vector(train_dat$Yi2)
  delta2 <- as.vector(train_dat$Di2)

  train_dat$y1L <- train_dat$y1U <- y1

  train_dat$y1U[which(delta1 == 0)] <- Inf

  train_dat$y2L <- train_dat$y2U <- y2

  train_dat$y2U[which(delta2 == 0)] <- Inf

  train_dat$LT <- rep(0, dim(train_dat)[1])

  startValues <- initiate.startValues_AFT(

    form, train_dat, model = myModel, nChain = 2)

  fit <- BayesID_AFT(form, train_dat, model = myModel, hyperParams,

    startValues, mcmcParams)

  theta <- summary(fit)$theta[1]
  
  k1 <- exp(summary(fit)$h0[1,1]);  a1 <- exp(summary(fit)$h0[2,1])
  k2 <- exp(summary(fit)$h0[1,4]);  a2 <- exp(summary(fit)$h0[2,4])
  k3 <- exp(summary(fit)$h0[1,7]);  a3 <- exp(summary(fit)$h0[2,7])

  p1 <- length(summary(fit)$coef[,1])
  p2 <- length(summary(fit)$coef[,1])
  p3 <- length(summary(fit)$coef[,1])

  beta1 <- summary(fit)$coef[,1]
  beta2 <- summary(fit)$coef[,4]
  beta3 <- summary(fit)$coef[,7]

  # Linear predictors
  
  h1 <- as.numeric(X_test %*% beta1)
  h2 <- as.numeric(X_test %*% beta2)
  h3 <- as.numeric(X_test %*% beta3)

  # Compute survival probability predictions at each t_pred

  preds <- sapply(t_pred, function(t) {
    
    (1 + theta * k1 * t^a1 * exp(h1) + theta * k2 * t^a2 * exp(h2))^(-1/theta)
  })

  # Compute the bbs() on the test fold
  
  res        <- bbs1(test_dat, preds, t_pred)
  bs_mat[,k] <- res$bs
  ibs_vec[k] <- res$ibs
}

# Point-wise mean & 95% CI of BBS

mean_bs <- rowMeans(bs_mat)
sd_bs   <- apply(bs_mat, 1, sd)
se_bs   <- sd_bs / sqrt(K)
t_q     <- qt(0.975, df = K-1)

ci_lower <- mean_bs - t_q * se_bs
ci_upper <- mean_bs + t_q * se_bs

# Assemble for plotting

plot_df <- data.frame(
  time  = t_pred,
  mean  = mean_bs,
  lower = ci_lower,
  upper = ci_upper
)

# Example plot

plot(plot_df$time, plot_df$mean, 
     
  type = "l", lwd = 2, xlab = "Time (yrs)", ylab = "BBS")

polygon(c(t_pred, rev(t_pred)), c(ci_lower, rev(ci_upper)),
        
  col = rgb(0,0,1,0.2), border = NA)

lines(plot_df$time, plot_df$mean, lwd = 2, col = "blue")

# Integrated BBS & 95% CI

mean_ibs <- mean(ibs_vec)
se_ibs   <- sd(ibs_vec) / sqrt(K)
ci_ibs   <- mean_ibs + c(-1,1) * t_q * se_ibs

cat(sprintf("iBBS = %.4f (95%% CI: %.4f – %.4f)\n",
            
  mean_ibs, ci_ibs[1], ci_ibs[2]))

#--- METHOD 03 -----------------------------------------------------------------

#- Prepare folds and prediction grid

K      <- 5
folds  <- createFolds(data$Di1, k = K)
t_pred <- seq(0, 5, length.out = 100)

# Storage

bs_mat  <- matrix(NA, nrow = length(t_pred), ncol = K)
ibs_vec <- numeric(K)

for (k in seq_along(folds)) {
  
  # Split train/test
  
  test_idx  <- folds[[k]]
  train_dat <- data[-test_idx, ]
  test_dat  <- data[ test_idx, ]
  X_train   <- X[-test_idx,-c(1,25)]
  X_test    <- X[test_idx,-c(1,25)]
  
  y1     <- as.vector(train_dat$Yi1)
  delta1 <- as.vector(train_dat$Di1)
  y2     <- as.vector(train_dat$Yi2)
  delta2 <- as.vector(train_dat$Di2)
  
  V <- y1
  
  a <- train_dat$Ti1 == train_dat$Ci
  
  a[is.na(a)] <- FALSE
  
  train_dat$Ci <- ifelse(a, train_dat$Ci+1e-5, train_dat$Ci)
  
  b <- data$Ti1 == data$Ti2
  
  b[is.na(b)] <- FALSE
  
  train_dat$Ti2 <- ifelse(b, train_dat$Ti2+1e-5, train_dat$Ti2)
  
  delta1 <- train_dat$Yi1 == train_dat$Ti1
  
  delta1[is.na(delta1)] <- FALSE
  
  delta2 <- train_dat$Yi1 == train_dat$Ti2
  
  delta2[is.na(delta2)] <- FALSE
  
  W <- ifelse(!delta1, 0, train_dat$Yi2)
  
  delta3 <- (W == train_dat$Ti2)
  
  delta3[is.na(delta3)] <- FALSE
  
  #- Covariates
  
  Z12 <- Z13 <- Z23 <- X_train
  
  #-  Fit Markov Model on Training Folds
  
  fit <- estimate_all(z12 = Z12, z13 = Z13, z23 = Z23, 
                      
    obs_times1 = V, obs_times2 = W,
                         
    delta1 = delta1, delta2 = delta2, delta3 = delta3
  )
  
  beta1 <- log(summary(fit)$coef[,1])
  beta2 <- log(summary(fit)$coef[,4])
  beta3 <- log(summary(fit)$coef[,7])
  
  h1 <- as.matrix(dat_sim[, paste0("X", 1:12)]) %*% beta1
  h2 <- as.matrix(dat_sim[, paste0("X", 1:12)]) %*% beta2
  h3 <- as.matrix(dat_sim[, paste0("X", 1:12)]) %*% beta3
  
  theta <- c(summary(fit)$theta[1])
  
  mu1  <- c(summary(fit)$h0[1,1], n_sim)
  sig1 <- c(summary(fit)$h0[2,1], n_sim)
  
  mu2  <- c(summary(fit)$h0[1,4], n_sim)
  sig2 <- c(summary(fit)$h0[2,4], n_sim)
  
  mu3  <- c(summary(fit)$h0[1,7], n_sim)
  sig3 <- c(summary(fit)$h0[2,7], n_sim)
  
  # Compute survival probability predictions at each t_pred
  
  preds <- sapply(t_pred, function(t) {
    
    (1 + theta *mu1*t^sig1 * exp(h1) + theta*mu2*t^sig2 * exp(h2))^(-1/theta)
  })
  
  # Compute the bbs() on the test fold
  
  res        <- bbs1(test_dat, preds, t_pred)
  bs_mat[,k] <- res$bs
  ibs_vec[k] <- res$ibs
}

# Point-wise mean & 95% CI of BBS

mean_bs <- rowMeans(bs_mat)
sd_bs   <- apply(bs_mat, 1, sd)
se_bs   <- sd_bs / sqrt(K)
t_q     <- qt(0.975, df = K-1)

ci_lower <- mean_bs - t_q * se_bs
ci_upper <- mean_bs + t_q * se_bs

# Assemble for plotting

plot_df <- data.frame(
  time  = t_pred,
  mean  = mean_bs,
  lower = ci_lower,
  upper = ci_upper
)

# Example plot

plot(plot_df$time, plot_df$mean, 
     
     type = "l", lwd = 2, xlab = "Time (yrs)", ylab = "BBS")

polygon(c(t_pred, rev(t_pred)), c(ci_lower, rev(ci_upper)),
        
        col = rgb(0,0,1,0.2), border = NA)

lines(plot_df$time, plot_df$mean, lwd = 2, col = "blue")

# Integrated BBS & 95% CI

mean_ibs <- mean(ibs_vec)
se_ibs   <- sd(ibs_vec) / sqrt(K)
ci_ibs   <- mean_ibs + c(-1,1) * t_q * se_ibs

cat(sprintf("iBBS = %.4f (95%% CI: %.4f – %.4f)\n",
            
            mean_ibs, ci_ibs[1], ci_ibs[2]))

#--- METHOD 04 -----------------------------------------------------------------

#- Prepare folds and prediction grid

K      <- 5
folds  <- createFolds(data$Di1, k = K)
t_pred <- seq(0, 5, length.out = 100)

# Storage

bs_mat  <- matrix(NA, nrow = length(t_pred), ncol = K)
ibs_vec <- numeric(K)

for (k in seq_along(folds)) {
  
  # Split train/test
  
  test_idx  <- folds[[k]]
  train_dat <- data[-test_idx, ]
  test_dat  <- data[ test_idx, ]
  X_train   <- X[-test_idx,-c(1,25)]
  X_test    <- X[test_idx,-c(1,25)]
  
  V <- train_dat[, "Yi1"]
  W <- train_dat[, "Yi2"]
  
  delta1 <- train_dat$Yi1 == train_dat$Ti1
  delta2 <- train_dat$Yi1 == train_dat$Ti2
  delta3 <- train_dat$Yi2 == train_dat$Ti2
  
  #- Covariates
  
  Z12 <- Z13 <- Z23 <- X_train
  
  fit <- estimation_with_frailty(X01 = Z12, X02 = Z13, X12 = Z23, V = V,
                                    
    W = W, delta1 = delta1, delta2 = delta2, delta3 = delta3, B = 0,
                                    
    print = F)
  
  beta1 <- fit$ests[ 2:13, 1]
  beta2 <- fit$ests[14:25, 1]
  beta3 <- fit$ests[26:37, 1]
  
  h1 <- as.matrix(dat_sim[, paste0("X", 1:12)]) %*% beta1
  h2 <- as.matrix(dat_sim[, paste0("X", 1:12)]) %*% beta2
  h3 <- as.matrix(dat_sim[, paste0("X", 1:12)]) %*% beta3
  
  theta <- c(fit$ests[1, 1])
  
  lam01 <- fit$H0s[[1]]
  lam02 <- fit$H0s[[2]]
  lam03 <- fit$H0s[[3]]
  
  
  # Compute survival probability predictions at each t_pred
  
  preds <- sapply(t_pred, function(t) {
    
    (1 + theta * lam01 * exp(h1) + theta * lam02 * exp(h2))^(-1/theta)
  })
  
  # Compute the bbs() on the test fold
  
  res        <- bbs1(test_dat, preds, t_pred)
  bs_mat[,k] <- res$bs
  ibs_vec[k] <- res$ibs
}

# Point-wise mean & 95% CI of BBS

mean_bs <- rowMeans(bs_mat)
sd_bs   <- apply(bs_mat, 1, sd)
se_bs   <- sd_bs / sqrt(K)
t_q     <- qt(0.975, df = K-1)

ci_lower <- mean_bs - t_q * se_bs
ci_upper <- mean_bs + t_q * se_bs

# Assemble for plotting

plot_df <- data.frame(
  time  = t_pred,
  mean  = mean_bs,
  lower = ci_lower,
  upper = ci_upper
)

# Example plot

plot(plot_df$time, plot_df$mean, 
     
     type = "l", lwd = 2, xlab = "Time (yrs)", ylab = "BBS")

polygon(c(t_pred, rev(t_pred)), c(ci_lower, rev(ci_upper)),
        
        col = rgb(0,0,1,0.2), border = NA)

lines(plot_df$time, plot_df$mean, lwd = 2, col = "blue")

# Integrated BBS & 95% CI

mean_ibs <- mean(ibs_vec)
se_ibs   <- sd(ibs_vec) / sqrt(K)
ci_ibs   <- mean_ibs + c(-1,1) * t_q * se_ibs

cat(sprintf("iBBS = %.4f (95%% CI: %.4f – %.4f)\n",
            
            mean_ibs, ci_ibs[1], ci_ibs[2]))

#=== END =======================================================================