#===============================================================================
#
#  PROGRAM: BLCS_TREE.R
#
#  AUTHORS: Stephen Salerno, Zhilin Zhang
#
#  PURPOSE: Real-data analysis of the BLSC data with the tree-based method
#
#  INPUT:   BLCS_Data_Updated_February2023.xlxs
#
#  NOTES: 
#
#  We cannot provide the real data here, however we provide the codes here for
#  review and reproducibility with similarly structured data. This script
#  implements the real data analysis for the tree-based method.
#
#===============================================================================

#=== INITIALIZATION ============================================================

#--- LOAD NECESSARY PACKAGES ---------------------------------------------------

library(pacman)

p_load(gtsummary, survival, survcomp, survminer, ggsci, patchwork,
      
  readxl, lubridate, here, tidyverse, SemiCompRisks, dplyr, update = F)

library(dplyr, warn.conflicts = FALSE)

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
  
  dplyr::select(ptime, pevent, dtime, devent,
                
    Yi1, Yi2, YiS, Di1, Di2, DiS, 
    
    agedx, sex, race, ethnic, education, bmi, smk, pkyrs, trt, 
    
    ctype, cstage, copd, asthma, egfr, kras, everything())

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
  
  dplyr::select(Yi1, Yi2, YiS, Di1, Di2, DiS, Ti1, Ti2, Ci, stage_splt, cstage, trt,
         
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
    
    asthma = asthma |> fct_na_value_to_level("Unknown")
  ) |>
  
  dplyr::select(Yi1, Yi2, Di1, Di2, stage_splt, ctype, trt,
         
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

library(randomForestSRC)

library(caret)

rhs <- ~ agedx + sexFemale + sexUnknown + raceOther + raceAsian + 
  `raceBlack/African American` + raceUnknown + ethnicHispanic + ethnicUnknown + 
  copdYes + copdUnknown + asthmaYes + asthmaUnknown + stage_splt1 + trtOther + 
  trtRadiation + trtSurgery + trtUnknown + `smkFormer smoker` +
  `smkCurrent smoker` + `smkSmoker, status unknown` + pkyrs + egfrYes + 
  `egfrNot Tested` + krasYes

K      <- 5

folds  <- createFolds(data$Di1, k = K)

t_pred <- seq(0, 5, length.out = 100)

bs_mat   <- matrix(NA, nrow = length(t_pred), ncol = K)
ibs_vec  <- numeric(K)

cind <- numeric(K)

for (k in seq_along(folds)) {
  
  test_idx  <- folds[[k]]
  train_dat <- data[-test_idx, ]
  test_dat  <- data[ test_idx, ]
  X_train <- X[-test_idx,-c(1,27)]
  X_test <- X[test_idx,-c(1,27)]
  
  form1 <- update(rhs, Surv(Yi1, Di1) ~ .)
  rsf_h1 <- rfsrc(
    formula = form1,
    data = train_dat,
    ntree = 500,
    splitrule = "logrank"
  )
  
  form2 <- update(rhs, Surv(Yi2, Di2) ~ .)
  rsf_h2 <- rfsrc(
    formula = form2,
    data = train_dat,
    ntree = 500,
    splitrule = "logrank"
  )
  
  pr1 <- predict(rsf_h1, newdata = test_dat) 
  pr2 <- predict(rsf_h2, newdata = test_dat)  
  
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
      cbind(Hg[,1,drop=FALSE], 
            Hg[, -1, drop=FALSE] 
            - Hg[, -ncol(Hg), drop=FALSE])
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
  
  preds <- cr_surv_from_trees(pr1, pr2, t_pred)$S
  
  res        <- bbs1(test_dat, preds, t_pred)
  bs_mat[,k] <- res$bs
  ibs_vec[k] <- res$ibs
  
  test_dat$event_all <- ifelse(test_dat$Di1 == 1 | test_dat$Di2 == 1, 1, 0)
  test_dat$time_all  <- pmin(test_dat$Yi1, test_dat$Yi2, na.rm = TRUE)
  
  risk_score <- rowMeans(1 - preds)
  
  cind[k] <- survConcordance(
    
    Surv(time_all, event_all) ~ risk_score, data = test_dat)
}

t_q <- qt(0.975, df = K-1)

mean(unlist(cind))
mean(unlist(cind)) - t_q * sd(unlist(cind))
mean(unlist(cind)) + t_q * sd(unlist(cind))

mean(ibs_vec)
mean(ibs_vec) - t_q * sd(ibs_vec)
mean(ibs_vec) + t_q * sd(ibs_vec)

#=== END =======================================================================
