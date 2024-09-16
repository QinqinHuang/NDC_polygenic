#----------------------------------------------------------------------
# 2024-09-13 last modified
#
# Test associations between:
# Diagnostic status and factors, with or without correcting for PGS
# Diagnostic status and PGS, with or without correcting for factors
# PGS and factors affecting diagnosis
#----------------------------------------------------------------------
library(data.table)
library(dplyr)
library(foreach)


## read input data
# a table that contains PGS and phenotypes (diagnostic status and factors associated with it)
phenoPGS = fread("phenotype_table_NDDGBR_probands_2023_07_28.txt")

# recode some phenotypes
phenoPGS$id.dd.ordered = factor(phenoPGS$id.dd.ordered, levels = c("no", "mild", "moderate", "unknown", "severe"))
phenoPGS$only_patient_affected = factor(phenoPGS$only_patient_affected, levels = c("No", "Unknown", "Yes"))
# make a new variable to compare severe ID against mild and moderate
phenoPGS$id.dd_severe_vs_m = ""
phenoPGS[id.dd == "severe"]$id.dd_severe_vs_m = "severe"
phenoPGS[id.dd %in% c("mild", "moderate")]$id.dd_severe_vs_m = "mild or moderate"
phenoPGS[id.dd_severe_vs_m == ""]$id.dd_severe_vs_m = NA
phenoPGS$id.dd_severe_vs_m = factor(phenoPGS$id.dd_severe_vs_m, levels = c("mild or moderate", "severe"))
# any affected family members
phenoPGS$affected_any = ifelse(phenoPGS$only_patient_affected == "No", yes = 1, no = 0)
# treat unknown as NA
phenoPGS[only_patient_affected == "Unknown", affected_any := NA]


# N=7549 unrelated probands regardless of trio status
# Note we didn't remove Scottish samples
NDDGBR = phenoPGS

# trios N=2497
NDDGBR_trio = NDDGBR[!is.na(mum_EA) & !is.na(dad_EA) & mother_rel_to_remove == 0 & father_rel_to_remove == 0]


#----- (0) function to get summary stats from fitted regression models -----
getcoef = function(mymodel) {
  mycoeff = summary(mymodel)$coefficients
  returndd = as.data.table(mycoeff)
  returndd = returndd[-1, -3]
  names(returndd) = c("estimate", "SE", "p")
  returnddS[, OR := signif(exp(estimate), 4)]
  returnddS[, OR_L := signif(exp(estimate - 1.96*SE), 4)]
  returnddS[, OR_U := signif(exp(estimate + 1.96*SE), 4)]
  returndd$variable = rownames(mycoeff)[-1]
  returndd$N = length(mymodel$residuals)
  # if the predictor is a binary variable
  if(length(unique(mymodel$model[,1])) == 2) {
    returndd$N0 = table(mymodel$model[,1])[1]
    returndd$N1 = table(mymodel$model[,1])[2]
  } else {returndd$N1 = returndd$N0 = NA}
  return(returndd)
}



#----- (1) diagnostic status ~ x and/or PGS -----
# compare three models: diagnosed ~ PGS, diagnosed ~ x, and diagnosed ~ x + PGS
# list of factors of interest
factors_investigated = c("trio", "proband_sex", "affected_any", "prematurity", "id.dd_severe_vs_m", "FROH_0.0625", "maternal_diabetes")

## controlling for proband's PGS
diagnosed_vs_childPGS = foreach(pp = c("EA", "CP", "EANonCog", "SCZ", "NDD"), .combine = rbind) %do% {
  # current PGS
  mydd = NDDGBR
  mydd$PGS = mydd[, get(paste0("proband_",pp))]
  mydd = mydd[!is.na(PGS)]
  
  # (1) diagnosed ~ PGS
  PGSmodel = getcoef(glm(formula = "diagnosed ~ PGS", 
                         data = mydd, family = binomial))
  PGSmodel$model = "diagnosed ~ PGS"
  
  # (2) diagnosed ~ x and (3) diagnosed ~ PGS + x
  summ_diagnosed = foreach(vv = factors_investigated, .combine = rbind) %do% {
    model1 = getcoef(glm(formula = as.formula(paste0("diagnosed ~ ", vv)), 
                         data = mydd, family = binomial))
    model1$model = paste0("diagnosed ~ ", vv)
    model2 = getcoef(glm(formula = as.formula(paste0("diagnosed ~ PGS + ", vv)), 
                         data = mydd, family = binomial))
    model2$model = paste0("diagnosed ~ PGS + ", vv)
    return(rbind(model1, model2))
  }
  
  # (4) fitting multiple factors
  model_full = getcoef(glm(formula = diagnosed ~ PGS + trio + prematurity, data = mydd, family = binomial))
  model_full$model = "diagnosed ~ PGS + trio + prematurity"
  
  d = rbind(PGSmodel, summ_diagnosed, model_full)
  d$PGS = pp
  return(d)
}



#----- (2) PGS ~ factor -----
##### proband_PGS ~ x #####
probandPGS = foreach(pp = c("EA", "CP", "EANonCog", "SCZ", "NDD"), .combine = rbind) %do% {
  # current PGS
  mydd = NDDGBR
  mydd$PGS = mydd[, get(paste0("proband_",pp))]
  mydd = mydd[!is.na(PGS)]
  
  summ_PGS = foreach(vv = factors_investigated, .combine = rbind) %do% {
    # linear regression
    returndd = getcoef(lm(formula = paste0("PGS ~ ", vv), data = mydd[!is.na(mydd[, get(vv)])]))
    return(returndd)
  }
  summ_PGS$PGS = pp
  
  return(summ_PGS)
}

##### mother_PGS ~ x #####
# remove mothers due to relatedness
mumPGS = foreach(pp = c("EA", "CP", "EANonCog", "SCZ", "NDD"), .combine = rbind) %do% {
  # current PGS
  mydd = NDDGBR_trio
  mydd$PGS = mydd[, get(paste0("mum_",pp))]
  mydd = mydd[!is.na(PGS)]
  
  summ_PGS = foreach(vv = factors_investigated, .combine = rbind) %do% {
    # linear regression
    returndd = getcoef(lm(formula = paste0("PGS ~ ", vv), data = mydd[!is.na(mydd[, get(vv)])]))
    return(returndd)
  }
  summ_PGS$PGS = pp
  
  return(summ_PGS)
}

##### father_PGS ~ x ####
# remove fathers due to relatedness
dadPGS = foreach(pp = c("EA", "CP", "EANonCog", "SCZ", "NDD"), .combine = rbind) %do% {
  # current PGS
  mydd = NDDGBR_trio
  mydd$PGS = mydd[, get(paste0("dad_",pp))]
  mydd = mydd[!is.na(PGS)]
  
  summ_PGS = foreach(vv = factors_investigated, .combine = rbind) %do% {
    # linear regression
    returndd = getcoef(lm(formula = paste0("PGS ~ ", vv), data = mydd[!is.na(mydd[, get(vv)])]))
    return(returndd)
  }
  summ_PGS$PGS = pp
  
  return(summ_PGS)
}

probandPGS$model = "proband PGS ~ x"
mumPGS$model = "mother PGS ~ x"
dadPGS$model = "father PGS ~ x"

results = rbind(probandPGS, mumPGS, dadPGS)

