#------------------------------------------------------------
# 2024-09-13 last modified
#
# (1) regress out PCs and normalisation
# (2) calculate average PGS in various subgroups
# (3) two-sample t-tests comparing subgroups
#------------------------------------------------------------
library(data.table)
library(dplyr)
library(foreach)
source("/lustre/scratch123/hgi/projects/ddd/users/qh1/pgs_paper/scripts/functions__modified.R")


#----- (1) adjusting for PCs and standardising using controls -----
mymodel_20PCs = as.formula(paste0("PGS ~ ", paste0("PC", 1:20, collapse = " + ")))

# read raw PGS in GEL
rawPGSGEL = fread("/re_gecip/paediatrics/QQHuang/GEL_heritability/PGS_common_SNPs/LDpred_PGS_GEL.txt")

# linear regression coefficients trained using all array samples outside of GEL RE
lmweights = fread("/re_gecip/RR397/DDD_genotype/qh1/PGS/LDpred_PGS_postQC_lm_20PCs_gcta_coefficients.txt")

# read PCs - projected on to PC space calculated using unrelated array samples
GEL_gcta = fread("/re_gecip/RR397/DDD_genotype/qh1/merged_array_joint_PCA/gcta_pca_GEL_projected.proj.eigenvec.txt")

## adjust for 20 PCs
PGSPCs = merge(rawPGSGEL, select(GEL_gcta, -FID), by = "IID")
nothing = foreach(whichtrait = c("EA", "EANonCog", "EAsib", "CP", "SCZ", "NDD"), .combine = cbind) %do% {
  # raw PGS
  PGSPCs$PGS = PGSPCs[, get(paste0(whichtrait, "_raw"))]
  
  # calculate the residuals
  myweights = lmweights[-1, get(whichtrait)]
  PGSPCs[, adjPGS := PGS - (lmweights[, get(whichtrait)][1] + PC1 * myweights[1] + PC2 * myweights[2] + PC3 * myweights[3] + PC4 * myweights[4] + PC5 * myweights[5] + PC6 * myweights[6] + PC7 * myweights[7] + PC8 * myweights[8] + PC9 * myweights[9] + PC10 * myweights[10] + PC11 * myweights[11] + PC12 * myweights[12] + PC13 * myweights[13] + PC14 * myweights[14] + PC15 * myweights[15] + PC16 * myweights[16] + PC17 * myweights[17] + PC18 * myweights[18] + PC19 * myweights[19] + PC20 * myweights[20])]
  
  # update variable names
  names(PGSPCs)[ncol(PGSPCs)] = paste0(whichtrait, "_resid")
  PGSPCs$PGS = NULL
  
  return(NULL)
}


## PGS normalisation
# controls from GEL 
GELctrl = PGSPCs[ctrl==1]
GELctrl = select(GELctrl, ends_with("_resid"))
# controls from UKHLS
UKHLS = PGSarray_gcta[cohort == "UKHLS"]
UKHLS = select(UKHLS, ends_with("_resid"))
# concatenate controls
PGSctrl = rbind(UKHLS, GELctrl)

# calculate mean and sd
mean_EA = mean(PGSctrl$EA_resid)
mean_EAsib = mean(PGSctrl$EAsib_resid)
mean_EANonCog = mean(PGSctrl$EANonCog_resid)
mean_CP = mean(PGSctrl$CP_resid)
mean_SCZ = mean(PGSctrl$SCZ_resid)
mean_NDD = mean(GELctrl$NDD_resid)

sd_EA = sd(PGSctrl$EA_resid)
sd_EAsib = sd(PGSctrl$EAsib_resid)
sd_EANonCog = sd(PGSctrl$EANonCog_resid)
sd_CP = sd(PGSctrl$CP_resid)
sd_SCZ = sd(PGSctrl$SCZ_resid)
sd_NDD = sd(GELctrl$NDD_resid)

# Standardise PGS so that UKHLS and GEL controls have mean 0 and SD 1
# for the NDD PGS, UKHLS was excluded
PGSPCs$EA_GELUKHLS = (PGSPCs$EA_resid - mean_EA)/sd_EA
PGSPCs$EAsib_GELUKHLS = (PGSPCs$EAsib_resid - mean_EAsib)/sd_EAsib
PGSPCs$EANonCog_GELUKHLS = (PGSPCs$EANonCog_resid - mean_EANonCog)/sd_EANonCog
PGSPCs$CP_GELUKHLS = (PGSPCs$CP_resid - mean_CP)/sd_CP
PGSPCs$SCZ_GELUKHLS = (PGSPCs$SCZ_resid - mean_SCZ)/sd_SCZ
PGSPCs$NDD_GELUKHLS = (PGSPCs$NDD_resid - mean_NDD)/sd_NDD



#----- (2) calculate average PGS in various subgroups -----
## Function to calculate 95% CIs for normally distributed data
ci_95_normal <- function(group) {
  
  mean <- mean(group)
  stdev <- sd(group)
  sterr <- stdev/(sqrt(length(group)))
  
  upper_CI <- mean + 1.96*(sterr)
  lower_CI <- mean - 1.96*(sterr)
  
  CIs <- c(lower_CI, upper_CI)
  
  return(CIs)
}

## my function to return mean PGS and 95 CI in a given table
list_PGS = paste0(c("EA", "CP", "EANonCog", "SCZ", "NDD"), "_GELUKHLS")
mywrapper = function(mydata, cohortlabel, name_subset) {
  return(
    foreach(pgscname = list_PGS, .combine = rbind) %do% {
      # PGS
      myPGS = mydata[, get(pgscname)]
      myPGS = myPGS[!is.na(myPGS)]
      return(
        data.table(cohort = cohortlabel,
                   subgroup = name_subset,
                   PGS = pgscname,
                   N = length(myPGS),
                   mean = mean(myPGS),
                   sd = sd(myPGS),
                   CI_L = ci_95_normal(myPGS)[1],
                   CI_U = ci_95_normal(myPGS)[2])
      )
    }
  )
}

## load three PGS tables for controls, GEL samples, and DDD samples
#PGSctrl
#PGSGEL
#PGSDDD

# a. PGS mean and CI of controls 
summstats_ctrl = mywrapper(PGSctrl, cohortlabel = "control", name_subset = "control")

# b. PGS mean and CI of DDD/GEL subsets
# list of subsets: column names; each subgroup of interest is indicated by 1/0 in one column 
list_subsets = grep("NDD_trios|parents", names(PGSDDD), value = T)
summstats = foreach(mysubset = list_subsets, .combine = rbind) %do% {
  # DDD
  df_ddd = select(PGSDDD[which(PGSDDD[, get(mysubset)] == 1)], ends_with("_GELUKHLS"))
  # GEL
  df_gel = select(PGSGEL[which(PGSGEL[, get(mysubset)] == 1)], ends_with("_GELUKHLS"))
  return(
    rbind(
      mywrapper(df_ddd, cohortlabel = "DDD", name_subset = mysubset),
      # GEL
      mywrapper(df_gel, cohortlabel = "GEL", name_subset = mysubset),
      # DDD+GEL combined
      mywrapper(rbind(df_ddd, df_gel), cohortlabel = "DDD+GEL", name_subset = mysubset)
    )
  )
}



#----- (3) two-sample t-tests comparing subgroups ------
# read mean and sd of PGS in all subgroups
meanPGS = fread("all_mean_sd_PGS.csv")

## split the table into NDD probands, parents and controls
mean_NDD = meanPGS[type == "NDD"]
mean_parent = meanPGS[type == "parent"]
mean_ctrl = meanPGS[type == "control"]
mean_ctrl = mean_ctrl[N!=0] #remove the NDD PGS in UKHLS


## my function to generate t-test summary stats ####
# c1/c2: NDD subgroups are from either DDD, GEL or DDD+GEL; keep it NA for control subgroups
wrapper_test_meansd = function(meantable, g1, g2, c1=NA, c2=NA, pgs) {
  d1 = meantable[subgroup == g1 & PGS == pgs]
  d2 = meantable[subgroup == g2 & PGS == pgs]
  if(!is.na(c1)) { d1 = d1[cohort == c1] }
  if(!is.na(c2)) { d2 = d2[cohort == c1] }
  if(nrow(d1) == 0) {
    cat("no group:", g1, "with PGS", pgs, "\n")
    return(NULL)
  }
  if(nrow(d2) == 0) {
    cat("no group:", g2, "with PGS", pgs, "\n")
    return(NULL)
  }
  if(nrow(d1) >1 | nrow(d2) >1) {
    cat("Error -", g1, g2, pgs); return(NULL)
  }
  dd = data.table(PGS = pgs, cohort = c1, group1 = g1, group2 = g2,
                  N1 = d1$N, N2 = d2$N, mean1 = d1$mean, mean2 = d2$mean)
  dd = cbind(dd, t_test_summ(d1$mean, d2$mean, d1$sd, d2$sd, d1$N, d2$N))
  return(dd)
}

# e.g. to run the function to compare undiagnosed and diagnosed probands:
t_results = wrapper_test_meansd(mean_ctrl, g1="undiagnosed", g2="diagnosed", pgs="EA")

# Bonferroni correction
t_results$Significance = ""
t_results[p_ttest < 0.05]$Significance = "*"
t_results[p_ttest < 0.05/5]$Significance = "**"





