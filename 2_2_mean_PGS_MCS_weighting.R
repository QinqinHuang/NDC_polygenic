#------------------------------------------------------------
# 2024-09-13 last modified
#
# calculate reweighted average PGS in MCS unrelated children 
# and trio children
#------------------------------------------------------------
library(data.table)
library(dplyr)
library(foreach)
source("/lustre/scratch123/hgi/projects/ddd/users/qh1/pgs_paper/scripts/functions__modified.R")


#----- input files -----
## MCS weights
MCSweights = fread("qinqin_weights.csv")
MCSweights$FID = gsub(" ", "", MCSweights$FID)
MCSweights_children = MCSweights[analaysis == "ALL", .(FID, PRS, weight)]
MCSweights_triochildren = MCSweights[analaysis == "TRIO", .(FID, PRS, weight)]

# birth cohort master file
master_birthcohort = fread("/lustre/scratch123/hgi/projects/ddd/users/qh1/pgs_paper/data/alspac_mcs_master_file.txt")

# PGS
PGSarray = fread("/lustre/scratch123/hgi/projects/ddd/users/qh1/PGS_common_SNPs/PC_adjusted_PGS/LDpred_PGS_postQC_adjusted_20PCs_gcta_scaled_GELUKHLS.txt")


#---- prepare data -----
PGSarray = select(PGSarray, IID, FID, ends_with("GELUKHLS"))
names(PGSarray) = gsub("_GELUKHLS", "", names(PGSarray))

# add all MCS PGS 
# 6036 unrelated children
PGS_MCSchildren = PGSarray[IID %in% master_birthcohort[cohort == "MCS" & unrelated_popctrl == 1]$IID]
length(unique(PGS_MCSchildren$FID)) & all(MCSweights_children$FID %in% PGS_MCSchildren$FID)
MCSweights_children = merge(MCSweights_children, PGS_MCSchildren, by = "FID")

# 2498 trio children
PGS_MCStriochildren = PGSarray[IID %in% master_birthcohort[cohort == "MCS" & trio_unrelated == 1]$IID]
length(unique(PGS_MCStriochildren$FID)) & all(MCSweights_triochildren$FID %in% PGS_MCStriochildren$FID)
MCSweights_triochildren = merge(MCSweights_triochildren, PGS_MCStriochildren, by = "FID")


#----- mean PGS after weighting -----
meanweighted = foreach(pp = c("EA", "CP", "EANonCog", "SCZ", "NDD"), .combine = rbind) %do% {
  MCSweights_children$PGS = MCSweights_children[, get(pp)]
  MCSweights_triochildren$PGS = MCSweights_triochildren[, get(pp)]
  
  myfit1 = summary(lm(PGS ~ 1, weights = weight, data = MCSweights_children))$coefficients
  myfit2 = summary(lm(PGS ~ 1, weights = weight, data = MCSweights_triochildren))$coefficients
  
  return(data.table(
    subgroup = c("MCS_children_weighted", "MCS_triochildren_weighted"),
    N = c(nrow(MCSweights_children), nrow(MCSweights_triochildren)),
    PGS = pp,
    mean = c(myfit1[1], myfit2[1]),
    se = c(myfit1[2], myfit2[2])
  ))
}

meanweighted[, sd := se * sqrt(N)]
meanweighted[, CI_L := mean - 1.96*se]
meanweighted[, CI_U := mean + 1.96*se]



# --> 
write.csv(meanweighted, file = "mean_sd_PGS_MCS_afterweighting.csv", row.names = F)






