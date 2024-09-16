#----------------------------------------------------------------------
# 2024-09-13 last modified
#
# Aim: to calculate weighted and unweighted correlations in MCS between
# (1) PC adjusted RVBS in one parent and PGS in the other parent
# (2) PC adjusted RVBS and PGS in the same parent
# (3) PC adjusted RVBS in the child and child's PGS
#
# in these trios:
# NDD trios where parents are unaffected
# NDD undiagnosed trios where parents are unaffected
# NDD trios with a de novo diagnosis and unaffected parents
#
# for parents, using the same weight as their child
#----------------------------------------------------------------------
library(data.table)
library(foreach)
library(dplyr)


## input data
# MCS weights
MCSweights = fread("/lustre/scratch123/hgi/projects/ddd/users/qh1/pgs_paper/analysis_MCS_weighting/qinqin_weights.csv")
MCSweights$FID = gsub(" ", "", MCSweights$FID)
MCSweights_children = MCSweights[analaysis == "ALL", .(FID, PRS, weight)]
MCSweights_triochildren = MCSweights[analaysis == "TRIO", .(FID, PRS, weight)]

# RVBS in MCS trios
RVBS_triosamples = fread("/lustre/scratch123/hgi/projects/ddd/users/qh1/pgs_paper/data/rare_variants/birth_cohort_MCS/RVBS_adj20PCs_MCS_trios.txt")

# add weights to my RVBS tables
# 2246 out of 2297 trios have a valid weight
table(unique(RVBS_triosamples$FID) %in% MCSweights_triochildren$FID)
RVBS_triosamples_wt = merge(RVBS_triosamples, MCSweights_triochildren, by = "FID")


#----- function to calculate correlation in children -----
## a function to calculate correlation within the same person
cor_wrapper_sameperson = function(df, myvar, myPGS, myflag = NA) {
  # normalise column name
  df$RVBS = df[, get(myvar)]
  df$PGS = df[, get(myPGS)]
  df = df[!is.na(PGS)]
  
  ## unweighted correlation using lm
  t1 = summary(lm(scale(RVBS) ~ scale(PGS), data = df))$coefficients
  d1 = data.table(N = nrow(df), subgroup = myflag, weighted = F, 
                  method = "lm; scaled",
                  var = myvar, PGS = myPGS, 
                  cor = t1["scale(PGS)", "Estimate"], 
                  stat = t1["scale(PGS)", 3], p = t1["scale(PGS)", 4])
  
  ## weighted correlation using lm
  t2 = summary(lm(scale(RVBS) ~ scale(PGS), data = df, weights = weight))$coefficients
  d2 = data.table(N = nrow(df), subgroup = myflag, weighted = T, 
                  method = "lm; scaled",
                  var = myvar, PGS = myPGS, 
                  cor = t2["scale(PGS)", "Estimate"], 
                  stat = t2["scale(PGS)", 3], p = t2["scale(PGS)", 4])
  
  return(rbind(d1, d2))
}


## a function to calculate correlation cross parents
cor_wrapper = function(df, myvar, myPGS, myflag = NA) {
  # normalise column name
  df$RVBS = df[, get(myvar)]
  df$PGS = df[, get(myPGS)]
  df = df[!is.na(PGS)]
  
  # mother 
  df_mother = df[persontype == "mother", .(FID, RVBS, PGS, weight)]
  # father 
  df_father = df[persontype == "father", .(FID, RVBS, PGS, weight)]
  df_father = df_father[match(df_mother$FID, FID)]
  
  df_mother$partnerPGS = df_father$PGS
  df_father$partnerPGS = df_mother$PGS
  
  df = rbind(df_mother, df_father)
  
  ## unweighted correlation using lm
  t1 = summary(lm(scale(RVBS) ~ scale(partnerPGS), data = df))$coefficients
  d1 = data.table(N = nrow(df), subgroup = myflag, weighted = F, 
                    method = "lm; scaled",
                    var = myvar, PGS = myPGS, 
                    cor = t1["scale(partnerPGS)", "Estimate"], 
                    stat = t1["scale(partnerPGS)", 3], p = t1["scale(partnerPGS)", 4])
  
  ## weighted correlation
  t2 = summary(lm(scale(RVBS) ~ scale(partnerPGS), data = df, weights = weight))$coefficients
  d2 = data.table(N = nrow(df), subgroup = myflag, weighted = T, 
                  method = "lm; scaled",
                  var = myvar, PGS = myPGS, cor = t2["scale(partnerPGS)", "Estimate"], stat = t2["scale(partnerPGS)", 3], p = t2["scale(partnerPGS)", 4])
  
  return(rbind(d1, d2))
}




#----- calculate weighted and unweighted correlations -----
## correlations in MCS children 
corPGS = foreach(ii = c("EA", "EANonCog", "CP", "SCZ", "NDD"), .combine = rbind) %do% {
  foreach(jj = c("PTV_ddg2p_mono_lof_adj20PCs", "both_ddg2p_mono_lof_adj20PCs", "synonymous_ddg2p_mono_lof_adj20PCs", "PTV_pLI_adj20PCs", "both_pLI_adj20PCs", "synonymous_pLI_adj20PCs"), .combine = rbind) %do% {
    rbind(
      # within trio children
      cor_wrapper_sameperson(RVBS_triosamples_wt[persontype == "child"], myvar = jj, myPGS = ii, myflag = "MCS trio children"),
      
      # within trio parents
      cor_wrapper_sameperson(RVBS_triosamples_wt[persontype != "child"], myvar = jj, myPGS = ii, myflag = "MCS trio parents"),
      
      # cross parental correlation
      cor_wrapper(RVBS_triosamples_wt[persontype != "child"], myvar = jj, myPGS = ii, myflag = "MCS cross-parental")
    )
  }
} 

corPGS$sig = ifelse(corPGS$p < 0.05, yes = "*", no = "")
corPGS[p < 0.05/(5*2*3), sig := "**"]

