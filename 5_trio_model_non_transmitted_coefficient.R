#------------------------------------------------------------
# 2024-09-13 last modified
#
# regression models
#   NDD status ~ child PRS
#   NDD status ~ child PRS + dad PRS + mum PRS 
#
#// Main analysis - the mega analysis
# GEL+DDD NDD trios vs GEL+MCS+ASLAPC trios
#// Sensitivity analysis
# (1) GEL only analysis
# (2) with each control separately
#   - GEL+DDD vs GEL control trios
#   - GEL+DDD vs MCS control trios
#   - GEL+DDD vs ALSPAC control trios
#------------------------------------------------------------
library(data.table)
library(dplyr)
library(foreach)


## my wrapper to get summary stats from fitted models
get_stats = function(myfit) {
  fit_stats = summary(myfit)$coefficients
  predictorvar = rownames(fit_stats)
  fit_stats = as.data.table(fit_stats)[-1,]
  names(fit_stats) = c("beta", "se", "zvalue", "p")
  fit_stats$variable = predictorvar[-1]
  fit_stats$N = nrow(myfit$data)
  fit_stats$Ncases = sum(myfit$data$ndd)
  return(fit_stats)
}

## my function to run all regression models
run_glm = function(mydatatable, whichPGS, scalectrl = "GELUKHLS", whichanalysis) {
  
  # remove DDD CE samples in the GWAS when testing the NDD PGS
  if(whichPGS == "NDD") {
    mydatatable = mydatatable[cohort != "DDD_GSA"]
  }
  
  if(nrow(mydatatable[ndd == 0]) < 2 & nrow(mydatatable[ndd == 1]) < 2) {
    return(NULL)
  }
  
  # formulas
  formula_child = "ndd ~ child_PGS"
  formula_trio = "ndd ~ child_PGS + mother_PGS + father_PGS"
  
  ## normalise predictor names
  mydatatable$child_PGS = mydatatable[, get(paste0("child_", whichPGS, "_", scalectrl))]
  mydatatable$mother_PGS = mydatatable[, get(paste0("mother_", whichPGS, "_", scalectrl))]
  mydatatable$father_PGS = mydatatable[, get(paste0("father_", whichPGS, "_", scalectrl))]
  

  ## proband only model: NDD ~ child PGS
  to_return = get_stats(glm(formula_child, data = mydatatable, family = binomial))
  to_return = cbind(data.table(PGS = whichPGS, model = "child", cohort_cov = 0),
                    to_return)
  
  ## trio model: NDD ~ child PGS + father PGS + mother PGS
  to_return2 = get_stats(glm(formula_trio, data = mydatatable, family = binomial))
  to_return2 = cbind(data.table(PGS = whichPGS, model = "child+mum+dad", cohort_cov = 0),
                     to_return2)
  to_return = rbind(to_return, to_return2)
  
  to_return$subgroup = whichanalysis
  
  return(to_return)
}



## run regressions
# df is a table for undiagnosed trios with unaffected parents; each row is a trio family
summstats_all = rbind(
  # EA
  run_glm(mydatatable = df, whichPGS = "EA", whichanalysis = "mega"),
  run_glm(mydatatable = df[cohort == "GEL"], whichPGS = "EA", whichanalysis = "GEL only"),
  run_glm(mydatatable = df[ndd == 1 | (ndd == 0 & cohort == "GEL")], whichPGS = "EA", whichanalysis = "NDD vs GEL control"),
  run_glm(mydatatable = df[ndd == 1 | cohort == "MCS"], whichPGS = "EA", whichanalysis = "NDD vs MCS"),
  run_glm(mydatatable = df[ndd == 1 | cohort == "ALSPAC"], whichPGS = "EA", whichanalysis = "NDD vs ALSPAC"),
  
  # EA-NonCog
  run_glm(mydatatable = df, whichPGS = "EANonCog", whichanalysis = "mega"),
  run_glm(mydatatable = df[cohort == "GEL"], whichPGS = "EANonCog", whichanalysis = "GEL only"),
  run_glm(mydatatable = df[ndd == 1 | (ndd == 0 & cohort == "GEL")], whichPGS = "EANonCog", whichanalysis = "NDD vs GEL control"),
  run_glm(mydatatable = df[ndd == 1 | cohort == "MCS"], whichPGS = "EANonCog", whichanalysis = "NDD vs MCS"),
  run_glm(mydatatable = df[ndd == 1 | cohort == "ALSPAC"], whichPGS = "EANonCog", whichanalysis = "NDD vs ALSPAC"),
  
  # CP
  run_glm(mydatatable = df, whichPGS = "CP", whichanalysis = "mega"),
  run_glm(mydatatable = df[cohort == "GEL"], whichPGS = "CP", whichanalysis = "GEL only"),
  run_glm(mydatatable = df[ndd == 1 | (ndd == 0 & cohort == "GEL")], whichPGS = "CP", whichanalysis = "NDD vs GEL control"),
  run_glm(mydatatable = df[ndd == 1 | cohort == "MCS"], whichPGS = "CP", whichanalysis = "NDD vs MCS"),
  run_glm(mydatatable = df[ndd == 1 | cohort == "ALSPAC"], whichPGS = "CP", whichanalysis = "NDD vs ALSPAC"),
  
  # SCZ
  run_glm(mydatatable = df, whichPGS = "SCZ", whichanalysis = "mega"),
  run_glm(mydatatable = df[cohort == "GEL"], whichPGS = "SCZ", whichanalysis = "GEL only"),
  run_glm(mydatatable = df[ndd == 1 | (ndd == 0 & cohort == "GEL")], whichPGS = "SCZ", whichanalysis = "NDD vs GEL control"),
  run_glm(mydatatable = df[ndd == 1 | cohort == "MCS"], whichPGS = "SCZ", whichanalysis = "NDD vs MCS"),
  run_glm(mydatatable = df[ndd == 1 | cohort == "ALSPAC"], whichPGS = "SCZ", whichanalysis = "NDD vs ALSPAC"),
  
  # NDD
  run_glm(mydatatable = df, whichPGS = "NDD", whichanalysis = "mega"),
  run_glm(mydatatable = df[cohort == "GEL"], whichPGS = "NDD", whichanalysis = "GEL only"),
  run_glm(mydatatable = df[ndd == 1 | (ndd == 0 & cohort == "GEL")], whichPGS = "NDD", whichanalysis = "NDD vs GEL control"),
  run_glm(mydatatable = df[ndd == 1 | cohort == "MCS"], whichPGS = "NDD", whichanalysis = "NDD vs MCS"),
  run_glm(mydatatable = df[ndd == 1 | cohort == "ALSPAC"], whichPGS = "NDD", whichanalysis = "NDD vs ALSPAC")
)



## trio model correcting for prematurity:
formula_trio_cov = "ndd ~ child_PGS + mother_PGS + father_PGS + prematurity"


## trio model correcting for RVBS:
formula_trio_cov = "ndd ~ child_PGS + mother_PGS + father_PGS + child_rare + mother_rare + father_rare"


## two-sided z-score test to compare if non-transmitted coefficients are significantly different
z_test = function(estimate1, estimate2, se1, se2) {
  zscore = (estimate1 - estimate2) / sqrt(se1^2 + se2^2)
  # two-sided p-value
  twosidedp = 2*pnorm(abs(zscore), lower.tail=FALSE)
  return(twosidedp)
}

