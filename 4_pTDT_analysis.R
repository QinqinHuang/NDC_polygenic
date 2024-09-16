#------------------------------------------------------------
# 2024-09-13 last modified
#
# run pTDT in trios with unaffected parents
# probands: undiagnosed, diagnosed, and these excluding
# autistic probands; 
# in undiagonosed male and female probands separately
#------------------------------------------------------------
library(data.table)
library(dplyr)
library(foreach)
source("/re_gecip/paediatrics/QQHuang/pgs_paper/scripts/functions__modified.R")



#----- function to run pTDT in a given sample -----
pTDT = function(myNDDtrios, PGSname, flag, cohort) {
  
  myNDDtrios$child_PGS = myNDDtrios[, get(paste0("child_", PGSname))]
  myNDDtrios$mother_PGS = myNDDtrios[, get(paste0("mother_", PGSname))]
  myNDDtrios$father_PGS = myNDDtrios[, get(paste0("father_", PGSname))]
  myNDDtrios = myNDDtrios[!is.na(child_PGS)]
  
  # Get the mean of the parents
  myNDDtrios[, mid_parent_PGS := (mother_PGS + father_PGS)/2]
  
  # Take the difference of the proband PRS from the mean parent
  # and standardise by the SD of the mean parent distribution
  myNDDtrios[, child_dev_PGS := (child_PGS - mid_parent_PGS)/sd(myNDDtrios$mid_parent_PGS)]
  
  # Conduct a two-sided one-sample t-test of the child deviation
  mytest = t.test(myNDDtrios$child_dev_PGS, mu = 0, alternative = "two.sided")
  
  # Put results in a table
  results_table = data.table(
    group = flag, analysis = cohort, PGS = PGSname, N_trios = nrow(myNDDtrios), 
    diff_mean = mytest$estimate,
    lower_CI = signif((mytest$conf.int)[1], digits = 4),
    upper_CI = signif((mytest$conf.int)[2], digits = 4),
    se = mytest$stderr,
    p_value = mytest$p.value
  )
  
  return(results_table)
}


#----- run pTDT -----
# trios: a table of trios with unaffected parents; each row is a trio family

results = foreach(ii = c("EA","CP","EANonCog","SCZ","NDD"), .combine = rbind) %do% {
  rbind(
    # Main analysis: undiagnosed probands, unaffected parents
    pTDT(trios[ndd == 1 & undiagnosed == 1 & unaffected_parents == 1], ii, "undiagnosed", cohort = "DDD+GEL"),
    # diagnosed probands, unaffected parents
    pTDT(trios[ndd == 1 & undiagnosed == 0 & unaffected_parents == 1], ii, "diagnosed", cohort = "DDD+GEL"),
    # excluding autistic, undiagnosed probands, unaffected parents
    pTDT(trios[ndd == 1 & undiagnosed == 1 & unaffected_parents == 1 & autistic == 0], ii, "undiagnosed_no_aut", cohort = "DDD+GEL"),
    # excluding autistic, diagnosed probands, unaffected parents
    pTDT(trios[ndd == 1 & undiagnosed == 0 & unaffected_parents == 1 & autistic == 0], ii, "diagnosed_no_aut", cohort = "DDD+GEL"),
    
    # female undiagnosed probands with unaffected parents
    pTDT(trios[ndd == 1 & undiagnosed == 1 & unaffected_parents == 1 & sex == "F"], ii, "undiagnosed_female", cohort = "DDD+GEL"),
    # male undiagnosed probands with unaffected parents
    pTDT(trios[ndd == 1 & undiagnosed == 1 & unaffected_parents == 1 & sex == "M"], ii, "undiagnosed_male", cohort = "DDD+GEL")
  )
}




