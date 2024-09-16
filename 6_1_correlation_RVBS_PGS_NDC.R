#----------------------------------------------------------------------
# 2024-09-13 last modified
#
# Aim: to calculate correlations between
# (1) PC adjusted RVBS in one parent and PGS in the other parent
# (2) PC adjusted RVBS and PGS in the same parent
# (3) PC adjusted RVBS in the child and child's PGS
#
# in these trios:
# NDD trios where parents are unaffected
# NDD undiagnosed trios where parents are unaffected
# NDD trios with a de novo diagnosis and unaffected parents
#----------------------------------------------------------------------
library(data.table)
library(foreach)
library(dplyr)


# input table for RVBS residuals (adjusted for 20 PCs) and PGS
# each row is a trio family
RVBS = fread("/re_gecip/paediatrics/QQHuang/pgs_paper/data/rare_variants/RVBS_adj20PCs_trios_DDD2390_GEL2174_masked_FMISSING0.05.txt")


# make a long table; each row is an individual
RVBS_mother = select(RVBS, ends_with("stable_id"), undiagnosed, unaffected_parents, diagnosed_denovo, starts_with("mother"))
RVBS_mother$type = "mother"
names(RVBS_mother) = gsub("mother_", "", names(RVBS_mother))
names(RVBS_mother) = gsub("father_", "partner_", names(RVBS_mother))

RVBS_father = select(RVBS, ends_with("stable_id"), undiagnosed, unaffected_parents, diagnosed_denovo, starts_with("father"))
RVBS_father$type = "father"
names(RVBS_father) = gsub("father_", "", names(RVBS_father))
names(RVBS_father) = gsub("mother_", "partner_", names(RVBS_father))

if(identical(RVBS_mother$partner_stable_id, RVBS_father$stable_id)) {
  partner = RVBS_father[, .(EA, EANonCog, CP, SCZ, NDD)]
  names(partner) = paste0("partner_", names(partner))
  RVBS_mother = cbind(RVBS_mother, partner)
}
if(identical(RVBS_father$partner_stable_id, RVBS_mother$stable_id)) {
  partner = RVBS_mother[, .(EA, EANonCog, CP, SCZ, NDD)]
  names(partner) = paste0("partner_", names(partner))
  RVBS_father = cbind(RVBS_father, partner)
}

RVBS_con = rbind(RVBS_mother, RVBS_father)

RVBS_child = select(RVBS, ends_with("stable_id"), undiagnosed, unaffected_parents, diagnosed_denovo, starts_with("child"))
RVBS_child$type = "child"
names(RVBS_child) = gsub("child_", "", names(RVBS_child))
RVBS_child[, ':=' (child_stable_id = stable_id, partner_stable_id = NA, mother_stable_id = NULL, father_stable_id = NULL, partner_EA = NA, partner_EANonCog = NA, partner_CP = NA, partner_SCZ = NA, partner_NDD = NA)]

RVBS_con = rbind(RVBS_con, RVBS_child)



#----- correlation between RVBS and PGS or partner's PGS -----
## a function to calculate correlation in the parents
cor_wrapper = function(df, myvar, myPGS, flag = NA) {
  # normalise column name
  df$RVBS = df[, get(myvar)]
  df$PGS = df[, get(myPGS)]
  df$partner_PGS = df[, get(paste0("partner_", myPGS))]
  df = df[!is.na(PGS)]
  
  t1 = cor.test(df$RVBS, df$PGS)
  d1 = data.table(subset = flag, N = nrow(df), var = myvar, PGS = myPGS, person = "same parent",
                  cor = t1$estimate, stat = t1$statistic, p = t1$p.value)
  t2 = cor.test(df$RVBS, df$partner_PGS)
  d2 = data.table(subset = flag, N = nrow(df), var = myvar, PGS = myPGS, person = "partner",
                  cor = t2$estimate, stat = t2$statistic, p = t2$p.value)
  
  return(rbind(d1, d2))
}

## a function to calculate correlation in the probands
cor_wrapper_probands = function(df, myvar, myPGS, flag = NA) {
  # normalise column name
  df$RVBS = df[, get(myvar)]
  df$PGS = df[, get(myPGS)]
  df = df[!is.na(PGS)]
  
  t1 = cor.test(df$RVBS, df$PGS)
  d1 = data.table(subset = flag, N = nrow(df), var = myvar, PGS = myPGS, person = "child",
                  cor = t1$estimate, stat = t1$statistic, p = t1$p.value)
  return(d1)
}


## calculate correlation in various subsets
corPGS = foreach(ii = c("EA", "EANonCog", "CP", "SCZ", "NDD"), .combine = rbind) %do% {
  foreach(jj = c("PTV_ddg2p_mono_lof_adj20PCs", "both_ddg2p_mono_lof_adj20PCs", "synonymous_ddg2p_mono_lof_adj20PCs", "PTV_pLI_adj20PCs", "both_pLI_adj20PCs", "synonymous_pLI_adj20PCs"), .combine = rbind) %do% {
    rbind(
      # NDD trios with unaffected parents
      cor_wrapper(RVBS_con[unaffected_parents == 1 & type != "child"], myvar = jj, myPGS = ii, flag = "NDD unaffected parents"),
      cor_wrapper_probands(RVBS_con[unaffected_parents == 1 & type == "child"], myvar = jj, myPGS = ii, flag = "NDD unaffected parents"),
      # NDD trios, undiganosed probands who have unaffected parents
      cor_wrapper(RVBS_con[unaffected_parents == 1 & undiagnosed == 1 & type != "child"], myvar = jj, myPGS = ii, flag = "undiagnosed NDD unaffected parents"),
      cor_wrapper_probands(RVBS_con[unaffected_parents == 1 & undiagnosed == 1 & type == "child"], myvar = jj, myPGS = ii, flag = "undiagnosed NDD unaffected parents"),
      # NDD trios, probands with a de novo diagnosed and unaffected parents
      cor_wrapper(RVBS_con[unaffected_parents == 1 & undiagnosed == 0 & diagnosed_denovo == 1 & type != "child"], myvar = jj, myPGS = ii, flag = "De novo diagnoses unaffected parents"),
      cor_wrapper_probands(RVBS_con[unaffected_parents == 1 & undiagnosed == 0 & diagnosed_denovo == 1 & type == "child"], myvar = jj, myPGS = ii, flag = "De novo diagnoses unaffected parents")
      
    )
  }
} 



