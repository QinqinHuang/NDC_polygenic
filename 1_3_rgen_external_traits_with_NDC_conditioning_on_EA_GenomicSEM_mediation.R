#-------------------------------------------------------
# 2024-09-13 last modified
# 
# (1) run GenomicSEM to estimate rgen between EA and NonEA
# latent factors and a latent variable for an external trait
# (2) mediation analysis to estimate the contribution of the EA
# and NonEA latent factors to the genetic correlation
# between NDCs and an external trait
#
# input GWAS: EA, NDD, and a third trait
#-------------------------------------------------------
library(GenomicSEM)

# reference data
path_hm3 = "/lustre/scratch126/humgen/teams/martin/users/qh1/Downloaded_data/LDScores/eur_w_ld_chr/w_hm3.noMHC.snplist"
dir_LD = "/lustre/scratch126/humgen/teams/martin/users/qh1/Downloaded_data/LDScores/eur_w_ld_chr/"
# prepared by the authors: 
#$wget https://utexas.box.com/s/vkd36n197m8klbaio3yzoxsee6sxo11v
path_1kGMAF = "/lustre/scratch126/humgen/teams/martin/users/qh1/Downloaded_data/GenomicSEM_reference_data/reference.1000G.maf.0.005.txt"



#----- (1) function to run rgen -----
# based on github: https://github.com/PerlineDemange/non-cognitive/blob/master/GenomicSEM/Genetic%20correlations/Without%20using%20SNP%20effects/function_rG_woSNP.R

LDSCoutput_rG_woSNP <- function(
    data1="/lustre/scratch123/hgi/projects/ddd/users/qh1/GWAS_by_subtraction/GWAS_munged/EA.sumstats.gz", 
    data2="/lustre/scratch123/hgi/projects/ddd/users/qh1/GWAS_by_subtraction/GWAS_munged/NDD.sumstats.gz", 
    data3, sample.prev_trait, population.prev_trait, trait.name){
  
  traits <- c(data1, data2, data3)
  sample.prev <- c(NA,0.5,sample.prev_trait) 
  population.prev <- c(NA,0.01,population.prev_trait) 
  ld<-dir_LD
  wld <- dir_LD
  trait.names<-c("EA", "NDD", "T")
  LDSCoutput<- ldsc(traits, sample.prev, population.prev, ld, wld, trait.names=c("EA", "NDD", "T"))
  
  modelrg<-'E=~NA*NDD + start(0.4)*EA
  NE=~NA*NDD
  LT=~NA*T
  NE~~1*NE
  E~~1*E
  LT~~1*LT
  EA~~0*EA
  NDD~~0*NDD
  T~~0*T
  E~~0*NE
  EA~~0*NDD
  EA~~0*T
  NDD~~0*T
  '
  output<-usermodel(LDSCoutput,estimation="DWLS",model=modelrg)
  saveRDS(output, file = paste0("NDD_minus_EA_with", trait.name, ".RDS"))  
  
  cog <- output$results[5,] #line of the results for the C~~LT
  noncog <- output$results[6,] #line of the results for the NC~~LT
  results <- cbind(trait.name, cog, noncog)
  results
}




#----- (2) calculate rgen conditioning on EA -----
externaltraits = rbind(
  LDSCoutput_rG_woSNP(data3 = "/lustre/scratch123/hgi/projects/ddd/users/qh1/GWAS_by_subtraction/GWAS_munged/CP.sumstats.gz", sample.prev_trait = NA, population.prev_trait = NA, trait.name = "CP"),
  LDSCoutput_rG_woSNP(data3 = "/lustre/scratch123/hgi/projects/ddd/users/qh1/GWAS_by_subtraction/GWAS_munged/NonCogEA.sumstats.gz", sample.prev_trait = NA, population.prev_trait = NA, trait.name = "NonCogEA"),
  LDSCoutput_rG_woSNP(data3 = "/lustre/scratch123/hgi/projects/ddd/users/qh1/GWAS_by_subtraction/GWAS_munged/ADHD.sumstats.gz", sample.prev_trait = 0.5, population.prev_trait = 0.05, trait.name = "ADHD"),
  LDSCoutput_rG_woSNP(data3 = "/lustre/scratch123/hgi/projects/ddd/users/qh1/GWAS_by_subtraction/GWAS_munged/SCZ.sumstats.gz", sample.prev_trait = 0.4175, population.prev_trait = 0.004, trait.name = "SCZ"),
  LDSCoutput_rG_woSNP(data3 = "/lustre/scratch123/hgi/projects/ddd/users/qh1/GWAS_by_subtraction/GWAS_munged/MDD.sumstats.gz", sample.prev_trait = 0.3414, population.prev_trait = 0.05, trait.name = "Depression"),
  LDSCoutput_rG_woSNP(data3 = "/lustre/scratch123/hgi/projects/ddd/users/qh1/GWAS_by_subtraction/GWAS_munged/TS.sumstats.gz", sample.prev_trait = 0.5, population.prev_trait = 0.003, trait.name = "TS"),
  LDSCoutput_rG_woSNP(data3 = "/lustre/scratch123/hgi/projects/ddd/users/qh1/GWAS_by_subtraction/GWAS_munged/PTSD.sumstats.gz", sample.prev_trait = 0.1329, population.prev_trait = 0.039, trait.name = "PTSD"),
  LDSCoutput_rG_woSNP(data3 = "/lustre/scratch123/hgi/projects/ddd/users/qh1/GWAS_by_subtraction/GWAS_munged/EatingDisorder.sumstats.gz", sample.prev_trait = 0.2343, population.prev_trait = 0.03, trait.name = "EatingDisorder"),
  LDSCoutput_rG_woSNP(data3 = "/lustre/scratch123/hgi/projects/ddd/users/qh1/GWAS_by_subtraction/GWAS_munged/Bipolar.sumstats.gz", sample.prev_trait = 0.1014, population.prev_trait = 0.01, trait.name = "Bipolar"),
  LDSCoutput_rG_woSNP(data3 = "/lustre/scratch123/hgi/projects/ddd/users/qh1/GWAS_by_subtraction/GWAS_munged/Epilepsy.sumstats.gz", sample.prev_trait = 0.5, population.prev_trait = 0.002, trait.name = "Epilepsy"),
  LDSCoutput_rG_woSNP(data3 = "/lustre/scratch123/hgi/projects/ddd/users/qh1/GWAS_by_subtraction/GWAS_munged/OCD.sumstats.gz", sample.prev_trait = 0.5, population.prev_trait = 0.01, trait.name = "OCD"),
  LDSCoutput_rG_woSNP(data3 = "/lustre/scratch123/hgi/projects/ddd/users/qh1/GWAS_by_subtraction/GWAS_munged/ASD.sumstats.gz", sample.prev_trait = 0.3966, population.prev_trait = 0.01, trait.name = "ASD")
)

write.csv(externaltraits, "GenomicSEM_rgen_woSNP_NDD_minus_EA_with_othertraits.csv", row.names = F)



#----- (3) calculate rgen with prenatal risk factor conditioning on EA -----
prenatalrisk = rbind(
  LDSCoutput_rG_woSNP(data3 = "/lustre/scratch123/hgi/projects/ddd/users/qh1/GWAS_by_subtraction/GWAS_munged/Preterm.sumstats.gz", sample.prev_trait = 0.06609, population.prev_trait = 0.07, trait.name = "Preterm"),
  LDSCoutput_rG_woSNP(data3 = "/lustre/scratch123/hgi/projects/ddd/users/qh1/GWAS_by_subtraction/GWAS_munged/SmkInit.sumstats.gz", sample.prev_trait = 0.4888, population.prev_trait = 0.133, trait.name = "SmkInit"),
  LDSCoutput_rG_woSNP(data3 = "/lustre/scratch123/hgi/projects/ddd/users/qh1/GWAS_by_subtraction/GWAS_munged/DrinksWk.sumstats.gz", sample.prev_trait = NA, population.prev_trait = NA, trait.name = "DrinksWk"),
  LDSCoutput_rG_woSNP(data3 = "/lustre/scratch123/hgi/projects/ddd/users/qh1/GWAS_by_subtraction/GWAS_munged/SleepApnea_MTAG.sumstats.gz", sample.prev_trait = 0.5, population.prev_trait = 0.04, trait.name = "SleepApnea_MTAG"),
  LDSCoutput_rG_woSNP(data3 = "/lustre/scratch123/hgi/projects/ddd/users/qh1/GWAS_by_subtraction/GWAS_munged/GestationlHTN.sumstats.gz", sample.prev_trait = 0.5, population.prev_trait = 0.1, trait.name = "GestationlHTN")
)

write.csv(prenatalrisk, "GenomicSEM_rgen_woSNP_NDD_minus_EA_with_prenatal.csv", row.names = F)



#----- (4) function to run the mediation analysis -----
# modified from: https://github.com/PerlineDemange/non-cognitive/blob/master/GenomicSEM/Genetic%20correlations/Decomposition%20of%20the%20correlation/function_rg_mediation.R

LDSCoutput_rG_mediation <- function(data1="/lustre/scratch123/hgi/projects/ddd/users/qh1/GWAS_by_subtraction/GWAS_munged/EA.sumstats.gz", data2="/lustre/scratch123/hgi/projects/ddd/users/qh1/GWAS_by_subtraction/GWAS_munged/NDD.sumstats.gz", data3, sample.prev_trait, population.prev_trait, trait.name){
  
  traits <- c(data1, data2, data3)
  sample.prev <- c(NA, 0.5, sample.prev_trait) 
  population.prev <- c(NA, 0.01, population.prev_trait) 
  ld <- dir_LD
  wld <- dir_LD
  trait.names <- c("EA", "NDD", "T")
  LDSCoutput <- ldsc(traits, sample.prev, population.prev, ld, wld, trait.names = c("EA", "NDD", "T"))
  
  #replace EA by NDD, replace CP by EA, replace C by E, and NC by NE
  model <- 'E=~NA*NDD + start(0.4)*EA + a*NDD
  NE=~NA*NDD + b*NDD
  LT=~NA*T + e*T
  NE~~1*NE
  E~~1*E
  LT~~1*LT
  EA~~0*EA
  NDD~~0*NDD
  T~~0*T
  E~~0*NE
  EA~~0*NDD
  EA~~0*T
  NDD~~0*T
  E~~c*LT
  NE~~d*LT
  EC:=a*c*e
  ENC:=b*d*e
  total:=a*c*e + b*d*e'
  
  output <- usermodel(LDSCoutput, estimation="DWLS", model=model)
  saveRDS(output, file = paste0("mediation_NDD_minus_EA_with_", trait.name, ".RDS"))
  
  cog <- output$results[7,] #line of the results for the  EC :=       a*c*e 
  noncog <- output$results[8,] #line of the results for the  ENC :=       b*d*e 
  total <- output$results[9,] #line of the results for total := a*c*e+b*d*e
  results <- cbind(trait.name, cog, noncog, total)
  results 
}


#----- (5) run the decomposition analysis of rgen -----
rg_mediation_results = rbind(
  LDSCoutput_rG_mediation(data3 = "/lustre/scratch123/hgi/projects/ddd/users/qh1/GWAS_by_subtraction/GWAS_munged/CP.sumstats.gz", sample.prev_trait = NA, population.prev_trait = NA, trait.name = "CP"),
  LDSCoutput_rG_mediation(data3 = "/lustre/scratch123/hgi/projects/ddd/users/qh1/GWAS_by_subtraction/GWAS_munged/NonCogEA.sumstats.gz", sample.prev_trait = NA, population.prev_trait = NA, trait.name = "NonCogEA"),
  LDSCoutput_rG_mediation(data3 = "/lustre/scratch123/hgi/projects/ddd/users/qh1/GWAS_by_subtraction/GWAS_munged/ADHD.sumstats.gz", sample.prev_trait = 0.5, population.prev_trait = 0.05, trait.name = "ADHD"),
  LDSCoutput_rG_mediation(data3 = "/lustre/scratch123/hgi/projects/ddd/users/qh1/GWAS_by_subtraction/GWAS_munged/SCZ.sumstats.gz", sample.prev_trait = 0.4175, population.prev_trait = 0.004, trait.name = "SCZ"),
  LDSCoutput_rG_mediation(data3 = "/lustre/scratch123/hgi/projects/ddd/users/qh1/GWAS_by_subtraction/GWAS_munged/MDD.sumstats.gz", sample.prev_trait = 0.3414, population.prev_trait = 0.05, trait.name = "Depression"),
  LDSCoutput_rG_mediation(data3 = "/lustre/scratch123/hgi/projects/ddd/users/qh1/GWAS_by_subtraction/GWAS_munged/TS.sumstats.gz", sample.prev_trait = 0.5, population.prev_trait = 0.003, trait.name = "TS"),
  LDSCoutput_rG_mediation(data3 = "/lustre/scratch123/hgi/projects/ddd/users/qh1/GWAS_by_subtraction/GWAS_munged/PTSD.sumstats.gz", sample.prev_trait = 0.1329, population.prev_trait = 0.039, trait.name = "PTSD"),
  LDSCoutput_rG_mediation(data3 = "/lustre/scratch123/hgi/projects/ddd/users/qh1/GWAS_by_subtraction/GWAS_munged/Bipolar.sumstats.gz", sample.prev_trait = 0.1014, population.prev_trait = 0.01, trait.name = "Bipolar"),
  LDSCoutput_rG_mediation(data3 = "/lustre/scratch123/hgi/projects/ddd/users/qh1/GWAS_by_subtraction/GWAS_munged/EatingDisorder.sumstats.gz", sample.prev_trait = 0.2343, population.prev_trait = 0.03, trait.name = "EatingDisorder"),
  LDSCoutput_rG_mediation(data3 = "/lustre/scratch123/hgi/projects/ddd/users/qh1/GWAS_by_subtraction/GWAS_munged/OCD.sumstats.gz", sample.prev_trait = 0.5, population.prev_trait = 0.01, trait.name = "OCD"),
  LDSCoutput_rG_mediation(data3 = "/lustre/scratch123/hgi/projects/ddd/users/qh1/GWAS_by_subtraction/GWAS_munged/Epilepsy.sumstats.gz", sample.prev_trait = 0.5, population.prev_trait = 0.002, trait.name = "Epilepsy"),
  LDSCoutput_rG_mediation(data3 = "/lustre/scratch123/hgi/projects/ddd/users/qh1/GWAS_by_subtraction/GWAS_munged/ASD.sumstats.gz", sample.prev_trait = 0.3966, population.prev_trait = 0.01, trait.name = "ASD"),
  
  # prenatal risk factor
  LDSCoutput_rG_mediation(data3 = "/lustre/scratch123/hgi/projects/ddd/users/qh1/GWAS_by_subtraction/GWAS_munged/SmkInit.sumstats.gz", sample.prev_trait = 0.4888, population.prev_trait = 0.133, trait.name = "SmkInit"),
  LDSCoutput_rG_mediation(data3 = "/lustre/scratch123/hgi/projects/ddd/users/qh1/GWAS_by_subtraction/GWAS_munged/Preterm.sumstats.gz", sample.prev_trait = 0.06609, population.prev_trait = 0.07, trait.name = "Preterm"),
  LDSCoutput_rG_mediation(data3 = "/lustre/scratch123/hgi/projects/ddd/users/qh1/GWAS_by_subtraction/GWAS_munged/GestationlHTN.sumstats.gz", sample.prev_trait = 0.5, population.prev_trait = 0.1, trait.name = "GestationlHTN"),
  LDSCoutput_rG_mediation(data3 = "/lustre/scratch123/hgi/projects/ddd/users/qh1/GWAS_by_subtraction/GWAS_munged/SleepApnea_MTAG.sumstats.gz", sample.prev_trait = 0.5, population.prev_trait = 0.04, trait.name = "SleepApnea_MTAG"),
  LDSCoutput_rG_mediation(data3 = "/lustre/scratch123/hgi/projects/ddd/users/qh1/GWAS_by_subtraction/GWAS_munged/DrinksWk.sumstats.gz", sample.prev_trait = NA, population.prev_trait = NA, trait.name = "DrinksWk")
)


write.csv(rg_mediation_results, "GenomicSEM_rgen_mediation_woSNP_NDD_minus_EA_with_othertraits.csv", row.names = F)




