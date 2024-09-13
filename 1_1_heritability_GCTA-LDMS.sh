#------------------------------------------------------------
# 2024-09-13 last modified
# 
# run GCTA (GREML-LDMS) to estimate heritability using DDD
# imputed data
#------------------------------------------------------------

# (1) calculate LD scores 
# 16253 samples including 9270 UKHLS controls and 6983 cases before excluding Scottish samples
bsub -J GCTA -o bjob_LDscore.out -e bjob_LDscore.error \
-R "select[mem>36000] rusage[mem=36000]" -M 36000 -n4 -R "span[hosts=1]" -q long \
"gcta --bfile /lustre/scratch123/hgi/projects/ddd_genotype/qh1/DDD_CoreExome_HRC_from_Mari/\
Well_imputed_INFO0.8_MAF0.01/HRC_GWAS_INFO_0.8_MAF_0.01_excldup_hwe6_diffmissing \
--thread-num 4 --ld-score-region 200 --out test"


# (2) stratifying SNPs by individual LD scores (this is done in R)
library(data.table)
lds_seg = fread("test.score.ld")

medianscore = median(lds_seg$ldscore_SNP)

fwrite(lds_seg[ldscore_SNP <= medianscore, .(SNP)], "snp_LD_group1.txt", row.names=F, col.names=F)
fwrite(lds_seg[ldscore_SNP > medianscore, .(SNP)], "snp_LD_group2.txt", row.names=F, col.names=F)


# (3) computing GRMs using the stratified SNPs 
# 2 LD groups and 3 MAF groups (0.01-0.05, 0.05-0.1, ≥0.1)
# LD group1
bsub -J LD1 -o bjob_GRM_LD1MAF1.out -e bjob_GRM_LD1MAF1.error \
-n4 -R "span[hosts=1]" -R "select[mem>4000] rusage[mem=4000]" -M 4000 \
"gcta --bfile /lustre/scratch123/hgi/projects/ddd_genotype/qh1/DDD_CoreExome_HRC_from_Mari/\
Well_imputed_INFO0.8_MAF0.01/HRC_GWAS_INFO_0.8_MAF_0.01_excldup_hwe6_diffmissing \
--thread-num 4 --extract snp_LD_group1.txt --make-grm --maf 0.01 --max-maf 0.05 \
--out DDD_GWAS_MAF_0.01_0.05_LDgroup1"

bsub -J LD1 -o bjob_GRM_LD1MAF2.out -e bjob_GRM_LD1MAF2.error \
-n4 -R "span[hosts=1]" -R "select[mem>4000] rusage[mem=4000]" -M 4000 \
"gcta --bfile /lustre/scratch123/hgi/projects/ddd_genotype/qh1/DDD_CoreExome_HRC_from_Mari/\
Well_imputed_INFO0.8_MAF0.01/HRC_GWAS_INFO_0.8_MAF_0.01_excldup_hwe6_diffmissing \
--thread-num 4 --extract snp_LD_group1.txt --make-grm --maf 0.05 --max-maf 0.1 \
--out DDD_GWAS_MAF_0.05_0.1_LDgroup1"

bsub -J LD1 -o bjob_GRM_LD1MAF3.out -e bjob_GRM_LD1MAF3.error \
-n4 -R "span[hosts=1]" -R "select[mem>10000] rusage[mem=10000]" -M 10000 \
"gcta --bfile /lustre/scratch123/hgi/projects/ddd_genotype/qh1/DDD_CoreExome_HRC_from_Mari/\
Well_imputed_INFO0.8_MAF0.01/HRC_GWAS_INFO_0.8_MAF_0.01_excldup_hwe6_diffmissing \
--thread-num 4 --extract snp_LD_group1.txt --make-grm --maf 0.1 \
--out DDD_GWAS_MAF_0.1_LDgroup1"

# LD group2
bsub -J LD2 -o bjob_GRM_LD2MAF1.out -e bjob_GRM_LD2MAF1.error \
-n4 -R "span[hosts=1]" -R "select[mem>4000] rusage[mem=4000]" -M 4000 \
"gcta --bfile /lustre/scratch123/hgi/projects/ddd_genotype/qh1/DDD_CoreExome_HRC_from_Mari/\
Well_imputed_INFO0.8_MAF0.01/HRC_GWAS_INFO_0.8_MAF_0.01_excldup_hwe6_diffmissing \
--thread-num 4 --extract snp_LD_group2.txt --make-grm --maf 0.01 --max-maf 0.05 \
--out DDD_GWAS_MAF_0.01_0.05_LDgroup2"

bsub -J LD2 -o bjob_GRM_LD2MAF2.out -e bjob_GRM_LD2MAF2.error \
-n4 -R "span[hosts=1]" -R "select[mem>4000] rusage[mem=4000]" -M 4000 \
"gcta --bfile /lustre/scratch123/hgi/projects/ddd_genotype/qh1/DDD_CoreExome_HRC_from_Mari/\
Well_imputed_INFO0.8_MAF0.01/HRC_GWAS_INFO_0.8_MAF_0.01_excldup_hwe6_diffmissing \
--thread-num 4 --extract snp_LD_group2.txt --make-grm --maf 0.05 --max-maf 0.1 \
--out DDD_GWAS_MAF_0.05_0.1_LDgroup2"

bsub -J LD2 -o bjob_GRM_LD2MAF3.out -e bjob_GRM_LD2MAF3.error \
-n4 -R "span[hosts=1]" -R "select[mem>10000] rusage[mem=10000]" -M 10000 \
"gcta --bfile /lustre/scratch123/hgi/projects/ddd_genotype/qh1/DDD_CoreExome_HRC_from_Mari/\
Well_imputed_INFO0.8_MAF0.01/HRC_GWAS_INFO_0.8_MAF_0.01_excldup_hwe6_diffmissing \
--thread-num 4 --extract snp_LD_group2.txt --make-grm --maf 0.1 \
--out DDD_GWAS_MAF_0.1_LDgroup2"



# make a file "multi_GRMs.txt" that contains GRM names
#DDD_GWAS_MAF_0.01_0.05_LDgroup1
#DDD_GWAS_MAF_0.05_0.1_LDgroup1
#DDD_GWAS_MAF_0.1_LDgroup1
#DDD_GWAS_MAF_0.01_0.05_LDgroup2
#DDD_GWAS_MAF_0.05_0.1_LDgroup2
#DDD_GWAS_MAF_0.1_LDgroup2



# (4) run gcta to calculate heritability
bsub -J gcta -o bjob_reml.out -e bjob_reml.error -q long \
-n4 -R "span[hosts=1]" -R "select[mem>20000] rusage[mem=20000]" -M 20000 \
"gcta --reml --mgrm multi_GRMs.txt --thread-num 4 \
--pheno test.pheno.txt --prevalence 0.01 \
--covar test.covar.discrete.txt --qcovar test.covar.quantitative.txt \
--out DDD_GWAS_reml_h2"
