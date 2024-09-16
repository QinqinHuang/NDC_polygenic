# Dissecting the role of common variants in rare neurodevelopmental conditions

Preprint: https://www.medrxiv.org/content/10.1101/2024.03.05.24303772v1

This repository contains code to perform the analyses for this project.


## Heritability and genetic correlations
1_1_heritability_GCTA-LDMS.sh: calculates heritability of neurodevelopmental conditions (NDCs) using LD score regression.

1_2_heritability_LDSC_meta_GWAS.sh: calculates genetic correlations between NDCs and brain-related traits and conditions using LD score regression.

1_3_rgen_external_traits_with_NDC_conditioning_on_EA_GenomicSEM_mediation.R: calculates genetic correlations conditioning on educational attainment (EA) using GenomicSEM and the proportion of the genetic correlation explained by the latent genetic component of NDCs that is correlated with EA. 


## Testing differences in polygenic scores between subgroups
2_1_standardisation_PGS_t_test_comparing_subgroups.R: standardises polygenic scores (PGS) using control individuals, calculates average PGS in subgroups, and performs two-sample t-tests to compare subgroups.

2_2_mean_PGS_MCS_weighting.R: applies statistical weights adjusting for sampling and non-response bias in MCS to calculate weighted average PGS that are likely to be representative of the full UK population.

## PGS and monogenic diagnosis
3_PGS_genetic_diagnosis_factors.R: performs regression models to test for associations between PGSs and diagnostic status, and whether associations between the chance of getting a diagnosis and factors affecting it might be confounding, or be confounded by PGS.

## Polygenic transmission disequilibrium test (pTDT)
4_pTDT_analysis.R: performs pTDT to assess if unaffected parents over-transmit polygenic risk to undiagnosed probands.

## Association with non-transmitted alleles
5_trio_model_non_transmitted_coefficient.R: runs the proband-only model and the trio model, with the latter assessing whether parental alleles not transmitted to the child are correlated with the child's risk of NDCs.

## Correlation between common and rare variant risk
6_1_correlation_RVBS_PGS_NDC.R: calculates the correlation between the rare variant burden score (RVBS) in one parent and the PGS in the other parent, as well as the correlation within the same parent or in the same proband. 

6_2_correlation_RVBS_PGS_MCS_weighted.R: calculates similar correlations in a control cohort, with and without applying the weights to correct for sampling and non-response bias in MCS.

