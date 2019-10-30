#-----------------------------------------
Short summary of scripts/directories within lg-ukbiobank/analysis/amaia/phenotypes/ukb25465_ukb25468/
Written: 27.02.2019
By: amacar
#-----------------------------------------


#-----------------------------------------
# General
#-----------------------------------------
# input files created by: lg-ukbiobank/analysis/amaia/primary_analysis/*_sqc_fam_combine.R

# ukb25465_ukb25468_selectCols_subsetSamples.R
Select columns of interest and subset the sample, to e.g. subset containing imaging data.
If run using knitr generates an html report, saved into: ../reports


#-----------------------------------------
# Imaging
#-----------------------------------------

# ukb25465_ukb25468_extractPhenotypes_imagingSubset.R
Extraction of phenotypes and covariates, to be used in downstream (genetic) analyses - using the imaging subset.

# ukb25465_ukb25468_LM_residuals_imagingSubset.R
Use output from *extractPhenotypes_imagingSubset.R to run LM with covariates, and extract to residuals

# ukb25465_ukb25468_LM_residuals_imagingSubset_perSex.R
Stratified analysis per sex for specific region of interest (to define as parameter).
Use output from *extractPhenotypes_imagingSubset.R to run LM with covariates, and extract to residuals

# ukb25465_ukb25468_totalBV_LM_residuals_imagingSubset.R
For total brain volume phenotype:
Use output from *extractPhenotypes_imagingSubset.R to run LM with covariates, and extract to residuals

#-----------------------------------------
# Planum Temporale
#-----------------------------------------
## Exploration of phenotypes
Perform exploratory analyses on the phenotypes. Correlation across AIs, relation to h2 estimates, factor analysis, etc.
# ukb25465_ukb25468_PT_TBV_imagingSubset.R
# ukb25465_ukb25468_PTadj_totalBV_LM_residuals_imagingSubset.R


#-----------------------------------------
# Cognitive phenotypes
#-----------------------------------------
# ukb25465_ukb25468_extractPhenotypesCogn_all.R
Extraction of cognitive phenotypes and covariates, to be used in downstream (genetic) analyses - using the full dataset.
It also creates subsets of data: all, imaging_T1, non_imaging_T1.

