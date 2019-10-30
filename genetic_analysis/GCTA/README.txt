#01-02-2019
#By: Amaia Carrion Castillo

# Scripts to run GCTA (#http://cnsgenomics.com/software/gcta/) to calculate SNP heritability and genetic correlations from genotyping data.
# These scripts were used for the UKB PT project.
# They come with no guarantees.
# Note: Paths to intput files may have been re-structured.
# In case of questions, feel free to contact me: amaia.carrioncastillO@mpi.nl

#----------------------#
# GCTA-REML
#----------------------#
GCTA_REML_calUKBv2_PT_N18057.sh
# For PT phenotypes (AI, L and R)
# h2 + genetic correlations
# get genetic correlations between the phenotypes.

GCTA_REML_calUKBv2_PT_perSex_N18057.sh
# For PT phenotypes (AI, L and R), stratified per sex.
# h2 + genetic correlations
# get genetic correlations between the phenotypes.

GCTA_REML_calUKBv2_PTadjh2_N18057.sh
# For PT phenotypes (AI, L and R), adjusted for total brain volume.
# h2 + genetic correlations
# get genetic correlations between the phenotypes.
#----------------------#

#----------------------#
# GCTA-REML output
#----------------------#
# R scripts to parse REML output.

Plot_hsq_results_PT_N18057.R
# Generates plots for PT reml results.
# Note: requires to run an additional script (for all volumes), which will generate input files for this script:
## /data/workspaces/lag/workspaces/lg-ukbiobank/analysis/amaia/genetic_analysis/release_v2/imagingT1_N18057/Plot_hsq_results_N18057.R

Plot_parse_hsq_results_PT_TBV_18057.R
# Parses output and generates plots for the PTadjTBV reml results.

#----------------------#
# GCTA-GREML Power Calculator
#----------------------#
GCTA_PowerCalculator_gencor_GWASes.R
# GCTA-GREML Power Calculator for specific analyses (http://cnsgenomics.com/shiny/gctaPower/)