#01-02-2019
#By: Amaia Carrion Castillo

# Scripts to run ldsc (https://github.com/bulik/ldsc) and calculate SNP heritability and genetic correlations from summary statistics.
# These scripts were used for the UKB PT project.
# They come with no guarantees.
# In case of questions, feel free to contact me: amaia.carrioncastillO@mpi.nl

LDSC_GenCor_GWASes_bgenie_PT_N18057.sh
# For PT phenotypes (AI, L and R)
# preliminaries + h2 running
# get genetic correlations between my phenotypes and publicly available summary statistics (which were already prepared for input)

LDSC_GenCor_GWASes_bgenie_PT_perSex_N18057.sh
# For PT phenotypes (AI, L and R), stratified per sex.
# preliminaries + h2 running
# get genetic correlations between my phenotypes and publicly available summary statistics (which were already prepared for input)

LDSC_GenCor_GWASes_bgenie_PTadjTBV_N18057.sh
# For PT phenotypes (AI, L and R), adjusted for total brain volume.
# preliminaries + h2 running
# get genetic correlations between my phenotypes and publicly available summary statistics (which were already prepared for input)
