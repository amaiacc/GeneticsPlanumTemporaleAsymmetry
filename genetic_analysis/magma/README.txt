#01-02-2019
#By: Amaia Carrion Castillo

# Scripts to run MAGMA (https://ctg.cncr.nl/software/magma) for gene and gene-set enrichment analysis.
# These scripts were used for the UKB PT project.
# They come with no guarantees.
# In case of questions, feel free to contact me: amaia.carrioncastillO@mpi.nl


magma_template.sh
# Template script to run magma that can be used to run in the grid.

magma_run_PT.sh
# Script I used to run magma_template.sh with specific settings for the PT project.

magma_summarize_output_PT.R
# R script to read magma output, generate some plots and tables of interest.