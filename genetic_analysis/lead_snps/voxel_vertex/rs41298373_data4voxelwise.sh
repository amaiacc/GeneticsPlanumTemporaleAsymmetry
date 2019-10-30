# Prepare data to perform a voxelwise association between rs41298373 and AI within PT

#----------------------------------------------------------------------
# Define paths and directories
#----------------------------------------------------------------------
resource_dir=/data/workspaces/lag/shared_spaces/Resource_DB/
analysis_dir=/data/workspaces/lag/workspaces/lg-ukbiobank/analysis/amaia/genetic_analysis/PT/
primary_dir=/data/workspaces/lag/workspaces/lg-ukbiobank/working_data/amaia/genetic_data/release_v2/cal/QC/
working_dir=/data/workspaces/lag/workspaces/lg-ukbiobank/working_data/amaia/genetic_data/PT/
UKB_phenos_dir=/data/workspaces/lag/workspaces/lg-ukbiobank/working_data/amaia/pheno_files/genetic_v2/ukb21288_ukb21293/summary_phenotypes/
vbm_dir=/data/
#----------------------------------------------------------------------
# define pheno/cov file
pheno_file=${UKB_phenos_dir}ukb21288_ukb21293_imaging_noBioCovsVolumePhenotypes_residuals_wHeader.table

# define snp of interest
snp=rs41298373
chr=10

mkdir -p ${working_dir}${snp}
cd ${working_dir}${snp}
# extract genotype data for this snp using plink
plink -bfile ${primary_dir}ukb_cal_snpQC_imagingT1_N12245_clean_${chr} --snp ${snp} --recode --out ${working_dir}${snp}/ukb_cal_snpQC_imagingT1_N12245_${snp}


# get list of subject IDs for which VBM data is available
cp /data/clusterfs/lag/projects/lg-ukbiobank/working_data/vbm_sids.txt ${working_dir}${snp}/


# run this script to generate summary csv with genotype data for this snp, and all the covariates
Rscript ${analysis_dir}/rs41298373_data4voxelwise.R

# Xiangzhen has run the association of the on the VBM
##bbc-design -behcsv rs41298373_covariates_4voxelwise.csv -contvar rs41298373add sex age zage2 PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10 -catevar array -odt freesurfer