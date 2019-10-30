# Prepare data to perform a voxelwise association between rs41298373/rs74201660 and AI within PT


pheno_root="ukb25465_ukb25468"
subset_name="imagingT1_N18057"
#----------------------------------------------------------------------
# Define paths and directories
#----------------------------------------------------------------------
resource_dir=/data/workspaces/lag/shared_spaces/Resource_DB/
analysis_dir=/data/workspaces/lag/workspaces/lg-ukbiobank/analysis/amaia/genetic_analysis/PT/${subset_name}/
geno_dir=/data/workspaces/lag/workspaces/lg-ukbiobank/working_data/amaia/genetic_data/PT/${subset_name}/lead_snps/genotypes/
working_dir=/data/workspaces/lag/workspaces/lg-ukbiobank/working_data/amaia/genetic_data/PT/${subset_name}/lead_snps/voxel_vertex/
UKB_phenos_dir=/data/workspaces/lag/workspaces/lg-ukbiobank/working_data/amaia/pheno_files/genetic_v2/${pheno_root}/summary_phenotypes/
# vbm_dir=/data/clusterfs/lag/users/xiakon/dat/ ## to be updated!


# create working_dir
mkdir -p ${working_dir}

# create list of subjects for which VBM is available
# ls /data/clusterfs/lag/users/xiakon/dat/* > ls /data/clusterfs/lag/projects/lg-ukbiobank/working_data/vbm_sids.txt

# get list of subject IDs for which VBM data is available
cp /data/clusterfs/lag/projects/lg-ukbiobank/working_data/vbm_sids.txt ${working_dir}/

#----------------------------------------------------------------------
# define pheno/cov file
pheno_file=${UKB_phenos_dir}${pheno_root}_imaging_noBioCovsVolumePhenotypes_residuals_wHeader.table



# run this script to generate summary csv with genotype data for this snp, and all the covariates
Rscript ${analysis_dir}/rs41298373_data4voxelwise.R

# Xiangzhen has run the association of the on the VBM
##bbc-design -behcsv rs41298373_covariates_4voxelwise.csv -contvar rs41298373add sex age zage2 PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10 -catevar array -odt freesurfer