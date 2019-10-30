# define clusters from VBM analysis
#----------------------------------------------------------------------
# Define paths and directories
#----------------------------------------------------------------------
vbm_dir=/data/workspaces/lag/workspaces/lg-ukbiobank/analysis/xiangzhen/vol_asym/analysis/
working_dir=/data/workspaces/lag/workspaces/lg-ukbiobank/working_data/amaia/genetic_data/PT/${subset_name}/lead_snps/voxel_vertex/
subset_name=imagingT1_N18057
# extract clusters
for snp in rs41298373 rs7420166
 do
 # get all negative associations (directions is reversed, given SNP coding)
 cluster -i ${vbm_dir}${snp}/glmout_ai/rs.contrast/z.nii.gz  -t 3 --mm > ${working_dir}/${snp}_rs.contrast_neg_clusters.txt
 # get all positive associations (directions is reversed, given SNP coding)
 cluster -i ${vbm_dir}${snp}/glmout_ai/rs_neg.contrast/z.nii.gz  -t 3 --mm > ${working_dir}/${snp}_rs_neg.contrast_pos_clusters.txt
 done

 
# to convert a nifti into negative:
##  fslmaths z.nii.gz -mul -1 ~/z_minus.nii.gz
