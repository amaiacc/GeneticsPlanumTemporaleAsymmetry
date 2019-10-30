# get cluster sizes

# Xiangzhen ran the VBM analsyis and association with SNPs of interst

# use 3dClustSim to get clusters of associated voxels
## https://afni.nimh.nih.gov/pub/dist/doc/program_help/3dClustSim.html
##3dClustSim computes a cluster-size threshold for a given voxel-wise p-value threshold, 
## such that the probability of anything surviving the dual thresholds 
## is at some given level (specified by the '-athr' option).

PATH=/usr/local/apps/afni/:${PATH}
subset_name="imagingT1_N18057"
working_dir=/data/workspaces/lag/workspaces/lg-ukbiobank/working_data/amaia/genetic_data/PT/${subset_name}/lead_snps/voxel_vertex/
vbm_dir=/data/workspaces/lag/workspaces/lg-ukbiobank/analysis/xiangzhen/vol_asym/analysis/

cd ${working_dir}
for snp in rs7420166 rs41298373
do
 zfile=${vbm_dir}/${snp}/glmout_ai/rs.contrast/z.nii.gz
 3dClustSim -mask ${zfile} -fwhm 6 -cmd 3dClustSim_${snp}.cmd
done
