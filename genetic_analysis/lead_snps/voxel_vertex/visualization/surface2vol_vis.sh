#
working_dir=/data/workspaces/lag/workspaces/lg-ukbiobank/working_data/amaia/genetic_data/PT/imagingT1_N18057/lead_snps/voxel_vertex/
data_dir=/data/workspaces/lag/workspaces/lg-ukbiobank/analysis/xiangzhen/vol_asym/analysis/

# "Results for each SNP are in the corresponding folder."
## ./glmout_ai/rs.contrast/z.nii.gz
## .//glmout_vbm/rs.contrast/z.nii.gz


# use mri_vol2surf to visualize VBM results
## https://surfer.nmr.mgh.harvard.edu/fswiki/mri_vol2surf
## mri_vol2surf [<options>] --src inputfile --out outpufile --srcreg registrationfile --hemi hemisphere 
snp=rs41298373
cd ${data_dir}${snp}/glmout_ai/rs.contrast/
cd ${working_dir}


mri_vol2surf --src ${data_dir}${snp}/glmout_ai/rs.contrast/  --src_type analyze --hemi lh --out ${working_dir}${snp}-lh.img


# Antonietta:
#I use surfice (https://www.nitrc.org/projects/surfice/)
#To convert mgh files to gii # ask Xiangzhen how does he do this?
#mris_convert -c ${data_dir}${snp}/glmout_ai/rs.contrast/z.mgh $FREESURFER_HOME/subjects/fsaverage_sym/surf/rh.white ${working_dir}${snp}.gii

# I already have the gii.gz, so I can use surfice directly
## https://www.nitrc.org/plugins/mwiki/index.php/surfice:MainPage#Overlay_a_normalized_SPM_or_FSL_statistical_map
# 1-    Launch Surf Ice
# 2-    Choose File/Open and select the image "BrainMesh_ICBM152.gii" - this is a template brain which is distributed with BrainNet.
# 3-    Choose Overlays/AddOverlay and select the image "motor_4t95.nii.gz" - this is a SPM map showing brain activation during a finger tapping task.
## choose files within ${data_dir}${snp}/glmout_ai/rs.contrast/
# 4-    A new panel appears in the toolbar named "Overlays": set the "min" value to be 4.95, to reflect the fact that in this case SPM reported that T values greater than 4.95 survive p <0.05 when corrected for multiple comparisons.
# 4-    Your image should look like the one shown above. 