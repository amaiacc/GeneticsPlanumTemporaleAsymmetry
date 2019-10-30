base_data_path='/homes_unix/pepe/workspace/VBM_SNP_for_Amaia/data' 
base_code_dir='/homes_unix/pepe/workspace/VBM_SNP_for_Amaia/code
base_out_dir ='/homes_unix/pepe/workspace/VBM_SNP_for_Amaia/out_plots

MNI152_T1_2mm_nii_gz='/srv/shares/softs/fsl/data/standard/MNI152_T1_2mm.nii.gz'

declare -a model_name=('rs7420166' 'rs41298373')

my_template='fsaverage_sym'

export SUBJECTS_DIR=${base_data_path}/SUBJECTS_DIR
mkdir -p  ${SUBJECTS_DIR}




###########################################################
# Plotting settings
###########################################################
my_surf='pial';
my_hemi='lh' 
min_threshold='3'  
max_threshold='7'  
surf_overlay='pial'




###########################################################
# Upsampling and change image orientation of the FSL image
###########################################################

cp ${MNI152_T1_2mm_nii_gz} ${base_data_path}
cp -r ${FREESURFER_HOME}/subjects/${my_template} ${SUBJECTS_DIR}

mri_convert ${MNI152_T1_2mm_nii_gz} ${base_data_path}/MNI152_T1_2mm.mgz --upsample 2 --in_orientation LAS --out_orientation LIA

mri_robust_register --mov ${base_data_path}/MNI152_T1_2mm.mgz \
--dst ${SUBJECTS_DIR}/${my_template}/mri/T1.mgz  \
--lta ${base_data_path}/MNI152_T1_2mm_to_fsaverage_sym.lta \
--mapmov ${base_data_path}/MNI152_T1_2mm_to_fsaverage_sym.mgz \
--satit




###########################################################
# 
###########################################################
for current_model_name in "${model_name[@]}"; do

current_vol_filename=${base_data_path}/${current_model_name}'/z.nii.gz'

if [ -f ${current_vol_filename} ]; then
echo processing file ${current_vol_filename} 

cp -r ${FREESURFER_HOME}/subjects/${my_template} ${SUBJECTS_DIR}
mv ${SUBJECTS_DIR}/${my_template} ${SUBJECTS_DIR}/fake_subject
cp -r ${FREESURFER_HOME}/subjects/${my_template} ${SUBJECTS_DIR}

current_vol_filename_nii_gz=${SUBJECTS_DIR}/fake_subject/mri/${current_model_name}_z.nii.gz
current_vol_filename_mgz=${SUBJECTS_DIR}/fake_subject/mri/${current_model_name}_z.mgz
current_surf_filename=${SUBJECTS_DIR}/fake_subject/surf/${current_model_name}_z_pial.mgh
current_temp_filename_mgz=${SUBJECTS_DIR}/${my_template}/mri/${current_model_name}_z.mgz

cp ${current_vol_filename} ${current_vol_filename_nii_gz} 
mri_convert ${current_vol_filename_nii_gz} ${current_vol_filename_mgz} --upsample 2 --in_orientation LAS --out_orientation LIA
cp ${current_vol_filename_mgz} ${current_temp_filename_mgz} 

## assign values from a volume to each surface vertex 
mri_vol2surf --src ${current_vol_filename_mgz} \
--out ${current_surf_filename} \
--reg ${base_data_path}/MNI152_T1_2mm_to_fsaverage_sym.lta \
--hemi 'lh' \
--surf 'pial' \
--sd ${SUBJECTS_DIR} \
--srcsubject fake_subject 


matlab -nosplash -nodisplay -nodesktop -r "plot_individual_maps_AP('${base_data_path}','${SUBJECTS_DIR}','${base_code_dir}','${base_out_dir}', '${current_model_name},'${file_name}','${my_template}','${my_surf}','fake_subject','lh', '${min_threshold}', '${max_threshold}', '${surf_overlay}','${my_surf}'); exit"

fi
done
whos
