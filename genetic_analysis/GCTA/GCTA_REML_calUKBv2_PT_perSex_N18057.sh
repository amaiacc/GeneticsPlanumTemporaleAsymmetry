#!/bin/bash

# Compute heritability estimates using GCTA
# #http://cnsgenomics.com/software/gcta/

# to run after GCTA_calUKBv2_imaging.sh (to generate GRMs and some test REMLs, also in: GCTA_REML.sh)

# Phenotypes: imaging derived phenotypes (IDPs) provided by UKB and AI of bilateral volumes, DTI tracts and DTI FA skeleton
# residualized for: c("assessment_centre","age","sex",paste("PC",1:10,sep=""),"array")

#----------------------------------------------------------------------
# Notes
#----------------------------------------------------------------------
## will need to include further technical and biological confounds, from Elliot et al. bioRxiv: 
## 		technical: scanner x,y,z, position, motion, height, volumetric scaling factor 
## 		biological: height, diastolic and systolic pressure
## others: ICV ???

## but, collider bias? keep in mind
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# Define paths and directories
#----------------------------------------------------------------------
resource_dir=/data/workspaces/lag/shared_spaces/Resource_DB/
analysis_dir=/data/workspaces/lag/workspaces/lg-ukbiobank/analysis/amaia/genetic_analysis/release_v2/
primary_dir=/data/workspaces/lag/workspaces/lg-ukbiobank/working_data/amaia/genetic_data/release_v2/cal/QC/
working_dir=/data/workspaces/lag/workspaces/lg-ukbiobank/working_data/amaia/genetic_data/release_v2/cal/h2/
working_dir2=/data/clusterfs/lag/users/amacar/ukb/
UKB_phenos_dir=/data/workspaces/lag/workspaces/lg-ukbiobank/working_data/amaia/pheno_files/genetic_v2/
# if necessary, run script to create templates for GREML analyses:
#bash ${analysis_dir}/GCTA_REML_generate_templates.sh

#----------------------------------------------------------------------
# define root that will be contained within generated input and output files; parameter to change if input files change
pheno_root=ukb25465_ukb25468
covs_name=noBioCovs_noAssessmentC
subset_name=imagingT1_N18057
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# copy the data (GRMs) to sge2 directory first, otherwise it won't run in the grid
#----------------------------------------------------------------------
mkdir -p ${working_dir2}/grm/ ${working_dir2}/reml ${working_dir2}/pheno ${working_dir2}/logs
mkdir -p ${working_dir}/grm/ ${working_dir}/reml ${working_dir}/pheno ${working_dir}/logs
# copy genetic data to working_dir2 in sge, if not present yet!
if [ ! -f ${working_dir2}/grm/ukb_cal_snpQC_${subset_name}_clean_rm025_adj.grm.bin ] ; then cp ${working_dir}/grm/ukb_cal_snpQC_${subset_name}_clean_rm025_adj.grm* ${working_dir2}/grm/ ; fi
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# Define phenotypes and covariate files: Run $analysis_dir/ukb25465_ukb25468_extractPhenotypes_imagingSubset.R and $analysis_dir/ukb25465_ukb25468_LM_residuals_imagingSubset.R
#----------------------------------------------------------------------
cp ${UKB_phenos_dir}/${pheno_root}/summary_phenotypes/${pheno_root}_imaging_${covs_name}_${key}_*table ${working_dir2}/pheno/
#----------------------------------------------------------------------
root=ukb_cal_snpQC_${subset_name}_clean_rm025_adj
pheno_root=ukb25465_ukb25468
covs_name=noBioCovs_noAssessmentC
key=Planum_Temporale

#----------------------------------------------------------------------
# Run sex-stratified analyses for PT
#--------------------------------------------
cd ${analysis_dir}/sge_jobs/gcta/

for sex in males females
 do
 type=Volume_${sex}
 pheno_file=$(ls ${working_dir2}/pheno/${pheno_root}_imaging_${covs_name}_${key}_${sex}_Phenotypes_residuals_wHeader.table)
 head -n 1 ${pheno_file} | sed 's/\t/\n/g'> ${key}_${sex}_Phenotypes_residuals.header
 
 for phen in $(grep residuals ${key}_${sex}_Phenotypes_residuals.header | grep ${key})
  do
   phen=$(echo ${phen} | sed 's/residuals_//g')
   # only run if output does not exist
   if [ ! -f ${working_dir2}/reml/${pheno_root}_${covs_name}_reml_${phen}_${type}.hsq ] && [ ! -f ${working_dir}/reml/${pheno_root}_${covs_name}_reml_${phen}_${type}.hsq ]
   then
    echo ${root} ${pheno_root} ${covs_name} ${phen} ${type} ${pheno_file}
    qsub ${analysis_dir}/GCTA_REML_calUKBv2_template.sh ${root} ${pheno_root} ${covs_name} ${phen} ${type} ${pheno_file}
   fi
 done
 rm ${key}_${sex}_Phenotypes_residuals.header
done
 
# correlations between left and right
listLR=${working_dir2}pheno/${pheno_root}_imaging_${covs_name}_nonAI_LeftRight_phenotypes_regions.list

for sex in males females
 do
 type=Volume_${sex}
 pheno_file=$(ls ${working_dir2}/pheno/${pheno_root}_imaging_${covs_name}_${key}_${sex}_Phenotypes_residuals_wHeader.table)
  
 for listLR in $(ls ${working_dir2}pheno/${pheno_root}_imaging_${covs_name}_nonAI_LeftRight_*regions.list )
 do
 while read phen
  do
  grep -v -e ${key} <(echo ${phen})
  if [ $? -ne 0 ]
   then 
   if [ ! -f ${working_dir2}/reml/${root}_${pheno_root}_${covs_name}_reml_covars_${phen}_${type}_left_right_diff1.hsq ] && [ ! -f ${working_dir}/reml/${root}_${pheno_root}_${covs_name}_reml_covars_${phen}_${type}_left_right_diff1.hsq ]
   then
   echo Run bivar REML: ${phen}
   echo ${root} ${pheno_root} ${covs_name} ${phen} ${type} ${pheno_file}
    #qsub -p -10 ${analysis_dir}/GCTA_REMLbivar_calUKBv2_template.sh ${root} ${pheno_root} ${covs_name} ${phen} ${type} ${pheno_file}
   fi
   fi
 done < ${listLR} 
 done
done
#----------------------------------------------------------------------

#--------------------------------------------
# Run bivariate analysis for different phenotypes: PT AI, L, R, TBV
#--------------------------------------------
cd ${analysis_dir}/sge_jobs/gcta/

# already defined by: GCTA_REML_calUKBv2_PT_N18057.sh
grep ${key} ${working_dir2}/pheno/PT_TBV.list > ${working_dir2}/pheno/PT.list
list=${working_dir2}/pheno/PT.list

for sex in males females
 do
 echo Sex: ${sex}
 type=Volume_${sex}
 
 while read phen1
 do
 while read phen2
  do
  if [ ${phen1} != ${phen2} ]
  then
  phenotype=${phen1}_${phen2}_${type}
  rphenotype=${phen2}_${phen1}_${type}
  if [ ! -f ${working_dir2}/reml/${root}_${pheno_root}_reml_covars_${phenotype}_bivar_diff0.hsq ] && [ ! -f ${working_dir}/reml/${root}_${pheno_root}_reml_${covs_name}_${phenotype}_bivar_diff0.hsq ] \
  && [ ! -f ${working_dir2}/reml/${root}_${pheno_root}_reml_covars_${rphenotype}_bivar_diff0.hsq ] && [ ! -f ${working_dir}/reml/${root}_${pheno_root}_reml_${covs_name}_${rphenotype}_bivar_diff0.hsq ] \
  && [ ! -f ${working_dir2}'logs/'${root}'_'${pheno_root}'_reml_'${phenotype}'_bivar_diff0_residuals.log' ] && [ ! -f ${working_dir2}'logs/'${root}'_'${pheno_root}'_reml_'${rphenotype}'_bivar_diff0_residuals.log' ]
  then
   if [[ ${phen1} == *"Planum"* ]]; then type1=_${key}_${sex}; covs_name1=${covs_name} ; else type1=_totalBV_${sex}; covs_name1=${covs_name0};  fi
   pheno_file1=$(ls ${working_dir2}/pheno/${pheno_root}_imaging_${covs_name1}*${type1}*.table)
   if [[ ${phen2} == *"Planum"* ]]; then type2=_${key}_${sex}; covs_name2=${covs_name} ; else type2=_totalBV_${sex}; covs_name2=${covs_name0};  fi
   pheno_file2=$(ls ${working_dir2}/pheno/${pheno_root}_imaging_${covs_name2}*${type2}*.table)
   #
   echo Run bivar REML: ${phenotype}
   #echo ${pheno_root} ${phen1} ${phen2} ${pheno_file1} ${pheno_file2} ${type}
   qsub ${analysis_dir}/GCTA_REMLbivar2_calUKBv2_template.sh ${root} ${pheno_root} ${phen1} ${phen2} ${pheno_file1} ${pheno_file2} ${type}
  fi
  fi
 done < ${list}
 sleep 10s
done < ${list}

done

# rename outputs to reflect the covariates that phenotypes have been adjusted for
rename covars ${covs_name} *males*diff0.hsq
#----------------------------------------------------------------------


#--------------------------------------------
# move output to working_dir, remove from clusterfs
#--------------------------------------------
mv ${working_dir2}/reml/${root}_${pheno_root}_reml_${covs_name}*males* ${working_dir}reml/
# clean all in $working_dir2 # clusterfs
mv ${working_dir2}/reml/*males*table ${working_dir}/reml
mv ${working_dir2}/logs/*reml* ${working_dir}/logs

#--------------------------------------------

#----------------------------------------------------------------------
# organize output
#--------------------------------------------
# extract summary tables for h2 and rho
cd ${working_dir}/reml/

########################################
## h2 ##						########
########################################

# volumes
hsq_file=$(ls -1 *males*hsq | grep $pheno_root | grep -v null | sort | uniq | grep diff -v)
if [ $(echo ${hsq_file} | wc | awk '{print $2}') -ne 0 ]
 then
  paste <(grep "^n" -H ${hsq_file} ) <(grep "/" ${hsq_file} | awk '{print $2,$3}') <(grep "Pval" ${hsq_file} | awk '{print $2}') > hsq_summary_AI_hemisVols_v2cal_${pheno_root}_imaging_${covs_name}_PT_perSex.table
fi

########################################
## rho ##						########
########################################
# extract genetic correlation values
# and with pvalue

# volumes
hsq_file=$(ls -1 *males*_left_right*hsq | grep $pheno_root | grep diff1)
if [ $(echo ${hsq_file} | wc | awk '{print $2}') -ne 0 ]
 then
 paste <(grep "^rG" -H ${hsq_file} ) <(grep "Pval" ${hsq_file} | awk '{print $2}')> hsq_summary_LRrG_pval_hemisVol_v2cal_${pheno_root}_imaging_${covs_name}_PT_perSex.table
fi
#----------------------------------------------------------------------

