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
cp ${UKB_phenos_dir}/${pheno_root}/summary_phenotypes/${pheno_root}_imaging_${covs_name}*table ${working_dir2}/pheno/

#----------------------------------------------------------------------
# Define list of phenotypes, to run!
## I think this is not necessary, only if want to run many phenotypes... but well
#--------------------------------------------
grep -e left -e right ${UKB_phenos_dir}${pheno_root}_imaging_${covs_name}_nonAI_VolumePhenotypes.list > ${UKB_phenos_dir}${pheno_root}_imaging_${covs_name}_nonAI_LeftRight_Volumephenotypes.list
#grep -e left -e right ${UKB_phenos_dir}${pheno_root}_imaging_${covs_name}_nonAI_FAskeletonPhenotypes.list > ${UKB_phenos_dir}${pheno_root}_imaging_${covs_name}_nonAI_LeftRight_FAskeleton.list
#grep -e left -e right ${UKB_phenos_dir}${pheno_root}_imaging_${covs_name}_nonAI_tractsPhenotypes.list > ${UKB_phenos_dir}${pheno_root}_imaging_${covs_name}_nonAI_LeftRight_tracts.list

# create lists of non-lateral measures
grep -e left -e right -v ${UKB_phenos_dir}${pheno_root}_imaging_${covs_name}_nonAI_VolumePhenotypes.list > ${UKB_phenos_dir}${pheno_root}_imaging_${covs_name}_nonLateral_Volumephenotypes.list
#grep -e left -e right -v ${UKB_phenos_dir}${pheno_root}_imaging_${covs_name}_nonAI_FAskeletonPhenotypes.list > ${UKB_phenos_dir}${pheno_root}_imaging_${covs_name}_nonLateral_FAskeleton.list
#grep -e left -e right -v ${UKB_phenos_dir}${pheno_root}_imaging_${covs_name}_nonAI_tractsPhenotypes.list > ${UKB_phenos_dir}${pheno_root}_imaging_${covs_name}_nonLateral_tracts.list

all_lists=$(echo ${pheno_root}_imaging_${covs_name}_AI_VolumePhenotypes.list \
 ${pheno_root}_imaging_${covs_name}_nonAI_LeftRight_Volumephenotypes.list ${pheno_root}_imaging_${covs_name}_nonLateral_Volumephenotypes.list \
 ${pheno_root}_imaging_${covs_name}_AI_FAskeletonPhenotypes.list \
 ${pheno_root}_imaging_${covs_name}_nonAI_LeftRight_FAskeleton.list ${pheno_root}_imaging_${covs_name}_nonLateral_FAskeleton.list \
 ${pheno_root}_imaging_${covs_name}_AI_tractsPhenotypes.list ${pheno_root}_imaging_${covs_name}_nonAI_LeftRight_tracts.list \
 ${pheno_root}_imaging_${covs_name}_nonLateral_tracts.list \
 )

# only vols
all_lists=$(echo ${pheno_root}_imaging_${covs_name}_AI_VolumePhenotypes.list ${pheno_root}_imaging_${covs_name}_nonAI_LeftRight_Volumephenotypes.list)

# for every region with left and right, create input pheno files to run bivariate analyses reml
# get region list without hemisphere info
sed 's/_left//g' ${UKB_phenos_dir}${pheno_root}_imaging_${covs_name}_nonAI_LeftRight_Volumephenotypes.list | sed 's/_right//g' | grep -v AI | sort | uniq > ${UKB_phenos_dir}${pheno_root}_imaging_${covs_name}_nonAI_LeftRight_Volumephenotypes_regions.list # 72
#sed 's/_left//g' ${UKB_phenos_dir}${pheno_root}_imaging_${covs_name}_nonAI_LeftRight_FAskeleton.list | sed 's/_right//g' | grep -v AI | sort | uniq > ${UKB_phenos_dir}${pheno_root}_imaging_${covs_name}_nonAI_LeftRight_FAskeleton_regions.list # 189
#sed 's/_left//g' ${UKB_phenos_dir}${pheno_root}_imaging_${covs_name}_nonAI_LeftRight_tracts.list | sed 's/_right//g' | grep -v AI | sort | uniq > ${UKB_phenos_dir}${pheno_root}_imaging_${covs_name}_nonAI_LeftRight_tracts_regions.list # 108
# how many of each type
wc ${UKB_phenos_dir}${pheno_root}_imaging_${covs_name}_nonAI_LeftRight_*regions.list
cp ${UKB_phenos_dir}${pheno_root}_imaging_${covs_name}_nonAI_LeftRight_*regions.list ${working_dir2}pheno/ # 72, 189, 108
#----------------------------------------------------------------------


#----------------------------------------------------------------------
# Run REML for phenotypes within list / key phenotype
#----------------------------------------------------------------------
# run loop, or first just one to test
cd ${analysis_dir}/sge_jobs/gcta/

root=ukb_cal_snpQC_${subset_name}_clean_rm025_adj
pheno_root=ukb25465_ukb25468
covs_name=noBioCovs_noAssessmentC
key=Planum_Temporale

for list_file in ${all_lists}
 do
  list=${UKB_phenos_dir}${list_file}
  type=$(echo $list_file | sed "s/${pheno_root}_imaging_${covs_name}_//g" | sed 's/.list//g' | sed 's/LeftRight_//g' | awk -F'_' '{print $2}' | sed 's/Phenotypes//g' | sed 's/phenotypes//g' )
  pheno_file=$(ls ${working_dir2}/pheno/${pheno_root}_imaging_${covs_name}*${type}*residuals_wHeader.table)

  #while read phen
  for phen in $(grep ${key} ${list})
   do
    #only submit job if output does not exist, neither in working_dir nor working_dir2
    if [ ! -f ${working_dir2}/reml/${pheno_root}_${covs_name}_reml_${phen}_${type}.hsq ] && [ ! -f ${working_dir}/reml/${pheno_root}_${covs_name}_reml_${phen}_${type}.hsq ]
     then
     echo ${phen}
      qsub ${analysis_dir}/GCTA_REML_calUKBv2_template.sh ${root} ${pheno_root} ${covs_name} ${phen} ${type} ${pheno_file}
    fi
  done
  #done < ${list}
 done

#--------------------------------------------
# Bivariate analysis, gen cor Left and Right, diff to 1
#--------------------------------------------
cd ${analysis_dir}/sge_jobs/gcta/
listLR=${working_dir2}pheno/${pheno_root}_imaging_${covs_name}_nonAI_LeftRight_phenotypes_regions.list
pheno_file=$(ls ${working_dir2}/pheno/${pheno_root}_imaging_${covs_name}*${type}*residuals_wHeader.table)

for listLR in $(ls ${working_dir2}pheno/${pheno_root}_imaging_${covs_name}_nonAI_LeftRight_*regions.list )
 do
 while read phen
  do
  grep -v -e ${key} <(echo ${phen})
  if [ $? -ne 0 ]
   then 
   if [ ! -f ${working_dir2}/reml/${root}_${pheno_root}_reml_covars_${phen}_left_right_diff1.hsq ] && [ ! -f ${working_dir}/reml/${root}_${pheno_root}_reml_covars_${phen}_left_right_diff1.hsq ]
   then
   echo Run bivar REML: ${phen}
    #qsub -p -10 ${analysis_dir}/GCTA_REMLbivar_calUKBv2_template.sh ${root} ${pheno_root} ${covs_name} ${phen} ${type} ${pheno_file}
   fi
   fi
 done < ${listLR} 
done

#--------------------------------------------
# Run bivariate analysis for different phenotypes: PT AI, L, R, TBV
#--------------------------------------------
cd ${analysis_dir}/sge_jobs/gcta/

#pheno_file1=${working_dir2}/pheno/${pheno_root}_imaging_${covs_name}_Planum_Temporale_all_Phenotypes_residuals_wHeader.table # should be the same
pheno_file1=${working_dir2}/pheno/${pheno_root}_imaging_${covs_name}VolumePhenotypes_residuals_wHeader.table
pheno_file2=${working_dir2}/pheno/${pheno_root}_imaging_${covs_name}_totalBVPhenotypes_residuals_wHeader.table

# Define list of phenotypes, to run!
head -n 1 $pheno_file1 | sed 's/\t/\n/g' | grep -e residuals | grep -e Planum_Temporale > ${working_dir2}/pheno/PT_TBV.list
head -n 1 $pheno_file2 | sed 's/\t/\n/g' | grep -e residuals | grep -e totalBV >> ${working_dir2}/pheno/PT_TBV.list
sed -i 's/residuals_//g' ${working_dir2}/pheno/PT_TBV.list
list=${working_dir2}/pheno/PT_TBV.list

while read phen1
do
 while read phen2
  do
  if [ ${phen1} != ${phen2} ]
  then
   if [[ ${phen1} == *"Planum"* ]]; then type1=Volume; else type1=totalBV;  fi
   pheno_file1=$(ls ${working_dir2}/pheno/${pheno_root}_imaging_${covs_name}*${type1}*.table)
   if [[ ${phen2} == *"Planum"* ]]; then type2=Volume; else type2=totalBV;  fi
   pheno_file2=$(ls ${working_dir2}/pheno/${pheno_root}_imaging_${covs_name}*${type2}*.table)
   #
   phenotype=${phen1}_${phen2}_${type1}
   rphenotype=${phen2}_${phen1}_${type2}
   if [ ! -f ${working_dir2}/reml/${root}_${pheno_root}_reml_covars_${phenotype}_bivar_diff0.hsq ] && [ ! -f ${working_dir}/reml/${root}_${pheno_root}_reml_${covs_name}_${phenotype}_bivar_diff0.hsq ] \
   && [ ! -f ${working_dir2}/reml/${root}_${pheno_root}_reml_covars_${rphenotype}_bivar_diff0.hsq ] && [ ! -f ${working_dir}/reml/${root}_${pheno_root}_reml_${covs_name}_${rphenotype}_bivar_diff0.hsq ] \
   && [ ! -f ${working_dir2}'logs/'${root}'_'${pheno_root}'_reml_'${phenotype}'_bivar_diff0_residuals.log' ] && [ ! -f ${working_dir2}'logs/'${root}'_'${pheno_root}'_reml_'${rphenotype}'_bivar_diff0_residuals.log' ]
   then
   echo Run bivar REML: ${phenotype}
   echo ${pheno_file1}
   echo ${pheno_file2}
   #echo ${pheno_root} ${phen1} ${phen2} ${pheno_file1} ${pheno_file2} ${type1}
   qsub ${analysis_dir}/GCTA_REMLbivar2_calUKBv2_template.sh ${root} ${pheno_root} ${phen1} ${phen2} ${pheno_file1} ${pheno_file2} ${type1}
  fi
  fi
 done < ${list}
 sleep 10s
done < ${list}

# rename outputs to reflect the covariates that phenotypes have been adjusted for
rename covars ${covs_name} *diff0.hsq
#----------------------------------------------------------------------

#--------------------------------------------
# move output to working_dir, remove from clusterfs
#--------------------------------------------
mv ${working_dir2}/reml/*${pheno_root}*${covs_name}* ${working_dir}reml/
#--------------------------------------------
# clean all in $working_dir2 # clusterfs
#--------------------------------------------
cd ${working_dir2}
mv ${working_dir2}/reml/*table ${working_dir}/reml
mv ${working_dir2}/logs/*reml* ${working_dir}/logs


#--------------------------------------------
# organize output
#--------------------------------------------
# extract summary tables for h2 and rho
cd ${working_dir}/reml/

########################################
## h2 ##						########
########################################

# volumes
hsq_file=$(ls -1 *AI_VOLUME*hsq *Volume*hsq | grep $pheno_root | grep -v null | sort | uniq | grep diff -v)
if [ $(echo ${hsq_file} | wc | awk '{print $2}') -ne 0 ]
 then
  paste <(grep "^n" -H ${hsq_file} ) <(grep "/" ${hsq_file} | awk '{print $2,$3}') <(grep "Pval" ${hsq_file} | awk '{print $2}') > hsq_summary_AI_hemisVols_v2cal_${pheno_root}_imaging_${covs_name}.table
fi

########################################
## rho ##						########
########################################
# extract genetic correlation values
# and with pvalue

# volumes
hsq_file=$(ls -1 *${pheno_root}*Volume*_left_right*hsq | grep $pheno_root | grep diff1)
if [ $(echo ${hsq_file} | wc | awk '{print $2}') -ne 0 ]
 then
 paste <(grep "^rG" -H ${hsq_file} ) <(grep "Pval" ${hsq_file} | awk '{print $2}')> hsq_summary_LRrG_pval_hemisVol_v2cal_${pheno_root}_imaging_${covs_name}.table
fi

# PT and totalBV
hsq_file=$(ls -1 *${pheno_root}_reml_${covs_name}_*bivar_diff0.hsq | grep -v adjTBV | grep -e totalBV -e grey_white -e Planum_Temporale)
paste <(grep "^rG" -H ${hsq_file} ) <(grep "Pval" ${hsq_file} | awk '{print $2}')> hsq_summary_PT_TBV_rG_pval_v2cal_${pheno_root}_imaging_${covs_name}.table

