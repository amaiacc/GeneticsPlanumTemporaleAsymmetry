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

#----------------------------------------------------------------------
# define root that will be contained within generated input and output files; parameter to change if input files change
pheno_root=ukb25465_ukb25468
covs_name=noBioCovs_noAssessmentC_TBV
covs_name00=noBioCovs_noAssessmentC_TBV # including TBV as covariate in lm, from which residuals were taken from
covs_name0=noBioCovs_noAssessmentC
subset_name=imagingT1_N18057
region=Planum_Temporale
root=ukb_cal_snpQC_${subset_name}_clean_rm025_adj
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# copy the data (GRMs) to sge2 directory first, otherwise it won't run in the grid
mkdir -p ${working_dir2}/grm/ ${working_dir2}/reml ${working_dir2}/pheno ${working_dir2}/logs
mkdir -p ${working_dir}/grm/ ${working_dir}/reml ${working_dir}/pheno ${working_dir}/logs
# copy genetic data to working_dir2 in sge, if not present yet!
if [ ! -f ${working_dir2}/grm/ukb_cal_snpQC_${subset_name}_clean_rm025_adj.grm.bin ]
 then 
 echo Copy GRM files to clusterfs
 cp ${working_dir}/grm/ukb_cal_snpQC_${subset_name}_clean_rm025_adj.grm* ${working_dir2}/grm/
 fi

#----------------------------------------------------------------------
# Define phenotypes and covariate files:
#----------------------------------------------------------------------
pheno_file1=${UKB_phenos_dir}/${pheno_root}/summary_phenotypes/${pheno_root}_imaging_${covs_name}_PTPhenotypes_residuals_wHeader.table
pheno_file2=${UKB_phenos_dir}/${pheno_root}/summary_phenotypes/${pheno_root}_imaging_${covs_name0}_totalBVPhenotypes_residuals_wHeader.table

cp ${pheno_file1} ${working_dir2}/pheno/
cp ${pheno_file2} ${working_dir2}/pheno/


# Define list of phenotypes, to run!
head -n 1 $pheno_file1 | sed 's/\t/\n/g' | grep -e residuals | grep -e Planum_Temporale > ${working_dir2}/pheno/PT_TBV.list
head -n 1 $pheno_file2 | sed 's/\t/\n/g' | grep -e residuals | grep -e totalBV >> ${working_dir2}/pheno/PT_TBV.list
sed -i 's/residuals_//g' ${working_dir2}/pheno/PT_TBV.list
list=${working_dir2}/pheno/PT_TBV.list

#--------------------------------------------
# Run REML for phenotypes within list 
## scripts created within: GCTA_REML_calUKBv2_IDPs_2.sh
#--------------------------------------------
# run loop, or first just one to test
cd ${working_dir2}/sge_jobs/gcta/

#--------------------------------------------
# h2
#--------------------------------------------
while read phen
 do
  if [[ ${phen} == *"Planum"* ]]; then
   type=PT
   covs_name=${covs_name00}
  else
   type=_totalBV
   covs_name=${covs_name0}
  fi
  pheno_file=$(ls ${working_dir2}/pheno/${pheno_root}_imaging_${covs_name}*${type}*.table)
  if [ ! -f ${working_dir2}/reml/${pheno_root}_${covs_name}_reml_${phen}_${type}.hsq ] && [ ! -f ${working_dir}/reml/${pheno_root}_${covs_name}_reml_${phen}_${type}.hsq ]
  then
   echo Submit job for: ${root} ${pheno_root} ${covs_name} ${phen} ${type} ${pheno_file}
   qsub ${analysis_dir}/GCTA_REML_calUKBv2_template.sh ${root} ${pheno_root} ${covs_name} ${phen} ${type} ${pheno_file}
  fi
  done < ${list}

#--------------------------------------------
# Run for different phenotypes
#--------------------------------------------
# run loop, or first just one to test
#cd ${analysis_dir}/sge_jobs/gcta/
cd ${working_dir2}/sge_jobs/gcta/

type=adjTBV

while read phen1
do
 while read phen2
  do
  if [ ${phen1} != ${phen2} ]
  then
  phenotype=${phen1}_${phen2}_${type}
  rphenotype=${phen2}_${phen1}_${type}
  if [ ! -f ${working_dir2}/reml/${root}_${pheno_root}_reml_covars_${phenotype}_bivar_diff0.hsq ] && [ ! -f ${working_dir}/reml/${root}_${pheno_root}_reml_${covs_name00}_${phenotype}_bivar_diff0.hsq ] \
  && [ ! -f ${working_dir2}/reml/${root}_${pheno_root}_reml_covars_${rphenotype}_bivar_diff0.hsq ] && [ ! -f ${working_dir}/reml/${root}_${pheno_root}_reml_${covs_name00}_${rphenotype}_bivar_diff0.hsq ] \
  && [ ! -f ${working_dir2}'logs/'${root}'_'${pheno_root}'_reml_'${phenotype}'_bivar_diff0_residuals.log' ] && [ ! -f ${working_dir2}'logs/'${root}'_'${pheno_root}'_reml_'${rphenotype}'_bivar_diff0_residuals.log' ]
  then
   if [[ ${phen1} == *"Planum"* ]]; then type1=_PT; covs_name1=${covs_name00} ; else type1=_totalBV; covs_name1=${covs_name0};  fi
   pheno_file1=$(ls ${working_dir2}/pheno/${pheno_root}_imaging_${covs_name1}*${type1}*.table)
   if [[ ${phen2} == *"Planum"* ]]; then type2=_PT; covs_name2=${covs_name00} ; else type2=_totalBV; covs_name2=${covs_name0};  fi
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

# rename outputs to reflect that phenotypes have been adjTBV
rename covars ${covs_name00} *${type}*diff0.hsq

#--------------------------------------------
# move output to working_dir, remove from clusterfs
#--------------------------------------------
mv ${working_dir2}/reml/*adjTBV*diff0.hsq ${working_dir}reml/
mv ${working_dir2}/reml/*adjTBV*.hsq ${working_dir}reml/
mv ${working_dir2}/reml/*totalBV*.hsq ${working_dir}reml/
#--------------------------------------------


#--------------------------------------------
# organize output
#--------------------------------------------
# extract summary tables for h2 and rho
cd ${working_dir}/reml/
pheno_root=ukb25465_ukb25468
covs_name00=noBioCovs_noAssessmentC_TBV # including TBV as covariate in lm, from which residuals were taken from
covs_name0=noBioCovs_noAssessmentC # without including biological covariates
subset_name=imagingT1_N18057
region=Planum_Temporale
root=ukb_cal_snpQC_${subset_name}_clean_rm025_adj
#
covs_name=${covs_name00}

########################################
## h2 ##						########
########################################

# totalBV
hsq_file=$(ls -1 *.hsq | grep -e BV -e grey_white | grep ${pheno_root} | grep -v diff )
paste <(grep "^n" -H ${hsq_file} ) <(grep "/" ${hsq_file} | awk '{print $2,$3}') <(grep "Pval" ${hsq_file} | awk '{print $2}') > hsq_summary_TBV_v2cal_${pheno_root}_imaging_${covs_name}.table

########################################
## rho ##						########
########################################
# extract genetic correlation values
# and with pvalue

# PT and totalBV
hsq_file=$(ls -1 *${root}_${pheno_root}_reml_${covs_name00}_*bivar_diff0.hsq | grep -e totalBV -e grey_white -e Planum_Temporale)
paste <(grep "^rG" -H ${hsq_file} ) <(grep "Pval" ${hsq_file} | awk '{print $2}')> hsq_summary_PT_TBV_rG_pval_v2cal_${pheno_root}_imaging_${covs_name}.table

#--------------------------------------------
# clean all in $working_dir2 # clusterfs
#--------------------------------------------
mv ${working_dir2}/reml/*table ${working_dir}/reml/
mv ${working_dir2}/logs/*reml* ${working_dir}/logs/

