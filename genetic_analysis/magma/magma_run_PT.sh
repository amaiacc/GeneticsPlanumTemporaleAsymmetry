#!/bin/bash
#----------------------------------------------------------------------
# Define paths and directories
#----------------------------------------------------------------------
resource_dir=/data/workspaces/lag/shared_spaces/Resource_DB/
analysis_dir=/data/workspaces/lag/workspaces/lg-ukbiobank/analysis/amaia/genetic_analysis/PT/
magma_dir=${resource_dir}/magma_v1.06b/
#--------
# define variables, again, just for the easy of copy-pasting...
pheno_root=ukb25465_ukb25468
subset_name=imagingT1_N18057
region=Planum_Temporale
type=Volume
#--------
#primary_dir=/data/workspaces/lag/workspaces/lg-ukbiobank/working_data/amaia/genetic_data/release_v3/imp/bgenie/output/PT/clean/
#working_dir=/data/workspaces/lag/workspaces/lg-ukbiobank/working_data/amaia/genetic_data/PT/magma/
primary_dir=/data/workspaces/lag/workspaces/lg-ukbiobank/working_data/amaia/genetic_data/release_v3/imp/bgenie/output/PT/${subset_name}/clean/
working_dir=/data/workspaces/lag/workspaces/lg-ukbiobank/working_data/amaia/genetic_data/PT/${subset_name}/magma/
#--------

# define pheno root, for different rounds of analysis
#pheno_root=ukb21288_ukb21293
pheno_root=ukb25465_ukb25468
#----------------------------------------------------------------------
mkdir -p ${working_dir}
# select phenotypes
cd ${primary_dir}

key=AI_VOLUME_Planum_Temporale
phens=$(ls *${key}*1e-07HWEp_0.7INFO_0.001MAF.txt.gz | grep ${key} )
# edit template to create jobs per phenotype
mkdir -p ${analysis_dir}/sge_jobs/magma
cd ${analysis_dir}/sge_jobs/magma
up_win=35
down_win=10

for p in ${phens}
 do
 phen=$(echo $p | sed "s/${pheno_root}//g" | sed 's/_CHRall_1e-07HWEp_0.7INFO_0.001MAF.txt.gz//g' | sed 's/_$//g')
 # check that output dir does not exist; alternatively could condition submission of job to existence of output file... but since there are many...
 if [ ! -d ${working_dir}/${phen}_${pheno_root}/${up_win}_${down_win}kb.log ]
 then
  echo ${pheno_root} ${phen}
  qsub ${analysis_dir}/magma_template.sh ${pheno_root} ${phen} ${up_win} ${down_win} ${subset_name}
 fi
 done

 
