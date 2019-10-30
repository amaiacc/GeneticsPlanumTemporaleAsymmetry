#!/bin/bash

resource_dir=/data/workspaces/lag/shared_spaces/Resource_DB/
ldscores_dir=${resource_dir}/LDscores/
## info on LD score formats/ estimation
# https://github.com/bulik/ldsc/wiki/LD-File-Formats
# https://github.com/bulik/ldsc/wiki/LD-Score-Estimation-Tutorial
ldscores_dir2=${working_dir2}/LDscores/

#--------------------------
# define variables, again, just for the easy of copy-pasting...
subset_name=imagingT1_N18057
working_dir=/data/workspaces/lag/workspaces/lg-ukbiobank/working_data/amaia/genetic_data/PT/${subset_name}/ldsc/
analysis_dir=/data/workspaces/lag/workspaces/lg-ukbiobank/analysis/amaia/genetic_analysis/PT/${subset_name}/ldsc/
working_dir2=/data/clusterfs/lag/users/amacar/ukb/
#--------------------------


#--------------------------
#cts_name=Cahoy
#cts_name=GTEx_brain
cts_name=Multi_tissue_gene_expr
# copy ldscores_dir annotations to clusterfs
if [ ! -d ${working_dir2}/LDscores/Phase3/ ]
then
 mkdir -p ${working_dir2}/LDscores/Phase3/
 cd ${working_dir2}/LDscores/Phase3/
 cp -R ${ldscores_dir}/Phase3/1000G_Phase3_weights_hm3_no_MHC ./
 cp -R ${ldscores_dir}/Phase3/1000G_Phase3_frq ./
 cp -R ${ldscores_dir}/Phase3/baseline_v1.1 ./
 #cp -R ${ldscores_dir}/Phase3/baselineLD_v1.1/ ./
 cp -R ${ldscores_dir}/Phase3/baselineLD_v2.2/ ./
fi
#--------
if [ ! -f ${working_dir2}/LDscores/Phase3/${cts_name}.ldcts ]
then
 cd ${working_dir2}/LDscores/Phase3/
 cp -R ${ldscores_dir}/Phase3/${cts_name}.ldcts ./
 cp -R ${ldscores_dir}/Phase3/${cts_name}_1000Gv3_ldscores/ ./
fi
#--------------------------

#--------------------------
# Run partitioned heritability
#--------------------------
# sumstats generated by : LDSC_GenCor_GWASes_bgenie_PT*_N18057.sh
## sumstats format: https://github.com/bulik/ldsc/wiki/Summary-Statistics-File-Format

## use the grid to run it
#--------
# define phenotypes
cd ${working_dir}
phenos=$(ls *sumstats.gz | sed 's/.sumstats.gz//g')
#
mkdir -p ${working_dir2}/ldsc/
cd ${working_dir2}/ldsc/
cp ${analysis_dir}LDSC_h2_partitioned_grid.sh ./

for p in ${phenos}
 do
 #if [ ! -f ${working_dir}/partitioned/${p}_h2_baselineLD.log ]
 #then
 # if [ ! -f ${p}_h2_baselineLD.log ]
 # then
  if [ ! -f ${p}.sumstats.gz ]
   then
   echo Copy sumstats
   cp ${working_dir}${p}.sumstats.gz ./
  fi
  echo Run ${p}
  #qsub LDSC_h2_partitioned_grid.sh ${p} ${working_dir2}/ldsc/ ${working_dir2}/LDscores/ ${cts_name}
  bash LDSC_h2_partitioned_grid.sh ${p} ${working_dir2}/ldsc/ ${working_dir2}/LDscores/ ${cts_name}
 #fi
 #fi
 done
#--------
# move results to workspace
mkdir -p ${working_dir}/partitioned/
mv *Planum_Temporale*log ${working_dir}/partitioned/
mv *Planum_Temporale*results* ${working_dir}/partitioned/
#--------

#--------
# Clean
#--------
# remove sumstats
rm *Planum_Temporale*sumstats.gz 
rm *.o*
rm *.e*
rm *.sh
# clean also the LDscores
cd ${working_dir2}/LDscores/Phase3/
rm 1000G_Phase3_weights_hm3_no_MHC/*
rm 1000G_Phase3_frq/*
rm baseline_v1.1/*
rm baselineLD_v1.1/*
rmdir 1000G_Phase3_weights_hm3_no_MHC 1000G_Phase3_frq baseline_v1.1 baselineLD_v1.1
# remove directories
cd ${working_dir2}
rmdir LDscores/Phase3/
rmdir ldsc 



#----------------------------------------------------
# ** work in progress **
#----------------------------------------------------
# for cell-type specific analyses

  if [ ! -f ${p}_h2_baselineLD.log ]
   then
   #---------------
    echo '-----------> Baseline LD model (Finucane et al. 2018)'
    ldsc.py \
    --h2 ${p}.sumstats.gz \
    --ref-ld-chr ${ldscores_dir}/Phase3/baselineLD_v1.1/baselineLD. \
    --w-ld-chr ${ldscores_dir}/Phase3/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
    --overlap-annot \
    --frqfile-chr ${ldscores_dir}/Phase3/1000G_Phase3_frq/1000G.EUR.QC. \
    --print-coefficients \
    --out ${p}_h2_baselineLD
   #---------------
  fi

