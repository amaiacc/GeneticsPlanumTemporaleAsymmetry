#!/bin/bash

# Check LD structure around rs41298373 and other lead SNPs, within the UKB (samples that were used for the analyses)

subset_name=imagingT1_N18057
#----------------------------------------------------------------------
# Define paths and directories
#----------------------------------------------------------------------
PATH=/home/winhomes/amacar/programs/HaploView/:$PATH
#PATH=/data/workspaces/lag/shared_spaces/Resource_DB/:$PATH
# use qctool within /usr/local/bin/qctool
PATH=/data/workspaces/lag/shared_spaces/Resource_DB/src/tabix/:${PATH}
haploview_dir=/home/winhomes/amacar/programs/HaploView/
resource_dir=/data/workspaces/lag/shared_spaces/Resource_DB/
analysis_dir=/data/workspaces/lag/workspaces/lg-ukbiobank/analysis/amaia/genetic_analysis/PT/${subset_name}/
UKB=/data/workspaces/lag/workspaces/lg-ukbiobank/working_data/amaia/genetic_data/
UKB_cal=${UKB}release_v2/cal/QC/
UKB_imp=${UKB}release_v3/imp/subset_${subset_name}/
assoc_dir=/data/workspaces/lag/workspaces/lg-ukbiobank/working_data/amaia/genetic_data/release_v3/imp/bgenie/output/PT/${subset_name}/clean/clumped/
working_dir=/data/workspaces/lag/workspaces/lg-ukbiobank/working_data/amaia/genetic_data/PT/${subset_name}/lead_snps/LD/

mkdir -p ${working_dir}
cd ${working_dir}

# lead SNPs identified by: clumping_results_impUKBv3_imaging_PT_N18057.sh
list=${assoc_dir}/all.lead.snps

# as in FUMA: http://fuma.ctglab.nl/tutorial#refpanel
# --r2 --ld-window 99999 --ld-window-r2 0.05

while read line
 do  
 snp=$(echo ${line} | awk '{print $1}')
 chr=$(echo ${line} | awk '{print $2}' | awk -F":" '{print $1}' )
 pos=$(echo ${line} | awk '{print $2}' | awk -F":" '{print $2}' )
 if [ "$snp" != "rsid" ]
 then
 echo Check if ${snp} was genotyped
 grep ${snp} ${UKB_cal}ukb_cal_snpQC_${subset_name}_clean.bim
 for range_kb in 200 2 20 50
  do  
   if [ ! -f ${snp}_imp_${range_kb}kb_HVfile.info ]
   then    
    echo '------------------'    
    min=$((${pos}-${range_kb}*1000))
    max=$((${pos}+${range_kb}*1000))
    echo '------------------'    
    if [ ${#chr} == 2 ];	then range=${chr}:${min}-${max}; else  range=0${chr}:${min}-${max}; fi
    echo Range that will be extracted from bgen file: ${range}
    echo Convert to plink files using qctool
    qctool -g ${UKB_imp}ukb_imp_chr${chr}_v3_${subset_name}.bgen -s ${UKB_imp}ukb_imp_chr${chr}_v3_${subset_name}.sample -missing-code 0     -threshold 0.9     -incl-range ${range}     -ofiletype binary_ped -og ukb_imp_chr${chr}_${min}_${max}_v3_${subset_name}
    echo '------------------'
    plink --bfile ukb_imp_chr${chr}_${min}_${max}_v3_${subset_name} \
     --r2 --ld-snp ${snp} \
     --ld-window-kb 1000 --ld-window 99999 --ld-window-r2 0 \
     --out LD_${r2}r2_${window}kb_${snp}_imp
    echo Recode in plink to 12 and convert to HaploView format
    plink -bfile ukb_imp_chr${chr}_${min}_${max}_v3_${subset_name} --recode 12 --make-bed --out ukb_imp_chr${chr}_${min}_${max}_v3_${subset_name}_recode12
    plink -file ukb_imp_chr${chr}_${min}_${max}_v3_${subset_name}_recode12 --recodeHV --out ${snp}_imp_${range_kb}kb_HVfile
    echo '------------------'
   fi
  done
  fi
 done < ${list}

 
 # use haploview to create figure:
 ## in windows, java memory problems depending on the ped file size
 ## format: linkage format
 ## data file: *HVfile.ped
 ## locus information file: *HVfile.info
 
# for all lead SNPs: get list of SNPs with r2>0.6
for f in $(ls LD_r2_kb_*_imp.ld)
 do
 f2=$(echo ${f} | sed 's/r2/0.6r2/g')
 awk '{ if ($7 >= 0.6) print $0 }' ${f} > ${f2}
 done
