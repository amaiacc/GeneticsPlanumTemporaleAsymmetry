# create input to query GTEx v7, using the eQTL calculator
#----------------------------------------------------------------------
subset_name=imagingT1_N18057
lead_snp=rs7420166
#
primary_dir=/data/workspaces/lag/shared_spaces/Resource_DB/OBrien2018_supplementary/
working_dir=/data/workspaces/lag/workspaces/lg-ukbiobank/working_data/amaia/genetic_data/PT/GTEx/
ld_dir=/data/workspaces/lag/workspaces/lg-ukbiobank/working_data/amaia/genetic_data/PT/${subset_name}/lead_snps/LD/
fuma_dir=/data/workspaces/lag/workspaces/lg-ukbiobank/working_data/amaia/genetic_data/PT/${subset_name}/FUMA/FUMA_job36766/
mkdir -p ${working_dir}
cd ${working_dir}
#----------------------------------------------------------------------
transcripts=$(awk -F"," '{print $6}' ${fuma_dir}/blood_FDR_eqtl.csv | sort | grep -v gene | uniq | sed 's/"//g')
ld_files=$(ls ${ld_dir}/*0.6r2*rs7420166*.ld)
snps=$(cat ${ld_files} | awk '{print $6}' | sort | uniq | grep -v SNP_B)
# create input file
touch ${lead_snp}_input4GTEx_eQTL_calculator.txt
for t in ${transcripts}
 do
 for s in ${snps}
  do
  echo ${s}","${t}",Brain_Cortex" >> ${lead_snp}_input4GTEx_eQTL_calculator.txt
  echo ${s}","${t}",Brain_Frontal_Cortex_BA9"  >> ${lead_snp}_input4GTEx_eQTL_calculator.txt
  done
 
 done
 
wc ${lead_snp}_input4GTEx_eQTL_calculator.txt # 264

# (1) Input this content into: https://gtexportal.org/home/testyourown
# (2) Download output as: rs7420166_GTEx_eQTL_calculator.csv
# (3) Run GTEx_chr2_snps.R
