# Check stats and info for all genome-wide significant snps

subset_name=imagingT1_N18057
#----------------------------------------------------------------------
# Define paths and directories
#----------------------------------------------------------------------
primary_dir=/data/workspaces/lag/workspaces/lg-ukbiobank/primary_data/genetic_data/release_v3_March2018/data/mfi/
hrc_info=/data/workspaces/lag/workspaces/lg-ukbiobank/working_data/genetic_data/snp_QC/
snp_stats_dir=/data/workspaces/lag/workspaces/lg-ukbiobank/working_data/amaia/genetic_data/release_v3/imp/subset_${subset_name}/
working_dir=/data/workspaces/lag/workspaces/lg-ukbiobank/working_data/amaia/genetic_data/PT/${subset_name}/lead_snps/snp_stats/
assoc_dir=/data/workspaces/lag/workspaces/lg-ukbiobank/working_data/amaia/genetic_data/release_v3/imp/bgenie/output/PT/${subset_name}/clean/clumped/
#----------------------------------------------------------------------
mkdir -p ${working_dir}


#---------------------------------
# check info from original mfi files for the whole dataset
# header info for mfi files: http://biobank.ctsu.ox.ac.uk/crystal/refer.cgi?id=531
#---------------------------------
# Imputation MAF+info
#===================
#This file lists the minor allele frequency and info score for each of the markers in the imputed data, calculated using QCTOOL. The order of markers in these files is not guaranteed to be the same as the BGEN files.
#Alternate_id
#RS_id
#Position
#Allele1
#Allele2
#Minor Allele
#MAF
#Info score
#---------------------------------
cd ${primary_dir}
echo "Alternate_id RS_id Position Allele1 Allele2 MinorAllele MAF InfoScore" > ${working_dir}/header_mfi.txt
touch ${working_dir}/lookup_v2_NOT_HRCsnps.txt ${working_dir}/lookup_v2_HRCsnps.txt 
#---------------------------------
# define SNP and chr to check

# all genome-wide significant SNPs identified by: clumping_results_impUKBv3_imaging_PT_N18057.sh
list=${assoc_dir}/all.GWS.snps

while read line
 do 
 snp=$(echo ${line} | awk '{print $1}')
 chr=$(echo ${line} | awk '{print $2}' | awk -F":" '{print $1}' )
 pos=$(echo ${line} | awk '{print $2}' | awk -F":" '{print $2}' )
 if [ "$snp" != "rsid" ]
 then
  cd ${primary_dir}
  echo '------------------'  
  echo 'SNP: ' ${snp}
  echo 'Extract MFI info provided by UKB v3'
  grep ${snp} *chr${chr}_v3.txt > ${working_dir}/${snp}_mfi.txt
  echo 'Extract info from QCtools snp-stats'
  cd ${snp_stats_dir}
  grep rsid *chr${chr}_*_snp-stats.txt > ${working_dir}/${snp}_snp-stats.txt
  grep ${snp} *chr${chr}_*_snp-stats.txt >> ${working_dir}/${snp}_snp-stats.txt
  echo 'Check if SNPs are part of HRC or not'
  cd ${hrc_info}
  grep ${snp} ukb_mfi_chr${chr}_v2_NOT_HRCsnps.list >> ${working_dir}/lookup_v2_NOT_HRCsnps.txt
  grep ${snp} ukb_mfi_chr${chr}_v2_HRCsnps.list >> ${working_dir}/lookup_v2_HRCsnps.txt
  echo '------------------'
 fi
 done < ${list}

#---------------------------------
# combine all
#---------------------------------
cd ${working_dir}
cat *snp-stats.txt | sort -n | uniq > ${working_dir}/lookup_snp-stats.txt
cat *_mfi.txt | sort -n | uniq > ${working_dir}/lookup_snp_mfi.txt


