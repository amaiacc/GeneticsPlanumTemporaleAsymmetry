# Extract genotypes for lead SNPs and check whether there is an interaction with sex
## SNPs around the signal already extracted from the bgen file when checking LD: LD_around_leadSNPs_withinUKB.sh
subset_name=imagingT1_N18057

#----------------------------------------------------------------------
# Define paths and directories
#----------------------------------------------------------------------
working_dir=/data/workspaces/lag/workspaces/lg-ukbiobank/working_data/amaia/genetic_data/PT/${subset_name}/lead_snps/genotypes/
ld_dir=/data/workspaces/lag/workspaces/lg-ukbiobank/working_data/amaia/genetic_data/PT/${subset_name}/lead_snps/LD/
assoc_dir=/data/workspaces/lag/workspaces/lg-ukbiobank/working_data/amaia/genetic_data/release_v3/imp/bgenie/output/PT/${subset_name}/clean/clumped/
list=${assoc_dir}/all.lead.snps
#----------------------------------------------------------------------
mkdir -p  ${working_dir}
cd ${working_dir}

# extract genotype data for this snp using plink
# plink files: ukb_imp_chr${chr}_${min}_${max}_v3_imagingT1_N12245
range_kb=2

while read line
 do 
 snp=$(echo ${line} | awk '{print $1}')
 chr=$(echo ${line} | awk '{print $2}' | awk -F":" '{print $1}' )
 pos=$(echo ${line} | awk '{print $2}' | awk -F":" '{print $2}' )
 if [ "$snp" != "rsid" ]
 then
  echo '------------------'  
  echo 'Extract genotypes for: ' ${snp}
  min=$((${pos}-${range_kb}*1000))
  max=$((${pos}+${range_kb}*1000))
  # extract genotype for snp
  plink --bfile ${ld_dir}ukb_imp_chr${chr}_${min}_${max}_v3_${subset_name} --snp ${snp} --recode --out ${working_dir}/ukb_imp_v3_${subset_name}_${snp}
  echo '------------------'  
 fi
 done < ${list}
