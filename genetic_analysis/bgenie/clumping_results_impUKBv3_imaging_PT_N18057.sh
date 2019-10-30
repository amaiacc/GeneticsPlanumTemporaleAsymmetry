#!/bin/bash

#----------------------------------------------------------------------
resource_dir=/data/workspaces/lag/shared_spaces/Resource_DB/
KGPref=${resource_dir}/1KG_phase1_all
ld_ref_eur=${KGPref}/1kg_phase1_eur
#----------------------------------------------------------------------
# define variables
pheno_root=ukb25465_ukb25468
subset_name=imagingT1_N18057
region=Planum_Temporale
type=Volume
# paths to directories
analysis_dir=/data/workspaces/lag/workspaces/lg-ukbiobank/analysis/amaia/genetic_analysis/PT/bgenie/
assoc_dir=/data/workspaces/lag/workspaces/lg-ukbiobank/working_data/amaia/genetic_data/release_v3/imp/bgenie/output/PT/${subset_name}/clean/
working_dir=${assoc_dir}/clumped/
#----------------------------------------------------------------------
mkdir -p ${working_dir}

#----------------------------------------------------------------------
# Use plink for clumping
# follow plink webpage, and fuma as indication
# http://fuma.ctglab.nl/tutorial#snp2gene
# http://zzz.bwh.harvard.edu/plink/clump.shtml
# --clump-p1 0.0001            Significance threshold for index SNPs
# --clump-p2 0.01              Secondary significance threshold for clumped SNPs
# --clump-r2 0.50              LD threshold for clumping
# --clump-kb 250               Physical distance threshold for clumping

cd ${assoc_dir}
for f in $(ls *gz | grep short)
 do
 p=$(echo ${f} | awk -F"_CHR" '{print $1}')
 if [ ! -f ${working_dir}/${p}.assoc ]
 then
  echo Format files and get genome-wide significant SNPs ${p}
  zcat ${f} > ${p}.tmp
  # create new variable as chr:pos, and save file with: rsid, SNP, P
  awk '(NR==1){print "rsid SNP P"}(NR>1){SNP=$1":"$3; print $2,SNP,$10}' ${p}.tmp > ${working_dir}/${p}.assoc
  # get all genome-wide significant SNPs
  head -n 1 ${working_dir}/${p}.assoc > ${working_dir}/${p}.GWS.snps
  awk '$3<5e-8 {print}' ${working_dir}/${p}.assoc >> ${working_dir}/${p}.GWS.snps
  # clean intermediate files
  rm ${p}.tmp
 fi
 
 if [ ! -f ${working_dir}/${p}_clumped.clumped ]
 then
 echo Run clumping for ${p}
  plink --bfile ${ld_ref_eur} \
  --clump ${working_dir}/${p}.assoc \
  --clump-kb 200 --clump-p1 5e-08 --clump-p2 5e-02 --clump-r2 0.4 --noweb --out ${working_dir}/${p}_clumped
 fi
done

# list all genome-wide association signals together, for checkups
cd ${working_dir}

cat *GWS.snps | awk '{print $1,$2}' | sort -r | uniq > all.GWS.snps
cat *.clumped | awk '{print $3}' | sort -r | uniq > pos.lead.snps
# add rsid to the lead snps
awk 'NR==FNR{a[$1]=$1;next}($2) in a{print $0,a[$1]}' pos.lead.snps all.GWS.snps > all.lead.snps
rm pos.lead.snps
