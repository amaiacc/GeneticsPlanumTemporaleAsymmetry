#!/bin/bash

PATH=/usr/local/apps/ldsc/:${PATH}
resource_dir=/data/workspaces/lag/shared_spaces/Resource_DB/
#ukb_neale=/data/workspaces/lag/workspaces/lg-ukbiobank/working_data/amaia/downloads/NealeLab/assoc/
ldscores_dir=${resource_dir}/LDscores/
## info on LD score formats/ estimation
# https://github.com/bulik/ldsc/wiki/LD-File-Formats
# https://github.com/bulik/ldsc/wiki/LD-Score-Estimation-Tutorial
ldscores_dir=${resource_dir}/LDscores/
#

## sumstats format: https://github.com/bulik/ldsc/wiki/Summary-Statistics-File-Format
# required columns:
##    SNP -- SNP identifier (e.g., rs number)
##  N -- sample size (which may vary from SNP to SNP).
##  Z -- z-score. Sign with respect to A1 (warning, possible gotcha)
##  A1 -- first allele (effect allele)
##  A2-- second allele (other allele)
# Note that ldsc filters out all variants that are not SNPs and strand-ambiguous SNPs.

#--------
# define variables, again, just for the easy of copy-pasting...
pheno_root=ukb25465_ukb25468
subset_name=imagingT1_N18057
region=Planum_Temporale
key=totalBV
type=Volume
#--------

assoc_dir=/data/workspaces/lag/workspaces/lg-ukbiobank/working_data/amaia/genetic_data/release_v3/imp/bgenie/output/PT/${subset_name}/clean/
working_dir=/data/workspaces/lag/workspaces/lg-ukbiobank/working_data/amaia/genetic_data/PT/${subset_name}/ldsc/
#--------
pheno_file=/data/workspaces/lag/workspaces/lg-ukbiobank/working_data/amaia/genetic_data/release_v3/imp/bgenie/${pheno_root}_${type}_sample_${subset_name}_${region}_phenos4BGENIE.table
mkdir -p ${working_dir}
cd ${working_dir}
#--------
head -n 1 ${pheno_file} | sed 's/ /\n/g' > ${working_dir}pheno_header.table
# define phenotypes:
grep ${region} pheno_header.table
phenos=$(grep ${region} pheno_header.table | grep -v ukb | sed 's/residuals_//g')
phensLR=$(echo $phenos | sed "s/ /\\n/g" | grep -v AI | sed 's/_left//g' | sed 's/_right//g' |  sort | uniq)


# BGENIE: 'the regression model we code the first and second alleles as 0 and 1 respectively, so the beta coefficient refers to the effect of having an extra copy of the second allele.'
# a_0 The character code for the first allele (a string).
# a_1 The character code for the second allele (a string). --> effect allele
## https://jmarchini.org/bgenie-usage/

# ignore beta column, otherwise this error comes up for totalBV
## https://groups.google.com/forum/#!topic/ldsc_users/RLbVw3e_PU0


for p in ${phenos}
 do
 echo ${p}
 if [ ! -f ${p}.adjTBV.sumstats.gz ]
 then
  p1=${region}_${key}
  p2=$(echo $p | sed "s/$region/$p1/g")
  p_file=${assoc_dir}${p2}_CHRall_1e-07HWEp_0.7INFO_0.001MAF.txt.gz
  # define N, should get it from pheno file...
  n=$(grep -nr residuals_$p$ pheno_header.table | awk -F':' '{print $1}')
  N1=$(awk -v c=${n} '{print $c}' ${pheno_file} | grep -v -e '\-999' -e ${p} | wc -l) # need to get the right number!
  
  if [ -f ${p_file} ]
   then
   echo 'Run munge_sumstats to create input files for LDSC - ' ${p}
   munge_sumstats.py \
    --sumstats ${p_file} \
    --N ${N1} \
    --p P \
    --snp rsid \
    --info INFO.UKB \
    --out ${p}.adjTBV \
    --a1 a_1 --a2 a_0 \
    --merge-alleles ${ldscores_dir}/w_hm3.snplist
  #---------------
    echo 'Estimate heritability for' ${p}
    ldsc.py \
    --h2 ${p}.adjTBV.sumstats.gz \
    --ref-ld-chr ${ldscores_dir}/eur_w_ld_chr/ \
    --w-ld-chr ${ldscores_dir}/eur_w_ld_chr/ \
    --out ${p}.adjTBV_h2
   #---------------
   fi
  fi
 done

#------------------------------------
# Run genetic correlations L and R
#------------------------------------
for p in ${phensLR}
do
if [ ! -f ${p}.adjTBV_LR.log ]
then
p1=${p}_left
p2=${p}_right
p1_file=${p1}.adjTBV.sumstats.gz
p2_file=${p2}.adjTBV.sumstats.gz
if [ -f ${p1_file} ] && [ -f ${p2_file} ]
then
echo 'Run LDSC to calculate genetic correlation between left and right for ' ${p}
ldsc.py \
--rg ${p1_file},${p2_file} \
--ref-ld-chr ${ldscores_dir}/eur_w_ld_chr/ \
--w-ld-chr ${ldscores_dir}/eur_w_ld_chr/ \
--out ${p}.adjTBV_LR
fi
fi
done
#------------------------------------

#------------------------------------
# GenCor with GWAS summary statistics
#------------------------------------ 
gwas_dir=/data/workspaces/lag/workspaces/lg-multilateral/working/Data/GWAS_sumstats/
# define phenotypes
#p2=Sniekers2017
p3=Ripke2014_SCZ
#p4=Okbay2016_EA
#p5=Wray2018_MDD
p6=Grove_ASD
p7=Demontis_ADHDeur
p8=Savage2018_IQ
#p9=2017_OCD
#p10=Sklar2012_BIP
#p10=Lambert2013_AD # skip because too few SNPs went into analysis
p11=Lee2018_EA3
#p12=Lee2018_CP
# and corresponding sumstats files
## created by:P:\workspaces\lg-multilateral\working\Analysis\amaia\GWAS_sumstats\GWAS_sumstats_ldsc.sh
#p2_file=${gwas_dir}/ldsc/${p2}.sumstats.gz
p3_file=${gwas_dir}/ldsc/${p3}.sumstats.gz
#p4_file=${gwas_dir}/ldsc/${p4}.sumstats.gz
#p5_file=${gwas_dir}/ldsc/${p5}.sumstats.gz
p6_file=${gwas_dir}/ldsc/${p6}.sumstats.gz
p7_file=${gwas_dir}/ldsc/${p7}.sumstats.gz
p8_file=${gwas_dir}/ldsc/${p8}.sumstats.gz
#p9_file=${gwas_dir}/ldsc/${p9}.sumstats.gz
#p10_file=${gwas_dir}/ldsc/${p10}.sumstats.gz
p11_file=${gwas_dir}/ldsc/${p11}.sumstats.gz
#p12_file=${gwas_dir}/ldsc/${p12}.sumstats.gz
#--------

#--------
cd ${working_dir}

for p1 in ${phenos}
 do
 p1_file=${working_dir}${p1}.adjTBV.sumstats.gz
 if [ ! -f ${p1}.adjTBV_GWASsumstats.log ]
 then
 #  --rg ${p1_file},${p2_file},${p3_file},${p4_file},${p5_file},${p6_file},${p7_file},${p8_file},${p9_file},${p10_file},${p11_file},${p12_file}
  echo Run using univariate model
  ldsc.py \
    --rg ${p1_file},${p3_file},${p6_file},${p7_file},${p8_file},${p11_file} \
   --ref-ld-chr ${ldscores_dir}/eur_w_ld_chr/ \
   --w-ld-chr ${ldscores_dir}/eur_w_ld_chr/ \
   --out ${p1}.adjTBV_GWASsumstats
  fi
# if [ ! -f ./partitioned_h2/${p1}_GWASsumstats_baselineLD.log ]
# then
# echo Run using baselineLD model from Gazal et al.
# ldsc.py \
# --rg ${p1_file},${p2_file},${p3_file},${p4_file},${p5_file},${p6_file},${p7_file} \
# --ref-ld-chr ${ldscores_dir}/Phase3/baselineLD_v1.1/baselineLD. \
# --w-ld-chr ${ldscores_dir}/Phase3/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
# --overlap-annot \
# --frqfile-chr ${ldscores_dir}/Phase3/1000G_Phase3_frq/1000G.EUR.QC. \
# --out ./partitioned_h2/${p1}_GWASsumstats_baselineLD
# fi
 done

#------------------------------------
# create summary tables
#------------------------------------
## Genetic correlations
# get line containing header of table, plus next line after match (A1)
grep -A11 p2 *log | grep -v Analysis | grep -e p2 -e sumstats | grep -v -e 'log-$' -e GWAS -e ALSPAC > summary_p1_p2_gencor_ldsc.table
grep -A11 p2 *GWAS*log | grep -v Analysis | grep -e p2 -e sumstats | grep -v -e 'log-$' -e 's$' -e 'p1  $' > summary_p1_p2_GWASgencor_ldsc.table

## h2 estimates
# columns to get:
##Total Observed scale h2: 0.0103 (0.0027)
##Lambda GC: 1.0466
##Mean Chi^2: 1.0469
##Intercept: 1.0093 (0.006)
##Ratio: 0.1984 (0.1285)
h2_file=$(ls -1 *h2.log )
paste <(echo 'File') <(echo 'h2 (se)') <(echo 'Lambda GC') <(echo 'Mean Chi^2') <(echo 'Intercept (se)') <(echo 'Ratio (se)')> summary_h2_ldsc.table
paste <(ls -1 ${h2_file}) \
      <(grep "h2:" ${h2_file} | awk '{print $5,$6}' ) \
      <(grep "Lambda GC:" ${h2_file} | awk '{print $3}') \
      <(grep "Mean Chi^2" ${h2_file} | awk '{print $3}') \
      <(grep "Intercept" ${h2_file} | awk '{print $2,$3}') \
      <(grep "Ratio" ${h2_file} | awk '{print $2,$3}') >> summary_h2_ldsc.table


#------------------------------------
# GenCor with ALSPAC measures
#------------------------------------ 
alspac_dir=/data/workspaces/lag/workspaces/lg-multilateral/working/Data/ALSPAC/ldsc/
# define phenotypes
p2=f7ws076
p3=f8ws110
p4=f8sl040
p5=f8sl105
p6=fg5830
# and corresponding sumstats files
p2_file=${alspac_dir}${p2}.sumstats.gz
p3_file=${alspac_dir}${p3}.sumstats.gz
p4_file=${alspac_dir}${p4}.sumstats.gz
p5_file=${alspac_dir}${p5}.sumstats.gz
p6_file=${alspac_dir}${p6}.sumstats.gz

#--------
# re-define phens of interest to run assoc with sumstats from publicly available phenos
#phens2=$(echo $phens | sed 's/ /\n/g' | grep -e $p -e $p2)
phens2=$(echo $h2phenosLR $h2phenos)
#--------
cd ${working_dir}

for p1 in ${phenos}
 do
 p1_file=${working_dir}${p1}.adjTBV.sumstats.gz
# if [ ! -f ${p1}.adjTBV_ALSPACsumstats.log ]
# then
  echo Run using univariate model
  ldsc.py \
   --rg ${p1_file},${p2_file},${p3_file},${p4_file},${p5_file},${p6_file} \
   --ref-ld-chr ${ldscores_dir}/eur_w_ld_chr/ \
   --w-ld-chr ${ldscores_dir}/eur_w_ld_chr/ \
   --out ${p1}.adjTBV_ALSPACsumstats_int0 --intercept-gencov 0,0,0,0,0,0
#  fi
 done
grep -A5 p2 *ALSPAC*log | grep -v Analysis | grep -e p2 -e sumstats | grep -v  -e 'log-$' -e ':p1' > summary_p1_p2_ALSPACgencor_ldsc.table



#------------------------------------
# clean and organize all files
#------------------------------------
#mkdir -p  ${working_dir}/logs/
#mv  ${working_dir}*log  ${working_dir}/logs/
