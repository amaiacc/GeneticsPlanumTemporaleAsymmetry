#!/bin/bash

PATH=/usr/local/apps/ldsc/:${PATH}
resource_dir=/data/workspaces/lag/shared_spaces/Resource_DB/
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
covs_name=noBioCovs_noAssessmentC
region=Planum_Temporale
type=Volume
#--------
assoc_dir=/data/workspaces/lag/workspaces/lg-ukbiobank/working_data/amaia/genetic_data/release_v3/imp/bgenie/output/PT/${subset_name}/clean/
working_dir=/data/workspaces/lag/workspaces/lg-ukbiobank/working_data/amaia/genetic_data/PT/${subset_name}/ldsc/
UKB_phenos_dir=/data/workspaces/lag/workspaces/lg-ukbiobank/working_data/amaia/pheno_files/genetic_v2/${pheno_root}/summary_phenotypes/
working_dir2=/data/clusterfs/lag/users/amacar/ukb/
#--------
# file for genetic correlations
gwas_dir=/data/workspaces/lag/workspaces/lg-multilateral/working/Data/GWAS_sumstats/ldsc/
# define phenotypes
#p2=Sniekers2017
p3=Ripke2014_SCZ
p4=Okbay2016_EA
#p5=Wray2018_MDD
p6=Grove_ASD
p7=Demontis_ADHDeur
p8=Savage2018_IQ
#p9=2017_OCD
#p10=Sklar2012_BIP
#p10=Lambert2013_AD # skip because too few SNPs went into analysis
p11=Lee2018_EA3
#p12=Lee2018_CP
#  and corresponding sumstats files
## created by:P:\workspaces\lg-multilateral\working\Analysis\amaia\GWAS_sumstats\GWAS_sumstats_ldsc.sh
#p2_file=${gwas_dir}/${p2}.sumstats.gz
p3_file=${gwas_dir}/${p3}.sumstats.gz
#p4_file=${gwas_dir}/${p4}.sumstats.gz
#p5_file=${gwas_dir}/${p5}.sumstats.gz
p6_file=${gwas_dir}/${p6}.sumstats.gz
p7_file=${gwas_dir}/${p7}.sumstats.gz
p8_file=${gwas_dir}/${p8}.sumstats.gz
#p9_file=${gwas_dir}/${p9}.sumstats.gz
#p10_file=${gwas_dir}/${p10}.sumstats.gz
p11_file=${gwas_dir}/${p11}.sumstats.gz
#p12_file=${gwas_dir}/${p12}.sumstats.gz
#--------

mkdir -p ${working_dir}
cd ${working_dir}

for sex in males females
 do
 type=Volume_${sex}
 pheno_file=$(ls ${UKB_phenos_dir}/${pheno_root}_imaging_${covs_name}_${region}_${sex}_Phenotypes_residuals_wHeader.table)
 head -n 1 ${pheno_file} | sed 's/\t/\n/g'> ${key}_${sex}_Phenotypes_residuals.header
 # define phenotypes:
 grep ${region} pheno_header.table
 phenos=$(grep ${region} pheno_header.table | grep -v ukb | sed 's/residuals_//g')
 phensLR=$(echo $phenos | sed "s/ /\\n/g" | grep -v ukb | grep -v AI | sed 's/_left//g' | sed 's/_right//g' |  sort | uniq)

 # BGENIE: 'the regression model we code the first and second alleles as 0 and 1 respectively, so the beta coefficient refers to the effect of having an extra copy of the second allele.'
 # a_0 The character code for the first allele (a string).
 # a_1 The character code for the second allele (a string). --> effect allele
 ## https://jmarchini.org/bgenie-usage/

 for p in ${phenos}
  do
  echo ${p}
  ps=${p}_${sex}
  if [ ! -f ${ps}.sumstats.gz ]
  then
   # define N, should get it from pheno file...
   n=$(grep -nr $p$ pheno_header.table | awk -F':' '{print $1}')
   N1=$(awk -v c=${n} '{print $c}' ${pheno_file} | grep -v -e '\-999' -e ${p} | wc -l) # need to get the right number!
   #p_file=${assoc_dir}${p}_CHRall_1e-07HWEp_0.7INFO_0.001MAF_short.txt.gz
   p_file=${assoc_dir}${ps}_CHRall_1e-07HWEp_0.7INFO_0.001MAF.txt.gz
   if [ -f ${p_file} ]
    then
    echo 'Run munge_sumstats to create input files for LDSC - ' ${ps}
    munge_sumstats.py \
     --sumstats ${p_file} \
     --N ${N1} \
     --p P \
     --snp rsid \
     --info INFO.UKB \
     --out ${ps} \
     --a1 a_1 --a2 a_0 \
     --merge-alleles ${ldscores_dir}/w_hm3.snplist
    #---------------
     echo 'Estimate heritability for' ${ps}
     ldsc.py \
     --h2 ${ps}.sumstats.gz \
     --ref-ld-chr ${ldscores_dir}/eur_w_ld_chr/ \
     --w-ld-chr ${ldscores_dir}/eur_w_ld_chr/ \
     --out ${ps}_h2
    #---------------
    fi
   fi
  done

 #------------------------------------
 # Run genetic correlations L and R
 #------------------------------------
 for p in ${phensLR}
 do
 ps=${p}_${sex}
 if [ ! -f ${ps}_LR.log ]
  then
  p1=${p}_left_${sex}
  p2=${p}_right_${sex}
  p1_file=${p1}.sumstats.gz
  p2_file=${p2}.sumstats.gz
  if [ -f ${p1_file} ] && [ -f ${p2_file} ]
   then
    echo 'Run LDSC to calculate genetic correlation between left and right for ' ${p}
    ldsc.py \
    --rg ${p1_file},${p2_file} \
    --ref-ld-chr ${ldscores_dir}/eur_w_ld_chr/ \
    --w-ld-chr ${ldscores_dir}/eur_w_ld_chr/ \
    --out ${ps}_LR
   fi;fi
  done
 #------------------------------------

 #------------------------------------
 # GenCor with GWAS summary statistics
 #------------------------------------ 

 for p1 in ${phenos}
  do
  p1_file=${working_dir}${p1}_${sex}.sumstats.gz
  if [ ! -f ${p1}_${sex}_GWASsumstats.log ]
  then
  echo Run using univariate model
  ldsc.py \
   --rg ${p1_file},${p3_file},${p6_file},${p7_file},${p8_file},${p11_file} \
   --ref-ld-chr ${ldscores_dir}/eur_w_ld_chr/ \
   --w-ld-chr ${ldscores_dir}/eur_w_ld_chr/ \
   --out ${p1}_${sex}_GWASsumstats
  fi
 done
 
 # sex specific sumstats for ADHD and EA (Okbay2016)
 p4s_file=${gwas_dir}/${p4}_${sex}.sumstats.gz
 p7s_file=${gwas_dir}/${p7}_${sex}.sumstats.gz
 for p1 in ${phenos}
  do
  p1_file=${working_dir}${p1}_${sex}.sumstats.gz
  if [ ! -f ${p1}_${sex}_2_GWASsumstats.log ]
  then
  echo Run using univariate model
  ldsc.py \
   --rg ${p1_file},${p4s_file},${p7s_file} \
   --ref-ld-chr ${ldscores_dir}/eur_w_ld_chr/ \
   --w-ld-chr ${ldscores_dir}/eur_w_ld_chr/ \
   --out ${p1}_${sex}_2_GWASsumstats
  fi
 done
 
 
done


#   --rg ${p1_file},${p2_file},${p3_file},${p4_file},${p5_file},${p6_file},${p7_file},${p8_file},${p9_file},${p10_file},${p11_file},${p12_file} 

#------------------------------------
# create summary tables * hemen nago *
#------------------------------------
## Genetic correlations
# get line containing header of table, plus next x lines after match (A1)
grep -A11 p2 *log | grep -v Analysis | grep -e p2 -e sumstats | grep -v -e 'log-$' -e 'time elapsed' > summary_p1_p2_gencor_ldsc.table
grep -A11 p2 *GWAS*log | grep -v Analysis | grep -e p2 -e sumstats | grep -v -e 'log-$' -e 'time elapsed' > summary_p1_p2_GWASgencor_ldsc.table

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

