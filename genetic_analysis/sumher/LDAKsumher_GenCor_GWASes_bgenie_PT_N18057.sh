#!/bin/bash
PATH=/data/workspaces/lag/shared_spaces/Resource_DB/ldak5.linux_/:${PATH}
#cd /data/workspaces/lag/shared_spaces/Resource_DB/ldak5.linux_/
#ln -s ldak5.linux ldak
resource_dir=/data/workspaces/lag/shared_spaces/Resource_DB/
ldscores_dir=${resource_dir}/LDscores/
ldak_dir=${resource_dir}/ldak5/
#--------
# define variables, again, just for the easy of copy-pasting...
pheno_root=ukb25465_ukb25468
subset_name=imagingT1_N18057
region=Planum_Temporale
type=Volume
#--------
assoc_dir=/data/workspaces/lag/workspaces/lg-ukbiobank/working_data/amaia/genetic_data/release_v3/imp/bgenie/output/PT/${subset_name}/clean/
working_dir=/data/workspaces/lag/workspaces/lg-ukbiobank/working_data/amaia/genetic_data/PT/${subset_name}/sumher/GenCor_GWASes/
gwas_dir=/data/workspaces/lag/workspaces/lg-multilateral/working/Data/GWAS_sumstats/sumher/
#--------
pheno_file=/data/workspaces/lag/workspaces/lg-ukbiobank/working_data/amaia/genetic_data/release_v3/imp/bgenie/${pheno_root}_${type}_sample_${subset_name}_${region}_phenos4BGENIE.table
mkdir -p ${working_dir}
cd ${working_dir}
#---------------------------------------------------------------------- 
# Define phenotypes
#--------------------------------------
# PT phenos
#--------------------------------------
head -n 1 ${pheno_file} | sed 's/ /\n/g' > ${working_dir}pheno_header.table
# define phenotypes:
grep ${region} pheno_header.table
phenos=$(grep ${region} pheno_header.table | grep -v ukb | sed 's/residuals_//g')
#--------------------------------------
#--------------------------------------
# donwnloaded and organized the summary statistics from publicly available GWASes
#--------------------------------------
# files are formatted for LDAK/sumher (see LDAK_sumher_IDPs.sh)
gwas_dir=/data/workspaces/lag/workspaces/lg-multilateral/working/Data/GWAS_sumstats/sumher/
#------------------------------------
# GenCor with GWAS summary statistics
#------------------------------------ 
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
#---------------------------------------------------------------------- 

#----------------------------------------------------------------------
## Reference panel: done, as in LDAK_sumher_IDPs_bgenie.sh
# Preliminaries: formatting input file per phenotype, already done in: LDAKsumher_bgenie_PT_N18057.sh
mkdir -p ${working_dir}
cd ${working_dir}
# combine exclusion and nonamb markers
wc ../h2phenos*
#

awk '(NR==FNR){arr[$1];next}($1 in arr)' ../h2phenos_PT.nonamb ${gwas_dir}/all.nonamb > h2phenos_PT_gwas.nonamb
echo '' >> h2phenos_PT_gwas.nonamb # add newline at the end of file, otherwise get error, last snps cannot be correctly read
cat ../h2phenos_PT.exclude ${gwas_dir}/all.exclude > h2phenos_PT_gwas.exclude
wc h2phenos_PT_gwas*
# 140482   140482  1494654 h2phenos_PT_gwas.exclude
# 3847337  3847336 40237245 h2phenos_PT_gwas.nonamb
#--------

#----------------------------------------------------
# Get weightings
#----------------------------------------------------
# The basic way to compute weightings is using, for LDAK model
echo 'Calculate weightings'
for j in {1..22}; do
 ldak --cut-weights corsect${j} --bfile ${ldak_dir}ref_1kg/ref --extract h2phenos_PT_gwas.nonamb --exclude h2phenos_PT_gwas.exclude --chr $j
 ldak --calc-weights-all corsect${j} --bfile ${ldak_dir}ref_1kg/ref --extract h2phenos_PT_gwas.nonamb --exclude h2phenos_PT_gwas.exclude --chr $j
done
mkdir corsect
cat corsect*/weights.short > corsect/weights.short
mv corsect* ./corsect/
#----------------------------------------------------
# Tagging File
#----------------------------------------------------
# calculate a tagging file under different models
# (speed up by dividing by chromosome)
#---------------
## LDAK Model
#---------------
for j in {1..22}; do
 if [ ! -f corldak${j}.tagging ]
 then
 ldak5.linux --calc-tagging corldak${j} --bfile ${ldak_dir}ref_1kg/ref --weights corsect/weights.short --power -0.25 --extract h2phenos_PT_gwas.nonamb --exclude h2phenos_PT_gwas.exclude --window-cm 1 --chr ${j}
 fi
done
rm list.txt; for j in {1..22}; do echo "corldak$j.tagging" >> list.txt; done
ldak5.linux --join-tagging corldak --taglist list.txt
rm list.txt

#---------------
## GCTA/LDSC Model
# To instead create a tagging file assuming the LDSC Model, replace --weights sumsect/weights.short --power -0.25 with --ignore-weights YES --power -1.
#---------------
for j in {1..22}; do
 if [ ! -f corgcta${j}.tagging ]
 then
  ldak5.linux --calc-tagging corgcta${j} --bfile ${ldak_dir}ref_1kg/ref --ignore-weights YES --power -1 --extract h2phenos_PT_gwas.nonamb --exclude h2phenos_PT_gwas.exclude --window-cm 1 --chr ${j}
 fi
done
rm list.txt; for j in {1..22}; do echo "corgcta$j.tagging" >> list.txt; done
ldak5.linux --join-tagging corgcta --taglist list.txt
rm list.txt

#----------------------------------------------------------------------
# Run genetic correlations
#----------------------------------------------------------------------
#ldak_tagging_dir=${working_dir}

mkdir -p ${working_dir}/
cd ${working_dir}

# we regress the two sets of summary statistics onto the tagging file (there will likely be a few inconsistent alleles, so add --check-sums NO to ignore).
## p1 - each of the target phenotypes, AIs
## p2 - GWAS sums
for p in ${phenos}
 do
 echo ${p}
 mkdir -p ${p}
 p_file=${working_dir}/../${p}.txt
 for p_gwas in ${p3} ${p6} ${p7} ${p8} ${p11} 
 do
 p_gwas_file=${gwas_dir}${p_gwas}.txt
  if [ -f ${p_gwas_file} ]
  then
   echo '-----------------------------------------'
   echo Run cor for ${p} and ${p_gwas}
   if [ ! -f ${p}/${p}_${p_gwas}_corldak.cors ]
   then
   echo Using LDAK model
   ldak --sum-cors ${p}/${p}_${p_gwas}_corldak --tagfile corldak.tagging --summary ${p_file} --summary2 ${p_gwas_file} --genomic-control YES --check-sums NO > ${p}/${p}_${p_gwas}_corldak.log
   ldak --sum-cors ${p}/${p}_${p_gwas}_corldak_noGC --tagfile corldak.tagging --summary ${p_file} --summary2 ${p_gwas_file} --check-sums NO > ${p}/${p}_${p_gwas}_corldak_noGC.log
   fi
   if [ ! -f ${p}/${p}_${p_gwas}_corgcta.cors ]
   then
   echo Using GCTA model
   ldak --sum-cors ${p}/${p}_${p_gwas}_corgcta --tagfile corgcta.tagging --summary ${p_file} --summary2 ${p_gwas_file} --genomic-control YES --check-sums NO > ${p}/${p}_${p_gwas}_corgcta.log
   ldak --sum-cors ${p}/${p}_${p_gwas}_corgcta_noGC --tagfile corgcta.tagging --summary ${p_file} --summary2 ${p_gwas_file} --check-sums NO > ${p}/${p}_${p_gwas}_corgcta_noGC.log
   fi
  fi
  done
done

#----------------------------------------------------
# summary file for genetic correlations
#----------------------------------------------------
# extract h2 estimates from cors files
paste <(echo 'File') <(echo 'rg') <(echo 'rg_SD')> summary_rg_estimates_sumher.table
#Component Heritability Her_SD Size Mega_Intensity Int_SD
paste <(ls -1 *Planum*/*cors) <(grep Cor_ALL *Planum*/*cors | awk -F ":Cor_ALL " '{print $2}' ) | sed 's/ /\t/g' >> summary_rg_estimates_sumher.table
#----------------------------------------------------------------------
