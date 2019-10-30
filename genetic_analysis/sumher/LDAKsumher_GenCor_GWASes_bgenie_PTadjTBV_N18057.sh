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
key=totalBV
type=Volume
#--------
assoc_dir=/data/workspaces/lag/workspaces/lg-ukbiobank/working_data/amaia/genetic_data/release_v3/imp/bgenie/output/PT/${subset_name}/clean/
working_dir=/data/workspaces/lag/workspaces/lg-ukbiobank/working_data/amaia/genetic_data/PT/${subset_name}/sumher/GenCor_GWASes/
gwas_dir=/data/workspaces/lag/workspaces/lg-multilateral/working/Data/GWAS_sumstats/sumher/
#--------
pheno_file=/data/workspaces/lag/workspaces/lg-ukbiobank/working_data/amaia/genetic_data/release_v3/imp/bgenie/${pheno_root}_${type}_sample_${subset_name}_${region}_${key}_phenos4BGENIE.table
mkdir -p ${working_dir}
cd ${working_dir}
#--------
head -n 1 ${pheno_file} | sed 's/ /\n/g' > ${working_dir}pheno_header.table
# define phenotypes:
grep ${region} pheno_header.table
phenos=$(grep ${region} pheno_header.table | grep -v -e ukb -e AIfa -e AIabs | sed 's/residuals_//g' | sort | uniq )
#---------------------------------------------------------------------- 

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
# Run genetic correlations
#----------------------------------------------------------------------
ldak_tagging_dir=${working_dir}

mkdir -p ${working_dir}/
cd ${working_dir}

# p.adjTBV as phenotype

# we regress the two sets of summary statistics onto the tagging file (there will likely be a few inconsistent alleles, so add --check-sums NO to ignore).
## p1 - each of the target phenotypes, AIs
## p2 - GWAS sums
for p in ${phenos}
 do
 echo ${p}
 p=${p}.adjTBV
 mkdir -p ${p}
 p_file=${working_dir}/../${p}.txt
 for p_gwas in ${p3} ${p6} ${p7} ${p8} ${p11} # ${p2} ${p3} ${p4} ${p5} ${p6} ${p7} ${p8} ${p9} ${p10} ${p11} ${p12}
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
