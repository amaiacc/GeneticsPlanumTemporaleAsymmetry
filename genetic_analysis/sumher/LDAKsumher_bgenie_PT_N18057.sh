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
working_dir=/data/workspaces/lag/workspaces/lg-ukbiobank/working_data/amaia/genetic_data/PT/${subset_name}/sumher/
gwas_dir=/data/workspaces/lag/workspaces/lg-multilateral/working/Data/GWAS_sumstats/sumher/
#--------
pheno_file=/data/workspaces/lag/workspaces/lg-ukbiobank/working_data/amaia/genetic_data/release_v3/imp/bgenie/${pheno_root}_${type}_sample_${subset_name}_${region}_phenos4BGENIE.table
mkdir -p ${working_dir}
cd ${working_dir}
#--------
head -n 1 ${pheno_file} | sed 's/ /\n/g' > ${working_dir}pheno_header.table
# define phenotypes:
grep ${region} pheno_header.table
phenos=$(grep ${region} pheno_header.table | grep -v ukb | sed 's/residuals_//g')
#---------------------------------------------------------------------- 

#----------------------------------------------------------------------
# Preliminaries: 
#----------------------------------------------------------------------
## Reference panel: done, as in LDAK_sumher_IDPs_bgenie.sh

# create input file and run h2s
for p in ${phenos}
 do
 #p_file=${assoc_dir}${p}_CHRall_0.7INFO_0.001MAF.txt.gz
 p_file=${assoc_dir}${p}_CHRall_1e-07HWEp_0.7INFO_0.001MAF.txt.gz
 # define N, should get it from pheno file...
 n=$(grep -nr $p$ pheno_header.table | awk -F':' '{print $1}')
 N1=$(awk -v c=${n} '{print $c}' ${pheno_file} | grep -v -e '\-999' -e ${p} | wc -l) # need to get the right number!
 if [ -f ${p_file} ] && [ ! -f ${p}.txt ] 
 then
 #gunzip -c $p_file | head
 # filer on INFO --> # info=INFO.UKB
 # direction=beta
 # stat=T
 ## check cols
 gunzip -c $p_file | head | awk '{print $36,$1,$2,$5,$4,$15,$17,$39}'
 # filter and format file
 gunzip -c $p_file | awk '$36>0.95 {print}' | awk -v N=$N1 '(NR>1){snp=$2;a1=$5;a2=$4;dir=$15;stat=($17)^2;p=$39;n=N}(NR==1)\
 {print "Predictor A1 A2 Direction Stat P n"}(NR>1 && (a1=="A"||a1=="C"||a1=="G"||a1=="T") \
 && (a2=="A"||a2=="C"||a2=="G"||a2=="T")){print snp, a1, a2, dir, stat, p, n}' > ${p}.txt
  # check if there are duplicate predictor names
   awk < ${p}.txt '{print $1}' | sort | uniq -d | head
   # remove duplicates (except one), using magic command
   mv ${p}.txt ${p}2.txt
   awk '!seen[$1]++' ${p}2.txt > ${p}.txt
   # To create a list of predictors for which we have summary statistics, use
   awk < ${p}.txt '(NR>1){print $1}' > ${p}.predictors
   # To reduce these to predictors with non-ambiguous alleles, use
   awk <  ${p}.txt '( ($2=="A"&&$3=="C") || ($2=="A"&&$3=="G") || ($2=="C"&&$3=="A") || ($2=="C"&&$3=="T") || ($2=="G"&&$3=="A") || ($2=="G"&&$3=="T") || ($2=="T"&&$3=="C") || ($2=="T"&&$3=="G") ){print $1}' > ${p}.nonamb
   # To identify which (reference panel) predictors are in the MHC, use
   awk < ${ldak_dir}ref_1kg/ref.bim '($1==6 && $4>25000000 && $4<34000000){print $2}' > mhc.snps
   # The variance explained by a predictor is Stat/(Stat+n), so we can create a list of large-effect predictors using; here $5 is stat=T, $7 is N
   awk '(NR>1){print $5}' ${p}.txt | sort -n | tail # check which are the largest T values
   awk < ${p}.txt '(NR>1 && $5>$7/99){print $1}' > ${p}.big
   # There are no large-effect predictors ... but were there some, we could then identify which (reference panel) predictors are in LD with these using the command
   ##ldak.out --remove-tags ${p} --bfile ref --top-preds ${p}.big --window-cm 1 --min-cor .1 ##${p}.out 
   cat mhc.snps > ${p}.exclude # 
  fi
 done

# remove back-up files
ls *2.txt
rm *2.txt

# check whether the exclude files are very different
wc *.exclude
comm -23 *left.exclude *right.exclude | wc # suppress common -> nothing left
comm -23 AI*.exclude *right.exclude | wc

# since they are the same, we can use just use one of these to create the weights and tagging
cat *.exclude | sort | uniq > h2phenos_PT.exclude
wc h2phenos_PT.exclude

# Check nonamb and exclude files across PT phenotypes
wc *.nonamb
awk '(NR==FNR){arr[$1];next}($1 in arr)' AI_VOLUME_Planum_Temporale.nonamb Volume_of_grey_matter_in_Planum_Temporale_right.nonamb > tmp.nonamb
awk '(NR==FNR){arr[$1];next}($1 in arr)' tmp.nonamb Volume_of_grey_matter_in_Planum_Temporale_left.nonamb > h2phenos_PT.nonamb
rm tmp.nonamb
wc *.nonamb


#----------------------------------------------------
# Heritability model, get weightings
#----------------------------------------------------
# The basic way to compute weightings is using, for LDAK model
ldak --cut-weights sumsect --bfile ${ldak_dir}ref_1kg/ref --extract h2phenos_PT.nonamb --exclude h2phenos_PT.exclude
ldak --calc-weights-all sumsect --bfile ${ldak_dir}ref_1kg/ref --extract h2phenos_PT.nonamb --exclude h2phenos_PT.exclude
#----------------------------------------------------
# Tagging File
#----------------------------------------------------
# calculate a tagging file under different models
#---------------
## LDAK Model
#---------------
#ldak --calc-tagging sumldak --bfile ${ldak_dir}ref_1kg/ref --weights sumsect/weights.short  --window-cm 1 --extract ${p}.predictors --exclude ${p}.exclude --power -0.25
for j in {1..22}; do
 if [ ! -f sumldak${j}.tagging ]
 then
 ldak5.linux --calc-tagging sumldak${j} --bfile ${ldak_dir}ref_1kg/ref --weights sumsect/weights.short --power -0.25 --extract ${p}.predictors --exclude ${p}.exclude --window-cm 1 --chr ${j}
 fi
done
rm list.txt; for j in {1..22}; do echo "sumldak$j.tagging" >> list.txt; done
ldak5.linux --join-tagging sumldak --taglist list.txt
rm list.txt

#---------------
## GCTA/LDSC Model
# To instead create a tagging file assuming the LDSC Model, replace --weights sumsect/weights.short --power -0.25 with --ignore-weights YES --power -1.
#---------------
for j in {1..22}; do
 if [ ! -f sumgcta${j}.tagging ]
 then
  ldak5.linux --calc-tagging sumgcta${j} --bfile ${ldak_dir}ref_1kg/ref --ignore-weights YES --power -1 --extract ${p}.predictors --exclude ${p}.exclude --window-cm 1 --chr ${j}
 fi
done
rm list.txt; for j in {1..22}; do echo "sumgcta$j.tagging" >> list.txt; done
ldak5.linux --join-tagging sumgcta --taglist list.txt
rm list.txt

#---------------
# clean intermediate files
for j in {1..9}; do rm sumgcta${j}*; done
for j in {1..9}; do rm sumldak${j}*; done
#----------------------------------------------------

#----------------------------------------------------
# Run heritabilitites
#----------------------------------------------------
for p in ${phenos} #${h2phenosLR} #${h2phenos}
 do
 echo ${p}
 mkdir -p ${p}
 #if [ ! -f ${p}.ldak.log ]
  #then
  echo 'Estimate SNP heritability using LDAK and GCTA models'
   # Once you have a Tagging File, you can estimate SNP heritabililty:
   ldak5.linux --sum-hers ${p}.ldak --tagfile sumldak.tagging --summary ${p}.txt --check-sums NO > ${p}.ldak.log
   ldak5.linux --sum-hers ${p}.gcta --tagfile sumgcta.tagging --summary ${p}.txt --check-sums NO > ${p}.gcta.log
  #fi
done


#----------------------------------------------------
# extract h2 estimates from hers files
paste <(echo 'File') <(echo 'h2') <(echo 'h2_SD') <(echo 'Size') <(echo 'Mega_Intensity') <(echo 'Int_SD')> summary_h2_estimates_sumher.table
#Component Heritability Her_SD Size Mega_Intensity Int_SD
paste <(ls -1 *hers) <(grep Her_ALL *hers | awk -F ":Her_ALL " '{print $2}' ) | sed 's/ /\t/g' >> summary_h2_estimates_sumher.table
#----------------------------------------------------------------------

