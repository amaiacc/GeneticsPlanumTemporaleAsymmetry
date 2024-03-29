#!/bin/bash

#----------------------------------------------------------------------
# Define paths and directories
#----------------------------------------------------------------------
# BGENIE v1.2
PATH=/data/workspaces/lag/shared_spaces/Resource_DB/bgenie/v1.2/:$PATH
#----------------------------------------------------------------------
resource_dir=/data/workspaces/lag/shared_spaces/Resource_DB/
analysis_dir=/data/workspaces/lag/workspaces/lg-ukbiobank/analysis/amaia/genetic_analysis/release_v3/
working_dir=/data/workspaces/lag/workspaces/lg-ukbiobank/working_data/amaia/genetic_data/release_v3/imp/bgenie/
working_dir2=/data/clusterfs/lag/users/amacar/ukb/
UKB=/data/workspaces/lag/workspaces/lg-ukbiobank/primary_data/genetic_data/
UKB_cal=/data/workspaces/lag/workspaces/lg-ukbiobank/working_data/amaia/genetic_data/release_v2/cal/QC/
UKB_imp=${UKB}release_v3_March2018/data/imp/
UKB_phenos_dir=/data/workspaces/lag/workspaces/lg-ukbiobank/working_data/amaia/genetic_data/release_v3/imp/bgenie/
#----------------------------------------------------------------------
mkdir -p ${working_dir}/output ${analysis_dir}

#----------------------------------------------------------------------
# define variables
pheno_root=ukb25465_ukb25468
subset_name=imagingT1_N18057
region=Planum_Temporale
type=Volume
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# copy all necessary data to to clusfterfs
# bgenie
mkdir -p ${working_dir2}/bin/bgenie/v1.2/
if [ ! -f ${working_dir2}/bin/bgenie/v1.2/bgenie ]
 then
  cp ${resource_dir}/bgenie/v1.2 ${working_dir2}/bin/bgenie/v1.2/
 fi
# imputed data
for i in {1..22} X
 do
  if [ ! -f ${working_dir2}/input/imp/ukb_imp_chr${i}_v3_${subset_name}.bgen ]
  then
   echo Copying bgen file for chromosome ${i}
   cp ${UKB_imp}ukb_imp_chr${i}_v3_${subset_name}.bgen ${working_dir2}/input/imp/
   date +"%D %R"
   echo -----------------------------------------------
  fi
 done

cp ${UKB_imp}*chr1_*.sample ${working_dir2}/input/imp/
cp ${UKB_imp}*chrX_*.sample ${working_dir2}/input/imp/

#----------------------------------------------------------------------
# Create template job to subset bgen files for imaging dataset using QCtools v.2.0
# already done:
## ${analysis_dir}qctool_subset_samples.sh

cd ${working_dir2}/sge_jobs/
cp ${analysis_dir}qctool_subset_samples.sh ${working_dir2}/sge_jobs/
# sampleQC_v2_imaging_not_related.samples

# define file with samples to include, i.e. subset
cp ${UKB_cal}/${subset_name}_samples.list ${working_dir2}/input/
incl_file=${working_dir2}/input/${subset_name}_samples.list

#----------------------------------------------------------------------
# Run QCtool subsetting, snp-stats and sample file generation: already run association_bgenie_impUKBv3_imaging_PT.sh
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# Create phenotype file, matching sample file
#----------------------------------------------------------------------
# need to create phenotype files in which the order of the phenotypes is the same as in the sample file; 
## and all the unknowns are set to -999

# sample file, output of qctool for one chromosome
# phenotype files... generated by ukb21288_ukb21293_LM_residuals_imagingSubset.R
## ${UKB_phenos_dir}/ukb21288_ukb21293_imaging_noBioCovsVolumePhenotypes_residuals_wHeader.table

# create a pheno file that matches sample file, for an autosome and chr22
Rscript ${analysis_dir}/association_phenotypes4BGENIE_PT_all_subsets.R ${working_dir2}/input/imp/ukb_imp_chr22_v3_${subset_name}.sample "A" ${pheno_root} ${subset_name}
Rscript ${analysis_dir}/association_phenotypes4BGENIE_PT_all_subsets.R ${working_dir2}/input/imp/ukb_imp_chrX_v3_${subset_name}.sample "X" ${pheno_root} ${subset_name}

# these are residualized phenotypes, so no need to specify covariates when running BGENIE

# copy phenotype data to clusterfs
cp ${UKB_phenos_dir}/*${subset_name}*${region}*phenos4BGENIE.table ${working_dir2}/pheno/
#----------------------------------------------------------------------
# BGENIE
# https://jmarchini.org/bgenie-usage/
#----------------------------------------------------------------------
# important note, for the phenotype file: 
## The subjects must be presented in the same order that they are found in the BGEN file. 
## If a SAMPLE file has been provided to accompany the BGEN files, then this order can be found by examining the SAMPLE file.

# Create script to run association analysis using BGENIE
# within clusterfs
## chr and pheno file as parameters
##${analysis_dir}bgenie_pheno_type_single.sh

## chr and pheno file as parameters + region as well!
echo "#----------------------------------------------------------------------
#!/bin/sh
#$ -N BGENIE_assoc_UKBv3__TYPE__chr_CHR_
#$ -cwd
#$ -q single.q
#$ -S /bin/bash
#$ -M amaia.carrioncastillo@mpi.nl
#$ -m beas
#----------------------------------------------------------------------
PATH=${working_dir2}/bin/bgenie/v1.2/:\$PATH
#----------------------------------------------------------------------
working_dir2=/data/clusterfs/lag/users/amacar/ukb/
bgenie_dir=${working_dir2}/bin/bgenie/v1.2/
UKB_imp=\${working_dir2}/input/imp/

## parameters to edit
#pheno_root=${pheno_root} # should be defined before creating template
#subset_name=
#pheno_type=_TYPE_
#i=_CHR_

#alternatively, set parameters as command line arguments:
pheno_root=\$1
subset_name=\$2 # e.g.
pheno_type=\$3 # e.g. Volume
i=\$4 # chr
region=\$5 # region specific analyses, e.g. Planum Temporale

if [ $i == "X" ];  then   sample=sampleX;   else   sample=sample;   fi

# define pheno file
pheno_file=\$(ls \${working_dir2}/pheno/\${pheno_root}*\${pheno_type}*\${sample}*\${subset_name}*\${region}*table)
head -n1 \$pheno_file | sed 's/ /\n/g' | grep res > \${working_dir2}/pheno/\${pheno_root}_\${subset_name}_\${pheno_type}_\${region}.cov

# define bgen file
bgen_file=\$(ls \${UKB_imp}ukb_imp_chr\${i}_v3*\${subset_name}*.bgen)


#----------------------------------------------------------------------
# timestamp
date +\"%D %R\"

if [ ! -f \${working_dir2}/output/\${pheno_root}_\${subset_name}_bgenie_\${pheno_type}_chr\${i}.out.gz ]
 then 
 # start
 echo 'Run BGENIE in ' \${subset_name} ' for chr' \${i} '; phenotypes: ' \${pheno_type}
 bgenie --bgen \${bgen_file} \
  --thread 8 \
  --pvals \
  --pheno \${pheno_file} \
  --include_pheno \${working_dir2}/pheno/\${pheno_root}_\${subset_name}_\${pheno_type}_\${region}.cov \
  --out \${working_dir2}/output/\${pheno_root}_\${subset_name}_bgenie_\${pheno_type}_\${region}_chr\${i}.out > \${working_dir2}/logs/\${pheno_root}_\${subset_name}_bgenie_\${pheno_type}_\${region}_chr\${i}.log

 else
 echo 'Output file already exists, skip.'
 fi

# timestamp
date +\"%D %R\" " > ${analysis_dir}bgenie_pheno_type_region_single.sh



# specify which phenotypes should be included in include_pheno,test just one chr, one phenotype
cp ${analysis_dir}bgenie_pheno_type_*.sh ${working_dir2}/sge_jobs/


# Run BGENIE
cd ${working_dir2}/sge_jobs/
chr=X
qsub bgenie_pheno_type_region_single.sh ${pheno_root} ${subset_name} ${type} ${chr} ${region}

for chr in {11..22} {1..9} # {1..20} 22 X #{1..6} # #{17..19} # 8 10 11 12..14 16 9 # 21 22
 do
  echo ${pheno_root} ${type} ${region} ${chr}
  qsub bgenie_pheno_type_region_single.sh ${pheno_root} ${subset_name} ${type} ${chr} ${region}
 done

# move jobs that have finished
#i=10 # 22 and 21
mv ${working_dir2}/output/*${subset_name}*chr*.out.gz ${working_dir}/output/PT/
for i in ${1..22} X
 do
 mv ${working_dir2}/output/*${subset_name}*chr${i}*.out.gz ${working_dir}/output/PT/
 mv ${working_dir2}/logs/*${subset_name}*bgenie*chr${i}*.log ${working_dir}/logs/
done

# clean also sge job logs
mv ${working_dir2}/sge_jobs/BGENIE_assoc_UKBv3__TYPE__chr_CHR_.* -lh ${working_dir}/output/
mv ${working_dir2}/sge_jobs/QCtool_UKBv3.* -lh ${working_dir}/output/


# run bgenie_output_filter.sh to combine with QCed SNPs (snps2keep, from: qctool_subset_imaging_N18057_imp.sh)
## bash ${analysis_dir}/bgenie_output_filter.sh
# next, run: BGENIE2Manhattan_PT_N18057.R to generate Manhattan plots, and clean data for downstream analyses

cd ${working_dir}/output/PT/${subset_name}/clean
# compress files
for f in $(ls *INFO*MAF*.txt)
 do
 gzip ${f}
 done

# still too big for FUMA (max 600Mb), so remove unnecessary cols
## all cols:
# chr rsid pos a_0 a_1 af info ID2 Chr pos rsid a_0 a_1 hrc beta se t -log10p ukb9246_ukb10785_beta ukb9246_ukb10785_se ukb9246_ukb10785_t ukb9246_ukb10785-log10p ukb21288_ukb21293_beta ukb21288_ukb21293_se ukb21288_ukb21293_t ukb21288_ukb21293-log10p ukb25465_ukb25468_beta ukb25465_ukb25468_se ukb25465_ukb25468_t ukb25465_ukb25468-log10p ukb9246_ukb21293_beta ukb9246_ukb21293_se ukb9246_ukb21293_t ukb9246_ukb21293-log10p MAF.UKB INFO.UKB info.imaging.QCtool HW_exact_p_value.imaging.QCtool P P_z
## target cols:
# chr rsid pos a_0 a_1 af beta se t P INFO.UKB HW_exact_p_value.imaging.QCtool

zless AI_VOLUME_Planum_Temporale_CHRall_1e-07HWEp_0.7INFO_0.001MAF.txt.gz 
zcat AI_VOLUME_Planum_Temporale_CHRall_1e-07HWEp_0.7INFO_0.001MAF.txt.gz | awk '{print $1,$2,$3,$4,$5,$6,$15,$16,$17,$39,$36,$38}' | head
zcat AI_VOLUME_Planum_Temporale_CHRall_1e-07HWEp_0.7INFO_0.001MAF.txt.gz | awk '{print $1,$2,$3,$4,$5,$6,$15,$16,$17,$39,$36,$38}' | gzip > AI_VOLUME_Planum_Temporale_CHRall_1e-07HWEp_0.7INFO_0.001MAF_short.txt.gz 
zcat Volume_of_grey_matter_in_Planum_Temporale_left_CHRall_1e-07HWEp_0.7INFO_0.001MAF.txt.gz  | awk '{print $1,$2,$3,$4,$5,$6,$15,$16,$17,$39,$36,$38}' | gzip > Volume_of_grey_matter_in_Planum_Temporale_left_CHRall_1e-07HWEp_0.7INFO_0.001MAF_short.txt.gz 
zcat Volume_of_grey_matter_in_Planum_Temporale_right_CHRall_1e-07HWEp_0.7INFO_0.001MAF.txt.gz | awk '{print $1,$2,$3,$4,$5,$6,$15,$16,$17,$39,$36,$38}' | gzip > Volume_of_grey_matter_in_Planum_Temporale_right_CHRall_1e-07HWEp_0.7INFO_0.001MAF_short.txt.gz 

