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
#----------------------------------------------------------------------
mkdir -p ${working_dir}/output ${analysis_dir}

#----------------------------------------------------------------------
# define variables
pheno_root=ukb25465_ukb25468
subset_name=imagingT1_N18057
region=Planum_Temporale
type=Volume
#----------------------------------------------------------------------

UKB=/data/workspaces/lag/workspaces/lg-ukbiobank/primary_data/genetic_data/
UKB_cal=/data/workspaces/lag/workspaces/lg-ukbiobank/working_data/amaia/genetic_data/release_v2/cal/QC/
UKB_phenos_dir=/data/workspaces/lag/workspaces/lg-ukbiobank/working_data/amaia/genetic_data/release_v3/imp/bgenie/
UKB_imp=/data/workspaces/lag/workspaces/lg-ukbiobank/working_data/amaia/genetic_data/release_v3/imp/subset_${subset_name}/

# template for association, males and females

## chr and pheno file as parameters + region as well!
echo "#----------------------------------------------------------------------
#!/bin/sh
#$ -N BGENIE_assoc_UKBv3_chr_per_sex
#$ -cwd
#$ -q single.q
#$ -S /bin/bash
#$ -M amaia.carrioncastillo@mpi.nl
#$ -m beas
#----------------------------------------------------------------------
working_dir2=/data/clusterfs/lag/users/amacar/ukb/
bgenie_dir=\${working_dir2}/bin/bgenie/v1.2/
UKB_imp=\${working_dir2}/input/imp/

PATH=\${working_dir2}/bin/bgenie/v1.2/:\$PATH
#----------------------------------------------------------------------
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
sex=\$6
snp=\$7

if [ \$i == "X" ];  then   sample=sampleX;   else   sample=sample;   fi

# if snp specified, include is as an argument...
echo \${snp}
if [ ! \${snp}='' ]; then rsid=\$(echo '--rsid ' \${snp}); else echo No SNP defined, will run whole chromosome \${chr}; fi

#rsid=\$(echo '--rsid ' \${snp})
echo This is the rsid parameter for bgenie: 
echo \${rsid}

# define pheno file
pheno_file=\$(ls \${working_dir2}/pheno/\${pheno_root}*\${pheno_type}*\${sample}*\${subset_name}*\${region}*_\${sex}*table)
head -n1 \$pheno_file | sed 's/ /\n/g' | grep res > \${working_dir2}/pheno/\${pheno_root}_\${subset_name}_\${pheno_type}_\${region}_\${sex}.cov

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
  --pheno \${pheno_file} \${rsid} \
  --include_pheno \${working_dir2}/pheno/\${pheno_root}_\${subset_name}_\${pheno_type}_\${region}_\${sex}.cov \
  --out \${working_dir2}/output/\${pheno_root}_\${subset_name}_bgenie_\${pheno_type}_\${region}_\${sex}_chr\${i}.out > \${working_dir2}/logs/\${pheno_root}_\${subset_name}_bgenie_\${pheno_type}_\${region}_\${sex}_chr\${i}\${snp}.log

 else
 echo 'Output file already exists, skip.'
 fi

# timestamp
date +\"%D %R\" " > ${analysis_dir}bgenie_pheno_type_region_perSex_SNP_single.sh


# specify which phenotypes should be included in include_pheno,test just one chr, one phenotype
cp ${analysis_dir}bgenie_pheno_type_region_perSex_SNP_single.sh ${working_dir2}/sge_jobs/


# after association_bgenie_impUKBv3_imaging_PT.sh, check how the top hits perform in the sex-specific analyses
# 1- create residuals (wihtout sex as covariate): ukb21288_ukb21293_LM_residuals_imagingSubset_perSex.R
for sex in males females
 do
 # 2- make input file for bgenie (samples in the same order as sample file)
 # create a pheno file that matches sample file, for an autosome and chr22
 # these are residualized phenotypes, so no need to specify covariates when running BGENIE
 Rscript ${analysis_dir}/association_phenotypes4BGENIE_PT_perSex.R ${working_dir2}/input/imp/ukb_imp_chr10_v3_${subset_name}.sample "A" ${pheno_root} ${subset_name} ${sex}
 Rscript ${analysis_dir}/association_phenotypes4BGENIE_PT_perSex.R ${working_dir2}/input/imp/ukb_imp_chrX_v3_${subset_name}.sample "X" ${pheno_root} ${subset_name} ${sex}
done
# copy phenotype data to clusterfs
ls -lh ${UKB_phenos_dir}/*${subset_name}*${region}*males*phenos4BGENIE.table
cp ${UKB_phenos_dir}/*${subset_name}*${region}*males*phenos4BGENIE.table ${working_dir2}/pheno/



#----------------------------------------------------------------------
# BGENIE
# https://jmarchini.org/bgenie-usage/
#----------------------------------------------------------------------

# copy all necessary data to to clusfterfs
# bgenie
if [ ! -d ${working_dir2}/bin/bgenie/v1.2/ ]
then
 mkdir -p ${working_dir2}/bin/bgenie/v1.2/
 cp ${resource_dir}/bgenie/v1.2 ${working_dir2}/bin/bgenie/v1.2/
fi
## imputed data
for i in {1..22} X
 do
 if [ ! -f ${working_dir2}/input/imp/ukb_imp_chr${i}_v3_${subset_name}.bgen ]
 then
 # ls ${UKB_imp}*chr${i}_*${subset_name}.bgen -lh
  cp ${UKB_imp}*chr${i}_*${subset_name}.bgen ${working_dir2}/input/imp/
 fi
 done

cp ${UKB_imp}*chr1_*.sample ${working_dir2}/input/imp/
cp ${UKB_imp}*chrX_*.sample ${working_dir2}/input/imp/
#----------------------------------------------------------------------

# 3- Run BGENIE
cd ${working_dir2}/sge_jobs/
# define variables, again, just for the easy of copy-pasting...
pheno_root=ukb25465_ukb25468
subset_name=imagingT1_N18057
region=Planum_Temporale
type=Volume
#----------------------------------------------------------------------
# run all chrs, per sex
#----------------------------------------------------------------------
cd ${working_dir2}/sge_jobs/
for chr in {1..22} X
do
 for sex in males females
  do
  if [ ! -f ${working_dir2}/output/${pheno_root}_${subset_name}_bgenie_${type}_${region}_${sex}_chr${chr}.out.gz ]
  then
   echo Submit: ${pheno_root} ${subset_name} ${type} ${chr} ${region} ${sex}
   qsub ${analysis_dir}bgenie_pheno_type_region_perSex_SNP_single.sh ${pheno_root} ${subset_name} ${type} ${chr} ${region} ${sex}
  fi
 done
done
# single SNP, to check
chr=10
snp=rs41298373
for sex in males females
 do
 #qsub ${analysis_dir}bgenie_pheno_type_region_perSex_SNP_single.sh ${pheno_root} ${subset_name} ${type} ${chr} ${region} ${sex} ${snp}
done
#----------------------------------------------------------------------
# clean intermediate files from clusterfs
mv ${working_dir2}/output/*${subset_name}*chr*.out.gz ${working_dir}/output/PT/
mv ${working_dir2}/logs/*${subset_name}*bgenie*chr*.log ${working_dir}/logs/
# clean also sge job logs
mv ${working_dir2}/sge_jobs/BGENIE_* ${analysis_dir}/sge_jobs/bgenie/
mv ${working_dir2}/sge_jobs/bgenie/BGENIE_* ${analysis_dir}/sge_jobs/bgenie/
mv ${working_dir2}/sge_jobs/QCtool_UKBv3.* ${analysis_dir}/sge_jobs/qctool/
#----------------------------------------------------------------------
# run bgenie_output_filter.sh to combine with QCed SNPs (snps2keep, from: qctool_subset_imaging_N18057_imp.sh)
## bash ${analysis_dir}/bgenie_output_filter.sh


# next, run: BGENIE_output_PT_perSex_N18057.R to compare results from males and females
cd ${working_dir}/output/PT/${subset_name}/clean
for f in $(ls *males*INFO*MAF*.txt)
 do
 if [ ! -f ${f}.gz ]; then  gzip ${f}; fi
 done

# still too big for FUMA (max 600Mb), so remove unnecessary cols
zless AI_VOLUME_Planum_Temporale_males_CHRall_1e-07HWEp_0.7INFO_0.001MAF.txt.gz 
zcat AI_VOLUME_Planum_Temporale_males_CHRall_1e-07HWEp_0.7INFO_0.001MAF.txt.gz | awk '{print $1,$2,$3,$4,$5,$6,$15,$16,$17,$23,$20,$22}' | head
zcat AI_VOLUME_Planum_Temporale_males_CHRall_1e-07HWEp_0.7INFO_0.001MAF.txt.gz | awk '{print $1,$2,$3,$4,$5,$6,$15,$16,$17,$23,$20,$22}' | gzip > AI_VOLUME_Planum_Temporale_males_CHRall_1e-07HWEp_0.7INFO_0.001MAF_short.txt.gz
zcat AI_VOLUME_Planum_Temporale_females_CHRall_1e-07HWEp_0.7INFO_0.001MAF.txt.gz | awk '{print $1,$2,$3,$4,$5,$6,$15,$16,$17,$23,$20,$22}' | gzip > AI_VOLUME_Planum_Temporale_females_CHRall_1e-07HWEp_0.7INFO_0.001MAF_short.txt.gz
