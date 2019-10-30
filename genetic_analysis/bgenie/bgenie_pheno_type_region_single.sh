#----------------------------------------------------------------------
#!/bin/sh
#$ -N BGENIE_assoc_UKBv3__TYPE__chr_CHR_
#$ -cwd
#$ -q single.q
#$ -S /bin/bash
#$ -M amaia.carrioncastillo@mpi.nl
#$ -m beas
#----------------------------------------------------------------------
PATH=/data/clusterfs/lag/users/amacar/ukb//bin/bgenie/v1.2/:$PATH
#----------------------------------------------------------------------
working_dir2=/data/clusterfs/lag/users/amacar/ukb/
bgenie_dir=/data/clusterfs/lag/users/amacar/ukb//bin/bgenie/v1.2/
UKB_imp=${working_dir2}/input/imp/

## parameters to edit
#pheno_root=ukb25465_ukb25468 # should be defined before creating template
#subset_name=
#pheno_type=_TYPE_
#i=_CHR_

#alternatively, set parameters as command line arguments:
pheno_root=$1
subset_name=$2 # e.g.
pheno_type=$3 # e.g. Volume
i=$4 # chr
region=$5 # region specific analyses, e.g. Planum Temporale

echo ${pheno_root} ${subset_name} ${type} ${chr} ${region}


if [ $i == X ];  then   sample=sampleX;   else   sample=sample;   fi

# define pheno file
pheno_file=$(ls ${working_dir2}/pheno/${pheno_root}*${pheno_type}*${sample}_*${subset_name}*${region}*table)
head -n 10 $pheno_file 
head -n1 $pheno_file | sed 's/ /\n/g' | grep res > ${working_dir2}/pheno/${pheno_root}_${subset_name}_${pheno_type}_${region}.cov

# define bgen file
bgen_file=$(ls ${UKB_imp}ukb_imp_chr${i}_v3*${subset_name}*.bgen)
echo ${bgen_file}

#----------------------------------------------------------------------
# timestamp
date +"%D %R"

if [ ! -f ${working_dir2}/output/${pheno_root}_${subset_name}_bgenie_${pheno_type}_${region}_chr${i}.out.gz ]
 then 
 # start
 echo 'Run BGENIE in '${subset_name} ' for chr' ${i} '; phenotypes: ' ${pheno_type} '; key: ' ${region}
 bgenie --bgen ${bgen_file}   --thread 8   --pvals   --pheno ${pheno_file}   --include_pheno ${working_dir2}/pheno/${pheno_root}_${subset_name}_${pheno_type}_${region}.cov   --out ${working_dir2}/output/${pheno_root}_${subset_name}_bgenie_${pheno_type}_${region}_chr${i}.out > ${working_dir2}/logs/${pheno_root}_${subset_name}_bgenie_${pheno_type}_${region}_chr${i}.log

 else
 echo 'Output file already exists, skip.'
 fi

# timestamp
date +"%D %R" 
