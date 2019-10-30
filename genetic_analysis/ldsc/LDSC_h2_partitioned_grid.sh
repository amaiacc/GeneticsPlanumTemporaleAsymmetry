#!/bin/sh
#$ -N ldsc_partitionedh2
#$ -cwd
#$ -q single.q
#$ -S /bin/bash
#$ -M amaia.carrioncastillo@mpi.nl
#$ -m beas
#----------------------------------------------------------------------
ldsc_dir=/usr/local/apps/ldsc/
#--------------------------
p=$1
working_dir=$2
ldscores_dir=$3
cts_name=$4
#--------------------------

echo 'Estimate partitioned heritability for' ${p}
#---------------
if [ ! -f ${working_dir}${p}_h2_baseline.log ]
   then
#    echo '-----------> Baseline model (Finucane et al. 2015)'
#    ${ldsc_dir}ldsc.py \
    --h2 ${working_dir}${p}.sumstats.gz \
    --ref-ld-chr ${ldscores_dir}/Phase3/baseline_v1.1/baseline. \
    --w-ld-chr ${ldscores_dir}/Phase3/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
    --overlap-annot \
    --frqfile-chr ${ldscores_dir}/Phase3/1000G_Phase3_frq/1000G.EUR.QC. \
    --print-coefficients \
    --out ${working_dir}${p}_h2_baseline
fi
#---------------
if [ ! -f ${working_dir}${p}_h2_baselineLD.log ]
then
#    echo '-----------> Baseline LD model (Gazal et al. 2017)'
#    ${ldsc_dir}ldsc.py \
    --h2 ${working_dir}${p}.sumstats.gz \
    --ref-ld-chr ${ldscores_dir}/Phase3/baselineLD_v2.2/baselineLD. \
    --w-ld-chr ${ldscores_dir}/Phase3/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
    --overlap-annot \
    --frqfile-chr ${ldscores_dir}/Phase3/1000G_Phase3_frq/1000G.EUR.QC. \
    --print-coefficients \
    --out ${working_dir}${p}_h2_baselineLD 
fi
#---------------

if [ ! -f ${working_dir}${p}_h2_${cts_name}.log ]
then
cd ${ldscores_dir}/Phase3/
    echo '-----------> Cell-type analyses (Finucane et al. 2018)'
    echo ${cts_name}
    echo '-----------'
    ${ldsc_dir}ldsc.py \
    --h2-cts ${working_dir}${p}.sumstats.gz \
    --ref-ld-chr ${ldscores_dir}/Phase3/baseline_v1.1/baseline. \
    --w-ld-chr ${ldscores_dir}/Phase3/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
    --ref-ld-chr-cts ${ldscores_dir}/Phase3/${cts_name}.ldcts \
    --out ${working_dir}${p}_h2_${cts_name} \
    --print-coefficients 
fi

#---------------