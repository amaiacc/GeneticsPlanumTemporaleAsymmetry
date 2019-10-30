#!/bin/sh
#$ -N magma__PHEN_
#$ -cwd
#$ -q single.q
#$ -S /bin/bash
#$ -M amaia.carrioncastillo@mpi.nl
#$ -m beas
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# Template to run MAGMA on summary statistics
#----------------------------------------------------------------------
# parameters to define input files
## one example, hard defined
#pheno_root=ukb21288_ukb21293
#phen=AI_VOLUME_Planum_Temporale
#pheno_root=_PHENOROOT_
#phen=_PHEN_
pheno_root=${1}
phen=${2}
up_win=${3} #35
down_win=${4} #10
subset_name=${5}

# if Windows not specified, set to 0,0
if [ -z ${up_win} ]; then up_win=0; fi
if [ -z ${down_win} ]; then down_win=0; fi
#---------------------------------------------------------------------- 

#----------------------------------------------------------------------
# Define paths and directories
#----------------------------------------------------------------------
# magma
PATH=/data/workspaces/lag/shared_spaces/Resource_DB/magma_v1.06b/:${PATH}
#----------------------------------------------------------------------
resource_dir=/data/workspaces/lag/shared_spaces/Resource_DB/
analysis_dir=/data/workspaces/lag/workspaces/lg-ukbiobank/analysis/amaia/genetic_analysis/PT/
magma_dir=${resource_dir}/magma_v1.06b/
#--------
#primary_dir=/data/workspaces/lag/workspaces/lg-ukbiobank/working_data/amaia/genetic_data/release_v3/imp/bgenie/output/PT/clean/
#working_dir=/data/workspaces/lag/workspaces/lg-ukbiobank/working_data/amaia/genetic_data/PT/magma/
primary_dir=/data/workspaces/lag/workspaces/lg-ukbiobank/working_data/amaia/genetic_data/release_v3/imp/bgenie/output/PT/${subset_name}/clean/
working_dir=/data/workspaces/lag/workspaces/lg-ukbiobank/working_data/amaia/genetic_data/PT/${subset_name}/magma/
#--------
out_prefix=${phen}_CHRall_1e-07HWEp_0.7INFO_0.001MAF.txt #_${pheno_root}
input_file=$(ls ${primary_dir}/${out_prefix}*.gz)
# get N from phenotype file
type=Volume
region=Planum_Temporale
pheno_file=$(ls /data/workspaces/lag/workspaces/lg-ukbiobank/working_data/amaia/genetic_data/release_v3/imp/bgenie/${pheno_root}_*_sample_${subset_name}_${region}_phenos4BGENIE.table)
head -n 1 ${pheno_file} | sed 's/ /\n/g' > ${working_dir}pheno_header.table
n=$(grep -nr $phen$ ${working_dir}pheno_header.table | awk -F':' '{print $1}')
N=$(awk -v c=${n} '{print $c}' ${pheno_file} | grep -v -e '\-999' -e ${phen} | wc -l) # need to get the right number!
rm ${working_dir}pheno_header.table
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# create dir
mkdir -p ${working_dir}
cd ${working_dir}
# create directory for this magma run and move there
mkdir -p ${working_dir}/${phen}
cd ${working_dir}/${phen}

#----------------------------------------------------------------------
# Preliminaries: create correct SNP_loc and input files
#----------------------------------------------------------------------
# unzip file and extract just the first three cols
## SNP CHR BP P
if [ ! -f ${out_prefix}.SNP.loc ]; then
 zcat $input_file | head -n 1 | sed 's/ /\n/g' > header.input.txt
 snp=$(grep '^rsid$' -nr header.input.txt  | head -n 1| awk -F':' '{print $1}')
 chr=$(grep '^chr$' -nr header.input.txt | head -n 1 | awk -F':' '{print $1}')
 bp=$(grep '^pos$' -nr header.input.txt | head -n 1 | awk -F':' '{print $1}')
 p=$(grep '^P$' -nr header.input.txt | head -n 1 | awk -F':' '{print $1}')
# define output file name
 zcat ${input_file} | awk -v snp=${snp} -v chr=${chr} -v bp=${bp} -v p=${p} '{print $snp,$chr,$bp,$p}' > ${out_prefix}.SNP.loc
fi

snploc=${out_prefix}.SNP.loc
synonyms=${magma_dir}dbsnp/dbsnp147.synonyms
geneloc=${magma_dir}/NCBI/NCBI37.3.gene.loc


# get subset of gene sets, only GO, KEGG, REACTOME

#----------------------------------------------------------------------
# Annotate: get definition of genes (i.e. SNPs within genes)
#----------------------------------------------------------------------
# default, without window; to include window around genes:-annotate window=5,1.5’
#up_win=35
#down_win=10
magma --annotate window=${up_win},${down_win} --snp-loc ${snploc} --gene-loc ${geneloc} --out ${up_win}_${down_win}kb


#----------------------------------------------------------------------
# Gene analysis on SNP p-value data
#----------------------------------------------------------------------
# the genotype data specified by --bfile is used to specify the reference data used to estimate LD between SNPs
magma --bfile ${magma_dir}/g1000/g1000_eur synonyms=${synonyms} \
 --gene-annot ${up_win}_${down_win}kb.genes.annot \
 --pval ${snploc} use=rsid,P N=${N} \
 --gene-model snp-wise=mean \
 --out ${up_win}_${down_win}kb_gene

#----------------------------------------------------------------------
# Gene-set analysis: GO sets from MSIGDB
# use latest version (6.1), and v5.2 for reference with FUMA
#----------------------------------------------------------------------
for ver in 6.1 # 5.2 
 do 
 for c in 5 # 2
 do
  setfile=$(echo ${magma_dir}/msigdb/msigdb_v${ver}_GMTs/c${c}.all.v${ver}.entrez.gmt)
  magma --gene-results ${up_win}_${down_win}kb_gene.genes.raw --set-annot ${setfile} --out sets_${up_win}_${down_win}kb_c${c}_v${ver}
 done
 done
 




