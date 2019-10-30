PATH=/data/workspaces/lag/shared_spaces/Resource_DB/inrich.v.1.1/:${PATH}

# load modules
module load plink/1.9b6 
#
resource_dir=/data/workspaces/lag/shared_spaces/Resource_DB/
inrich_dir=/data/workspaces/lag/shared_spaces/Resource_DB/inrich.v.1.1/
KGPref=${resource_dir}/1KG_phase1_all
ld_ref_eur=${KGPref}/1kg_phase1_eur
#
primary_dir=/data/workspaces/lag/workspaces/lg-ukbiobank/working_data/amaia/genetic_data/release_v3/imp/bgenie/output/PT/${subset_name}/clean/
magma_dir=/data/workspaces/lag/workspaces/lg-ukbiobank/working_data/amaia/genetic_data/PT/${subset_name}/magma/
working_dir=/data/workspaces/lag/workspaces/lg-ukbiobank/working_data/amaia/genetic_data/PT/${subset_name}/inrich/
mkdir -p ${working_dir}
#--------
cd ${working_dir}
subset_name=imagingT1_N18057
phen=AI_VOLUME_Planum_Temporale
out_prefix=${phen}_CHRall_1e-07HWEp_0.7INFO_0.001MAF.txt #_${pheno_root}
magma_file=$(ls ${magma_dir}/${phen}/${out_prefix}*.loc) # use input for magma, simplified: rsid chr pos P
# create new variable as chr:pos, and save file with: rsid, SNP, P
awk '(NR==1){print "rsid SNP P"}(NR>1){SNP=$2":"$3; print $1,SNP,$4}' ${magma_file} > ${phen}.assoc

#-------------------------------------------------------
# Create input files for INRICH
#-------------------------------------------------------
#--------
# Reference SNP File (-m): List of reference SNPs examined in the association study
#Col 1: Chromosome
#Col 2: SNP Base Pair Position
awk '{ print $2 }' ${phen}.assoc | sed 's/:/ /g' > tmp
grep -v SNP tmp > ${phen}.snp.map
rm tmp
#--------
# Associated Interval File (-a): List of LD-independent associated interval regions for enrichment tests
#Col 1  Chromosome
#Col 2: Interval Start Base Pair Position
#Col 3: Interval End Base Pair Position

# run clumping in plink to define intervals, use the same parameters as Tulio did for his 2015 paper (Guadalupe et al. 2015, Cortex)
## e.g.
##Options in effect:
##	--noweb
##	--bfile BIG
##	--clump ../SHIP2_females.csv
##	--clump-p1 0.001
##	--clump-p2 0.05
##	--clump-r2 0.5
##	--clump-range empty.txt
##	--clump-range-border 20
##	--out SHIP2_females

# use 1KG eur for LD reference
plink --bfile ${ld_ref_eur} \
      --clump ${phen}.assoc \
      --clump-p1 0.001 \
      --clump-p2 0.05 \
      --clump-r2 0.5 \
      --clump-range ${inrich_dir}/ref/entrez_gene.hg19.map \
      --clump-range-border 20 \
      --out ${phen}

awk 'NR>1 {print $5 }' ${phen}.clumped.ranges | sed 's/:/ /g' | sed 's/\.\./ /g' | sed 's/chr//g' > ${phen}.p.0.001.int

#-------------------------------------------------------
# Run INRICH
#-------------------------------------------------------
inrich -g ${inrich_dir}/ref/entrez_gene.hg19.map -m ${phen}.snp.map -t ${inrich_dir}/sets/go.set -a AI_VOLUME_Planum_Temporale.p.0.001.int

# using parameters and candidate gene set from Tulio Guadalupe's runs
inrich -g ${inrich_dir}/ref/entrez_gene.hg19.map \
        -m ${phen}.snp.map  -a ${phen}.p.0.001.int \
        -t andro_estro_steroid_autis.set \
        -i 10 -j 200 -z 5 -w 100000 -d 0.1 -r 10000 -q 5000 -p 1 \
        -o ${phen}_Guadalupe2015sets

#-------------------------------------------------------


# try to replicate the results from Guadalupe et al., using the same data from Guadalupe but the same reference genes as in the UKB run
guadalupe_dir=/data/workspaces/lag/workspaces/lg-ukbiobank/analysis/amaia/genetic_analysis/PT/Guadalupe2015/enrichment/input_files/

# BIG males
inrich -g ${inrich_dir}/ref/entrez_gene.hg19.map \
        -m ${guadalupe_dir}BIGsnps.map  -a ${guadalupe_dir}BIG_males.intervals \
        -t ${guadalupe_dir}andro_estro_steroid_autis.set \
        -i 10 -j 200 -z 5 -w 100000 -d 0.1 -r 10000 -q 5000 -p 1 \
        -f 2058903252 \
        -o BIG_males_seed
