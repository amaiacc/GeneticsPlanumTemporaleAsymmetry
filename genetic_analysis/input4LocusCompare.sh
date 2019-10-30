# create input file for locus compare:
# http://locuscompare.com/
## format:
# rsid	pval

working_dir=/data/workspaces/lag/workspaces/lg-ukbiobank/working_data/amaia/genetic_data/PT/imagingT1_N18057/
cd ${working_dir}
mkdir -p locuscompare
cd locuscompare

phenotype=AI_VOLUME_Planum_Temporale


# chr 10
chr=10
snp=rs41298373

# chr 2
#chr=2
#snp=rs7420166

file=${working_dir}/locuszoom/${phenotype}_CHRall_1e-07HWEp_0.7INFO_0.001MAF_chr${chr}.input4locuszoom

# get line containing lead SNP
n=$(grep ${snp} ${file} -nr | awk -F":" '{print $1}')
X=$((${n}-5000))
Y=$((${n}+5000))
cat <(echo 'rsid pval') <(sed 's/ /\t/g' ${file} | head -n "$Y" | tail -n +"$X") > ${phenotype}_CHRall_1e-07HWEp_0.7INFO_0.001MAF_chr${chr}_${snp}.input4locuscompare