# subset data from OBrien et al to all the genes/transcripts/SNPs of interest
#----------------------------------------------------------------------
subset_name=imagingT1_N18057
#
primary_dir=/data/workspaces/lag/shared_spaces/Resource_DB/OBrien2018_supplementary/
working_dir=/data/workspaces/lag/workspaces/lg-ukbiobank/working_data/amaia/genetic_data/PT/OBrien2018/
ld_dir=/data/workspaces/lag/workspaces/lg-ukbiobank/working_data/amaia/genetic_data/PT/${subset_name}/lead_snps/LD/
fuma_dir=/data/workspaces/lag/workspaces/lg-ukbiobank/working_data/amaia/genetic_data/PT/${subset_name}/FUMA/FUMA_job36766/
mkdir -p ${working_dir}
cd ${working_dir}
#----------------------------------------------------------------------
# get header from all files to query:
## gene level
zcat ${primary_dir}expression_gene.bed.gz | head -n 1 > ${working_dir}/expression_genes_interest.bed
zcat ${primary_dir}all_eqtls_gene.txt.gz | head -n 1 > ${working_dir}/all_eqtls_gene_interest.txt
zcat ${primary_dir}top_eqtls_gene.txt.gz | head -n 1 > ${working_dir}/top_eqtls_gene_interest.txt
## transcript level
zcat ${primary_dir}expression_transcript.bed.gz | head -n 1 > ${working_dir}/expression_transcript_interest.bed
zcat ${primary_dir}all_eqtls_transcript.txt.gz | head -n 1 > ${working_dir}/all_eqtls_transcript_interest.txt
zcat ${primary_dir}top_eqtls_transcript.txt.gz | head -n 1 > ${working_dir}/top_eqtls_transcript_interest.txt
## snp level
touch ${working_dir}/snp_positions_interest.bed
zcat ${primary_dir}all_eqtls_gene.txt.gz | head -n 1 > ${working_dir}/all_eqtls_snp_interest.txt
#zcat ${primary_dir}snp_positions.bed.gz | head -n 1 > ${working_dir}/snp_positions_interest.bed
#----------------------------------------------------------------------


#----------------------------------------------------------------------
## get Gencode IDs
# ITIH5: ENSG00000123243
# BOK: ENSG00000176720
# BOK-AS1: ENSG00000234235
# DTYMK: ENSG00000168393
# ING5: ENSG00000168395
transcripts='ENSG00000123243 ENSG00000176720 ENSG00000168393 ENSG00000168395 ENSG00000234235'
transcripts=$(awk -F"," '{print $6}' ${fuma_dir}/blood_FDR_eqtl.csv | sort | grep -v gene | uniq | sed 's/"//g')

#----------------------------------------------------------------------

for t in ENSG00000123243 ${transcripts}
 do
 ## gene level
 zgrep ${t} ${primary_dir}expression_gene.bed.gz >> ${working_dir}/expression_genes_interest.bed
 zgrep ${t} ${primary_dir}all_eqtls_gene.txt.gz >> ${working_dir}/all_eqtls_gene_interest.txt
 zgrep ${t} ${primary_dir}top_eqtls_gene.txt.gz >> ${working_dir}/top_eqtls_gene_interest.txt
 ## transcript level
 zgrep ${t} ${primary_dir}expression_transcript.bed.gz >> ${working_dir}/expression_transcript_interest.bed
 zgrep ${t} ${primary_dir}all_eqtls_transcript.txt.gz >> ${working_dir}/all_eqtls_transcript_interest.txt
 zgrep ${t} ${primary_dir}top_eqtls_transcript.txt.gz >> ${working_dir}/top_eqtls_transcript_interest.txt
 done

#----------------------------------------------------------------------
## snp level
# get snps: r2>0.6 with lead snps, as defined by: LD_around_leadSNPs_withinUKB.sh
ld_files=$(ls ${ld_dir}/*0.6r2*.ld)
cat ${ld_files} | awk '{print $6}' | sort | uniq | grep -v SNP_B > leadSNPs_0.6r2.list
 zcat ${primary_dir}all_eqtls_gene.txt.gz | head -n 1 > ${working_dir}/all_eqtls_snp_interest.txt
while read s
 do
  echo ${s}
  #zgrep ${s} ${primary_dir}snp_positions.bed.gz >> ${working_dir}/snp_positions_interest.bed
  zgrep ${s} ${primary_dir}all_eqtls_gene.txt.gz >> ${working_dir}/all_eqtls_snp_interest.txt
 done < leadSNPs_0.6r2.list
 #----------------------------------------------------------------------