# locuszoom 
PATH=/data/workspaces/lag/shared_spaces/Resource_DB/locuszoom/bin/:${PATH}

subset_name=imagingT1_N18057
# 
ld_dir=/data/workspaces/lag/workspaces/lg-ukbiobank/working_data/amaia/genetic_data/PT/LD/
assoc_dir=/data/workspaces/lag/workspaces/lg-ukbiobank/working_data/amaia/genetic_data/release_v3/imp/bgenie/output/PT/${subset_name}/clean/
list=${assoc_dir}/clumped/all.lead.snps
working_dir=/data/workspaces/lag/workspaces/lg-ukbiobank/working_data/amaia/genetic_data/PT/${subset_name}/locuszoom/


mkdir -p ${working_dir}
cd ${working_dir}

# create input files for locuszoom: for AI, L and R
## association results in metal format, i.e. MarkerName P-value
#zless ${assoc_dir}AI_VOLUME_Planum_Temporale_CHRall_1e-07HWEp_0.7INFO_0.001MAF.txt.gz
## MarkerName=rsid
## P-value=P
#--------------------

dist=200

while read line
 do 
 snp=$(echo ${line} | awk '{print $1}')
 chr=$(echo ${line} | awk '{print $2}' | awk -F":" '{print $1}' )
 pos=$(echo ${line} | awk '{print $2}' | awk -F":" '{print $2}' )
 if [ "$snp" != "rsid" ]
 then
 echo '------------------'  
 echo 'Format data for locuszoom (metal format): '${snp}', chr'${chr}
 zgrep -e '^'${chr}' ' -e chr ${assoc_dir}AI_VOLUME_Planum_Temporale_CHRall_1e-07HWEp_0.7INFO_0.001MAF_short.txt.gz | awk -v OFS='\t'  '{print $2,$10}' > AI_VOLUME_Planum_Temporale_CHRall_1e-07HWEp_0.7INFO_0.001MAF_chr${chr}.input4locuszoom
 zgrep -e '^'${chr}' ' -e chr ${assoc_dir}AI_VOLUME_Planum_Temporale_males_CHRall_1e-07HWEp_0.7INFO_0.001MAF_short.txt.gz | awk -v OFS='\t'  '{print $2,$10}' > AI_VOLUME_Planum_Temporale_females_CHRall_1e-07HWEp_0.7INFO_0.001MAF_chr${chr}.input4locuszoom
 zgrep -e '^'${chr}' ' -e chr ${assoc_dir}Volume_of_grey_matter_in_Planum_Temporale_left_CHRall_1e-07HWEp_0.7INFO_0.001MAF_short.txt.gz | awk -v OFS='\t'  '{print $2,$10}' > Volume_of_grey_matter_in_Planum_Temporale_left_CHRall_1e-07HWEp_0.7INFO_0.001MAF_chr${chr}.input4locuszoom
 zgrep -e '^'${chr}' ' -e chr ${assoc_dir}Volume_of_grey_matter_in_Planum_Temporale_right_CHRall_1e-07HWEp_0.7INFO_0.001MAF_short.txt.gz | awk -v OFS='\t'  '{print $2,$10}' > Volume_of_grey_matter_in_Planum_Temporale_right_CHRall_1e-07HWEp_0.7INFO_0.001MAF_chr${chr}.input4locuszoom
 echo '------------------'
 # first plot, conditional on existing chr specific db
 if [ ! -d AI_VOLUME_Planum_Temporale_${dist}kb_noRecomb_${snp} ]
 then
  if [ -f ld_cache_chr${chr}.db ]
  then
  locuszoom --metal AI_VOLUME_Planum_Temporale_CHRall_1e-07HWEp_0.7INFO_0.001MAF_chr${chr}.input4locuszoom \
   --markercol rsid --pvalcol P --refsnp ${snp} --flank ${dist}kb \
   --pop EUR --build hg19 --source 1000G_March2012 \
   --prefix AI_VOLUME_Planum_Temporale_${dist}kb_noRecomb --no-date showRecomb=F \
   --cache ld_cache_chr${chr}.db
  else 
  locuszoom --metal AI_VOLUME_Planum_Temporale_CHRall_1e-07HWEp_0.7INFO_0.001MAF_chr${chr}.input4locuszoom \
   --markercol rsid --pvalcol P --refsnp ${snp} --flank ${dist}kb \
   --pop EUR --build hg19 --source 1000G_March2012 \
   --prefix AI_VOLUME_Planum_Temporale_${dist}kb_noRecomb --no-date showRecomb=F
  mv ld_cache.db ld_cache_chr${chr}.db
  fi
 fi
 # the rest using the chr specific db
 for p in AI_VOLUME_Planum_Temporale_males AI_VOLUME_Planum_Temporale_females Volume_of_grey_matter_in_Planum_Temporale_left Volume_of_grey_matter_in_Planum_Temporale_right
 do
  if [ ! -d ${p}_${dist}kb_noRecomb_${snp} ]
  then
  locuszoom --metal ${p}_CHRall_1e-07HWEp_0.7INFO_0.001MAF_chr${chr}.input4locuszoom \
   --markercol rsid --pvalcol P --refsnp ${snp} --flank ${dist}kb \
   --pop EUR --build hg19 --source 1000G_March2012 \
   --prefix ${p}_${dist}kb_noRecomb   --no-date showRecomb=F
  fi
 done
 fi
 done < ${list}
