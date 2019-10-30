# Lookup GW significant SNPs in across GWASes

subset_name=imagingT1_N18057
#----------------------------------------------------------------------
# Define paths and directories
#----------------------------------------------------------------------
working_dir=/data/workspaces/lag/workspaces/lg-ukbiobank/working_data/amaia/genetic_data/PT/${subset_name}/lead_snps/lookups/
ld_dir=/data/workspaces/lag/workspaces/lg-ukbiobank/working_data/amaia/genetic_data/PT/${subset_name}/lead_snps/LD/
assoc_dir=/data/workspaces/lag/workspaces/lg-ukbiobank/working_data/amaia/genetic_data/release_v3/imp/bgenie/output/PT/${subset_name}/clean/clumped/
data_dir=/data/workspaces/lag/workspaces/lg-ukbiobank/working_data/amaia/genetic_data/release_v3/imp/bgenie/output/PT/${subset_name}/clean/
fs_dir=/data/workspaces/lag/workspaces/lg-ukbiobank/working_data/amaia/genetic_data/release_v3/imp/bgenie/output/multilateral/clean/
gwas_dir=/data/workspaces/lag/workspaces/lg-multilateral/working/Data/GWAS_sumstats/downloaded_data/

#----------------------------------------------------------------------
mkdir -p ${working_dir}
mkdir -p ${working_dir}/PT_checks/
# all genome-wide significant SNPs identified by: clumping_results_impUKBv3_imaging_PT_N18057.sh
list=${assoc_dir}/all.lead.snps
#------------------------------------------------------------------
# Check within PT GWAS results
#------------------------------------------------------------------
cd ${data_dir}

while read line
 do 
 snp=$(echo ${line} | awk '{print $1}')
 chr=$(echo ${line} | awk '{print $2}' | awk -F":" '{print $1}' )
 pos=$(echo ${line} | awk '{print $2}' | awk -F":" '{print $2}' )
 if [ "$snp" != "rsid" ]
 then
  echo '------------------'  
  echo 'Look-ups for SNP: ' ${snp}
  ld_file=$(ls $ld_dir/LD_r2_*_${snp}_imp.ld)
  if [ -f ${ld_file} ]
  then
   ld_snps=$(awk '{if ( $7 >= 0.6) print $6}' $ld_file | grep -v SNP)
   grep_snps=$(echo $ld_snps | sed 's/ / -e /g')
  else
   grep_snps=${snp}
  fi
  # loop through all sumstats to get relevant lines
  for file in $(ls *CHRall_1e-07HWEp_0.7INFO_0.001MAF.txt.gz)
  do
  p=$(echo ${file} | sed 's/_CHRall_1e-07HWEp_0.7INFO_0.001MAF.txt.gz//g' )
   if [ ! -f ${working_dir}/PT_checks/lookup_${snp}_r2_0.6_${p}.txt ]
   then
    echo ${p}
    zcat ${file} | head -n 1 > ${working_dir}/PT_checks/lookup_${snp}_r2_0.6_${p}.txt
    zgrep -e ${grep_snps} ${file}  >> ${working_dir}/PT_checks/lookup_${snp}_r2_0.6_${p}.txt
   fi
  done
  echo '------------------'
 fi
 done < ${list}

 
#------------------------------------------------------------------
# Check within FS GWAS results
#------------------------------------------------------------------
cd ${fs_dir}
mkdir -p ${working_dir}/FS_checks/
while read line
 do 
 snp=$(echo ${line} | awk '{print $1}')
 chr=$(echo ${line} | awk '{print $2}' | awk -F":" '{print $1}' )
 pos=$(echo ${line} | awk '{print $2}' | awk -F":" '{print $2}' )
 if [ "$snp" != "rsid" ]
 then
  echo '------------------'  
  echo 'Look-ups for SNP: ' ${snp}
  ld_file=$(ls $ld_dir/LD_r2_*_${snp}_imp.ld)
  if [ -f ${ld_file} ]
  then
   ld_snps=$(awk '{if ( $7 >= 0.6) print $6}' $ld_file | grep -v SNP)
   grep_snps=$(echo $ld_snps | sed 's/ / -e /g')
  else
   grep_snps=${snp}
  fi
  # loop through all sumstats to get relevant lines
  for file in $(ls *CHRall_1e-07HWEp_0.7INFO_0.001MAF.txt.gz)
  do
  p=$(echo ${file} | sed 's/_CHRall_1e-07HWEp_0.7INFO_0.001MAF.txt.gz//g' )
   if [ ! -f ${working_dir}/PT_checks/lookup_${snp}_r2_0.6_${p}.txt ]
    then
    echo ${p}
    cat ${file} | head -n 1 > ${working_dir}/FS_checks/lookup_FS_${snp}_r2_0.6_${p}.txt
    grep -e ${grep_snps} ${file}  >> ${working_dir}/FS_checks/lookup_FS_${snp}_r2_0.6_${p}.txt
   fi
  done
  echo '------------------'
 fi
 done < ${list}


 

#------------------------------------------------------------------
# Check within GWAS results of other traits (public sumstats + other projects)
#------------------------------------------------------------------
mkdir -p ${working_dir}/GWAS_checks/

while read line
 do 
 snp=$(echo ${line} | awk '{print $1}')
 chr=$(echo ${line} | awk '{print $2}' | awk -F":" '{print $1}' )
 pos=$(echo ${line} | awk '{print $2}' | awk -F":" '{print $2}' )
 if [ "$snp" != "rsid" ]
 then
  echo '------------------'  
  mkdir -p ${working_dir}/GWAS_checks/${snp}/
  echo 'Look-ups for SNP: ' ${snp}
  ld_file=$(ls $ld_dir/LD_r2_*_${snp}_imp.ld)
  if [ -f ${ld_file} ]
  then
   ld_snps=$(awk '{if ( $7 >= 0.6) print $6}' $ld_file | grep -v SNP)
   grep_snps=$(echo $ld_snps | sed 's/ / -e /g')
  else
   grep_snps=${snp}
  fi
  # for each GWAS results...
  # EA3
  cd ${gwas_dir}Lee2018_EA3
  zgrep MarkerName GWAS_EA_excl23andMe.txt > ${working_dir}/GWAS_checks/${snp}/lookup_${snp}_Lee_EA3
  zgrep -e ${grep_snps} *txt >> ${working_dir}/GWAS_checks/${snp}/lookup_${snp}_Lee_EA3
  # intelligence 2018
  cd ${gwas_dir}SavageJansen_IntMeta_sumstats/sumstats
  zgrep SNP SavageJansen_2018_intelligence_metaanalysis.txt > ${working_dir}/GWAS_checks/${snp}/lookup_${snp}_Savage_IQ
  zgrep -e ${grep_snps} *txt >> ${working_dir}/GWAS_checks/${snp}/lookup_${snp}_Savage_IQ
  # ADHD
  cd ${gwas_dir}PGC/ADHD
  zgrep SNP  adhd_jul2017.gz > ${working_dir}/GWAS_checks/${snp}/lookup_${snp}_PGC_ADHD
  zgrep -e ${grep_snps} *gz >> ${working_dir}/GWAS_checks/${snp}/lookup_${snp}_PGC_ADHD
  # ASD
  cd ${gwas_dir}PGC/ASD
  zgrep SNP  iPSYCH-PGC_ASD_Nov2017.gz  > ${working_dir}/GWAS_checks/${snp}/lookup_${snp}_PGC_ASD
  zgrep -e ${grep_snps} *gz >> ${working_dir}/GWAS_checks/${snp}/lookup_${snp}_PGC_ASD
  # SCZ2
  cd ${gwas_dir}PGC/SCZ2
  zgrep snpid ckqny.scz2snpres.gz > ${working_dir}/GWAS_checks/${snp}/lookup_${snp}_PGC_SCZ
  zgrep -e ${grep_snps} *gz >> ${working_dir}/GWAS_checks/${snp}/lookup_${snp}_PGC_SCZ
  # crossdisorder
  cd ${gwas_dir}PGC/PGC_cross_full
  grep snpid pgc.cross.full.2013-03.txt > ${working_dir}/GWAS_checks/${snp}/lookup_${snp}_PGC_cross
  grep -e ${grep_snps} *txt >> ${working_dir}/GWAS_checks/${snp}/lookup_${snp}_PGC_cross
  # other GWAS results
  # ALSPAC 
  cd /data/workspaces/lag/workspaces/lg-multilateral/working/Data/ALSPAC/sumstats/
  grep SNP ALSPAC.f7ws076.r.language.txt > ${working_dir}/GWAS_checks/${snp}/lookup_${snp}_ALSPAC
  grep -e ${grep_snps} *.txt >> ${working_dir}/GWAS_checks/${snp}/lookup_${snp}_ALSPAC
  # Tulio's Planum Temporale GWAS (Guadalupe 2015)
  cd /data/workspaces/lag/workspaces/lg-multilateral/working/Data/Guadalupe_BIG_PT_GWAS/meta_analysis/
  grep MarkerName meta_PT_AI1.results > ${working_dir}/GWAS_checks/${snp}/lookup_${snp}_Guadalupe_meta_aiPT
  grep -e ${grep_snps} *results >> ${working_dir}/GWAS_checks/${snp}/lookup_${snp}_Guadalupe_aiPT
  grep -e ${grep_snps} */*.tbl >> ${working_dir}/GWAS_checks/${snp}/lookup_${snp}_Guadalupe_aiPT
  chrPos=$(zgrep $snp $data_dir/AI_VOLUME_Planum_Temporale_CHRall_1e-07HWEp_0.7INFO_0.001MAF.txt.gz | awk '{print $2,$3}' | sed 's/ /:/g')
  cd /data/workspaces/lag/workspaces/lg-multilateral/working/Data/Guadalupe_BIG_PT_GWAS/GWAS/
  grep CHR ./SHIP_T/zres_Planum_Index_females-SHIP-T-females.assoc.linear > ${working_dir}/GWAS_checks/${snp}/lookup_${snp}_Guadalupe_SHIP_aiPT
  grep ${chrPos} ./SHIP*/*linear >> ${working_dir}/GWAS_checks/${snp}/lookup_${snp}_Guadalupe_SHIP_aiPT
  grep SNP ./BIG/BIG_females_SPM8_NEW.results > ${working_dir}/GWAS_checks/${snp}/lookup_${snp}_Guadalupe_BIG_aiPT
  echo '------------------'
  
 fi 
done < ${list}



#------------------------------------------------------------------
# combine all in a somewhat coherent way
#------------------------------------------------------------------
${working_dir}/GWAS_checks/

#lookups in GWASes
cat rs*/*SCZ* | sed 's/ /\t/g' | grep -v -e males > lookups_SCZ.txt
cat rs*/*ASD* | sed 's/ /\t/g'  | grep -v -e males  > lookups_ASD.txt
cat rs*/*ADHD* | sed 's/ /\t/g' | grep -v -e males  > lookups_ADHD.txt
cat rs*/*EA3* | sed 's/ /\t/g' | grep -v -e males  > lookups_EA.txt
cat rs*/*Savage_IQ* | sed 's/ /\t/g' | grep -v -e males  > lookups_IQ.txt

# lookupus in Guadalupe et al.
cat rs*/*Guadalupe* | sed 's/ /\t/g' > lookups_Guadalupe2015.txt
grep linear lookups_Guadalupe2015.txt > lookups_Guadalupe2015_SNPs.txt

#------------------------------------------------------------------
# combine output with: lookups_leadSNPs_PT.R
#------------------------------------------------------------------

