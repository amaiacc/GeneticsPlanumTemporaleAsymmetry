# Check eQTL effects of lead SNPs
## data downloaded from GTEx
library(dplyr)
#----------------------------------------------------------------------
# define working_dirs
if (Sys.info()['sysname']=='Windows') {dir="P://workspaces/"} else {dir="/data/workspaces/lag/workspaces/"}
#----------------------------------------------------------------------
subset_name="imagingT1_N18057"
ld_dir=paste(dir,"lg-ukbiobank/working_data/amaia/genetic_data/PT/",subset_name,"/lead_snps/LD/",sep="")
working_dir=paste(dir,"lg-ukbiobank/working_data/amaia/genetic_data/PT/xQTL/queries/",sep="")
setwd(working_dir)
options(stringsAsFactors = FALSE)
#----------------------------------------------------------------------
# read ld files
ld_files<-paste(ld_dir,list.files(ld_dir,pattern="r2.*_rs*.*.ld"),sep="")
for (f in ld_files){
  t<-read.table(f,header=T,strip.white = T)
  if (f==ld_files[1]){ld<-t} else {ld<-rbind(ld,t)}
  rm(t)
}
rm(f,ld_files)
w<-which(duplicated(ld))
if (length(w)>0){
  ld<-ld[-w,]
}
rm(w)
ld<-subset(ld,R2>=0.6)

#---------------
# read GTEx output from: GTEx_eQTL_calculator_input.sh
for (f in list.files(pattern="eSTLs_query.csv")){
  t<-read.csv(f)
  if (f==list.files(pattern="eSTLs_query.csv")[1]) { data<-t
  } else { data<-rbind(data,t)}
  rm(t)
}
rm(f)
#----------------------------------------------------------------------
# combine
data_ld<-merge(data,ld,by.x=c("SNPid","SNPchromosome","SNPpos"),by.y=c("SNP_B","CHR_B","BP_B"))
table(data_ld$Gene.Symbol,data_ld$Tissue)
# subset correcting for number of tested genes
genes<-unique(data$featureName) # 11
#----------------------------------------------------------------------
# wide format
data_wide<-data_ld[,c("SNPid","SNPchromosome","SNPpos","R2","featureName","pValue")] %>% spread(featureName,pValue)
# select only SNPxgene combinations where Pval for either tissue is < 0.05/length(genes)
thr=0.05/11
w<-unlist(apply(data_wide[,genes],2,function(x){
   which(x<thr)
  }))
data_wide[w,]
# sort by chr position
data_wide<-data_wide[order(data_wide$SNPpos),]
write.csv(data_wide,"xQTL_DLPC_eQTLs_wide.csv",row.names = FALSE)
