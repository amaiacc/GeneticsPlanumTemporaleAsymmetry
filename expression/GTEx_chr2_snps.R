# Check eQTL effects of lead SNPs
## data downloaded from GTEx
library(dplyr)
#----------------------------------------------------------------------
# define working_dirs
if (Sys.info()['sysname']=='Windows') {dir="P://workspaces/"} else {dir="/data/workspaces/lag/workspaces/"}
#----------------------------------------------------------------------
subset_name="imagingT1_N18057"
ld_dir=paste(dir,"lg-ukbiobank/working_data/amaia/genetic_data/PT/",subset_name,"/lead_snps/LD/",sep="")
working_dir=paste(dir,"lg-ukbiobank/working_data/amaia/genetic_data/PT/GTEx/",sep="")
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
ld<-subset(ld,R2>=0.4)
#---------------
# read GTEx output from: GTEx_eQTL_calculator_input.sh
for (f in list.files(pattern="GTEx*.*csv")){
  t<-read.csv(f)
  if (f==list.files(pattern="GTEx*.*csv")[1]) { data<-t
  } else { data<-rbind(data,t)}
  rm(t)
}
rm(f)
# remove lines that contain no info
data<-data[grep("not sufficiently expressed|N/A",data$Gene.Symbol,invert=TRUE),]
#----------------------------------------------------------------------
# combine
data_ld<-merge(data,ld,by.x="SNP",by.y="SNP_B",all.x=TRUE)
# remove duplicated rows
w<-which(duplicated(data_ld))
if (length(w)>0){
  data_ld<-data_ld[-w,]
}
rm(w)
table(data_ld$Gene.Symbol,data_ld$Tissue)
# subset correcting for number of tested genes
genes<-unique(data$Gene.Symbol) # 11
# subset only to cortex
data_c<-subset(data_ld,Tissue=="Brain - Cortex")

# wide format
data_wide<-data_c[,c("SNP","CHR_B","BP_B","R2","Gene.Symbol","Tissue","P.Value")] %>% spread(Gene.Symbol,P.Value)
# select only SNPxgene combinations where Pval for either tissue is < 0.05/length(genes)
thr=0.05/length(genes)
w<-unlist(apply(data_wide[,genes],2,function(x){
   which(x<thr)
  }))
data_wide[w,]
# sort by chr position
data_wide<-data_wide[order(data_wide$BP_B),]
write.csv(data_wide,"GTEx_BrainCortex_wide.csv",row.names = FALSE)
