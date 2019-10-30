library(tidyr)
library(ggplot2)
options(stringsAsFactors = FALSE)
#----------------------------------------------------------------------
# define working_dirs
if (Sys.info()['sysname']=='Windows') {dir="P://workspaces/"} else {dir="/data/workspaces/lag/workspaces/"}
#----------------------------------------------------------------------
subset_name="imagingT1_N18057"
primary_dir=paste(dir,"../shared_spaces/Resource_DB/OBrien2018_supplementary/",sep="")
ld_dir=paste(dir,"lg-ukbiobank/working_data/amaia/genetic_data/PT/",subset_name,"/lead_snps/LD/",sep="")
fuma_dir=paste(dir,"lg-ukbiobank/working_data/amaia/genetic_data/PT/",subset_name,"/FUMA/FUMA_job36828_1KGref/",sep="")
working_dir=paste(dir,"lg-ukbiobank/working_data/amaia/genetic_data/PT/OBrien2018/",sep="")
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
# gene symbol gencode match


#----------------------------------------------------------------------
# read info
## sample
sample_info<-read.table(paste(primary_dir,"covariates.txt",sep=""),header=TRUE)
snp_info<-read.table(paste(primary_dir,"snp_positions.bed.gz",sep=""),header=TRUE)
#----------------------------------------------------------------------
## genes
gene_names<- read.csv(paste(fuma_dir,"blood_FDR_eqtl_symbol_gencode.csv",sep=""),header=TRUE)
eqtls_genes<- read.table(paste(working_dir,"all_eqtls_gene_interest.txt",sep=""),header=TRUE,comment.char = "")
colnames(eqtls_genes)<-gsub("X\\.","",colnames(eqtls_genes))
expr_genes<-read.table(paste(working_dir,"expression_genes_interest.bed",sep=""),header=TRUE,comment.char = "")
colnames(expr_genes)<-gsub("X\\.","",colnames(expr_genes))
#----------------------------------------------------------------------
## snps
eqtls_snps<-read.table(paste(working_dir,"all_eqtls_snp_interest.txt",sep=""),header=TRUE,comment.char = "")
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# check eQTL effects: snp x gene
# subset only to genes of interst
w<-which(eqtls_snps$gene_id %in% unique(eqtls_genes$gene_id))
eqtls_snps<-eqtls_snps[w,]
subset(eqtls_snps,pval_nominal<0.1)
merge()
data_ld<-merge(eqtls_snps,ld,by.x=c("variant_id"),by.y=c("SNP_B"))
data_ld<-merge(data_ld,gene_names,by.x="gene_id",by.y="gene")
# convert to wide
data_ld_wide<-data_ld[,c("variant_id","CHR_B","BP_B","R2","symbol","pval_nominal")] %>% spread(symbol,pval_nominal)
write.csv(data_ld_wide,"OBrien2018_fetalBrain_eQTLs_wide.csv",row.names = FALSE)
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# transform
# from wide to long format: 
expr_genes_long <- gather(expr_genes,key="Sample",value="norm_count", X12545:X13144, factor_key=TRUE)
expr_genes_long$Sample<-gsub("^X","",expr_genes_long$Sample)
table(expr_genes_long$Sample)
# combine with sample info:
expr_genes_long<-merge(sample_info,expr_genes_long,by="Sample",all.y=TRUE)
# add gene names
expr_genes_long<-merge(gene_names,expr_genes_long,by="ID",all=TRUE)
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# plot expression of genes of interest across time:
# data-points per timepoint/gene

gene_time<-ggplot(data=expr_genes_long,aes(x=PCW,y=norm_count)) + geom_point() + 
  facet_grid(Sex~gene) + theme_bw() +  ggtitle("Expression in fetal brain\nO'Brien et al. (2018)")

# get highest and lowest expression values and timepoints per gene:
summary<-do.call("rbind",lapply(unique(expr_genes_long$gene),function(g) {
  tmp<-subset(expr_genes_long,gene==g)
  min_v<-min(tmp$norm_count)
  min_t<-tmp$PCW[which(tmp$norm_count==min_v)]
  max_v<-max(tmp$norm_count)
  max_t<-tmp$PCW[which(tmp$norm_count==max_v)]
  d<-as.data.frame(rbind(cbind("min",min_v,min_t),cbind("max",max_v,max_t)))
  colnames(d)<-c("val","norm_count","PCW")
  d$gene<-g
  return(d) 
  }   ))

summary<-as.data.frame(summary)
