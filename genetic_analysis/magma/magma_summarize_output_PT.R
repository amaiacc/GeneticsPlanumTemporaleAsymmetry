#----------------------------------------------------------------------
# plot and visualize results from MAGMA runs
#----------------------------------------------------------------------
if("qqman" %in% rownames(installed.packages()) == FALSE) {install.packages("qqman")}
if("gridExtra" %in% rownames(installed.packages()) == FALSE) {install.packages("gridExtra")}
library(qqman)
library(gridExtra)
library(ggplot2)
library(pander)
library("biomaRt")
ensembl <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl") #) #, )# listAttributes(ensembl)
options(stringsAsFactors = FALSE) 
#----------------------------------------------------------------------
# define working_dirs
if (Sys.info()['sysname']=='Windows') {dir="P://workspaces/"} else {dir="/data/workspaces/lag/workspaces/"}
subset_name="imagingT1_N18057"
magma_dir=paste(dir,"lg-ukbiobank/working_data/amaia/genetic_data/PT/",subset_name,"/magma/",sep="")
# working_dir=paste(magma_dir,"bgenie_info/AI_VOLUME_Planum_Temporale_ukb21288_ukb21293/",sep="")
# working_dir=paste(magma_dir,"multi_snp_wise_AI_VOLUME_Planum_Temporale_CHRall_1e-07HWEp_0.7INFO_0.001MAF.txt.gz/",sep="")
working_dir=paste(magma_dir,"AI_VOLUME_Planum_Temporale/",sep="")
setwd(working_dir)
## read definitions of gene-sets used by magma
misg<-"../shared_spaces/Resource_DB//magma_v1.06b//msigdb/msigdb_v6.1_GMTs/"
c5<-read.table(paste(dir,misg,"c5.all.v6.1.entrez.gmt",sep=""),header=FALSE,sep=" ")
# c2<-read.table(paste(dir,misg,"c2.all.v6.1.entrez.gmt",sep=""),header=FALSE,sep=" ")

#----------------------------------------------------------------------
# Read data for PT
#----------------------------------------------------------------------
# s<-read.table(paste(p,"/sets_c5_v6.1.sets.out",sep=""),header=TRUE)
# g<-read.table(paste(p,"/magma.genes.out",sep=""),header=TRUE)
p="AI_VOLUME_Planum_Temporale"
gene_def="35_10kb" #"35_10kb"
pt_gene<-read.table(paste(gene_def,"_gene.genes.out",sep=""),header=TRUE)
pt_set5<-read.table(paste("sets_",gene_def,"_c5_v6.1.sets.out",sep=""),header=TRUE)
gene_names<- getBM(attributes = c("entrezgene","external_gene_name","hgnc_symbol"),filters="entrezgene", values = pt_gene$GENE, mart = ensembl )
pt_gene<-merge(pt_gene,gene_names,by.x="GENE",by.y="entrezgene",all.x=TRUE)
# # read previous round output
# pt_gene_p<-read.table(paste(working_dir3,gene_def,"_gene.genes.out",sep=""),header=TRUE)
# pt_set5_p<-read.table(paste(working_dir3,"sets_",gene_def,"_c5_v6.1.sets.out",sep=""),header=TRUE)
# # read output from fuma
# fuma_genes<-read.table(paste(working_dir_fuma,"magma.genes.out",sep=""),header=TRUE)
## combine
# genes<-merge(pt_gene,pt_gene_p,by=c("GENE","CHR","START","STOP"),suffixes=c("",".bgenie_info"),all=TRUE)
# genes<-merge(genes,fuma_genes,by.x=c("CHR","hgnc_symbol"),by.y=c("CHR","SYMBOL"),suffixes=c("",".FUMA"),all=TRUE)
# sets<-merge(pt_set5,pt_set5_p,by=c("SET","FULL_NAME"),suffixes=c("",".bgenie_info"),all=TRUE)
#----------------------------------------------------------------------


#----------------------------------------------------------------------
genes<-pt_gene
pval="P"

# pval="P_SNPWISE_MEAN"
# pval="P_SNPWISE_TOP1"
#pval="P_JOINT"
#
genes$Q<-p.adjust(genes[,pval],method="fdr")
genes2label=subset(genes,Q<0.05)$hgnc_symbol
subset(genes,Q<0.05)
t_g<-0.05/NROW(genes)
# genes[which(genes[,pval]<=t_g),]

# code chr as numeric for plot
genes$CHR[genes$CHR=="X"]<-"23"
genes$CHR<-as.numeric(genes$CHR)

# plot and save
png(paste("ManhattanGene_",gene_def,"_",p,".png",sep=""),width=490,height=443)
manhattan(genes,chr="CHR",bp="START",p=pval,snp="hgnc_symbol",
          annotatePval=t_g, #highlight=genes2label,
          genomewideline =-log10(t_g),main=gsub("_"," ",p))
dev.off()

png(paste("QQplotGene_",gene_def,"_",p,".png",sep=""))
qq(genes[,pval],main=paste("Q-Q plot of gene-based GWAS p-values\n",gsub("_"," ",p),sep=""))
dev.off()
#----------------------------------------------------------------------
## sets
sets<-pt_set5

sets<-subset(sets,NGENES>10)
t_s<-0.05/NROW(sets)
sets$p.adj<-sets$P*NROW(sets)

sets$Q<-p.adjust(sets$P,method="fdr")
sets$Padj<-p.adjust(sets$P,method="bonferroni")
top_sets<-subset(sets,P<t_s)$FULL_NAME
top_sets_genes<-lapply(top_sets,function(s){
  x<-unlist(strsplit(c5[grep(s,c5$V1),],"\t"))[-c(1:2)]
  x_names<-getBM(attributes = c("entrezgene","external_gene_name","hgnc_symbol"),filters="entrezgene", values = x, mart = ensembl )
  x_genes<-merge(x_names,genes)
  return(x_genes)
})
names(top_sets_genes)<-top_sets

paste(sort(top_sets_genes[[1]]$hgnc_symbol),sep=", ",1,collapse="")


head(sets[with(sets, order(P)),])

# check Guadalupe et al's pathways
guadalupe_gos<-c( 'Steroid hormone receptor activity', 'Steroid binding', 'Steroid biosynthetic process', 'Androgen biosynthetic process', 
'Steroid metabolic process', 'Androgen metabolic process', 'Estrogen metabolic process', 'Steroid hydroxylase activity',
'Estrogen receptor binding', 'Steroid hormone receptor signallingpathway', 'Estrogen receptor signalling pathway',
'Androgen receptor signalling pathway', 'Response to progesterone stimulus', 'Response to estrogen stimulus',
'Response to steroid hormone stimulus' ,'Androgen receptor binding')
guadalupe_gos<-gsub(" ","_",toupper(guadalupe_gos))

guadalupe_keys<-c("")
Guadalupe2015_candidates<-sets[grep(paste(guadalupe_gos,collapse="|"),sets$FULL_NAME),]
top_Guadalupe<-subset(Guadalupe2015_candidates,P<0.05)$FULL_NAME
top_Guadalupe_genes<-lapply(top_Guadalupe,function(s){
  x<-unlist(strsplit(c5[grep(s,c5$V1),],"\t"))[-c(1:2)]
  x_names<-getBM(attributes = c("entrezgene","external_gene_name","hgnc_symbol"),filters="entrezgene", values = x, mart = ensembl )
  x_genes<-merge(x_names,genes)
  return(x_genes)
})
names(top_Guadalupe_genes)<-top_Guadalupe

library(xtable)
sink(paste("MAGMA_sets_",gene_def,"_Guadalupe2015.tex",sep=""))
print(xtable(Guadalupe2015_candidates[,c("FULL_NAME","NGENES","BETA","BETA_STD","SE","P")]),include.rownames = FALSE)
sink()

write.csv(Guadalupe2015_candidates[,c("FULL_NAME","NGENES","BETA","BETA_STD","SE","P")],file=paste("MAGMA_sets_",gene_def,"_Guadalupe2015.csv",sep=""),row.names = FALSE,quote=FALSE)
write.csv(sets[,c("FULL_NAME","NGENES","BETA","BETA_STD","SE","P")],file=paste("MAGMA_sets_",gene_def,".csv",sep=""),row.names = FALSE,quote=FALSE)

#----------------------------------------------------------------------