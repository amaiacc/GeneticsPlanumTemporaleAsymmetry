# Read FUMA's eQTL results and use Blood datasets to define genes of interest within a locus.

subset_name="imagingT1_N18057"
# define snps of interest:
lead_snp<-"rs7420166"
#----------------------------------------------------------------------
# Set working directory, paths and files
#----------------------------------------------------------------------
# Define directories, dependent on  system
if (Sys.info()['sysname']=='Windows') {dir="P://workspaces/"} else {dir="/data/workspaces/lag/workspaces/"}
ld_dir=paste(dir,"lg-ukbiobank/working_data/amaia/genetic_data/PT/",subset_name,"/lead_snps/LD/",sep="")
fuma_dir=paste(dir,"lg-ukbiobank/working_data/amaia/genetic_data/PT/",subset_name,"/FUMA/FUMA_job36828_1KGref/",sep="")
setwd(fuma_dir)
#----------------------------------------------------------------------
# read eqtl data, all nominal associations from FUMA
eqtl<-read.table("eqtl.txt",header=TRUE)
eqtl$CHR<-sapply(strsplit(eqtl$uniqID,":"),"[[",1)
eqtl$BP<-sapply(strsplit(eqtl$uniqID,":"),"[[",2)
table(eqtl$db)
# subset blood tissues to define genes of interest
eqtl_blood<-subset(eqtl,db=="eQTLGen"|db=="BIOSQTL")
eqtl_blood_fdr<-subset(eqtl_blood,FDR<0.05)
table(eqtl_blood_fdr$uniqID,eqtl_blood_fdr$symbol)
# all these snps are in high LD, so most of them detect the same signals
#----------------------------------------------------------------------
# read LD data, to get r2 values with respect to the lead SNP
ld<-read.table(paste(ld_dir,list.files(ld_dir,pattern=paste("0.6r2*.*",lead_snp,sep="")),sep=""),header=TRUE)
#
eqtl_blood_fdr_ld<-merge(eqtl_blood_fdr,ld,by.x=c("CHR","BP"),by.y=c("CHR_B","BP_B"))
table(eqtl_blood_fdr_ld$R2)
eqtl_blood_fdr_ld_single<-subset(eqtl_blood_fdr_ld,R2==max(eqtl_blood_fdr_ld$R2))
write.csv(eqtl_blood_fdr_ld,"blood_FDR_eqtl.csv",row.names = FALSE)
# get genes to check in other tissues
gene_symbols<-eqtl_blood_fdr_ld[,c("gene","symbol")]
gene_symbols<-gene_symbols[-which(duplicated(gene_symbols)),]
write.csv(gene_symbols,"blood_FDR_eqtl_symbol_gencode.csv",row.names = FALSE,quote=FALSE)
#----------------------------------------------------------------------
eqtl_genes<-eqtl[eqtl$symbol %in% genes2check,]
#----------------------------------------------------------------------
# convert to wide format
head(eqtl_blood_fdr_ld)
# two dbs are included (Westra2013 and x), flag them to convert to wide
snps<-eqtl_blood_fdr_ld$SNP_B
genes<-eqtl_blood_fdr_ld$symbol
eqtl_blood_fdr_ld$flag<-NA
for (s in snps){
  for (g in genes){
    t<-subset(eqtl_blood_fdr_ld,symbol==g&SNP_B==s)
    eqtl_blood_fdr_ld[which(eqtl_blood_fdr_ld$symbol==g&eqtl_blood_fdr_ld$SNP_B==s),"flag"]<-1:NROW(t)
    rm(t)
  }
}
rm(s,g)
## add
eqtl_fdr_wide<-eqtl_blood_fdr_ld[,c("SNP_B","CHR","BP","R2","RiskIncAllele","db","flag","symbol","tissue","FDR")] %>% spread(symbol,FDR)
# sort by chromosome position
eqtl_fdr_wide<-eqtl_fdr_wide[order(eqtl_fdr_wide$BP),]

write.csv(eqtl_fdr_wide,"Blood_BIOSQTL_FDR_eqtl_wide.csv",row.names = FALSE)
