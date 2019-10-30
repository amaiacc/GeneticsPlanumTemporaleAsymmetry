## eQTL analysis
# query data from Braineac (http://www.braineac.org/)

library(ggplot2)
library(grid);library(gridExtra)
options(stringsAsFactors = FALSE)

subset_name="imagingT1_N18057"
#----------------------------------------------------------------------
# Set working directory, paths and files
#----------------------------------------------------------------------
if (Sys.info()['sysname']=='Windows') {dir="P://workspaces/"} else {dir="/data/workspaces/lag/workspaces/"}
data_dir=paste(dir,"/lg-ukbiobank/working_data/amaia/genetic_data/PT/braineac/data/",sep="")
out_dir=paste(dir,"/lg-ukbiobank/working_data/amaia/genetic_data/PT/braineac/output/",sep="")
ld_dir=paste(dir,"/lg-ukbiobank/working_data/amaia/genetic_data/PT/",subset_name,"/lead_snps/LD/",sep="")
fuma_dir=paste(dir,"lg-ukbiobank/working_data/amaia/genetic_data/PT/",subset_name,"/FUMA/FUMA_job36828_1KGref/",sep="")
# read braineac data for genes of interest, defined by eQTLs within blood
setwd(data_dir)
#----------------------------------------------------------------------

#----------------------------------------------------------------------
tissue_int="TCTX" # temporal cortex as tissue of interest
# define genes of interest: # or from file...
blood_eqtls<-read.csv(paste(fuma_dir,"blood_FDR_eqtl.csv",sep=""))
genes<-unique(blood_eqtls$symbol)
# define snps of interest:
lead_snp<-"rs7420166"
snps<-read.table(paste(ld_dir,list.files(ld_dir,pattern=paste("0.6r2*.*",lead_snp,sep="")),sep=""),header=TRUE)$SNP_B
#----------------------------------------------------------------------
# Run per gene: read data and generate plots
for (gene in genes){
  if (length(list.files(pattern=gene))>0){
  #----------------------------------------------------------------------
  # load Braineac data for gene of interest
  gene_dir<-list.files(pattern=gene)
  rda_file<-paste(gene_dir,list.files(gene_dir,pattern="rda"),sep="/")
  load(rda_file)
  rm(gene_dir,rda_file)
  # combine all expression data into a data.frame: # each col is a probe, each row is a sample
  gene_expr<-as.data.frame(do.call("rbind",
                     lapply(names(expr),function(x){
                        tmp<-as.data.frame(t(expr[[x]]))
                        tmp$tissue<-x
                        tmp$sample<-row.names(tmp)
                        return(tmp)
                      })
              ))
  
  #------------------------------
  # find snps of interest
  w<-match(snps,markers.info$rsid)
  if (length(w)>0){
    w1<-w[!is.na(w)]
    w2<-snps[which(is.na(w))]
  }
  if (length(w2)>0){
   print("The following SNPs were not found:")
   print(w2)
  }
  # subset to snps of interest
  if (length(w1)>0){
    snps_info<-markers.info[w1,]
  }
  rm(w1,w2,w)
  snps_genos<-as.data.frame(t(markers[rownames(snps_info),]))
  # add additive genotypes:
  snps_genos[,snps_info$rsid]<-apply(snps_genos,2,function(x) as.factor(round(x)))
  # add sample as column
  snps_genos$sample<-row.names(snps_genos)
  #------------------------------
  # Merge expresion data and genotypes for snps of interest
  t<-merge(snps_genos,gene_expr,by=c("sample"),all=TRUE)
  # get name for total transcipt column
  t_col<-colnames(t)[grep("^t",colnames(t))]
  t_col<-t_col[grep("tissue",t_col,invert=TRUE)]
  # subset to tissues of interest:
  ## aveALL and temporal cortex:
  tissue_int<-"TCTX"
  t2<-subset(t,tissue=="aveALL"|tissue==tissue_int)
  #------------------------------
  # make plot per snp
  # sapply(snps, function(snp){
  for (snp in snps_info$rsid){
    # get alleles for snp
    a1<-subset(snps_info,rsid==snp)$Al1
    a2<-subset(snps_info,rsid==snp)$Al2
    if (subset(snps_info,rsid==snp)$Freq1>0.5) {
      genotypes<-c(paste(rep(a2,2),collapse=""),paste(a2,a1,collapse=""),paste(rep(a1,2),collapse=""))
    } else {
      genotypes<-c(paste(rep(a1,2),collapse=""),paste(a1,a2,collapse=""),paste(rep(a2,2),collapse=""))
    }
    
    # run lm
    formula1<-paste(t_col,"~ tissue +",snp," + tissue*",snp,collapse="")
    formula2<-paste(t_col,"~ ",snp,collapse="")
    lm_eqtl<-do.call("lm",list (as.formula(formula1),data=t))
    lm1_eqtl<-do.call("lm",list (as.formula(formula2),data=subset(t,tissue=="aveALL")))
    lm2_eqtl<-do.call("lm",list (as.formula(formula2),data=subset(t,tissue==tissue_int)))
    # extract info to add into plot as text
    lm_pvals<-c(anova(lm1_eqtl)[snp,"Pr(>F)"],anova(lm2_eqtl)[snp,"Pr(>F)"])
    lm_pvals<-paste("Pr(F) = ",format(lm_pvals,digits=3))
    lm_F<-c(anova(lm1_eqtl)[snp,"F value"],anova(lm2_eqtl)[snp,"F value"])
    lm_F<-paste("F value = ",format(lm_F,digits=3),sep="")
    summary_lm<-as.data.frame(cbind(pval=lm_pvals,F=lm_F,tissue=c("aveALL",tissue_int)))
    summary<-rbind(cbind(anova(lm1_eqtl),tissue="aveALL"),cbind(anova(lm1_eqtl),tissue=tissue_int))
    summary$snp<-snp
    summary$gene<-gene
    # get minimum value, to print text in plot
    min<-min(t[,t_col],na.rm=TRUE)
    # plot
    p<-ggplot(t2,aes_string(x=snp,y=t_col)) + 
      geom_boxplot(varwidth = TRUE) +
      # geom_text(data=cors_ai, aes(label=paste("r=", cor, sep="")), x=0.05, y=(0.22),size=5,color="red") +
      geom_text(data=summary_lm, aes(label=F),x=1,y=min+0.5) +
      geom_text(data=summary_lm, aes(label=pval),x=1,y=min+0.4) +
      scale_x_discrete(labels=genotypes) + 
      scale_y_continuous(name=expression("Expression level - log"[2])) +
      labs(title=paste(gene," (",t_col,")",sep="")) +
      facet_wrap(~tissue) + theme_classic()
    
    # 
    assign(paste("plot_eQTL_",snp,sep=""),p)
    assign(paste("summary_",gene,"_",snp,sep=""),summary)
    # clean intermediate objects
    rm(a1,a2,genotypes)
    rm(formula1,formula2,lm_eqtl,lm1_eqtl,lm2_eqtl,lm_pvals,lm_F,summary_lm)
    rm(min)
    rm(p,summary)
    
  }
  rm(snp)
  #------------------------------
  # arrange and save plots
  lst=lapply(ls(pattern="plot_eQTL"),get)
  png(file=paste(out_dir,"/plots/",gene,"_eQTL.png",sep=""),width=200*length(snps_info$rsid),height=150*length(snps_info$rsid))
  do.call(grid.arrange,lst)
  dev.off()
  rm(list=ls(pattern="plot_eQTL"))
  rm(lst)
  #------------------------------
  # clean other intermediate objects
  rm(expr,expr.info,markers.info,markers)
  rm(snps_info,snps_geno)
  rm(gene_expr,snps_genos)
  rm(t,t2)
  }
}
#------------------------------
# save
# combine all summary tables, all genes
lst_s<-lapply(ls(pattern="summary_"),get)
summary_all<-do.call(rbind,lst_s)
summary_all$p<- summary_all$`Pr(>F)`
subset(summary_all,p<0.1)
