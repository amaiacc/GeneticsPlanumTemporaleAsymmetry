#----------------------------------------------------------------------
# Load libraries for plotting
#----------------------------------------------------------------------
library(ggplot2)
library(grid)
library(ggbeeswarm)
library(grid.Extra); library(grid)
library(lme4);library(fmsb)

# allen brain packages
# devtools::install_github('oganm/allenBrain')
library(allenBrain)
## if needed
# install.packages("devtools")
## main package
library(devtools)
install.packages('oro.nifti')
# install_github('aaronjfisher/ggBrain',build_vignettes=TRUE)
## to access help pages
library(ggBrain)
help(package=ggBrain)
#
# function to make first letter of a string capital
simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  s2<- paste(toupper(substring(s, 1,1)), substring(s, 2),  sep="", collapse=" ")
  return(s2)
}

library(dplyr)
library(tidyr)

#----------------------------------------------------------------------
# white background 
#----------------------------------------------------------------------
t1<-theme(
  plot.background = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.border = element_rect(colour='grey', fill=NA),
  panel.background = element_blank(),
  axis.line = element_line(size=.4)
)
mytheme2<-theme_bw() + theme(panel.spacing = unit(0, "lines"), 
                             strip.background = element_rect(fill="white"), strip.text = element_text(size=16), 
                             # axis.title.x=element_blank(), axis.ticks.x =element_blank(), # axis.text.x =element_blank(),
                             # axis.text.x = element_text(size=16,colour="black",hjust=1,vjust=.5),
                             # axis.title.y=element_text(size=16),axis.text.y=element_text(size=16,colour="black"),
                             title=element_text(size=16),
                             axis.title=element_text(size=16),axis.text=element_text(size=16,colour="black"),
                             legend.text=element_text(size=16), legend.title =element_text(size=16) ) 
cols=c("#d95f02","#1f78b4","#33a02c","#000000") # '#1b9e77','#d95f02','#7570b3 # "#a6cee3"  ,"#7570b3"
#----------------------------------------------------------------------
args<-commandArgs(TRUE)
# From: P:/lg-dyslexia-exomes/working/analysis/Gene_expression_AllenBrainAtlas
#----------------------------------------------------------------------
# Set working directory, paths and files
#----------------------------------------------------------------------
# Define directories, dependent on  system
if (Sys.info()['sysname']=='Windows') {dir="P://workspaces/"} else {dir="/data/workspaces/lag/workspaces/"}
setwd(paste(dir,"/lg-dyslexia-exomes/working/",sep="")) # this is where I initally donwloaded the data
allen_path<-paste(getwd(),'/expression_data/AllenBrainAtlas/genes',sep="")
output_plots_path<-paste(dir,'lg-ukbiobank/working_data/amaia/genetic_data/PT/allen_brain/',sep="")
output_tables_path<-paste(dir,'lg-ukbiobank/working_data/amaia/genetic_data/PT/allen_brain/',sep="")
#----------------------------------------------------------------------
probes<-read.csv(paste(dir,"/lg-dyslexia-exomes/working/expression_data/AllenBrainAtlas/H0351.1009/probes.csv",sep=""))
# a<-paste(dir,"/lg-dyslexia-exomes/working/working_data/expression/genes_ITIH5.txt",sep="")
a<-paste(dir,"/lg-dyslexia-exomes/working/working_data/expression/genes_PT.txt",sep="")
# cand_genes<-read.table(a,stringsAsFactors=FALSE)
# cand_genes<-c("ITIH5","SLC35E2A","NADK","TMEM52","BOK","C19orf12","PPP1R14A","AC011479.2") # SPINT2" coded as AC011479.2
cand_genes<-c("ITIH5","BOK","BOK-AS1","ING5","DTYMK","AC114730.11") 

#----------------------------------------------------------------------
# Read generated table (including all genes in cand_genes) from candidate_genes_allen_expression.R
#----------------------------------------------------------------------
date<-"2019-02-20"# "2018-11-30" #"2018-10-26"# "2018-06-08"
table_name<-paste(getwd(),"/working_data/expression/",  tail(unlist(strsplit(a,"/")),n=1),"_microarrayExp_AllenBrainAtlas_",date,".csv", sep="")
# table_name<-"working_data/expression/genes_ITIH5.txt_microarrayExp_AllenBrainAtlas_2018-06-08.csv"
all_genes<-read.csv(table_name,quote = "")
table(all_genes$Gene_name)
# PAC 
# Contains a present/absent flag which indicates whether the probe's
#     expression is well above background.  It is set to 1 when both of the
#     following conditions are met.
#     
#         1) The 2-sided t-test p-value is lower than 0.01, (indicating the mean
#            signal of the probe's expression is significantly different from the
# corresponding background).
#         2) The difference between the background subtracted signal and the
# background is significant (> 2.6 * background standard deviation).
#----------------------------------------------------------------------
all_genes_exp<-subset(all_genes,PAC==1)

# average expression of gene per probe and sample
all_genes_exp2<- all_genes_exp  %>% group_by(Gene_name,Structure_acronym,Structure_info) %>%  summarise(avgExp = mean(Expression_PAC)) 
all_genes_avgExp_list<-lapply(as.character(unique(all_genes_exp2$Gene_name)),function(x){
  g<-subset(all_genes_exp2,Gene_name==x)
  maxAcr<-as.character(g$Structure_acronym[which(g$avgExp==max(g$avgExp))])
  maxReg<-as.character(g$Structure_info[which(g$avgExp==max(g$avgExp))])
  minAcr<-as.character(g$Structure_acronym[which(g$avgExp==min(g$avgExp))])
  minReg<-as.character(g$Structure_info[which(g$avgExp==min(g$avgExp))])
  s<-cbind(Gene=x,Region_maxExpr=maxReg,Acr_maxExpr=maxAcr,Region_MinExpr=minReg,Acr_MinExpr=minAcr,t(summary(g$avgExp)))
  return(s)
  })
all_genes_avgExp<-do.call("rbind",all_genes_avgExp_list)

write.csv(all_genes_avgExp,file=paste(output_tables_path,"all_genes_avgExp.csv",sep="/"),row.names = FALSE)

#----------------------------------------------------------------------
# linear model for exp?
#----------------------------------------------------------------------
for (ngene in c("ITIH5","BOK","DTYMK")){
  gene<-subset(all_genes,Gene_name==ngene&PAC==1) # if PAC=0 -> not diff to backaground
  # lm_gene<-lm(Expression~ Probe_id + Subject + Structure_acronym + Hemis, data=subset(gene))
  # lm_gene2<-lm(Expression~ Probe_id + Subject + Structure_acronym*Hemis, data=subset(gene))
  # anova(lm_gene,lm_gene2)
  # 
  # lmm_gene<-lmer(Expression~  + (1|Subject) + Probe_id + Structure_acronym + Hemis, data=subset(gene))
  # summary(lmm_gene)
  # lmm_gene2<-lmer(Expression~ + (1|Subject) + Probe_id + Structure_acronym*Hemis, data=subset(gene))
  # 
  # anova(lmm_gene,lmm_gene2)
  # rm(lmm_gene,lmm_gene2)
  # 
  # # only temporal regions
  # temp<-subset(gene,Region=="Temporal lobe")
  # lm_gene_tmp0<-lm(Expression ~ Probe_id + Subject + Structure_acronym + Hemis, data=temp)
  # lm_gene_tmp<-lm(Expression ~ Probe_id + Subject + Structure_acronym + Hemis, data=temp)
  # lm_gene_tmp2<-lm(Expression ~ Probe_id + Subject + Structure_acronym*Hemis, data=temp)
  # anova(lm_gene_tmp,lm_gene_tmp2)
  # 
  # lmm_gene_tmp<-lmer(Expression~  + (1|Subject) + Probe_id + Structure_acronym + Hemis, data=temp)
  # lmm_gene_tmp2<-lmer(Expression~ + (1|Subject) + Probe_id + Structure_acronym*Hemis, data=temp)
  # anova(lmm_gene_tmp,lmm_gene_tmp2)
  # # clean all models
  # rm(lm_gene_tmp0, lm_gene_tmp,lm_gene_tmp2)
  # rm(lmm_gene_tmp,lmm_gene_tmp2)
  # rm(temp)
  #----------------------------------------------------------------------
  # Plot all regions
  #----------------------------------------------------------------------
  for (r in unique(gene$Region)){
    g_r<-ggplot(data=subset(gene,Region==r),aes(x=as.factor(Probe_id),y=Expression_PAC,colour=factor(Hemis))) + #geom_quasirandom(dodge.width=1)  +  
      geom_boxplot(alpha = 0.2,outlier.shape = NA) +
      facet_wrap(~ Structure_acronym,scales="free_x") + mytheme2 + #coord_cartesian(ylim=c(3,10)) + 
      labs(x=NULL,y="Expression * PAC", title= r ) +
      theme(legend.position="none") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      scale_color_manual(values=cols[2:3]) +
      theme(legend.position="bottom",legend.title=element_blank()) 
    assign(gsub(" ","_",paste(ngene,"_",r,sep="")),g_r)
    rm(g_r)
  }
  
  #----------------------------------------------------------------------
  # candidate region: Planum Temporale
  #----------------------------------------------------------------------
  # # all regions
  # x0<-lmer(Expression~  (1|Subject) + Probe_id + Structure_acronym , data=gene)
  # x1<-lmer(Expression~  (1|Subject) + Probe_id + Structure_acronym + Hemis, data=gene)
  # anova(x0,x1)
  # # check interaction with region, to see if effect in PT is specific
  # x2<-lmer(Expression~  (1|Subject) + Probe_id + Structure_acronym * Hemis, data=gene)
  # anova(x1,x2)
  # pt only
  pt<-subset(gene,Structure_acronym=="PLT")
  z0<-lmer(Expression~  (1|Subject) + Probe_id , data=pt)
  z<-lmer(Expression ~  (1|Subject) + Probe_id  + Hemis, data=pt)
  anova(z0,z)
  summary(z)
  rm(z0,z)
  chisq_v<-round(anova(z0,z)$`Chisq`[2],digits=2)
  p_v<-format(anova(z0,z)$`Pr(>Chisq)`[2],digits=4)
  note<-paste("Chisq(1) =",chisq_v,"; p-value=",p_v,sep="")
  #
  pt$sample<-gsub("sampleH0351.","",pt$Subject)
  pt$Probe_id<-gsub("probe","",pt$Probe_id)
  # combine with probes, to get correct name
  pt<-merge(pt,probes,by.x=c("Probe_id"),by.y="probe_id",stringAsFactors=FALSE,all.x=TRUE)
  
  # plots for planum temporale
  g<-ggplot(data=pt,aes(x=as.factor(probe_name),y=Expression_PAC,colour=factor(Hemis))) + 
    geom_boxplot(fill="white",outlier.colour = NA, position = position_dodge(width=0.9)) +
    geom_quasirandom(dodge.width=0.9, aes(colour=factor(Hemis)),size=3,alpha=0.4)  +
    # geom_text(label=note,position=2,size=3) +???
    # facet_grid(.~ Structure_acronym,scales="free_x") + 
    mytheme2 + 
    coord_cartesian(ylim=c(0,10)) + 
    labs(x=NULL,y="Expression * PAC", title= ngene ) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_color_manual(values=cols[2:3]) +
    theme(legend.position="bottom",legend.title=element_blank()) 
  # g_pt<-ggplot(data=pt,aes(x=as.factor(probe_name),y=Expression_PAC,colour=factor(Hemis))) + 
  #   geom_boxplot(fill="white",outlier.colour = NA, position = position_dodge(width=0.9),alpha=0.5) +
  #   geom_text(aes(label=sample),position=position_jitterdodge(dodge.width=1),size=3,alpha=1,angle=20) +
  #   # facet_grid(.~ Structure_acronym,scales="free_x") + 
  #   mytheme2 + 
  #   coord_cartesian(ylim=c(0,10)) + 
  #   labs(x=NULL,y="Expression * PAC", title= ngene ) +
  #   theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  #   scale_color_manual(values=cols[2:3]) +
  #   theme(legend.position="bottom",legend.title=element_blank())
  
  assign(paste(ngene,"_g",sep=""),g)
  assign(paste(ngene,"_g_pt",sep=""),g_pt)
  # save pt data for gene
  assign(paste(ngene,"_pt_data"),pt)
  
  ggsave(g,file=paste(output_plots_path,ngene,"_PT_summary.png",sep=""),height=6,width=6)
  ggsave(g_pt,file=paste(output_plots_path,ngene,"_PT_summary_samples.png",sep=""),height=6,width=8)
  
  # clean
  rm(g,g_pt,pt)
  rm(gene)
}

rm(ngene)



#----------------------------------------------------------------------
# Plot temporal lobe
#----------------------------------------------------------------------
ITIH5_Temporal_lobe
BOK_Temporal_lobe
DTYMK_Temporal_lobe

# PT 
library(gridExtra)
ITIH5_BOK_pt<-grid.arrange(ITIH5_g + theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) ,
             BOK_g + theme(axis.text.x = element_text(angle = 0, hjust = 0.5)),
             ncol=2,widths=c(3,2))

ITIH5_BOK_DTYMK_pt<-grid.arrange(ITIH5_g + theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) ,
                           BOK_g + theme(axis.text.x = element_text(angle = 0, hjust = 0.5)),
                           DTYMK_g + theme(axis.text.x = element_text(angle = 0, hjust = 0.5)),
                           ncol=3,widths=c(3,2,2))

ggsave(ITIH5_BOK_pt,file=paste(output_plots_path,"ITIH5_BOK","_PT_summary.png",sep=""),width=15,height=7)
ggsave(ITIH5_BOK_DTYMK_pt,file=paste(output_plots_path,"ITIH5_BOK_DTYMK","_PT_summary.png",sep=""),width=20,height=7)

