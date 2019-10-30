# summarize and plot genetic correlations from the Cognitive and PT analyses, using GCTA, LDSC and SUMHER
library(ggpubr)
library(ggplot2)
library(tidyr)
library(dplyr)
library(reshape2)
library(grid)
library(gridExtra)
library(GGally)
if (!require("gplots")) {
  install.packages("gplots", dependencies = TRUE)
  library(gplots)
}
if (!require("RColorBrewer")) {
  install.packages("RColorBrewer", dependencies = TRUE)
  library(RColorBrewer)
}
if (!require("corrplot")) {
  install.packages("corrplot", dependencies = TRUE)
  library(corrplot)
}

# function to make first letter of a string capital
simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  s2<- paste(toupper(substring(s, 1,1)), substring(s, 2),  sep="", collapse=" ")
  return(s2)
}

options(stringsAsFactors = FALSE)
#------------------------------------------
# plots, some general parameters
#------------------------------------------
colrs=c("#d95f02","#1f78b4","#33a02c","#000000") # '#1b9e77','#d95f02','#7570b3 # "#a6cee3"  ,"#7570b3"
mytheme<-theme_bw() + theme(panel.spacing = unit(0, "lines"), 
                            strip.background = element_rect(fill="white"), strip.text = element_text(size=16), 
                            # axis.title.x=element_blank(), axis.ticks.x =element_blank(), # axis.text.x =element_blank(),
                            # axis.text.x = element_text(size=16,colour="black",hjust=1,vjust=.5),
                            # axis.title.y=element_text(size=16),axis.text.y=element_text(size=16,colour="black"),
                            title=element_text(size=16),
                            axis.title=element_text(size=16),axis.text=element_text(size=16,colour="black"),
                            legend.text=element_text(size=16), legend.title =element_text(size=16) ) 

mytheme2<-theme_bw() + theme(panel.spacing = unit(0, "lines"), 
                            strip.background = element_rect(fill="white"), strip.text = element_text(size=16), 
                            title=element_text(size=16),
                            axis.title=element_text(size=24),axis.text=element_text(size=24,colour="black"),
                            axis.title.x=element_blank(),
                            legend.text=element_text(size=24), legend.title=element_blank(), legend.position="bottom") 
#extract legend
#https://github.com/hadley/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}


#----------------------------------------------------------------------
# define working_dirs
if (Sys.info()['sysname']=='Windows') {dir="P://workspaces/"} else {dir="/data/workspaces/lag/workspaces/"}
#----------------------------------------------------------------------
working_dir=paste(dir,"lg-ukbiobank/working_data/amaia/genetic_data/cognitive/ldsc/",sep="")
#
subset_name="imagingT1_N18057"
root="ukb25465_ukb25468"
pattern_run="_noBioCovs_noAssessmentC"
pattern=paste(root,"imaging",pattern_run,"_",sep="")
#----------------------------------------------------------------------
setwd(working_dir)
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# read summary result form LDSC run (from BGENIE output)  
#----------------------------------------------------------------------
## heritability
tmp<-read.table(file=paste(working_dir,"summary_h2_ldsc.table",sep=""),sep="\t",header=TRUE)
tmp$stat<-"h2"
tmp<-subset(tmp,h2..se.!="")
tmp$h2<-sapply(strsplit(tmp$h2..se.," "),"[[",1)
tmp$h2.SE<-gsub("\\(|\\)","",sapply(strsplit(tmp$h2..se.," "),"[[",2))
tmp$Intercept<-sapply(strsplit(tmp$Intercept..se.," "),"[[",1)
tmp$Intercept.SE<-gsub("\\(|\\)","",sapply(strsplit(tmp$Intercept..se.," "),"[[",2))
tmp$Ratio<-sapply(strsplit(tmp$Ratio..se.," "),"[[",1)
tmp$Ratio.SE<-gsub("\\(|\\)","",sapply(strsplit(tmp$Ratio..se.," "),"[[",2))
tmp<-tmp[,-grep("\\.\\.",colnames(tmp))]
tmp$phenotype<-gsub("_h2.log","",tmp$File)
tmp$phenotype<-gsub("_max","Max",tmp$phenotype)
tmp$pheno<- sapply(strsplit(tmp$phenotype,"_"),"[[",1)
tmp$measure<-gsub("EduYears|FluidInt","",tmp$pheno)
#
tmp$subset<- sapply(strsplit(tmp$phenotype,"_"),"[[",2)
#
tmp$sample<-"total"
tmp$sample[grep("_female",tmp$phenotype)]<-"females"
tmp$sample[grep("_male",tmp$phenotype)]<-"males"
tmp$sample<-factor(tmp$sample,levels=c("total","females","males"))
#
tmp$program<-"LDSC"
tmp$sumstats<-"BGENIE"
#
tmp$estimate<-as.numeric(tmp$h2)
tmp$se<-as.numeric(tmp$h2.SE)
# save
ldsc_h2<-tmp; rm(tmp)
ldsc_h2$measure
h2_p2<-ggplot(data=subset(ldsc_h2), aes(x=subset,y=estimate,width=0.7,alpha=sample,fill=measure)) + geom_bar(aes(y=estimate),position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=estimate-(1.96*se), ymax=estimate+(1.96*se),width=.5), position=position_dodge(.7)) +
  mytheme +
  facet_grid(.~pheno) +
  scale_alpha_discrete( range=c(1,0.4), na.value = 0) + 
  ylab(bquote('Estimate ('*h^2*')')) +
  labs(title="") + 
  scale_fill_hue(l = 50, c = 20) + guides(fill = FALSE)
NULL
#----------------------------------------------------------------------
## genetic correlation
# with GWAS sumstats from publicly available data
tmp<-read.table(paste(working_dir,"summary_p1_p2_GWASgencor_ldsc.table",sep=""))
colnames(tmp)<-c("p1",tmp[1,-1])
tmp<-subset(tmp,p2!="p2")
tmp$p1<-gsub("/data/workspaces/lag/workspaces/lg-multilateral/working/Data/ALSPAC/ldsc/|/data/workspaces/lag/workspaces/lg-ukbiobank/working_data/amaia/genetic_data/release_v3/imp/ldsc/","",gsub("//ldsc/","/ldsc/",tmp$p1))
tmp$p2<-gsub("/data/workspaces/lag/workspaces/lg-multilateral/working/Data/ALSPAC/ldsc/|/data/workspaces/lag/workspaces/lg-ukbiobank/working_data/amaia/genetic_data/release_v3/imp/ldsc/","",gsub("//ldsc/","/ldsc/",tmp$p2))
tmp$p1<-gsub("/data/workspaces/lag/workspaces/lg-ukbiobank/working_data/amaia/genetic_data/release_v3/imp/imagingT1_N18057/ldsc/","",tmp$p1)
tmp$p2<-gsub("/data/workspaces/lag/workspaces/lg-ukbiobank/working_data/amaia/genetic_data/release_v3/imp/imagingT1_N18057/ldsc/","",tmp$p2)
tmp$p1<-gsub("/data/workspaces/lag/workspaces/lg-multilateral/working/Data/GWAS_sumstats/ldsc/","",tmp$p1)
tmp$p2<-gsub("/data/workspaces/lag/workspaces/lg-multilateral/working/Data/GWAS_sumstats/ldsc/","",tmp$p2)
tmp$p1<-gsub("_GWASsumstats.log","",tmp$p1)
tmp$p2<-gsub("_GWASsumstats.log","",tmp$p2)
tmp$p1<-sapply(strsplit(tmp$p1,":|-"),"[[",1)
#
tmp$p1<-gsub("_max","Max",tmp$p1)
tmp$subset<- sapply(strsplit(tmp$p1,"_"),"[[",2)
#
tmp$sample<-"total"
tmp$sample[grep("_female",tmp$p1)]<-"females"
tmp$sample[grep("_male",tmp$p1)]<-"males"
tmp$sample<-factor(tmp$sample,levels=c("total","females","males"))
#
tmp$program<-"LDSC"
tmp$sumstats<-"BGENIE"
#
tmp$estimate<-as.numeric(tmp$rg)
tmp$se<-as.numeric(tmp$se)
#
ldsc_gencor<-tmp; rm(tmp)
ldsc_gencor$p1<-gsub(".sumstats.gz","",ldsc_gencor$p1)
ldsc_gencor$p2<-gsub(".sumstats.gz","",ldsc_gencor$p2)
ldsc_gencor$p2[ldsc_gencor$p2=="Sniekers2017"]<-paste("Sniekers2017","IQ",sep="_")
# for every phenotype, make first letter capital
ldsc_gencor$pheno1<-gsub("_"," ",gsub(":","\n",sapply(ldsc_gencor$p1,simpleCap)))
ldsc_gencor$pheno1<-sapply(strsplit(ldsc_gencor$p1,"_"),"[[",1)
ldsc_gencor$pheno2<-gsub("_"," ",gsub(":","\n",sapply(ldsc_gencor$p2,simpleCap)))
ldsc_gencor$pheno2<-sapply(strsplit(ldsc_gencor$pheno2," "),"[[",2)
ldsc_gencor$pheno2_cat<-"Psychiatric"
w<-which(ldsc_gencor$pheno2=="EA"|ldsc_gencor$pheno2=="EA3"|
           ldsc_gencor$pheno2=="CP"|ldsc_gencor$pheno2=="IQ")
ldsc_gencor$pheno2_cat[w]<-"Cognitive"
table(ldsc_gencor$pheno2_cat)
# save
write.csv(ldsc_gencor,file=paste(working_dir,"summary_GWASgencor_ldsc.csv",sep=""),row.names = FALSE,quote = FALSE)



# with GWAS sumstats from publicly available data
pt_dir="/data/workspaces/lag/workspaces/lg-ukbiobank/working_data/amaia/genetic_data/PT/imagingT1_N18057/ldsc/"
tmp<-read.table(paste(working_dir,"summary_p1_p2_gencor_ldsc.table",sep=""))
colnames(tmp)<-c("p1",tmp[1,-1])
tmp<-subset(tmp,p2!="p2")
tmp$p1<-gsub("/data/workspaces/lag/workspaces/lg-multilateral/working/Data/ALSPAC/ldsc/|/data/workspaces/lag/workspaces/lg-ukbiobank/working_data/amaia/genetic_data/release_v3/imp/ldsc/","",gsub("//ldsc/","/ldsc/",tmp$p1))
tmp$p2<-gsub("/data/workspaces/lag/workspaces/lg-multilateral/working/Data/ALSPAC/ldsc/|/data/workspaces/lag/workspaces/lg-ukbiobank/working_data/amaia/genetic_data/release_v3/imp/ldsc/","",gsub("//ldsc/","/ldsc/",tmp$p2))
tmp$p1<-gsub("/data/workspaces/lag/workspaces/lg-ukbiobank/working_data/amaia/genetic_data/release_v3/imp/imagingT1_N18057/ldsc/","",tmp$p1)
tmp$p2<-gsub("/data/workspaces/lag/workspaces/lg-ukbiobank/working_data/amaia/genetic_data/release_v3/imp/imagingT1_N18057/ldsc/","",tmp$p2)
tmp$p1<-gsub(pt_dir,"",gsub("/data/workspaces/lag/workspaces/lg-multilateral/working/Data/GWAS_sumstats/ldsc/","",tmp$p1))
tmp$p2<-gsub(pt_dir,"",gsub("/data/workspaces/lag/workspaces/lg-multilateral/working/Data/GWAS_sumstats/ldsc/","",tmp$p2))
tmp$p1<-gsub("_GWASsumstats.log|_PTsumstats.log","",tmp$p1)
tmp$p2<-gsub("_GWASsumstats.log|_PTsumstats.log","",tmp$p2)
tmp$p1<-sapply(strsplit(tmp$p1,":|-"),"[[",1)
#
tmp$p1<-gsub("_max","Max",tmp$p1)
tmp$subset<- sapply(strsplit(tmp$p1,"_"),"[[",2)
#
tmp$sample<-"total"
tmp$sample[grep("_female",tmp$p1)]<-"females"
tmp$sample[grep("_male",tmp$p1)]<-"males"
tmp$sample<-factor(tmp$sample,levels=c("total","females","males"))
#
tmp$program<-"LDSC"
tmp$sumstats<-"BGENIE"
#
tmp$estimate<-as.numeric(tmp$rg)
tmp$se<-as.numeric(tmp$se)
#
ldsc_gencor<-tmp; rm(tmp)
ldsc_gencor$p1<-gsub(".sumstats.gz","",ldsc_gencor$p1)
ldsc_gencor$p2<-gsub(".sumstats.gz","",ldsc_gencor$p2)
ldsc_gencor$p2[ldsc_gencor$p2=="Sniekers2017"]<-paste("Sniekers2017","IQ",sep="_")
ldsc_gencor$p2<-gsub("_females|_males","",ldsc_gencor$p2)
# for every phenotype, make first letter capital
ldsc_gencor$pheno1<-gsub("_"," ",gsub(":","\n",sapply(ldsc_gencor$p1,simpleCap)))
ldsc_gencor$pheno1<-sapply(strsplit(ldsc_gencor$p1,"_"),"[[",1)
ldsc_gencor$pheno2<-gsub("_"," ",gsub(":","\n",sapply(ldsc_gencor$p2,simpleCap)))
ldsc_gencor$pheno2<-sapply(strsplit(ldsc_gencor$pheno2," "),"[[",2)
ldsc_gencor$pheno2_cat<-"Psychiatric"
w<-which(ldsc_gencor$pheno2=="EA"|ldsc_gencor$pheno2=="EA3"|
           ldsc_gencor$pheno2=="CP"|ldsc_gencor$pheno2=="IQ")
ldsc_gencor$pheno2_cat[w]<-"Cognitive"
ldsc_gencor$pheno2_cat[grep("Planum_Temporale",ldsc_gencor$p2)]<-"PT"
table(ldsc_gencor$pheno2_cat)
rm(w)
#
ldsc_gencor$pheno2[ldsc_gencor$pheno2_cat=="PT"]<- gsub("_.adj",".adj",gsub("_$|^_","",gsub("Planum_Temporale|Volume_of_grey_matter_in_|VOLUME_","",ldsc_gencor$p2[ldsc_gencor$pheno2_cat=="PT"])))
# save
# write.csv(ldsc_gencor,file=paste(working_dir,"summary_GWASgencor_all_cogn_ldsc.csv",sep=""),row.names = FALSE,quote = FALSE)


# plot this, select relevant phenos
t_gwas<-subset(ldsc_gencor,pheno2_cat!="PT" &
               p2!="Okbay2016_EA" & p2!="Sniekers2017_IQ" & p2!="2017_OCD" & p2!="Lee2018_CP" &
                 p2!="Sklar2012_BIP" & p2!="Wray2018_MDD")
t_gwas<-t_gwas[with(t_gwas,order(pheno2_cat)),]
levs=c("IQ","EA3","ASD","ADHDeur","SCZ")
levs2=c("Intelligence","EA","ASD","ADHD","SCZ")
t_gwas$pheno2<-factor(t_gwas$pheno2,levels=levs)
levels(t_gwas$pheno2)
levels(t_gwas$pheno2)<-levs2

ldsc_gencorGWAS_barplot<-ggplot(data=t_gwas,color="black",aes(x=subset,y=estimate,width=0.7,fill=pheno1,alpha=sample)) +  #fill=program
  geom_bar(,position=position_dodge(), stat="identity") +
  # 95% CI
  geom_errorbar(aes(ymin=estimate-(1.96*se), ymax=estimate+(1.96*se),width=.5), position=position_dodge(.7)) +
  facet_grid(pheno1~pheno2) + #,scales="free_y"
  scale_alpha_discrete( range=c(1,0.4), na.value = 0) + 
  mytheme +
  # theme(axis.text.x = element_blank(), axis.title.x=element_blank(),legend.title=element_blank()) + 
  ylab(bquote('Estimate ('*rho*')')) +
  scale_fill_hue(l = 50, c = 20) + guides(fill = FALSE) +
  labs(title="") 

# plot this, select relevant phenos
t_gwas2<-subset(ldsc_gencor,pheno2_cat=="PT")
t_gwas2<-t_gwas2[with(t_gwas2,order(pheno2_cat)),]
t_gwas2$run<-""
t_gwas2$run[grep("adjTBV",t_gwas2$pheno2)]<-"adjTBV"
t_gwas2$pheno2<-gsub(".adjTBV","",t_gwas2$pheno2)
levs=c("AI","left","right")
levs2=c("AI","L","R")
t_gwas2$pheno2<-factor(t_gwas2$pheno2,levels=levs)
levels(t_gwas2$pheno2)
levels(t_gwas2$pheno2)<-levs2

ldsc_gencorPT_barplot<-ggplot(data=subset(t_gwas2,run==""),color="black",aes(x=subset,y=estimate,width=0.7,fill=pheno2,alpha=sample)) +  #fill=program
  geom_bar(,position=position_dodge(), stat="identity") +
  # 95% CI
  geom_errorbar(aes(ymin=estimate-(1.96*se), ymax=estimate+(1.96*se),width=.5), position=position_dodge(.7)) +
  facet_grid(pheno1~pheno2) + #,scales="free_y"
  scale_alpha_discrete( range=c(1,0.4), na.value = 0) + 
  mytheme +
  theme(axis.title.x=element_blank(),legend.title=element_blank()) + 
  ylab(bquote('Estimate ('*rho*')')) +
  # scale_fill_hue(l = 50, c = 20) + guides(fill = FALSE) +
  scale_fill_manual(values=colrs) +
  labs(title="") 



ldsc_gencorPT_EA_barplot<-ggplot(data=subset(t_gwas2,pheno1=="EduYearsMax"),color="black",aes(x=program,y=estimate,width=0.7,fill=pheno2,alpha=sample)) +  #fill=program
  geom_bar(,position=position_dodge(), stat="identity") +
  # 95% CI
  geom_errorbar(aes(ymin=estimate-(1.96*se), ymax=estimate+(1.96*se),width=.5), position=position_dodge(.7)) +
  facet_grid(run~subset) + #,scales="free_y"
  scale_alpha_discrete( range=c(1,0.4), na.value = 0) + 
  mytheme +
  theme(axis.title.x=element_blank(),legend.title=element_blank()) +
  ylab(bquote('Estimate ('*rho*')')) +
  scale_fill_manual(values=colrs) +
  labs(title="Sumstats") 

ldsc_gencorPT_FI_barplot<-ggplot(data=subset(t_gwas2,pheno1!="EduYearsMax"),color="black",aes(x=program,y=estimate,width=0.7,fill=pheno2,alpha=sample)) +  #fill=program
  geom_bar(,position=position_dodge(), stat="identity") +
  # 95% CI
  geom_errorbar(aes(ymin=estimate-(1.96*se), ymax=estimate+(1.96*se),width=.5), position=position_dodge(.7)) +
  facet_grid(run~subset) + #,scales="free_y"
  scale_alpha_discrete( range=c(1,0.4), na.value = 0) + 
  mytheme +
  theme(axis.title.x=element_blank(),legend.title=element_blank()) +
  ylab(bquote('Estimate ('*rho*')')) +
  scale_fill_manual(values=colrs) +
  labs(title="FluidInt ~ PT") 



#----------------------------------------------------------------------


