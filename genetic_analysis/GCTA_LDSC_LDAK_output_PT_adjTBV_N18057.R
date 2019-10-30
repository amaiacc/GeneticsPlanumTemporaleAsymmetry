# summarize and plot genetic correlations from the PT analyses, using GCTA, LDSC and SUMHER
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
if (Sys.info()['sysname']=='Windows') {
  # dir="P://workspaces/"
  dir="\\\\data/lag/workspaces/"
  } else {dir="/data/workspaces/lag/workspaces/"}
#----------------------------------------------------------------------
working_dir=paste(dir,"lg-ukbiobank/working_data/amaia/",sep="")
#
subset_name="imagingT1_N18057"
out_dir=paste(paste(dir,"lg-ukbiobank/working_data/amaia/genetic_data/PT/",subset_name,"/",sep=""))
# gcta results
working_dir_gcta=paste(dir,"lg-ukbiobank/working_data/amaia/genetic_data/release_v2/cal/h2/",subset_name,"/output/PT/",sep="")
# sumstats from bgenie
working_dir_ldsc2=paste(out_dir,"/ldsc/",sep="")
working_dir_ldak2=paste(out_dir,"/sumher/",sep="")
#

root="ukb25465_ukb25468"
pattern_run="_noBioCovs_noAssessmentC"
pattern=paste(root,"imaging",pattern_run,"_",sep="")
#
pheno_dir=paste(working_dir,"/pheno_files/genetic_v2/",root,"/summary_phenotypes/volume_PT",pattern_run,"/",sep="")

#----------------------------------------------------------------------
setwd(working_dir)
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# read summary result form GCTA run:
#----------------------------------------------------------------------
## heritability and genetic correlation
h2_files<-list.files(working_dir_gcta,pattern=paste(root,"_imaging",pattern_run,"*.*.csv",sep=""))
GCTAsummary1<-read.csv(paste(working_dir_gcta,h2_files[1],sep=""))
GCTAsummary1<-subset(GCTAsummary1,region!="")
GCTAsummary1$adj<-""
# adjTBV
GCTAsummary2<-read.csv(paste(working_dir_gcta,h2_files[2],sep=""))
GCTAsummary2$phenotype<-gsub("_PT","",GCTAsummary2$phenotype)
GCTAsummary2$region<-GCTAsummary2$measure<-""
GCTAsummary2$region[grep("Planum_Temporale",GCTAsummary2$file)]<-"Planum Temporale"
GCTAsummary2$region[grep("totalBV$",GCTAsummary2$phenotype)]<-"totalBV"
GCTAsummary2$measure[grep("left",GCTAsummary2$file)]<-"L"
GCTAsummary2$measure[grep("right",GCTAsummary2$file)]<-"R"
GCTAsummary2$measure[grep("AI",GCTAsummary2$file)]<-"AI"
GCTAsummary2<-subset(GCTAsummary2,region!="")
GCTAsummary2$adj<-"adjTBV"
# # combine
GCTAsummary<-merge(GCTAsummary1,GCTAsummary2,all=TRUE,stringAsFactors=FALSE)
# GCTAsummary<-GCTAsummary1
GCTAsummary$pheno<-gsub("_"," ",gsub(":","\n",sapply(GCTAsummary$region,simpleCap)))
GCTAsummary$program<-"GCTA"
GCTAsummary$sumstats<-'-'
GCTAsummary$region<-gsub("_"," ",gsub(":","\n",sapply(GCTAsummary$region,simpleCap)))
GCTAsummary<-subset(GCTAsummary,stat=="h2")
GCTAsummary$region[grep("Cerebellum",GCTAsummary$file)]<-paste(GCTAsummary$region[grep("Cerebellum",GCTAsummary$file)]," Cerebellum",sep="")


rm(GCTAsummary1,GCTAsummary2)
#----------------------------------------------------------------------


#----------------------------------------------------------------------
# read summary result form LDAK, from BGENIE output
#----------------------------------------------------------------------
## heritability
tmp<-read.table(file=paste(working_dir_ldak2,"summary_h2_estimates_sumher.table",sep=""),sep="\t",header=TRUE)
tmp$phenotype<-gsub(".gcta|.ldak|.hers|.adjTBV","",tmp$File)
tmp$phenotype<-gsub("Weighted_mean","WeightedMean",tmp$phenotype)
tmp$region<-gsub("tract_","",tmp$phenotype)
tmp$type<- sapply(strsplit(tmp$phenotype,"_"),"[[",2)
tmp$pheno<-gsub("_"," ",gsub(":","\n",sapply(tmp$region,simpleCap)))
tmp$measure<-"NotLat"
tmp$measure[grep("AI",tmp$region)]<-"AI"
tmp$measure[grep("_left|left",tmp$region)]<-"L"
tmp$measure[grep("_right|right",tmp$region)]<-"R"
table(tmp$measure)
tmp$stat<-"h2"
# clean names
tmp$region<-gsub(paste(paste(unique(tmp$type),collapse="_|"),"_",sep=""),"",tmp$region)
tmp$region<-gsub(paste(paste(unique(tmp$type),collapse="_|"),"_",sep=""),"",gsub("_"," ",gsub("_left|_right|AI_","",tmp$region)))
tmp$region<-gsub("Weighted","",gsub("Mean in | on skeleton|Volume grey matter in ", "", tmp$region))
tmp$region<-gsub("_"," ",gsub(":","\n",sapply(tmp$region,simpleCap)))
tmp$pheno<-gsub("_"," ",gsub(":","\n",sapply(tmp$region,simpleCap)))
# 
tmp$model<-sapply(strsplit(gsub(".adjTBV","",tmp$File),"\\."),"[[",2)
tmp$program<-paste("SumHer","-",tmp$model,sep="")
tmp$sumstats<-'BGENIE'
# extra cols to match GCTA format
tmp$volume<-""
tmp$volume[grep("volume",tolower(tmp$File))]<-"GM"
#tmp$type[tmp$volume=="GM"&]<-"Cortex" # need to fix this!
tmp$type[grep("Cerebell|Pallid|Thalam|Hippocampus|Ventral Striatum|Putam|Caudate|Pallidum|Amygdala",tmp$region)]<-"Subcortical"
tmp$type[grep("Cerebell|Pallid|Thalam|Hippocampus|Ventral Striatum|Putam|Caudate|Pallidum|Amygdala",tmp$region,invert=TRUE)]<-"Cortex"
# calculate se from sd?
#tmp$se<-tmp$h2_SD # *kontuz hemen!*

# save
LDAKsummary_h2<-tmp; rm(tmp)


# genetic correlation from LDAK
tmp<-read.table(file=paste(working_dir_ldak2,"/GenCor_GWASes/summary_rg_estimates_sumher.table",sep=""),sep="\t",header=TRUE)
tmp$phenotype<-gsub(".gcta|.ldak|.hers","",tmp$File)
tmp$phenotype<-gsub("Weighted_mean","WeightedMean",tmp$phenotype)
tmp$region<-sapply(strsplit(tmp$phenotype,"\\/"),"[[",1)
tmp$pheno<-gsub("_"," ",gsub(":","\n",sapply(tmp$region,simpleCap)))
tmp$measure<-"NotLat"
tmp$measure[grep("AI",tmp$region)]<-"AI"
tmp$measure[grep("_left|left",tmp$region)]<-"L"
tmp$measure[grep("_right|right",tmp$region)]<-"R"
table(tmp$measure)
tmp$stat<-"rg"
# clean names
tmp$region<-gsub(paste(paste(unique(tmp$type),collapse="_|"),"_",sep=""),"",tmp$region)
tmp$region<-gsub(paste(paste(unique(tmp$type),collapse="_|"),"_",sep=""),"",gsub("_"," ",gsub("_left|_right|AI_","",tmp$region)))
tmp$region<-gsub("Weighted","",gsub("Mean in | on skeleton|Volume grey matter in ", "", tmp$region))
tmp$region<-gsub("_"," ",gsub(":","\n",sapply(tmp$region,simpleCap)))
tmp$pheno<-gsub("_"," ",gsub(":","\n",sapply(tmp$region,simpleCap)))
# 
tmp$model<-gsub("\\.","",sapply(strsplit(tmp$File,"cor"),"[[",2))
tmp$program<-paste("SumHer","-",tmp$model,sep="")
tmp$GC<-"yes"
tmp$GC[grep("noGC",tmp$program)]<-"no"
tmp$program<-gsub("_noGC","",tmp$program)
#
tmp$sumstats<-'BGENIE'
# extra cols to match GCTA format
tmp$volume<-""
tmp$volume[grep("volume",tolower(tmp$File))]<-"GM"
#tmp$type[tmp$volume=="GM"&]<-"Cortex" # need to fix this!
# tmp$type[grep("Cerebell|Pallid|Thalam|Hippocampus|Ventral Striatum|Putam|Caudate|Pallidum|Amygdala",tmp$region)]<-"Subcortical"
# tmp$type[grep("Cerebell|Pallid|Thalam|Hippocampus|Ventral Striatum|Putam|Caudate|Pallidum|Amygdala",tmp$region,invert=TRUE)]<-"Cortex"

# calculate se from sd?
LDAKsummary_gencor<-tmp; rm(tmp)
LDAKsummary_gencor$p1<-sapply(strsplit(LDAKsummary_gencor$File,"/"),"[[",1)
LDAKsummary_gencor$p2<-gsub("_corgcta|_corldak|.cors|AI_VOLUME_Planum_Temporale_|Volume_of_grey_matter_in_Planum_Temporale_|right_|left_|left.adjTBV_|right.adjTBV_|Planum_Temporale.adjTBV_","",sapply(strsplit(LDAKsummary_gencor$File,"/"),"[[",2))
LDAKsummary_gencor$p2<-gsub("_noGC","",LDAKsummary_gencor$p2)
LDAKsummary_gencor$p2[LDAKsummary_gencor$p2=="Sniekers2017"]<-paste("Sniekers2017","IQ",sep="_")
table(LDAKsummary_gencor$p2)
# for every phenotype, make first letter capital
LDAKsummary_gencor$pheno1<-gsub("_"," ",gsub(":","\n",sapply(LDAKsummary_gencor$p1,simpleCap)))
LDAKsummary_gencor$pheno2<-gsub("_"," ",gsub(":","\n",sapply(LDAKsummary_gencor$p2,simpleCap)))
LDAKsummary_gencor$pheno2<-gsub("AI VOLUME |Volume of grey matter in | left| right","",LDAKsummary_gencor$pheno2)
LDAKsummary_gencor$pheno2<-sapply(strsplit(LDAKsummary_gencor$pheno2," "),"[[",2)
#
LDAKsummary_gencor$pheno2_cat<-"Psychiatric"
w<-which(LDAKsummary_gencor$pheno2=="EA"|LDAKsummary_gencor$pheno2=="EA3"|
           LDAKsummary_gencor$pheno2=="CP"|LDAKsummary_gencor$pheno2=="IQ")
LDAKsummary_gencor$pheno2_cat[w]<-"Cognitive"
table(LDAKsummary_gencor$pheno2_cat)
rm(w)
#
LDAKsummary_gencor$region<-gsub("AI VOLUME |Volume of grey matter in | left| right","",LDAKsummary_gencor$pheno1)
LDAKsummary_gencor$pheno2<-gsub("AI VOLUME |Volume of grey matter in |left |right ","",LDAKsummary_gencor$pheno2)
LDAKsummary_gencor$measure<-"NotLat"
LDAKsummary_gencor$measure[grep("AI",LDAKsummary_gencor$pheno1)]<-"AI"
LDAKsummary_gencor$measure[grep("left",LDAKsummary_gencor$pheno1)]<-"L"
LDAKsummary_gencor$measure[grep("right",LDAKsummary_gencor$pheno1)]<-"R"
LDAKsummary_gencor$adj<-""
LDAKsummary_gencor$adj[grep("adjTBV",LDAKsummary_gencor$pheno1)]<-"adjTBV"
LDAKsummary_gencor$rg<-as.numeric(LDAKsummary_gencor$rg)
LDAKsummary_gencor$se<-as.numeric(LDAKsummary_gencor$rg_SD)
LDAKsummary_gencor$p<-pchisq((LDAKsummary_gencor$rg/LDAKsummary_gencor$se)^2,1,lower.tail = FALSE)
LDAKsummary_gencor$p<-as.numeric(LDAKsummary_gencor$p)
LDAKsummary_gencor$sig<-""
LDAKsummary_gencor$sig[LDAKsummary_gencor$p<0.1&LDAKsummary_gencor$p>=0.05]<-"."
LDAKsummary_gencor$sig[LDAKsummary_gencor$p<0.05]<-"*"
write.csv(LDAKsummary_gencor,file=paste(working_dir_ldak2,"summary_GWASgencor_sumher.csv",sep=""),row.names = FALSE,quote = FALSE)
#LDAKsummary_gencor$sig[LDAKsummary_gencor$p<0.01]<-"**"
LDAKsummary_gencor$sig<-factor(LDAKsummary_gencor$sig,levels=c("",".","*"))
# plot this, select relevant phenos
table(LDAKsummary_gencor$p2)
t_gwas<-subset(LDAKsummary_gencor,
               p2!="Okbay2016_EA" & p2!="Sniekers2017_IQ" & p2!="2017_OCD" & p2!="Lee2018_CP" &
                p2!="Sklar2012_BIP" & p2!="Wray2018_MDD")
t_gwas<-t_gwas[with(t_gwas,order(pheno2_cat)),]
table(t_gwas$p2)
levs=c("IQ","EA3","ASD","ADHDeur","SCZ")
levs2=c("Intelligence","EA","ASD","ADHD","SCZ")

t_gwas$pheno2<-factor(t_gwas$pheno2,levels=levs)
levels(t_gwas$pheno2)<-levs2

ldak_gencor_barplot_gc<-ggplot(data=t_gwas,color="black",aes(x=region,y=rg,width=0.7,fill=measure,alpha=GC)) +  #fill=program,,alpha=sig
  geom_bar(aes(y=rg),position=position_dodge(), stat="identity") +
  # geom_errorbar(data=t2,aes(ymin=CI95lower, ymax=CI95upper),
  # width=.2, position=position_dodge(.9)) + 
  scale_alpha_discrete( range=c(0.3,1), na.value = 0) +
  # geom_errorbar(data=t_gwas,aes(ymin=rg-se, ymax=rg+se,width=.5), position=position_dodge(.7)) +
  geom_errorbar(data=t_gwas,aes(ymin=rg-(1.96*se), ymax=rg+(1.96*se),width=.5), position=position_dodge(.7)) +
  # coord_cartesian(ylim = c(-0.20,0.10))+
  scale_fill_manual(values=colrs) +
  facet_grid(program~pheno2) + #,scales="free_y"
  mytheme + theme(axis.text.x = element_blank(), axis.title.x=element_blank(),legend.title=element_blank()) + 
  ylab(bquote('Estimate ('*rho*')')) +
  labs(title="SumHer") 

ldak_gencor_barplot<-ggplot(data=subset(t_gwas,GC=="yes"&adj==""),color="black",aes(x=region,y=rg,width=0.7,fill=measure)) +  #fill=program,,alpha=sig
  geom_bar(aes(y=rg),position=position_dodge(), stat="identity") +
  geom_errorbar(data=subset(t_gwas,GC=="yes"&adj==""),aes(ymin=rg-(1.96*se), ymax=rg+(1.96*se),width=.5), position=position_dodge(.7)) +
  coord_cartesian(ylim = c(-0.25,0.50))+
  scale_fill_manual(values=colrs) +
  facet_grid(program~pheno2) + #,scales="free_y"
  mytheme + theme(axis.text.x = element_blank(), axis.title.x=element_blank(),legend.title=element_blank()) + 
  ylab(bquote('Estimate ('*rho*')')) +
  labs(title="SumHer") 
  # labs(title="Genetic correlation (SumHer) with disorders") 
ldak_gencor_barplot
ggsave(ldak_gencor_barplot,file=paste(working_dir_ldak2,"summary_GWASgencor_sumher_GCTAandLDAK.pdf",sep=""),width=10,height=10)

ldak_gencor_barplot3<-ggplot(data=subset(t_gwas,GC=="yes"),color="black",aes(x=region,y=rg,width=0.7,fill=measure,alpha=adj)) +  #fill=program,,alpha=sig
  geom_bar(aes(y=rg),position=position_dodge(), stat="identity") +
  scale_alpha_discrete( range=c(1,0.5), na.value = 0) +
  geom_errorbar(data=subset(t_gwas,GC=="yes"),aes(ymin=rg-(1.96*se), ymax=rg+(1.96*se),width=.5), position=position_dodge(.7)) +
  coord_cartesian(ylim = c(-0.25,0.50))+
  scale_fill_manual(values=colrs) +
  facet_grid(program~pheno2) + #,scales="free_y"
  mytheme + theme(axis.text.x = element_blank(), axis.title.x=element_blank(),legend.title=element_blank()) + 
  ylab(bquote('Estimate ('*rho*')')) +
  labs(title="SumHer")
ldak_gencor_barplot3
ggsave(ldak_gencor_barplot3,file=paste(working_dir_ldak2,"summary_GWASgencor_sumher_GCTAandLDAK_adjTBV.pdf",sep=""),width=10,height=10)


ldak_gencor_barplot2<-ggplot(data=subset(t_gwas,model=="ldak"&GC=="yes"&adj==""),color="black",aes(x=region,y=rg,width=0.7,fill=measure)) +  #fill=program,,alpha=sig
  geom_bar(aes(y=rg),position=position_dodge(), stat="identity") +
  # geom_errorbar(data=t2,aes(ymin=CI95lower, ymax=CI95upper),
  # width=.2, position=position_dodge(.9)) + 
  # scale_alpha_discrete( range=c(0.3,1), na.value = 0) +
  # geom_errorbar(data=subset(t_gwas,model=="ldak"),aes(ymin=rg-se, ymax=rg+se,width=.5), position=position_dodge(.7)) +
  geom_errorbar(data=subset(t_gwas,model=="ldak"&GC=="yes"&adj==""),aes(ymin=rg-(1.96*se), ymax=rg+(1.96*se),width=.5), position=position_dodge(.7)) +
  coord_cartesian(ylim = c(-0.25,0.50))+
  scale_fill_manual(values=colrs) +
  facet_grid(.~pheno2) + #,scales="free_y"
  mytheme + theme(axis.text.x = element_blank(), axis.title.x=element_blank(),legend.title=element_blank()) + 
  ylab(bquote('Estimate ('*rho*')')) +
  labs(title="SumHer-ldak") 
ldak_gencor_barplot2
ggsave(ldak_gencor_barplot2,file=paste(working_dir_ldak2,"summary_GWASgencor_sumher_LDAK.pdf",sep=""),width=10,height=10)


#----------------------------------------------------------------------
# read summary result form LDSC run (from BGENIE output)  
#----------------------------------------------------------------------
## heritability
tmp<-read.table(file=paste(working_dir_ldsc2,"summary_h2_ldsc.table",sep=""),sep="\t",header=TRUE)
tmp<-subset(tmp,h2..se.!="")
tmp$h2<-sapply(strsplit(tmp$h2..se.," "),"[[",1)
tmp$h2.SE<-gsub("\\(|\\)","",sapply(strsplit(tmp$h2..se.," "),"[[",2))
tmp$Intercept<-sapply(strsplit(tmp$Intercept..se.," "),"[[",1)
tmp$Intercept.SE<-gsub("\\(|\\)","",sapply(strsplit(tmp$Intercept..se.," "),"[[",2))
tmp$Ratio<-sapply(strsplit(tmp$Ratio..se.," "),"[[",1)
tmp$Ratio.SE<-gsub("\\(|\\)","",sapply(strsplit(tmp$Ratio..se.," "),"[[",2))
tmp<-tmp[,-grep("\\.\\.",colnames(tmp))]
tmp$phenotype<-gsub("^handedness_|_h2.log","",tmp$File)
tmp$phenotype<-gsub("Weighted_mean","WeightedMean",tmp$phenotype)
tmp$region<-gsub("tract_","",tmp$phenotype)
tmp$type<- sapply(strsplit(tmp$phenotype,"_"),"[[",2)
tmp$pheno<-gsub("_"," ",gsub(":","\n",sapply(tmp$region,simpleCap)))
tmp$measure<-"NotLat"
tmp$measure[grep("AI",tmp$region)]<-"AI"
tmp$measure[grep("_left",tmp$region)]<-"L"
tmp$measure[grep("_right",tmp$region)]<-"R"
table(tmp$measure)
tmp$stat<-"h2"
# clean names
tmp$region<-gsub(paste(paste(unique(tmp$type),collapse="_|"),"_",sep=""),"",tmp$region)
tmp$region<-gsub(paste(paste(unique(tmp$type),collapse="_|"),"_",sep=""),"",gsub("_"," ",gsub("_left|_right|AI_","",tmp$region)))
tmp$region<-gsub("Weighted|Volume grey matter in ","",gsub("Mean in | on skeleton", "", tmp$region))
tmp$region<-gsub("_"," ",gsub(":","\n",sapply(tmp$region,simpleCap)))
tmp$pheno<-gsub("_"," ",gsub(":","\n",sapply(tmp$region,simpleCap)))
# 
tmp$program<-"LDSC"
tmp$sumstats<-"BGENIE"
# extra cols to match GCTA format
tmp$volume<-""
tmp$volume[grep("volume",tolower(tmp$File))]<-"GM"
# #tmp$type[tmp$volume=="GM"&]<-"Cortex" # need to fix this!
# tmp$type[grep("Cerebell|Pallid|Thalam|Hippocampus|Ventral Striatum|Putam|Caudate|Pallidum|Amygdala",tmp$region)]<-"Subcortical"
# tmp$type[grep("Cerebell|Pallid|Thalam|Hippocampus|Ventral Striatum|Putam|Caudate|Pallidum|Amygdala",tmp$region,invert=TRUE)]<-"Cortex"
# save
LDSCsummary_h2<-tmp; rm(tmp)


## genetic correlation
tmp<-read.table(paste(working_dir_ldsc2,"summary_p1_p2_GWASgencor_ldsc.table",sep=""))
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
tmp$program<-"LDSC"
tmp$sumstats<-"BGENIE"
LDSCsummary2_gencor<-tmp; rm(tmp)
LDSCsummary2_gencor$p1<-gsub(".sumstats.gz","",LDSCsummary2_gencor$p1)
LDSCsummary2_gencor$p2<-gsub(".sumstats.gz","",LDSCsummary2_gencor$p2)
LDSCsummary2_gencor$p2[LDSCsummary2_gencor$p2=="Sniekers2017"]<-paste("Sniekers2017","IQ",sep="_")
# for every phenotype, make first letter capital
LDSCsummary2_gencor$pheno1<-gsub("_"," ",gsub(":","\n",sapply(LDSCsummary2_gencor$p1,simpleCap)))
LDSCsummary2_gencor$pheno2<-gsub("_"," ",gsub(":","\n",sapply(LDSCsummary2_gencor$p2,simpleCap)))
LDSCsummary2_gencor$pheno2<-sapply(strsplit(LDSCsummary2_gencor$pheno2," "),"[[",2)
LDSCsummary2_gencor$pheno2_cat<-"Psychiatric"
w<-which(LDSCsummary2_gencor$pheno2=="EA"|LDSCsummary2_gencor$pheno2=="EA3"|
           LDSCsummary2_gencor$pheno2=="CP"|LDSCsummary2_gencor$pheno2=="IQ")
LDSCsummary2_gencor$pheno2_cat[w]<-"Cognitive"
table(LDSCsummary2_gencor$pheno2_cat)
#
LDSCsummary2_gencor$region<-gsub("AI VOLUME |Volume of grey matter in | left| right","",LDSCsummary2_gencor$pheno1)
LDSCsummary2_gencor$measure<-"NotLat"
LDSCsummary2_gencor$measure[grep("AI",LDSCsummary2_gencor$pheno1)]<-"AI"
LDSCsummary2_gencor$measure[grep("left",LDSCsummary2_gencor$pheno1)]<-"L"
LDSCsummary2_gencor$measure[grep("right",LDSCsummary2_gencor$pheno1)]<-"R"
# LDSCsummary2_gencor$measure[grep("adjTBV",LDSCsummary2_gencor$pheno1)]<-paste(LDSCsummary2_gencor$measure[grep("adjTBV",LDSCsummary2_gencor$pheno1)],"adjTBV",sep="")
LDSCsummary2_gencor$adj<-""
LDSCsummary2_gencor$adj[grep("adjTBV",LDSCsummary2_gencor$pheno1)]<-"adjTBV"
LDSCsummary2_gencor$rg<-as.numeric(LDSCsummary2_gencor$rg)
LDSCsummary2_gencor$se<-as.numeric(LDSCsummary2_gencor$se)
LDSCsummary2_gencor$p<-as.numeric(LDSCsummary2_gencor$p)
LDSCsummary2_gencor$sig<-""
LDSCsummary2_gencor$sig[LDSCsummary2_gencor$p<0.05]<-"*"
LDSCsummary2_gencor$sig<-factor(LDSCsummary2_gencor$sig,levels=c("","*"))

# save
write.csv(LDSCsummary2_gencor,file=paste(working_dir_ldsc2,"summary_GWASgencor_ldsc.csv",sep=""),row.names = FALSE,quote = FALSE)

# plot this, select relevant phenos
table(LDSCsummary2_gencor$p2)
t_gwas<-subset(LDSCsummary2_gencor,
               p2!="Okbay2016_EA" & p2!="Sniekers2017_IQ" & p2!="2017_OCD" & p2!="Lee2018_CP" &
                 p2!="Sklar2012_BIP" & p2!="Wray2018_MDD")
t_gwas<-t_gwas[with(t_gwas,order(pheno2_cat)),]
levs=c("IQ","EA3","ASD","ADHDeur","SCZ")
levs2=c("Intelligence","EA","ASD","ADHD","SCZ")
t_gwas$pheno2<-factor(t_gwas$pheno2,levels=levs)
levels(t_gwas$pheno2)
levels(t_gwas$pheno2)<-levs2
# t_gwas<-subset(LDSCsummary2_gencor,pheno2=="Grove ASD"|pheno2=="Ripke2014 SCZ")
ldsc_gencor_barplot<-ggplot(data=subset(t_gwas,adj==""),color="black",aes(x=region,y=rg,width=0.7,fill=measure)) +  #fill=program
  geom_bar(aes(y=rg),position=position_dodge(), stat="identity") +
  # 95% CI
  geom_errorbar(data=subset(t_gwas,adj==""),aes(ymin=rg-(1.96*se), ymax=rg+(1.96*se),width=.5), position=position_dodge(.7)) +
  coord_cartesian(ylim = c(-0.25,0.50))+
  # scale_alpha_discrete( range=c(1,0.5), na.value = 0) +
  scale_fill_manual(values=colrs) +
  facet_grid(~pheno2) + #,scales="free_y"
  mytheme + theme(axis.text.x = element_blank(), axis.title.x=element_blank(),legend.title=element_blank()) + 
  ylab(bquote('Estimate ('*rho*')')) +
  labs(title="LDSC") 
ldsc_gencor_barplot +
geom_hline(yintercept = 0.10,color="gray")

ldsc_gencor_barplot_adjTBV<-ggplot(data=t_gwas,color="black",aes(x=region,y=rg,width=0.7,fill=measure,alpha=adj)) +  #fill=program
  geom_bar(aes(y=rg),position=position_dodge(), stat="identity") +
  # 95% CI
  geom_errorbar(data=t_gwas,aes(ymin=rg-(1.96*se), ymax=rg+(1.96*se),width=.5), position=position_dodge(.7)) +
  coord_cartesian(ylim = c(-0.25,0.50))+
  scale_alpha_discrete( range=c(1,0.5), na.value = 0) +
  scale_fill_manual(values=colrs) +
  facet_grid(program~pheno2) + #,scales="free_y
  mytheme + theme(axis.text.x = element_blank(), axis.title.x=element_blank(),legend.title=element_blank()) + 
  ylab(bquote('Estimate ('*rho*')')) +
  labs(title="LDSC") 

ggsave(ldsc_gencor_barplot,file=paste(working_dir_ldsc2,"summary_GWASgencor_ldsc.pdf",sep=""),width=15,height=10)
ggsave(ldsc_gencor_barplot_adjTBV,file=paste(working_dir_ldsc2,"summary_adjTBV_GWASgencor_ldsc.pdf",sep=""),width=15,height=10)


ldsc_gencor_barplot<-ggplot(data=subset(t_gwas,adj==""),color="black",aes(x=region,y=rg,width=0.7,fill=measure)) +  #fill=program
  geom_bar(aes(y=rg),position=position_dodge(), stat="identity") +
  # 95% CI
  geom_errorbar(data=subset(t_gwas,adj==""),aes(ymin=rg-(1.96*se), ymax=rg+(1.96*se),width=.5), position=position_dodge(.7)) +
  coord_cartesian(ylim = c(-0.25,0.50))+
  # scale_alpha_discrete( range=c(1,0.5), na.value = 0) +
  scale_fill_manual(values=colrs) +
  facet_grid(program~pheno2) + #,scales="free_y"
  mytheme + theme(axis.text.x = element_blank(), axis.title.x=element_blank(),legend.title=element_blank()) + 
  ylab(bquote('Estimate ('*rho*')')) +
  labs(title="") 

#----------------------------------------------------------------------

#----------------------------------------------------------------------
# combine all gencor results, and save to out_dir
#----------------------------------------------------------------------
# extract legend
mylegend<-g_legend(ldak_gencor_barplot2+theme(legend.position="bottom"))

pdf(file=paste(out_dir,"summary_GWASgencor_ldsc_SUMER_gcta_ldak_comparison.pdf",sep=""),width=7,height=10)
# grid.arrange(ldsc_gencor_barplot+theme(legend.position="none"),ldak_gencor_barplot+theme(legend.position="bottom"),heights=c(1.5,2.5))
grid.arrange(ldsc_gencor_barplot+theme(legend.position="none"),
             ldak_gencor_barplot+theme(legend.position="none")+labs(title=""),
             mylegend,heights=c(1.5,2.5,0.2))
dev.off()


mylegend2<-g_legend(ldak_gencor_barplot3+theme(legend.position="bottom"))

pdf(file=paste(out_dir,"summary_GWASgencor_adjTBV_ldsc_SUMER_gcta_ldak_comparison.pdf",sep=""),width=7,height=10)
# grid.arrange(ldsc_gencor_barplot+theme(legend.position="none"),ldak_gencor_barplot+theme(legend.position="bottom"),heights=c(1.5,2.5))
grid.arrange(ldsc_gencor_barplot_adjTBV+theme(legend.position="none") +labs(title=""),
             ldak_gencor_barplot3+theme(legend.position="none")+labs(title=""),
             mylegend2,heights=c(1.5,2.5,0.2))
dev.off()



pdf(file=paste(out_dir,"summary_GWASgencor_ldsc_SUMER_ldak_comparison.pdf",sep=""),width=7,height=8)
grid.arrange(ldsc_gencor_barplot+theme(legend.position="none"),ldak_gencor_barplot2+theme(legend.position="none"),mylegend,heights=c(1,1,0.2))
dev.off()


# combine LDSC and LDAK results
cols<-c("program","pheno1","pheno2","rg","se","p")
summary_GWASgencor<-rbind(
  subset(LDSCsummary2_gencor, p2!="Okbay2016_EA" & p2!="Sniekers2017_IQ" & p2!="2017_OCD" & p2!="Sklar2012_BIP"  & p2!="Sklar2012_BIP" & p2!="Lee2018_CP" & p2!="Wray2018_MDD")[,cols],
  subset(LDAKsummary_gencor,(GC=="yes")&
           (p2!="Okbay2016_EA" & p2!="Sniekers2017_IQ" & p2!="2017_OCD" & p2!="Sklar2012_BIP"  & p2!="Sklar2012_BIP" & p2!="Lee2018_CP" & p2!="Wray2018_MDD")
           )[,cols]
)
# summary_GWASgencor$program<-gsub("-ldak|-gcta","",summary_GWASgencor$program)
summary_GWASgencor$pheno2<-gsub("3|eur","",summary_GWASgencor$pheno2)
summary_GWASgencor$pheno1<-gsub(" VOLUME Planum Temporale|Volume of grey matter in Planum Temporale ","",summary_GWASgencor$pheno1)
summary_GWASgencor$pheno1<-gsub("left","L",summary_GWASgencor$pheno1)
summary_GWASgencor$pheno1<-gsub("right","R",summary_GWASgencor$pheno1)

# wide format
summary_GWASgencor_wide<-summary_GWASgencor %>% 
  gather(variable,value,-(program:pheno2)) %>% 
  unite(temp, program, variable) %>% spread(temp, value)

write.csv(summary_GWASgencor_wide,
          file=paste(out_dir,"summary_GWASgencor_PT_LDSC_SUMHER.csv",sep=""),row.names = FALSE,quote=FALSE)
# round numbers
summary_GWASgencor$rg<-as.character(round(summary_GWASgencor$rg,digits=3))
summary_GWASgencor$se<-as.character(round(summary_GWASgencor$se,digits=2))
summary_GWASgencor$p2<-as.character(round(summary_GWASgencor$p,digits=2))
summary_GWASgencor$p2[summary_GWASgencor$p2=="0"]<-as.character(format(summary_GWASgencor$p,digits=2,scientific=TRUE))[summary_GWASgencor$p2=="0"]
summary_GWASgencor$p<-summary_GWASgencor$p2
summary_GWASgencor$rg_se<-paste(summary_GWASgencor$rg," (",summary_GWASgencor$se,")",sep="")
summary_GWASgencor_wide<-summary_GWASgencor %>% select(-rg,-se) %>% 
  gather(variable,value,-(program:pheno2)) %>% 
  unite(temp, program, variable) %>% spread(temp, value)
colnames(summary_GWASgencor_wide)
summary_GWASgencor_wide<-summary_GWASgencor_wide[,c("pheno1","pheno2","LDSC_rg_se","LDSC_p",
                                                    "SumHer-gcta_rg_se","SumHer-gcta_p",
                                                    "SumHer-ldak_rg_se","SumHer-ldak_p"
                                                    )]
sink(paste(out_dir,"summary_GWASgencor_PT_LDSC_SUMHER.tex"))
print(xtable(summary_GWASgencor_wide),include.rownames = FALSE)
sink()

write.csv(summary_GWASgencor_wide,
          file=paste(out_dir,"summary_GWASgencor_PT_LDSC_SUMHER_formated.csv",sep=""),row.names = FALSE,quote=FALSE)
#----------------------------------------------------------------------
# Compare results: GCTA vs LDSC vs LDAK(sumher) / sumstats from BGENIE
#----------------------------------------------------------------------
# H2
h2cols<-c("phenotype","region","measure","program","sumstats","h2","stat","type","volume") 
summary_h2<-merge(subset(GCTAsummary,sample=="Total"),LDSCsummary_h2,by.x=c(h2cols,"se","file","pheno"),by.y=c(h2cols,"h2.SE","File","pheno"),all=TRUE)
summary_h2<-merge(summary_h2,LDAKsummary_h2,by.x=c(h2cols,"se","file","pheno"),by.y=c(h2cols,"h2_SD","File","pheno"),all=TRUE)
summary_h2$adj<-""
summary_h2$adj[grep("adj|TBV",summary_h2$file)]<-"adjTBV"
summary_h2$region<-gsub(".adjTBV","",summary_h2$region)
summary_h2$program<-gsub(".adjTBV","",summary_h2$program)
summary_h2<-summary_h2[,c("region","measure","adj","program","sumstats","h2","se")] 
if (sum(duplicated(summary_h2))>0){summary_h2<-summary_h2[-which(duplicated(summary_h2)),]}
# wide format for program
summary_h2_wide <- summary_h2 %>%  
                  gather(variable,value,-(region:sumstats)) %>%
                  unite(temp, adj:sumstats, variable) %>% spread(temp, value)
colnames(summary_h2_wide)<-gsub("^_","",colnames(summary_h2_wide))
# subset(summary_h2,region=="Paracingulate Gyrus")[,c("region","DTI","volume","type","measure","program","h2","se")] %>%
# plot comparison
#------------------------------------------
t<-subset(summary_h2_wide,!is.na(`GCTA_-_h2`)&!is.na(`LDSC_BGENIE_h2`))
t[,colnames(t)[grep("h2|se",colnames(t))]]<-
  apply(t[,colnames(t)[grep("h2|se",colnames(t))]],2,as.numeric)
#(GCTA_h2,LDSC_h2,`LDSC-baselineLD_h2`)) %>%  gather(key=program,value=se,c(GCTA_se,LDSC_se,LDSC-baselineLD_se))
t2<- t %>%  gather(key=program1,value=h2,contains('_h2')) %>%  gather(key=program2,value=se,contains("_se"))
w<-which(gsub("_h2","",t2$program1)==gsub("_se","",t2$program2))
t2$adj<-""
t2$adj[grep("adjTBV",t2$program1)]<-"adjTBV"
t2$program1<-gsub("adjTBV_","",t2$program1)
t2$program<-sapply(strsplit(t2$program1,"_"),"[[",1)
t2$sumstats<-sapply(strsplit(t2$program1,"_"),"[[",2)
t2<-t2[w,c("region","measure","h2","se","program","sumstats","adj")];rm(w)
# compute 95CI: = (mean + (1.96 x SE)) to (mean - (1.96 x SE))
t2$CI95upper<-t2$h2+(1.96*t2$se)
t2$CI95lower<-t2$h2-(1.96*t2$se)
# to plot only subset that contains all rounds
t3<-subset(t,!is.na(LDSC_BGENIE_h2)) %>%  gather(key=program1,value=h2,contains('_h2')) %>%  gather(key=program2,value=se,contains("_se"))
w<-which(gsub("_h2","",t3$program1)==gsub("_se","",t3$program2))
t3$adj<-""
t3$adj[grep("adjTBV",t3$program1)]<-"adjTBV"
t3$program1<-gsub("adjTBV_","",t3$program1)
t3$program<-sapply(strsplit(t3$program1,"_"),"[[",1)
t3$sumstats<-sapply(strsplit(t3$program1,"_"),"[[",2)
t3<-t3[w,c("region","measure","h2","se","program","sumstats","adj")];rm(w)
# compute 95CI: = (mean + (1.96 x SE)) to (mean - (1.96 x SE))
t3$CI95upper<-t3$h2+(1.96*t3$se)
t3$CI95lower<-t3$h2-(1.96*t3$se)
# The palette with black:
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

t_comp_barplot2<-ggplot(data=t2,color="black",aes(x=measure,y=h2,width=0.7,fill=program,alpha=adj)) +  #fill=program
  geom_bar(aes(y=h2),position=position_dodge(), stat="identity") +
  # geom_errorbar(data=t2,aes(ymin=CI95lower, ymax=CI95upper),
                # width=.2, position=position_dodge(.9)) + 
  # geom_errorbar(data=t2,aes(ymin=h2-se, ymax=h2+se,width=.5), position=position_dodge(.7)) +
  geom_errorbar(data=t2,aes(ymin=h2-(1.96*se), ymax=h2+(1.96*se),width=.5), position=position_dodge(.7)) +
  scale_alpha_discrete(range=c(1,0.5)) +
  # facet_grid(adj~.) + #,scales="free_y"
  mytheme + theme(axis.text.x = element_text(angle = 0, hjust = 1, vjust=0.5), axis.title.x=element_blank()) +
  ylab(bquote('Estimate ('*h^2*')')) +
  # labs(title="Heritability estimates") + 
  mytheme2 + scale_fill_manual(values=cbbPalette ) + labs(title="")

t_comp_barplot1<-ggplot(data=subset(t2,adj==""),color="black",aes(x=measure,y=h2,width=0.7,fill=program)) +  #fill=program
  geom_bar(aes(y=h2),position=position_dodge(), stat="identity") +
  geom_errorbar(data=subset(t2,adj==""),aes(ymin=h2-(1.96*se), ymax=h2+(1.96*se),width=.5), position=position_dodge(.7)) +
  scale_alpha_discrete(range=c(0.5,1)) +
  # facet_grid(adj~.) + #,scales="free_y"
  # mytheme + 
  coord_cartesian(ylim = c(0, 0.6)) +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, vjust=0.5), axis.title.x=element_blank()) +
  ylab(bquote('Estimate ('*h^2*')')) +
  # labs(title="Heritability estimates") + 
  mytheme2 + scale_fill_manual(values=cbbPalette ) + labs(title="")

# t_comp_barplot1

t_comp_barplot_adj<-ggplot(data=subset(t2,adj=="adjTBV"),color="black",aes(x=measure,y=h2,width=0.7,fill=program,alpha=0.5)) +  #fill=program
  geom_bar(aes(y=h2),position=position_dodge(), stat="identity") +
  geom_errorbar(data=subset(t2,adj=="adjTBV"),aes(ymin=h2-(1.96*se), ymax=h2+(1.96*se),width=.5), position=position_dodge(.7)) +
  # scale_alpha_discrete(range=c(1,0.5)) +
  # facet_grid(adj~.) + #,scales="free_y"
  # mytheme + 
  coord_cartesian(ylim = c(0, 0.6)) +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, vjust=0.5), axis.title.x=element_blank()) +
  ylab(bquote('Estimate ('*h^2*')')) +
  # labs(title="Heritability estimates") + 
  mytheme2 + scale_fill_manual(values=cbbPalette ) + labs(title="adjTBV")

legend<- get_legend(t_comp_barplot1 + theme(legend.position="bottom"))
combined<-plot_grid(t_comp_barplot1 +  theme(legend.position="none"),
          t_comp_barplot_adj +  theme(legend.position="none"),
          labels=c("A","B"))
t_comp_barplot_combined_l <- plot_grid( combined, legend, ncol = 1, rel_heights = c(1.5, .1))

# To use for fills, add
ggsave(t_comp_barplot1,file=paste(out_dir,"h2_estimates_comparison_barplot.pdf",sep=""),width=10,height=5)
ggsave(t_comp_barplot_adjTBV,file=paste(out_dir,"h2_estimates_comparison_barplot_adjTBV.pdf",sep=""),width=10,height=5)
ggsave(t_comp_barplot2,file=paste(out_dir,"h2_estimates_comparison_barplot_adjTBV_1plot.pdf",sep=""),width=12,height=7)
ggsave(t_comp_barplot_combined_l ,file=paste(out_dir,"h2_estimates_comparison_barplot_adjTBV_AB.pdf",sep=""),width=12,height=7)

# clean
rm(t,t_comp_barplot1,t_comp_barplot2,t_comp_barplot_adjTBV,t_comp_barplot_combined_l)


write.csv(t3,paste(out_dir,"/estimates4power.csv",sep=""),row.names=FALSE)
# t4<-subset(t3,program!="SUMHER-gcta")
# t5<-data.frame(do.call("rbind",tapply(t4$h2,t4$program,summary)))
# rm(t4,t5)

