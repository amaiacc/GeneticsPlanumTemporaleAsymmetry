# read data from PT analysis, and rs41298373 genotypes, and create input files for
#---------------------------------------------------------------------
# Plots
#---------------------------------------------------------------------
library(GGally);library(ggplot2)
library(gridExtra);library(grid)
library(plyr)
source("https://gist.githubusercontent.com/benmarwick/2a1bb0133ff568cbe28d/raw/fb53bd97121f7f9ce947837ef1a4c65a73bffb3f/geom_flat_violin.R")
cols=c("#d95f02","#1f78b4","#33a02c","#000000") # '#1b9e77','#d95f02','#7570b3 # "#a6cee3"  ,"#7570b3"
cols2=c('#1b9e77','#d95f02','#7570b3' )
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
                             title=element_text(size=24),
                             axis.title=element_text(size=24),axis.text=element_text(size=24,colour="black"),
                             axis.title.x=element_blank(),
                             legend.text=element_text(size=24), legend.title=element_blank(), legend.position="bottom") 
# function to make first letter of a string capital
simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  s2<- paste(toupper(substring(s, 1,1)), substring(s, 2),  sep="", collapse=" ")
  return(s2)
}
lb <- function(x) {mean(x,na.rm=TRUE) - sd(x,na.rm=TRUE)}
ub <- function(x) {mean(x,na.rm=TRUE) + sd(x,na.rm=TRUE)}
#---------------------------------------------------------------------
pheno_root="ukb25465_ukb25468"
subset_name="imagingT1_N18057"
# Define directories, dependent on  system
if (Sys.info()['sysname']=='Windows') {dir="P://workspaces/"} else {dir="/data/workspaces/lag/workspaces/"}
pheno_dir=paste(dir,"/lg-ukbiobank/working_data/amaia/pheno_files/genetic_v2/",pheno_root,"/summary_phenotypes/",sep="")
working_dir=paste(dir,"/lg-ukbiobank/working_data/amaia/genetic_data/PT/",subset_name,"/lead_snps/genotypes/",sep="")
setwd(working_dir)

# read pheno data, which also include covariates
phenos<-read.table(paste(pheno_dir,"/",pheno_root,"_imaging_noBioCovs_noAssessmentCVolumePhenotypes_residuals_wHeader.table",sep=""),header=TRUE)
phenos<-phenos[,c(1:52,grep("Planum_Temporale",colnames(phenos)))] # all covariates, note that handedness is all NA, because it was converted to factor and then info lost when trying to recode...
phenos$sex<-as.factor(phenos$sex)
# remove outliers if |value|>mean+4*sd
ps<-colnames(phenos)[grep("Planum_Temporale",colnames(phenos))]
ps<-ps[grep("residuals",ps,invert=TRUE)]
pheno_names<-ps[grep("AI_|left$|right$",ps)]
t=4
pt2<-data.frame(apply(phenos[,pheno_names],2,function(x){
  m<-mean(x,na.rm=T)
  sd<-sd(x,na.rm=T)
  min<-m-t*sd
  max<-m+t*sd
  x[which(x<min|x>max)]<-NA
  return(x)
})
)
phenos[,pheno_names]<-pt2
rm(pt2,ps,t)


snps<-gsub(".ped|_","",gsub(subset_name,"",gsub("ukb_imp_v3_","",list.files(working_dir,pattern="ped"))))
# run sex interaction models per snp
for (snp in snps){
  # read genotype data for snp
  ped<-read.table(paste(working_dir,"/ukb_imp_v3_",subset_name,"_",snp,".ped",sep=""),stringsAsFactors = FALSE)
  colnames(ped)<-c("ID1","ID2","mID","pID","sex","pheno","A1","A2")
  ped$ID1<-ped$ID2
  ped$geno<-paste(ped$A1,ped$A2,sep="")
  ped$snp<-snp
  alleles<-names(sort(table(c(ped$A1,ped$A2)))) # minor, major
  if (length(which(alleles=="0"))>0){
    alleles<-alleles[-(which(alleles=="0"))]
  }
  # recode for the additive model
  ped$genoAdd<-NA
  ped$genoAdd[ped$geno=="00"]<-NA
  ped$genoAdd[ped$geno==paste(rep(alleles[2],2),collapse="")]<-0 # two major alleles
  ped$genoAdd[ped$geno==paste(alleles[2],alleles[1],collapse="",sep="")]<-1 # heterozygote
  ped$genoAdd[ped$geno==paste(alleles[1],alleles[2],collapse="",sep="")]<-1 # heterozygote
  ped$genoAdd[ped$geno==paste(rep(alleles[1],2),collapse="")]<-2 # two minor alleles
  # check
  table(ped$geno)
  table(ped$genoAdd)
  # merge with phenos
  phenos_ped<-merge(ped[,c("ID1","ID2","genoAdd","geno","snp")],phenos,stringsAsFactors=FALSE)
  # general lm
  effect<-lm(AI_VOLUME_Planum_Temporale~
               age+zage2+sex+
               PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+
               array+
               Scanner_lateral_X_brain_position+Scanner_transverse_Y_brain_position+Scanner_longitudinal_Z_brain_position +
               genoAdd, 
             data=phenos_ped)
  
  interact<-lm(AI_VOLUME_Planum_Temporale~
                 genoAdd*sex+
                 age+zage2+sex+
                 PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+
                 array+
                 Scanner_lateral_X_brain_position+Scanner_transverse_Y_brain_position+Scanner_longitudinal_Z_brain_position, 
               data=phenos_ped)
  # save
  sink(file=paste(working_dir,"/LM_AI_PT_",snp,"_sex_interaction.txt",sep=""))
  print(summary(effect))
  print(anova(effect,interact))
  print(summary(interact))
  sink()
  
    # clean
  rm(ped,alleles,phenos_ped,effect,interact)
  
}

