# read data from PT analysis, and rs41298373 genotypes, and create input files for
#---------------------------------------------------------------------
# Plots
#---------------------------------------------------------------------
library(GGally);library(ggplot2)
library(gridExtra);library(grid)
library(plyr)
library(tidyr)
library(lme4); library(lmerTest)
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
# Define directories, dependent on  system
if (Sys.info()['sysname']=='Windows') {
  # dir="P://workspaces/"
  dir="\\\\data/lag/workspaces/"
} else {dir="/data/workspaces/lag/workspaces/"}
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

phenos$PT_LplusR<-phenos$Volume_of_grey_matter_in_Planum_Temporale_left+phenos$Volume_of_grey_matter_in_Planum_Temporale_right

# convert phenos to long format:
LR_cols<-colnames(phenos)[grep("left|right",colnames(phenos))]
LR_cols<-LR_cols[grep("residuals",LR_cols,invert=TRUE)]
phenos_wide <- phenos %>% 
                    select_at(vars(-contains("residuals"))) %>% 
                    gather_("side", "measure", LR_cols)


#
snps<-gsub(".ped|_","",gsub(subset_name,"",gsub("ukb_imp_v3_","",list.files(working_dir,pattern="ped"))))
# run side interaction models per snp
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
  phenos_ped<-merge(ped[,c("ID1","ID2","genoAdd","geno","snp")],phenos_wide,stringsAsFactors=FALSE)
  phenos_ped<-subset(phenos_ped,!is.na(genoAdd))
  # general lm
  effect1<-lmer(measure~
               (1|ID1) +
               age+zage2+sex+
               PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+
               array+
               Scanner_lateral_X_brain_position+Scanner_transverse_Y_brain_position+Scanner_longitudinal_Z_brain_position +
               side + genoAdd, 
             data=phenos_ped)
  
  interact1<-lmer(measure~
                   (1|ID1) +
                   age+zage2+sex+
                   PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+
                   array+
                   Scanner_lateral_X_brain_position+Scanner_transverse_Y_brain_position+Scanner_longitudinal_Z_brain_position +
                   genoAdd*side,
                 data=phenos_ped)
  
  anova(effect1,interact1)
  
  #
  effect2<-lmer(measure~
                 (1|ID1) +
                 age+zage2+sex+
                 PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+
                 array+
                 Scanner_lateral_X_brain_position+Scanner_transverse_Y_brain_position+Scanner_longitudinal_Z_brain_position +
                 totalBV +
                 side + genoAdd, 
               data=phenos_ped)
  
  interact2<-lmer(measure~
                   (1|ID1) +
                   age+zage2+sex+
                   PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+
                   array+
                   Scanner_lateral_X_brain_position+Scanner_transverse_Y_brain_position+Scanner_longitudinal_Z_brain_position +
                  totalBV + 
                  genoAdd*side,
                 data=phenos_ped)
  anova(effect2,interact2)
  
  #
  effect3<-lmer(measure~
                  (1|ID1) +
                  age+zage2+sex+
                  PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+
                  array+
                  Scanner_lateral_X_brain_position+Scanner_transverse_Y_brain_position+Scanner_longitudinal_Z_brain_position +
                  PT_LplusR +
                  side + genoAdd, 
                data=phenos_ped)
  
  interact3<-lmer(measure~
                    (1|ID1) +
                    age+zage2+sex+
                    PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+
                    array+
                    Scanner_lateral_X_brain_position+Scanner_transverse_Y_brain_position+Scanner_longitudinal_Z_brain_position +
                    PT_LplusR + 
                    genoAdd*side,
                  data=phenos_ped)
  anova(effect3,interact3)
  
  #
  
  # save
  sink(file=paste(working_dir,"/LMER_PT_",snp,"_side_interaction.txt",sep=""))
  
  cat("\n------------------------------------\nSNP: ",snp,"\n------------------------------------\n")
  
  
    cat("\n------------------------------------\nModel 1\n------------------------------------\n")
    cat("\n-----------Base model-----------\n")
    print(summary(effect1))
    cat("\n-----------Anova-----------\n")
    print(anova(effect1,interact1))
    cat("\n-----------Test model-----------\n")
    print(summary(interact1))
    cat("\n------------------------------------\nModel 2: adjusted for TBV\n------------------------------------\n")
    cat("\n-----------Base model-----------\n")
    print(summary(effect2))
    cat("\n-----------Anova-----------\n")
    print(anova(effect2,interact2))
    cat("\n-----------Test model-----------\n")
    print(summary(interact2))
    cat("\n------------------------------------\nModel 3: adjusted for (L+R)\n------------------------------------\n")
    cat("\n-----------Base model-----------\n")
    print(summary(effect3))
    cat("\n-----------Anova-----------\n")
    print(anova(effect3,interact3))
    cat("\n-----------Test model-----------\n")
    print(summary(interact3))
    
  sink()
  
  # clean
  rm(ped,alleles,phenos_ped)
  rm(effect1,interact1,effect2,interact2,effect3,interact3)
  
}

