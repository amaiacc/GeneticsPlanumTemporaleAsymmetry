# read data from PT analysis, and rs41298373 genotypes, and create input files for


# Define directories, dependent on  system
if (Sys.info()['sysname']=='Windows') {dir="P://workspaces/"} else {dir="/data/workspaces/lag/workspaces/"}
pheno_dir=paste(dir,"/lg-ukbiobank/working_data/amaia/pheno_files/genetic_v2/ukb21288_ukb21293/summary_phenotypes/",sep="")
working_dir=paste(dir,"/lg-ukbiobank/working_data/amaia/genetic_data/PT/rs41298373/",sep="")
setwd(working_dir)

# read pheno data, which also include covariates
phenos<-read.table(paste(pheno_dir,"ukb21288_ukb21293_imaging_noBioCovs_noAssessmentCVolumePhenotypes_residuals_wHeader.table",sep="/"),header=TRUE)
phenos<-phenos[,c(1:43,grep("Planum_Temporale",colnames(phenos)))] # all covariates, note that handedness is all NA, because it was converted to factor and then info lost when trying to recode...

# list of phenotype IDs for which VBM is available
vbm_ids<-read.table("vbm_sids.txt"); colnames(vbm_ids)<-"ID1"
vbm_ids$VBM<-1
# read genotype data for rs41298373
ped<-read.table("ukb_cal_snpQC_imagingT1_N12245_rs41298373.ped")
colnames(ped)<-c("ID1","ID2","mID","pID","sex","pheno","A1","A2")
ped$rs41298373<-paste(ped$A1,ped$A2,sep="")
# recode for the additive model
ped$rs41298373add<-NA
ped$rs41298373add[ped$rs41298373=="00"]<-NA
ped$rs41298373add[ped$rs41298373=="GG"]<-0
ped$rs41298373add[ped$rs41298373=="AG"]<-1
ped$rs41298373add[ped$rs41298373=="AA"]<-2

phenos_ped<-merge(ped,phenos,stringsAsFactors=FALSE)


# general lm
effect<-lm(AI_VOLUME_Planum_Temporale~
             age+zage2+sex+
             PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+
             array+
             Scanner_lateral_X_brain_position+Scanner_transverse_Y_brain_position+Scanner_longitudinal_Z_brain_position +
             rs41298373add, 
           data=phenos_ped)

interact<-lm(AI_VOLUME_Planum_Temporale~
              rs41298373add*sex+
             age+zage2+sex+
             PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+
             array+
             Scanner_lateral_X_brain_position+Scanner_transverse_Y_brain_position+Scanner_longitudinal_Z_brain_position, 
           data=phenos_ped)
anova(effect,interact)

# check whether any of the SNP is associated with any of the scanner parameters 
x<-lm(phenos_ped$Scanner_lateral_X_brain_position~phenos_ped$rs41298373add)
y<-lm(phenos_ped$Scanner_transverse_Y_brain_position~phenos_ped$rs41298373add)
z<-lm(phenos_ped$Scanner_longitudinal_Z_brain_position~phenos_ped$rs41298373add)
tbv<-lm(phenos_ped$Volume_of_brain_grey_white~phenos_ped$rs41298373add)






# define covs
covs<-c("array","age","zage2",paste("PC",1:10,sep=""),
        "Scanner_lateral_X_brain_position","Scanner_transverse_Y_brain_position","Scanner_longitudinal_Z_brain_position")

# combine ped file with the phenotype file and VBM IDs
phenos2<-merge(phenos,vbm_ids,by="ID1",all=TRUE)
phenos_ped<-merge(ped,phenos2,by=c("ID1","ID2","sex")); rm(phenos2)
sum(is.na(phenos_ped$rs41298373add))
sum(is.na(phenos_ped$AI_VOLUME_Planum_Temporale))
sum(is.na(phenos_ped$VBM))
phenos_ped<-subset(phenos_ped,!(is.na(rs41298373add))) # remove if no genotype data avail.
table(is.na(phenos_ped$rs41298373add))
# check how many NA's there are per column
apply(phenos_ped,2,function(x) sum(is.na(x)))
# double check covariates
apply(phenos_ped[,covs],2,function(x) sum(is.na(x)))
phenos_ped<-subset(phenos_ped,!is.na(phenos_ped$Scanner_lateral_X_brain_position))

# subset for vbm
phenos_ped4voxelwise<-subset(phenos_ped,!is.na(VBM)) # only samples that did not fail VBM
# select cols to save
phenos_ped4voxelwise<-phenos_ped4voxelwise[,c("ID1","sex","rs41298373add","rs41298373",covs)]
colnames(phenos_ped4voxelwise)<-c("SID","sex","rs41298373add","rs41298373",covs)
write.csv(phenos_ped4voxelwise,file="rs41298373_covariates_4voxelwise.csv",row.names = FALSE,quote=FALSE)


# select cols to save for vertex-wise analyses
phenos_ped4vertexwise<-phenos_ped[,c("ID1","sex","rs41298373add","rs41298373",covs)]
colnames(phenos_ped4vertexwise)<-c("SID","sex","rs41298373add","rs41298373",covs)
write.csv(phenos_ped4vertexwise,file="rs41298373_covariates_4vertexwise.csv",row.names = FALSE,quote=FALSE)


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
# check correlation between AI and lateral volumes
p<-"Planum_Temporale"
# select columns, and clean colnames
pt<-phenos[,c(1,grep("Planum_Temporale",colnames(phenos)))]
colnames(pt)<-gsub("_VOLUME_Planum_Temporale|Volume_of_grey_matter_in_Planum_Temporale_","",colnames(pt))
pt<-pt[,colnames(pt)[grep("ID|left|right|AI$",colnames(pt))]]

# plot, all values, including outliers
ggpairs(pt[, -1], 
        upper = list(continuous = wrap("cor", size = 10)), 
        lower = list(continuous = "smooth")) + theme_bw()
# remove outliers if |value|>mean+4*sd
t=4
pt2<-data.frame(apply(pt[,-1],2,function(x){
  m<-mean(x,na.rm=T)
  sd<-sd(x,na.rm=T)
  min<-m-t*sd
  max<-m+t*sd
  x[which(x<min|x>max)]<-NA
  return(x)
})
)
pt2<-cbind(ID1=pt[,1],pt2)
p2<-gsub("_"," ",p)
## These parameters are applied to each plot we create
pdf(file=paste("../correlations_",p,"_phenotypes.pdf",sep=""),width=5,height=5)
ggpairs(pt2[,-1], 
        upper = list(continuous = wrap("cor", size = 8)), 
        lower = list(continuous = "smooth") #,line_color="red",alpha = 0.7
        ) + 
  theme_bw() + ggtitle(p2)
dev.off()

# check distribution of PT phenotypes, dependent on rs41298373 genotype
snp="rs41298373"
phenos_ped2<-merge(ped,pt2,by=c("ID1"))
pheno_names<-colnames(phenos_ped2)[grep("AI$|left$|right$",colnames(phenos_ped2))]
for (pheno in pheno_names){
  pheno2=simpleCap(gsub("^_|_$","",gsub(p,"",gsub("_VOLUME_|Volume_of_grey_matter_in_","",pheno))))
  print(tapply(phenos_ped2[,pheno],phenos_ped2[,snp],function(x) mean(x,na.rm=TRUE)))
  t<-phenos_ped2[,c(snp,pheno)];colnames(t)<-c("snp","pheno")
  t<-subset(t,!is.na(pheno)&snp!="00")
  sumld =  ddply(t, ~snp, summarise, mean = mean(pheno,na.rm=TRUE), 
                 median = median(pheno,na.rm=TRUE), 
                 lower = lb(pheno), upper = ub(pheno))
  
  # plot
  dist_snp<-ggplot(data=t,aes(y=pheno,x=snp,fill=snp)) +
    geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8) +
    geom_point(aes(y =pheno, color = snp)) + #, position = position_jitter(width = .15), size = .5, alpha = 0.8
    geom_boxplot(width = .1, guides = FALSE, outlier.shape = NA, alpha = 0.5) +
    geom_point(data = sumld, aes_string(x = "snp", y = "mean"), position = position_nudge(x = 0.3), size = 2.5) +
    geom_errorbar(data = sumld, aes_string(x = "snp",ymin = "lower", ymax = "upper", y = "mean"), position = position_nudge(x = 0.3), width = 0) +
    guides(fill = FALSE) + guides(color = FALSE) +
    scale_color_manual(values=cols2) +
    scale_fill_manual(values=cols2) +
    coord_flip() + ggtitle(pheno2) +
    mytheme2 +
    theme(axis.title.x = element_blank(),axis.title.y = element_blank())
    
  pname<-paste("distr_snp",p,pheno2,sep="_")
  assign(pname,dist_snp)
  rm(pheno2,dist_snp,pname,sumld,t)
}

# combine into one plot
list_p<-lapply(ls(pattern=paste("distr_snp*.*",p,sep="")),get)
g<-grid.arrange(arrangeGrob(grobs=list_p), top=textGrob(paste("Effect of ", snp, " on\n",p2, " grey matter volume phenotypes",sep=""),gp=gpar(fontsize=24,font=2)))
ggsave(file=paste(snp,"_effects_on_",p,".pdf",sep=""),g,width=10,height=10)

