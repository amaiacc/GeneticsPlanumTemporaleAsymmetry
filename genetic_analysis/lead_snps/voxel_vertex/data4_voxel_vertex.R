# read data from PT analysis, and genotypes for lead SNPs
# and create input files for VBM and vertex-wise analyses

options(stringsAsFactors = FALSE)
#
pheno_root="ukb25465_ukb25468"
subset_name="imagingT1_N18057"
# Define directories, dependent on  system
if (Sys.info()['sysname']=='Windows') {dir="P://workspaces/"} else {dir="/data/workspaces/lag/workspaces/"}
pheno_dir=paste(dir,"/lg-ukbiobank/working_data/amaia/pheno_files/genetic_v2/",pheno_root,"/summary_phenotypes/",sep="")
geno_dir=paste(dir,"/lg-ukbiobank/working_data/amaia/genetic_data/PT/",subset_name,"/lead_snps/genotypes/",sep="")
working_dir=paste(dir,"/lg-ukbiobank/working_data/amaia/genetic_data/PT/",subset_name,"/lead_snps/voxel_vertex/",sep="")
vbm_dir=paste(dir,"/lg-ukbiobank/analysis/xiangzhen/vol_asym/analysis/",sep="")
setwd(working_dir)

# read pheno data, which also include covariates
phenos<-read.table(paste(pheno_dir,pheno_root,"_imaging_noBioCovs_noAssessmentCVolumePhenotypes_residuals_wHeader.table",sep=""),header=TRUE)

n<-grep("Volumetric_scaling_from_T1_head_image_to_standard_space",colnames(phenos))[1]
phenos<-phenos[,c(1:n,grep("Planum_Temporale",colnames(phenos)))] # all covariates
rm(n)

# define covariates as used for the LM pre-GWAS:
scanner_position<-colnames(phenos)[grep("Scanner",colnames(phenos))]
covs<-c("age","zage2","sex",paste("PC",1:10,sep=""),"array",scanner_position ) 

# list of phenotype IDs for which VBM is available
vbm_ids<-read.table("vbm_sids.txt"); colnames(vbm_ids)<-"ID1"
vbm_ids$VBM<-1

# combine phenos and VBM flag
phenos<-merge(phenos,vbm_ids,all.x=TRUE,by="ID1")

# read genotype data for lead SNPs

for (p in list.files(geno_dir,pattern="ped")) {
  ped<-read.table(paste(geno_dir,p,sep="/"))
  colnames(ped)<-c("ID1","ID2","mID","pID","sex","pheno","A1","A2")  
  # get snp
  snp<-sapply(strsplit(gsub(".ped","",p),"_"),"[[",6)
  # create genotype column
  ped[,snp]<-paste(ped$A1,ped$A2,sep="")
  # get possible genotypes
  genos<-names(sort(table(ped[,snp])))
  w<-which(genos=="00")
  if (length(w)>0){
    genos<-genos[-w]    
  }
  rm(w)
  # code genotypes as additive
  ped[,paste(snp,"add",sep="")]<-NA
  ## homozygous minor
  ped[ped[,snp]==genos[1],paste(snp,"add",sep="")]<-2
  ## heterozygous
  ped[ped[,snp]==genos[2],paste(snp,"add",sep="")]<-1
  ## homozyous major
  ped[ped[,snp]==genos[3],paste(snp,"add",sep="")]<-0
  
  # remove alleles, otherwise merging won't work
  ped<-ped[,c("ID1","ID2",snp,paste(snp,"add",sep=""))]
  
  # combine
  if (p==list.files(geno_dir,pattern="ped")[1]){ ped_all<-ped
  } else { ped_all<-merge(ped_all,ped,stringAsFactors=FALSE,all=TRUE)}
  # clean
  rm(ped,snp,genos)
}

# make ID1==ID2
ped_all$ID1<-ped_all$ID2

# combine phenos and genotypes
phenos_ped<-merge(ped_all,phenos,stringsAsFactors=FALSE)

# define snps  of interest for VBM and vertex-wise analyses
snps<-c("rs41298373","rs7420166")
cols_snps<-colnames(phenos_ped)[grep(paste(snps,collapse="|"),colnames(phenos_ped))]


# subset for vbm:
## only samples that contain all covariates
phenos_ped4voxelwise<-subset(phenos_ped,!is.na(Scanner_lateral_X_brain_position)) # 5 removed
# select cols to save
phenos_ped4voxelwise<-phenos_ped4voxelwise[,c("ID1","sex",cols_snps,covs)]
colnames(phenos_ped4voxelwise)<-c("SID","sex",cols_snps,covs)
# # only samples that did not fail VBM
# phenos_ped4voxelwise<-subset(phenos_ped,!is.na(VBM)) 

# for each SNP, only samples that have genotype data
for (s in snps){
  cols<-c("SID","sex",cols_snps[grep(s,cols_snps)],covs)
  tmp<-phenos_ped4voxelwise[which(phenos_ped4voxelwise[,s]!="00"),cols]
  write.csv(tmp,file=paste(s,"_covariates_4voxelwise.csv",sep=""),row.names = FALSE,quote=FALSE)
  rm(tmp,cols)
}



# select cols to save for vertex-wise analyses
phenos_ped4vertexwise<-phenos_ped[,c("ID1","sex",cols_snps,covs)]
colnames(phenos_ped4vertexwise)<-c("SID","sex",cols_snps,covs)
write.csv(phenos_ped4vertexwise,file="leadSNPs_covariates_4vertexwise.csv",row.names = FALSE,quote=FALSE)


#---------------------------------------------------------------------
# Plots and LM
# Further exploration of SNP effects
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

# check how many NA's there are per column
apply(phenos_ped[,covs],2,function(x) sum(is.na(x)))
phenos_ped2<-subset(phenos_ped,!is.na(phenos_ped$Scanner_lateral_X_brain_position))

# check correlation between AI and lateral volumes
p<-"Planum_Temporale"
# select columns, and clean colnames
pt<-phenos[,c(1,grep("Planum_Temporale",colnames(phenos)))]
colnames(pt)<-gsub("_VOLUME_Planum_Temporale|Volume_of_grey_matter_in_Planum_Temporale_","",colnames(pt))
pt<-pt[,colnames(pt)[grep("ID|left|right|AI$",colnames(pt))]]

# plot, all values, including outliers
ggpairs(pt[,-c(1,grep("residuals",colnames(pt)))], 
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
pdf(file=paste("../../correlations_",p,"_phenotypes.pdf",sep=""),width=5,height=5)
ggpairs(pt2[,-c(1,grep("residuals",colnames(pt2)))], 
        upper = list(continuous = wrap("cor", size = 8)), 
        lower = list(continuous = "smooth") #,line_color="red",alpha = 0.7
) + 
  theme_bw() + ggtitle(p2)
dev.off()

# check distribution of PT phenotypes, dependent on genotypes
pheno_names<-colnames(phenos_ped2)[grep("AI_|left$|right$",colnames(phenos_ped2))]
pheno_names<-pheno_names[grep(p,pheno_names)]

## per SNP
snps2<-colnames(phenos_ped2)[grep("^rs",colnames(phenos_ped2))]
snps2<-snps2[grep("add$",snps2,invert=TRUE)]
for (s in snps2){
  for (pheno in pheno_names){
    pheno2=simpleCap(gsub("_"," ",gsub("__","_",gsub("^_|_$","",gsub(p,"",gsub("_VOLUME_|Volume_of_grey_matter_in_","",pheno))))))
    print(tapply(phenos_ped2[,pheno],phenos_ped2[,s],function(x) mean(x,na.rm=TRUE)))
    t<-phenos_ped2[,c(s,pheno)];colnames(t)<-c("snp","pheno")
    # remove samples if no pheno or genotype data avail.
    t<-subset(t,!is.na(pheno)&snp!="00")
    sumld =  ddply(t, ~snp, summarise, mean = mean(pheno,na.rm=TRUE), 
                   median = median(pheno,na.rm=TRUE), 
                   lower = lb(pheno), upper = ub(pheno))
    
    # plot
    dist_snp<-ggplot(data=t,aes(y=pheno,x=snp,fill=snp)) +
      geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .5) +
      geom_point(aes(y =pheno, color = snp)) + #, position = position_jitter(width = .15), size = .5, alpha = 0.8
      geom_boxplot(width = .1, outlier.shape = NA, alpha = 0.2) +
      geom_point(data = sumld, aes_string(x = "snp", y = "mean"), position = position_nudge(x = 0.3), size = 2.5) +
      geom_errorbar(data = sumld, aes_string(x = "snp",ymin = "lower", ymax = "upper", y = "mean"), position = position_nudge(x = 0.3), width = 0) +
      guides(fill = FALSE) + guides(color = FALSE) +
      scale_color_manual(values=cols2) +
      scale_fill_manual(values=cols2) +
      coord_flip() + ggtitle(pheno2) +
      mytheme2 +
      theme(axis.title.x = element_blank(),axis.title.y = element_blank())
    
    pname<-gsub(" ","_",paste("distr",s,p,pheno2,sep="_"))
    assign(pname,dist_snp)
    rm(pheno2,dist_snp,pname,sumld,t)
  }
  rm(pheno)
  
  # combine into one plot
  l<-ls(pattern=paste("distr_",s,"*.*",p,sep=""))
  list_p<-lapply( l[grep("Residuals",l,invert=TRUE)],get)
  g<-grid.arrange(arrangeGrob(grobs=list_p), top=textGrob(paste("Effect of ", s, " on\n",p2, " grey matter volume phenotypes",sep=""),gp=gpar(fontsize=24,font=2)))
  ggsave(file=paste("../",s,"_effects_on_",p,".pdf",sep=""),g,width=10,height=10)
  # residuals
  list_pr<-lapply(l[grep("Residuals",l)],get)
  gr<-grid.arrange(arrangeGrob(grobs=list_pr), top=textGrob(paste("Effect of ", s, " on\nresidualized ",p2, "\ngrey matter volume phenotypes",sep=""),gp=gpar(fontsize=24,font=2)))
  ggsave(file=paste("../",s,"_effects_on_residuals_",p,".pdf",sep=""),gr,width=10,height=10)
  # clean
  rm(l,list_p,g,list_pr,gr)
  rm(list=ls(pattern=paste("distr_",s,"*.*",p,sep="")))

}

# check effect of covariates/interactions
model_formula=

for (s in snps2){
  
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
  
  interact<-lm(AI_VOLUME_Planum_Temporale~
                 rs41298373add*sex+
                 age+zage2+sex+
                 PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+
                 array+
                 Scanner_lateral_X_brain_position+Scanner_transverse_Y_brain_position+Scanner_longitudinal_Z_brain_position, 
               data=phenos_ped)
  
  
  # check whether any of the SNP is associated with any of the scanner parameters 
  x<-lm(phenos_ped$Scanner_lateral_X_brain_position~phenos_ped$rs41298373add)
  y<-lm(phenos_ped$Scanner_transverse_Y_brain_position~phenos_ped$rs41298373add)
  z<-lm(phenos_ped$Scanner_longitudinal_Z_brain_position~phenos_ped$rs41298373add)
  tbv<-lm(phenos_ped$Volume_of_brain_grey_white~phenos_ped$rs41298373add)
  h<-lm(phenos_ped$handedness~phenos_ped$rs41298373add)
  
}