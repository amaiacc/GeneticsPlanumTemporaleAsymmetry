#' ---
#' title: ukb25465_ukb25468 combined - extraction of cognitive phenotypes and covariates - for cognitive variables
#' author: amacar
#' date: "Edited from ukb25465_ukb25468_extractPhenotypes_imagingSubset.R; . Modified: `r Sys.time()`"
#' output:
#'   html_document:
#'     toc: true
#'     toc_depth: 2
#'     theme: "flatly"
#'     highlight: "textmate"
#'   pdf_document:
#'     keep_tex: true
#' ---

library(pander)
library(ggplot2)
library(car)
library(dplyr)
library(tidyr)
library(grid)
library(gridExtra)
library(knitr)
library("GGally")

opts_chunk$set(include=TRUE, warnings=FALSE, echo=FALSE, results = "asis", tidy=TRUE, width=50, fig.width=8, fig.height=6)

panderOptions('knitr.auto.asis',FALSE)
# opts_chunk$set(dev="pdf", 
#                dev.args=list(type="cairo"),
#                dpi=96)
set.seed(50)

release="ukb25465_ukb25468"
# Define directories, dependent on  system
if (Sys.info()['sysname']=='Windows') {dir="P://workspaces/"} else {dir="/data/workspaces/lag/workspaces/"}
primary_dir=paste(dir,"lg-ukbiobank/primary_data",sep="")
working_dir=paste(dir,"lg-ukbiobank/working_data/amaia/",sep="")
out_dir=paste(working_dir,"/pheno_files/genetic_v2/",release,sep="")
out_plots_dir=paste(working_dir,"/pheno_files/genetic_v2/",release,"/plots/",sep="")
# Set working directory
opts_knit$set(root.dir = working_dir)
# create ouput directories if they don't already exist
dir.create(file.path(out_dir), showWarnings = FALSE)
dir.create(file.path(out_plots_dir), showWarnings = FALSE)

# get functions for some plots, and relatedness removal
source(paste(working_dir,"/../../analysis/amaia/phenotypes/AIfunctions.R",sep=""))
source(paste(working_dir,"/../../analysis/amaia/phenotypes/relatedness_functions.R",sep=""))
# Start
setwd(working_dir)
input_file=paste(working_dir,"/demographic/",release,"_sqc_fam_selectedPhenotypes_allSamples.csv",sep="")
data<-read.csv(input_file,header=TRUE,stringsAsFactors = FALSE) #nrows=3
colnames_old<-colnames(data)
# edit colnames, remove points, etc
# clean colnames to match fields...
colnames(data)<-gsub("_.left.","_left",gsub("_.right.","_right",colnames(data)))
colnames(data)<-gsub("_\\.|\\+|\\.","_",colnames(data))
colnames(data)<-gsub("^f_","f.",colnames(data))
colnames(data)<-gsub("_$","",gsub("__","_",colnames(data)))

# info about variables, colnames etc
data_cols<-read.csv(gsub("imagingSamples|allSamples","variables",input_file),header=TRUE,stringsAsFactors = FALSE) #nrows=3
#' QC information
sqc_info<-read.csv(paste(primary_dir,"/genetic_data/release_May2017/data/sqc_info.csv",sep=""),sep=",",quote = "\"",header=FALSE)
#' Flag samples to keep after QC (genetic_data)
#' QC-ed samples, from genetic analysis:
sampleQC<-read.table("/data/workspaces/lag/workspaces/lg-ukbiobank/working_data/genetic_data/sample_QC/ukb9246_BritishUKB.samples2keep")
sampleQC$sample2keep<-1
data<-merge(data,sampleQC[,c(-2)],by.x="f.eid",by.y="V1",all.x=TRUE)
data$sample2keep[is.na(data$sample2keep)]<-0
#' White British ancestry
pander(table(data$in_white_British_ancestry_subset)) # 409703 white British ancestry
#' Genetic QC
pander(table(data$sample2keep)) # 408256 sample2keep, white british + genetic sample QC
pander(table(sample2keep=data$sample2keep,whiteBritish=data$in_white_British_ancestry_subset)) 
#' There are 1447 whiteBritish that are not kept due to QC

## Flag samples used in the imaging analysis
imaging_samples<-read.table(paste(release,"_sampleQC_v2_imaging_not_related.samples",sep=""),header=FALSE)
colnames(imaging_samples)<-c("f.eid","ID2")
# add flag for imaging samples after QC:
imaging_samples$imaging_T1_QC<-1

# merge with data, to get flag
data<-merge(data,imaging_samples,by="f.eid",all=TRUE)
data$imaging_T1_QC[is.na(data$imaging_T1_QC)]<-0
table(data$imaging_T1_QC) # 18057 QCed imaging T1
rm(imaging_samples)

#' Flag individuals that contain phenotypes of interest:
#' Imaging/ Fluid Intelligence/ Education
## fluid intelligence score: f.20016_0_0, Fluid_intelligence_score
## education: f.6138_0_0
## check missingness
table(is.na(data$f.20016_0_0))
table(is.na(data$f.6138_0_0))
table(is.na(data$f.20016_0_0)&is.na(data$f.6138_0_0))
table(!is.na(data$f.20016_0_0)&is.na(data$f.6138_0_0))
table(!is.na(data$f.20016_0_0)&!is.na(data$f.6138_0_0))
#' So all the non NA for intelligence has data for EA: 165486
#' Create a flag for CognPhenos: 1 data present/0 not. Be inclusive, i.e. include all the EA
data$Cogn_info<-0
data$Cogn_info[!is.na(data$f.20016_0_0)|!is.na(data$f.6138_0_0)]<-1

## Relatedness

#' If relatedness filtering based on cognitive phenotypes has not been run yet, run and save files
#' Else, read already saved files to filter based on those.
#' Since the filtering is a random process, avoid running this step several times, since it may end up excluding different individuals
#' Although I tried to avoid this by defining the seed in the beginning of the script: set.seed(50)

#' Take into account the samples that were already removed when selecting unrelated individuals for the imaging subset.
#' Otherwise, there will be samples removed that had T1 based on the random process (and the number of individuals with both imaging and cognitive variables will be ~2000 samples smaller).
if (file.exists(paste(release,"_sampleQC_v2_allCogn_not_related.samples",sep=""))==FALSE){
  #' Check relatedness. Relationship pairs, if Kinship coeff > 0.0442
  rel_pairs<-read.table(paste(primary_dir,"/genetic_data/release_May2017/data/ukb1606_rel_s488366.dat",sep=""),stringsAsFactors = FALSE,header=TRUE)

  # include only pairs if both are within present after QC and contain phenotypes of interest
  w1<-which(rel_pairs$ID1 %in% data$f.eid[data$sample2keep==1&data$Cogn_info==1])
  w2<-which(rel_pairs$ID2 %in% data$f.eid[data$sample2keep==1&data$Cogn_info==1])
  w<-intersect(w1,w2) # 91095 related individuals
  rel_pairs_all<-rel_pairs[w,]
  table((rel_pairs_all$ID1 %in% data$f.eid) & (rel_pairs_all$ID2 %in% data$f.eid) )
  samples2check<-unique(c(rel_pairs_all$ID1,rel_pairs_all$ID2)) # 127386 unique samples
  rm(w1,w2,w)
  #' Kinship coefficient within imaging subset
  hist(rel_pairs_all$Kinship,breaks=100)
  summary(rel_pairs_all$Kinship)
  table(rel_pairs_all$Kinship> 0.0625) # 58142 third degree relatives
  # after removing -> 33008 individuals
  
  # create new variable to flag samples with some related sample
  data$sampleRel<-0
  data$sampleRel[which(data$f.eid %in% samples2check)]<-1
  table(data$sample2keep)
  pander(table(sample2keep=data$sample2keep,relatedness=data$sampleRel))
  ## check types of relatedness
  #' pairs of MZ twins:
  dim(subset(rel_pairs_all,Kinship>0.4)) # 157
  subset(data,f.eid=="3429006"|f.eid=="3047939")[,c("f.eid","imaging","imaging_T1","sample2keep","sampleRel")]
  #' pairs of 1st degree relatives
  dim(subset(rel_pairs_all,Kinship>0.15&Kinship<0.4))[1] # 24768 pairs of 1st degree relatives
  #' pairs of 2nd degree relatives
  dim(subset(rel_pairs_all,Kinship>0.06&Kinship<0.15))[1] # 37865 pairs
  #' Need to remove: 1 relative per pair
  # data contains 407814 samples where sample2keep==1 and at least one of the two cognitive phenotypes is present
  data1<-related_clean(sdata=subset(data,sample2keep==1&Cogn_info==1),rel_info=rel_pairs_all)
  table(data1$imaging_excl_sampleRel)
  table(data1$RelRemove1)
  table(is.na(data1$RelRemove1))
  table(is.na(data1$RelRemove),data1$imaging_T1)
  
  #' Number of samples to remove due to relatedness
  #' 73851 to remove based on relatedness
  #' 333963 clean QCed samples
  all_samples2exclude<-subset(data1,RelRemove1==1)[,c("f.eid","f.eid")]# 73851 samples to exclude due to relatedness
  # save file, for reference:
  write.table(all_samples2exclude,file=paste(release,"_related_allCogn_exclude.samples",sep=""),sep="\t",quote=FALSE,col.names = FALSE,row.names = FALSE)
  #' Save list of samples belonging to the "white British ancestry subset"
  samples2keep<-subset(data1,sample2keep==1&is.na(RelRemove1))[,c("f.eid","f.eid")] # 333963
  write.table(samples2keep,file=paste(release,"_sampleQC_v2_allCogn_not_related.samples",sep=""),sep="\t",quote=FALSE,col.names = FALSE,row.names = FALSE)
  rm(samples2keep,all_samples2exclude)
  
  #----------------------------
  ## create new flag: for non-imaging subset: i.e. including only non-imaging samples into the relatedness removal
  # include only pairs if both are within present after QC and contain phenotypes of interest
  # and if they have not used in imaging subset
  w1<-which(rel_pairs$ID1 %in% data$f.eid[data$imaging_T1_QC==0&(data$sample2keep==1&data$Cogn_info==1)])
  w2<-which(rel_pairs$ID2 %in% data$f.eid[data$imaging_T1_QC==0&(data$sample2keep==1&data$Cogn_info==1)])
  w<-intersect(w1,w2) # 83684 related individuals
  rel_pairs_nonT1<-rel_pairs[w,]
  table((rel_pairs_nonT1$ID1 %in% data$f.eid) & (rel_pairs_nonT1$ID2 %in% data$f.eid) )
  samples2check_nonT1<-unique(c(rel_pairs_nonT1$ID1,rel_pairs_nonT1$ID2)) # 127386 unique samples
  rm(w1,w2,w)
  # create new variable to flag samples with some related sample
  data$sampleRel_nonT1<-0
  data$sampleRel_nonT1[which(data$f.eid %in% samples2check_nonT1)]<-1
  table(data$sampleRel_nonT1,data$sampleRel)
  # clean relatedness, for non imaging subset # ***hemen nago***
  data_non_imaging<-subset(data,imaging_T1_QC==0&(sample2keep==1&Cogn_info==1))[,c("f.eid","imaging_T1_QC")]
  data1_non_imaging<-related_clean(sdata=data_non_imaging,rel_info=rel_pairs_nonT1)
  table( data1_non_imaging$RelRemove1)
  colnames(data1_non_imaging)<-c("f.eid","imaging_T1_QC","RelRemove1_nonT1")
  rm(data_non_imaging)
  
  #' Number of samples to remove due to relatedness
  all_samples2exclude_non_imaging<-subset(data1_non_imaging,RelRemove1_nonT1==1)[,c("f.eid","f.eid")]# 73851 samples to exclude due to relatedness
  # save file, for reference:
  write.table(all_samples2exclude_non_imaging,file=paste(release,"_related_allCogn_nonT1_exclude.samples",sep=""),sep="\t",quote=FALSE,col.names = FALSE,row.names = FALSE)
  #' Save list of samples belonging to the "white British ancestry subset"
  samples2keep_non_imaging<-subset(data1_non_imaging,is.na(RelRemove1_nonT1))[,c("f.eid","f.eid")] # 333963
  write.table(samples2keep_non_imaging,file=paste(release,"_sampleQC_v2_allCogn_nonT1_not_related.samples",sep=""),sep="\t",quote=FALSE,col.names = FALSE,row.names = FALSE)
  rm(samples2keep_non_imaging,all_samples2exclude_non_imaging)
  
  
  #-----------------------------------------
  # combine with data1
  data2<-merge(data1,data1_non_imaging,all=TRUE)
  # do not subset data to exclude related individuals:
  # because not the same subsets will be for: imaging_T1; imaging_nonT1; all
  table(data2$imaging_T1_QC) # 18057-88 (i.e. 88 imaging individuals that did not have cognitive phenos)
  table(is.na(data2$RelRemove1)) # considering all possible samples
  table(is.na(data2$RelRemove1),is.na(data2$RelRemove1_nonT1)) # considering only non T1 samples
  # create clean flags for the three subsets of the data
  ## imaging_T1_QC: already exists, and it does not include related individuals 
  ## non_imaging_T1_QC
  data2$all_QC<-data2$non_imaging_T1_QC<-0 
  # flag individuals after QC
  ## for total sample
  w<-which(is.na(data2$RelRemove1))
  data2$all_QC[w]<-1
  rm(w)
  ## for nonimaging sample
  w<-which(data2$imaging_T1_QC==0&is.na(data2$RelRemove1_nonT1))
  data2$non_imaging_T1_QC[w]<-1
  rm(w)
  
  #----------------------------
  #' Save data2 as data, and remove intermediate files
  data<-subset(data2,all_QC==1|non_imaging_T1_QC==1|imaging_T1_QC==1)
  rm(data1,data1_non_imaging,data2)
  #----------------------------
  rm(rel_pairs,rel_pairs_all,rel_pairs_nonT1)
  
} else {
  # imaging_T1_QC will already be there, so need to define:
  ## all_QC
  samples2keep_all<-read.table(paste(release,"_sampleQC_v2_allCogn_not_related.samples",sep=""))[,1]
  samples2keep_all<-data.frame(samples2keep_all); colnames(samples2keep_all)<-c("f.eid")
  samples2keep_all$all_QC<-1
  ## non_imaging_T1_QC
  samples2keep_nonT1<-read.table(paste(release,"_sampleQC_v2_allCogn_nonT1_not_related.samples",sep=""))[,1]
  samples2keep_nonT1<-data.frame(samples2keep_nonT1); colnames(samples2keep_nonT1)<-c("f.eid")
  samples2keep_nonT1$non_imaging_T1_QC<-1
  # combine both
  samples2keep2<-merge(samples2keep_all,samples2keep_nonT1,all=TRUE)
  samples2keep2$all_QC[which(is.na(samples2keep2$all_QC))]<-0
  samples2keep2$non_imaging_T1_QC[which(is.na(samples2keep2$non_imaging_T1_QC))]<-0
  #
  data2<-merge(data,samples2keep2,all=TRUE)
  rm(samples2keep_all,samples2keep_nonT1,samples2keep2)
}


#----------------------------------------------------------
# Phenotypes of interest: 
## fluid intelligence score: f.20016_0_0, Fluid_intelligence_score
## education: f.6138_0_0


#' check how many of the kepts samples contain imaging data as well
table(data$imaging_T1_QC) #
table(imaging=data$imaging_T1_QC,intelligence=!is.na(data$f.20016_0_0)) 
table(imaging=data$imaging_T1_QC,intellignece=!is.na(data$f.20016_2_0)) #imaging visit for fluid intelligence
#' 16238 samples with imaging and with intelligence
table(imaging=data$imaging_T1_QC,EA=!is.na(data$f.6138_0_0)) 
table(imaging=data$imaging_T1_QC,EA=!is.na(data$f.6138_2_0))
#' 17969 samples with imaging and with educational data, only 17847 with data from imaging visit

#' Follow Lee et al. (2018)/ Ge et al. (2018)
#' to map from Qualifications to EduYears
#' code EduYears
#' UK Biobank | UK Biobank qualifications | ISCED  level| US years of coding schooling
# coding for education, to impute years of education
# 1	College or University degree - 
# 2	A levels/AS levels or equivalent
# 3	O levels/GCSEs or equivalent
# 4	CSEs or equivalent
# 5	NVQ or HND or HNC or equivalent
# 6	Other professional qualifications eg: nursing, teaching
# -7 None of the above 1 7
# -3 Prefer not to answer - -


#' There are three entries for Qualifications: f.6138.0.0, f.6138.1.0, f.6138.2.0
# work with 0.0 or highest?
summary(data$f.6138_0_0) # NO NA's
table(is.na(data$f.6138_0_0)) 
head(data[,grep("f.6138.|Qualific",colnames(data))])

#' Create EduYears variable
data$EduYears02<-data$EduYears02<-data$EduYears00<-NA
# 0.0
data$EduYears00[data$f.6138_0_0==-7]<-7 # ISCED 1
data$EduYears00[data$f.6138_0_0==6]<-15 # ISCED 4
data$EduYears00[data$f.6138_0_0==5]<-19 # ISCED 5 (19)
data$EduYears00[data$f.6138_0_0==4]<-10 # ISCED 2
data$EduYears00[data$f.6138_0_0==3]<-10 # ISCED 2
data$EduYears00[data$f.6138_0_0==2]<-13 # ISCED 3
data$EduYears00[data$f.6138_0_0==1]<-20 # ISCED 5 (20)
# 1.0
data$EduYears10[data$f.6138_1_0==-7]<-7 # ISCED 1
data$EduYears10[data$f.6138_1_0==6]<-15 # ISCED 4
data$EduYears10[data$f.6138_1_0==5]<-19 # ISCED 5 (19)
data$EduYears10[data$f.6138_1_0==4]<-10 # ISCED 2
data$EduYears10[data$f.6138_1_0==3]<-10 # ISCED 2
data$EduYears10[data$f.6138_1_0==2]<-13 # ISCED 3
data$EduYears10[data$f.6138_1_0==1]<-20 # ISCED 5 (20)
# 2.0
data$EduYears20[data$f.6138_2_0==-7]<-7 # ISCED 1
data$EduYears20[data$f.6138_2_0==6]<-15 # ISCED 4
data$EduYears20[data$f.6138_2_0==5]<-19 # ISCED 5 (19)
data$EduYears20[data$f.6138_2_0==4]<-10 # ISCED 2
data$EduYears20[data$f.6138_2_0==3]<-10 # ISCED 2
data$EduYears20[data$f.6138_2_0==2]<-13 # ISCED 3
data$EduYears20[data$f.6138_2_0==1]<-20 # ISCED 5 (20)
#
table(data$EduYears00)

#' Combine them into a single one
#' Assign highest category to respondents who selected multiple options
data$EduYears_collapsed<-apply(data[,c("EduYears00", "EduYears10", "EduYears20")],1,function(x){
  x2<-x[!is.na(x)]
  y<-paste(x2,collapse=",")
  return(y)
})
data$EduYears_max<-apply(data[,c("EduYears00", "EduYears10", "EduYears20")],1,function(x){
  if (sum(!is.na(x))>0){
  m<-max(x,na.rm=TRUE)
  } else { m<-NA }
  return(m)
})

table(data$EduYears_collapsed)
table(data$EduYears_max)
plot(data$EduYears00,data$EduYears_max)
hist(data$EduYears00)
hist(data$EduYears_max) 
#' check summary stats
#' From Lee et al. (2018): mean(SD) 14.87(5.11), N=442,183, did not exclude based on relatedness -> used BOLT-LMM for association
#' From Ge et al. (2018): mean(SD) 14.8 (5.1), N=332,613
summary(data$EduYears00); sd(data$EduYears00,na.rm=TRUE)
summary(data$EduYears_max); sd(data$EduYears_max,na.rm=TRUE)

hist(data$EduYears00) 
#' This is very non-normal...

table(data$sex[!is.na(data$EduYears00)])
summary(data$f.21003_0_0[!is.na(data$EduYears00)]) # age at instance 00

#----------------------------------------------------------
#' Fluid intelligence score
#' From Ge et al.: mean(SD) 6.2(2.1), N=108,147
table(is.na(data$f.20016_0_0))
hist(data$f.20016_0_0)
summary(data$f.20016_0_0); sd(data$f.20016_0_0,na.rm=TRUE)
#----------------------------------------------------------
# relationship between Fluid Intelligence and EduYears
boxplot(data$f.20016_0_0~data$EduYears00)
#----------------------------------------------------------
#' Define columns containing imaging variables of interest
phenos_all<-colnames(data)[grep("f.20016_|f.6138_|Edu|intelligence|Qual",colnames(data))]
# phenos<-colnames(data)[grep("f.20016_0_0|f.6138_0_0|^EduYears00",colnames(data))]
#' Define phenotypes of interest
phenos<-c("EduYears00","EduYears_max","f.20016_0_0","f.20016_2_0","Fluid_intelligence_score") 
# where 00 indicates first assessment only, max=maximum across instances, and Fluid_intelligence_score is the average across non NA instances
# where 02 indicates assesment in 
#----------------------------------------------------------
#' ## Define 'default covariates'  
#' 
#' PCs, available from release  
# check visually
# plot(data[,"PC1"],data[,"PC2"])
# plot(data[data$in_white_British_ancestry_subset==1,"PC1"],data[data$in_white_British_ancestry_subset==1,"PC2"])
pcs<-colnames(data)[grep("^PC",colnames(data))]
##if 1-40 PCs were to be included, edit here
#' Extract all 40PCs, will define how many to include within the lm in next steps
maxPC=40
# define first 10 PCs
pcs<-pcs[1:maxPC]
print("is.na PC1'")
table(is.na(data[,pcs[2]])) 
#' sample QC has been completed, so all contain PCs

#' Age at assesment, instance 2 --> age at imaging visit
age1<-colnames(data)[grep("^Age_when",colnames(data))] # this is average age across three visits
age<-"f.21003_0_0" # age at initial visit
age_imaging<-"f.21003_2_0"
print("is.na age")
table(is.na(data[,age])) # no missing data
table(is.na(data[,age_imaging]),data$imaging_T1_QC) # no missing data for the imaging subset
# compute age2, as recommended: (age-meanAge)2
meanAge_all<-mean(data[(data$all_QC==1),age],na.rm = TRUE) # 56.86607
meanAge_imagingT1<-mean(data[(data$imaging_T1_QC==1),age_imaging],na.rm=TRUE) # 62.77739
meanAge_non_imagingT1<-mean(data[(data$non_imaging_T1_QC==1),age]) # 56.95803
# creat zage2 variables per group
data$zage2_all<-data$zage2_imagingT1<-data$zage2_non_imagingT1<-NA
# fill in those variables
data$zage2_all[data$all_QC==1]<-(data[data$all_QC==1,age]-meanAge_all)^2 # age2: (age-meanAge)^2
data$zage2_imagingT1[data$imaging_T1_QC==1]<-(data[data$imaging_T1_QC==1,age_imaging]-meanAge_imagingT1)^2
data$zage2_non_imagingT1[data$non_imaging_T1_QC==1]<-(data[data$non_imaging_T1_QC==1,age]-meanAge_non_imagingT1)^2
rm(meanAge_all,meanAge_imagingT1,meanAge_non_imagingT1)
zage2_cols<-colnames(data)[grep("zage2",colnames(data))]

#' Assesment center
colnames(data)[grep("UK_Biobank_assessment|f.54",colnames(data))] # missing, need to add!
# UK_Biobank_assessment_centre
assesmentC<-"f.54_0_0"
assesmentC_imaging<-"f.54_2_0"
#' Sex  
sex<-"sex"
# boxplot(data[,age]~data[,sex]) # males are older than females on average
print("is.na sex")
pander(table(is.na(data[,sex]))) 
#' Some missing data, same as people without Genetic info.  
#' Select Genotyping array and batch.  
array<-colnames(data)[grep("array",colnames(data))]
batch<-colnames(data)[grep("^batch",colnames(data))]

#' Define binary and quantitative covariates
bin_covs<-c(sex,assesmentC,assesmentC_imaging,array,batch)
quant_covs<-c(age,age_imaging,zage2_cols,pcs)

#' Visually inspect covariates to see whether there are outliers # if so, should blank?!
quant_covs_hist_out<-lapply(quant_covs,function(x){histo_out_AIs(data2=subset(data,sample2keep==1),x,thr=6)})
# save plot
hists<-sapply(quant_covs_hist_out,"[[",2)
pdf(paste(out_plots_dir,"quant_covariates_outliers_Cogn_all.pdf",sep=""),onefile = TRUE)
for (h in 1:length(hists)){ grid.arrange(hists[[h]]) }
dev.off()
rm(hists,quant_covs_hist_out)

#----------------------------------------------------------
# select phenotypes and samples2keep
flags4subsets<-c("all_QC","imaging_T1_QC","non_imaging_T1_QC")
data_s<-data[,c("f.eid","f.eid",bin_covs,quant_covs,phenos,flags4subsets)]
# number of samples per subset:
table(data_s$all_QC) # 333963
table(data_s$imaging_T1_QC) # 17969
table(data_s$non_imaging_T1_QC) # 321776

# new colnames
phenos2<-c("EduYears00","EduYears_max","FluidInt00","FluidInt_imaging","FluidInt_mean")
bin_covs2<-c("sex","assessment_centre","assessment_centre_imaging","array","batch")
quant_covs2<-c("age","age_imaging",zage2_cols,pcs)
if (
    length(colnames(data_s))==
    length(c("ID1","ID2",bin_covs2,quant_covs2,phenos2,flags4subsets))
) { 
  colnames(data_s)<-c("ID1","ID2",bin_covs2,quant_covs2,phenos2,flags4subsets)
}

# define covariate order as well
## per subset
covs_order_all<-c("assessment_centre","age","zage2_all","sex",pcs,"array")
covs_order_imagingT1<-c("assessment_centre_imaging","age_imaging","zage2_imagingT1",pcs,"array")
covs_order_non_imagingT1<-c("assessment_centre","age","zage2_non_imagingT1",pcs,"array")

#' Convert binary covariates into factors
data_s[,bin_covs2]<-apply(data_s[,bin_covs2],2,function(x) as.factor(x))

#' Save files: imaging and nonimaging subsets
write.table(data_s,file=paste(out_dir,"_PhenotypesCogn_Covs_flags4subsets.table",col.names=TRUE,row.names=FALSE))

#' Clean: remove other objects that are not necessary for further steps.
rm(data,data_cols)
rm(sampleQC,sqc_info,samples2check)
rm(age,age1,age_imaging,array,assesmentC,assesmentC_imaging,sex,batch)
rm(zage2_cols)
rm(bin_covs,bin_covs2,quant_covs,quant_covs2)
rm(phenos,phenos2,phenos_all)
rm(colnames_old)
rm(h)

#' Next: run ukb25465_ukb25468_LM_residuals_Cogn_all.R to run LMs and generate residualized phenotypes.
#' and ukb25465_ukb25468_LM_residuals_Cogn_imaging.R to run LMs and generate residualized phenotypes for the imaging subset.

# save image
save.image(paste(out_dir,"/",release,"_extractPhenotypesCogn_PhenoCovs_subsets.RData",sep=""))
