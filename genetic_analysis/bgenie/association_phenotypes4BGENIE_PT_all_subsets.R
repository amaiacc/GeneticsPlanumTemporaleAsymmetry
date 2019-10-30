# association_phenotypes4BGENIE.R
# date: 13.10.2017; edited: 09.05.2018; edited 23.05.2018, edited for PT 08.06.2018, 11.01.2019
# create phenotype input files for BGENIE (which do not have sample IDs, and need to match the order within the BGEN file)
# genetic data: v3; imputed data, 
## includes chromosome X, and XY, but sample files are different for these chromosomes, so new input files will be required
library(data.table)

#------------------------------------
# Get command line arguments
args = commandArgs(trailingOnly=TRUE)

# args=c("/data/clusterfs/lag/users/amacar/ukb/input/imp/ukb_imp_chr10_v3_imagingT1_N18057.sample",
#        "10",
#        "ukb25465_ukb25468",
#        "imagingT1_N18057")

# get sample file as argument
sample_file=args[1]
chr_type=args[2] # A for autosomes, X or XY
pheno_root=args[3] # ukb phenotype batch names
subset_name=args[4] # subset of data that will be included,, in case not samples are in; should reflect sample_file info
#------------------------------------


options(stringsAsFactors = FALSE)
#------------------------------------
# Set directories, dependent on  system
if (Sys.info()['sysname']=='Windows') {dir="P://workspaces/"} else {dir="/data/workspaces/lag/workspaces/"}
imaging_dir=paste(dir,"lg-ukbiobank/working_data/amaia/pheno_files/genetic_v2/",pheno_root,"/summary_phenotypes/",sep="")
out_dir_base=paste(dir,"lg-ukbiobank/working_data/amaia/genetic_data/",sep="")
out_dir=paste(out_dir_base,"release_v3/imp/bgenie/",sep="")
dir.create(file.path(out_dir_base,"release_v3/imp/bgenie"),showWarnings = FALSE)

# define different subsets
ukb9246_ukb10785<-read.csv(paste(dir,"lg-ukbiobank/working_data/amaia/demographic/",
                        "ukb9246_ukb10785_sqc_fam_selectedPhenotypes_imagingSamples.csv",sep=""))
ukb9246_ukb10785$ukb9246_ukb10785<-1
ukb9246_ukb10785<-ukb9246_ukb10785[,c("f.eid","ukb9246_ukb10785")]
colnames(ukb9246_ukb10785)<-c("ID1","ukb9246_ukb10785")
# 
ukb21288_ukb21293<-read.csv(paste(dir,"lg-ukbiobank/working_data/amaia/demographic/",
                                  "ukb21288_ukb21293_sqc_fam_selectedPhenotypes_imagingSamples.csv",sep=""))
ukb21288_ukb21293$ukb21288_ukb21293<-1
ukb21288_ukb21293<-ukb21288_ukb21293[,c("f.eid","ukb21288_ukb21293")]
colnames(ukb21288_ukb21293)<-c("ID1","ukb21288_ukb21293")
#
ukb25465_ukb25468<-read.csv(paste(dir,"lg-ukbiobank/working_data/amaia/demographic/",
                                  "ukb25465_ukb25468_sqc_fam_selectedPhenotypes_imagingSamples.csv",sep=""))
ukb25465_ukb25468$ukb25465_ukb25468<-1
ukb25465_ukb25468<-ukb25465_ukb25468[,c("f.eid","ukb25465_ukb25468")]
colnames(ukb25465_ukb25468)<-c("ID1","ukb25465_ukb25468")

#
all_subsets<-merge(ukb9246_ukb10785,ukb21288_ukb21293,all=TRUE)
all_subsets<-merge(all_subsets,ukb25465_ukb25468,all=TRUE)
# replace NAs with 0s
all_subsets[is.na(all_subsets)]<-0
# for each subset only flag the new added samples
all_subsets$ukb21288_ukb21293[all_subsets$ukb9246_ukb10785==1]<-0
all_subsets$ukb25465_ukb25468[all_subsets$ukb9246_ukb10785==1|all_subsets$ukb21288_ukb21293==1]<-0
# combine first two subsets into one
all_subsets$ukb9246_ukb21293<-0 # combined v2, previous analysis = all
all_subsets$ukb9246_ukb21293[all_subsets$ukb9246_ukb10785==1|all_subsets$ukb21288_ukb21293==1]<-1
# counts per subset
apply(all_subsets[,-1],2,table)
# clean subsets
rm(ukb9246_ukb10785,ukb21288_ukb21293,ukb25465_ukb25468)
#------------------------------------
# Sample files,  read
#------------------------------------

## read
sample<-read.table(sample_file,header=TRUE) # actually, use v3: they have the same number of individuals for the autosomes
# some checks
table(sample$missing)
#table(sample$sex)


# remove first line which is only 0
sample[1,]
sample<-sample[-1,1:2]

#------------------------------------
# Phenotype file for imaging phenotypes, output from: 
## ukb21288_ukb21293_LM_residuals_imagingSubset.R
#------------------------------------
covs_name="_noBioCovs_noAssessmentC"
region="Planum_Temporale"
type="Volume"
#------------------------------------
  # define files
  pheno_imaging_file<-paste(imaging_dir,list.files(imaging_dir,pattern=paste("imaging",covs_name,type,"Phenotypes_residuals_wHeader.table",sep="")),sep="")
  out_file=paste(out_dir,pheno_root,"_",type,"_sample_",subset_name,"_",region,"_phenos4BGENIE.table",sep="")
  #  read data
  pheno_imaging<-read.table(pheno_imaging_file,header=TRUE)
  # add col to mark wether samples belong to previous release
  pheno_imaging<-merge(pheno_imaging,all_subsets,all.x=TRUE,by="ID1")
  subset_cols<-colnames(pheno_imaging)[grep("ukb",colnames(pheno_imaging))]
  apply(pheno_imaging[,grep("ukb",colnames(pheno_imaging))],2,table)
  # get phenos:
  phenos<-colnames(pheno_imaging)[grep(paste("residuals*.*",region,sep=""),colnames(pheno_imaging))]
  # create subsets of the data
  ## if no subset is specified: total N will be available
  for (s in subset_cols){
    # copy all first
    pheno_imaging[,paste(phenos,"_",s,sep="")]<-pheno_imaging[,phenos]
    # then, blank inds not present within each subset
    w<-which(pheno_imaging[,s]==0)
    pheno_imaging[w,paste(phenos,"_",s,sep="")]<-NA
    rm(w)
  }
  rm(s)
  # check number of non NAs per phenotype
  apply(pheno_imaging[,grep("ukb",colnames(pheno_imaging))],2,function (x) table(!is.na(x)))
  
  # combine pheno files
  if (length(grep("X",chr_type))>0) {
        out_file2<-gsub("sample",paste("sample",chr_type,sep=""),out_file)
        } else {out_file2<-out_file }
  if (file.exists(out_file2)==FALSE){
      # make sure that the pheno file matches the order from the sample
      pheno_samples<-merge(sample,pheno_imaging,by.x=c("ID_1","ID_2"),by.y=c("ID1","ID2"),all.x=TRUE,stringsAsFactors=FALSE) #,
      pheno_samples$array<-as.character(pheno_samples$array)
      pheno_samples$batch<-as.character(pheno_samples$batch)
      pheno_cols<-colnames(pheno_samples)[grep(paste(phenos,collapse="|"),colnames(pheno_samples))]
      pheno_samples<-pheno_samples[,c("ID_1",pheno_cols)]
      rm(pheno_cols)
      # check if order is the same
      table(sample$ID_1==pheno_samples$ID_1)
      pheno_samples<-pheno_samples[match(sample$ID_1,pheno_samples$ID_1),]
      table(sample$ID_1==pheno_samples$ID_1)
      
      # missing data is NA, need to specify this when running BGENIE, replace to -999
      pheno_samples[is.na(pheno_samples)]<-(-999)
      
      # or alternatively use the --miss parameter to code for missingness (NA)
      
      # save file with phenotypes and covariates
      # include just one ID, to be able to double check
      phenos2<-colnames(pheno_samples)[grep(region,colnames(pheno_samples))]
      write.table(pheno_samples[,c("ID_1",phenos2)],file=out_file2,row.names=FALSE,quote=FALSE)
      rm(phenos2)      
    rm(pheno_samples)
    }
    # clean intermediate files
  rm(out_file2)


rm(out_file,pheno_imaging_file,pheno_imaging,phenos)
rm(all_subsets,subset_cols)
  

