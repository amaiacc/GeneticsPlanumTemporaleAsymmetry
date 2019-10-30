# association_phenotypes4BGENIE.R
# date: 13.10.2017; edited: 09.05.2018; edited 23.05.2018, edited for PT 08.06.2018
# create phenotype input files for BGENIE (which do not have sample IDs, and need to match the order within the BGEN file)
# genetic data: v3; imputed data, 
## includes chromosome X, and XY, but sample files are different for these chromosomes, so new input files will be required
library(data.table)

#------------------------------------
# Get command line arguments
args = commandArgs(trailingOnly=TRUE)

# args=c("/data/clusterfs/lag/users/amacar/ukb/input/imp/ukb_imp_chrX_v3_imagingT1_N12245.sample",
#        "X",
#        "ukb21288_ukb21293",
#        "imagingT1_N12245")

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

ukb9246<-read.csv(paste(dir,"lg-ukbiobank/working_data/amaia/demographic/",
                        "ukb9246_sqc_fam_selectedPhenotypes_imagingSamples.csv",sep=""))
ukb9246$ukb9246<-1
ukb9246<-ukb9246[,c("f.eid","ukb9246")]
colnames(ukb9246)<-c("ID1","ukb9246")
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
  pheno_imaging<-merge(pheno_imaging,ukb9246,all=TRUE,by="ID1")
  pheno_imaging$ukb9246[is.na(pheno_imaging$ukb9246)]<-0
  table(pheno_imaging$ukb9246)
  # get phenos:
  phenos<-colnames(pheno_imaging)[grep(paste("residuals*.*",region,sep=""),colnames(pheno_imaging))]
  # combine pheno files
  if (length(grep("X",chr_type))>0) {
        out_file2<-gsub("sample",paste("sample",chr_type,sep=""),out_file)
        } else {out_file2<-out_file }
  if (file.exists(out_file2)==FALSE){
      # make sure that the pheno file matches the order from the sample
      pheno_samples<-merge(sample,pheno_imaging,by.x=c("ID_1","ID_2"),by.y=c("ID1","ID2"),all.x=TRUE,stringsAsFactors=FALSE) #,
      pheno_samples$array<-as.character(pheno_samples$array)
      pheno_samples$batch<-as.character(pheno_samples$batch)
      pheno_samples<-pheno_samples[,c("ID_1","ukb9246",phenos)]
      # check if order is the same
      table(sample$ID_1==pheno_samples$ID_1)
      pheno_samples<-pheno_samples[match(sample$ID_1,pheno_samples$ID_1),]
      table(sample$ID_1==pheno_samples$ID_1)
      
      # create subsets of the data, given: ukb9456 or ukb21288_ukb21293
      pheno_samples[,paste(phenos,"_ukb9246",sep="")]<-pheno_samples[,phenos]
      pheno_samples[,paste(phenos,"_",pheno_root,sep="")]<-pheno_samples[,phenos]
      # blank inds not present within each subset
      w<-which(pheno_samples$ukb9246==1)
      pheno_samples[w,paste(phenos,"_ukb9246",sep="")]<-NA
      pheno_samples[-w,paste(phenos,"_",pheno_root,sep="")]<-NA
      
      # missing data is NA, need to specify this when running BGENIE, replace to -999
      pheno_samples[is.na(pheno_samples)]<-(-999)
      
      # or alternatively use the --miss parameter to code for missingness (NA)
      
      # save file with phenotypes and covariates
      # include just one ID, to be able to double check
      phenos2<-colnames(pheno_samples)[grep(region,colnames(pheno_samples))]
      write.table(pheno_samples[,c("ID_1",phenos2)],file=out_file2,row.names=FALSE,quote=FALSE)
      
    # rm(pheno_samples)
    }
    # clean intermediate files
  rm(out_file2)


rm(out_file,pheno_imaging_file,pheno_imaging,phenos)
  

