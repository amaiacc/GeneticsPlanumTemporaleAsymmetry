if("qqman" %in% rownames(installed.packages()) == FALSE) {install.packages("qqman")}
if("gridExtra" %in% rownames(installed.packages()) == FALSE) {install.packages("gridExtra")}
library(qqman)
library(gridExtra); library(grid)
library(ggplot2)
library("data.table")
library(dplyr);library(tidyr)
options(stringsAsFactors = FALSE) #,digits=22
# set number of digits to 22, default is 7, but this affects calculation of lambda, which always becomes 1

#------------------------------------------
# plots, some general parameters
#------------------------------------------
cols=c("#d95f02","#1f78b4","#33a02c","#000000") # '#1b9e77','#d95f02','#7570b3 # "#a6cee3"  ,"#7570b3"
mytheme<-theme_bw() + theme(panel.spacing = unit(0, "lines"), 
                            strip.background = element_rect(fill="white"), strip.text = element_text(size=16), 
                            # axis.title.x=element_blank(), axis.ticks.x =element_blank(), # axis.text.x =element_blank(),
                            # axis.text.x = element_text(size=16,colour="black",hjust=1,vjust=.5),
                            # axis.title.y=element_text(size=16),axis.text.y=element_text(size=16,colour="black"),
                            axis.text.x = element_text(angle = 45, hjust = 1),
                            title=element_text(size=16),
                            axis.title=element_text(size=16),axis.text=element_text(size=16,colour="black"),
                            legend.text=element_text(size=16), legend.title =element_text(size=16) ) 


#---------------------------------------------------#
# args=commandArgs()
# file=args[1]
# p_col=args[2]
#---------------------------------------------------#

if (Sys.info()['sysname']=='Windows') {dir="P://workspaces/"} else {dir="/data/workspaces/lag/workspaces/"}
subset_name="imagingT1_N18057"

primary_dir=paste(dir,"/lg-ukbiobank/working_data/amaia/genetic_data/release_v3/imp/bgenie/output/PT/filtered/",sep="")
working_dir=paste(dir,"/lg-ukbiobank/working_data/amaia/genetic_data/release_v3/imp/bgenie/output/PT/",subset_name,"/",sep="")
# read filtered snps based on HWE<1e-07 and INFO<0.7 (INFO from QCtool - total sample, not BGENIE!)
qc_dir=paste(dir,"/lg-ukbiobank/working_data/amaia/genetic_data/release_v3/imp/subset_",subset_name,"/QC/",sep="")
#
if (!file.exists(working_dir)){
  dir.create(file.path(working_dir))
}

setwd(working_dir) 

# create ./plots ./clean directories if they do not exist
if (!file.exists("plots")){
  dir.create(file.path("plots"))
}

if (!file.exists("clean")){
  dir.create(file.path("clean"))
}

# Parse arguments:
## chr, phenotype
#------------------------------------
# Get command line arguments
# args = commandArgs(trailingOnly=TRUE)

# arguments for test round
region="Planum_Temporale"
phenos2test<-c(paste("AI_VOLUME",region,sep="_"),
               paste("Volume_of_grey_matter_in",region,c("left","right"),sep="_")             )

for (phenotype in phenos2test){
  
  args<-c("all",phenotype,
          "ukb25465_ukb25468",
          subset_name,
          "Volume",
          region)
  
  chr=args[1] #either all, or a specific chr name
  pheno=args[2] # phenotype name
  pheno_root=args[3] # ukb phenotype batch names
  # subset_name=args[4] # subset of data that will be included,, in case not samples are in; should reflect sample_file info
  type=args[5]
  region=args[6]
  #------------------------------------
  # define thresholds
  maf_thr=0.001 # args[6]
  info_thr=0.7 # args[7]
  hwe_thr=1e-7
  filter_pattern=paste("_",hwe_thr,"HWEp_",info_thr,"INFO_",maf_thr,"MAF",sep="")
  #------------------------------------
  rm(args)
  
  # define phenotype name
  name0=gsub("_"," ",gsub("residuals_","",pheno))
  #------------------------------------
  for (sex in c("females","males")){ 
    name=paste(name0,sex,sep=" ")
    print("Check if Manhattan plot (or some other output file) exists, only proceed if it does not")
    if (file.exists(paste("plots/Manhattan_",gsub(" ","_",name),"_CHR",chr,".png",sep=""))==FALSE){
      
    
    if (chr!="all") {c=paste(chr,"\\.",sep="")} else {c=".*."} # just one chr!
    # define files to parse, will be either one (if one chr, or all)
    pattern1=paste(pheno_root,subset_name,"bgenie",type,region,sex,"chr",sep="_")
    pattern4files<-paste(pattern1,c,sep="")
    files<-list.files(primary_dir,pattern=pattern4files)
    rm(pattern4files,c)
      
    cat('Info filter threshold: ',info_thr,'\n')
    cat('Frequency filter threshold:',maf_thr,'\n')
    cat('HWE p-value filter threshold:',hwe_thr,'\n')
    # function to parse files, read one at a time, and subset only cols for phenotype of interest
    for (file in files) {
      # define chromosome
      chrom=gsub(".out.gz","",gsub(filter_pattern,"",gsub(pattern1,"",file)))
        #-----------------------------------
        # start processing file
        # cat("Opening conection for chromosome ",chrom,"\nFile: ",file,'\n',sep="")
        zz=gzfile(paste(primary_dir,file,sep=""),'rt')  
        # get header first, and define columns of interest
        cat('Reading header of file...\n')
        header=read.table(zz,header=T,skipNul=TRUE,nrows=1)
        
        # define general cols
        snp_cols<-grep("residuals",colnames(header),invert=TRUE)
        pheno_cols<-grep(gsub("_totalBV","",pheno),colnames(header))
        
        # read file
        cat(paste('Reading file...\n'))
        cat(as.character(Sys.time()),'\n')
        d0=fread(paste('zcat ',paste(primary_dir,file,sep=""),sep=""),header=T,select=c(snp_cols,pheno_cols))
        cat(as.character(Sys.time()),'\n')
        cat('Number of SNPs in choromsome ',chrom,':',NROW(d0),'\n')
        # check filters:
        # maf
        cat('Number of variants with af<',maf_thr,'\n')
        cat(table(d0$af<=maf_thr))
        cat('Number of variants with af>',(1-maf_thr),'\n')
        cat(table(d0$af>=(1-maf_thr)))
        cat('Number of filtered out based on frequency (TRUE) for chromosome',chrom,'\n')
        cat(table(d0$af<=maf_thr|d0$af>=(1-maf_thr)))
        # subset on maf, calculated by BGENIE
        d1<-subset(d0,(af>=maf_thr&af<=(1-maf_thr)))
        cat('Number of variants with HWE pval <',hwe_thr,'\n')
        cat(table(d1$hwe<hwe_thr))
        cat('Number of variants with INFO.qctool  <',info_thr,'\n')
        cat(table(d1$info.qctool<info_thr))
        # double check INFO field from qctool, and hwe
        cat('Filter out variants with UKB.INFO<', info_thr,' or HWE pval <', hwe_thr)
        d2<-subset(d1,INFO.UKB>info_thr&HW_exact_p_value.imaging.QCtool>hwe_thr) 
        # some extra checks for chr X
        if (chrom=="X"){
          d2<-subset(d1,INFO.UKB>info_thr&HW_exact_p_value.imaging.QCtool>hwe_thr) 
          d2<-d2[,colnames(d),with=FALSE]
          w=which(colnames(d2)=="female")
          if(length(w)>0) {
            d2<-d2[,-w,with=FALSE]
          }
          rm(w)
        }
        cat('Number of SNPs in choromsome ',chrom,' after filtering:',NROW(d2),'\n')
        # variants filtered out because of hwe / INFO
        
        # save 1/d2 into all chromosomes
        if (!exists("d")) {d<-d2} else {d<-rbind(d,d2)}
        
        # clean
        rm(chrom,zz,header,snp_cols,pheno_cols)
        rm(d0,d1,d2)
      } ; rm(file,files)
      # select columns to keep, otherwise there are too many
      cols<-colnames(d)[grep(".QCtool|.UKB|.BGENIE|.HRC",colnames(d),invert=TRUE)]
      cols2<-colnames(d)[grep(".QCtool|.UKB",colnames(d),invert=FALSE)]
      cols2<-cols2[grep("hw_exact_p_value.|^info\\.|maf\\.",tolower(cols2))]
      #
      data<-d[,c(cols,cols2),with=FALSE]
      rm(cols,cols2)
      rm(snp_qc,d)
      # check SNPs per chr
      dim(data)
      table(data$chr)
      #-----------------------------------------------------------
      # rename phenotype columns
      colnames(data)<-gsub("^_","",gsub(pheno,"",gsub("residuals_","",colnames(data))))
      # get p from -logP or beta!
      data$P<-10^(-data$`-log10p`)
      data$P_z<-2*pnorm(-abs(data$beta/data$se))
      plot(data$P,data$P_z)
      p_col="P"
      # subset
      head(data)
      # subset to plot in Manhattan
      data2<-subset(data,P<0.01)
      data_sig<-subset(data2,P<5e-07)
      write.csv(data_sig,file=paste("clean/",gsub(" |\n","_",name),"_CHR",chr,"_Pmin7.csv",sep=""),row.names = FALSE)
      write.table(data,file=paste("clean/",gsub(" |\n","_",name),"_CHR",chr,"_",hwe_thr,"HWEp_",info_thr,"INFO_",maf_thr,"MAF.txt",sep=""),row.names = FALSE,col.names=TRUE,quote=FALSE)
      #---------------------------------------------------#
      # Plots
      #---------------------------------------------------#
      # For p-values, calculate chi-squared statistic
      
      chisq <- qchisq(1-data$P,1)
      mchisq<-summary(chisq)[c("Mean","Median","Max.")]
      
      # Calculate lambda gc (??gc)
      print("Calculate lambda")
      lambda <- round(median(chisq,na.rm=TRUE)/qchisq(0.5,1),digits=5)
      print(lambda)
      
      # save lambda values into file
      t<-cbind(name=name,
               lambda_all=lambda,
               mean_chisq_all=mchisq[1],median_chisq_all=mchisq[2],max_chisq_all=mchisq[3]
      )
      t<-as.data.frame(t)
      write.csv(t,file=paste("clean/",gsub(" |\n","_",name),"_CHR",chr,"_lambdas_chisqStats.csv",sep=""),row.names = FALSE)
      rm(t)
      
      pmin=ceiling(abs(log10(min(data2$P))))
      data2$chrom<-data2$chr
      data2$chrom[data2$chr=="X"]<-23
      data2$chrom<-as.numeric(data2$chrom)
      print("Plot\n")
      png(paste("plots/Manhattan_",gsub(" |\n","_",name),"_CHR",chr,".png",sep=""),width=1000)#, height=1000,
      # pdf(paste("plots/Manhattan_",gsub(" |\n","_",name),"_CHR",chr,".pdf",sep=""),width=15) #, height=1000,
      manhattan(data2,chr="chrom",bp="pos",p="P",snp="rsid",ylim=c(2,pmin),main=name
                # suggestiveline= (-log10(5e-8)),genomewideline= (-log10(1e-11))
      )
      
      dev.off()
      
      # pdf(paste("plots/QQplot_",gsub(" |\n","_",name),"_CHR",chr,".pdf",sep=""))
      
      png(paste("plots/QQplot_",gsub(" |\n","_",name),"_CHR",chr,".png",sep=""))
      qq(data$P,main=paste("Q-Q plot of GWAS p-values\n",name,sep=""))
      text(x=3,y=0.5,bquote(lambda==  .(lambda)))
      dev.off()
      print("Finished, hopefully without errors")
      # clean
      rm(data,data2)
      rm(data_sig)
      rm(pmin,lambda,mchisq,chisq)
      
    } else {
      print("Manhattan plot already exists, skip and continue")
    }
    
    # some cleaning...
    rm(name)
    rm(pattern1)
    #
    
    print(Sys.time())
    print("---------------------------------------------")
    #
  }
  rm(chr,pheno,pheno_root,type,info_thr,maf_thr,hwe_thr,p_col)
  gc()
  }
