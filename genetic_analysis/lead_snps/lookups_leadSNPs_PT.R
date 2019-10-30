# After: lookups_leadSNPs.sh
# Get stats for lead SNPs and SNPs in LD (r2>0.6) ## lookups_leadSNPs.sh

library(dplyr);library(tidyverse)
options(stringsAsFactors = FALSE)
#----------------------------------------------------------------------
# define working_dirs
# define working_dir
if (Sys.info()['sysname']=='Windows') {
  # dir="P://workspaces/"
  dir="\\\\data/lag/workspaces/"
} else {dir="/data/workspaces/lag/workspaces/"}
#----------------------------------------------------------------------
subset_name="imagingT1_N18057"
ld_dir=paste(dir,"lg-ukbiobank/working_data/amaia/genetic_data/PT/",subset_name,"/lead_snps/LD/",sep="")
working_dir=paste(dir,"lg-ukbiobank/working_data/amaia/genetic_data/PT/",subset_name,"/lead_snps/lookups/PT_checks/",sep="")
working_dir2=paste(dir,"lg-ukbiobank/working_data/amaia/genetic_data/PT/",subset_name,"/lead_snps/lookups/GWAS_checks/",sep="")
#----------------------------------------------------------------------
# read ld files
ld_files<-paste(ld_dir,list.files(ld_dir,pattern="r2.*_rs*.*.ld"),sep="")
for (f in ld_files){
  t<-read.table(f,header=T,strip.white = T)
  if (f==ld_files[1]){ld<-t} else {ld<-rbind(ld,t)}
  rm(t)
}
rm(f,ld_files)
ld<-subset(ld,R2>=0.6)
#----------------------------------------------------------------------
# read stats for these snps
files<-paste(working_dir,list.files(working_dir,pattern="lookup.*.txt"),sep="")
for (f in files){
  t<-read.table(f,sep=" ",header=T)
  if (NROW(t)>0){
    t$f<-f
    if (length(grep(":",t$chr)) >0 ) {
     file_i<-unique(unlist(lapply(strsplit(t$chr,":"),length)))
     if ( length(file_i)==1 & file_i==2 ){
      t$assoc_file<-sapply(strsplit(t$chr,":"),"[[",1)
      t$chr<-sapply(strsplit(t$chr,":"),"[[",2)
     }
    }
    t$leadSNP<-sapply(strsplit(sapply(strsplit(f,"lookup_"),"[[",2),"_"),"[[",1)
    # t<-read.table(f,sep=" ",nrow=24,header=T)
    if (f==files[1]){d<-t} else {d<-merge(d,t,all=TRUE)}
  }
  rm(t)
}
rm(f,files)
# extract relevant info
table(d$chr)
d$file<-sapply(strsplit(d$f,"PT_checks/"),"[[",2)
d$phenotype<-gsub(".txt","",sapply(strsplit(d$file,"_r2_0.6_"),"[[",2))
d$phenotype<-gsub("VOLUME_Planum_Temporale_|_VOLUME_Planum_Temporale|Volume_of_grey_matter_in_Planum_Temporale_","",d$phenotype)
# # fix: there were some AI flagged as totalBV
# subset(d,assoc_file=="AI_VOLUME_Planum_Temporale_totalBV_CHRall_1e-07HWEp_0.7INFO_0.001MAF.txt.gz")$phenotype
# d<-subset(d,assoc_file!="AI_VOLUME_Planum_Temporale_totalBV_CHRall_1e-07HWEp_0.7INFO_0.001MAF.txt.gz"|is.na(assoc_file))

# define sample
d$sample<-"total"
d$sample[grep("_females",d$phenotype)]<-"females"
d$sample[grep("_males",d$phenotype)]<-"males"

d$region<-NA
d$region[grep("PTHG",d$file)]<-"PTHG"
d$region[grep("Planum_Temporale",d$file)]<-"PT"
d$region[grep("lookup_rs7420166_r2_0.6_totalBV.txt",d$file)]<-"TBV"
table(d$region)

d$measure<-"NotLat"
d$measure[grep("AI_|AI$",d$phenotype)]<-"AI"
d$measure[grep("AIdiff_|AIdiff$",d$phenotype)]<-"(L-R)"
d$measure[grep("_left|left",d$phenotype)]<-"L"
d$measure[grep("_right|right",d$phenotype)]<-"R"
d$measure[d$phenotype=="totalBV"]<-"TBV"
table(d$measure)

# add cols for adj
d$adj<-"-"
d$adj[grep("adjTBV|_TBV|_totalBV",d$f)]<-"TBV"
d$adj[grep("adjLplusR",d$f)]<-"(L+R)"
table(d$adj)
d$adj[d$measure=="TBV"]<-"-" # because TBV was not corrected for TBv
d$adj<-factor(d$adj,levels=c("-","TBV","(L+R)"))

# merge with ld
d_ld<-merge(d,ld,by.x=c("rsid","chr","pos"),by.y=c("SNP_B","CHR_B","BP_B"),all.x=TRUE)
# define order of cols to keep
cols<-c(colnames(d_ld)[1:6],"leadSNP","INFO.UKB","HW_exact_p_value.imaging.QCtool","SNP_A","R2","region","measure","adj","sample","beta","se","t","P")

# remove duplicated rows
w<-which(duplicated(d_ld[,cols]))
if (length(w)>0){d_ld<-d_ld[-w,]}
rm(w)

head(d_ld)

# convert to wide format
d_ld_wide<-d_ld[,cols] %>% gather(variable,value,-(rsid:sample)) %>%
  unite(temp, adj, variable) %>% spread(temp, value)
rm(cols)
colnames(d_ld_wide)<-gsub("-_","",colnames(d_ld_wide))
# order columns of wide
stats=c("t","P") #"beta","se",
w_cols<-c(colnames(d_ld_wide)[1:14],
          paste(rep(c("","TBV_","(L+R)_"),each=length(stats)),stats,sep="")
          )

# sort by chr, pos
d_ld_wide<-d_ld_wide[,c(w_cols)] %>% arrange(region,chr,pos,rsid,sample,measure)
# clean, some SNPs got in because of partial grep matches...
subset(d_ld_wide,is.na(R2)&P>0.1&measure=="AI")
d_ld_wide<-subset(d_ld_wide,!(is.na(R2)))

# save
write.csv(d_ld_wide,file=paste(working_dir,"leadSNPs_r20.6_PT_stats.csv"),quote=FALSE,row.names = FALSE)
# rm(d_ld_wide,d_ld,d)
#----------------------------------------------------------------------
# read lookups from public/ other GWASes
#----------------------------------------------------------------------
files<-paste(working_dir2,list.files(working_dir2,pattern="lookups_*.*.txt"),sep="")
files<-files[grep("EA|IQ|ASD|ADHD|SCZ",files)]
for (f in files){
  t<-read.table(f,sep="\t",header=T)
  w<-grep("chr",tolower(colnames(t)))
  w2<-which(t[,w]==colnames(t)[w])
  p<-gsub("lookups_|.txt","",gsub(working_dir2,"",f))
  t2<-t[-w2,]
  t2$PHENO<-p
  if (length(grep(":",t2[,w]))==0) {t2$FILE<-p
  } else {
    t2$FILE<-sapply(strsplit(t2[,1],":"),"[[",1)
    t2$CHR<-gsub("chr","",sapply(strsplit(t2[,w],":"),"[[",2))
    
  }
  assign(p,t2)
  # clean
  rm(p,t,t2)
  rm(w,w2)
}
rm(f,files)

# some extra tweaks
ADHD$CHR<-sapply(strsplit(ADHD$CHR,":"),"[[",1)
SCZ$SNP<-SCZ$snpid
colnames(SCZ)<-toupper(colnames(SCZ))
EA$P<-EA$Pval
EA$BP<-EA$POS
EA$SNP<-sapply(strsplit(EA$MarkerName,":"),"[[",2)
EA$FILE<-sapply(strsplit(EA$MarkerName,":"),"[[",1)
IQ$FILE<-sapply(strsplit(IQ$SNP,":"),"[[",1)
IQ$SNP<-sapply(strsplit(IQ$SNP,":"),"[[",2)
# combine
all<-merge(ASD,SCZ,all=TRUE,stringsAsFactors=FALSE)
all<-merge(all,ADHD,all=TRUE,stringsAsFactors=FALSE)
all<-merge(all,EA,all=TRUE,stringsAsFactors=FALSE)
all<-merge(all,IQ,all=TRUE,stringsAsFactors=FALSE)
#
cols<-c("SNP","CHR","POS","PHENO","FILE","P","INFO","minINFO","OR","Beta","stdBeta","SE","Zscore")
all_s<-all[,cols]
# combine with LD info
all_s_ld<-merge(all_s,ld,by.x=c("CHR","SNP"),by.y=c("CHR_B","SNP_B"),all=TRUE)
all_s_ld<-all_s_ld[,c(colnames(all_s),"SNP_A","R2")]
write.csv(all_s_ld,paste(working_dir2,"lookup_leadSNPs_GWASes.csv",sep=""),row.names = FALSE)

#----------------------------------------------------------------------
f<-paste(working_dir2,list.files(working_dir2,pattern="lookups_*.*Guadalupe2015_SNPs.txt"),sep="")
t<-read.table(f,sep=" ",header=F,strip.white = T)
l<-strsplit(t$V1,"\t")
l2<-lapply(l,function(x){x[-which(unlist(x)=="")]})
t2<-data.frame(do.call("rbind",l2))
#BETA							SE						L95						U95									STAT												P	
colnames(t2)<-c("file","CHR","SNP","POS","A1","TEST","NMISS","BETA","SE","L95","U95","STAT","P")
# combine with LD info
g<-merge(t2,ld,by.x=c("CHR","POS"),by.y=c("CHR_B","BP_B"),all=TRUE)
write.csv(g,paste(working_dir2,"lookup_leadSNPs_Guadalupe2015.csv",sep=""),row.names = FALSE)
