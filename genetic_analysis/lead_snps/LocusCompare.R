#devtools::install_github("boxiangliu/locuscomparer")
library(locuscomparer)

subset="imagingT1_N18057"
#
if (Sys.info()['sysname']=='Windows') {dir="P://"} else {dir="/data/workspaces/lag/"}
working_dir=paste(dir,"workspaces/lg-ukbiobank/working_data/amaia/genetic_data/PT/",subset,"/locuscompare/",sep="")
setwd(working_dir)


chr2_fn<-"AI_VOLUME_Planum_Temporale_CHRall_1e-07HWEp_0.7INFO_0.001MAF_chr2_rs7420166.input4locuscompare"
chr10_fn<-"AI_VOLUME_Planum_Temporale_CHRall_1e-07HWEp_0.7INFO_0.001MAF_chr10_rs41298373.input4locuscompare"

locuscompare(in_fn1 = chr2_fn, in_fn2 = eqtl_fn, title1 = 'CAD GWAS', title2 = 'Coronary Artery eQTL')


d<-read.table(eqtl_fn,header=TRUE)
d$chr<-"10"
make_locuszoom(d)
