library(ggpubr)
library(ggplot2)
library(reshape2)
library(grid)
library(gridExtra)
library(tidyr)
#----------------------------------------------------------------------
#------------------------------------------
# plots, some general parameters
#------------------------------------------
cols=c("#d95f02","#1f78b4","#33a02c","#000000") # '#1b9e77','#d95f02','#7570b3 # "#a6cee3"  ,"#7570b3"
mytheme<-theme_bw() + theme(panel.spacing = unit(0, "lines"), 
                            strip.background = element_rect(fill="white"), strip.text = element_text(size=16), 
                            axis.title.x=element_blank(), axis.ticks.x =element_blank(), # axis.text.x =element_blank(),
                            # axis.text.x = element_text(angle = 90, size=16,colour="black",hjust=1,vjust=.5),
                            axis.text.x = element_text(size=16,colour="black",hjust=1),
                            axis.title.y=element_text(size=18),axis.text.y=element_text(size=16,colour="black"),legend.text=element_text(size=16)) 

meas_labels<-unlist(c(bquote(''*h^2*'(AI)'),bquote(''*h^2*'(L)'),bquote(''*h^2*'(R)'),bquote(''*rho*'(L,R)')))
#----------------------------------------------------------------------
options(stringsAsFactors = FALSE)
subset_name="imagingT1_N18057"
# define pattern for this run
root="ukb25465_ukb25468"
pattern_run= "noBioCovs_noAssessmentC"
pattern_run2=paste("_",root,"_imaging_",pattern_run,sep="")
pattern_run3=paste("_",root,"_imaging_",pattern_run,"_TBV",sep="")

# define working_dir
if (Sys.info()['sysname']=='Windows') {dir="P://workspaces/"} else {dir="/data/workspaces/lag/workspaces/"}
reml_dir=paste(dir,"lg-ukbiobank/working_data/amaia/genetic_data/release_v2/cal/h2/",subset_name,"/reml/",sep="")
out_dir=paste(dir,"lg-ukbiobank/working_data/amaia/genetic_data/release_v2/cal/h2/",subset_name,"/output/PT/",sep="")
#
dir.create(file.path(out_dir), showWarnings = FALSE)
# out_dir="K://written/tex/multilateral/presentations/figures/"
setwd(out_dir)
#----------------------------------------------------------------------
tbv_h2<-read.table(paste(reml_dir,"hsq_summary_TBV_v2cal",pattern_run3,".table",sep=""),stringsAsFactors = FALSE)
colnames(tbv_h2)<-c("file","n","h2","se","pval")
tbv_h2$file<-gsub(":n","",tbv_h2$file)
tbv_h2$phenotype<-gsub(".hsq","",sapply(strsplit(tbv_h2$file,"_reml_"),"[[",2)) #_residuals_
tbv_h2$phenotype<-gsub("_TBV|residuals_","",gsub(pattern_run,"",tbv_h2$phenotype))
tbv_h2$measure<-"NotLat"
tbv_h2$stat<-"h2"
tbv_h2$sample<-"Total"
tbv_h2$region<-gsub("_PT|_Volume|VOLUME_|Volume_of_|grey_matter_in_|ukb_cal_snpQC_imagingT1_clean_rm025_adj_|reml_|residuals_|.hsq|:n","",tbv_h2$phenotype)
tbv_h2$measure<-"NotLat"
tbv_h2$measure[grep("AI",tbv_h2$region)]<-"AI"
tbv_h2$measure[grep("_left",tbv_h2$region)]<-"L"
tbv_h2$measure[grep("_right",tbv_h2$region)]<-"R"
table(tbv_h2$measure)
tbv_h2$region<-gsub("_"," ",gsub("_left|_right|AI_","",tbv_h2$region))
tbv_h2$lobe[grep("Temporal",tbv_h2$region)]<-"Temporal"
tbv_h2$stat<-"h2"
tbv_h2$estimate<-tbv_h2$h2
# save 
write.csv(tbv_h2,  file=paste(out_dir,"gcta_estimates_TBV",pattern_run2,".csv",sep=""),row.names = FALSE,quote=TRUE)
#-------------------------
summary_rG1<-read.table(paste(reml_dir,"hsq_summary_PT_TBV_rG_pval_v2cal",pattern_run2,".table",sep=""),stringsAsFactors = FALSE) # will need to add significance!
summary_rG2<-read.table(paste(reml_dir,"hsq_summary_PT_TBV_rG_pval_v2cal",pattern_run3,".table",sep=""),stringsAsFactors = FALSE) # will need to add significance!
summary_rG<-rbind(summary_rG1,summary_rG2);rm(summary_rG1,summary_rG2)
colnames(summary_rG)<-c("file","rG","se","pval")
summary_rG$run<-pattern_run
summary_rG$run[grep(paste(pattern_run,"_TBV",sep=""),summary_rG$file)]<-paste(pattern_run,"_TBV",sep="")
#
summary_rG$file2<-gsub("^_|^_TBV_|_bivar_diff0.hsq:rG","",sapply(strsplit(summary_rG$file,pattern_run),"[[",2))
# get both phenotypes from file name...
summary_rG$phenos<-gsub("^_|_$","",gsub("Volume_of_grey_matter_in_Planum_Temporale|_VOLUME_Planum_Temporale","",summary_rG$file2))
summary_rG$phenos<-gsub("Volume_of_brain_grey_white","totalBV",summary_rG$phenos)
summary_rG$phenos<-gsub("_normalised_for_head_size","adjHeadSize",summary_rG$phenos)
summary_rG$phenos<-gsub("__","_",summary_rG$phenos)
summary_rG$p1<-sapply(strsplit(summary_rG$phenos,"_"),"[[",1)
summary_rG$p2<-sapply(strsplit(summary_rG$phenos,"_"),"[[",2)
# define sample: all, males or females
summary_rG$sample<-"total"
summary_rG$sample[grep("_males",summary_rG$phenos)]<-"males"
summary_rG$sample[grep("_females",summary_rG$phenos)]<-"females"
# order pheno cols, so that we can remove duplicates
summary_rG[,c("p1","p2")]<-t(apply(summary_rG[,c("p1","p2")], 1, sort))
# delete duplicated
w<-which(duplicated(summary_rG[,c("rG","se","p1","p2","pval","sample","run")]))
summary_rG[w,]
if (length(w)>0){summary_rG<-summary_rG[-w,]}

# convert to wide format, for summary table
# wide format
summary_rG_wide<-summary_rG[,c("run","sample","p1","p2","rG","se","pval")] %>% 
  gather(variable,value,-(run:p2)) %>% 
  unite(temp, run, variable) %>% spread(temp, value)
colnames(summary_rG_wide)<-gsub("noBioCovs_noAssessmentC_","",colnames(summary_rG_wide))

# if present in TBV run, adjust name of L, R and AI
w<-which(summary_rG$run==paste(pattern_run,"_TBV",sep=""))
summary_rG$p1[w]<-gsub("$","AdjTBV",summary_rG$p1[w])
summary_rG$p2[w]<-gsub("$","AdjTBV",summary_rG$p2[w])
summary_rG$p1[w]<-gsub("totalBVAdjTBV","totalBV",summary_rG$p1[w])
summary_rG$p2[w]<-gsub("totalBVAdjTBV","totalBV",summary_rG$p2[w])
rm(w)
#
summary_rG$measure<-"rG"
summary_rG$stat<-"rG"

#-------------------------
# save 
write.csv(summary_rG,  file=paste(out_dir,"rG_PT_TBV",pattern_run3,".csv",sep=""),row.names = FALSE,quote=TRUE)
write.csv(summary_rG_wide,  file=paste(out_dir,"rG_PT_TBV",pattern_run2,"_wide.csv",sep=""),row.names = FALSE,quote=TRUE)
#-------------------------


#----------------------------------------------------------------------
# read results after adjusting for TBV as well
h2_pt_tbv<-read.csv(paste(out_dir,"gcta_estimates","_TBV",pattern_run2,".csv",sep=""))
h2_pt_tbv$sample<-"Total"
h2_pt_tbv$estimate<-h2_pt_tbv$h2
h2_pt_tbv$region<-gsub("_Volume|VOLUME_|Volume_of_|grey_matter_in_|ukb_cal_snpQC_imagingT1_clean_rm025_adj_|reml_|residuals_|.hsq|:n","",h2_pt_tbv$phenotype)
h2_pt_tbv$measure<-"NotLat"
h2_pt_tbv$measure[grep("AI",h2_pt_tbv$region)]<-"AI"
h2_pt_tbv$measure[grep("_left",h2_pt_tbv$region)]<-"L"
h2_pt_tbv$measure[grep("_right",h2_pt_tbv$region)]<-"R"
rg_pt_tbv<-read.csv(paste(out_dir,"rG_PT_TBV",pattern_run2,"_TBV.csv",sep=""))

pt_tbv<-merge(tbv_h2,summary_rG,all=TRUE,stringsAsFactors=FALSE)
pt_tbv$adj<-"TBV"
pt_tbv$p1<-gsub("AdjTBV","",pt_tbv$p1)
pt_tbv$p2<-gsub("AdjTBV","",pt_tbv$p2)

# read output from PT results (without adjusting for TBV)
pt<-read.csv(paste(out_dir,"gcta_estimates_PT",pattern_run2,".csv",sep=""))
pt1<-merge(pt,pt_tbv,all=T)
pt1$adj[is.na(pt1$adj)]<-""

# plot heritabilities from h2
gcta_h2_barplot_tbv<-ggplot(data=subset(pt1,stat=="h2"&sample=="Total"&measure!="NotLat"),color="black",aes(x=measure,y=estimate,width=0.7,fill=measure,alpha=as.factor(adj))) + 
  geom_bar(aes(y=estimate),position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=estimate-se, ymax=estimate+se,width=.5), position=position_dodge(.7)) +
  mytheme + 
  scale_fill_manual(values=cols,name = "") + coord_cartesian(ylim = c(0,1)) +
  scale_alpha_discrete( range=c(1,0.5), na.value = 0,name="") +
  ylab(bquote('Estimate ('*h^2*')')) +
  theme(axis.text.x = element_text(angle = 0)) # + theme(legend.position="none") 


png(paste(out_dir,"gcta_estimates_PT_",pattern_run3,".png",sep=""))
gcta_h2_barplot_tbv
dev.off()

