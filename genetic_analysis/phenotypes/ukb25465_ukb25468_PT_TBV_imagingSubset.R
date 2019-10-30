library(GGally);library(ggplot2)
library(gridExtra);library(grid)
library("dplyr")
library(tidyr)
if (!require("gplots")) {
  install.packages("gplots", dependencies = TRUE)
  library(gplots)
}
if (!require("RColorBrewer")) {
  install.packages("RColorBrewer", dependencies = TRUE)
  library(RColorBrewer)
}
if (!require("corrplot")) {
  install.packages("corrplot", dependencies = TRUE)
  library(corrplot)
}
library(reshape2)

options(digits=7)
#------------------------------------------
# plots, some general parameters
#------------------------------------------
cols=c("#d95f02","#1f78b4","#33a02c","#000000") # '#1b9e77','#d95f02','#7570b3 # "#a6cee3"  ,"#7570b3"
mytheme<-theme_bw() + theme(panel.spacing = unit(0, "lines"), 
                            strip.background = element_rect(fill="white"), strip.text = element_text(size=16), 
                            # axis.title.x=element_blank(), axis.ticks.x =element_blank(), # axis.text.x =element_blank(),
                            # axis.text.x = element_text(size=16,colour="black",hjust=1,vjust=.5),
                            # axis.title.y=element_text(size=16),axis.text.y=element_text(size=16,colour="black"),
                            title=element_text(size=16),
                            axis.title=element_text(size=16),axis.text=element_text(size=16,colour="black"),
                            legend.text=element_text(size=16), legend.title =element_text(size=16) ) 
#------------------------------------------
# some color palettes
col1 <- colorRampPalette(c("#7F0000", "red", "#FF7F00", "yellow", "white",
                           "cyan", "#007FFF", "blue","#00007F"))
col2 <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                           "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                           "#4393C3", "#2166AC", "#053061"))
col3 <- colorRampPalette(c("red", "white", "blue"))
col4 <- colorRampPalette(c("#7F0000", "red", "#FF7F00", "yellow", "#7FFF7F",
                           "cyan", "#007FFF", "blue", "#00007F"))
wb <- c("white", "black")
#------------------------------------------

#-------------------------------------------
# Define directories, dependent on  system
if (Sys.info()['sysname']=='Windows') {dir="P://workspaces/"} else {dir="/data/workspaces/lag/workspaces/"}
primary_dir=paste(dir,"lg-ukbiobank/primary_data/",sep="")
working_dir=paste(dir,"lg-ukbiobank/working_data/amaia/pheno_files/",sep="")
h2_dir=paste(dir,"lg-ukbiobank/working_data/amaia/genetic_data/release_v2/cal/h2/output/PT/",sep="")
#
subset_name="imagingT1_N18057"
root="ukb25465_ukb25468"
pattern_run="_noBioCovs_noAssessmentC"
pattern_run2="_noBioCovs_noAssessmentC_TBV"
#
h2_dir=paste(dir,"lg-ukbiobank/working_data/amaia/genetic_data/release_v2/cal/h2/",subset_name,"/output/PT/",sep="")
#
lm_dir=paste(working_dir,"/genetic_v2/",root,"/lm/",sep="")
out_dir=paste(working_dir,"/genetic_v2/",root,"/summary_phenotypes/PT/",sep="")
out_plots_dir=paste(out_dir,"plots/",sep="")
#--------------------------------------------

#--------------------------------------------
# Read genetic correlation and h2 info
tBV_h2<-read.csv(paste(h2_dir,"gcta_estimates_TBV_",root,"_imaging",pattern_run,".csv",sep=""))
# tBV_rg<-read.csv(paste(h2_dir,"rG_PT_tBV_ukb21288_ukb21293_imaging_noBioCovs_noAssessmentC.csv",sep=""))
tBV_rg<-read.csv(paste(h2_dir,"rG_PT_TBV_",root,"_imaging_noBioCovs_noAssessmentC_TBV.csv",sep=""))
PT_gen<-read.csv(paste(h2_dir,"gcta_estimates_PT_",root,"_imaging",pattern_run,".csv",sep=""))
#
PT_rg<-subset(PT_gen,stat=="rG")
PT_rg$p1<-"left"; PT_rg$p2<-"right"
# combine rg
rg<-merge(PT_rg,tBV_rg,all=TRUE)
rg$run<-""
w<-c(grep("AdjTBV",rg$p1),grep("AdjTBV",rg$p2))
if (length(w)>0){
 rg$run[w]<-"adjTBV"
}
rm(w)
# match sample names
rg$sample<-tolower(rg$sample)
rg$sample[rg$sample=="all"]<-"total"
table(rg$sample)
# subset to only total or sex-specific results
rg_per_sex<-subset(rg,sample!="total")
rg<-subset(rg,sample=="total")

#
w<-which(duplicated(rg[,c("p1","p2","rG","se","pval")]))
if (length(w)>0){
rg[w,c("run","p1","p2","rG","se","pval")]
rg<-rg[-which(rg$p1=="left"&rg$p2=="right"&rg$run==""&rg$rG=="0.785719"),]
}
rm(w)

#--------------------------------------------

#--------------------------------------------
# Read phenotype data data
#--------------------------------------------
setwd(out_dir)
# read data, used for the genetic analyses
all_data<-read.table(paste("../",root,"_imaging",pattern_run,"VolumePhenotypes_residuals_wHeader.table",sep=""),header=TRUE,sep="\t")
pt2<-read.table(paste("../",root,"_imaging",pattern_run2,"_PTPhenotypes_residuals_wHeader.table",sep=""),header=TRUE,sep="\t") # adjusted
# flag columns with residuals to indicate that they were adjusted for TBV as well
colnames(pt2)<-gsub("residuals","residuals_adjTBV",colnames(pt2))
# select only residuals, because the number of digits from the different files does match, which creates mismatches when merging using all columns 
pt2<-pt2[,c(1:2,grep("res",colnames(pt2)))]

#
p="Planum_Temporale"
pt_cols<-grep(p,colnames(all_data))
n=grep("^totalBV",colnames(all_data))-1
pt<-all_data[c(1:n,pt_cols)]#;rm(pt_cols)
rm(all_data)
# combine pt and pt2
pt_all<-merge(pt,pt2,all=TRUE) # three columns added: residuals_adjTBV, per phenotype AI, L and R
rm(pt,pt2)


# tBV data
tBV_data<-read.table(paste("../",root,"_imaging",pattern_run,"_totalBVPhenotypes_residuals_wHeader.table",sep=""),header=TRUE,sep="\t")
plot(tBV_data$totalBV,tBV_data$Volume_of_brain_grey_white)
# combine
pt_tBV<-merge(pt_all,tBV_data[,c(1:2,grep("res|totalBV",colnames(tBV_data)))],stringsAsFactors=FALSE,all=TRUE)

# clean
# rm(pt_all,tBV_data)
#--------------------------------------------
# summary stats for phenotypes
pt_measures<-c("AI_VOLUME_Planum_Temporale","Volume_of_grey_matter_in_Planum_Temporale_left","Volume_of_grey_matter_in_Planum_Temporale_right")
tBV_measures<-c("totalBV") #,"Volume_of_grey_matter_normalised_for_head_size"
measures<-c(pt_measures,tBV_measures)
measures2<-c("AI","left","right","totalBV","totalBVadjHeadSize")
measures3<-c("AI","left","right","totalBV")
# define residuals
res_measures<-colnames(pt_tBV)[grep("residuals",colnames(pt_tBV))]
# compute summary per phenotype
summary_table<-as.data.frame(t(apply(pt_tBV[,c(measures,res_measures)],2,function(x) summary(x,digits=10))))
summary_table$SD<-apply(pt_tBV[,c(measures,res_measures)],2,function(x) sd(x,na.rm=TRUE) )
summary_table$N<-apply(pt_tBV[,c(measures,res_measures)],2,function(x) sum(!is.na(x)))
summary_table$measure<-rownames(summary_table)
write.csv(summary_table,file="summary_stats_PT_TBV_raw_residuals.csv",quote=FALSE,row.names=FALSE)
# summary stats for covariates
scanner_position<-colnames(pt_tBV)[grep("Scanner",colnames(pt_all))]
covs<-c("age","zage2","sex",paste("PC",1:10,sep=""),"array",scanner_position ) 
## quantitative covariates
covs_q<-c("age","zage2",paste("PC",1:10,sep=""),scanner_position )
summary_table_covs<-as.data.frame(t(apply(pt_all[,covs_q],2,function(x) summary(x)[1:6] ) ))
# summary_table_covs$`NA's`[which(summary_table_covs$`NA's`==summary_table_covs$Min.)]<-0
summary_table_covs$SD<-apply(pt_all[,covs_q],2,function(x) sd(x,na.rm=TRUE) )
summary_table_covs$N<-apply(pt_all[,covs_q],2,function(x) sum(!is.na(x)))
summary_table_covs$measure<-rownames(summary_table_covs)
write.csv(summary_table_covs,file="summary_stats_quant_covs.csv",quote=FALSE,row.names=FALSE)
## factors
covs_f<-c("sex","array")
summary_table_covsf<-as.data.frame(apply(pt_all[,covs_f],2,table))
write.csv(summary_table_covsf,file="summary_stats_bin_covs.csv",quote=FALSE,row.names=FALSE)
#--------------------------------------------
# Compute phenotypic correlations
phenos4cor<-pt_tBV[,measures];colnames(phenos4cor)<-measures3
rphenos4cor<-pt_tBV[,paste("residuals_",measures,sep="")];colnames(rphenos4cor)<-paste("residuals_",measures3,sep="")
r2phenos4cor<-pt_tBV[,gsub("adjTBV_totalBV","totalBV",paste("residuals_adjTBV_",measures,sep=""))];colnames(r2phenos4cor)<-gsub("adjTBV_totalBV","totalBV",paste("residuals_adjTBV_",measures3,sep=""))
# check if TBV is still a significant predictor of the residuals adjTBV 
x<-lm(pt_tBV$residuals_adjTBV_Volume_of_grey_matter_in_Planum_Temporale_left~pt_tBV$totalBV)
summary(x)
rm(x)
# correlation matrix
pcor<-cor(phenos4cor, use= "pairwise.complete.obs")
rpcor<-cor(rphenos4cor, use= "pairwise.complete.obs")
r2pcor<-cor(r2phenos4cor, use= "pairwise.complete.obs")
# create matrix for significance
pcor_sig<-cor.mtest(phenos4cor)$p
colnames(pcor_sig)<-row.names(pcor_sig)<-colnames(pcor)
#
rpcor_sig<-cor.mtest(rphenos4cor)$p
colnames(rpcor_sig)<-row.names(rpcor_sig)<-colnames(rpcor)
#
r2pcor_sig<-cor.mtest(r2phenos4cor)$p
colnames(r2pcor_sig)<-row.names(r2pcor_sig)<-colnames(r2pcor)

# some plots
png(file="./plots/PT_tBV_ggpairs.png")
ggpairs(phenos4cor[,measures3])
dev.off()


png(file="./plots/residuals_PT_tBV_ggpairs.png")
ggpairs(rphenos4cor)
dev.off()


png(file="./plots/residuals_adjTBV_PT_tBV_ggpairs.png")
ggpairs(r2phenos4cor)
dev.off()

#--------------------------------------------
# convert to long format, to save later
#--------------------------------------------
## raw measures, residuals and residuals adjTBV
for (t in c("pcor","rpcor","r2pcor")){
  tmp<-get(t)
  stmp<-get(paste(t,"_sig",sep=""))
  ## correlation 
  pcor_df<-as.data.frame(tmp)
  pcor_df$region1<-rownames(pcor_df)
  rownames(pcor_df)<-1:NROW(pcor_df)
  pcor_df<-pcor_df %>% gather(key="region2",value="r",-region1)
  ## significance
  pcor_sig_df<-data.frame(stmp)
  pcor_sig_df$region1<-rownames(pcor_sig_df)
  rownames(pcor_sig_df)<-1:NROW(pcor_sig_df)
  pcor_sig_df<-pcor_sig_df%>% gather(key="region2",value="pval",-region1)
  colnames(pcor_sig_df)<-c("region1","region2","pval")
  ## combine
  out<-merge(pcor_df,pcor_sig_df,by=c("region1","region2"))
  ## sort region1 and region2, to remove duplicated items
  sort_df2<- t(apply(out[,1:2], 1, sort))
  out<-out[!duplicated(sort_df2),]
  out<-out[-which(out$region1==out$region2),]
  ##
  assign(paste("df_",t,sep=""),out)
  ## clean intermediate objects
  rm(pcor_df,pcor_sig_df,sort_df2)
  rm(tmp,stmp,out)
}

# combine all three
all_cor_sig<-rbind(df_pcor,df_rpcor,df_r2pcor)
rm(df_pcor,df_rpcor,df_r2pcor)
#--------------------------------------------

#--------------------------------------------
# Phenotypes
summary3m<-pcor
summary3pmat<-pcor_sig
# duplicate, for now only gen. corr.; and significance matrix
summary3m[is.na(summary3m)]<-0
summary3pmat[is.na(summary3pmat)]<-1

# rG
# make duplicate to have complete matrix
rg2<-rg[,c("p2","p1","rG","pval","run")]; colnames(rg2)<-c("p1","p2","rG","pval","run")
rg3_all<-rbind(rg[,c("p2","p1","rG","pval","run")],rg2[,c("p2","p1","rG","pval","run")])
#
rg3<-subset(rg3_all,run=="")
rm(rg2)
# convert rG to matrix format
rg_df<-acast(rg3[,c("p2","p1","rG")], p1~p2, value.var="rG",fun.aggregate = mean)
rg_sig_df<-acast(rg3[,c("p1","p2","pval")], p1~p2, value.var="pval",fun.aggregate = mean)

# fill in lower triangle with genetic correlations
summary3m[lower.tri(summary3m)]<- NA
summary3m[lower.tri(summary3m)]<- rg_df[lower.tri(rg_df)]
summary3pmat[lower.tri(summary3pmat)]<- NA
summary3pmat[lower.tri(summary3pmat)]<- rg_sig_df[lower.tri(rg_sig_df)]

pdf(file=paste(h2_dir,"corplot_pheno_geno_PT_tBV.pdf",sep=""),height=5,width=15)
corrplot(summary3m, 
         method="color",
         tl.col="black", tl.srt=45, #Text label color and rotation
         # Combine with significance
         p.mat=summary3pmat, sig.level = 0.05, insig = "blank",
         addCoef.col = "black", # Add coefficient of correlation
         # hide correlation coefficient on the principal diagonal
         diag=FALSE,
         mar=c(0,0,2,1) ,
         # type="lower",
         # order="hclust",
         tl.cex=1,number.cex=0.75,
         # xlim=c(-0.3,0.3),ylim=c(-0.3,0.3),
         # cl.lim=c(-1,1),col= col3(100),
         title="Phenotypic (raw) and genotypic correlations\n PT volumes and global brain volume measures"
         # title=""
         )
dev.off()

# adjusting for TBV
rg3<-subset(rg3_all,run=="adjTBV")

# convert rG to matrix format
rg_df<-acast(rg3[,c("p2","p1","rG")], p1~p2, value.var="rG",fun.aggregate = mean)
rg_sig_df<-acast(rg3[,c("p1","p2","pval")], p1~p2, value.var="pval",fun.aggregate = mean)
rg_df[upper.tri(rg_df)]<- NA
rg_sig_df[upper.tri(rg_sig_df)]<- NA

summary3m_adj<-rbind(summary3m,rg_df)
summary3pmat_adj<-rbind(summary3pmat,rg_sig_df)

# # fill in lower triangle with genetic correlations
# summary3m[lower.tri(summary3m)]<- NA
# summary3m[lower.tri(summary3m)]<- rg_df[lower.tri(rg_df)]
# summary3pmat[lower.tri(summary3pmat)]<- NA
# summary3pmat[lower.tri(summary3pmat)]<- rg_sig_df[lower.tri(rg_sig_df)]

rownames(summary3m_adj)<-gsub("TBV","totalBV",gsub("right","R",gsub("left","L",gsub("Adj","_adj",rownames(summary3m_adj)))))
colnames(summary3m_adj)<-gsub("TBV","totalBV",gsub("right","R",gsub("left","L",gsub("Adj","_adj",colnames(summary3m_adj)))))

pdf(file=paste(h2_dir,"corplot_pheno_geno_PTadjTBV_tBV.pdf",sep=""),height=5,width=15)
corrplot(summary3m_adj, 
         method="color",
         tl.col="black", tl.srt=45, #Text label color and rotation
         # Combine with significance
         p.mat=summary3pmat_adj, sig.level = 0.05, insig = "blank",
         addCoef.col = "black", # Add coefficient of correlation
         # hide correlation coefficient on the principal diagonal
         diag=FALSE,
         na.label="-",
         mar=c(0,0,2,1) ,
         # type="lower",
         # order="hclust",
         tl.cex=1,number.cex=0.75,
         # xlim=c(-0.3,0.3),ylim=c(-0.3,0.3),
         # cl.lim=c(-1,1),col= col3(100),
         title="Phenotypic (raw) and genotypic correlations (w/wo totalBV adjustment)\n PT volumes and global brain volume measures"
         # title=""
)
dev.off()

# use residuals for phenotypic correlations
# Phenotypes
## replace upper triangles for correlations of residualized phenotypes
## with and without totalBV adjustment
summary3m_r<-summary3m
summary3m_r[upper.tri(summary3m_r)]<- rpcor[upper.tri(rpcor)]
summary3m_r_sig<-summary3pmat
summary3m_r_sig[upper.tri(summary3m_r_sig)]<- rpcor_sig[upper.tri(rpcor_sig)]
# adjTBV
summary3m_r2<-r2pcor
summary3m_r2[lower.tri(summary3m_r2)]<- rg_df[lower.tri(rg_df)]
summary3m_r2_sig<-r2pcor_sig
summary3m_r2_sig[lower.tri(summary3m_r2_sig)]<- rg_sig_df[lower.tri(rg_sig_df)]

rownames(summary3m_r2)<-colnames(summary3m_r2)<-gsub("residuals_","",colnames(summary3m_r2))

# add adjusted
summary3m_adj2<-rbind(summary3m_r,summary3m_r2)
summary3pmat_adj2<-rbind(summary3m_r_sig,summary3m_r2_sig)

# clean
# rm(summary3m_r,summary3m_r2)
# rm(summary3m_r_sig,summary3m_r2_sig)

#
pdf(file=paste(h2_dir,"corplot_respheno_geno_PTadjTBV_tBV.pdf",sep=""),height=5,width=15)
corrplot(summary3m_r, 
         type="lower",
         method="color",
         tl.col="black", tl.srt=45, #Text label color and rotation
         # Combine with significance
         p.mat=summary3m_r_sig, sig.level = 0.05, insig = "blank",
         addCoef.col = "black", # Add coefficient of correlation
         # hide correlation coefficient on the principal diagonal
         diag=FALSE,
         na.label="-",
         mar=c(0,0,2,1) ,
         # type="lower",
         # order="hclust",
         tl.cex=1,number.cex=0.75,
         # xlim=c(-0.3,0.3),ylim=c(-0.3,0.3),
         # cl.lim=c(-1,1),col= col3(100),
         title="Phenotypic (residuals) and genotypic correlations (without totalBV adjustment)\n PT volumes and global brain volume measures"
         # title=""
)

corrplot(summary3m_r2, 
         method="color",
         tl.col="black", tl.srt=45, #Text label color and rotation
         # Combine with significance
         p.mat=summary3m_r2_sig, sig.level = 0.05, insig = "blank",
         addCoef.col = "black", # Add coefficient of correlation
         # hide correlation coefficient on the principal diagonal
         diag=FALSE,
         na.label="-",
         mar=c(0,0,2,1) ,
         # type="lower",
         # order="hclust",
         tl.cex=1,number.cex=0.75,
         # xlim=c(-0.3,0.3),ylim=c(-0.3,0.3),
         # cl.lim=c(-1,1),col= col3(100),
         title="Phenotypic (residuals) and genotypic correlations (with totalBV adjustment)\n PT volumes and global brain volume measures"
         # title=""
)


dev.off()




#------------------------------------------



# combine rg and pheno r
rg$p1<-gsub("AdjTBV","",rg$p1)
rg$p2<-gsub("AdjTBV","",rg$p2)
rg_wide<-rg[,c("run","p1","p2","rG","se","pval")] %>% 
  gather(variable,value,-(run:p2)) %>% 
  unite(temp, run, variable) %>% spread(temp, value)
colnames(rg_wide)<-gsub("^_","",colnames(rg_wide))

summary_all<-merge(rg_wide,all_cor_sig[,c("region1","region2","r","pval")], by.x=c("p1","p2"),by.y=c("region1","region2"),all=TRUE,suffixes=c(".rG",".rPheno"))
# save files as tables as well
write.csv(summary_all,file=paste(h2_dir,"rG_phenoCor_PTadjTBV_TBV_summary.csv",sep=""),row.names=FALSE)
#------------------------------------------

