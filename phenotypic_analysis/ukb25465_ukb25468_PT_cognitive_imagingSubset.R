library(GGally);library(ggplot2)
library(gridExtra);library(grid)
library("dplyr")
library(cowplot) 
library(corrplot)
library(PerformanceAnalytics)
library(plyr)
library(Hmisc)
options(stringsAsFactors = FALSE)


#------------------------------------------
# functions
#------------------------------------------
# to adjust a column by another
adj_by_col<-function(x,y,data){
  t<-data.frame(cbind(x=data[,x],y=data[,y]))
  lmt<-lm(x~y,data=t)
  # which are NA
  t$x_adjy<-NA
  w<- which(is.na(t$x)|is.na(t$y))
  if (length(w)>0){
    t$x_adjy[-w]<-lmt$residuals
  } else {
    t$x_adjy<-lmt$residuals
  }
  rm(w)
  return(t$x_adjy)
}


# ++++++++++++++++++++++++++++
# flattenCorrMatrix
# ++++++++++++++++++++++++++++
# from: http://www.sthda.com/english/wiki/correlation-matrix-a-quick-start-guide-to-analyze-format-and-visualize-a-correlation-matrix-using-r-software
# cormat : matrix of the correlation coefficients
# pmat : matrix of the correlation p-values
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}

#------------------------------------------
# plots, some general parameters
#------------------------------------------
cols=c("#d95f02","#1f78b4","#33a02c","#000000") # '#1b9e77','#d95f02','#7570b3 # "#a6cee3"  ,"#7570b3"
mytheme<-theme_bw()+ theme(panel.spacing = unit(0, "lines"), 
                            strip.background = element_rect(fill="white"), strip.text = element_text(size=16), 
                            # axis.title.x=element_blank(), axis.ticks.x =element_blank(), # axis.text.x =element_blank(),
                            # axis.text.x = element_text(size=16,colour="black",hjust=1,vjust=.5),
                            # axis.title.y=element_text(size=16),axis.text.y=element_text(size=16,colour="black"),
                            title=element_text(size=16),
                            axis.title=element_text(size=16),axis.text=element_text(size=16,colour="black"),
                            legend.text=element_text(size=16), legend.title =element_text(size=16) ) 
#--------------------------------------------
root="ukb25465_ukb25468"
pattern_run="_noBioCovs_noAssessmentC"
#
# Define directories, dependent on  system
if (Sys.info()['sysname']=='Windows') {
  # dir="P://workspaces/"
  dir="\\\\data/lag/workspaces/"
} else {dir="/data/workspaces/lag/workspaces/"}
primary_dir=paste(dir,"lg-ukbiobank/primary_data/",sep="")
working_dir=paste(dir,"lg-ukbiobank/working_data/amaia/pheno_files/",sep="")
out_dir=paste(working_dir,"/genetic_v2/",root,"/summary_phenotypes/PT/",sep="")
h2_dir=paste(dir,"lg-ukbiobank/working_data/amaia/pheno_files/genetic_data/release_v2/cal/h2/cognitive/",sep="")
#--------------------------------------------
# Read phenotype data data
#--------------------------------------------
setwd(out_dir)
# read data
vol_data<-read.table(paste("../",root,"_imaging",pattern_run,"VolumePhenotypes_residuals_wHeader.table",sep=""),header=TRUE,sep="\t")
p="Planum_Temporale|Heschls"
pt_cols<-grep(p,colnames(vol_data))
n_col=grep("^Volumetric_scaling_from_T1_head_image_to_standard_space",colnames(vol_data))
pt_data<-vol_data[c(1:n_col,pt_cols)];rm(pt_cols)
rm(vol_data,n_col,p)
#------------------------------------
# read data containing icd diagnoses
all_data<-read.table(paste(working_dir,"../demographic/",root,"_sqc_fam_selectedPhenotypes_imagingSamples.csv",sep=""),header=TRUE,sep=",")
icd_cols<-colnames(all_data)[grep("ICD10",colnames(all_data))]
icd10_data<-all_data[,c("f.eid",icd_cols)]
rm(all_data)
## check counts of each category within main and secondary icd10 cols
icd10_codes<-c("R48","F80","F90", # reading
               paste("F",20:29,sep=""), # scz + types
               "Q89.3","Q24.0","Q24.8" # situs related
)
summary_icd10_counts<-sapply(icd10_codes,function(x){
  length(unique(as.vector(unlist(apply(icd10_data[,icd_cols],2,function(c){grep(x,c) } )))))
}  )
# flag the 8individuals with ICD10 code F20 (scz)
icd10_F20_IDs<-icd10_data$f.eid[unique(as.vector(unlist(apply(icd10_data[,icd_cols],2,function(c){grep("F20",c) } ))))]
# check that the selection is right
icd10_data[match(icd10_F20_IDs,icd10_data$f.eid),icd_cols]
# make a different flag columns per SCZ (F20) and Speech disorders (F80)
icd10_data$SCZ<-icd10_data$Speech<-0
icd10_data$SCZ[match(icd10_F20_IDs,icd10_data$f.eid)]<-1
icd10_data$Speech[unique(as.vector(unlist(apply(icd10_data[,icd_cols],2,function(c){grep("F80",c) } ))))]<-1
rm(icd_cols,icd10_codes)
#------------------------------------
# cognitive data
cogn_data<-read.table(paste("../",root,"_Cogn_all_noBioCovsPhenotypes_residuals_wHeader.table",sep=""),header=TRUE,sep="\t")
cogn_data<-subset(cogn_data,imaging_T1_QC==1)
cogn_cols<-grep("_T1|PC|Assessment",colnames(cogn_data))


#------------------------------------
# combine
cogn_icd10<-merge(cogn_data[,c(1,cogn_cols)],icd10_data,all=TRUE,by.x="ID1",by.y="f.eid")
pt_cogn<-merge(pt_data,cogn_icd10,all.x=TRUE)
pt_cogn$sex<-as.factor(pt_cogn$sex)
pt_cogn$handedness<-as.factor(pt_cogn$handedness)
levels(pt_cogn$handedness)<-c("RH","LH")
rm(icd10_data,cogn_data,pt_data,cogn_icd10,cogn_cols)
#------------------------------------

#------------------------------------
# create new variables related to PT
## PT
pt_cogn$AIdiff_VOLUME_Planum_Temporale <- pt_cogn$Volume_of_grey_matter_in_Planum_Temporale_left- pt_cogn$Volume_of_grey_matter_in_Planum_Temporale_right
pt_cogn$PT_LplusR<- pt_cogn$Volume_of_grey_matter_in_Planum_Temporale_right+pt_cogn$Volume_of_grey_matter_in_Planum_Temporale_left

## PTHG: Planum + Heschl's
pt_cogn$PTHG_left <- pt_cogn$Volume_of_grey_matter_in_Planum_Temporale_left+pt_cogn$Volume_of_grey_matter_in_Heschls_Gyrus_includes_H1_and_H2_left
pt_cogn$PTHG_right <- pt_cogn$Volume_of_grey_matter_in_Planum_Temporale_right+pt_cogn$Volume_of_grey_matter_in_Heschls_Gyrus_includes_H1_and_H2_right
pt_cogn$AIdiff_PTHG <- pt_cogn$PTHG_left-pt_cogn$PTHG_right
pt_cogn$AI_PTHG <- (pt_cogn$PTHG_left-pt_cogn$PTHG_right)/((pt_cogn$PTHG_left+pt_cogn$PTHG_right)/2)
pt_cogn$PTHG_LplusR<- pt_cogn$PTHG_left+pt_cogn$PTHG_right

# adjusted measures
x<-colnames(pt_cogn)[grep("Planum|PT",colnames(pt_cogn))]
cols2adj<-x[grep("fa|abs|residuals|plus|adj",x,invert=TRUE)]
rm(x)
pt_cols<-cols2adj[grep("Planum",cols2adj)]
pthg_cols<-cols2adj[grep("PTHG",cols2adj)]
# all adjusted for TBV
pt_cogn[,paste(cols2adj,"_adjTBV",sep="")]<-as.data.frame(do.call("cbind",lapply(cols2adj,function(x) adj_by_col(x,"totalBV",data=pt_cogn))))
# PT adjusted for PT L+R
pt_cogn[,paste(pt_cols,"_adjPTplus",sep="")] <-  as.data.frame(do.call("cbind",lapply(pt_cols,function(x) adj_by_col(x,"PT_LplusR",data=pt_cogn))))
# PTHG adjusted for PTHG L+R
pt_cogn[,paste(pthg_cols,"_adjPTHGplus",sep="")] <-  as.data.frame(do.call("cbind",lapply(pthg_cols,function(x) adj_by_col(x,"PTHG_LplusR",data=pt_cogn))))
# clean
rm(pt_cols,pthg_cols,cols2adj)


#------------------------------------
# Plots for rebuttal
p="Planum_Temporale"
#------------------------------------

#------------------------------------
# mean of AI for LH and RH
mu <- ddply(pt_cogn, "handedness", summarise, grp.mean=mean(AI_VOLUME_Planum_Temporale,na.rm=TRUE))
head(mu)
pt_ai_hand_hist<- ggplot(data=subset(pt_cogn,!is.na(handedness)),aes(x=AI_VOLUME_Planum_Temporale,linetype=handedness,color=handedness)) +
  # geom_rect(aes(xmin = 0, xmax = 1, ymin = 0, ymax = 5),
  #           fill = "grey", alpha = 0.05) +
  geom_density() +
  geom_vline(data=mu, aes(xintercept=grp.mean, linetype=handedness,color=handedness)) +
  geom_vline(xintercept=0, color="blue") +
  scale_color_manual(values=c("darkgreen","darkgrey")) +
  scale_linetype_manual(values=c("dashed","solid")) +
  xlab("AI - PT") +
  mytheme +
  NULL
rm(mu)

pt_ai_hand_hist

# mean of the SCZ, and speech (expressive) disorder individuals
ddply(pt_cogn, "SCZ", summarise, grp.mean=mean(AI_VOLUME_Planum_Temporale,na.rm=TRUE))
ddply(pt_cogn, "Speech", summarise, grp.mean=mean(AI_VOLUME_Planum_Temporale,na.rm=TRUE))
pt_ai_icd10_hist<- ggplot(data=subset(pt_cogn,!is.na(SCZ)),aes(x=AI_VOLUME_Planum_Temporale)) +
  # geom_rect(aes(xmin = 0, xmax = 1, ymin = 0, ymax = 5),
  #           fill = "grey", alpha = 0.05) +
  geom_density() +
  geom_point(data=subset(pt_cogn,SCZ==1),aes(y=0.8),color="red",alpha=0.3,size=3) +
  geom_point(data=subset(pt_cogn,Speech==1),aes(y=0.7),color="darkgreen",alpha=0.3, size=3) +
  geom_vline(xintercept=0, color="blue") +
  xlab("AI - PT") +
  geom_text(x=0.8,y=1,label="ICD10 codes",color="black",fontface=2) +
  geom_text(x=0.8,y=0.8,label="F20 (SCZ)",color="red") +
  geom_text(x=0.8,y=0.7,label="F80.1 (Expressive language disorder)",color="darkgreen") +
  mytheme +
  NULL
# save plots
ggsave(pt_ai_hand_hist,file=paste(out_dir,"/plots/AI_VOLUME_",p,"_hand_histograms",".pdf",sep=""),width=10,height=7)
ggsave(pt_ai_icd10_hist,file=paste(out_dir,"/plots/AI_VOLUME_",p,"_icd10_histograms",".pdf",sep=""),width=10,height=7)
rm(list=ls(pattern="dens|scatter|hist"))

#------------------------------------
# Distribution of EduYears and FluidInt

EduYears_hist<- ggplot(data=pt_cogn,aes(x=EduYears_max_T1)) +
  geom_histogram() +
  # geom_density() +
  xlab("EduYears") +
  mytheme +
  NULL

rEduYears_hist<- ggplot(data=pt_cogn,aes(x=residuals_EduYears_max_T1)) +
  geom_histogram() +
  # geom_density() +
  xlab("EduYears") +
  mytheme +
  NULL


FluidInt_hist<- ggplot(data=pt_cogn,aes(x=FluidInt00_T1)) +
  geom_histogram() +
  # geom_density() +
  xlab("FluidInt") +
  mytheme +
  NULL
rFluidInt_hist<- ggplot(data=pt_cogn,aes(x=residuals_FluidInt00_T1)) +
  geom_histogram() +
  # geom_density() +
  xlab("FluidInt - residuals") +
  mytheme +
  NULL

plot_grid(EduYears_hist,FluidInt_hist,
          # rEduYears_hist,rFluidInt_hist,
          ncol=2)

#------------------------------------
# Scatterplot of L vs R volumes
pt_LvsR_scatter<- ggplot(data=pt_cogn,aes(x=pt_cogn$Volume_of_grey_matter_in_Planum_Temporale_left,
                                          y=pt_cogn$Volume_of_grey_matter_in_Planum_Temporale_right)) +
  geom_abline(intercept=0,slope=1,color="blue") +
  geom_point(alpha=0.2) +
  # coord_cartesian(xlim=c(600,4000),ylim=c(600,3000)) +
  geom_density_2d() +
  labs(x="Left",y="Right") + 
  # mytheme
  NULL
# pt_LvsR_scatter
pt_L_dens<- ggplot(data=pt_cogn,aes(x=pt_cogn$Volume_of_grey_matter_in_Planum_Temporale_left)) +
  geom_density() +
  # coord_cartesian(xlim=c(600,4000)) +
  # labs(x="Left") + 
  theme(axis.title.x=element_blank(),
        axis.text=element_blank(),
        axis.line=element_blank(),
        axis.ticks=element_blank()) +
  NULL

pt_R_dens<- ggplot(data=pt_cogn,aes(x=pt_cogn$Volume_of_grey_matter_in_Planum_Temporale_right)) +
  geom_density() +
  # labs(x="Right") +
  coord_flip() +
  # coord_flip(xlim =c(600,3000) ) +
  theme(axis.title.y=element_blank(),
          axis.text=element_blank(),
          axis.line=element_blank(),
          axis.ticks=element_blank()) +
  NULL


dist1=pt_L_dens + theme(plot.margin = unit(c(0.5, 0, 0, 0.7), "cm"))
scatter1=pt_LvsR_scatter + theme(plot.margin = unit(c(0, 0, 0.5, 0.5), "cm"))
dist2=pt_R_dens + theme(plot.margin = unit(c(0, 0.5, 0.5, 0), "cm"))

# Combine all plots together and crush graph density with rel_heights
first_col = plot_grid(dist1, scatter1, ncol = 1, rel_heights = c(1, 3))
second_col = plot_grid(NULL, dist2, ncol = 1, rel_heights = c(1, 3))
pt_LvsR_scatter_density = plot_grid(first_col, second_col, ncol = 2, rel_widths = c(3, 1))
rm(first_col,second_col,dist1,dist2,scatter1)

# define boolean R>L
pt_cogn$AIdiff_PT_neg<-pt_cogn$AIdiff_VOLUME_Planum_Temporale<0
pt_LvsR_adjTBV_scatter<- ggplot(data=pt_cogn, #subset(pt_cogn,!is.na(AIdiff_PT_neg)),
                                aes(x=pt_cogn$Volume_of_grey_matter_in_Planum_Temporale_left_adjTBV,
                                    y=pt_cogn$Volume_of_grey_matter_in_Planum_Temporale_right_adjTBV)) +
                                    # color=as.factor(handedness))) +
  geom_abline(intercept=0,slope=1,color="blue") +
  geom_point(alpha=0.2, aes(color=AIdiff_PT_neg)) +
  geom_density_2d() +
  labs(x="Left - adjTBV",y="Right - adjTBV") + 
  scale_color_manual(values=c("darkgrey","red"),name="Right > Left",labels=c("FALSE","TRUE","")) +
  theme(legend.justification=c(0, 1), legend.position=c(0.8, 0.5)) +
  # mytheme
  NULL
# pt_LvsR_adjTBV_scatter


pt_LadjTBV_dens<- ggplot(data=pt_cogn,aes(x=pt_cogn$Volume_of_grey_matter_in_Planum_Temporale_left_adjTBV)) +
  geom_density() +
  # coord_cartesian(xlim=c(600,4000)) +
  # labs(x="Left") + 
  theme(axis.title.x=element_blank(),
        axis.text=element_blank(),
        axis.line=element_blank(),
        axis.ticks=element_blank()) +
  NULL

pt_RadjTBV_dens<- ggplot(data=pt_cogn,aes(x=pt_cogn$Volume_of_grey_matter_in_Planum_Temporale_right_adjTBV)) +
  geom_density() +
  # labs(x="Right") +
  coord_flip() +
  # coord_flip(xlim =c(600,3000) ) +
  theme(axis.title.y=element_blank(),
        axis.text=element_blank(),
        axis.line=element_blank(),
        axis.ticks=element_blank()) +
  NULL


dist1=pt_LadjTBV_dens + theme(plot.margin = unit(c(0.5, 0, 0, 0.7), "cm"))
scatter1=pt_LvsR_adjTBV_scatter + theme(plot.margin = unit(c(0, 0, 0.5, 0.5), "cm"))
dist2=pt_RadjTBV_dens + theme(plot.margin = unit(c(0, 0.5, 0.5, 0), "cm"))

# Combine all plots together and crush graph density with rel_heights
first_col = plot_grid(dist1, scatter1, ncol = 1, rel_heights = c(1, 3))
second_col = plot_grid(NULL, dist2, ncol = 1, rel_heights = c(1, 3))
pt_LvsR_adjTBV_scatter_density = plot_grid(first_col, second_col, ncol = 2, rel_widths = c(3, 1))
rm(first_col,second_col,dist1,dist2,scatter1)

pt_LvsR_all<-plot_grid(pt_LvsR_scatter_density,pt_LvsR_adjTBV_scatter_density,ncol=1) #,labels=c("A","B"))


# save plots
ggsave(pt_LvsR_scatter,file=paste(out_dir,"/plots/VOLUME_",p,"_LvsR_scatter",".png",sep=""))
ggsave(pt_LvsR_scatter_density,file=paste(out_dir,"/plots/VOLUME_",p,"_LvsR_scatter_densities",".png",sep=""))
ggsave(pt_LvsR_all,file=paste(out_dir,"/plots/VOLUME_",p,"_LvsR_adjTBV_scatter_densities",".png",sep=""),width = 10, height=8)
rm(list=ls(pattern="dens|scatter|hist"))

#------------------------------------
# Relation of L and R PT with age
pt_L_age_scatter<- ggplot(data=pt_cogn,aes(x=age,y=pt_cogn$Volume_of_grey_matter_in_Planum_Temporale_left,color=sex)) +
  geom_point(alpha=0.2) +
  # geom_density_2d() +
  geom_smooth(method = "lm", se = FALSE) +
  scale_color_manual(values=c("orange","darkgreen"),labels=c("male","female")) +
  labs(y="Left",x="age (years)") + 
  # theme(axis.title.x=element_blank(),
  #       axis.text=element_blank(),
  #       axis.line=element_blank(),
  #       axis.ticks=element_blank()) +
  NULL

pt_LadjTBV_age_scatter<- ggplot(data=pt_cogn,aes(x=age,y=pt_cogn$Volume_of_grey_matter_in_Planum_Temporale_left_adjTBV,color=sex)) +
  geom_point(alpha=0.2) +
  # geom_density_2d() +
  geom_smooth(method = "lm", se = FALSE) +
  scale_color_manual(values=c("orange","darkgreen"),labels=c("male","female")) +
  labs(y="Left - adjTBV",x="age (years)") + 
  # theme(axis.title.x=element_blank(),
  #       axis.text=element_blank(),
  #       axis.line=element_blank(),
  #       axis.ticks=element_blank()) +
  NULL

pt_R_age_scatter<- ggplot(data=pt_cogn,aes(x=age,y=pt_cogn$Volume_of_grey_matter_in_Planum_Temporale_right,color=sex)) +
  geom_point(alpha=0.2) +
  # geom_density_2d() +
  geom_smooth(method = "lm", se = FALSE) +
  scale_color_manual(values=c("orange","darkgreen"),labels=c("male","female")) +
  labs(y="Right",x="age (years)") + 
  # theme(axis.title.x=element_blank(),
  #       axis.text=element_blank(),
  #       axis.line=element_blank(),
  #       axis.ticks=element_blank()) +
  NULL

pt_RadjTBV_age_scatter<- ggplot(data=pt_cogn,aes(x=age,y=pt_cogn$Volume_of_grey_matter_in_Planum_Temporale_right_adjTBV,color=sex)) +
  geom_point(alpha=0.2) +
  # geom_density_2d() +
  geom_smooth(method = "lm", se = FALSE) +
  scale_color_manual(values=c("orange","darkgreen"),labels=c("male","female")) +
  labs(y="Right - adjTBV",x="age (years)") + 
  # theme(axis.title.x=element_blank(),
  #       axis.text=element_blank(),
  #       axis.line=element_blank(),
  #       axis.ticks=element_blank()) +
  NULL

pt_AI_age_scatter<- ggplot(data=pt_cogn,aes(x=age,y=AI_VOLUME_Planum_Temporale,color=sex)) +
  geom_point(alpha=0.2) +
  # geom_density_2d() +
  geom_smooth(method = "lm", se = FALSE) +
  scale_color_manual(values=c("orange","darkgreen"),labels=c("male","female")) +
  labs(y="AI",x="age (years)") + 
  # theme(axis.title.x=element_blank(),
  #       axis.text=element_blank(),
  #       axis.line=element_blank(),
  #       axis.ticks=element_blank()) +
  NULL

pt_AIadjTBV_age_scatter<- ggplot(data=pt_cogn,aes(x=age,y=AI_VOLUME_Planum_Temporale_adjTBV,color=sex)) +
  geom_point(alpha=0.2) +
  # geom_density_2d() +
  geom_smooth(method = "lm", se = FALSE) +
  scale_color_manual(values=c("orange","darkgreen"),labels=c("male","female")) +
  labs(y="AI - adjTBV",x="age (years)") + 
  # theme(axis.title.x=element_blank(),
  #       axis.text=element_blank(),
  #       axis.line=element_blank(),
  #       axis.ticks=element_blank()) +
  NULL


legend_b <- get_legend(pt_L_age_scatter + theme(legend.position="bottom"))
combined<-plot_grid(pt_L_age_scatter + theme(legend.position="none"),
          pt_R_age_scatter + theme(legend.position="none"),
          pt_AI_age_scatter + theme(legend.position="none"),
          pt_LadjTBV_age_scatter + theme(legend.position="none"),
          pt_RadjTBV_age_scatter + theme(legend.position="none"),
          pt_AIadjTBV_age_scatter + theme(legend.position="none"),
          ncol=3,
          labels=c("A","B","C","D","E","F"))
LR_age_scatter <- plot_grid( combined, legend_b, ncol = 1, rel_heights = c(1, .1))
# save plots and clean
ggsave(LR_age_scatter,file=paste(out_dir,"/plots/VOLUME_",p,"_LRAI_age_adjTBV_scatter",".png",sep=""),width = 10, height=8)
rm(list=ls(pattern="dens|scatter|hist"))
rm(legend_b,combined,p)
#--------------------------------------------


#--------------------------------------------
# Phenotypic relationship between PT measures and cognitive measures/disorders
#--------------------------------------------

# for disorders, N is too small
table(pt_cogn$SCZ,pt_cogn$Speech)
ddply(pt_cogn, "SCZ", summarise, grp.mean=mean(AI_VOLUME_Planum_Temporale,na.rm=TRUE))
ddply(pt_cogn, "Speech", summarise, grp.mean=mean(AI_VOLUME_Planum_Temporale,na.rm=TRUE))
#

# define columns to assess
pt_cols<-colnames(pt_cogn)[grep("Planum|PT|PTHG",colnames(pt_cogn))]
cogn_cols<-colnames(pt_cogn)[grep("max*_T1$|imaging*_T1$",colnames(pt_cogn))]
cols<-c(pt_cols,cogn_cols)
cols<-cols[grep("residuals|abs|fa|neg|LplusR",cols,invert=TRUE)]
rm(pt_cols,cogn_cols)

#
r0<-rcorr(as.matrix(pt_cogn[,c(cols)]),type="pearson")

# visualize some of these correlations
names<-colnames(r0$r)
pt_ai_w<-grep("AI*.*Planum",names)
pthg_ai_w<-grep("AI*.*PTHG",names)
pt_w<-grep("Planum",names)
pthg_w<-grep("PTHG",names)
cogn_w<-grep("Edu|Fluid",names)

# relationship between different measures of asymmetry
## pt
tr<-r0$r[pt_ai_w,pt_ai_w]
colnames(tr)<-gsub("PTplus","(L+R)",gsub("AIdiff","(L-R)",gsub("_VOLUME_Planum_Temporale","",colnames(tr))))
rownames(tr)<-colnames(tr)
corrplot(tr, type = "upper", method="number",
         tl.col = "black", tl.srt = 45)
rm(tr)

## pthg
corrplot(r0$r[pthg_ai_w,pthg_ai_w], type = "upper", method="number",
         tl.col = "black", tl.srt = 45)

# checke PT in relation to cogn variables
tr<-r0$r[pt_w,cogn_w]
rownames(tr)<-gsub("PTplus","(L+R)",gsub("AIdiff","(L-R)",gsub("_Planum_Temporale|_VOLUME|Volume_of_grey_matter_in_","",rownames(tr))))

corrplot(tr, type = "upper", method="number",
         tl.col = "black", tl.srt = 45)


## pthg
r0$r[pthg_ai_w,pthg_ai_w]
#pt measures and cognitive measures
r0$r[pt_w,cogn_w]
#pthg measures and cognitive measures
r0$r[pthg_w,cogn_w]

# convert into table
r1<-flattenCorrMatrix(r0$r, r0$P)
colnames(r1)<-c("p1","p2","cor","p")
# duplicate to get rows and columns interchanged (and then filter)
r2<-r1[,c("p2","p1","cor","p")]
colnames(r2)<-colnames(r1)
# combine
r<-rbind(r1,r2)
# get correlations across AI measures
r_pt<-r[intersect(grep("AI*.*Planum",r$p1),grep("AI*.*Planum",r$p2)),]
r_pthg<-r[intersect(grep("AI*.*PTHG",r$p1),grep("AI*.*PTHG",r$p2)),]

# filter
w1<-grep("PT|Planum",r$p1)
w2<-grep("PT|Planum",r$p2,invert=TRUE)
w<-intersect(w1,w2)
r<-r[w,]
# clean intermediate objects
rm(r1,r2,w1,w2,w)


# create new variables to organize output
r$phenotype<-r$measure<-r$adj<-NA
r$phenotype[grep("PTHG",r$p1)]<-"PTHG"
r$phenotype[grep("PTHG",r$p1,invert=TRUE)]<-"PT"
r$measure[grep("left",r$p1)]<-"L"
r$measure[grep("right",r$p1)]<-"R"
r$measure[grep("AIdiff",r$p1)]<-"AIdiff"
r$measure[grep("AIdiff|left|right",r$p1,invert=TRUE)]<-"AI"
r$adj[grep("adj",r$p1,invert=TRUE)]<-""
r$adj[grep("adjTBV",r$p1,)]<-"TBV"
r$adj[grep("plus",r$p1)]<-"SumRegion"

r$measure<-factor(r$measure,levels=c("L","R","AI","AIdiff"))
levels(r$measure)<-c("L","R","AI","L-R")

# organize and convert to wide format 
cols_order<-c("p2","phenotype","measure",paste(rep(c("","TBV_","SumRegion_"),each=2),rep(c("cor","p"),times=3),sep=""))

# all
r_wide<-r %>% select(cor,p,adj,measure,p2,phenotype) %>%
  gather(variable, value, cor:p)  %>%
  unite(temp, adj, variable) %>%
  spread(temp, value) %>%
  arrange(phenotype,p2,measure)
colnames(r_wide)<-gsub("^_","",colnames(r_wide))
r_wide<-r_wide[,cols_order]



p="Planum_Temporale"
write.csv(r_wide,file=paste(out_dir,p,"_correlation_matrix_EA_FI",".csv",sep=""),row.names=FALSE)

#--------------------------------------------
library("olsrr")

pt_cols<-colnames(pt_cogn)[grep("Planum|PT|PTHG",colnames(pt_cogn))]
pt_cols<-pt_cols[grep("residuals|abs|fa|neg|LplusR|adj",pt_cols,invert=TRUE)]


maxPC=10 # or 40...
scanner_position=c("Scanner_lateral_X_brain_position","Scanner_transverse_Y_brain_position","Scanner_longitudinal_Z_brain_position")
covs_order<-c("age","zage2","sex",
              # paste("PC",1:maxPC,sep=""),"array",
              scanner_position )
cogn_cols=c("EduYears_max_T1","FluidInt_imaging_T1")

lm_effect_cogn<-function(p,covs_order,cogn_cols){
  
  model_formula<-paste(p,"~",paste(c(covs_order),collapse=" + "))
  lm0<-do.call("lm",list (as.formula(model_formula),data=pt_cogn))
  # # evaluate the effect of extra covariates
  # lm1<-update(lm0,  . ~ . + 
  #               Volumetric_scaling_from_T1_head_image_to_standard_space +
  #               Volume_of_brain_grey_white +
  #               Volume_of_brain_stem_4th_ventricle +
  #               # Mean_tfMRI_head_motion_averaged_across_space_and_time_points + Mean_rfMRI_head_motion_averaged_across_space_and_time_points +
  #               height
  # )
  #--------------------------------------------
  # Assess effect of EA and FI
  #--------------------------------------------
  # EA
  # p2="EduYears_max_T1"
  
  
  l<-lapply(cogn_cols,function(p2){
    f1<-paste(model_formula, "+", p2)
    f2<-paste(model_formula, "+ totalBV +", p2)
    lm1<-do.call("lm",list (as.formula(f1),data=pt_cogn))
    lm2<-do.call("lm",list (as.formula(f2),data=pt_cogn))
    
    # F test from anova
    t1<-data.frame(anova(lm1))[p2,]
    t2<-data.frame(anova(lm2))[p2,]
    
    # # partial correlations
    # data.frame(ols_correlations(lm1)[p2,])
    # data.frame(ols_correlations(lm2)[p2,])
    
    
    # combine output into single summary table for all models
    summary_t<-data.frame(cbind(
      model=c("m1","m2"),
      N=c(length(resid(lm1)),length(resid(lm2))),
      rsq= c(summary(lm1)$r.squared,  summary(lm2)$r.squared),
      rbind(
        cbind(t1,data.frame(ols_correlations(lm1)[p2,])),
        cbind(t2,data.frame(ols_correlations(lm2)[p2,]))
      )
    ))
    summary_t$p1<-p
    summary_t$p2<-p2
    rownames(summary_t)<-1:2
    
    return(summary_t)
    
  })
  summary<-do.call("rbind",l)
  return(summary)
}


lms_summary<-do.call("rbind",lapply(pt_cols,function(p){lm_effect_cogn(p,covs_order,cogn_cols=c("EduYears_max_T1","FluidInt_imaging_T1"))}))
# create new variables to organize output
lms_summary$phenotype<-lms_summary$measure<-NA
lms_summary$phenotype[grep("PTHG",lms_summary$p1)]<-"PTHG"
lms_summary$phenotype[grep("PTHG",lms_summary$p1,invert=TRUE)]<-"PT"
lms_summary$measure[grep("left",lms_summary$p1)]<-"L"
lms_summary$measure[grep("right",lms_summary$p1)]<-"R"
lms_summary$measure[grep("AIdiff",lms_summary$p1)]<-"AIdiff"
lms_summary$measure[grep("AIdiff|left|right",lms_summary$p1,invert=TRUE)]<-"AI"

lms_summary$measure<-factor(lms_summary$measure,levels=c("L","R","AI","AIdiff"))
levels(lms_summary$measure)<-c("L","R","AI","L-R")


# all
lms_wide<-  lms_summary %>% select(-Df,-Sum.Sq,-Mean.Sq,-Zero.order,-Part) %>%
  gather(variable, value, N:Partial)  %>%
  unite(temp, model, variable) %>%
  spread(temp, value)  %>% 
  arrange(phenotype,p2,measure)

lms_wide<-lms_wide[,
                   c("phenotype","p2","measure",c(paste(rep(c("m1","m2"),each=5),rep(c("N","rsq","F.value","Pr..F.","Partial"),times=2),sep="_")))
                   ]


p="Planum_Temporale"
write.csv(lms_wide,file=paste(out_dir,p,"_lm_EA_FI",".csv",sep=""),row.names=FALSE)

# plots to illustrate these relationships
cols=c("#d95f02","#1f78b4","#33a02c","#000000")
# PT and EA (neg with AI, pos with R, nothing with L)
p_AI_EA<-ggplot(data=pt_cogn,aes(x=EduYears_max_T1,y=AI_VOLUME_Planum_Temporale)) + 
  geom_smooth(span=0.3,method='lm',colour="#d95f02") +
  # geom_violin(width=1.4) +
  # geom_boxplot(width=0.1, color="grey", alpha=0.2) +
  labs(x="EduYears",y="AI") +
  NULL
p_AIadjTBV_EA<-ggplot(data=pt_cogn,aes(x=EduYears_max_T1,y=AI_VOLUME_Planum_Temporale_adjTBV)) + 
  geom_smooth(span=0.3,method='lm',colour="#d95f02",linetype="dashed") +
  # geom_violin(width=1.4) +
  # geom_boxplot(width=0.1, color="grey", alpha=0.2) +
  labs(x="EduYears",y="AI - adj TBV") +
  NULL

p_AIadjPTp_EA<-ggplot(data=pt_cogn,aes(x=EduYears_max_T1,y=AI_VOLUME_Planum_Temporale_adjPTplus)) + 
  geom_smooth(span=0.3,method='lm',colour="#d95f02",linetype="dotted") +
  # geom_violin(width=1.4) +
  # geom_boxplot(width=0.1, color="grey", alpha=0.2) +
  labs(x="EduYears",y="AI - adj (L+R)") +
  NULL
#
p_AIdiff_EA<-ggplot(data=pt_cogn,aes(x=EduYears_max_T1,y=AIdiff_VOLUME_Planum_Temporale)) + 
  geom_smooth(span=0.3,method='lm',colour="#d95f02",alpha=0.3) +
  # geom_violin(width=1.4) +
  # geom_boxplot(width=0.1, color="grey", alpha=0.2) +
  labs(x="EduYears",y="(L-R)") +
  NULL
p_AIdiffadjTBV_EA<-ggplot(data=pt_cogn,aes(x=EduYears_max_T1,y=AIdiff_VOLUME_Planum_Temporale_adjTBV)) + 
  geom_smooth(span=0.3,method='lm',colour="#d95f02",linetype="dashed",alpha=0.3) +
  # geom_violin(width=1.4) +
  # geom_boxplot(width=0.1, color="grey", alpha=0.2) +
  labs(x="EduYears",y="(L-R) - adj TBV") +
  NULL

p_AIdiffadjPTp_EA<-ggplot(data=pt_cogn,aes(x=EduYears_max_T1,y=AIdiff_VOLUME_Planum_Temporale_adjPTplus)) + 
  geom_smooth(span=0.3,method='lm',colour="#d95f02",linetype="dotted",alpha=0.3) +
  # geom_violin(width=1.4) +
  # geom_boxplot(width=0.1, color="grey", alpha=0.2) +
  labs(x="EduYears",y="(L-R) - adj (L+R)") +
  NULL

#
p_RL_EA<-ggplot(data=pt_cogn,aes(x=EduYears_max_T1)) + 
  geom_smooth(span=0.3,method='lm',colour="#1f78b4", aes(y=Volume_of_grey_matter_in_Planum_Temporale_left)) +
  geom_smooth(span=0.3,method='lm',colour="#33a02c", aes(y=Volume_of_grey_matter_in_Planum_Temporale_right)) +
  # geom_violin(width=1.4) +
  # geom_boxplot(width=0.1, color="grey", alpha=0.2) +
  labs(x="EduYears",y="PT") +
  NULL

p_RLadjTBV_EA<-ggplot(data=pt_cogn,aes(x=EduYears_max_T1)) + 
  geom_smooth(span=0.3,method='lm',colour="#1f78b4", linetype="dashed", aes(y=Volume_of_grey_matter_in_Planum_Temporale_left_adjTBV)) +
  geom_smooth(span=0.3,method='lm',colour="#33a02c", linetype="dashed", aes(y=Volume_of_grey_matter_in_Planum_Temporale_right_adjTBV)) +
  # geom_violin(width=1.4) +
  # geom_boxplot(width=0.1, color="grey", alpha=0.2) +
  labs(x="EduYears",y="PT - adj TBV") +
  NULL

p_RLadjPTp_EA<-ggplot(data=pt_cogn,aes(x=EduYears_max_T1)) + 
  geom_smooth(span=0.3,method='lm',colour="#1f78b4", linetype="dotted", aes(y=Volume_of_grey_matter_in_Planum_Temporale_left_adjPTplus)) +
  geom_smooth(span=0.3,method='lm',colour="#33a02c", linetype="dotted", aes(y=Volume_of_grey_matter_in_Planum_Temporale_right_adjPTplus)) +
  # geom_violin(width=1.4) +
  # geom_boxplot(width=0.1, color="grey", alpha=0.2) +
  labs(x="EduYears",y="PT - adj (L+R)") +
  NULL

p_PT_EA_c<-plot_grid(p_RL_EA,p_RLadjTBV_EA,p_RLadjPTp_EA,
                     p_AI_EA,p_AIadjTBV_EA,p_AIadjPTp_EA,
                     p_AIdiff_EA,p_AIdiffadjTBV_EA,p_AIdiffadjPTp_EA,ncol=3)


p_rAI_EA<-ggplot(data=pt_cogn,aes(x=EduYears_max_T1,y=residuals_AI_VOLUME_Planum_Temporale)) + 
  geom_smooth(span=0.3,method='lm',colour="#d95f02") +
  # geom_violin(width=1.4) +
  # geom_boxplot(width=0.1, color="grey", alpha=0.2) +
  labs(x="EduYears",y="AI -residuals") +
  NULL

#
p_rRL_EA<-ggplot(data=pt_cogn,aes(x=EduYears_max_T1)) + 
  geom_smooth(span=0.3,method='lm',colour="#1f78b4", aes(y=residuals_Volume_of_grey_matter_in_Planum_Temporale_left)) +
  geom_smooth(span=0.3,method='lm',colour="#33a02c", aes(y=residuals_Volume_of_grey_matter_in_Planum_Temporale_right)) +
  # geom_violin(width=1.4) +
  # geom_boxplot(width=0.1, color="grey", alpha=0.2) +
  labs(x="EduYears",y="PT - residuals") +
  NULL

p_PT_EA_c3<-plot_grid(p_rRL_EA, p_rAI_EA,
                     ncol=2)

rm(p_RL_EA,p_RLadjTBV_EA,p_RLadjPTp_EA,
   p_AI_EA,p_AIadjTBV_EA,p_AIadjPTp_EA,
   p_AIdiff_EA,p_AIdiffadjTBV_EA,p_AIdiffadjPTp_EA)





# PT and FI (nothing with AI, pos with L and R)
p_AI_FI<-ggplot(data=pt_cogn,aes(x=FluidInt_imaging_T1,y=AI_VOLUME_Planum_Temporale)) + 
  geom_smooth(span=0.3,method='lm',colour="#d95f02") +
  # geom_violin(width=1.4) +
  # geom_boxplot(width=0.1, color="grey", alpha=0.2) +
  labs(x="Fluid Int",y="AI") +
  NULL
p_AIadjTBV_FI<-ggplot(data=pt_cogn,aes(x=FluidInt_imaging_T1,y=AI_VOLUME_Planum_Temporale_adjTBV)) + 
  geom_smooth(span=0.3,method='lm',colour="#d95f02",linetype="dashed") +
  # geom_violin(width=1.4) +
  # geom_boxplot(width=0.1, color="grey", alpha=0.2) +
  labs(x="Fluid Int",y="AI - adj TBV") +
  NULL

p_AIadjPTp_FI<-ggplot(data=pt_cogn,aes(x=FluidInt_imaging_T1,y=AI_VOLUME_Planum_Temporale_adjPTplus)) + 
  geom_smooth(span=0.3,method='lm',colour="#d95f02",linetype="dotted") +
  # geom_violin(width=1.4) +
  # geom_boxplot(width=0.1, color="grey", alpha=0.2) +
  labs(x="Fluid Int",y="AI - adj (L+R)") +
  NULL
#
p_AIdiff_FI<-ggplot(data=pt_cogn,aes(x=FluidInt_imaging_T1,y=AIdiff_VOLUME_Planum_Temporale)) + 
  geom_smooth(span=0.3,method='lm',colour="#d95f02",alpha=0.3) +
  # geom_violin(width=1.4) +
  # geom_boxplot(width=0.1, color="grey", alpha=0.2) +
  labs(x="Fluid Int",y="(L-R)") +
  NULL
p_AIdiffadjTBV_FI<-ggplot(data=pt_cogn,aes(x=FluidInt_imaging_T1,y=AIdiff_VOLUME_Planum_Temporale_adjTBV)) + 
  geom_smooth(span=0.3,method='lm',colour="#d95f02",linetype="dashed",alpha=0.3) +
  # geom_violin(width=1.4) +
  # geom_boxplot(width=0.1, color="grey", alpha=0.2) +
  labs(x="Fluid Int",y="(L-R) - adj TBV") +
  NULL

p_AIdiffadjPTp_FI<-ggplot(data=pt_cogn,aes(x=FluidInt_imaging_T1,y=AIdiff_VOLUME_Planum_Temporale_adjPTplus)) + 
  geom_smooth(span=0.3,method='lm',colour="#d95f02",linetype="dotted",alpha=0.3) +
  # geom_violin(width=1.4) +
  # geom_boxplot(width=0.1, color="grey", alpha=0.2) +
  labs(x="Fluid Int",y="(L-R) - adj (L+R)") +
  NULL

#
p_RL_FI<-ggplot(data=pt_cogn,aes(x=FluidInt_imaging_T1)) + 
  geom_smooth(span=0.3,method='lm',colour="#1f78b4", aes(y=Volume_of_grey_matter_in_Planum_Temporale_left)) +
  geom_smooth(span=0.3,method='lm',colour="#33a02c", aes(y=Volume_of_grey_matter_in_Planum_Temporale_right)) +
  # geom_violin(width=1.4) +
  # geom_boxplot(width=0.1, color="grey", alpha=0.2) +
  labs(x="Fluid Int",y="PT") +
  NULL

p_RLadjTBV_FI<-ggplot(data=pt_cogn,aes(x=FluidInt_imaging_T1)) + 
  geom_smooth(span=0.3,method='lm',colour="#1f78b4", linetype="dashed", aes(y=Volume_of_grey_matter_in_Planum_Temporale_left_adjTBV)) +
  geom_smooth(span=0.3,method='lm',colour="#33a02c", linetype="dashed", aes(y=Volume_of_grey_matter_in_Planum_Temporale_right_adjTBV)) +
  # geom_violin(width=1.4) +
  # geom_boxplot(width=0.1, color="grey", alpha=0.2) +
  labs(x="Fluid Int",y="PT - adj TBV") +
  NULL

p_RLadjPTp_FI<-ggplot(data=pt_cogn,aes(x=FluidInt_imaging_T1)) + 
  geom_smooth(span=0.3,method='lm',colour="#1f78b4", linetype="dotted", aes(y=Volume_of_grey_matter_in_Planum_Temporale_left_adjPTplus)) +
  geom_smooth(span=0.3,method='lm',colour="#33a02c", linetype="dotted", aes(y=Volume_of_grey_matter_in_Planum_Temporale_right_adjPTplus)) +
  # geom_violin(width=1.4) +
  # geom_boxplot(width=0.1, color="grey", alpha=0.2) +
  labs(x="Fluid Int",y="PT - adj (L+R)") +
  NULL

p_PT_FI_c<-plot_grid(p_RL_FI,p_RLadjTBV_FI,p_RLadjPTp_FI,
                     p_AI_FI,p_AIadjTBV_FI,p_AIadjPTp_FI,
                     p_AIdiff_FI,p_AIdiffadjTBV_FI,p_AIdiffadjPTp_FI,
                     ncol=3)

rm(p_RL_FI,p_RLadjTBV_FI,p_RLadjPTp_FI,
   p_AI_FI,p_AIadjTBV_FI,p_AIadjPTp_FI,
   p_AIdiff_FI,p_AIdiffadjTBV_FI,p_AIdiffadjPTp_FI)

# assess effect with handedness via regression
lms_summary<-do.call("rbind",lapply(pt_cols,function(p){lm_effect_cogn(p,covs_order,cogn_cols=c("handedness"))}))
# create new variables to organize output
lms_summary$phenotype<-lms_summary$measure<-NA
lms_summary$phenotype[grep("PTHG",lms_summary$p1)]<-"PTHG"
lms_summary$phenotype[grep("PTHG",lms_summary$p1,invert=TRUE)]<-"PT"
lms_summary$measure[grep("left",lms_summary$p1)]<-"L"
lms_summary$measure[grep("right",lms_summary$p1)]<-"R"
lms_summary$measure[grep("AIdiff",lms_summary$p1)]<-"AIdiff"
lms_summary$measure[grep("AIdiff|left|right",lms_summary$p1,invert=TRUE)]<-"AI"

lms_summary$measure<-factor(lms_summary$measure,levels=c("L","R","AI","AIdiff"))
levels(lms_summary$measure)<-c("L","R","AI","L-R")


# all
lms_wide<-  lms_summary %>% select(-Df,-Sum.Sq,-Mean.Sq,-Zero.order,-Part) %>%
  gather(variable, value, N:Partial)  %>%
  unite(temp, model, variable) %>%
  spread(temp, value)  %>% 
  arrange(phenotype,p2,measure)

lms_wide<-lms_wide[,
                   c("phenotype","p2","measure",c(paste(rep(c("m1","m2"),each=5),rep(c("N","rsq","F.value","Pr..F.","Partial"),times=2),sep="_")))
                   ]


p="Planum_Temporale"
write.csv(lms_wide,file=paste(out_dir,p,"_lm_handedness",".csv",sep=""),row.names=FALSE)
