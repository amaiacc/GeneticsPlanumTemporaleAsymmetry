
#----------------------------------------------------------------------
# define working_dirs
if (Sys.info()['sysname']=='Windows') {dir="P://workspaces/"} else {dir="/data/workspaces/lag/workspaces/"}
#----------------------------------------------------------------------
library(ggplot2)
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

mytheme2<-theme_bw() + theme(panel.spacing = unit(0, "lines"), 
                             strip.background = element_rect(fill="white"), strip.text = element_text(size=16), 
                             title=element_text(size=16),
                             axis.title=element_text(size=24),axis.text=element_text(size=24,colour="black"),
                             axis.title.x=element_blank(),
                             legend.text=element_text(size=24), legend.title=element_blank(), legend.position="bottom") 



#--------------------------------------------
# Define directories, dependent on  system
if (Sys.info()['sysname']=='Windows') {dir="P://workspaces/"} else {dir="/data/workspaces/lag/workspaces/"}
primary_dir=paste(dir,"lg-ukbiobank/primary_data/",sep="")
working_dir=paste(dir,"lg-ukbiobank/working_data/amaia/",sep="")
#
subset_name="imagingT1_N18057"
root="ukb25465_ukb25468"
pattern_run="_noBioCovs_noAssessmentC"

#
h2_dir=paste(dir,"lg-ukbiobank/working_data/amaia/genetic_data/release_v2/cal/h2/",subset_name,"/output/",sep="")
#
lm_dir=paste(working_dir,"/pheno_files/genetic_v2/",root,"/lm/",sep="")
out_dir=paste(working_dir,"/pheno_files/genetic_v2/",root,"/summary_phenotypes/PT/",sep="")
out_plots_dir=paste(out_dir,"plots/",sep="")
# create directory for plots and files
dir.create(file.path(out_dir), showWarnings = FALSE)
dir.create(file.path(out_plots_dir), showWarnings = FALSE)
#--------------------------------------------
# Set working directory, and define general options
opts_knit$set(root.dir = out_dir)
options(stringsAsFactors = FALSE)
#--------------------------------------------
# Read data
#--------------------------------------------
setwd(out_dir) # just in case it's not knitted
# read data
all_data<-read.table(paste("../",root,"_imaging",pattern_run,"VolumePhenotypes_residuals_wHeader.table",sep=""),header=TRUE,sep="\t")
p="Planum_Temporale"
pt_cols<-grep(p,colnames(all_data))
pt<-all_data[c(1:43,pt_cols)];rm(pt_cols)
#--------------------------------------------
# histogram with shade
pt_ai<- ggplot(data=pt,aes(x=AI_VOLUME_Planum_Temporale)) +
        geom_rect(aes(xmin = 0, xmax = 1, ymin = -Inf, ymax = Inf),
            fill = "grey", alpha = 0.8) +
        geom_histogram() + xlab("AI") + 
        coord_cartesian(xlim = c(-0.5,1)) +
        geom_vline(xintercept=0,color="black") +
        geom_vline(xintercept=mean(pt$AI_VOLUME_Planum_Temporale,na.rm=T),color="red") +
       mytheme + theme(axis.title.y=element_blank())

pt_rai<- ggplot(data=pt,aes(x=residuals_AI_VOLUME_Planum_Temporale)) +
  # geom_rect(aes(xmin = -1, xmax = 0, ymin = -Inf, ymax = Inf),fill = "blue", alpha = 0.03) +
  geom_histogram() + xlab("AI - residuals") + 
  coord_cartesian(xlim = c(-1,1))+
  mytheme + theme(axis.title.y=element_blank())
# lateral volumes, 

pt_l<- ggplot(data=pt,aes(x=Volume_of_grey_matter_in_Planum_Temporale_left)) +
  geom_histogram() + xlab("Left") +
  geom_vline(xintercept=mean(pt$Volume_of_grey_matter_in_Planum_Temporale_left,na.rm=T),color="red") +
  mytheme + theme(axis.title.y=element_blank())

pt_rl<- ggplot(data=pt,aes(x=residuals_Volume_of_grey_matter_in_Planum_Temporale_left)) +
  # geom_rect(aes(xmin = -1, xmax = 0, ymin = -Inf, ymax = Inf),fill = "blue", alpha = 0.03) +
  geom_histogram() + xlab("Left - residuals") + 
  mytheme + theme(axis.title.y=element_blank())

pt_r<- ggplot(data=pt,aes(x=Volume_of_grey_matter_in_Planum_Temporale_right)) +
  geom_histogram() + xlab("Right") +
  geom_vline(xintercept=mean(pt$Volume_of_grey_matter_in_Planum_Temporale_right,na.rm=T),color="red") +
  mytheme + theme(axis.title.y=element_blank())

pt_rr<- ggplot(data=pt,aes(x=residuals_Volume_of_grey_matter_in_Planum_Temporale_right)) +
  # geom_rect(aes(xmin = -1, xmax = 0, ymin = -Inf, ymax = Inf),fill = "blue", alpha = 0.03) +
  geom_histogram() + xlab("Right - residuals ") + 
  mytheme + theme(axis.title.y=element_blank())



ggsave(grid.arrange(pt_l,pt_r,pt_ai,
                    pt_rl,pt_rr,pt_rai,
                    ncol=3),file=paste(out_plots_dir,p,"_histograms_shade",pattern_run,".pdf",sep=""),width=10,height=7)

#--------------------------------------------
pt$sex<-as.factor(pt$sex)
levels(pt$sex)<-c("males","females")
pt_ai_sex<- ggplot(data=pt,aes(y=AI_VOLUME_Planum_Temporale,x=sex)) +
  # geom_rect(aes(xmin = -Inf, xmax = Inf,ymin = 0, ymax = 1),  fill = "grey", alpha = 0.8) +
  geom_violin() +
  xlab("") + ylab("AI PT") +
  coord_cartesian(ylim = c(-0.5,1)) +
  # geom_hline(yintercept=0,color="black") +
  # geom_vline(xintercept=mean(pt$AI_VOLUME_Planum_Temporale,na.rm=T),color="red") +
  mytheme 
pt_ai_sex2<-pt_ai_sex + 
    geom_dotplot(binaxis='y', stackdir='center', dotsize=0.01) +
    stat_summary(fun.y=mean, geom="point", size=0.05, color="red")
    # stat_summary(fun.data=mean, geom="pointrange", color="red",size=0.05)


pt_ai_sex3<-pt_ai_sex + stat_summary(fun.data="mean_sdl", 
               geom="crossbar", width=0.05 )


library(plyr)
mu <- ddply(pt, "sex", summarise, grp.mean=mean(AI_VOLUME_Planum_Temporale,na.rm=TRUE))
head(mu)

pt_ai_sex_hist<- ggplot(data=pt,aes(x=AI_VOLUME_Planum_Temporale,linetype=sex)) +
  # geom_rect(aes(xmin = 0, xmax = 1, ymin = 0, ymax = 5),
  #           fill = "grey", alpha = 0.05) +
        geom_density() +
        geom_vline(data=mu, aes(xintercept=grp.mean, linetype=sex)) +
        geom_vline(xintercept=0, color="blue") +
       
        scale_linetype_manual(values=c("dashed","solid")) +
        xlab("AI - PT") +
        mytheme
        
ggsave(pt_ai_sex_hist,file=paste(out_plots_dir,"AI_VOLUME_",p,"_sex_histograms",pattern_run,".pdf",sep=""),width=10,height=7)

lm(AI_VOLUME_Planum_Temporale~sex,data=pt)

lm_base<-lm(AI_VOLUME_Planum_Temporale~age+zage2+sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+
              array+
              Scanner_lateral_X_brain_position+
              Scanner_transverse_Y_brain_position+
              Scanner_longitudinal_Z_brain_position+
              totalBV,data=subset(all_data,!is.na(handedness)))

# check effect of handedness
all_data$handedness<-as.factor(all_data$handedness)
all_data$sex<-as.factor(all_data$sex)
lm_hand<-update(lm_base,.~.+handedness)
summary(lm_hand)
anova(lm_base,lm_hand)
        