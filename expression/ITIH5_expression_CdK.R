
#----------------------------------------------------------------------
# define working_dirs
if (Sys.info()['sysname']=='Windows') {dir="P://workspaces/"} else {dir="/data/workspaces/lag/workspaces/"}
#----------------------------------------------------------------------
library(ggplot2);library(gridExtra);library(grid)
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


cols=c("#1f78b4","#33a02c","#000000") # '#1b9e77','#d95f02','#7570b3 # "#a6cee3"  ,"#7570b3"
#--------------------------------------------
# Define directories, dependent on  system
if (Sys.info()['sysname']=='Windows') {dir="P://workspaces/"} else {dir="/data/workspaces/lag/workspaces/"}
primary_dir=paste(dir,"lg-ukbiobank/primary_data/",sep="")
working_dir=paste(dir,"lg-ukbiobank/working_data/amaia/genetic_data/PT/",sep="")
setwd(working_dir)
#--------------------------------------------
itih5<-read.csv("itih5_for_amaia_CdK.csv")
# some harmonization of variables
itih5$sex[grep("F",itih5$sex)]<-"F"
itih5$sex[grep("M",itih5$sex)]<-"M"
itih5$sex<-as.factor(itih5$sex)


itih5$age<-itih5$age.pc.weeks
itih5$age.pc.weeks<-as.numeric(itih5$age)
itih5$age.years<-NA
itih5$age.years[grep("Y",itih5$age)]<-as.numeric(gsub("Y","",itih5$age[grep("Y",itih5$age)]))


age.pc_plot<-ggplot(data=itih5) + geom_point(aes(x=age.pc.weeks,y=log_cpm_ITIH5,shape=sex,color=side)) + facet_grid(Dataset~.) +
        scale_color_manual(values=cols) +
        mytheme    + theme(legend.position="none")
   
age.y_plot<-ggplot(data=subset(itih5,Dataset=="Kang")) + geom_point(aes(x=age.years,y=log_cpm_ITIH5,shape=sex,color=side)) + facet_grid(Dataset~.) +
  scale_color_manual(values=cols) +
  mytheme    + theme(legend.position="bottom")


age_plot<-grid.arrange(age.pc_plot,age.y_plot,heights=c(3,1.5),top=textGrob( "ITIH5 expression in temporal areas",gp=gpar(fontsize=16,font=2))
                        )
ggsave(age_plot,file="ITIH5_expression_brain_temporal_CdK.png",height = 8,width=5)




