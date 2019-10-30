library(ggpubr)
library(ggplot2)
library(reshape2)
library(grid)
library(gridExtra)
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
# define working_dir
if (Sys.info()['sysname']=='Windows') {dir="P://workspaces/"} else {dir="/data/workspaces/lag/workspaces/"}
# dir="\\\\data/lag/workspaces/"
out_dir=paste(dir,"lg-ukbiobank/working_data/amaia/genetic_data/release_v2/cal/h2/",subset_name,"/output/PT/",sep="")
#
dir.create(file.path(out_dir), showWarnings = FALSE)
# out_dir="K://written/tex/multilateral/presentations/figures/"
setwd(out_dir)

# define pattern for this run
root="ukb25465_ukb25468"
pattern_run= "noBioCovs_noAssessmentC"
pattern_run2=paste("_",root,"_imaging_",pattern_run,sep="")
#----------------------------------------------------------------------
summary<-read.csv(paste(out_dir,"/../summary_AIvolumes_h2_REML",pattern_run2,".csv",sep=""))
#----------------------------------------------------------------------
pt<-subset(summary,region=="Planum Temporale")
pt$estimate<-pt$h2; pt$estimate[is.na(pt$h2)]<- pt$rG[is.na(pt$h2)]
pt$estimate<-as.numeric(pt$estimate)
pt$se<-as.numeric(pt$se)
pt$measure<-gsub("rG","L,R",pt$measure)
rm(summary)

#----------------------------------------------------------------------
# plot heritabilities from h2
gcta_h2_barplot<-ggplot(data=subset(pt,stat=="h2"&sample=="Total"),color="black",aes(x=measure,y=estimate,width=0.7,fill=measure)) + 
  geom_bar(aes(y=estimate),position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=estimate-se, ymax=estimate+se,width=.5), position=position_dodge(.7)) +
  mytheme + 
  scale_fill_manual(values=cols,name = "") + coord_cartesian(ylim = c(0,1)) +
  ylab(bquote('Estimate ('*h^2*')')) +
  theme(axis.text.x = element_text(angle = 0)) + theme(legend.position="none") 
  
gcta_rg_barplot<-ggplot(data=subset(pt,stat=="rG"&sample=="Total"),color="black",aes(x=measure,y=estimate,width=0.7,fill=measure)) + 
  geom_bar(aes(y=estimate),position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=estimate-se, ymax=estimate+se,width=.5), position=position_dodge(.7)) +
  mytheme + 
  scale_fill_manual(values=cols[4],name = "")  + coord_cartesian(ylim = c(0,1)) +
  ylab(bquote('Estimate ('*rho*')')) +
  theme(legend.position="none") 
#----------------------------------------------------------------------
png(paste(out_dir,"gcta_estimates_PT",pattern_run2,".png",sep=""))
gcta_estimates_barplot<-grid.arrange(gcta_h2_barplot,gcta_rg_barplot,ncol=2,widths=c(3,1.25), top = textGrob("Planum Temporale",gp=gpar(fontsize=20,font=3)))
dev.off()
write.csv(pt,paste(out_dir,"gcta_estimates_PT",pattern_run2,".csv",sep=""),row.names=FALSE)
#----------------------------------------------------------------------


# per sex
# plot heritabilities from h2
gcta_h2_barplot_sex<-ggplot(data=subset(pt,stat=="h2"),color="black",aes(x=measure,y=estimate,width=0.7,fill=measure,alpha=sample)) + 
  geom_bar(aes(y=estimate),position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=estimate-se, ymax=estimate+se,width=.5), position=position_dodge(.7)) +
  mytheme + 
  scale_fill_manual(values=cols,name = "") + coord_cartesian(ylim = c(0,1)) +
  scale_alpha_discrete( range=c(0.4,1), na.value = 0) +
  ylab(bquote('Estimate ('*h^2*')')) +
  theme(axis.text.x = element_text(angle = 0)) # + theme(legend.position="none") 

gcta_rg_barplot_sex<-ggplot(data=subset(pt,stat=="rG"),color="black",aes(x=measure,y=estimate,width=0.7,fill=measure,alpha=sample)) +  
  geom_bar(aes(y=estimate),position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=estimate-se, ymax=estimate+se,width=.5), position=position_dodge(.7)) +
  mytheme + 
  scale_fill_manual(values=cols[4],name = "")  + coord_cartesian(ylim = c(0,1)) +
  scale_alpha_discrete( range=c(0.4,1), na.value = 0) +
  ylab(bquote('Estimate ('*rho*')')) +
  theme(legend.position="none") 

png(paste(out_dir,"gcta_estimates_PT_sex",pattern_run2,".png",sep=""))
gcta_estimates_barplot_sex<-grid.arrange(gcta_h2_barplot_sex,gcta_rg_barplot_sex,ncol=2,widths=c(3,1.25), top = textGrob("Planum Temporale",gp=gpar(fontsize=20,font=3)))
dev.off()
