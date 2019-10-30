#' ---
#' title: ukb25465_ukb25468 combined - Runs LM/plots/checks/ and saves residuals for cognitive traits (EA/FI) in all subsets (all, imaging, non-imaging)
#' author: amacar
#' date: "Edited from ukb9246_extractPhenotypes_imagingSubset.R; . Modified: `r Sys.time()`"
#' output:
#'   html_document:
#'     toc: true
#'     toc_depth: 2
#'     theme: "flatly"
#'     highlight: "textmate"
#'   pdf_document:
#'     keep_tex: true
#' ---

## NOTE: it does not include biological covariates (totalBV, height and BMI) that had been included previously.


library(pander)
library(ggplot2)
library(car)
library(dplyr)
library(tidyr)
library(grid)
library(gridExtra)
library(knitr)
library("GGally")

opts_chunk$set(include=TRUE, warnings=FALSE, echo=FALSE, results = "asis", tidy=TRUE, width=50, fig.width=8, fig.height=6)

panderOptions('knitr.auto.asis',FALSE)
# opts_chunk$set(dev="pdf", 
#                dev.args=list(type="cairo"),
#                dpi=96)

set.seed(50)

release="ukb25465_ukb25468"
# Define directories, dependent on  system
# if (Sys.info()['sysname']=='Windows') {dir="P://workspaces/"} else {dir="/data/workspaces/lag/workspaces/"}
primary_dir=paste(dir,"lg-ukbiobank/primary_data",sep="")
working_dir=paste(dir,"lg-ukbiobank/working_data/amaia/",sep="")
out_dir=paste(working_dir,"/pheno_files/genetic_v2/",release,"/",sep="")
out_lmplots_dir=paste(out_dir,"/lm_plots/",sep="")
out_lm_dir=paste(out_dir,"/lm/",sep="")

# create ouput directories if they don't already exist
dir.create(file.path(out_dir), showWarnings = FALSE)
dir.create(file.path(out_lmplots_dir), showWarnings = FALSE)
dir.create(file.path(out_lm_dir), showWarnings = FALSE)

# Set working directory
opts_knit$set(root.dir = out_dir)

# Define options
maxPC=40 #40
covs_name="_noBioCovs"
# t, threshold to define outliers: +-t*SD
t=6

#' Start!
setwd(out_dir)

#' check if RData/data file exists, load, otherwise print warning
rdata<-paste(out_dir,release,"_extractPhenotypesCogn_PhenoCovs_subsets.RData",sep="")
# fdata<-paste(out_dir,"pheno_files/genetic_v2/","ukb21288_ukb21293_imaging_extracted_PhenoCovs.table",sep="")
if (file.exists(rdata)) {
load(rdata) 
  } else {
  cat("Warning!!\nThe data you tried to load does not exist.\nYou may have to run ukb21288_ukb21293_extractPhenotypesCogn_all.R")
  stop()
  }

#' Define all possible covariates
pcs<-paste("PC",1:maxPC,sep="")
zage2_cols<-colnames(data_s)[grep("zage2",colnames(data_s))]
#
bin_covs<-c("sex","assessment_centre","assessment_centre_imaging","array","batch")
quant_covs<-c("age","age_imaging",zage2_cols,pcs)

#' Define phenotypes of interest
phenos<-c("EduYears00","EduYears_max","FluidInt00","FluidInt_imaging","FluidInt_mean")

#' For each phenotype, create new variables that contain QC subsets for: 
#' all, imaging (T1) and nonimaging (nonT1)
#' define indices per subset
w_all<-which(data_s$all_QC==1)
w_T1<-which(data_s$imaging_T1_QC==1)
w_nonT1<-which(data_s$non_imaging_T1_QC==1)
for (p in phenos){
  # create new cols; filled in with total data for phenoytpe x
  data_s[,paste(p,c("_all","_nonT1","_T1"),sep="")]<-data_s[,p]
  # blank relevant cols
  data_s[-w_all,paste(p,"_all",sep="")]<-NA
  data_s[-w_T1,paste(p,"_T1",sep="")]<-NA
  data_s[-w_nonT1,paste(p,"_nonT1",sep="")]<-NA
}
rm(p)
rm(w_all,w_T1,w_nonT1)

#' Check missingness pattern per phenotype
w<-grep(paste(phenos,collapse="|"),colnames(data_s))
phenos_all<-colnames(data_s)[w]
do.call("rbind",lapply(data_s[,phenos_all],function(x) table(is.na(x))))
rm(w)


#' Check mean
w<-grep(paste(phenos,collapse="|"),colnames(data_s))
phenos_all<-colnames(data_s)[w]
do.call("rbind",lapply(data_s[,phenos_all],function(x) mean(x,na.rm=T)))
rm(w)


#--------------------------------------
#' ## Regress out covariates from all the phenotypes, and save residuals. 
#--------------------------------------
#' Define order of covariates to include within lm, which will be dependent on the subset
#' Drop batch, too many levels # ,"batch"
#' Include assessement center as well

# define covariate order as well
## per subset
covs_order_all<-c("assessment_centre","age","zage2_all","sex",pcs,"array")
covs_order_T1<-c("assessment_centre_imaging","age_imaging","zage2_imagingT1","sex",pcs,"array")
covs_order_nonT1<-c("assessment_centre","age","zage2_non_imagingT1","sex",pcs,"array")

# run all LMs, but for imaging: only the imaging will be used? maximize N
for (s in c("_all","_nonT1","_T1")){
  covs_order<-get(paste("covs_order",s,sep=""))
  #' Define covariates of interest, to run LM
  quant_covs2<-quant_covs[quant_covs %in% covs_order]
  bin_covs2<-bin_covs[bin_covs %in% covs_order]
  # define age cols
  age<-covs_order[grep("^age",covs_order)]
  zage2<-covs_order[grep("^zage2",covs_order)]
    
  for (phen in phenos) {
    print(paste("Dependent variable: ",phen,"\n",sep=""))
    print(paste("subset: ",gsub("_","",s),"\n ",sep=""))
    # define column to test
    p=paste(phen,s,sep="")
    #----------------------
    # for each phenotype, first check if output exists and only proceed if it does not
    lm_file=paste(out_lm_dir,p,"_lm",covs_name,".txt",sep="")
    if (file.exists(lm_file)==FALSE) {
      #----------------------
      print('Detect and blank outliers"\n  ')
      m<-mean(data_s[,p],na.rm=TRUE)
      sd<-sd(data_s[,p],na.rm=TRUE)
      out<-which( data_s[,p] > m+t*sd |data_s[,p] < m-t*sd )
      print('Remove outliers: ')
      print(length(out))
      print('')
      if (length(out)>0){
        data_s[out,p]<-NA
      }
      rm(m,sd,out)
      #----------------------
      #print('Regress out covariates and save residuals"\n  ')
      print(paste('Run a linear model to get residuals for ',p,'\n  ',sep=""))
      # model_formula<-paste(p,"~",paste(c(bin_covs2,quant_covs2),collapse=" + "))
      model_formula<-paste(p,"~",paste(c(covs_order),collapse=" + "))
      model_formula2<-paste(model_formula,"+ sex*",age," + sex*",zage2,sep="")
      if (class(data_s[,p])=="numeric" | class(data_s[,p])=="integer") {
        lm<-do.call("lm",list (as.formula(model_formula),data=data_s))
        lm2<-do.call("lm",list (as.formula(model_formula2),data=data_s))
      } else if (class(data_s[,p])=="factor") {
        lm<-do.call("glm",list (as.formula(model_formula),data=data_s,family=binomial(link='logit')))
        lm2<-do.call("glm",list (as.formula(model_formula2),data=data_s,family=binomial(link='logit')))
      }
      pander(anova(lm),round=1000)
      pander(anova(lm2),round=1000)
      #----------
      print('Fit check\n  ')
      pdf(file=paste(out_lmplots_dir,p,"_lm_plots",covs_name,".pdf",sep=""))
      # par(mfrow=c(2,2))
      # png(file=paste(out_lmplots_dir,p,"_lm_plots",covs_name,".png",sep=""))
      plot(lm)
      dev.off()
      pdf(file=paste(out_lmplots_dir,p,"_lm2_plots",covs_name,".pdf",sep=""))
      # png(file=paste(out_lmplots_dir,p,"_lm2_plots",covs_name,".png",sep=""))
      plot(lm2)
      dev.off()
      par(mfrow=c(1,1))
      #----------
      # lm and lm2
      # save model + summary model estimates in txt file as well
      df<-as.data.frame(anova(lm)); df$variable<-rownames(df)
      df<-df[,c(NCOL(df),1:(NCOL(df)-1))]
      df2<-as.data.frame(summary(lm)$coefficients); df2$variable<-rownames(df2)
      df2<-df2[,c(NCOL(df2),1:(NCOL(df2)-1))]
      write.table(df,file=paste(out_lm_dir,p,"_lm",covs_name,".txt",sep=""),sep="\t",row.names=FALSE,quote=FALSE)
      write.table(df2,file=paste(out_lm_dir,p,"_coefficients_lm",covs_name,".txt",sep=""),sep="\t",row.names=FALSE,quote=FALSE)
      rm(df,df2)
      #----------
      df<-as.data.frame(anova(lm2)); df$variable<-rownames(df)
      df<-df[,c(NCOL(df),1:(NCOL(df)-1))]
      df2<-as.data.frame(summary(lm2)$coefficients); df2$variable<-rownames(df2)
      df2<-df2[,c(NCOL(df2),1:(NCOL(df2)-1))]
      write.table(df,file=paste(out_lm_dir,p,"_lm2",covs_name,".txt",sep=""),sep="\t",row.names=FALSE,quote=FALSE)
      write.table(df2,file=paste(out_lm_dir,p,"_coefficients_lm2",covs_name,".txt",sep=""),sep="\t",row.names=FALSE,quote=FALSE)
      rm(df,df2)
      #----------
      tmp_df<-data_s[,c("ID1","ID2",p,bin_covs2,quant_covs2)]
      covNA<-names(which(apply(tmp_df,1,function(x) sum(is.na(x)))>0))
      w<-match(covNA,rownames(tmp_df)); rm(covNA)
      if (length(w)>0) { tmp_df<-tmp_df[-w,] }
      tmp_df$predicted=predict(lm)
      tmp_df$residuals=residuals(lm)
      tmp_df$predicted2=predict(lm2)
      tmp_df$residuals2=residuals(lm2)
      #---------- 
      # combine residuals with data_s
      data_s[,paste(c("residuals_","residuals2_"),p,sep="")]<-NA
      if (length(w)>0) { 
        data_s[-w,paste("residuals_",p,sep="")]<-tmp_df$residuals
        data_s[-w,paste("residuals2_",p,sep="")]<-tmp_df$residuals2
      } else { 
        data_s[,paste("residuals_",p,sep="")]<-tmp_df$residuals 
        data_s[,paste("residuals2_",p,sep="")]<-tmp_df$residuals2
      }
      rm(w)
      #----------
      # skip these plots, they take forever...
      # or try to make them better, i.e. not plot every single dot!
      #----------
      # some more plots, save them as png instead of pdf, otherwise too many points
      if (class(data_s[,p])!="factor") {
        print('More model plots\n  ')
        # plot phenotype in relation to binary covariates, boxplots
        print(paste('Plot ',p,' in relation to binary covariates, boxplots\n  ',sep=""))
        for (cov in bin_covs2) {
          if (class(data_s[,p])=="numeric"| class(data_s[,p])=="integer") {
            plot<-ggplot(data=data_s) + geom_boxplot(aes_string(x=cov,y=p),varwidth = TRUE) + theme_bw()
          } else if (class(data_s[,p])=="factor"){
            plot<-ggplot(data=subset(data_s[-which(is.na(data_s[,p])),],!is.na(sex))) +
              geom_bar(aes_string(x=cov,fill=p),position="fill") +
              scale_fill_grey(start=.3,end=.7) + theme_bw()
          }
          assign(paste(p,cov,"_plot",sep="_"),plot) ; rm(plot)
        }
        # pander(anova(lm),round=100)
        ncol=round(length(bin_covs2)/2)
        list_p<-lapply(ls(pattern=paste(p,"*.*plot",sep="")),get)
        g<-grid.arrange(arrangeGrob(grobs=list_p))
        # ggsave(file=paste(out_lmplots_dir,p,"_plot_binary_covs",covs_name,".pdf",sep=""),g,width=10,height=10)
        ggsave(file=paste(out_lmplots_dir,p,"_plot_binary_covs",covs_name,".png",sep=""),g,width=10,height=10)
        rm(cov,ncol,list_p,g)
        rm(list=ls(pattern="_plot"))
        
        ## skipping this: takes too much time and it's not very informative
        # # plot phenotype in relation to quantitative covariates, scatterplot
        # print(paste('Plot ',p,' in relation to quantitative covariates\nScatterplot with lm and loess smooth',sep=""))
        # for (cov in quant_covs2) {
        #   if (class(data_s[,p])=="numeric"| class(data_s[,p])=="integer") {
        #     plot<-ggplot(data=data_s, aes_string(x=cov,y=p)) + 
        #       stat_binhex(bins=nrow(data_s)/250,alpha=0.8) +
        #       # geom_point() +
        #       geom_smooth(span=0.3,method='loess',color="black",linetype="dashed",se=FALSE) + 
        #       geom_smooth(span=0.3,method='lm',colour='red') +
        #       scale_fill_gradient(low="#AAAAFF",high="#000080") +
        #       theme_bw()
        #   } else if (class(data_s[,p])=="factor"){
        #     plot<-ggplot(data=data_s, aes_string(x=p,y=cov)) +
        #       geom_boxplot(varwidth = TRUE) +
        #       theme_bw()
        #   }
        #   assign(paste(p,cov,"plot",sep="_"),plot) ; rm(plot)
        # }
        # ncol=round(length(quant_covs2)/2)
        # list_p<-lapply(paste(p,quant_covs2,"plot",sep="_"),get)
        # g<-grid.arrange(arrangeGrob(grobs=list_p))
        # ggsave(file=paste(out_lmplots_dir,p,"_plot_quant_covs",covs_name,".pdf",sep=""),g,width=10,height=10)
        # # ggsave(file=paste(out_lmplots_dir,p,"_plot_quant_covs",covs_name,".png",sep=""),g,width=10,height=10)
        # rm(cov,ncol,list_p,g)
        # rm(list=ls(pattern="_plot"))
        
        rm(list=ls(pattern="^p_"))
        if (class(data_s[,p])=="numeric"| class(data_s[,p])=="integer") {
          # plot histogram of original phenotype and residuals
          p_hist<-ggplot(data=tmp_df, aes_string(x=p)) + geom_histogram() + geom_vline(xintercept = mean(tmp_df[,p])-sd(tmp_df[,p])*t) + geom_vline(xintercept = mean(tmp_df[,p])+sd(tmp_df[,p])*t) + theme_bw()
          p_rhist<-ggplot(data=tmp_df, aes_string(x="residuals")) + geom_histogram() + geom_vline(xintercept = mean(tmp_df[,p])-sd(tmp_df$residuals)*t) + geom_vline(xintercept = mean(tmp_df[,p])+sd(tmp_df$residuals)*t) + theme_bw()
        } else if (class(data_s[,p])=="factor") {
          p_hist<-ggplot(data=tmp_df, aes_string(x=p)) + geom_bar() + theme_bw()
          p_rhist<-ggplot(data=tmp_df, aes_string(x="residuals")) + geom_histogram() + theme_bw()
        }
        list_p<-lapply(ls(pattern="^p_"),get)
        g<-grid.arrange(arrangeGrob(grobs=list_p))
        ggsave(file=paste(out_lmplots_dir,p,"_histograms",covs_name,".pdf",sep=""),g,width=10,height=10)
        rm(list=ls(pattern="^p_"))
        # clean
        rm(lm,tmp_df,file,lm_file)
      }
    }
  }
}
#' Save raw phenotypes with covariates + list of phenotypes
#' note: these raw phenotypes have been cleaned for outliers (i.e. blanked to NA)

# save volumes
f<-paste(out_dir,paste("/summary_phenotypes/",release,"_Cogn_all",covs_name,"Phenotypes_residuals_wHeader.table",sep=""),sep="")
write.table(data_s,f,sep="\t",quote=FALSE,row.names = FALSE,col.names=TRUE)
write.table(colnames(data_s),gsub("_wHeader.table",".list",f),col.names = FALSE,row.names = FALSE,quote=FALSE)



# clean all
# rm(data_s)
# save image
save.image(paste(out_dir,release,"_LM_residuals_Cogn_all",covs_name,".RData",sep=""))

