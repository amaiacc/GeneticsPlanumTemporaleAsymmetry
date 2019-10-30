#' ---
#' title: ukb25465_ukb25468 combined - Runs LM/plots/checks/ and saves residuals
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

#' NOTE: it does not include biological covariates (totalBV, height and BMI) that had been included previously.


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

release="ukb25465_ukb25468"

# Define directories, dependent on  system
if (Sys.info()['sysname']=='Windows') {dir="P://workspaces/"} else {dir="/data/workspaces/lag/workspaces/"}
primary_dir=paste(dir,"lg-ukbiobank/primary_data/",sep="")
working_dir=paste(dir,"lg-ukbiobank/working_data/amaia/",sep="")
out_dir=paste(working_dir,"/pheno_files/genetic_v2/",release,"/",sep="")
out_lmplots_dir=paste(working_dir,"/pheno_files/genetic_v2/",release,"/lm_plots/",sep="")
out_lm_dir=paste(working_dir,"/pheno_files/genetic_v2/",release,"/lm/",sep="")
# create ouput directories if they don't already exist
dir.create(file.path(out_dir), showWarnings = FALSE)
dir.create(file.path(out_lmplots_dir), showWarnings = FALSE)
dir.create(file.path(out_lm_dir), showWarnings = FALSE)

# Set working directory
opts_knit$set(root.dir = working_dir)


#' Start!
setwd(working_dir)

#' check if RData/data file exists, load, otherwise print warning
rdata<-paste(out_dir,release,"_extractPhenotypes_imaging_extracted_PhenoCovs.RData",sep="")
# fdata<-paste(out_dir,release,"_imaging_covsImagingPhenotypes_wHeader.table",sep="")
if (file.exists(rdata)) {
load(rdata) 
  } else {
  cat("Warning!!\nThe data you tried to load does not exist.\nYou may have to run ukb21288_ukb21293_extractPhenotypes_imagingSubset.R")
  stop()
}


#--------------------------------------
#' ## Regress out covariates from all the phenotypes, and save residuals. 
#--------------------------------------
#' Define order of covariates to include within lm
#' Drop batch, too many levels # ,"batch"
#' Drop assessement center as well, not relevant for imaging phenos, all assessment took place in the same place
maxPC=10 # or 40...
covs_order<-c("age","zage2","sex",paste("PC",1:maxPC,sep=""),"array",scanner_position ) #, "totalBV","height","BMI") "assessment_centre", 
covs_name="_noBioCovs_noAssessmentC"
# if we want to include height and BMI as biological covariates, just in case
# covs_order<-c("assessment_centre","age","zage2","sex",paste("PC",1:10,sep=""),"array",scanner_position, "totalBV","height","BMI")
# covs_name=""

quant_covs2<-quant_covs2[quant_covs2 %in% covs_order]
bin_covs2<-bin_covs2[bin_covs2 %in% covs_order]

#' Define phenotypes of interest, to run LM
phenos3_AI<-phenos2[grep("AI_",phenos2)] # select only AIs
phenos3_AIfa<-phenos2[grep("AIfa_",phenos2)] # select only AIs
phenos3_nonAI<-phenos2[grep("AI",phenos2,invert=TRUE)] # 851, select only non AIs; i.e. left/right or non-symmetric
phenos3_nonAI<-phenos3_nonAI[-(1:4)]
#----------------------------------------------------------
# all phenotypes to run lms and get residuals from:
# 72 volumes + 1 total volume = 73
phenos3_AI_1<-phenos3_AI[grep("VOLUME|Volume",phenos3_AI)]
phenos3_nonAI_1<-phenos3_nonAI[grep("VOLUME|Volume",phenos3_nonAI)]
# 189 DTI (FA skeleton)
phenos3_AI_2<-phenos3_AI[grep("FA_skeleton",phenos3_AI)]
phenos3_nonAI_2<-phenos3_nonAI[grep("FA_skeleton",phenos3_nonAI)]
# 108 DTI (tracts)
phenos3_AI_3<-phenos3_AI[grep("tract_",phenos3_AI)][grep("FA_skeleton",phenos3_AI[grep("tract_",phenos3_AI)],invert=TRUE)]
phenos3_nonAI_3<-phenos3_nonAI[grep("tract_",phenos3_nonAI)][grep("FA_skeleton",phenos3_nonAI[grep("tract_",phenos3_nonAI)],invert=TRUE)]
#----------------------------------------------------------
# threshold to define outliers --> 4 sd deviations from mean
t=4
length(c(phenos3_AI,phenos3_nonAI)) # 1218 measures (T1 volume + dti (FA skeleton and tracts)); AI,left, right)

pheno_names<-c("Volume","FAskeleton","tracts")
# run three rounds, one per phenotype type:
for (i in 2:3) {
  phen_round<-paste(c("phenos3_AI","phenos3_nonAI"),i,sep="_")
  for (p in c(get(phen_round[1]),get(phen_round[2]))) {
    print(paste("Dependent variable: ",p,"\n  ",sep=""))
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
      data_s[out,p]<-NA
      rm(m,sd,out)
      #----------------------
      #print('Regress out covariates and save residuals"\n  ')
      print(paste('Run a linear model to get residuals for ',p,'\n  ',sep=""))
      # model_formula<-paste(p,"~",paste(c(bin_covs2,quant_covs2),collapse=" + "))
      model_formula<-paste(p,"~",paste(c(covs_order),collapse=" + "))
      if (class(data_s[,p])=="numeric" | class(data_s[,p])=="integer") {
        lm<-do.call("lm",list (as.formula(model_formula),data=data_s))
      } else if (class(data_s[,p])=="factor") {
        lm<-do.call("glm",list (as.formula(model_formula),data=data_s,family=binomial(link='logit')))
      }
      pander(anova(lm),round=1000)
      #----------
      print('Fit check\n  ')
      # par(mfrow=c(2,2))
      pdf(file=paste(out_lmplots_dir,p,"_lm_plots",covs_name,".pdf",sep=""))
      plot(lm)
      dev.off()
      par(mfrow=c(1,1))
      #----------
      # save model + summary model estimates in txt file as well
      df<-as.data.frame(anova(lm)); df$variable<-rownames(df)
      df<-df[,c(NCOL(df),1:(NCOL(df)-1))]
      df2<-as.data.frame(summary(lm)$coefficients); df2$variable<-rownames(df2)
      df2<-df2[,c(NCOL(df2),1:(NCOL(df2)-1))]
      write.table(df,file=paste(out_lm_dir,p,"_lm",covs_name,".txt",sep=""),sep="\t",row.names=FALSE,quote=FALSE)
      write.table(df2,file=paste(out_lm_dir,p,"_coefficients_lm",covs_name,".txt",sep=""),sep="\t",row.names=FALSE,quote=FALSE)
      rm(df,df2)
      #----------
      tmp_df<-data_s[,c("ID1","ID2",p,bin_covs2,quant_covs2)]
      covNA<-names(which(apply(tmp_df,1,function(x) sum(is.na(x)))>0))
      w<-match(covNA,rownames(tmp_df)); rm(covNA)
      if (length(w)>0) { tmp_df<-tmp_df[-w,] }
      tmp_df$predicted=predict(lm)
      tmp_df$residuals=residuals(lm)
      #---------- 
      m<-mean(tmp_df$residuals,na.rm = TRUE)
      sd<-sd(tmp_df$residuals,na.rm = TRUE)
      res_out<-which(tmp_df$residuals > m+t*sd | tmp_df$residuals < m-t*sd )
      print('Number of outliers for residuals, deviating more than 4SD from mean:')
      print(length(res_out))
      if (length(res_out) > 0) {
        tmp_df$residuals[res_out]<-NA
      }
      rm(m,sd,res_out)
      # combine residuals with data_s
      data_s[,paste("residuals_",p,sep="")]<-NA
      if (length(w)>0) { data_s[-w,paste("residuals_",p,sep="")]<-tmp_df$residuals
      } else { data_s[,paste("residuals_",p,sep="")]<-tmp_df$residuals }
      rm(w)
      #----------
      # some more plots
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
        ggsave(file=paste(out_lmplots_dir,p,"_plot_binary_covs",covs_name,".pdf",sep=""),g,width=10,height=10)
        rm(cov,ncol,list_p,g)
        rm(list=ls(pattern="_plot"))
        # plot phenotype in relation to quantitative covariates, scatterplot
        print(paste('Plot ',p,' in relation to quantitative covariates\nScatterplot with lm and loess smooth',sep=""))
        for (cov in quant_covs2) {
          if (class(data_s[,p])=="numeric"| class(data_s[,p])=="integer") {
            plot<-ggplot(data=data_s, aes_string(x=cov,y=p)) + geom_point() +
              geom_smooth(span=0.3,method='loess') + geom_smooth(span=0.3,method='lm',colour='red') +
              theme_bw()
          } else if (class(data_s[,p])=="factor"){
            plot<-ggplot(data=data_s, aes_string(x=p,y=cov)) +
              geom_boxplot(varwidth = TRUE) +
              theme_bw()
          }
          assign(paste(p,cov,"plot",sep="_"),plot) ; rm(plot)
        }
        ncol=round(length(quant_covs2)/2)
        list_p<-lapply(paste(p,quant_covs2,"plot",sep="_"),get)
        g<-grid.arrange(arrangeGrob(grobs=list_p))
        ggsave(file=paste(out_lmplots_dir,p,"_plot_quant_covs",covs_name,".pdf",sep=""),g,width=10,height=10)
        rm(cov,ncol,list_p,g)
        rm(list=ls(pattern="_plot"))
        
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
        ## skip this, no need to save 100% extra files...
        # print('Get residuals and save for downstream analyses\n')
        # file=paste(working_dir,"/pheno_files/genetic_v2/ukb21288_ukb21293_imaging",covs_name,"_residuals.",p,"_wHeader.table",sep="")
        # write.table(tmp_df[,c("ID1","ID2","residuals")],file=file,sep="\t",quote=FALSE,row.names = FALSE,col.names=TRUE)
        # clean
        rm(lm,tmp_df,file,lm_file)
      }
    }
  }
  # save list of phenotypes
  # AI
  write.table(get(paste("phenos3_AI",i,sep="_")),paste("pheno_files/genetic_v2/",release,"_imaging",covs_name,"_AI_",pheno_names[i],"Phenotypes.list",sep=""),col.names = FALSE,row.names = FALSE,quote=FALSE)
  ## non AI
  write.table(get(paste("phenos3_nonAI",i,sep="_")),paste("pheno_files/genetic_v2/",release,"_imaging",covs_name,"_nonAI_",pheno_names[i],"Phenotypes.list",sep=""),col.names = FALSE,row.names = FALSE,quote=FALSE)
  # clean
  rm(phen_round)
}

#' Save raw phenotypes with covariates + list of phenotypes
#' note: these raw phenotypes have been cleaned for outliers (i.e. blanked to NA)

# define columns: general (non imaging covariates/factors), volumes (raw + residualized), dti (raw + residualized)
other_cols<-colnames(data_s)[grep("FA_skeleton|tract|VOLUME|Volume",invert=TRUE,colnames(data_s))]
volume_cols_raw_res<-colnames(data_s)[grep("VOLUME|Volume",colnames(data_s))]
dti_cols_raw_res<-colnames(data_s)[grep("FA_skeleton",colnames(data_s))]
dti_tracts_cols_raw_res<-colnames(data_s)[grep("_tract_",colnames(data_s))][grep("FA_skeleton",colnames(data_s)[grep("_tract_",colnames(data_s))],invert=TRUE)]

# # save volumes
# data_s_vols<-data_s[,c(other_cols,volume_cols_raw_res)] # 656 variables
# f_vols<-paste(out_dir,paste("/summary_phenotypes/",release,"_imaging",covs_name,"VolumePhenotypes_residuals_wHeader.table",sep=""),sep="")
# write.table(data_s_vols,f_vols,sep="\t",quote=FALSE,row.names = FALSE,col.names=TRUE)
# write.table(colnames(data_s_vols),gsub("_wHeader.table",".list",f_vols),col.names = FALSE,row.names = FALSE,quote=FALSE)
# rm(f_vols,data_s_vols)

# save dti
data_s_dti<-data_s[,c(other_cols,dti_cols_raw_res)] # 1655 variables
f_dti_fa<-paste(out_dir,paste("/summary_phenotypes/",release,"_imaging",covs_name,"FAskeletonPhenotypes_residuals_wHeader.table",sep=""),sep="")
write.table(data_s_dti,f_dti_fa,sep="\t",quote=FALSE,row.names = FALSE,col.names=TRUE)
write.table(colnames(data_s_dti),gsub("_wHeader.table",".list",f_dti_fa),col.names = FALSE,row.names = FALSE,quote=FALSE)
rm(f_dti_fa,data_s_dti)

# save tracts
data_s_dti_tracts<-data_s[,c(other_cols,dti_tracts_cols_raw_res)] # 1324 variables
f_dti_tracts<-paste(out_dir,paste("/summary_phenotypes/",release,"_imaging",covs_name,"tractsPhenotypes_residuals_wHeader.table",sep=""),sep="")
write.table(data_s_dti_tracts,f_dti_tracts,sep="\t",quote=FALSE,row.names = FALSE,col.names=TRUE)
write.table(colnames(data_s_dti_tracts),gsub("_wHeader.table",".list",f_dti_tracts),col.names = FALSE,row.names = FALSE,quote=FALSE)
rm(f_dti_tracts,data_s_dti_tracts)

# clean all
# rm(data_s)
# save image
save.image(paste(out_dir,release,"_LM_residuals_imagingSubset",covs_name,".RData",sep=""))

