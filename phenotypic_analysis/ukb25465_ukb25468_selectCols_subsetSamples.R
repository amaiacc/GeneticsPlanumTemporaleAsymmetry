#' ---
#' title: 'ukb25465_ukb25468_sqc - select columns for variables of interest and subset sample (to imaging subset), also save all'
#' author: amacar
#' date: "Created: 2018-04-06. Modified: `r Sys.time()`"
#' output:
#'   html_document:
#'     toc: true
#'     toc_depth: 2
#'     theme: "flatly"
#'     highlight: "textmate"
#'   pdf_document:
#'     keep_tex: true
#' ---
library(knitr)
opts_chunk$set(include=TRUE, echo=TRUE, error = FALSE, warning = FALSE, results = "asis", tidy=TRUE, width=50, fig.width=8, fig.height=6)
# opts_chunk$set(dev="pdf", 
#                dev.args=list(type="cairo"),
#                dpi=96)

# Set working directory
# Set directories, dependent on  system
if (Sys.info()['sysname']=='Windows') {dir="P://workspaces/"} else {dir="/data/workspaces/lag/workspaces/"}
working_dir=paste(dir,"lg-ukbiobank/working_data/amaia",sep="")
primary_dir=paste(dir,"lg-ukbiobank/primary_data",sep="")
#setwd(primary_dir)
# 
#' abstract: " Subset to variables and samples of interest."
opts_knit$set(root.dir = working_dir)

library(rmarkdown)
library(pander)
library(Cairo)
library(ggplot2)
require(grid)
library(gridExtra)
if (!require("RColorBrewer")) {
  install.packages("RColorBrewer")
  library(RColorBrewer)
}
if (!require("gridExtra")) {
  install.packages("gridExtra")
  library(gridExtra)
}

if (!require("pander")) {
  install.packages("pander")
  library(pander)
}

#+ setup, include=TRUE
# echo=FALSE, warning=FALSE, message=FALSE, 

#options(bitmapType='cairo')
options(stringsAsFactors = FALSE)
#options(bitmapType='cairo',device=cairo_pdf,stringsAsFactors = FALSE)
panderOptions('table.split.table', Inf)
panderOptions('digits', 2)
panderOptions('round', 2)
panderOptions('keep.trailing.zeros', TRUE)
panderOptions('table.split.cells', 10)

# Some functions
#extract legend
# http://stackoverflow.com/questions/13649473/add-a-common-legend-for-combined-ggplots
grid_arrange_shared_legend <- function(...) {
  plots <- list(...)
  g <- ggplotGrob(plots[[1]] + theme(legend.position="bottom"))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  grid.arrange(
    do.call(arrangeGrob, lapply(plots, function(x)
      x + theme(legend.position="none"))),
    legend,
    ncol = length(plots),
    heights = unit.c(unit(1, "npc") - lheight, lheight))
}

# Select subset of variables to work on them

# small function to detect which are the column index of variables of interest - but not all are included, since not all are part of field.tsv
get_column_index<-function(string, data=field, data2_cols=colnames(bd_sqc_fam), info="info2", type="type"){
  # get field name in lower case, to match regardless of case
  string2=tolower(string)
  w<-grep(string2,tolower(data[,info]))
  #
  data_code<-data[,type][w]
  data_info<-data[w,c("info1","info2")]
  data_code_f<-paste("^f\\.",data_code,"\\.",sep="")
  data_cols_index<-unlist(lapply(data_code_f, function(x)
    if (length(grep(x,data2_cols) )> 0) { w<-paste(grep(x,data2_cols),collapse=";")} else { w<-NA}
  ))
  data_colname<-unlist( lapply(data_cols_index, function(x)
    if (is.na(x))  { w2<-NA} else  { w2<- paste(data2_cols[as.numeric(unlist(strsplit(x,";")))] ,collapse=";") }   ))   #    data2_cols[as.numeric(data_cols_index)]
  out<-cbind(data_info,data_code,data_cols_index,data_colname)
  return(out)
}

# Summary statistics and plots - covariates 
# to get summary stats and basic plot per variable
summary_func<-function(label,data2=data,vars2=vars){
  w<-vars2$data_colname[vars2$info2==label]
  if (is.numeric(data[,w])==TRUE) {  x<-summary(data[,w])} else { table(data[,w])}
  return(x)
}

summary_per_col<-function(data=data,column="") {
  x<-data[,c]
  y<-subset(data,imaging==1)[,c]
  w<-  vars$info2[match(c,vars$data_colname)]
  if (length(w)>0) {
    print(w)
    if (class(x)=="integer"| class(x)=="numeric") {
      print("Total sample")
      print(summary(x))
      print(sd(x),na.rm=TRUE)
      
    }
    if (length(grep("factor",class(x))>0)) {
      print(table(x)  ) / sum(!is.na(x))
    }
  } }

## Start!

# List of available fields
field<-read.csv(paste(primary_dir,"ukb_field.tsv",sep="/"),sep="\t",header=FALSE,stringsAsFactors=FALSE)
colnames(field)<-c("V1","type","info1","info2")

#' Load data with all phenotypes and sample QC:
load(file=paste(working_dir,"/ukb25465_ukb25468_sqc_fam.RData",sep=""))
# define working dir again, because it's different to the one in the RData
working_dir=paste(dir,"lg-ukbiobank/working_data/amaia",sep="")
samples2exclude<-read.csv(paste(working_dir,"/w16066_20180503.csv",sep=""),header=FALSE) #
#' Get list of samples included in the imaging analyses. These samples contain the three imaging data that we need (T1, T2_FLAIR, rsMRI).
imaging_samples<-read.table(paste(dir,"lg-ukbiobank/working_data/imaging_data/subject_to_process/pruned_list_subjects_to_process_FINAL.txt",sep=""),stringsAsFactors = FALSE)
imaging_samples$imaging<-1

#' Combine phenotypic variables, samples2exclude and imaging sample info.
# exclusion
bd_sqc_fam$exclude<-0
w<-which(bd_sqc_fam$f.eid %in% samples2exclude$V1)
bd_sqc_fam$exclude[w]<-1
# imaging, these only flag imaging samples that have the three imaging modalities (T1, T2-flair, rfMRI)
bd_sqc_fam$imaging<-0
# check, just for one
w<-which(bd_sqc_fam$f.eid %in% imaging_samples$V1) # index of samples with imaging data within bd
#' Flag imaging ok samples (i.e. T1, T2-FLAIR and fMRI data available)
if (length(w)>0) {bd_sqc_fam$imaging[w]<-1 }
#' Are all samples flagged as imaging, listed within the list of imaging samples?
sum(bd_sqc_fam$f.eid[bd_sqc_fam$imaging==1] %in% imaging_samples$V1) ==NROW(imaging_samples)
rm(w)

#' Make a second imaging_T1 flag: for all samples that have a non-missing value for data field f.25000
#' f.25000 = Volumetric scaling from T1 head image to standard space
w<-which(!is.na(bd_sqc_fam$f.25000.2.0)) # index of samples with imaging data within bd
bd_sqc_fam$imaging_T1<-0
bd_sqc_fam$imaging_T1[w]<-1
rm(w)


#' Some summary table of counts: total dataset, imaging subset and genetic subset
t<- table(bd_sqc_fam$imaging,bd_sqc_fam$in.white.British.ancestry.subset)
rownames(t)<-c("Imaging - NO", "Imaging - YES")
colnames(t)<-c("white.British - NO", "white.British - YES")
pander(t); rm(t)

t<- table(bd_sqc_fam$imaging_T1,bd_sqc_fam$in.white.British.ancestry.subset)
rownames(t)<-c("Imaging T1 - NO", "Imaging T1 - YES")
colnames(t)<-c("white.British - NO", "white.British - YES")
pander(t); rm(t)


#' Summary table of counts: total dataset, imaging subset and genetic subset
t<- table(bd_sqc_fam$imaging,is.na(bd_sqc_fam$PC1))
rownames(t)<-c("Imaging - NO", "Imaging - YES")
colnames(t)<-c("Genetic PC1 - NO", "Genetic PC1 - YES")
pander(t); rm(t)

#' Summary table of counts: total dataset, imaging subset (T1) and genetic subset
t<- table(bd_sqc_fam$imaging_T1,!is.na(bd_sqc_fam$PC1))
rownames(t)<-c("Imaging T1 - NO", "Imaging T1 - YES")
colnames(t)<-c("Genetic PC1 - NO", "Genetic PC1 - YES")
pander(t); rm(t)

t<- table(bd_sqc_fam$imaging_T1,bd_sqc_fam$in.white.British.ancestry.subset)
rownames(t)<-c("Imaging T1 - NO", "Imaging T1 - YES")
colnames(t)<-c("Genetic PC1 - NO", "Genetic PC1 - YES")
pander(t); rm(t)

#' Summary table of counts: imaging subset (T1) vs imaging subset (T1,T2-Flair,rfMRI)
t<- table(bd_sqc_fam$imaging,bd_sqc_fam$imaging_T1)
rownames(t)<-c("Imaging  - NO", "Imaging - YES")
colnames(t)<-c("Imaging (T1) - NO", "Imaging (T1) - YES")
pander(t); rm(t)

#' Identify variables of interest, including: baseline characteristics, age, sex, handedness, early 

# Baseline characteristics
var_baseline<-get_column_index(string="Baseline",info="info1") # Baseline characteristics
var_age<-get_column_index(string="Age at recruitment|Age when attended assessment centre|Year of birth")
var_date<-get_column_index(string="^Date |^Month|^Year")
var_sex<-get_column_index(string="^Sex|^Gender") # Not a perfect search, not all may be related
var_hand<-get_column_index(string="hand |handed")
# var_birth=get_column_index(string="birth")
var_size=get_column_index(string="^Body mass|^Weight| height$|")
# early life factors
var_early<-get_column_index(string="Early",info="info1") #   Early life factors
# Genotype related variables
var_genomics<-get_column_index(string="Genomics",info="info1")
var_genotype<-get_column_index(string="^Genotype|genotype quality",info="info2")
var_UKBiL<-get_column_index(string="UKBiL")
var_UKBassessment<-get_column_index(string="UK Biobank",info="info2")
# Ethnicity
var_ethnic<-get_column_index(string="ethnic")
# Phenotypes
var_t1<-get_column_index(string="T1",info="info1") # Not a perfect search, not all may be related - not all are included, e.g. grey matter volumes are not contained in ukb_field.tsv
var_dti<-get_column_index(string="dMRI",info="info1") # Not a perfect search, not all may be related
var_head<-get_column_index(string="head")
var_brain<-get_column_index(string="brain") # Not a perfect search, not all may be related
# Illnesses - neurological
var_mental<-get_column_index(string="^Mental| mental |psychiatric") # Not a perfect search, not all may be related
var_mental2<-get_column_index(string="psychiatric",info="info1") # Not a perfect search, not all may be related
var_diag<-get_column_index(string="diagnos",info="info1") # Not a perfect search, not all may be related
var_illness<-get_column_index(string="illness") # Not a perfect search, not all may be related
var_medical<-get_column_index(string="medical")
## cognitive fields: 20018, 10137, 20016, 20023, 4282, 10146, 10144, 10612, 10610
# field[match(as.character(c(20018, 10137, 20016, 20023, 4282, 10146, 10144, 10612, 10610)),field$type),]
var_cognitive<-get_column_index(string="Prospective memory result|Number of incorrect matches|Fluid intelligence score|Mean time to correctly identify matches|Maximum digits remembered correctly|Pattern of lights as remembered|Time taken to complete lights test|Number of words beginning with 'S'|Word count|Time to complete round",info="info2") #   
# other
var_socio<-get_column_index(string="socio",info="info1") #  
var_educ<-get_column_index(string=" educ|^educ",info="info1") #   
var_hearing<-get_column_index(string="hearing$",info="info1") #   
var_smoking<-get_column_index(string="smoking status")
# var_lang<-get_column_index(string="speech|language",info="info2") #   


# Define variables of interest in a list
## also icnluding gemographic variables from Antonietta's script (extract_demographic_info_new_new.sh)
## which are:
vars_demo_antonietta<-c("f.eid", ## SUBJECT.ID
                        ## HANDEDNESS
                        unlist(strsplit(var_hand$data_colname[var_hand$data_code!="38"],";")), # handedness
                        ## GENDER
                        unlist(strsplit(var_baseline$data_colname[var_baseline$data_code=="31"],";")), # gender
                        unlist(strsplit(var_genomics$data_colname[var_genomics$data_code=="22001"],";")), # genetic sex
                        ## YEAR AND MONTH OF BIRTH
                        unlist(strsplit(var_baseline$data_colname[var_baseline$data_code=="34"|var_baseline$data_code=="52"],";")), # year and month of birth
                        ## AGE
                        unlist(strsplit(var_age$data_colname[var_age$data_code=="21003"],";")), # age at asessment
                        ## DATE
                        unlist(strsplit(var_date$data_colname[var_date$data_code=="53"],";")), # date of asessment
                        ## HEAD MOTION
                        unlist(strsplit(var_head$data_colname[var_head$data_code=="25741"],";")), 
                        ## HEARING TESTS
                        unlist(strsplit(var_hearing$data_colname[var_hearing$data_code=="2247"],";")), 
                        ## COGNITIVE SCORES
                        unlist(strsplit(var_cognitive$data_colname[!is.na(var_cognitive$data_colname)],";")), 
                        ## SOCIODEMOGRAPHIC SCORES AND EDUCATION LEVEL
                        unlist(strsplit(var_educ$data_colname[!is.na(var_educ$data_colname)],";")), 
                        ## SMOKING STATUS
                        var_smoking$data_colname, # not available, it's NA
                        ## ETHNICITY AND COUNTRY OF BIRTH
                        unlist(strsplit(var_ethnic$data_colname[!is.na(var_ethnic$data_colname)],";")), # ethnicity
                        unlist(strsplit(var_early$data_colname[var_early$data_code=="20115"|var_early$data_code=="1647"|var_early$data_code=="130"|var_early$data_code=="129"],";")),  # country of birth
                        ## EARLY LIFE FACTORS
                        unlist(strsplit(var_early$data_colname[var_early$data_code=="1777"|var_early$data_code=="1787"|var_early$data_code=="1677"|var_early$data_code=="20022"],";")),  # country of birth
                        ## Illnesses - neurological
                        unlist(strsplit(var_mental$data_colname[!is.na(var_mental$data_colname)],";")), # mental categories
                        unlist(strsplit(var_illness$data_colname[var_illness$data_code=="2188"|var_illness$data_code=="20002"],";")), # long-standing illness, disability or infirmity / self-reported (some are not available)
                        unlist(strsplit(var_illness$data_colname[var_diag$data_code=="41202"|var_diag$data_code=="41204"],";")), # ICD10 main and secondary
                        unlist(strsplit(var_illness$data_colname[var_diag$data_code=="41203"|var_diag$data_code=="41205"],";")), # ICD9 main and secondary                         
                        ## VOLUME OF GRAY MATTER AND WHITE MATTER
                        unlist(strsplit(var_t1$data_colname[var_t1$data_code=="25006"|
                                                                 var_t1$data_code=="25005"|
                                                                 var_t1$data_code=="25009"|
                                                                 var_t1$data_code=="25010"|
                                                                 var_t1$data_code=="25025"],";")),
                        
                        # all vols from 2578-25920, not present within the fields - get colnames
                        colnames(bd_sqc_fam)[grep(paste(c(25782:25920),collapse="|"),colnames(bd_sqc_fam))]
                        )
table(is.na(vars_demo_antonietta)) # four variables are missing: smoking, some non-cancer illness codes
vars_demo_antonietta<-vars_demo_antonietta[!is.na(vars_demo_antonietta)] 
#' Demographic variables selected from Antonietta's list: 447 column IDs total
#' ls /data/workspaces/lag/workspaces/imaging_data/demographic/*9246_FINAL.txt | wc -l # 193 files ## where does this mismatch come from?

#' Other biological variables, to be used as imaging covariates, from ukb10785:
var_blood<-get_column_index(string="blood pressure")
# BMD: bone mineral density
var_bmd<-get_column_index(string="BMD")

# all variables of interest - combine into a dataframe
vars<-data.frame()
for (x in objects()[grep("var_",objects())] ) {
  tmp<-get(x);tmp$var_type<-x
  vars<-rbind(vars,tmp); rm(tmp)
  } 
rm(x)
vars<-subset(vars,!is.na(data_colname))
# clean if duplicated
w<-which(duplicated(vars[-6])) # col 6 is "var_type" which could make the cols nonduplciated if the same variable has been accessed several times
if (length(w)>0) {
vars<-vars[-which(duplicated(vars[-6])),]
}
#' Create new colnames for these variables, i.e. by using info2, and replacing spaces with underscore
vars$colnames_new<-gsub(" ","_",vars$info2)
vars$var_within_antonietta<-NA
# flag variables that are within Antonietta's list
vars$var_within_antonietta[grep(paste(vars_demo_antonietta,collapse="|"),vars$data_colname)]<-1
# clean previous tables, now the var_type info is coded within the vars datframe, so should be able to subset
# rm(list = objects()[grep("var_",objects())] )

#' These are not within field (should download latest field file), so just name them...
# f.25756-25759: brain position variables
# f.25742: tfMRI head motion
vars_scanner<-paste("f.",c(25756:25759,25742),".2.0",sep="")

#' Summary table of all selected variables: 125 data types
length(unlist(strsplit(vars$data_colname,";"))) #'2014 columns total
cols_total<-unique(unlist(strsplit(vars$data_colname,";")))
#' This does not include all brain measurements, so combine with vars selected by Antonietta:
cols_total2<-unique(c(cols_total,vars_demo_antonietta,vars_scanner)) #' 1446
#' Other variables of interest, from the sqc and fam, and variables that I created (e.g. imaging)
cols_sqc_fam<-colnames(bd_sqc_fam)[grep("^f\\.",colnames(bd_sqc_fam),invert=TRUE)]

#' Total number of variables:
cols_all<-c(cols_total2,cols_sqc_fam)
length(cols_all)
#' 2233, including dti # prev round: 1471

#' For each data type that has more than one column --> combine across columns by:
#' 1) for numeric values --> compute mean across variables
#' 2) for categorical values --> collapse across each instance; doesn't make sense always but...
#' 
# note, only to variables coded as f.xxx
f.vars<-cols_all[grep("^f",cols_all)]
f.vars_per_col_table<-table( sapply(strsplit(gsub("f.","",f.vars), "\\."),"[[",1) )
f.vars_across_cols<-paste("f.",names(which(f.vars_per_col_table>1)),sep="")
# 
bd_sqc_fam[,f.vars_across_cols]<-sapply(f.vars_across_cols,function(x) {
    tmp<-bd_sqc_fam[,grep(paste(x,".",sep=""),colnames(bd_sqc_fam),fixed=TRUE)]
    cols<-colnames(tmp)
    var_class<-unique(apply(tmp,2,class))
    # create new variable - empty
    tmp[,x]<-NA
    # if variable type is integer or numeric, then compute mean
    if (var_class=="integer"|var_class=="numeric") {
        tmp[,x]<-apply(tmp[,cols],1, function(z) {
            mean(unique(z[!is.na(z)]),na.rm=TRUE)  } )
      } else if (var_class=="factor"|var_class=="character") {
        # if variable type is factor or character, collapse across unique values
        tmp[,x]<-apply(tmp[,cols],1, function(z) {
          paste(unique(z[!is.na(z)]),collapse=";")    } ) } 
  return(tmp[,x]) 
  } )


# select variables of interest
cols_all2<-c(cols_all,f.vars_across_cols)
# Subset with selected variables, columns
bd_sqc_fam_selected<-bd_sqc_fam[,cols_all2]
#' Rename colnames
#' by replacing the data_code by the content of colnames_new within vars
vars$length<-unlist(lapply(strsplit(vars$data_colname,";"),length))
vars$colname2replace<-NA
vars$colname2replace[vars$length==1]<-vars$data_colname[vars$length==1]
vars$colname2replace[vars$length>1]<-paste("f.",vars$data_code[vars$length>1],sep="")

backup_colnames<-colnames(bd_sqc_fam_selected)
new_colnames<-backup_colnames
index2replace<-match(vars$colname2replace,colnames(bd_sqc_fam_selected))
content2replace<-vars$colnames_new
new_colnames[index2replace]<-content2replace
# replace column names
colnames(bd_sqc_fam_selected)<-new_colnames

# subset imaging dataset, based on imaging_T1 (which contains all samples within imaging==1)
bd_sqc_fam_selected_imaging<-subset(bd_sqc_fam_selected,imaging_T1==1)

# clean up
rm(backup_colnames,new_colnames,index2replace,content2replace)
#' Save tables:
#' vars --> as selected_variables_summary.table
var_file<-paste(working_dir,"/demographic/","ukb25465_ukb25468_sqc_fam_selectedPhenotypes_variables.csv",sep="")
write.csv(vars,file=var_file,quote=TRUE,row.names=FALSE)
# all samples, selected phenotypes
file1=paste(working_dir,"/demographic/","ukb25465_ukb25468_sqc_fam_selectedPhenotypes_allSamples.csv",sep="")
write.csv(bd_sqc_fam_selected,file=file1,row.names=FALSE)
# imaging samples, selected phenotypes
file2=paste(working_dir,"/demographic/","ukb25465_ukb25468_sqc_fam_selectedPhenotypes_imagingSamples.csv",sep="")
write.csv(bd_sqc_fam_selected_imaging,file=file2,row.names=FALSE)
# # clean more
# rm(var_file,file1,file2)


# Save one file per variable, as Antonietta did. Make sure that this makes sense to her.
