options(stringsAsFactors = FALSE)
library(reshape)
library(ggplot2)
colrs=c("#d95f02","#1f78b4","#33a02c","#000000") # '#1b9e77','#d95f02','#7570b3 # "#a6cee3"  ,"#7570b3"
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
#----------------------------------------------------------------------
subset_name="imagingT1_N18057"
root="ukb25465_ukb25468"
pattern_run="_noBioCovs_noAssessmentC"
pattern=paste(root,"imaging",pattern_run,"_",sep="")
#----------------------------------------------------------------------
# define working_dirs
if (Sys.info()['sysname']=='Windows') {dir="P://workspaces/"} else {dir="/data/workspaces/lag/workspaces/"}
working_dir=paste(paste(dir,"lg-ukbiobank/working_data/amaia/genetic_data/PT/",subset_name,"/ldsc/partitioned/",sep=""))
#----------------------------------------------------------------------
# read data
setwd(working_dir)
files<-list.files(pattern="h2_*.*results$")
for (f in files){
  t<-read.table(f,sep="\t",header=T)
  t$file<-f
  if (f==files[1]){
    h2part<-t
  } else {h2part<-merge(h2part,t,all=TRUE)}
  rm(t)
}
rm(f,files)


#
files<-list.files(pattern="h2_*.*cell_type_results.txt")
for (f in files){
  t<-read.table(f,sep="\t",header=T)
  t$file<-f
  if (f==files[1]){
    h2celltype<-t
  } else {h2celltype<-merge(h2celltype,t,all=TRUE)}
  rm(t)
}


#----------------------------------------------------------------------
# define phenotype and baselineLD
#----------------------------------------------------------------------
table(h2part$file)
h2part$model<-gsub(".results","",sapply(strsplit(h2part$file,"_h2_"),"[[",2))
h2part$phenotype<-sapply(strsplit(h2part$file,"_h2_"),"[[",1)
h2part$sample<-"total"
h2part$sample[grep("_females",h2part$phenotype)]<-"females"
h2part$sample[grep("_males",h2part$phenotype)]<-"males"
h2part$adj<-""
h2part$adj[grep("adjTBV",h2part$phenotype)]<-"adjTBV"
# clean pheno name
h2part$pheno<-gsub("_","",gsub("_females|_males|.adjTBV","",gsub("^_","",gsub("Volume_of_grey_matter_in_|_VOLUME_|Planum_Temporale","",h2part$phenotype))))

# compute -log(Pvalue) for the enrichment, to plot
h2part$logP<-(-log10(h2part$Enrichment_p))
# flag if significant
h2part$sig<-"ns"
h2part$sig[h2part$Enrichment_p<0.05]<-"0.001<p<0.05"
h2part$sig[h2part$Enrichment_p<0.05/52]<-"p<0.001"
h2part$sig<-factor(h2part$sig,levels=c("ns","0.001<p<0.05","p<0.001"))
# clean category names
h2part$category<-gsub(".bed","",gsub("L2_0","",h2part$Category))
h2part$category<-gsub("\\.extend.500|\\.flanking.500","",h2part$category)
# flag if it's extend or not
h2part$extend500<-0
h2part$extend500[grep("extend.500|flanking.500",h2part$Category)]<-1
h2part$extend500<-as.factor(h2part$extend500)
# define category types...
h2part$category<-gsub("non_","non",h2part$category)
#
h2part$cat_source<-NA
h2part$cat_source[grep("UCSC",h2part$category)]<-"UCSC"
h2part$cat_source[grep("Hoffman",h2part$category)]<-"Hoffman"
h2part$cat_source[grep("ENCODE",h2part$category)]<-"ENCODE"
h2part$cat_source[grep("Trynka",h2part$category)]<-"Trynka"
h2part$cat_source[grep("Andersson",h2part$category)]<-"Andersson"
h2part$cat_source[grep("Vahedi",h2part$category)]<-"Vahedi"
h2part$cat_source[grep("Hnisz",h2part$category)]<-"Hnisz"
h2part$cat_source[grep("PGC2",h2part$category)]<-"PGC2"
h2part$cat_source[grep("LindbladToh",h2part$category)]<-"LindbladToh"
# v1.1 LD
h2part$cat_source[grep("GTEx",h2part$category)]<-"GTEx"
h2part$cat_source[grep("GERP",h2part$category)]<-"GERP"
h2part$cat_source[grep("MAF",h2part$category)]<-"MAF"
h2part$cat_source[grep("CpG",h2part$category)]<-"CpG"
h2part$cat_source[grep("Recomb_Rate",h2part$category)]<-"Recomb_Rate"
h2part$cat_source[grep("Backgrd_Selection",h2part$category)]<-"Backgrd_Selection"
h2part$cat_source[grep("Nucleotide_Diversity",h2part$category)]<-"Nucleotide_Diversity"
# v2.2 LD
h2part$cat_source[grep("Villar",h2part$category)]<-"Villar"
h2part$cat_source[grep("BLUEPRINT",h2part$category)]<-"BLUEPRINT"
h2part$cat_source[grep("phastCons46way",h2part$category)]<-"phastCons46way"
h2part$cat_source[grep("BivFlnk",h2part$category)]<-"Roadmap"
h2part$cat_source[grep("Ancient_Sequence_Age",h2part$category)]<-"Hujoel"

# cat type
h2part$cat_type<-sapply(strsplit(h2part$category,"_"),"[[",1)
h2part$cat_type[grep("MAF",h2part$category)]<-"MAF"
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# edit celltype results file to have relevant columns
h2celltype$model<-gsub(".cell_type_results.txt","",sapply(strsplit(h2celltype$file,"_h2_"),"[[",2))
h2celltype$phenotype<-sapply(strsplit(h2celltype$file,"_h2_"),"[[",1)
h2celltype$sample<-"total"
h2celltype$sample[grep("_females",h2celltype$phenotype)]<-"females"
h2celltype$sample[grep("_males",h2celltype$phenotype)]<-"males"
h2celltype$adj<-""
h2celltype$adj[grep("adjTBV",h2celltype$phenotype)]<-"adjTBV"
# clean pheno name
h2celltype$pheno<-gsub("_","",gsub("_females|_males|.adjTBV","",gsub("^_","",gsub("Volume_of_grey_matter_in_|_VOLUME_|Planum_Temporale","",h2celltype$phenotype))))
# compute -log(Pvalue) for the enrichment, to plot
h2celltype$logP<-(-log10(h2celltype$Coefficient_P_value))
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# Plots
#----------------------------------------------------------------------
## Enrichment of functional categories
plot_enrichment<-function(data,m,p="AI") {
  t<-subset(data,model==m)
  n<-length(unique(t$Category))
  # sort by phenotype of interest
  levs<-t$category[t$pheno==p][order(t$Enrichment[t$pheno==p],decreasing = FALSE)]
  # convert category into factor
  t$category<-factor(t$category,levels=levs)
  #
  max<-max(t$Enrichment)+10
  min<-min(t$Enrichment)+10

  if (max>100) { max<-100}
  if (min<100) { min<-(-100)}
  
  plot_e<-ggplot(t,aes(fill=cat_type,alpha=sig)) + 
    geom_bar(position = 'dodge',stat='identity',aes(x=category,y=Enrichment))  +
    geom_errorbar(aes(x=category,ymin=Enrichment-Enrichment_std_error, ymax=Enrichment+Enrichment_std_error,width=.5, color=cat_type), position=position_dodge(1),alpha=0.7) +
    facet_grid(.~pheno,scales="free") +
    geom_hline(yintercept=1) +
    mytheme +
    theme(legend.position="bottom") +
    labs(title=gsub("_"," ",m),x="",y="Enrichment\nPr(h2)/Pr(SNPs)",fill="") +
    # theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_alpha_discrete(range=c(0.2,1)) +
    guides(color = FALSE,alpha=guide_legend(show=FALSE)) +
    coord_flip(ylim = c(min,max), expand = TRUE) +
    # coord_flip() +
    NULL
  #
  return(plot_e) 
  
}

h2part_p<-subset(h2part,category!="base"&sample=="males"&adj==""&extend500=="0")
plot_enrichment(data=h2part_p,m="baseline",p="AI") 
plot_enrichment(data=h2part_p,m="baselineLD",p="AI")


#----------------------------------------------------------------------
## Enrichment of cell types

plot_celltype_coef<-function(data,m) {
  t<-subset(data,model==m)
  n<-length(unique(t$Name))
  
  plot_celltype_t<-ggplot(t,aes(x=Name,y=logP,color=Name)) + 
    geom_jitter(size=3) +
    # geom_bar(position = 'dodge',stat='identity')  + 
    facet_grid(~pheno) +
    mytheme +
    # theme(axis.text.x=element_text(angle=90,hjust=1)) +
    # coord_flip() +
    theme(legend.position="bottom",
          axis.text.x = element_blank()) +
    labs(title=gsub("_"," ",m),x="",y="-log10(P)",color="") +
    geom_hline(yintercept = -log10(0.05),linetype="dashed") +
    geom_hline(yintercept = -log10(0.05/n),linetype="solid") +
    NULL
  return(plot_celltype_t)
}


h2celltype_p<-subset(h2celltype,sample=="total"&adj=="")

for (m in unique(h2celltype_p$model)){
  plot_t<-plot_celltype_coef(data=h2celltype_p,m=m)
  assign(paste("plot_celltype_",m,sep=""),plot_t)
  rm(plot_t)
}
#m="Cahoy"


#----------------------------------------------------------------------
