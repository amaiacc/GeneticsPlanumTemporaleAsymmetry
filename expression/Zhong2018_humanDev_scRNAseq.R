#############################################
# adapted from: Zhong_etal_plot_script.R
##title: "Zhong et al scRNAseq expression data plotting"
##author: "Else Eising"
##date: "9 April 2019"
##############################################

rm(list=ls())
##---------------------------------------------------
# Set directories, dependent on  system
if (Sys.info()['sysname']=='Windows') {dir="P://"} else {dir="/data/workspaces/lag/"}
working_dir=paste(dir,"shared_spaces/Resource_DB/scRNAseq/Zhong_scRNAseq_neurodevelopment/",sep="")
setwd(working_dir)
##---------------------------------------------------

##---------------------------------------------------
# Step 1: get libraries loaded
##---------------------------------------------------

library(ggplot2)
# library(ggthemes)
library(plyr)
library(reshape2)

#### Steps to follow:
# 1) load and organize the data
# 2) Rename ages and brain structures
# 3) Make a subset of the data you are interested in
# 4) Plot results

##---------------------------------------------------
# Step 2: load and organize expression data
##---------------------------------------------------
TPM <- read.table("GSE104276_all_pfc_2394_UMI_TPM_NOERCC.txt", header=TRUE)
samples = read.table("GSE104276_Sample_Info.txt", header=TRUE, sep="\t")
samples$class<-NA
table(samples$cell_type)
samples$class[samples$cell_type=="Neurons"|samples$cell_type=="Stem cells"|samples$cell_type=="GABAergic neurons"]<-"Neuronal"
samples$class[(samples$cell_type!="Neurons"&samples$cell_type!="Stem cells"&samples$cell_type!="GABAergic neurons")]<-"Glial"
table(samples$class)
#pseudotime info


#subclusters
astrocytes <- read.table("GSE104276_PFC_astrocyte_rf_ident.txt", header=TRUE)
microglia <- read.table("GSE104276_PFC_microglia_rf_ident.txt")
opc <- read.table("GSE104276_PFC_OPC_rf_c4_ident.txt")
en <- read.table("GSE104276_PFC_ExcitatoryNeurons_rf_c7.txt")
#
npc <- read.table("GSE104276_PFC_NPCs_rf_c9_monocle_ident.txt")
int <- read.table("GSE104276_PFC_interneurons_rf_c8_monocle.txt")
neuroglia <- read.table("GSE104276_PFC_NPCneuglia_monocle_ident.txt")
neuroglia$sample<-rownames(neuroglia)
#
sub1<-rbind(int,npc)
sub1$sample<-rownames(sub1)
sub2<- rbind(astrocytes, microglia,  en, opc)
sub2$sample<-rownames(sub2)
sub3<-merge(sub2,neuroglia,all=TRUE)
sub<-merge(sub1,sub3,all=TRUE,stringsAsFactors=FALSE)
#
samples2 <- merge(samples,neuroglia,all.x=TRUE)
rm(sub1,sub2,sub3,sub)
#


# keep data from samples in sample file
TPM2 <- TPM[,which(colnames(TPM) %in% samples2$sample)]

# log-transform TPM
TPM3 <- log((TPM2/10) + 1)

# remove low-expressed genes
zeros <- rowSums(TPM3 < 1)
maxzeros <- ncol(TPM3) - 3
TPM4 <- TPM3[which(zeros <= maxzeros),]

# transpose and merge with sample info
# TPM5 <- as.data.frame(t(TPM4))
# TPM5$sample <- rownames(TPM5)
# TPM6 <- join(samples2, TPM5, by="sample")

##---------------------------------------------------
# Step 3: select data of your genes of interest (GOI)
##---------------------------------------------------

# GOI <- c("PAX6","SOX2","OLIG1","EOMES","DAAM1") 
GOI <- c("ITIH5","BOK","DTYMK")

GOI_data0 <- TPM4[GOI,]
GOI_data1 <- as.data.frame(t(GOI_data0))
GOI_data1$sample <- rownames(GOI_data1)
GOI_data2 <- join(samples2[,c("sample","week","cell_type","Pseudotime","State","class")], GOI_data1, by="sample")
print(GOI_data2[1:10,])


##---------------------------------------------------
# Step 4: Use melt to get the data in a "long" shape
##---------------------------------------------------

GOI_data3 <- melt(GOI_data2, 
                  id.vars = c("sample","week","cell_type","Pseudotime","State","class"), 
                  variable.name = "Gene", 
                  value.name = "TPM")
GOI_data3$TPM <- as.numeric(GOI_data3$TPM) #must make sure R knows these are numeric values
summary(GOI_data3$Gene)
hist(GOI_data3$TPM)


##---------------------------------------------------
# Step 5: Plot
##---------------------------------------------------

ggplot(subset(GOI_data3,!is.na(cell_type)),aes(x = Pseudotime, y = TPM)) +
			geom_jitter(aes(colour = cell_type),
			            shape=16, size = 2, alpha=0.5,  width=0.3) +
      geom_smooth(aes(colour=cell_type), method="auto",  se=FALSE, size=2) +
      facet_grid(. ~ Gene, scales = "free_y") +
      theme_classic() +
			theme(strip.background = element_rect(colour = "white", fill = FALSE)) +
			theme(panel.border = element_rect(colour = "black", fill=FALSE, size=1)) +
            # theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
            theme(strip.text.y = element_text(angle = 0, size = 12))+
			theme(legend.position = "bottom") +
            labs(y = "log10(TPM)", x = "Pseudotime") +
  NULL

