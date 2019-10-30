#############################################

#----------------------------------------------------------------------
# Set working directory
#----------------------------------------------------------------------
if (Sys.info()['sysname']=='Windows') {dir="P://workspaces/"} else {dir="/data/workspaces/lag/workspaces/"}
working_dir=paste(dir,"lg-ukbiobank/working_data/amaia/genetic_data/PT/brainspan/",sep="")
setwd(working_dir)
#----------------------------------------------------------------------

#----------------------------------------------------------------------
library(reshape2)
library(ggplot2)
library(ggthemes)
require(grid)
library("GGally")
library(limma)
library(tidyr)
library(dplyr)
library(RColorBrewer)
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# Read brainspan data
#----------------------------------------------------------------------
bs_dir=paste(dir,"lg-dyslexia-exomes/working/expression_data/BrainSpan/genes_matrix_csv/",sep="")
bs_row = read.csv(paste(bs_dir,"rows_metadata.csv",sep=""))
bs_col = read.csv(paste(bs_dir,"columns_metadata.csv",sep=""))
brainspan <- read.csv(paste(bs_dir,"expression_matrix.csv",sep=""), header=FALSE, row.names = 1)
brainspan <- cbind(brainspan, bs_row)
rownames(brainspan)<-as.character(brainspan$gene_symbol)
# brainspan_data1<-as.data.frame(t(brainspan))

#----------------------------------------------------------------------

#----------------------------------------------------------------------
# Format data: select specific genes and brain regions
#----------------------------------------------------------------------
#Steps:
#1. Define Genes of interest
#2. Subset the brainspan data by that GOI list (x rows by 524/529 cols)
#3. Transpose the data, so each gene is a column
#4. Use cbind to attach some of the Brainspan column metadata (bs_col) to the GOI data
#5. Use melt from the reshape2 package to turn the data into its “Long” form, rather than “Wide”. This is needed for ggplot2.
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# define genes of interest
GOI<-c("SLC35E2A","NADK","TMEM52","BOK","C19orf12","ITIH5","SPINT2","PPP1R14A",
       "BOK-AS1","ING5","DTYMK","AC114730.11"
       ) #paste("ITIH",1:5,sep="")
# GOI<-"ITIH5"
#-----------------------
# Get indices for the genes of interest on the brainspain data, which rows do they correspond to?
a<-unlist(sapply(GOI,function(x) which(brainspan$gene_symbol==x)))
genes<-names(a) # this will be the actual genes, because some may be missing
# Subset and format
GOI_data0 <- brainspan[a,]  # select only rows for the genes of interest
GOI_data1 <- as.data.frame(t(GOI_data0)) # transpose: columns= genes, rows=regions
genes<-gsub("-","_",genes)
colnames(GOI_data1) <- genes
GOI_data2 <- data.frame(GOI_data1[grep("V",rownames(GOI_data1)),]) # remove metadata for genes: rows 525:529, make it more informative: all the region rows are defined as V, so select those
colnames(GOI_data2) <- genes
GOI_data3 <- cbind(GOI_data2,bs_col[,c("age","structure_acronym","structure_name")]) # add: age, structure acronym and structure name
GOI_data3[, genes] <-apply(GOI_data3[, genes],2,as.numeric)
GOI_data3_early<-GOI_data3[grepl("pcw|mos|^1 yrs|^2 yrs",GOI_data3$age),]
# clean: remove all intermediate objects
rm(GOI_data0,GOI_data1,GOI_data2,a)

# Check correlation across genes (all structures and time-points)
cor_genes<-ggpairs(GOI_data3[, genes],title = "Relation of neural expression across genes \nfor all structures and ages")
# ggsave(paste("plots/",name,"_p_cor_genes",".png",sep="") )
# ggsave("cases_genes_cor.png")
cor_genes_early<-ggpairs(GOI_data3_early[, genes],columns=genes,title = "Relation of neural expression across genes \nfrom 8 pcw till 2 years of age, and across all structures")
cor_genes_early2<-ggpairs(GOI_data3_early,aes(colour=structure_acronym ),columns=genes,title = "Relation of neural expression across genes \nfrom 8 pcw till 2 years of age, and across all structures")

#  plot mds
plotMDS(subset(GOI_data3)[,genes])

# dist<-dist(GOI_data3[,genes])
# fit <- cmdscale(dist, eig = TRUE, k = 2)
# x <- fit$points[, 1]
# y <- fit$points[, 2]
# plot(x, y, pch = 19, xlim = range(x) + c(0, 600))
# text(x, y, pos = 4, labels = genes)
# rm(dist,fit,x,y)

#-----------------------
# wide to long format
# Melt from will now set up the data so that the GOIs are stacked (from wide to long), with all the other column data repeated. 
# If you have 2 GOIs, the number of rows will double. This makes it easier to make plots in ggplot2, where we might want to show the expression data across time and have a line for each GOI on the same graph to compare the patterns.
GOI_data4 <- melt(GOI_data3, 
                  id.vars = c("age","structure_acronym","structure_name"), 
                  variable.name = "Gene", 
                  value.name = "FPKM")
GOI_data4$FPKM <- as.numeric(GOI_data4$FPKM)

# Number of rows per gene, one per region
summary(GOI_data4$Gene)

# find number of samples that have an expression level above the threshold in any tissue included in Brainspan
# filter out if less than 1 FPKM
threshold <- 1
GOI_data4<-subset(GOI_data4,FPKM>=threshold)
summary(GOI_data4$Gene)
# dist<-dist(GOI_data4[, c("Gene","FPKM")])
# fit <- cmdscale(dist, eig = TRUE, k = 2)
# x <- fit$points[, 1]
# y <- fit$points[, 2]
# plot(x, y, pch = 19, xlim = range(x) + c(0, 600))
# text(x, y, pos = 4, labels = genes)
#-----------------------

structures <- unique(GOI_data4$structure_name)
ages <- unique(bs_col$age)
GOI_data4$age <- factor(GOI_data4$age,levels=ages)
pwc<-grep("pcw",GOI_data4$age) # index of rows that contain entries before birth

# GOI_data4 <- mutate(GOI_data4, structure_name_wrap= str_wrap(structure_name, width = 15))
# str_exclude<-names(which(table(GOI_data4$structure_acronym)<length(GOI)*10))
# w<- which(GOI_data4$structure_acronym %in% str_exclude)
# if (w>0){GOI_data4<-GOI_data4[-w,]}
# rm(w)


# average per time-point
table(GOI_data4$age,GOI_data4$Gene,GOI_data4$structure_acronym)

GOI_data5 <- GOI_data4  %>% group_by(Gene,age,structure_acronym,structure_name) %>% summarise(avgFPKM = mean(FPKM)) 
table(GOI_data5$age,GOI_data5$Gene,GOI_data5$structure_acronym)

# additional variables for plotting
GOI_data5$area<-""
GOI_data5$area[grep("cortex",GOI_data5$structure_name)]<-"Cortex"
GOI_data5$area[grep("cerebell|rhombic lip",GOI_data5$structure_name)]<-"Cerebellum"
GOI_data5$area[grep("hippocamp|amygdal|anglionic eminence|striatum|caudate|thalamus",GOI_data5$structure_name)]<-"Subcortical"
# GOI_data5$area[grep("thalamus",GOI_data5$structure_name)]<-"Thalamus"
# GOI_data5$area[grep("striatum|caudate",GOI_data5$structure_name)]<-"Basal ganglia"
# GOI_data5$area[grep("ganglionic eminence",GOI_data5$structure_name)]<-"Ganglionic eminence"
table(GOI_data5$area)
#
GOI_data5$logAvgFPKM<-(log2(GOI_data5$avgFPKM))
# create variable that contains age info in years
table(GOI_data5$age)
GOI_data5$ageY<-NA
GOI_data5$ageY[grep("yrs",GOI_data5$age)]<-as.numeric(gsub(" yrs","",GOI_data5$age)[grep("yrs",GOI_data5$age)])
GOI_data5$ageY[grep("mos",GOI_data5$age)]<-as.numeric(gsub(" mos","",GOI_data5$age)[grep("mos",GOI_data5$age)])/12
GOI_data5$ageY[grep("pcw",GOI_data5$age)]<-(-as.numeric(gsub(" pcw","",GOI_data5$age)[grep("pcw",GOI_data5$age)])/(12*4.3))
GOI_data5$ageT<-as.numeric(GOI_data5$age)




# subset of data
GOI_data8 <- GOI_data5 %>% 
  filter(avgFPKM>threshold,
         grepl("pcw|mos|^1 yrs|^2 yrs",age)
         )


#-----------------------
# Plots
#-----------------------
# p_all_time<-ggplot(subset(GOI_data5), aes(x = age,y = logAvgFPKM,colour=Gene)) +
#             geom_line(aes(group=Gene)) +
#             geom_point() +
#             theme_classic()+
# 			geom_hline(yintercept=0, colour="black") +
#             theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))+
#             theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
#             theme(strip.text.y = element_text(angle = 0))+
#   facet_grid(Gene ~ ., scales = "free_y") +
#   # facet_wrap(~structure_name) +
#             labs(title = paste("Developmental gene expression - Brainspan"),name,sep="\n") +
#   scale_x_discrete(breaks = levels(GOI_data8$age)[seq(1, length(levels(GOI_data8$age)), by = 4)] ) + 
#   theme(legend.position="bottom", axis.text=element_text(size=10,colour="black"), axis.title=element_text(size=14,face="bold"), strip.text = element_text(size=13) ) +
#   scale_color_manual(values=c("purple","darkgrey","darkturquoise","green","red")) +
#   ylab(bquote(''*log[2]*' (average FPKM)')) +
#   NULL


p_all_time<-ggplot(subset(GOI_data5), aes(x = age, y = logAvgFPKM)) +
  geom_point(alpha=0.3) +   geom_line(aes(group=structure_name),alpha=0.3) +
  geom_smooth(color="red", method = "loess",se = TRUE, aes(x=ageT,y=logAvgFPKM)) +
  geom_vline(xintercept = 14,color="black",linetype="dashed",alpha=0.8) +
  theme_classic()+geom_hline(yintercept=0, colour="black") +  
  theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  theme(strip.text.y = element_text(angle = 0))+
  facet_grid(.~Gene, scales = "free_y") +
  # facet_wrap(~structure_name) +
  labs(title = paste("Developmental gene expression - Brainspan")) +
  scale_x_discrete(breaks = levels(GOI_data5$age)[seq(1, length(levels(GOI_data5$age)), by = 4)] ) + 
  theme(legend.position="bottom", axis.text=element_text(size=10,colour="black"), axis.title=element_text(size=14,face="bold"), strip.text = element_text(size=13) ) +
  # scale_color_manual(values=c("purple","darkgrey","darkturquoise","green","red")) +
  ylab(bquote(''*log[2]*' (average FPKM)')) +
  NULL

p_all_time_c<-ggplot(subset(GOI_data5), aes(x = age, y = logAvgFPKM, color=area)) +
  geom_point(alpha=0.1) +  
  # geom_line(aes(group=structure_name),alpha=0.2) +
  geom_smooth( method = "loess",se = FALSE, aes(x=ageT,y=logAvgFPKM,color=area)) +
  geom_vline(xintercept = 14,color="black",linetype="dashed",alpha=0.8) +
  theme_classic()+geom_hline(yintercept=0, colour="black") +  
  theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  theme(strip.text.y = element_text(angle = 0))+
  facet_grid(.~Gene, scales = "free_y") +
  # facet_wrap(~structure_name) +
  labs(title = paste("Developmental gene expression - Brainspan")) +
  scale_x_discrete(breaks = levels(GOI_data5$age)[seq(1, length(levels(GOI_data5$age)), by = 4)] ) + 
  theme(legend.position="bottom", axis.text=element_text(size=10,colour="black"), axis.title=element_text(size=14,face="bold"), strip.text = element_text(size=13) ) +
  scale_color_manual(values=c("darkgreen","blue","red","darkturquoise","darkgrey")) +
  ylab(bquote(''*log[2]*' (average FPKM)')) +
  NULL



ggsave(p_all_time,file=paste("GOI_brainspan_p_all",".pdf",sep=""),width=15,height=10)
ggsave(p_all_time_c,file=paste("GOI_brainspan_p_regions",".pdf",sep=""),width=15,height=10)


# per gene
# name="ITIH5"
for (name in GOI){
  tmp_p<-ggplot(subset(GOI_data5,Gene==name), aes(x = age, y = logAvgFPKM)) +
    geom_point(alpha=0.3) +   geom_line(aes(group=structure_name),alpha=0.3) +
    geom_smooth(color="red", method = "loess",se = TRUE, aes(x=ageT,y=logAvgFPKM)) +
    geom_vline(xintercept = 14,color="black",linetype="dashed",alpha=0.8) +
    theme_classic()+geom_hline(yintercept=0, colour="black") +  
    theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
    theme(strip.text.y = element_text(angle = 0))+
    # facet_grid(Gene ~ ., scales = "free_y") +
    # facet_wrap(~structure_name) +
    labs(title = paste("Developmental gene expression - Brainspan",name,sep="\n")) +
    scale_x_discrete(breaks = levels(GOI_data8$age)[seq(1, length(levels(GOI_data8$age)), by = 4)] ) + 
    theme(legend.position="bottom", axis.text=element_text(size=10,colour="black"), axis.title=element_text(size=14,face="bold"), strip.text = element_text(size=13) ) +
    # scale_color_manual(values=c("purple","darkgrey","darkturquoise","green","red")) +
    ylab(bquote(''*log[2]*' (average FPKM)')) +
    NULL
  ggsave(tmp_p,file=paste("GOI_brainspan_p_all_",name,".pdf",sep=""),width=5,height=5)
  assign(paste("p_",name,sep=""),tmp_p)
  
  tmp_p2<-ggplot(subset(GOI_data5,Gene==name),  aes(x = age, y = logAvgFPKM, color=area)) +
    geom_point(alpha=0.2) +  
    # geom_line(aes(group=structure_name),alpha=0.2) +
    geom_smooth( method = "loess",se = FALSE, aes(x=ageT,y=logAvgFPKM,color=area)) +
    geom_vline(xintercept = 14,color="black",linetype="dashed",alpha=0.8) +
    theme_classic()+geom_hline(yintercept=0, colour="black") +  
    theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
    theme(strip.text.y = element_text(angle = 0))+
    # facet_grid(.~Gene, scales = "free_y") +
    # facet_wrap(~structure_name) +
    labs(title = paste("Developmental gene expression - Brainspan")) +
    scale_x_discrete(breaks = levels(GOI_data5$age)[seq(1, length(levels(GOI_data5$age)), by = 4)] ) + 
    theme(legend.position="bottom", axis.text=element_text(size=10,colour="black"), axis.title=element_text(size=14,face="bold"), strip.text = element_text(size=13) ) +
    scale_color_manual(values=c("darkgreen","blue","red","darkturquoise","darkgrey")) +
    ylab(bquote(''*log[2]*' (average FPKM)')) +
    NULL
  ggsave(tmp_p2,file=paste("GOI_brainspan_p_regions_",name,".pdf",sep=""),width=5,height=5)
  assign(paste("p_regions_",name,sep=""),tmp_p2)
  
}


# only BOK, ITIH5 and # |Gene=="PPP1R14A"
t<-subset(GOI_data5,Gene=="ITIH5"|Gene=="BOK"|Gene=="DTYMK")
t$Gene<-factor(t$Gene,levels=c("ITIH5","BOK","DTYMK"))
p_all_time_c3<-ggplot(t, aes(x = age, y = logAvgFPKM, color=area)) +
 geom_point(alpha=0.1) +  
  # geom_line(aes(group=structure_name),alpha=0.2) +
  geom_smooth( method = "loess",se = FALSE, aes(x=ageT,y=logAvgFPKM,color=area)) +
  geom_vline(xintercept = 14,color="darkgrey",linetype="dashed") +
  theme_classic()+
  geom_hline(yintercept=0, colour="black") +
  theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  theme(strip.text.y = element_text(angle = 0))+
  facet_grid(~Gene) +
  theme(strip.background = element_rect(colour="white")) +
  # facet_wrap(~structure_name) +
  labs(title = paste("Developmental gene expression - Brainspan")) +
  scale_x_discrete(breaks = levels(GOI_data5$age)[seq(1, length(levels(GOI_data5$age)), by = 4)] ) + 
  theme(legend.position="bottom", axis.text=element_text(size=10,colour="black"), axis.title=element_text(size=14,face="bold"), strip.text = element_text(size=13) ) +
  scale_color_manual(values=c("darkgreen","blue","red","darkturquoise","darkgrey")) +
  ylab(bquote(''*log[2]*' (average FPKM)')) +
  
  NULL

#|Gene=="PPP1R14A"
p_all_time_c4<-ggplot(t, aes(x = age, y = logAvgFPKM, color=structure_name)) +
  geom_point(alpha=0.3) +  
  geom_line(aes(group=structure_name),alpha=0.2) +
  # geom_smooth( method = "loess",se = FALSE, aes(x=ageT,y=logAvgFPKM,color=area)) +
  geom_vline(xintercept = 14,color="darkgrey",linetype="dashed") +
  facet_grid(~Gene) +
  # facet_wrap(~structure_name) +
  labs(title = paste("Developmental gene expression - Brainspan")) +
  scale_x_discrete(breaks = levels(GOI_data5$age)[seq(1, length(levels(GOI_data5$age)), by = 4)] ) + 
  theme_classic()+geom_hline(yintercept=0, colour="black") +  
  theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  theme(strip.text.y = element_text(angle = 0))+
  theme(legend.position="bottom", axis.text=element_text(size=10,colour="black"), axis.title=element_text(size=14,face="bold"), strip.text = element_text(size=13) ) +
  # scale_color_manual(values=c("darkgreen","blue","red","darkturquoise","darkgrey")) +
  ylab(bquote(''*log[2]*' (average FPKM)')) +
  NULL

ggsave(p_all_time_c3,file=paste("ITIH5_BOK_DTYMK_brainspan_p_regions.png",sep=""),width=10,height=5)


p_all_time_tmp<-ggplot(t[grep("tempo",t$structure_name),], 
                            aes(x = age, y = logAvgFPKM, color=structure_name)) +
  geom_point(alpha=0.8) +  
  # geom_line(aes(group=structure_name),alpha=0.2) +
  geom_smooth( method = "loess",se = FALSE, aes(x=ageT,y=logAvgFPKM,color=structure_name)) +
  geom_vline(xintercept = 14,color="darkgrey",linetype="dashed") +
  theme_classic()+geom_hline(yintercept=0, colour="black") +  
  theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  theme(strip.text.y = element_text(angle = 0))+
  facet_grid(~Gene) +
  # facet_wrap(~structure_name) +
  labs(title = paste("Developmental gene expression - Brainspan Temporal")) +
  scale_x_discrete(breaks = levels(GOI_data5$age)[seq(1, length(levels(GOI_data5$age)), by = 4)] ) + 
  theme(legend.position="bottom", axis.text=element_text(size=10,colour="black"), axis.title=element_text(size=14,face="bold"), strip.text = element_text(size=13) ) +
  scale_color_manual(values=brewer.pal(9, "Blues")[c(4,6,9)]) +
  # scale_color_brewer(palette="Blues") +
  ylab(bquote(''*log[2]*' (average FPKM)')) +
  NULL

ggsave(p_all_time_tmp,file=paste("ITIH5_BOK_DTYMK_brainspan_p_temporal_regions.png",sep=""),width=10,height=5)



# include all possible genes with eQTL effects in 2q37
genes<-subset(GOI_data5,Gene=="ITIH5"|Gene=="BOK"|Gene=="BOK_AS1"|Gene=="ING5"|Gene=="DTYMK"|Gene=="AC114730.11")
genes$Gene<-factor(genes$Gene,levels=c("BOK","BOK_AS1","DTYMK","ING5","AC114730.11","ITIH5"))
levels(genes$Gene)<-c("BOK","BOK-AS1","DTYMK","ING5","AC114730.11","ITIH5")
p_all_time_10_2<-ggplot(genes, 
                      aes(x = age, y = logAvgFPKM, color=area)) +
  geom_point(alpha=0.1) +  
  # geom_line(aes(group=structure_name),alpha=0.2) +
  geom_smooth( method = "loess",se = FALSE, aes(x=ageT,y=logAvgFPKM,color=area)) +
  geom_vline(xintercept = 14,color="darkgrey",linetype="dashed") +
  theme_classic()+geom_hline(yintercept=0, colour="black") +  
  theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  theme(strip.text.y = element_text(angle = 0))+
  facet_grid(~Gene) +
  # facet_wrap(~structure_name) +
  labs(title = paste("Developmental gene expression - Brainspan")) +
  scale_x_discrete(breaks = levels(GOI_data5$age)[seq(1, length(levels(GOI_data5$age)), by = 4)] ) + 
  theme(legend.position="bottom", axis.text=element_text(size=10,colour="black"), axis.title=element_text(size=14,face="bold"), strip.text = element_text(size=13) ) +
  scale_color_manual(values=c("darkgreen","blue","red","darkturquoise","darkgrey")) +
  ylab(bquote(''*log[2]*' (average FPKM)')) +
  NULL

ggsave(p_all_time_10_2,file=paste("ITIH5_genes2q37_brainspan_p_regions.png",sep=""),width=15,height=5)


p_all_time_10_2_str<-ggplot(genes,
                      aes(x = age, y = logAvgFPKM, color=structure_name)) +
  geom_point(alpha=0.3) +  
  geom_line(aes(group=structure_name),alpha=0.2) +
  # geom_smooth( method = "loess",se = FALSE, aes(x=ageT,y=logAvgFPKM,color=area)) +
  geom_vline(xintercept = 14,color="darkgrey",linetype="dashed") +
  theme_classic()+geom_hline(yintercept=0, colour="black") +  
  theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  theme(strip.text.y = element_text(angle = 0))+
  facet_grid(~Gene) +
  # facet_wrap(~structure_name) +
  labs(title = paste("Developmental gene expression - Brainspan")) +
  scale_x_discrete(breaks = levels(GOI_data5$age)[seq(1, length(levels(GOI_data5$age)), by = 4)] ) + 
  theme(legend.position="bottom", axis.text=element_text(size=10,colour="black"), axis.title=element_text(size=14,face="bold"), strip.text = element_text(size=13) ) +
  # scale_color_manual(values=c("darkgreen","blue","red","darkturquoise","darkgrey")) +
  ylab(bquote(''*log[2]*' (average FPKM)')) +
  NULL



# temporal regions only
p_all_time_10_2_tmp<-ggplot(genes[grep("tempo",genes$structure_name),], 
                        aes(x = age, y = logAvgFPKM, color=structure_name)) +
  geom_point(alpha=0.8) +  
  # geom_line(aes(group=structure_name),alpha=0.2) +
  geom_smooth( method = "loess",se = FALSE, aes(x=ageT,y=logAvgFPKM,color=structure_name)) +
  geom_vline(xintercept = 14,color="darkgrey",linetype="dashed") +
  theme_classic()+geom_hline(yintercept=0, colour="black") +  
  theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  theme(strip.text.y = element_text(angle = 0))+
  facet_grid(~Gene) +
  # facet_wrap(~structure_name) +
  labs(title = paste("Developmental gene expression - Brainspan Temporal")) +
  scale_x_discrete(breaks = levels(GOI_data5$age)[seq(1, length(levels(GOI_data5$age)), by = 4)] ) + 
  theme(legend.position="bottom", axis.text=element_text(size=10,colour="black"), axis.title=element_text(size=14,face="bold"), strip.text = element_text(size=13) ) +
  scale_color_manual(values=brewer.pal(9, "Blues")[c(4,6,9)]) +
  # scale_color_brewer(palette="Blues") +
  ylab(bquote(''*log[2]*' (average FPKM)')) +
  NULL

ggsave(p_all_time_10_2_tmp,file=paste("ITIH5_genes2q37_brainspan_p_temp_regions.png",sep=""),width=15,height=5)
