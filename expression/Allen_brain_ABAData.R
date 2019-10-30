source("https://bioconductor.org/biocLite.R")
biocLite("ABAData")
biocLite("ABAEnrichment")
require(ABAData);require(ABAEnrichment)

## get name and superstructures of brain region 4679
get_name(4679)
get_superstructures(4679)
get_name(get_superstructures(4679))
## require averaged gene expression data (microarray) from adult human brain regions
data(dataset_adult)
## look at first lines
head(dataset_adult)
## -----------------------------------------------------------------------------------------------------------
## require averaged gene expression data (RNA-seq) for 5 age categories
data(dataset_5_stages)
## look at first lines
head(dataset_5_stages)

## -----------------------------------------------------------------------------------------------------------
## require developmental effect score for genes in brain regions
data(dataset_dev_effect)
## look at first lines
head(dataset_dev_effect)

## -----------------------------------------------------------------------------------------------------------
# gene of interest
i<-paste("ITIH",1:5,sep="")
c<-paste("^",paste(i,collapse="$|^"),"$",sep="")
# get subsets
d1<-dataset_adult[grep(c,dataset_adult$hgnc_symbol),]
d2<-dataset_5_stages[grep(c,dataset_5_stages$hgnc_symbol),]
d3<-dataset_dev_effect[grep(c,dataset_dev_effect$hgnc_symbol),]

# 
d1$structure_name<-get_name(d1$structure)
d2$structure_name<-get_name(d2$structure)
d3$structure_name<-get_name(d3$structure)
regions_d1<-paste("Allen:",unique(d1$structure),sep="")
regions_d2<-paste("Allen:",unique(d2$structure),sep="")
regions_d3<-paste("Allen:",unique(d3$structure),sep="")


# some random plots
ggplot(data=d1) + geom_point(aes(x=structure,y=signal,color=as.factor(structure))) + facet_wrap(~hgnc_symbol,ncol=1)
ggplot(data=d2) + geom_line(aes(x=age_category,y=signal,color=structure_name)) + facet_wrap(~hgnc_symbol,ncol=1)
ggplot(data=d3) + geom_point(aes(x=structure,y=signal,color=structure_name))+ facet_wrap(~hgnc_symbol,ncol=1)

# use package functions to plot
plot_expression(structure_ids=regions_d1, 
                gene_ids=unique(d1$ensembl_gene_id), dataset="adult", dendro=TRUE)

plot_expression(structure_ids=regions_d2, 
                gene_ids=unique(d1$ensembl_gene_id), dataset="5_stages", dendro=TRUE)

plot_expression(structure_ids=regions_d3, 
                gene_ids=unique(d1$ensembl_gene_id), dataset="dev_effect", dendro=TRUE)

