# Checking expression of genes Allen Brain Atlas
cd /data/workspaces/lag/workspaces/lg-dyslexia-exomes/working/working_data/

#------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------
## PREPARE SAMPLES 
#------------------------------------------------------------------------------------------
# Define path for the downloaded Allen Brain Atlas data
# Also, the anmes of the six brains (which will be in different directories within $allen_path
#------------------------------------------------------------------------------------------
allen_path='/data/workspaces/lag/workspaces/lg-dyslexia-exomes/working/expression_data/AllenBrainAtlas'
allen_samples=(H0351.1009  H0351.1012  H0351.1015  H0351.1016  H0351.2001  H0351.2002)

#------------------------------------------------------------------------------------------
# Each of the allen_samples has the same structure, which is explained in the Readme.txt file.
# In short:
ls $allen_path/${allen_samples[0]}
#	1	MicroarrayExpression.csv (probe,samples)									Probe's expression above background
#	2	Probes.csv									Metadata for  1. Name of the gene in column gene_symbol
#	3	PACall.csv									Present/Absent flags for expression different to background
#	4	Ontology.csv								The ontology on brain structures used for sampling. Brain areas etc
#	5	SampleAnnot.csv								The samples are listed in the same order as the columns in MicroarrayExpression.csv
#------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------
## Save order of brain structures for each samples to file $allen_path/$i/structure_acronym_names
## SampleAnnot.csv $5 contains acronym name for the expression pattern levels, when transposed, it should be the header name of the MicroarrayExpression.csv
#------------------------------------------------------------------------------------------
#for i in "${allen_samples[@]}"
# do
#  awk -F, '{print $5}' $allen_path/$i/SampleAnnot.csv | paste -sd ','| sed s/structure_acronym/probe_id/g > $allen_path/$i/structure_acronym_names
# done
## Information about the hemisphere embedded in the structure_names
#for i in "${allen_samples[@]}"
# do
#   awk -F'","|,' '{print $6,$7,$8}' $allen_path/$i/SampleAnnot.csv | paste -sd ','| sed s/'structure_name polygon_id mri_voxel_x'/probe_id/g > $allen_path/$i/structure_names
# done
#------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------
## SELECT GENES
#------------------------------------------------------------------------------------------
# Define file containing genes of interest, and convert to unix format
#------------------------------------------------------------------------------------------
cd /data/workspaces/lag/workspaces/lg-dyslexia-exomes/working/working_data/
# check just in one sample, since all samples have same probe info
wc $allen_path/H0351.*/Probes.csv
i=H0351.1009
#gene="ITIH5"
# SPINT2 has name AC011479.2 in the data: checked by grepping one of the probes from: http://human.brain-map.org/microarray/search/show?exact_match=false&search_term=SPINT2&search_type=gene
for gene in DTYMK BOK-AS1 ING5 AC114730.11 #AC011479.2 PPP1R14A # ITIH5 BOK C19orf12 NADK TMEM52 SLC35E2A
 do 
 #------------------------------------------------------------------------------------------
 # probes per gene into gene folder
 #------------------------------------------------------------------------------------------
 if [ ! -d "$allen_path/genes/$gene" ]; then
  mkdir -p $allen_path/genes/$gene
  grep '"'$gene'"' $allen_path/$i/Probes.csv | awk -F , '{OFS=",";print $1}'  > $allen_path/genes/$gene/probes_$gene
 fi

 #------------------------------------------------------------------------------------------
 # for each gene --> for each sample, find the expression of the probes for the given gene ($gene)
 #------------------------------------------------------------------------------------------
 #if [ -z `ls $allen_path/genes/$gene/*.csv` ]
 #then  
 # check that the csv files have not been generated for this gene yet
 for i in "${allen_samples[@]}"; do
  while read j; do
   awk 'BEGIN{FS=OFS=","}{$1="";sub(",","")}1' $allen_path/$i/structure_acronym_names > $allen_path/genes/$gene/exp'_'$gene'_probe'$j'_sample'$i.csv
   awk 'BEGIN{FS=OFS=","}{$1="";sub(",","")}1' $allen_path/$i/structure_names >> $allen_path/genes/$gene/exp'_'$gene'_probe'$j'_sample'$i.csv
   grep ^$j $allen_path/$i/MicroarrayExpression.csv | awk 'BEGIN{FS=OFS=","}{$1="";sub(",","")}1'  >> $allen_path/genes/$gene/exp'_'$gene'_probe'$j'_sample'$i.csv
   grep ^$j $allen_path/$i/PACall.csv | awk 'BEGIN{FS=OFS=","}{$1="";sub(",","")}1'  >> $allen_path/genes/$gene/exp'_'$gene'_probe'$j'_sample'$i.csv
  done < $allen_path/genes/$gene/probes_$gene
 done
 #fi
 done

#------------------------------------------------------------------------------------------
# Create directories for the plots of each gene
#------------------------------------------------------------------------------------------
plots='/data/workspaces/lag/workspaces/lg-dyslexia-exomes/working/working_data/expression/allen_plots'
mkdir -p $plots/$gene

#------------------------------------------------------------------------------------------
# Run Rscript candidate_genes_allen_expression to get output table with info of all genes
#------------------------------------------------------------------------------------------
cand_genes='/data/workspaces/lag/workspaces/lg-dyslexia-exomes/working/working_data/expression/genes_PT.txt'
echo ITIH5 BOK C19orf12 NADK TMEM52 PPP1R14A AC011479.2 DTYMK BOK-AS1 ING5 AC114730.11 | sed 's/ /\n/g' > $cand_genes # SLC35E2A exclude because there is no expression info #AC011479.2
r_genes='/data/workspaces/lag/workspaces/lg-dyslexia-exomes/working/analysis/AllenBrainAtlas/candidate_genes_allen_expression.R'
# first argument: path to candidate list file
Rscript --verbose $r_genes $cand_genes

# get values and plots: /data/workspaces/lag/workspaces/lg-ukbiobank/analysis/amaia/genetic_analysis/PT/candidate_genes_allen_expression_byregions_PT.R