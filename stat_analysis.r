###############################################################################
########################      STATISTICAL ANALYSIS      #######################
###############################################################################


#### Description
# This script performs a statistical analysis of the count matrix obtained by the workflow. It uses the DESeq2 and FactoMineR packages.
# It takes the count file as input and returns various files as output. The files which are returned are a PCA graph, the table of the 
# results of the differential expression analysis, a plot of the counts of the gene most differentially expressed between the two groups,
# the heatmap of the normalized counts of the most variable genes and the name of these. genes in a text file.


#### Library Loading
install.packages("FactoMineR")  # Installation of the package to perform the PCA
library(DESeq2)                 # Loading package for differential analysis
library(FactoMineR)             # Loading packahe for the PCA


#### Data Loading
# Counts Matrix
count_file = read.table("output.counts", header=T)    # Loading count matrix file 
counts = count_file[,7:14]
rownames(counts) = count_file[,1]

parserexpcolnames <- function(nomamodifier){     # remove ".bam" in colnames function
  nomamodifier <- substring(nomamodifier, 1, 9)
  return(nomamodifier)
}

colnames(counts) = sapply(colnames(counts), parserexpcolnames)   # remove ".bam" in colnames

# Annotation Phenotype
list_phenotype = c()   # list containing phenotype of each SRR
for (i in 1:dim(counts)[2]){
  if((colnames(counts)[i] == 'SRR628582') | (colnames(counts)[i] == 'SRR628583') | (colnames(counts)[i] == 'SRR628589')){
    list_phenotype = c(list_phenotype, "Mutated")
  }
  else{
    list_phenotype = c(list_phenotype, "Wild")
  }
}
phenotype = as.data.frame(matrix(list_phenotype, ncol = 1))
rownames(phenotype) = colnames(counts)
colnames(phenotype) = "Type"



#### Analysis

## PCA

res_pca = PCA(t(counts), graph= FALSE)              # Calculation of PCA
pdf("PCA_GraphOfIndividuals.pdf")
plot(res_pca, axes = c(1, 2), choix = c("ind"))     # Calculation is made in 5 dimensions, and the display in 2 
dev.off()


## Differencial Analysis

counts = counts[rowSums(counts)>0,]   # Removal of all unexpressed genes
counts=as.matrix(counts)              # Converting the object's class


dds = DESeqDataSetFromMatrix(countData = counts, colData = phenotype, design = ~ Type)  # Object construction
dds = DESeq(dds)     # performs differential expression analysis
res = results(dds)   # differential expression analysis
res_mat = as.matrix(res)
write.table(res_mat, "DESeq_results.txt", sep="\t")  # output file (table format txt) resulting of the differential expression analysis

pdf("plot_counts.pdf")
plotCounts(dds, gene=which.min(res$padj), intgroup="Type")  # output file (pdf) showing the most variable gene bewteen the 2 groups (wild and mutated)
dev.off()


select = order(res$padj, decreasing = FALSE)[1:30]  # Sorting and selection of the most variable genes

pdf("heatmap_MostVariableGenes.pdf")
heatmap(counts(dds, normalized=T)[select,], cexRow = 0.5, cexCol = 0.8, Rowv = NA)   # output file (pdf) showing normalized counts of the most variable genes for each individual
dev.off()

MostVariableGenes = rownames(counts[select,])  # Retrieving names names of the most variable genes

write.table(MostVariableGenes, "MostVariableGenes.txt", sep=";")   # output file (list format txt) containing names of the most variable genes