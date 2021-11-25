## Analysis Script

library(DESeq2)

# Annotations

#count_file = read.table("counts.txt", header = TRUE)
count_file = read.table("output.counts", header=T)
counts = count_file[,7:14]
rownames(counts) = count_file[,1]

list_names = c("SRR628584.bam", "SRR628587.bam", "SRR628586.bam", "SRR628583.bam", "SRR628585.bam", "SRR628588.bam", "SRR628582.bam", "SRR628589.bam")
list_phenotype = c("mutated", "wild", "wild", "mutated", "wild", "wild","mutated", "wild")

phenotype = data.frame(list_phenotype[1:8])
colnames(phenotype) = "type"
rownames(phenotype) = list_names[1:8]



## supression des gènes "non exprimés"
datexp= counts[rowSums(counts)>0,]

## formatage de la matrice de comptage
dataexp=as.matrix(datexp)

## visualisation des log comptages
hist(log(dataexp+1),col="blue")



# construction de l'objet 
dds <- DESeqDataSetFromMatrix(dataexp, DataFrame(phenotype), ~ type)
dds <- estimateSizeFactors(dds)
dds <- dds[ rowSums(counts(dds)) > 0, ]
rld <- vst(dds, blind=TRUE)
d.vst <- assay(rld)


# Récupération des 1000 gènes qui varient le plus #
gvar <- apply(d.vst, 1, var)
mostvargenes <- order(gvar, decreasing=TRUE)[1:3000]
d.vst_1000genes <- d.vst[mostvargenes,]


pdf("heatmap.pdf")
heatmap(d.vst_1000genes)
dev.off()

write.table(d.vst_1000genes, "d.vst_1000genes.txt", sep="\t")
