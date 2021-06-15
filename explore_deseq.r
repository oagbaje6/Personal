#!/usr/local/bin/Rscript

#Loading data sample and assay data
#sampleFile <- "/courses/example-data/csv/toy_sample_table.csv"
#assayFile <- "/courses/example-data/csv/toy_assay_table.csv"

sampleFile <- "/courses/example-data/csv/LUAD_sample_table.1.csv"
assayFile <- "/courses/example-data/csv/LUAD_assay_table.1.csv"


sampleTable <- read.table(sampleFile, row.names = 1, header = TRUE, sep = ",")
assayTable <-read.csv(assayFile, row.names=1)


library("DESeq2")
dds <- DESeqDataSetFromMatrix(countData = assayTable, colData = sampleTable, design = ~sampleGroup)


# Removing rows of the DESeqDataSet that have no counts, or only a single count
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 1, ]
nrow(dds)


# Stabilizing variance with regularized-logarithm transformation

rld <- rlog(dds, blind=FALSE)
head(assay(rld), 3)

#Sample distances
sampleDists <- dist( t( assay(rld) ) )

# Visualize distances with pheatmap
library("pheatmap")
library("RColorBrewer")
sampleDistMatrix <- as.matrix( sampleDists )
#rownames(sampleDistMatrix) <- paste( rld$sampleGroup, rld$cell, sep="-" )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists, col=colors)

#PCA plot
library(ggplot2)
plotPCA(rld, intgroup = c("sampleGroup")) + geom_text(aes(label=name),vjust=2)
