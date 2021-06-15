#!/usr/local/bin/Rscript

#sampleFile <- "/courses/example-data/csv/toy_sample_table.csv"
#assayFile <- "/courses/example-data/csv/toy_assay_table.csv"

sampleFile <- "/courses/example-data/csv/LUAD_sample_table.1.csv"
assayFile <- "/courses/example-data/csv/LUAD_assay_table.1.csv"

sampleTable <- read.table(sampleFile, row.names = 1, header = TRUE, sep = ",")
assayTable <-read.csv(assayFile, row.names=1)

library("DESeq2")
dds <- DESeqDataSetFromMatrix(countData = assayTable, colData = sampleTable, design = ~sampleGroup)

nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 1, ]
nrow(dds)

# Running the differential expression pipeline
dds <- DESeq(dds)

# Building the results table
res <- results(dds, contrast=c("sampleGroup", "normal", "tumor"))

# Summarizing the results
summary(res)

# Changing padj to filter more significant genes
res.05 <- results(dds, alpha=.05)

res.05 <- na.omit(res.05)

filtered_res.05 <- res.05[res.05$padj < .05, ]
ordered_res.05 <- filtered_res.05[order(filtered_res.05$padj),]

write.csv(ordered_res.05, file = "genes_with_padj_under_0.05.csv")

table(res.05$padj < .05)


resLFC1 <- results(dds, lfcThreshold=1)
table(resLFC1$padj < 0.1)
