#!/usr/local/bin/Rscript
# RNA seq analysis pipeline

# prep
sampleTable <- read.csv("/courses/example-data/csv/sample_table.reduced.csv",row.names=1)
fileNames <- file.path("/courses/example-data/bam/", paste0(sampleTable$Run, "_all.bam"))

library("Rsamtools")

# Locations of bam files
bamFiles <- BamFileList(fileNames, yieldSize=2000000)

library("GenomicFeatures")

# Generating text database with GFF annotations
txdb <- makeTxDbFromGFF("/courses/example-data/gtf/Homo_sapiens.GRCh37.75.gtf", format="gtf",circ_seqs=character())
ebg <- exonsBy(txdb, by="gene")

library("GenomicAlignments")
library("BiocParallel")
register(SerialParam())

# Resolves overlapping reads
se <- summarizeOverlaps(features=ebg, reads=bamFiles, mode="Union", singleEnd=FALSE, ignore.strand=TRUE, fragments=TRUE )
colData(se) <- DataFrame(sampleTable)

# Changing factor levels to focus on untreated
se$dex <- relevel(se$dex, "untrt")


# PT2
library("DESeq2")

# Directs input data into range standardized experiment
dds <- DESeqDataSet(se, design = ~ dex)
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 1, ]
nrow(dds)

# Transforms counts to log scale, returns top three values
rld <- rlog(dds, blind=FALSE)
head(assay(rld), 3)
sampleDists <- dist( t( assay(rld) ) )

library("pheatmap")
library("RColorBrewer")

sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( rld$dex, rld$cell, sep="-" )
colnames(sampleDistMatrix) <- NULL

# colors for heatmap
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists, col=colors)

# Preparing PCA map
plotPCA(rld, intgroup = c("dex", "cell"))


# PT 3
library("DESeq2")

# Analysis
dds <- DESeqDataSet(se, design = ~ cell + dex)
dds <- DESeq(dds)

# Results, metadata, and summary
res <- results(dds)
mcols(res, use.names=TRUE)
summary(res)

# Filtering for significance
res.05 <- results(dds, alpha=.05)
table(res.05$padj < .05)
resLFC1 <- results(dds, lfcThreshold=1)
table(resLFC1$padj < 0.1)


# PT 4
# Results and comparison
results(dds, contrast=c("cell", "N052611", "N080611"))


# PT 5
# summing p-values, removing NA
sum(res$pvalue < 0.05, na.rm=TRUE)
sum(!is.na(res$pvalue))
sum(res$padj < 0.1, na.rm=TRUE)

# Subsetting p-values and ordering top  values
resSig <- subset(res, padj < 0.1)
head(resSig[ order(resSig$log2FoldChange), ])
head(resSig[ order(resSig$log2FoldChange, decreasing=TRUE), ])


# PT 6
# Choosing significant p-values and setting dex as interesting group
topGene <- rownames(res)[which.min(res$padj)]
plotCounts(dds, gene=topGene, intgroup=c("dex"))

# PT 7
library("AnnotationDbi")
library("org.Hs.eg.db")

columns(org.Hs.eg.db)

# Extracts data from database
res$symbol <- mapIds(org.Hs.eg.db, keys=row.names(res),column="SYMBOL", keytype="ENSEMBL",multiVals="first")
res$entrez <- mapIds(org.Hs.eg.db, keys=row.names(res), column="ENTREZID", keytype="ENSEMBL", multiVals="first")

# Rearranges p-values, selects first 100 valuees
resOrdered <- res[order(res$padj),]
head(resOrdered)
resOrderedDF <- as.data.frame(resOrdered)[1:100,]

# Saves results as csv
write.csv(resOrderedDF, file="results.csv")
