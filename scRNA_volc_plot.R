require(ggplot2)
library(dplyr)
library(data.table)
library(gplots)
library(tidyr)


##Highlight genes that have an absolute fold change > 2 and a p-value < Bonferroni cut-off

# Enter the file location
deseq_res <- fread("DESEq2_all_results.txt")

# Mark genes that are significant with paj < 0.05
deseq_res$threshold <- as.factor(deseq_res$padj < 0.05)

# Set Range to squish data beyond some value
deseq_res$log_q_inv <- scales::squish(sapply(deseq_res$padj, function (x) -log10(x)), range = c(0, 2))

# Construct the plot
g = ggplot(data=deseq_res, aes(x=log2FoldChange, y=log_q_inv)) +
  geom_point(aes(colour = threshold), alpha=0.4, size=2) + # The size affects the point-size
  theme(legend.position = "none", panel.background = element_blank(), legend.key.size = unit(15, "mm"), axis.text = element_text(size = 20)) +
  xlim(c(min(deseq_res$log2FoldChange)-1, max(deseq_res$log2FoldChange)+1)) +
  ylim(c(0, 2)) +
  xlab("") + ylab("") + scale_colour_manual(values = c("TRUE" = "Red", "FALSE" = "BLUE" ))
png("Volcano_plot.png")
g
dev.off()
