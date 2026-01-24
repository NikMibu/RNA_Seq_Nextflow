# ══════════════════════════════════════════════════════════
# PCA aus featureCounts Output
# ══════════════════════════════════════════════════════════

library(DESeq2)
library(ggplot2)
library(tidyverse)

# 1. Counts einlesen (einzelne featureCounts Dateien zusammenführen)
count_files <- list.files(
    path = "/mnt/d/RNA_seq/results/counts",
    pattern = "SRR.*\\.counts\\.txt$",
    full.names = TRUE
)

# Erste Datei für Gene IDs
first <- read.table(count_files[1], header=TRUE, row.names=1, skip=1, comment.char="")
gene_ids <- rownames(first)

# Count-Matrix erstellen
count_matrix <- matrix(0, nrow=length(gene_ids), ncol=length(count_files))
rownames(count_matrix) <- gene_ids
sample_names <- gsub("\\.counts\\.txt", "", basename(count_files))
colnames(count_matrix) <- sample_names

# Alle Counts einlesen
for (i in seq_along(count_files)) {
    temp <- read.table(count_files[i], header=TRUE, row.names=1, skip=1, comment.char="")
    count_matrix[, i] <- temp[gene_ids, 6]  # Spalte 6 enthält die Counts
}

counts <- as.data.frame(count_matrix)

# 2. Sample-Metadaten erstellen
sample_info <- data.frame(
    sample = colnames(counts),
    condition = c(
        "SF3B1_mutant", "SF3B1_mutant", "SF3B1_mutant",
        "SF3B1_wildtype", "SF3B1_wildtype", "SF3B1_wildtype",
        "SF3B1_wildtype", "SF3B1_wildtype"
    )
)
rownames(sample_info) <- sample_info$sample

# 3. DESeq2 Objekt erstellen (nur für Normalisierung)
dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = sample_info,
    design = ~ condition
)

# Filtern: Gene mit mindestens 10 Reads in mindestens 3 Samples
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep, ]

# 4. Variance Stabilizing Transformation (für PCA)
vsd <- vst(dds, blind = TRUE)

# 5. PCA berechnen
pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

# 6. Plot erstellen
ggplot(pca_data, aes(x = PC1, y = PC2, color = condition, label = name)) +
    geom_point(size = 4) +
    geom_text(vjust = -1, size = 3) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    scale_color_manual(values = c(
        "SF3B1_mutant" = "coral",
        "SF3B1_wildtype" = "steelblue"
    )) +
    theme_minimal(base_size = 14) +
    ggtitle("PCA: SF3B1 Wildtype vs. Mutant") +
    theme(legend.position = "bottom")

ggsave("/mnt/d/RNA_seq/results/PCA_SF3B1.pdf", width = 8, height = 6)