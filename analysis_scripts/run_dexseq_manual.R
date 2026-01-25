#!/usr/bin/env Rscript

library(DEXSeq)
library(tidyverse)

setwd("results/dexseq")

# Target genes
target_genes <- c("ABCC5", "CRNDE", "UQCC", "GUSBP11", "ANKHD1", "ADAM12")

# Read samples
samples <- read.csv("samplesheet.csv")
samples <- data.frame(
    condition = factor(samples$condition)
)
rownames(samples) <- read.csv("samplesheet.csv")$sample_id

cat("Sample data:\n")
print(samples)

# Read count files manually
count_files <- paste0(rownames(samples), ".txt")
cat("\nReading count files...\n")

# Read first file to get exon IDs - use quote="" to preserve quotes
first_file <- read.table(count_files[1], header=FALSE, stringsAsFactors=FALSE, sep="\t", quote="")
exon_ids <- first_file[,1]
cat("Number of exons:", length(exon_ids), "\n")
cat("Example exon IDs:\n")
print(head(exon_ids, 3))

# Create count matrix (will set rownames after parsing)
count_matrix <- matrix(0, nrow=length(exon_ids), ncol=length(count_files))
colnames(count_matrix) <- rownames(samples)

# Fill count matrix
for (i in seq_along(count_files)) {
    cat("Reading", count_files[i], "...\n")
    counts <- read.table(count_files[i], header=FALSE, stringsAsFactors=FALSE, sep="\t", quote="")
    # Verify exon IDs match
    if (!all(counts[,1] == exon_ids)) {
        stop("Exon IDs don't match in ", count_files[i])
    }
    count_matrix[,i] <- counts[,2]
}

cat("\nCount matrix dimensions:", dim(count_matrix), "\n")

# Parse exon IDs to get gene and exon info
cat("\nParsing exon IDs...\n")
exon_info <- str_match(exon_ids, '^"([^"]+)":"([^"]+)"$')
gene_ids <- exon_info[,2]
exon_nums <- exon_info[,3]

# Check for parsing failures
if (any(is.na(gene_ids))) {
    cat("WARNING: Some exon IDs could not be parsed:\n")
    failed_ids <- exon_ids[is.na(gene_ids)]
    print(head(failed_ids, 10))
    stop("Exon ID parsing failed")
}

# Create unique exon IDs
unique_exon_ids <- paste0(gene_ids, ":E", exon_nums)

# Check for duplicates
if (any(duplicated(unique_exon_ids))) {
    cat("WARNING: Duplicate exon IDs found:\n")
    dups <- unique_exon_ids[duplicated(unique_exon_ids)]
    print(head(dups, 10))
    # Make them unique by adding a counter
    unique_exon_ids <- make.unique(unique_exon_ids, sep="_")
}

# Create feature data
feature_data <- data.frame(
    geneID = gene_ids,
    exonID = unique_exon_ids,
    stringsAsFactors = FALSE
)
rownames(feature_data) <- unique_exon_ids

cat("Feature data dimensions:", dim(feature_data), "\n")
cat("First few rows:\n")
print(head(feature_data))

# Set rownames for count matrix using unique exon IDs
rownames(count_matrix) <- unique_exon_ids

# Create DEXSeqDataSet manually
cat("\nCreating DEXSeqDataSet...\n")
dxd <- DEXSeqDataSet(
    countData = count_matrix,
    sampleData = samples,
    design = ~ sample + exon + condition:exon,
    featureID = feature_data$exonID,
    groupID = feature_data$geneID
)

cat("DEXSeqDataSet created successfully!\n")
print(dxd)

# Run DEXSeq analysis
cat("\nEstimating size factors...\n")
dxd <- estimateSizeFactors(dxd)

cat("Estimating dispersions...\n")
dxd <- estimateDispersions(dxd)

cat("Testing for DEU...\n")
dxd <- testForDEU(dxd)

cat("Estimating fold changes...\n")
dxd <- estimateExonFoldChanges(dxd, fitExpToVar = "condition")

# Save results
cat("Saving results...\n")
results <- DEXSeqResults(dxd)
write.csv(as.data.frame(results), "dexseq_results.csv")

# Plot target genes
cat("Generating plots...\n")
for (gene in target_genes) {
    pdf(paste0(gene, "_splicing.pdf"))
    tryCatch({
        plotDEXSeq(results, gene, displayTranscripts = TRUE)
    }, error = function(e) {
        plot.new()
        text(0.5, 0.5, paste("Gene", gene, "not found"))
    })
    dev.off()
}

# HTML report
cat("Generating HTML report...\n")
DEXSeqHTML(results, FDR = 0.1, path = "dexseq_report.html")

cat("\nDEXSeq analysis complete!\n")
cat("Results saved in: results/dexseq/\n")
