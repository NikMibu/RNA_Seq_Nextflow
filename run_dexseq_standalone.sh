#!/bin/bash

# Standalone DEXSeq analysis using existing BAM files
# Run from project root: bash run_dexseq_standalone.sh

set -e

echo "Starting DEXSeq analysis with existing data..."

# Create output directory
mkdir -p results/dexseq
cd results/dexseq

# Copy BAM files from results
echo "Linking BAM files..."
for bam in ../hisat2/*.bam; do
    ln -sf "$bam" .
done

# Copy GTF and samplesheet
cp ../../reference/annotation/gencode.v43.annotation.gtf .
cp ../../data/samplesheet.csv .

# Run R script
echo "Running DEXSeq analysis..."
Rscript - <<'EOF'
library(DEXSeq)
library(tidyverse)

# Target genes from paper
target_genes <- c("ABCC5", "CRNDE", "UQCC", 
                  "GUSBP11", "ANKHD1", "ADAM12")

# Prepare GTF for DEXSeq
dexseq_scripts <- file.path(Sys.getenv("CONDA_PREFIX"), "lib/R/library/DEXSeq/python_scripts")
system(paste("python", file.path(dexseq_scripts, "dexseq_prepare_annotation.py"), "gencode.v43.annotation.gtf", "dexseq.gff"))

# Count reads per exon for each BAM
bam_files <- list.files(pattern = "\\.bam$")

for (bam in bam_files) {
    sample_id <- gsub("\\.Aligned.*", "", bam)
    sample_id <- gsub("\\.bam$", "", sample_id)
    count_file <- paste0(sample_id, ".txt")
    cmd <- paste(
        "python", file.path(dexseq_scripts, "dexseq_count.py"),
        "-p yes -r pos -s no -f bam",
        "dexseq.gff", bam, count_file
    )
    cat("Processing", bam, "...\n")
    system(cmd)
    
    # Remove meta lines (they cause parsing errors)
    system(paste("grep -v '^_' ", count_file, "> temp.txt && mv temp.txt", count_file))
}

# Read sample info
samples <- read.csv("samplesheet.csv")
count_files <- paste0(samples$sample_id, ".txt")

# Create DEXSeq dataset
cat("Creating DEXSeq dataset...\n")
dxd <- DEXSeqDataSetFromHTSeq(
    count_files,
    sampleData = samples,
    design = ~ sample + exon + condition:exon,
    flattenedfile = "dexseq.gff"
)

# Run DEXSeq
cat("Running DEXSeq analysis...\n")
dxd <- estimateSizeFactors(dxd)
dxd <- estimateDispersions(dxd)
dxd <- testForDEU(dxd)
dxd <- estimateExonFoldChanges(dxd, fitExpToVar = "condition")

# Results
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

cat("DEXSeq analysis complete!\n")
cat("Results saved in: results/dexseq/\n")
EOF

echo ""
echo "Done! Check results/dexseq/ for output files:"
echo "  - dexseq_results.csv"
echo "  - *_splicing.pdf"
echo "  - dexseq_report.html"
