#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
 * Pipeline parameters
 */
params.input = "data/samplesheet.csv"
params.outdir = "results"
params.hisat2_index = "reference/hisat2_index"
params.gtf = "reference/annotation/gencode.v43.annotation.gtf"

log.info """\
    RNA-SEQ PIPELINE
    ================
    input      : ${params.input}
    outdir     : ${params.outdir}
    hisat2_index : ${params.hisat2_index}
    gtf        : ${params.gtf}
    """
    .stripIndent()

/*
 * Main workflow
 */
workflow {
    
    // Convert paths to file objects
    hisat2_index_path = file(params.hisat2_index)
    gtf_file = file(params.gtf)
    samplesheet_file = file(params.input)
    
    // Validate paths
    if (!hisat2_index_path.exists()) {
        error "HISAT2 index not found: ${params.hisat2_index}"
    }
    if (!gtf_file.exists()) {
        error "GTF file not found: ${params.gtf}"
    }
    if (!samplesheet_file.exists()) {
        error "Samplesheet not found: ${params.input}"
    }
    
    // Read samplesheet
    Channel
        .fromPath(params.input)
        .splitCsv(header: true)
        .map { row -> 
            tuple(
                row.sample_id, 
                row.condition,
                file(row.fastq_1), 
                file(row.fastq_2)
            )
        }
        .set { samples_ch }
    
    // Quality control
    FASTQC(samples_ch)
    
    // Alignment with HISAT2
    HISAT2_ALIGN(samples_ch, hisat2_index_path)
    
    // Count reads per gene (runs once per sample)
    FEATURECOUNTS(HISAT2_ALIGN.out.bam, gtf_file)
    
    // DEXSeq for differential splicing
    DEXSEQ_ANALYSIS(
        HISAT2_ALIGN.out.bam.map { sample_id, condition, bam -> bam }.collect(),
        gtf_file,
        samplesheet_file
    )
}

/*
 * Processes
 */

process FASTQC {
    tag "$sample_id"
    publishDir "${params.outdir}/fastqc", mode: 'copy'
    
    input:
    tuple val(sample_id), val(condition), path(read1), path(read2)
    
    output:
    path "*.html"
    path "*.zip"
    
    script:
    """
    fastqc -t 2 ${read1} ${read2}
    """
}

process HISAT2_ALIGN {
    tag "$sample_id"
    publishDir "${params.outdir}/hisat2", mode: 'copy'
    
    input:
    tuple val(sample_id), val(condition), path(read1), path(read2)
    path index_dir
    
    output:
    tuple val(sample_id), val(condition), path("${sample_id}.bam"), emit: bam
    path "${sample_id}.log"
    
    script:
    """
    hisat2 -p ${task.cpus} \\
           -x ${index_dir}/genome \\
           -1 ${read1} \\
           -2 ${read2} \\
           --summary-file ${sample_id}.log \\
           | samtools sort -o ${sample_id}.bam -

    samtools index ${sample_id}.bam
    """
}

process FEATURECOUNTS {
    tag "$sample_id"
    publishDir "${params.outdir}/counts", mode: 'copy'
    
    input:
    tuple val(sample_id), val(condition), path(bam)
    path gtf
    
    output:
    path "${sample_id}.counts.txt"
    
    script:
    """
    featureCounts -p -T ${task.cpus} \\
                  -a ${gtf} \\
                  -o ${sample_id}.counts.txt \\
                  ${bam}
    """
}

process DEXSEQ_ANALYSIS {
    publishDir "${params.outdir}/dexseq", mode: 'copy'
    
    input:
    path bams
    path gtf
    path samplesheet
    
    output:
    path "dexseq_results.csv"
    path "*.pdf"
    path "dexseq_report.html"
    
    script:
    """
    #!/usr/bin/env Rscript
    
    library(DEXSeq)
    library(tidyverse)
    
    # Target genes from paper
    target_genes <- c("ABCC5", "CRNDE", "UQCC", 
                      "GUSBP11", "ANKHD1", "ADAM12")
    
    # Prepare GTF for DEXSeq
    dexseq_scripts <- file.path(Sys.getenv("CONDA_PREFIX"), "lib/R/library/DEXSeq/python_scripts")
    system(paste("python", file.path(dexseq_scripts, "dexseq_prepare_annotation.py"), "${gtf}", "dexseq.gff"))
    
    # Count reads per exon for each BAM
    bam_files <- list.files(pattern = "\\\\.bam\$")
    
    for (bam in bam_files) {
        sample_id <- gsub("\\\\.Aligned.*", "", bam)
        sample_id <- gsub("\\\\.bam\$", "", sample_id)
        count_file <- paste0(sample_id, ".txt")
        cmd <- paste(
            "python", file.path(dexseq_scripts, "dexseq_count.py"),
            "-p yes -r pos -s no -f bam",
            "dexseq.gff", bam, count_file
        )
        system(cmd)
        
        # Remove meta lines (they cause parsing errors)
        system(paste("grep -v '^_' ", count_file, "> temp.txt && mv temp.txt", count_file))
    }
    
    # Read sample info
    samples <- read.csv("${samplesheet}")
    count_files <- paste0(samples\$sample_id, ".txt")
    
    # Create DEXSeq dataset
    dxd <- DEXSeqDataSetFromHTSeq(
        count_files,
        sampleData = samples,
        design = ~ sample + exon + condition:exon,
        flattenedfile = "dexseq.gff"
    )
    
    # Run DEXSeq
    dxd <- estimateSizeFactors(dxd)
    dxd <- estimateDispersions(dxd)
    dxd <- testForDEU(dxd)
    dxd <- estimateExonFoldChanges(dxd, fitExpToVar = "condition")
    
    # Results
    results <- DEXSeqResults(dxd)
    write.csv(as.data.frame(results), "dexseq_results.csv")
    
    # Plot target genes
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
    DEXSeqHTML(results, FDR = 0.1, path = "dexseq_report.html")
    
    cat("DEXSeq analysis complete!\\n")
    """
}

workflow.onComplete {
    log.info """\
        Pipeline completed!
        -------------------
        Status   : ${workflow.success ? 'SUCCESS' : 'FAILED'}
        Duration : ${workflow.duration}
        Results  : ${params.outdir}
        """
        .stripIndent()
}

