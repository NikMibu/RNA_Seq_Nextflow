#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
 * Pipeline parameters
 */
params.input = "data/samplesheet.csv"
params.outdir = "results"
params.hisat2_index = "reference/hisat2_index"
params.gtf = "reference/annotation/gencode.v43.annotation.gtf"
params.rmats_path = System.getenv("RMATS_PATH") ?: "/mnt/d/RNA_seq/rmats-turbo/rmats.py"

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
    
    // rMATS for differential splicing
    RMATS(
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

process RMATS {
    tag "rMATS_analysis"
    publishDir "${params.outdir}/rmats", mode: 'copy'
    cpus 8
    memory 16.GB
    
    input:
    path bams
    path gtf
    path samplesheet
    
    output:
    path "rmats_output/*"
    path "*.pdf"
    
    script:
    """
    #!/bin/bash
    set -e

    RMATS_CMD="python ${params.rmats_path}"
    GTF_PATH=\$(realpath ${gtf})

    # Read sample lists
    grep "SF3B1_mutant" ${samplesheet} | cut -d',' -f1 > mutant_samples.txt
    grep "SF3B1_wildtype" ${samplesheet} | cut -d',' -f1 > wildtype_samples.txt

    # Build BAM lists (direct matching, no subshell)
    > b1.txt
    while read sample; do
        for bam in \${sample}.bam; do
            if [ -f "\$bam" ]; then
                realpath "\$bam" >> b1.txt
            fi
        done
    done < mutant_samples.txt

    > b2.txt
    while read sample; do
        for bam in \${sample}.bam; do
            if [ -f "\$bam" ]; then
                realpath "\$bam" >> b2.txt
            fi
        done
    done < wildtype_samples.txt

    # Convert to comma-separated
    paste -sd, b1.txt > b1_list.txt
    paste -sd, b2.txt > b2_list.txt

    # Debug
    echo "b1_list.txt:"
    cat b1_list.txt
    echo "b2_list.txt:"
    cat b2_list.txt

    # Run rMATS
    \${RMATS_CMD} --b1 b1_list.txt \\
                --b2 b2_list.txt \\
                --gtf \${GTF_PATH} \\
                --od rmats_output \\
                --tmp rmats_tmp \\
                -t paired \\
                --readLength 202 \\
                --nthread ${task.cpus} \\
                --libType fr-unstranded \\
                --variable-read-length \\
                --allow-clipping
    
    # Plot target genes
    Rscript - <<'RSCRIPT'
    library(ggplot2)
    library(dplyr)
    
    target_genes <- c("ABCC5", "CRNDE", "UQCC", "GUSBP11", "ANKHD1", "ADAM12")
    
    # Read rMATS results
    for (event_type in c("RI", "SE", "A5SS", "A3SS", "MXE")) {
        file <- paste0("rmats_output/", event_type, ".MATS.JC.txt")
        if (file.exists(file)) {
            data <- read.table(file, header=TRUE, sep="\\t", stringsAsFactors=FALSE)
            
            # Filter for target genes
            data_target <- data %>% 
                filter(geneSymbol %in% target_genes, FDR < 0.05)
            
            if (nrow(data_target) > 0) {
                pdf(paste0(event_type, "_significant_events.pdf"))
                print(
                    ggplot(data_target, aes(x=geneSymbol, y=IncLevelDifference)) +
                    geom_bar(stat="identity") +
                    theme_minimal() +
                    labs(title=paste(event_type, "- Significant Events"),
                         x="Gene", y="Inclusion Level Difference") +
                    theme(axis.text.x = element_text(angle=45, hjust=1))
                )
                dev.off()
            }
        }
    }
    
    cat("rMATS analysis complete!\\n")
    RSCRIPT
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

