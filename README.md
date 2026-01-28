# Scientific context of the Project
## What is uveal melanoma?
Uveal melanoma is a rare but aggressive cancer of the eye. More precisely what “uveal” means:
The uvea is the middle, pigmented layer of the eye, and it has three parts:
- Iris (colored part at the front)
- Ciliary body
- Choroid (vascular layer at the back of the eye)
  
Uveal melanoma arises from melanocytes (pigment-producing cells) located in one of these uveal tissues—most commonly the choroid.

### Key characteristics
It is the most common primary intraocular cancer in adults distinct from cutaneous (skin) melanoma:
- Different mutation spectrum
- Different biology
- Different clinical behavior
Often diagnosed by eye examination rather than biopsy

### Genetics & molecular biology (why it matters for our project)
Uveal melanoma is genetically quite simple but very specific:

Early “driver” mutations (usually mutually exclusive):
- GNAQ
- GNA11

Progression / prognostic mutations:
- BAP1 → associated with metastasis (poor prognosis)
- SF3B1 → associated with late-onset metastasis
- EIF1AX → generally better prognosis
  
This is why SF3B1 is so interesting:
it encodes a core splicing factor, making uveal melanoma a natural model to study mutation-driven splicing dysregulation.

### Clinical relevance
- Primary tumor can often be treated locally (radiotherapy, surgery)
- ~50% of patients develop metastases, mainly to the liver
- Once metastatic, prognosis is poor
- Molecular profiling is used for risk stratification

## What is SF3B1?
SF3B1 is a core component of the RNA splicing machinery, and one of the most frequently mutated splicing factors in human cancer.

### What SF3B1 does (normal function)
SF3B1 stands for splicing factor 3B subunit 1.
It is part of the U2 small nuclear ribonucleoprotein (U2 snRNP), a key complex of the spliceosome.
In simple terms:
- Genes are transcribed into pre-mRNA containing exons and introns
- The spliceosome removes introns and joins exons
- SF3B1 helps the spliceosome recognize the correct 3′ splice site
  
Specifically, SF3B1:
- Binds near the branch point sequence
- Stabilizes U2 snRNP binding
- Ensures accurate selection of intron–exon boundaries.
  
Without SF3B1, splicing would be inaccurate or fail.

### SF3B1 mutations in cancer
SF3B1 is not randomly mutated. Cancer-associated mutations are:
- Recurrent
- Heterozygous
- Missense
- Clustered at specific residues
  
In uveal melanoma:
- Hotspot mutation: R625 (codon 625)
- Found in ~15–25% of tumors
- Associated with late-onset metastasis
  
In other cancers:
- Myelodysplastic syndromes (K700E hotspot)
- Chronic lymphocytic leukemia
- Breast cancer
- Pancreatic cancer
  
The same protein, but different hotspots in different cancers.

### What SF3B1 mutations do to splicing
Mutant SF3B1 does not shut down splicing. Instead, it causes systematic errors:
- Preferential use of cryptic / alternative 3′ splice sites
- Typically 10–30 nucleotides upstream of the canonical site
- Leads to:
  - Exon truncation
  - Frameshifts
  - Premature stop codons
  - Altered protein isoforms
    
This explains why:
- Total gene expression may look normal
- But isoform composition is altered
  
This is exactly why early RNA-seq analyses often missed the effect.

### Why SF3B1 is central to our RNA-seq project
SF3B1 is:
- A direct molecular link between mutation and transcriptome phenotype
- One of the cleanest splicing-factor mutation models in cancer
- A perfect test case for:
  - Junction-level analysis
  - PSI-based methods
  - Annotation-aware vs annotation-free tools
  
This explains the discrepancy:
- Paper [1]: splicing-aware reanalysis → clear differences
- Paper [2]: gene-level focus → no splicing difference


# Experimental design
### Number of samples
- After quality control, 8 RNA-seq samples were retained for analysis, including n₁ SF3B1-mutant and n₂ SF3B1–wild-type tumors.
### SF3B1-mutated vs wild-type
- Samples were classified as SF3B1-mutant if they carried a recurrent missense mutation at codon 625; all other samples were considered SF3B1–wild-type.
### RNA-seq type (critical for splicing)
- RNA sequencing was performed using paired-end reads. Library preparation was (stranded / unstranded). These characteristics were taken into account when selecting splicing analysis tools.
### Any confounders we control for (batch, sex, etc.)
- Potential confounders including sequencing batch and patient sex were evaluated. Where metadata was available, batch was included as a covariate in downstream analyses; otherwise, exploratory analyses were used to assess its impact.


# Aim of the analysis
### What exactly are we trying to demonstrate with this analysis, and why is it worth doing now?
- The aim of this analysis is to reassess the impact of SF3B1 mutations on alternative splicing in uveal melanoma using modern RNA-seq processing and splicing-aware analysis methods. Specifically, we seek to reproduce previously reported SF3B1-associated splicing events, evaluate the robustness of these findings with current tools, and characterize the nature of the resulting splicing alterations.



# RNA-Seq pipeline for uveal melanoma

Nextflow pipeline for analyzing RNA-Seq data with a focus on differential splicing analysis in SF3B1 mutations.

## Overview

This pipeline performs the following steps:
1. **Quality Control** - FastQC for FASTQ-Files
2. **Alignment** - HISAT2 for read mapping
3. **Quantification** - featureCounts for gene-level counts
4. **Splicing analysis** - rMATS for differential splicing events

## Requirements

- Micromamba/Conda
- Nextflow
- ~100 GB free storage space

## Installation

### 1. Environment setup

Create the Conda/micromamba environment with all the necessary tools:

```bash
micromamba env create -f environment.yml
```

Activate the environment:

```bash
micromamba activate uveal-melanoma
```

Test the installation:

```bash
nextflow -version
fastqc --version
hisat2 --version
```

### 2. Download files

**Important:** Every scripts has to be executed from the project root directory!

#### FASTQ-files (RNA-Seq data)

```bash
./scripts/01_download_sra.sh
```

This downloads 8 samples (~4 GB compressed) and automatically creates `data/samplesheet.csv`.

#### Reference genome and annotation

```bash
./scripts/02_download_reference.sh
```

Downloads:
- Human reference genome (GRCh38)
- GENCODE annotation (v43)

#### HISAT2 index creation

```bash
./scripts/03_build_hisat_index.sh
```

Creates HISAT2 index (~8GB). The index build takes ~30-60 minutes.

#### rMATS installation

```bash
./scripts/04_install_rmats.sh
```

installs rMATS-turbo from GitHub and creates a Symlink. The dependencies (Cython, GSL, GCC, etc.) have to be already installed over the `environment.yml`.

**Note:** Building takes ~5-10 minutes.

**Important - Configure rMATS Path:**

The script installs rMATS to `~/rmats-turbo/`. The pipeline needs to know where rMATS is located:

**Option 1 (recommended):** Set environment variable:
```bash
export RMATS_PATH=~/rmats-turbo/rmats.py
```

**Option 2:** Adjust the default path in `nextflow.config`:
```groovy
params {
    rmats_path = "${System.getProperty('user.home')}/rmats-turbo/rmats.py"
}
```

## Run pipeline

### Standard run

```bash
nextflow run main.nf
```

### With resume (after error/interruption)

```bash
nextflow run main.nf -resume
```

Nextflow caches successfully completed jobs and only restarts failed ones.

## Configuration

Adjustments in `nextflow.config`:

```groovy
params {
    input = 'data/samplesheet.csv'
    outdir = 'results'
    hisat2_index = 'reference/hisat2_index'
    gtf = 'reference/annotation/gencode.v43.annotation.gtf'
}

process {
    cpus = 4          // Standard-CPUs
    memory = 16.GB    // Standard-RAM
    
    withName: HISAT2_ALIGN {
        cpus = 8
        memory = 12.GB
    }
}
```

## Results

After a successful run, you will find the results in `results/`:

```
results/
├── fastqc/              # QC reports (HTML)
├── hisat2/              # BAM-Files and alignment logs
├── counts/              # featureCounts tables
└── rmats/
    ├── rmats_output/            # rMATS results
    │   ├── RI.MATS.JC.txt       # Retained introns
    │   ├── SE.MATS.JC.txt       # Skipped exons
    │   ├── A3SS.MATS.JC.txt     # Alternative 3' splice sites
    │   ├── A5SS.MATS.JC.txt     # Alternative 5' splice sites
    │   └── MXE.MATS.JC.txt      # Mutually exclusive exons
    └── *_significant_events.pdf # Plots for significant events
```

### Important Files

- `RI.MATS.JC.txt` - Retained intron events (z.B. ABCC5)
- `SE.MATS.JC.txt` - Skipped exon events
- `A3SS.MATS.JC.txt` - Alternative 3' splice site events (z.B. CRNDE)
- `A5SS.MATS.JC.txt` - Alternative 5' splice site events
- `MXE.MATS.JC.txt` - Mutually exclusive exon events
- `*_significant_events.pdf` - Plots for significant events (FDR < 0.05) in target genes:
  - ABCC5, CRNDE, UQCC
  - GUSBP11, ANKHD1, ADAM12

## Further analyses

Additional analyses can be performed after the pipeline. The scripts are located in `analysis_scripts/`:

### 1. Detailed rMATS plots

Creates advanced visualizations of splicing events:

```bash
cd results/rmats/rmats_output
Rscript ../../../analysis_scripts/improved_plots_rMATS.R
```

**Output:**
- `signifikante_events_count.pdf` - Overview of all significant events per gene
- `signifikante_events_table.txt` - Table with all events
- `GENE_events.pdf` - Detailed plots for each affected gene (ABCC5, CRNDE, ANKHD1)

### 2. PCA analysis

Principal component analysis of gene expression:

```bash
Rscript analysis_scripts/PCA.R
```

**Output:** `results/PCA_SF3B1.pdf` - PCA-plot shows separation between SF3B1 mutants and wild type

### 3. Venn diagram

Comparison of the genes found with literature (Furney et al.):

```bash
Rscript analysis_scripts/venn_diagramm.R
```

**Output:** `results/Venn_Splicing_Genes.pdf` - Overlaps between own analysis and the paper

**Hinweis:** The scripts use absolute paths. If necessary, the `setwd()` or paths must be adjusted.

## Troubleshooting

### Pipeline does not start

**Issue:** `hisat2: command not found`

**Solution:** Activate conda environment:
```bash
micromamba activate uveal-melanoma
```

### Resume not working

**Issue:** Nextflow restarts pipeline from the beginning

**Solution:** 
- Use `-resume` flag
- Check whether the `.nextflow/` directory exists
- If problems occur: `rm -rf .nextflow/cache_backup/`

### Memory issues

**Issue:** `OutOfMemoryError` or jobs are terminated

**Solution:** Reduce parallelism in `nextflow.config`:
```groovy
process {
    cpus = 2
    memory = 8.GB
}
```

### Damaged FASTQ files

**Issue:** `gzip: unexpected end of file`

**Solution:** Download the file again:
```bash
rm data/fastq/SRR628XXX_1.fastq.gz
wget -P data/fastq/ ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR628/SRR628XXX/SRR628XXX_1.fastq.gz
```

## Project structure

```
RNA_seq/
├── main.nf                  # Nextflow pipeline
├── nextflow.config          # Configurations
├── environment.yml          # Conda environment
├── README.md               # This file
├── scripts/
│   ├── 01_download_sra.sh
│   ├── 02_download_reference.sh
│   ├── 03_build_hisat_index.sh
│   └── 04_install_rmats.sh
├── data/
│   ├── samplesheet.csv     # Sample metadaten
│   └── fastq/              # FASTQ files
├── reference/
│   ├── genome/             # Referenz genome
│   ├── annotation/         # GTF annotation
│   └── hisat2_index/       # HISAT2 index
├── results/                # Pipeline outputs
└── work/                   # Nextflow work directory (temp)
```

## Samples

The pipeline analyzes 8 uveal melanoma samples:

| Sample ID  | Condition       |
|-----------|-----------------|
| SRR628582 | SF3B1_mutant   |
| SRR628583 | SF3B1_mutant   |
| SRR628584 | SF3B1_mutant   |
| SRR628585 | SF3B1_wildtype |
| SRR628586 | SF3B1_wildtype |
| SRR628587 | SF3B1_wildtype |
| SRR628588 | SF3B1_wildtype |
| SRR628589 | SF3B1_wildtype |

## Resource requirements

### Minimum requirements
- 4 CPUs
- 16 GB RAM
- 100 GB storage

### Recommended
- 8+ CPUs
- 32 GB RAM
- 150 GB storage

### Runtime
- Complete pipeline: ~2-3 hours (8 CPUs)
- With resume: Depends on failed jobs

## References

- **HISAT2:** Kim et al. (2019) - Graph-based genome alignment
- **featureCounts:** Liao et al. (2014) - Read summarization
- **rMATS:** Shen et al. (2014) - rMATS: robust and flexible detection of differential alternative splicing
- **Nextflow:** Di Tommaso et al. (2017) - Workflow management


