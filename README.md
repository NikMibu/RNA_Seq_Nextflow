# Scientific context of the Project
## What is uveal melanoma?
Uveal melanoma is a rare but aggressive cancer of the eye. More precisely what “uveal” means:
The uvea is the middle, pigmented layer of the eye, and it has three parts:
- Iris (colored part at the front)
- Ciliary body
- Choroid (vascular layer at the back of the eye)
  
Uveal melanoma arises from melanocytes (pigment-producing cells) located in one of these uveal tissues—most commonly the choroid.

### Key characteristics
It is the most common primary intraocular cancer in adults
Distinct from cutaneous (skin) melanoma:
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
SF3B1 stands for Splicing Factor 3B Subunit 1.
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
  
In uveal melanoma
- Hotspot mutation: R625 (codon 625)
- Found in ~15–25% of tumors
- Associated with late-onset metastasis
  
In other cancers
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
  - junction-level analysis
  - PSI-based methods
  - annotation-aware vs annotation-free tools
  
This explains the discrepancy:
- Paper [1]: splicing-aware reanalysis → clear differences
- Paper [2]: gene-level focus → no splicing difference


# Experimental design
### Number of samples
- After quality control, 8 RNA-seq samples were retained for analysis, including n₁ SF3B1-mutant and n₂ SF3B1–wild-type tumors.
### SF3B1-mutated vs wild-type
- Samples were classified as SF3B1-mutant if they carried a recurrent missense mutation at codon 625; all other samples were considered SF3B1–wild-type.
### RNA-seq type (critical for splicing)
- RNA sequencing was performed using paired-end reads of approximately X bp. Library preparation was (stranded / unstranded). These characteristics were taken into account when selecting splicing analysis tools.
### Any confounders we control for (batch, sex, etc.)
- Potential confounders including sequencing batch and patient sex were evaluated. Where metadata was available, batch was included as a covariate in downstream analyses; otherwise, exploratory analyses were used to assess its impact.


# Aim of the analysis
### What exactly are we trying to demonstrate with this analysis, and why is it worth doing now?
- The aim of this analysis is to reassess the impact of SF3B1 mutations on alternative splicing in uveal melanoma using modern RNA-seq processing and splicing-aware analysis methods. Specifically, we seek to reproduce previously reported SF3B1-associated splicing events, evaluate the robustness of these findings with current tools, and characterize the nature of the resulting splicing alterations.



# RNA-Seq Pipeline für Uveal Melanoma

Nextflow-Pipeline zur Analyse von RNA-Seq Daten mit Fokus auf differentielle Splicing-Analyse bei SF3B1-Mutationen.

## Übersicht

Diese Pipeline führt folgende Schritte aus:
1. **Quality Control** - FastQC für FASTQ-Dateien
2. **Alignment** - HISAT2 für Read-Mapping
3. **Quantifizierung** - featureCounts für Gen-Level Counts
4. **Splicing-Analyse** - rMATS für differentielle Splicing-Events

## Voraussetzungen

- Micromamba/Conda
- Nextflow
- ~100 GB freier Speicherplatz

## Installation

### 1. Environment Setup

Erstelle die Conda/micromamba-Umgebung mit allen benötigten Tools:

```bash
micromamba env create -f environment.yml
```

Aktiviere die Umgebung:

```bash
micromamba activate uveal-melanoma
```

Teste die Installation:

```bash
nextflow -version
fastqc --version
hisat2 --version
```

### 2. Daten herunterladen

**Wichtig:** Alle Scripts müssen aus dem Projekt-Root-Verzeichnis ausgeführt werden!

#### FASTQ-Dateien (RNA-Seq Daten)

```bash
./scripts/01_download_sra.sh
```

Dies lädt 8 Samples herunter (~4 GB komprimiert) und erstellt automatisch `data/samplesheet.csv`.

#### Referenz-Genom und Annotation

```bash
./scripts/02_download_reference.sh
```

Downloads:
- Human Referenz-Genom (GRCh38)
- GENCODE Annotation (v43)

#### HISAT2-Index erstellen

```bash
./scripts/03_build_hisat_index.sh
```

Erstellt den HISAT2-Index (~8 GB). Der Index-Build dauert ~30-60 Minuten.

#### rMATS installieren

```bash
./scripts/04_install_rmats.sh
```

Installiert rMATS-turbo von GitHub und erstellt einen Symlink. Die Dependencies (Cython, GSL, GCC, etc.) müssen bereits über `environment.yml` installiert sein.

**Hinweis:** Der Build dauert ~5-10 Minuten.

**Wichtig - rMATS Path konfigurieren:**

Das Script installiert rMATS nach `~/rmats-turbo/`. Die Pipeline muss wissen, wo rMATS liegt:

**Option 1 (empfohlen):** Environment Variable setzen:
```bash
export RMATS_PATH=~/rmats-turbo/rmats.py
```

**Option 2:** In `nextflow.config` den Default-Path anpassen:
```groovy
params {
    rmats_path = "${System.getProperty('user.home')}/rmats-turbo/rmats.py"
}
```

## Pipeline ausführen

### Standard-Run

```bash
nextflow run main.nf
```

### Mit Resume (nach Fehler/Unterbrechung)

```bash
nextflow run main.nf -resume
```

Nextflow cached erfolgreich abgeschlossene Jobs und startet nur fehlgeschlagene neu.

## Konfiguration

Anpassungen in `nextflow.config`:

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

## Ergebnisse

Nach erfolgreichem Run findest du die Ergebnisse in `results/`:

```
results/
├── fastqc/              # QC-Reports (HTML)
├── hisat2/              # BAM-Dateien und Alignment-Logs
├── counts/              # featureCounts Tabellen
└── rmats/
    ├── rmats_output/            # rMATS Ergebnisse
    │   ├── RI.MATS.JC.txt       # Retained Introns
    │   ├── SE.MATS.JC.txt       # Skipped Exons
    │   ├── A3SS.MATS.JC.txt     # Alternative 3' Splice Sites
    │   ├── A5SS.MATS.JC.txt     # Alternative 5' Splice Sites
    │   └── MXE.MATS.JC.txt      # Mutually Exclusive Exons
    └── *_significant_events.pdf # Plots für signifikante Events
```

### Wichtige Dateien

- `RI.MATS.JC.txt` - Retained Intron Events (z.B. ABCC5)
- `SE.MATS.JC.txt` - Skipped Exon Events
- `A3SS.MATS.JC.txt` - Alternative 3' Splice Site Events (z.B. CRNDE)
- `A5SS.MATS.JC.txt` - Alternative 5' Splice Site Events
- `MXE.MATS.JC.txt` - Mutually Exclusive Exon Events
- `*_significant_events.pdf` - Plots für signifikante Events (FDR < 0.05) in Target-Genen:
  - ABCC5, CRNDE, UQCC
  - GUSBP11, ANKHD1, ADAM12

## Weiterführende Analysen

Nach der Pipeline können zusätzliche Analysen durchgeführt werden. Die Scripts befinden sich in `analysis_scripts/`:

### 1. Detaillierte rMATS Plots

Erstellt erweiterte Visualisierungen der Splicing-Events:

```bash
cd results/rmats/rmats_output
Rscript ../../../analysis_scripts/improved_plots_rMATS.R
```

**Output:**
- `signifikante_events_count.pdf` - Übersicht aller signifikanten Events pro Gen
- `signifikante_events_table.txt` - Tabelle mit allen Events
- `GENE_events.pdf` - Detaillierte Plots für jedes betroffene Gen (ABCC5, CRNDE, ANKHD1)

### 2. PCA-Analyse

Principal Component Analysis der Genexpression:

```bash
Rscript analysis_scripts/PCA.R
```

**Output:** `results/PCA_SF3B1.pdf` - PCA-Plot zeigt Separation zwischen SF3B1-Mutanten und Wildtyp

### 3. Venn-Diagramm

Vergleich der gefundenen Gene mit Literatur (Furney et al.):

```bash
Rscript analysis_scripts/venn_diagramm.R
```

**Output:** `results/Venn_Splicing_Genes.pdf` - Überlappung zwischen eigener Analyse und Paper

**Hinweis:** Die Scripts verwenden absolute Pfade. Bei Bedarf müssen die `setwd()` oder Pfade angepasst werden.

## Troubleshooting

### Pipeline startet nicht

**Problem:** `hisat2: command not found`

**Lösung:** Aktiviere die Conda-Umgebung:
```bash
micromamba activate uveal-melanoma
```

### Resume funktioniert nicht

**Problem:** Nextflow startet Pipeline von vorne

**Lösung:** 
- Verwende `-resume` Flag
- Prüfe, ob `.nextflow/` Verzeichnis existiert
- Bei Problemen: `rm -rf .nextflow/cache_backup/`

### Speicherprobleme

**Problem:** `OutOfMemoryError` oder Jobs werden gekillt

**Lösung:** Reduziere Parallelität in `nextflow.config`:
```groovy
process {
    cpus = 2
    memory = 8.GB
}
```

### Beschädigte FASTQ-Dateien

**Problem:** `gzip: unexpected end of file`

**Lösung:** Datei neu herunterladen:
```bash
rm data/fastq/SRR628XXX_1.fastq.gz
wget -P data/fastq/ ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR628/SRR628XXX/SRR628XXX_1.fastq.gz
```

## Projekt-Struktur

```
RNA_seq/
├── main.nf                  # Nextflow-Pipeline
├── nextflow.config          # Konfiguration
├── environment.yml          # Conda-Environment
├── README.md               # Diese Datei
├── scripts/
│   ├── 01_download_sra.sh
│   ├── 02_download_reference.sh
│   ├── 03_build_hisat_index.sh
│   └── 04_install_rmats.sh
├── data/
│   ├── samplesheet.csv     # Sample-Metadaten
│   └── fastq/              # FASTQ-Dateien
├── reference/
│   ├── genome/             # Referenz-Genom
│   ├── annotation/         # GTF-Annotation
│   └── hisat2_index/       # HISAT2-Index
├── results/                # Pipeline-Outputs
└── work/                   # Nextflow Work-Dir (temp)
```

## Samples

Die Pipeline analysiert 8 Uveal Melanoma Samples:

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

## Ressourcen-Anforderungen

### Minimale Anforderungen
- 4 CPUs
- 16 GB RAM
- 100 GB Speicher

### Empfohlen
- 8+ CPUs
- 32 GB RAM
- 150 GB Speicher

### Laufzeit
- Komplette Pipeline: ~2-3 Stunden (8 CPUs)
- Mit Resume: Abhängig von fehlgeschlagenen Jobs

## Referenzen

- **HISAT2:** Kim et al. (2019) - Graph-based genome alignment
- **featureCounts:** Liao et al. (2014) - Read summarization
- **rMATS:** Shen et al. (2014) - rMATS: robust and flexible detection of differential alternative splicing
- **Nextflow:** Di Tommaso et al. (2017) - Workflow management


