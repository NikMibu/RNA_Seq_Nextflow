
# RNA-Seq Pipeline für Uveal Melanoma

Nextflow-Pipeline zur Analyse von RNA-Seq Daten mit Fokus auf differentielle Splicing-Analyse bei SF3B1-Mutationen.

## Übersicht

Diese Pipeline führt folgende Schritte aus:
1. **Quality Control** - FastQC für FASTQ-Dateien
2. **Alignment** - HISAT2 für Read-Mapping
3. **Quantifizierung** - featureCounts für Gen-Level Counts
4. **Splicing-Analyse** - DEXSeq für differentielle Exon-Nutzung (von AI Erstellt)

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


Downloads:
- Human Referenz-Genom (GRCh38)
- GENCODE Annotation (v43)
- Erstellt HISAT2-Index (~8 GB)

**Hinweis:** Der Index-Build dauert ~30-60 Minuten.

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
└── dexseq/
    ├── dexseq_results.csv       # Differentielle Exon-Nutzung
    ├── dexseq_report.html       # Interaktiver Report
    └── *_splicing.pdf           # Plots für Target-Gene
```

### Wichtige Dateien

- `dexseq_results.csv` - Vollständige DEXSeq-Ergebnisse mit p-Werten
- `dexseq_report.html` - Interaktiver HTML-Report
- `GENE_splicing.pdf` - Splicing-Plots für Target-Gene:
  - ABCC5, CRNDE, UQCC
  - GUSBP11, ANKHD1, ADAM12

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
│   └── 02_download_reference.sh
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
- **DEXSeq:** Anders et al. (2012) - Differential exon usage
- **Nextflow:** Di Tommaso et al. (2017) - Workflow management

=======
# RNA-Seq Pipeline für Uveal Melanoma

Nextflow-Pipeline zur Analyse von RNA-Seq Daten mit Fokus auf differentielle Splicing-Analyse bei SF3B1-Mutationen.

## Übersicht

Diese Pipeline führt folgende Schritte aus:
1. **Quality Control** - FastQC für FASTQ-Dateien
2. **Alignment** - HISAT2 für Read-Mapping
3. **Quantifizierung** - featureCounts für Gen-Level Counts
4. **Splicing-Analyse** - DEXSeq für differentielle Exon-Nutzung

## Voraussetzungen

- Micromamba/Conda
- Nextflow
- ~100 GB freier Speicherplatz

## Installation

### 1. Environment Setup

Erstelle die Conda-Umgebung mit allen benötigten Tools:

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

#### FASTQ-Dateien (RNA-Seq Daten)

```bash
cd scripts
./01_download_sra.sh
```

Dies lädt 8 Samples herunter (~4 GB komprimiert) und erstellt automatisch `data/samplesheet.csv`.

#### Referenz-Genom und Annotation

```bash
./02_download_reference.sh
```

Downloads:
- Human Referenz-Genom (GRCh38)
- GENCODE Annotation (v43)
- Erstellt HISAT2-Index (~8 GB)

**Hinweis:** Der Index-Build dauert ~30-60 Minuten.

## Pipeline ausführen

### Standard-Run

```bash
nextflow run main.nf
```

### Mit Resume (nach Fehler/Unterbrechung)

```bash
nextflow run main.nf --resume
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
└── dexseq/
    ├── dexseq_results.csv       # Differentielle Exon-Nutzung
    ├── dexseq_report.html       # Interaktiver Report
    └── *_splicing.pdf           # Plots für Target-Gene
```

### Wichtige Dateien

- `dexseq_results.csv` - Vollständige DEXSeq-Ergebnisse mit p-Werten
- `dexseq_report.html` - Interaktiver HTML-Report
- `GENE_splicing.pdf` - Splicing-Plots für Target-Gene:
  - ABCC5, CRNDE, UQCC
  - GUSBP11, ANKHD1, ADAM12

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
- Verwende `--resume` Flag
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
│   └── 02_download_reference.sh
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
- **DEXSeq:** Anders et al. (2012) - Differential exon usage
- **Nextflow:** Di Tommaso et al. (2017) - Workflow management


