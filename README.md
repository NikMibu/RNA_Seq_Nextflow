
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


