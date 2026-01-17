#!/bin/bash
# download_reference.sh

# Wir brauchen hier KEIN conda activate, da wget und gunzip System-Tools sind.

echo "=== Downloading Reference Genome & Annotation ==="

# Ordnerstruktur erstellen
mkdir -p reference/{genome,annotation,star_index}

# In den Ordner wechseln (und abbrechen, falls das nicht klappt)
cd reference || exit

# -------------------------------------------------------
# 1. Genome FASTA
# -------------------------------------------------------
if [[ ! -f "genome/GRCh38.primary_assembly.genome.fa" ]]; then
    echo "→ Downloading genome..."
    
    # -c erlaubt Fortsetzen bei Abbruch
    wget -c -P genome/ \
        https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/GRCh38.primary_assembly.genome.fa.gz
    
    # Entpacken (nur wenn Download erfolgreich war)
    if [[ -f "genome/GRCh38.primary_assembly.genome.fa.gz" ]]; then
        gunzip genome/GRCh38.primary_assembly.genome.fa.gz
        echo "✓ Genome downloaded and unzipped"
    else
        echo "❌ Fehler beim Download des Genoms"
    fi
else
    echo "✓ Genome already exists"
fi

# -------------------------------------------------------
# 2. GTF Annotation
# -------------------------------------------------------
if [[ ! -f "annotation/gencode.v43.annotation.gtf" ]]; then
    echo "→ Downloading annotation..."
    
    wget -c -P annotation/ \
        https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/gencode.v43.annotation.gtf.gz
    
    if [[ -f "annotation/gencode.v43.annotation.gtf.gz" ]]; then
        gunzip annotation/gencode.v43.annotation.gtf.gz
        echo "✓ Annotation downloaded and unzipped"
    else
        echo "❌ Fehler beim Download der Annotation"
    fi
else
    echo "✓ Annotation already exists"
fi

# Zurück zum Hauptordner
cd ..

echo "✓ Reference download complete!"