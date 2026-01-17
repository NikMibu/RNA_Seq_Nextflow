#!/bin/bash
# download_hardcoded.sh
# 1. Lädt RNA-Seq Samples herunter
# 2. Erstellt am Ende automatisch das samplesheet.csv

mkdir -p data/fastq

# Liste aller RNA-Seq IDs
IDS="SRR628588 SRR628589 SRR628585 SRR628584 SRR628583 SRR628586 SRR628587 SRR628582"

echo "=== Starte Download für RNA-Seq Samples ==="

for id in $IDS; do
    echo "---------------------------------------"
    echo "Bearbeite: $id"
    
    URL="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR628/$id"
    
    # Datei 1
    if [ ! -f "data/fastq/${id}_1.fastq.gz" ]; then
        wget -P data/fastq/ "$URL/${id}_1.fastq.gz"
    else
        echo "✓ ${id}_1 schon da."
    fi

    # Datei 2
    if [ ! -f "data/fastq/${id}_2.fastq.gz" ]; then
        wget -P data/fastq/ "$URL/${id}_2.fastq.gz"
    else
        echo "✓ ${id}_2 schon da."
    fi
done

echo "---------------------------------------"
echo "=== Erstelle data/samplesheet.csv ==="

# Hier wird das Samplesheet geschrieben
cat > data/samplesheet.csv << 'EOF'
sample_id,condition,fastq_1,fastq_2
SRR628582,SF3B1_mutant,data/fastq/SRR628582_1.fastq.gz,data/fastq/SRR628582_2.fastq.gz
SRR628583,SF3B1_mutant,data/fastq/SRR628583_1.fastq.gz,data/fastq/SRR628583_2.fastq.gz
SRR628584,SF3B1_mutant,data/fastq/SRR628584_1.fastq.gz,data/fastq/SRR628584_2.fastq.gz
SRR628585,SF3B1_wildtype,data/fastq/SRR628585_1.fastq.gz,data/fastq/SRR628585_2.fastq.gz
SRR628586,SF3B1_wildtype,data/fastq/SRR628586_1.fastq.gz,data/fastq/SRR628586_2.fastq.gz
SRR628587,SF3B1_wildtype,data/fastq/SRR628587_1.fastq.gz,data/fastq/SRR628587_2.fastq.gz
SRR628588,SF3B1_wildtype,data/fastq/SRR628588_1.fastq.gz,data/fastq/SRR628588_2.fastq.gz
SRR628589,SF3B1_wildtype,data/fastq/SRR628589_1.fastq.gz,data/fastq/SRR628589_2.fastq.gz
EOF

echo "✓ Samplesheet erstellt: data/samplesheet.csv"
echo "✓ Alles erledigt!"
