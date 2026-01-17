#!/bin/bash
# build_hisat_index.sh

mkdir -p reference/hisat2_index

echo "Baue HISAT2 Index ..."

# hisat2-build [fasta] [output_prefix]
hisat2-build -p 4 \
    reference/genome/GRCh38.primary_assembly.genome.fa \
    reference/hisat2_index/genome

echo "✓ Fertig!"