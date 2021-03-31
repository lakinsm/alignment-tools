#!/usr/bin/env bash

# Wrapper script for reorient_contigs_to_reference.py for nucleotide alignment using BLASTN
# Must have ncbi-blast toolkit installed

query=$1
reference=$2
output=$3

makeblastdb -dbtype nucl -in $reference
blastn -db $reference -query $query -out $output -out