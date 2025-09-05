#!/bin/bash
# Main pipeline to check allergenicity with sliding windows and BLAST

# ---- USER INPUT ----
QUERY_FASTA="subtilisin.fasta"           # your protein sequence
WINDOW_SCRIPT="generate_windows.py"         # Python script to create sliding windows
BLAST_DB_FASTA="COMPARE2025-FastA-Seq-01-27-2025.fasta"
BLAST_DB_NAME="allergen_db"
WINDOWS_FASTA="sliding_windows.fasta"
BLAST_RESULTS="blast_results.tsv"
POSITIVE_WINDOWS="positive_windows.fasta"
POSITIVE_HITS="positive_hits.csv"

# ---- 1. Generate sliding windows using Python ----
echo "Generating sliding windows..."
python3 $WINDOW_SCRIPT $QUERY_FASTA $WINDOWS_FASTA

# ---- 2. Create BLAST database ----
if [ ! -f "${BLAST_DB_NAME}.pin" ]; then
    echo "Creating BLAST database..."
    makeblastdb -in $BLAST_DB_FASTA -dbtype prot -out $BLAST_DB_NAME
else
    echo "BLAST database already exists."
fi

# ---- 3. Run BLASTP ----
echo "Running BLASTP..."
blastp -query $WINDOWS_FASTA -db $BLAST_DB_NAME -out $BLAST_RESULTS \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"

# ---- 4. Filter BLAST results with Python ----
echo "Filtering positive hits..."
python3 parse_blast.py $WINDOWS_FASTA $BLAST_RESULTS $POSITIVE_WINDOWS $POSITIVE_HITS $BLAST_DB_FASTA

echo "Pipeline finished!"
