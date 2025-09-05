#!/usr/bin/env python3
import sys
import pandas as pd
from Bio import SeqIO

if len(sys.argv) < 6:
    print("Usage: parse_blast.py <query_fasta> <blast_results.tsv> <positive_windows.fasta> <positive_hits.csv> <db_fasta>")
    sys.exit(1)

query_fasta = sys.argv[1]
blast_file = sys.argv[2]
output_fasta = sys.argv[3]
output_csv = sys.argv[4]
db_fasta = sys.argv[5]

print("Parsing BLAST results...")

# Load BLAST tabular results
cols = ["qseqid","sseqid","pident","length","mismatch","gapopen",
        "qstart","qend","sstart","send","evalue","bitscore"]
df = pd.read_csv(blast_file, sep="\t", names=cols)

# Keep only alignments with identity >=35% AND alignment length == 80
df = df[(df["pident"] >= 35) & (df["length"] == 80)]

if df.empty:
    print("No positive hits found.")
    sys.exit(0)

# Map subject IDs to descriptions from db_fasta
desc_map = {}
for record in SeqIO.parse(db_fasta, "fasta"):
    desc_map[record.id] = record.description

df["sdescription"] = df["sseqid"].map(desc_map)

# Save all positive hits
df.to_csv(output_csv, index=False)
print(f"Positive hits saved to {output_csv}")

# Deduplicate by sdescription, keeping row with highest pident
cleaned = df.loc[df.groupby("sdescription")["pident"].idxmax()]
cleaned = cleaned.sort_values(by="pident", ascending=False)

clean_csv = output_csv.replace(".csv", "_clean.csv")
cleaned.to_csv(clean_csv, index=False)
print(f"Cleaned hits saved to {clean_csv}")

# Extract corresponding query sequences (positive windows)
positive_ids = set(df["qseqid"])
with open(output_fasta, "w") as out_f:
    for record in SeqIO.parse(query_fasta, "fasta"):
        if record.id in positive_ids:
            SeqIO.write(record, out_f, "fasta")

print(f"Positive window sequences written to {output_fasta}")
