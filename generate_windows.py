import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

query_fasta = sys.argv[1]
output_fasta = sys.argv[2]
window_size = 80

records_out = []

for record in SeqIO.parse(query_fasta, "fasta"):
    seq_len = len(record.seq)
    for i in range(seq_len - window_size + 1):
        window_seq = record.seq[i:i+window_size]
        window_id = f"{record.id}_window{i+1}"
        records_out.append(SeqRecord(window_seq, id=window_id, description=""))

from Bio import SeqIO
SeqIO.write(records_out, output_fasta, "fasta")
print(f"Generated {len(records_out)} sliding windows of size {window_size}.")
