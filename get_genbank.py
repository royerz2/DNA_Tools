from Bio import Entrez
from Bio import SeqIO

nucleotide_id = "M77789.2"
Entrez.email = "A.N.Other@example.com"
with Entrez.efetch(db="nucleotide", rettype="gb", retmode="text", id=nucleotide_id) as handle:
    seq_record = SeqIO.read(handle, "gb")  # using "gb" as an alias for "genbank"

print("%s with %i features" % (seq_record.id, len(seq_record.features)))

