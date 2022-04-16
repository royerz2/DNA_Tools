from Bio import SeqIO, Entrez, Restriction
import SequenceCompilerLibrary as SeqC

with Entrez.efetch(db="nucleotide", rettype="gb", retmode="text", id='X71080.1') as handle:
    seq_record = SeqIO.read(handle, "gb")  # using "gb" as an alias for "genbank"

restriction1 = 'XbaI'
restriction2 = 'XbaI'



fw_primer, rev_primer, fw_assembly, rev_assembly = SeqC.design_primers(seq_record.seq, restriction1, restriction2)

pcr = SeqC.pcr_simulation('GGAACC'+fw_primer, 'GGAACC'+rev_primer, seq_record.seq)

print(pcr.figure())
