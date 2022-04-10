from Bio import SeqIO, Entrez, Restriction
from Bio.Seq import Seq
from pydna.dseqrecord import Dseqrecord
from pydna.design import primer_design

cdna_seq = Seq('ACACTGATCGATCGATCGATCGATCGATCGTACGTACGTACG')
cdna_dseq = Dseqrecord(cdna_seq, str(cdna_seq.complement), linear=True, circular=False)

ampl = primer_design(cdna_dseq)

fw_seq = ampl.forward_primer
rev_seq = ampl.reverse_primer[::-1]

codon = 'AAA'
anti_codon = 'TTT'

# TODO fix overhang addition data structures
fw_assembly = codon.lower() + fw_seq.seq + getattr(getattr(Restriction, 'XbaI'), 'site').lower()
rev_assembly = anti_codon.lower() + rev_seq.seq + getattr(getattr(Restriction, 'XbaI'), 'site').lower()

print(type(str(fw_assembly)))
