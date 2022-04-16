import SequenceCompilerLibrary as SeqLib
from Bio import Entrez, SeqIO, Restriction
from Bio.Seq import Seq
from pydna.amplify import pcr
from pydna.dseqrecord import Dseqrecord
from pydna.design import primer_design

codons = ['ATT', 'ATC', 'ATA', 'CTT', 'CTC', 'CTA', 'CTG', 'TTA', 'TTG', 'GTT', 'GTC', 'GTA', 'GTG', 'TTT', 'TTC',
          'ATG', 'TGT', 'TGC', 'GCT', 'GCC', 'GCA', 'GCG', 'GGT', 'GGC', 'GGA', 'GGG', 'CCT', 'CCC', 'CCA', 'CCG',
          'ACT', 'ACC', 'ACA', 'ACG', 'TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC', 'TAT', 'TAC', 'TGG', 'CAA', 'CAG',
          'AAT', 'AAC', 'CAT', 'CAC', 'GAA', 'GAG', 'GAT', 'GAC', 'AAA', 'AAG', 'CGT', 'CGC', 'CGA', 'CGG', 'AGA',
          'AGG', 'TAA', 'TAG', 'TGA']

with Entrez.efetch(db="nucleotide", rettype="gb", retmode="text", id='X71080.1') as handle:
    seq_record = SeqIO.read(handle, "gb")  # using "gb" as an alias for "genbank"


def design_primers(dna_sequence, restriction1, restriction2):
    dna_Seq = Seq(dna_sequence)
    dna_Dseq = Dseqrecord(dna_Seq, str(dna_Seq.complement), linear=True, circular=False)
    ampl = primer_design(dna_Dseq)

    fw_primer = ampl.forward_primer
    rev_primer = ampl.reverse_primer

    for codon in codons:  # Check for hairpin structure in primer.
        anti_codon = str(Seq(codon).complement())

        # check if codon in ForwardPrimer and ReversePrimer
        if fw_primer.seq.find(codon) == -1 and rev_primer.seq.find(anti_codon) == -1:
            break

    # Assembly primers with the adition of extra codon and restriction sites selected
    fw_assembly = codon.lower() + fw_primer.seq + getattr(getattr(Restriction,
                                                                  restriction1),
                                                          'site').lower()

    rev_assembly = anti_codon.lower() + rev_primer.seq + getattr(getattr(Restriction,
                                                                         restriction2),
                                                                 'site').lower()

    # Create insert object with the properties to return and display.
    return fw_primer, rev_primer, fw_assembly, rev_assembly


print(design_primers(str(seq_record.seq), 'XbaI', 'XbaI'))
