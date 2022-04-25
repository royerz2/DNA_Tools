import urllib.error
from Bio import SeqIO, Entrez, Restriction
from Bio.Seq import Seq
from PyQt5.QtWidgets import QTableWidgetItem
from pydna.assembly import Assembly
from pydna.dseqrecord import Dseqrecord
from pydna.design import primer_design
import random
from PyQt5 import QtWidgets, uic
import sys
from pydna.amplify import pcr

verbose = False  # Set to true to see output as they are imported from
codons = ['ATT', 'ATC', 'ATA', 'CTT', 'CTC', 'CTA', 'CTG', 'TTA', 'TTG', 'GTT', 'GTC', 'GTA', 'GTG', 'TTT', 'TTC',
          'ATG', 'TGT', 'TGC', 'GCT', 'GCC', 'GCA', 'GCG', 'GGT', 'GGC', 'GGA', 'GGG', 'CCT', 'CCC', 'CCA', 'CCG',
          'ACT', 'ACC', 'ACA', 'ACG', 'TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC', 'TAT', 'TAC', 'TGG', 'CAA', 'CAG',
          'AAT', 'AAC', 'CAT', 'CAC', 'GAA', 'GAG', 'GAT', 'GAC', 'AAA', 'AAG', 'CGT', 'CGC', 'CGA', 'CGG', 'AGA',
          'AGG', 'TAA', 'TAG', 'TGA']


class Insert:
    def __init__(self, entrez_id, aa_sequence, dna, mrna, fw_primer, rev_primer, fw_assembly, rev_assembly):
        self.entrez_id = entrez_id
        self.aa_sequence = aa_sequence
        self.dna = dna
        self.mrna = mrna
        self.fw_primer = fw_primer
        self.rev_primer = rev_primer
        self.fw_assembly = fw_assembly
        self.rev_assembly = rev_assembly


def aa_to_codon(n):
    # Switch that converts amino acids to respective codons with codon optimization weights.
    switch = {
        "A": (random.choices(["GCU", "GCC", "GCA", "GCG"], weights=(0.19, 0.25, 0.22, 0.34), k=1)),
        "I": (random.choices(["AUU", "AUC", "AUA"], weights=(0.47, 0.46, 0.07), k=1)),
        "L": (random.choices(["UUA", "UUG", "CUU", "CUC", "CUA", "CUG"],
                             weights=(0.11, 0.11, 0.10, 0.10, 0.03, 0.55), k=1)),
        "M": ["AUG"],
        "V": (random.choices(["GUU", "GUC", "GUA", "GUG"], weights=(0.29, 0.20, 0.17, 0.34), k=1)),
        "F": (random.choices(["UUU", "UUC"], weights=(0.51, 0.49), k=1)),
        "W": ["UUG"],
        "Y": (random.choices(["UAU", "UAC"], weights=(0.53, 0.47), k=1)),
        "N": (random.choices(["AAU", "AAC"], weights=(0.39, 0.61), k=1)),
        "C": (random.choices(["UGU", "UGC"], weights=(0.43, 0.57), k=1)),
        "Q": (random.choices(["CAA", "CAG"], weights=(0.31, 0.69), k=1)),
        "S": (random.choices(["UCU", "UCC", "UCA", "UCG"], weights=(0.19, 0.17, 0.12, 0.13), k=1)),
        "T": (random.choices(["ACU", "ACC", "ACA", "ACG"], weights=(0.16, 0.47, 0.13, 0.24), k=1)),
        "D": (random.choices(["GAU", "GAC"], weights=(0.59, 0.41), k=1)),
        "E": (random.choices(["GAA", "GAG"], weights=(0.7, 0.3), k=1)),
        "R": (random.choices(["CCU", "CGC", "CGA", "CGG"], weights=(0.42, 0.37, 0.05, 0.08), k=1)),
        "H": (random.choices(["AAU", "AAC"], weights=(0.39, 0.61), k=1)),
        "K": (random.choices(["AAA", "AAG"], weights=(0.76, 0.24), k=1)),
        "G": (random.choices(["GGU", "GGC", "GGA", "GGG"], weights=(0.38, 0.40, 0.09, 0.13), k=1)),
        "P": (random.choices(["CCU", "CCC", "CCA", "CCG"], weights=(0.16, 0.10, 0.2, 0.54), k=1)),
    }

    codon_elem = switch[n]
    codon_str = f"{codon_elem[0]}"

    return codon_str


def entrez_fetch_protein(protein_id, restriction1, restriction2):
    try:
        # Get AA sequence from Entrez database.
        with Entrez.efetch(db="protein", id=protein_id, rettype="fasta") as handle:
            record = SeqIO.read(handle, "fasta")

        aa_sequence = record.seq

        # Reverse translate Protein to RNA
        codon_list = list(map(aa_to_codon, aa_sequence))  # Use RT to get codon list from AAs

        mrna = ''.join(codon_list)  # Merge codons to get mRNA strand.

        cdna = mrna.replace("U", "T")  # Reverse translate mRNA to Crick (anti-coding) strand.
        cdna_seq = Seq(cdna)  # Convert cDNA sequence to Seq data structure.

        fw_primer, rev_primer, fw_assembly, rev_assembly = design_primers(cdna, restriction1, restriction2)

        # Create insert object with the properties to return and display.
        insert = Insert(str(protein_id), str(aa_sequence), str(cdna), str(mrna), str(fw_primer.seq),
                        str(rev_primer.seq), str(fw_assembly), str(rev_assembly))

        # Set status for notifying user of method success.
        return insert

    except urllib.error.HTTPError:  # If cannot find AA sequence, fill all columns "n/a".
        print("HTTP Error.")

    except Exception as e:
        print(e)


def entrez_fetch_nucleotide(nucleotide_id, restriction1, restriction2):
    try:
        # Get AA sequence from Entrez database.
        with Entrez.efetch(db="nucleotide", rettype="gb", retmode="text", id=nucleotide_id) as handle:
            record = SeqIO.read(handle, "gb")  # using "gb" as an alias for "genbank"

        # Get DNA sequence as a string
        dna_sequence = record.seq

        # DNA to RNA
        mrna = Seq.transcribe(dna_sequence)

        # mRNA to AA_Sequence
        mrna.translate()

        # Get primers for the dna and its assembly
        fw_primer, rev_primer, fw_assembly, rev_assembly = design_primers(dna_sequence, restriction1, restriction2)

        # Create insert object with the properties to return and display.
        insert = Insert(str(nucleotide_id), 'Not Translated', str(dna_sequence), str(mrna), str(fw_primer.seq),
                        str(rev_primer.seq), str(fw_assembly), str(rev_assembly))

        # Set status for notifying user of method success.

        return insert

    except urllib.error.HTTPError:  # If cannot find AA sequence, fill all columns "n/a".
        if verbose:
            print("HTTP Error.")

    except Exception as e:
        print(e)


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
    fw_assembly = codon.lower() + getattr(getattr(Restriction, restriction1), 'site').lower() + fw_primer.seq

    rev_assembly = anti_codon.lower() + getattr(getattr(Restriction, restriction2), 'site').lower() + rev_primer.seq

    # Create insert object with the properties to return and display.
    return fw_primer, rev_primer, fw_assembly, rev_assembly


def pcr_simulation(forward_primer, reverse_primer, dna):
    pcr_prod = pcr(forward_primer, reverse_primer, dna)

    return pcr_prod


def assembly_analysis(fragments_csv, linear):
    fragments_tuple = tuple(fragments_csv.split(','))
    assembly = Assembly(fragments_tuple, limit=14)

    if linear:
        assembly_result = assembly.assemble_linear()
    else:
        assembly_result = assembly.assemble_circular()

    return assembly, assembly_result


def plasmid_analysis(plasmid_id, enzyme=None, digest=False):
    with Entrez.efetch(db="nucleotide", rettype="gb", retmode="text", id=plasmid_id) as handle:
        plasmid = SeqIO.read(handle, "gb")  # using "gb" as an alias for "genbank"

    # If user wants digestion data, digest the plasmid with the specified
    if plasmid.linear:
        print("This DNA is not a plasmid.")
    elif digest:
        # Check if enzyme exists in database
        try:
            restriction_object = getattr(Restriction, enzyme)
        except Exception as e:
            print("Restriction enzyme does not exist.")

        restriction_site = restriction_object.site

        # Check if plasmid has digestion site:
        getattr(plasmid.seq).has(restriction_site)

        linear_plasmid = plasmid.linearize(enzyme)
        promoter_location = None

        # Here, the restriction sites right after the promoter will be determined and omitted.

        advised_restrictions = None

    return linear_plasmid, promoter_location, advised_restrictions


if __name__ == "__main__":
    print('This file is not made to run.')
