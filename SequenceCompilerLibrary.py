import urllib.error
from Bio import SeqIO, Entrez, Restriction, SeqFeature
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


def prepare_insert(insert_id, insert_type, insert_restriction_5, insert_restriction_3):
    if insert_type == 'nucleotide':
        newInsert = entrez_fetch_nucleotide(insert_id, insert_restriction_5, insert_restriction_3)

    elif insert_type == 'protein':
        newInsert = entrez_fetch_protein(insert_id, insert_restriction_5, insert_restriction_3)

    else:
        raise IndexError('This function can only accept nucleotide or protein type inserts.')

    # Run PCR with the assembly primers to introduce restriction sites we need.
    pcr_prod = pcr(newInsert.fw_assembly, newInsert.rev_assembly, newInsert.dna)

    # Convert pcr product sequence to Dseqrecord object.
    pcr_prod_dseqrecord = Dseqrecord(pcr_prod.seq, linear=True, circular=False)

    if insert_restriction_3 == insert_restriction_5:
        # Cut 5' and 3' ends.
        restricted_pcr_prod_dseqrecords = pcr_prod_dseqrecord.cut(getattr(Restriction, insert_restriction_5))
        restricted_pcr_prod_dseqrecord = restricted_pcr_prod_dseqrecords[1]

    else:
        # Cut 5' end.
        five_restricted_inserts = pcr_prod_dseqrecord.cut(getattr(Restriction, insert_restriction_5))
        five_restricted_insert = five_restricted_inserts[1]

        # Cut 3' end.
        all_restricted_inserts = five_restricted_insert.cut(getattr(Restriction, insert_restriction_3))
        restricted_pcr_prod_dseqrecord = all_restricted_inserts[1]

    return restricted_pcr_prod_dseqrecord


def gibson_assembly(plasmid_id='L37382.1',
                    plasmid_digestion='HindIII',
                    insert_ids=(('X71080.1', 'HindIII', 'nucleotide'),
                                ('X71080.1', 'XbaI', 'protein'))):

    # Fetch and Digest plasmid
    with Entrez.efetch(db='nucleotide', rettype='gb', retmode='text', id=plasmid_id) as handle:
        seq_record = SeqIO.read(handle, 'gb')

    dseq_record = Dseqrecord(seq_record, linear=False, circular=True)
    restricted_plasmid = dseq_record.linearize(getattr(Restriction, plasmid_digestion))

    if len(insert_ids) == 1:
        single_insert = True
    else:
        single_insert = False

    # if there is only one insert the logic chooses to insert digestion sites to the insert homogeneously.
    if single_insert:
        insert_id = insert_ids[0][0]
        insert_type = insert_ids[0][2]
        insert_restriction_5 = plasmid_digestion
        insert_restriction_3 = plasmid_digestion

        restricted_pcr_prod_dseqrecord = prepare_insert(insert_id,
                                                        insert_type,
                                                        insert_restriction_5,
                                                        insert_restriction_3)

        assemblyResult = (restricted_plasmid + restricted_pcr_prod_dseqrecord).looped()

    # If there is only one insert the logic chooses to insert digestion sites to the insert heterogeneously with next
    # restriction site as the first of the next one.
    else:
        ampliconsList = []
        restricted_ampliconsList = []
        linear_assembly = restricted_plasmid
        i = 0
        # If in the first insert position, digest 5' with plasmid restricting enzyme, 3' with tuple spesifcied enzyme.
        if i == 0:
            insert_id = insert_ids[0][0]
            insert_type = insert_ids[0][2]
            insert_restriction_5 = plasmid_digestion
            insert_restriction_3 = insert_ids[i + 1][1]

            restricted_pcr_prod_dseqrecord = prepare_insert(insert_id,
                                                            insert_type,
                                                            insert_restriction_5,
                                                            insert_restriction_3)

            # Add results to the respective lists.
            restricted_ampliconsList.append(restricted_pcr_prod_dseqrecord)
            print(restricted_plasmid, restricted_pcr_prod_dseqrecord)
            linear_assembly = restricted_plasmid + restricted_pcr_prod_dseqrecord

        # Elif the last insert position, digest 5' with tuple spesifcied enzyme, 3' with plasmid restricting enzyme.
        elif i == len(insert_ids)-1:
            insert_id = insert_ids[0][0]
            insert_type = insert_ids[0][2]
            insert_restriction_5 = insert_ids[i][1]
            insert_restriction_3 = plasmid_digestion

            restricted_pcr_prod_dseqrecord = prepare_insert(insert_id,
                                                            insert_type,
                                                            insert_restriction_5,
                                                            insert_restriction_3)

            restricted_ampliconsList.append(restricted_pcr_prod_dseqrecord)
            linear_assembly = restricted_plasmid + restricted_pcr_prod_dseqrecord

        # Else in a mid insert position, digest 5' and 3' with tuple spesifcied enzyme.
        else:
            insert_id = insert_ids[0][0]
            insert_type = insert_ids[0][2]
            insert_restriction_5 = insert_ids[i][1]
            insert_restriction_3 = insert_ids[i + 1][1]
            restricted_pcr_prod_dseqrecord = prepare_insert(insert_id,
                                                            insert_type,
                                                            insert_restriction_5,
                                                            insert_restriction_3)

            restricted_ampliconsList.append(restricted_pcr_prod_dseqrecord)
            linear_assembly = restricted_plasmid + restricted_pcr_prod_dseqrecord

        assemblyResult = linear_assembly.looped()

    print('Doing circular assembly.')
    print(restricted_plasmid)

    return assemblyResult


def plasmid_analysis(plasmid_id):
    with Entrez.efetch(db="nucleotide", rettype="gb", retmode="text", id=plasmid_id) as handle:
        plasmid = SeqIO.read(handle, "gb")  # using "gb" as an alias for "genbank"

    features = plasmid.features




if __name__ == "__main__":
    plasmid_analysis('L37382.1')