import urllib.error
import pandas as pd
from Bio import SeqIO, Entrez
from Bio.Seq import Seq
from openpyxl import workbook, load_workbook
from pydna.dseq import Dseq
from pydna.dseqrecord import Dseqrecord
from pydna.design import primer_design
import random
from tqdm import tqdm

# Parameters
ProteinsSequenced = True  # False activates Entrez sequence fetching
SeqRevTranscribed = False  # False activates reverse transcription
filename = "VBPO.xlsx"
Entrez.email = "royerz.2002@gmail.com"
verbose = False


def aa_to_codon(n):
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


def id_to_aa_sequence(protein_id):
    try:
        # Get AA sequence from Entrez database.
        handle = Entrez.efetch(db="protein", id=protein_id, rettype="fasta")
        record = SeqIO.read(handle, "fasta")

        aa_sequence = record.seq
        if verbose: print(f"AA Sequence: {aa_sequence}")
        ws[f'D{i + 2}'] = f"{aa_sequence}"

        return aa_sequence

    except urllib.error.HTTPError:  # If cannot find AA sequence, fill all columns in excel "n/a".
        if verbose: print(f"HTTP Error at position {i}.")
        ws[f'C{i + 2}'] = "n/a"
        ws[f'D{i + 2}'] = "n/a"
        ws[f'E{i + 2}'] = "n/a"
        ws[f'F{i + 2}'] = "n/a"
        ws[f'G{i + 2}'] = "n/a"
        ws[f'H{i + 2}'] = "n/a"

    except Exception as e:  # Handle unexpected EOF and continue loop.
        print(e)


def aa_to_mrna(aa_sequence):
    try:
        # Reverse translate Protein to RNA
        codon_list = list(map(aa_to_codon, aa_sequence))  # Use RT to get codon list from AAs)
        if verbose: print(f"Codons: {codon_list}")

        mrna = ''.join(codon_list)  # Merge codons to get mRNA strand.
        if verbose: print("RNA: " + mrna)

        return mrna

    except Exception as e:
        print(e)


def mrna_to_cdna_seq(mrna):
    try:
        cdna = mrna.replace("U", "T")  # Reverse translate mRNA to Crick (anti-coding) strand.
        cdna_seq = Seq(cdna)  # Convert cDNA sequence to Seq data structure.
        if verbose: print("cDNA (5-3): " + cdna_seq)
        if verbose: print("cDNA (3-5): " + cdna_seq.complement())
        ws[f'E{i + 2}'] = f"{cdna_seq.complement()}"  # Send cDNA Watson (coding) strand to excel sheet.
        ws[f'F{i + 2}'] = f"{cdna_seq}"  # Send cDNA Crick (anti-coding) strand to excel sheet.
        return cdna_seq

    except Exception as e:
        print(e)


def cdna_to_primer(cdna_seq):
    try:
        # Design primers.
        cdna_dseq = Dseqrecord(cdna_seq, str(cdna_seq.complement)[::-1], linear=True, circular=False)
        ampl = primer_design(cdna_dseq)
        fp = ampl.forward_primer
        rp = ampl.reverse_primer
        ws[f'G{i + 2}'] = f"{str(fp.seq)}"  # Send forward primer to excel sheet.
        ws[f'H{i + 2}'] = f"{str(rp.seq)}"  # Send reverse primer to excel sheet.

        primers = (fp, rp)

        return primers

    except Exception as e:
        print(e)


def gibson_assembly(plasmid, *inserts):
    print(f'Plasmid: {plasmid}')
    i = 0
    for insert in inserts:

        print(f'Insert {i}: {insert}.')
        i += 1



wb = load_workbook(filename)
sheets = wb.sheetnames
ws = wb[sheets[0]]

dataframe = pd.DataFrame(ws.values, columns=["prot", "spp", "code", "aa", "cDNA(5-3)",
                                             "cDNA(3-5)", "primer", "rev_primer"])
dataframe.drop(index=dataframe.index[0], axis=0, inplace=True)

print("Starting protein fetching, reverse translation and primer design.")

code_list = dataframe.code.tolist()

if verbose: print(code_list)
if verbose: print(len(code_list))

for i in tqdm(range(len(code_list))):
    aa_sequence = id_to_aa_sequence(code_list[i])
    mrna = aa_to_mrna(aa_sequence)
    cdna_seq = mrna_to_cdna_seq(mrna)
    primers = cdna_to_primer(cdna_seq)

wb.save(filename=filename)
