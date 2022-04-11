import urllib.error
from Bio import SeqIO, Entrez, Restriction
from Bio.Seq import Seq
from PyQt5.QtWidgets import QTableWidgetItem
from pydna.dseqrecord import Dseqrecord
from pydna.design import primer_design
import random
from PyQt5 import QtWidgets, uic
import sys


def entrez_fetch_nucleotide(self, nucleotideID, restriction1, restriction2):  # TODO connect to ui logic
    try:
        # Get AA sequence from Entrez database.
        handle = Entrez.efetch(db="nucleotide", id=nucleotideID, rettype="fasta")
        record = SeqIO.read(handle, "fasta")

        dna_sequence = record.seq

        # DNA to RNA

        mrna = Seq.transcribe(dna_sequence)

        # mRNA to AA_Sequence

        mrna.translate()

        # Design primers.

        cdna_dseq = Dseqrecord(cdna_seq, str(cdna_seq.complement), linear=True, circular=False)
        ampl = primer_design(cdna_dseq)

        fw_seq = ampl.forward_primer
        rev_seq = ampl.reverse_primer[::-1]

        for codon in codon_list:  # Check for hairpin structure in primer.
            codon = codon
            anti_codon = str(Seq(codon).complement())

            # check if codon in ForwardPrimer and ReversePrimer
            if fw_seq.seq.find(codon) == -1 and rev_seq.seq.find(anti_codon) == -1:
                break

        # Assembly primers with the adition of extra codon and restriction sites selected
        fw_assembly = codon.lower() + fw_seq.seq + getattr(getattr(Restriction,
                                                                   restriction1),
                                                           'site').lower()

        rev_assembly = anti_codon.lower() + rev_seq.seq + getattr(getattr(Restriction,
                                                                          restriction2),
                                                                  'site').lower()

        # Create insert object with the properties to return and display.
        insert = Insert(str(nucleotideID), str(aa_sequence), str(dna_sequence), str(mrna), str(fw_seq.seq),
                        str(rev_seq.seq), str(fw_assembly), str(rev_assembly))

        # Set status for notifying user of method success.

        return insert

    except urllib.error.HTTPError:  # If cannot find AA sequence, fill all columns "n/a".
        if verbose:
            print("HTTP Error.")

    except Exception as e:
        print(e)
