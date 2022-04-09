import urllib.error
import pandas as pd
from Bio import SeqIO, Entrez, Restriction
from Bio.Seq import Seq
from openpyxl import workbook, load_workbook
from pydna.dseqrecord import Dseqrecord
from pydna.design import primer_design
import random
from tqdm import tqdm
import sqlite3
from PyQt5 import QtCore, QtGui, QtWidgets, Qt, uic
from PyQt5.QtCore import pyqtSignal, Qt
import sys

# Parameters
verbose = False
codons = ['ATT', 'ATC', 'ATA', 'CTT', 'CTC', 'CTA', 'CTG', 'TTA', 'TTG', 'GTT', 'GTC', 'GTA', 'GTG', 'TTT',
          'TTC', 'ATG', 'TGT', 'TGC', 'GCT', 'GCC', 'GCA', 'GCG', 'GGT', 'GGC', 'GGA', 'GGG', 'CCT', 'CCC',
          'CCA', 'CCG', 'ACT', 'ACC', 'ACA', 'ACG', 'TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC', 'TAT', 'TAC',
          'TGG', 'CAA', 'CAG', 'AAT', 'AAC', 'CAT', 'CAC', 'GAA', 'GAG', 'GAT', 'GAC', 'AAA', 'AAG', 'CGT',
          'CGC', 'CGA', 'CGG', 'AGA', 'AGG', 'TAA', 'TAG', 'TGA']

qtCreatorFile = "SeqC.ui"

Ui_MainWindow, QtBaseClass = uic.loadUiType(qtCreatorFile)


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


class MyWindow(QtWidgets.QMainWindow, Ui_MainWindow):

    def __init__(self, parent=None):
        super(MyWindow, self).__init__(parent)

        QtWidgets.QMainWindow.__init__(self)
        Ui_MainWindow.__init__(self)
        self.setupUi(self)
        self.tabWidget.setCurrentIndex(0)
        self.setFixedSize(1161, 781)
        self.setStyleSheet = 'QPushButton{border-radius:40px;}'
        Entrez.email = self.emailSetting.text()

        #  On Click Actions for List View entrez Scraping Page
        self.addtoList.clicked.connect(self.list_results)

        #  On click actions for table view entrez scraping page.
        self.addtoTable.clicked.connect(self.table_results)

    # On click methods for buttons.
    def table_results(self):
        pass

    def list_results(self):
        print(self.entrez_fetch_protein())
        self.entrez_fetch_protein()
        self.idListWidget.addItem(self.idLineEdit.text())
        self.aaListWidget.addItem(str(self.aa_sequence))
        self.rnaListWidget.addItem(self.mrna)
        self.dnaListWidget.addItem(self.cdna)
        self.fwListWidget.addItem(self.fwString)
        self.revListWidget.addItem(self.revString)
        self.assemblyFwListWidget.addItem(self.fwRestricting)
        self.assemblyRevListWidget.addItem(self.revRestricting)


    def entrez_fetch_protein(self):
        try:
            # Get AA sequence from Entrez database.
            proteinID = self.idLineEdit.text()
            handle = Entrez.efetch(db="protein", id=proteinID, rettype="fasta")
            record = SeqIO.read(handle, "fasta")

            aa_sequence = record.seq
            if verbose: print(f"AA Sequence: {aa_sequence}")

            self.aa_sequence = aa_sequence

        except urllib.error.HTTPError:  # If cannot find AA sequence, fill all columns in excel "n/a".
            if verbose: print(f"HTTP Error.")


        except Exception as e:  # Handle unexpected EOF and continue loop.
            print(e)

        try:
            # Reverse translate Protein to RNA
            codon_list = list(map(aa_to_codon, aa_sequence))  # Use RT to get codon list from AAs)
            if verbose: print(f"Codons: {codon_list}")

            mrna = ''.join(codon_list)  # Merge codons to get mRNA strand.
            if verbose: print("RNA: " + mrna)

            self.mrna = mrna

            cdna = mrna.replace("U", "T")  # Reverse translate mRNA to Crick (anti-coding) strand.
            cdna_seq = Seq(cdna)  # Convert cDNA sequence to Seq data structure.
            if verbose: print("cDNA (5-3): " + cdna_seq)
            if verbose: print("cDNA (3-5): " + cdna_seq.complement())

            self.cdna = cdna

            # Design primers.
            cdna_dseq = Dseqrecord(cdna_seq, str(cdna_seq.complement), linear=True, circular=False)
            ampl = primer_design(cdna_dseq)

            fw_seq = ampl.forward_primer
            rev_seq = ampl.reverse_primer[::-1]

            self.fwString = str(fw_seq)
            self.revString = str(rev_seq)

            # Check for hairpin structure in primer.

            for codon in codon_list:
                codon = codon
                anti_codon = str(Seq(codon).complement())

                # check if codon in ForwardPrimer and ReversePrimer
                if self.fwString.find(codon) == -1 and self.revString.find(anti_codon) == -1:
                    break

            self.fwRestricting = codon + str(self.fwString) + str(getattr(Restriction, self.restrictionLine_1.text()))
            self.revRestricting = anti_codon + self.revString + str(getattr(Restriction, self.restrictionLine_2.text()))

        except Exception as e:
            print(e)



if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)
    window = MyWindow()
    window.show()
    sys.exit(app.exec_())
