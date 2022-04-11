import urllib.error
from Bio import SeqIO, Entrez, Restriction
from Bio.Seq import Seq
from PyQt5.QtWidgets import QTableWidgetItem
from pydna.dseqrecord import Dseqrecord
from pydna.design import primer_design
import random
from PyQt5 import QtWidgets, uic
import sys

# Parameters
verbose = False  # Set to true to see output as they are imported from
codons = ['ATT', 'ATC', 'ATA', 'CTT', 'CTC', 'CTA', 'CTG', 'TTA', 'TTG', 'GTT', 'GTC', 'GTA', 'GTG', 'TTT', 'TTC',
          'ATG', 'TGT', 'TGC', 'GCT', 'GCC', 'GCA', 'GCG', 'GGT', 'GGC', 'GGA', 'GGG', 'CCT', 'CCC', 'CCA', 'CCG',
          'ACT', 'ACC', 'ACA', 'ACG', 'TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC', 'TAT', 'TAC', 'TGG', 'CAA', 'CAG',
          'AAT', 'AAC', 'CAT', 'CAC', 'GAA', 'GAG', 'GAT', 'GAC', 'AAA', 'AAG', 'CGT', 'CGC', 'CGA', 'CGG', 'AGA',
          'AGG', 'TAA', 'TAG', 'TGA']

qtCreatorFile = "SeqC.ui"

Ui_MainWindow, QtBaseClass = uic.loadUiType(qtCreatorFile)


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

        #  On click actions and settings for table view entrez scraping page.
        self.entrezTableWidget.setRowCount(1)  # Row count
        self.entrezTableWidget.setColumnCount(8)  # Column count

        self.addtoTable.clicked.connect(self.table_results)
        self.entrezTableWidget.setItem(0, 0, QTableWidgetItem("ID"))
        self.entrezTableWidget.setItem(0, 1, QTableWidgetItem("DNA"))
        self.entrezTableWidget.setItem(0, 2, QTableWidgetItem("RNA"))
        self.entrezTableWidget.setItem(0, 3, QTableWidgetItem("Protein"))
        self.entrezTableWidget.setItem(0, 4, QTableWidgetItem("Fw Primer"))
        self.entrezTableWidget.setItem(0, 5, QTableWidgetItem("Rev Primer"))
        self.entrezTableWidget.setItem(0, 6, QTableWidgetItem("Fw Assembly"))
        self.entrezTableWidget.setItem(0, 7, QTableWidgetItem("Rev Assembly"))

    # On click methods for buttons.
    def table_results(self):
        last_row = self.entrezTableWidget.rowCount()
        row_toAdd = last_row + 1

        # Increase row number to inhabit incoming data in the following lines.
        self.entrezTableWidget.setRowCount(row_toAdd)

        if self.proteinRadio.isChecked():
            resultInsert = self.entrez_fetch_protein(self.idLineEdit.text(), self.restrictionLine1.text(),
                                                     self.restrictionLine2.text())
        elif self.nucleotideRadio.isChecked():
            resultInsert = self.entrez_fetch_nucleotide(self.idLineEdit.text(), self.restrictionLine1.text(),
                                                        self.restrictionLine2.text())
        else:
            raise Exception



        self.entrezTableWidget.setItem(last_row, 1, QTableWidgetItem(resultInsert.dna))
        self.entrezTableWidget.setItem(last_row, 2, QTableWidgetItem(resultInsert.mrna))
        self.entrezTableWidget.setItem(last_row, 3, QTableWidgetItem(resultInsert.aa_sequence))
        self.entrezTableWidget.setItem(last_row, 4, QTableWidgetItem(resultInsert.fw_primer))
        self.entrezTableWidget.setItem(last_row, 5, QTableWidgetItem(resultInsert.rev_primer))
        self.entrezTableWidget.setItem(last_row, 6, QTableWidgetItem(resultInsert.fw_assembly))
        self.entrezTableWidget.setItem(last_row, 7, QTableWidgetItem(resultInsert.rev_assembly))
        # Added last so that line index continuity is kept when exception occurs
        self.entrezTableWidget.setItem(last_row, 0, QTableWidgetItem(self.idLineEdit.text()))

        resultInsert = None

    def entrez_fetch_protein(self, proteinID, restriction1, restriction2):
        try:
            # Get AA sequence from Entrez database.
            handle = Entrez.efetch(db="protein", id=proteinID, rettype="fasta")
            record = SeqIO.read(handle, "fasta")

            aa_sequence = record.seq
            if verbose:
                print(f"AA Sequence: {aa_sequence}")

            # Reverse translate Protein to RNA
            codon_list = list(map(aa_to_codon, aa_sequence))  # Use RT to get codon list from AAs

            mrna = ''.join(codon_list)  # Merge codons to get mRNA strand.

            if verbose:
                print("RNA: " + mrna)
                print(f"Codons: {codon_list}")

            mrna = mrna

            cdna = mrna.replace("U", "T")  # Reverse translate mRNA to Crick (anti-coding) strand.
            cdna_seq = Seq(cdna)  # Convert cDNA sequence to Seq data structure.

            if verbose:
                print("cDNA (5-3): " + cdna_seq)
                print("cDNA (3-5): " + cdna_seq.complement())

            fw_primer, rev_primer, fw_assembly, rev_assembly= design_primers(cdna, restriction1, restriction2)

            # Create insert object with the properties to return and display.
            insert = Insert(str(proteinID), str(aa_sequence), str(cdna), str(mrna), str(fw_primer.seq),
                            str(rev_primer.seq), str(fw_assembly), str(rev_assembly))

            # Set status for notifying user of method success.
            print(insert.dna)
            return insert

        except urllib.error.HTTPError:  # If cannot find AA sequence, fill all columns "n/a".
            if verbose:
                print("HTTP Error.")

        except Exception as e:
            print(e)

    def entrez_fetch_nucleotide(self, nucleotideID, restriction1, restriction2):
        try:
            # Get AA sequence from Entrez database.
            handle = Entrez.efetch(db="nucleotide", id=nucleotideID, rettype="fasta")
            record = SeqIO.read(handle, "fasta")

            # Get DNA sequence as a string
            dna_sequence = record.seq

            # DNA to RNA
            mrna = Seq.transcribe(dna_sequence)

            # mRNA to AA_Sequence
            mrna.translate()

            # Get primers for the dna and its assembly
            fw_primer, rev_primer, fw_assembly, rev_assembly = design_primers(dna_sequence, restriction1, restriction2)

            # Create insert object with the properties to return and display.
            insert = Insert(str(nucleotideID), 'Not Translated', str(dna_sequence), str(mrna), str(fw_primer.seq),
                            str(rev_primer.seq), str(fw_assembly), str(rev_assembly))

            # Set status for notifying user of method success.

            return insert

        except urllib.error.HTTPError:  # If cannot find AA sequence, fill all columns "n/a".
            if verbose:
                print("HTTP Error.")

        except Exception as e:
            print(e)


def design_primers(dna_sequence, restriction1, restriction2):
    try:
        dna_Seq = Seq(dna_sequence)
        dna_Dseq = Dseqrecord(dna_Seq, str(dna_Seq.complement), linear=True, circular=False)
        ampl = primer_design(dna_Dseq)

        fw_primer = ampl.forward_primer
        rev_primer = ampl.reverse_primer[::-1]

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

    except Exception as e:
        print(e)


if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)
    window = MyWindow()
    window.show()
    sys.exit(app.exec_())
