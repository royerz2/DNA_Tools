import urllib.error
from Bio import SeqIO, Entrez, Restriction
from Bio.Seq import Seq
from PyQt5.QtWidgets import QTableWidgetItem
from pydna.dseqrecord import Dseqrecord
from pydna.design import primer_design
import random
from PyQt5 import QtWidgets, uic
import sys
import SequenceCompilerLibrary as SeqLib
import sqlite3

# DataBase
conn = sqlite3.connect('seqc.db')

c = conn.cursor()

c.execute('''CREATE TABLE IF NOT EXISTS "Inserts" (
            "id" TEXT NOT NULL,
            "dna_sequence" TEXT NOT NULL,
            "rna_sequence"	TEXT NOT NULL,
            "amino_acid" TEXT NOT NULL,
            "forward_primer" TEXT NOT NULL,
            "reverse_primer" TEXT NOT NULL,
            "assembly_fw" TEXT NOT NULL,
            "assembly_rev" TEXT NOT NULL,
            "type" TEXT NOT NULL
            )''')

c.execute('''CREATE TABLE IF NOT EXISTS "Plasmids" (
            "id" TEXT NOT NULL,
            "dna_sequence" TEXT NOT NULL,
            "forward_primer" TEXT NOT NULL,
            "reverse_primer" TEXT NOT NULL,
            "assembly_fw" TEXT NOT NULL,
            "assembly_rev" TEXT NOT NULL,
            "type" TEXT NOT NULL
            
            )''')

c.execute('''CREATE TABLE IF NOT EXISTS "Assemblies" (
            "id" TEXT NOT NULL,
            "part_ids" TEXT NOT NULL,
            "type" TEXT NOT NULL,
            "digestion" TEXT,
            "fw_primer" TEXT,
            "rev_primer" TEXT
            )''')

# Get UI file.
qtCreatorFile = "SeqC.ui"

# Create Window and QtBase instances.
Ui_MainWindow, QtBaseClass = uic.loadUiType(qtCreatorFile)


# The Class that encapsulates the UI elements, inheriting from QtMainWindow.
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

        # On click actions for PCR simulation tab.

        self.getPcr.clicked.connect(self.pcr_results)

    # On click methods for buttons.

    # On click method for inhabiting the table with protein or nucleotide data.
    def table_results(self):
        last_row = self.entrezTableWidget.rowCount()  # Get the last row index of the table
        row_toAdd = last_row + 1  # Add 1 to last row for indexing.

        # Increase row number to inhabit incoming data in the following lines.
        self.entrezTableWidget.setRowCount(row_toAdd)

        # Get protein data using function in SequenceCompiler library
        if self.proteinRadio.isChecked():
            resultInsert = SeqLib.entrez_fetch_protein(self.idLineEdit.text(), self.restrictionLine1.text(),
                                                       self.restrictionLine2.text())

        # Get nucleotide data using function in SequenceCompiler library
        elif self.nucleotideRadio.isChecked():
            resultInsert = SeqLib.entrez_fetch_nucleotide(self.idLineEdit.text(), self.restrictionLine1.text(),
                                                          self.restrictionLine2.text())
        else:
            raise Exception

        # Add values to the next row of the table
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

    # On click method for running PCR simulation.
    def pcr_results(self):
        if self.automaticPcrRadio.isChecked():
            print('Getting values from Entrez database.')
            with Entrez.efetch(db="nucleotide", rettype="gb", retmode="text", id=self.nucleotidePcr.text()) as handle:
                seq_record = SeqIO.read(handle, "gb")  # using "gb" as an alias for "genbank"ss
            dna_Seq = Seq(seq_record.seq)
            dna_Dseq = Dseqrecord(dna_Seq, str(dna_Seq.complement), linear=True, circular=False)

            ampl = primer_design(dna_Dseq)

            dna_sequence = seq_record.seq
            fw_primer = ampl.forward_primer
            rev_primer = ampl.reverse_primer

            pcr_prod = SeqLib.pcr_simulation(fw_primer, rev_primer, dna_sequence)
            print(pcr_prod.figure())

            pcr_success = True

            # Print results to the pcrResult QTextEdit on front end.
            self.pcrResult.setPlainText(str(pcr_prod.figure()) +
                                        '\n\n\n' +
                                        str(pcr_prod.program()))

        else:
            fw_primer = self.fwPcr.text()
            rev_primer = self.revPcr.text()
            dna_sequence = self.dnaPcr.toPlainText()

        pcr_prod = SeqLib.pcr_simulation(fw_primer, rev_primer, dna_sequence)
        print(pcr_prod.figure())

        pcr_success = True

        # Print results to the pcrResult QTextEdit on front end.
        self.pcrResult.setPlainText(str(pcr_prod.figure()) +
                                    '\n\n\n' +
                                    str(pcr_prod.program()))


    # On click method for running assembly simulation.
    def assembly_results(self):
        pass

    # On click method for running plasmid analysis.
    def plasmid_results(self):
        pass


if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)
    window = MyWindow()
    window.show()
    sys.exit(app.exec_())
