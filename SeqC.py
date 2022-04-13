import urllib.error
from Bio import SeqIO, Entrez, Restriction
from Bio.Seq import Seq
from PyQt5.QtWidgets import QTableWidgetItem
from pydna.dseqrecord import Dseqrecord
from pydna.design import primer_design
import random
from PyQt5 import QtWidgets, uic
import sys
import SequenceCompiler as SeqC

# Parameters
qtCreatorFile = "SeqC.ui"

Ui_MainWindow, QtBaseClass = uic.loadUiType(qtCreatorFile)


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

    # On click method for inhabitating the table with protein or nucleotide data.
    def table_results(self):
        last_row = self.entrezTableWidget.rowCount()
        row_toAdd = last_row + 1

        # Increase row number to inhabit incoming data in the following lines.
        self.entrezTableWidget.setRowCount(row_toAdd)

        # Get protein data using function in SequenceCompiler library
        if self.proteinRadio.isChecked():
            resultInsert = SeqC.entrez_fetch_protein(self.idLineEdit.text(), self.restrictionLine1.text(),
                                                     self.restrictionLine2.text())

        # Get nucleotide data using function in SequenceCompiler library
        elif self.nucleotideRadio.isChecked():
            resultInsert = SeqC.entrez_fetch_nucleotide(self.idLineEdit.text(), self.restrictionLine1.text(),
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

    # On click method for running PCR simulation.
    def pcr_results(self):
        pass

    def assembly_results(self):
        pass

    def plasmid_results(self):
        pass


if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)
    window = MyWindow()
    window.show()
    sys.exit(app.exec_())
