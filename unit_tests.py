import unittest
import gibson_asembly as gibby
from Bio import Entrez, SeqIO
import SequenceCompilerLibrary as SeqLib


class RefSeqScrapingTests(unittest.TestCase):
    def test_protein_fetch(self):
        protein_id = 'CAA50398.1'
        restrictions = 'XbaI'
        resultInsert = SeqLib.entrez_fetch_protein(protein_id, 'XbaI', 'XbaI')
        self.assertTrue(resultInsert is not None)

    def test_nucleotide_fetch(self):
        nucleotide_id = 'X71080.1'
        restrictions = 'XbaI'
        resultInsert = SeqLib.entrez_fetch_nucleotide(nucleotide_id, 'XbaI', 'XbaI')
        self.assertTrue(resultInsert is not None)


class AssemblyTests(unittest.TestCase):
    def test_single_gibson_nucleotide(self):
        plasmid_id = 'L37382.1'
        plasmid_digestion = 'XbaI'
        insert_ids = [('X71080.1', 'XbaI', 'nucleotide')]
        vector = gibby.assembly_analysis(plasmid_id=plasmid_id,
                                         plasmid_digestion=plasmid_digestion,
                                         insert_ids=insert_ids)
        self.assertTrue(vector.circular, 'Insert did not circularize!')  # Check if circular assembly happens.
        # Expected vector size after looping can also be asserted for higher level of precision.

    def test_single_gibson_protein(self):
        plasmid_id = 'L37382.1'
        plasmid_digestion = 'HindIII'
        insert_ids = [('CAA50398.1', 'XbaI', 'protein')]
        vector = gibby.assembly_analysis(plasmid_id=plasmid_id,
                                         plasmid_digestion=plasmid_digestion,
                                         insert_ids=insert_ids)
        self.assertTrue(vector.circular, 'Insert did not circularize!')  # Check if circular assembly happens.
        # Expected vector size after looping can also be asserted for higher level of precision.

    def test_multiple_gibson(self):
        plasmid_id = 'L37382.1'
        plasmid_digestion = 'XbaI'
        insert_ids = [('X71080.1', 'XbaI', 'nucleotide'),
                      ('X71080.1', 'XbaI', 'nucleotide'),
                      ('CAA50398.1', 'HindIII', 'protein')]

        vector = gibby.assembly_analysis(plasmid_id=plasmid_id,
                                         plasmid_digestion=plasmid_digestion,
                                         insert_ids=insert_ids)
        self.assertTrue(vector.circular, 'Insert did not circularize!')
        # Expected vector size after looping can also be asserted for higher level of precision.


class PcrSimulationTests:
    pass


class UITests:
    pass


class PlasmidAnalysisTests:
    pass


class PrimerDesignTests:
    pass


if __name__ == '__main__':
    AssemblyTests
    RefSeqScrapingTests
