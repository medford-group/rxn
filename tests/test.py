import unittest

import networkx as nx
import networkx.algorithms.isomorphism as iso
from networkx.readwrite import json_graph

import sys
sys.path.insert(0, '../')

from rxn.data_structures import MolGraph, RxnGraph
from rxn.rxn_network import scissions


class MolGraphTestCase(unittest.TestCase):
    """Tests custom methods of MolGraph data structure"""

    def generate_ONNO_test_molecules(self):

        ONNO = MolGraph()
        ONNO.generate('ONNO', data_format='smi')

        ONNO_1 = ONNO.copy()
        ONNO_1.remove_node(sorted(ONNO.nodes())[0])

        ONNO_2 = ONNO.copy()
        ONNO_2.remove_node(sorted(ONNO.nodes())[1])

        ONNO_3 = ONNO.copy()
        ONNO_3.remove_node(sorted(ONNO.nodes())[2])

        ONNO_4 = nx.relabel.convert_node_labels_to_integers(ONNO)

        return ONNO, ONNO_1, ONNO_2, ONNO_3, ONNO_4

    def test_repr_CCO(self):
        # generate ethanol test case using smiles format
        EtOH = MolGraph()
        EtOH.generate('CCO', data_format='smi')
        self.assertEqual(eval(repr(EtOH)), EtOH)

    def test_eq_ONNO_isomers(self):
        ONNO, ONNO_1, ONNO_2, ONNO_3, ONNO_4 = self.generate_ONNO_test_molecules()

        self.assertNotEqual(ONNO_1, ONNO_2)
        self.assertNotEqual(ONNO_2, ONNO_3)

        self.assertEqual(ONNO_3, ONNO_1)  # test isomer equivalency

    def test_eq_ONNO_labels(self):
        ONNO, ONNO_1, ONNO_2, ONNO_3, ONNO_4 = self.generate_ONNO_test_molecules()

        self.assertEqual(ONNO, ONNO_4)  # test node label invariance


class RxnGraphTestCase(unittest.TestCase):
    """Tests custom methods of MolGraph data structure"""

    def generate_ONNO_test_rxns(self):
        ONNO = MolGraph()
        ONNO.generate('ONNO', 'smi')

        H = MolGraph()
        H.generate('[H]', 'smi')

        G1 = ONNO.copy()
        G1.remove_node(sorted(ONNO.nodes())[0])

        G2 = ONNO.copy()
        G2.remove_node(sorted(ONNO.nodes())[1])

        G3 = ONNO.copy()
        G3.remove_node(sorted(ONNO.nodes())[2])

        G4 = nx.relabel.convert_node_labels_to_integers(G3)

        rxn_1 = [[ONNO], [G1, H]]
        rxn_2 = [[ONNO], [G2, H]]
        rxn_3 = [[ONNO], [G3, H]]
        rxn_4 = [[ONNO], [G4, H]]

        R1 = RxnGraph()
        R1.from_rxn_list([rxn_1])

        R2 = RxnGraph()
        R2.from_rxn_list([rxn_2])

        R3 = RxnGraph()
        R3.from_rxn_list([rxn_3])

        R4 = RxnGraph()
        R4.from_rxn_list([rxn_4])

        return R1, R2, R3, R4

    def test_eq_ONNO_isomers(self):

        R1, R2, R3, R4 = self.generate_ONNO_test_rxns()

        self.assertNotEqual(R1, R2)
        self.assertNotEqual(R2, R3)

        self.assertEqual(R1, R3)

    def test_eq_ONNO_labels(self):
        R1, R2, R3, R4 = self.generate_ONNO_test_rxns()

        self.assertEqual(R3, R4)


class ScissionsTestCase(unittest.TestCase):
    """Tests scissions function."""

    def test_self_consistency(self):
        """Test that running scissions multiple times does not find new reactions"""

        ONNO = MolGraph().generate('ONNO', 'smi')

        all_rxns = scissions(ONNO)
        new_rxns = scissions(ONNO, rxns=all_rxns)

        self.assertEqual(all_rxns, new_rxns)


if __name__ == '__main__':
    unittest.main()
