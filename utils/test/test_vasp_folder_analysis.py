import unittest
from pymatgen.io.vasp.outputs import Vasprun
from pymatgen import Element
import os
import sys
sys.path.insert(0, '../')
from vasp_folder_analysis import VasprunAnalyzer

class TestVasprunAnalyzer(unittest.TestCase):
    def setUp(self):
        self.vr = Vasprun('/Users/yao/Google Drive/data/2116/solidsolution/Cs2Ag1Sb1Cl6/vasprun.xml')
        self.va = VasprunAnalyzer(self.vr)

    def test_get_pd_with_open_element(self):
        # gcpds, min_chempots, avg_chempots, max_chempots = self.va.get_pd_with_open_element(Element('Ag'))
        # print(min_chempots)
        # print(avg_chempots)
        # print(max_chempots)
        pass

    def test_get_stable_and_unstable_with_open_element(self):
        #self.va.get_stable_and_unstable_with_open_element(Element('Ag'), 'Cs2Ag1Bi1Cl6_Ag_open')
        pass

    def test_get_chempot_range(self):
        #print(self.va.get_chempot_range(Element('Ag')))
        pass

    def test_get_equilibrium_reaction_energy(self):
        pass
        #print(self.va.get_equilibrium_reaction_energy())

    def test_get_stable_range(self):
        print(self.va.get_stable_range(Element('Ag')))


if __name__ == '__main__':
    unittest.main()