import unittest
from pymatgen.io.vasp.inputs import Poscar
from SolidSolutionMaker import SolidSolutionMaker
import os

class TestSolidSolutionMaker(unittest.TestCase):
    def setUp(self):
        self.struc1 = Poscar.from_file('test/POSCAR_Br').structure
        self.struc2 = Poscar.from_file('test/POSCAR_Cl').structure
        self.a1 = self.struc1.lattice.a
        self.a2 = self.struc2.lattice.a

    def test_get_strain(self):
        ssm = SolidSolutionMaker(self.struc1, self.struc2, 1)
        self.assertEqual(ssm.get_strain(), 0)
        ssm = SolidSolutionMaker(self.struc1, self.struc2, 0)
        self.assertAlmostEqual(ssm.get_strain(), (self.a2-self.a1)/self.a1, delta=1e-6)

    def test_apply_strain(self):
        ssm = SolidSolutionMaker(self.struc1, self.struc2, 1)
        self.assertAlmostEqual(ssm.apply_strain().lattice.a, self.a1, delta=1e-6)
        ssm = SolidSolutionMaker(self.struc1, self.struc2, 0) 
        self.assertAlmostEqual(ssm.apply_strain().lattice.a, self.a2, delta=1e-6) 

    def test_get_supercell(self):
        ssm = SolidSolutionMaker(self.struc1, self.struc2, 1, [2,2,1])
        ssm.get_supercell()

    def test_get_mixing_elements(self):
        ssm = SolidSolutionMaker(self.struc1, self.struc2, 1)
        self.assertEqual(ssm.get_mixing_elements()[0].symbol, 'Br')
        self.assertEqual(ssm.get_mixing_elements()[1].symbol, 'Cl')

    def test_get_random_supercell(self):
        ssm = SolidSolutionMaker(self.struc1, self.struc2, 0.1)
        ssm.get_random_supercell()
        ssm = SolidSolutionMaker(self.struc1, self.struc2, 0.3)
        ssm.get_random_supercell()
        ssm = SolidSolutionMaker(self.struc1, self.struc2, 1)
        ssm.get_random_supercell()
        
if __name__ == '__main__':
    unittest.main()
