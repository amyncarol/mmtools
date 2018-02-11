import unittest
from pymatgen.io.vasp.inputs import Poscar
from mmtools.structure.strain import StrainMaker

class TestStrainMaker(unittest.TestCase):
	def setUp(self):
		self.struc = Poscar.from_file('../data/POSCAR_Br').structure

	def test_get_strained_structure(self):
		sm = StrainMaker(self.struc, 0)
		new_struc = sm.get_strained_structure()
		self.assertAlmostEqual(new_struc.lattice.a, 8.123452, delta=1e-5)
		self.assertAlmostEqual(self.struc.lattice.a, 8.123452, delta=1e-5)
		sm = StrainMaker(self.struc, 0.5)
		new_struc = sm.get_strained_structure()
		self.assertAlmostEqual(new_struc.lattice.a, 12.185178, delta=1e-5)
		self.assertAlmostEqual(self.struc.lattice.a, 8.123452, delta=1e-5)
		sm = StrainMaker(self.struc, -0.1)
		new_struc = sm.get_strained_structure()
		self.assertAlmostEqual(new_struc.lattice.a, 7.3111068, delta=1e-5)
		self.assertAlmostEqual(self.struc.lattice.a, 8.123452, delta=1e-5)



if __name__ == '__main__':
    unittest.main()


