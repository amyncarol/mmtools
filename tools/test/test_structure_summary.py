import unittest
from mmtools.tools.structure_summary import print_info
from pymatgen.io.vasp.inputs import Poscar


class Test(unittest.TestCase):
	def setUp(self):
		self.struc = Poscar.from_file('../../data/POSCAR_Br').structure

	def test_print_info(self):
		print_info(self.struc)

if __name__ == '__main__':
    unittest.main()