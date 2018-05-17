import unittest
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.io.xyz import XYZ
from mmtools.structure.voronoi import *

class TestColorVonoroi(unittest.TestCase):
        def setUp(self):
                self.struc = Poscar.from_file("/Users/yao/Google Drive/mmtools/data/POSCAR.mp-27544_Cs3Bi2Br9.vasp").structure
                self.cv = ColorVonoroi(self.struc)

        def test_expand_structure(self):
            molecule = ColorVonoroi.expand_structure(self.struc)
            xyz = XYZ(molecule)
            xyz.write_file("/Users/yao/Google Drive/mmtools/data/expaned_structure.xyz")
            #print(molecule)

        def test_get_polyhedron_hull_volume_area(self):
            for i in range(self.cv.n):
                _, vol, area = self.cv.get_polyhedron_hull_volume_area(i)
                print(vol, area)
            
        def test_get_polyhedron_faces(self):
            for i in range(self.cv.n):
                ridge_shape = self.cv.get_polyhedron_faces(i)
                print('for atom {}'.format(i))
                for j in ridge_shape:
                    print(j[0], j[1])
                print('-------------')

        def test_get_distinct_atoms(self):
            print('----------')
            print('the distinct atoms are:')
            print(self.cv.get_distinct_atoms())
            print('----------')

        def test_plot_polyhedron_faces(self):
            self.cv.plot_polyhedron_faces(0)


if __name__ == '__main__':
    unittest.main()