import sys
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.core.structure import Structure
from pymatgen.io.vasp.sets import MPRelaxSet
from numpy.linalg import norm
import numpy as np

class IonVacuumFileWriter(object):
    """
    write vasp files ready for running using MPRelaxSet
    """
    def __init__(self, struc1, center_index, neighboring_index=None, vacuum=20, max_distance=4, charge=0):
        """
        Args:
            struc1: the original crystal structure
            center_index: the index of the center atom in original structure (starting from 1, as used by vesta)
            neighboring_index: a list of index for neighboring atoms (starting from 1, as used by vesta)
            vacuum: dimension of the vacuum(cubic)
            max_distance: the max distance between center and neighboring atoms
            charge: the net charge of the ion
        """
        self.struc1 = struc1
        self.charge = charge
        self.center_index = center_index
        self.neighboring_index = neighboring_index
        self.vacuum = vacuum
        self.max_distance = max_distance
        self.struc2 = self.generate_structure()

    def generate_structure(self):
        """
        select the center atom and several neighboring atoms, generate
        a structure with the cluster in the center surrounded by vacuum

        Returns:
            structure with cluster in vacuum
        """
        sites = self.struc1.sites
        lattice_matrix = self.struc1.lattice.matrix
        cart_coords = []
        species = []
        new_lattice = [[self.vacuum, 0, 0], [0, self.vacuum, 0], [0, 0, self.vacuum]]

        ##the center atom
        center_coords = sites[self.center_index-1].coords
        cart_coords.append(center_coords)
        species.append(sites[self.center_index-1].specie)

        ##the neighboring atoms
        translation_vectors = []
        for a in [-lattice_matrix[0, :], np.zeros(3), lattice_matrix[0, :]]:
        	for b in [-lattice_matrix[1, :], np.zeros(3), lattice_matrix[1, :]]:
        		for c in [-lattice_matrix[2, :], np.zeros(3), lattice_matrix[2, :]]:
        			translation_vectors.append(a+b+c)
 
        if self.neighboring_index != None:
            for index in self.neighboring_index:
                for vector in translation_vectors:
                        translated_coords = sites[index-1].coords+vector
                        if norm(translated_coords-center_coords) < self.max_distance:
                            cart_coords.append(translated_coords)
                species.append(sites[index-1].specie)
        
        new_struc = Structure(new_lattice, species, cart_coords, to_unit_cell=False, coords_are_cartesian=True)
        return new_struc

    def write_vasp_files(self, wd, neutral=False):
        """
        Args:
            neutral: whether the ion is neutral or not
            wd: working directory
        """
        ##neutral 
        mpset_neutral = MPRelaxSet(self.struc2, user_incar_settings={'EDIFF': 1e-5, 'EDIFFG': -0.01, 'ALGO': 'F', 'ISMEAR': 0, \
            'ISIF': 2})
        if neutral or self.charge == 0:
            mpset_neutral.write_input(wd)
        else:
            ##charged
            nelect = mpset_neutral.nelect - self.charge
            mpset = MPRelaxSet(self.struc2,  user_incar_settings={'EDIFF': 1e-5, 'EDIFFG': -0.01, 'ALGO': 'F', 'ISMEAR': 0, \
                'ISIF': 2, 
                'NELECT': nelect, 'IDIPOL': 4, 'DIPOL': self._get_dipol(), 'LDIPOL': True})
            mpset.write_input(wd)
      

    def _get_dipol(self):
        """
        get the position where the charge is in the supercell in direct coords
        """
        frac_coords = np.array([i.frac_coords for i in self.struc2.sites])
        center_mass = frac_coords.mean(axis=0)
        return "{} {} {}".format(center_mass[0], center_mass[1], center_mass[2])


if __name__ == '__main__':
    ##PbBr3
    # POSCAR_input = '/Users/yao/Google Drive/data/lattice_energy_testset/mp-567629_CsPbBr3/POSCAR.vasp'
    # wd = '/Users/yao/Google Drive/data/lattice_energy_testset/PbBr3_ion'
    # wd_neutral = '/Users/yao/Google Drive/data/lattice_energy_testset/PbBr3_neutral'
    # center_index = 5
    # neighboring_index = [9, 15, 18]
    # charge = -1

    ##InCl3
    POSCAR_input = '/Users/yao/Google Drive/data/lattice_energy_testset/Cs2Ag1In1Cl6/POSCAR'
    wd = '/Users/yao/Google Drive/data/lattice_energy_testset/InCl3_ion'
    center_index = 3
    neighboring_index = [6, 7, 8]
    charge = 0

    ##InCl4
    # POSCAR_input = '/Users/yao/Google Drive/data/lattice_energy_testset/Cs2Ag1In1Cl6/POSCAR'
    # wd = '/Users/yao/Google Drive/data/lattice_energy_testset/InCl4_ion'
    # center_index = 3
    # neighboring_index = [5, 6, 7, 8]
    # charge = -1

    # ##InCl5
    # POSCAR_input = '/Users/yao/Google Drive/data/lattice_energy_testset/Cs2Ag1In1Cl6/POSCAR'
    # wd = '/Users/yao/Google Drive/data/lattice_energy_testset/InCl5_ion'
    # center_index = 3
    # neighboring_index = [5, 6, 7, 8, 9]
    # charge = -2


    poscar = Poscar.from_file(POSCAR_input)
    struc = poscar.structure
    filewriter = IonVacuumFileWriter(struc1=struc, center_index=center_index, neighboring_index=neighboring_index, charge=charge)
    filewriter.write_vasp_files(wd = wd, neutral=False)


