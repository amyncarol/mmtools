from recalculate_mp import FileWriter
import copy
import os
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.io.vasp.sets import MPRelaxSet
from pymatgen.io.vasp.outputs import Vasprun

class PhononFileWriter(object):
    def __init__(self, struc, wd, scaling_matrix, encut):
        """
        scaling matrix:
        a) A full 3x3 scaling matrix defining the linear combination the old lattice vectors. E.g., [[2,1,0],[0,3,0],[0,0, 1]] generates a new structure with lattice vectors a’ = 2a + b, b’ = 3b, c’ = c where a, b, and c are the lattice vectors of the original structure.
        b) An sequence of three scaling factors. E.g., [2, 1, 1] specifies that the supercell should have dimensions 2a x b x c.
        c) A number, which simply scales all lattice vectors by the same factor.
        """
        self.scaling_matrix = scaling_matrix
        self.encut = encut
        self.wd = wd
        self.struc = struc

    def generate_structure(self, structure):
        supercell = copy.deepcopy(structure)
        supercell.make_supercell(self.scaling_matrix)
        return supercell

    def write_vasp_files_relax(self):
        """
        should relax the structure again in tight criteria before running dfpt calculation
        """
        mpset = MPRelaxSet(self.struc, user_incar_settings={'GGA': 'PS', \
                                                            'ENCUT': self.encut, 'EDIFF': 1e-8, 'ISMEAR': 0,\
                                                            'EDIFFG': -0.001,\
                                                            'PREC': 'Accurate', 'LREAL': False, 'ADDGRID': True}, \
                 user_kpoints_settings={'reciprocal_density': 250}, force_gamma=True)
        mpset.write_input(os.path.join(self.wd, 'relax'))


    def write_vasp_files_dfpt(self):
        """
        should run only after tight, proper relaxation
        """
        try:
            #relaxed_structure = Poscar.from_file(os.path.join(self.wd, 'relax/CONTCAR')).structure
            relaxed_structure = self.struc
            supercell = self.generate_structure(relaxed_structure)
            mpset = MPRelaxSet(supercell, user_incar_settings={'GGA': 'PS', \
                                                        'ENCUT': self.encut, 'EDIFF': 1e-8, 'ALGO': 'N', 'ISMEAR': 0, \
                                                        'IBRION': 8, 'NSW': 1,\
                                                        'PREC': 'Accurate', 'LREAL': False, 'ADDGRID': True, \
                                                        'LWAVE': False, 'LCHARG': False}, \
                user_kpoints_settings={'reciprocal_density': 250}, force_gamma=True)
            mpset.write_input(os.path.join(self.wd, 'dfpt'))
        except:
            print('make sure you have relaxed the structure in tight criteria(ediff 1e-8, ediffg 1e-3')


if __name__ == '__main__':
    structure = Poscar.from_file('/Users/yao/Google Drive/data/113/MAPbI3/phonon/MAPbI3_walsh.vasp').structure
    ph = PhononFileWriter(structure, '/Users/yao/Google Drive/data/113/MAPbI3/phonon/', 2, 700)
    #ph.write_vasp_files_relax()

    ph.write_vasp_files_dfpt()
