import sys
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.io.vasp.sets import MPRelaxSet, MPHSERelaxSet, MPStaticSet
from recalculate_mp import FileWriter
from pymatgen.io.xyz import XYZ
import glob

class MoleculeVacuumFileWriter(FileWriter):
    """
    this class takes a molecule structure and generate vasp input files for it
    """

    def __init__(self, struc, wd, cube_length):
        self.cube_length = cube_length
        super().__init__(struc, wd)

    def generate_structure(self):
        """
        takes the input molecule structure and output a structure for poscar
        """
        a = b = c = self.cube_length
        return self.struc.get_boxed_structure(a, b, c)

    def write_vasp_files(self):
        """
        Writes:
            a folder that includes vasp input file for a single structure, this is for relaxation
        """
        mpset = MPRelaxSet(self.out_struc, user_incar_settings={'EDIFF': 1e-5, 'EDIFFG': 0.001, 'ALGO': 'F', 'ISMEAR': 0, 'ISIF': 2, 'IBRION': 2, 'NSW': 500}, \
            user_kpoints_settings={'reciprocal_density': 1})
        mpset.write_input(self.wd)

    def write_vasp_files_static(self):
        """
        Writes:
            a folder that includes vasp input file for a single structure
        """
        mpset = MPStaticSet(self.out_struc, reciprocal_density=1, user_incar_settings={'EDIFF': 1e-5, 'ALGO': 'F', 'ISMEAR': 0})
        mpset.write_input(self.wd)

def xyz2struc(xyz_file):
    """
    read xyz file and return molecule structure
    """
    return XYZ.from_file(xyz_file).molecule

if __name__=='__main__':
    folder = '/Users/yao/Google Drive/data/2116/solidsolution/structures_exp_guess/modfiles'
    xyz_list = glob.glob(folder+'/*.xyz')
    for xyz in xyz_list:
        writer = MoleculeVacuumFileWriter(xyz2struc(xyz), xyz[:-4], 15)
        writer.write_vasp_files()
    

