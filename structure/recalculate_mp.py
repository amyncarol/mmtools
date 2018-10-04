from abc import ABC, abstractmethod
import sys
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.io.vasp.sets import MPRelaxSet, MPHSERelaxSet
from pymatgen.io.cif import CifParser


class FileWriter(ABC):
    """
    the abstract base class for generating a single structure and writes vasp input files, 
    all subclasses should follow this template
    """
    def __init__(self, struc, wd):
        """
        Args:
            struc: the input structure to generate output structure from
            wd: the directory to write vasp input files
        """
        self.struc = struc
        self.wd = wd
        self.out_struc = self.generate_structure()
        super().__init__()

    @abstractmethod
    def generate_structure(self):
        """
        Returns:
            output structure
        """
        pass

    @abstractmethod
    def write_vasp_files(self):
        """
        Writes:
            a folder that includes vasp input file for a single structure
        """
        pass

class RecalculateMPFileWriter(FileWriter):
    """
    recalculate with higher k-point density
    """

    def generate_structure(self):
        return self.struc

    def write_vasp_files(self):
        mpset = MPRelaxSet(self.out_struc, user_incar_settings={'EDIFF': 1e-5, 'EDIFFG': -0.01, 'ALGO': 'F', 'ISMEAR': 0}, \
            user_kpoints_settings={'reciprocal_density': 250})
        mpset.write_input(self.wd)

class RecalculateHSEFileWriter(FileWriter):
    """
    recalculate with HSE 
    """
    def generate_structure(self, primitive=True):
        return self.struc

    def write_vasp_files(self):
        mpset = MPHSERelaxSet(self.out_struc, user_incar_settings={'EDIFF': 1e-5, 'EDIFFG': -0.01, 'ALGO': 'F', 'ISMEAR': 0}, \
            user_kpoints_settings={'reciprocal_density': 250})
        mpset.write_input(self.wd)


if __name__=='__main__':
    # poscar_input = "/Users/yao/Google Drive/data/lattice_energy_testset/mp-22862_NaCl/POSCAR"
    # wd = "/Users/yao/Google Drive/data/lattice_energy_testset/mp-22862_NaCl_recal"

    # poscar_input = "/Users/yao/Google Drive/data/lattice_energy_testset/mp-22865_CsCl/POSCAR"
    # wd = "/Users/yao/Google Drive/data/lattice_energy_testset/mp-22865_CsCl_recal"

    # poscar_input = "/Users/yao/Google Drive/data/phonopy/NaCl/NaCl-vasp-dfpt-exercise/POSCAR"
    # wd = "/Users/yao/Google Drive/data/phonopy/NaCl/NaCl-vasp-dfpt-exercise"
    # struc = Poscar.from_file(poscar_input).structure

    cif_input = "/Users/yao/Google Drive/data/2116/solidsolution/structures_exp_guess/inorganic/EntryWithCollCode30657_Cs2AgCl3.cif"
    wd = "/Users/yao/Google Drive/data/2116/solidsolution/structures_exp_guess/inorganic/Cs2AgCl3"
    struc = CifParser(cif_input).get_structures()[0]
    filewriter = RecalculateMPFileWriter(struc, wd)
    filewriter.write_vasp_files()



