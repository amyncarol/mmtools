from copy import deepcopy
from pymatgen.io.vasp.sets import MPRelaxSet
from pymatgen.io.vasp.inputs import Poscar
from shutil import copyfile, move, rmtree
from glob import glob
import os
import numpy as np

class StrainMaker(object):
    """
    This class makes a strained structure or strained structures from the original structure
    """
    def __init__(self, struc):
        """
        Args:

        struc: a strucutre object from pymatgen

        """
        self.struc = struc

    def iso_strain(self, strain):
        """
        return a structure with isotropic strain

        Args:

        strain: a list of float numbers indicating the percent strain applied, 
        can be positive(tensile) or negative(compress)

        Returns:

        A list of Structure objects corresponding to each strain

        """
        struc_new = deepcopy(self.struc)
        struc_new.apply_strain(strain)
        return struc_new

    def biaxial_strain(self, strain, direction):
        """
        return a structure with biaxial strain

        Args:

        strain: a list of float numbers indicating the percent strain applied, 
        can be positive(tensile) or negative(compress)

        direction: strain applied to a plane, and direction is the normal direction of that plane

        Returns:
        
        A list of Structure objects corresponding to each strain
        """
        pass

class StrainFileWriter(object):
    """
    this writes the vasp files for a series of strain calculation
    """
    def __init__(self, struc, strain_list, wd, copy_input):
        """
        wd: working directory
        copy_input: if true, copy the inputs from wd
        """
        self.strain_list = strain_list
        self.wd = wd
        self.sm = StrainMaker(struc)
        self.copy_input = copy_input

    def write_vasp_files(self):
        #clean the folder
        for path in glob(self.wd+'/strain_*'):
            rmtree(path)

        #prepare inputs
        for strain in self.strain_list:
            struc_new = self.sm.iso_strain(strain)
            subfolder = self.wd + '/strain_' + '{0:.3f}'.format(strain)

            if self.copy_input:
                os.mkdir(subfolder)
                print(subfolder)
                copyfile(self.wd+'/KPOINTS', subfolder+'/KPOINTS')
                copyfile(self.wd+'/INCAR', subfolder+'/INCAR')
                copyfile(self.wd+'/POTCAR', subfolder+'/POTCAR')
                Poscar(struc_new).write_file(subfolder+'/POSCAR')

            else:
                mpset = MPRelaxSet(struc_new,  user_incar_settings={'EDIFF': 1e-5, 'EDIFFG': -0.01, 'ALGO': 'N', 'ISMEAR': 0, 'ISIF': 2},\
                user_kpoints_settings={"reciprocal_density": 100}, force_gamma=True)

                mpset.write_input(subfolder)

if __name__=='__main__':
    wd = '/Users/yao/Google Drive/data/2116/InTl/tetragonal_exp/10atom_strain_HSE'
    struc = Poscar.from_file(os.path.join(wd, 'POSCAR')).structure
    strain_list = np.linspace(-0.06, 0, 8)
    sfw = StrainFileWriter(struc, strain_list, wd, copy_input = True)
    sfw.write_vasp_files()




        