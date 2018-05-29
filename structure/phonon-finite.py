from recalculate_mp import FileWriter
import copy
import os
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.io.vasp.sets import MPRelaxSet, MPStaticSet
from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from numpy.linalg import inv
import numpy as np
from pymatgen.symmetry.bandstructure  import HighSymmKpath
from pymatgen.io.phonopy import *
from pymatgen.phonon.plotter import PhononBSPlotter
import subprocess
from glob import glob

class PhononFiniteFileWriter(object):
    """
    This implements the workflow of phonon calculation using finite displacement. 

    Using Vasp and Phonopy.

    1) Create supercell with displacements
    2) Perform vasp static calculations for the displaced structures
    3) Band calculation using phonopy 

    Folder structure:
        - wd
          POSCAR-starting
          INCAR-relax
          INCAR-finite
          POSCAR (the relaxed structure)
          disp.yaml
          POSCAR-001
          POSCAR-002
          ....

          - relax
            POSCAR
            CONTCAR

          - disp-001
            POSCAR-001
            POSCAR

          - disp-002
            POSCAR-002
            POSCAR

          - disp-003
            POSCAR-003
            POSCAR

          .....

    """
    def __init__(self, wd, scaling_matrix, struc=None):
        if struc:
            self.struc = struc  ## the unrelaxed structure
        self.wd = wd
        self.scaling_matrix = scaling_matrix

    def generate_structure(self):
        p = subprocess.Popen(['cp', 'relax/CONTCAR', 'POSCAR'], cwd = self.wd)
        p.wait()
        scaling_matrix = ' '.join([str(i) for i in self.scaling_matrix])
        p = subprocess.Popen(['phonopy', '-d', '--dim='+scaling_matrix], cwd = self.wd)
        p.wait()

    def write_vasp_files_disp(self, incar_path):
        try: 
            poscar_list = glob(self.wd+'/POSCAR-[0-9][0-9][0-9]')
            print('{} displacement poscar in total'.format(len(poscar_list)))

            for poscar in poscar_list:
                folder = 'disp-'+poscar.split('-')[-1]
                structure = Poscar.from_file(poscar).structure
                mpset = MPStaticSet(structure, user_kpoints_settings={'reciprocal_density': 250}, force_gamma=True)
                mpset.write_input(os.path.join(self.wd, folder))
                p = subprocess.Popen(['cp', incar_path, 'INCAR'], cwd=os.path.join(self.wd, folder))
                p.wait()
        except:
            print('No displacement poscars')

if __name__=='__main__':
    phonon_path = '/Users/yao/Google Drive/data/113/MAPbI3/phonon/finite_1st_try'
    ph = PhononFiniteFileWriter(phonon_path, [2, 2, 2])
    #ph.generate_structure()
    incar_path = os.path.join(phonon_path, 'INCAR-finite')
    ph.write_vasp_files_disp(incar_path)

