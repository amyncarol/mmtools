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

class PhononFileWriter(object):
    """
    This filewriter follows the phonon calculation workflow detailed in 
    https://atztogo.github.io/phonopy/vasp-dfpt.html

    i.e. using vasp dfpt to perform phonon calculation. Note that vasp uses supercell-gamma approach with 
    dfpt, while QE uses unitcell-generic q points approach

    1) Relax the unitcell to get equilibrium structure (tight criteria)
    2) Perform DFPT calculation with supercell
    3) Phonopy: calculate force constant (needs vasprun.xml from dfpt calculation)
    4) Prepare band.conf
    5) Phonopy: calculate phonon band structure (needs unitcell poscar and band.conf)
    6) Pymatgen: Plot bands


    Folder structure:
        - wd
          INCAR-dfpt
          INCAR-relax
          POSCAR-starting

          - relax
            POSCAR
            CONTCAR

          - dfpt-1

          - dfpt-2 
    """
    def __init__(self, struc, wd, dfpt_folder, scaling_matrix, encut):
        """
        struc: input structure
        wd: working directory
        scaling matrix: E.g., [2, 1, 1] specifies that the supercell should have dimensions 2a x b x c.
        encut: encut
        """
        self.scaling_matrix = scaling_matrix
        self.encut = encut
        self.wd = wd
        self.struc = struc
        self.dfpt_folder = dfpt_folder


    def generate_structure(self, structure):
        supercell = copy.deepcopy(structure)
        supercell.make_supercell(self.scaling_matrix)
        return supercell

    def write_vasp_files_relax(self, incar_path):
        """
        should relax the structure again in tight criteria before running dfpt calculation

        pymatgen inputset may not work very well, so just copy incar from a source
        """
        mpset = MPRelaxSet(self.struc, user_kpoints_settings={'reciprocal_density': 250}, force_gamma=True)
        mpset.write_input(os.path.join(self.wd, 'relax'))
  
        p = subprocess.Popen(['cp', incar_path, 'INCAR'], cwd=os.path.join(self.wd, 'relax'))
        p.wait()

    def write_vasp_files_dfpt(self, incar_path):
        """
        should run only after tight, proper relaxation

        Write vasp files to a folder named dfpt

        pymatgen inputset may not work very well, so just copy incar from a source
        """
        try:
            relaxed_structure = Poscar.from_file(os.path.join(self.wd, 'relax/CONTCAR')).structure
            supercell = self.generate_structure(relaxed_structure)
            mpset = MPStaticSet(supercell, user_kpoints_settings={'reciprocal_density': 250}, force_gamma=True)
            mpset.write_input(os.path.join(self.wd, self.dfpt_folder))

            p = subprocess.Popen(['cp', incar_path, 'INCAR'], cwd=os.path.join(self.wd, self.dfpt_folder))
            p.wait()

        except:
            print('make sure you have relaxed the structure in tight criteria(ediff 1e-8, ediffg 1e-3')

    def write_band_conf(self):
        """
        prepare the input file(band.conf) for phonon band structure
        """
        primi_structure = self.struc.get_primitive_structure()
        primi_lattice = primi_structure.lattice.matrix.T
        lattice = self.struc.lattice.matrix.T
        transformation = np.round(inv(lattice) @ primi_lattice, 5).flatten().tolist()
        t_string = [str(i) for i in transformation]

        kpath = HighSymmKpath(self.struc.get_primitive_structure()).kpath  ##the primitive structure???
        print(primi_structure.get_space_group_info())
        
        kpath_string = 'BAND ='
        print('Please indicate the k-path, disconnect paths are indicated by \',\', we have the following kpoints:')
        print(kpath['kpoints'])
        a = input('Enter your input:')
        self.kpath_name = a
        a_split = a.split(',')
        for i, item in enumerate(a_split):
            for point in item.strip().split(' '):
                kpath_string += ' '
                kpath_string += ' '.join([str(i) for i in kpath['kpoints'][point].tolist()])
            if i < len(a_split)-1:
                kpath_string += ','

        a = a.replace('\Gamma', '$\Gamma$')
        label_string = 'BAND_LABELS = ' + a
       
        with open(os.path.join(self.wd, self.dfpt_folder+'/band.conf'), 'w') as f:
            f.write('DIM = {} {} {}\n'.format(self.scaling_matrix[0], self.scaling_matrix[1], self.scaling_matrix[2]))
            f.write('PRIMITIVE_AXIS = ' + ' '.join(t_string) + '\n')
            f.write(kpath_string + '\n')
            f.write(label_string + '\n')
            f.write('FORCE_CONSTANTS = READ\n')

    def plot_phonon_bands(self, filename = None, ylim=None, units='thz'):
        bs = get_ph_bs_symm_line(os.path.join(self.wd, self.dfpt_folder+'/band.yaml'))
        plotter = PhononBSPlotter(bs)
        if filename==None:
            filename = os.path.join(self.wd, self.dfpt_folder+'/'+self.kpath_name+'.eps')
        plotter.save_plot(filename, ylim=ylim, units=units)

    def get_contcar_file(self):
        return os.path.join(self.wd, 'relax/CONTCAR')


if __name__ == '__main__':
    ## 0) inputs 
    phonon_path = '/Users/yao/Google Drive/data/113/MAPbI3/phonon/1st_try/'
    dfpt_folder = 'dfpt-ps'
    poscar_file = os.path.join(phonon_path, 'POSCAR-starting')  ##before relax
    dfpt_incar = os.path.join(phonon_path, 'INCAR-dfpt') #INCAR for dfpt calculation
    relax_incar = os.path.join(phonon_path, 'INCAR-relax') #INCAR for relax calculation

    structure = Poscar.from_file(poscar_file).structure
    ph = PhononFileWriter(structure, phonon_path, dfpt_folder, [2, 2, 2], 500)
    contcar_file = ph.get_contcar_file()

    ## 1) relax the structure
    #ph.write_vasp_files_relax(relax_incar)

    ## 2) prepare dfpt vasp files
    #ph.write_vasp_files_dfpt(dfpt_incar)

    # 3) create force constant
    p = subprocess.Popen(['phonopy', '--fc', 'vasprun.xml'], cwd=os.path.join(phonon_path, dfpt_folder))
    p.wait()

    # 4) prepare band.conf
    ph.write_band_conf()

    # 5) calculate bands 
    p = subprocess.Popen(['phonopy', '-c', contcar_file, 'band.conf'], cwd=os.path.join(phonon_path, dfpt_folder))
    p.wait()

    # 6) plot bands
    ph.plot_phonon_bands(ylim=[-2, 7], units='mev')
    #ph.plot_phonon_bands(filename=os.path.join(phonon_path, dfpt_folder+'/1.eps'))
