from pymatgen.io.vasp.inputs import Poscar
from copy import deepcopy
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
import random
from pymatgen.io.vasp.sets import MPRelaxSet
import numpy as np

class SolidSolutionMaker:
    """
    make solid solution between two materials, which have the same structure but different 
    elements and lattice constant, given a certain alloying ratio
    """
    def __init__(self, struc1, struc2, percent, scaling_list=[1,1,1]):
        self.struc1 = struc1
        self.struc2 = struc2
        self.percent = percent
        self.scaling_list = scaling_list

    def get_strain(self):
        """
        get strain relative to the first structure
        ratio: alloying percent of struc1 out of total
        """
        V1 = self.struc1.volume
        V2 = self.struc2.volume
        V = self.percent*V1 + (1-self.percent)*V2
        return (V/V1)**(1/3)-1

    def apply_strain(self):
        """
        get the new structure with the right volume but elements of the first structure
        """
        strain = self.get_strain()
        struc_new = deepcopy(self.struc1)
        struc_new.apply_strain(strain)
        return struc_new

    def get_supercell(self):
        struc = self.apply_strain()
        c_struc = SpacegroupAnalyzer(struc).get_conventional_standard_structure()
        c_struc.make_supercell(self.scaling_list)
        print('Your supercell has {} atoms in total\n'.format(len(c_struc.sites)))
        return c_struc

    def get_mixing_elements(self):
    	s1 = set(self.struc1.species)
    	s2 = set(self.struc2.species)
    	e1 = (s1-s2).pop()
    	e2 = (s2-s1).pop()
    	return e1, e2

    def get_random_supercell(self):
    	"""
    	returns the solid solution Structure object and the actual mixing percent for struc1
    	"""
    	struc = self.get_supercell()
    	sites = [i for i in range(len(struc.sites)) if struc.sites[i].specie==self.get_mixing_elements()[0]]
    	n = int((1-self.percent)*len(sites))
    	print('Now you will replace {} {} atoms with {} atoms, the actual mixing percent for {} is {}\n'. \
    		format(n, self.get_mixing_elements()[0], self.get_mixing_elements()[1], \
    			self.get_mixing_elements()[0], 1-float(n)/len(sites)))
    	random_choose = random.sample(sites, n)
    	for i in random_choose:
    		struc.replace(i, self.get_mixing_elements()[1])
    	return struc, 1-float(n)/len(sites)

class SolidSolutionFileWriter:
    """
    write vasp files ready for running using MPRelaxSet
    """
    def __init__(self, struc1, struc2, wd, percent_list=np.linspace(0.05, 0.95,5), scaling_list=[1,1,1]):
        self.struc1 = struc1
        self.struc2 = struc2
        self.percent_list = percent_list
        self.scaling_list = scaling_list
        self.wd = wd
    def write_vasp_files(self):
        for i in self.percent_list:
            ssm = SolidSolutionMaker(self.struc1, self.struc2, i, scaling_list=self.scaling_list)
            e1, e2 = ssm.get_mixing_elements()
            struc, percent = ssm.get_random_supercell()
            mpset = MPRelaxSet(struc,  user_incar_settings={'EDIFF': 1e-5, 'EDIFFG': -0.01, 'ALGO': 'N', 'ISMEAR': 0})
            mpset.write_input(self.wd+e1.symbol+'_{0:.3f}'.format(percent)+e2.symbol+'_{0:.3f}'.format(1-percent))

    def get_true_percent(self):
        true_percent = []
        for i in self.percent_list:
            ssm = SolidSolutionMaker(self.struc1, self.struc2, i, scaling_list=self.scaling_list)
            e1, e2 = ssm.get_mixing_elements()
            struc, percent = ssm.get_random_supercell()
            true_percent.append(percent)
        return true_percent








    