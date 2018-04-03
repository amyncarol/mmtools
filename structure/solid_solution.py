from pymatgen.io.vasp.inputs import Poscar
from copy import deepcopy
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
import random
from pymatgen.io.vasp.sets import MPRelaxSet
import numpy as np
from pymatgen.analysis.structure_matcher import *

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

    def can_match(self):
    	"""
    	judge whether struc1 and struc2 match or not. Only consider if their structures are compatible.

    	Returns:
    		boolean
    	"""
    	matcher = StructureMatcher(comparator=FrameworkComparator())
    	return matcher.fit(self.struc1, self.struc2)

    def can_mix(self, radii_diff = 0.2):
    	"""
    	judge whether struc1 and struc2 can mix or not. 

    	1\ Consider structure matching. 

        2\ We only allow mixing at one set of symmetrically equivalent sites,
         at these sites, the mixing element should have similar atomic radii.
         At all other site, the element should be the same.

        Args: 
            radii_diff(Angstrom): the difference in atomic radii should be less than radii_diff. 
    
    	Returns:
    		boolean
    	"""
    	if not self.can_match():
    		return False
    	else:
            sym_struc1 = SpacegroupAnalyzer(self.struc1).get_symmetrized_structure()
            sym_struc2 = SpacegroupAnalyzer(self.struc2).get_symmetrized_structure()

            sites_searched = []
            num_slight_diff = 0
            for site1, site2 in zip(sym_struc1.sites, sym_struc2.sites):
                if site1 not in sites_searched:
                    if site1.specie.symbol != site2.specie.symbol:
                        abs_diff = abs(site1.specie.atomic_radius - site2.specie.atomic_radius)
                        #print(site1.specie.atomic_radius, site2.specie.atomic_radius)
                        #print(abs_diff)
                        if abs_diff < radii_diff:
                            num_slight_diff += 1
                        else:
                            return False
                sites_searched += sym_struc1.find_equivalent_sites(site1)
            if num_slight_diff == 1:
                return True
            else:
                return False

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

    def get_supercell(self, isprint=True):
        struc = self.apply_strain()
        c_struc = SpacegroupAnalyzer(struc).get_conventional_standard_structure()
        c_struc.make_supercell(self.scaling_list)
        if isprint:
        	print('Your supercell has {} atoms in total\n'.format(len(c_struc.sites)))
        return c_struc

    def get_mixing_elements(self):
    	s1 = set(self.struc1.species)
    	s2 = set(self.struc2.species)
    	e1 = (s1-s2).pop()
    	e2 = (s2-s1).pop()
    	return e1, e2

    def get_random_supercell(self, isprint=True):
    	"""
    	returns the solid solution Structure object and the actual mixing percent for struc1
    	"""
    	struc = self.get_supercell(isprint=isprint)
    	sites = [i for i in range(len(struc.sites)) if struc.sites[i].specie==self.get_mixing_elements()[0]]
    	n = int((1-self.percent)*len(sites))
    	if isprint:
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
            if percent == 1.0 or percent == 0.0:
                pass
            else:
                mpset = MPRelaxSet(struc,  user_incar_settings={'EDIFF': 1e-5, 'EDIFFG': -0.01, 'ALGO': 'F', 'ISMEAR': 0})
                mpset.write_input(self.wd+'_'+e1.symbol+'_{0:.3f}'.format(percent)+'_'+e2.symbol+'_{0:.3f}'.format(1-percent))

    def get_true_percent(self):
        true_percent = []
        for i in self.percent_list:
            ssm = SolidSolutionMaker(self.struc1, self.struc2, i, scaling_list=self.scaling_list)
            e1, e2 = ssm.get_mixing_elements()
            struc, percent = ssm.get_random_supercell()
            true_percent.append(percent)
        return true_percent








    