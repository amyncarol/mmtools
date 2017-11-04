from copy import deepcopy
from pymatgen.io.vasp.sets import MPRelaxSet
from pymatgen.io.vasp.inputs import Kpoints
class StrainMaker(object):
	"""
	this class makes a strained structure from the original structure
	"""
	def __init__(self, struc, strain):
		"""
		struc: a strucutre object from pymatgen
		strain: the percent strain applied
		"""
		self.struc = struc
		self.strain = strain

	def get_strained_structure(self):
		struc_new = deepcopy(self.struc)
		struc_new.apply_strain(self.strain)
		return struc_new

class StrainFileWriter(object):
	"""
	this writes the vasp files for a series of strain calculation
	"""
	def __init__(self, struc, strain_list, wd):
		"""
		wd: working directory
		"""
		self.struc = struc
		self.strain_list = strain_list
		self.wd = wd

	def write_vasp_files(self):
		for i in self.strain_list:

			sm = StrainMaker(self.struc, i)

			struc_new = sm.get_strained_structure()

			mpset = MPRelaxSet(struc_new,  user_incar_settings={'EDIFF': 1e-5, 'EDIFFG': -0.01, 'ALGO': 'N', 'ISMEAR': 0, 'ISIF': 2},\
				user_kpoints_settings={"reciprocal_density": 100}, force_gamma=True)

			mpset.write_input(self.wd + '/strain_' + '{0:.3f}'.format(i))


		