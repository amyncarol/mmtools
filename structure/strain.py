from copy import deepcopy
from pymatgen.io.vasp.sets import MPRelaxSet

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
		pass
		struc_new = deepcopy(self.struc)
		struc_new.apply_strain(self.strain)
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





		