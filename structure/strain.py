from copy import deepcopy

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





		