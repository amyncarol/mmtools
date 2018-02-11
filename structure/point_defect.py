class PointDefect(object):
	"""
	From the input structure, this class creates structures with point defects, which includes:

	1\ substitution
	2\ antisite
	3\ vacancy

	"""
	def __init__(self, struc, supercell_multiples):
		"""
		Args:

			struc: a strucutre object from pymatgen

			supercell: a three-element list that indicates the multiples on lattice vector a, b, c
			i.e. [1, 2, 2]

		Attributes:

			supercell: the Structure object for the supercell

			conventional_structure: the conventional structure that are corresponding to this structure,
			some times they are different, and in some cases in may be better to use the conventional 
			structure

		"""
		self.struc = struc
		self.supercell = 
		self.conventional_structure = 

	def substitution(self, element1, element2):
		"""
		Create substitution point defects by replacing element1 with element2
		Ideally, this should only return the defected structure that are symmetrically distinct from each other.

		Args:
			element1: a string represents the element that is to be replaced
			element2: a string represents the doping element

		Returs:
			A list of structures with substitution point defect.
		"""
		pass

	def antisite(self, element1, element2):
		"""
		Create antisite point defects by exchanging element1 and element2
		Ideally, this should only return the defected structure that are symmetrically distinct from each other.

		Args:
			element1: a string represents the element
			element2: a string represents the other element

		Returs:
			A list of structures with antisite point defect.
		"""
		pass

	def vacancy(self, element):
		"""
		Create vacancy point defects by removing one of the atoms
		Ideally, this should only return the defected structure that are symmetrically distinct from each other.

		Args:
			element1: a string represents the element to be removed

		Returs:
			A list of structures with vacancy point defect.
		"""
		pass



