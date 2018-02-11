
class BondChanger(object):
	"""
	This class includes methods that changes the some bond lengths of a structure, which includes:

	1\ Keep symmtry, and change bond length, but not bond direction. Why we do this? For example, in band structure analysis, 
	sometimes we want to see hwo does a change in bond length affect the band structure; but we don't 
	want to change the symmetry since that will change the band structure symmetry as well.

	2\ Keep symmetry, and rotate the buliding block(tetrahedra, octahedra, etc). 

	3\ Don't worry about symmetry, stretch a bond or compress a bond.

	"""
	def __init__(self, struc):
		"""
		Args:

			struc: a strucutre object from pymatgen

		Attributes:

			spacegroup: the space group symbol of this structure
			
			free_lattice_parameters: a dictionary of free lattice parameters. For example, {'a': 4.00, 'c': 3.44},
			which means this structure have only two free lattice parameters, which are a and c. We don't consider 
			angles, since usually we don't change angles. 

			free_coordinates: a nested dictionary of free coordinates. For example, {2: {x: 3.00, y: 4.00}}, 
			which means no.2 site(numbering starts from 0) has free coordinates whose change don't change the space group and we thus can 
			change them, its x and y coordinates are free, but z coordinate is not.
		
		"""
		self.struc = struc
		self.spacegroup = 
		self.free_lattice_parameters = 
		self.free_coordinates = 

	def symmetry_stretch(self, bond_to_change, bonds_to_keep, lattice_constants, stretching_percent):
		"""
		Keep symmtry, and change bond length, but not bond direction(of the bonds we care about, other bonds may change 
		directions). 

		Args:
			bond_to_change: a tuple (site1, site2), the bond(and all equivalent bonds) we want to stretch or compress.

				site1: the number of the site on one side of the bond, this site usually has a high symmetry than site2

				site2: the number of the site on the other side of the bond, this site usually has a lower symmetry and has 
				free coordinates that we can tune to stretch the bond.

			bonds_to_keep: a list of tuples [(site3, site4), (site5, site6)], the bonds(and all equivalent bonds) 
			whose length we want to keep unchanged.

			lattice_constants: a list of lattice constants that can be changed to stretch to bond. For example, ['a', 'c']

			stretching_percent: a list of stretching percentage applied to the bond, can be positive or negative.

		Returns:

			A list of structures with stretched bonds, corresponding to each stretching percent.


		
		"""
		pass

	def symmetry_rotation(self, bond_to_rotate, normal_direction, angles):

		"""
		Keep symmetry, and rotate the buliding block(tetrahedra, octahedra, etc). 
		
		Args:
			bond_to_rotate: a tuple (site1, site2), the bond(and all equivalent bonds) we want to rotate by some degree.

				site1: the number of the site on one side of the bond, this site usually has a high symmetry than site2

				site2: the number of the site on the other side of the bond, this site usually has a lower symmetry and has 
				free coordinates that we can tune to rotate the bond.

			normal_direction: the normal direction of the bond rotation in crystal coordinates. For example: [0, 0, 1].

			angles: a list of angles in degree by which we want to rotate.

		Returns:

			A list of structures with rotated bonds, corresponding to each angle.
		


		"""
		pass



	def bond_stretch(self, sites, bond_length, stretching_percent):
	"""
	Given a list of atom sites and a bond_length, stretch the bond by moving the site, return a list of structures.

	Args:
		sites: a list of atom sites(pymatgen sites numbering starts from 0). Usually the sites are symmetrically 
		equivalent sites.

		bond_length: the bond_length of targeted bond

		stretching_percent: a list of percentage by which we should stretch the bonds. Can be positive or negative.

	Returns:
		A list of structures with different stretching_percent.
	"""

		pass



