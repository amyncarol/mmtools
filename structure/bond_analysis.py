def bond_diff_analysis(s1, s2, atoms, bond_length, threshold = 0.03):
	"""
	Given two structures with slight difference(usually two different structures from two 
	relaxation), the purpose is to compare which bonds are longer which are shorter. Print 
	all the information. 

	Args:
		atoms: a list of atom numbers(remember in pymatgen numbering starts from 0). We need
		to find all bonds around these atoms. This is by define a length, and find all 
		neighbors insides the circle of this length. Both the neighbor positions and actual 
		bond lengths can be obtained. 

		bond_length(angstrom): the length within which all neighbors need to be found. Typically this is 
		slight larger than your target bond length so that you can find all streched or compressed 
		bonds.

		threshold(angstrom): only print the bond difference that are larger than some threshold. We don't 
		consider very small difference as real difference.

	Returns:
		None

	Prints:
		The number of found bonds in both structures, should make sure they are equal so that we don't 
		miss any bonds. 

		Should print all bonds(site-site pairs), and their length in structure1, and the difference in 
		them between the two slightly different structures. 
	"""

	pass