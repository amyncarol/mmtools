from pymatgen.io.vasp.inputs import Poscar
import sys
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

def get_structure(filename, filetype = 'POSCAR'):
	"""
	Read either a poscar file or quantum espresso 'in' file and return the structure object.
	
	Args:
		filename(str) : filename containing Poscar or quantum espresso 'in' file.
		filetype(str) : can be "POSCAR", "QE".

	Returns:
		Structure object.
	"""

	if filetype == 'POSCAR':
		poscar = Poscar.from_file(filename)
		return poscar.structure
	if filetype == 'QE':
		print('QE input to be implemented')

def print_info(struc):
	"""
	Read a Structure object and print structure information.

	Args: 
	    struc : pymatgen Structure object.

	Returns:
		None.
	"""

	print(struc)

	print('\n')

	sga = SpacegroupAnalyzer(struc)

	print('The space group: {} {}\n'.format(sga.get_space_group_number(), sga.get_space_group_symbol()))


	equi_sites = sga.get_symmetrized_structure().equivalent_sites
	for sites in equi_sites:
		print('the symmtry equivalent sites are \n{}\n'.format(sites))

	sym_operations = sga.get_symmetry_operations()
	print('there are {} symmetry operations: \n {} \n'.format(len(sym_operations), sym_operations))

print('\nPlease input filename and filetype. \n For example:\n python structure_summary.py POSCAR_Br POSCAR \n python structure_summary.py in QE\n')
if len(sys.argv)>2:
	structure = get_structure(str(sys.argv[1]), str(sys.argv[2]))
else:
	structure = get_structure(str(sys.argv[1]))
print_info(structure)

