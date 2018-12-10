from pymatgen.io.vasp.inputs import Poscar
from pymatgen.core.lattice import Lattice
from pymatgen.core.sites import Site, PeriodicSite
from pymatgen.core.structure import Structure

def double_unit_cell(struc):
	"""
	doubling the size of the unit cell

	use case: given a perovskite structure, return the corresponding double perovskite structure

	Args: 
		struc: the input structure

	Returns:
	 	the doubled new structure
	"""
	lattice = struc.lattice.matrix
	vector_a = lattice[0] + lattice[1]
	vector_b = lattice[1] + lattice[2]
	vector_c = lattice[0] + lattice[2]
	new_lattice = Lattice([vector_a, vector_b, vector_c])
	
	#the atoms
	new_sites = []
	for site in struc.sites:
		new_sites.append(PeriodicSite(site.species_and_occu, site.coords, new_lattice, coords_are_cartesian=True))
		new_sites.append(PeriodicSite(site.species_and_occu, site.coords+lattice[0], new_lattice, coords_are_cartesian=True))

	species = []
	coords = []
	for site in new_sites:
		species.append(site.specie)
		coords.append(site.coords)

	new_struc = Structure(new_lattice, species, coords, coords_are_cartesian=True)
	return new_struc

if __name__=='__main__':
	folder = '/Users/yao/Google Drive/data/2116/InTl/tetragonal_exp/'
	struc = Poscar.from_file(folder+'POSCAR.vasp').structure
	new_struc = double_unit_cell(struc)
	Poscar(new_struc).write_file(folder+'POSCAR_dp.vasp')



