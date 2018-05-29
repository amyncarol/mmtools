from pymatgen.io.vasp.inputs import Poscar
from scipy.stats import special_ortho_group
import numpy as np
from numpy.linalg import norm
import copy
from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure, Molecule
from pymatgen.io.vasp.sets import MPRelaxSet
from pymatgen.io.xyz import XYZ
import os

class RandomOrientFileWriter(object):
	def __init__(self, struc, ion_indices, center_index, scaling_matrix, wd):
		"""
		Args:
			struc: the unit cell
			ion_indices: a list of indices of atoms in the molecule to be rotated
			center_index: provide a index of an atom in the molecule to make molecule the center of the unit cell
			scaling_matrix: i.e. [2, 2, 2],  a 2*2*2 supercell
			wd: the working directory
		"""
		self.struc = self.re_center(struc, center_index)
		self.ion_indices = ion_indices
		self.scaling_matrix = scaling_matrix
		self.molecule, self.others = self.seperate_ion()
		self.wd = wd
		
	def re_center(self, struc, center_index):
		"""
		translate all atoms to make the molecule in the center of the unit cell
		"""
		center_frac_coords = struc.sites[center_index].frac_coords
		struc.translate_sites([i for i in range(len(struc.sites))], vector = np.array([0.5, 0.5, 0.5])-center_frac_coords)
		return struc

	def seperate_ion(self):
		"""
		given a structure and a list of atom indices(of a molecule), return cartesian coordinates of the molecule
		ions and cartesian coordinates of the remaining aomts. 
		"""
		sites = self.struc.sites
		molecule = [(sites[i].specie, sites[i].coords) for i in self.ion_indices]
		others = [(sites[i].specie, sites[i].coords) for i in range(len(sites)) if i not in self.ion_indices]
		return molecule, others

	def generate_structure(self):
		"""
		generate a supercell structure with random oriented molecules
		"""
		molecule_array = np.array([atom[1] for atom in self.molecule])
		others_array = np.array([atom[1] for atom in self.others])
		center = np.mean(molecule_array, axis = 0)
		lattice = self.struc.lattice.matrix
		a, b, c = lattice[0], lattice[1], lattice[2]
		lattice = np.array([a*self.scaling_matrix[0], b*self.scaling_matrix[1], c*self.scaling_matrix[2]])
		lattice = Lattice(lattice)
		species = []
		coords = []
		for i in range(self.scaling_matrix[0]):
			for j in range(self.scaling_matrix[1]):
				for k in range(self.scaling_matrix[2]):
					rotation_matrix = special_ortho_group.rvs(3)
					rotated_molecule = (rotation_matrix @ (molecule_array-center).T).T + center
					for m, atom in enumerate(self.molecule):
						species.append(atom[0])
						coords.append(rotated_molecule[m] + a*i + b*j + c*k)
					for m, atom in enumerate(self.others):
						species.append(atom[0])
						coords.append(others_array[m] + a*i + b*j + c*k)
		return Structure(lattice, species, coords, coords_are_cartesian=True)
						
	def write_vasp_files(self):
		struc = self.generate_structure()
		mpset = MPRelaxSet(struc, user_incar_settings={'GGA': 'PS', 'EDIFF': 1e-5, 'EDIFFG': -0.01, 'ALGO': 'F', 'ISMEAR': 0}, \
            user_kpoints_settings={'reciprocal_density': 100}, force_gamma=True)
		mpset.write_input(self.wd)


if __name__ == '__main__':
	struc = Poscar.from_file("/Users/yao/Google Drive/data/113/MAPbI3/phonon/MAPbI3_walsh.vasp").structure
	for j in range(10):
		fw = RandomOrientFileWriter(struc, [i for i in range(8)], 0, [2, 2, 2], '/Users/yao/Google Drive/data/113/MAPbI3/phonon/rotation'+str(j))
		fw.write_vasp_files()

	
