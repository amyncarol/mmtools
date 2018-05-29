import re
import numpy as np
from numpy.linalg import norm
from pymatgen.io.vasp.inputs import Poscar
import copy
def get_eigenvalue_eigenvector(outcar_file):
	"""
	return a list of eigenvalues and eigenvector of vasp dfpt phonon calculation
	"""
	eigenvalues = []
	with open(outcar_file, 'r') as f:
		lines = f.readlines()
		for l in lines:
			match = re.search('[\s]+f[\s]+=[\s]+', l)
			if match:
				freq = float(l.strip().split(' ')[-2])
				eigenvalues.append(freq)
		for l in lines:
			match = re.search('f/i', l)
			if match:
				freq = -float(l.strip().split(' ')[-2])
				eigenvalues.append(freq)

		n_modes = len(eigenvalues)
		n_atoms = int(n_modes/3)

		eigenvectors = []

		for i, l in enumerate(lines):
			match1 = re.search('[\s]+f[\s]+=[\s]+', l)
			match2 = re.search('f/i', l)
			if match1 or match2:
				eigenvector = []
				data_block = lines[i+2: i+2+n_atoms]
				for item in data_block:
					alist = item.strip().split()
					eigenvector.append([float(alist[-3]), float(alist[-2]), float(alist[-1])])
				eigenvector = np.array(eigenvector)
				eigenvectors.append(eigenvector)
	return eigenvalues, eigenvectors

def which_q_point(supercell, eigenvector, scaling_matrix):
	"""
	given a eigenvector and the supercell, judge which q point is this eigenvector

	scaling_matrix: i.e.[2, 2, 2]

	Warning: this is not perfect, this is a special case for [2, 2, 2]
	"""
	n = scaling_matrix[0]*scaling_matrix[1]*scaling_matrix[2] ##the number of unit cells
	atom_set = supercell.sites[:n] ##assume the first n atoms are in different unit cells
	relative_position = []
	for atom in atom_set:
		relative_position.append(atom.frac_coords-atom_set[0].frac_coords)

	q_point = np.zeros(3)
	for i in range(n):
		if (np.abs(eigenvector[i]+eigenvector[0])<1e-4).all() and abs(np.sum(relative_position[i])-0.5)<1e-4:
			q_point += relative_position[i]
	return q_point

def convert_to_label(q_point):
	if abs(np.sum(q_point)-0)<1e-4:
		return 'G'
	elif abs(np.sum(q_point)-0.5)<1e-4:
		return 'X'
	elif abs(np.sum(q_point)-1)<1e-4:
		return 'M'
	elif abs(np.sum(q_point)-1.5)<1e-4:
		return 'R'

def is_accoustic(supercell, eigenvector, scaling_matrix):
	"""
	given an eigenvector, judge it is optic or acoustic

	Not sure how to judge this.....
	"""
	n = scaling_matrix[0]*scaling_matrix[1]*scaling_matrix[2] ##the number of unit cells
	n_atoms = eigenvector.shape[0]
	weights = []
	for i, site in enumerate(supercell.sites):
		if i%n == 0:
			weights.append(site.specie.atomic_mass)
	weights = np.array(weights)

	displacement_sum = weights @ eigenvector[0:n_atoms:n]
	#displacement_sum = np.sum(eigenvector[0:n_atoms:n], axis=0)
	print(displacement_sum)
	# if not (np.abs(displacement_sum)<1e-4).all():
	# 	print(displacement_sum)
	# 	return True
	# else:
	# 	return False

def generate_distortion(supercell, eigenvector):
	"""
	given a supercell and phonon eigenvector, generate the distorted structures
	"""
	new_cell = copy.deepcopy(supercell)
	for i in range(len(new_cell.sites)):
		new_cell.translate_sites(i, eigenvector[i], frac_coords=False)

	weights = []
	for site in supercell.sites:
		weights.append(site.specie.atomic_mass)
	weights = np.array(weights)
	normal_mode_coords = np.sqrt(weights.T) @ norm(eigenvector, axis=1)
	print('the normal mode coords is {}, not sure if this is right though...'.format(normal_mode_coords))
	return new_cell


	
if __name__ == '__main__':
	outcar_file = '/Users/yao/Google Drive/data/113/MAPbI3/phonon/1st_try/dfpt-ps/OUTCAR'
	eigenvalues, eigenvectors = get_eigenvalue_eigenvector(outcar_file)
	supercell = Poscar.from_file('/Users/yao/Google Drive/data/113/MAPbI3/phonon/1st_try/dfpt-ps/POSCAR').structure

	## the phonon modes summary
	# G = []
	# M = []
	# X = []
	# R = []
	# for i in range(100):
	# 	q_point = which_q_point(supercell, eigenvectors[-i-1], [2, 2, 2])
	# 	q_label = convert_to_label(q_point)
	# 	if q_label == 'G':
	# 		G.append(eigenvalues[-i-1])
	# 	elif q_label == 'M':
	# 		M.append(eigenvalues[-i-1])
	# 	elif q_label == 'X':
	# 		X.append(eigenvalues[-i-1])
	# 	elif q_label == 'R':
	# 		R.append(eigenvalues[-i-1])	
		
	# print('G:')	
	# print(G)
	# print('M:')	
	# print(M)
	# print('X:')	
	# print(X)
	# print('R:')	
	# print(R)

	distorted_cell = generate_distortion(supercell, eigenvectors[-1])
	Poscar(distorted_cell).write_file('/Users/yao/Google Drive/data/113/MAPbI3/phonon/1st_try/dfpt-ps-distorted-RTA1/POSCAR')


		
		
	




	