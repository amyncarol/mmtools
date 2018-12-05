from mmtools.io.vasp.phonon_dfpt import *

def main():
	"""
	an example usage case for phonon_dfpt
	"""
	outcar_file = '/Users/yao/Google Drive/data/113/MAPbI3/phonon/deutoration/dfpt-ps/OUTCAR'
	eigenvalues, eigenvectors = get_eigenvalue_eigenvector(outcar_file)
	supercell = Poscar.from_file('/Users/yao/Google Drive/data/113/MAPbI3/phonon/deutoration/dfpt-ps/POSCAR').structure

	# the phonon modes summary
	G = []
	M = []
	X = []
	R = []

	## index starts from 1
	G_index = []
	M_index = []
	X_index = []
	R_index = []

	for i in range(40):
		q_point = which_quasi_q_point(supercell, eigenvectors[-i-1], [2, 2, 2])
		q_label = convert_to_label(q_point)
		if q_label == 'G':
			G.append(eigenvalues[-i-1])
			G_index.append(len(eigenvectors)-i)
		elif q_label == 'M':
			M.append(eigenvalues[-i-1])
			M_index.append(len(eigenvectors)-i)
		elif q_label == 'X':
			X.append(eigenvalues[-i-1])
			X_index.append(len(eigenvectors)-i)
		elif q_label == 'R':
			R.append(eigenvalues[-i-1])
			R_index.append(len(eigenvectors)-i)

	print('G:')
	print(G)
	#print(G_index)
	print('M:')
	print(M)
	#print(M_index)
	print('X:')
	print(X)
	#print(X_index)
	print('R:')
	print(R)
	#print(R_index)

def compare_calculations():
	"""
	here I want to compare the modes of two calculations, the hypothesis is that there is not much 
	difference in terms of phonon mode. And I want to make sure that the vibrations for each mode for
	the two calculations is nearly the same, by only comparing the largest vibration in each mode. 
	"""
	outcar_file = '/Users/yao/Google Drive/data/113/MAPbI3/phonon/1st_try/dfpt-ps/OUTCAR'
	eigenvalues1, eigenvectors1 = get_eigenvalue_eigenvector(outcar_file)
	outcar_file = '/Users/yao/Google Drive/data/113/MAPbI3/phonon/deutoration/dfpt-ps/OUTCAR'
	eigenvalues2, eigenvectors2 = get_eigenvalue_eigenvector(outcar_file)
	for i in range(len(eigenvectors1)-1, len(eigenvectors1)-41, -1):
		# atom1, _,  amplitude1 = largest_vibration(eigenvectors1[i])
		# atom2, _,  amplitude2 = largest_vibration(eigenvectors2[i])

		print('the {}th mode'.format(i+1))
		print(eigenvalues1[i], largest_vibration(eigenvectors1[i]))
		print(eigenvalues2[i], largest_vibration(eigenvectors2[i]))

def compare_calculations2():
	"""
	here I want to test if v1 is not orthogonal to v2, if they not orthogonal, then v1 == v2
	"""
	outcar_file = '/Users/yao/Google Drive/data/113/MAPbI3/phonon/1st_try/dfpt-ps/OUTCAR'
	eigenvalues1, eigenvectors1 = get_eigenvalue_eigenvector(outcar_file)
	supercell1 = Poscar.from_file('/Users/yao/Google Drive/data/113/MAPbI3/phonon/1st_try/dfpt-ps/POSCAR').structure

	outcar_file = '/Users/yao/Google Drive/data/113/MAPbI3/phonon/deutoration/dfpt-ps/OUTCAR'
	eigenvalues2, eigenvectors2 = get_eigenvalue_eigenvector(outcar_file)
	supercell2 = Poscar.from_file('/Users/yao/Google Drive/data/113/MAPbI3/phonon/deutoration/dfpt-ps/POSCAR').structure

	miss_matched_i = []
	for i in range(len(eigenvalues1)-40, len(eigenvalues1)):
		print('-----------')
		v1 = eigenvectors1[i]
		for j in range(len(eigenvalues2)):
			v2 = eigenvectors2[j]
			if orthogonal(eigenvectors1[i], eigenvectors2[j], supercell1) != 0:
				if i==j:
					break
		else:
			for j in range(len(eigenvalues2)):
				v2 = eigenvectors2[j]
				if orthogonal(eigenvectors1[i], eigenvectors2[j], supercell1) != 0:
					q1 = which_quasi_q_point(supercell1, v1, [2, 2, 2])
					q2 = which_quasi_q_point(supercell2, v2, [2, 2, 2])
					print(i, j, convert_to_label(q1), convert_to_label(q2), eigenvalues1[i]-eigenvalues2[j])


if __name__ == '__main__':
    #main()
    compare_calculations2()