import unittest
import numpy as np
from mmtools.io.vasp.phonon_dfpt import *
from math import pi
from numpy.linalg import norm
from numpy import (array, dot, arccos)
from numpy.testing import assert_almost_equal, assert_equal
class TestPhononDfpt(unittest.TestCase):		
	def test_angle(self):
		a = np.array([0, 0, 1])
		b = np.array([0, 1, 0])
		self.assertEqual(angle_between_vector(a, b), pi/2)

		a = np.array([1, 1, 0])
		b = np.array([0, 1, 0])
		self.assertAlmostEqual(angle_between_vector(a, b), pi/4)

		a = np.array([1, 1, 0])
		b = np.array([0, -1, 0])
		self.assertAlmostEqual(angle_between_vector(a, b), 3*pi/4)

		v1 = np.array([0.018856,  0.000117, -0.010526])
		v2 = np.array([0.018856,  0.000117, -0.010526])
		self.assertAlmostEqual(angle_between_vector(v1, v2), 0)

	def test_which_q_point_perfect_supercell(self):
		"""
		test with the perfect cell
		"""
		print('test perfect cell')
		outcar_file = '/Users/yao/Google Drive/data/113/MAPbI3/phonon/1st_try/dfpt-ps/OUTCAR'
		eigenvalues, eigenvectors = get_eigenvalue_eigenvector(outcar_file)
		supercell = Poscar.from_file('/Users/yao/Google Drive/data/113/MAPbI3/phonon/1st_try/dfpt-ps/POSCAR').structure

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

		for i in range(20):
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

		assert_almost_equal(np.array(G), np.array([-0.081009, -0.070038, -0.057882, 1.790555]))
		assert_almost_equal(np.array(M), np.array([-1.891642, 0.443186, 1.838702, 2.082854, 2.585056]))
		assert_almost_equal(np.array(X), np.array([0.889369, 1.867747, 1.908073, 2.008787, 2.089795, 2.368314, 2.395276, 2.403152, 2.522396]))
		assert_almost_equal(np.array(R), np.array([-2.217779, -0.96519]))

		assert_equal(np.array(G_index), np.array([285, 284, 283, 280]))
		assert_equal(np.array(M_index), np.array([287, 282, 279, 275, 269]))
		assert_equal(np.array(X_index), np.array([281, 278, 277, 276, 274, 273, 272, 271, 270]))
		assert_equal(np.array(R_index), np.array([288, 286]))

	# def test_which_q_point_distorted_supercell(self):
	# 	"""
	# 	test with the distorted cell, the problem in this case is not clearly defined
	# 	"""
	# 	print('test distorted cell')
	# 	outcar_file = '/Users/yao/Google Drive/data/113/MAPbI3/phonon/1st_try/RTA1-dfpt-ps/smear0_07_expansion1_4/OUTCAR'
	# 	eigenvalues, eigenvectors = get_eigenvalue_eigenvector(outcar_file)
	# 	supercell = Poscar.from_file('/Users/yao/Google Drive/data/113/MAPbI3/phonon/1st_try/dfpt-ps/POSCAR').structure

	# 	# the phonon modes summary
	# 	G = []
	# 	M = []
	# 	X = []
	# 	R = []

	# 	## index starts from 1
	# 	G_index = []
	# 	M_index = []
	# 	X_index = []
	# 	R_index = []

	# 	for i in range(100):
	# 		q_point = which_quasi_q_point(supercell, eigenvectors[-i-1], [2, 2, 2])
	# 		q_label = convert_to_label(q_point)
	# 		if q_label == 'G':
	# 			G.append(eigenvalues[-i-1])
	# 			G_index.append(len(eigenvectors)-i)
	# 		elif q_label == 'M':
	# 			M.append(eigenvalues[-i-1])
	# 			M_index.append(len(eigenvectors)-i)
	# 		elif q_label == 'X':
	# 			X.append(eigenvalues[-i-1])
	# 			X_index.append(len(eigenvectors)-i)
	# 		elif q_label == 'R':
	# 			R.append(eigenvalues[-i-1])
	# 			R_index.append(len(eigenvectors)-i)

	# 	print('G:')
	# 	print(G)
	# 	#print(G_index)
	# 	print('M:')
	# 	print(M)
	# 	#print(M_index)
	# 	print('X:')
	# 	print(X)
	# 	#print(X_index)
	# 	print('R:')
	# 	print(R)
	# 	#print(R_index)

	#this is not right!!!
	# def test_is_transverse(self):
	# 	outcar_file = '/Users/yao/Google Drive/data/113/MAPbI3/phonon/1st_try/dfpt-ps/OUTCAR'
	# 	eigenvalues, eigenvectors = get_eigenvalue_eigenvector(outcar_file)
	# 	supercell = Poscar.from_file('/Users/yao/Google Drive/data/113/MAPbI3/phonon/1st_try/dfpt-ps/POSCAR').structure

	# 	# the phonon modes summary
	# 	G = []
	# 	M = []
	# 	X = []
	# 	R = []

	# 	for i in range(20):
	# 		q_point = which_quasi_q_point(supercell, eigenvectors[-i-1], [2, 2, 2])
	# 		q_label = convert_to_label(q_point)
	# 		trans = is_transverse(supercell, eigenvectors[-i-1], q_point)
	# 		if q_label == 'G':
	# 			G.append(trans)
	# 		elif q_label == 'M':
	# 			M.append(trans)
	# 		elif q_label == 'X':
	# 			X.append(trans)
	# 		elif q_label == 'R':
	# 			R.append(trans)

	# 	print('G:')
	# 	print(G)

	# 	print('M:')
	# 	print(M)

	# 	print('X:')
	# 	print(X)

	# 	print('R:')
	# 	print(R)




if __name__ == '__main__':
    unittest.main()