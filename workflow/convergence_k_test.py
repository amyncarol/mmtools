"""
given a folder of INCAR, POSCAR, POTCAR 
generate the several folders for different k mesh densities
"""
from pymatgen.io.vasp.inputs import Poscar, Kpoints
import os
from shutil import copyfile, move, rmtree
from glob import glob

def convergence_k_test(folder, start_k_density, ):
	"""
	Args:
		folder: the working folder under which has the input files
		start_k_density: the start k density we use to generate KPOINTS files

	Writes:
		several folders with same vasp inputs except for the KPOINTS
	"""
	#clean the folder
	for path in glob(folder+'/kpoints-*'):
		rmtree(path)

	#if CONTCAR is input instead of POSCAR
	if os.path.exists(folder + '/CONTCAR'):
		move(folder + '/CONTCAR', folder + '/POSCAR')
		
	#generate 5 k-points file and create folders
	struc = Poscar.from_file(folder + '/POSCAR').structure
	k_set = set()
	i = 0
	while len(k_set) < 5:
		kppa = start_k_density + i*150
		i += 1
		kpoints = Kpoints.automatic_gamma_density(struc, kppa)
		if kpoints.kpts[0][0] in k_set:
			continue
		else:
			k_set.add(kpoints.kpts[0][0])			
			subfolder = folder+'/kpoints-'+str(kpoints.kpts[0][0])
			os.mkdir(subfolder)
			print(subfolder)
			copyfile(folder+'/INCAR', subfolder+'/INCAR')
			copyfile(folder+'/POSCAR', subfolder+'/POSCAR')
			copyfile(folder+'/POTCAR', subfolder+'/POTCAR')
			kpoints.write_file(subfolder+'/KPOINTS')

if __name__=='__main__':
	for folder in glob('/Users/yao/Google Drive/data/2116/data/convergence/*'):
		print(folder)
		convergence_k_test(folder, 10)
	