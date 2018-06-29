from pymatgen.io.vasp.inputs import Poscar
import os

def poscar2cartesian():
	"""
	read POSCAR in the current directory and 
	output POSCAR_cartesian in the same directory
	"""
	dir_path = os.path.dirname(os.path.realpath(__file__))
	poscar = Poscar.from_file(os.path.join(dir_path, 'POSCAR'))
	poscar_string = poscar.get_string(direct=False)
	with open(os.path.join(dir_path, 'POSCAR_cartesian'), 'w') as f:
		f.write(poscar_string)


if __name__=='__main__':
	poscar2cartesian()

