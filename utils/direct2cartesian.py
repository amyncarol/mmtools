import sys
from pymatgen.io.vasp.inputs import Poscar
def direct2cartesian(POSCAR_input, POSCAR_output):
	"""
	given a poscar file in direct coords, output a poscar file in cartesian coords
	"""
	poscar = Poscar.from_file(POSCAR_input)
	poscar.write_file(POSCAR_output, direct = False)

if __name__ == '__main__':
	print('usage: direct2cartesian.py POSCAR_input POSCAR_output')
	poscar_input = sys.argv[1]
	poscar_output = sys.argv[2]
	direct2cartesian(poscar_input, poscar_output)
