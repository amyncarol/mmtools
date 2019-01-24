from pymatgen.symmetry.bandstructure import HighSymmKpath 
from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.io.vasp.inputs import Kpoints, Poscar

def generateBsKpoints(structure, file, divisions = 15, path = None):
	ibz = HighSymmKpath(structure, symprec=0.01, angle_tolerance=5)
	print(ibz.kpath['kpoints'])

	if not path:
		kpoints = []
		labels = []
		for path in ibz.kpath['path']:
			k, l = helper(path, ibz)
			kpoints += k
			labels += l
	else:
		kpoints, labels = helper(path, ibz)

	kp = Kpoints("Line_mode KPOINTS file",
                       style=Kpoints.supported_modes.Line_mode,
                       coord_type="Reciprocal",
                       kpts=kpoints,
                       labels=labels,
                       num_kpts=int(divisions))
	kp.write_file(file)

def helper(path, ibz):
	kpoints = []
	labels = []

	kpoints.append(ibz.kpath["kpoints"][path[0]])
	labels.append(path[0])
	for i in range(1, len(path) - 1):
	    kpoints.append(ibz.kpath["kpoints"][path[i]])
	    labels.append(path[i])
	    kpoints.append(ibz.kpath["kpoints"][path[i]])
	    labels.append(path[i])
	kpoints.append(ibz.kpath["kpoints"][path[-1]])
	labels.append(path[-1])
	return kpoints, labels

if __name__=='__main__':
	structure = Poscar.from_file('/Users/yao/Google Drive/data/2116/InTl/tetragonal_exp/10atom_strain_PBE/POSCAR').structure
	file = '/Users/yao/Google Drive/data/2116/InTl/tetragonal_exp/KPOINTS'
	path = [u'X', u'\\Gamma', u'W', u'\\Gamma', u'L']
	generateBsKpoints(structure, file, divisions = 15, path = path)

