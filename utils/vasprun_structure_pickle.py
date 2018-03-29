from pymatgen.io.vasp.outputs import Vasprun
import os
import pickle

def store_structures(folder, filename):
	"""
	read in vasprun files, convert to structure files and pickle them

	Args:
		folder: under with the subfolder with vasprun.xml files in it
				i.e. folder/subfolder/vasprun.xml
		filename: filename to store the pickle(in folder)
	"""
	structure_list = []
	for subfolder in os.listdir(folder):
		vasprun_file = os.path.join(folder, subfolder+'/vasprun.xml')
		if os.path.exists(vasprun_file):
			struc = Vasprun(vasprun_file).final_structure
			structure_list.append(struc)

	with open(os.path.join(folder, filename), 'wb') as f:
		pickle.dump(structure_list, f, -1)


def load_structures(filepath):
	"""
	read the pickled structure_list in file as given by filepath

	Returns: 
		a list of structures
	"""
	with open(filepath, 'rb') as f:
		structure_list = pickle.load(f)

	#print(structure_list)

if __name__ == '__main__':
	store_structures('/Users/yao/Google Drive/mmtools/data/sample_vasp_calculation', 'structure.pkl')
	load_structures('/Users/yao/Google Drive/mmtools/data/sample_vasp_calculation/structure.pkl')