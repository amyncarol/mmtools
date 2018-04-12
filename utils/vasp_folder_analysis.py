import os
from multiprocessing import Pool

from pymatgen import MPRester
from pymatgen.entries.compatibility import MaterialsProjectCompatibility
from pymatgen.analysis.phase_diagram import PhaseDiagram
from pymatgen.io.vasp.outputs import Vasprun

class VaspFolderAnalyzerSerial(object):
	"""
	Given a folder containing subfolders of vasp results, analyze all subfolders
	"""
	def __init__(self, folder):
		"""
		Args:
			folder: the parent folder that contains all subfolders

		Attributes:
			subfolders: a list of subfolders
			vasprun_paths: a list of absolute paths to the vasprun files
		"""
		self.folder = folder
		self.subfolders = [i for i in os.listdir(self.folder) if os.path.isdir(os.path.join(self.folder, i))]
		self.vasprun_paths = [os.path.join(os.path.join(folder, i), 'vasprun.xml') for i in self.subfolders]

	def write_formation_energyies(self, filename = 'formation_energy'):
		"""		
		Writes: 
			each line is of the following format:
				folder_name formation_energy
		"""
		with open(os.path.join(self.folder, filename), 'w') as f:
			for (subfolder, vasprun_file) in zip(self.subfolders, self.vasprun_paths):
				try: 
					vasprun = Vasprun(vasprun_file)
					energy = VasprunAnalyzer(vasprun).get_formation_energy()
					f.write(subfolder+' '+str(energy)+'\n')
				except:
					print('something wrong with {}'.format(subfolder))

class VasprunAnalyzer(object):
	"""
	analyze various materials properits from a vasprun.xml file
	"""
	def __init__(self, vasprun):
		"""
		Args:

		vasprun: a pymatgen vasprun object
		"""
		self.vasprun = vasprun

	def get_pd(self):
		"""
		get the phase diagram object

		Returns:
			pd: the phase diagram object
			entry: entry contained in vasprun
		"""
		entry = self.vasprun.get_computed_entry()
		compat = MaterialsProjectCompatibility()
		entry = compat.process_entry(entry)

		el = [specie.symbol for specie in entry.composition.keys()]
		with MPRester(api_key="64JmsIV32c8lUaxu") as mpr:
			entries = mpr.get_entries_in_chemsys(el)
		entries.append(entry)
		pd = PhaseDiagram(entries)
		return pd, entry

	def get_formation_energy(self):
		"""
		Returns:
			the formation energy per atom in eV/atom
		"""
		pd, entry = self.get_pd()
		return pd.get_form_energy_per_atom(entry)

	def get_energy_above_hull(self):
		"""
		phases in the phase diagram are retrieved from materials projects

		Returns:
			decomp_reduce: the lists of decomposition phases
			hull: the energy above hull per atom in eV/atom

		"""
		pd, entry = self.get_pd()
		(decomp, hull) = pd.get_decomp_and_e_above_hull(entry)
		decomp_reduced = [compound.composition.reduced_formula for compound in decomp]
		return decomp_reduced, hull

	def get_absolute_energy(self):
		"""
		Returns:
			the absolute energy per atom in eV/atom
		"""
		pass

	def get_band_gap(self):
		"""
		Returns:
			the band gap in eV
		"""
		pass

	def get_effective_mass(self, outcar):
		"""
		Arg: 
			outcar: a pymatgen outcar object

		Returns:
			the effective mass
		"""
		pass

if __name__ == '__main__':
	folder_ana = VaspFolderAnalyzerSerial('/Users/yao/Google Drive/data/2116/solidsolution/solid_solution_ml/complete')
	folder_ana.write_formation_energyies()
