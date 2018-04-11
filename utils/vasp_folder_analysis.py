import os
from multiprocessing import Pool

from pymatgen import MPRester
from pymatgen.entries.compatibility import MaterialsProjectCompatibility
from pymatgen.phasediagram.maker import PhaseDiagram
from pymatgen.phasediagram.analyzer import PDAnalyzer
from pymatgen.io.vasp.outputs import Vasprun

class VaspFolderAnalyzer(object):
	"""
	Given a folder containing subfolders of vasp results, use multithreading to analyze all subfolders
	"""
	def __init__(self, folder, n_threads=4):
		"""
		Args:
			folder: the parent folder that contains all subfolders
			n_threads: the number of threads to use

		Attributes:
			subfolders: a list of subfolders
			vasprun_paths: a list of absolute paths to the vasprun files
		"""
		self.n_threads = n_threads
		self.folder = folder
		self.subfolders = [i for i in os.listdir(self.folder) if os.path.isdir(os.path.join(self.folder, i))]
		self.vasprun_paths = [os.path.join(os.path.join(folder, i), 'vasprun.xml') for i in self.subfolders]

	def get_vaspruns(self):
		"""
		get the list of vasprun objects

		Returns: 
			a list of vasprun objects
		"""
		p = Pool(self.n_threads)
		return p.map(Vasprun, self.vasprun_paths)

	@staticmethod
	def get_formation_energy(vasprun):
		"""
		wrapper

		Args:
			vasprun: a vasprun object

		Returns:
			formation energy per atom
		"""
		va = VasprunAnalyzer(vasprun)
		return va.get_formation_energy()

	def get_formation_energyies(self):
		"""
		get a list of lines which are to be writen
		
		Returns: 
			a list of lines, each line is of the following format:
				folder_name formation_energy
		"""
		vasprun_list = self.get_vaspruns()
		p = Pool(self.n_threads)
		formation_energy_list = p.map(self.get_formation_energy, vasprun_list)
		return formation_energy_list
	
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
		pda = PDAnalyzer(pd)
		(decomp, hull) = pda.get_decomp_and_e_above_hull(entry)
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
	folder_ana = VaspFolderAnalyzer('/work/05018/tg843171/stampede2/Ab/2116/solid_solution_ml/complete/', 8)
	print(folder_ana.get_formation_energyies())
