import os
from multiprocessing import Pool

from pymatgen import MPRester
from pymatgen.entries.compatibility import MaterialsProjectCompatibility
from pymatgen.analysis.phase_diagram import PhaseDiagram, GrandPotentialPhaseDiagram
from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.core.periodic_table import Element

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
            entries: all entries used to construct the phase diagram
        """
        entry = self.vasprun.get_computed_entry()
        compat = MaterialsProjectCompatibility()
        entry = compat.process_entry(entry)

        el = [specie.symbol for specie in entry.composition.keys()]
        with MPRester(api_key="64JmsIV32c8lUaxu") as mpr:
            entries = mpr.get_entries_in_chemsys(el)
        entries.append(entry)
        pd = PhaseDiagram(entries)
        return pd, entry, entries

    def get_formation_energy(self):
        """
        Returns:
            the formation energy per atom in eV/atom
        """
        pd, entry, _ = self.get_pd()
        return pd.get_form_energy_per_atom(entry)

    def get_energy_above_hull(self):
        """
        phases in the phase diagram are retrieved from materials projects

        Returns:
            decomp_reduce: the lists of decomposition phases
            hull: the energy above hull per atom in eV/atom

        """
        pd, entry, _ = self.get_pd()
        (decomp, hull) = pd.get_decomp_and_e_above_hull(entry)
        decomp_reduced = [compound.composition.reduced_formula for compound in decomp]
        return decomp_reduced, hull

    def get_pd_with_open_element(self, open_el):
        """
        Return a list of grand canonical phase diagrams with one open element

        Args:
            open_el: the open element, a pymaten Element object

        Returns:
            a list of grand canonical phase diagrams, 
            a list of corresponding min_chempots, 
            a list of corresponding avg_chempots, 
            a list of corresponding max_chempots, 
        """
        pd, entry, entries = self.get_pd()
        chempots = pd.get_transition_chempots(open_el)

        gcpds = []
        min_chempots = []
        max_chempots = []
        avg_chempots = []

        for i in range(len(chempots)):
            if i == len(chempots)-1:
                avg_chempot = chempots[i] - 0.1
                min_chempot = None
            else:
                avg_chempot = 0.5*(chempots[i]+chempots[i+1])
                min_chempot = chempots[i+1]

            gcpds.append(GrandPotentialPhaseDiagram(entries, {open_el: avg_chempot}, pd.elements))
            min_chempots.append(min_chempot)
            max_chempots.append(chempots[i])
            avg_chempots.append(avg_chempot)
            
        return gcpds, min_chempots, avg_chempots, max_chempots

    def get_stable_and_unstable_with_open_element(self, open_el, filename):
        """
        given open element, print all stable and unstable phases under different chemical potentials for the open element

        Args:
            open_el: the open element, a pymaten Element object
            filename: the path to the file that stores the information

        Writes:
            a file contains all information
        """

        pds, min_chempots, avg_chempots, max_chempots = self.get_pd_with_open_element(open_el)

        with open(filename, 'w') as f:
            for pd, min_chempot, avg_chempot, max_chempot in zip(pds, min_chempots, avg_chempots, max_chempots):

                f.write("-------Chempot range: {} to {}, calculated at {}--------\n".format(min_chempot, max_chempot, avg_chempot))

                f.write('--------Stable Entries (formula, materials_id)--------\n')
                for e in pd.stable_entries:
                    f.write("{}, {}\n".format(e.original_comp.reduced_formula, e.entry_id))

                f.write('---------Unstable Entries (formula, materials_id, e_above_hull (eV/atom), decomposes_to)--------\n')
                for e in pd.unstable_entries:
                    decomp, e_above_hull = pd.get_decomp_and_e_above_hull(e)
                    pretty_decomp = [("{}:{}".format(k.original_comp.reduced_formula, k.entry_id), round(v, 2)) for k, v in decomp.items()]
                    f.write("{}, {}, {}, {}\n".format(e.original_comp.reduced_formula, e.entry_id, "%.3f" % e_above_hull, pretty_decomp))
                f.write('\n')
        
    def get_chempot_range(self, open_el):
        """
        Not working

        given a open element, get the chemical potential range for this entry

        Args:
            open_el: the open element, a pymaten Element object
        Returns:
            a list of simplexes
        """

        pd, entry, _ = self.get_pd()
        entry_dict = pd.get_chempot_range_map([open_el], referenced=True)
        return entry_dict[entry]


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
    #folder_ana = VaspFolderAnalyzerSerial('/Users/yao/Google Drive/mmtools/data/sample_vasp_calculation')
    #folder_ana = VaspFolderAnalyzerSerial('/Users/yao/Google Drive/data/2116/solidsolution/solid_solution_ml/complete')
    #folder_ana.write_formation_energyies()
    pass
