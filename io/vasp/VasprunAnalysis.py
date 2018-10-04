from pymatgen import MPRester
from pymatgen.entries.compatibility import MaterialsProjectCompatibility
from pymatgen.phasediagram.analyzer import PDAnalyzer
from pymatgen.phasediagram.maker import PhaseDiagram
from pymatgen.io.vasp.outputs import Vasprun

API_KEY = "64JmsIV32c8lUaxu"
class VasprunAnalysis():
    """
    can perform all custermized analysis of a single vasprun.xml file

    Includes:
        phase diagram analysis, energy above hull, formation energy.....

        will add more functionality in the future
    """
    def __init__(self, vasprun_file):
        self.vasprun = Vasprun(vasprun_file)
        self.pd, self.entry = self.get_pd()

    def get_pd(self):
        """
        get the phase diagram object for this compound

        Returns:
            phase diagram, entry
        """
        #make MP compatible entry from vasprun
        entry = self.vasprun.get_computed_entry()
        compat = MaterialsProjectCompatibility()
        entry = compat.process_entry(entry)

        el = [specie.symbol for specie in entry.composition.keys()]
        with MPRester(api_key=API_KEY) as mpr:
            entries = mpr.get_entries_in_chemsys(el)
        entries.append(entry)
        pd = PhaseDiagram(entries)
        return pd, entry

    def get_e_above_hull(self):
        """ 
        Get e_above_hull for this compound

        Args: 
            allow_negative: whether to calculate negative energy above hull for stable compound

        Returns:
            decomposition, energy above hull
        """
        pda = PDAnalyzer(self.pd)
        (decomp, hull) = pda.get_decomp_and_e_above_hull(self.entry)
        decomp = [compound.composition.reduced_formula for compound in decomp]
        return (decomp, hull)

    def get_formation_energy(self):
        """
        get formation energy per atom for this compound

        Returns: 
            formation energy per atom
        """
        return self.pd.get_form_energy_per_atom(self.entry)


    def get_equilibrium_reaction_energy(self):
        """
        Adapted from pymatgen. Only work if entry is stable(hull = 0)

        Provides the reaction energy of a stable entry from the neighboring
        equilibrium stable entries (also known as the inverse distance to
        hull).

        Returns:
            Equilibrium reaction energy of entry. Stable entries should have
            equilibrium reaction energy <= 0.
        """
        if self.entry not in self.pd.stable_entries:
            raise ValueError("Equilibrium reaction energy is available only "
                             "for stable entries.")
        entries = [e for e in self.pd.stable_entries if e != self.entry] #all stable entries without this stable entry
        modpd = PhaseDiagram(entries, self.pd.elements)
        analyzer = PDAnalyzer(modpd)
        (decomp, hull) = analyzer.get_decomp_and_e_above_hull(self.entry,
                                                    allow_negative=True)
        decomp = [compound.composition.reduced_formula for compound in decomp]
        return (decomp, hull)


if __name__=='__main__':
    vasprun_file = '/Users/yao/Google Drive/data/2116/solidsolution/Cs2Ag1In1Cl6/vasprun.xml'
    va = VasprunAnalysis(vasprun_file)
    print(va.get_equilibrium_reaction_energy())


