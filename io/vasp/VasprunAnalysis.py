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

        Band gap....

        will add more functionality in the future
    """
    def __init__(self, vasprun_file, cal_pd = False):
        self.vasprun = Vasprun(vasprun_file)
        if cal_pd:
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

    def get_eg(self):
        """
        return the band gap(by me), band gap(by pymatgen) and direct band gap

        two ways of calculating band gap to test the correctness of pymatgen algorithm

        Returns:
            eg1: band gap calculated by me
            eg2: band gap calculated by pymatgen algorithm
            eg_direct: direct band gap
        """
        eg1 = self.vasprun.eigenvalue_band_properties[0]
        efermi = self.vasprun.efermi
        vbm = self.vasprun.eigenvalue_band_properties[2]

        eig = [self.vasprun.eigenvalues[i] for i in self.vasprun.eigenvalues]
        occu = [i.tolist() for i in eig]
        occu_flat = [i[1] for sublist in occu for subsublist in sublist for i in subsublist]
        occu_set = set(occu_flat)
        
        if efermi < vbm:
            eg1 = 0
        elif len(occu_set) > 2:
            eg1 = 0

        bs = self.vasprun.get_band_structure()
        eg2 = bs.get_band_gap()['energy']
        eg_direct = bs.get_direct_band_gap()

        print('eg1 = {}, eg2 = {}, eg_direct = {}'.format(eg1, eg2, eg_direct))
        print('check if eg1 = eg2, it should')
        print('transition: {}'.format(bs.get_band_gap()['transition']))

        return eg1, eg2, eg_direct


if __name__=='__main__':
    # vasprun_file = '/Users/yao/Google Drive/data/2116/solidsolution/Cs2Ag1In1Cl6/vasprun.xml'
    # va = VasprunAnalysis(vasprun_file, cal_pd = True)
    # print(va.get_equilibrium_reaction_energy())

    #test get_eg
    # vasprun_file = '/Users/yao/Google Drive/data/2116/HSE_SOC/Cs2In1In1Cl6/vasprun.xml'
    # va = VasprunAnalysis(vasprun_file)
    # print(va.get_eg())

    vasprun_file = '/Users/yao/Google Drive/data/2116/InTl/tetragonal_exp/HSEsol_relax_lattice/vasprun.xml'
    va = VasprunAnalysis(vasprun_file)
    print(va.get_eg())

