###density-density dielectric function, vasp5.4.4
import re
import os
import numpy as np
import pandas as pd
from scipy.integrate import simps, trapz
from math import exp, sqrt, pi, log
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.io import FortranFile
from numpy.linalg import norm
from pymatgen.io.vasp.outputs import Outcar, Vasprun
from pymatgen.electronic_structure.core import Spin

h = 6.62607004e-34
c = 299792458.0
e = 1.60217662e-19
kB = 1.38064852e-23
T = 300
sun_solid_angle = 6.807e-5
path = os.path.abspath(__file__)
dir_path = os.path.dirname(path)

class OpticalAbsorption:
    def __init__(self, filename, waveder_file_ip, vasprun_file_static, nedos, method = 'ip', scissor = 0):
        """
        filename: OUTCAR file
        method: can be 'ip' or 'rpa'
        """
        self.filename = filename
        self.waveder_file_ip = waveder_file_ip
        self.vasprun_file_static = vasprun_file_static
        self.nedos = nedos
        self.method = method
        self.sun = SunSpectrum()
        self.scissor = scissor

    def get_dielect(self):
        """
        LOPTICS parsing:
        returns energy and imaginery part of epsilon, real part of epsilon

        tested
        """
        if self.method == 'ip':
            return self.get_dielect_ip()
        elif self.method == 'rpa':
            return self.get_dielect_rpa()
        else:
            raise Exception('only support method = ip or rpa')


    def get_dielect_ip(self):
        """
        tested
        """
        i = -4
        data = np.zeros((self.nedos, 3))
        with open(self.filename) as f:
            for line in f:
                if 'IMAGINARY' in line.split() and 'density-density' in line.split():
                    i = -3
                if i in [-3, -2, -1]:
                    i += 1 
                    continue
                if i == self.nedos:
                    break
                if i >= 0:
                    data[i][0] = float(line.split()[0])
                    data[i][1] = float(line.split()[1])
                    i += 1

        i = -4
        with open(self.filename) as f:
            for line in f:
                if 'REAL' in line.split() and 'density-density' in line.split():
                    i = -3
                if i in [-3, -2, -1]:
                    i += 1 
                    continue
                if i == self.nedos:
                    break
                if i >= 0:
                    data[i][2] = float(line.split()[1])
                    i += 1             
        return np.transpose(data)

    def get_dielect_rpa(self):
        """
        tested
        """
        rpa_block = False
        data = np.zeros((self.nedos, 3))
        i = 0
        with open(self.filename) as f:
            for line in f:
                if 'DIELECTRIC' in line.split() and 'RPA' in line.split():
                    rpa_block = True
                if rpa_block and 'dielectric' in line.split() and line.strip().split()[3] == 'dielectric':
                    data[i][0] = float(line.split()[0])
                    data[i][1] = float(line.split()[2])
                    data[i][2] = float(line.split()[1])
                    i += 1
                if 'XI_TO_W:' in line.split():
                    rpa_block = False
                    print(rpa_block)
        return np.transpose(data)

    def get_absorb(self):
        """
        absorption coefficient unit: cm-1
        returns energy and absorption coefficient

        energy in eV

        tested
        """ 
        data = self.get_dielect()
        energy = data[0]+self.scissor
        imag_iso = data[1]
        real_iso = data[2]
        constant = 4*pi*e/(h*c)*0.01
        norm_iso = np.zeros(self.nedos)
        k_iso = np.zeros(self.nedos)
        absorb_iso = np.zeros(self.nedos)
        for i in range(self.nedos):
            norm_iso[i] = sqrt(imag_iso[i]**2+real_iso[i]**2)
            k_iso[i] = sqrt(2)/2*sqrt(norm_iso[i]-real_iso[i])
            absorb_iso[i] = energy[i]*k_iso[i]*constant

        ##we should make sure energy start from 0
        if self.scissor != 0:
            energy_old = np.copy(energy)
            absorb_iso_old = np.copy(absorb_iso)
            energy = np.zeros(self.nedos+100)
            absorb_iso = np.zeros(self.nedos+100)
            energy[:100] = np.linspace(0, self.scissor, 100)
            energy[100:] = energy_old
            absorb_iso[100:] = absorb_iso_old
        return energy, absorb_iso

    def get_absorb_onset(self, is_scissor, threshold = 1e4):
        """
        return the absorption onset where absorb_iso > 1e4

        Args:
            is_scissor: whether we use scissor operation here
        """
        energy, absorb_iso = self.get_absorb()
        if not is_scissor:
            energy = energy - self.scissor
        for i in range(len(energy)):
            if absorb_iso[i] >= threshold:
                return energy[i]

    def get_absorptivity(self, thickness):
        """
        thickness (unit: cm)
        absorptivity (no unit)
        return absorptivity a(E) given thickness and energy

        formula reference: Yu and Zunger, PRL 108, 068701 (2012)
        In SLME, L is the thickness of the thin film with a
        zero-reflectivity front surface and unity-reflectivity back surface.

        Returns:
            photon energies are in nm

        tested
        """
        photon_eV, absorb_iso = self.get_absorb()
        a_iso = 1-np.exp(-2*absorb_iso*thickness)
        photon_nm = self.sun.eV_to_nm(photon_eV)
        return photon_nm[1:], a_iso[1:]

    def get_sc_current(self, thickness):
        """
        calculalte the short circuit current density for this material with certain thickness
        unit: A/m^2

        formula reference: Yu and Zunger, PRL 108, 068701 (2012)

        Note: all photon energies are in nm here

        tested: the magnitude is correct
        """
        photon_sun, photon_flux = self.sun.get_sun_photon_flux() #energy range 0~4 eV
        photon_material, absorptivity = self.get_absorptivity(thickness) #energy range 0~60 eV
        f_absorb = interp1d(photon_material, absorptivity)
        absorptivity_new = f_absorb(photon_sun)
        #plt.plot(photon_sun, absorptivity_new*photon_flux)
        #plt.show()
        return e*simps(absorptivity_new*photon_flux, photon_sun)

    def get_fraction_radiative(self):
        if self.method != 'ip':
            raise Exception('Not implemented yet!')
        else: 
            mea = MatrixElementAnalysis(self.waveder_file_ip, self.filename, self.vasprun_file_static)
            E_da_g = mea.get_dipole_allow_eg()
            E_g = mea.get_band_gap()
            E_delta = E_da_g - E_g
            return exp(-E_delta*e/kB/T)

    def get_reverse_saturation_current_radiative(self, thickness):
        """
        calculalte the radiative component of reverse saturation current for this material with certain thickness

        formula reference: Yu and Zunger, PRL 108, 068701 (2012)

        unit: A/m^2

        tested: the magnitude is correct
        """
        photon_material, absorptivity = self.get_absorptivity(thickness) #energy range 0~60 eV
        f_absorb = interp1d(photon_material, absorptivity)
        photon_nm = np.linspace(min(photon_material), 4000, 1000)
        bb_photon_flux = self.sun.get_blackbody_photon_flux(photon_nm)
        absorptivity_new = f_absorb(photon_nm)

        # plt.plot(photon_nm, bb_photon_flux*absorptivity_new)
        # plt.show()
        return e*pi*simps(absorptivity_new*bb_photon_flux, photon_nm)

    def get_reverse_saturation_current(self, thickness):
        """
        unit: A/m^2

        tested: the magnitude is correct
        """
        return self.get_reverse_saturation_current_radiative(thickness)/self.get_fraction_radiative()

    def get_J(self, V, thickness):
        """
        Note: the formula in Yu and Zunger, PRL 108, 068701 (2012) is wrong

        Correct Ref: IEEE TRANSACTIONS ON ELECTRON DEVICES, VOL. ED-31, NO. 5 , MAY 1984
        equation (6)

        J: current density
        """
        J_sc = self.get_sc_current(thickness)
        J_0 = self.get_reverse_saturation_current(thickness)
        return J_sc - J_0 * np.exp(e*V/kB/T)

    def get_V_oc(self, thickness):
        """
        unit: Volt

        tested: the magnitude is correct
        """
        J_sc = self.get_sc_current(thickness)
        J_0 = self.get_reverse_saturation_current(thickness)
        return log(J_sc/J_0)*kB*T/e

    def get_max_power_density(self, thickness):
        """
        power density: power per area: unit: W/m^2

        tested: the magnitude is correct
        """
        V_oc = self.get_V_oc(thickness)
        V = np.linspace(0, V_oc, num=1000)
        J = self.get_J(V, thickness)
        plt.plot(V, J)
        plt.xlabel('Voltage(V)')
        plt.ylabel('Current Density(A/m^2)')
        plt.show()
        return np.max(V*J)

    def get_max_efficiency(self, thickness):
        """
        tested: the magnitude is correct
        """
        P_max = self.get_max_power_density(thickness)
        P_sun = self.sun.get_sun_power()
        return P_max/P_sun

class SunSpectrum:
    def __init__(self, filename = os.path.join(dir_path, "ASTMG173.csv")):
        """
        https://rredc.nrel.gov/solar//spectra/am1.5/ASTMG173/ASTMG173.html
        """
        a = pd.read_csv(filename)
        a = a.convert_objects(convert_numeric=True)
        #photon in nm
        self.photon_nm = a.iloc[1:,0].values
        ##the photon energy in J
        self.photon_J = SunSpectrum.nm_to_J(self.photon_nm)
        ##the spectral irradiance in W/(m^2*nm)
        self.etr = a.iloc[1:,1].values
        self.tilt = a.iloc[1:,2].values
        self.direct = a.iloc[1:,3].values

    @staticmethod
    def nm_to_eV(a):
        """
        a is a vector
        convert photon energy from nm to eV

        tested
        """
        return h*c/(a*1e-9)/e

    @staticmethod
    def eV_to_nm(a):
        """
        a is a vector
        convert photon energy from eV to nm
        """
        return h*c/(a*1e-9)/e

    @staticmethod
    def nm_to_J(wavelength):
        """
        a is a vector
        convert photon energy from nm to J

        tested
        """
        return h*c/(wavelength*1e-9)

    def get_sun_photon_flux(self, spectrum='tilt'):
        """
        photon flux(number of photons per second per unit area) for a photon energy interval
        unit: # of photons/(sec*m^2*nm)

        Note! Since sun spectrum is continuous, we can only have photon flux for an energy interval instead of a exact energy

        https://www.pveducation.org/pvcdrom/properties-of-sunlight/photon-flux
        https://www.pveducation.org/pvcdrom/properties-of-sunlight/spectral-irradiance
        https://www.pveducation.org/pvcdrom/properties-of-sunlight/radiant-power-density

        Args:
            spectrum: can be 'tilt' or 'direct'

        Returns:
            photon flux from sun
        """
        if spectrum == 'tilt':
            spectral_irradiance = self.tilt
        elif spectrum == 'direct':
            spectral_irradiance = self.direct
        else:
            raise Exception('spectrum: can be tilt or direct')

        ##the photon_flux in # of photons/(sec*m^2*nm)
        photon_flux = spectral_irradiance/self.photon_J
        return self.photon_nm, photon_flux

    def get_sun_irradiance(self):
        """
        spectral irradiance: It gives the power density at a particular wavelength. 
        The units of spectral irradiance are in W m-2 nm-1

        https://www.pveducation.org/pvcdrom/properties-of-sunlight/spectral-irradiance

        tested
        """        
        return self.photon_nm, self.etr, self.tilt, self.direct

    def get_sun_power(self, spectrum='tilt'):
        """
        sun power: unit: W m-2
        """
        if spectrum == 'tilt':
            spectral_irradiance = self.tilt
        elif spectrum == 'direct':
            spectral_irradiance = self.direct
        else:
            raise Exception('spectrum: can be tilt or direct')

        return trapz(spectral_irradiance, self.photon_nm)

    @staticmethod
    def get_blackbody_photon_flux(photon_nm): #photon in nm
        """
        blackbody photon flux: unit: # of photons/(sec*m^2*nm)

        Ref: www.spectralcalc.com, Calculating Blackbody Radiance
        Ref: IEEE TRANSACTIONS ON ELECTRON DEVICES, VOL. ED-31, NO. 5 , MAY 1984
        The refs are consistent with each other

        tested
        """
        photon_J = SunSpectrum.nm_to_J(photon_nm)
        return SunSpectrum.get_blackbody_spectral_irradiance(photon_nm)/photon_J

    @staticmethod
    def get_blackbody_spectral_irradiance(photon_nm): #photon in nm
        """
        blackbody photon spectral irradiance: unit: W/(m^2*nm)

        Ref: www.spectralcalc.com, Calculating Blackbody Radiance

        tested
        """
        photon_J = SunSpectrum.nm_to_J(photon_nm)
        return 2*1e36*h*c**2/photon_nm**5*np.exp(-photon_J/kB/T)

class Waveder:
    """
    Class for reading a WAVEDER file.
    The LOPTICS tag produces a WAVEDER file.
    The WAVEDER contains the derivative of the orbitals with respect to k.
    Author: Kamal Choudhary, NIST

    Args:
        filename: Name of file containing WAVEDER.
    """

    def __init__(self, filename):
        with FortranFile(filename, "r") as f:
            val = (f.read_reals(dtype=np.int32))
            nbands = int(val[0])
            nelect = int(val[1])
            nk = int(val[2])
            ispin = int(val[3])
            self.nodes_in_dielectric_function = (f.read_reals(dtype=np.float))
            self.wplasmon = (f.read_reals(dtype=np.float))

            self.cder = np.array((f.read_reals(dtype=np.float)))
            #self.cder_data = self.cder.reshape((nbands, nelect, nk, ispin, 3)) #pymatgen version
            self.cder_data = self.cder.reshape((nbands, nelect, nk, 6))
            self._nkpoints = nk
            self.ispin = ispin
            self._nelect = nelect
            self._nbands = nbands

    @property
    def nbands(self):
        """
        Returns the number of bands in the calculation
        """
        return self._nbands

    @property
    def nkpoints(self):
        """
        Returns the number of k-points in the calculation
        """
        return self._nkpoints

    @property
    def nelect(self):
        """
        Returns the number of electrons in the calculation
        """
        return self._nelect


class MatrixElementAnalysis:
    def __init__(self, waveder_file_ip, outcar_file_ip,  vasprun_file_static):
        """
        waveder and outcar files from the linear optical independent particle approximation
        vasprun file from the associated static calculation
        """
        waveder = Waveder(waveder_file_ip)
        outcar = Outcar(outcar_file_ip)
        vasprun = Vasprun(vasprun_file_static)

        M_norm = norm(waveder.cder_data, axis = -1)*(1e-10/e)
        self.M_squared = M_norm*M_norm

        nelect = round(outcar.nelect)
        self.vbm_band = nelect//2-1
        self.cbm_band = nelect//2

        self.bs = vasprun.get_band_structure()
        self.nbands = self.bs.nb_bands
        self.kpoints = self.bs.kpoints
        self.nkpoints = waveder.nkpoints

    def get_dipole_allow_eg(self):
        ##get the smallest direct dipole allowed band gap
        best_v_band = None
        best_c_band = None
        best_kpoint_index = None
        best_kpoint = None
        best_strength = None
        best_gap = 100

        for v in range(self.vbm_band, -1, -1):
            for c in range(self.cbm_band, self.nbands):
                for k in range(self.nkpoints):
                    gap = self.bs.bands[Spin.up][c][k]-self.bs.bands[Spin.up][v][k]
                    strength = self.M_squared[c][v][k] 
                    if strength >= 1e-3 and gap < best_gap:
                        best_v_band = v
                        best_c_band = c
                        best_kpoint_index = k
                        best_kpoint = self.kpoints[k].frac_coords
                        best_strength = strength
                        best_gap = gap
        print('The smallest direct dipole allowed transition:')
        print('from valence band {} to conduction band {}'.format(best_v_band, best_c_band))
        print('at kpoint {}:{}'.format(best_kpoint_index, best_kpoint))
        print('Has transition matrix squared {}'.format(best_strength))
        print('Eda_g (dipole allowed gap) is {}'.format(best_gap))
        print('\n')

        ##the information for the smallest direct band gap
        direct_band_gap = self.bs.get_direct_band_gap_dict()[Spin.up]
        v = direct_band_gap['band_indices'][0]
        c = direct_band_gap['band_indices'][1]
        k = direct_band_gap['kpoint_index']
        gap = direct_band_gap['value']
        print('Compare with direct band gap:')
        print('from valence band {} to conduction band {}'.format(v, c))
        print('at kpoint {}:{}'.format(k, self.kpoints[k].frac_coords))
        print('Has transition matrix squared {}'.format(self.M_squared[c][v][k]))
        print('Eda_g (dipole allowed gap) is {}'.format(gap))
        print('\n')

        return best_gap

    def get_band_gap(self):
        gap = self.bs.get_band_gap()['energy']
        print('The band gap is {}'.format(gap))
        print('\n')
        return gap


if __name__=='__main__':
    # plot blackbody spectrum and sun spectrum
    sun = SunSpectrum()
    photon_nm_sun, etr, tilt, direct = sun.get_sun_irradiance()
    photon_nm = np.linspace(0, 2500, 3000)
    irradiance = sun.get_blackbody_spectral_irradiance(photon_nm)*sun_solid_angle
    plt.plot(photon_nm_sun, etr, label='Etr')
    plt.plot(photon_nm_sun, tilt, label='Global Tilt')
    plt.plot(photon_nm_sun, direct, label='Direct+circumsolar')
    plt.plot(photon_nm, irradiance, label='5523K blackbody')
    plt.legend()
    plt.xlabel('Photon nm')
    plt.ylabel('Spectral irradiance W m-2 nm-1')
    plt.show()




