#This is to use the 216paper env and pymatgen4.5.7 with python3.5, I really hate hwo pymatgen updated!!!!!
#and matplotlib 2.0.1
from pymatgen.io.vasp.outputs import Vasprun
import matplotlib.pyplot as plt
import matplotlib
import os
import numpy as np
from matplotlib import gridspec
import re
from pymatgen.core.periodic_table import Element
from pymatgen.electronic_structure.core import Spin
import matplotlib.ticker as ticker
matplotlib.rcParams.update({'font.size': 20})

def get_k_dist(kpoints):
    dist = np.zeros(len(kpoints))
    for i in range(len(kpoints)-1):
        dist[i+1]=dist[i]+np.linalg.norm(kpoints[i]-kpoints[i+1])
    return dist
    
def get_special_k(kpoints):
    k_index = [0]
    for i in range(len(kpoints)-1):
        if np.linalg.norm(kpoints[i]-kpoints[i+1]) < 1e-8:
            k_index.append(i)
    k_index.append(len(kpoints)-1)
    return get_k_dist(kpoints)[k_index]

def get_k_label(kpoints, labels):
    k_index = [0]
    for i in range(len(kpoints)-1):
        if np.linalg.norm(kpoints[i]-kpoints[i+1]) < 1e-8:
            k_index.append(i)
    k_index.append(len(kpoints)-1)
    return labels[k_index]

def scissor_shift(energy, scissor):
    import copy
    energy_copy = copy.deepcopy(energy)
    for i in range(energy.shape[0]):
        if energy[i] > 0:
            energy_copy[i] = energy[i]+scissor
    return energy_copy

##This is just for bands with up spin(in a spin-polarized calculation for non-magnetic materials)
def draw_bands(bs, ylim):
    kpoints = [i.cart_coords for i in bs.kpoints]
    kpoints_dist = get_k_dist(kpoints)
    label = get_k_label(kpoints, np.array([i.label for i in bs.kpoints]))
    bands = [bs.bands[i] for i in bs.bands]
    bands_up = np.array(bands[0]) - bs.get_vbm()['energy']
    print(bands_up.shape)
    for i in range(len(bands_up)):
        plt.plot(kpoints_dist, scissor_shift(bands_up[i], scissor),'b')
    plt.plot([0, kpoints_dist[-1]], [0, 0], 'k--')
    for i in get_special_k(kpoints):
        plt.plot([i, i], ylim, 'k')
    plt.xlim([0, kpoints_dist[-1]])
    plt.ylim(ylim)
    label_loc = get_special_k(kpoints)
    label[label=='\\Gamma'] = r'$\Gamma$'
    plt.xticks(label_loc, label)
    plt.ylabel('Energy (eV)')

## This is just for SOC bands ksplit calculations
def draw_bands_ksplit(kpoints_obj, bands, vbm, ylim):
    kpoints = [i.cart_coords for i in kpoints_obj]
    kpoints_dist = get_k_dist(kpoints)
    label = get_k_label(kpoints, np.array([i.label for i in kpoints_obj]))
    bands = bands - vbm
    for i in range(len(bands)):
        plt.plot(kpoints_dist, bands[i],'b')
    plt.plot([0, kpoints_dist[-1]], [0, 0], 'k--')
    for i in get_special_k(kpoints):
        plt.plot([i, i], ylim, 'k')
    plt.xlim([0, kpoints_dist[-1]])
    plt.ylim(ylim)
    label_loc = get_special_k(kpoints)
    label[label=='\\Gamma'] = r'$\Gamma$'
    plt.xticks(label_loc, label)
    plt.ylabel('Energy (eV)')

#draw total DOS only
def draw_total_dos(v, xlim, ylim=[-8, 12], scissor=0, spinpolarized=True):
    vbm = v.eigenvalue_band_properties[2]
    comp_dos = v.complete_dos
    elements = [i for i in comp_dos.get_element_dos()]
    for element in elements:
        el_dos = comp_dos.get_element_dos()[element]
        p = plt.plot(el_dos.get_densities(spin=Spin.up), scissor_shift(el_dos.energies-vbm, scissor), label = element)
        if spinpolarized:
            plt.plot(-np.array(el_dos.get_densities(spin=Spin.down)), scissor_shift(el_dos.energies-vbm, scissor), color=p[0].get_color())

        
    tdos = v.tdos
    plt.plot(tdos.get_densities(spin=Spin.up), scissor_shift(el_dos.energies-vbm, scissor), 'k', label = 'Total')
    if spinpolarized:
        plt.plot(-np.array(tdos.get_densities(spin=Spin.down)), scissor_shift(el_dos.energies-vbm, scissor), 'k')

    plt.legend(prop={'size': 15})
    plt.ylim(ylim)
    plt.xlim(xlim)
    plt.plot(xlim, [0, 0], 'k--')
    plt.xlabel('DOS (states/eV)')
    
    #plt.yticks([])

#draw pDOS only
def draw_pdos(v, element_str, xlim, ylim=[-8, 12], scissor=0, spinpolarized=True):
    vbm = v.eigenvalue_band_properties[2]
    comp_dos = v.complete_dos
    elements = [i for i in comp_dos.get_element_dos()]
    #element = elements[element_no]
    element = Element(element_str)
    pdos = comp_dos.get_element_spd_dos(element)

    orbitals = [i for i in pdos]
    for o in orbitals:
        p_dos = pdos[o]
        p = plt.plot(p_dos.get_densities(Spin.up), scissor_shift(p_dos.energies-vbm, scissor), label = element.symbol + '-' + o.name)
        if spinpolarized:
            plt.plot(-np.array(p_dos.get_densities(Spin.down)), scissor_shift(p_dos.energies-vbm, scissor), color=p[0].get_color())

    plt.legend(prop={'size': 15})
    plt.ylim(ylim)
    plt.xlim(xlim)
    plt.plot(xlim, [0, 0], 'k--')
    plt.xlabel('DOS (states/eV)')
    plt.yticks([])

def get_sites(v):
    struc = v.structures[-1]
    sites = [i for i in struc.sites]
    return sites

def draw_site_dos(v, site_no, xlim, ylim=[-8, 12], scissor=0, spinpolarized=True):
    vbm = v.eigenvalue_band_properties[2]
    comp_dos = v.complete_dos
    sites = get_sites(v)
    site = sites[site_no-1]
    pdos = comp_dos.get_site_spd_dos(site)

    print(site.specie.symbol+' '+str(site_no))
    print(site.frac_coords)

    orbitals = [i for i in pdos]
    for o in orbitals:
        p_dos = pdos[o]
        p = plt.plot(p_dos.get_densities(Spin.up), scissor_shift(p_dos.energies-vbm, scissor), label = site.specie.symbol + str(site_no) + '-' + o.name)
        if spinpolarized:
            plt.plot(-np.array(p_dos.get_densities(Spin.down)), scissor_shift(p_dos.energies-vbm, scissor), color=p[0].get_color())

    plt.legend(prop={'size': 15})
    plt.ylim(ylim)
    plt.xlim(xlim)
    plt.plot(xlim, [0, 0], 'k--')
    plt.xlabel('DOS (states/eV)')
    plt.yticks([])



def add_delete_nbands(nbands, bands):
    nbands_2 = bands.shape[0]
    if nbands_2 < nbands:
        bands = np.concatenate((bands, np.ones((nbands-nbands_2, bands.shape[1]))*100))
    if nbands_2 > nbands:
        bands = bands[nbands, :]
    return bands