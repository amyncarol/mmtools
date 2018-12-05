from pymatgen.io.vasp.inputs import Poscar, Kpoints
from pymatgen.io.vasp.outputs import Vasprun
import os
from shutil import copyfile, move, rmtree
from glob import glob
from xml.etree.ElementTree import ParseError

def convergence_k_test(folder, start_k_density):
    """
    given a folder of INCAR, POSCAR, POTCAR
    generate the several folders for different k mesh densities

    Args:
        folder: the working folder under which has the input files
        start_k_density: the start k density we use to generate KPOINTS files

    Writes:
        several folders with same vasp inputs except for the KPOINTS
    """
    #clean the folder
    for path in glob(folder+'/kpoints-*'):
        rmtree(path)

    #if CONTCAR is input instead of POSCAR
    if os.path.exists(folder + '/CONTCAR'):
        move(folder + '/CONTCAR', folder + '/POSCAR')
        
    #generate 5 k-points file and create folders
    struc = Poscar.from_file(folder + '/POSCAR').structure
    k_set = set()
    i = 0
    while len(k_set) < 5:
        kppa = start_k_density + i*150
        i += 1
        kpoints = Kpoints.automatic_gamma_density(struc, kppa)
        if kpoints.kpts[0][0] in k_set:
            continue
        else:
            k_set.add(kpoints.kpts[0][0])           
            subfolder = folder+'/kpoints-'+str(kpoints.kpts[0][0])
            os.mkdir(subfolder)
            print(subfolder)
            copyfile(folder+'/INCAR', subfolder+'/INCAR')
            copyfile(folder+'/POSCAR', subfolder+'/POSCAR')
            copyfile(folder+'/POTCAR', subfolder+'/POTCAR')
            kpoints.write_file(subfolder+'/KPOINTS')

def convergence_encut_test(folder, start_encut):
    """
    given a folder of INCAR, POSCAR, POTCAR, KPOINTS
    generate the several folders for different encut

    Args:
        folder: the working folder under which has the input files
        start_encut: the start encut we use to generate INCAR files

    Writes:
        several folders with same vasp inputs except for encut
    """

    #clean the folder
    for path in glob(folder+'/*-*'):
        # if os.path.isfile(path):
        #   os.remove(path)
        # if os.path.isdir(path):
        rmtree(path)

    #if CONTCAR is input instead of POSCAR
    if os.path.exists(folder + '/CONTCAR'):
        move(folder + '/CONTCAR', folder + '/POSCAR')
        
    #generate 5 k-points file and create folders
    for i in range(5):
        encut = start_encut + i*40      
        subfolder = folder+'/encut-'+str(encut)
        os.mkdir(subfolder)
        print(subfolder)
        copyfile(folder+'/KPOINTS', subfolder+'/KPOINTS')
        copyfile(folder+'/POSCAR', subfolder+'/POSCAR')
        copyfile(folder+'/POTCAR', subfolder+'/POTCAR')

        #handle INCAR
        lines = []
        with open(folder+'/INCAR', 'r') as f:
            for line in f:
                if line.startswith('ENCUT'):
                    line = 'ENCUT = ' + str(encut) + '\n'
                lines.append(line)

        with open(subfolder+'/INCAR', 'w') as f:
            for line in lines:
                f.write(line)

def convergence_report(folder, k_test = True):
    """
    print the report for convergence test, only works for cubic structure now

    Args:
        folder: in this folder there are several folders with different k meshes
        k_test: whether it is k test or encut test
    """
    tuple_list = []
    for path in glob(folder+'/*-*'):
        try:
            vasprun = Vasprun(path+'/vasprun.xml')
            lattice_constant_a = vasprun.final_structure.lattice.a 
            energy = vasprun.final_energy
            if k_test:
                parameter = vasprun.kpoints.kpts[0][0]
            else:
                parameter = vasprun.incar.get("ENCUT")
            tuple_list.append((parameter, energy, lattice_constant_a))
        except FileNotFoundError:
            print('{}/vasprun.xml not found'.format(path))
        except ParseError:
            print('cannot parse {}/vasprun.xml'.format(path))

    tuple_list.sort(key=lambda x:x[0])

    if k_test:
        string = 'k'
    else:
        string = 'encut'

    for tup in tuple_list:
        print('for {} = {}, the total energy is {} from convergence, and the lattice is {} Angstrom from convergence'\
                .format(string, tup[0], tup[1]-tuple_list[-1][1], tup[2]-tuple_list[-1][2]))

if __name__=='__main__':
    # #generate k test inputs
    # for folder in glob('/Users/yao/Google Drive/data/2116/data/convergence_k/*'):
    #   print(folder)
    #   convergence_k_test(folder, 10)

    # #generate encut test inputs
    # for folder in glob('/Users/yao/Google Drive/data/2116/data/convergence_encut/*'):
    #   print(folder)
    #   convergence_encut_test(folder, 440)

    #convergence k report
    for folder in glob('/Users/yao/Google Drive/data/2116/data/convergence_k/*'):
        print(folder)
        convergence_report(folder, True)
        print('----------------')

    #convergence encut report
    for folder in glob('/Users/yao/Google Drive/data/2116/data/convergence_encut/*'):
        print(folder)
        convergence_report(folder, False)
        print('----------------')


    
