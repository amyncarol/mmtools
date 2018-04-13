"""
This is for generating the inputs for cgcnn code writen by Tian Xie

The structure of the root_dir should be:

root_dir
├── id_prop.csv
├── atom_init.json
├── id0.cif
├── id1.cif
├── ...

id_prop.csv: a CSV file with two columns. The first column recodes a unique ID for each crystal, 
and the second column recodes the value of target property. If you want to predict material properties
 with predict.py, you can put any number in the second column. (The second column is still needed.)

atom_init.json: a JSON file that stores the initialization vector for each element. An example of 
atom_init.json is cgcnn/data/sample-regression/atom_init.json, which should be good for most applications.

ID.cif: a CIF file that recodes the crystal structure, where ID is the unique ID for the crystal.
"""

import os
import time
import pandas as pd
from pymatgen.io.cif import CifWriter

from models.data.solid_solution.generate_2116_solid_solution_vasp_input import get_structure_dict, get_pymatgen_formula

def generate_dataset(energy_file, structure_file, root_dir):
    """
    given energy_file and structure_file, generate the all .cif files and id_prop.csv file

    Args:

    energy_file: the file that contains the energy value for each structure

        Cs2Ag1Al1F6 -3.06382545971 eV
        Cs2Ag1As1F6 -2.46506153321 eV
        Cs2Ag1Au1Br6 -0.897941594 eV

    structure_file: the structure.pkl file

    root_dir: folder we should write all files

    Writes: 

    cifs and id_prop.csv defined above
    """
    structure_dict = get_structure_dict(structure_file)
    
    energies = []

    with open(energy_file, 'r') as f:

        start_time = time.time()
        count = 0

        lines = f.readlines()
        for line in lines:
            count += 1
            if count%10 == 0:
                print("Processed {} items, using {} sec.".format(count, time.time()-start_time))

            compound = line.split(' ')[0]
            formula = get_pymatgen_formula(compound)

            if formula not in structure_dict:
                print("something wrong with {}".format(formula))

            else:
                ##ID for this structure
                ID = len(energies)

                ##energy
                energies.append(float(line.split(' ')[1]))
               
                ##the structure
                structure = structure_dict[formula]
                cifwriter = CifWriter(structure)
                cifwriter.write_file(os.path.join(root_dir, '{}.cif'.format(ID)))

    #pd.DataFrame(energies).to_csv(os.path.join(root_dir, 'id_prop.csv'))
    pd.DataFrame(energies, columns=['formation_energy']).to_csv(os.path.join(root_dir, 'id_prop.csv'), header=False)

if __name__=='__main__':
    generate_dataset('/Users/yao/Google Drive/data/2116/solidsolution/solid_solution_ml/complete/cnn_learning/formation_energy', \
        '/Users/yao/Google Drive/data/2116/solidsolution/solid_solution_ml/complete/cnn_learning/structure.pkl', \
        '/Users/yao/Google Drive/data/2116/solidsolution/solid_solution_ml/complete/cnn_learning')

