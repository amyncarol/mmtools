from pymatgen.io.vasp.inputs import Poscar
from pymatgen.io.xyz import XYZ
import os

def largercube(cube_length):
    """
    read POSCAR in the current directory and 
    output another POSCAR with larger cube size
    """
    dir_path = os.path.dirname(os.path.realpath(__file__))
    poscar = Poscar.from_file(os.path.join(dir_path, 'POSCAR'))
    structure = poscar.structure
    XYZ(structure).write_file(os.path.join(dir_path, 'POSCAR.xyz'))
    molecule = XYZ.from_file(os.path.join(dir_path, 'POSCAR.xyz')).molecule
    a = b = c = cube_length
    structure = molecule.get_boxed_structure(a, b, c)
    Poscar(structure).write_file(os.path.join(dir_path, 'POSCAR_larger.vasp'))

if __name__=='__main__':
    largercube(20)