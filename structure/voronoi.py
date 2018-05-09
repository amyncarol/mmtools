from pymatgen.io.vasp.inputs import Poscar
from scipy.spatial import Voronoi, ConvexHull
import numpy as np
from pymatgen.analysis.structure_analyzer import VoronoiConnectivity
import copy

class ColorVonoroi(object):
    """
    given a structure, find the vonoroi polyhedra for each atom and color each face of the polyhedra
    according to the neighbor.

    I try to use this to visualize the difference between polymorphs. Different polymorphs must have 
    different shapes of vonoroi polyhedra. The ultimate goal is to find a way to quantify the difference 
    using polyhedra.

    Another way to visualize the difference is to use octahedra as building blocks and look at the atom 
    sharing scheme, vertice sharing, edge sharing, face sharing, I hope there is a method to quantify the 
    different sharing schemes for different polymorphs. 

    Another way is close packing, which I don't quite like since it has less emphasis on 3D bonding. 

    Complete this, just for fun.
    """
    def __init__(self, structure):
        self.structure = structure
        self.primitive = structure.get_primitive_structure()
        self.supercell = copy.deepcopy(self.primitive)
        self.supercell.make_supercell(3)
        self.cart_coords = self.supercell.cart_coords
        self.vor = Voronoi(self.cart_coords)
        self.elements = [i.symbol for i in set(self.primitive.species)]
        self.color = ['r', 'g', 'b', 'c', 'y']
        self.color_dict = {}
        for i, e in enumerate(self.elements):
            self.color_dict[e] = self.color[i]
        #print(self.supercell)

    def get_polyhedron_volume_area(self, isite):
        """
        Returns:
             (volume, area) for site i
        """
        region_index = self.vor.point_region[isite]
        vertices_index = self.vor.regions[region_index]
        if -1 in vertices_index:
            #print("boudary atoms, give it up")
            return 0, 0
        else:
            region_vertices = np.zeros((len(vertices_index), 3))
            for i in range(len(vertices_index)):
                region_vertices[i, :] = self.vor.vertices[vertices_index[i]]
            hull = ConvexHull(region_vertices)
            return round(hull.volume, 3), round(hull.area, 3)

    def get_distinct_atoms(self):
        """
        return the distinct atom indices together with its voronoi volume and area

        Returns:
            {atom_index: (volume, area)}
        """
        isites = [i for i in range(len(self.cart_coords))]
        atom_dict = {}
        for i in isites:
            volume, area = self.get_polyhedron_volume_area(i)
            if volume != 0 and (volume, area) not in atom_dict.values():
                atom_dict[i] = (volume, area)
        return atom_dict

    def get_polyhedron_faces(self, isite):
        """
        Returns:
            a list of faces
            each face is tuple (element in the mirror, a list of vertices of the face)
        """
        ridge_points = self.vor.ridge_points
        ridge_vertices = self.vor.ridge_vertices
        ridge_index = []
        ridge_shape = []
        for i in range(len(ridge_points)):
            if isite in ridge_points[i]:
                ridge_index.append(i)
        for i in ridge_index:
            other_point_index = [p for p in ridge_points[i] if p not in [isite]][0]
            other_element_symbol = self.supercell.sites[other_point_index].specie.symbol
            vertices = [self.vor.vertices[j] for j in ridge_vertices[i]]
            ridge_shape.append((other_element_symbol, vertices))
        return ridge_shape

    def plot_polyhedron_faces(self, isite):
        pass

if __name__ == '__main__':
    struc = Poscar.from_file("/Users/yao/Downloads/POSCAR.mp-614013_CsSnI3.vasp").structure
    cv = ColorVonoroi(struc)
    print(cv.get_distinct_atoms())

    print(cv.get_polyhedron_faces(40))

    #print(cv.primitive)
    #print(len(cv.cart_coords))

# vc = VoronoiConnectivity(struc)
# print(vc.connectivity_array[1])
