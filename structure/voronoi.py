from pymatgen.io.vasp.inputs import Poscar
from scipy.spatial import Voronoi, ConvexHull
import numpy as np
from pymatgen.analysis.structure_analyzer import VoronoiConnectivity
import copy
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from pymatgen.core.structure import Molecule
from numpy.linalg import norm
from mayavi.mlab import *

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
        self.structure = structure.get_primitive_structure()    
        self.n = len(self.structure.sites)  #number of atoms in structure  
        self.expanded_structure = ColorVonoroi.expand_structure(self.structure)
        self.cart_coords = [site.coords for site in self.expanded_structure]

        elements = list(set(self.structure.species))
        self.elements = [i.symbol for i in elements]
        electronegativity = [i.X for i in elements]
        self.elements = [x for _, x in sorted(zip(electronegativity, self.elements), reverse=True)]
        self.color = [(1, 0, 0), (1, 1, 0), (0, 1, 0), (0, 1, 1), (0, 0, 1)]
        self.color_dict = {}
        for i, e in enumerate(self.elements):
            self.color_dict[e] = self.color[i]

        self.vor = Voronoi(self.cart_coords)

    @staticmethod
    def expand_structure(structure, radius=10):
        """expand the structure to include more atoms. 
        For each atom in the original structure, ensure 
        that around it a raidus of 10 angstrom there are atoms in 
        the structure. Keep the index for atoms in the original structure

        Return:
            a expanded structure with index for the original structure fixed
        """
        species = structure.species
        coords = structure.cart_coords
        molecule = Molecule(species, coords)

        a = structure.lattice.matrix[0]
        b = structure.lattice.matrix[1]
        c = structure.lattice.matrix[2]
        ni = max(int(10/norm(a)), 1)
        nj = max(int(10/norm(b)), 1)
        nk = max(int(10/norm(c)), 1)
        
        for i in range(-ni, ni+1):
            for j in range(-nj, nj+1):
                for k in range(-nk, nk+1):
                    if i == 0 and j ==0 and k==0:
                        continue
                    for site in structure.sites:
                        specie = site.specie
                        coords = site.coords + i*a+j*b+k*c
                        molecule.append(specie, coords)
        return molecule

    def get_polyhedron_hull_volume_area(self, isite):
        """
        Returns:
             (region_vertices, volume, area) for site i
        """
        region_index = self.vor.point_region[isite]
        vertices_index = self.vor.regions[region_index]
        if -1 in vertices_index:
            #print("boudary atoms, give it up")
            return 0, 0, 0
        else:
            region_vertices = np.zeros((len(vertices_index), 3))
            for i in range(len(vertices_index)):
                region_vertices[i, :] = self.vor.vertices[vertices_index[i]]
            hull = ConvexHull(region_vertices)
            return region_vertices, round(hull.volume, 2), round(hull.area, 2)

    def get_distinct_atoms(self):
        """
        return the distinct atom indices together with its voronoi volume and area

        Returns:
            {atom_index: (volume, area)}
        """
        atom_dict = {}
        for i in range(self.n):
            hull, volume, area = self.get_polyhedron_hull_volume_area(i)
            if volume != 0 and (volume, area) not in atom_dict.values():
                atom_dict[i] = (volume, area)
        return atom_dict

    def get_polyhedron_faces(self, isite):
        """
        Returns:
            a list of faces
            each face is tuple (index of element, element in the mirror, a list of vertices of the face)
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
            other_element_symbol = self.expanded_structure.sites[other_point_index].specie.symbol
            vertices = np.array([self.vor.vertices[j] for j in ridge_vertices[i]])
            ridge_shape.append((other_point_index, other_element_symbol, vertices))
        return ridge_shape

    @staticmethod
    def in_hull(p, hull):
        """
        Test if points in `p` are in `hull`

        `p` should be a `NxK` coordinates of `N` points in `K` dimensions
        `hull` is either a scipy.spatial.Delaunay object or the `MxK` array of the 
        coordinates of `M` points in `K`dimensions for which Delaunay triangulation
        will be computed
        It returns a boolean array where True values indicate points that lie in the given convex hull
        """
        from scipy.spatial import Delaunay
        if not isinstance(hull, Delaunay):
            hull = Delaunay(hull)

        return hull.find_simplex(p)>=0

    def _get_x_y_z_triangle(self, isite):
        """
        get x, y, z and triangle to plot triangle mesh (to plot using mayavi)

        Returns:
            x: list of the x coords of vertices
            y: list of the y coords of vertices
            z: list of the z coords of vertices
            triangle: list of tuples (_, _, _) which indicates indices of vertices of triangle
        """
        ridge_points = self.vor.ridge_points
        ridge_vertices = self.vor.ridge_vertices
        ridge_index = []
        ridge_shape = []
        for i in range(len(ridge_points)):
            if isite in ridge_points[i]:
                ridge_index.append(i)

        triangle = {e:[] for e in self.elements}
        for i in ridge_index:
            other_point_index = [p for p in ridge_points[i] if p not in [isite]][0]
            other_element_symbol = self.expanded_structure.sites[other_point_index].specie.symbol
            vertices_index = [j for j in ridge_vertices[i]]
            polyn = len(vertices_index)
            for (p, q) in zip(range(1, polyn-1), range(2, polyn)):
               triangle[other_element_symbol].append((vertices_index[0], vertices_index[p], vertices_index[q]))
            #triangle[other_element_symbol].append((vertices_index[0], vertices_index[1], vertices_index[2]))
        x, y, z = self.vor.vertices[:, 0], self.vor.vertices[:, 1], self.vor.vertices[:, 2]
        return x, y, z, triangle

    @staticmethod
    def cut_large_number(x):
        return (abs(x)<1e6) * x + (x>=1e6) * 1e6 + (x<=-1e6) * (-1e6)

    def plot_polyhedra(self, isite):
        """
        plot the voronoi polyhedra for isite using mayavi
        """
        _, vol, area = self.get_polyhedron_hull_volume_area(isite)

        x, y, z, triangle = self._get_x_y_z_triangle(isite)

        ##overflow here, mayavi cannot handle, manually cut the large number to 1e6
        x = ColorVonoroi.cut_large_number(x)
        y = ColorVonoroi.cut_large_number(y)
        z = ColorVonoroi.cut_large_number(z)
        
        i = 0.05
        for e in triangle:
            if triangle[e]:
                triangular_mesh(x, y, z, triangle[e], color=self.color_dict[e])
                text(0.05, i, e, color=self.color_dict[e], width=0.08)
                i+=0.12
        title('isite: {}, center atom: {}'.format(isite, self.expanded_structure[isite].specie.symbol))
        text(0.3, 0.05, 'volume:{}, \n area:{}'.format(vol, area), width=0.4)
        show()

    
if __name__ == '__main__':
    pass
