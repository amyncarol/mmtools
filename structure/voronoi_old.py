from pymatgen.io.vasp.inputs import Poscar
from scipy.spatial import Voronoi, ConvexHull
import numpy as np
from pymatgen.analysis.structure_analyzer import VoronoiConnectivity
import copy
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from pymatgen.core.structure import Molecule
from numpy.linalg import norm

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

        self.elements = [i.symbol for i in set(self.structure.species)]
        self.color = ['r', 'g', 'b', 'c', 'y']
        self.color_dict = {}
        for i, e in enumerate(self.elements):
            self.color_dict[e] = self.color[i]
        print(self.color_dict)

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


    #################Plotting#####################################################################
    @staticmethod
    def _plot_vertices(ridge_shape, ax):
        for face in ridge_shape:
            points = face[2]
            ax.scatter(xs = points[:, 0], ys = points[:, 1], zs=points[:, 2])

    # @staticmethod
    # def _plot_surface(shape):
    #     fig = plt.figure()
    #     ax = fig.add_subplot(111, projection='3d')
    #     for face in shape:
    #         points = face[1]
    #         ax.scatter(xs = points[:, 0], ys = points[:, 1], zs=points[:, 2])
    #     plt.show()

    @staticmethod
    def _get_plane_equation(coords):
        """
        At least three points coords should be provide

        Returns:
            a function that takes x, y and output z
        """
        cofficient = np.linalg.solve(coords[:3, :], np.ones((3, 1)))

        index = np.argmax(np.abs(cofficient))

        if index == 0:
            def _equation(y, z):
                return (1-cofficient[1]*y-cofficient[2]*z)/cofficient[0]
        elif index == 1:
            def _equation(x, z):
                return (1-cofficient[0]*x-cofficient[2]*z)/cofficient[1]
        elif index == 2:
            def _equation(x, y):
                return (1-cofficient[0]*x-cofficient[1]*y)/cofficient[2]
        return index, _equation


    def plot_plane(self, coords, ax, hull, element):
        index, equation = ColorVonoroi._get_plane_equation(coords)
        X = np.arange(min(coords[:, 0]), max(coords[:, 0])+0.25, 0.25)
        Y = np.arange(min(coords[:, 1]), max(coords[:, 1])+0.25, 0.25)
        Z = np.arange(min(coords[:, 2]), max(coords[:, 2])+0.25, 0.25)
        if index == 0:
            Y, Z = np.meshgrid(Y, Z)
            X = equation(Y, Z)
        elif index == 1:
            X, Z = np.meshgrid(X, Z)
            Y = equation(X, Z)
        elif index == 2:
            X, Y = np.meshgrid(X, Y)
            Z = equation(X, Y)

        points = np.array([X.flatten(), Y.flatten(), Z.flatten()]).T
        hull_index = ColorVonoroi.in_hull(points, hull)
        hull_index = hull_index.reshape(X.shape[0], -1)
        X[hull_index != True] = None
        Y[hull_index != True] = None
        Z[hull_index != True] = None
        # X = np.ma.masked_where(hull_index != True, X)
        # Y = np.ma.masked_where(hull_index != True, Y)
        # Z = np.ma.masked_where(hull_index != True, Z)
        color = self.color_dict[element]
        ax.plot_surface(X, Y, Z, color=color, alpha=0.5)

    def plot_polyhedron_faces(self, isite):
        center_atom = self.expanded_structure.sites[isite].specie.symbol
        hull, _, _ = self.get_polyhedron_hull_volume_area(isite)
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ridge_shape = self.get_polyhedron_faces(isite)
        for face in ridge_shape:
            vertices = face[2]
            element = face[1]
            self.plot_plane(vertices, ax, hull, element)
        ColorVonoroi._plot_vertices(ridge_shape, ax)
        plt.show()
        

if __name__ == '__main__':
    pass
