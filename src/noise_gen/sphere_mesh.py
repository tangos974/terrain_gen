# pylint: disable=unused-import

import random

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import meshio
import noise
import numpy as np
from scipy.spatial import ConvexHull  # pylint: disable=no-name-in-module
from scipy.spatial import KDTree  # pylint: disable=no-name-in-module
from scipy.spatial import SphericalVoronoi


class MeshFractionalSphere:

    def __init__(
        self,
        radius: float,
        num_points: int,
        theta_fraction: float,
        phi_fraction: float,
        center: np.ndarray = np.array([0, 0, 0]),
    ):
        """
        Initializes a MeshFractionalSphere object with the given parameters.

        Args:
            radius (float): The radius of the sphere.
            num_points (int): The number of points to generate on the sphere.
            theta_fraction (float): The fraction of the sphere to include in the theta direction.
            phi_fraction (float): The fraction of the sphere to include in the phi direction.
            center (np.ndarray): The center of the sphere.

        Attributes:
            radius (float): The radius of the sphere.
            num_points (int): The number of points to generate on the sphere.

            theta_fraction (float): The fraction of the sphere to include in the theta direction.
            0 < fraction <= 1
            phi_fraction (float): The fraction of the sphere to include in the phi direction.
            0 < fraction <= 1
            Examples:
                # 1/8 of sphere: phi_fraction = 1/4, theta_fraction = 1/2
                # 1/4 of sphere: phi_fraction = 1/2, theta_fraction = 1/2
                # 1/2 of sphere: phi_fraction = 1/2, theta_fraction = 1.0
                # full sphere  : phi_fraction = 1.0, theta_fraction = 1.0

            coords (numpy.ndarray): The generated points on the sphere.
            topo (numpy.ndarray): The triangular surface mesh of the sphere.

        Returns:
            None
        """

        self.radius = radius
        self.num_points = num_points
        self.theta_fraction = theta_fraction
        self.phi_fraction = phi_fraction
        self.center = center

        # Generate coordinates
        self.coords = self.generate_points()

    def generate_points(self) -> np.ndarray:
        """
        Generates points on a sphere using the Fibonacci lattice method.

        Returns:
        --------
        np.ndarray : The generated points on the sphere.
        Shape is (num_points, 3).
        """
        points = self._fibonacci_sphere(self.num_points, self.radius)

        # Shift points to the specified center
        points += self.center

        return points

    def generate_spherical_mesh(self):
        """
        Generates a triangular mesh of the sphere.

        Returns:
        --------
        numpy.ndarray: The triangular surface mesh of the sphere.
        Shape is (num_triangles, 3).
        """
        points = self.coords

        # Shift points back to the origin for calculations
        shifted_points = points - self.center

        # Filter points based on theta_fraction and phi_fraction
        theta_limit = np.pi * self.theta_fraction
        phi_limit = np.pi * self.phi_fraction

        filtered_points = []
        for p in shifted_points:
            x, y, z = p
            theta = np.arccos(z / self.radius)
            phi = np.arctan2(y, x)
            if 0 <= theta <= theta_limit and -phi_limit <= phi <= phi_limit:
                filtered_points.append(p)

        coords = np.array(filtered_points) + self.center

        # Perform ConvexHull to get triangular surface mesh
        hull = ConvexHull(coords)
        print(hull.simplices, hull.simplices.shape, type(hull.simplices))
        triangles = hull.simplices

        return triangles

    def generate_spherical_voronoi_regions_mesh(self):
        sv = SphericalVoronoi(self.coords, self.radius, self.center)
        sv.sort_vertices_of_regions()
        print(sv.regions)

        voronoi_vertices = sv.vertices
        voronoi_regions = sv.regions
        cells = []

        for region in voronoi_regions:
            if len(region) >= 3:  # Only consider regions with 3 or more vertices
                cells.append(region)

        return voronoi_vertices, cells

    def _fibonacci_sphere(self, samples: int = 1000, radius: float = 1) -> np.ndarray:
        """
        Generates points on a sphere using the Fibonacci lattice method.

        Parameters:
        -----------
        samples : int, optional
            The number of points to generate (default is 1000).

        radius : float, optional
            The radius of the sphere (default is 1).

        Returns:
        --------
        np.ndarray : The generated points on the sphere.
        Shape is (samples, 3).
        """
        points = []
        phi = np.pi * (3.0 - np.sqrt(5.0))  # golden angle in radians

        for i in range(samples):
            y = 1 - (i / float(samples - 1)) * 2  # y goes from 1 to -1
            r = np.sqrt(1 - y * y)  # radius at y

            theta = phi * i  # golden angle increment

            x = np.cos(theta) * r
            z = np.sin(theta) * r

            points.append([x * radius, y * radius, z * radius])

        return np.array(points)

    def to_vtk(self):
        """
        Write mesh to VTK
        """
        # Create mesh
        mesh = meshio.Mesh(self.coords, [("triangle", self.generate_spherical_mesh())])
        return mesh

    def to_obj_triangles(self):
        """
        Write mesh to OBJ, ensuring triangular faces.
        """
        # Create mesh
        mesh = meshio.Mesh(self.coords, [("triangle", self.generate_spherical_mesh())])
        return mesh

    def to_obj_voronoi(self):
        voronoi_vertices, voronoi_cells = self.generate_spherical_voronoi_regions_mesh()
        cells = []
        for region in voronoi_cells:
            cells.append([region])
        voronoi_cells = np.array(cells, dtype=object)
        mesh = meshio.Mesh(voronoi_vertices, [("polygon", voronoi_cells)])
        return mesh


def main_fibo():

    radius = 8.4567
    num_points = 10
    center = np.array([1.0, 20.0, 3.0])

    phi_fraction = 1.00
    theta_fraction = 1.00
    sphere = MeshFractionalSphere(
        radius, num_points, theta_fraction, phi_fraction, center
    )

    # Write into obj
    # sphere.to_obj_triangles().write("mesh_sphere.obj")
    sphere.to_obj_voronoi().write("mesh_sphere_voronoi.vtk")

    sphere.generate_spherical_voronoi_regions_mesh()

    # Write into vtk
    # sphere.to_vtk().write("mesh_sphere.vtk")


if __name__ == "__main__":
    main_fibo()
