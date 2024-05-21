import random

import meshio
import numpy as np
from scipy.spatial import ConvexHull  # pylint: disable=no-name-in-module


class MeshSimpleSphere:
    """
    The mesh class
    """

    def __init__(self, radius, num_theta, num_phi, theta_fraction, phi_fraction):
        """
        Constructor
        """
        self.radius = radius
        self.num_theta = num_theta
        self.num_phi = num_phi
        self.theta_fraction = theta_fraction
        self.phi_fraction = phi_fraction

        # Generate mesh
        self.coords, self.topo = self.generate_mesh()

    def generate_mesh(self):
        """
        Mesh sphere
        """
        # Create the mesh grid (0 <= theta <= pi, 0 <= phi <= 2*pi)
        theta = np.linspace(0.0, 1.0 * np.pi * self.theta_fraction, self.num_theta)
        phi = np.linspace(0.0, 2.0 * np.pi * self.phi_fraction, self.num_phi)

        # Generate the nodal coordinates
        th, ph = np.meshgrid(theta, phi, indexing="ij")
        x = self.radius * np.sin(th) * np.cos(ph)
        y = self.radius * np.sin(th) * np.sin(ph)
        z = self.radius * np.cos(th)

        # Reshape the arrays
        coords = np.vstack([x.ravel(), y.ravel(), z.ravel()]).T

        # Generate the elemental connectivity array manually
        topo = []
        for i in range(self.num_theta - 1):
            for j in range(self.num_phi - 1):
                n0 = i * self.num_phi + j
                n1 = n0 + 1
                n2 = n0 + self.num_phi
                n3 = n2 + 1

                # Add two triangles for each quadrilateral face
                topo.append([n0, n1, n3])
                topo.append([n0, n3, n2])

        return coords, np.array(topo)

    def to_vtk(self):
        """
        Write mesh to VTK
        """
        # Create mesh
        mesh = meshio.Mesh(self.coords, [("tetra", self.topo)])
        return mesh

    def to_obj(self):
        """
        Write mesh to OBJ, ensuring triangular faces.
        """
        # Create mesh
        mesh = meshio.Mesh(self.coords, [("triangle", self.topo)])
        return mesh


def main_simple():
    """
    main entry point
    """
    # Example usage:
    radius = 5
    num_theta = 200
    num_phi = 200

    phi_fraction = 1.00  # Fraction of mesh is phi direction, 0 < fraction <= 1
    theta_fraction = 1.00  # fraction of mesh in theta direction, 0 < fraction <= 1
    mesh = MeshSimpleSphere(radius, num_theta, num_phi, theta_fraction, phi_fraction)

    # Write into obj
    mesh.to_obj().write("mesh_sphere.obj")


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

        # Generate mesh
        self.topo = self.generate_mesh()

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

    def generate_mesh(self):
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
        triangles = hull.simplices

        return triangles

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
        mesh = meshio.Mesh(self.coords, [("triangle", self.topo)])
        return mesh

    def to_obj(self):
        """
        Write mesh to OBJ, ensuring triangular faces.
        """
        # Create mesh
        mesh = meshio.Mesh(self.coords, [("triangle", self.topo)])
        return mesh


def main_fibo():

    radius = 8.4567
    num_points = 100000
    center = np.array([1.0, 2.0, 3.0])

    phi_fraction = 1.00
    theta_fraction = 1.00
    sphere = MeshFractionalSphere(
        radius, num_points, theta_fraction, phi_fraction, center
    )

    # Write into obj
    sphere.to_obj().write("mesh_sphere.obj")

    # Write into vtk
    # sphere.to_vtk().write("mesh_sphere.vtk")


if __name__ == "__main__":
    main_fibo()
