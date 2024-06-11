# pylint: disable=unused-import

import meshio
import numpy as np
from scipy.spatial import SphericalVoronoi


class MeshSimpleSphere:
    """
    The mesh class
    """

    def __init__(
        self, radius, num_theta, num_phi, theta_fraction, phi_fraction, center=(0, 0, 0)
    ):
        """
        Constructor
        """
        self.radius = radius
        self.num_theta = num_theta
        self.num_phi = num_phi
        self.theta_fraction = theta_fraction
        self.phi_fraction = phi_fraction
        self.center = center

        # Generate mesh
        self.coords, self.topo = self.generate_mesh()

    def generate_mesh(self):
        points = self.coords - self.center
        theta_limit = np.pi * self.theta_fraction
        phi_limit = np.pi * self.phi_fraction

        filtered_points = []
        for p in points:
            x, y, z = p
            theta = np.arccos(z / self.radius)
            phi = np.arctan2(y, x)
            if 0 <= theta <= theta_limit and -phi_limit <= phi <= phi_limit:
                filtered_points.append(p)

        coords = np.array(filtered_points) + self.center
        sv = SphericalVoronoi(coords, self.radius, self.center)
        # sv.sort_vertices_of_regions()

        triangles = []
        for region in sv.regions:
            for i in range(1, len(region) - 1):
                triangles.append([region[0], region[i], region[i + 1]])

        return np.array(triangles)

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


if __name__ == "__main__":
    main_simple()
