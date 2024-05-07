# planet_pyvista.py
import numpy as np
import pyvista as pv


class Planet:
    def __init__(self, resolution=10):
        self.resolution = resolution
        self.mode = "surface"  # Modes: 'points', 'surface'
        self.mesh = self.generate_mesh()
        self.plotter = pv.Plotter()
        self.add_mesh()

    def generate_mesh(self):
        # Create a sphere mesh using the specified resolution
        sphere = pv.Sphere(
            theta_resolution=self.resolution, phi_resolution=self.resolution
        )
        return sphere

    def add_mesh(self):
        # Add the mesh to the plotter with the current mode
        if self.mode == "points":
            self.plotter.add_mesh(
                self.mesh, render_points_as_spheres=True, point_size=10
            )
        else:
            self.plotter.add_mesh(self.mesh)

    def toggle_points(self):
        # Toggle between 'surface' and 'points' mode
        self.mode = "points" if self.mode == "surface" else "surface"
        self.plotter.clear()
        self.add_mesh()

    def show(self):
        self.plotter.add_text(
            "Press T to toggle points/surface mode", position="upper_left", font_size=12
        )
        self.plotter.add_key_event("t", self.toggle_points)
        self.plotter.show()


if __name__ == "__main__":
    planet = Planet(resolution=20)
    planet.show()
