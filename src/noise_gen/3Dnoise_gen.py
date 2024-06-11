import math
import random

import noise
import numpy as np
from opensimplex import OpenSimplex
from PIL import Image, ImageDraw
from sphere_mesh import MeshFractionalSphere

# https://loady.one/blog/terrain_mesh.html
# https://github.com/rbehrou/EasyFEM/blob/master/meshFractionalSphere/meshFractionalSphere.py
MAP_SIZE = (512, 512, 512)
SCALE = 256
EXPO_HEIGHT = 2
COLORS = {
    "grass": (34, 139, 34),
    "forest": (0, 100, 0),
    "sand": (238, 214, 175),
    "water": (65, 105, 225),
    "rock": (139, 137, 137),
    "snow": (255, 250, 250),
}


def create_sphere(
    cx: float, cy: float, cz: float, r: float, resolution: int = 360
) -> np.ndarray:
    """
    Create a sphere with center (cx, cy, cz) and radius r.
    https://stackoverflow.com/questions/61048426/python-generating-3d-sphere-in-numpy
    Parameters:
        cx (float): The x-coordinate of the center of the sphere.
        cy (float): The y-coordinate of the center of the sphere.
        cz (float): The z-coordinate of the center of the sphere.
        r (float): The radius of the sphere.
        resolution (int, optional): The resolution of the sphere. Defaults to 360.

    Returns:
        np.ndarray: A 3D array representing the coordinates of the sphere.
    """

    phi = np.linspace(0, 2 * np.pi, 2 * resolution)
    theta = np.linspace(0, np.pi, resolution)

    theta, phi = np.meshgrid(theta, phi)

    r_xy = r * np.sin(theta)
    x = cx + np.cos(phi) * r_xy
    y = cy + np.sin(phi) * r_xy
    z = cz + r * np.cos(theta)

    return np.stack([x, y, z])


test = create_sphere(0, 0, 0, 10, 1000)
print(test.shape)

def generate_sphere_vertices(sphere_points: np.ndarray) -> :