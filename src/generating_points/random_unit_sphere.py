import numpy as np
import pandas as pd


def generate_n_poins_on_sphere(n):
    """Generate n random points on the surface of a unit sphere."""

    # Generate random numbers
    theta = 2 * np.pi * np.random.rand(n)
    phi = np.arccos(1 - 2 * np.random.rand(n))

    # Spherical to Cartesian coordinates conversion
    x = np.sin(phi) * np.cos(theta)
    y = np.sin(phi) * np.sin(theta)
    z = np.cos(phi)

    # Create DataFrame
    data = pd.DataFrame({"Theta": theta, "Phi": phi, "x": x, "y": y, "z": z})

    #plot
    fig = plt.figure()
