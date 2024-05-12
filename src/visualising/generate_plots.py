import matplotlib.pyplot as plt
import numpy as np


def plot_unit_sphere(ax, number_of_points=100):
    u = np.linspace(0, 2 * np.pi, number_of_points)
    v = np.linspace(0, np.pi, number_of_points)
    x = np.outer(np.cos(u), np.sin(v))
    y = np.outer(np.sin(u), np.sin(v))
    z = np.outer(np.ones(np.size(u)), np.cos(v))
    ax.plot_surface(x, y, z, color="b", alpha=0.2)



if __name__ == "__main__":
    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    plot_unit_sphere(ax)
    plt.savefig("unit_sphere.png")
