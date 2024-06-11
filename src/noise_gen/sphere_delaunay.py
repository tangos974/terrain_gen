import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from scipy.spatial import Delaunay


def spherical_to_cartesian(lat_deg, lon_deg):
    lat_rad = np.deg2rad(lat_deg)
    lon_rad = np.deg2rad(lon_deg)
    x = np.cos(lat_rad) * np.cos(lon_rad)
    y = np.cos(lat_rad) * np.sin(lon_rad)
    z = np.sin(lat_rad)
    return x, y, z


def generate_fibonacci_sphere(N, jitter):
    points = []
    s = 3.6 / np.sqrt(N)
    dz = 2.0 / N
    for k in range(N):
        z = 1 - dz / 2 - k * dz
        r = np.sqrt(1 - z * z)
        lat_deg = np.rad2deg(np.arcsin(z))
        lon_deg = (k * s / r) % (2 * np.pi)
        if jitter > 0:
            lat_deg += jitter * (np.random.rand() - 0.5)
            lon_deg += jitter * (np.random.rand() - 0.5)
        points.append(spherical_to_cartesian(lat_deg, lon_deg))
    return np.array(points)


def add_south_pole(points):
    south_pole = np.array([0, 0, -1])
    return np.vstack([points, south_pole])


def generate_delaunay(points):
    projected_points = points[:, :2] / (1 - points[:, 2, np.newaxis])
    delaunay = Delaunay(projected_points)
    return delaunay


def circumcenter(a, b, c):
    ac = c - a
    ab = b - a
    abXac = np.cross(ab, ac)
    numerator = np.cross(abXac, ab) * np.sum(ac**2) + np.cross(ac, abXac) * np.sum(
        ab**2
    )
    to_circumcenter = numerator / (2 * np.sum(abXac**2))
    return a + to_circumcenter


def generate_voronoi(points, delaunay):
    centers = []
    for simplex in delaunay.simplices:
        a, b, c = points[simplex]
        centers.append(circumcenter(a, b, c))
    return np.array(centers)


def plot(points, delaunay, voronoi_centers):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    ax.scatter(points[:, 0], points[:, 1], points[:, 2], color="blue", s=1)

    for simplex in delaunay.simplices:
        triangle = points[simplex]
        poly3d = Poly3DCollection([triangle], alpha=0.3, linewidths=1)
        poly3d.set_edgecolor("r")
        ax.add_collection3d(poly3d)

    ax.scatter(
        voronoi_centers[:, 0],
        voronoi_centers[:, 1],
        voronoi_centers[:, 2],
        color="green",
        s=5,
    )

    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    ax.view_init(elev=30, azim=30)  # Adjusting the view for better perspective
    plt.show()


def main(N=1000, jitter=0.0):
    points = generate_fibonacci_sphere(N, jitter)
    points_with_south_pole = add_south_pole(points)
    delaunay = generate_delaunay(points_with_south_pole)
    voronoi_centers = generate_voronoi(points_with_south_pole, delaunay)
    plot(points_with_south_pole, delaunay, voronoi_centers)


if __name__ == "__main__":
    main()
