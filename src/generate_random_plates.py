import pygplates
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import numpy as np
from shapely.geometry import Polygon


def generate_random_plate(num_points: int) -> Polygon:
    """Generate random spherical polygons representing tectonic plates."""
    points = []
    for _ in range(num_points):
        lat = np.random.uniform(-90, 90)
        lon = np.random.uniform(-180, 180)
        points.append((lon, lat))  # Create tuples for longitude and latitude
    polygon = Polygon(points)  # Create a Shapely Polygon
    return polygon


def main(number_of_plates=5, points_per_plate=10):
    """Generate a map with random tectonic plates."""
    fig = plt.figure(figsize=(10, 5))
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.coastlines()

    # Randomly generate tectonic plates
    for _ in range(number_of_plates):
        plate = generate_random_plate(points_per_plate)
        # Plot the plate
        ax.add_geometries(
            [plate],
            ccrs.PlateCarree(),
            edgecolor="black",
            facecolor="none",
        )

    plt.savefig("tectonic_plates.png")  # Save the plot as a PNG file


if __name__ == "__main__":
    main()
