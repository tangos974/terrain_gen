# main.py
from ursina import *

from planet import Planet


def toggle_points():
    global point_mode
    point_mode = not point_mode
    mode = "point" if point_mode else "triangle"
    print(f"Toggling to mode: {mode}")  # Print the mode to debug
    for mesh in planet.meshes:
        print(f"Mesh: {mesh}")  # Print the mesh to debug
        if mesh:
            mesh.mode = mode
        else:
            print("Error: Mesh is None")  # Catch None mesh issues


if __name__ == "__main__":
    app = Ursina(fullscreen=True)  # Keeping the Ursina app instance named `app`

    planet = Planet(resolution=10)
    EditorCamera()

    point_mode = False
    toggle_button = Button(
        text="Toggle Points",
        parent=planet,
        color=color.azure,
        scale=(0.15, 0.05),
        position=(0, -0.4),
    )
    # toggle_button.on_click = toggle_points

    app.run()
