# planet.py

from typing import List

from ursina import *

from ursina.terrain_face import TerrainFace


class Planet(Entity):
    def __init__(self, resolution=10, **kwargs):
        super().__init__(**kwargs)
        self.resolution = resolution
        self.terrain_faces: List[TerrainFace] = []
        self.meshes: List[Mesh] = [Mesh() for _ in range(6)]
        self.initialize()

    def initialize(self):
        directions = [
            Vec3(0, 1, 0),
            Vec3(0, -1, 0),
            Vec3(-1, 0, 0),
            Vec3(1, 0, 0),
            Vec3(0, 0, 1),
            Vec3(0, 0, -1),
        ]

        for i, direction in enumerate(directions):
            mesh = Mesh()
            self.terrain_faces.append(
                TerrainFace(mesh, self.resolution, local_up=direction)
            )
            self.meshes[i] = mesh
            face_entity = Entity(model=self.meshes[i], parent=self, double_sided=True)

        for i, mesh in enumerate(self.meshes):
            self.meshes[i] = self.terrain_faces[i].construct_mesh()
            face_entity = Entity(model=self.meshes[i], parent=self, double_sided=True)

    def generate_mesh(self):
        for i in range(6):
            self.meshes[i] = self.terrain_faces[i].construct_mesh()
