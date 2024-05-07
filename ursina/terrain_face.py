# terrain_face.py

from typing import List, Tuple

from ursina import *


class TerrainFace:
    def __init__(self, mesh: Mesh, resolution: int, local_up: Vec3) -> None:
        self.mesh: Mesh = mesh
        self.resolution: int = resolution
        self.local_up: Vec3 = local_up

        self.axis_a: Vec3 = Vec3(local_up.y, local_up.z, local_up.x)
        self.axis_b: Vec3 = Vec3.cross(local_up, self.axis_a)

    def construct_mesh(self) -> None:
        vertices: List[Vec3] = []
        triangles: List[int] = []

        for y in range(self.resolution):
            for x in range(self.resolution):
                i: int = x + y * self.resolution
                percent: Vec2 = Vec2(x, y) / (self.resolution - 1)
                point_on_unit_cube: Vec3 = (
                    self.local_up
                    + (self.axis_a * (percent.x - 0.5) * 2)
                    + (self.axis_b * (percent.y - 0.5) * 2)
                )
                point_on_unit_sphere: Vec3 = point_on_unit_cube.normalized()
                vertices.append(point_on_unit_sphere)

                if x != self.resolution - 1 and y != self.resolution - 1:
                    triangles.extend(
                        [
                            i,
                            i + self.resolution + 1,
                            i + self.resolution,
                            i,
                            i + 1,
                            i + self.resolution + 1,
                        ]
                    )

        self.mesh.vertices = vertices
        self.mesh.triangles = triangles
        self.mesh.uvs = [Vec2(0, 0) for _ in vertices]
        self.mesh.generate_normals()
