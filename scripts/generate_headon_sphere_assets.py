from __future__ import annotations

import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
ASSET_DIR = ROOT / "assets" / "headon_spheres"
MODEL_DIR = ASSET_DIR / "models"
RATIO_ASSET_DIR = ROOT / "assets" / "headon_spheres_mass_ratio"
RATIO_MODEL_DIR = RATIO_ASSET_DIR / "models"


def write_uv_sphere_obj(path: Path, radius: float, slices: int, stacks: int) -> None:
    vertices = []
    faces = []

    for stack in range(stacks + 1):
        v = stack / stacks
        phi = math.pi * v
        sin_phi = math.sin(phi)
        cos_phi = math.cos(phi)
        for slice_idx in range(slices):
            u = slice_idx / slices
            theta = 2.0 * math.pi * u
            x = radius * sin_phi * math.cos(theta)
            y = radius * cos_phi
            z = radius * sin_phi * math.sin(theta)
            vertices.append((x, y, z))

    def vid(stack: int, slice_idx: int) -> int:
        return stack * slices + (slice_idx % slices)

    for stack in range(stacks):
        for slice_idx in range(slices):
            a = vid(stack, slice_idx)
            b = vid(stack, slice_idx + 1)
            c = vid(stack + 1, slice_idx)
            d = vid(stack + 1, slice_idx + 1)
            if stack != 0:
                faces.append((a + 1, c + 1, b + 1))
            if stack != stacks - 1:
                faces.append((b + 1, c + 1, d + 1))

    with path.open("w", encoding="ascii", newline="\n") as f:
        f.write("# Generated UV sphere for head-on equal sphere collision benchmark\n")
        for x, y, z in vertices:
            f.write(f"v {x:.9f} {y:.9f} {z:.9f}\n")
        for i, j, k in faces:
            f.write(f"f {i} {j} {k}\n")


def main() -> None:
    radius = 0.05
    slices = 64
    stacks = 32

    MODEL_DIR.mkdir(parents=True, exist_ok=True)
    RATIO_MODEL_DIR.mkdir(parents=True, exist_ok=True)
    write_uv_sphere_obj(MODEL_DIR / "ball_sphere.obj", radius, slices, stacks)
    write_uv_sphere_obj(RATIO_MODEL_DIR / "ball_sphere.obj", radius, slices, stacks)

    equal_mass_model_json = {
        "sphere_radius": radius,
        "sphere_a_density": 1000.0,
        "sphere_b_density": 1000.0,
        "sphere_mesh": "models/ball_sphere.obj",
        "sphere_a_init_pos": [-0.15, 0.0, 0.0],
        "sphere_b_init_pos": [0.15, 0.0, 0.0],
        "sphere_a_init_vel": [1.0, 0.0, 0.0],
        "sphere_b_init_vel": [0.0, 0.0, 0.0],
        "gravity_y": 0.0,
        "friction": 0.0,
        "restitution": 1.0,
        "step_size": 5.0e-4,
        "total_time": 0.5,
    }

    with (ASSET_DIR / "headon_spheres_model.json").open("w", encoding="utf-8", newline="\n") as f:
        json.dump(equal_mass_model_json, f, indent=2)
        f.write("\n")

    mass_ratio_model_json = {
        "sphere_radius": radius,
        "sphere_a_density": 1000.0,
        "sphere_b_density": 2000.0,
        "sphere_mesh": "models/ball_sphere.obj",
        "sphere_a_init_pos": [-0.15, 0.0, 0.0],
        "sphere_b_init_pos": [0.15, 0.0, 0.0],
        "sphere_a_init_vel": [0.0, 0.0, 0.0],
        "sphere_b_init_vel": [-1.0, 0.0, 0.0],
        "gravity_y": 0.0,
        "friction": 0.0,
        "restitution": 1.0,
        "step_size": 5.0e-4,
        "total_time": 0.5,
        "analytic_solution": {
            "description": "1D elastic head-on collision, equal radii, m_B = 2 m_A, B moving left with speed 1 m/s",
            "sphere_a_postimpact_vx": -4.0 / 3.0,
            "sphere_b_postimpact_vx": -1.0 / 3.0
        }
    }

    with (RATIO_ASSET_DIR / "headon_spheres_mass_ratio_model.json").open("w", encoding="utf-8", newline="\n") as f:
        json.dump(mass_ratio_model_json, f, indent=2)
        f.write("\n")


if __name__ == "__main__":
    main()
