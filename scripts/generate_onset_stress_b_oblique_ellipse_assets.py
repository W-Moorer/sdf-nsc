import json
import math
from pathlib import Path

from generate_onset_stress_b_assets import make_cylinder, write_obj


ROOT = Path(__file__).resolve().parents[1]
ASSET_DIR = ROOT / "assets" / "onset_stress_b_oblique_ellipse"
MODEL_DIR = ASSET_DIR / "models"


def make_extruded_ellipse(a_x, a_y, half_thickness, segments, center_x=0.0, center_y=0.0):
    vertices = []
    faces = []

    for z_sign in (-1.0, 1.0):
        z = z_sign * half_thickness
        for s in range(segments):
            theta = 2.0 * math.pi * s / segments
            x = center_x + a_x * math.cos(theta)
            y = center_y + a_y * math.sin(theta)
            vertices.append((x, y, z))

    bottom_center_index = len(vertices) + 1
    vertices.append((center_x, center_y, -half_thickness))
    top_center_index = len(vertices) + 1
    vertices.append((center_x, center_y, half_thickness))

    bottom_start = 1
    top_start = 1 + segments

    for s in range(segments):
        b0 = bottom_start + s
        b1 = bottom_start + ((s + 1) % segments)
        t0 = top_start + s
        t1 = top_start + ((s + 1) % segments)
        faces.append((b0, b1, t1))
        faces.append((b0, t1, t0))
        faces.append((bottom_center_index, b1, b0))
        faces.append((top_center_index, t0, t1))

    return vertices, faces


def ellipse_support_contact_y(cam_radius, cam_eccentricity, phase, motor_speed, t, follower_x, ellipse_ax, ellipse_ay):
    theta = phase + motor_speed * t
    cx = cam_eccentricity * math.cos(theta)
    cy = cam_eccentricity * math.sin(theta)
    dx = follower_x - cx
    if abs(dx) > ellipse_ax:
        return None
    s = max(1.0 - (dx * dx) / (ellipse_ax * ellipse_ax), 0.0)
    return cy + cam_radius + ellipse_ay * math.sqrt(s)


def main():
    params = {
        "cam_radius": 0.030,
        "cam_eccentricity": 0.006,
        "cam_thickness": 0.020,
        "follower_ellipse_ax": 0.020,
        "follower_ellipse_ay": 0.008,
        "follower_thickness": 0.018,
        "cam_segments": 160,
        "follower_segments": 160,
        "phase": math.pi,
        "cam_init_pos": [0.0, 0.0, 0.0],
        "motor_joint_pos": [0.0, 0.0, 0.0],
        "gravity_y": 0.0,
        "density": 7800.0,
        "friction": 0.0,
        "restitution": 0.0,
        "motor_speed": -2.0,
        "step_size": 0.001,
        "total_time": 0.25,
        "target_onset_time": 0.15,
        "dynamics_substeps": 8,
        "follower_x_offset": 0.012,
        "variant": "oblique_ellipse",
    }

    cam_half_t = 0.5 * params["cam_thickness"]
    follower_half_t = 0.5 * params["follower_thickness"]
    center_x = params["cam_eccentricity"] * math.cos(params["phase"])
    center_y = params["cam_eccentricity"] * math.sin(params["phase"])

    cam_vertices, cam_faces = make_cylinder(
        params["cam_radius"],
        cam_half_t,
        params["cam_segments"],
        center_x=center_x,
        center_y=center_y,
        center_z=0.0,
    )
    follower_vertices, follower_faces = make_extruded_ellipse(
        params["follower_ellipse_ax"],
        params["follower_ellipse_ay"],
        follower_half_t,
        params["follower_segments"],
        center_x=0.0,
        center_y=0.0,
    )

    follower_x = params["follower_x_offset"]
    onset_y = ellipse_support_contact_y(
        params["cam_radius"],
        params["cam_eccentricity"],
        params["phase"],
        params["motor_speed"],
        params["target_onset_time"],
        follower_x,
        params["follower_ellipse_ax"],
        params["follower_ellipse_ay"],
    )
    params["follower_init_pos"] = [follower_x, onset_y, 0.0]
    params["preload_anchor_pos"] = [follower_x, 0.0, 0.0]
    params["preload_rest_length"] = onset_y
    params["preload_stiffness"] = 670.0
    params["preload_damping"] = 47.0

    write_obj(MODEL_DIR / "onset_cam.obj", cam_vertices, cam_faces)
    write_obj(MODEL_DIR / "ellipse_follower.obj", follower_vertices, follower_faces)

    metadata_path = ASSET_DIR / "onset_stress_b_oblique_ellipse_model.json"
    metadata_path.parent.mkdir(parents=True, exist_ok=True)
    with metadata_path.open("w", encoding="utf-8", newline="\n") as f:
        json.dump(params, f, indent=2)
        f.write("\n")

    print(f"Wrote {MODEL_DIR / 'onset_cam.obj'}")
    print(f"Wrote {MODEL_DIR / 'ellipse_follower.obj'}")
    print(f"Wrote {metadata_path}")


if __name__ == "__main__":
    main()
