import json
import math
from pathlib import Path

from generate_onset_stress_b_assets import write_obj, make_cylinder


ROOT = Path(__file__).resolve().parents[1]
ASSET_DIR = ROOT / "assets" / "onset_stress_b_oblique_shoe"
MODEL_DIR = ASSET_DIR / "models"


def sample_rounded_rect_contour(half_width, half_height, radius, arc_segments):
    if radius <= 0.0 or radius >= half_width or radius >= half_height:
        raise ValueError("radius must be positive and smaller than both half extents")

    pts = []

    def add_arc(cx, cy, start_deg, end_deg):
        for i in range(arc_segments + 1):
            t = i / arc_segments
            ang = math.radians(start_deg + (end_deg - start_deg) * t)
            pts.append((cx + radius * math.cos(ang), cy + radius * math.sin(ang)))

    top_y = half_height
    bot_y = -half_height
    left_x = -half_width
    right_x = half_width
    flat_x = half_width - radius
    flat_y = half_height - radius

    pts.append((-flat_x, top_y))
    pts.append((flat_x, top_y))
    add_arc(flat_x, flat_y, 90.0, 0.0)
    pts.append((right_x, -flat_y))
    add_arc(flat_x, -flat_y, 0.0, -90.0)
    pts.append((-flat_x, bot_y))
    add_arc(-flat_x, -flat_y, -90.0, -180.0)
    pts.append((left_x, flat_y))
    add_arc(-flat_x, flat_y, 180.0, 90.0)

    cleaned = []
    for p in pts:
        if not cleaned or abs(cleaned[-1][0] - p[0]) > 1.0e-12 or abs(cleaned[-1][1] - p[1]) > 1.0e-12:
            cleaned.append(p)
    if abs(cleaned[0][0] - cleaned[-1][0]) < 1.0e-12 and abs(cleaned[0][1] - cleaned[-1][1]) < 1.0e-12:
        cleaned.pop()
    return cleaned


def extrude_polygon(points_xy, half_thickness):
    vertices = []
    faces = []
    n = len(points_xy)
    for z_sign in (-1.0, 1.0):
        z = z_sign * half_thickness
        for x, y in points_xy:
            vertices.append((x, y, z))

    bottom_center_index = len(vertices) + 1
    vertices.append((0.0, 0.0, -half_thickness))
    top_center_index = len(vertices) + 1
    vertices.append((0.0, 0.0, half_thickness))

    bottom_start = 1
    top_start = 1 + n
    for i in range(n):
        b0 = bottom_start + i
        b1 = bottom_start + ((i + 1) % n)
        t0 = top_start + i
        t1 = top_start + ((i + 1) % n)
        faces.append((b0, b1, t1))
        faces.append((b0, t1, t0))
        faces.append((bottom_center_index, b1, b0))
        faces.append((top_center_index, t0, t1))

    return vertices, faces


def rounded_rect_top_support(dx, half_width, half_height, radius):
    ax = abs(dx)
    if ax > half_width:
        return None
    flat_limit = half_width - radius
    if ax <= flat_limit:
        return half_height
    local = ax - flat_limit
    return (half_height - radius) + math.sqrt(max(radius * radius - local * local, 0.0))


def shoe_support_contact_y(cam_radius, cam_eccentricity, phase, motor_speed, t, follower_x, half_width, half_height, radius):
    theta = phase + motor_speed * t
    cx = cam_eccentricity * math.cos(theta)
    cy = cam_eccentricity * math.sin(theta)
    dx = follower_x - cx
    top = rounded_rect_top_support(dx, half_width, half_height, radius)
    if top is None:
        return None
    return cy + cam_radius + top


def main():
    params = {
        "cam_radius": 0.030,
        "cam_eccentricity": 0.006,
        "cam_thickness": 0.020,
        "shoe_half_width": 0.020,
        "shoe_half_height": 0.010,
        "shoe_corner_radius": 0.004,
        "shoe_thickness": 0.018,
        "cam_segments": 160,
        "shoe_arc_segments": 48,
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
        "variant": "oblique_shoe",
    }

    cam_half_t = 0.5 * params["cam_thickness"]
    shoe_half_t = 0.5 * params["shoe_thickness"]
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
    contour = sample_rounded_rect_contour(
        params["shoe_half_width"],
        params["shoe_half_height"],
        params["shoe_corner_radius"],
        params["shoe_arc_segments"],
    )
    shoe_vertices, shoe_faces = extrude_polygon(contour, shoe_half_t)

    follower_x = params["follower_x_offset"]
    onset_y = shoe_support_contact_y(
        params["cam_radius"],
        params["cam_eccentricity"],
        params["phase"],
        params["motor_speed"],
        params["target_onset_time"],
        follower_x,
        params["shoe_half_width"],
        params["shoe_half_height"],
        params["shoe_corner_radius"],
    )
    if onset_y is None:
        raise RuntimeError("target onset lies outside the shoe top-support domain; adjust offset or width")
    params["follower_init_pos"] = [follower_x, onset_y, 0.0]
    params["preload_anchor_pos"] = [follower_x, 0.0, 0.0]
    params["preload_rest_length"] = onset_y
    params["preload_stiffness"] = 670.0
    params["preload_damping"] = 47.0

    write_obj(MODEL_DIR / "onset_cam.obj", cam_vertices, cam_faces)
    write_obj(MODEL_DIR / "shoe_follower.obj", shoe_vertices, shoe_faces)

    metadata_path = ASSET_DIR / "onset_stress_b_oblique_shoe_model.json"
    metadata_path.parent.mkdir(parents=True, exist_ok=True)
    with metadata_path.open("w", encoding="utf-8", newline="\n") as f:
        json.dump(params, f, indent=2)
        f.write("\n")

    print(f"Wrote {MODEL_DIR / 'onset_cam.obj'}")
    print(f"Wrote {MODEL_DIR / 'shoe_follower.obj'}")
    print(f"Wrote {metadata_path}")


if __name__ == "__main__":
    main()
