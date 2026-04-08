import json
import math
from pathlib import Path

from generate_onset_stress_b_assets import make_cylinder, write_obj


ROOT = Path(__file__).resolve().parents[1]
ASSET_DIR = ROOT / "assets" / "onset_stress_b_oblique"
MODEL_DIR = ASSET_DIR / "models"


def main():
    params = {
        "cam_radius": 0.030,
        "cam_eccentricity": 0.006,
        "cam_thickness": 0.020,
        "roller_radius": 0.010,
        "roller_thickness": 0.018,
        "cam_segments": 160,
        "roller_segments": 128,
        "phase": math.pi,
        "cam_init_pos": [0.0, 0.0, 0.0],
        "motor_joint_pos": [0.0, 0.0, 0.0],
        "gravity_y": 0.0,
        "density": 7800.0,
        "friction": 0.0,
        "restitution": 0.0,
        "motor_speed": -2.0,
        "step_size": 0.001,
        "total_time": 0.45,
        "target_onset_time": 0.15,
        "dynamics_substeps": 32,
        "follower_x_offset": 0.012,
        "variant": "oblique_onset",
    }

    cam_half_t = 0.5 * params["cam_thickness"]
    roller_half_t = 0.5 * params["roller_thickness"]
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
    roller_vertices, roller_faces = make_cylinder(
        params["roller_radius"],
        roller_half_t,
        params["roller_segments"],
        center_x=0.0,
        center_y=0.0,
        center_z=0.0,
    )

    onset_theta = params["phase"] + params["motor_speed"] * params["target_onset_time"]
    reach = params["cam_radius"] + params["roller_radius"]
    follower_x = params["follower_x_offset"]
    lateral_gap = follower_x - params["cam_eccentricity"] * math.cos(onset_theta)
    onset_y = params["cam_eccentricity"] * math.sin(onset_theta) + math.sqrt(
        max(reach * reach - lateral_gap * lateral_gap, 0.0)
    )
    params["follower_init_pos"] = [follower_x, onset_y, 0.0]
    params["preload_anchor_pos"] = [follower_x, 0.0, 0.0]
    params["preload_rest_length"] = onset_y
    params["preload_stiffness"] = 670.0
    params["preload_damping"] = 47.0

    write_obj(MODEL_DIR / "onset_cam.obj", cam_vertices, cam_faces)
    write_obj(MODEL_DIR / "roller_follower.obj", roller_vertices, roller_faces)

    metadata_path = ASSET_DIR / "onset_stress_b_oblique_model.json"
    metadata_path.parent.mkdir(parents=True, exist_ok=True)
    with metadata_path.open("w", encoding="utf-8", newline="\n") as f:
        json.dump(params, f, indent=2)
        f.write("\n")

    print(f"Wrote {MODEL_DIR / 'onset_cam.obj'}")
    print(f"Wrote {MODEL_DIR / 'roller_follower.obj'}")
    print(f"Wrote {metadata_path}")


if __name__ == "__main__":
    main()
