import json
import math
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
ASSET_DIR = ROOT / "assets" / "rev_joint_clearance"
RAW_MODEL_DIR = ASSET_DIR / "models"
FINAL_MODEL_DIR = ASSET_DIR / "models"
RMD_PATH = ASSET_DIR / "rev_clearance_joint.rmd"


def read_obj_vertices_faces(path: Path):
    vertices = []
    faces = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith("v "):
            _, xs, ys, zs = line.split()
            vertices.append((float(xs), float(ys), float(zs)))
        elif line.startswith("f "):
            _, a, b, c = line.split()
            faces.append((int(a), int(b), int(c)))
    return vertices, faces


def write_obj(path: Path, vertices, faces):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as f:
        f.write("# normalized from RecurDyn rev_clearance_joint.rmd\n")
        for x, y, z in vertices:
            f.write(f"v {x:.12f} {y:.12f} {z:.12f}\n")
        for i, j, k in faces:
            f.write(f"f {i} {j} {k}\n")


def parse_rmd_float(text: str) -> float:
    s = text.strip().replace("D", "E").replace("d", "e")
    if s.endswith(("E", "e")):
        s += "0"
    return float(s)


def parse_markers(lines):
    markers = {}
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if not line.startswith("MARKER /"):
            i += 1
            continue

        marker_id = int(line.split("/")[1].strip())
        qp = None
        reuler = None
        i += 1
        while i < len(lines):
            stripped = lines[i].strip()
            if stripped.startswith(("MARKER /", "PART /", "GGEOM /", "JOINT /", "GGEOMCONTACT /", "!================================")):
                break
            if "QP =" in stripped:
                vals = stripped.split("=", 1)[1].split(",")
                qp = tuple(parse_rmd_float(v) for v in vals[:3])
            elif "REULER =" in stripped:
                vals = stripped.split("=", 1)[1].split(",")
                reuler = tuple(parse_rmd_float(v) for v in vals[:3])
            i += 1
        markers[marker_id] = {
            "qp": qp if qp is not None else (0.0, 0.0, 0.0),
            "reuler": reuler if reuler is not None else (0.0, 0.0, 0.0),
        }
    return markers


def parse_ggeom_rm(lines):
    out = {}
    i = 0
    current_name = None
    while i < len(lines):
        line = lines[i].strip()
        if not line.startswith("GGEOM /"):
            i += 1
            continue
        current_name = None
        rm = None
        i += 1
        while i < len(lines):
            stripped = lines[i].strip()
            if stripped.startswith(("GGEOM /", "PART /", "MARKER /", "JOINT /", "GGEOMCONTACT /", "!================================")):
                break
            if "NAME =" in stripped:
                match = re.search(r"NAME\s*=\s*'([^']+)'", stripped)
                if match:
                    current_name = match.group(1)
            elif "RM =" in stripped:
                rm = int(stripped.split("=")[1].strip())
            i += 1
        if current_name and rm is not None:
            out[current_name] = rm
    return out


def matmul3(a, b):
    return tuple(
        tuple(sum(a[r][k] * b[k][c] for k in range(3)) for c in range(3))
        for r in range(3)
    )


def matvec3(a, v):
    return (
        a[0][0] * v[0] + a[0][1] * v[1] + a[0][2] * v[2],
        a[1][0] * v[0] + a[1][1] * v[1] + a[1][2] * v[2],
        a[2][0] * v[0] + a[2][1] * v[1] + a[2][2] * v[2],
    )


def transpose3(a):
    return (
        (a[0][0], a[1][0], a[2][0]),
        (a[0][1], a[1][1], a[2][1]),
        (a[0][2], a[1][2], a[2][2]),
    )


def reuler_to_matrix(reuler):
    ax, ay, az = reuler
    cx, sx = math.cos(ax), math.sin(ax)
    cy, sy = math.cos(ay), math.sin(ay)
    cz, sz = math.cos(az), math.sin(az)
    rx = ((1.0, 0.0, 0.0), (0.0, cx, -sx), (0.0, sx, cx))
    ry = ((cy, 0.0, sy), (0.0, 1.0, 0.0), (-sy, 0.0, cy))
    rz = ((cz, -sz, 0.0), (sz, cz, 0.0), (0.0, 0.0, 1.0))
    # RecurDyn surface RM markers on this model match the part-frame transform p_part = R^T * p_raw + q,
    # with R built from intrinsic XYZ REULER rotations.
    return transpose3(matmul3(matmul3(rx, ry), rz))


def apply_surface_marker(vertices, qp, reuler):
    r = reuler_to_matrix(reuler)
    out = []
    for p in vertices:
        rp = matvec3(r, p)
        out.append((rp[0] + qp[0], rp[1] + qp[1], rp[2] + qp[2]))
    return out


def main():
    raw_body1 = RAW_MODEL_DIR / "body1_subtract1.obj"
    raw_body3 = RAW_MODEL_DIR / "body3_cylinder1.obj"
    if not raw_body1.exists() or not raw_body3.exists():
        raise SystemExit(
            "Raw OBJ files are missing. Run scripts/extract_recurdyn_rmd_ggeom_to_obj.py first."
        )

    body1_vertices, body1_faces = read_obj_vertices_faces(raw_body1)
    body3_vertices, body3_faces = read_obj_vertices_faces(raw_body3)

    lines = RMD_PATH.read_text(encoding="utf-8-sig").splitlines()
    markers = parse_markers(lines)
    ggeom_rm = parse_ggeom_rm(lines)

    body1_marker = markers[ggeom_rm["##GSURFACE##_Body1.Subtract1"]]
    body3_marker = markers[ggeom_rm["##GSURFACE##_Body3.Cylinder1"]]

    body1_final = apply_surface_marker(body1_vertices, body1_marker["qp"], body1_marker["reuler"])
    body3_final = apply_surface_marker(body3_vertices, body3_marker["qp"], body3_marker["reuler"])

    body1_out = FINAL_MODEL_DIR / "body1_subtract1_centered.obj"
    body3_out = FINAL_MODEL_DIR / "body3_cylinder1_centered.obj"
    write_obj(body1_out, body1_final, body1_faces)
    write_obj(body3_out, body3_final, body3_faces)

    model_json = {
        "body1_mesh_path": "assets/rev_joint_clearance/models/body1_subtract1_centered.obj",
        "body3_mesh_path": "assets/rev_joint_clearance/models/body3_cylinder1_centered.obj",
        "surface_rm_note": "Meshes are restored to the part frame using the RecurDyn surface RM marker transform.",
        "body1_init_pos": [0.0, 0.0, 0.0],
        "body3_init_pos": [0.0, 0.0, 0.0],
        "body2_cm_offset": [1.64388630568497e-06, 8.35588256883334e-05, 2.65983598394373],
        "body3_mass": 12984.280977103,
        "body3_inertia_xx": [1314.65844893167, 7971.80750823384, 7971.80750823384],
        "body3_inertia_xy": [0.0, 3.02128692434659e-13, 0.0],
        "body2_mass": 271.0,
        "body2_inertia_xx": [202.0, 359.0, 203.0],
        "body2_inertia_xy": [0.244945783962547, -0.0295154838928402, -1.39749450062735],
        "gravity_y": -9.80665,
        "friction": 0.0,
        "restitution": 0.0,
        "contact_compliance": 1.0e-9,
        "contact_compliance_t": 1.0e-9,
        "contact_damping_f": 1.0e-5,
        "collision_envelope": 2.0e-3,
        "step_size": 0.001,
        "total_time": 3.0,
        "dynamics_substeps": 1,
        "env_prefix": "SPCC_REVCLR",
        "body2_reference_csv": "assets/rev_joint_clearance/data/body2.csv",
        "body3_reference_csv": "assets/rev_joint_clearance/data/body3.csv",
        "body2_ideal_csv": "assets/rev_joint_clearance/data/body2_ideal.csv",
    }
    model_path = ASSET_DIR / "rev_joint_clearance_model.json"
    model_path.write_text(json.dumps(model_json, indent=2), encoding="utf-8", newline="\n")

    print(f"Wrote {body1_out}")
    print(f"Wrote {body3_out}")
    print(f"Wrote {model_path}")


if __name__ == "__main__":
    main()
