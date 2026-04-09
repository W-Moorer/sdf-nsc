import argparse
import json
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


GGEOM_TARGETS = {
    "##GSURFACE##_Body1.Subtract1": "body1_subtract1.obj",
    "##GSURFACE##_Body3.Cylinder1": "body3_cylinder1.obj",
}


def parse_ggeom_blocks(lines):
    blocks = []
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if not line.startswith("GGEOM /"):
            i += 1
            continue

        block = {
            "id": int(line.split("/")[1].strip()),
            "name": None,
            "no_patch": None,
            "no_node": None,
            "patches": [],
            "nodes": [],
        }
        i += 1

        while i < len(lines):
            raw = lines[i].rstrip("\n")
            stripped = raw.strip()

            if stripped.startswith("GGEOM /") or stripped.startswith("PART /") or stripped.startswith("MARKER /") or stripped.startswith("JOINT /") or stripped.startswith("GGEOMCONTACT /") or stripped.startswith("!================================"):
                break

            if "NAME =" in stripped:
                match = re.search(r"NAME\s*=\s*'([^']+)'", stripped)
                if match:
                    block["name"] = match.group(1)
            elif "NO_PATCH" in stripped:
                block["no_patch"] = int(stripped.split("=")[1].strip())
            elif "NO_NODE" in stripped:
                block["no_node"] = int(stripped.split("=")[1].strip())
            elif stripped.startswith(", PATCHES"):
                i += 1
                while i < len(lines) and len(block["patches"]) < block["no_patch"]:
                    patch_line = lines[i].strip().lstrip(",").strip()
                    if not patch_line:
                        i += 1
                        continue
                    parts = [p.strip() for p in patch_line.split(",")]
                    if len(parts) >= 4:
                        block["patches"].append(tuple(int(v) for v in parts[1:4]))
                    i += 1
                continue
            elif stripped.startswith(", NODES"):
                i += 1
                while i < len(lines) and len(block["nodes"]) < block["no_node"]:
                    node_line = lines[i].strip().lstrip(",").strip()
                    if not node_line:
                        i += 1
                        continue
                    parts = [p.strip() for p in node_line.split(",")]
                    if len(parts) >= 3:
                        block["nodes"].append(tuple(float(v) for v in parts[:3]))
                    i += 1
                continue

            i += 1

        if block["name"] and block["patches"] and block["nodes"]:
            blocks.append(block)

    return blocks


def write_obj(path, nodes, faces):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as f:
        f.write("# extracted from RecurDyn RMD\n")
        for x, y, z in nodes:
            f.write(f"v {x:.12f} {y:.12f} {z:.12f}\n")
        for i, j, k in faces:
            f.write(f"f {i} {j} {k}\n")


def bbox(nodes):
    xs = [p[0] for p in nodes]
    ys = [p[1] for p in nodes]
    zs = [p[2] for p in nodes]
    return {
        "min": [min(xs), min(ys), min(zs)],
        "max": [max(xs), max(ys), max(zs)],
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", default=str(ROOT / "assets" / "rev_joint_clearance" / "rev_clearance_joint.rmd"))
    parser.add_argument("--output-dir", default=str(ROOT / "assets" / "rev_joint_clearance" / "models"))
    args = parser.parse_args()

    input_path = Path(args.input)
    output_dir = Path(args.output_dir)
    lines = input_path.read_text(encoding="utf-8-sig").splitlines()
    blocks = parse_ggeom_blocks(lines)
    by_name = {b["name"]: b for b in blocks}

    summary = {}
    for name, out_name in GGEOM_TARGETS.items():
        if name not in by_name:
            raise RuntimeError(f"Missing GGEOM block: {name}")
        block = by_name[name]
        out_path = output_dir / out_name
        write_obj(out_path, block["nodes"], block["patches"])
        summary[name] = {
            "obj": str(out_path.relative_to(ROOT)).replace("\\", "/"),
            "id": block["id"],
            "no_patch": block["no_patch"],
            "no_node": block["no_node"],
            "bbox": bbox(block["nodes"]),
        }
        print(f"Wrote {out_path}")

    summary_path = output_dir / "extract_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8", newline="\n")
    print(f"Wrote {summary_path}")


if __name__ == "__main__":
    main()
