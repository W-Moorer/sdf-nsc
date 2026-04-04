#!/usr/bin/env python3

import argparse
import json
import shutil
import subprocess
from pathlib import Path


ALGORITHMS = ("mesh", "sdf_1st", "sdf_2nd")


def parse_dt_list(value: str) -> list[str]:
    parts = [item.strip() for item in value.split(",")]
    out = [item for item in parts if item]
    if not out:
        raise argparse.ArgumentTypeError("dt list must not be empty")
    return out


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Export full VTK frame sequences for cam contact algorithm sweeps."
    )
    parser.add_argument(
        "--output-root",
        default="data/outputs/cam_vtk_exports",
        help="Root directory for algorithm/dt subfolders.",
    )
    parser.add_argument(
        "--dt-list",
        type=parse_dt_list,
        default=["0.005", "0.01"],
        help="Comma-separated list of step sizes to export.",
    )
    parser.add_argument("--T", default="0.22", help="Total simulation time passed to baseline_cam_nsc.")
    parser.add_argument("--speed", default=None, help="Optional motor speed override.")
    parser.add_argument("--vtk-stride", default="1", help="VTK export stride.")
    parser.add_argument(
        "--binary",
        default="./_build/project/baseline_cam_nsc",
        help="Path to the baseline_cam_nsc binary, relative to the repo root.",
    )
    parser.add_argument(
        "--clean",
        action="store_true",
        help="Delete the output root before running the sweep.",
    )
    args = parser.parse_args()

    repo_root = Path(__file__).resolve().parents[1]
    output_root = (repo_root / args.output_root).resolve()
    binary_path = (repo_root / args.binary).resolve()

    if not binary_path.exists():
        raise SystemExit(f"Binary not found: {binary_path}")

    if args.clean and output_root.exists():
        shutil.rmtree(output_root)
    output_root.mkdir(parents=True, exist_ok=True)

    runs = []
    for algorithm in ALGORITHMS:
        for dt in args.dt_list:
            dt_folder = f"dt_{dt}"
            run_dir = output_root / algorithm / dt_folder
            run_dir.mkdir(parents=True, exist_ok=True)

            csv_path = run_dir / "result.csv"
            vtk_dir = run_dir / "frames"
            vtk_dir.mkdir(parents=True, exist_ok=True)

            command = [
                str(binary_path),
                "--contact-algorithm",
                algorithm,
                "--dt",
                dt,
                "--T",
                args.T,
                "--output",
                str(csv_path.relative_to(repo_root)),
                "--vtk-dir",
                str(vtk_dir.relative_to(repo_root)),
                "--vtk-stride",
                args.vtk_stride,
            ]
            if args.speed is not None:
                command.extend(["--speed", args.speed])

            print(f"[RUN] algorithm={algorithm} dt={dt}")
            print("      " + " ".join(command))
            subprocess.run(command, cwd=repo_root, check=True)

            runs.append(
                {
                    "algorithm": algorithm,
                    "dt": dt,
                    "output_csv": str(csv_path.relative_to(repo_root)),
                    "vtk_dir": str(vtk_dir.relative_to(repo_root)),
                }
            )

    manifest_path = output_root / "manifest.json"
    manifest_path.write_text(json.dumps({"runs": runs}, indent=2), encoding="utf-8")
    print(f"[DONE] Wrote manifest: {manifest_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
