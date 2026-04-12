#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
from pathlib import Path


CASE_PRESETS = {
    "eccentric_roller": {
        "binary": "baseline_eccentric_roller_nsc",
        "extra_args": [],
        "description": "Analytical eccentric disk / roller follower benchmark.",
    },
    "onset_stress_b": {
        "binary": "baseline_onset_stress_b_nsc",
        "extra_args": [],
        "description": "Onset Stress Benchmark Version B.",
    },
    "cam_small_dt": {
        "binary": "baseline_cam_nsc",
        "extra_args": ["--dt", "0.001"],
        "description": "Industrial cam benchmark with the smaller paper timestep.",
    },
    "cam_large_dt": {
        "binary": "baseline_cam_nsc",
        "extra_args": ["--dt", "0.005"],
        "description": "Industrial cam benchmark with the larger paper timestep.",
    },
    "twin_gear": {
        "binary": "baseline_simple_gear_nsc",
        "extra_args": [],
        "description": "Twin gear meshing benchmark.",
    },
}


def parse_case_list(value: str) -> list[str]:
    raw_items = [item.strip() for item in value.split(",") if item.strip()]
    if not raw_items:
        raise argparse.ArgumentTypeError("case list must not be empty")
    if len(raw_items) == 1 and raw_items[0] == "all":
        return list(CASE_PRESETS.keys())

    unknown = [item for item in raw_items if item not in CASE_PRESETS]
    if unknown:
        raise argparse.ArgumentTypeError(
            f"unknown case preset(s): {', '.join(unknown)}; choose from: {', '.join(CASE_PRESETS)} or all"
        )
    return raw_items


def resolve_binary(binary_dir: Path, binary_name: str) -> Path:
    candidate = (binary_dir / binary_name).resolve()
    if candidate.exists():
        return candidate

    exe_candidate = (binary_dir / f"{binary_name}.exe").resolve()
    if exe_candidate.exists():
        return exe_candidate

    raise FileNotFoundError(f"Binary not found: {candidate} or {exe_candidate}")


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Export VTK frame sequences plus PVD manifests for the follow-up paper benchmarks."
    )
    parser.add_argument(
        "--cases",
        type=parse_case_list,
        default=list(CASE_PRESETS.keys()),
        help="Comma-separated case preset list, or 'all'.",
    )
    parser.add_argument(
        "--output-root",
        default="data/outputs/followup_case_vtk_exports",
        help="Root directory for exported case folders.",
    )
    parser.add_argument(
        "--binary-dir",
        default="build",
        help="Directory containing the compiled benchmark binaries.",
    )
    parser.add_argument(
        "--contact-algorithm",
        default="sdf_2nd",
        help="Contact algorithm passed to the benchmark apps.",
    )
    parser.add_argument("--vtk-stride", default="1", help="VTK export stride.")
    parser.add_argument("--T", default=None, help="Optional total-time override for every selected case.")
    parser.add_argument("--speed", default=None, help="Optional motor speed override for every selected case.")
    parser.add_argument(
        "--clean",
        action="store_true",
        help="Delete the output root before running exports.",
    )
    args = parser.parse_args()

    repo_root = Path(__file__).resolve().parents[1]
    output_root = (repo_root / args.output_root).resolve()
    binary_dir = (repo_root / args.binary_dir).resolve()

    if args.clean and output_root.exists():
        shutil.rmtree(output_root)
    output_root.mkdir(parents=True, exist_ok=True)

    runs: list[dict[str, str]] = []
    for case_name in args.cases:
        preset = CASE_PRESETS[case_name]
        binary_path = resolve_binary(binary_dir, preset["binary"])

        run_dir = output_root / case_name / args.contact_algorithm
        run_dir.mkdir(parents=True, exist_ok=True)
        vtk_dir = run_dir / "frames"
        vtk_dir.mkdir(parents=True, exist_ok=True)
        csv_path = run_dir / "result.csv"

        command = [
            str(binary_path),
            *preset["extra_args"],
            "--contact-algorithm",
            args.contact_algorithm,
            "--output",
            str(csv_path.relative_to(repo_root)),
            "--vtk-dir",
            str(vtk_dir.relative_to(repo_root)),
            "--vtk-stride",
            args.vtk_stride,
        ]
        if args.T is not None:
            command.extend(["--T", args.T])
        if args.speed is not None:
            command.extend(["--speed", args.speed])

        print(f"[RUN] case={case_name} algorithm={args.contact_algorithm}")
        print("      " + " ".join(command))
        subprocess.run(command, cwd=repo_root, check=True)

        pvd_files = sorted(vtk_dir.glob("*.pvd"))
        vtk_series_files = sorted(vtk_dir.glob("*.vtk.series"))
        runs.append(
            {
                "case": case_name,
                "description": preset["description"],
                "binary": str(binary_path.relative_to(repo_root)),
                "output_csv": str(csv_path.relative_to(repo_root)),
                "vtk_dir": str(vtk_dir.relative_to(repo_root)),
                "pvd": str(pvd_files[0].relative_to(repo_root)) if pvd_files else "",
                "vtk_series": str(vtk_series_files[0].relative_to(repo_root)) if vtk_series_files else "",
            }
        )

    manifest_path = output_root / "manifest.json"
    manifest_path.write_text(
        json.dumps(
            {
                "contact_algorithm": args.contact_algorithm,
                "vtk_stride": args.vtk_stride,
                "runs": runs,
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    print(f"[DONE] Wrote manifest: {manifest_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
