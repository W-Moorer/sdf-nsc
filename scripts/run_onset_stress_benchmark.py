import csv
import html
import json
import math
import re
import subprocess
import sys
import time
from pathlib import Path

import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parents[1]
ASSET_JSON = ROOT / "assets" / "onset_stress" / "onset_stress_model.json"
OUTPUT_DIR = ROOT / "data" / "outputs"
DOCS_DIR = ROOT / "docs"
FIG_DIR = DOCS_DIR / "figures"
REPORT_PATH = DOCS_DIR / "onset_stress_report.html"


def to_wsl_path(path: Path) -> str:
    drive = path.drive[:-1].lower()
    tail = path.as_posix()[2:]
    return f"/mnt/{drive}{tail}"


def ensure_assets():
    subprocess.run([sys.executable, str(ROOT / "scripts" / "generate_onset_stress_assets.py")], check=True)


def load_csv(path: Path):
    rows = []
    with path.open("r", encoding="utf-8") as f:
        reader = csv.reader(f)
        next(reader)
        for row in reader:
            if not row:
                continue
            rows.append({"t": float(row[0]), "pos": float(row[10]), "vel": float(row[13])})
    return rows


def piecewise_reference(params, times):
    radius = params["cam_radius"]
    ecc = params["cam_eccentricity"]
    roller = params["roller_radius"]
    omega = params["motor_speed"]
    phase = params["phase"]
    y0 = params["follower_init_pos"][1]
    reach = radius + roller
    out = []
    onset_time = None
    for t in times:
        theta = phase + omega * t
        cx = ecc * math.cos(theta)
        cy = ecc * math.sin(theta)
        root = math.sqrt(max(reach * reach - cx * cx, 0.0))
        env_pos = cy + root
        dpos_dtheta = ecc * math.cos(theta) + (ecc * ecc * math.cos(theta) * math.sin(theta)) / max(root, 1.0e-12)
        env_vel = dpos_dtheta * omega
        if env_pos >= y0:
            if onset_time is None:
                onset_time = t
            out.append({"t": t, "pos": env_pos, "vel": env_vel})
        else:
            out.append({"t": t, "pos": y0, "vel": 0.0})
    return out, onset_time


def compute_metrics(rows, ref_rows):
    n = min(len(rows), len(ref_rows))
    pos_err = [rows[i]["pos"] - ref_rows[i]["pos"] for i in range(n)]
    vel_err = [rows[i]["vel"] - ref_rows[i]["vel"] for i in range(n)]
    return {
        "pos_rmse": math.sqrt(sum(e * e for e in pos_err) / n),
        "vel_rmse": math.sqrt(sum(e * e for e in vel_err) / n),
        "vel_mae": sum(abs(e) for e in vel_err) / n,
        "vel_max_abs": max(abs(e) for e in vel_err),
    }


def detect_motion_onset(rows, y0):
    for row in rows:
        if abs(row["pos"] - y0) > 1.0e-6 or abs(row["vel"]) > 1.0e-4:
            return row["t"]
    return None


def run_case(algorithm: str, output_path: Path):
    root_wsl = to_wsl_path(ROOT)
    output_wsl = to_wsl_path(output_path)
    cmd = (
        f"cd {root_wsl} && "
        f"./_build/project/baseline_onset_stress_nsc "
        f"--contact-algorithm {algorithm} "
        f"--output {output_wsl}"
    )
    t0 = time.perf_counter()
    proc = subprocess.run(
        ["wsl", "bash", "-lc", cmd],
        cwd=str(ROOT),
        text=True,
        capture_output=True,
        check=True,
    )
    elapsed = time.perf_counter() - t0
    match = re.search(r"\[PERF\]\s+Total Run Time:\s+([0-9.]+)\s+s", proc.stdout)
    runtime = float(match.group(1)) if match else elapsed
    return runtime, proc.stdout


def make_plots(ref_rows, series_map, onset_time):
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    times = [r["t"] for r in ref_rows]
    ref_pos = [r["pos"] for r in ref_rows]
    ref_vel = [r["vel"] for r in ref_rows]

    plt.figure(figsize=(10, 5.0))
    plt.plot(times, ref_pos, "k--", linewidth=2.0, label="Piecewise analytic reference")
    for label, rows, color in series_map:
        plt.plot(times, [r["pos"] for r in rows], color=color, linewidth=1.5, label=label)
    if onset_time is not None:
        plt.axvline(onset_time, color="#666666", linestyle=":", linewidth=1.1)
    plt.xlabel("Time (s)")
    plt.ylabel("Follower Y (m)")
    plt.title("Onset Stress Case: Position")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    pos_png = FIG_DIR / "onset_stress_position.png"
    plt.savefig(pos_png, dpi=180)
    plt.close()

    plt.figure(figsize=(10, 5.0))
    plt.plot(times, ref_vel, "k--", linewidth=2.0, label="Piecewise analytic reference")
    for label, rows, color in series_map:
        plt.plot(times, [r["vel"] for r in rows], color=color, linewidth=1.5, label=label)
    if onset_time is not None:
        plt.axvline(onset_time, color="#666666", linestyle=":", linewidth=1.1)
    plt.xlim(0.0, 0.3)
    plt.xlabel("Time (s)")
    plt.ylabel("Follower Y velocity (m/s)")
    plt.title("Onset Stress Case: Velocity")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    vel_png = FIG_DIR / "onset_stress_velocity.png"
    plt.savefig(vel_png, dpi=180)
    plt.close()

    return pos_png, vel_png


def write_report(params, results, onset_time, pos_png: Path, vel_png: Path):
    rows_html = []
    order = ["mesh", "sdf_1st", "sdf_2nd"]
    labels = {
        "mesh": "Chrono NSC native mesh",
        "sdf_1st": "SDF 1st-order",
        "sdf_2nd": "SDF 2nd-order",
    }
    for key in order:
        item = results[key]
        rows_html.append(
            "<tr>"
            f"<td>{html.escape(labels[key])}</td>"
            f"<td>{item['runtime']:.3f}</td>"
            f"<td>{item['metrics']['pos_rmse']:.6f}</td>"
            f"<td>{item['metrics']['vel_rmse']:.6f}</td>"
            f"<td>{item['metrics']['vel_max_abs']:.6f}</td>"
            f"<td>{item['motion_onset'] if item['motion_onset'] is not None else 'n/a'}</td>"
            "</tr>"
        )

    report = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8" />
  <title>Onset Stress Benchmark</title>
  <style>
    body {{
      font-family: "Times New Roman", Times, serif;
      margin: 24px auto;
      max-width: 1120px;
      line-height: 1.45;
      color: #111;
      padding: 0 18px 40px;
      background: #faf8f3;
    }}
    h1, h2 {{ margin: 0 0 10px; }}
    p {{ margin: 8px 0 16px; }}
    .cards {{
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(280px, 1fr));
      gap: 14px;
      margin: 18px 0 22px;
    }}
    .card {{
      background: white;
      border: 1px solid #ddd6c8;
      border-radius: 10px;
      padding: 14px 16px;
    }}
    table {{
      width: 100%;
      border-collapse: collapse;
      background: white;
      margin: 14px 0 24px;
    }}
    th, td {{
      border: 1px solid #d9d1c3;
      padding: 8px 10px;
      text-align: right;
    }}
    th:first-child, td:first-child {{ text-align: left; }}
    .fig {{
      background: white;
      border: 1px solid #ddd6c8;
      border-radius: 10px;
      padding: 12px;
      margin: 0 0 18px;
    }}
    img {{
      width: 100%;
      height: auto;
      display: block;
    }}
    code {{
      background: #f2eee6;
      padding: 1px 4px;
      border-radius: 4px;
    }}
  </style>
</head>
<body>
  <h1>Onset Stress Case (Version A)</h1>
  <p>
    This version-A benchmark keeps a smooth eccentric disk / roller geometry but starts from a positive initial
    clearance. The follower remains at constant height before first contact, then enters unilateral contact and
    continues sliding while the cam envelope rises. It is intended to isolate manifold onset behavior before
    introducing a spring-preloaded version B.
  </p>
  <div class="cards">
    <div class="card">
      <h2>Geometry</h2>
      <p>Cam radius: {params['cam_radius']:.3f} m<br/>
      Eccentricity: {params['cam_eccentricity']:.3f} m<br/>
      Roller radius: {params['roller_radius']:.3f} m</p>
    </div>
    <div class="card">
      <h2>Setup</h2>
      <p>Motor speed: {params['motor_speed']:.3f} rad/s<br/>
      Step size: {params['step_size']:.4f} s<br/>
      Total time: {params['total_time']:.3f} s<br/>
      Gravity: {params['gravity_y']:.2f} m/s²</p>
    </div>
    <div class="card">
      <h2>Target Onset</h2>
      <p>Configured onset target: {params['target_onset_time']:.3f} s<br/>
      Piecewise reference onset: {onset_time if onset_time is not None else 'n/a'} s</p>
    </div>
  </div>
  <table>
    <thead>
      <tr>
        <th>Method</th>
        <th>Runtime (s)</th>
        <th>Pos RMSE (m)</th>
        <th>Vel RMSE (m/s)</th>
        <th>Vel Max Abs (m/s)</th>
        <th>Detected motion onset (s)</th>
      </tr>
    </thead>
    <tbody>
      {"".join(rows_html)}
    </tbody>
  </table>
  <div class="fig">
    <h2>Position</h2>
    <img src="{pos_png.relative_to(DOCS_DIR).as_posix()}" alt="Onset stress position" />
  </div>
  <div class="fig">
    <h2>Velocity</h2>
    <img src="{vel_png.relative_to(DOCS_DIR).as_posix()}" alt="Onset stress velocity" />
  </div>
</body>
</html>
"""
    REPORT_PATH.write_text(report, encoding="utf-8", newline="\n")


def main():
    ensure_assets()
    params = json.loads(ASSET_JSON.read_text(encoding="utf-8"))
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    outputs = {
        "mesh": OUTPUT_DIR / "onset_stress_mesh.csv",
        "sdf_1st": OUTPUT_DIR / "onset_stress_sdf1.csv",
        "sdf_2nd": OUTPUT_DIR / "onset_stress_sdf2.csv",
    }
    results = {}
    for algorithm, output_path in outputs.items():
        print(f"[RUN] {algorithm} -> {output_path}")
        runtime, stdout = run_case(algorithm, output_path)
        rows = load_csv(output_path)
        ref_rows, onset_time = piecewise_reference(params, [row["t"] for row in rows])
        results[algorithm] = {
            "runtime": runtime,
            "stdout": stdout,
            "rows": rows,
            "metrics": compute_metrics(rows, ref_rows),
            "motion_onset": detect_motion_onset(rows, params["follower_init_pos"][1]),
        }

    ref_rows, onset_time = piecewise_reference(params, [row["t"] for row in results["mesh"]["rows"]])
    pos_png, vel_png = make_plots(
        ref_rows,
        [
            ("Chrono mesh", results["mesh"]["rows"], "#1f77b4"),
            ("SDF 1st-order", results["sdf_1st"]["rows"], "#d62728"),
            ("SDF 2nd-order", results["sdf_2nd"]["rows"], "#2ca02c"),
        ],
        onset_time,
    )
    write_report(params, results, onset_time, pos_png, vel_png)

    for algorithm in ("mesh", "sdf_1st", "sdf_2nd"):
        metrics = results[algorithm]["metrics"]
        print(
            f"[SUMMARY] {algorithm}: runtime={results[algorithm]['runtime']:.3f}s "
            f"pos_rmse={metrics['pos_rmse']:.6f} vel_rmse={metrics['vel_rmse']:.6f} "
            f"vel_max={metrics['vel_max_abs']:.6f} onset={results[algorithm]['motion_onset']}"
        )
    print(f"[REPORT] {REPORT_PATH}")


if __name__ == "__main__":
    main()
