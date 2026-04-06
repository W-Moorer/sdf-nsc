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
ASSET_JSON = ROOT / "assets" / "eccentric_roller" / "eccentric_roller_model.json"
OUTPUT_DIR = ROOT / "data" / "outputs"
DOCS_DIR = ROOT / "docs"
FIG_DIR = DOCS_DIR / "figures"
REPORT_PATH = DOCS_DIR / "eccentric_roller_report.html"


def to_wsl_path(path: Path) -> str:
    drive = path.drive[:-1].lower()
    tail = path.as_posix()[2:]
    return f"/mnt/{drive}{tail}"


def ensure_assets():
    subprocess.run([sys.executable, str(ROOT / "scripts" / "generate_eccentric_roller_assets.py")], check=True)


def load_csv(path: Path):
    rows = []
    with path.open("r", encoding="utf-8") as f:
        reader = csv.reader(f)
        next(reader)
        for row in reader:
            if not row:
                continue
            rows.append(
                {
                    "t": float(row[0]),
                    "pos": float(row[10]),
                    "vel": float(row[13]),
                }
            )
    return rows


def analytic_reference(params, times):
    radius = params["cam_radius"]
    ecc = params["cam_eccentricity"]
    roller = params["roller_radius"]
    omega = params["motor_speed"]
    reach = radius + roller
    out = []
    for t in times:
        theta = omega * t
        cos_t = math.cos(theta)
        sin_t = math.sin(theta)
        cx = ecc * cos_t
        cy = ecc * sin_t
        root = math.sqrt(max(reach * reach - cx * cx, 0.0))
        pos = cy + root
        dpos_dtheta = ecc * cos_t + (ecc * ecc * cos_t * sin_t) / max(root, 1.0e-12)
        vel = dpos_dtheta * omega
        out.append({"t": t, "pos": pos, "vel": vel})
    return out


def compute_metrics(rows, ref_rows):
    n = min(len(rows), len(ref_rows))
    pos_err = [rows[i]["pos"] - ref_rows[i]["pos"] for i in range(n)]
    vel_err = [rows[i]["vel"] - ref_rows[i]["vel"] for i in range(n)]
    return {
        "pos_rmse": math.sqrt(sum(e * e for e in pos_err) / n),
        "pos_mae": sum(abs(e) for e in pos_err) / n,
        "vel_rmse": math.sqrt(sum(e * e for e in vel_err) / n),
        "vel_mae": sum(abs(e) for e in vel_err) / n,
        "vel_max_abs": max(abs(e) for e in vel_err),
    }


def run_case(algorithm: str, output_path: Path):
    root_wsl = to_wsl_path(ROOT)
    output_wsl = to_wsl_path(output_path)
    cmd = (
        f"cd {root_wsl} && "
        f"./_build/project/baseline_eccentric_roller_nsc "
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


def make_plots(ref_rows, series_map):
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    times = [r["t"] for r in ref_rows]
    ref_pos = [r["pos"] for r in ref_rows]
    ref_vel = [r["vel"] for r in ref_rows]

    plt.figure(figsize=(10, 5.2))
    plt.plot(times, ref_pos, "k--", linewidth=2.0, label="Analytic reference")
    for label, rows, color in series_map:
        plt.plot(times, [r["pos"] for r in rows], color=color, linewidth=1.6, label=label)
    plt.xlabel("Time (s)")
    plt.ylabel("Follower Y (m)")
    plt.title("Eccentric Disk / Roller Follower: Position")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    pos_png = FIG_DIR / "eccentric_roller_position.png"
    plt.savefig(pos_png, dpi=180)
    plt.close()

    plt.figure(figsize=(10, 5.2))
    plt.plot(times, ref_vel, "k--", linewidth=2.0, label="Analytic reference")
    for label, rows, color in series_map:
        plt.plot(times, [r["vel"] for r in rows], color=color, linewidth=1.6, label=label)
    plt.xlabel("Time (s)")
    plt.ylabel("Follower Y velocity (m/s)")
    plt.title("Eccentric Disk / Roller Follower: Velocity")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    vel_png = FIG_DIR / "eccentric_roller_velocity.png"
    plt.savefig(vel_png, dpi=180)
    plt.close()

    plt.figure(figsize=(10, 4.8))
    for label, rows, color in series_map:
        errs = [rows[i]["vel"] - ref_rows[i]["vel"] for i in range(min(len(rows), len(ref_rows)))]
        plt.plot(times[: len(errs)], errs, color=color, linewidth=1.4, label=f"{label} - analytic")
    plt.axhline(0.0, color="k", linewidth=1.0, alpha=0.4)
    plt.xlabel("Time (s)")
    plt.ylabel("Velocity error (m/s)")
    plt.title("Velocity Error Against Analytic Reference")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    err_png = FIG_DIR / "eccentric_roller_velocity_error.png"
    plt.savefig(err_png, dpi=180)
    plt.close()

    return pos_png, vel_png, err_png


def write_report(params, results, pos_png: Path, vel_png: Path, err_png: Path):
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
            "</tr>"
        )

    report = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8" />
  <title>Eccentric Roller Benchmark</title>
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
  <h1>Eccentric Disk / Roller Follower Benchmark</h1>
  <p>
    This benchmark compares Chrono native mesh contact against the current SDF 1st-order and 2nd-order
    contact pipelines on a procedurally generated eccentric disk cam and vertically constrained roller follower.
    The analytic reference is the ideal no-separation geometry
    <code>y(t)=e sin(&omega;t)+sqrt((R+r)^2-(e cos(&omega;t))^2)</code>.
  </p>

  <div class="cards">
    <div class="card">
      <h2>Geometry</h2>
      <p>Cam radius: {params['cam_radius']:.3f} m<br/>
      Eccentricity: {params['cam_eccentricity']:.3f} m<br/>
      Roller radius: {params['roller_radius']:.3f} m</p>
    </div>
    <div class="card">
      <h2>Dynamics</h2>
      <p>Motor speed: {params['motor_speed']:.3f} rad/s<br/>
      Step size: {params['step_size']:.4f} s<br/>
      Total time: {params['total_time']:.6f} s</p>
    </div>
    <div class="card">
      <h2>Interpretation</h2>
      <p>
        Native mesh is the Chrono baseline. The two SDF modes reuse the current sliding-patch contact path and
        can be read directly against both the baseline and the analytic envelope.
      </p>
    </div>
    <div class="card">
      <h2>Dedicated SDF 2nd Defaults</h2>
      <p>
        This benchmark uses a dedicated <code>SPCC_ECC</code> tuning profile with
        <code>max_active_keep=1</code>, <code>cluster_radius=1.5e-3</code>,
        <code>local_scan_radius=6e-3</code>, <code>avg_point=false</code>,
        and <code>dynamics_substeps=2</code>.
      </p>
    </div>
  </div>

  <h2>Metrics</h2>
  <table>
    <thead>
      <tr>
        <th>Method</th>
        <th>Runtime (s)</th>
        <th>Pos RMSE (m)</th>
        <th>Vel RMSE (m/s)</th>
        <th>Vel Max Abs (m/s)</th>
      </tr>
    </thead>
    <tbody>
      {"".join(rows_html)}
    </tbody>
  </table>

  <div class="fig">
    <h2>Position</h2>
    <img src="{pos_png.relative_to(DOCS_DIR).as_posix()}" alt="Position comparison" />
  </div>
  <div class="fig">
    <h2>Velocity</h2>
    <img src="{vel_png.relative_to(DOCS_DIR).as_posix()}" alt="Velocity comparison" />
  </div>
  <div class="fig">
    <h2>Velocity Error</h2>
    <img src="{err_png.relative_to(DOCS_DIR).as_posix()}" alt="Velocity error comparison" />
  </div>
</body>
</html>
"""
    REPORT_PATH.write_text(report, encoding="utf-8", newline="\n")


def main():
    ensure_assets()
    params = json.loads(ASSET_JSON.read_text(encoding="utf-8"))

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    results = {}
    outputs = {
        "mesh": OUTPUT_DIR / "eccentric_roller_mesh.csv",
        "sdf_1st": OUTPUT_DIR / "eccentric_roller_sdf1.csv",
        "sdf_2nd": OUTPUT_DIR / "eccentric_roller_sdf2.csv",
    }

    for algorithm, output_path in outputs.items():
        print(f"[RUN] {algorithm} -> {output_path}")
        runtime, stdout = run_case(algorithm, output_path)
        results[algorithm] = {"runtime": runtime, "stdout": stdout, "rows": load_csv(output_path)}

    ref_rows = analytic_reference(params, [row["t"] for row in results["mesh"]["rows"]])
    for algorithm in outputs:
        results[algorithm]["metrics"] = compute_metrics(results[algorithm]["rows"], ref_rows)

    pos_png, vel_png, err_png = make_plots(
        ref_rows,
        [
            ("Chrono mesh", results["mesh"]["rows"], "#1f77b4"),
            ("SDF 1st-order", results["sdf_1st"]["rows"], "#d62728"),
            ("SDF 2nd-order", results["sdf_2nd"]["rows"], "#2ca02c"),
        ],
    )
    write_report(params, results, pos_png, vel_png, err_png)

    for algorithm in ("mesh", "sdf_1st", "sdf_2nd"):
        metrics = results[algorithm]["metrics"]
        print(
            f"[SUMMARY] {algorithm}: runtime={results[algorithm]['runtime']:.3f}s "
            f"pos_rmse={metrics['pos_rmse']:.6f} vel_rmse={metrics['vel_rmse']:.6f} "
            f"vel_max={metrics['vel_max_abs']:.6f}"
        )
    print(f"[REPORT] {REPORT_PATH}")


if __name__ == "__main__":
    main()
