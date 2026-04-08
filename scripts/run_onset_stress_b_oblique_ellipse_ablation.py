import csv
import html
import json
import math
import re
import subprocess
import sys
from pathlib import Path

import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parents[1]
ASSET_JSON = ROOT / "assets" / "onset_stress_b_oblique_ellipse" / "onset_stress_b_oblique_ellipse_model.json"
OUTPUT_DIR = ROOT / "data" / "outputs"
DOCS_DIR = ROOT / "docs"
FIG_DIR = DOCS_DIR / "figures"
REPORT_PATH = DOCS_DIR / "onset_stress_b_oblique_ellipse_ablation_report.html"

FIT_RE = re.compile(r"fit_attempted=(\d+)\s+fit_applied=(\d+)\s+fit_reject=(\d+)")
PERF_RE = re.compile(r"\[PERF\]\s+Total Run Time:\s+([0-9.]+)\s+s")


def to_wsl_path(path: Path) -> str:
    drive = path.drive[:-1].lower()
    tail = path.as_posix()[2:]
    return f"/mnt/{drive}{tail}"


def ensure_assets():
    subprocess.run([sys.executable, str(ROOT / "scripts" / "generate_onset_stress_b_oblique_ellipse_assets.py")], check=True)


def load_csv(path: Path):
    rows = []
    with path.open("r", encoding="utf-8") as f:
        reader = csv.reader(f)
        next(reader)
        for row in reader:
            if row:
                rows.append({"t": float(row[0]), "pos": float(row[10]), "vel": float(row[13])})
    return rows


def ellipse_reference(params, times):
    cam_radius = params["cam_radius"]
    ecc = params["cam_eccentricity"]
    phase = params["phase"]
    omega = params["motor_speed"]
    follower_x = params["follower_init_pos"][0]
    follower_y0 = params["follower_init_pos"][1]
    ellipse_ax = params["follower_ellipse_ax"]
    ellipse_ay = params["follower_ellipse_ay"]

    out = []
    onset_time = None
    for t in times:
        theta = phase + omega * t
        cx = ecc * math.cos(theta)
        cy = ecc * math.sin(theta)
        dx = follower_x - cx
        if abs(dx) <= ellipse_ax:
            s = max(1.0 - (dx * dx) / (ellipse_ax * ellipse_ax), 0.0)
            env_pos = cy + cam_radius + ellipse_ay * math.sqrt(s)
            droot_dtheta = 0.0
            if s > 1.0e-12:
                droot_dtheta = (
                    ellipse_ay * (-dx) * ecc * math.sin(theta) /
                    (ellipse_ax * ellipse_ax * math.sqrt(s))
                )
            env_vel = (ecc * math.cos(theta) + droot_dtheta) * omega
        else:
            env_pos = -1.0e9
            env_vel = 0.0

        if env_pos >= follower_y0:
            if onset_time is None:
                onset_time = t
            out.append({"t": t, "pos": env_pos, "vel": env_vel})
        else:
            out.append({"t": t, "pos": follower_y0, "vel": 0.0})
    return out, onset_time


def compute_metrics(rows, ref_rows, onset_time):
    n = min(len(rows), len(ref_rows))
    pos_err = [rows[i]["pos"] - ref_rows[i]["pos"] for i in range(n)]
    vel_err = [rows[i]["vel"] - ref_rows[i]["vel"] for i in range(n)]
    onset_lo = (onset_time if onset_time is not None else 0.15) - 0.03
    onset_hi = (onset_time if onset_time is not None else 0.15) + 0.06
    onset_idx = [i for i in range(n) if onset_lo <= rows[i]["t"] <= onset_hi]
    return {
        "pos_rmse": math.sqrt(sum(e * e for e in pos_err) / n),
        "vel_rmse": math.sqrt(sum(e * e for e in vel_err) / n),
        "vel_max_abs": max(abs(e) for e in vel_err),
        "onset_vel_rmse": math.sqrt(sum(vel_err[i] * vel_err[i] for i in onset_idx) / max(len(onset_idx), 1)),
        "onset_vel_max_abs": max((abs(vel_err[i]) for i in onset_idx), default=0.0),
    }


def detect_motion_onset(rows, y0):
    for row in rows:
        if abs(row["pos"] - y0) > 1.0e-6 or abs(row["vel"]) > 1.0e-4:
            return row["t"]
    return None


def parse_debug_stats(stdout: str):
    entries = []
    for line in stdout.splitlines():
        match = FIT_RE.search(line)
        if match:
            entries.append(
                {
                    "attempted": int(match.group(1)),
                    "applied": int(match.group(2)),
                    "reject": int(match.group(3)),
                    "line": line.strip(),
                }
            )
    nonzero = [entry for entry in entries if entry["attempted"] > 0 or entry["applied"] > 0 or entry["reject"] > 0]
    return {
        "max_attempted": max((entry["attempted"] for entry in entries), default=0),
        "max_applied": max((entry["applied"] for entry in entries), default=0),
        "max_reject": max((entry["reject"] for entry in entries), default=0),
        "first_nonzero_line": nonzero[0]["line"] if nonzero else "none",
    }


def run_case(name: str, params, extra_env: dict[str, str], output_name: str):
    output_path = OUTPUT_DIR / output_name
    env = {
        "SPCC_DEBUG_CONTACT_STATS": "1",
        "SPCC_ONSET_B_DYNAMICS_SUBSTEPS": "8",
        "SPCC_ONSET_B_SURFACE_RES": "0.0005",
        "SPCC_ONSET_B_VOXEL_SIZE": "0.0002",
        "SPCC_ONSET_B_AVG_POINT": "1",
        "SPCC_ONSET_B_CLUSTER_RADIUS": "0.005",
        "SPCC_ONSET_B_LOCAL_SCAN_RADIUS": "0.01",
        "SPCC_ONSET_B_LOCAL_FIT_ONSET_STEPS": "1",
        "SPCC_ONSET_B_LOCAL_FIT_MIN_CLUSTER_SIZE": "2",
        "SPCC_ONSET_B_LOCAL_FIT_PATH_SAMPLES": "5",
        "SPCC_ONSET_B_LOCAL_FIT_MAX_SHIFT_RATIO": "1.0",
        "SPCC_ONSET_B_LOCAL_FIT_BLEND": "1.0",
        "SPCC_ONSET_B_LOCAL_FIT_REJECT_POSITIVE_PHI": "0.0",
        "SPCC_ONSET_B_ONSET_GATE_CURRENT_PHI_MAX": "0.0",
        "SPCC_ONSET_B_SINGLE_POINT_LOCAL_FIT_PATH_SAMPLES": "0",
    }
    env.update(extra_env)
    cmd = "cd {} && {} ./_build/project/baseline_onset_stress_b_nsc --cam-mesh {} --follower-mesh {} --contact-algorithm sdf_2nd --output {} --dt {} --T {} --speed {} --spring-k {} --spring-c {} --spring-rest-length {} --spring-anchor-y {} --spring-anchor-x {} --follower-x {} --follower-y {}".format(
        to_wsl_path(ROOT),
        " ".join([f"{k}={v}" for k, v in env.items()]),
        to_wsl_path(ROOT / "assets" / "onset_stress_b_oblique_ellipse" / "models" / "onset_cam.obj"),
        to_wsl_path(ROOT / "assets" / "onset_stress_b_oblique_ellipse" / "models" / "ellipse_follower.obj"),
        to_wsl_path(output_path),
        params["step_size"],
        params["total_time"],
        params["motor_speed"],
        params["preload_stiffness"],
        params["preload_damping"],
        params["preload_rest_length"],
        params["preload_anchor_pos"][1],
        params["preload_anchor_pos"][0],
        params["follower_init_pos"][0],
        params["follower_init_pos"][1],
    )
    proc = subprocess.run(["wsl", "bash", "-lc", cmd], cwd=str(ROOT), text=True, capture_output=True, check=True)
    runtime_match = PERF_RE.search(proc.stdout)
    runtime = float(runtime_match.group(1)) if runtime_match else float("nan")
    rows = load_csv(output_path)
    ref_rows, onset_time = ellipse_reference(params, [row["t"] for row in rows])
    return {
        "name": name,
        "runtime": runtime,
        "rows": rows,
        "ref_rows": ref_rows,
        "onset_time": onset_time,
        "motion_onset": detect_motion_onset(rows, params["follower_init_pos"][1]),
        "metrics": compute_metrics(rows, ref_rows, onset_time),
        "debug": parse_debug_stats(proc.stdout),
        "output_path": output_path,
    }


def make_plot(results, onset_time):
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    times = [row["t"] for row in results[0]["ref_rows"]]
    ref_vel = [row["vel"] for row in results[0]["ref_rows"]]
    plt.figure(figsize=(10, 5.1))
    plt.plot(times, ref_vel, "k--", linewidth=2.0, label="Numerical smooth-reference envelope")
    styles = [
        ("SDF 2nd-order", "#2ca02c", "-"),
        ("SDF 2nd-order without onset gate", "#ff7f0e", ":"),
        ("SDF 2nd-order without local fit", "#d62728", "--"),
    ]
    for result, (label, color, style) in zip(results, styles):
        plt.plot(times, [row["vel"] for row in result["rows"]], color=color, linestyle=style, linewidth=1.6, label=label)
    if onset_time is not None:
        plt.axvspan(onset_time - 0.03, onset_time + 0.06, color="#f1e4c6", alpha=0.35, zorder=0)
        plt.axvline(onset_time, color="#666666", linestyle=":", linewidth=1.1)
    plt.xlabel("Time (s)")
    plt.ylabel("Follower Y velocity (m/s)")
    plt.title("Oblique Ellipse Follower: Mechanism Ablation")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    out = FIG_DIR / "onset_stress_b_oblique_ellipse_ablation_velocity.png"
    plt.savefig(out, dpi=180)
    plt.close()
    return out


def write_report(results, plot_path: Path):
    labels = {
        "default": "SDF 2nd-order",
        "no_gate": "SDF 2nd-order without onset gate",
        "no_fit": "SDF 2nd-order without local fit",
    }
    rows_html = []
    notes = []
    for result in results:
        m = result["metrics"]
        d = result["debug"]
        rows_html.append(
            "<tr>"
            f"<td>{html.escape(labels[result['name']])}</td>"
            f"<td>{result['runtime']:.3f}</td>"
            f"<td>{m['vel_rmse']:.6f}</td>"
            f"<td>{m['vel_max_abs']:.6f}</td>"
            f"<td>{m['onset_vel_rmse']:.6f}</td>"
            f"<td>{m['onset_vel_max_abs']:.6f}</td>"
            f"<td>{result['motion_onset'] if result['motion_onset'] is not None else 'n/a'}</td>"
            f"<td>{d['max_attempted']}</td>"
            f"<td>{d['max_applied']}</td>"
            f"<td>{d['max_reject']}</td>"
            "</tr>"
        )
        notes.append(
            f"<li><strong>{html.escape(labels[result['name']])}</strong>: "
            f"<code>{html.escape(d['first_nonzero_line'])}</code></li>"
        )

    report = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8" />
  <title>Oblique Ellipse Follower Ablation</title>
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
      font-family: Consolas, "Courier New", monospace;
      font-size: 0.92em;
    }}
  </style>
</head>
<body>
  <h1>Oblique Ellipse Follower: Mechanism Ablation</h1>
  <p>
    This stronger oblique benchmark replaces the circular roller with a horizontally elongated smooth ellipse.
    The wider contact footprint makes patch-level local fit materially affect the onset trajectory. Under this geometry,
    removing local fit measurably increases the onset-window peak error, while the current onset gate still mainly changes
    internal candidate rejection statistics rather than the final trajectory.
  </p>
  <table>
    <thead>
      <tr>
        <th>Method</th>
        <th>Runtime (s)</th>
        <th>Vel RMSE (m/s)</th>
        <th>Vel Max Abs (m/s)</th>
        <th>Onset-window Vel RMSE (m/s)</th>
        <th>Onset-window Vel Max Abs (m/s)</th>
        <th>Detected motion onset (s)</th>
        <th>Max fit attempted</th>
        <th>Max fit applied</th>
        <th>Max positive-gap reject</th>
      </tr>
    </thead>
    <tbody>
      {"".join(rows_html)}
    </tbody>
  </table>
  <div class="fig">
    <h2>Velocity Comparison</h2>
    <img src="{plot_path.relative_to(DOCS_DIR).as_posix()}" alt="Oblique ellipse follower ablation" />
  </div>
  <h2>First Nonzero Debug Lines</h2>
  <ul>
    {"".join(notes)}
  </ul>
</body>
</html>
"""
    REPORT_PATH.write_text(report, encoding="utf-8", newline="\n")


def main():
    ensure_assets()
    params = json.loads(ASSET_JSON.read_text(encoding="utf-8"))
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    cases = [
        ("default", {}, "onset_stress_b_oblique_ellipse_default.csv"),
        (
            "no_gate",
            {
                "SPCC_ONSET_B_LOCAL_FIT_REJECT_POSITIVE_PHI": "-1",
                "SPCC_ONSET_B_ONSET_GATE_CURRENT_PHI_MAX": "-1",
            },
            "onset_stress_b_oblique_ellipse_no_gate.csv",
        ),
        (
            "no_fit",
            {
                "SPCC_ONSET_B_LOCAL_FIT_ONSET_STEPS": "0",
                "SPCC_ONSET_B_ONSET_GATE_CURRENT_PHI_MAX": "-1",
            },
            "onset_stress_b_oblique_ellipse_no_fit.csv",
        ),
    ]
    results = [run_case(name, params, env, output_name) for name, env, output_name in cases]
    plot_path = make_plot(results, results[0]["onset_time"])
    write_report(results, plot_path)
    print(f"[DONE] Wrote {REPORT_PATH}")


if __name__ == "__main__":
    main()
