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
ASSET_JSON = ROOT / "assets" / "onset_stress_b_oblique" / "onset_stress_b_oblique_model.json"
OUTPUT_DIR = ROOT / "data" / "outputs"
DOCS_DIR = ROOT / "docs"
FIG_DIR = DOCS_DIR / "figures"
REPORT_PATH = DOCS_DIR / "onset_stress_b_oblique_ablation_report.html"

FIT_RE = re.compile(r"fit_attempted=(\d+)\s+fit_applied=(\d+)\s+fit_reject=(\d+)")
PERF_RE = re.compile(r"\[PERF\]\s+Total Run Time:\s+([0-9.]+)\s+s")


def to_wsl_path(path: Path) -> str:
    drive = path.drive[:-1].lower()
    tail = path.as_posix()[2:]
    return f"/mnt/{drive}{tail}"


def ensure_assets():
    subprocess.run([sys.executable, str(ROOT / "scripts" / "generate_onset_stress_b_oblique_assets.py")], check=True)


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
    x_f = params["follower_init_pos"][0]
    y0 = params["follower_init_pos"][1]
    reach = radius + roller
    out = []
    onset_time = None
    for t in times:
        theta = phase + omega * t
        cx = ecc * math.cos(theta)
        cy = ecc * math.sin(theta)
        dx = x_f - cx
        root = math.sqrt(max(reach * reach - dx * dx, 0.0))
        env_pos = cy + root
        droot_dtheta = -(dx * ecc * math.sin(theta)) / max(root, 1.0e-12)
        dpos_dtheta = ecc * math.cos(theta) + droot_dtheta
        env_vel = dpos_dtheta * omega
        if env_pos >= y0:
            if onset_time is None:
                onset_time = t
            out.append({"t": t, "pos": env_pos, "vel": env_vel})
        else:
            out.append({"t": t, "pos": y0, "vel": 0.0})
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


def build_common_args(params, output_path: Path):
    return [
        "./_build/project/baseline_onset_stress_b_nsc",
        "--contact-algorithm",
        "sdf_2nd",
        "--output",
        to_wsl_path(output_path),
        "--dt",
        str(params["step_size"]),
        "--T",
        str(params["total_time"]),
        "--speed",
        str(params["motor_speed"]),
        "--spring-k",
        str(params["preload_stiffness"]),
        "--spring-c",
        str(params["preload_damping"]),
        "--spring-rest-length",
        str(params["preload_rest_length"]),
        "--spring-anchor-y",
        str(params["preload_anchor_pos"][1]),
        "--spring-anchor-x",
        str(params["preload_anchor_pos"][0]),
        "--follower-x",
        str(params["follower_init_pos"][0]),
        "--follower-y",
        str(params["follower_init_pos"][1]),
    ]


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
        "num_stat_lines": len(entries),
        "num_nonzero_lines": len(nonzero),
        "max_attempted": max((entry["attempted"] for entry in entries), default=0),
        "max_applied": max((entry["applied"] for entry in entries), default=0),
        "max_reject": max((entry["reject"] for entry in entries), default=0),
        "first_nonzero_line": nonzero[0]["line"] if nonzero else "none",
    }


def run_case(name: str, params, output_path: Path, extra_env: dict[str, str]):
    root_wsl = to_wsl_path(ROOT)
    cmd = "cd {} && {} {}".format(
        root_wsl,
        " ".join([f"{k}={v}" for k, v in {"SPCC_DEBUG_CONTACT_STATS": "1", **extra_env}.items()]),
        " ".join(build_common_args(params, output_path)),
    )
    proc = subprocess.run(
        ["wsl", "bash", "-lc", cmd],
        cwd=str(ROOT),
        text=True,
        capture_output=True,
        check=True,
    )
    match = PERF_RE.search(proc.stdout)
    runtime = float(match.group(1)) if match else float("nan")
    rows = load_csv(output_path)
    ref_rows, onset_time = piecewise_reference(params, [row["t"] for row in rows])
    return {
        "name": name,
        "runtime": runtime,
        "rows": rows,
        "ref_rows": ref_rows,
        "onset_time": onset_time,
        "motion_onset": detect_motion_onset(rows, params["follower_init_pos"][1]),
        "metrics": compute_metrics(rows, ref_rows, onset_time),
        "debug": parse_debug_stats(proc.stdout),
    }


def make_plot(results, onset_time):
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    times = [row["t"] for row in results[0]["ref_rows"]]
    ref_vel = [row["vel"] for row in results[0]["ref_rows"]]
    plt.figure(figsize=(10, 5.0))
    plt.plot(times, ref_vel, "k--", linewidth=2.0, label="Piecewise analytic reference")
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
    plt.title("Oblique Onset Stress B: Mechanism Ablation")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    out = FIG_DIR / "onset_stress_b_oblique_ablation_velocity.png"
    plt.savefig(out, dpi=180)
    plt.close()
    return out


def write_report(results, plot_path: Path):
    rows_html = []
    labels = {
        "default": "SDF 2nd-order",
        "no_gate": "SDF 2nd-order without onset gate",
        "no_fit": "SDF 2nd-order without local fit",
    }
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

    notes = []
    for result in results:
        d = result["debug"]
        notes.append(
            f"<li><strong>{html.escape(labels[result['name']])}</strong>: "
            f"<code>{html.escape(d['first_nonzero_line'])}</code></li>"
        )

    report = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8" />
  <title>Oblique Onset Stress B Ablation</title>
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
  <h1>Oblique Onset Stress B: Mechanism Ablation</h1>
  <p>
    This report compares the default second-order SDF configuration against two mechanism ablations on the same
    oblique-onset geometry: removing the positive-gap onset gate and removing the local-fit path correction. The
    current operating point is intentionally kept aligned with the existing Version B solver defaults; therefore the
    benchmark is used to confirm that the oblique onset genuinely exercises local-fit logic, even if the net trajectory
    difference remains small under this smooth single-roller geometry.
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
    <img src="{plot_path.relative_to(DOCS_DIR).as_posix()}" alt="Oblique onset ablation velocity" />
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
        ("default", {}, OUTPUT_DIR / "onset_stress_b_oblique_ablation_default.csv"),
        (
            "no_gate",
            {"SPCC_ONSET_B_LOCAL_FIT_REJECT_POSITIVE_PHI": "-1"},
            OUTPUT_DIR / "onset_stress_b_oblique_ablation_no_gate.csv",
        ),
        (
            "no_fit",
            {
                "SPCC_ONSET_B_LOCAL_FIT_ONSET_STEPS": "0",
                "SPCC_ONSET_B_SINGLE_POINT_LOCAL_FIT_PATH_SAMPLES": "0",
            },
            OUTPUT_DIR / "onset_stress_b_oblique_ablation_no_fit.csv",
        ),
    ]
    results = [run_case(name, params, output_path, env) for name, env, output_path in cases]
    plot_path = make_plot(results, results[0]["onset_time"])
    write_report(results, plot_path)
    print(f"[DONE] Wrote {REPORT_PATH}")


if __name__ == "__main__":
    main()
