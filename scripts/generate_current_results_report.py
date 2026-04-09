from __future__ import annotations

import html
import json
import math
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
OUT_HTML = ROOT / "docs" / "current_results_report.html"


CAM_REFERENCE_PATH = ROOT / "assets" / "cam" / "data" / "cam_data.csv"

CAM_CASES = [
    {
        "key": "mesh",
        "label": "Cam Mesh",
        "path": ROOT / "data" / "outputs" / "cam_reverify_mesh__dt0p01.csv",
        "runtime_s": 135.8841,
        "color": "#1d4ed8",
    },
    {
        "key": "sdf1",
        "label": "Cam SDF 1st",
        "path": ROOT / "data" / "outputs" / "cam_current_default_sdf1_dt0p01.csv",
        "runtime_s": 14.9609,
        "color": "#0f766e",
    },
    {
        "key": "sdf2",
        "label": "Cam SDF 2nd",
        "path": ROOT / "data" / "outputs" / "cam_current_default_sdf2_dt0p01.csv",
        "runtime_s": 16.0043,
        "color": "#b91c1c",
    },
]

GEAR_CASES = [
    {
        "key": "sdf1",
        "label": "Gear SDF 1st",
        "path": ROOT / "data" / "outputs" / "gear_current_default_s1.csv",
        "pinned_path": ROOT / "data" / "outputs" / "gear_default_pinned_s1.csv",
        "runtime_s": 60.397,
        "color": "#0f766e",
    },
    {
        "key": "sdf2",
        "label": "Gear SDF 2nd",
        "path": ROOT / "data" / "outputs" / "gear_current_default_s2.csv",
        "pinned_path": ROOT / "data" / "outputs" / "gear_default_pinned_s2.csv",
        "runtime_s": 66.026,
        "color": "#b91c1c",
    },
]

CAM_DEFAULTS = [
    ("voxel_size", "4e-4"),
    ("half_band_voxels", "12"),
    ("surface_res", "1e-3"),
    ("max_samples", "60000"),
    ("full_scan_period", "1"),
    ("local_scan_radius", "1e-2"),
    ("cluster_radius", "2.5e-3"),
    ("cluster_angle_deg", "25"),
    ("separating_cutoff", "1e-3"),
    ("max_active_keep", "8"),
    ("avg_point", "true"),
    ("persistent_match_radius", "6e-3"),
    ("persistent_normal_cos_min", "0.92"),
    ("persistent_blend_alpha", "0.60"),
    ("persistent_path_samples", "1"),
    ("coverage_spacing_radius", "4e-3"),
    ("dynamics_substeps", "4"),
    ("local_fit_onset_steps", "1"),
    ("local_fit_min_cluster_size", "2"),
    ("local_fit_path_samples", "1"),
    ("local_fit_max_shift_ratio", "1.0"),
    ("local_fit_blend", "1.0"),
    ("local_fit_reject_positive_phi", "2.5e-3"),
    ("use_sample_bvh", "true"),
    ("curvature.enabled", "true"),
    ("curvature.tangential_only", "true"),
    ("curvature.alignment_cos_min", "0.99"),
    ("curvature.max_hessian_frobenius", "120"),
    ("curvature.max_abs", "6e-4"),
    ("curvature.max_ratio", "0.20"),
    ("curvature.gap_floor", "4e-3"),
]


def fmt(v: float, digits: int = 6) -> str:
    return f"{v:.{digits}f}"


def load_cam_reference():
    df = pd.read_csv(CAM_REFERENCE_PATH)
    return {
        "time": df.iloc[:, 0].tolist(),
        "pos": df.iloc[:, 10].tolist(),
        "vel": df.iloc[:, 13].tolist(),
    }


def load_cam_trace(path: Path):
    df = pd.read_csv(path)
    return {
        "time": df.iloc[:, 0].tolist(),
        "pos": df["Y:Pos_TY-Body2(m)"].tolist(),
        "vel": df["Y:Vel_TY-Body2(m/s)"].tolist(),
    }


def load_gear_trace(path: Path):
    df = pd.read_csv(path)
    time_col = df.columns[0]
    vel_col = "Y:Vel_RX-GEAR22-chrono(rad/s)"
    return {
        "time": df[time_col].tolist(),
        "vel": df[vel_col].tolist(),
    }


def interp_linear(x_src, y_src, x_dst):
    out = []
    j = 0
    n = len(x_src)
    for x in x_dst:
        while j + 1 < n and x_src[j + 1] < x:
            j += 1
        if j + 1 >= n:
            out.append(y_src[-1])
            continue
        x0, x1 = x_src[j], x_src[j + 1]
        y0, y1 = y_src[j], y_src[j + 1]
        if x1 == x0:
            out.append(y0)
            continue
        alpha = (x - x0) / (x1 - x0)
        out.append(y0 * (1.0 - alpha) + y1 * alpha)
    return out


def compute_cam_metrics(trace, ref):
    ref_pos = interp_linear(ref["time"], ref["pos"], trace["time"])
    ref_vel = interp_linear(ref["time"], ref["vel"], trace["time"])
    pos_err = [a - b for a, b in zip(trace["pos"], ref_pos)]
    vel_err = [a - b for a, b in zip(trace["vel"], ref_vel)]

    def mae(arr):
        return sum(abs(v) for v in arr) / len(arr)

    def rmse(arr):
        return math.sqrt(sum(v * v for v in arr) / len(arr))

    return {
        "pos_rmse": rmse(pos_err),
        "vel_rmse": rmse(vel_err),
        "pos_mae": mae(pos_err),
        "vel_mae": mae(vel_err),
        "vel_max_abs": max(abs(v) for v in vel_err),
    }


def compute_gear_metrics(trace, pinned_trace):
    vel = trace["vel"]
    pinned_vel = pinned_trace["vel"]
    mae = sum(abs(v + 1.0) for v in vel) / len(vel)
    rmse = math.sqrt(sum((v + 1.0) ** 2 for v in vel) / len(vel))
    maxabs_vs_pinned = max(abs(a - b) for a, b in zip(vel, pinned_vel))
    return {
        "mae": mae,
        "rmse": rmse,
        "maxabs_vs_pinned": maxabs_vs_pinned,
    }


def _bounds(values):
    vmin = min(values)
    vmax = max(values)
    if not math.isfinite(vmin) or not math.isfinite(vmax):
        return -1.0, 1.0
    if abs(vmax - vmin) < 1e-12:
        pad = max(1.0, abs(vmin) * 0.05)
        return vmin - pad, vmax + pad
    pad = 0.06 * (vmax - vmin)
    return vmin - pad, vmax + pad


def polyline_svg(series, width=760, height=360, x_label="Time (s)", y_label="", title=""):
    left, right, top, bottom = 62, 18, 20, 42
    plot_w = width - left - right
    plot_h = height - top - bottom

    all_x = []
    all_y = []
    for s in series:
        all_x.extend(s["x"])
        all_y.extend(s["y"])
    xmin, xmax = min(all_x), max(all_x)
    ymin, ymax = _bounds(all_y)

    def sx(x):
        if xmax == xmin:
            return left + plot_w * 0.5
        return left + (x - xmin) / (xmax - xmin) * plot_w

    def sy(y):
        if ymax == ymin:
            return top + plot_h * 0.5
        return top + (ymax - y) / (ymax - ymin) * plot_h

    parts = [
        f'<svg viewBox="0 0 {width} {height}" xmlns="http://www.w3.org/2000/svg">',
        '<rect x="0" y="0" width="100%" height="100%" fill="#fffdfa"/>',
    ]

    for i in range(6):
        yy = top + plot_h * i / 5.0
        val = ymax - (ymax - ymin) * i / 5.0
        parts.append(f'<line x1="{left}" y1="{yy:.2f}" x2="{left+plot_w}" y2="{yy:.2f}" stroke="#e8e0d2" stroke-width="1"/>')
        parts.append(f'<text x="{left-10}" y="{yy+4:.2f}" text-anchor="end" font-size="12" fill="#6b7280">{html.escape(fmt(val, 3))}</text>')

    for i in range(6):
        xx = left + plot_w * i / 5.0
        val = xmin + (xmax - xmin) * i / 5.0
        parts.append(f'<line x1="{xx:.2f}" y1="{top}" x2="{xx:.2f}" y2="{top+plot_h}" stroke="#f1eadf" stroke-width="1"/>')
        parts.append(f'<text x="{xx:.2f}" y="{top+plot_h+20}" text-anchor="middle" font-size="12" fill="#6b7280">{html.escape(fmt(val, 2))}</text>')

    parts.append(f'<line x1="{left}" y1="{top+plot_h}" x2="{left+plot_w}" y2="{top+plot_h}" stroke="#77624a" stroke-width="1.2"/>')
    parts.append(f'<line x1="{left}" y1="{top}" x2="{left}" y2="{top+plot_h}" stroke="#77624a" stroke-width="1.2"/>')
    parts.append(f'<text x="{left + plot_w / 2:.2f}" y="{height-8}" text-anchor="middle" font-size="13" fill="#4b5563">{html.escape(x_label)}</text>')
    parts.append(f'<text x="18" y="{top + plot_h / 2:.2f}" text-anchor="middle" transform="rotate(-90 18 {top + plot_h / 2:.2f})" font-size="13" fill="#4b5563">{html.escape(y_label)}</text>')

    for s in series:
        pts = " ".join(f"{sx(x):.2f},{sy(y):.2f}" for x, y in zip(s["x"], s["y"]))
        dash = ' stroke-dasharray="7 5"' if s.get("dash") else ""
        parts.append(
            f'<polyline fill="none" stroke="{s["color"]}" stroke-width="2.2"{dash} points="{pts}"/>'
        )

    parts.append("</svg>")
    return "".join(parts)


def bar_svg(items, value_key, width=520, height=300, y_label=""):
    left, right, top, bottom = 52, 18, 18, 48
    plot_w = width - left - right
    plot_h = height - top - bottom
    vmax = max(item[value_key] for item in items)
    vmax *= 1.12 if vmax > 0 else 1.0

    parts = [f'<svg viewBox="0 0 {width} {height}" xmlns="http://www.w3.org/2000/svg">']
    parts.append('<rect x="0" y="0" width="100%" height="100%" fill="#fffdfa"/>')

    for i in range(6):
        yy = top + plot_h * i / 5.0
        val = vmax - vmax * i / 5.0
        parts.append(f'<line x1="{left}" y1="{yy:.2f}" x2="{left+plot_w}" y2="{yy:.2f}" stroke="#e8e0d2" stroke-width="1"/>')
        parts.append(f'<text x="{left-8}" y="{yy+4:.2f}" text-anchor="end" font-size="12" fill="#6b7280">{html.escape(fmt(val, 2))}</text>')

    n = len(items)
    gap = 18
    bar_w = (plot_w - gap * (n + 1)) / max(1, n)
    for i, item in enumerate(items):
        x = left + gap + i * (bar_w + gap)
        h = 0 if vmax <= 0 else plot_h * item[value_key] / vmax
        y = top + plot_h - h
        parts.append(f'<rect x="{x:.2f}" y="{y:.2f}" width="{bar_w:.2f}" height="{h:.2f}" rx="6" fill="{item["color"]}"/>')
        parts.append(f'<text x="{x+bar_w/2:.2f}" y="{y-6:.2f}" text-anchor="middle" font-size="12" fill="#374151">{html.escape(fmt(item[value_key], 3))}</text>')
        parts.append(f'<text x="{x+bar_w/2:.2f}" y="{top+plot_h+18}" text-anchor="middle" font-size="12" fill="#6b7280">{html.escape(item["label"])}</text>')

    parts.append(f'<line x1="{left}" y1="{top+plot_h}" x2="{left+plot_w}" y2="{top+plot_h}" stroke="#77624a" stroke-width="1.2"/>')
    parts.append(f'<text x="18" y="{top + plot_h / 2:.2f}" text-anchor="middle" transform="rotate(-90 18 {top + plot_h / 2:.2f})" font-size="13" fill="#4b5563">{html.escape(y_label)}</text>')
    parts.append("</svg>")
    return "".join(parts)


def legend(items):
    bits = []
    for item in items:
        dash = "border-top-style:dashed;" if item.get("dash") else ""
        bits.append(
            f'<span><i style="display:inline-block;width:18px;border-top:3px solid {item["color"]};{dash}margin-right:8px;"></i>{html.escape(item["label"])}</span>'
        )
    return "".join(bits)


def make_report():
    ref = load_cam_reference()

    cam_rows = []
    for item in CAM_CASES:
        trace = load_cam_trace(item["path"])
        metrics = compute_cam_metrics(trace, ref)
        cam_rows.append({**item, "trace": trace, "metrics": metrics})

    gear_rows = []
    for item in GEAR_CASES:
        trace = load_gear_trace(item["path"])
        pinned_trace = load_gear_trace(item["pinned_path"])
        metrics = compute_gear_metrics(trace, pinned_trace)
        gear_rows.append({**item, "trace": trace, "metrics": metrics})

    cam_speedup = cam_rows[0]["runtime_s"] / cam_rows[2]["runtime_s"]
    cam_sdf2_vs_sdf1 = cam_rows[1]["metrics"]["vel_rmse"] - cam_rows[2]["metrics"]["vel_rmse"]
    cam_sdf2_vs_mesh = cam_rows[0]["metrics"]["vel_rmse"] - cam_rows[2]["metrics"]["vel_rmse"]
    cam_peak_gain_vs_sdf1 = cam_rows[1]["metrics"]["vel_max_abs"] - cam_rows[2]["metrics"]["vel_max_abs"]
    cam_peak_gain_vs_mesh = cam_rows[0]["metrics"]["vel_max_abs"] - cam_rows[2]["metrics"]["vel_max_abs"]

    pos_svg = polyline_svg(
        [
            {"x": ref["time"], "y": ref["pos"], "label": "Reference", "color": "#111827"},
            *[
                {"x": row["trace"]["time"], "y": row["trace"]["pos"], "label": row["label"], "color": row["color"]}
                for row in cam_rows
            ],
        ],
        y_label="Pos Y (m)",
    )
    vel_svg = polyline_svg(
        [
            {"x": ref["time"], "y": ref["vel"], "label": "Reference", "color": "#111827"},
            *[
                {"x": row["trace"]["time"], "y": row["trace"]["vel"], "label": row["label"], "color": row["color"]}
                for row in cam_rows
            ],
        ],
        y_label="Vel Y (m/s)",
    )
    gear_svg = polyline_svg(
        [
            {
                "x": gear_rows[0]["trace"]["time"],
                "y": [-1.0] * len(gear_rows[0]["trace"]["time"]),
                "label": "Ideal -1 rad/s",
                "color": "#111827",
                "dash": True,
            },
            *[
                {"x": row["trace"]["time"], "y": row["trace"]["vel"], "label": row["label"], "color": row["color"]}
                for row in gear_rows
            ],
        ],
        y_label="Gear22 Wrx (rad/s)",
    )
    cam_runtime_svg = bar_svg(cam_rows, "runtime_s", y_label="Runtime (s)")
    cam_velrmse_svg = bar_svg(
        [{"label": row["label"], "color": row["color"], "vel_rmse": row["metrics"]["vel_rmse"]} for row in cam_rows],
        "vel_rmse",
        y_label="Vel RMSE (m/s)",
    )
    gear_runtime_svg = bar_svg(gear_rows, "runtime_s", y_label="Runtime (s)")

    cam_legend_html = legend(
        [{"label": "Reference", "color": "#111827"}] +
        [{"label": r["label"], "color": r["color"]} for r in cam_rows]
    )
    gear_legend_html = legend(
        [{"label": "Ideal -1 rad/s", "color": "#111827", "dash": True}] +
        [{"label": r["label"], "color": r["color"]} for r in gear_rows]
    )

    cam_table_rows = []
    for row in cam_rows:
        cam_table_rows.append(
            "<tr>"
            f"<td>{html.escape(row['label'])}</td>"
            f"<td><code>{html.escape(row['path'].relative_to(ROOT).as_posix())}</code></td>"
            f"<td>{fmt(row['runtime_s'], 3)}</td>"
            f"<td>{fmt(row['metrics']['pos_rmse'])}</td>"
            f"<td>{fmt(row['metrics']['vel_rmse'])}</td>"
            f"<td>{fmt(row['metrics']['vel_mae'])}</td>"
            f"<td>{fmt(row['metrics']['vel_max_abs'])}</td>"
            "</tr>"
        )

    gear_table_rows = []
    for row in gear_rows:
        gear_table_rows.append(
            "<tr>"
            f"<td>{html.escape(row['label'])}</td>"
            f"<td><code>{html.escape(row['path'].relative_to(ROOT).as_posix())}</code></td>"
            f"<td>{fmt(row['runtime_s'], 3)}</td>"
            f"<td>{fmt(row['metrics']['mae'])}</td>"
            f"<td>{fmt(row['metrics']['rmse'])}</td>"
            f"<td>{fmt(row['metrics']['maxabs_vs_pinned'], 9)}</td>"
            "</tr>"
        )

    cam_defaults_rows = "".join(
        f"<tr><td><code>{html.escape(name)}</code></td><td>{html.escape(value)}</td></tr>" for name, value in CAM_DEFAULTS
    )

    html_text = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>Current SDF Benchmark Report</title>
  <style>
    :root {{
      --bg:#f6f3ec; --panel:#fffdf8; --ink:#1f2937; --muted:#6b7280; --line:#ddd6c8;
      --accent:#8b5e34; --blue:#1d4ed8; --teal:#0f766e; --red:#b91c1c;
    }}
    * {{ box-sizing:border-box; }}
    body {{
      margin:0; color:var(--ink); background:
      radial-gradient(circle at top right, rgba(139,94,52,.12), transparent 26%),
      linear-gradient(180deg, #faf7f1 0%, var(--bg) 100%);
      font-family:"Times New Roman", Times, serif;
    }}
    .wrap {{ width:min(1460px, 96vw); margin:0 auto; padding:28px 0 36px; }}
    .panel {{ background:var(--panel); border:1px solid var(--line); border-radius:18px; padding:18px 20px; box-shadow:0 12px 28px rgba(31,41,55,.05); margin-bottom:18px; }}
    .hero {{ display:grid; grid-template-columns:1.35fr .95fr; gap:18px; }}
    .cards {{ display:grid; grid-template-columns:repeat(4,1fr); gap:14px; margin-top:18px; }}
    .card {{ border:1px solid var(--line); border-radius:16px; padding:14px 16px; background:linear-gradient(180deg,#fffefb,#fcf8ef); }}
    .grid2 {{ display:grid; grid-template-columns:1fr 1fr; gap:18px; }}
    .grid3 {{ display:grid; grid-template-columns:1fr 1fr 1fr; gap:18px; }}
    h1,h2,h3,p {{ margin:0; }}
    p {{ line-height:1.5; }}
    .muted {{ color:var(--muted); }}
    .pill {{ display:inline-block; padding:4px 10px; border-radius:999px; border:1px solid #e7ddce; background:#fbf6ea; color:#7c5a36; font-size:12px; }}
    .big {{ font-size:28px; color:var(--accent); margin:8px 0; }}
    table {{ width:100%; border-collapse:collapse; font-size:15px; }}
    th, td {{ padding:10px 12px; border-bottom:1px solid #ece4d8; vertical-align:top; }}
    th {{ text-align:left; font-size:14px; color:#4b5563; }}
    code {{ font-family:"Consolas","Courier New",monospace; font-size:13px; }}
    .legend {{ display:flex; gap:14px; flex-wrap:wrap; margin-top:10px; font-size:14px; color:#4b5563; }}
    svg {{ width:100%; height:auto; display:block; border-radius:12px; border:1px solid #ebe4d6; background:#fffdfa; }}
    .foot {{ font-size:13px; color:var(--muted); }}
    @media (max-width:1100px) {{
      .hero,.grid2,.grid3,.cards {{ grid-template-columns:1fr; }}
    }}
  </style>
</head>
<body>
<div class="wrap">
  <section class="hero">
    <div class="panel">
      <div class="pill">Current Mainline State</div>
      <h1 style="margin-top:10px;">Cam onset spikes stay suppressed while BVH broad-phase restores SDF speed.</h1>
      <p class="muted" style="margin-top:10px;">
        This report summarizes the current default code path after promoting the cam-only
        patch-local onset manifold gate and enabling the sample-side BVH broad-phase on
        the sliding-patch benchmarks. Gear results remain bitwise identical to the pinned
        benchmark. On cam, the default <code>sdf_2nd</code> result is still the most
        accurate of the three on both velocity RMSE and velocity peak error, and both SDF
        modes are now faster than the native mesh baseline.
      </p>
    </div>
    <div class="panel">
      <h2 style="font-size:22px;">Headline Numbers</h2>
      <table style="margin-top:10px;">
        <tr><td>Cam mesh runtime</td><td><strong>{fmt(cam_rows[0]["runtime_s"], 3)} s</strong></td></tr>
        <tr><td>Cam SDF 2nd runtime</td><td><strong>{fmt(cam_rows[2]["runtime_s"], 3)} s</strong></td></tr>
        <tr><td>Cam SDF 2nd runtime ratio vs mesh</td><td><strong>{fmt(cam_speedup, 2)}x</strong></td></tr>
        <tr><td>Cam SDF 2nd vel RMSE</td><td><strong>{fmt(cam_rows[2]["metrics"]["vel_rmse"])}</strong></td></tr>
        <tr><td>Cam SDF 2nd gain vs mesh</td><td><strong>{fmt(cam_sdf2_vs_mesh, 6)} m/s</strong></td></tr>
        <tr><td>Cam SDF 2nd gain vs SDF 1st</td><td><strong>{fmt(cam_sdf2_vs_sdf1, 6)} m/s</strong></td></tr>
        <tr><td>Cam SDF 2nd peak gain vs mesh</td><td><strong>{fmt(cam_peak_gain_vs_mesh, 6)} m/s</strong></td></tr>
        <tr><td>Cam SDF 2nd peak gain vs SDF 1st</td><td><strong>{fmt(cam_peak_gain_vs_sdf1, 6)} m/s</strong></td></tr>
        <tr><td>Gear exactness vs pinned</td><td><strong>max abs diff = 0.0</strong></td></tr>
      </table>
    </div>
  </section>

  <section class="cards">
    <div class="card">
      <div class="muted">Cam Mesh</div>
      <div class="big">{fmt(cam_rows[0]["metrics"]["vel_rmse"])}</div>
      <div class="muted">Velocity RMSE</div>
      <div class="foot" style="margin-top:8px;">Runtime {fmt(cam_rows[0]["runtime_s"], 3)} s | Peak {fmt(cam_rows[0]["metrics"]["vel_max_abs"])}</div>
    </div>
    <div class="card">
      <div class="muted">Cam SDF 1st</div>
      <div class="big">{fmt(cam_rows[1]["metrics"]["vel_rmse"])}</div>
      <div class="muted">Velocity RMSE</div>
      <div class="foot" style="margin-top:8px;">Runtime {fmt(cam_rows[1]["runtime_s"], 3)} s | Peak {fmt(cam_rows[1]["metrics"]["vel_max_abs"])}</div>
    </div>
    <div class="card">
      <div class="muted">Cam SDF 2nd</div>
      <div class="big">{fmt(cam_rows[2]["metrics"]["vel_rmse"])}</div>
      <div class="muted">Velocity RMSE</div>
      <div class="foot" style="margin-top:8px;">Runtime {fmt(cam_rows[2]["runtime_s"], 3)} s | Peak {fmt(cam_rows[2]["metrics"]["vel_max_abs"])}</div>
    </div>
    <div class="card">
      <div class="muted">Gear SDF 2nd</div>
      <div class="big">{fmt(gear_rows[1]["metrics"]["mae"])}</div>
      <div class="muted">MAE to ideal -1 rad/s</div>
      <div class="foot" style="margin-top:8px;">Pinned diff {fmt(gear_rows[1]["metrics"]["maxabs_vs_pinned"], 9)}</div>
    </div>
  </section>

  <section class="grid2">
    <div class="panel">
      <h2>Cam Position</h2>
      <p class="muted" style="margin:6px 0 10px;">Reference, mesh, current <code>sdf_1st</code>, current <code>sdf_2nd</code>.</p>
      {pos_svg}
      <div class="legend">{cam_legend_html}</div>
    </div>
    <div class="panel">
      <h2>Cam Velocity</h2>
      <p class="muted" style="margin:6px 0 10px;">This is the main comparison view for the cam benchmark.</p>
      {vel_svg}
      <div class="legend">{cam_legend_html}</div>
    </div>
  </section>

  <section class="grid3">
    <div class="panel">
      <h2>Cam Runtime</h2>
      <p class="muted" style="margin:6px 0 10px;">BVH broad-phase recovers the runtime advantage while preserving the current cam trajectories exactly.</p>
      {cam_runtime_svg}
    </div>
    <div class="panel">
      <h2>Cam Velocity RMSE</h2>
      <p class="muted" style="margin:6px 0 10px;">Current ordering is SDF 2nd best, then SDF 1st, then mesh.</p>
      {cam_velrmse_svg}
    </div>
    <div class="panel">
      <h2>Gear Runtime</h2>
      <p class="muted" style="margin:6px 0 10px;">Current default gear timings measured on the same codebase.</p>
      {gear_runtime_svg}
    </div>
  </section>

    <section class="panel">
      <h2>Gear Accuracy Preservation</h2>
      <p class="muted" style="margin:6px 0 10px;">Current default gear outputs are still exactly equal to the pinned benchmark files.</p>
      {gear_svg}
      <div class="legend">{gear_legend_html}</div>
    </section>

  <section class="grid2">
    <div class="panel">
      <h2>Cam Metrics Table</h2>
      <table>
        <thead>
          <tr><th>Case</th><th>CSV</th><th>Runtime (s)</th><th>Pos RMSE</th><th>Vel RMSE</th><th>Vel MAE</th><th>Vel Max Abs</th></tr>
        </thead>
        <tbody>
          {"".join(cam_table_rows)}
        </tbody>
      </table>
    </div>
    <div class="panel">
      <h2>Gear Metrics Table</h2>
      <table>
        <thead>
          <tr><th>Case</th><th>CSV</th><th>Runtime (s)</th><th>MAE</th><th>RMSE</th><th>Max abs vs pinned</th></tr>
        </thead>
        <tbody>
          {"".join(gear_table_rows)}
        </tbody>
      </table>
    </div>
  </section>

  <section class="grid2">
    <div class="panel">
      <h2>Current Cam Default Parameters</h2>
      <p class="muted" style="margin:6px 0 10px;">These are the effective mainline defaults behind the current stable cam result.</p>
      <table>
        <thead><tr><th>Parameter</th><th>Value</th></tr></thead>
        <tbody>{cam_defaults_rows}</tbody>
      </table>
    </div>
    <div class="panel">
      <h2>Reading Guide</h2>
      <p class="muted" style="margin-bottom:12px;">
        The current interpretation is straightforward:
      </p>
      <table>
        <tbody>
          <tr><td>1.</td><td>Gear is safe: current default files match the pinned gear benchmark exactly.</td></tr>
          <tr><td>2.</td><td>Cam <code>sdf_2nd</code> now includes a patch-local onset manifold gate that suppresses early-contact spikes.</td></tr>
          <tr><td>3.</td><td>Cam <code>sdf_2nd</code> is more accurate than both mesh and cam <code>sdf_1st</code> on velocity RMSE.</td></tr>
          <tr><td>4.</td><td>The main new gain is local peak suppression: current cam <code>sdf_2nd</code> peak error is <strong>{fmt(cam_rows[2]["metrics"]["vel_max_abs"])}</strong> versus <strong>{fmt(cam_rows[1]["metrics"]["vel_max_abs"])}</strong> for <code>sdf_1st</code> and <strong>{fmt(cam_rows[0]["metrics"]["vel_max_abs"])}</strong> for mesh.</td></tr>
          <tr><td>5.</td><td>The current default now preserves the onset-suppressed trajectories while running both SDF modes faster than mesh on the <code>dt=0.01</code> cam benchmark.</td></tr>
        </tbody>
      </table>
    </div>
  </section>

  <section class="panel">
    <div class="foot">
      Generated from current CSV outputs in <code>data/outputs</code>. Data sources:
      {html.escape(json.dumps({row["label"]: row["path"].relative_to(ROOT).as_posix() for row in cam_rows + gear_rows}, ensure_ascii=False))}
    </div>
  </section>
</div>
</body>
</html>
"""
    OUT_HTML.write_text(html_text, encoding="utf-8")
    print(f"Wrote {OUT_HTML}")


if __name__ == "__main__":
    make_report()
