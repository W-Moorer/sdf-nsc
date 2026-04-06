from __future__ import annotations

import json
import math
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
OUT_HTML = ROOT / "docs" / "cam_sdf_analysis_report.html"


def load_reference():
    ref = pd.read_csv(ROOT / "assets" / "cam" / "data" / "cam_data.csv")
    return {
        "time": ref.iloc[:, 0].tolist(),
        "pos": ref.iloc[:, 10].tolist(),
        "vel": ref.iloc[:, 13].tolist(),
    }


def load_trace(path: Path):
    df = pd.read_csv(path)
    cols = set(df.columns)
    if "Y:Pos_TY-Body2(m)" in cols:
        time_col = df.columns[0]
        pos_col = "Y:Pos_TY-Body2(m)"
        vel_col = "Y:Vel_TY-Body2(m/s)"
    elif {"t", "py", "vy"}.issubset(cols):
        time_col = "t"
        pos_col = "py"
        vel_col = "vy"
    else:
        raise RuntimeError(f"Unsupported columns in {path}: {list(df.columns)}")
    return {
        "time": df[time_col].tolist(),
        "pos": df[pos_col].tolist(),
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


def compute_metrics(trace, ref):
    rp = interp_linear(ref["time"], ref["pos"], trace["time"])
    rv = interp_linear(ref["time"], ref["vel"], trace["time"])
    pos_err = [a - b for a, b in zip(trace["pos"], rp)]
    vel_err = [a - b for a, b in zip(trace["vel"], rv)]

    def mae(arr):
        return sum(abs(v) for v in arr) / len(arr)

    def rmse(arr):
        return math.sqrt(sum(v * v for v in arr) / len(arr))

    return {
        "pos_mae": mae(pos_err),
        "pos_rmse": rmse(pos_err),
        "vel_mae": mae(vel_err),
        "vel_rmse": rmse(vel_err),
    }


def obj_bounds(path: Path):
    mins = [float("inf")] * 3
    maxs = [float("-inf")] * 3
    vertices = 0
    with path.open("r", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if not line.startswith("v "):
                continue
            parts = line.split()
            xyz = [float(parts[1]), float(parts[2]), float(parts[3])]
            for i, v in enumerate(xyz):
                mins[i] = min(mins[i], v)
                maxs[i] = max(maxs[i], v)
            vertices += 1
    size = [maxs[i] - mins[i] for i in range(3)]
    return {"vertices": vertices, "mins": mins, "maxs": maxs, "size": size}


def estimate_cam_follower_samples():
    path = ROOT / "assets" / "cam" / "models" / "cam_body2.obj"
    vertices = []
    faces = []
    with path.open("r", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if line.startswith("v "):
                p = line.split()
                vertices.append((float(p[1]), float(p[2]), float(p[3])))
            elif line.startswith("f "):
                idx = []
                for token in line.split()[1:4]:
                    idx.append(int(token.split("/")[0]) - 1)
                faces.append(tuple(idx))

    def sub(a, b):
        return (a[0] - b[0], a[1] - b[1], a[2] - b[2])

    def cross(a, b):
        return (
            a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0],
        )

    def norm(a):
        return math.sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2])

    out = []
    for res in [5e-4, 7.5e-4, 1e-3, 1.25e-3, 1.5e-3, 2e-3]:
        samples = len(vertices)
        for i, j, k in faces:
            v0, v1, v2 = vertices[i], vertices[j], vertices[k]
            area = 0.5 * norm(cross(sub(v1, v0), sub(v2, v0)))
            target = min(500, max(0, int(area / (res * res * 0.5))))
            if target > 0:
                steps = max(2, int(math.sqrt(target * 2)))
                samples += (steps - 1) * (steps - 2) // 2
        out.append({"surface_res": res, "samples": samples})
    return out


def make_report():
    ref = load_reference()
    traces = [
        {
            "key": "reference",
            "label": "Reference",
            "kind": "reference",
            "path": "assets/cam/data/cam_data.csv",
            "color": "#0f172a",
            "trace": ref,
            "metrics": {"pos_mae": 0.0, "pos_rmse": 0.0, "vel_mae": 0.0, "vel_rmse": 0.0},
            "runtime_s": None,
        },
        {
            "key": "hist_sdf2",
            "label": "Historical SDF 2nd dt=0.01",
            "kind": "historical",
            "path": "data/outputs/baseline_cam_nsc_hessian_dt0.01.csv",
            "color": "#16a34a",
            "runtime_s": None,
        },
        {
            "key": "current_mesh",
            "label": "Current Mesh dt=0.01",
            "kind": "current_mesh",
            "path": "data/outputs/cam_reverify_mesh__dt0p01.csv",
            "color": "#2563eb",
            "runtime_s": 120.34,
        },
        {
            "key": "best_current_sdf2",
            "label": "Current best stable SDF 2nd",
            "kind": "current_sdf2",
            "path": "data/outputs/cam_sweep_f.csv",
            "color": "#dc2626",
            "runtime_s": 31.26,
        },
        {
            "key": "current_sdf1_probe",
            "label": "Current SDF 1st probe",
            "kind": "current_sdf1",
            "path": "data/outputs/cam_probe_sdf1_best_dt0p01.csv",
            "color": "#a855f7",
            "runtime_s": 34.56,
        },
    ]

    for item in traces[1:]:
        item["trace"] = load_trace(ROOT / item["path"])
        item["metrics"] = compute_metrics(item["trace"], ref)

    sweep_runtime = {
        "cam_sweep_a.csv": 21.32,
        "cam_sweep_b.csv": 16.29,
        "cam_sweep_c.csv": 36.63,
        "cam_sweep_d.csv": 13.13,
        "cam_sweep_e.csv": 82.63,
        "cam_sweep_f.csv": 31.26,
        "cam_sweep_g.csv": 32.56,
        "cam_sweep_h.csv": 66.96,
        "cam_sweep_i.csv": 74.75,
        "cam_sweep_j.csv": 49.83,
    }
    sweep_points = []
    for fname, runtime in sweep_runtime.items():
        trace = load_trace(ROOT / "data" / "outputs" / fname)
        metrics = compute_metrics(trace, ref)
        sweep_points.append(
            {
                "name": fname.replace("cam_sweep_", "").replace(".csv", ""),
                "runtime_s": runtime,
                "pos_rmse": metrics["pos_rmse"],
                "vel_rmse": metrics["vel_rmse"],
            }
        )

    param_rows = [
        {
            "name": "SPCC_CAM_VOXEL_SIZE",
            "current": "1e-4",
            "recommended": "3e-4 to 5e-4, start at 4e-4",
            "why": "Current cam body is about 0.42 m wide. 1e-4 is memory-heavy; 4e-4 was also the paper-scale value for the successful cam story.",
        },
        {
            "name": "SPCC_CAM_HALF_BAND_VOXELS",
            "current": "50",
            "recommended": "12 to 16 when voxel_size=4e-4",
            "why": "Keeps the physical half-band near 4.8 to 6.4 mm, close to the current 5 mm, but with far fewer active voxels.",
        },
        {
            "name": "SPCC_CAM_SURFACE_RES",
            "current": "5e-4",
            "recommended": "1.0e-3 to 1.25e-3",
            "why": "At 5e-4 the follower produces about 447,936 samples. At 1.0e-3 it drops to about 72,192; at 1.25e-3 it is about 38,976.",
        },
        {
            "name": "SPCC_CAM_MAX_SAMPLES",
            "current": "0 (uncapped)",
            "recommended": "4e4 to 8e4, start at 6e4",
            "why": "Cam should not run uncapped while gear is explicitly capped. This is the simplest runtime guardrail.",
        },
        {
            "name": "SPCC_CAM_FULL_SCAN_PERIOD",
            "current": "10",
            "recommended": "1 to 2",
            "why": "The cam rotates at about pi rad/s, so a contact point around r~0.2 m travels roughly 6 mm per 0.01 s step. Ten-step reuse is too stale.",
        },
        {
            "name": "SPCC_CAM_LOCAL_SCAN_RADIUS",
            "current": "8e-4",
            "recommended": "8e-3 to 1.2e-2, start at 1e-2",
            "why": "0.8 mm is gear-scale. Cam contact migration per step is already on the order of millimeters.",
        },
        {
            "name": "SPCC_CAM_CLUSTER_RADIUS",
            "current": "3e-4",
            "recommended": "2e-3 to 4e-3, start at 2.5e-3",
            "why": "0.3 mm over-fragments an extended sliding contact manifold. Cam needs patch grouping at the millimeter scale.",
        },
        {
            "name": "SPCC_CAM_CLUSTER_ANGLE_DEG",
            "current": "30",
            "recommended": "20 to 35, start at 25 or 30",
            "why": "This knob is not the main failure source. Keep it moderate and tune only after scan radius and sampling are fixed.",
        },
        {
            "name": "SPCC_CAM_SEPARATING_CUTOFF",
            "current": "1e-4",
            "recommended": "5e-4 to 1e-3",
            "why": "Current cutoff is too permissive for cam follower velocities that reach roughly 1.48 m/s in magnitude.",
        },
        {
            "name": "SPCC_CAM_AVG_POINT",
            "current": "0",
            "recommended": "1",
            "why": "Cam behaves more like sliding distributed contact than the compact gear tooth patch. Averaging the patch point is more defensible here.",
        },
        {
            "name": "SPCC_CAM_DELTA_ON / DELTA_OFF",
            "current": "4e-3 / 5e-3",
            "recommended": "Keep in the 1e-3 range, start with current values",
            "why": "These are not the dominant failure knobs. The bigger issues are scan staleness, sample count, and manifold compression.",
        },
        {
            "name": "SPCC_CAM_HOLD_STEPS",
            "current": "10",
            "recommended": "1 to 3, start at 2",
            "why": "Long hold makes stale patches linger while the cam keeps sliding underneath them.",
        },
        {
            "name": "SPCC_CAM_MAX_CONTACTS",
            "current": "4",
            "recommended": "8 to 12 after the cap is really enforced",
            "why": "Right now this knob is effectively dead: it is stored in Configure(), but BuildActiveSet() never trims to max_active_keep_.",
        },
        {
            "name": "SPCC_CAM_DIRECT_PHI_HESSIAN",
            "current": "0",
            "recommended": "Keep at 0",
            "why": "The current evidence points to activation/manifold issues first, not to a need for direct phi Hessians.",
        },
    ]

    summary_cards = [
        {
            "title": "Historical good SDF 2nd",
            "subtitle": "baseline_cam_nsc_hessian_dt0.01.csv",
            "value": "vel RMSE 0.0233",
            "detail": "pos RMSE 0.00111, shows the method previously worked well on cam.",
        },
        {
            "title": "Current mesh reference",
            "subtitle": "cam_reverify_mesh__dt0p01.csv",
            "value": "vel RMSE 0.0576",
            "detail": "Stable and reproducible. This is the current trustworthy cam baseline.",
        },
        {
            "title": "Current best stable SDF 2nd",
            "subtitle": "cam_sweep_f.csv",
            "value": "vel RMSE 0.5013",
            "detail": "Fast enough to run, but far worse than both historical SDF2 and current mesh.",
        },
        {
            "title": "Current SDF 1st probe",
            "subtitle": "cam_probe_sdf1_best_dt0p01.csv",
            "value": "vel RMSE 1.8951",
            "detail": "Shows the current cam path is not just a Hessian issue; manifold activation is already off.",
        },
    ]

    root_causes = [
        {
            "title": "Cam defaults moved into a gear-sized local scan regime",
            "detail": "Current local scan defaults are FULL_SCAN_PERIOD=10 and LOCAL_SCAN_RADIUS=8e-4. That makes sense for the small gear, but not for a cam whose contact patch can move about 6 mm in one 0.01 s step.",
            "refs": [
                "src/backend/chrono/spcc/ContactActivation.cpp:224-225",
                "assets/cam/cam_model.json: speed_rad_s = 3.1415926",
            ],
        },
        {
            "title": "Cam still builds a very heavy SDF and a huge follower sample cloud",
            "detail": "Current cam defaults are voxel_size=1e-4, half_band=50, surface_res=5e-4, and max_samples=0. That is a large narrow-band build plus an uncapped follower cloud of about 447,936 samples.",
            "refs": [
                "src/backend/chrono/ChronoRigidSystemNSC.cpp:464-487",
                "src/backend/chrono/ChronoRigidSystemNSC.cpp:44-114",
            ],
        },
        {
            "title": "Patch geometry is now compressed too aggressively for cam",
            "detail": "Current BuildActiveSet first queries only phi+grad, then clusters candidates, then asks for Hessian only on the deepest sample. That is efficient for gear, but cam uses extended sliding contact and loses too much manifold information under this compression.",
            "refs": [
                "src/backend/chrono/spcc/ContactActivation.cpp:287",
                "src/backend/chrono/spcc/ContactActivation.cpp:420-428",
                "src/backend/chrono/spcc/ContactActivation.cpp:458-466",
            ],
        },
        {
            "title": "The patch rules favor deepest-point contact, not sliding-contact continuity",
            "detail": "Normals are averaged, but phi and Hessian stay attached to the deepest sample. That works on compact tooth contacts more often than on cam-style rolling/sliding manifolds.",
            "refs": [
                "src/backend/chrono/spcc/ContactActivation.cpp:397-418",
                "src/backend/chrono/spcc/ContactActivation.cpp:420-434",
            ],
        },
        {
            "title": "One declared cam knob is currently not honored",
            "detail": "SPCC_CAM_MAX_CONTACTS is passed into Configure(), but BuildActiveSet never trims output to max_active_keep_. The code stores the number but does not enforce it.",
            "refs": [
                "src/backend/chrono/spcc/ContactActivation.cpp:123-130",
                "src/backend/chrono/spcc/ContactActivation.cpp:480-481",
            ],
        },
    ]

    return {
        "traces": traces,
        "sweep_points": sweep_points,
        "summary_cards": summary_cards,
        "param_rows": param_rows,
        "root_causes": root_causes,
        "sample_estimates": estimate_cam_follower_samples(),
        "cam_bounds": obj_bounds(ROOT / "assets" / "cam" / "models" / "cam_body1.obj"),
        "follower_bounds": obj_bounds(ROOT / "assets" / "cam" / "models" / "cam_body2.obj"),
        "gear_bounds": obj_bounds(ROOT / "assets" / "simple_gear" / "model" / "GEAR21.obj"),
    }


HTML_HEAD = """<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>Cam SDF Analysis Report</title>
<style>
:root{--bg:#f5f2e8;--panel:#fffdf7;--ink:#1f2937;--muted:#6b7280;--line:#d6d0c3;--accent:#8b5e34}
*{box-sizing:border-box}body{margin:0;font-family:"Times New Roman",Times,serif;color:var(--ink);background:radial-gradient(circle at top right,rgba(139,94,52,.10),transparent 28%),linear-gradient(180deg,#f8f4ea 0%,var(--bg) 100%)}
.wrap{width:min(1400px,96vw);margin:0 auto;padding:28px 0 40px}.panel{background:var(--panel);border:1px solid var(--line);border-radius:18px;padding:18px 20px;box-shadow:0 12px 28px rgba(31,41,55,.06);margin-bottom:18px}
.hero{display:grid;grid-template-columns:1.35fr .95fr;gap:18px}.cards{display:grid;grid-template-columns:repeat(4,1fr);gap:14px}.card{padding:14px 16px;border-radius:16px;border:1px solid var(--line);background:linear-gradient(180deg,#fffefb,#fcf8ef)}
.grid2{display:grid;grid-template-columns:1fr 1fr;gap:18px}.plot-head{display:flex;justify-content:space-between;align-items:baseline;gap:12px;margin-bottom:10px}.plot-sub{color:var(--muted);font-size:14px}
h1,h2,h3,p{margin:0}p{line-height:1.5}svg{width:100%;height:auto;display:block;border-radius:12px;background:#fffdfa;border:1px solid #ebe4d6}.legend{display:flex;gap:14px;flex-wrap:wrap;margin-top:10px;font-size:14px}.legend span{display:inline-flex;align-items:center;gap:8px}.sw{width:18px;height:3px;border-radius:999px;display:inline-block}
table{width:100%;border-collapse:collapse;font-size:15px}th,td{padding:10px 12px;border-bottom:1px solid #e9e1d4;vertical-align:top}th{text-align:left;color:#4b5563;font-size:14px}
.pill{display:inline-block;padding:4px 10px;border-radius:999px;font-size:12px;border:1px solid #e5dccd;background:#fbf6ea;color:#7c5a36}.note,.refs,.muted{color:var(--muted);font-size:14px}.cause{padding:14px 0;border-bottom:1px solid #e9e1d4}.cause:last-child{border-bottom:0}.cause h3{font-size:19px;margin-bottom:6px}
code{font-family:"Consolas","Courier New",monospace}@media (max-width:1100px){.hero,.grid2,.cards{grid-template-columns:1fr}}
</style></head><body><div class="wrap">"""


HTML_TAIL = """</div><script>
const DATA=__DATA__;
function legend(id,traces){document.getElementById(id).innerHTML=traces.map(t=>`<span><i class="sw" style="background:${t.color}"></i>${t.label}</span>`).join('')}
function linePlot(id,traces,key,yLabel){const svg=document.getElementById(id),W=640,H=360,m={l:58,r:18,t:16,b:42};const xs=traces.flatMap(t=>t.trace.time),ys=traces.flatMap(t=>t.trace[key]);const xmin=Math.min(...xs),xmax=Math.max(...xs),ymin=Math.min(...ys),ymax=Math.max(...ys),pad=(ymax-ymin||1)*.08,y0=ymin-pad,y1=ymax+pad;const sx=x=>m.l+(x-xmin)/(xmax-xmin)*(W-m.l-m.r),sy=y=>H-m.b-(y-y0)/(y1-y0)*(H-m.t-m.b);let s=[];for(let i=0;i<=6;i++){let xv=xmin+(xmax-xmin)*i/6,x=sx(xv);s.push(`<line x1="${x}" y1="${m.t}" x2="${x}" y2="${H-m.b}" stroke="#ece5d8"/>`);s.push(`<text x="${x}" y="${H-14}" text-anchor="middle" font-size="12" fill="#6b7280">${xv.toFixed(1)}</text>`)}for(let i=0;i<=5;i++){let yv=y0+(y1-y0)*i/5,y=sy(yv);s.push(`<line x1="${m.l}" y1="${y}" x2="${W-m.r}" y2="${y}" stroke="#ece5d8"/>`);s.push(`<text x="${m.l-8}" y="${y+4}" text-anchor="end" font-size="12" fill="#6b7280">${yv.toFixed(3)}</text>`)}s.push(`<line x1="${m.l}" y1="${m.t}" x2="${m.l}" y2="${H-m.b}" stroke="#8b8171"/>`);s.push(`<line x1="${m.l}" y1="${H-m.b}" x2="${W-m.r}" y2="${H-m.b}" stroke="#8b8171"/>`);s.push(`<text x="${W/2}" y="${H-2}" text-anchor="middle" font-size="13" fill="#6b7280">time (s)</text>`);s.push(`<text x="18" y="${H/2}" text-anchor="middle" font-size="13" fill="#6b7280" transform="rotate(-90 18 ${H/2})">${yLabel}</text>`);traces.forEach(t=>{const d=t.trace.time.map((x,i)=>`${i===0?'M':'L'} ${sx(x).toFixed(2)} ${sy(t.trace[key][i]).toFixed(2)}`).join(' ');s.push(`<path d="${d}" fill="none" stroke="${t.color}" stroke-width="${t.kind==='reference'?2.8:2.2}" stroke-linecap="round" stroke-linejoin="round"/>`)});svg.innerHTML=s.join('')}
function scatter(id,pts){const svg=document.getElementById(id),W=760,H=360,m={l:62,r:18,t:16,b:44};const xs=pts.map(p=>p.runtime_s),ys=pts.map(p=>p.vel_rmse),xmin=0,xmax=Math.max(...xs)*1.08,ymin=0,ymax=Math.max(...ys)*1.08,sx=x=>m.l+(x-xmin)/(xmax-xmin)*(W-m.l-m.r),sy=y=>H-m.b-(y-ymin)/(ymax-ymin)*(H-m.t-m.b);let s=[];for(let i=0;i<=6;i++){let xv=xmin+(xmax-xmin)*i/6,x=sx(xv);s.push(`<line x1="${x}" y1="${m.t}" x2="${x}" y2="${H-m.b}" stroke="#ece5d8"/>`);s.push(`<text x="${x}" y="${H-14}" text-anchor="middle" font-size="12" fill="#6b7280">${xv.toFixed(0)}</text>`)}for(let i=0;i<=5;i++){let yv=ymin+(ymax-ymin)*i/5,y=sy(yv);s.push(`<line x1="${m.l}" y1="${y}" x2="${W-m.r}" y2="${y}" stroke="#ece5d8"/>`);s.push(`<text x="${m.l-8}" y="${y+4}" text-anchor="end" font-size="12" fill="#6b7280">${yv.toFixed(2)}</text>`)}s.push(`<line x1="${m.l}" y1="${m.t}" x2="${m.l}" y2="${H-m.b}" stroke="#8b8171"/>`);s.push(`<line x1="${m.l}" y1="${H-m.b}" x2="${W-m.r}" y2="${H-m.b}" stroke="#8b8171"/>`);s.push(`<text x="${W/2}" y="${H-2}" text-anchor="middle" font-size="13" fill="#6b7280">runtime (s)</text>`);s.push(`<text x="18" y="${H/2}" text-anchor="middle" font-size="13" fill="#6b7280" transform="rotate(-90 18 ${H/2})">velocity RMSE</text>`);pts.forEach(p=>{const x=sx(p.runtime_s),y=sy(p.vel_rmse);s.push(`<circle cx="${x}" cy="${y}" r="5.5" fill="#8b5e34" stroke="#fffdfa" stroke-width="2"/>`);s.push(`<text x="${x+8}" y="${y-8}" font-size="12" fill="#4b5563">${p.name}</text>`)});svg.innerHTML=s.join('')}
linePlot('plot-pos',DATA.traces,'pos','Y:Pos_TY (m)');linePlot('plot-vel',DATA.traces,'vel','Y:Vel_TY (m/s)');legend('legend-pos',DATA.traces);legend('legend-vel',DATA.traces);scatter('plot-scatter',DATA.sweep_points);
</script></body></html>"""


def build_html(report):
    cards_html = "".join(
        f"""
        <div class="card">
          <div class="muted" style="margin-bottom:8px;">{c['title']}</div>
          <div style="font-size:28px; color:#8b5e34; margin-bottom:8px;">{c['value']}</div>
          <div class="muted" style="margin-bottom:6px;">{c['subtitle']}</div>
          <div style="font-size:14px; color:#4b5563;">{c['detail']}</div>
        </div>
        """
        for c in report["summary_cards"]
    )
    param_rows_html = "".join(
        f"<tr><td><code>{r['name']}</code></td><td><code>{r['current']}</code></td><td><strong>{r['recommended']}</strong></td><td>{r['why']}</td></tr>"
        for r in report["param_rows"]
    )
    sample_rows_html = "".join(
        f"<tr><td><code>{row['surface_res']:.6f}</code></td><td>{row['samples']:,}</td></tr>"
        for row in report["sample_estimates"]
    )
    root_causes_html = "".join(
        f"""
        <div class="cause">
          <h3>{cause['title']}</h3>
          <p>{cause['detail']}</p>
          <div class="refs">{'<br>'.join(cause['refs'])}</div>
        </div>
        """
        for cause in report["root_causes"]
    )

    cam_bbox = " x ".join(f"{v:.3f} m" for v in report["cam_bounds"]["size"])
    follower_bbox = " x ".join(f"{v:.3f} m" for v in report["follower_bounds"]["size"])
    gear_bbox = " x ".join(f"{v*1e-3:.4f} m" for v in report["gear_bounds"]["size"])

    body = f"""
    <section class="hero">
      <div class="panel">
        <div class="pill">Cam SDF report</div>
        <h1 style="margin-top:10px;">Cam defaults should move away from gear-scale heuristics</h1>
        <p style="color:#6b7280; margin-top:10px;">
          This report only uses the current source tree and existing CSV outputs. The main conclusion is that the current
          cam failure is not explained by the Hessian formula alone. The larger problem is that the activation and manifold
          logic was optimized around the small gear benchmark, while cam operates at a much larger spatial and kinematic scale.
        </p>
      </div>
      <div class="panel">
        <h2 style="font-size:22px;">Executive summary</h2>
        <p style="margin-top:10px; color:#6b7280;">Historical cam SDF 2nd-order was strong, but the current code path is now biased toward compact gear patches and stale local scans.</p>
        <p style="margin-top:10px; color:#6b7280;">The most important default shifts are voxel_size from <code>1e-4</code> toward <code>4e-4</code>, surface_res from <code>5e-4</code> toward <code>1e-3</code>, FULL_SCAN_PERIOD from <code>10</code> toward <code>1-2</code>, and LOCAL_SCAN_RADIUS from <code>8e-4</code> toward <code>1e-2</code>.</p>
      </div>
    </section>
    <section class="cards">{cards_html}</section>
    <section class="grid2">
      <div class="panel"><div class="plot-head"><h2>Follower Position</h2><div class="plot-sub">Y:Pos_TY against reference</div></div><svg id="plot-pos" viewBox="0 0 640 360"></svg><div class="legend" id="legend-pos"></div></div>
      <div class="panel"><div class="plot-head"><h2>Follower Velocity</h2><div class="plot-sub">Y:Vel_TY against reference</div></div><svg id="plot-vel" viewBox="0 0 640 360"></svg><div class="legend" id="legend-vel"></div></div>
    </section>
    <section class="panel"><div class="plot-head"><h2>Runtime vs Velocity Error</h2><div class="plot-sub">Current SDF 2nd sweep points, dt=0.01</div></div><svg id="plot-scatter" viewBox="0 0 760 360"></svg><div class="note" style="margin-top:10px;">Historical good SDF 2nd is excluded from the scatter because only the output CSV remains; its wall time was not recovered from current logs.</div></section>
    <section class="grid2">
      <div class="panel"><div class="plot-head"><h2>Recommended SPCC_CAM_* default magnitudes</h2><div class="plot-sub">Analysis recommendations, not yet re-validated as a final tuned preset</div></div><table><thead><tr><th>Variable</th><th>Current</th><th>Recommended start</th><th>Why</th></tr></thead><tbody>{param_rows_html}</tbody></table></div>
      <div class="panel"><div class="plot-head"><h2>Geometry and sampling scale</h2><div class="plot-sub">Why gear defaults do not transfer</div></div>
      <table><thead><tr><th>Item</th><th>Value</th><th>Interpretation</th></tr></thead><tbody>
      <tr><td>Cam bbox</td><td>{cam_bbox}</td><td>About 0.42 m wide. This is a large object for a 1e-4 SDF.</td></tr>
      <tr><td>Follower bbox</td><td>{follower_bbox}</td><td>About 0.24 m wide, so the contact manifold is not a tiny point feature.</td></tr>
      <tr><td>Gear bbox after 1e-3 scale</td><td>{gear_bbox}</td><td>Only millimeter scale. This is why tight local radii worked there.</td></tr>
      <tr><td>Cam physical half-band now</td><td>50 * 1e-4 = 5 mm</td><td>Reasonable thickness, but at a very expensive voxel count.</td></tr>
      <tr><td>Suggested half-band</td><td>12 * 4e-4 = 4.8 mm</td><td>Almost the same physical band, much cheaper build.</td></tr>
      <tr><td>Contact travel per step</td><td>about 6 mm at dt=0.01</td><td>Based on r~0.2 m and omega~pi rad/s. Much larger than the current 0.8 mm local scan radius.</td></tr>
      </tbody></table><div style="height:14px"></div><table><thead><tr><th>Follower surface_res</th><th>Estimated samples</th></tr></thead><tbody>{sample_rows_html}</tbody></table></div>
    </section>
    <section class="panel"><div class="plot-head"><h2>What likely made cam worse in the code</h2><div class="plot-sub">Compared against the earlier project import state and the surviving historical good CSV</div></div>{root_causes_html}</section>
    """
    return (HTML_HEAD + body + HTML_TAIL).replace(
        "__DATA__", json.dumps({"traces": report["traces"], "sweep_points": report["sweep_points"]}, ensure_ascii=False)
    )


def main():
    report = make_report()
    OUT_HTML.parent.mkdir(parents=True, exist_ok=True)
    OUT_HTML.write_text(build_html(report), encoding="utf-8")
    print(f"Wrote {OUT_HTML}")


if __name__ == "__main__":
    main()
