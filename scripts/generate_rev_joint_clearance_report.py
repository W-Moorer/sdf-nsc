from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
ASSET_DIR = ROOT / "assets" / "rev_joint_clearance"
OUTPUT_DIR = ROOT / "data" / "outputs"
DOC_DIR = ROOT / "docs"
FIG_DIR = DOC_DIR / "figures"


OURS = {
    "Mesh": OUTPUT_DIR / "rev_joint_clearance_mesh.csv",
    "SDF 1st": OUTPUT_DIR / "rev_joint_clearance_sdf1.csv",
    "SDF 2nd": OUTPUT_DIR / "rev_joint_clearance_sdf2.csv",
}
PLOT_METHODS = ("SDF 1st", "SDF 2nd")

REF_BODY2 = ASSET_DIR / "data" / "body2.csv"
REF_BODY3 = ASSET_DIR / "data" / "body3.csv"
REF_BODY2_IDEAL = ASSET_DIR / "data" / "body2_ideal.csv"
SUMMARY_JSON = ASSET_DIR / "models" / "extract_summary.json"
MODEL_JSON = ASSET_DIR / "rev_joint_clearance_model.json"
REPORT_HTML = DOC_DIR / "rev_joint_clearance_report.html"


plt.rcParams.update(
    {
        "font.family": "serif",
        "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
        "axes.labelsize": 11,
        "axes.titlesize": 12,
        "legend.fontsize": 9,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
    }
)


def _load_recurdyn_csv(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    time_col = next(c for c in df.columns if c.startswith("X:Pos_TX"))
    out = {"Time": df[time_col].to_numpy(float)}
    for col in df.columns:
        if col.startswith("Y:"):
            out[col[2:]] = df[col].to_numpy(float)
    return pd.DataFrame(out)


def _rmse(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.sqrt(np.mean((a - b) ** 2)))


def _mae(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.mean(np.abs(a - b)))


def _maxabs(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.max(np.abs(a - b)))


def _interp(ref_t: np.ndarray, ref_y: np.ndarray, t: np.ndarray) -> np.ndarray:
    return np.interp(t, ref_t, ref_y)


def _unwrap_angle_series(raw: np.ndarray) -> np.ndarray:
    return np.unwrap(raw)


def _derive_swing_angle_from_positions(
    body2_y: np.ndarray,
    body2_z: np.ndarray,
    body3_y: np.ndarray,
    body3_z: np.ndarray,
    offset_y: float,
    offset_z: float,
) -> np.ndarray:
    ref_angle = np.arctan2(offset_z, offset_y)
    cur_angle = np.arctan2(body2_z - body3_z, body2_y - body3_y)
    return _unwrap_angle_series(cur_angle - ref_angle)


def _central_diff(t: np.ndarray, y: np.ndarray) -> np.ndarray:
    return np.gradient(y, t)


def _plot_body2(
    body2_ref: pd.DataFrame,
    body2_ideal: pd.DataFrame,
    ours: dict[str, pd.DataFrame],
    out_path: Path,
) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(12.5, 7.6), constrained_layout=True)
    ref_t = body2_ref["Time"].to_numpy(float)
    ideal_t = body2_ideal["Time"].to_numpy(float)
    time = next(iter(ours.values()))["Time"].to_numpy(float)

    series = [
        ("Body2_Pos_TY", "Pos_TY-Body2-rev_clearance_joint(m)", "Pos_TY-Body2-rev_clearance_joint_ideal(m)", "Body2 TY Position (m)"),
        ("Body2_Pos_TZ", "Pos_TZ-Body2-rev_clearance_joint(m)", "Pos_TZ-Body2-rev_clearance_joint_ideal(m)", "Body2 TZ Position (m)"),
        ("Body2_Vel_TY", "Vel_TY-Body2-rev_clearance_joint(m/s)", "Vel_TY-Body2-rev_clearance_joint_ideal(m/s)", "Body2 TY Velocity (m/s)"),
        ("Body2_Vel_TZ", "Vel_TZ-Body2-rev_clearance_joint(m/s)", "Vel_TZ-Body2-rev_clearance_joint_ideal(m/s)", "Body2 TZ Velocity (m/s)"),
    ]

    colors = {"SDF 1st": "#d62728", "SDF 2nd": "#2ca02c"}
    for ax, (our_col, ref_col, ideal_col, ylabel) in zip(axes.ravel(), series):
        ax.plot(ref_t, body2_ref[ref_col].to_numpy(float), color="black", linewidth=1.6, label="RecurDyn clearance")
        ax.plot(
            ideal_t,
            body2_ideal[ideal_col].to_numpy(float),
            color="#555555",
            linestyle="--",
            linewidth=1.4,
            label="RecurDyn ideal revolute",
        )
        for name in PLOT_METHODS:
            df = ours[name]
            ax.plot(time, df[our_col].to_numpy(float), color=colors[name], linewidth=1.2, label=name)
        ax.set_xlabel("Time (s)")
        ax.set_ylabel(ylabel)
        ax.grid(True, linewidth=0.3, alpha=0.5)

    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=4, frameon=False, bbox_to_anchor=(0.5, 1.02))
    fig.savefig(out_path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def _plot_body3(body3_ref: pd.DataFrame, ours: dict[str, pd.DataFrame], out_path: Path) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(12.5, 7.6), constrained_layout=True)
    ref_t = body3_ref["Time"].to_numpy(float)
    time = next(iter(ours.values()))["Time"].to_numpy(float)
    series = [
        ("Body3_Pos_TY", "Pos_TY-Body3-rev_clearance_joint(m)", "Body3 TY Position (m)"),
        ("Body3_Pos_TZ", "Pos_TZ-Body3-rev_clearance_joint(m)", "Body3 TZ Position (m)"),
        ("Body3_Vel_TY", "Vel_TY-Body3-rev_clearance_joint(m/s)", "Body3 TY Velocity (m/s)"),
        ("Body3_Vel_TZ", "Vel_TZ-Body3-rev_clearance_joint(m/s)", "Body3 TZ Velocity (m/s)"),
    ]
    colors = {"SDF 1st": "#d62728", "SDF 2nd": "#2ca02c"}
    for ax, (our_col, ref_col, ylabel) in zip(axes.ravel(), series):
        ax.plot(ref_t, body3_ref[ref_col].to_numpy(float), color="black", linewidth=1.6, label="RecurDyn clearance")
        for name in PLOT_METHODS:
            df = ours[name]
            ax.plot(time, df[our_col].to_numpy(float), color=colors[name], linewidth=1.2, label=name)
        ax.set_xlabel("Time (s)")
        ax.set_ylabel(ylabel)
        ax.grid(True, linewidth=0.3, alpha=0.5)

    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=4, frameon=False, bbox_to_anchor=(0.5, 1.02))
    fig.savefig(out_path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def _plot_body3_rotation(
    t_ref: np.ndarray,
    angle_ref: np.ndarray,
    angvel_ref: np.ndarray,
    t_ideal: np.ndarray,
    angle_ideal: np.ndarray,
    angvel_ideal: np.ndarray,
    ours: dict[str, pd.DataFrame],
    out_path: Path,
) -> None:
    fig, axes = plt.subplots(2, 1, figsize=(12.5, 6.8), constrained_layout=True)
    colors = {"SDF 1st": "#d62728", "SDF 2nd": "#2ca02c"}
    time = next(iter(ours.values()))["Time"].to_numpy(float)

    axes[0].plot(t_ref, angle_ref, color="black", linewidth=1.6, label="RecurDyn clearance")
    axes[0].plot(t_ideal, angle_ideal, color="#555555", linestyle="--", linewidth=1.4, label="Ideal revolute")
    for name in PLOT_METHODS:
        df = ours[name]
        axes[0].plot(time, df["Body3_Ang_X"].to_numpy(float), color=colors[name], linewidth=1.2, label=name)
    axes[0].set_xlabel("Time (s)")
    axes[0].set_ylabel("Body3 Swing Angle about X (rad)")
    axes[0].grid(True, linewidth=0.3, alpha=0.5)

    axes[1].plot(t_ref, angvel_ref, color="black", linewidth=1.6, label="RecurDyn clearance")
    axes[1].plot(t_ideal, angvel_ideal, color="#555555", linestyle="--", linewidth=1.4, label="Ideal revolute")
    for name in PLOT_METHODS:
        df = ours[name]
        axes[1].plot(time, df["Body3_AngVel_X"].to_numpy(float), color=colors[name], linewidth=1.2, label=name)
    axes[1].set_xlabel("Time (s)")
    axes[1].set_ylabel("Body3 Angular Velocity about X (rad/s)")
    axes[1].grid(True, linewidth=0.3, alpha=0.5)

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=4, frameon=False, bbox_to_anchor=(0.5, 1.02))
    fig.savefig(out_path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def _collect_metrics(
    ours: dict[str, pd.DataFrame],
    body2_ref: pd.DataFrame,
    body3_ref: pd.DataFrame,
    body2_ideal: pd.DataFrame,
    body3_angle_ref: np.ndarray,
    body3_angvel_ref: np.ndarray,
    body3_angle_ideal: np.ndarray,
    body3_angvel_ideal: np.ndarray,
) -> dict[str, dict[str, dict[str, float]]]:
    out: dict[str, dict[str, dict[str, float]]] = {}
    ref_t = body2_ref["Time"].to_numpy(float)
    body3_t = body3_ref["Time"].to_numpy(float)
    ideal_t = body2_ideal["Time"].to_numpy(float)

    for name, df in ours.items():
        t = df["Time"].to_numpy(float)
        result: dict[str, dict[str, float]] = {}

        for title, ref_df, ref_time, cols in [
            (
                "Body2 vs RecurDyn clearance",
                body2_ref,
                ref_t,
                [
                    ("Body2_Pos_TY", "Pos_TY-Body2-rev_clearance_joint(m)"),
                    ("Body2_Pos_TZ", "Pos_TZ-Body2-rev_clearance_joint(m)"),
                    ("Body2_Vel_TY", "Vel_TY-Body2-rev_clearance_joint(m/s)"),
                    ("Body2_Vel_TZ", "Vel_TZ-Body2-rev_clearance_joint(m/s)"),
                ],
            ),
            (
                "Body2 vs ideal revolute",
                body2_ideal,
                ideal_t,
                [
                    ("Body2_Pos_TY", "Pos_TY-Body2-rev_clearance_joint_ideal(m)"),
                    ("Body2_Pos_TZ", "Pos_TZ-Body2-rev_clearance_joint_ideal(m)"),
                    ("Body2_Vel_TY", "Vel_TY-Body2-rev_clearance_joint_ideal(m/s)"),
                    ("Body2_Vel_TZ", "Vel_TZ-Body2-rev_clearance_joint_ideal(m/s)"),
                ],
            ),
            (
                "Body3 vs RecurDyn clearance",
                body3_ref,
                body3_t,
                [
                    ("Body3_Pos_TY", "Pos_TY-Body3-rev_clearance_joint(m)"),
                    ("Body3_Pos_TZ", "Pos_TZ-Body3-rev_clearance_joint(m)"),
                    ("Body3_Vel_TY", "Vel_TY-Body3-rev_clearance_joint(m/s)"),
                    ("Body3_Vel_TZ", "Vel_TZ-Body3-rev_clearance_joint(m/s)"),
                ],
            ),
        ]:
            section: dict[str, float] = {}
            for our_col, ref_col in cols:
                y = df[our_col].to_numpy(float)
                y_ref = _interp(ref_time, ref_df[ref_col].to_numpy(float), t)
                section[f"{our_col} RMSE"] = _rmse(y, y_ref)
                section[f"{our_col} MAE"] = _mae(y, y_ref)
                section[f"{our_col} MaxAbs"] = _maxabs(y, y_ref)
            result[title] = section

        rotation_clearance = {}
        y = df["Body3_Ang_X"].to_numpy(float)
        y_ref = _interp(body3_t, body3_angle_ref, t)
        rotation_clearance["Body3_Ang_X RMSE"] = _rmse(y, y_ref)
        rotation_clearance["Body3_Ang_X MAE"] = _mae(y, y_ref)
        rotation_clearance["Body3_Ang_X MaxAbs"] = _maxabs(y, y_ref)
        w = df["Body3_AngVel_X"].to_numpy(float)
        w_ref = _interp(body3_t, body3_angvel_ref, t)
        rotation_clearance["Body3_AngVel_X RMSE"] = _rmse(w, w_ref)
        rotation_clearance["Body3_AngVel_X MAE"] = _mae(w, w_ref)
        rotation_clearance["Body3_AngVel_X MaxAbs"] = _maxabs(w, w_ref)
        result["Body3 rotation vs RecurDyn clearance"] = rotation_clearance

        rotation_ideal = {}
        y_ideal = _interp(ideal_t, body3_angle_ideal, t)
        rotation_ideal["Body3_Ang_X RMSE"] = _rmse(y, y_ideal)
        rotation_ideal["Body3_Ang_X MAE"] = _mae(y, y_ideal)
        rotation_ideal["Body3_Ang_X MaxAbs"] = _maxabs(y, y_ideal)
        w_ideal = _interp(ideal_t, body3_angvel_ideal, t)
        rotation_ideal["Body3_AngVel_X RMSE"] = _rmse(w, w_ideal)
        rotation_ideal["Body3_AngVel_X MAE"] = _mae(w, w_ideal)
        rotation_ideal["Body3_AngVel_X MaxAbs"] = _maxabs(w, w_ideal)
        result["Body3 rotation vs ideal revolute"] = rotation_ideal
        out[name] = result

    return out


def _metrics_table_html(title: str, section: str, metrics: dict[str, dict[str, dict[str, float]]]) -> str:
    keys = [
        ("Body2_Pos_TY RMSE", "TY pos RMSE"),
        ("Body2_Pos_TZ RMSE", "TZ pos RMSE"),
        ("Body2_Vel_TY RMSE", "TY vel RMSE"),
        ("Body2_Vel_TZ RMSE", "TZ vel RMSE"),
    ]
    if "Body3" in section:
        keys = [
            ("Body3_Pos_TY RMSE", "TY pos RMSE"),
            ("Body3_Pos_TZ RMSE", "TZ pos RMSE"),
            ("Body3_Vel_TY RMSE", "TY vel RMSE"),
            ("Body3_Vel_TZ RMSE", "TZ vel RMSE"),
        ]
    if "rotation" in section:
        keys = [
            ("Body3_Ang_X RMSE", "Angle RMSE"),
            ("Body3_AngVel_X RMSE", "Angular vel RMSE"),
        ]

    rows = []
    for name in ("Mesh", "SDF 1st", "SDF 2nd"):
        row = f"<tr><th>{name}</th>"
        sec = metrics[name][section]
        for key, _label in keys:
            row += f"<td>{sec[key]:.6f}</td>"
        row += "</tr>"
        rows.append(row)

    head = "".join(f"<th>{label}</th>" for _key, label in keys)
    return f"""
    <section>
      <h3>{title}</h3>
      <table>
        <thead><tr><th>Method</th>{head}</tr></thead>
        <tbody>
          {''.join(rows)}
        </tbody>
      </table>
    </section>
    """


def main() -> None:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    DOC_DIR.mkdir(parents=True, exist_ok=True)

    ours = {name: pd.read_csv(path) for name, path in OURS.items()}
    body2_ref = _load_recurdyn_csv(REF_BODY2)
    body3_ref = _load_recurdyn_csv(REF_BODY3)
    body2_ideal = _load_recurdyn_csv(REF_BODY2_IDEAL)
    summary = json.loads(SUMMARY_JSON.read_text(encoding="utf-8"))
    model = json.loads(MODEL_JSON.read_text(encoding="utf-8"))
    offset_y = float(model["body2_cm_offset"][1])
    offset_z = float(model["body2_cm_offset"][2])

    body2_fig = FIG_DIR / "rev_joint_clearance_body2_compare.png"
    body3_fig = FIG_DIR / "rev_joint_clearance_body3_compare.png"
    body3_rot_fig = FIG_DIR / "rev_joint_clearance_body3_rotation_compare.png"
    _plot_body2(body2_ref, body2_ideal, ours, body2_fig)
    _plot_body3(body3_ref, ours, body3_fig)

    body3_ref_angle = _derive_swing_angle_from_positions(
        body2_ref["Pos_TY-Body2-rev_clearance_joint(m)"].to_numpy(float),
        body2_ref["Pos_TZ-Body2-rev_clearance_joint(m)"].to_numpy(float),
        body3_ref["Pos_TY-Body3-rev_clearance_joint(m)"].to_numpy(float),
        body3_ref["Pos_TZ-Body3-rev_clearance_joint(m)"].to_numpy(float),
        offset_y,
        offset_z,
    )
    body3_ref_angvel = _central_diff(body3_ref["Time"].to_numpy(float), body3_ref_angle)

    body3_ideal_angle = _derive_swing_angle_from_positions(
        body2_ideal["Pos_TY-Body2-rev_clearance_joint_ideal(m)"].to_numpy(float),
        body2_ideal["Pos_TZ-Body2-rev_clearance_joint_ideal(m)"].to_numpy(float),
        np.zeros(len(body2_ideal), dtype=float),
        np.zeros(len(body2_ideal), dtype=float),
        offset_y,
        offset_z,
    )
    body3_ideal_angvel = _central_diff(body2_ideal["Time"].to_numpy(float), body3_ideal_angle)

    _plot_body3_rotation(
        body3_ref["Time"].to_numpy(float),
        body3_ref_angle,
        body3_ref_angvel,
        body2_ideal["Time"].to_numpy(float),
        body3_ideal_angle,
        body3_ideal_angvel,
        ours,
        body3_rot_fig,
    )

    metrics = _collect_metrics(
        ours,
        body2_ref,
        body3_ref,
        body2_ideal,
        body3_ref_angle,
        body3_ref_angvel,
        body3_ideal_angle,
        body3_ideal_angvel,
    )

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8" />
  <title>Revolute Clearance Joint Comparison</title>
  <style>
    body {{
      margin: 24px auto;
      max-width: 1100px;
      font-family: "Times New Roman", Times, serif;
      color: #1b1b1b;
      line-height: 1.45;
      padding: 0 16px 40px;
    }}
    h1, h2, h3 {{ margin-bottom: 0.3em; }}
    p {{ margin-top: 0.3em; }}
    img {{ width: 100%; border: 1px solid #d8d8d8; }}
    .figure {{ margin: 18px 0 28px; }}
    table {{
      border-collapse: collapse;
      width: 100%;
      margin: 8px 0 24px;
      font-size: 15px;
    }}
    th, td {{
      border: 1px solid #cfcfcf;
      padding: 8px 10px;
      text-align: center;
    }}
    th {{
      background: #f2f2f2;
    }}
    .note {{
      background: #f7f8fb;
      border-left: 4px solid #8ca0d7;
      padding: 10px 14px;
      margin: 18px 0;
    }}
    code {{ font-family: Consolas, monospace; font-size: 0.95em; }}
  </style>
</head>
<body>
  <h1>Revolute Clearance Joint Benchmark</h1>
  <p>
    This report compares the extracted RecurDyn revolute-clearance model under the current project implementations
    (<code>mesh</code>, <code>sdf_1st</code>, and <code>sdf_2nd</code>) against the reference clearance trajectories
    (<code>body2.csv</code>, <code>body3.csv</code>) and the ideal revolute-joint reference
    (<code>body2_ideal.csv</code>).
  </p>
  <div class="note">
    <strong>Extracted OBJ assets</strong><br/>
    Body1 housing surface: <code>{summary["##GSURFACE##_Body1.Subtract1"]["obj"]}</code><br/>
    Body3 pin surface: <code>{summary["##GSURFACE##_Body3.Cylinder1"]["obj"]}</code>
  </div>

  <h2>Body2 Comparison</h2>
  <p>
    Body2 is the payload body rigidly attached to the clearance pin. The plot overlays the current project trajectories
    with both the RecurDyn clearance simulation and the ideal revolute-joint reference.
  </p>
  <div class="figure">
    <img src="figures/{body2_fig.name}" alt="Body2 comparison" />
  </div>

  <h2>Body3 Comparison</h2>
  <p>
    Body3 is the clearance pin. RecurDyn provides the clearance trajectory for this body; there is no separate ideal
    body3 reference because the ideal case replaces the clearance joint with a pure revolute joint.
  </p>
  <div class="figure">
    <img src="figures/{body3_fig.name}" alt="Body3 comparison" />
  </div>

  <h2>Body3 Rotation Response</h2>
  <p>
    Because Body2 is rigidly attached to Body3, the relative vector from Body3 to Body2.CM defines the effective
    swing angle of the pin-payload assembly. The curves below compare that inferred rotation from RecurDyn with the
    directly exported project-side angular response.
  </p>
  <div class="figure">
    <img src="figures/{body3_rot_fig.name}" alt="Body3 rotation comparison" />
  </div>

  <h2>Error Summary</h2>
  {_metrics_table_html("Body2 vs RecurDyn clearance", "Body2 vs RecurDyn clearance", metrics)}
  {_metrics_table_html("Body2 vs ideal revolute", "Body2 vs ideal revolute", metrics)}
  {_metrics_table_html("Body3 vs RecurDyn clearance", "Body3 vs RecurDyn clearance", metrics)}
  {_metrics_table_html("Body3 rotation vs RecurDyn clearance", "Body3 rotation vs RecurDyn clearance", metrics)}
  {_metrics_table_html("Body3 rotation vs ideal revolute", "Body3 rotation vs ideal revolute", metrics)}

  <h2>Current Reading</h2>
  <p>
    The position comparison shows that Body3 already tracks the RecurDyn clearance orbit much more closely than Body2.
    The added rotation-response comparison helps isolate the remaining source of Body2 error: if the effective swing
    angle of the Body3-Body2 assembly diverges from the RecurDyn clearance rotation, then Body2.CM will be displaced
    even when the Body3 translation itself is already close.
  </p>
  <p>
    This means the benchmark is executable end-to-end and the extracted OBJ assets are usable, but the remaining
    mismatch should now be interpreted primarily as a rotational-response mismatch of the pin-payload assembly rather
    than a simple output-point or mesh-registration error.
  </p>
</body>
</html>
"""
    REPORT_HTML.write_text(html, encoding="utf-8")
    print(f"Wrote {body2_fig}")
    print(f"Wrote {body3_fig}")
    print(f"Wrote {body3_rot_fig}")
    print(f"Wrote {REPORT_HTML}")


if __name__ == "__main__":
    main()
