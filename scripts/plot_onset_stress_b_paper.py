import csv
import json
import math
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parents[1]
OUTPUT_DIR = ROOT / "data" / "outputs"
FIG_DIR = ROOT / "papers" / "figures"


def load_series(path: Path):
    times, pos, vel = [], [], []
    with path.open("r", encoding="utf-8") as f:
        reader = csv.reader(f)
        next(reader)
        for row in reader:
            if not row:
                continue
            times.append(float(row[0]))
            pos.append(float(row[10]))
            vel.append(float(row[13]))
    return times, pos, vel


def piecewise_reference(times, params):
    radius = params["cam_radius"]
    ecc = params["cam_eccentricity"]
    roller = params["roller_radius"]
    omega = params["motor_speed"]
    phase = params["phase"]
    y0 = params["follower_init_pos"][1]
    reach = radius + roller

    pos_ref, vel_ref = [], []
    for t in times:
        theta = phase + omega * t
        cx = ecc * math.cos(theta)
        cy = ecc * math.sin(theta)
        root = math.sqrt(max(reach * reach - cx * cx, 0.0))
        env_pos = cy + root
        dpos_dtheta = ecc * math.cos(theta) + (ecc * ecc * math.cos(theta) * math.sin(theta)) / max(root, 1.0e-12)
        env_vel = dpos_dtheta * omega
        if env_pos >= y0:
            pos_ref.append(env_pos)
            vel_ref.append(env_vel)
        else:
            pos_ref.append(y0)
            vel_ref.append(0.0)
    return pos_ref, vel_ref


def main():
    params = json.loads((ROOT / "assets" / "onset_stress_b" / "onset_stress_b_model.json").read_text(encoding="utf-8"))
    mesh_t, mesh_p, mesh_v = load_series(OUTPUT_DIR / "onset_stress_b_mesh.csv")
    sdf1_t, sdf1_p, sdf1_v = load_series(OUTPUT_DIR / "onset_stress_b_sdf1.csv")
    sdf2_t, sdf2_p, sdf2_v = load_series(OUTPUT_DIR / "onset_stress_b_sdf2.csv")
    ref_p, ref_v = piecewise_reference(mesh_t, params)

    mpl.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Times", "Nimbus Roman No9 L", "DejaVu Serif"],
            "mathtext.fontset": "stix",
            "font.size": 10,
            "axes.labelsize": 10,
            "axes.titlesize": 10,
            "legend.fontsize": 9,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )

    fig, axes = plt.subplots(
        3,
        1,
        figsize=(7.2, 7.1),
        sharex=True,
        gridspec_kw={"height_ratios": [1.0, 1.0, 0.8]},
    )

    colors = {
        "mesh": "#1f77b4",
        "sdf1": "#d62728",
        "sdf2": "#2ca02c",
        "ref": "#111111",
    }

    ax = axes[0]
    ax.plot(mesh_t, ref_p, linestyle="--", color=colors["ref"], linewidth=1.8, label="Piecewise analytic reference")
    ax.plot(mesh_t, mesh_p, color=colors["mesh"], linewidth=1.5, label="Native mesh")
    ax.plot(sdf1_t, sdf1_p, color=colors["sdf1"], linewidth=1.5, linestyle="--", dashes=(6, 3), label="SDF 1st-order")
    ax.plot(sdf2_t, sdf2_p, color=colors["sdf2"], linewidth=1.6, label="SDF 2nd-order")
    ax.set_ylabel("Position Y (m)")
    ax.set_title("Onset Stress Benchmark (Version B)")
    ax.grid(True, alpha=0.25)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, 1.02), ncol=2, frameon=True, framealpha=0.92)

    ax = axes[1]
    ax.plot(mesh_t, ref_v, linestyle="--", color=colors["ref"], linewidth=1.8)
    ax.plot(mesh_t, mesh_v, color=colors["mesh"], linewidth=1.5)
    ax.plot(sdf1_t, sdf1_v, color=colors["sdf1"], linewidth=1.5, linestyle="--", dashes=(6, 3))
    ax.plot(sdf2_t, sdf2_v, color=colors["sdf2"], linewidth=1.6)
    ax.set_ylabel("Velocity Y (m/s)")
    ax.grid(True, alpha=0.25)

    ax = axes[2]
    mesh_err = [mesh_v[i] - ref_v[i] for i in range(len(mesh_t))]
    sdf1_err = [sdf1_v[i] - ref_v[i] for i in range(len(sdf1_t))]
    sdf2_err = [sdf2_v[i] - ref_v[i] for i in range(len(sdf2_t))]
    ax.plot(mesh_t, mesh_err, color=colors["mesh"], linewidth=1.5, label="Mesh - reference")
    ax.plot(sdf1_t, sdf1_err, color=colors["sdf1"], linewidth=1.5, linestyle="--", dashes=(6, 3), label="1st - reference")
    ax.plot(sdf2_t, sdf2_err, color=colors["sdf2"], linewidth=1.6, label="2nd - reference")
    ax.axhline(0.0, color="#666666", linewidth=0.9, alpha=0.7)
    ax.set_ylabel("Vel. error (m/s)")
    ax.set_xlabel("Time (s)")
    ax.grid(True, alpha=0.25)
    ax.legend(loc="lower right", ncol=1, frameon=True, framealpha=0.92)

    fig.subplots_adjust(left=0.11, right=0.98, top=0.95, bottom=0.09, hspace=0.14)

    FIG_DIR.mkdir(parents=True, exist_ok=True)
    png_path = FIG_DIR / "onset_stress_benchmark.png"
    pdf_path = FIG_DIR / "onset_stress_benchmark.pdf"
    fig.savefig(png_path, dpi=220)
    fig.savefig(pdf_path)
    plt.close(fig)
    print(png_path)
    print(pdf_path)


if __name__ == "__main__":
    main()
