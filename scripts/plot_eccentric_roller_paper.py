import csv
import json
import math
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset


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


def analytic_reference(times, params):
    radius = params["cam_radius"]
    ecc = params["cam_eccentricity"]
    roller = params["roller_radius"]
    omega = params["motor_speed"]
    reach = radius + roller

    pos_ref, vel_ref = [], []
    for t in times:
        theta = omega * t
        cx = ecc * math.cos(theta)
        cy = ecc * math.sin(theta)
        root = math.sqrt(max(reach * reach - cx * cx, 0.0))
        pos = cy + root
        dpos_dtheta = ecc * math.cos(theta) + (ecc * ecc * math.cos(theta) * math.sin(theta)) / max(root, 1.0e-12)
        vel = dpos_dtheta * omega
        pos_ref.append(pos)
        vel_ref.append(vel)
    return pos_ref, vel_ref


def main():
    params = json.loads((ROOT / "assets" / "eccentric_roller" / "eccentric_roller_model.json").read_text(encoding="utf-8"))
    mesh_t, mesh_p, mesh_v = load_series(OUTPUT_DIR / "eccentric_roller_mesh.csv")
    sdf1_t, sdf1_p, sdf1_v = load_series(OUTPUT_DIR / "eccentric_roller_sdf1.csv")
    sdf2_t, sdf2_p, sdf2_v = load_series(OUTPUT_DIR / "eccentric_roller_sdf2.csv")
    ref_p, ref_v = analytic_reference(mesh_t, params)

    mpl.rcParams.update({
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
    })

    fig, axes = plt.subplots(
        3,
        1,
        figsize=(7.2, 7.0),
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
    ax.plot(mesh_t, ref_p, linestyle="--", color=colors["ref"], linewidth=1.8, label="Analytic reference")
    ax.plot(mesh_t, mesh_p, color=colors["mesh"], linewidth=1.5, label="Native mesh")
    ax.plot(sdf1_t, sdf1_p, color=colors["sdf1"], linewidth=1.5, linestyle="--", dashes=(6, 3),
            label="SDF 1st-order")
    ax.plot(sdf2_t, sdf2_p, color=colors["sdf2"], linewidth=1.6, label="SDF 2nd-order")
    ax.set_ylabel("Position Y (m)")
    ax.set_title("Eccentric Disk / Roller Follower Benchmark")
    ax.grid(True, alpha=0.25)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, 1.02), ncol=2, frameon=True, framealpha=0.92)

    ax = axes[1]
    ax.plot(mesh_t, ref_v, linestyle="--", color=colors["ref"], linewidth=1.8)
    ax.plot(mesh_t, mesh_v, color=colors["mesh"], linewidth=1.5)
    ax.plot(sdf1_t, sdf1_v, color=colors["sdf1"], linewidth=1.5, linestyle="--", dashes=(6, 3))
    ax.plot(sdf2_t, sdf2_v, color=colors["sdf2"], linewidth=1.6)
    ax.set_ylabel("Velocity Y (m/s)")
    ax.grid(True, alpha=0.25)

    zoom_x0, zoom_x1 = 0.67, 0.79
    zoom_indices = [i for i, t in enumerate(mesh_t) if zoom_x0 <= t <= zoom_x1]
    zoom_vals = (
        [ref_v[i] for i in zoom_indices]
        + [mesh_v[i] for i in zoom_indices]
        + [sdf1_v[i] for i in zoom_indices]
        + [sdf2_v[i] for i in zoom_indices]
    )
    zoom_y0 = min(zoom_vals) - 0.002
    zoom_y1 = max(zoom_vals) + 0.002

    axins = inset_axes(ax, width="37%", height="46%", loc="lower right", borderpad=2.0)
    axins.plot(mesh_t, ref_v, linestyle="--", color=colors["ref"], linewidth=1.2)
    axins.plot(mesh_t, mesh_v, color=colors["mesh"], linewidth=1.1)
    axins.plot(sdf1_t, sdf1_v, color=colors["sdf1"], linewidth=1.2, linestyle="--", dashes=(6, 3))
    axins.plot(sdf2_t, sdf2_v, color=colors["sdf2"], linewidth=1.2)
    axins.set_xlim(zoom_x0, zoom_x1)
    axins.set_ylim(zoom_y0, zoom_y1)
    axins.set_xticks([0.68, 0.73, 0.78])
    axins.set_yticks([])
    axins.grid(True, alpha=0.2)
    mark_inset(ax, axins, loc1=2, loc2=4, fc="none", ec="#666666", lw=0.8)

    ax = axes[2]
    mesh_err = [mesh_v[i] - ref_v[i] for i in range(len(mesh_t))]
    sdf1_err = [sdf1_v[i] - ref_v[i] for i in range(len(sdf1_t))]
    sdf2_err = [sdf2_v[i] - ref_v[i] for i in range(len(sdf2_t))]
    ax.plot(mesh_t, mesh_err, color=colors["mesh"], linewidth=1.5, label="Mesh - analytic")
    ax.plot(sdf1_t, sdf1_err, color=colors["sdf1"], linewidth=1.5, linestyle="--", dashes=(6, 3),
            label="1st - analytic")
    ax.plot(sdf2_t, sdf2_err, color=colors["sdf2"], linewidth=1.6, label="2nd - analytic")
    ax.axhline(0.0, color="#666666", linewidth=0.9, alpha=0.7)
    ax.set_ylabel("Vel. error (m/s)")
    ax.set_xlabel("Time (s)")
    ax.grid(True, alpha=0.25)
    ax.legend(loc="upper right", ncol=3, frameon=True, framealpha=0.92)

    fig.subplots_adjust(left=0.11, right=0.98, top=0.95, bottom=0.09, hspace=0.14)

    FIG_DIR.mkdir(parents=True, exist_ok=True)
    png_path = FIG_DIR / "eccentric_roller_benchmark.png"
    pdf_path = FIG_DIR / "eccentric_roller_benchmark.pdf"
    fig.savefig(png_path, dpi=220)
    fig.savefig(pdf_path)
    plt.close(fig)
    print(png_path)
    print(pdf_path)


if __name__ == "__main__":
    main()
