import csv
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parents[1]
REF_PATH = ROOT / "assets" / "cam" / "data" / "cam_data.csv"
FIG_DIR = ROOT / "papers" / "paper1" / "figures" / "generated"


def pick_existing(candidates):
    for candidate in candidates:
        if candidate.exists():
            return candidate
    raise FileNotFoundError("Missing required benchmark trace: " + ", ".join(str(c) for c in candidates))


def load_trace(path: Path):
    with path.open("r", encoding="utf-8", errors="replace", newline="") as f:
        rows = list(csv.reader(f))
    data = [[float(x) for x in row] for row in rows[1:] if row]
    return {
        "time": [row[0] for row in data],
        "pos_y": [row[10] for row in data],
        "vel_y": [row[13] for row in data],
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
        a = (x - x0) / (x1 - x0)
        out.append(y0 * (1.0 - a) + y1 * a)
    return out


def rmse(a, b):
    return (sum((x - y) ** 2 for x, y in zip(a, b)) / len(a)) ** 0.5


def configure_plotting():
    mpl.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Times", "Nimbus Roman No9 L", "DejaVu Serif"],
            "mathtext.fontset": "stix",
            "font.size": 10,
            "axes.labelsize": 10,
            "axes.titlesize": 10,
            "legend.fontsize": 8.8,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def main():
    configure_plotting()

    mesh_path = pick_existing(
        [
            ROOT / "data" / "outputs" / "paper_cam_mesh_dt0p01.csv",
            ROOT / "data" / "outputs" / "cam_reverify_mesh__dt0p01.csv",
        ]
    )
    sdf1_path = pick_existing(
        [
            ROOT / "data" / "outputs" / "paper_cam_sdf1_dt0p01.csv",
            ROOT / "data" / "outputs" / "cam_current_default_sdf1_dt0p01.csv",
            ROOT / "data" / "outputs" / "cam_finalcheck_sdf1_dt0p01.csv",
        ]
    )

    ref = load_trace(REF_PATH)
    mesh = load_trace(mesh_path)
    sdf1 = load_trace(sdf1_path)

    ref_pos_mesh = interp_linear(ref["time"], ref["pos_y"], mesh["time"])
    ref_vel_mesh = interp_linear(ref["time"], ref["vel_y"], mesh["time"])
    ref_pos_sdf1 = interp_linear(ref["time"], ref["pos_y"], sdf1["time"])
    ref_vel_sdf1 = interp_linear(ref["time"], ref["vel_y"], sdf1["time"])

    mesh_pos_rmse = rmse(mesh["pos_y"], ref_pos_mesh)
    mesh_vel_rmse = rmse(mesh["vel_y"], ref_vel_mesh)
    sdf1_pos_rmse = rmse(sdf1["pos_y"], ref_pos_sdf1)
    sdf1_vel_rmse = rmse(sdf1["vel_y"], ref_vel_sdf1)

    colors = {"ref": "#111111", "mesh": "#1f77b4", "sdf1": "#d62728"}
    fig, (ax_pos, ax_vel) = plt.subplots(2, 1, figsize=(6.8, 5.8), sharex=True)

    ax_pos.plot(ref["time"], ref["pos_y"], linestyle="--", color=colors["ref"], linewidth=1.8, label="Reference")
    ax_pos.plot(
        mesh["time"],
        mesh["pos_y"],
        color=colors["mesh"],
        linewidth=1.4,
        label=f"Native mesh (RMSE = {mesh_pos_rmse:.4f} m)",
    )
    ax_pos.plot(
        sdf1["time"],
        sdf1["pos_y"],
        color=colors["sdf1"],
        linewidth=1.5,
        label=f"1st-order SDF (RMSE = {sdf1_pos_rmse:.4f} m)",
    )
    ax_pos.set_ylabel("Follower position Y (m)")
    ax_pos.grid(True, alpha=0.25)
    ax_pos.legend(loc="upper right", frameon=True, framealpha=0.92)

    ax_vel.plot(ref["time"], ref["vel_y"], linestyle="--", color=colors["ref"], linewidth=1.8, label="Reference")
    ax_vel.plot(
        mesh["time"],
        mesh["vel_y"],
        color=colors["mesh"],
        linewidth=1.4,
        label=f"Native mesh (RMSE = {mesh_vel_rmse:.4f} m/s)",
    )
    ax_vel.plot(
        sdf1["time"],
        sdf1["vel_y"],
        color=colors["sdf1"],
        linewidth=1.5,
        label=f"1st-order SDF (RMSE = {sdf1_vel_rmse:.4f} m/s)",
    )
    ax_vel.set_xlabel("Time (s)")
    ax_vel.set_ylabel("Follower velocity Y (m/s)")
    ax_vel.grid(True, alpha=0.25)
    ax_vel.legend(loc="lower right", frameon=True, framealpha=0.92)

    fig.subplots_adjust(left=0.12, right=0.98, top=0.98, bottom=0.10, hspace=0.18)
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG_DIR / "cam_benchmark_dt0p01.pdf")
    fig.savefig(FIG_DIR / "cam_benchmark_dt0p01.png", dpi=220)
    plt.close(fig)


if __name__ == "__main__":
    main()
