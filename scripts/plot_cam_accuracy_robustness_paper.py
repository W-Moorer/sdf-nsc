import csv
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
OUTPUT_DIR = ROOT / "data" / "outputs"
FIG_DIR = ROOT / "papers" / "figures"
REF_PATH = ROOT / "assets" / "cam" / "data" / "cam_data.csv"


def load_reference():
    df = pd.read_csv(REF_PATH)
    return {
        "time": df.iloc[:, 0].tolist(),
        "vel": df.iloc[:, 13].tolist(),
    }


def load_trace(path: Path):
    times, vel = [], []
    with path.open("r", encoding="utf-8") as f:
        reader = csv.reader(f)
        next(reader)
        for row in reader:
            if not row:
                continue
            times.append(float(row[0]))
            vel.append(float(row[13]))
    return times, vel


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


def rmse(trace_v, ref_v):
    n = len(trace_v)
    return (sum((a - b) ** 2 for a, b in zip(trace_v, ref_v)) / n) ** 0.5


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


def make_figure(
    out_pdf: Path,
    out_png: Path,
    title: str,
    mesh_path: Path,
    sdf1_path: Path,
    sdf2_path: Path,
):
    ref = load_reference()
    t_mesh, v_mesh = load_trace(mesh_path)
    t_sdf1, v_sdf1 = load_trace(sdf1_path)
    t_sdf2, v_sdf2 = load_trace(sdf2_path)

    ref_mesh = interp_linear(ref["time"], ref["vel"], t_mesh)
    ref_sdf1 = interp_linear(ref["time"], ref["vel"], t_sdf1)
    ref_sdf2 = interp_linear(ref["time"], ref["vel"], t_sdf2)

    mesh_rmse = rmse(v_mesh, ref_mesh)
    sdf1_rmse = rmse(v_sdf1, ref_sdf1)
    sdf2_rmse = rmse(v_sdf2, ref_sdf2)

    fig, ax = plt.subplots(figsize=(6.8, 3.8))
    colors = {
        "ref": "#111111",
        "mesh": "#1f77b4",
        "sdf1": "#d62728",
        "sdf2": "#2ca02c",
    }

    ax.plot(ref["time"], ref["vel"], linestyle="--", color=colors["ref"], linewidth=1.8, label="Reference")
    ax.plot(t_mesh, v_mesh, color=colors["mesh"], linewidth=1.5, label=f"Native mesh (RMSE = {mesh_rmse:.3f})")
    ax.plot(
        t_sdf1,
        v_sdf1,
        color=colors["sdf1"],
        linewidth=1.5,
        linestyle="--",
        dashes=(6, 3),
        label=f"1st-order SDF (RMSE = {sdf1_rmse:.3f})",
    )
    ax.plot(t_sdf2, v_sdf2, color=colors["sdf2"], linewidth=1.6, label=f"2nd-order SDF (RMSE = {sdf2_rmse:.3f})")
    ax.set_title(title)
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Follower velocity Y (m/s)")
    ax.grid(True, alpha=0.25)
    ax.legend(loc="lower right", frameon=True, framealpha=0.92)
    fig.subplots_adjust(left=0.12, right=0.98, top=0.90, bottom=0.16)
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_pdf)
    fig.savefig(out_png, dpi=220)
    plt.close(fig)


def main():
    configure_plotting()
    make_figure(
        FIG_DIR / "fig_accuracy.pdf",
        FIG_DIR / "fig_accuracy.png",
        r"Industrial cam benchmark at $\Delta t = 0.005$ s",
        OUTPUT_DIR / "cam_reverify_mesh__dt0p005.csv",
        OUTPUT_DIR / "cam_current_default_sdf1_dt0p005.csv",
        OUTPUT_DIR / "cam_current_default_sdf2_dt0p005.csv",
    )
    make_figure(
        FIG_DIR / "fig_robustness.pdf",
        FIG_DIR / "fig_robustness.png",
        r"Industrial cam benchmark at $\Delta t = 0.01$ s",
        OUTPUT_DIR / "cam_reverify_mesh__dt0p01.csv",
        OUTPUT_DIR / "cam_current_default_sdf1_dt0p01.csv",
        OUTPUT_DIR / "cam_current_default_sdf2_dt0p01.csv",
    )


if __name__ == "__main__":
    main()
