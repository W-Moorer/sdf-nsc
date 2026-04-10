import csv
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parents[1]
OUTPUT_DIR = ROOT / "data" / "outputs"
FIG_DIR = ROOT / "papers" / "paper1" / "figures" / "generated"


def load_series(path: Path):
    times, vel_a, vel_b = [], [], []
    with path.open("r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            times.append(float(row["Time"]))
            vel_a.append(float(row["SphereA_Vel_X"]))
            vel_b.append(float(row["SphereB_Vel_X"]))
    return times, vel_a, vel_b


def analytic_reference(times, impact_time=0.2):
    ref_a = [1.0 if t < impact_time else 0.0 for t in times]
    ref_b = [0.0 if t < impact_time else 1.0 for t in times]
    return ref_a, ref_b


def rmse(values, reference):
    n = max(1, len(values))
    return (sum((a - b) ** 2 for a, b in zip(values, reference)) / n) ** 0.5


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

    mesh_t, mesh_a, mesh_b = load_series(OUTPUT_DIR / "headon_spheres_mesh.csv")
    sdf1_t, sdf1_a, sdf1_b = load_series(OUTPUT_DIR / "headon_spheres_sdf1.csv")

    ref_a, ref_b = analytic_reference(mesh_t)

    mesh_rmse_a = rmse(mesh_a, ref_a)
    sdf1_rmse_a = rmse(sdf1_a, ref_a)
    mesh_rmse_b = rmse(mesh_b, ref_b)
    sdf1_rmse_b = rmse(sdf1_b, ref_b)

    colors = {
        "ref": "#111111",
        "mesh": "#1f77b4",
        "sdf1": "#d62728",
    }

    fig, axes = plt.subplots(2, 1, figsize=(7.2, 5.2), sharex=True)

    ax = axes[0]
    ax.plot(mesh_t, ref_a, linestyle="--", color=colors["ref"], linewidth=1.8, label="Analytic reference")
    ax.plot(mesh_t, mesh_a, color=colors["mesh"], linewidth=1.5,
            label=f"Native mesh (RMSE = {mesh_rmse_a:.3e})")
    ax.plot(sdf1_t, sdf1_a, color=colors["sdf1"], linewidth=1.5, linestyle="--", dashes=(6, 3),
            label=f"SDF 1st-order (RMSE = {sdf1_rmse_a:.3e})")
    ax.set_ylabel("Sphere A velocity X (m/s)")
    ax.set_title("Head-on Equal-Sphere Impact Benchmark")
    ax.grid(True, alpha=0.25)
    ax.legend(loc="upper right", frameon=True, framealpha=0.92)

    ax = axes[1]
    ax.plot(mesh_t, ref_b, linestyle="--", color=colors["ref"], linewidth=1.8, label="Analytic reference")
    ax.plot(mesh_t, mesh_b, color=colors["mesh"], linewidth=1.5,
            label=f"Native mesh (RMSE = {mesh_rmse_b:.3e})")
    ax.plot(sdf1_t, sdf1_b, color=colors["sdf1"], linewidth=1.5, linestyle="--", dashes=(6, 3),
            label=f"SDF 1st-order (RMSE = {sdf1_rmse_b:.3e})")
    ax.set_ylabel("Sphere B velocity X (m/s)")
    ax.set_xlabel("Time (s)")
    ax.grid(True, alpha=0.25)
    ax.legend(loc="lower right", frameon=True, framealpha=0.92)

    fig.subplots_adjust(left=0.12, right=0.98, top=0.92, bottom=0.12, hspace=0.16)
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    pdf_path = FIG_DIR / "headon_spheres_benchmark.pdf"
    png_path = FIG_DIR / "headon_spheres_benchmark.png"
    fig.savefig(pdf_path)
    fig.savefig(png_path, dpi=220)
    plt.close(fig)
    print(pdf_path)
    print(png_path)


if __name__ == "__main__":
    main()
