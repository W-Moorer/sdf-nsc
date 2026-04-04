from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent.parent
FIG_DIR = ROOT / "papers" / "figures"
DATA_DIR = ROOT / "data" / "outputs"
REF_PATH = ROOT / "assets" / "simple_gear" / "data" / "Gear22.csv"
FIRST_PATH = DATA_DIR / "gear_tricubic_1st.csv"
SECOND_PATH = DATA_DIR / "gear_tricubic_2nd.csv"


plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["DejaVu Sans", "Arial"]
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42
plt.rcParams["font.size"] = 11
plt.rcParams["mathtext.fontset"] = "dejavusans"


def load_csv(path: Path, wrx_col_idx: int) -> tuple[np.ndarray, np.ndarray]:
    if not path.exists():
        return np.array([]), np.array([])

    lines = path.read_text(encoding="utf-8", errors="ignore").splitlines()
    start = 0
    for i, line in enumerate(lines):
        if "Time" in line or "Y:Vel_RX" in line or "Y:Vel" in line:
            start = i + 1
            break

    times = []
    wrxs = []
    for line in lines[start:]:
        parts = line.strip().split(",")
        if len(parts) <= wrx_col_idx:
            continue
        try:
            times.append(float(parts[0]))
            wrxs.append(float(parts[wrx_col_idx]))
        except ValueError:
            continue
    return np.asarray(times), np.asarray(wrxs)


def calc_mae(y: np.ndarray) -> float:
    if len(y) == 0:
        return 0.0
    return float(np.mean(np.abs(y - (-1.0))))


t_ref, wrx_ref = load_csv(REF_PATH, 6)
t_1, wrx_1 = load_csv(FIRST_PATH, 3)
t_2, wrx_2 = load_csv(SECOND_PATH, 3)

mae_ref = calc_mae(wrx_ref)
mae_1 = calc_mae(wrx_1)
mae_2 = calc_mae(wrx_2)

fig, ax = plt.subplots(figsize=(8, 6))

ax.axhline(y=-1.0, color="gray", linestyle="-.", linewidth=1.0, alpha=0.65)

if len(t_ref) > 0:
    ax.plot(
        t_ref,
        wrx_ref,
        color="black",
        linewidth=1.6,
        label=f"Commercial Industrial Software (MAE: {mae_ref:.3f})",
    )
if len(t_1) > 0:
    ax.plot(
        t_1,
        wrx_1,
        color="#2f3fff",
        linewidth=1.4,
        alpha=0.98,
        label=f"1st-order SDF (MAE: {mae_1:.3f})",
    )
if len(t_2) > 0:
    ax.plot(
        t_2,
        wrx_2,
        color="#ff4a45",
        linewidth=1.3,
        alpha=0.98,
        label=f"2nd-order SDF (MAE: {mae_2:.3f})",
    )

ax.text(0.347, -0.93, "Theoretical Solution (-1)", color="gray", fontsize=13, alpha=0.85)
ax.set_xlabel("Time (s)", fontsize=14)
ax.set_ylabel(r"Driven Gear Velocity $\omega_{rx}$ (rad/s)", fontsize=14)
ax.set_xlim([0.0, 0.5])
ax.set_ylim([-3.5, 1.8])
ax.grid(True, linestyle="--", linewidth=0.8, alpha=0.45)
ax.legend(loc="upper right", framealpha=1.0, edgecolor="black", fontsize=12)

plt.tight_layout()
FIG_DIR.mkdir(parents=True, exist_ok=True)
pdf_path = FIG_DIR / "gear_error.pdf"
png_path = FIG_DIR / "gear_error_render.png"
plt.savefig(pdf_path, format="pdf", dpi=300)
plt.savefig(png_path, format="png", dpi=150)
print(f"Saved {pdf_path}")
print(f"Saved {png_path}")
