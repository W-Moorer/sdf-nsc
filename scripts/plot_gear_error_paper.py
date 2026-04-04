from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import font_manager


ROOT = Path(__file__).resolve().parent.parent
FIG_DIR = ROOT / "papers" / "figures"
DATA_DIR = ROOT / "data" / "outputs"
REF_PATH = ROOT / "assets" / "simple_gear" / "data" / "Gear22.csv"
FIRST_PATH = DATA_DIR / "gear_default_pinned_s1.csv"
SECOND_PATH = DATA_DIR / "gear_default_pinned_s2.csv"


def configure_times_font() -> None:
    font_candidates = [
        Path("C:/Windows/Fonts/times.ttf"),
        Path("C:/Windows/Fonts/timesbd.ttf"),
        Path("C:/Windows/Fonts/timesi.ttf"),
        Path("C:/Windows/Fonts/timesbi.ttf"),
        Path("/mnt/c/Windows/Fonts/times.ttf"),
        Path("/mnt/c/Windows/Fonts/timesbd.ttf"),
        Path("/mnt/c/Windows/Fonts/timesi.ttf"),
        Path("/mnt/c/Windows/Fonts/timesbi.ttf"),
    ]

    for font_path in font_candidates:
        if font_path.exists():
            try:
                font_manager.fontManager.addfont(str(font_path))
            except Exception:
                pass

    plt.rcParams["font.family"] = "serif"
    plt.rcParams["font.serif"] = ["Times New Roman", "Times", "DejaVu Serif"]
    plt.rcParams["pdf.fonttype"] = 42
    plt.rcParams["ps.fonttype"] = 42
    plt.rcParams["font.size"] = 11
    plt.rcParams["mathtext.fontset"] = "stix"
    plt.rcParams["axes.unicode_minus"] = False


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


def calc_mae(values: np.ndarray) -> float:
    if len(values) == 0:
        return 0.0
    return float(np.mean(np.abs(values + 1.0)))


def main() -> None:
    configure_times_font()

    t_ref, wrx_ref = load_csv(REF_PATH, 6)
    t_1, wrx_1 = load_csv(FIRST_PATH, 3)
    t_2, wrx_2 = load_csv(SECOND_PATH, 3)

    mae_ref = calc_mae(wrx_ref)
    mae_1 = calc_mae(wrx_1)
    mae_2 = calc_mae(wrx_2)

    all_values = np.concatenate(
        [arr for arr in (wrx_ref, wrx_1, wrx_2) if len(arr) > 0]
    )
    y_min = float(np.min(all_values))
    y_max = float(np.max(all_values))
    y_span = max(1e-6, y_max - y_min)

    fig, ax = plt.subplots(figsize=(6.8, 4.1))

    ax.axhline(
        y=-1.0,
        color="#666666",
        linestyle=(0, (5, 2.5)),
        linewidth=1.0,
        alpha=0.85,
        label="Theoretical solution",
        zorder=1,
    )

    if len(t_ref) > 0:
        ax.plot(
            t_ref,
            wrx_ref,
            color="black",
            linewidth=1.8,
            label=f"Commercial software (MAE = {mae_ref:.3f})",
            zorder=4,
        )
    if len(t_1) > 0:
        ax.plot(
            t_1,
            wrx_1,
            color="#1f4fd6",
            linewidth=1.7,
            label=f"1st-order SDF (MAE = {mae_1:.3f})",
            zorder=3,
        )
    if len(t_2) > 0:
        ax.plot(
            t_2,
            wrx_2,
            color="#d43f2f",
            linewidth=1.7,
            label=f"2nd-order SDF (MAE = {mae_2:.3f})",
            zorder=2,
        )

    ax.set_xlabel("Time (s)", fontsize=13)
    ax.set_ylabel(r"Driven gear velocity $\omega_{rx}$ (rad/s)", fontsize=13)
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(y_min - 0.06 * y_span, y_max + 0.08 * y_span)
    ax.margins(x=0.0)
    ax.grid(True, linestyle="--", linewidth=0.6, alpha=0.28)

    for spine in ax.spines.values():
        spine.set_linewidth(0.9)

    ax.tick_params(axis="both", labelsize=11, width=0.9)
    ax.legend(
        loc="upper right",
        frameon=True,
        framealpha=0.96,
        facecolor="white",
        edgecolor="black",
        fontsize=10.5,
        borderpad=0.45,
        labelspacing=0.35,
        handlelength=2.2,
        handletextpad=0.6,
    )

    fig.subplots_adjust(left=0.12, right=0.99, bottom=0.15, top=0.98)
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    pdf_path = FIG_DIR / "gear_error_paper_default.pdf"
    png_path = FIG_DIR / "gear_error_paper_default.png"
    fig.savefig(pdf_path, dpi=300, bbox_inches="tight", pad_inches=0.02)
    fig.savefig(png_path, dpi=300, bbox_inches="tight", pad_inches=0.02)

    print(f"Saved {pdf_path}")
    print(f"Saved {png_path}")
    print(f"Commercial MAE: {mae_ref:.6f}")
    print(f"1st-order MAE: {mae_1:.6f}")
    print(f"2nd-order MAE: {mae_2:.6f}")


if __name__ == "__main__":
    main()
