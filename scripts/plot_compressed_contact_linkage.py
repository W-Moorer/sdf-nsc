from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path
from typing import Dict, List

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import font_manager as fm


ROOT = Path(__file__).resolve().parent.parent
DYNAMIC_CSV = ROOT / "data" / "generated" / "compressed_contact_dynamics.csv"
FIG_DIR = ROOT / "papers" / "paper1" / "figures" / "generated"
CORR_CSV = ROOT / "data" / "generated" / "compressed_contact_linkage_correlations.csv"

TERMES_FONT_FILES = [
    Path("/usr/share/texmf/fonts/opentype/public/tex-gyre/texgyretermes-regular.otf"),
    Path("/usr/share/texmf/fonts/opentype/public/tex-gyre/texgyretermes-bold.otf"),
    Path("/usr/share/texmf/fonts/opentype/public/tex-gyre/texgyretermes-italic.otf"),
    Path("/usr/share/texmf/fonts/opentype/public/tex-gyre/texgyretermes-bolditalic.otf"),
]

SCENARIO_ORDER = [
    ("tilted_plate_impact", "Tilted Plate Impact"),
    ("tilted_plate_slide", "Tilted Plate Slide"),
]

STATIC_DYNAMIC_PAIRS = [
    ("epsF", "linear_impulse_error", r"$\varepsilon_F$", "Linear Impulse Error"),
    ("epsM", "angular_impulse_error", r"$\varepsilon_M$", "Angular Impulse Error"),
    ("epsCoP", "ang_vel_error", r"$\varepsilon_{CoP}$", "Angular Velocity Error"),
    ("epsGap", "energy_drift_diff", r"$\varepsilon_{gap}$", "Energy Drift Difference"),
]

TEMPORAL_PAIRS = [
    ("temporal_hausdorff", "ang_vel_error", "Temporal Hausdorff", "Angular Velocity Error"),
    ("temporal_mean_drift", "angular_impulse_error", "Mean Support Drift", "Angular Impulse Error"),
    ("support_churn", "linear_impulse_error", "Support Churn", "Linear Impulse Error"),
]


for font_path in TERMES_FONT_FILES:
    if font_path.exists():
        fm.fontManager.addfont(str(font_path))

plt.rcParams["font.family"] = "TeX Gyre Termes"
plt.rcParams["font.serif"] = [
    "TeX Gyre Termes",
    "Tinos",
    "Nimbus Roman",
    "Nimbus Roman No9 L",
]
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42
plt.rcParams["font.size"] = 10
plt.rcParams["mathtext.fontset"] = "custom"
plt.rcParams["mathtext.rm"] = "TeX Gyre Termes"
plt.rcParams["mathtext.it"] = "TeX Gyre Termes:italic"
plt.rcParams["mathtext.bf"] = "TeX Gyre Termes:bold"


def pearson_corr(x: np.ndarray, y: np.ndarray) -> float:
    if x.size < 2 or y.size < 2:
        return 0.0
    if np.allclose(x, x[0]) or np.allclose(y, y[0]):
        return 0.0
    return float(np.corrcoef(x, y)[0, 1])


def load_dynamic_rows(path: Path) -> Dict[str, Dict[str, np.ndarray]]:
    rows_by_scenario: Dict[str, Dict[str, List[float]]] = defaultdict(lambda: defaultdict(list))
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            scenario = row["scenario"]
            for key, value in row.items():
                if key == "scenario":
                    continue
                rows_by_scenario[scenario][key].append(float(value))

    out: Dict[str, Dict[str, np.ndarray]] = {}
    for scenario, values in rows_by_scenario.items():
        out[scenario] = {key: np.asarray(series, dtype=float) for key, series in values.items()}
    return out


def apply_log_like_scale(ax: plt.Axes, values: np.ndarray, axis: str) -> None:
    positive = np.abs(values[np.abs(values) > 0.0])
    if positive.size == 0:
        return
    linthresh = max(1.0e-8, float(np.min(positive)))
    if axis == "x":
        ax.set_xscale("symlog", linthresh=linthresh)
    else:
        ax.set_yscale("symlog", linthresh=linthresh)


def apply_time_series_y_scale(ax: plt.Axes, values: np.ndarray) -> None:
    positive = np.abs(values[np.abs(values) > 0.0])
    if positive.size == 0:
        return
    vmax = float(np.max(positive))
    vmin = float(np.min(positive))
    if vmax / max(vmin, 1.0e-12) > 100.0:
        ax.set_yscale("symlog", linthresh=max(1.0e-8, vmin))


def save_figure(fig: plt.Figure, stem: str) -> None:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    pdf_path = FIG_DIR / f"{stem}.pdf"
    png_path = FIG_DIR / f"{stem}.png"
    fig.savefig(pdf_path, format="pdf", dpi=300, bbox_inches="tight")
    fig.savefig(png_path, format="png", dpi=180, bbox_inches="tight")
    print(f"saved {pdf_path}")
    print(f"saved {png_path}")


def plot_time_series(data: Dict[str, Dict[str, np.ndarray]]) -> None:
    fig, axes = plt.subplots(len(STATIC_DYNAMIC_PAIRS), len(SCENARIO_ORDER), figsize=(12.5, 10.5), sharex="col")
    static_color = "#2f55d4"
    dynamic_color = "#d94841"

    for col, (scenario, label) in enumerate(SCENARIO_ORDER):
        series = data[scenario]
        t = series["time"]
        contact_mask = (series["dense_contacts"] > 0.0) | (series["reduced_contacts"] > 0.0)

        for row, (static_key, dynamic_key, static_label, dynamic_label) in enumerate(STATIC_DYNAMIC_PAIRS):
            ax = axes[row, col]
            ax2 = ax.twinx()
            ax.plot(t, series[static_key], color=static_color, linewidth=1.5, label=static_label)
            ax2.plot(t, series[dynamic_key], color=dynamic_color, linewidth=1.3, label=dynamic_label)
            if np.any(contact_mask):
                ymin, ymax = ax.get_ylim()
                ax.fill_between(
                    t,
                    np.full_like(t, ymin),
                    np.full_like(t, ymax),
                    where=contact_mask,
                    color="#d9d9d9",
                    alpha=0.18,
                    step="mid",
                )
                ax.set_ylim(ymin, ymax)

            if row == 0:
                ax.set_title(label, fontsize=12, weight="bold")
            if col == 0:
                ax.set_ylabel(static_label, color=static_color)
            if col == len(SCENARIO_ORDER) - 1:
                ax2.set_ylabel(dynamic_label, color=dynamic_color)
            if row == len(STATIC_DYNAMIC_PAIRS) - 1:
                ax.set_xlabel("Time (s)")

            ax.grid(True, linestyle="--", linewidth=0.7, alpha=0.4)
            apply_time_series_y_scale(ax, series[static_key])
            apply_time_series_y_scale(ax2, series[dynamic_key])
            ax.tick_params(axis="y", colors=static_color)
            ax2.tick_params(axis="y", colors=dynamic_color)

    fig.suptitle("Instantaneous Static Errors and Dynamic Response", fontsize=14, weight="bold", y=0.995)
    fig.tight_layout()
    save_figure(fig, "compressed_contact_static_dynamic_timeseries")
    plt.close(fig)


def plot_correlation_scatter(data: Dict[str, Dict[str, np.ndarray]]) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(11.5, 9.0))
    colors = {
        "tilted_plate_impact": "#1f77b4",
        "tilted_plate_slide": "#ff7f0e",
    }

    corr_rows = []
    for ax, (static_key, dynamic_key, static_label, dynamic_label) in zip(axes.flat, STATIC_DYNAMIC_PAIRS):
        all_x = []
        all_y = []
        for scenario, scenario_label in SCENARIO_ORDER:
            series = data[scenario]
            x = series[static_key]
            y = series[dynamic_key]
            mask = (series["dense_contacts"] > 0.0) | (series["reduced_contacts"] > 0.0)
            x = x[mask]
            y = y[mask]
            all_x.append(x)
            all_y.append(y)
            ax.scatter(
                x,
                y,
                s=18,
                alpha=0.7,
                color=colors[scenario],
                label=scenario_label if static_key == "epsF" else None,
                edgecolors="none",
            )
            corr_rows.append((scenario, static_key, dynamic_key, pearson_corr(x, y)))

        x_all = np.concatenate(all_x) if all_x else np.asarray([])
        y_all = np.concatenate(all_y) if all_y else np.asarray([])
        corr_all = pearson_corr(x_all, y_all)
        corr_rows.append(("all", static_key, dynamic_key, corr_all))

        apply_log_like_scale(ax, x_all, "x")
        apply_log_like_scale(ax, y_all, "y")
        ax.set_xlabel(static_label)
        ax.set_ylabel(dynamic_label)
        ax.grid(True, linestyle="--", linewidth=0.7, alpha=0.35)
        ax.text(
            0.03,
            0.96,
            rf"$\rho={corr_all:.3f}$",
            transform=ax.transAxes,
            ha="left",
            va="top",
            bbox={"facecolor": "white", "edgecolor": "#cccccc", "alpha": 0.9, "pad": 3},
        )

    handles, labels = axes[0, 0].get_legend_handles_labels()
    if handles:
        fig.legend(handles, labels, loc="upper center", ncol=2, frameon=False, bbox_to_anchor=(0.5, 1.02))

    fig.suptitle("Static-to-Dynamic Error Correlation", fontsize=14, weight="bold", y=0.99)
    fig.tight_layout()
    save_figure(fig, "compressed_contact_static_dynamic_correlation")
    plt.close(fig)

    CORR_CSV.parent.mkdir(parents=True, exist_ok=True)
    with CORR_CSV.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["scenario", "x_metric", "y_metric", "pearson_r"])
        writer.writerows(corr_rows)
        for temporal_key, dynamic_key, _, _ in TEMPORAL_PAIRS:
            for scenario, _ in SCENARIO_ORDER:
                series = data[scenario]
                mask = (series["dense_contacts"] > 0.0) | (series["reduced_contacts"] > 0.0)
                x = series[temporal_key][mask]
                y = series[dynamic_key][mask]
                writer.writerow([scenario, temporal_key, dynamic_key, pearson_corr(x, y)])
            x_all = np.concatenate(
                [data[scenario][temporal_key][(data[scenario]["dense_contacts"] > 0.0) |
                                              (data[scenario]["reduced_contacts"] > 0.0)]
                 for scenario, _ in SCENARIO_ORDER]
            )
            y_all = np.concatenate(
                [data[scenario][dynamic_key][(data[scenario]["dense_contacts"] > 0.0) |
                                             (data[scenario]["reduced_contacts"] > 0.0)]
                 for scenario, _ in SCENARIO_ORDER]
            )
            writer.writerow(["all", temporal_key, dynamic_key, pearson_corr(x_all, y_all)])
    print(f"saved {CORR_CSV}")


def plot_temporal_coherence(data: Dict[str, Dict[str, np.ndarray]]) -> None:
    fig, axes = plt.subplots(len(TEMPORAL_PAIRS), len(SCENARIO_ORDER), figsize=(12.5, 8.8), sharex="col")
    temporal_color = "#2c7a4b"
    dynamic_color = "#7b3294"

    for col, (scenario, label) in enumerate(SCENARIO_ORDER):
        series = data[scenario]
        t = series["time"]
        contact_mask = (series["dense_contacts"] > 0.0) | (series["reduced_contacts"] > 0.0)

        for row, (temporal_key, dynamic_key, temporal_label, dynamic_label) in enumerate(TEMPORAL_PAIRS):
            ax = axes[row, col]
            ax2 = ax.twinx()
            ax.plot(t, series[temporal_key], color=temporal_color, linewidth=1.5, label=temporal_label)
            ax2.plot(t, series[dynamic_key], color=dynamic_color, linewidth=1.3, label=dynamic_label)
            if np.any(contact_mask):
                ymin, ymax = ax.get_ylim()
                ax.fill_between(
                    t,
                    np.full_like(t, ymin),
                    np.full_like(t, ymax),
                    where=contact_mask,
                    color="#d9d9d9",
                    alpha=0.18,
                    step="mid",
                )
                ax.set_ylim(ymin, ymax)

            if row == 0:
                ax.set_title(label, fontsize=12, weight="bold")
            if col == 0:
                ax.set_ylabel(temporal_label, color=temporal_color)
            if col == len(SCENARIO_ORDER) - 1:
                ax2.set_ylabel(dynamic_label, color=dynamic_color)
            if row == len(TEMPORAL_PAIRS) - 1:
                ax.set_xlabel("Time (s)")

            ax.grid(True, linestyle="--", linewidth=0.7, alpha=0.4)
            apply_time_series_y_scale(ax, series[temporal_key])
            apply_time_series_y_scale(ax2, series[dynamic_key])
            ax.tick_params(axis="y", colors=temporal_color)
            ax2.tick_params(axis="y", colors=dynamic_color)

    fig.suptitle("Temporal Coherence and Dynamic Stability", fontsize=14, weight="bold", y=0.995)
    fig.tight_layout()
    save_figure(fig, "compressed_contact_temporal_coherence")
    plt.close(fig)


def main() -> None:
    if not DYNAMIC_CSV.exists():
        raise FileNotFoundError(f"dynamic csv not found: {DYNAMIC_CSV}")

    data = load_dynamic_rows(DYNAMIC_CSV)
    missing = [scenario for scenario, _ in SCENARIO_ORDER if scenario not in data]
    if missing:
        raise RuntimeError(f"missing scenarios in dynamic csv: {missing}")

    plot_time_series(data)
    plot_correlation_scatter(data)
    plot_temporal_coherence(data)


if __name__ == "__main__":
    main()
