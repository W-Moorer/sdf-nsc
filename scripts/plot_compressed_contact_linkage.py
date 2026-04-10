from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Tuple

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
plt.rcParams["axes.labelsize"] = 8.5
plt.rcParams["axes.titlesize"] = 8.5
plt.rcParams["xtick.labelsize"] = 7.5
plt.rcParams["ytick.labelsize"] = 7.5
plt.rcParams["legend.fontsize"] = 7.5
plt.rcParams["mathtext.fontset"] = "custom"
plt.rcParams["mathtext.rm"] = "TeX Gyre Termes"
plt.rcParams["mathtext.it"] = "TeX Gyre Termes:italic"
plt.rcParams["mathtext.bf"] = "TeX Gyre Termes:bold"

FULL_WIDTH_IN = 7.05

STATIC_COLOR = "#1f4e79"
DYNAMIC_COLOR = "#b03a2e"
TEMPORAL_COLOR = "#2f6f44"
TEMPORAL_DYNAMIC_COLOR = "#6c3483"
CONTACT_SHADE = "#d9d9d9"
SCENARIO_STYLE = {
    "tilted_plate_impact": {"color": "#1f77b4", "marker": "o", "label": "Impact"},
    "tilted_plate_slide": {"color": "#d97706", "marker": "^", "label": "Slide"},
}

PANEL_LABELS = "abcdefghijklmnopqrstuvwxyz"


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
            if "variant" in row and row["variant"] != "full":
                continue
            scenario = row["scenario"]
            for key, value in row.items():
                if key in {"scenario", "variant"}:
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


def panel_label(index: int) -> str:
    return f"({PANEL_LABELS[index]})"


def decorate_panel(ax: plt.Axes, index: int) -> None:
    ax.text(
        0.02,
        0.98,
        panel_label(index),
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=8.0,
        fontweight="bold",
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.85, "pad": 1.8},
    )


def add_pair_tag(ax: plt.Axes, text: str) -> None:
    ax.text(
        0.5,
        0.98,
        text,
        transform=ax.transAxes,
        ha="center",
        va="top",
        fontsize=7.5,
        bbox={"facecolor": "white", "edgecolor": "#cccccc", "alpha": 0.9, "pad": 2.0},
    )


def add_contact_band(ax: plt.Axes, t: np.ndarray, contact_mask: np.ndarray) -> None:
    if not np.any(contact_mask):
        return
    ymin, ymax = ax.get_ylim()
    ax.fill_between(
        t,
        np.full_like(t, ymin),
        np.full_like(t, ymax),
        where=contact_mask,
        color=CONTACT_SHADE,
        alpha=0.18,
        step="mid",
        zorder=0,
    )
    ax.set_ylim(ymin, ymax)


def binned_median_curve(x: np.ndarray, y: np.ndarray, bins: int = 12) -> Tuple[np.ndarray, np.ndarray]:
    if x.size < bins:
        return np.asarray([]), np.asarray([])
    edges = np.quantile(x, np.linspace(0.0, 1.0, bins + 1))
    edges = np.unique(edges)
    if edges.size < 3:
        return np.asarray([]), np.asarray([])
    mids: List[float] = []
    meds: List[float] = []
    for left, right in zip(edges[:-1], edges[1:]):
        mask = (x >= left) & (x <= right if right == edges[-1] else x < right)
        if np.count_nonzero(mask) < 3:
            continue
        mids.append(float(np.median(x[mask])))
        meds.append(float(np.median(y[mask])))
    return np.asarray(mids), np.asarray(meds)


def plot_time_series(data: Dict[str, Dict[str, np.ndarray]]) -> None:
    fig, axes = plt.subplots(
        len(STATIC_DYNAMIC_PAIRS),
        len(SCENARIO_ORDER),
        figsize=(FULL_WIDTH_IN, 8.2),
        sharex="col",
    )
    panel_idx = 0

    for col, (scenario, label) in enumerate(SCENARIO_ORDER):
        series = data[scenario]
        t = series["time"]
        contact_mask = (series["dense_contacts"] > 0.0) | (series["reduced_contacts"] > 0.0)

        for row, (static_key, dynamic_key, static_label, dynamic_label) in enumerate(STATIC_DYNAMIC_PAIRS):
            ax = axes[row, col]
            ax2 = ax.twinx()
            ax.plot(t, series[static_key], color=STATIC_COLOR, linewidth=1.2)
            ax2.plot(t, series[dynamic_key], color=DYNAMIC_COLOR, linewidth=1.1, linestyle="--")
            add_contact_band(ax, t, contact_mask)
            decorate_panel(ax, panel_idx)
            add_pair_tag(ax, f"{static_label} / {dynamic_label}")
            panel_idx += 1

            if row == 0:
                ax.set_title(label, fontsize=9.0, fontweight="bold", pad=6)
            if col == 0:
                ax.set_ylabel(static_label, color=STATIC_COLOR)
            if col == len(SCENARIO_ORDER) - 1:
                ax2.set_ylabel(dynamic_label, color=DYNAMIC_COLOR)
            if row == len(STATIC_DYNAMIC_PAIRS) - 1:
                ax.set_xlabel("Time (s)")

            ax.grid(True, linestyle="--", linewidth=0.55, alpha=0.3)
            apply_time_series_y_scale(ax, series[static_key])
            apply_time_series_y_scale(ax2, series[dynamic_key])
            ax.tick_params(axis="y", colors=STATIC_COLOR, length=3)
            ax2.tick_params(axis="y", colors=DYNAMIC_COLOR, length=3)
            ax.tick_params(axis="x", length=3)

    line_static = plt.Line2D([], [], color=STATIC_COLOR, linewidth=1.2, label="Static error")
    line_dynamic = plt.Line2D([], [], color=DYNAMIC_COLOR, linewidth=1.1, linestyle="--", label="Dynamic response error")
    band = plt.Rectangle((0, 0), 1, 1, fc=CONTACT_SHADE, alpha=0.18, label="Active-contact interval")
    fig.legend(
        handles=[line_static, line_dynamic, band],
        loc="upper center",
        ncol=3,
        frameon=False,
        bbox_to_anchor=(0.5, 0.995),
        columnspacing=1.5,
        handlelength=2.0,
    )
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.965))
    save_figure(fig, "compressed_contact_static_dynamic_timeseries")
    plt.close(fig)


def plot_correlation_scatter(data: Dict[str, Dict[str, np.ndarray]]) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(FULL_WIDTH_IN, 5.75))
    panel_idx = 0

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
            style = SCENARIO_STYLE[scenario]
            ax.scatter(
                x,
                y,
                s=11,
                alpha=0.38,
                color=style["color"],
                marker=style["marker"],
                label=style["label"] if static_key == "epsF" else None,
                edgecolors="none",
            )
            corr_rows.append((scenario, static_key, dynamic_key, pearson_corr(x, y)))

        x_all = np.concatenate(all_x) if all_x else np.asarray([])
        y_all = np.concatenate(all_y) if all_y else np.asarray([])
        corr_all = pearson_corr(x_all, y_all)
        corr_rows.append(("all", static_key, dynamic_key, corr_all))
        trend_x, trend_y = binned_median_curve(x_all, y_all)

        apply_log_like_scale(ax, x_all, "x")
        apply_log_like_scale(ax, y_all, "y")
        if trend_x.size > 0:
            ax.plot(trend_x, trend_y, color="#222222", linewidth=1.2, label="Median trend" if static_key == "epsF" else None)
        ax.set_xlabel(static_label)
        ax.set_ylabel(dynamic_label)
        ax.set_title(f"{static_label} vs {dynamic_label}", pad=4)
        ax.grid(True, linestyle="--", linewidth=0.55, alpha=0.28)
        decorate_panel(ax, panel_idx)
        panel_idx += 1
        ax.text(
            0.97,
            0.96,
            rf"$\rho={corr_all:.3f}$",
            transform=ax.transAxes,
            ha="right",
            va="top",
            bbox={"facecolor": "white", "edgecolor": "#cccccc", "alpha": 0.9, "pad": 3},
        )

    handles, labels = axes[0, 0].get_legend_handles_labels()
    if handles:
        fig.legend(handles, labels, loc="upper center", ncol=3, frameon=False, bbox_to_anchor=(0.5, 0.995))

    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.955))
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
    fig, axes = plt.subplots(
        len(TEMPORAL_PAIRS),
        len(SCENARIO_ORDER),
        figsize=(FULL_WIDTH_IN, 6.2),
        sharex="col",
    )
    panel_idx = 0

    for col, (scenario, label) in enumerate(SCENARIO_ORDER):
        series = data[scenario]
        t = series["time"]
        contact_mask = (series["dense_contacts"] > 0.0) | (series["reduced_contacts"] > 0.0)

        for row, (temporal_key, dynamic_key, temporal_label, dynamic_label) in enumerate(TEMPORAL_PAIRS):
            ax = axes[row, col]
            ax2 = ax.twinx()
            ax.plot(t, series[temporal_key], color=TEMPORAL_COLOR, linewidth=1.2)
            ax2.plot(t, series[dynamic_key], color=TEMPORAL_DYNAMIC_COLOR, linewidth=1.1, linestyle="--")
            add_contact_band(ax, t, contact_mask)
            decorate_panel(ax, panel_idx)
            add_pair_tag(ax, f"{temporal_label} / {dynamic_label}")
            panel_idx += 1

            if row == 0:
                ax.set_title(label, fontsize=9.0, fontweight="bold", pad=6)
            if col == 0:
                ax.set_ylabel(temporal_label, color=TEMPORAL_COLOR)
            if col == len(SCENARIO_ORDER) - 1:
                ax2.set_ylabel(dynamic_label, color=TEMPORAL_DYNAMIC_COLOR)
            if row == len(TEMPORAL_PAIRS) - 1:
                ax.set_xlabel("Time (s)")

            ax.grid(True, linestyle="--", linewidth=0.55, alpha=0.3)
            apply_time_series_y_scale(ax, series[temporal_key])
            apply_time_series_y_scale(ax2, series[dynamic_key])
            ax.tick_params(axis="y", colors=TEMPORAL_COLOR, length=3)
            ax2.tick_params(axis="y", colors=TEMPORAL_DYNAMIC_COLOR, length=3)
            ax.tick_params(axis="x", length=3)

    line_temporal = plt.Line2D([], [], color=TEMPORAL_COLOR, linewidth=1.2, label="Temporal-coherence metric")
    line_dynamic = plt.Line2D([], [], color=TEMPORAL_DYNAMIC_COLOR, linewidth=1.1, linestyle="--", label="Dynamic response error")
    band = plt.Rectangle((0, 0), 1, 1, fc=CONTACT_SHADE, alpha=0.18, label="Active-contact interval")
    fig.legend(
        handles=[line_temporal, line_dynamic, band],
        loc="upper center",
        ncol=3,
        frameon=False,
        bbox_to_anchor=(0.5, 0.995),
        columnspacing=1.4,
        handlelength=2.0,
    )
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.955))
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
