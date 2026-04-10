from __future__ import annotations

import csv
from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import font_manager as fm
from matplotlib.colors import TwoSlopeNorm


ROOT = Path(__file__).resolve().parent.parent
STATIC_CSV = ROOT / "data" / "generated" / "compressed_contact_validation.csv"
DYNAMIC_CSV = ROOT / "data" / "generated" / "compressed_contact_dynamics.csv"
FIG_DIR = ROOT / "papers" / "paper1" / "figures" / "generated"
SECTION_DIR = ROOT / "papers" / "paper1" / "sections" / "generated"

TERMES_FONT_FILES = [
    Path("/usr/share/texmf/fonts/opentype/public/tex-gyre/texgyretermes-regular.otf"),
    Path("/usr/share/texmf/fonts/opentype/public/tex-gyre/texgyretermes-bold.otf"),
    Path("/usr/share/texmf/fonts/opentype/public/tex-gyre/texgyretermes-italic.otf"),
    Path("/usr/share/texmf/fonts/opentype/public/tex-gyre/texgyretermes-bolditalic.otf"),
]

VARIANT_ORDER = [
    "full",
    "fixed4",
    "single_patch",
    "no_dense_micro",
    "no_eC",
    "no_sentinel",
    "no_impulse_transport",
    "no_reinj_accept",
]

VARIANT_LABELS = {
    "full": "Full",
    "fixed4": "Fixed-4",
    "single_patch": "Single-patch",
    "no_dense_micro": "No dense micro",
    "no_eC": "No $e_{\\mathcal{C}}$",
    "no_sentinel": "No sentinel",
    "no_impulse_transport": "No transport",
    "no_reinj_accept": "No reinj. accept",
}

SCENARIO_ORDER = [
    ("tilted_plate_impact", "Tilted Plate Impact"),
    ("tilted_plate_slide", "Tilted Plate Slide"),
]

DYNAMIC_METRICS = [
    ("ang_vel_error", "Max\nang. vel."),
    ("linear_impulse_error", "Max\nlin. imp."),
    ("angular_impulse_error", "Max\nang. imp."),
    ("energy_drift_diff", "Max\nenergy drift"),
]

FULL_WIDTH_IN = 7.05


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
plt.rcParams["font.size"] = 9
plt.rcParams["axes.labelsize"] = 8.5
plt.rcParams["axes.titlesize"] = 9
plt.rcParams["xtick.labelsize"] = 7.5
plt.rcParams["ytick.labelsize"] = 7.5
plt.rcParams["mathtext.fontset"] = "custom"
plt.rcParams["mathtext.rm"] = "TeX Gyre Termes"
plt.rcParams["mathtext.it"] = "TeX Gyre Termes:italic"
plt.rcParams["mathtext.bf"] = "TeX Gyre Termes:bold"


def load_static_rows() -> List[Dict[str, str]]:
    with STATIC_CSV.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def load_dynamic_rows() -> List[Dict[str, str]]:
    with DYNAMIC_CSV.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def format_decimal(value: float, digits: int = 4) -> str:
    return f"{value:.{digits}f}"


def format_small(value: float) -> str:
    abs_value = abs(value)
    if abs_value == 0.0:
        return "0"
    if abs_value < 1.0e-3 or abs_value >= 1.0e3:
        mantissa, exponent = f"{value:.2e}".split("e")
        return rf"${mantissa}\times 10^{{{int(exponent)}}}$"
    if abs_value < 1.0e-2:
        return f"{value:.5f}"
    if abs_value < 1.0e-1:
        return f"{value:.4f}"
    return f"{value:.3f}"


def summarize_static(rows: List[Dict[str, str]]) -> Dict[str, Dict[str, float | str]]:
    summary: Dict[str, Dict[str, float | str]] = {}
    for variant in VARIANT_ORDER:
        vrows = [row for row in rows if row["variant"] == variant]
        if not vrows:
            continue
        pass_count = sum(1 for row in vrows if row["pass"] == "1")
        check = f"{pass_count}/{len(vrows)}"
        if variant == "no_sentinel":
            check = "n/a"
        summary[variant] = {
            "check": check,
            "avg_ratio": float(np.mean([float(row["compression_ratio"]) for row in vrows])),
            "max_epsM": max(float(row["epsilon_M"]) for row in vrows),
            "max_epsCoP": max(float(row["epsilon_CoP"]) for row in vrows),
            "max_epsGap": max(float(row["epsilon_gap"]) for row in vrows),
        }
    return summary


def summarize_dynamic(rows: List[Dict[str, str]]) -> Dict[Tuple[str, str], Dict[str, float]]:
    summary: Dict[Tuple[str, str], Dict[str, float]] = {}
    for variant in VARIANT_ORDER:
        for scenario, _ in SCENARIO_ORDER:
            srows = [row for row in rows if row["variant"] == variant and row["scenario"] == scenario]
            if not srows:
                continue
            summary[(variant, scenario)] = {
                metric: max(float(row[metric]) for row in srows)
                for metric, _ in DYNAMIC_METRICS
            }
    return summary


def save_figure(fig: plt.Figure, stem: str) -> None:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    pdf_path = FIG_DIR / f"{stem}.pdf"
    png_path = FIG_DIR / f"{stem}.png"
    fig.savefig(pdf_path, format="pdf", dpi=300, bbox_inches="tight")
    fig.savefig(png_path, format="png", dpi=180, bbox_inches="tight")
    print(f"saved {pdf_path}")
    print(f"saved {png_path}")


def plot_dynamic_variant_heatmap(summary: Dict[Tuple[str, str], Dict[str, float]]) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(FULL_WIDTH_IN, 3.7))
    norm = TwoSlopeNorm(vmin=0.8, vcenter=1.0, vmax=1.8)
    cmap = plt.get_cmap("RdBu_r")

    last_im = None
    for idx, (scenario, title) in enumerate(SCENARIO_ORDER):
        ax = axes[idx]
        matrix = []
        for variant in VARIANT_ORDER:
            baseline = summary[("full", scenario)]
            metrics = summary[(variant, scenario)]
            matrix.append([metrics[key] / baseline[key] for key, _ in DYNAMIC_METRICS])
        data = np.asarray(matrix, dtype=float)
        last_im = ax.imshow(data, cmap=cmap, norm=norm, aspect="auto")

        ax.set_title(title, pad=6, fontweight="bold")
        ax.set_xticks(np.arange(len(DYNAMIC_METRICS)))
        ax.set_xticklabels([label for _, label in DYNAMIC_METRICS])
        ax.set_yticks(np.arange(len(VARIANT_ORDER)))
        ax.set_yticklabels([VARIANT_LABELS[variant] for variant in VARIANT_ORDER])
        if idx == 1:
            ax.set_yticklabels([])

        for row in range(len(VARIANT_ORDER)):
            for col in range(len(DYNAMIC_METRICS)):
                value = data[row, col]
                text_color = "white" if value > 1.24 else "#222222"
                ax.text(
                    col,
                    row,
                    f"{value:.2f}x",
                    ha="center",
                    va="center",
                    color=text_color,
                    fontsize=7.2,
                )

        ax.axhline(2.5, color="#444444", linewidth=0.8)
        ax.set_xticks(np.arange(-0.5, len(DYNAMIC_METRICS), 1), minor=True)
        ax.set_yticks(np.arange(-0.5, len(VARIANT_ORDER), 1), minor=True)
        ax.grid(which="minor", color="white", linewidth=0.7)
        ax.tick_params(which="minor", bottom=False, left=False)

    fig.subplots_adjust(left=0.16, right=0.84, bottom=0.14, top=0.90, wspace=0.08)
    cax = fig.add_axes([0.86, 0.18, 0.02, 0.64])
    cbar = fig.colorbar(last_im, cax=cax)
    cbar.set_label("Ratio to full (lower is better)")
    save_figure(fig, "compressed_contact_variant_dynamics")
    plt.close(fig)


def write_static_table(summary: Dict[str, Dict[str, float | str]]) -> None:
    SECTION_DIR.mkdir(parents=True, exist_ok=True)
    path = SECTION_DIR / "compressed_contact_variant_static_table.tex"
    lines = [
        r"\begin{tabular}{@{}lcccc@{}}",
        r"\toprule",
        r"Variant & Static check & Avg.\ $N_{\mathrm{red}}/N_{\mathrm{dense}}$ & Max $\varepsilon_M$ & Max $\varepsilon_{\mathrm{CoP}}$ / Max $\varepsilon_g$ \\",
        r"\midrule",
    ]
    for variant in VARIANT_ORDER:
        row = summary[variant]
        lines.append(
            " & ".join(
                [
                    VARIANT_LABELS[variant],
                    str(row["check"]),
                    format_decimal(float(row["avg_ratio"]), 4),
                    format_small(float(row["max_epsM"])),
                    f"{format_small(float(row['max_epsCoP']))} / {format_small(float(row['max_epsGap']))}",
                ]
            )
            + r" \\"
        )
    lines.extend(
        [
            r"\bottomrule",
            r"\end{tabular}",
        ]
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"saved {path}")


def write_dynamic_table(summary: Dict[Tuple[str, str], Dict[str, float]]) -> None:
    SECTION_DIR.mkdir(parents=True, exist_ok=True)
    path = SECTION_DIR / "compressed_contact_variant_dynamic_table.tex"
    lines = [
        r"\begin{tabular}{@{}lcccc@{}}",
        r"\toprule",
        r"Variant & Impact max ang.-vel. & Impact max lin.-imp. & Slide max ang.-vel. & Slide max lin.-imp. \\",
        r"\midrule",
    ]
    for variant in VARIANT_ORDER:
        impact = summary[(variant, "tilted_plate_impact")]
        slide = summary[(variant, "tilted_plate_slide")]
        lines.append(
            " & ".join(
                [
                    VARIANT_LABELS[variant],
                    format_decimal(impact["ang_vel_error"], 3),
                    format_decimal(impact["linear_impulse_error"], 4),
                    format_decimal(slide["ang_vel_error"], 3),
                    format_decimal(slide["linear_impulse_error"], 4),
                ]
            )
            + r" \\"
        )
    lines.extend(
        [
            r"\bottomrule",
            r"\end{tabular}",
        ]
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"saved {path}")


def main() -> None:
    if not STATIC_CSV.exists():
        raise FileNotFoundError(f"missing static csv: {STATIC_CSV}")
    if not DYNAMIC_CSV.exists():
        raise FileNotFoundError(f"missing dynamic csv: {DYNAMIC_CSV}")

    static_summary = summarize_static(load_static_rows())
    dynamic_summary = summarize_dynamic(load_dynamic_rows())
    plot_dynamic_variant_heatmap(dynamic_summary)
    write_static_table(static_summary)
    write_dynamic_table(dynamic_summary)


if __name__ == "__main__":
    main()
