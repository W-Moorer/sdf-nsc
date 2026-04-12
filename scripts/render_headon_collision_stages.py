#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch


ROOT = Path(__file__).resolve().parents[1]
FIG_DIR = ROOT / "papers" / "figures"

COLORS = {
    "ball_a_fill": "#dbe4ee",
    "ball_a_edge": "#475569",
    "ball_b_fill": "#fde6c8",
    "ball_b_edge": "#b45309",
    "text": "#111827",
    "accent": "#0f172a",
}

PANEL_WIDTH = 6.9
PANEL_HEIGHT = 2.8
BALL_RADIUS = 0.72
BALL_Y = 0.96
ARROW_Y = 1.92


def configure_matplotlib() -> None:
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Times", "Nimbus Roman", "DejaVu Serif"],
            "font.size": 11,
            "axes.titlesize": 13,
            "axes.titleweight": "semibold",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
        }
    )


def add_ball(ax, center, label: str, fill: str, edge: str) -> None:
    circle = Circle(center, radius=BALL_RADIUS, facecolor=fill, edgecolor=edge, linewidth=2.0)
    ax.add_patch(circle)
    ax.text(center[0], center[1], label, ha="center", va="center", fontsize=13, color=COLORS["text"], fontweight="semibold")


def add_velocity(ax, start, dx: float, text: str, color: str) -> None:
    if abs(dx) > 1e-9:
        arrow = FancyArrowPatch(
            posA=(start[0], start[1]),
            posB=(start[0] + dx, start[1]),
            arrowstyle="-|>",
            mutation_scale=15,
            linewidth=1.8,
            color=color,
        )
        ax.add_patch(arrow)
        mid_x = start[0] + 0.5 * dx
        label_y = start[1] + 0.14
    else:
        mid_x = start[0]
        label_y = start[1] + 0.12
    ax.text(mid_x, label_y, text, ha="center", va="bottom", fontsize=9.5, color=color, fontstyle="italic")


def add_contact_marker(ax, x: float, y: float) -> None:
    ax.scatter([x], [y], s=42, marker="*", color=COLORS["accent"], zorder=6)


def setup_panel(ax, title: str) -> None:
    ax.set_xlim(0.0, PANEL_WIDTH)
    ax.set_ylim(0.0, PANEL_HEIGHT)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_frame_on(False)
    ax.text(
        0.5,
        1.035,
        title,
        transform=ax.transAxes,
        ha="center",
        va="top",
        fontsize=10.8,
        color=COLORS["text"],
        clip_on=False,
    )


def build_equal_mass_figure() -> plt.Figure:
    fig, axes = plt.subplots(1, 3, figsize=(8.55, 2.65))
    fig.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.02, wspace=0.03)
    titles = ["Initial", "Just Contact", "Final"]
    for ax, title in zip(axes, titles):
        setup_panel(ax, title)

    # Initial
    add_ball(axes[0], (1.58, BALL_Y), "A", COLORS["ball_a_fill"], COLORS["ball_a_edge"])
    add_ball(axes[0], (5.32, BALL_Y), "B", COLORS["ball_b_fill"], COLORS["ball_b_edge"])
    add_velocity(axes[0], (0.70, ARROW_Y), 1.38, "vA = +1.0 m/s", COLORS["ball_a_edge"])
    add_velocity(axes[0], (5.32, ARROW_Y), 0.0, "vB = 0", COLORS["ball_b_edge"])

    # Impact
    add_ball(axes[1], (2.73, BALL_Y), "A", COLORS["ball_a_fill"], COLORS["ball_a_edge"])
    add_ball(axes[1], (4.17, BALL_Y), "B", COLORS["ball_b_fill"], COLORS["ball_b_edge"])
    add_velocity(axes[1], (1.62, ARROW_Y), 1.16, "vA = +1.0 m/s", COLORS["ball_a_edge"])
    add_velocity(axes[1], (4.17, ARROW_Y), 0.0, "vB = 0", COLORS["ball_b_edge"])
    add_contact_marker(axes[1], 3.45, BALL_Y)

    # Final
    add_ball(axes[2], (1.72, BALL_Y), "A", COLORS["ball_a_fill"], COLORS["ball_a_edge"])
    add_ball(axes[2], (4.98, BALL_Y), "B", COLORS["ball_b_fill"], COLORS["ball_b_edge"])
    add_velocity(axes[2], (1.72, ARROW_Y), 0.0, "vA = 0", COLORS["ball_a_edge"])
    add_velocity(axes[2], (3.80, ARROW_Y), 1.42, "vB = +1.0 m/s", COLORS["ball_b_edge"])
    return fig


def build_unequal_mass_figure() -> plt.Figure:
    fig, axes = plt.subplots(1, 3, figsize=(8.55, 2.65))
    fig.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.02, wspace=0.03)
    titles = ["Initial", "Just Contact", "Final"]
    for ax, title in zip(axes, titles):
        setup_panel(ax, title)

    # Initial
    add_ball(axes[0], (1.58, BALL_Y), "A", COLORS["ball_a_fill"], COLORS["ball_a_edge"])
    add_ball(axes[0], (5.32, BALL_Y), "B", COLORS["ball_b_fill"], COLORS["ball_b_edge"])
    add_velocity(axes[0], (1.58, ARROW_Y), 0.0, "vA = 0", COLORS["ball_a_edge"])
    add_velocity(axes[0], (6.20, ARROW_Y), -1.38, "vB = -1.0 m/s", COLORS["ball_b_edge"])

    # Impact
    add_ball(axes[1], (2.73, BALL_Y), "A", COLORS["ball_a_fill"], COLORS["ball_a_edge"])
    add_ball(axes[1], (4.17, BALL_Y), "B", COLORS["ball_b_fill"], COLORS["ball_b_edge"])
    add_velocity(axes[1], (2.73, ARROW_Y), 0.0, "vA = 0", COLORS["ball_a_edge"])
    add_velocity(axes[1], (5.12, ARROW_Y), -1.12, "vB = -1.0 m/s", COLORS["ball_b_edge"])
    add_contact_marker(axes[1], 3.45, BALL_Y)

    # Final
    add_ball(axes[2], (1.88, BALL_Y), "A", COLORS["ball_a_fill"], COLORS["ball_a_edge"])
    add_ball(axes[2], (4.82, BALL_Y), "B", COLORS["ball_b_fill"], COLORS["ball_b_edge"])
    add_velocity(axes[2], (3.12, ARROW_Y), -1.80, "vA = -4/3 m/s", COLORS["ball_a_edge"])
    add_velocity(axes[2], (5.76, ARROW_Y), -0.84, "vB = -1/3 m/s", COLORS["ball_b_edge"])
    return fig


def save_figure(fig: plt.Figure, stem: str) -> None:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    pdf_path = FIG_DIR / f"{stem}.pdf"
    svg_path = FIG_DIR / f"{stem}.svg"
    fig.savefig(pdf_path, bbox_inches="tight", pad_inches=0.01)
    fig.savefig(svg_path, bbox_inches="tight", pad_inches=0.01)
    print(f"Saved {pdf_path}")
    print(f"Saved {svg_path}")
    plt.close(fig)


def main() -> None:
    configure_matplotlib()
    save_figure(build_equal_mass_figure(), "headon_equal_spheres_states")
    save_figure(build_unequal_mass_figure(), "headon_unequal_mass_spheres_states")


if __name__ == "__main__":
    main()
