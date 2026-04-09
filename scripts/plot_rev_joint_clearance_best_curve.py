from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
ASSET_DIR = ROOT / "assets" / "rev_joint_clearance"
OUTPUT_DIR = ROOT / "docs" / "figures"

REF_BODY3 = ASSET_DIR / "data" / "body3.csv"
OURS = {
    "SDF 1st-order": ROOT / "data" / "outputs" / "rev_joint_clearance_sdf1.csv",
    "SDF 2nd-order": ROOT / "data" / "outputs" / "rev_joint_clearance_sdf2.csv",
}

PNG_OUT = OUTPUT_DIR / "rev_joint_clearance_best_curve.png"
PDF_OUT = OUTPUT_DIR / "rev_joint_clearance_best_curve.pdf"


def load_recurdyn_csv(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    time_col = next(c for c in df.columns if c.startswith("X:Pos_TX"))
    out = {"Time": df[time_col].astype(float)}
    for col in df.columns:
        if col.startswith("Y:"):
            out[col[2:]] = df[col].astype(float)
    return pd.DataFrame(out)


def main() -> None:
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
            "axes.labelsize": 11,
            "axes.titlesize": 12,
            "legend.fontsize": 9,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
        }
    )

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    body3_ref = load_recurdyn_csv(REF_BODY3)
    ours = {name: pd.read_csv(path) for name, path in OURS.items()}
    time = next(iter(ours.values()))["Time"].to_numpy(float)
    ref_t = body3_ref["Time"].to_numpy(float)

    fig, ax = plt.subplots(figsize=(8.6, 4.8), constrained_layout=True)

    ax.plot(
        ref_t,
        body3_ref["Pos_TZ-Body3-rev_clearance_joint(m)"].to_numpy(float),
        color="black",
        linewidth=2.0,
        label="RecurDyn clearance",
    )
    ax.plot(
        time,
        ours["SDF 2nd-order"]["Body3_Pos_TZ"].to_numpy(float),
        color="#2ca02c",
        linewidth=1.8,
        label="SDF 2nd-order",
    )
    ax.plot(
        time,
        ours["SDF 1st-order"]["Body3_Pos_TZ"].to_numpy(float),
        color="#d62728",
        linewidth=1.1,
        alpha=0.7,
        label="SDF 1st-order",
    )
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Body3 TZ Position (m)")
    ax.set_title("Best Current Comparison: Body3 TZ Position")
    ax.grid(True, linewidth=0.3, alpha=0.5)
    ax.legend(loc="best", frameon=True)

    fig.savefig(PNG_OUT, dpi=220)
    fig.savefig(PDF_OUT)
    plt.close(fig)


if __name__ == "__main__":
    main()
