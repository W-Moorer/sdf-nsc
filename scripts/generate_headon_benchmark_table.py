import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUTPUT_DIR = ROOT / "data" / "outputs"
TABLE_DIR = ROOT / "papers" / "paper1" / "sections" / "generated"


def load_series(path: Path):
    times, vel_a, vel_b = [], [], []
    with path.open("r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            times.append(float(row["Time"]))
            vel_a.append(float(row["SphereA_Vel_X"]))
            vel_b.append(float(row["SphereB_Vel_X"]))
    return times, vel_a, vel_b


def analytic_equal(times, impact_time=0.2):
    ref_a = [1.0 if t < impact_time else 0.0 for t in times]
    ref_b = [0.0 if t < impact_time else 1.0 for t in times]
    return ref_a, ref_b, 0.0, 1.0


def analytic_mass_ratio(times, impact_time=0.2):
    ref_a = [0.0 if t < impact_time else (-4.0 / 3.0) for t in times]
    ref_b = [-1.0 if t < impact_time else (-1.0 / 3.0) for t in times]
    return ref_a, ref_b, -4.0 / 3.0, -1.0 / 3.0


def rmse(values, reference):
    n = max(1, len(values))
    return (sum((a - b) ** 2 for a, b in zip(values, reference)) / n) ** 0.5


def summarize(case_name, csv_name, analytic_fn):
    times, vel_a, vel_b = load_series(OUTPUT_DIR / csv_name)
    ref_a, ref_b, final_a_ref, final_b_ref = analytic_fn(times)
    final_a = vel_a[-1]
    final_b = vel_b[-1]
    return {
        "case": case_name,
        "final_a_ref": final_a_ref,
        "final_b_ref": final_b_ref,
        "final_a": final_a,
        "final_b": final_b,
        "abs_err_a": abs(final_a - final_a_ref),
        "abs_err_b": abs(final_b - final_b_ref),
        "rmse_a": rmse(vel_a, ref_a),
        "rmse_b": rmse(vel_b, ref_b),
    }


def format_row(case_label, model_label, row):
    analytic_pair = f"$({row['final_a_ref']:.3f},\\ {row['final_b_ref']:.3f})$"
    final_pair = f"$({row['final_a']:.3f},\\ {row['final_b']:.3f})$"
    return (
        f"{case_label} & {model_label} & {analytic_pair} & {final_pair} & "
        f"{row['abs_err_a']:.2e} & {row['abs_err_b']:.2e} & "
        f"{row['rmse_a']:.2e} & {row['rmse_b']:.2e} \\\\"
    )


def main():
    equal_mesh = summarize("Equal mass", "headon_spheres_mesh.csv", analytic_equal)
    equal_sdf1 = summarize("Equal mass", "headon_spheres_sdf1.csv", analytic_equal)
    ratio_mesh = summarize("Mass ratio 1:2", "headon_spheres_massratio_mesh.csv", analytic_mass_ratio)
    ratio_sdf1 = summarize("Mass ratio 1:2", "headon_spheres_massratio_sdf1.csv", analytic_mass_ratio)

    lines = [
        "\\begin{tabular}{llcccccc}",
        "\\toprule",
        "Case & Model & Analytic $\\,(v_A^+,v_B^+)\\,$ & Final $\\,(v_A,v_B)\\,$ & $|\\Delta v_A|$ & $|\\Delta v_B|$ & RMSE$_A$ & RMSE$_B$ \\\\",
        "\\midrule",
        format_row("Equal mass", "Native mesh", equal_mesh),
        format_row("", "SDF 1st-order", equal_sdf1),
        "\\midrule",
        format_row("Mass ratio 1:2", "Native mesh", ratio_mesh),
        format_row("", "SDF 1st-order", ratio_sdf1),
        "\\bottomrule",
        "\\end{tabular}",
    ]

    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    out_path = TABLE_DIR / "headon_benchmark_table.tex"
    out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(out_path)


if __name__ == "__main__":
    main()
