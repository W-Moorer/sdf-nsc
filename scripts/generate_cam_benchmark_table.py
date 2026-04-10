import csv
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
REF_PATH = ROOT / "assets" / "cam" / "data" / "cam_data.csv"
OUT_PATH = ROOT / "papers" / "paper1" / "sections" / "generated" / "cam_benchmark_table.tex"


def pick_existing(candidates):
    for candidate in candidates:
        if candidate.exists():
            return candidate
    raise FileNotFoundError("Missing required benchmark trace: " + ", ".join(str(c) for c in candidates))


def load_trace(path: Path):
    with path.open("r", encoding="utf-8", errors="replace", newline="") as f:
        rows = list(csv.reader(f))
    data = [[float(x) for x in row] for row in rows[1:] if row]
    return {
        "time": [row[0] for row in data],
        "pos_y": [row[10] for row in data],
        "vel_y": [row[13] for row in data],
    }


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
        a = (x - x0) / (x1 - x0)
        out.append(y0 * (1.0 - a) + y1 * a)
    return out


def rmse(a, b):
    return math.sqrt(sum((x - y) ** 2 for x, y in zip(a, b)) / len(a))


def metrics(trace, ref):
    ref_pos = interp_linear(ref["time"], ref["pos_y"], trace["time"])
    ref_vel = interp_linear(ref["time"], ref["vel_y"], trace["time"])
    return {
        "pos_rmse": rmse(trace["pos_y"], ref_pos),
        "vel_rmse": rmse(trace["vel_y"], ref_vel),
        "pos_max": max(abs(a - b) for a, b in zip(trace["pos_y"], ref_pos)),
        "vel_max": max(abs(a - b) for a, b in zip(trace["vel_y"], ref_vel)),
        "final_pos": trace["pos_y"][-1],
        "final_vel": trace["vel_y"][-1],
        "ref_final_pos": ref_pos[-1],
        "ref_final_vel": ref_vel[-1],
    }


def main():
    ref = load_trace(REF_PATH)
    mesh = load_trace(
        pick_existing(
            [
                ROOT / "data" / "outputs" / "paper_cam_mesh_dt0p01.csv",
                ROOT / "data" / "outputs" / "cam_reverify_mesh__dt0p01.csv",
            ]
        )
    )
    sdf1 = load_trace(
        pick_existing(
            [
                ROOT / "data" / "outputs" / "paper_cam_sdf1_dt0p01.csv",
                ROOT / "data" / "outputs" / "cam_current_default_sdf1_dt0p01.csv",
                ROOT / "data" / "outputs" / "cam_finalcheck_sdf1_dt0p01.csv",
            ]
        )
    )

    rows = [("Native mesh", metrics(mesh, ref)), ("1st-order SDF", metrics(sdf1, ref))]
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    with OUT_PATH.open("w", encoding="utf-8", newline="\n") as f:
        f.write("\\begin{tabular}{lcccccc}\n")
        f.write("\\toprule\n")
        f.write("Model & Pos RMSE & Vel RMSE & Peak $|\\Delta v|$ & Final $y$ & Final $v_y$ & Reference $(y, v_y)$ \\\\\n")
        f.write("\\midrule\n")
        for label, m in rows:
            ref_pair = f"({m['ref_final_pos']:.3f},\\ {m['ref_final_vel']:.3f})"
            f.write(
                f"{label} & {m['pos_rmse']:.3e} & {m['vel_rmse']:.3e} & {m['vel_max']:.3e} & "
                f"{m['final_pos']:.3f} & {m['final_vel']:.3f} & ${ref_pair}$ \\\\\n"
            )
        f.write("\\bottomrule\n")
        f.write("\\end{tabular}\n")


if __name__ == "__main__":
    main()
