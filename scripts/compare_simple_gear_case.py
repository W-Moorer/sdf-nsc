#!/usr/bin/env python3
import argparse
import csv
import json
import math
import os
import re
from pathlib import Path


def _find_column_index(header, patterns):
    if not patterns:
        return None
    compiled = [re.compile(p, re.IGNORECASE) for p in patterns]
    for idx, name in enumerate(header):
        if all(p.search(name) for p in compiled):
            return idx
    return None


def _resolve_indices(header, column_map):
    t_idx = _find_column_index(header, column_map.get("time_col_patterns"))
    wrx_idx = _find_column_index(header, column_map.get("wrx_col_patterns"))
    wry_idx = _find_column_index(header, column_map.get("wry_col_patterns"))
    wrz_idx = _find_column_index(header, column_map.get("wrz_col_patterns"))

    if t_idx is None:
        t_idx = column_map.get("time_col_idx", 0)
    if wrx_idx is None:
        wrx_idx = column_map.get("wrx_col_idx", 3)
    if wry_idx is None:
        wry_idx = column_map.get("wry_col_idx", 4)
    if wrz_idx is None:
        wrz_idx = column_map.get("wrz_col_idx", 5)

    return t_idx, wrx_idx, wry_idx, wrz_idx


def read_gear_csv(filepath, column_map, source_name):
    data = {"Time": [], "Wrx": [], "Wry": [], "Wrz": []}

    with open(filepath, "r", encoding="utf-8", errors="ignore") as f:
        reader = csv.reader(f)
        header = next(reader)
        t_idx, wrx_idx, wry_idx, wrz_idx = _resolve_indices(header, column_map)

        print(f"Resolved {source_name} columns:")
        print(f" - Time: [{t_idx}] {header[t_idx]}")
        print(f" - Wrx : [{wrx_idx}] {header[wrx_idx]}")
        print(f" - Wry : [{wry_idx}] {header[wry_idx]}")
        print(f" - Wrz : [{wrz_idx}] {header[wrz_idx]}")

        for row in reader:
            if not row:
                continue
            try:
                data["Time"].append(float(row[t_idx]))
                data["Wrx"].append(float(row[wrx_idx]))
                data["Wry"].append(float(row[wry_idx]))
                data["Wrz"].append(float(row[wrz_idx]))
            except (ValueError, IndexError):
                continue

    return data


def interpolate(t, times, values):
    if not times:
        return 0.0
    if t <= times[0]:
        return values[0]
    if t >= times[-1]:
        return values[-1]

    for i in range(len(times) - 1):
        if times[i] <= t <= times[i + 1]:
            t0, t1 = times[i], times[i + 1]
            v0, v1 = values[i], values[i + 1]
            if t1 == t0:
                return v0
            r = (t - t0) / (t1 - t0)
            return v0 + r * (v1 - v0)

    return values[-1]


def numerical_derivative(times, values):
    if not times or not values or len(times) != len(values):
        return []
    if len(times) == 1:
        return [0.0]

    deriv = [0.0] * len(values)

    # Forward difference at start.
    dt0 = times[1] - times[0]
    deriv[0] = (values[1] - values[0]) / dt0 if abs(dt0) > 1e-15 else 0.0

    # Central difference for interior points.
    for i in range(1, len(values) - 1):
        dt = times[i + 1] - times[i - 1]
        deriv[i] = (values[i + 1] - values[i - 1]) / dt if abs(dt) > 1e-15 else 0.0

    # Backward difference at end.
    dtn = times[-1] - times[-2]
    deriv[-1] = (values[-1] - values[-2]) / dtn if abs(dtn) > 1e-15 else 0.0

    return deriv


def compare_all_rows(sim_data, ref_data):
    stats = {
        "rows_ref": len(ref_data["Time"]),
        "rows_compared": 0,
        "rows_extrapolated_boundary": 0,
        "wrx_rmse": 0.0,
        "wry_rmse": 0.0,
        "wrz_rmse": 0.0,
        "wrx_max": 0.0,
        "wry_max": 0.0,
        "wrz_max": 0.0,
        "rows": [],
    }

    if not ref_data["Time"] or not sim_data["Time"]:
        return stats

    sim_t0 = sim_data["Time"][0]
    sim_t1 = sim_data["Time"][-1]

    for i, t in enumerate(ref_data["Time"]):
        if t < sim_t0 or t > sim_t1:
            stats["rows_extrapolated_boundary"] += 1

        s_wrx = interpolate(t, sim_data["Time"], sim_data["Wrx"])
        s_wry = interpolate(t, sim_data["Time"], sim_data["Wry"])
        s_wrz = interpolate(t, sim_data["Time"], sim_data["Wrz"])

        r_wrx = ref_data["Wrx"][i]
        r_wry = ref_data["Wry"][i]
        r_wrz = ref_data["Wrz"][i]

        e_wrx = abs(s_wrx - r_wrx)
        e_wry = abs(s_wry - r_wry)
        e_wrz = abs(s_wrz - r_wrz)

        stats["wrx_rmse"] += e_wrx * e_wrx
        stats["wry_rmse"] += e_wry * e_wry
        stats["wrz_rmse"] += e_wrz * e_wrz
        stats["wrx_max"] = max(stats["wrx_max"], e_wrx)
        stats["wry_max"] = max(stats["wry_max"], e_wry)
        stats["wrz_max"] = max(stats["wrz_max"], e_wrz)
        stats["rows_compared"] += 1

        stats["rows"].append([
            t,
            r_wrx,
            s_wrx,
            e_wrx,
            r_wry,
            s_wry,
            e_wry,
            r_wrz,
            s_wrz,
            e_wrz,
        ])

    if stats["rows_compared"] > 0:
        n = stats["rows_compared"]
        stats["wrx_rmse"] = math.sqrt(stats["wrx_rmse"] / n)
        stats["wry_rmse"] = math.sqrt(stats["wry_rmse"] / n)
        stats["wrz_rmse"] = math.sqrt(stats["wrz_rmse"] / n)

    return stats


def save_pointwise_csv(rows, out_plot):
    out_csv = str(Path(out_plot).with_suffix("")) + "_pointwise.csv"
    with open(out_csv, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow([
            "Time_Ref",
            "Ref_Wrx",
            "Sim_Wrx_Interp",
            "AbsErr_Wrx",
            "Ref_Wry",
            "Sim_Wry_Interp",
            "AbsErr_Wry",
            "Ref_Wrz",
            "Sim_Wrz_Interp",
            "AbsErr_Wrz",
        ])
        w.writerows(rows)
    return out_csv


def set_times_font(plt):
    from matplotlib import font_manager

    times_paths = [
        r"C:/Windows/Fonts/times.ttf",
        r"C:/Windows/Fonts/timesbd.ttf",
        r"C:/Windows/Fonts/timesi.ttf",
        r"C:/Windows/Fonts/timesbi.ttf",
        "/mnt/c/Windows/Fonts/times.ttf",
        "/mnt/c/Windows/Fonts/timesbd.ttf",
        "/mnt/c/Windows/Fonts/timesi.ttf",
        "/mnt/c/Windows/Fonts/timesbi.ttf",
    ]

    for fp in times_paths:
        if os.path.exists(fp):
            try:
                font_manager.fontManager.addfont(fp)
            except Exception:
                pass

    plt.rcParams["font.family"] = "serif"
    plt.rcParams["font.serif"] = ["Times New Roman", "Times", "DejaVu Serif"]
    plt.rcParams["axes.unicode_minus"] = False


def process(sim_file, ref_file, out_plot, column_map):
    sim_data = read_gear_csv(sim_file, column_map, "simulation")
    ref_data = read_gear_csv(ref_file, column_map, "reference")

    # Gear22 reference provides angular velocity only; derive angular acceleration numerically.
    sim_awrx = numerical_derivative(sim_data["Time"], sim_data["Wrx"])
    sim_awry = numerical_derivative(sim_data["Time"], sim_data["Wry"])
    sim_awrz = numerical_derivative(sim_data["Time"], sim_data["Wrz"])
    ref_awrx = numerical_derivative(ref_data["Time"], ref_data["Wrx"])
    ref_awry = numerical_derivative(ref_data["Time"], ref_data["Wry"])
    ref_awrz = numerical_derivative(ref_data["Time"], ref_data["Wrz"])

    stats = compare_all_rows(sim_data, ref_data)

    print("\n==============================")
    print(" SIMPLE GEAR ERROR METRICS")
    print("==============================")
    print(" Basis: compare all reference rows (sim interpolated at ref time)")
    print(f" Ref Total Rows     : {stats['rows_ref']}")
    print(f" Compared Rows      : {stats['rows_compared']}")
    print(f" Boundary-Extrapolated Rows: {stats['rows_extrapolated_boundary']}")
    print(f" Wrx RMSE           : {stats['wrx_rmse']:.6f} rad/s")
    print(f" Wrx Max Abs Error  : {stats['wrx_max']:.6f} rad/s")
    print(f" Wry RMSE           : {stats['wry_rmse']:.6e} rad/s")
    print(f" Wry Max Abs Error  : {stats['wry_max']:.6e} rad/s")
    print(f" Wrz RMSE           : {stats['wrz_rmse']:.6e} rad/s")
    print(f" Wrz Max Abs Error  : {stats['wrz_max']:.6e} rad/s")
    print("==============================\n")

    pointwise_csv = save_pointwise_csv(stats["rows"], out_plot)
    print(f"Pointwise comparison CSV saved to: {pointwise_csv}")

    try:
        import matplotlib.pyplot as plt

        set_times_font(plt)

        plt.figure(figsize=(10, 16))

        plt.subplot(6, 1, 1)
        plt.plot(sim_data["Time"], sim_data["Wrx"], "b-", label="Sim Wrx")
        plt.plot(ref_data["Time"], ref_data["Wrx"], "r--", label="Ref Wrx")
        plt.title("Gear22 Angular Velocity RX Comparison")
        plt.ylabel("Angular Velocity (rad/s)")
        plt.legend()
        plt.grid(True)

        plt.subplot(6, 1, 2)
        plt.plot(sim_data["Time"], sim_data["Wry"], "b-", label="Sim Wry")
        plt.plot(ref_data["Time"], ref_data["Wry"], "r--", label="Ref Wry")
        plt.title("Gear22 Angular Velocity RY Comparison")
        plt.ylabel("Angular Velocity (rad/s)")
        plt.legend()
        plt.grid(True)

        plt.subplot(6, 1, 3)
        plt.plot(sim_data["Time"], sim_data["Wrz"], "b-", label="Sim Wrz")
        plt.plot(ref_data["Time"], ref_data["Wrz"], "r--", label="Ref Wrz")
        plt.title("Gear22 Angular Velocity RZ Comparison")
        plt.ylabel("Angular Velocity (rad/s)")
        plt.legend()
        plt.grid(True)

        plt.subplot(6, 1, 4)
        plt.plot(sim_data["Time"], sim_awrx, "b-", label="Sim Awrx (dWrx/dt)")
        plt.plot(ref_data["Time"], ref_awrx, "r--", label="Ref Awrx (dWrx/dt)")
        plt.title("Gear22 Angular Acceleration RX Comparison (Numerical Derivative)")
        plt.ylabel("Angular Accel (rad/s^2)")
        plt.legend()
        plt.grid(True)

        plt.subplot(6, 1, 5)
        plt.plot(sim_data["Time"], sim_awry, "b-", label="Sim Awry (dWry/dt)")
        plt.plot(ref_data["Time"], ref_awry, "r--", label="Ref Awry (dWry/dt)")
        plt.title("Gear22 Angular Acceleration RY Comparison (Numerical Derivative)")
        plt.ylabel("Angular Accel (rad/s^2)")
        plt.legend()
        plt.grid(True)

        plt.subplot(6, 1, 6)
        plt.plot(sim_data["Time"], sim_awrz, "b-", label="Sim Awrz (dWrz/dt)")
        plt.plot(ref_data["Time"], ref_awrz, "r--", label="Ref Awrz (dWrz/dt)")
        plt.title("Gear22 Angular Acceleration RZ Comparison (Numerical Derivative)")
        plt.xlabel("Time (s)")
        plt.ylabel("Angular Accel (rad/s^2)")
        plt.legend()
        plt.grid(True)

        plt.tight_layout()
        plt.savefig(out_plot, dpi=150)
        print(f"Comparison plot saved to: {out_plot}")
    except ImportError:
        print("Warning: matplotlib not installed, skip plotting.")


def main():
    parser = argparse.ArgumentParser(description="Compare simple gear Chrono result against RecurDyn Gear22.csv")
    parser.add_argument("--sim", default="data/outputs/baseline_simple_gear_nsc.csv", help="Path to simulation CSV")
    parser.add_argument("--ref", default="assets/simple_gear/data/Gear22.csv", help="Path to reference CSV")
    parser.add_argument("--map", default="data/reference/column_maps/simple_gear.json", help="Path to column map JSON")
    parser.add_argument("--out-plot", default="data/outputs/plots/simple_gear_comparison.png", help="Output plot path")
    args = parser.parse_args()

    with open(args.map, "r", encoding="utf-8") as f:
        column_map = json.load(f)

    process(args.sim, args.ref, args.out_plot, column_map)


if __name__ == "__main__":
    main()
