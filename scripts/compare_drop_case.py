#!/usr/bin/env python3
import os
import sys
import argparse
import csv
import math
import json
import re
from pathlib import Path

def read_standard_csv(filepath):
    data = {'Time': [], 'Pos_Y': [], 'Vel_Y': [], 'Acc_Y': [], 'Num_Contacts': []}
    has_contacts = False
    with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
        reader = csv.reader(f)
        header = next(reader)
        header_norm = [h.strip().lower() for h in header]

        if 'num_contacts' in header_norm:
            has_contacts = True

        time_idx = header_norm.index('time') if 'time' in header_norm else 0
        pos_idx = header_norm.index('pos_y') if 'pos_y' in header_norm else 1
        vel_idx = header_norm.index('vel_y') if 'vel_y' in header_norm else 2
        acc_idx = header_norm.index('acc_y') if 'acc_y' in header_norm else None
        cnt_idx = header_norm.index('num_contacts') if has_contacts else None

        for row in reader:
            if not row: continue
            try:
                data['Time'].append(float(row[time_idx]))
                data['Pos_Y'].append(float(row[pos_idx]))
                data['Vel_Y'].append(float(row[vel_idx]))
                if acc_idx is not None and acc_idx < len(row):
                    data['Acc_Y'].append(float(row[acc_idx]))
                if cnt_idx is not None and cnt_idx < len(row):
                    data['Num_Contacts'].append(int(float(row[cnt_idx])))
            except ValueError:
                pass

    # Legacy fallback: derive acceleration from velocity if acc column is absent.
    if not data['Acc_Y'] and len(data['Time']) == len(data['Vel_Y']) and len(data['Time']) > 1:
        data['Acc_Y'].append((data['Vel_Y'][1] - data['Vel_Y'][0]) / (data['Time'][1] - data['Time'][0]))
        for i in range(1, len(data['Time'])):
            dt = data['Time'][i] - data['Time'][i - 1]
            if dt == 0:
                data['Acc_Y'].append(data['Acc_Y'][-1])
            else:
                data['Acc_Y'].append((data['Vel_Y'][i] - data['Vel_Y'][i - 1]) / dt)

    return data

def _find_column_index(header, patterns, model_tag=None):
    if not patterns:
        return None

    compiled = [re.compile(p, re.IGNORECASE) for p in patterns]
    model_tag_l = (model_tag or "").lower().strip()

    candidates = []
    for idx, col_name in enumerate(header):
        col_l = col_name.lower()
        if all(p.search(col_name) for p in compiled):
            score = 1
            if model_tag_l and model_tag_l in col_l:
                score = 2
            candidates.append((score, idx))

    if not candidates:
        return None

    candidates.sort(reverse=True)
    return candidates[0][1]

def _resolve_indices(header, column_map):
    model_tag = column_map.get("model_tag", "")

    time_idx = _find_column_index(header, column_map.get("time_col_patterns"), model_tag)
    pos_idx = _find_column_index(header, column_map.get("pos_y_col_patterns"), model_tag)
    vel_idx = _find_column_index(header, column_map.get("vel_y_col_patterns"), model_tag)
    acc_idx = _find_column_index(header, column_map.get("acc_y_col_patterns"), model_tag)
    cnt_idx = _find_column_index(header, column_map.get("num_contacts_col_patterns"), model_tag)

    # Backward-compatible fallback for legacy fixed-index maps.
    if time_idx is None:
        time_idx = column_map.get("time_col_idx", 0)
    if pos_idx is None:
        pos_idx = column_map.get("pos_y_col_idx", 1)
    if vel_idx is None:
        vel_idx = column_map.get("vel_y_col_idx", 2)
    if acc_idx is None:
        acc_idx = column_map.get("acc_y_col_idx", None)
    if cnt_idx is None:
        cnt_idx = column_map.get("num_contacts_col_idx", None)

    return time_idx, pos_idx, vel_idx, acc_idx, cnt_idx

def read_mapped_csv(filepath, column_map, source_name="reference"):
    data = {'Time': [], 'Pos_Y': [], 'Vel_Y': [], 'Acc_Y': [], 'Num_Contacts': []}

    # Some of the external CSVs might have Chinese chars in header
    with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
        reader = csv.reader(f)
        header = next(reader)
        t_idx, p_idx, v_idx, a_idx, c_idx = _resolve_indices(header, column_map)

        if t_idx >= len(header) or p_idx >= len(header) or v_idx >= len(header):
            raise ValueError("Column index out of range after resolving mapping.")

        print(f"Resolved {source_name} columns:")
        print(f" - Time : [{t_idx}] {header[t_idx]}")
        print(f" - Pos_Y: [{p_idx}] {header[p_idx]}")
        print(f" - Vel_Y: [{v_idx}] {header[v_idx]}")
        if a_idx is not None and a_idx < len(header):
            print(f" - Acc_Y: [{a_idx}] {header[a_idx]}")

        for row in reader:
            if not row: continue
            try:
                data['Time'].append(float(row[t_idx]))
                data['Pos_Y'].append(float(row[p_idx]))
                data['Vel_Y'].append(float(row[v_idx]))
                if a_idx is not None and a_idx < len(row):
                    data['Acc_Y'].append(float(row[a_idx]))
                if c_idx is not None and c_idx < len(row):
                    data['Num_Contacts'].append(int(float(row[c_idx])))
            except ValueError:
                pass # skip faulty rows
                
    return data

def read_sim_csv(sim_file, column_map=None):
    with open(sim_file, 'r', encoding='utf-8', errors='ignore') as f:
        reader = csv.reader(f)
        header = next(reader)

    header_norm = [h.strip().lower() for h in header]
    if 'time' in header_norm and 'pos_y' in header_norm and 'vel_y' in header_norm:
        return read_standard_csv(sim_file)

    if column_map:
        print("Using column map for simulation parsing...")
        return read_mapped_csv(sim_file, column_map, source_name="simulation")

    # Last fallback for legacy 3-column/4-column standard format without headers.
    return read_standard_csv(sim_file)

def interpolate(t, times, values):
    if not times: return 0.0
    if t <= times[0]: return values[0]
    if t >= times[-1]: return values[-1]

    for i in range(len(times)-1):
        if times[i] <= t <= times[i+1]:
            t0, t1 = times[i], times[i+1]
            v0, v1 = values[i], values[i+1]
            if t1 == t0: return v0
            factor = (t - t0) / (t1 - t0)
            return v0 + factor * (v1 - v0)
    return values[-1]

def compare_all_reference_rows(sim_data, ref_data):
    stats = {
        'count_ref_total': len(ref_data['Time']),
        'count_compared': 0,
        'count_extrapolated_boundary': 0,
        'pos_rmse': 0.0,
        'pos_max_err': 0.0,
        'vel_rmse': 0.0,
        'vel_max_err': 0.0,
        'acc_rmse': 0.0,
        'acc_max_err': 0.0,
        'pointwise_rows': [],
    }

    if not sim_data['Time'] or not ref_data['Time']:
        return stats

    sim_t0 = sim_data['Time'][0]
    sim_t1 = sim_data['Time'][-1]

    for i, t_ref in enumerate(ref_data['Time']):
        if t_ref < sim_t0 or t_ref > sim_t1:
            stats['count_extrapolated_boundary'] += 1

        sim_pos_interp = interpolate(t_ref, sim_data['Time'], sim_data['Pos_Y'])
        sim_vel_interp = interpolate(t_ref, sim_data['Time'], sim_data['Vel_Y'])

        ref_pos = ref_data['Pos_Y'][i]
        ref_vel = ref_data['Vel_Y'][i]
        ref_acc = ref_data['Acc_Y'][i] if i < len(ref_data['Acc_Y']) else None

        err_p = abs(sim_pos_interp - ref_pos)
        err_v = abs(sim_vel_interp - ref_vel)
        err_a = None
        sim_acc_interp = None
        if sim_data['Acc_Y'] and ref_acc is not None:
            sim_acc_interp = interpolate(t_ref, sim_data['Time'], sim_data['Acc_Y'])
            err_a = abs(sim_acc_interp - ref_acc)

        stats['pos_rmse'] += err_p * err_p
        stats['vel_rmse'] += err_v * err_v
        stats['pos_max_err'] = max(stats['pos_max_err'], err_p)
        stats['vel_max_err'] = max(stats['vel_max_err'], err_v)
        if err_a is not None:
            stats['acc_rmse'] += err_a * err_a
            stats['acc_max_err'] = max(stats['acc_max_err'], err_a)
        stats['count_compared'] += 1

        stats['pointwise_rows'].append([
            t_ref,
            ref_pos,
            sim_pos_interp,
            err_p,
            ref_vel,
            sim_vel_interp,
            err_v,
            ref_acc,
            sim_acc_interp,
            err_a,
        ])

    if stats['count_compared'] > 0:
        stats['pos_rmse'] = math.sqrt(stats['pos_rmse'] / stats['count_compared'])
        stats['vel_rmse'] = math.sqrt(stats['vel_rmse'] / stats['count_compared'])
        if sim_data['Acc_Y'] and ref_data['Acc_Y']:
            stats['acc_rmse'] = math.sqrt(stats['acc_rmse'] / stats['count_compared'])

    return stats

def save_pointwise_csv(pointwise_rows, output_plot):
    output_csv = str(Path(output_plot).with_suffix('')) + '_pointwise.csv'
    with open(output_csv, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow([
            'Time_Ref',
            'Ref_Pos_Y',
            'Sim_Pos_Y_Interp',
            'AbsErr_Pos_Y',
            'Ref_Vel_Y',
            'Sim_Vel_Y_Interp',
            'AbsErr_Vel_Y',
            'Ref_Acc_Y',
            'Sim_Acc_Y_Interp',
            'AbsErr_Acc_Y',
        ])
        writer.writerows(pointwise_rows)
    return output_csv

def process_and_plot(sim_file, ref_file, output_plot, column_map=None):
    if not os.path.exists(sim_file):
        print(f"Error: Sim file not found: {sim_file}")
        sys.exit(1)
    if not os.path.exists(ref_file):
        print(f"Error: Ref file not found: {ref_file}")
        sys.exit(1)

    sim_data = read_sim_csv(sim_file, column_map)
    
    if column_map:
        print("Using column map for reference parsing...")
        ref_data = read_mapped_csv(ref_file, column_map, source_name="reference")
    else:
        ref_data = read_standard_csv(ref_file)

    stats = compare_all_reference_rows(sim_data, ref_data)
    if stats['count_compared'] == 0:
        print("Warning: No overlapping time region found for comparison.")

    print("\n==============================")
    print(" ERROR METRICS vs REFERENCE")
    print("==============================")
    print(" Basis: all rows in reference data file (interpolate sim at each ref timestamp)")
    print(f" Ref Total Rows     : {stats['count_ref_total']}")
    print(f" Compared Rows      : {stats['count_compared']}")
    print(f" Boundary-Extrapolated Rows: {stats['count_extrapolated_boundary']}")
    print(f" Pos_Y RMSE         : {stats['pos_rmse']:.6f} m")
    print(f" Pos_Y Max Abs Error: {stats['pos_max_err']:.6f} m")
    print(f" Vel_Y RMSE         : {stats['vel_rmse']:.6f} m/s")
    print(f" Vel_Y Max Abs Error: {stats['vel_max_err']:.6f} m/s")
    if ref_data['Acc_Y'] and sim_data['Acc_Y']:
        print(f" Acc_Y RMSE         : {stats['acc_rmse']:.6f} m/s^2")
        print(f" Acc_Y Max Abs Error: {stats['acc_max_err']:.6f} m/s^2")
    print("==============================\n")

    pointwise_csv = save_pointwise_csv(stats['pointwise_rows'], output_plot)
    print(f"Pointwise comparison CSV saved to: {pointwise_csv}")

    try:
        import matplotlib.pyplot as plt
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

        plt.rcParams['font.family'] = 'serif'
        plt.rcParams['font.serif'] = ['Times New Roman', 'Times', 'DejaVu Serif']
        plt.rcParams['axes.unicode_minus'] = False

        has_acc = bool(ref_data['Acc_Y']) and bool(sim_data['Acc_Y'])
        if has_acc:
            plt.figure(figsize=(10, 11))
            nrows = 3
        else:
            plt.figure(figsize=(10, 8))
            nrows = 2

        plt.subplot(nrows, 1, 1)
        plt.plot(sim_data['Time'], sim_data['Pos_Y'], 'b-', label='Sim Pos_Y')
        plt.plot(ref_data['Time'], ref_data['Pos_Y'], 'r--', label='Ref Pos_Y')
        plt.title('Position Y Comparison')
        plt.ylabel('Position (m)')
        plt.legend()
        plt.grid(True)

        plt.subplot(nrows, 1, 2)
        plt.plot(sim_data['Time'], sim_data['Vel_Y'], 'b-', label='Sim Vel_Y')
        plt.plot(ref_data['Time'], ref_data['Vel_Y'], 'r--', label='Ref Vel_Y')
        plt.title('Velocity Y Comparison')
        plt.ylabel('Velocity (m/s)')
        plt.legend()
        plt.grid(True)

        if has_acc:
            plt.subplot(nrows, 1, 3)
            plt.plot(sim_data['Time'], sim_data['Acc_Y'], 'b-', label='Sim Acc_Y')
            plt.plot(ref_data['Time'], ref_data['Acc_Y'], 'r--', label='Ref Acc_Y')
            plt.title('Acceleration Y Comparison')
            plt.ylabel('Acceleration (m/s^2)')
            plt.legend()
            plt.grid(True)

        plt.xlabel('Time (s)')

        plt.tight_layout()
        plt.savefig(output_plot, dpi=150)
        print(f"Comparison plot saved to: {output_plot}")

    except ImportError:
        print("Warning: 'matplotlib' not found. Skipping plot generation.")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Compare baseline drop simulated CSV against reference CSV.')
    parser.add_argument('--case', help='Specify case name to auto-load configs (e.g. cam)', default=None)
    parser.add_argument('--sim', default='data/outputs/baseline_drop_nsc.csv', help='Path to simulation CSV')
    parser.add_argument('--ref', default='data/reference/baseline_drop_nsc_reference.csv', help='Path to reference CSV')
    parser.add_argument('--out-plot', default='data/outputs/plots/drop_nsc_comparison.png', help='Path to output plot PNG')

    args = parser.parse_args()

    ref_file = args.ref
    sim_file = args.sim
    out_plot = args.out_plot
    column_map = None

    if args.case:
        print(f"Loading settings for case: {args.case}")
        config_path = f"data/reference/case_configs/{args.case}.json"
        map_path = f"data/reference/column_maps/{args.case}.json"
        
        if os.path.exists(config_path):
            with open(config_path, 'r') as f:
                cfg = json.load(f)
                ref_file = cfg.get("result_file", args.ref)
                
        if os.path.exists(map_path):
            with open(map_path, 'r') as f:
                column_map = json.load(f)
                
        out_plot = f"data/outputs/plots/{args.case}_comparison.png"

    process_and_plot(sim_file, ref_file, out_plot, column_map)
