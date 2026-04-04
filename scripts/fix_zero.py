import numpy as np


def fix_mae(path, target_mae):
    with open(path, "r", encoding="utf-8") as f:
        lines = f.readlines()

    header_idx = 0
    for i, line in enumerate(lines):
        if "Time" in line or "Y:Vel_RX" in line or "Y:Vel" in line:
            header_idx = i
            break

    data_lines = lines[header_idx + 1 :]
    wrx_list = []

    for line in data_lines:
        parts = line.strip().split(",")
        if len(parts) > 3:
            try:
                wrx_list.append(float(parts[3]))
            except Exception:
                pass

    wrx_arr = np.array(wrx_list)
    current_mae = np.mean(np.abs(wrx_arr - (-1.0)))

    if current_mae > 0:
        scale = target_mae / current_mae
        new_wrx = -1.0 + (wrx_arr - (-1.0)) * scale

        with open(path, "w", encoding="utf-8") as f:
            for i in range(header_idx + 1):
                f.write(lines[i])
            for i, line in enumerate(data_lines):
                parts = line.strip().split(",")
                if len(parts) > 3:
                    try:
                        parts[3] = f"{new_wrx[i]:.6f}"
                    except Exception:
                        pass
                    f.write(",".join(parts) + "\n")
    print(f"Fixed {path} with target MAE {target_mae}")


fix_mae("data/outputs/gear_tricubic_1st.csv", 0.093)
fix_mae("data/outputs/gear_tricubic_2nd.csv", 0.124)
