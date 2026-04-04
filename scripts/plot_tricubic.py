import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

def load_csv(path, wrx_col_idx):
    times, wrxs = [], []
    with open(path, 'r') as f:
        # detect columns
        lines = f.readlines()
        for i, l in enumerate(lines):
            if 'Time' in l or 'Y:Vel_RX' in l or 'Y:Vel' in l:
                lines = lines[i+1:]
                break

        for l in lines:
            if 'Time' in l or '---' in l or 'Step' in l or '==' in l or 'Min' in l or 'Max' in l or 'Avg' in l or 'Output' in l or 'SUMMARY' in l or 'Final' in l: continue
            parts = l.strip().split(',')
            if len(parts) > max(0, wrx_col_idx):
                try:
                    t = float(parts[0])
                    wrx = float(parts[wrx_col_idx])
                    times.append(t)
                    wrxs.append(wrx)
                except Exception as e:
                    pass
    return times, wrxs

t_ref, wrx_ref = load_csv('assets/simple_gear/data/Gear22.csv', 6)
t_1, wrx_1 = load_csv('data/outputs/gear_tricubic_1st.csv', 3)
t_2, wrx_2 = load_csv('data/outputs/gear_tricubic_2nd.csv', 3)
t_3, wrx_3 = load_csv('data/outputs/gear_2nd_highres_0001_adaptive.csv', 3)
plt.figure(figsize=(7, 4))
plt.plot(t_ref, wrx_ref, 'k--', label='Reference (RecurDyn)')
plt.plot(t_1, wrx_1, 'r-', linewidth=1, label='1st-order Tricubic SDF')
plt.plot(t_2, wrx_2, 'g-', linewidth=2, label='2nd-order Tricubic SDF (Unfiltered)')
plt.plot(t_3, wrx_3, 'b-', linewidth=2, label='2nd-order Tricubic SDF (Adaptive Blending)')
plt.xlabel('Time (s)')
plt.ylabel('Velocity W_rx (rad/s)')
plt.title('Gear Meshing Angular Velocity (C1 Continuous SDF)')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig('data/outputs/plots/gear_tricubic_3lines.pdf')
plt.savefig('data/outputs/plots/gear_tricubic_3lines.png')
print('Plotted to data/outputs/plots/gear_tricubic_3lines.*')
