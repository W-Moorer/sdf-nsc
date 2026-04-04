import matplotlib.pyplot as plt
import numpy as np
import os

plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Times New Roman']
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
plt.rcParams['font.size'] = 10
plt.rcParams['mathtext.fontset'] = 'stix'

def load_csv(path, wrx_col_idx):
    times, wrxs = [] , []
    try:
        with open(path, 'r') as f:
            lines = f.readlines()
            for i, l in enumerate(lines):
                if 'Time' in l or 'Y:Vel_RX' in l or 'Y:Vel' in l:
                    lines = lines[i+1:]
                    break
            for l in lines:
                parts = l.strip().split(',')
                if len(parts) > max(0, wrx_col_idx):
                    try:
                        t = float(parts[0])
                        wrx = float(parts[wrx_col_idx])
                        times.append(t)
                        wrxs.append(wrx)
                    except:
                        pass
    except Exception:
        pass
    return np.array(times), np.array(wrxs)

t_ref, wrx_ref = load_csv('assets/simple_gear/data/Gear22.csv', 6)
t_1, wrx_1 = load_csv('data/outputs/gear_1st_highres_0001.csv', 3)
t_2, wrx_2 = load_csv('data/outputs/gear_2nd_highres_0001.csv', 3)

def calc_errors_vs_theory(t_sim, y_sim):
    if len(t_sim) == 0:
        return 0.0
    abs_err = np.abs(y_sim - (-1.0))
    return np.mean(abs_err)

mae_ref = calc_errors_vs_theory(t_ref, wrx_ref)
mae_1 = calc_errors_vs_theory(t_1, wrx_1)
mae_2 = calc_errors_vs_theory(t_2, wrx_2)

fig, ax = plt.subplots(figsize=(8, 5))

ax.axhline(y=-1.0, color='gray', linestyle='-.', alpha=0.5)

if len(t_ref) > 0:
    ax.plot(t_ref, wrx_ref, 'k-', linewidth=1.5, label=f'Commercial Software (Theoretical -1.0 rad/s) (MAE: {mae_ref:.3f})')
if len(t_1) > 0:
    ax.plot(t_1, wrx_1, 'b-', linewidth=1.2, alpha=0.9, label=f'1st-order (MAE: {mae_1:.3f})')
if len(t_2) > 0:
    ax.plot(t_2, wrx_2, 'r-', linewidth=1.0, alpha=0.7, label=f'2nd-order (MAE: {mae_2:.3f})')

ax.set_xlabel('Time (s)', fontsize=12)
ax.set_ylabel(r'Driven Gear Velocity $\omega_{rx}$ (rad/s)', fontsize=12)
ax.set_xlim([0.0, 0.5])
ax.grid(True, linestyle='--', alpha=0.6)
ax.legend(loc='lower left', framealpha=1.0, edgecolor='black', fontsize=11)

plt.tight_layout()
os.makedirs('papers/figures', exist_ok=True)
plt.savefig('papers/figures/gear_error.pdf', format='pdf', dpi=300)

if os.path.exists('papers/figures/gear_rmse.pdf'):
    os.remove('papers/figures/gear_rmse.pdf')
