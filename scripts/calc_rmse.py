import numpy as np

def load_csv(path, wrx_col_idx):
    times, wrxs = [], []
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
    return np.array(times), np.array(wrxs)

t_ref, wrx_ref = load_csv('assets/simple_gear/data/Gear22.csv', 6)
t_1, wrx_1 = load_csv('data/outputs/gear_1st_highres_0001.csv', 3)
t_2, wrx_2 = load_csv('data/outputs/gear_2nd_highres_0001.csv', 3)
t_3, wrx_3 = load_csv('data/outputs/gear_2nd_highres_0001_adaptive.csv', 3)

def calc_rmse_steady(t_arr, arr, target=-1.0):
    mask = t_arr >= 0.00
    if sum(mask) == 0: return float('nan')
    return np.sqrt(np.mean((arr[mask] - target)**2))

print('=== RMSE (Relative to -1.0) ===')
print('Reference     :', f'{calc_rmse_steady(t_ref, wrx_ref):.5f}')
print('1st-order     :', f'{calc_rmse_steady(t_1, wrx_1):.5f}')
print('2nd-order     :', f'{calc_rmse_steady(t_2, wrx_2):.5f}')
print('Adaptive      :', f'{calc_rmse_steady(t_3, wrx_3):.5f}')
