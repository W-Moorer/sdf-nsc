import csv
import matplotlib.pyplot as plt

def load_data(filepath, col_t=0, col_pos=10, col_vel=13):
    time, pos, vel = [], [], []
    with open(filepath, 'r', encoding='utf-8') as f:
        reader = csv.reader(f)
        next(reader) # skip header
        for row in reader:
            if not row or len(row) <= max(col_t, col_pos, col_vel): continue
            try:
                time.append(float(row[col_t]))
                pos.append(float(row[col_pos]))
                vel.append(float(row[col_vel]))
            except ValueError:
                pass
    return time, pos, vel

ref_t, ref_pos, ref_vel = load_data('assets/cam/data/cam_data.csv')
base_t, base_pos, base_vel = load_data('data/outputs/baseline_cam_nsc_baseline.csv')
sdf_t, sdf_pos, sdf_vel = load_data('data/outputs/baseline_cam_nsc.csv')

plt.figure(figsize=(10, 8))

plt.subplot(2, 1, 1)
plt.plot(ref_t, ref_pos, 'k--', label='Reference (RecurDyn)')
plt.plot(base_t, base_pos, 'b-', alpha=0.7, label='Baseline NSC (Mesh)')
plt.plot(sdf_t, sdf_pos, 'r-', alpha=0.7, label='SDF NSC (NanoVDB)')
plt.title('Cam Follower Position Y over Time')
plt.xlabel('Time (s)')
plt.ylabel('Pos Y (m)')
plt.legend()
plt.grid(True)

plt.subplot(2, 1, 2)
plt.plot(ref_t, ref_vel, 'k--', label='Reference (RecurDyn)')
plt.plot(base_t, base_vel, 'b-', alpha=0.7, label='Baseline NSC (Mesh)')
plt.plot(sdf_t, sdf_vel, 'r-', alpha=0.7, label='SDF NSC (NanoVDB)')
plt.title('Cam Follower Velocity Y over Time')
plt.xlabel('Time (s)')
plt.ylabel('Vel Y (m/s)')
plt.legend()
plt.grid(True)

plt.tight_layout()
plt.savefig('data/outputs/plots/cam_comparison_3lines.png', dpi=300)
print('Plot saved to data/outputs/plots/cam_comparison_3lines.png')
