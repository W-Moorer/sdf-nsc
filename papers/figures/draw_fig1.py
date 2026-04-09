import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.ndimage import distance_transform_edt
import matplotlib.patches as patches
from matplotlib.lines import Line2D

def generate_academic_sdf_diagram():
    # ==========================================
    # 0. 全局字体与样式设置 (Times New Roman)
    # ==========================================
    plt.rcParams['font.family'] = 'serif'
    plt.rcParams['font.serif'] = ['Times New Roman'] + plt.rcParams['font.serif']
    plt.rcParams['mathtext.fontset'] = 'stix' 

    # ==========================================
    # 1. 参数设置
    # ==========================================
    N_high = 300       
    N_low = 18         
    tau = 0.28         
    
    x = np.linspace(-1.5, 1.5, N_high)
    y = np.linspace(-1.5, 1.5, N_high)
    X, Y = np.meshgrid(x, y)

    # 生成基础几何形状
    theta = np.arctan2(Y, X)
    r_shape = 0.8 + 0.15 * np.sin(4 * theta) + 0.1 * np.cos(3 * theta)
    r_grid = np.sqrt(X**2 + Y**2)
    inside_mask = r_grid < r_shape

    # 计算精确的 SDF (内部负，外部正)
    dist_inside = distance_transform_edt(inside_mask)
    dist_outside = distance_transform_edt(~inside_mask)
    dx = x[1] - x[0]
    sdf = (dist_outside - dist_inside) * dx

    # ==========================================
    # 2. 布局设计 (保持紧凑度)
    # ==========================================
    fig = plt.figure(figsize=(12, 6.5))
    gs = gridspec.GridSpec(2, 2, width_ratios=[1.3, 1], height_ratios=[1, 1.4], wspace=0.05, hspace=0.1)
    
    ax_main = fig.add_subplot(gs[:, 0])      
    ax_legend = fig.add_subplot(gs[0, 1])    
    ax_voxel = fig.add_subplot(gs[1, 1])     

    # ==========================================
    # 3. 左侧主图 (连续 SDF)
    # ==========================================
    ax_main.set_xticks(np.linspace(-1.5, 1.5, 13))
    ax_main.set_yticks(np.linspace(-1.5, 1.5, 13))
    ax_main.grid(True, color='gray', linestyle='-', linewidth=0.5, alpha=0.3, zorder=1)

    cmap = plt.cm.RdBu
    cf = ax_main.contourf(X, Y, sdf, levels=np.linspace(-1.5, 1.5, 60), 
                          cmap=cmap, vmin=-1.5, vmax=1.5, alpha=0.6, zorder=2)
    
    ax_main.contour(X, Y, sdf, levels=[0], colors='black', linewidths=2, zorder=4)
    ax_main.contour(X, Y, sdf, levels=[-tau, tau], colors='black', linestyles='dashed', linewidths=1.2, zorder=4)

    # 【修改点1 & 2】：删除所有数字，精准放置 d=0 和 +-tau 标签
    # 设置文本底框，遮盖线条，提高清晰度
    bbox_props = dict(boxstyle="round,pad=0.1", fc="white", ec="none", alpha=0.8)
    
    # 将标签精确放置在右下角线条上
    ax_main.text(0.51, -0.51, r'$d=0$', fontsize=14, zorder=5, bbox=bbox_props, ha='center', va='center')
    ax_main.text(0.74, -0.74, r'$d=+\tau$', fontsize=13, zorder=5, bbox=bbox_props, ha='center', va='center')
    ax_main.text(0.30, -0.30, r'$d=-\tau$', fontsize=13, zorder=5, bbox=bbox_props, ha='center', va='center')

    # 指示线：Narrow-Band Area
    ax_main.annotate('Narrow-Band Area', 
                     xy=(-0.45, 0.95), xycoords='data',    
                     xytext=(-1.2, 1.4), textcoords='data', 
                     arrowprops=dict(arrowstyle="-|>", color='black', lw=1, connectionstyle="arc3,rad=-0.2"),
                     fontsize=14, ha='center')

    ax_main.set_aspect('equal')
    ax_main.set_xlim(-1.6, 1.6)
    ax_main.set_ylim(-1.6, 1.6)
    ax_main.set_xticklabels([]) 
    ax_main.set_yticklabels([])
    ax_main.tick_params(length=0) 

    # ==========================================
    # 4. 右上角：图例区
    # ==========================================
    ax_legend.axis('off') 
    
    legend_elements = [
        Line2D([0], [0], color='black', lw=2, label=r'Geometric Boundary ($d=0$)'),
        Line2D([0], [0], color='black', lw=1.2, linestyle='--', label=r'Narrow-Band Thresholds ($\pm\tau$)'),
        patches.Patch(facecolor='#a6cee3', edgecolor='none', label='Signed Distance Visualization'),
        patches.Patch(facecolor='darkgray', edgecolor='black', hatch='xxxx', label='Stored/Active Voxel'),
        patches.Patch(facecolor='none', edgecolor='lightgray', label='Pruned/Not Stored Voxel')
    ]
    
    # 【修改点3】：上移标题 (y 从 0.85 提升到 0.98)，拉开与图例的距离
    ax_legend.text(0.5, 0.98, "Sparse Narrow-Band SDF\n(SNB-SDF)", 
                   fontsize=16, fontweight='bold', ha='center', va='top', transform=ax_legend.transAxes)
    
    # 保持图例主体在偏下方
    ax_legend.legend(handles=legend_elements, loc='center', frameon=False, 
                     fontsize=13, borderaxespad=0., handlelength=2.5, bbox_to_anchor=(0.5, 0.4))

    # ==========================================
    # 5. 右下角：稀疏体素图与指示说明
    # ==========================================
    x_low = np.linspace(-1.5, 1.5, N_low)
    y_low = np.linspace(-1.5, 1.5, N_low)
    X_low, Y_low = np.meshgrid(x_low, y_low)
    
    theta_low = np.arctan2(Y_low, X_low)
    r_shape_low = 0.8 + 0.15 * np.sin(4 * theta_low) + 0.1 * np.cos(3 * theta_low)
    r_grid_low = np.sqrt(X_low**2 + Y_low**2)
    dist_in_low = distance_transform_edt(r_grid_low < r_shape_low)
    dist_out_low = distance_transform_edt(~(r_grid_low < r_shape_low))
    dx_low = x_low[1] - x_low[0]
    sdf_low = (dist_out_low - dist_in_low) * dx_low

    ax_voxel.contour(X, Y, sdf, levels=[0], colors='black', linewidths=1, zorder=5)

    grid_size = dx_low

    for i in range(N_low - 1):
        for j in range(N_low - 1):
            center_sdf = sdf_low[i, j]
            rect_x = x_low[j] - grid_size/2
            rect_y = y_low[i] - grid_size/2
            
            if abs(center_sdf) <= tau:
                rect = patches.Rectangle((rect_x, rect_y), grid_size, grid_size, 
                                         linewidth=0.8, edgecolor='black', facecolor='darkgray', 
                                         hatch='xxxx', alpha=0.8)
                ax_voxel.add_patch(rect)
            else:
                rect = patches.Rectangle((rect_x, rect_y), grid_size, grid_size, 
                                         linewidth=0.5, edgecolor='lightgray', facecolor='none')
                ax_voxel.add_patch(rect)

    ax_voxel.annotate("Stored Voxel\n(Active)", 
                      xy=(-0.8, -0.6), xycoords='data',
                      xytext=(-1.5, -1.3), textcoords='data',
                      arrowprops=dict(arrowstyle="-|>", color="black", connectionstyle="arc3,rad=-0.2"),
                      fontsize=12, ha='center', va='center')

    ax_voxel.annotate("Pruned Voxel\n(Empty)", 
                      xy=(0, -1.2), xycoords='data',
                      xytext=(0.8, -1.4), textcoords='data',
                      arrowprops=dict(arrowstyle="-|>", color="black", connectionstyle="arc3,rad=0.2"),
                      fontsize=12, ha='center', va='center')
    
    ax_voxel.annotate("Narrow-Band\nBoundary", 
                      xy=(0.6, 0.6), xycoords='data',
                      xytext=(1.2, 0.9), textcoords='data',
                      arrowprops=dict(arrowstyle="-|>", color="black", connectionstyle="arc3,rad=-0.3"),
                      fontsize=12, ha='left', va='center')

    ax_voxel.set_xlim(-1.6, 1.6)
    ax_voxel.set_ylim(-1.6, 1.6)
    ax_voxel.set_aspect('equal')
    ax_voxel.axis('off') 

    # ==========================================
    # 6. 保存与展示
    # ==========================================
    plt.tight_layout()
    plt.savefig('SNB_SDF_Diagram_Compact_V2.pdf', dpi=600, bbox_inches='tight', pad_inches=0.05)
    plt.show()

if __name__ == "__main__":
    generate_academic_sdf_diagram()