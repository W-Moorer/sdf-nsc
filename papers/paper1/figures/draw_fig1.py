import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.ndimage import distance_transform_edt
import matplotlib.patches as patches
from matplotlib.lines import Line2D

def generate_complex_2d_sdf_diagram():
    # 0. 字体设置
    plt.rcParams['font.family'] = 'serif'
    plt.rcParams['font.serif'] = ['Times New Roman'] + plt.rcParams['font.serif']
    plt.rcParams['mathtext.fontset'] = 'stix' 

    # 1. 参数设置
    N_high = 300       
    N_low = 18         
    tau = 0.28         
    
    x = np.linspace(-1.5, 1.5, N_high)
    y = np.linspace(-1.5, 1.5, N_high)
    X, Y = np.meshgrid(x, y)

    # ==========================================
    # 【核心修改区】生成更加复杂、高频的几何曲线
    # ==========================================
    theta = np.arctan2(Y, X)
    # 引入高频的正弦和余弦项，形成复杂的“变形虫”或“星形齿轮”状边界
    r_shape = 0.7 + 0.2 * np.sin(5 * theta) + 0.15 * np.cos(6 * theta) + 0.08 * np.sin(10 * theta)
    inside_mask = np.sqrt(X**2 + Y**2) < r_shape

    # 计算精确 SDF
    dist_inside = distance_transform_edt(inside_mask)
    dist_outside = distance_transform_edt(~inside_mask)
    dx = x[1] - x[0]
    sdf = (dist_outside - dist_inside) * dx

    # 2. 布局设计
    fig = plt.figure(figsize=(12, 6.5))
    gs = gridspec.GridSpec(2, 2, width_ratios=[1.3, 1], height_ratios=[1, 1.4], wspace=0.05, hspace=0.1)
    
    ax_main = fig.add_subplot(gs[:, 0])      
    ax_legend = fig.add_subplot(gs[0, 1])    
    ax_voxel = fig.add_subplot(gs[1, 1])     

    # ==========================================
    # 左侧主图：通过透明度极差，强力凸显“窄带”
    # ==========================================
    ax_main.set_xticks(np.linspace(-1.5, 1.5, 13))
    ax_main.set_yticks(np.linspace(-1.5, 1.5, 13))
    ax_main.grid(True, color='gray', linestyle='-', linewidth=0.5, alpha=0.3, zorder=1)
    cmap = plt.cm.RdBu

    # A. 绘制全局暗淡背景 (alpha=0.15)
    ax_main.contourf(X, Y, sdf, levels=np.linspace(-1.5, 1.5, 60), 
                     cmap=cmap, vmin=-1.5, vmax=1.5, alpha=0.15, zorder=2)
    
    # B. 掩膜提取窄带，绘制高亮窄带 (alpha=0.85)
    sdf_narrow = np.ma.masked_where(np.abs(sdf) > tau, sdf)
    ax_main.contourf(X, Y, sdf_narrow, levels=np.linspace(-tau, tau, 40), 
                     cmap=cmap, vmin=-1.5, vmax=1.5, alpha=0.85, zorder=3)

    # 线条
    ax_main.contour(X, Y, sdf, levels=[0], colors='black', linewidths=2, zorder=4)
    ax_main.contour(X, Y, sdf, levels=[-tau, tau], colors='black', linestyles='dashed', linewidths=1.2, zorder=4)

    # 标签 (带白底抗遮挡) - 坐标已根据新的复杂曲线进行了重新数学计算适配
    bbox_props = dict(boxstyle="round,pad=0.1", fc="white", ec="none", alpha=0.8)
    ax_main.text(0.54, -0.54, r'$d=0$', fontsize=14, zorder=5, bbox=bbox_props, ha='center', va='center')
    ax_main.text(0.74, -0.74, r'$d=+\tau$', fontsize=13, zorder=5, bbox=bbox_props, ha='center', va='center')
    ax_main.text(0.34, -0.34, r'$d=-\tau$', fontsize=13, zorder=5, bbox=bbox_props, ha='center', va='center')

    ax_main.annotate('Highlighted\nNarrow-Band', 
                     xy=(-0.55, 0.55), xycoords='data',    
                     xytext=(-1.3, 1.3), textcoords='data', 
                     arrowprops=dict(arrowstyle="-|>", color='black', lw=1.5, connectionstyle="arc3,rad=-0.2"),
                     fontsize=14, fontweight='bold', ha='center')

    ax_main.set_aspect('equal')
    ax_main.set_xlim(-1.6, 1.6)
    ax_main.set_ylim(-1.6, 1.6)
    ax_main.axis('off') 

    # ==========================================
    # 右上角：图例区
    # ==========================================
    ax_legend.axis('off') 
    legend_elements = [
        Line2D([0], [0], color='black', lw=2, label=r'Surface ($d=0$)'),
        Line2D([0], [0], color='black', lw=1.2, linestyle='--', label=r'Band Limits ($\pm\tau$)'),
        patches.Patch(facecolor='#1f77b4', edgecolor='black', hatch='////', alpha=0.8, label='Allocated Voxel (SDF Stored)'),
        patches.Patch(facecolor='none', edgecolor='lightgray', label='Empty Space (Memory Saved)')
    ]
    ax_legend.text(0.5, 0.98, "Sparse Narrow-Band SDF", 
                   fontsize=18, fontweight='bold', ha='center', va='top', transform=ax_legend.transAxes)
    ax_legend.legend(handles=legend_elements, loc='center', frameon=False, 
                     fontsize=13, borderaxespad=0., handlelength=2.5, bbox_to_anchor=(0.5, 0.4))

    # ==========================================
    # 右下角：高对比度稀疏体素图
    # ==========================================
    x_low = np.linspace(-1.5, 1.5, N_low)
    y_low = np.linspace(-1.5, 1.5, N_low)
    X_low, Y_low = np.meshgrid(x_low, y_low)
    
    r_grid_low = np.sqrt(X_low**2 + Y_low**2)
    theta_low = np.arctan2(Y_low, X_low)
    
    # 同步修改低分辨率网格的曲线生成函数
    r_shape_low = 0.7 + 0.2 * np.sin(5 * theta_low) + 0.15 * np.cos(6 * theta_low) + 0.08 * np.sin(10 * theta_low)
    dist_in_low = distance_transform_edt(r_grid_low < r_shape_low)
    dist_out_low = distance_transform_edt(~(r_grid_low < r_shape_low))
    sdf_low = (dist_out_low - dist_in_low) * (x_low[1] - x_low[0])

    ax_voxel.contour(X, Y, sdf, levels=[0], colors='black', linewidths=1.5, zorder=5)
    grid_size = x_low[1] - x_low[0]

    for i in range(N_low - 1):
        for j in range(N_low - 1):
            center_sdf = sdf_low[i, j]
            rect_x = x_low[j] - grid_size/2
            rect_y = y_low[i] - grid_size/2
            
            if abs(center_sdf) <= tau:
                rect = patches.Rectangle((rect_x, rect_y), grid_size, grid_size, 
                                         linewidth=1, edgecolor='black', facecolor='#1f77b4', 
                                         hatch='////', alpha=0.8)
                ax_voxel.add_patch(rect)
            else:
                rect = patches.Rectangle((rect_x, rect_y), grid_size, grid_size, 
                                         linewidth=0.3, edgecolor='#e0e0e0', facecolor='none')
                ax_voxel.add_patch(rect)

    ax_voxel.annotate("Allocated Data Structure\n(Sparse Memory)", 
                      xy=(-0.8, -0.6), xycoords='data',
                      xytext=(-1.5, -1.3), textcoords='data',
                      arrowprops=dict(arrowstyle="-|>", color="#1f77b4", lw=1.5, connectionstyle="arc3,rad=-0.2"),
                      fontsize=12, fontweight='bold', color='#1f77b4', ha='center', va='center')

    ax_voxel.annotate("Pruned\n(No Memory)", 
                      xy=(0, -1.2), xycoords='data',
                      xytext=(0.8, -1.4), textcoords='data',
                      arrowprops=dict(arrowstyle="-|>", color="gray", lw=1, connectionstyle="arc3,rad=0.2"),
                      fontsize=12, color='gray', ha='center', va='center')

    ax_voxel.set_xlim(-1.6, 1.6)
    ax_voxel.set_ylim(-1.6, 1.6)
    ax_voxel.set_aspect('equal')
    ax_voxel.axis('off') 

    plt.tight_layout()
    plt.savefig('SNB_SDF_Diagram.pdf', dpi=600, bbox_inches='tight', pad_inches=0.05)
    plt.show()

if __name__ == "__main__":
    generate_complex_2d_sdf_diagram()