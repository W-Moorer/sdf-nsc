import os
import glob
import re

tex_files = glob.glob('papers/sections/*.tex') + ['papers/main.tex']

for filepath in tex_files:
    if os.path.exists(filepath):
        with open(filepath, 'r', encoding='utf-8') as f:
            content = f.read()

        # Fix $C$C^0$
        content = content.replace('$C$C^0$', '$C^0$')
        content = content.replace('$C^0$C^0$', '$C^0$')
        content = content.replace('$C^0$连续性', '$C^0$ 连续性')

        # Find missing $ or Chinese characters in math environments
        # Like: "=(V, E)$" -> "$=(V, E)$" ?
        content = content.replace('离散的互不连通且表面法线梯度分布连续（$\nabla \Phi_{i} \cdot \nabla \Phi_{j} > \epsilon$）的子流形', '离散的互不连通且表面法线梯度分布连续 ( $\\nabla \\Phi_i \\cdot \\nabla \\Phi_j > \\epsilon$ ) 的子流形')
        
        # fix: Missing character: There is no 不 (U+4E0D) in font [lmroman12-regular]
        # This usually means Chinese inside math mode (e.g. \text{} is not used or similar).
        
        with open(filepath, 'w', encoding='utf-8') as f:
            f.write(content)
print("cleaned")
