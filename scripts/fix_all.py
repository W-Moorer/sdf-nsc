import os

def fix_file(filepath):
    if not os.path.exists(filepath):
        return
    with open(filepath, 'r', encoding='utf-8') as f:
        content = f.read()
    
    # 03_methodology.tex specific
    import re
    content = re.sub(r'\\subsection\{WNbQb\~gR:SNAmb_M\~v\\\(uirt.*?sumption\)\}', '', content, flags=re.DOTALL)
    content = content.replace(r'$G$G=(V, E)$', r'$G=(V, E)$')
    content = content.replace(r'结合 $ 之间的弱耦合接触', r'结合所有连通域之间的弱耦合接触')

    # 05_experiments.tex specific broken maths
    d = {
        r' \times 10^{-4}$ m': r' $4 \times 10^{-4}$ m',
        r'.15 \sim 0.20$ 秒': r'$0.15 \sim 0.20$ 秒',
        r'$ 面震荡': r'$y$ 面震荡',
        r'{curve}$ 项': r'$C_{curve}$ 项',
        r'=0.2$ s': r'$t=0.2$ s',
        r' \to -231.8$ m': r'$z \to -231.8$ m',
        r'^2$ 调节模块': r'$(\Delta t)^2$ 调节模块',
        r'^1$ 连续性': r'$C^1$ 连续性',
        r' \times \text{voxel\_size}$': r'$2 \times \text{voxel\_size}$',
        r' \times 3$ 张量': r'$3 \times 3$ 张量',
        r' \sim 5$ 倍扩容': r'$2 \sim 5$ 倍扩容',
        r'^0$ 连续性': r'$C^0$ 连续性',
        r' \mu\text{m}$ 分辨率': r'$2 \mu\text{m}$ 分辨率',
        r'^2$ 连续的光滑抛物面': r'$C^2$ 连续的光滑抛物面',
        r'^2$ 连续的': r'$C^2$ 连续的',
        r'^0$ 奇异几何体': r'$C^0$ 奇异几何体',
        r'^0$ 边界）处': r'$C^0$ 边界）处'
    }
    for k, v in d.items():
        content = content.replace(k, v)
    
    # Remove TV
    content = content.replace('，以及总变差（TV）', '')
    content = content.replace('和 TV（总变差）', '')
    content = content.replace('、TV突变', '')
    content = content.replace('总变差(TV)', '')
    content = content.replace('总变差（TV）', '')
    
    with open(filepath, 'w', encoding='utf-8') as f:
        f.write(content)

fix_file('../papers/sections/03_methodology.tex')
fix_file('../papers/sections/05_experiments.tex')
fix_file('../papers/sections/06_conclusion.tex')
