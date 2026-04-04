import os
import glob
import re

tex_files = glob.glob('papers/sections/*.tex') + ['papers/main.tex']

for filepath in tex_files:
    if os.path.exists(filepath):
        with open(filepath, 'r', encoding='utf-8') as f:
            content = f.read()

        # Fix missing $ before ^0
        content = re.sub(r'(\s|^)\^0\$', r' $C^0$', content)
        content = re.sub(r'(?<!\$)\^0\$', r'$C^0$', content)
        content = content.replace('^0连续性', '$C^0$连续性')
        content = content.replace('^0边界', '$C^0$边界')
        content = content.replace('^0 边界', '$C^0$ 边界')
        content = content.replace('^0 尖点', '$C^0$ 尖点')
        content = content.replace('^0奇异', '$C^0$奇异')
        content = content.replace('^0 奇异', '$C^0$ 奇异')
        content = content.replace('^0 级', '$C^0$ 级')
        content = content.replace('^0级', '$C^0$级')

        # Empty lines in math environments will crash latex
        def remove_blank_lines(match):
            eq_content = match.group(0)
            return re.sub(r'\n\s*\n+', '\n', eq_content)
        
        content = re.sub(r'\\begin\{equation\}(.*?)\\end\{equation\}', remove_blank_lines, content, flags=re.DOTALL)
        content = re.sub(r'\\begin\{align\}(.*?)\\end\{align\}', remove_blank_lines, content, flags=re.DOTALL)
        content = re.sub(r'\\begin\{cases\}(.*?)\\end\{cases\}', remove_blank_lines, content, flags=re.DOTALL)
        
        # In cases, usually you need a \& for alignment if it's not cases format but aligned.
        # But \begin{cases} handles it ok (e.g. `x & y \\ z & w`). But what if the user wrote two lines without `&` or just commas?
        # That's fine for simple equations, it will just typeset it sequentially. Or maybe it complains about not being able to find `. 

        # Fix missing `$` symbols or malformed math inside text
        content = content.replace(r'\ \times \text{Voxel\_Size}$', r'$n \times \text{Voxel\_Size}$')
        content = content.replace(r'\ \times 3$$', r'$3 \times 3$')
        content = content.replace(r'\ \times 3$', r'$3 \times 3$')
        content = content.replace(' =(V, E)$', '$G=(V, E)$')
        content = content.replace('=(V, E)$', '$G=(V, E)$')
        
        content = content.replace('第 $ 个接触', '第 $k$ 个接触')
        content = content.replace('区 $ 对应', '区 $k$ 对应')
        content = content.replace('形 $ 之上', '形 $k$ 之上')
        content = content.replace('长 \$ 直接', '长 $h$ 直接')

        with open(filepath, 'w', encoding='utf-8') as f:
            f.write(content)

print("done")
