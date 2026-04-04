import os
import glob
import re

tex_files = glob.glob('papers/sections/*.tex') + ['papers/main.tex']

for filepath in tex_files:
    if os.path.exists(filepath):
        with open(filepath, 'r', encoding='utf-8') as f:
            content = f.read()

        # Fix escaped dollar signs \$ -> $
        content = content.replace(r'\$', '$')
        
        # In 03_methodology.tex, fix the strange subsection title and corrupted references
        if '03_methodology' in filepath:
            # Fix garbage subsection title
            content = re.sub(r'\\subsection\{WNbQb~gR:SNAmb_M~v\(uirt\{eHh \(Topology-Aware Contact M  Manifold Subsumption\)\}', '', content)
            
            # Fix malformed inline math $ G=(V, E)$ instead of =(V, E)$ 
            content = content.replace(' =(V, E)$ ', ' =(V, E)$ ')
            content = content.replace('子流形 $ 之上', '子流形 $ 之上')
            content = content.replace('第 $ 个接触', '第 $ 个接触')
            content = content.replace('每个连通区 $ 对应', '每个连通区 $ 对应')

        # In 04_implementation.tex, there's \^0$ and \ \times 3$
        if '04_implementation.tex' in filepath:
            content = content.replace(r'\^0', r'^0')
            content = content.replace(r' \times 3', r' \times 3$')
            content = content.replace(r'\\$', r'$')
        
        # In 00_abstract, there may be ^0$ again
        if '00_abstract' in filepath:
            content = content.replace('^0$', '^0$')

        with open(filepath, 'w', encoding='utf-8') as f:
            f.write(content)

print('LaTeX math fixed.')
