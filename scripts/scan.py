with open('papers/sections/03_methodology.tex', encoding='utf-8') as f:
    lines = f.readlines()
for i, l in enumerate(lines):
    if 'Topology-Aware' in l or '方案' in l or 'Scheme' in l:
        print(f"Line {i+1}: {l.strip()}")
