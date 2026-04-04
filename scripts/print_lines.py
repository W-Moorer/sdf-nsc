with open('meth_lines.txt', 'w', encoding='utf-8') as f:
    with open('papers/sections/03_methodology.tex', encoding='utf-8') as fin:
        lines = fin.readlines()
        for i in range(215, min(240, len(lines))):
            f.write(f"L{i+1}: {lines[i]}")
