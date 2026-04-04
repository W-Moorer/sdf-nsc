with open('scan_out_utf8.txt', 'w', encoding='utf-8') as f:
    with open('papers/sections/03_methodology.tex', encoding='utf-8') as fin:
        for i, l in enumerate(fin.readlines()):
            if 'Topology-Aware' in l or '方案' in l or 'Scheme' in l:
                f.write(f"Line {i+1}: {l.strip()}\n")
