with open('scan_a.txt', 'w', encoding='utf-8') as f:
    import os
    for root, _, files in os.walk('papers'):
        for file in files:
            if not file.endswith('.tex'): continue
            path = os.path.join(root, file)
            with open(path, encoding='utf-8') as fin:
                for i, l in enumerate(fin.readlines()):
                    if 'A' in l and ('方案' in l or 'Scheme' in l):
                        f.write(f"Line {i+1} in {path}: {l.strip()}\n")
                    if 'Scheme A' in l or '方案 A' in l or '方案A' in l:
                        f.write(f"EXACT A!! Line {i+1} in {path}: {l.strip()}\n")
