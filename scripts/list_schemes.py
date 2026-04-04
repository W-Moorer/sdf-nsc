import os

output = []
for root, _, files in os.walk('papers'):
    for f in files:
        if f.endswith('.tex'):
            path = os.path.join(root, f)
            with open(path, encoding='utf-8') as fp:
                lines = fp.readlines()
                for i, l in enumerate(lines):
                    if '方案A' in l or '方案B' in l or '方案 A' in l or '方案 B' in l:
                        output.append(f"{path}:{i+1}:{l.strip()}")

with open('find_res.txt', 'w', encoding='utf-8') as f:
    f.write('\n'.join(output))

with open('find_res.txt', encoding='utf-8') as f:
    print(f.read())
