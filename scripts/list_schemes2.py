import os

output = []
for root, _, files in os.walk('papers/sections'):
    for f in files:
        if f.endswith('.tex'):
            path = os.path.join(root, f)
            with open(path, encoding='utf-8') as fp:
                lines = fp.readlines()
                for i, l in enumerate(lines):
                    if '方案' in l or 'Scheme' in l or 'scheme' in l:
                        output.append(f"{path}:{i+1}:{l.strip()}")

with open('find_res2.txt', 'w', encoding='utf-8') as f:
    f.write('\n'.join(output))

with open('find_res2.txt', encoding='utf-8') as f:
    print("RES:")
    print(f.read())
