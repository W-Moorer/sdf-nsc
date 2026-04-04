import re
with open('papers/sections/03_methodology.tex', encoding='utf-8') as f:
    lines = f.readlines()
    for i, l in enumerate(lines):
        print(f"{i+1:03d} {l.strip()}")
