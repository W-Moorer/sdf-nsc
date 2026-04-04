with open('papers/sections/03_methodology.tex', encoding='utf-8') as f:
    lines = f.read().splitlines()
    for i, l in enumerate(lines[105:135], 105):
        print(f"{i:03d} {l}")
