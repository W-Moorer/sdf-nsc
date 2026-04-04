lines = []
with open('papers/sections/03_methodology.tex', 'r', encoding='utf-8', errors='replace') as f:
    for line in f:
        if 'Topology-Aware Contact' in line and '(uirt' in line:
            pass
        elif 'WNbQb' in line:
            line = '\\\\subsection{Topology-Aware Contact Manifold Subsumption}\\n'
            lines.append(line)
        else:
            lines.append(line)
with open('papers/sections/03_methodology.tex', 'w', encoding='utf-8') as f:
    f.writelines(lines)
