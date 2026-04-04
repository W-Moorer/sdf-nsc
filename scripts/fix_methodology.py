import re
with open('papers/sections/03_methodology.tex', 'r', encoding='utf-8') as f:
    content = f.read()

# Replace the specific equation that has missing \& / missing $ errors:
old_eq = r'''\begin{equation}
\begin{cases}
\left\langle \mathbf{M} (\mathbf{v}^{l+1} - \mathbf{v}^l) - \Delta t \mathbf{f}_{ext} - \mathbf{J}_N^T \boldsymbol{\Lambda}_N, \delta \mathbf{v} - \mathbf{v}^{l+1} \right\rangle = 0 \\[8pt]
\left\langle \mathbf{J}_N \mathbf{v}^{l+1} + \frac{\alpha}{\Delta t} \mathbf{g}_{2nd}, \delta \boldsymbol{\Lambda}_N - \boldsymbol{\Lambda}_N \right\rangle \ge 0
\end{cases}
\end{equation}'''

new_eq = r'''\begin{equation}
\begin{cases}
\left\langle \mathbf{M} (\mathbf{v}^{l+1} - \mathbf{v}^l) - \Delta t \mathbf{f}_{ext} - \mathbf{J}_N^T \boldsymbol{\Lambda}_N, \delta \mathbf{v} - \mathbf{v}^{l+1} \right\rangle &= 0 \\[8pt]
\left\langle \mathbf{J}_N \mathbf{v}^{l+1} + \frac{\alpha}{\Delta t} \mathbf{g}_{2nd}, \delta \boldsymbol{\Lambda}_N - \boldsymbol{\Lambda}_N \right\rangle &\ge 0
\end{cases}
\end{equation}'''

# Cases actually doesn't use & before = directly like aligned. It uses & for the second column (condition).
# But if it complains about missing $, maybe the environment is broken or there is a blank line somewhere near it.
# Let's just use \begin{aligned} ... \end{aligned} since it works beautifully for two equations.
new_aligned = r'''\begin{equation}
\begin{aligned}
\left\langle \mathbf{M} (\mathbf{v}^{l+1} - \mathbf{v}^l) - \Delta t \mathbf{f}_{ext} - \mathbf{J}_N^T \boldsymbol{\Lambda}_N, \delta \mathbf{v} - \mathbf{v}^{l+1} \right\rangle &= 0 \\[8pt]
\left\langle \mathbf{J}_N \mathbf{v}^{l+1} + \frac{\alpha}{\Delta t} \mathbf{g}_{2nd}, \delta \boldsymbol{\Lambda}_N - \boldsymbol{\Lambda}_N \right\rangle &\ge 0
\end{aligned}
\end{equation}'''

# Replace normal cases block:
content = re.sub(r'\\begin\{equation\}\s*\\begin\{cases\}.*?\\end\{cases\}\s*\\end\{equation\}', new_aligned, content, flags=re.DOTALL)

with open('papers/sections/03_methodology.tex', 'w', encoding='utf-8') as f:
    f.write(content)

print("done")
