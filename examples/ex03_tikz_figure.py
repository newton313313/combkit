"""
範例 3：不用 manim，直接產生 LaTeX/TikZ 圖。
畫出五邊形的 5 種三角剖分，以及對應的 5 條 Dyck 路徑（同一欄是同一棵樹），依根／左／右上色。

    python ex03_tikz_figure.py        → 寫出 ex03_tikz_figure.tex
    xelatex ex03_tikz_figure.tex      → 得到 PDF（pdflatex 也可以）
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from combkit.counting import binary_trees, dyck_parts, tree_to_dyck, tree_to_triangulation
from combkit.tikz import dyck_prims, picture, place, tri_prims

lines = []
for k, t in enumerate(binary_trees(3)):
    lines += place(tri_prims(5, tree_to_triangulation(t, 0, 4), 0.6, tree=t), 1.7 * k, 0.0)
    lines += place(dyck_prims(tree_to_dyck(t), 0.2, dyck_parts(t)), 1.7 * k, -1.3)

doc = r"""\documentclass{article}
\usepackage{tikz}
\pagestyle{empty}
\colorlet{cA}{blue!75!black}\colorlet{cTop}{orange!90!black}
\colorlet{cLeft}{teal!80!black}\colorlet{cRight}{purple!80!black}\colorlet{cBad}{red!85!black}
\begin{document}
""" + picture(lines) + "\n\\end{document}\n"
out = Path(__file__).with_suffix(".tex")
out.write_text(doc, encoding="utf-8")
print(f"wrote {out.name}: {len(lines)} TikZ commands")
