"""
範例 9：不用 manim。高斯數「所有配對」，隨機矩陣數「不交叉配對」。

先算：6 個點的 15 種配對裡有 5 種不交叉（(2k-1)!! vs C_k）；GOE 的精確動差 (1/N) E Tr J^4 = 2 + 5/N + 5/N^2；
k=3 的 5 條樹形閉路徑 ↔ Dyck 路徑 ↔ 平面樹 ↔ 不交叉配對。全部 assert 之後再畫成 TikZ。

    python ex09_wick_moments.py        → 寫出 ex09_wick_moments.tex
    xelatex ex09_wick_moments.tex      → 得到 PDF
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from combkit.counting import (binary_to_plane, binary_trees, catalan, dyck_to_matching, matching_is_noncrossing,
                              perfect_matchings, tree_to_dyck)
from combkit.rmt_moments import (double_factorial_odd, goe_moment_poly, gue_moment_poly, tree_walks, walk_pairing,
                                 walk_to_dyck)
from combkit.tikz import dyck_prims, matching_prims, picture, place
from combkit.tikz_trees import plane_tree_prims

# ---- 先算、assert
allm = perfect_matchings(6)
assert len(allm) == double_factorial_odd(3) == 15
assert sum(matching_is_noncrossing(m) for m in allm) == catalan(3) == 5
assert goe_moment_poly(2) == {1: 2, 0: 5, -1: 5} and gue_moment_poly(2) == {1: 2, -1: 1}
walks = tree_walks(3)
for t, w in zip(binary_trees(3), walks):
    d = tree_to_dyck(t)
    assert walk_to_dyck(w) == d and walk_pairing(w) == sorted(dyck_to_matching(d))

# ---- 再畫：上排 15 種配對（不交叉的加底色），下面三排是 k=3 的 Dyck 路徑、平面樹、配對
lines = []
for idx, m in enumerate(allm):
    cx = idx * 1.25
    good = matching_is_noncrossing(m)
    prims = matching_prims(list(m), 6, 0.45)
    if not good:
        prims = [(p[0], p[1].replace("cA", "gray!60"), *p[2:]) if p[0] == "line" else p for p in prims]
    else:
        lines.append(rf"\fill[cA!8] ({cx:.2f},0) circle (0.58);")
    lines += place(prims, cx, 0)
for j, (t, w) in enumerate(zip(binary_trees(3), walks)):
    cx = 2.0 + j * 3.2
    d = tree_to_dyck(t)
    lines += place(dyck_prims(d, 0.28), cx, -1.6)
    lines += place(plane_tree_prims(binary_to_plane(t)), cx, -2.9)
    lines += place(matching_prims(dyck_to_matching(d), 6, 0.45), cx, -4.3)
    walk_tex = r"\,".join(map(str, w))
    lines.append(rf"\node[font=\footnotesize] at ({cx:.2f},-5.1) {{${walk_tex}$}};")

doc = r"""\documentclass{article}
\usepackage[margin=1cm,landscape]{geometry}
\usepackage{tikz}
\pagestyle{empty}
\colorlet{cA}{blue!75!black}\colorlet{cTop}{orange!90!black}
\colorlet{cLeft}{teal!80!black}\colorlet{cRight}{purple!80!black}\colorlet{cBad}{red!85!black}
\begin{document}
""" + picture(lines) + "\n\\end{document}\n"
out = Path(__file__).with_suffix(".tex")
out.write_text(doc, encoding="utf-8")
print(f"wrote {out.name}: {len(lines)} TikZ commands")
print("GOE (1/N) E Tr J^4 = 2 + 5/N + 5/N^2;  GUE = 2 + 1/N^2")
