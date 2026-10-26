r"""
tikz_trees.py —— 平面樹（有根、子樹有順序）的 TikZ 畫法（純 Python，不需要 manim）

平面樹用 combkit.counting 的表示法（binary_to_plane 的輸出）。版面：葉子依序等距排開，
父節點放在子節點的正中上方，根在最上面（橘色 cTop）。圖元格式與 combkit.tikz 相同，可以直接 place()。

    from combkit.counting import binary_trees, binary_to_plane
    from combkit.tikz import place, picture
    from combkit.tikz_trees import plane_tree_prims
    picture(place(plane_tree_prims(binary_to_plane(binary_trees(3)[0])), 0, 0))

來源：spglass 專案的 codes/make_solution_tikz.py（先修講義練習解答的圖 5、圖 7）。
"""
from __future__ import annotations

from combkit.counting import plane_tree_preorder
from combkit.tikz import node


def plane_tree_layout(pt, dx: float = 0.42, dy: float = 0.5, center: bool = True):
    """回傳 (preorder, pos)：preorder[v] = (父節點, 深度)，pos[v] = (x, y)；頂點按前序（第一次到訪）編號。
    center=True 時把根放在 x = 0。"""
    pre = plane_tree_preorder(pt)
    children = {i: [] for i in range(len(pre))}
    for i, (par, _) in enumerate(pre):
        if par is not None:
            children[par].append(i)
    xs: dict[int, float] = {}
    counter = [0]

    def lay(v):
        if not children[v]:
            xs[v] = counter[0]
            counter[0] += 1
        else:
            for c in children[v]:
                lay(c)
            xs[v] = sum(xs[c] for c in children[v]) / len(children[v])
    lay(0)
    x0 = xs[0] if center else 0.0
    pos = {v: ((xs[v] - x0) * dx, -pre[v][1] * dy) for v in sorted(xs)}
    return pre, pos


def plane_tree_prims(pt, dx: float = 0.42, dy: float = 0.5, labels=True, color: str = "cA",
                     label_style: str = "font=\\tiny,gray!70!black"):
    """平面樹的圖元。labels=True 標前序編號 0,1,2,…；也可以給一個字串或串列（例如 "ijk"）當頂點名稱；False 不標。"""
    pre, pos = plane_tree_layout(pt, dx, dy)
    prims = []
    for v, (par, _) in enumerate(pre):
        if par is not None:
            prims.append(("line", f"{color},line width=1.1pt", pos[par], pos[v]))
    for v, (x, y) in pos.items():
        prims.append(("dot", "cTop" if v == 0 else "black", (x, y), 0.06 if v == 0 else 0.045))
        if labels is True:
            prims.append(node((x + 0.17, y + 0.05), f"{v}", label_style))
        elif labels:
            prims.append(node((x + 0.22, y + 0.05), f"${labels[v]}$", label_style))
    return prims


if __name__ == "__main__":
    from combkit.counting import binary_to_plane, binary_trees, catalan, plane_tree_walk
    from combkit.tikz import picture, place

    for k in range(1, 6):
        trees = [binary_to_plane(t) for t in binary_trees(k)]
        assert len(trees) == catalan(k)
        for pt in trees:
            pre, pos = plane_tree_layout(pt)
            assert len(pos) == k + 1 and pos[0] == (0.0, 0.0)
            for v, (par, depth) in enumerate(pre):                    # 子節點在父節點下一層
                if par is not None:
                    assert abs(pos[v][1] - pos[par][1] + 0.5) < 1e-12
            assert sorted(set(plane_tree_walk(pt))) == list(range(k + 1))
    tex = picture(place(plane_tree_prims(binary_to_plane(binary_trees(3)[0]), labels="ijkl"), 0, 0))
    assert tex.startswith(r"\begin{tikzpicture}") and "$i$" in tex
    print("tikz_trees：平面樹版面與 TikZ 輸出全部驗證通過")
