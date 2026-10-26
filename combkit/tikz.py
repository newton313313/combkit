r"""
tikz.py —— 用 Python 產生 TikZ 程式碼（純 Python，不需要 manim）

做法：先在原點附近建立「圖元」（填色多邊形、線段、點、方塊、圓、文字節點），
算出外框後整體平移到指定位置，再輸出成 TikZ 指令。
卡塔蘭物件（三角剖分、二元樹、Dyck／格子路徑、圓上配對、括號）都由 counting.py 的樹產生，
畫法與 manim 端的 combkit.draw.catalan_draw 相同。

顏色名稱 cA（一般）、cTop（根）、cLeft（左）、cRight（右）、cBad（不合法）要在 LaTeX 導言區定義，例如：
    \colorlet{cA}{blue!75!black}  \colorlet{cTop}{orange!90!black}
    \colorlet{cLeft}{teal!80!black}  \colorlet{cRight}{purple!80!black}  \colorlet{cBad}{red!85!black}

來源：catalan930/notes/make_catalan_figs.py（講義第 1 節的圖就是用它產生的）。
"""
from __future__ import annotations

import math

from combkit.counting import (  # noqa: E402
    ballot_more_than_k_times,
    bertrand_strict_count,
    binary_trees,
    catalan,
    dyck_parts,
    dyck_to_ballot,
    dyck_to_brackets,
    dyck_to_lattice,
    dyck_to_matching,
    first_bad_step,
    lattice_ok,
    lattice_points,
    lattice_words,
    matching_parts,
    reflect_tail,
    tree_size,
    tree_to_dyck,
    tree_to_paren_tokens,
    tree_to_triangulation,
    triangulations,
)

COLOR = {"root": "cTop", "left": "cLeft", "right": "cRight", None: "cA"}
FILL = {"root": "cTop!38", "left": "cLeft!22", "right": "cRight!20"}


# =============================================================================
# 幾何：先在原點附近建立圖元，再整體平移
# =============================================================================
def f(x: float) -> str:
    s = f"{x:.3f}".rstrip("0").rstrip(".")
    return "0" if s in ("-0", "") else s


def pt(x: float, y: float) -> str:
    return f"({f(x)},{f(y)})"


def bbox(prims):
    xs, ys = [], []
    for pr in prims:
        kind = pr[0]
        if kind in ("fill", "poly"):
            pts = pr[2]
        elif kind == "line":
            pts = [pr[2], pr[3]]
        elif kind in ("dot", "square", "circle"):
            (x, y), r = pr[2], pr[3]
            pts = [(x - r, y - r), (x + r, y + r)]
        elif kind == "node":
            (x, y), (w, h) = pr[2], pr[4]
            pts = [(x - w, y - h), (x + w, y + h)]
        else:
            raise ValueError(kind)
        xs += [p[0] for p in pts]
        ys += [p[1] for p in pts]
    return min(xs), min(ys), max(xs), max(ys)


def emit(prims, dx: float = 0.0, dy: float = 0.0) -> list:
    def P(p):
        return pt(p[0] + dx, p[1] + dy)

    out = []
    for pr in prims:
        kind = pr[0]
        if kind == "fill":
            out.append(rf"\fill[{pr[1]}] " + " -- ".join(P(p) for p in pr[2]) + " -- cycle;")
        elif kind == "poly":
            out.append(rf"\draw[{pr[1]}] " + " -- ".join(P(p) for p in pr[2]) + " -- cycle;")
        elif kind == "line":
            out.append(rf"\draw[{pr[1]}] {P(pr[2])} -- {P(pr[3])};")
        elif kind == "dot":
            out.append(rf"\fill[{pr[1]}] {P(pr[2])} circle ({f(pr[3])});")
        elif kind == "square":
            (x, y), s = pr[2], pr[3]
            out.append(rf"\fill[{pr[1]}] {pt(x - s + dx, y - s + dy)} rectangle {pt(x + s + dx, y + s + dy)};")
        elif kind == "circle":
            out.append(rf"\draw[{pr[1]}] {P(pr[2])} circle ({f(pr[3])});")
        elif kind == "node":
            out.append(rf"\node[{pr[1]}] at {P(pr[2])} {{{pr[3]}}};")
    return out


def place(prims, cx: float, cy: float) -> list:
    """把圖元的外框中心移到 (cx, cy)。"""
    x0, y0, x1, y1 = bbox(prims)
    return emit(prims, cx - (x0 + x1) / 2, cy - (y0 + y1) / 2)


def node(p, text, opts="", size=(0.12, 0.12)):
    return ("node", opts, p, text, size)


def colored_math(chars, parts, wrap=lambda c: c) -> str:
    """把連續屬於同一部分的字元包成一組上色（數學模式）。"""
    out, i = [], 0
    while i < len(chars):
        j = i
        while j < len(chars) and parts[j] == parts[i]:
            j += 1
        seg = "".join(wrap(c) for c in chars[i:j])
        out.append(rf"{{\color{{{COLOR[parts[i]]}}}{seg}}}")
        i = j
    return "".join(out)


# =============================================================================
# 物件（與 codes/utils/catalan_draw.py 的畫法相同）
# =============================================================================
def polygon_vertices(m: int, r: float):
    """底邊（第 0 與第 m-1 點）水平放在最下面；第 0 點在左下，順時針到右下的第 m-1 點。"""
    out = []
    for k in range(m):
        a = -math.pi / 2 - math.pi / m - 2 * math.pi * k / m
        out.append((r * math.cos(a), r * math.sin(a)))
    return out


def tri_prims(m: int, diags, r: float, tree=None, diag_w=1.1, dot_r=0.03):
    """凸 m 邊形的三角剖分；給 tree 時把根三角形、左右子多邊形上色。"""
    v = polygon_vertices(m, r)
    prims, colors = [], {}
    if tree is not None:
        k = tree_size(tree[0]) + 1
        left, right = list(range(0, k + 1)), list(range(k, m))
        prims.append(("fill", FILL["root"], [v[0], v[k], v[m - 1]]))
        if len(left) >= 3:
            prims.append(("fill", FILL["left"], [v[i] for i in left]))
        if len(right) >= 3:
            prims.append(("fill", FILL["right"], [v[i] for i in right]))
        for a, b in diags:
            if (a, b) in ((0, k), (k, m - 1)):
                colors[(a, b)] = "cTop"
            elif b <= k:
                colors[(a, b)] = "cLeft"
            else:
                colors[(a, b)] = "cRight"
    prims.append(("poly", "gray!70,line width=0.6pt,line join=round", v))
    for a, b in diags:
        prims.append(("line", f"{colors.get((a, b), 'cA')},line width={diag_w}pt,line cap=round", v[a], v[b]))
    for p in v:
        prims.append(("dot", "black", p, dot_r))
    return prims


def tree_prims(tree, gap: float, level: float, colored=True, labels=None,
               node_r=0.05, leaf_s=0.04, edge_w=0.8):
    """
    滿二元樹，樹根在 (0, 0)。葉子依中序等距排列，內部節點放在它所有葉子的正中間
    （與 catalan_draw.BinaryTreeMob 相同）。labels 給字母時，葉子畫成字母。
    """
    leaves, nodes = [], []

    def rec(t, depth, part):
        if t is None:
            leaves.append({"x": len(leaves), "d": depth, "part": part})
            return ("leaf", len(leaves) - 1)
        left = rec(t[0], depth + 1, part)
        i = len(nodes)
        nodes.append({"d": depth, "part": part, "ch": None})
        right = rec(t[1], depth + 1, part)
        nodes[i]["ch"] = (left, right)
        return ("node", i)

    if tree is None:
        leaves.append({"x": 0, "d": 0, "part": "root"})
    else:
        lref = rec(tree[0], 1, "left")
        nodes.append({"d": 0, "part": "root", "ch": None})
        root = len(nodes) - 1
        rref = rec(tree[1], 1, "right")
        nodes[root]["ch"] = (lref, rref)

    def span(ref):
        kind, i = ref
        if kind == "leaf":
            return leaves[i]["x"], leaves[i]["x"]
        return span(nodes[i]["ch"][0])[0], span(nodes[i]["ch"][1])[1]

    x0 = (len(leaves) - 1) / 2

    def pos(ref):
        kind, i = ref
        if kind == "leaf":
            return (leaves[i]["x"] - x0) * gap, -leaves[i]["d"] * level
        lo, hi = span(ref)
        return ((lo + hi) / 2 - x0) * gap, -nodes[i]["d"] * level

    def col(part, default):
        return COLOR[part] if colored else default

    prims = []
    for i, nd in enumerate(nodes):
        start = pos(("node", i))
        for ref in nd["ch"]:
            end = pos(ref)
            if ref[0] == "leaf" and labels:
                dx, dy = end[0] - start[0], end[1] - start[1]
                length = math.hypot(dx, dy)
                end = (end[0] - 0.17 * dx / length, end[1] - 0.17 * dy / length)
            part = "root" if nd["part"] == "root" else nd["part"]
            prims.append(("line", f"{col(part, 'gray!70')},line width={edge_w}pt,line cap=round", start, end))
    for i, nd in enumerate(nodes):
        prims.append(("dot", col(nd["part"], "cA"), pos(("node", i)), node_r))
    for i, lf in enumerate(leaves):
        p = pos(("leaf", i))
        if labels:
            prims.append(node((p[0], p[1] - 0.03), f"${labels[i]}$",
                              f"font=\\small,inner sep=0pt,text={col(lf['part'], 'black')}", (0.1, 0.13)))
        else:
            prims.append(("square", col(lf["part"], "gray!60"), p, leaf_s))
    return prims


def dyck_prims(word: str, unit: float, parts=None, w=1.1, base=True):
    pts = [(0.0, 0.0)]
    for c in word:
        pts.append((pts[-1][0] + unit, pts[-1][1] + (unit if c == "U" else -unit)))
    prims = []
    if base:
        prims.append(("line", "gray!55,line width=0.4pt", (0.0, 0.0), (unit * max(len(word), 1), 0.0)))
    for j in range(len(word)):
        color = COLOR[parts[j]] if parts else "cA"
        prims.append(("line", f"{color},line width={w}pt,line cap=round", pts[j], pts[j + 1]))
    if not word:
        prims.append(("dot", "cA", (0.0, 0.0), 0.05))
    return prims


def lattice_prims(word: str, n: int, unit: float, color="cA", colors=None, w=1.1, ny=None):
    ny = n if ny is None else ny
    prims = []
    for i in range(n + 1):
        prims.append(("line", "gray!30,line width=0.3pt", (i * unit, 0.0), (i * unit, ny * unit)))
    for j in range(ny + 1):
        prims.append(("line", "gray!30,line width=0.3pt", (0.0, j * unit), (n * unit, j * unit)))
    prims.append(("line", "gray!60,dashed,line width=0.4pt", (0.0, 0.0), (n * unit, n * unit)))
    pts = [(x * unit, y * unit) for x, y in lattice_points(word)]
    for j in range(len(word)):
        c = colors[j] if colors else color
        prims.append(("line", f"{c},line width={w}pt,line cap=round", pts[j], pts[j + 1]))
    return prims


def circle_point(i: int, npts: int, r: float):
    a = math.radians(90 - 360 * i / npts)
    return r * math.cos(a), r * math.sin(a)


def matching_prims(pairs, npts: int, r: float, parts=None, root=False, w=1.1, labels=False):
    """圓周 npts 點（第 0 點在正上方、順時針）的配對；root=True 時第 0 點畫成橘色。"""
    prims = [("circle", "gray!65,line width=0.5pt", (0.0, 0.0), r)]
    for j, (a, b) in enumerate(pairs):
        color = COLOR[parts[j]] if parts else "cA"
        prims.append(("line", f"{color},line width={w}pt,line cap=round",
                      circle_point(a, npts, r), circle_point(b, npts, r)))
    for i in range(npts):
        prims.append(("dot", "black", circle_point(i, npts, r), 0.035))
    if root:
        prims.append(("dot", "cTop", circle_point(0, npts, r), 0.065))
    if labels:
        for i in range(npts):
            prims.append(node(circle_point(i, npts, r + 0.24), f"{i + 1}", "font=\\scriptsize"))
    return prims


def paren_math(tree, factors="abcd", colored=True) -> str:
    tokens, parts = tree_to_paren_tokens(tree, factors)
    if not colored:
        return "".join(r"\cdot " if t == "·" else t for t in tokens)
    out, i = [], 0
    while i < len(tokens):
        if parts[i] == "root":                      # 最後一次乘法：保留二元運算的間距
            out.append(r"\mathbin{\color{cTop}\cdot}")
            i += 1
            continue
        j = i
        while j < len(tokens) and parts[j] == parts[i]:
            j += 1
        seg = "".join(r"\cdot " if t == "·" else t for t in tokens[i:j])
        out.append(rf"{{\color{{{COLOR[parts[i]]}}}{seg}}}")
        i = j
    return "".join(out)


def word_math(word: str, parts, kind: str) -> str:
    """Dyck 字轉成括號字串或開票順序，逐字依根／左／右上色。"""
    text = dyck_to_brackets(word) if kind == "brackets" else dyck_to_ballot(word)
    wrap = (lambda c: c) if kind == "brackets" else (lambda c: rf"\mathrm{{{c}}}")
    return colored_math(text, parts, wrap)


def picture(lines, opts="") -> str:
    head = r"\begin{tikzpicture}" + (f"[{opts}]" if opts else "")
    return head + "\n" + "\n".join(lines) + "\n" + r"\end{tikzpicture}"
