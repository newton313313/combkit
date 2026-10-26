"""
catalan_draw.py —— 卡塔蘭基礎模組（場景 00a–00m）的圖形

所有物件都由 counting.py 的「樹」產生（None = 空，(L, R) = 一個內部節點），
所以同一棵樹畫出來的三角剖分、括號、二元樹、格子路徑、Dyck 路徑、配對彼此一一對應。

上色規則（與 Chu 的部分一致，見 styles.py）：
    一般的卡塔蘭物件：藍色 CATALAN_COLOR（＝ 2 點群組的顏色，也就是之後的 x_1）
    根（固定的那一個）：橘色 ROOT_COLOR
    左邊的子問題：青色 LEFT_COLOR（＝ 區域 R_1 的顏色）
    右邊的子問題：粉紅 RIGHT_COLOR（＝ 區域 R_2 的顏色）
    不合法（越界、交叉）：紅色 CROSSING_COLOR
"""
from manimlib import *

from combkit.counting import (
    dyck_parts,
    dyck_to_lattice,
    dyck_to_matching,
    matching_parts,
    tree_size,
    tree_to_dyck,
    tree_to_paren_tokens,
    tree_to_triangulation,
)
from combkit.draw.polygon_groups import make_configuration
from combkit.draw.styles import CROSSING_COLOR, REGION_COLORS, ROOT_COLOR, TYPE1_COLOR

CATALAN_COLOR = TYPE1_COLOR
LEFT_COLOR = REGION_COLORS[0]
RIGHT_COLOR = REGION_COLORS[1]
EDGE_COLOR = GREY_B
LEAF_COLOR = GREY_B
PART_COLORS = {"root": ROOT_COLOR, "left": LEFT_COLOR, "right": RIGHT_COLOR}


def part_color(part, default=CATALAN_COLOR):
    return PART_COLORS.get(part, default)


# =============================================================================
# 凸多邊形與三角剖分
# =============================================================================
def polygon_vertices(m: int, radius: float = 1.0, center=ORIGIN) -> np.ndarray:
    """
    正 m 邊形的頂點。底邊（第 0 與第 m-1 點）水平放在最下面；
    第 0 點在左下，沿左邊往上、經過上方，順時針到右下的第 m-1 點。
    """
    center = np.array(center, dtype=float)
    out = []
    for k in range(m):
        a = -PI / 2 - PI / m - TAU * k / m
        out.append(center + radius * np.array([np.cos(a), np.sin(a), 0.0]))
    return np.array(out)


class TriangulationMob(VGroup):
    """
    凸 m 邊形 + 對角線 + 頂點。子物件順序：[fills, outline, diagonals, dots]。
    屬性：.outline、.diagonals（每條有 .pair）、.dots、.fills、.m
    方法：.vertex(k)（物件移動後仍正確）、.region(indices, color, opacity)
    """

    def __init__(self, m: int, diagonals=(), radius: float = 1.0, center=ORIGIN,
                 edge_width: float = 3, diag_width: float = 4.5, dot_radius: float = 0.05,
                 diag_color=CATALAN_COLOR, edge_color=EDGE_COLOR, show_dots: bool = True):
        verts = polygon_vertices(m, radius, center)
        outline = Polygon(*verts, stroke_color=edge_color, stroke_width=edge_width)
        outline.set_fill(opacity=0)
        diag_mobs = VGroup()
        for a, b in diagonals:
            line = Line(verts[a], verts[b], stroke_color=diag_color, stroke_width=diag_width)
            line.pair = (a, b)
            diag_mobs.add(line)
        dots = VGroup(*[Dot(p, radius=dot_radius, fill_color=WHITE) for p in verts])
        if not show_dots:
            dots.set_opacity(0)
        dots.set_z_index(2)
        fills = VGroup()
        super().__init__(fills, outline, diag_mobs, dots)
        self.fills = fills
        self.outline = outline
        self.diagonals = diag_mobs
        self.dots = dots
        self.m = m
        self.radius = radius

    def vertex(self, k: int) -> np.ndarray:
        return self.dots[k].get_center()

    def region(self, indices, color, opacity: float = 0.35) -> Polygon:
        poly = Polygon(*[self.vertex(i) for i in indices], stroke_width=0)
        poly.set_fill(color, opacity=opacity)
        return poly

    def edge(self, a: int, b: int, color=ROOT_COLOR, width: float = 7) -> Line:
        return Line(self.vertex(a), self.vertex(b), stroke_color=color, stroke_width=width)


def triangulation_split(tree, m: int):
    """回傳 (k, left_indices, right_indices)：含底邊的三角形是 (0, k, m-1)。"""
    k = tree_size(tree[0]) + 1
    return k, list(range(0, k + 1)), list(range(k, m))


def triangulation_item(tree, radius: float = 1.0, center=ORIGIN, colored: bool = False,
                       **kwargs) -> TriangulationMob:
    """由樹畫出凸 (n+2) 邊形的三角剖分；colored=True 時把根三角形、左右子多邊形上色。"""
    m = tree_size(tree) + 2
    tm = TriangulationMob(m, tree_to_triangulation(tree, 0, m - 1), radius=radius, center=center, **kwargs)
    if colored and tree is not None:
        k, left, right = triangulation_split(tree, m)
        tm.fills.add(tm.region([0, k, m - 1], ROOT_COLOR, 0.55))
        if len(left) >= 3:
            tm.fills.add(tm.region(left, LEFT_COLOR, 0.4))
        if len(right) >= 3:
            tm.fills.add(tm.region(right, RIGHT_COLOR, 0.4))
        for d in tm.diagonals:
            a, b = d.pair
            if (a, b) in ((0, k), (k, m - 1)):
                d.set_stroke(ROOT_COLOR)
            elif b <= k:
                d.set_stroke(LEFT_COLOR)
            else:
                d.set_stroke(RIGHT_COLOR)
    return tm


# =============================================================================
# 滿二元樹
# =============================================================================
class BinaryTreeMob(VGroup):
    """
    滿二元樹：內部節點（實心點）、葉子（小方塊或字母）、邊。
    葉子依中序由左到右等距排列；內部節點放在它所有葉子的正中間；高度依深度。
    屬性（都依中序排列）：
        .nodes、.node_parts       內部節點與它屬於 root/left/right
        .leaves、.leaf_parts      葉子
        .edges、.edge_parts       邊（根連到兩個小孩的邊算 root）
        .node_pos、.leaf_pos      位置（建立當下）
    """

    def __init__(self, tree, root_pos=ORIGIN, leaf_gap: float = 0.6, level_height: float = 0.6,
                 node_radius: float = 0.08, leaf_labels=None, leaf_font: float = 36,
                 leaf_size: float = 0.1, edge_width: float = 3, node_color=CATALAN_COLOR,
                 edge_color=EDGE_COLOR, leaf_color=LEAF_COLOR, colored: bool = False,
                 root_color=None):
        root_pos = np.array(root_pos, dtype=float)
        leaves, nodes, edges = [], [], []

        def rec(t, depth, part):
            """回傳 ('leaf'|'node', index)。"""
            if t is None:
                leaves.append({"x": len(leaves), "depth": depth, "part": part})
                return ("leaf", len(leaves) - 1)
            left = rec(t[0], depth + 1, part)
            idx = len(nodes)
            nodes.append({"depth": depth, "part": part, "children": None})
            right = rec(t[1], depth + 1, part)
            nodes[idx]["children"] = (left, right)
            return ("node", idx)

        if tree is None:
            leaves.append({"x": 0, "depth": 0, "part": "root"})
        else:
            # 根自己是 root，左右子樹分別是 left / right
            l_ref = rec(tree[0], 1, "left")
            nodes.append({"depth": 0, "part": "root", "children": None})
            root_idx = len(nodes) - 1
            r_ref = rec(tree[1], 1, "right")
            nodes[root_idx]["children"] = (l_ref, r_ref)

        # 每個內部節點的 x = 它子樹葉子的中點
        def span(ref):
            kind, i = ref
            if kind == "leaf":
                return leaves[i]["x"], leaves[i]["x"]
            lo, _ = span(nodes[i]["children"][0])
            _, hi = span(nodes[i]["children"][1])
            return lo, hi

        n_leaves = len(leaves)
        x0 = (n_leaves - 1) / 2

        def pos_of(ref):
            kind, i = ref
            if kind == "leaf":
                x, d = leaves[i]["x"], leaves[i]["depth"]
            else:
                lo, hi = span(ref)
                x, d = (lo + hi) / 2, nodes[i]["depth"]
            return root_pos + np.array([(x - x0) * leaf_gap, -d * level_height, 0.0])

        # nodes 已經依中序建立（左子樹、自己、右子樹），正好對應括號裡由左到右的乘號。
        # 注意：中序不一定等於 x 座標的順序（偏左的樹裡，左子樹的節點可能畫在根的右邊）。
        node_pos = [pos_of(("node", i)) for i in range(len(nodes))]
        leaf_pos = [pos_of(("leaf", i)) for i in range(n_leaves)]

        # 字母葉子：邊停在字母上方一點，不穿過字母
        label_buff = 0.0055 * leaf_font if leaf_labels is not None else 0.0
        edge_mobs, edge_parts, edge_kinds = VGroup(), [], []
        for i, nd in enumerate(nodes):
            for ref in nd["children"]:
                start = node_pos[i]
                if ref[0] == "node":
                    end = node_pos[ref[1]]
                else:
                    end = leaf_pos[ref[1]]
                    end = end - label_buff * normalize(end - start)
                part = "root" if nd["part"] == "root" else nd["part"]
                line = Line(start, end, stroke_color=edge_color, stroke_width=edge_width)
                edge_mobs.add(line)
                edge_parts.append(part)
                edge_kinds.append(ref[0])

        node_mobs = VGroup(*[Dot(p, radius=node_radius, fill_color=node_color) for p in node_pos])
        if root_color is not None and nodes:
            for nd, dot in zip(nodes, node_mobs):
                if nd["part"] == "root":
                    dot.set_fill(root_color)
        leaf_mobs = VGroup()
        for i, p in enumerate(leaf_pos):
            if leaf_labels is not None:
                mob = Tex(leaf_labels[i], font_size=leaf_font, fill_color=WHITE)
                mob.move_to(p + 0.004 * leaf_font * DOWN)
            else:
                mob = Square(side_length=2 * leaf_size, stroke_width=0)
                mob.set_fill(leaf_color, opacity=1)
                mob.move_to(p)
            leaf_mobs.add(mob)

        super().__init__(edge_mobs, node_mobs, leaf_mobs)
        self.edges = edge_mobs
        self.edge_parts = edge_parts
        self.edge_kinds = edge_kinds          # 每條邊接到 'node' 還是 'leaf'
        self.nodes = node_mobs
        self.node_parts = [nd["part"] for nd in nodes]
        self.leaves = leaf_mobs
        self.leaf_parts = [lf["part"] for lf in leaves]
        self.node_pos = node_pos
        self.leaf_pos = leaf_pos

        if colored:
            self.color_parts()

    def color_parts(self):
        for mob, part in zip(self.nodes, self.node_parts):
            mob.set_fill(part_color(part))
        for mob, part in zip(self.edges, self.edge_parts):
            mob.set_stroke(part_color(part))
        for mob, part in zip(self.leaves, self.leaf_parts):
            mob.set_fill(part_color(part))
        return self


def tree_item(tree, scale: float = 1.0, colored: bool = False, **kwargs) -> BinaryTreeMob:
    kw = dict(leaf_gap=0.26 * scale, level_height=0.22 * scale, node_radius=0.05 * scale,
              leaf_size=0.035 * scale, edge_width=2.2)
    kw.update(kwargs)
    return BinaryTreeMob(tree, colored=colored, **kw)


# =============================================================================
# 完全括號化（Tex，每個符號一個 glyph）
# =============================================================================
def paren_mob(tree, factors="abcdefghij", font_size: float = 40, colored: bool = False,
              color=WHITE) -> Tex:
    """
    回傳 Tex，附加屬性：
        .glyphs  依閱讀順序排列的符號（與 .tokens 一一對應）
        .tokens  ['(', 'a', '·', ...]
        .parts   每個符號屬於 'root'（最後一次乘法）、'left'、'right'
    """
    tokens, parts = tree_to_paren_tokens(tree, factors)
    tex = Tex(" ".join(R"\cdot" if t == "·" else t for t in tokens), font_size=font_size)
    tex.set_fill(color)
    glyphs = sorted(tex.submobjects, key=lambda g: g.get_center()[0])
    assert len(glyphs) == len(tokens), (tokens, len(glyphs))
    tex.glyphs = glyphs
    tex.tokens = tokens
    tex.parts = parts
    if colored:
        for g, p in zip(glyphs, parts):
            g.set_fill(part_color(p))
    return tex


def catalan_product_tex(k: int, j: int, tail: str = "", font_size: float = 34) -> Tex:
    """
    「C_k C_j」＋ tail（例如 "=2"），C_k 用左邊的青色、C_j 用右邊的粉紅。
    依 glyph 位置上色（k = j 時 t2c 分不出左右）。k、j 只能是個位數。
    """
    tex = Tex(f"C_{{{k}}}C_{{{j}}}{tail}", font_size=font_size)
    glyphs = sorted(tex.submobjects, key=lambda g: g.get_center()[0])
    VGroup(*glyphs[0:2]).set_fill(LEFT_COLOR)
    VGroup(*glyphs[2:4]).set_fill(RIGHT_COLOR)
    return tex


def paren_item(tree, scale: float = 1.0, colored: bool = False, factors="abcd") -> Tex:
    return paren_mob(tree, factors, font_size=26 * scale, colored=colored)


def paren_depth_stages(tokens) -> list:
    """
    把括號依「由內而外」分批：回傳 [[(開括號位置, 閉括號位置), ...], ...]，
    第一批是最內層（巢狀最深）的括號。
    """
    stack, pairs = [], []
    for i, tok in enumerate(tokens):
        if tok == "(":
            stack.append(i)
        elif tok == ")":
            j = stack.pop()
            pairs.append((j, i, len(stack)))     # len(stack) = 外面還有幾層
    if not pairs:
        return []
    depths = sorted({d for _, _, d in pairs}, reverse=True)
    return [[(a, b) for a, b, d in pairs if d == depth] for depth in depths]


def non_paren_glyphs(mob: Tex) -> list:
    """括號以外的符號（因子與乘號），依閱讀順序。"""
    return [g for g, t in zip(mob.glyphs, mob.tokens) if t not in "()"]


def paren_build_animations(plain: Tex, target: Tex, run_time: float = 0.7):
    """
    從沒有括號的「a·b·c·d」長出某一種括號法：
        第 1 步：因子與乘號從 plain 飛到 target 的位置
        之後：括號由內而外一層一層出現
    回傳動畫串列的串列，依序 self.play(*step)。
    """
    src = sorted(plain.submobjects, key=lambda g: g.get_center()[0])
    dst = non_paren_glyphs(target)
    assert len(src) == len(dst)
    steps = [[TransformFromCopy(a, b) for a, b in zip(src, dst)]]
    for stage in paren_depth_stages(target.tokens):
        anims = []
        for a, b in stage:
            anims.append(GrowFromCenter(target.glyphs[a]))
            anims.append(GrowFromCenter(target.glyphs[b]))
        steps.append(anims)
    return steps


# =============================================================================
# 路徑：Dyck 路徑（↗ ↘）與格子路徑（→ ↑）
# =============================================================================
DYCK_STEPS = {"U": np.array([1.0, 1.0, 0.0]), "D": np.array([1.0, -1.0, 0.0])}
LATTICE_STEPS = {"R": np.array([1.0, 0.0, 0.0]), "U": np.array([0.0, 1.0, 0.0])}


class StepPathMob(VGroup):
    """
    一條由「步」組成的路徑。子物件：[steps, dots]。
    屬性：.steps（每步一條線段）、.dots（每個點一個圓點，含起點）、.word
    方法：.vertex(j)
    """

    def __init__(self, word, step_vectors: dict, unit: float = 0.5, origin=ORIGIN,
                 color=CATALAN_COLOR, colors=None, stroke_width: float = 5,
                 dot_radius: float = 0.05, show_dots: bool = True):
        origin = np.array(origin, dtype=float)
        pts = [origin]
        for c in word:
            pts.append(pts[-1] + unit * step_vectors[c])
        steps = VGroup()
        for j, c in enumerate(word):
            col = colors[j] if colors is not None else color
            steps.add(Line(pts[j], pts[j + 1], stroke_color=col, stroke_width=stroke_width))
        dots = VGroup(*[Dot(p, radius=dot_radius, fill_color=WHITE) for p in pts])
        if not show_dots:
            dots.set_opacity(0)
        dots.set_z_index(2)
        super().__init__(steps, dots)
        self.steps = steps
        self.dots = dots
        self.word = word
        self.unit = unit

    def vertex(self, j: int) -> np.ndarray:
        return self.dots[j].get_center()


def dyck_path_mob(word: str, unit: float = 0.5, origin=ORIGIN, **kwargs) -> StepPathMob:
    return StepPathMob(word, DYCK_STEPS, unit=unit, origin=origin, **kwargs)


def lattice_path_mob(word: str, unit: float = 0.5, origin=ORIGIN, **kwargs) -> StepPathMob:
    return StepPathMob(word, LATTICE_STEPS, unit=unit, origin=origin, **kwargs)


def lattice_grid(nx: int, ny: int, unit: float = 0.5, origin=ORIGIN, color=GREY_D,
                 stroke_width: float = 1.5) -> VGroup:
    """格子線：x = 0..nx、y = 0..ny。"""
    origin = np.array(origin, dtype=float)
    lines = VGroup()
    for i in range(nx + 1):
        lines.add(Line(origin + unit * np.array([i, 0, 0]), origin + unit * np.array([i, ny, 0]),
                       stroke_color=color, stroke_width=stroke_width))
    for j in range(ny + 1):
        lines.add(Line(origin + unit * np.array([0, j, 0]), origin + unit * np.array([nx, j, 0]),
                       stroke_color=color, stroke_width=stroke_width))
    return lines


def lattice_diagonal(n: int, unit: float = 0.5, origin=ORIGIN, color=GREY_A,
                     stroke_width: float = 2, offset: float = 0.0) -> DashedLine:
    """對角線 y = x + offset（從 x=0 畫到 y=n 或 x=n）。"""
    origin = np.array(origin, dtype=float)
    if offset >= 0:
        a, b = np.array([0, offset, 0]), np.array([n - offset, n, 0])
    else:
        a, b = np.array([-offset, 0, 0]), np.array([n, n + offset, 0])
    return DashedLine(origin + unit * a, origin + unit * b, dash_length=0.06,
                      stroke_color=color, stroke_width=stroke_width)


def lattice_thumbnail(word: str, n: int, unit: float = 0.25, origin=ORIGIN, valid=None,
                      stroke_width: float = 3.5) -> VGroup:
    """小圖：格子 + 對角線 + 路徑。valid=None 時路徑用一般藍色。屬性 .path、.grid、.diag。"""
    grid = lattice_grid(n, n, unit, origin, stroke_width=1)
    diag = lattice_diagonal(n, unit, origin, stroke_width=1.5)
    color = CATALAN_COLOR if valid in (None, True) else CROSSING_COLOR
    path = lattice_path_mob(word, unit, origin, color=color, stroke_width=stroke_width,
                            dot_radius=0.022, show_dots=False)
    group = VGroup(grid, diag, path)
    group.grid, group.diag, group.path = grid, diag, path
    return group


def dyck_item(tree, scale: float = 1.0, colored: bool = False, **kwargs) -> StepPathMob:
    word = tree_to_dyck(tree)
    colors = [part_color(p) for p in dyck_parts(tree)] if colored else None
    kw = dict(unit=0.15 * scale, stroke_width=3.2, dot_radius=0.025, show_dots=False)
    kw.update(kwargs)
    mob = dyck_path_mob(word, colors=colors, **kw)
    base = Line(mob.vertex(0), mob.vertex(len(word)), stroke_color=GREY_C, stroke_width=1.2)
    mob.add_to_back(base)
    mob.base = base
    return mob


def lattice_item(tree, scale: float = 1.0, colored: bool = False) -> VGroup:
    word = dyck_to_lattice(tree_to_dyck(tree))
    n = len(word) // 2
    thumb = lattice_thumbnail(word, n, unit=0.26 * scale, stroke_width=3.2)
    if colored:
        for step, part in zip(thumb.path.steps, dyck_parts(tree)):
            step.set_stroke(part_color(part))
    return thumb


# =============================================================================
# 圓周上的不交叉完美配對
# =============================================================================
def matching_item(tree, radius: float = 0.4, center=ORIGIN, colored: bool = False, **kwargs):
    """由樹畫出圓周 2n 點的不交叉完美配對；colored=True 時含第 0 點的弦是根。"""
    word = tree_to_dyck(tree)
    pairs = dyck_to_matching(word)
    colors = [part_color(p) for p in matching_parts(tree)] if colored else None
    return make_configuration(len(word), pairs, radius=radius, center=center,
                              root=0 if colored else None, colors=colors, **kwargs)


# =============================================================================
# 統一入口：同一棵樹的六種樣子
# =============================================================================
FAMILIES = ["triangulation", "parens", "tree", "lattice", "dyck", "matching"]
FAMILY_NAMES = {
    "triangulation": "凸五邊形的三角剖分",
    "parens": "4 個因子的完全括號化",
    "tree": "3 個內部節點的滿二元樹",
    "lattice": "不越過 y=x 的格子路徑",
    "dyck": "半長 3 的 Dyck 路徑",
    "matching": "圓周 6 點的不交叉配對",
}


def catalan_item(family: str, tree, scale: float = 1.0, colored: bool = False) -> Mobject:
    """畫出同一棵樹在某一種物件下的樣子，中心放在原點（大小約 0.9 × scale）。"""
    if family == "triangulation":
        mob = triangulation_item(tree, radius=0.42 * scale, colored=colored, edge_width=2,
                                 diag_width=3, dot_radius=0.03)
    elif family == "parens":
        mob = paren_item(tree, scale=scale, colored=colored)
    elif family == "tree":
        mob = tree_item(tree, scale=scale, colored=colored)
    elif family == "lattice":
        mob = lattice_item(tree, scale=scale, colored=colored)
    elif family == "dyck":
        mob = dyck_item(tree, scale=scale, colored=colored)
    elif family == "matching":
        mob = matching_item(tree, radius=0.4 * scale, colored=colored)
    else:
        raise ValueError(family)
    mob.move_to(ORIGIN)
    return mob
