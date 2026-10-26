"""
mit_draw.py —— MIT 18.212 L1–4 補充場景（M01–M14）的共用圖形

數學資料一律由 counting.py 第 6 節算（點、格子、列、行都是 0 起算），這裡只負責畫：
    ArcDiagramMob     直線上 2n 個點 + 上方的半圓弧（不交叉／不巢狀配對）
    arc_curve         兩點之間的一條上方半圓弧
    YoungDiagramMob   Young 圖（英式：第 0 列在上、靠左），可填數字、標出鉤
    PlaneTreeMob      平面樹（頂點依前序編號，與 counting.plane_tree_preorder 相同）
    box_row           一排方塊（排列、堆疊、佇列、數字列）

顏色沿用 catalan_draw：一般物件藍 CATALAN_COLOR、根橘、左青、右粉紅、不合法紅。
"""
from manimlib import *

from combkit.draw.catalan_draw import (
    CATALAN_COLOR,
    EDGE_COLOR,
    LEFT_COLOR,
    PART_COLORS,
    RIGHT_COLOR,
    ROOT_COLOR,
    part_color,
)
from combkit.counting import hook_cells, plane_tree_preorder, young_cells
from combkit.draw.styles import CROSSING_COLOR

NEST_COLOR = PURPLE_B          # 巢狀的弧（交叉用 CROSSING_COLOR 紅色）
HOOK_COLOR = YELLOW            # 鉤
CORNER_COLOR = GREEN           # 角落
CELL_STROKE = GREY_A


def _label_mob(value, font_size: float, color=WHITE) -> Mobject:
    """數字或短字串 -> Tex（中文請自己用 styles.zh 再傳 Mobject 進來）。"""
    if isinstance(value, Mobject):
        return value
    return Tex(str(value), font_size=font_size, fill_color=color)


# =============================================================================
# 弧圖
# =============================================================================
def arc_curve(p, q, color=CATALAN_COLOR, stroke_width: float = 4, height_ratio: float = 1.0,
              below: bool = False) -> Arc:
    """
    p、q（同一條水平線上）之間的半圓弧，從 p 畫到 q（ShowCreation 的方向）。
    預設在線的上方；height_ratio < 1 把弧壓扁；below=True 畫在下方。
    """
    p, q = np.array(p, dtype=float), np.array(q, dtype=float)
    mid = (p + q) / 2
    radius = np.linalg.norm(q - p) / 2
    start = angle_of_vector(p - mid)
    sweep = -PI if p[0] < q[0] else PI          # 從左往右順時針、從右往左逆時針 -> 都經過上方
    if below:
        sweep = -sweep
    arc = Arc(start_angle=start, angle=sweep, radius=radius, arc_center=mid,
              stroke_color=color, stroke_width=stroke_width)
    if height_ratio != 1.0:
        arc.stretch(height_ratio, 1, about_point=mid)
    return arc


class ArcDiagramMob(VGroup):
    """
    直線上 n_points 個點（0 起算，由左到右，間距 gap，整體中心 center）與上方的弧。
    子物件順序：[baseline, arcs, dots, labels]。
    屬性：.baseline、.arcs（第 k 條對應 pairs[k]，有 .pair）、.dots、.labels、
          .pairs（唯讀，每次由 .arcs 算出，.copy() 後仍正確）
    方法：.point(i)（物件移動後仍正確）、.make_arc(i, j, ...)（在目前位置做一條新弧，不加入）、
          .add_arc(i, j, ...)、.color_arcs(colors)
    labels：True -> 點下方寫 1..n（one_based=False 時寫 0..n-1）；list -> 自訂；None/False -> 不寫。
    colors：每條弧的顏色串列（或 None，全部 color）。
    """

    def __init__(self, n_points: int, pairs=(), gap: float = 0.6, center=ORIGIN, color=CATALAN_COLOR,
                 colors=None, stroke_width: float = 4, height_ratio: float = 1.0, labels=True,
                 one_based: bool = True, label_font: float = 26, label_color=GREY_A,
                 label_buff: float = 0.18, dot_radius: float = 0.06, dot_color=WHITE,
                 show_line: bool = True, line_color=GREY_C, line_width: float = 2, line_margin: float = 0.3):
        center = np.array(center, dtype=float)
        pts = [center + (i - (n_points - 1) / 2) * gap * RIGHT for i in range(n_points)]
        baseline = Line(pts[0] + line_margin * LEFT, pts[-1] + line_margin * RIGHT,
                        stroke_color=line_color, stroke_width=line_width)
        if not show_line:
            baseline.set_stroke(opacity=0)
        dots = VGroup(*[Dot(p, radius=dot_radius, fill_color=dot_color) for p in pts])
        dots.set_z_index(2)
        label_mobs = VGroup()
        if labels:
            texts = labels if isinstance(labels, (list, tuple)) else \
                [str(i + 1 if one_based else i) for i in range(n_points)]
            for t in texts:
                label_mobs.add(_label_mob(t, label_font, label_color))
            tallest = max((lab.get_height() for lab in label_mobs), default=0)
            for p, lab in zip(pts, label_mobs):      # 底部對齊（字母高低不同時仍整齊）
                lab.move_to(p, aligned_edge=DOWN)
                lab.shift((label_buff + tallest) * DOWN)
        self.dots = dots
        self.height_ratio = height_ratio
        self.arc_width = stroke_width
        self.arc_color = color
        arcs = VGroup()
        for k, (i, j) in enumerate(tuple(pr) for pr in pairs):
            col = colors[k] if colors is not None else color
            arcs.add(self.make_arc(i, j, color=col))
        super().__init__(baseline, arcs, dots, label_mobs)
        self.baseline = baseline
        self.arcs = arcs
        self.labels = label_mobs
        self.gap = gap

    @property
    def pairs(self) -> list:
        """目前每條弧的端點（由 .arcs 算出，所以 .copy() 之後仍正確）。"""
        return [arc.pair for arc in self.arcs]

    def point(self, i: int) -> np.ndarray:
        return self.dots[i].get_center()

    def make_arc(self, i: int, j: int, color=None, stroke_width=None, height_ratio=None) -> Arc:
        arc = arc_curve(self.point(i), self.point(j), color=color or self.arc_color,
                        stroke_width=stroke_width or self.arc_width,
                        height_ratio=self.height_ratio if height_ratio is None else height_ratio)
        arc.pair = (i, j)
        return arc

    def add_arc(self, i: int, j: int, color=None) -> Arc:
        arc = self.make_arc(i, j, color=color)
        self.arcs.add(arc)
        return arc

    def color_arcs(self, colors):
        for arc, col in zip(self.arcs, colors):
            arc.set_stroke(col)
        return self


# =============================================================================
# Young 圖
# =============================================================================
class YoungDiagramMob(VGroup):
    """
    Young 圖（英式）：格子 (i, j) = 第 i 列、第 j 行（0 起算），第 0 列在最上面、每列靠左。
    top_left 給定時，格子 (0, 0) 的左上角放在那裡；否則整張圖的中心放在 center。
    子物件順序：[cells, entries]。
    屬性：.lam、.cells（逐列由左到右，與 counting.young_cells 相同順序）、.entries（VGroup）、.cell_size
    方法：
        .cell(i, j)                 那一格的 Square
        .cell_center(i, j)          格子中心（物件移動後仍正確）
        .cell_group(cells)          VGroup(那些格子)
        .fill_cells(cells, color, opacity)          直接改格子的填色
        .set_entries(table, ...)    依列填入數字或字串（table[i][j]；None 跳過），回傳新的 entries
        .entry(i, j)                目前那一格的文字（沒有則 None）
        .hook_highlight(i, j, ...)  回傳鉤的半透明覆蓋方塊 VGroup（第 0 個是格子本身），不加入圖中
    """

    def __init__(self, lam, cell_size: float = 0.6, center=ORIGIN, top_left=None, stroke_color=CELL_STROKE,
                 stroke_width: float = 3, fill_color=BLACK, fill_opacity: float = 0.0, entries=None,
                 font_size=None, entry_color=WHITE):
        self.lam = tuple(lam)
        self.cell_size = cell_size
        origin = np.zeros(3) if top_left is None else np.array(top_left, dtype=float)
        cells = VGroup()
        self._index = {}
        for k, (i, j) in enumerate(young_cells(self.lam)):
            sq = Square(side_length=cell_size, stroke_color=stroke_color, stroke_width=stroke_width)
            sq.set_fill(fill_color, opacity=fill_opacity)
            sq.move_to(origin + (j + 0.5) * cell_size * RIGHT + (i + 0.5) * cell_size * DOWN)
            sq.ij = (i, j)
            cells.add(sq)
            self._index[(i, j)] = k
        self.entries = VGroup()
        self._entry_index = {}                   # (i, j) -> 在 .entries 裡的位置（.copy() 之後仍正確）
        super().__init__(cells, self.entries)
        self.cells = cells
        if top_left is None:
            self.move_to(center)
        if entries is not None:
            self.set_entries(entries, font_size=font_size, color=entry_color)

    def cell(self, i: int, j: int) -> Square:
        return self.cells[self._index[(i, j)]]

    def cell_center(self, i: int, j: int) -> np.ndarray:
        return self.cell(i, j).get_center()

    def cell_group(self, cells) -> VGroup:
        return VGroup(*[self.cell(i, j) for i, j in cells])

    def fill_cells(self, cells, color, opacity: float = 0.5):
        for i, j in cells:
            self.cell(i, j).set_fill(color, opacity=opacity)
        return self

    def set_entries(self, table, font_size=None, color=WHITE) -> VGroup:
        """table[i][j] 是數字、字串或 Mobject；None 表示空白。會換掉原本的 entries。"""
        font_size = font_size or 48 * self.cell_size
        new = VGroup()
        index = {}
        for i, row in enumerate(table):
            for j, value in enumerate(row):
                if value is None:
                    continue
                mob = _label_mob(value, font_size, color)
                mob.move_to(self.cell_center(i, j))
                mob.set_z_index(1)               # 在 hook_highlight／填色的上面
                index[(i, j)] = len(new)
                new.add(mob)
        self.remove(self.entries)
        self.entries = new
        self._entry_index = index
        self.add(new)
        return new

    def entry(self, i: int, j: int):
        k = self._entry_index.get((i, j))
        return None if k is None else self.entries[k]

    def hook_highlight(self, i: int, j: int, color=HOOK_COLOR, opacity: float = 0.35,
                       cell_color=None, cell_opacity: float = 0.7) -> VGroup:
        """鉤（自己、同列右邊、同行下面）的覆蓋方塊；第 0 個是格子本身（顏色 cell_color，預設同 color）。"""
        out = VGroup()
        for k, (a, b) in enumerate(hook_cells(self.lam, i, j)):
            sq = Square(side_length=self.cell_size * 0.86, stroke_width=0)
            col = (cell_color or color) if k == 0 else color
            sq.set_fill(col, opacity=cell_opacity if k == 0 else opacity)
            sq.move_to(self.cell_center(a, b))
            out.add(sq)
        return out


# =============================================================================
# 平面樹
# =============================================================================
class PlaneTreeMob(VGroup):
    """
    平面樹（counting 第 6 節的 tuple 形式）。根在 root_pos、往下長；同一個頂點的小孩由左到右。
    葉子依 DFS 順序等距排列（間距 leaf_gap），內部頂點放在第一個與最後一個小孩的正中間。
    子物件順序：[edges, nodes, labels]。
    屬性：
        .nodes      頂點，依前序（= counting.plane_tree_preorder、plane_tree_walk 的編號；0 是根）
        .edges      邊；.edges[k] 連到頂點 k+1（從它的父親），所以 .edge_to(v) = .edges[v-1]
        .tree_parents、.tree_depths、.node_pos（建立當下的位置）、
        .labels（dict 頂點 -> 文字，每次由 .label_group 算出，.copy() 後仍正確）、.label_group（VGroup）
    方法：.vertex(v)（物件移動後仍正確）、.edge_to(v)
    labels：{頂點: 字串或 Mobject}，放在頂點右邊。
    """

    def __init__(self, pt, root_pos=ORIGIN, leaf_gap: float = 0.6, level_height: float = 0.7,
                 node_radius: float = 0.08, node_color=CATALAN_COLOR, root_color=None,
                 edge_color=EDGE_COLOR, edge_width: float = 3, labels=None, label_font: float = 28,
                 label_color=WHITE):
        root_pos = np.array(root_pos, dtype=float)
        pre = plane_tree_preorder(pt)
        children = [[] for _ in pre]
        for v, (parent, _) in enumerate(pre):
            if parent is not None:
                children[parent].append(v)
        xs = [0.0] * len(pre)
        next_leaf = [0]

        def place(v):
            if not children[v]:
                xs[v] = next_leaf[0]
                next_leaf[0] += 1
                return
            for c in children[v]:
                place(c)
            xs[v] = (xs[children[v][0]] + xs[children[v][-1]]) / 2

        place(0)
        node_pos = [root_pos + (xs[v] - xs[0]) * leaf_gap * RIGHT + d * level_height * DOWN
                    for v, (_, d) in enumerate(pre)]
        nodes = VGroup(*[Dot(p, radius=node_radius, fill_color=node_color) for p in node_pos])
        if root_color is not None:
            nodes[0].set_fill(root_color)
        edges = VGroup()
        for v in range(1, len(pre)):
            edges.add(Line(node_pos[pre[v][0]], node_pos[v], stroke_color=edge_color, stroke_width=edge_width))
        label_mobs = VGroup()
        for v, text in (labels or {}).items():
            lab = _label_mob(text, label_font, label_color)
            lab.next_to(node_pos[v], RIGHT, buff=0.12)
            lab.tree_vertex = v
            label_mobs.add(lab)
        super().__init__(edges, nodes, label_mobs)
        self.label_group = label_mobs
        self.edges = edges
        self.nodes = nodes
        self.tree_parents = [p for p, _ in pre]
        self.tree_depths = [d for _, d in pre]
        self.node_pos = node_pos

    @property
    def labels(self) -> dict:
        """{頂點: 文字 Mobject}（由 .label_group 算出，所以 .copy() 之後仍正確）。"""
        return {lab.tree_vertex: lab for lab in self.label_group}

    def vertex(self, v: int) -> np.ndarray:
        return self.nodes[v].get_center()

    def edge_to(self, v: int) -> Line:
        assert v >= 1
        return self.edges[v - 1]


# =============================================================================
# 一排方塊（排列、堆疊、佇列、數字）
# =============================================================================
def box_row(values, size: float = 0.55, buff: float = 0.08, direction=RIGHT, colors=None, color=GREY_A,
            fill_opacity: float = 0.15, font_size=None, text_color=WHITE, stroke_width: float = 2.5,
            center=ORIGIN) -> VGroup:
    """
    一排方塊，每個方塊 = VGroup(方框, 文字)，有 .value、.square、.label。
    values 裡的 None 是空方塊；colors = 每格的顏色（list 或 dict 位置 -> 顏色），其餘用 color。
    direction=UP 可做直立的堆疊（第 0 個在最下面）。整排的中心放在 center。
    """
    font_size = font_size or 56 * size
    row = VGroup()
    for k, v in enumerate(values):
        if isinstance(colors, dict):
            col = colors.get(k, color)
        elif colors is not None:
            col = colors[k]
        else:
            col = color
        sq = Square(side_length=size, stroke_color=col, stroke_width=stroke_width)
        sq.set_fill(col, opacity=fill_opacity)
        box = VGroup(sq)
        label = None
        if v is not None:
            label = _label_mob(v, font_size, text_color).move_to(sq)
            box.add(label)
        box.value, box.square, box.label = v, sq, label
        row.add(box)
    row.arrange(direction, buff=buff)
    row.move_to(center)
    return row


def perm_row(perm, highlight=(), highlight_color=YELLOW, **kwargs) -> VGroup:
    """排列（或任何數列）的方塊列；highlight 裡的位置（0 起算）用 highlight_color 標出。"""
    colors = {k: highlight_color for k in highlight}
    return box_row(list(perm), colors=colors, **kwargs)
