"""
schubert_draw.py —— Fulton 第 9–10 章（場景 F46–F55）的圖形（ManimGL）

    STAR_COLOR、ONE_COLOR、ZERO_COLOR、BOX_COLOR    星號（自由座標）＝青、1（樞紐）＝黃、0＝灰、方框外框＝灰
    CellMatrixMob(pattern, cs=0.42)   列梯形矩陣：pattern 是 '1' '*' '0' 的二維 list（schubert.flag_cell_pattern 等）；
                                      .cell(i, j)（0 起算）、.stars()、.ones()
    BoxDiagramMob(lam, r, n, cs=0.4)  r×n 方框裡的 Young 圖（λ 填色、其餘只有格線）；.cells 是 λ 的格子
    HasseMob(levels, covers, label=perm_str, dx=1.3, dy=1.0, font_size=26)
                                      Bruhat 序的 Hasse 圖（schubert.bruhat_levels / bruhat_covers）：.node[w]、.edge[(u, v)]
    PipeDreamMob(D, m, cs=0.55)       pipe dream：階梯 i+j ≤ m 的格子，十字（D 裡）或一對彎管；.cross[(i, j)]
    poly_tex(f, **kw)                 多項式 → Tex（用 schubert.poly_str）
    word_tex(w)                       置換的一行記法字串（給 Tex 用，例如 "4\\,1\\,3\\,5\\,2"）

所有物件建好時中心在原點，用 move_to / next_to 擺位置；版面檢查用 fulton_draw.assert_in_frame。
注意：同一個 play 裡不要同時動畫「整個物件」和「它的子物件」（例如 FadeIn(mat) 加 LaggedStartMap(FadeIn, mat.stars())），
子物件會消失；先 FadeIn(mat)，下一個 play 再用 Indicate / Flash 強調 mat.stars()。
"""
from __future__ import annotations

from manimlib import *

from combkit.schubert import perm_str, poly_str

STAR_COLOR = TEAL_B
ONE_COLOR = YELLOW
ZERO_COLOR = GREY_C
BOX_COLOR = GREY_B


class CellMatrixMob(VGroup):
    def __init__(self, pattern, cs: float = 0.42, font_size: float | None = None, grid: bool = True, **kwargs):
        super().__init__(**kwargs)
        fs = font_size or 60 * cs
        self.rows, self.cols = len(pattern), len(pattern[0])
        self.entries = VGroup()
        self._grid = {}
        for i, row in enumerate(pattern):
            for j, s in enumerate(row):
                pos = np.array([(j - (self.cols - 1) / 2) * cs, ((self.rows - 1) / 2 - i) * cs, 0.0])
                color = {"1": ONE_COLOR, "*": STAR_COLOR}.get(s, ZERO_COLOR)
                t = Tex(r"\ast" if s == "*" else s, font_size=fs).set_color(color).move_to(pos)   # Tex 的 color= 參數在這版無效
                t.symbol = s
                self._grid[(i, j)] = t
                self.entries.add(t)
        w, h = self.cols * cs, self.rows * cs
        self.brackets = VGroup(
            VGroup(Line(UP * h / 2, DOWN * h / 2), Line(UP * h / 2, UP * h / 2 + RIGHT * 0.1),
                   Line(DOWN * h / 2, DOWN * h / 2 + RIGHT * 0.1)).shift(LEFT * (w / 2 + 0.05)),
            VGroup(Line(UP * h / 2, DOWN * h / 2), Line(UP * h / 2, UP * h / 2 + LEFT * 0.1),
                   Line(DOWN * h / 2, DOWN * h / 2 + LEFT * 0.1)).shift(RIGHT * (w / 2 + 0.05)),
        ).set_stroke(WHITE, 2)
        self.add(self.brackets, self.entries)

    def cell(self, i, j):
        return self._grid[(i, j)]

    def stars(self):
        return VGroup(*(t for t in self.entries if t.symbol == "*"))

    def ones(self):
        return VGroup(*(t for t in self.entries if t.symbol == "1"))


class BoxDiagramMob(VGroup):
    def __init__(self, lam, r: int, n: int, cs: float = 0.4, color=BLUE_C, fill_opacity: float = 0.55, **kwargs):
        super().__init__(**kwargs)
        lam = tuple(lam) + (0,) * (r - len(lam))
        self.grid = VGroup()
        self.cells = VGroup()
        for i in range(r):
            for j in range(n):
                sq = Square(cs).move_to(np.array([(j - (n - 1) / 2) * cs, ((r - 1) / 2 - i) * cs, 0.0]))
                if j < lam[i]:
                    sq.set_stroke(WHITE, 2).set_fill(color, fill_opacity)
                    self.cells.add(sq)
                else:
                    sq.set_stroke(BOX_COLOR, 1).set_fill(BLACK, 0)
                    self.grid.add(sq)
        self.add(self.grid, self.cells)


class HasseMob(VGroup):
    def __init__(self, levels, covers, label=perm_str, dx: float = 1.3, dy: float = 1.0, font_size: float = 26,
                 **kwargs):
        super().__init__(**kwargs)
        self.node, self.edge = {}, {}
        nodes, edges = VGroup(), VGroup()
        for k, lv in enumerate(levels):
            for a, w in enumerate(lv):
                pos = np.array([(a - (len(lv) - 1) / 2) * dx, (k - (len(levels) - 1) / 2) * dy, 0.0])
                t = Tex(label(w), font_size=font_size).move_to(pos)
                self.node[w] = t
                nodes.add(t)
        for u, v in covers:
            ln = Line(self.node[u].get_top(), self.node[v].get_bottom(), buff=0.06).set_stroke(GREY_B, 1.5)
            self.edge[(u, v)] = ln
            edges.add(ln)
        self.nodes, self.edges = nodes, edges
        self.add(edges, nodes)


class PipeDreamMob(VGroup):
    """格子 (i, j)（1 起算，i 列 j 行，i+j ≤ m）。十字畫兩條線，其餘畫一對彎管（左上、右下的四分之一圓）。"""

    def __init__(self, D, m: int, cs: float = 0.55, cross_color=YELLOW, pipe_color=GREY_A, **kwargs):
        super().__init__(**kwargs)
        D = set(D)
        self.cross = {}
        self.frame_cells = VGroup()
        self.tiles = VGroup()
        for i in range(1, m + 1):
            for j in range(1, m + 2 - i):
                c = np.array([(j - (m + 1) / 2) * cs, ((m + 1) / 2 - i) * cs, 0.0])
                sq = Square(cs).move_to(c).set_stroke(GREY_D, 1)
                self.frame_cells.add(sq)
                tl = Arc(start_angle=-PI / 2, angle=PI / 2, radius=cs / 2,          # 連左邊與上邊（圓心在左上角）
                         arc_center=c + np.array([-cs / 2, cs / 2, 0]))
                if i + j > m:                                   # 階梯外緣的半格：只畫一個彎管
                    tile = tl.set_stroke(pipe_color, 3)
                elif (i, j) in D:
                    tile = VGroup(Line(c + LEFT * cs / 2, c + RIGHT * cs / 2), Line(c + UP * cs / 2, c + DOWN * cs / 2))
                    tile.set_stroke(cross_color, 4)
                    self.cross[(i, j)] = tile
                else:
                    br = Arc(start_angle=PI / 2, angle=PI / 2, radius=cs / 2,       # 連右邊與下邊（圓心在右下角）
                             arc_center=c + np.array([cs / 2, -cs / 2, 0]))
                    tile = VGroup(tl, br).set_stroke(pipe_color, 3)
                self.tiles.add(tile)
        self.add(self.frame_cells, self.tiles)


def poly_tex(f, **kwargs) -> Tex:
    return Tex(poly_str(f), **kwargs)


def word_tex(w) -> str:
    return r"\,".join(map(str, w))
