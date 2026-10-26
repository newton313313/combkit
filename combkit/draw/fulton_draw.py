"""
fulton_draw.py —— Fulton《Young Tableaux》第 2–4 章場景（F09–F19）共用的畫圖工具（ManimGL）

    版面檢查      assert_in_frame(*mobs)、assert_apart(a, b, gap)、X_MIN.. 常數（字幕在 y=-3.45、標題在左上）
    表            tableau_mob(T, cs, fs, ...)：填好數字的 YoungDiagramMob（.entry(i,j).tab_value = 數字）
                  SkewMob(S, cs, ...)：斜表／有洞的表（前面的 None = μ，畫成灰格；中間的 None = 洞，畫成虛線格）
    字            WordMob(w, ...)：一排字母方塊；.apply_move(scene, i, rule) 播放一步基本 Knuth 變換
                  word_tex(w, colors)：只有文字的字（Tex，字母之間有空隙）
    陣列、矩陣    TwoRowArray(top, bottom)：雙列陣列（.col(k)）；BallMatrixMob(A, numbering)：矩陣球圖
    其他          tuple_str(λ) → "(3,1)"；labeled(tex, mob)：左邊一個公式、右邊一個物件

格子 (i, j) 一律 0 起算（同 tableaux.py）；字幕裡給人看的列號、行號用 1 起算。
"""
from __future__ import annotations

from manimlib import *

from combkit.draw.mit_draw import CELL_STROKE, YoungDiagramMob
from combkit.draw.styles import MUTED_COLOR
from combkit.tableaux import KNUTH_RULES, shape, young_boxes

X_MIN, X_MAX, Y_MIN, Y_MAX = -6.8, 6.8, -3.0, 3.2
HOLE_COLOR = YELLOW
MU_COLOR = GREY_D
RULE_COLORS = {"K'": TEAL, "K'⁻¹": TEAL, "K''": PINK, "K''⁻¹": PINK}
RULE_TEX = {"K'": "K'", "K'⁻¹": "K'^{-1}", "K''": "K''", "K''⁻¹": "K''^{-1}"}


# =============================================================================
# 版面檢查
# =============================================================================
def assert_in_frame(*mobs):
    """物件都在畫面的安全範圍內（不碰到標題與字幕）。"""
    for m in mobs:
        assert m.get_left()[0] >= X_MIN - 1e-6 and m.get_right()[0] <= X_MAX + 1e-6, m
        assert m.get_bottom()[1] >= Y_MIN - 1e-6 and m.get_top()[1] <= Y_MAX + 1e-6, m


def assert_apart(a, b, gap=0.05):
    """兩個物件的外框不重疊（至少一個方向分開 gap）。"""
    apart_x = a.get_right()[0] + gap <= b.get_left()[0] or b.get_right()[0] + gap <= a.get_left()[0]
    apart_y = a.get_top()[1] + gap <= b.get_bottom()[1] or b.get_top()[1] + gap <= a.get_bottom()[1]
    assert apart_x or apart_y, (a.get_corner(UL).round(2), a.get_corner(DR).round(2),
                                b.get_corner(UL).round(2), b.get_corner(DR).round(2))


def tuple_str(lam) -> str:
    return "(" + ",".join(str(p) for p in lam) + ")"


def labeled(tex: str, mob: Mobject, font_size=36, buff=0.25) -> VGroup:
    """VGroup(公式, 物件)，左右排好；公式的 Tex 用 tex。"""
    return VGroup(Tex(tex, font_size=font_size), mob).arrange(RIGHT, buff=buff)


# =============================================================================
# 表與斜表
# =============================================================================
def tableau_mob(T, cs: float = 0.55, fs=None, top_left=None, center=ORIGIN, **kwargs) -> YoungDiagramMob:
    """填好數字的 Young 圖；每個數字物件記下 .tab_value。空表傳回只有一個點的 VGroup 不適用，請自行處理。"""
    kw = {"center": center} if top_left is None else {"top_left": np.array(top_left, dtype=float)}
    d = YoungDiagramMob(shape(T), cell_size=cs, entries=T, font_size=fs or 56 * cs, **kw, **kwargs)
    for i, j in young_boxes(shape(T)):
        d.entry(i, j).tab_value = T[i][j]
    return d


class SkewMob(VGroup):
    """
    斜表或有洞的表：S 的每一列可以有 None。列首連續的 None 是 μ 的格子（灰色，show_mu=False 時不畫），
    其他的 None 是「洞」（黃色虛線框）。子物件：[cells, entries]。
    .cell(i, j)、.entry(i, j)（沒有數字時 None）、.cell_center(i, j)、.boxes（有畫出的格子座標）
    """

    def __init__(self, S, cs: float = 0.55, fs=None, top_left=None, center=ORIGIN, show_mu: bool = True,
                 stroke_color=CELL_STROKE, entry_color=WHITE):
        super().__init__()
        self.cs = cs
        origin = np.zeros(3) if top_left is None else np.array(top_left, dtype=float)
        self.cells, self.entries = VGroup(), VGroup()
        self._cell, self._entry = {}, {}
        for i, row in enumerate(S):
            lead = next((j for j, v in enumerate(row) if v is not None), len(row))
            for j, v in enumerate(row):
                pos = origin + (j + 0.5) * cs * RIGHT + (i + 0.5) * cs * DOWN
                if v is None and j < lead:                       # μ 的格子
                    if not show_mu:
                        continue
                    sq = Square(side_length=cs, stroke_color=MU_COLOR, stroke_width=1.5)
                    sq.set_fill(MU_COLOR, opacity=0.35)
                elif v is None:                                  # 洞
                    sq = DashedVMobject(Square(side_length=cs * 0.92), num_dashes=16)
                    sq.set_stroke(HOLE_COLOR, width=2.5)
                else:
                    sq = Square(side_length=cs, stroke_color=stroke_color, stroke_width=3)
                    lab = Tex(str(v), font_size=fs or 56 * cs, fill_color=entry_color)
                    lab.move_to(pos).set_z_index(1)
                    lab.tab_value = v
                    self._entry[(i, j)] = lab
                    self.entries.add(lab)
                sq.move_to(pos)
                self._cell[(i, j)] = sq
                self.cells.add(sq)
        self.boxes = list(self._cell)
        self.add(self.cells, self.entries)
        if top_left is None:
            self.move_to(center)

    def cell(self, i, j):
        return self._cell[(i, j)]

    def entry(self, i, j):
        return self._entry.get((i, j))

    def cell_center(self, i, j):
        return self._cell[(i, j)].get_center()


# =============================================================================
# 字
# =============================================================================
def word_tex(w, colors=None, font_size=36, gap=r"\ ") -> Tex:
    """字的文字版（例如 5\\,6\\,4…）；colors = {位置: 顏色}。字母依序是 tex[k]（每個字母都是一位數時）。"""
    assert all(0 < a < 10 for a in w), "word_tex 只處理一位數字母"
    t = Tex(gap.join(str(a) for a in w), font_size=font_size)
    for k, col in (colors or {}).items():
        t[k].set_fill(col)
    return t


class WordMob(VGroup):
    """
    一排字母方塊（同 mit_draw.box_row 的樣式）。位置（slot）固定，方塊會移動：
        .tiles[k]    目前在第 k 個位置的方塊（VGroup(方框, 文字)，.value = 字母）
        .letters     目前的字（tuple）
        .slot(k)     第 k 個位置的中心
        .apply_move(scene, i, rule, run_time)   播放一步基本 Knuth 變換（框出三個字母、交換其中兩個）
        .bracket(i, n, label, color)            框住位置 i..i+n-1 的圓角框（加上下方標籤），不加入圖中
        .color_slots({k: 顏色})                 把位置 k 的方塊改色（回傳動畫 list）
    """

    def __init__(self, w, size: float = 0.55, buff: float = 0.08, colors=None, color=GREY_A,
                 fill_opacity: float = 0.15, font_size=None, center=ORIGIN):
        super().__init__()
        self.size = size
        font_size = font_size or 56 * size
        tiles = []
        for k, v in enumerate(w):
            col = (colors or {}).get(k, color) if isinstance(colors, dict) else (colors[k] if colors else color)
            sq = Square(side_length=size, stroke_color=col, stroke_width=2.5).set_fill(col, opacity=fill_opacity)
            lab = Tex(str(v), font_size=font_size).move_to(sq)
            tile = VGroup(sq, lab)
            tile.value, tile.square, tile.label = v, sq, lab
            tiles.append(tile)
        self.add(*tiles)
        self.arrange(RIGHT, buff=buff).move_to(center)
        self.tiles = list(tiles)
        self.step = size + buff

    @property
    def letters(self):
        return tuple(t.value for t in self.tiles)

    def slot(self, k) -> np.ndarray:
        return self.tiles[0].get_center() + k * self.step * RIGHT

    def bracket(self, i, n=3, label=None, color=YELLOW, font_size=26) -> VGroup:
        box = RoundedRectangle(width=n * self.step + 0.06, height=self.size + 0.22, corner_radius=0.12)
        box.set_stroke(color, width=3).move_to(self.slot(i + (n - 1) / 2))
        out = VGroup(box)
        if label is not None:
            out.add(Tex(label, font_size=font_size, fill_color=color).next_to(box, DOWN, buff=0.12))
        return out

    def swap_anims(self, a, b, run_time=0.6):
        ta, tb = self.tiles[a], self.tiles[b]
        pa, pb = self.slot(a), self.slot(b)
        self.tiles[a], self.tiles[b] = tb, ta
        return [ta.animate(path_arc=-PI / 2).move_to(pb), tb.animate(path_arc=-PI / 2).move_to(pa)]

    def apply_move(self, scene, i, rule, run_time=0.6, hold=0.25, show_label=True):
        """一步基本 Knuth 變換：K'／K'⁻¹ 交換三個字母中的後兩個，K''／K''⁻¹ 交換前兩個。回傳新字。"""
        assert rule in KNUTH_RULES
        color = RULE_COLORS[rule]
        br = self.bracket(i, 3, RULE_TEX[rule] if show_label else None, color=color)
        scene.play(FadeIn(br), run_time=0.3)
        a = i + 1 if rule.startswith("K'") and not rule.startswith("K''") else i
        scene.play(*self.swap_anims(a, a + 1), run_time=run_time)
        if hold:
            scene.wait(hold)
        scene.play(FadeOut(br), run_time=0.25)
        return self.letters

    def color_slots(self, colors: dict, opacity=0.3):
        return [self.tiles[k].square.animate.set_fill(c, opacity=opacity).set_stroke(c) for k, c in colors.items()]


# =============================================================================
# 雙列陣列與矩陣球
# =============================================================================
class TwoRowArray(VGroup):
    """雙列陣列（上列 top、下列 bottom），左右有圓括號。.col(k) = VGroup(上, 下)，.cols。"""

    def __init__(self, top, bottom, font_size=36, buff=0.32, center=ORIGIN, top_color=WHITE, bottom_color=WHITE):
        super().__init__()
        cols = VGroup()
        for u, v in zip(top, bottom):
            col = VGroup(Tex(str(u), font_size=font_size, fill_color=top_color),
                         Tex(str(v), font_size=font_size, fill_color=bottom_color)).arrange(DOWN, buff=0.28)
            cols.add(col)
        cols.arrange(RIGHT, buff=buff)
        h = cols.get_height() + 0.25
        lp = Tex("(", font_size=font_size).stretch_to_fit_height(h)
        rp = Tex(")", font_size=font_size).stretch_to_fit_height(h)
        lp.next_to(cols, LEFT, buff=0.12)
        rp.next_to(cols, RIGHT, buff=0.12)
        self.cols = cols
        self.add(lp, cols, rp)
        self.move_to(center)

    def col(self, k):
        return self.cols[k]


class BallMatrixMob(VGroup):
    """
    矩陣球圖：m×n 的格子，(i, j) 裡畫 a(i,j) 個球（由西北往東南斜排），numbering = tableaux.ball_numbering(A)
    給定時球上寫號碼。列號、行號（1 起算）畫在左邊與上面。
    .ball(i, j, t)：(i, j) 的第 t 個球（VGroup(圓, 號碼)，.number）；.balls：全部的球；.cell_center(i, j)
    """

    def __init__(self, A, numbering=None, cell: float = 0.95, ball_radius: float = 0.17, center=ORIGIN,
                 color=BLUE_C, font_size=22, show_labels=True):
        super().__init__()
        m, n = len(A), len(A[0])
        self.m, self.n, self.cell = m, n, cell
        grid = VGroup()
        for i in range(m):
            for j in range(n):
                sq = Square(side_length=cell, stroke_color=GREY_C, stroke_width=1.5)
                sq.move_to((j + 0.5) * cell * RIGHT + (i + 0.5) * cell * DOWN)
                grid.add(sq)
        self.grid = grid
        self._centers = {(i, j): grid[i * n + j].get_center() for i in range(m) for j in range(n)}
        labels = VGroup()
        if show_labels:
            for i in range(m):
                labels.add(Tex(str(i + 1), font_size=26, fill_color=MUTED_COLOR)
                           .next_to(grid[i * n], LEFT, buff=0.18))
            for j in range(n):
                labels.add(Tex(str(j + 1), font_size=26, fill_color=MUTED_COLOR)
                           .next_to(grid[j], UP, buff=0.14))
        self.labels = labels
        self.balls = VGroup()
        self._balls = {}
        for i in range(m):
            for j in range(n):
                a = A[i][j]
                for t in range(a):
                    off = (t - (a - 1) / 2) * min(0.3, 0.75 * cell / max(a, 1))
                    pos = self._centers[(i, j)] + off * (RIGHT + DOWN)
                    circ = Circle(radius=ball_radius).set_fill(color, opacity=0.85).set_stroke(WHITE, width=1.5)
                    circ.move_to(pos)
                    ball = VGroup(circ)
                    ball.number = None
                    if numbering is not None:
                        ball.number = numbering[(i, j)][t]
                        ball.add(Tex(str(ball.number), font_size=font_size, fill_color=BLACK).move_to(pos))
                    self._balls[(i, j, t)] = ball
                    self.balls.add(ball)
        self.add(grid, labels, self.balls)
        self.move_to(center)

    def ball(self, i, j, t=0):
        return self._balls[(i, j, t)]

    def cell_square(self, i, j):
        return self.grid[i * self.n + j]

    def cell_center(self, i, j):
        return self.cell_square(i, j).get_center()
