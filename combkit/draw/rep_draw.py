"""
rep_draw.py —— Fulton 第 7 章場景（F28–F41）共用的畫圖元件（ManimGL）

    顏色        ROW_COLOR（列、a_T）、COL_COLOR（行、b_T）、OP_COLOR（目前操作）、POS_COLOR（正號／標準基底）、NEG_COLOR（負號）
    TabloidMob(rows)            列 tabloid：只畫各列上下的橫線（同一列內順序不重要）；.key、.entry(i, j)
    ColumnTabloidMob(rows, s)   帶方向的行 tabloid：每行左右畫直線、前面標 s 的符號；.key、.sign
    signed_sum(terms)           [(係數, Mobject), ...] → 「a·X + b·Y − …」，係數 ±1 只畫符號；正綠、負紅紫
    MatrixMob(rows)             整數矩陣（括號自動拉高）；.entry(i, j)、.col(j)、.row(i)
    perm_tex(p)                 置換的循環記法 Tex 字串（恆等寫 e）；oneline_tex(p) → "[2,1,3]"

畫面上的數字一律由 combkit.rep_data 算出；本檔只負責畫。版面檢查用 combkit.draw.fulton_draw 的 assert_in_frame／assert_apart。
"""
from __future__ import annotations

from manimlib import *

from combkit.rep_data import cycles, tabloid_key

ROW_COLOR = BLUE_C
COL_COLOR = ORANGE
OP_COLOR = YELLOW
POS_COLOR = GREEN_C
NEG_COLOR = "#E05AA0"


def perm_tex(p) -> str:
    cs = cycles(p, keep_fixed=False)
    return "e" if not cs else "".join("(" + "\\,".join(map(str, c)) + ")" for c in cs)


def oneline_tex(p) -> str:
    return "[" + ",".join(map(str, p)) + "]"


class TabloidMob(VGroup):
    """列 tabloid {T}：數字排成列，每列上下各一條橫線（Fulton p.85 的畫法）。"""

    def __init__(self, rows, cs: float = 0.5, font_size=None, color=WHITE, line_color=GREY_B, center=ORIGIN):
        super().__init__()
        self.key = tabloid_key(rows)
        fs = font_size or 60 * cs
        self.entries = VGroup()
        self._entry = {}
        lines = VGroup()
        width = max(len(r) for r in self.key) * cs
        for i, r in enumerate(self.key):
            y = -i * cs
            for j, v in enumerate(r):
                t = Tex(str(v), font_size=fs, fill_color=color).move_to(np.array([(j + 0.5) * cs, y, 0.0]))
                self._entry[(i, j)] = t
                self.entries.add(t)
            lines.add(Line(np.array([0, y + cs / 2, 0]), np.array([len(r) * cs, y + cs / 2, 0]),
                           stroke_color=line_color, stroke_width=2))
        last = len(self.key) - 1
        lines.add(Line(np.array([0, -last * cs - cs / 2, 0]), np.array([len(self.key[-1]) * cs, -last * cs - cs / 2, 0]),
                       stroke_color=line_color, stroke_width=2))
        self.lines = lines
        self.add(lines, self.entries)
        self.move_to(center)
        assert width > 0

    def entry(self, i, j):
        return self._entry[(i, j)]


class ColumnTabloidMob(VGroup):
    """帶方向的行 tabloid [T]（Fulton p.95）：每行左右各一條直線；sign = ±1 畫在前面。"""

    def __init__(self, rows, sign_: int = 1, cs: float = 0.5, font_size=None, line_color=GREY_B, center=ORIGIN,
                 show_plus=False):
        super().__init__()
        self.rows = tuple(map(tuple, rows))
        self.sign = sign_
        fs = font_size or 60 * cs
        body = VGroup()
        ncol = len(self.rows[0])
        for j in range(ncol):
            col = [r[j] for r in self.rows if len(r) > j]
            for i, v in enumerate(col):
                body.add(Tex(str(v), font_size=fs).move_to(np.array([(j + 0.5) * cs, -i * cs, 0.0])))
            top, bot = cs / 2, -(len(col) - 0.5) * cs
            for x in (j * cs + 0.06, (j + 1) * cs - 0.06):
                body.add(Line(np.array([x, top, 0]), np.array([x, bot, 0]), stroke_color=line_color, stroke_width=2))
        self.body = body
        self.add(body)
        if sign_ < 0 or show_plus:
            s = Tex("-" if sign_ < 0 else "+", font_size=fs, fill_color=NEG_COLOR if sign_ < 0 else POS_COLOR)
            s.next_to(body, LEFT, buff=0.12)
            self.add(s)
        self.move_to(center)


def signed_sum(terms, font_size=40, buff=0.2, dot=False) -> VGroup:
    """terms = [(係數, Mobject), ...]；回傳 VGroup，.terms[k] = VGroup(係數記號, 物件)。"""
    out = VGroup()
    out.terms = []
    for k, (c, mob) in enumerate(terms):
        if c == 0:
            continue
        mag = abs(c)
        s = ("-" if c < 0 else ("+" if k else "")) + ("" if mag == 1 else str(mag))
        part = VGroup()
        if s:
            part.add(Tex(s + (r"\cdot" if dot and mag != 1 else ""), font_size=font_size,
                         fill_color=NEG_COLOR if c < 0 else POS_COLOR))
        part.add(mob)
        part.arrange(RIGHT, buff=0.12)
        out.terms.append(part)
        out.add(part)
    out.arrange(RIGHT, buff=buff)
    return out


class MatrixMob(VGroup):
    """整數矩陣（rows）；括號依高度拉長。.entry(i, j)、.col(j)、.row(i)。"""

    def __init__(self, rows, font_size=36, h_buff=0.75, v_buff=0.5, colors=None, center=ORIGIN):
        super().__init__()
        self.m, self.n = len(rows), len(rows[0])
        grid = VGroup()
        self._e = {}
        for i, r in enumerate(rows):
            for j, v in enumerate(r):
                col = (colors or {}).get((i, j), WHITE)
                t = Tex(str(v), font_size=font_size, fill_color=col).move_to(np.array([j * h_buff, -i * v_buff, 0.0]))
                self._e[(i, j)] = t
                grid.add(t)
        h = (self.m - 1) * v_buff + 0.5
        lb = Tex("(", font_size=font_size).stretch_to_fit_height(h).next_to(grid, LEFT, buff=0.12)
        rb = Tex(")", font_size=font_size).stretch_to_fit_height(h).next_to(grid, RIGHT, buff=0.12)
        self.grid = grid
        self.add(lb, grid, rb)
        self.move_to(center)

    def entry(self, i, j):
        return self._e[(i, j)]

    def col(self, j):
        return VGroup(*[self._e[(i, j)] for i in range(self.m)])

    def row(self, i):
        return VGroup(*[self._e[(i, j)] for j in range(self.n)])
