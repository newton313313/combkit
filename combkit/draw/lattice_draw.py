# 由 combkit/tools/sync_from_walk017.py 從 walk017 專案的 codes/utils/lattice_draw.py 複製；請改 walk017 那份再同步。
"""
lattice_draw.py —— 《Elements of the Random Walk》第 2 章（格子上的生成函數）的畫圖元件（ManimGL）。數學在 walk_gf.py（先算、再畫）。
本模組不 import walk_draw（顏色照抄同一套），以便原樣複製進 combkit；一維的圖、長條、plot_axes、表格、筆誤框、極點
仍用 walk_draw。

座標：格點 (x, y) 畫在 origin + x·cell·RIGHT + y·cell·UP。元件都用參數決定大小，不要事後 .scale()（線寬不會跟著縮）。
顏色：SITE_COLOR 白（中心格點 y⃗）、NEIGHBOR_COLOR 灰（鄰居 w⃗_i）、TRAIL_COLOR 黃（足跡）、
      POS_COLOR 橘（正值）、NEG_COLOR 藍（負值）、ZERO_COLOR 深灰、GRID_COLOR 格線。

    lat_point(x, y, cell, origin)               格點 (x, y) 的畫面座標
    lattice_lines(n, cell, origin)              −n..n 的方格線（Fig. 2.3 的底圖）；外框在 ±(n+½)
    lattice_dots(n, cell, origin, radius)       −n..n 的格點小圓點
    count_grid(values, cell, origin, font_size, zero_dots)
                                                {(x, y): 整數} → 每個格點寫上數字（0 不寫；zero_dots=True 畫小灰點）；
                                                grid_number(g, x, y) 取出那個數字（沒有就 None）
    neighbor_star(cell, origin, labels=True)    Fig. 2.3：中心 y⃗（白）＋四個鄰居 w⃗1..w⃗4（灰，左下右上依書上 w1 左、w2 下、
                                                w3 右、w4 上）＋從鄰居指向中心的箭頭 → VGroup(center, neighbors, arrows, labels)
    trail(points, cell, origin, color, stroke_width, ends=True)
                                                足跡折線（格點或浮點座標都可以；Fig. 2.1–2.2）；起點綠點、終點紅點
    heat_color(v, vmax)                         發散色階：−vmax 藍 → 0 深灰 → +vmax 橘黃（超出就截斷）
    heatmap(fn, x_range, y_range, n, width, height, origin, vmax)
                                                在矩形 [x0,x1]×[y0,y1] 上取 n×n 個色塊，色塊 (i, j) 的值 fn(x, y)；回傳 VGroup（列優先）
    heat_cells(values, cell, origin, vmax)      {(x, y): 值} → 每個格點一個色塊（機率擴散；值 0 的格點不畫）
    bz_frame(width, origin, labels=True)        布里淵區 [−π, π]² 的方框＋q_x、q_y 軸＋±π 標籤
    roots_of_unity(n, radius, origin, m=1)      單位圓與 e^{2πimk/n}（k = 0..n−1）的箭頭 → VGroup(circle, arrows)
    root_chain(n, m, length, origin)            同一組向量頭接尾：m ≢ 0 (mod n) 圍成封閉多邊形（和為 0），m = 0 排成一直線
    iso_point(x, y, z, unit, origin)            斜投影：x 向右、z 向上、y 往右上方（縮 0.45）
    lattice3d(kind, unit, origin)               立方體 [−1,1]³ 的淡框＋中心＋最近鄰（"sc" 面心 6 個、"bcc" 角 8 個、
                                                "fcc" 稜中點 12 個）＋鍵 → VGroup(cube, bonds, neighbors, center)
    branch_cut(axes, x0, x1, color)             plot_axes 上從 x0 到 x1 的鋸齒線（分支切割）

注意：Tex 的 color= 參數無效，一律 .set_color()；同一個 play 不要同時動畫父物件與它的子物件。
"""
import numpy as np
from manimlib import *

from combkit.draw.styles import zh

SITE_COLOR = WHITE
NEIGHBOR_COLOR = GREY_B
TRAIL_COLOR = YELLOW
POS_COLOR = "#F2A33A"
NEG_COLOR = "#4A90D9"
ZERO_COLOR = "#2A2A2A"
GRID_COLOR = GREY_D
START_COLOR = "#7CE38B"
END_COLOR = RED


# ---------------------------------------------------------------- 方格與格點
def lat_point(x: float, y: float, cell: float = 0.5, origin=ORIGIN) -> np.ndarray:
    return np.array(origin, dtype=float) + x * cell * RIGHT + y * cell * UP


def lattice_lines(n: int, cell: float = 0.5, origin=ORIGIN, color=GRID_COLOR, stroke_width: float = 1.5) -> VGroup:
    g = VGroup()
    for k in range(-n, n + 1):
        g.add(Line(lat_point(-n - 0.5, k, cell, origin), lat_point(n + 0.5, k, cell, origin)))
        g.add(Line(lat_point(k, -n - 0.5, cell, origin), lat_point(k, n + 0.5, cell, origin)))
    return g.set_stroke(color, stroke_width)


def lattice_dots(n: int, cell: float = 0.5, origin=ORIGIN, radius: float = 0.035, color=GREY_C) -> VGroup:
    return VGroup(*[Dot(lat_point(x, y, cell, origin), radius=radius).set_fill(color, 1)
                    for y in range(n, -n - 1, -1) for x in range(-n, n + 1)])


def count_grid(values: dict, cell: float = 0.5, origin=ORIGIN, font_size: float = 22, color=WHITE,
               zero_dots: bool = False, n: int | None = None) -> VGroup:
    """values：{(x, y): 整數}。數字的自訂屬性 .lat_xy = (x, y)。zero_dots=True 時，−n..n 裡值為 0 的格點畫小灰點。"""
    g = VGroup()
    for (x, y), v in sorted(values.items()):
        if v:
            t = Tex(str(v), font_size=font_size).set_color(color)
            t.move_to(lat_point(x, y, cell, origin))
            t.lat_xy = (x, y)
            g.add(t)
    if zero_dots:
        n = n if n is not None else max(max(abs(x), abs(y)) for x, y in values)
        for y in range(n, -n - 1, -1):
            for x in range(-n, n + 1):
                if not values.get((x, y)):
                    d = Dot(lat_point(x, y, cell, origin), radius=0.025).set_fill(GREY_D, 1)
                    d.lat_xy = (x, y)
                    g.add(d)
    return g


def grid_number(g: VGroup, x: int, y: int):
    for m in g:
        if getattr(m, "lat_xy", None) == (x, y) and isinstance(m, Tex):
            return m
    return None


def neighbor_star(cell: float = 0.9, origin=ORIGIN, labels: bool = True, font_size: float = 26) -> VGroup:
    c = lat_point(0, 0, cell, origin)
    center = Dot(c, radius=0.13).set_fill(SITE_COLOR, 1)
    nbr_xy = [(-1, 0), (0, -1), (1, 0), (0, 1)]                     # 書上 Fig. 2.3：w1 左、w2 下、w3 右、w4 上
    neighbors = VGroup(*[Dot(lat_point(x, y, cell, origin), radius=0.11).set_fill(NEIGHBOR_COLOR, 1) for x, y in nbr_xy])
    arrows = VGroup(*[Arrow(lat_point(x, y, cell, origin), c, buff=0.16, thickness=2.5).set_color(NEIGHBOR_COLOR)
                      for x, y in nbr_xy])
    labs = VGroup()
    if labels:
        labs.add(Tex(r"\vec y", font_size=font_size).set_color(SITE_COLOR).next_to(c, DR, buff=0.06))
        for k, (x, y) in enumerate(nbr_xy):
            t = Tex(rf"\vec w_{k + 1}", font_size=font_size).set_color(NEIGHBOR_COLOR)
            t.next_to(lat_point(x, y, cell, origin), DL if x <= 0 and y <= 0 else UR, buff=0.06)
            labs.add(t)
    return VGroup(center, neighbors, arrows, labs)


def trail(points, cell: float = 0.2, origin=ORIGIN, color=TRAIL_COLOR, stroke_width: float = 3,
          ends: bool = True, end_radius: float = 0.07) -> VGroup:
    pts = [lat_point(p[0], p[1], cell, origin) for p in points]
    line = VMobject().set_points_as_corners(pts).set_stroke(color, stroke_width)
    out = VGroup(line)
    if ends:
        out.add(Dot(pts[0], radius=end_radius).set_fill(START_COLOR, 1),
                Dot(pts[-1], radius=end_radius).set_fill(END_COLOR, 1))
    return out


# ---------------------------------------------------------------- 色塊（熱圖）
def heat_color(v: float, vmax: float = 1.0):
    t = max(-1.0, min(1.0, v / vmax))
    if t >= 0:
        return interpolate_color(ZERO_COLOR, POS_COLOR, t)
    return interpolate_color(ZERO_COLOR, NEG_COLOR, -t)


def heatmap(fn, x_range=(-PI, PI), y_range=(-PI, PI), n: int = 30, width: float = 4.0, height: float | None = None,
            origin=ORIGIN, vmax: float = 1.0, opacity: float = 1.0) -> VGroup:
    """色塊 (i, j)：i 是列（由上往下）、j 是行（由左往右）；中心值 fn(x_j, y_i)。"""
    height = height or width
    x0, x1 = x_range
    y0, y1 = y_range
    w, h = width / n, height / n
    org = np.array(origin, dtype=float)
    g = VGroup()
    for i in range(n):
        y = y1 - (i + 0.5) * (y1 - y0) / n
        for j in range(n):
            x = x0 + (j + 0.5) * (x1 - x0) / n
            sq = Rectangle(width=w * 1.02, height=h * 1.02).set_stroke(width=0)
            sq.set_fill(heat_color(fn(x, y), vmax), opacity)
            sq.move_to(org + (-width / 2 + (j + 0.5) * w) * RIGHT + (height / 2 - (i + 0.5) * h) * UP)
            g.add(sq)
    return g


def heat_cells(values: dict, cell: float = 0.2, origin=ORIGIN, vmax: float = 1.0) -> VGroup:
    g = VGroup()
    for (x, y), v in sorted(values.items()):
        if v:
            sq = Square(side_length=cell * 1.02).set_stroke(width=0).set_fill(heat_color(float(v), vmax), 1)
            sq.move_to(lat_point(x, y, cell, origin))
            sq.lat_xy = (x, y)
            g.add(sq)
    return g


def bz_frame(width: float = 4.0, origin=ORIGIN, labels: bool = True, color=GREY_B, font_size: float = 24) -> VGroup:
    org = np.array(origin, dtype=float)
    box = Square(side_length=width).move_to(org).set_stroke(color, 2)
    ext = 0.35
    xa = Arrow(org + (width / 2 + ext) * LEFT, org + (width / 2 + ext) * RIGHT, buff=0, thickness=1.5).set_color(color)
    ya = Arrow(org + (width / 2 + ext) * DOWN, org + (width / 2 + ext) * UP, buff=0, thickness=1.5).set_color(color)
    g = VGroup(box, xa, ya)
    if labels:
        g.add(Tex("q_x", font_size=font_size).set_color(color).next_to(xa.get_end(), DOWN, buff=0.1),
              Tex("q_y", font_size=font_size).set_color(color).next_to(ya.get_end(), LEFT, buff=0.1),
              Tex(r"\pi", font_size=font_size).set_color(color).next_to(org + width / 2 * RIGHT, DR, buff=0.06),
              Tex(r"-\pi", font_size=font_size).set_color(color).next_to(org + width / 2 * LEFT, DL, buff=0.06),
              Tex(r"\pi", font_size=font_size).set_color(color).next_to(org + width / 2 * UP, UL, buff=0.06),
              Tex(r"-\pi", font_size=font_size).set_color(color).next_to(org + width / 2 * DOWN, DL, buff=0.06))
    return g


# ---------------------------------------------------------------- 單位根
def roots_of_unity(n: int, radius: float = 1.2, origin=ORIGIN, m: int = 1, color=TEAL) -> VGroup:
    org = np.array(origin, dtype=float)
    circ = Circle(radius=radius).move_to(org).set_stroke(GREY_C, 1.5)
    arrows = VGroup()
    for k in range(n):
        ang = 2 * PI * m * k / n
        tip = org + radius * np.array([np.cos(ang), np.sin(ang), 0])
        arrows.add(Arrow(org, tip, buff=0, thickness=2.5).set_color(color))
    return VGroup(circ, arrows)


def root_chain(n: int, m: int = 1, length: float = 0.8, origin=ORIGIN, color=TEAL) -> VGroup:
    p = np.array(origin, dtype=float)
    g = VGroup()
    for k in range(n):
        ang = 2 * PI * m * k / n
        q = p + length * np.array([np.cos(ang), np.sin(ang), 0])
        g.add(Arrow(p, q, buff=0, thickness=2.5).set_color(color))
        p = q
    return g


# ---------------------------------------------------------------- 三維晶格（斜投影）
DEPTH = np.array([0.45 * np.cos(PI / 6), 0.45 * np.sin(PI / 6), 0.0])


def iso_point(x: float, y: float, z: float, unit: float = 1.0, origin=ORIGIN) -> np.ndarray:
    return np.array(origin, dtype=float) + unit * (x * RIGHT + z * UP + y * DEPTH)


def _lattice_nbrs(kind: str):
    if kind == "sc":
        return [(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)]
    if kind == "bcc":
        return [(i, j, k) for i in (1, -1) for j in (1, -1) for k in (1, -1)]
    out = []
    for a in (1, -1):
        for b in (1, -1):
            out += [(a, b, 0), (a, 0, b), (0, a, b)]
    return out


def lattice3d(kind: str = "sc", unit: float = 1.0, origin=ORIGIN, color=TEAL, cube_color=GREY_D) -> VGroup:
    corners = [(i, j, k) for i in (-1, 1) for j in (-1, 1) for k in (-1, 1)]
    cube = VGroup()
    for a in corners:
        for b in corners:
            if a < b and sum(abs(s - t) for s, t in zip(a, b)) == 2:
                cube.add(Line(iso_point(*a, unit, origin), iso_point(*b, unit, origin)))
    cube.set_stroke(cube_color, 1.5)
    c = iso_point(0, 0, 0, unit, origin)
    nb = _lattice_nbrs(kind)
    bonds = VGroup(*[Line(c, iso_point(*p, unit, origin)) for p in nb]).set_stroke(color, 2.5, 0.8)
    dots = VGroup(*[Dot(iso_point(*p, unit, origin), radius=0.075).set_fill(color, 1) for p in nb])
    center = Dot(c, radius=0.1).set_fill(SITE_COLOR, 1)
    return VGroup(cube, bonds, dots, center)


# ---------------------------------------------------------------- 複平面的分支切割
def branch_cut(axes: Axes, x0: float, x1: float, color=RED, amp: float = 0.06, teeth: int = 14,
               stroke_width: float = 3) -> VMobject:
    p0, p1 = axes.c2p(x0, 0), axes.c2p(x1, 0)
    pts = []
    for k in range(teeth * 2 + 1):
        t = k / (teeth * 2)
        pts.append(p0 + t * (p1 - p0) + (amp if k % 2 else -amp) * UP * (0 < k < teeth * 2))
    return VMobject().set_points_as_corners(pts).set_stroke(color, stroke_width)
