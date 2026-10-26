# 由 combkit/tools/sync_from_walk017.py 從 walk017 專案的 codes/utils/walk_draw.py 複製；請改 walk017 那份再同步。
"""
walk_draw.py —— 《Elements of the Random Walk》第 1 章的畫圖元件（ManimGL）。數學在 walk1d.py（先算、再畫）。

放置方式：元件都用「原點 origin＋單位長 ux、uy」直接排好，不要事後 .scale()（線寬不會跟著縮）。
顏色：HEAD_COLOR 金（正面／向右）、TAIL_COLOR 藍灰（反面／向左）、WALK_COLOR 黃（主例路線）、
      BAR_COLOR 藍（精確的長條）、GAUSS_COLOR 橘（高斯近似）、EXACT_COLOR 綠（精確值）、
      WRONG_COLOR 紅（書上的筆誤）、RIGHT_COLOR 綠（更正）、MUTED 灰。

    coin(face, radius=0.3)                      'H' 金色「正」、'T' 藍灰「反」→ VGroup(圓, 內圈, 字)
    coin_row(coins, radius=0.3, buff=0.12)      例如 coin_row("HTH")（往右排，置中在 ORIGIN）
    number_line(N, unit=0.5, label_every=2)     −N..N 的數線 VGroup(底線, 刻度, 數字)；line_point(nl, x) → 座標 x 的點
    walker_dot(color=WALK_COLOR, radius=0.12)   走路的人（圓點＋外圈）
    walk_axes(N, ymax, ux, uy, origin, ...)     走法圖的座標軸：橫軸步數 0..N、縱軸位置 −ymax..ymax，附刻度、中文軸名
    path_point(k, x, ux, uy, origin)            走法圖上「第 k 步、位置 x」的點
    walk_path(pos, ux, uy, origin, color, ...)  走法圖折線（pos = walk1d.positions(steps)，長 N+1）→ VGroup(折線, 圓點們)
    fan(pos_list, ux, uy, origin, ...)          很多條路線疊在一起（細、半透明）→ VGroup
    sqrt_envelope(Nmax, ux, uy, origin, c=1)    ±c·√N 兩條曲線 → VGroup(上, 下)
    bar_chart(values, ux, uy, origin, ...)      {x: 高} → 長條 VGroup(長條們, 標籤們)；底在 origin 的水平線、中心 x 在 origin + x·ux
    bar_of(chart, x)                            取出 x 那一根長條（沒有就 None）
    curve(fn, x0, x1, ux, uy, origin, ...)      y = fn(x) 的曲線（與 bar_chart 同一套座標）；dashed=True 成虛線
    plot_axes(x_range, y_range, width, height, x_ticks, y_ticks, ...)
                                                一般函數圖：ManimGL Axes（c2p、get_graph 照用），刻度數字用 Tex 自己加
    pascal_triangle(rows, dx, dy, origin)       rows = walk1d.pascal_rows(n)；C(N,y) 放在 origin + y·dx·RIGHT − N·dy·UP；
                                                pascal_entry(tri, N, y) 取出那個數字
    delta_grid(N, cell, origin)                 (1.12) 的 N×N 格子 VGroup(格子們)；對角格亮（Δ_j²）、非對角暗（Δ_jΔ_k）；
                                                grid_cell(g, j, k) 取格子（0 起算）
    tex_table(rows, col_widths, row_h, font_size)  Tex 表格（rows[0] 是表頭；中文表頭請直接放 Mobject）；table_cell(t, i, j)
    erratum_box(where, wrong_tex, right_tex, width=5.6)
                                                「書上筆誤」框：第一行 where（中文，例如「p. 11，(1.34)」），
                                                第二行 錯（紅色＋刪除線）→ 對（綠色）
    contour_circle(axes, r, color)              plot_axes 上繞原點的逆時針圓周（附箭頭）
    pole_marker(axes, z, label_tex, color)      × 號＋標籤

注意：Tex 的 color= 參數無效，一律 .set_color()；同一個 play 不要同時動畫父物件與它的子物件。
"""
import numpy as np
from manimlib import *

from combkit.draw.styles import zh

HEAD_COLOR = "#F2C14E"
TAIL_COLOR = "#7FA7C9"
WALK_COLOR = YELLOW
BAR_COLOR = BLUE_C
GAUSS_COLOR = ORANGE
EXACT_COLOR = "#7CE38B"
WRONG_COLOR = RED
RIGHT_COLOR = "#7CE38B"
MUTED = GREY_B


# ---------------------------------------------------------------- 硬幣與數線
def coin(face: str = "H", radius: float = 0.3) -> VGroup:
    color = HEAD_COLOR if face == "H" else TAIL_COLOR
    disk = Circle(radius=radius).set_fill(color, 0.95).set_stroke(WHITE, 2)
    rim = Circle(radius=radius * 0.8).set_stroke(BLACK, 1.5, 0.35)
    label = zh("正" if face == "H" else "反", font_size=radius * 72, color=BLACK)
    label.move_to(disk)
    return VGroup(disk, rim, label)


def coin_row(coins: str, radius: float = 0.3, buff: float = 0.12) -> VGroup:
    return VGroup(*[coin(c, radius) for c in coins]).arrange(RIGHT, buff=buff)


def number_line(N: int, unit: float = 0.5, label_every: int = 2, font_size: float = 22,
                color=MUTED, labels: bool = True) -> VGroup:
    base = Line((-N - 0.5) * unit * RIGHT, (N + 0.5) * unit * RIGHT).set_stroke(color, 2)
    ticks = VGroup(*[Line(0.08 * DOWN, 0.08 * UP).move_to(x * unit * RIGHT).set_stroke(color, 2)
                     for x in range(-N, N + 1)])
    nums = VGroup()
    if labels:
        for x in range(-N, N + 1):
            if x % label_every == 0:
                t = Tex(str(x), font_size=font_size).set_color(color)
                t.next_to(x * unit * RIGHT, DOWN, buff=0.16)
                nums.add(t)
    return VGroup(base, ticks, nums)


def line_point(nl: VGroup, x: float) -> np.ndarray:
    """數線 nl 上座標 x 的點（nl 可以移動）。"""
    ticks = nl[1]
    n = len(ticks)
    step = ticks[1].get_center() - ticks[0].get_center()
    return ticks[n // 2].get_center() + x * step


def walker_dot(color=WALK_COLOR, radius: float = 0.12) -> VGroup:
    dot = Dot(radius=radius).set_fill(color, 1)
    ring = Circle(radius=radius * 1.7).set_stroke(color, 2, 0.6)
    return VGroup(dot, ring)


# ---------------------------------------------------------------- 走法圖
def path_point(k: float, x: float, ux: float, uy: float, origin=ORIGIN) -> np.ndarray:
    return np.array(origin) + k * ux * RIGHT + x * uy * UP


def walk_axes(N: int, ymax: int, ux: float = 0.5, uy: float = 0.4, origin=ORIGIN,
              x_every: int | None = None, y_every: int | None = None, font_size: float = 20,
              x_label: str = "步數", y_label: str = "位置", color=MUTED) -> VGroup:
    origin = np.array(origin)
    x_every = x_every or max(1, N // 5)
    y_every = y_every or max(1, ymax // 2)
    xa = Arrow(origin, origin + (N + 0.9) * ux * RIGHT, buff=0, thickness=2.0).set_color(color)
    ya = Arrow(origin + (ymax + 0.6) * uy * DOWN, origin + (ymax + 0.9) * uy * UP, buff=0,
               thickness=2.0).set_color(color)
    ticks, nums = VGroup(), VGroup()
    for k in range(x_every, N + 1, x_every):
        p = path_point(k, 0, ux, uy, origin)
        ticks.add(Line(p + 0.06 * DOWN, p + 0.06 * UP).set_stroke(color, 2))
        nums.add(Tex(str(k), font_size=font_size).set_color(color).next_to(p, DOWN, buff=0.12))
    for y in range(-ymax, ymax + 1, y_every):
        if y == 0:
            continue
        p = path_point(0, y, ux, uy, origin)
        ticks.add(Line(p + 0.06 * LEFT, p + 0.06 * RIGHT).set_stroke(color, 2))
        nums.add(Tex(str(y), font_size=font_size).set_color(color).next_to(p, LEFT, buff=0.12))
    xl = zh(x_label, font_size=font_size, color=color).next_to(xa.get_end(), DOWN, buff=0.12)
    yl = zh(y_label, font_size=font_size, color=color).next_to(ya.get_end(), LEFT, buff=0.12)
    return VGroup(xa, ya, ticks, nums, xl, yl)


def walk_path(pos, ux: float = 0.5, uy: float = 0.4, origin=ORIGIN, color=WALK_COLOR,
              stroke_width: float = 4, dots: bool = True, dot_radius: float = 0.05) -> VGroup:
    pts = [path_point(k, x, ux, uy, origin) for k, x in enumerate(pos)]
    line = VMobject().set_points_as_corners(pts).set_stroke(color, stroke_width)
    ds = VGroup(*[Dot(p, radius=dot_radius).set_fill(color, 1) for p in pts]) if dots else VGroup()
    return VGroup(line, ds)


def fan(pos_list, ux: float = 0.05, uy: float = 0.2, origin=ORIGIN, color=BLUE_B,
        stroke_width: float = 1.0, opacity: float = 0.35) -> VGroup:
    out = VGroup()
    for pos in pos_list:
        pts = [path_point(k, x, ux, uy, origin) for k, x in enumerate(pos)]
        out.add(VMobject().set_points_as_corners(pts).set_stroke(color, stroke_width, opacity))
    return out


def sqrt_envelope(Nmax: int, ux: float, uy: float, origin=ORIGIN, c: float = 1.0, color=GAUSS_COLOR,
                  stroke_width: float = 3, samples: int = 120) -> VGroup:
    ks = np.linspace(0, Nmax, samples)
    up = VMobject().set_points_smoothly([path_point(k, c * np.sqrt(k), ux, uy, origin) for k in ks])
    dn = VMobject().set_points_smoothly([path_point(k, -c * np.sqrt(k), ux, uy, origin) for k in ks])
    return VGroup(up.set_stroke(color, stroke_width), dn.set_stroke(color, stroke_width))


# ---------------------------------------------------------------- 長條圖與曲線
def bar_chart(values: dict, ux: float = 0.4, uy: float = 1.0, origin=ORIGIN, color=BAR_COLOR,
              width: float = 0.8, fill_opacity: float = 0.75, labels: bool = False, label_fmt=str,
              font_size: float = 20, label_color=WHITE) -> VGroup:
    origin = np.array(origin)
    bars, labs = VGroup(), VGroup()
    for x in sorted(values):
        h = float(values[x]) * uy
        r = Rectangle(width=width * ux, height=max(h, 1e-3))
        r.set_fill(color, fill_opacity).set_stroke(color, 1.5)
        r.move_to(origin + x * ux * RIGHT + h / 2 * UP)
        r.bar_x = x                     # 自訂屬性（不和 ManimGL 撞名）
        bars.add(r)
        if labels:
            t = Tex(label_fmt(values[x]), font_size=font_size).set_color(label_color)
            t.next_to(r, UP, buff=0.08)
            labs.add(t)
    return VGroup(bars, labs)


def bar_of(chart: VGroup, x):
    for r in chart[0]:
        if getattr(r, "bar_x", None) == x:
            return r
    return None


def curve(fn, x0: float, x1: float, ux: float = 0.4, uy: float = 1.0, origin=ORIGIN, color=GAUSS_COLOR,
          stroke_width: float = 3, samples: int = 200, dashed: bool = False, num_dashes: int = 40) -> VMobject:
    xs = np.linspace(x0, x1, samples)
    pts = [np.array(origin) + x * ux * RIGHT + fn(x) * uy * UP for x in xs]
    c = VMobject().set_points_smoothly(pts).set_stroke(color, stroke_width)
    if dashed:
        c = DashedVMobject(c, num_dashes=num_dashes).set_stroke(color, stroke_width)
    return c


def plot_axes(x_range, y_range, width: float = 6.0, height: float = 3.5, x_ticks=(), y_ticks=(),
              x_tick_labels=None, y_tick_labels=None, font_size: float = 20, color=MUTED,
              x_label: str | None = None, y_label: str | None = None) -> Axes:
    """x_range=(a, b, step)。刻度數字用 Tex 字串（x_tick_labels 可給對應的 TeX，否則用 str(值)）。
    x_label／y_label 是 TeX（例如 "w"、"q"）。回傳的 Axes 可以 c2p、get_graph。"""
    ax = Axes(x_range=x_range, y_range=y_range, width=width, height=height,
              axis_config=dict(stroke_color=color, stroke_width=2, include_tip=True, tick_size=0.05))
    labs = VGroup()
    for k, v in enumerate(x_ticks):
        s = x_tick_labels[k] if x_tick_labels else str(v)
        labs.add(Tex(s, font_size=font_size).set_color(color).next_to(ax.c2p(v, y_range[0]), DOWN, buff=0.12))
    for k, v in enumerate(y_ticks):
        s = y_tick_labels[k] if y_tick_labels else str(v)
        labs.add(Tex(s, font_size=font_size).set_color(color).next_to(ax.c2p(x_range[0], v), LEFT, buff=0.12))
    if x_label:
        labs.add(Tex(x_label, font_size=font_size + 4).set_color(color)
                 .next_to(ax.x_axis.get_end(), DOWN, buff=0.15))
    if y_label:
        labs.add(Tex(y_label, font_size=font_size + 4).set_color(color)
                 .next_to(ax.y_axis.get_end(), LEFT, buff=0.15))
    ax.add(labs)
    return ax


# ---------------------------------------------------------------- Pascal 三角與 Δ 格子
def pascal_triangle(rows, dx: float = 0.32, dy: float = 0.5, origin=ORIGIN, font_size: float = 24,
                    color=WHITE) -> VGroup:
    origin = np.array(origin)
    tri = VGroup()
    for N, row in enumerate(rows):
        r = VGroup()
        for y in sorted(row):
            t = Tex(str(row[y]), font_size=font_size).set_color(color)
            t.move_to(origin + y * dx * RIGHT + N * dy * DOWN)
            r.add(t)
        tri.add(r)
    return tri


def pascal_entry(tri: VGroup, N: int, y: int):
    return tri[N][(y + N) // 2]


def delta_grid(N: int, cell: float = 0.55, origin=ORIGIN, diag_color=HEAD_COLOR, off_color=GREY_D) -> VGroup:
    origin = np.array(origin)
    g = VGroup()
    for j in range(N):
        for k in range(N):
            sq = Square(side_length=cell * 0.94)
            if j == k:
                sq.set_fill(diag_color, 0.55).set_stroke(diag_color, 2)
            else:
                sq.set_fill(off_color, 0.45).set_stroke(GREY_B, 1)
            sq.move_to(origin + k * cell * RIGHT + j * cell * DOWN)
            g.add(sq)
    return g


def grid_cell(g: VGroup, j: int, k: int):
    n = int(round(len(g) ** 0.5))
    return g[j * n + k]


# ---------------------------------------------------------------- 表格與筆誤框
def tex_table(rows, col_widths, row_h: float = 0.5, font_size: float = 26, header_color=YELLOW,
              line_color=GREY_C) -> VGroup:
    """rows：list of list，元素是 TeX 字串或 Mobject（中文表頭用 zh(...)）。回傳 VGroup(格子們, 橫線)，左上角在 ORIGIN。"""
    cells = VGroup()
    x_edges = np.cumsum([0] + list(col_widths))
    for i, row in enumerate(rows):
        for j, item in enumerate(row):
            m = item if isinstance(item, Mobject) else Tex(item, font_size=font_size)
            if i == 0 and not isinstance(item, Mobject):
                m.set_color(header_color)
            m.move_to(np.array([(x_edges[j] + x_edges[j + 1]) / 2, -(i + 0.5) * row_h, 0]))
            cells.add(m)
    total_w = x_edges[-1]
    rules = VGroup(Line(ORIGIN, total_w * RIGHT).set_stroke(line_color, 1.5),
                   Line(row_h * DOWN, row_h * DOWN + total_w * RIGHT).set_stroke(line_color, 1.5),
                   Line(len(rows) * row_h * DOWN, len(rows) * row_h * DOWN + total_w * RIGHT).set_stroke(line_color, 1.5))
    return VGroup(cells, rules)


def table_cell(t: VGroup, i: int, j: int, ncols: int):
    return t[0][i * ncols + j]


def erratum_box(where: str, wrong_tex: str, right_tex: str, width: float = 5.6, font_size: float = 30) -> VGroup:
    title = zh("書上筆誤　" + where, font_size=font_size * 0.8, color=WRONG_COLOR)
    wrong = Tex(wrong_tex, font_size=font_size).set_color(WRONG_COLOR)
    strike = Line(wrong.get_left(), wrong.get_right()).set_stroke(WRONG_COLOR, 3)
    arrow = Tex(r"\Rightarrow", font_size=font_size)
    right = Tex(right_tex, font_size=font_size).set_color(RIGHT_COLOR)
    row = VGroup(VGroup(wrong, strike), arrow, right).arrange(RIGHT, buff=0.3)
    body = VGroup(title, row).arrange(DOWN, buff=0.25)
    if body.get_width() > width - 0.4:
        body.set_width(width - 0.4)
    frame = RoundedRectangle(width=max(width, body.get_width() + 0.4), height=body.get_height() + 0.45,
                             corner_radius=0.12).set_stroke(WRONG_COLOR, 2).set_fill(BLACK, 0.6)
    frame.move_to(body)
    return VGroup(frame, body)


# ---------------------------------------------------------------- 複平面
def contour_circle(axes: Axes, r: float, color=TEAL, stroke_width: float = 3) -> VGroup:
    c0 = axes.c2p(0, 0)
    rad = np.linalg.norm(axes.c2p(r, 0) - c0)
    circ = Circle(radius=rad).move_to(c0).set_stroke(color, stroke_width)
    tip_pt = c0 + rad * np.array([np.cos(PI / 4), np.sin(PI / 4), 0])
    tip = ArrowTip(angle=PI / 4 + PI / 2).set_fill(color, 1).set_height(0.18).move_to(tip_pt)
    return VGroup(circ, tip)


def pole_marker(axes: Axes, z: complex, label_tex: str = "", color=RED, font_size: float = 26) -> VGroup:
    p = axes.c2p(z.real, z.imag)
    s = 0.09
    cross = VGroup(Line(p + s * (UL), p + s * (DR)), Line(p + s * (UR), p + s * (DL))).set_stroke(color, 4)
    g = VGroup(cross)
    if label_tex:
        g.add(Tex(label_tex, font_size=font_size).set_color(color).next_to(p, UP, buff=0.15))
    return g
