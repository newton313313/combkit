# 由 combkit/tools/sync_from_symosc.py 從 symosc 專案的 codes/utils/osc_draw.py 複製；請改 symosc 那份再同步。
"""
osc_draw.py —— 對稱群與耦合振盪子的畫圖元件（ManimGL）。數學在 symosc.py（先算、再畫）。

元件都用參數決定大小，不要事後 .scale()。注意：Tex 的 color= 無效，一律 .set_color()。

    tri_points(radius, center)          正三角形三個頂點（擺 1 在上，逆時針 1→2→3；轉動 r 是逆時針 120°）
    spring(p, q, coils, amp, color)     兩點之間的鋸齒彈簧
    osc_triangle(radius, center, labels=True)   三個錨點＋三根彈簧＋編號 → VGroup(springs, anchors, labels)
    masses_updater(tracker_fn, radius, center, amp, color)
                                        always_redraw 的質點＋彈簧：tracker_fn() 回傳三個位移（沿徑向）
    radial_bars(values, radius, center, scale, color)   每個頂點一根徑向的長條（正值向外、負值向內）
    phasor(z, at, scale, color)         一個複數畫成從 at 出發的箭頭（長度 |z|·scale、角度 arg z）；|z| 很小時畫小點
    phasors(zs, points, scale, color, zero_dot)   一組相量（zero_dot=False：零分量不畫）
    level_diagram(levels, x, y0, yscale, width, color)
                                        levels = [(值, 標籤 TeX, 簡併數)] → 水平線＋右側標籤（簡併數 > 1 時寫 ×n）
    q8_positions(r_in, r_out, center)   Q8 的 8 個元素的位置：內圈 1, i, −1, −i（×i 一圈），外圈 j, −k, −j, k
    q8_graph(pos, colors)               Cayley 圖：h 與 h·s 相連（s = ±i 藍、±j 綠、±k 橘）→ VGroup(edges_i, edges_j, edges_k, nodes, labels)
"""
import math

import numpy as np
from manimlib import *

from combkit.draw.styles import zh

C_I = "#4A90D9"
C_J = "#7CE38B"
C_K = "#F2A33A"
MASS_COLOR = YELLOW
SPRING_COLOR = GREY_B


def tri_points(radius: float = 1.5, center=ORIGIN):
    c = np.array(center, dtype=float)
    return [c + radius * np.array([math.cos(a), math.sin(a), 0]) for a in (math.pi / 2, math.pi / 2 + 2 * math.pi / 3,
                                                                           math.pi / 2 + 4 * math.pi / 3)]


def spring(p, q, coils: int = 7, amp: float = 0.08, color=SPRING_COLOR, stroke_width: float = 2.5):
    p, q = np.array(p, dtype=float), np.array(q, dtype=float)
    d = q - p
    L = np.linalg.norm(d)
    if L < 1e-6:
        return VMobject()
    u = d / L
    nrm = np.array([-u[1], u[0], 0])
    pts = [p, p + 0.12 * d]
    n = 2 * coils
    for s in range(1, n):
        t = 0.12 + 0.76 * s / n
        pts.append(p + t * d + amp * (1 if s % 2 else -1) * nrm)
    pts += [p + 0.88 * d, q]
    return VMobject().set_points_as_corners(pts).set_stroke(color, stroke_width)


def osc_triangle(radius: float = 1.5, center=ORIGIN, labels: bool = True, font_size: float = 26):
    P = tri_points(radius, center)
    springs = VGroup(*[spring(P[a], P[b]) for a, b in ((0, 1), (1, 2), (2, 0))])
    anchors = VGroup(*[Circle(radius=0.07).set_stroke(GREY_C, 2).move_to(p) for p in P])
    labs = VGroup()
    if labels:
        c = np.array(center, dtype=float)
        for j, p in enumerate(P):
            u = (p - c) / np.linalg.norm(p - c)
            tang = np.array([-u[1], u[0], 0])                       # 標號放在切線方向，不擋住沿徑向擺動的質點
            labs.add(Tex(str(j + 1), font_size=font_size).set_color(GREY_A).move_to(p + 0.4 * tang))
    return VGroup(springs, anchors, labs)


def masses_updater(tracker_fn, radius: float = 1.5, center=ORIGIN, amp: float = 0.45, color=MASS_COLOR):
    P = tri_points(radius, center)
    c = np.array(center, dtype=float)
    U = [(p - c) / np.linalg.norm(p - c) for p in P]

    def build():
        x = tracker_fn()
        Q = [p + amp * xi * u for p, xi, u in zip(P, x, U)]
        sp = VGroup(*[spring(Q[a], Q[b]) for a, b in ((0, 1), (1, 2), (2, 0))])
        dots = VGroup(*[Dot(q, radius=0.13).set_fill(color, 1).set_stroke(WHITE, 1.5) for q in Q])
        return VGroup(sp, dots)

    return always_redraw(build)


def radial_bars(values, radius: float = 1.5, center=ORIGIN, scale: float = 0.35, color=MASS_COLOR, width: float = 7):
    P = tri_points(radius, center)
    c = np.array(center, dtype=float)
    out = VGroup()
    for p, v in zip(P, values):
        u = (p - c) / np.linalg.norm(p - c)
        if abs(v) < 1e-9:
            out.add(Dot(p, radius=0.06).set_fill(GREY_B, 1))
        else:
            out.add(Line(p, p + scale * v * u).set_stroke(color if v > 0 else BLUE_C, width))
    return out


def phasor(z: complex, at, scale: float = 0.45, color=YELLOW, zero_dot: bool = True):
    at = np.array(at, dtype=float)
    if abs(z) < 1e-9:
        return Dot(at, radius=0.05).set_fill(GREY_B, 1 if zero_dot else 0)
    end = at + scale * np.array([z.real, z.imag, 0])
    return Arrow(at, end, buff=0, thickness=3.0).set_color(color)


def phasors(zs, points, scale: float = 0.45, color=YELLOW, zero_dot: bool = True):
    return VGroup(*[phasor(z, p, scale, color, zero_dot) for z, p in zip(zs, points)])


def level_diagram(levels, x: float = 0.0, y0: float = -2.0, yscale: float = 1.0, width: float = 1.2, color=WHITE,
                  font_size: float = 26, label_side=RIGHT):
    out = VGroup()
    for val, tex, deg in levels:
        y = y0 + yscale * val
        line = Line([x - width / 2, y, 0], [x + width / 2, y, 0]).set_stroke(color, 4)
        s = tex + (rf"\ \times{deg}" if deg > 1 else "")
        lab = Tex(s, font_size=font_size).next_to(line, label_side, buff=0.15)
        out.add(VGroup(line, lab))
    return out


def q8_positions(r_in: float = 0.9, r_out: float = 1.9, center=ORIGIN):
    c = np.array(center, dtype=float)
    inner = ["1", "i", "-1", "-i"]
    outer = ["j", "-k", "-j", "k"]
    pos = {}
    for t, n in enumerate(inner):
        a = math.pi / 2 - t * math.pi / 2
        pos[n] = c + r_in * np.array([math.cos(a), math.sin(a), 0])
    for t, n in enumerate(outer):
        a = math.pi / 2 - t * math.pi / 2 + math.pi / 4
        pos[n] = c + r_out * np.array([math.cos(a), math.sin(a), 0])
    return pos


def q8_graph(pos: dict, mul, font_size: float = 26, node_radius: float = 0.22):
    """mul(a, b) 是 Q8 的乘法（symosc.q8_mul）。"""
    groups = []
    for u, col in (("i", C_I), ("j", C_J), ("k", C_K)):
        seen, g = set(), VGroup()
        for h in pos:
            for s in (u, "-" + u):
                e = frozenset((h, mul(h, s)))
                if e in seen:
                    continue
                seen.add(e)
                a, b = tuple(e)
                g.add(Line(pos[a], pos[b], buff=node_radius).set_stroke(col, 2.5, opacity=0.8))
        groups.append(g)
    nodes = VGroup(*[Circle(radius=node_radius).set_fill(BLACK, 1).set_stroke(WHITE, 2).move_to(p) for p in pos.values()])
    labels = VGroup(*[Tex(n.replace("-", "-"), font_size=font_size).move_to(p) for n, p in pos.items()])
    return VGroup(*groups, nodes, labels)
