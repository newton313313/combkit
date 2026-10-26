"""
young_art.py —— Young 表的藝術桌布工具（來自 youngart：Fulton《Young Tableaux》全書 60 張 4K 桌布，codes/art/ykit.py）

從 art_kit 匯入一切（Plate、PAL、background、path、glow、dots、glow_dots、ramp、rng、render_png、contact_sheet…），
也匯入 combkit.young_extras 的計算（patience、slide_path、shadow_lines、rsk_rows、pipe_route、meet…），再加上：

    YPlate                 角落標籤 PLATE NN · CHAPTER I…X／APPENDIX A、B（chapter = 11、12 是附錄；改 PARTS 可換成別本書）
    cell_center(i, j, cs, origin)                英式座標（第 0 列在上）的格子中心
    cell_squares(cells, cs, origin, colors, opacity, stroke, sw, gap)   方格：cells = [(i, j), …]，colors 單色或逐格
    tab_squares(T, cs, origin, cmap, vmax, …)    表 T 依數字上色（cmap 內插成 vmax−vmin+1 色）
    outline_pts(lam, cs, origin)                 分割 λ 的外框（從左上角順時針，封閉）；[1:-1] 只剩外側的階梯
    perm_pts(w, size, center)                    排列 w 的點（第 i 個點在第 w(i) 高；y 往上）
    rot3(P, yaw, pitch) / proj(P, yaw, pitch, scale, center)   3D 旋轉與正交投影（回傳 2D 座標與深度，越大越近）
    shape_of(T)                                  表的形狀

注意：ramp() 回傳的是顏色物件，要給 lerp_hex 時先用 rgb_to_hex(color_to_rgb(c)) 轉成字串；
空的點串列不要交給 dots／glow_dots。
"""
from combkit.draw.art_kit import *  # noqa: F401,F403
from combkit.draw import art_kit as _ak
from combkit.young_extras import *  # noqa: F401,F403

PARTS = {**{i: f"CHAPTER {_ak.ROMAN[i]}" for i in range(1, 11)}, 11: "APPENDIX A", 12: "APPENDIX B"}


class YPlate(Plate):
    """與 art_kit.Plate 相同，只有右上角的標籤用 PARTS（章或附錄）。"""

    def caption(self):
        ink = self.P["ink"]
        tag = Text(f"PLATE {self.plate:02d}   ·   {PARTS[self.chapter]}", font=SANS, font_size=15)
        tag.set_fill(ink, 0.55)
        tag.to_corner(UR, buff=0.42)
        self.add(tag)
        if not self.title:
            return
        title = Text(self.title, font=SERIF, font_size=40, slant=ITALIC)
        title.set_fill(ink, 0.95)
        parts = [title]
        if self.sub_tex:
            s = Tex(self.sub_tex, font_size=26)
            s.set_fill(ink, 0.75)
            parts.append(s)
        elif self.sub:
            s = Text(self.sub, font=SANS, font_size=19)
            s.set_fill(ink, 0.7)
            parts.append(s)
        al = RIGHT if self.caption_corner[0] > 0 else LEFT
        block = VGroup(*parts).arrange(DOWN, buff=0.16, aligned_edge=al)
        rule = Line(LEFT, RIGHT).set_width(title.get_width())
        rule.set_stroke([self.C(0), self.C(2)], 1.2, 0.8)
        rule.next_to(title, UP, buff=0.14, aligned_edge=al)
        cap = VGroup(rule, block)
        cap.to_corner(self.caption_corner, buff=0.5)
        if self.panel:
            under = self.P["bg"][0 if self.P["light"] else 1]
            for b, op in ((0.5, 0.2), (0.34, 0.3), (0.18, 0.45)):
                self.add(RoundedRectangle(width=cap.get_width() + 2 * b, height=cap.get_height() + 2 * b,
                                          corner_radius=b).move_to(cap).set_fill(under, op).set_stroke(width=0))
        self.add(cap)


def shape_of(T):
    return tuple(len(r) for r in T)


def cell_center(i, j, cs, origin):
    return np.array([origin[0] + (j + 0.5) * cs, origin[1] - (i + 0.5) * cs, 0.0])


def cell_squares(cells, cs=0.3, origin=(0, 0), colors=WHITE, opacity=0.85, stroke=None, sw=1.0, gap=0.08,
                 stroke_opacity=0.9):
    g = VGroup()
    cols = colors if isinstance(colors, (list, tuple)) else [colors] * len(cells)
    for (i, j), c in zip(cells, cols):
        sq = Square(cs * (1 - gap)).move_to(cell_center(i, j, cs, origin))
        sq.set_fill(c, opacity)
        sq.set_stroke(stroke or c, sw, stroke_opacity if stroke else 0)
        g.add(sq)
    return g


def tab_squares(T, cs=0.3, origin=(0, 0), cmap=None, vmax=None, vmin=1, **kw):
    cells = [(i, j) for i, r in enumerate(T) for j in range(len(r))]
    vals = [T[i][j] for i, j in cells]
    vmax = vmax or max(vals)
    pal = ramp(cmap or [WHITE, GREY], vmax - vmin + 1)
    return cell_squares(cells, cs, origin, [pal[v - vmin] for v in vals], **kw)


def outline_pts(lam, cs=0.3, origin=(0, 0)):
    lam = [x for x in lam if x]
    x0, y0 = origin
    pts = [(x0, y0), (x0 + lam[0] * cs, y0)]
    for i, l in enumerate(lam):
        y = y0 - (i + 1) * cs
        pts.append((x0 + l * cs, y))
        nxt = lam[i + 1] if i + 1 < len(lam) else 0
        pts.append((x0 + nxt * cs, y))
    pts.append((x0, y0))
    return np.array(pts, float)


def perm_pts(w, size=1.0, center=ORIGIN):
    n = len(w)
    s = size / max(n - 1, 1)
    return np.array([[center[0] + (i - (n - 1) / 2) * s, center[1] + (w[i] - 1 - (n - 1) / 2) * s] for i in range(n)])


def rot3(P, yaw=0.6, pitch=0.35):
    P = np.asarray(P, float)
    cy, sy, cp, sp = np.cos(yaw), np.sin(yaw), np.cos(pitch), np.sin(pitch)
    Rz = np.array([[cy, -sy, 0], [sy, cy, 0], [0, 0, 1]])
    Rx = np.array([[1, 0, 0], [0, cp, -sp], [0, sp, cp]])
    return P @ Rz.T @ Rx.T


def proj(P, yaw=0.6, pitch=0.35, scale=1.0, center=ORIGIN):
    Q = rot3(P, yaw, pitch)
    xy = np.c_[Q[:, 0], Q[:, 2]] * scale + np.asarray(center)[:2]
    return xy, -Q[:, 1]


def turn3d(items, view, angle=0.7, T=T_B):
    """動畫：items = [(物件, 3D 折線), …]（物件是用 proj 畫出的 path／glow 群組或點雲），view = (yaw, pitch, scale, center)；
    每格把 yaw 加上 angle·smooth(α) 重新投影（youngart 42、48 的動態段）。"""
    yaw, pitch, sc, ctr = view

    def f(mob, a):
        for m, P3 in items:
            xy = to3(proj(P3, yaw + angle * smooth(a), pitch, sc, ctr)[0])
            if isinstance(m, DotCloud):
                m.set_points(xy)
            else:
                for sm in m.family_members_with_points():
                    sm.set_points_as_corners(xy)
    return UpdateFromAlphaFunc(items[0][0], f, run_time=T)
