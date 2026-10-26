"""
藝術桌布共用工具（walkart；combkit 的 combkit/draw/art_kit.py 是同一份）：
調色盤、背景、線條與光暈、角落文字、Plate 基底類別、輸出 PNG 與接觸表。

每張桌布是 Plate 的子類別，只要寫 art(self)；背景與文字由 construct() 自動加上。
render_png() 直接從相機讀畫面存 PNG（不用 manimgl -s；這台的 -s 存圖會壞），
呼叫前要先讓 manimlib 以「寫檔、不開視窗」初始化：import manimlib 之前把 sys.argv 設成 ["manimgl", "-w"]。
背景與暫存圖放在目前工作資料夾的 _artcache/。
"""
import hashlib
from pathlib import Path

import numpy as np
from PIL import Image

from manimlib import *

BG_DIR = Path.cwd() / "_artcache"
FW, FH = FRAME_WIDTH, FRAME_HEIGHT

# bg: (中心色, 邊緣色)；c: 主色由強到弱；ink: 文字色；light: 亮底
PAL = {
    "NOCTURNE": dict(bg=("#16204A", "#04060F"), c=["#F2C14E", "#F4EDE0", "#4FB3BF", "#E07A5F"], ink="#F4EDE0", light=False),
    "AURORA":   dict(bg=("#0E3B42", "#020D10"), c=["#7CF5C8", "#B48CFF", "#3FD0FF", "#F7F3C6"], ink="#E8FFF6", light=False),
    "EMBER":    dict(bg=("#3A1209", "#070202"), c=["#FF8C42", "#E8463A", "#FFD07A", "#FFF1D6"], ink="#FFE9D2", light=False),
    "PAPER":    dict(bg=("#FBF7EE", "#E6DCC8"), c=["#1F2A44", "#B5523B", "#7E9C83", "#C9A227"], ink="#1F2A44", light=True),
    "ORCHID":   dict(bg=("#3A1452", "#0B0312"), c=["#FF7AC8", "#C7A6FF", "#FFC2A1", "#FFF0F7"], ink="#FFEFF8", light=False),
    "GLACIER":  dict(bg=("#FFFFFF", "#D3E2EC"), c=["#0F3557", "#2F7FA8", "#7FB8D6", "#C0392B"], ink="#0F3557", light=True),
    "FOREST":   dict(bg=("#163626", "#030A06"), c=["#B7F171", "#6FBF73", "#F3E9C6", "#3FA7A0"], ink="#F3F7E6", light=False),
    "MONO":     dict(bg=("#1A1A1D", "#030304"), c=["#F2F2F5", "#A9A9B3", "#6B6B75", "#D4AF37"], ink="#ECECF0", light=False),
    "SUNSET":   dict(bg=("#6B2350", "#140A26"), c=["#FF9E7A", "#FFD27D", "#FF6FA3", "#FFF4E0"], ink="#FFF2E6", light=False),
    "OCEAN":    dict(bg=("#0A4A6E", "#010C16"), c=["#6EE7F2", "#F5D6A1", "#3FA9F5", "#E9FBFF"], ink="#E9FBFF", light=False),
}

ROMAN = ["", "I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X", "XI", "XII"]
SERIF = "Garamond"
SANS = "Corbel"


def rng(seed):
    return np.random.default_rng(seed)


def hexrgb(h):
    h = h.lstrip("#")
    return np.array([int(h[i:i + 2], 16) for i in (0, 2, 4)], float)


def lerp_hex(a, b, t):
    return rgb_to_hex(((1 - t) * hexrgb(a) + t * hexrgb(b)) / 255)


def ramp(colors, n):
    """n 個沿 colors 漸變的顏色。"""
    return color_gradient(list(colors), max(n, 2))[:n]


# ---------- 背景 ----------

def background(pal, center=(0.5, 0.5), spread=1.0, grain=0.005, seed=0, w=3840, h=2160):
    """徑向漸層＋暗角＋細微顆粒的 4K PNG（快取在 _bg/），回傳鋪滿畫面的 ImageMobject。"""
    BG_DIR.mkdir(exist_ok=True)
    c0, c1 = PAL[pal]["bg"]
    key = hashlib.md5(repr((c0, c1, center, spread, grain, seed, w, h)).encode()).hexdigest()[:10]
    path = BG_DIR / f"{pal}_{key}.png"
    if not path.exists():
        y, x = np.mgrid[0:h, 0:w].astype(np.float32)
        d = np.hypot((x / w - center[0]) * w / h, y / h - center[1]) / (0.95 * spread)
        t = np.clip(d, 0, 1) ** 1.35
        img = (1 - t)[..., None] * hexrgb(c0) + t[..., None] * hexrgb(c1)
        img += rng(seed).normal(0, 255 * grain, (h, w, 1))  # 顆粒兼抖色，避免色帶
        Image.fromarray(np.clip(img, 0, 255).astype(np.uint8)).save(path)
    return ImageMobject(str(path), height=FH).move_to(ORIGIN)


def separator():
    path = BG_DIR / "sep.png"
    if not path.exists():
        BG_DIR.mkdir(exist_ok=True)
        Image.new("RGBA", (2, 2), (0, 0, 0, 0)).save(path)
    return ImageMobject(str(path), height=0.01).move_to([-FW / 2 + 0.02, FH / 2 - 0.02, 0])


# ---------- 線與點 ----------

def to3(pts):
    pts = np.asarray(pts, float)
    if pts.shape[1] == 2:
        pts = np.hstack([pts, np.zeros((len(pts), 1))])
    return pts


def path(pts, color=WHITE, width=2.0, opacity=1.0, smooth=False):
    m = VMobject()
    if smooth:
        m.set_points_smoothly(to3(pts))
    else:
        m.set_points_as_corners(to3(pts))
    m.set_fill(opacity=0)
    m.set_stroke(color, width, opacity)
    return m


def gradient_path(pts, colors, width=2.0, opacity=1.0, smooth=False):
    m = path(pts, width=width, opacity=opacity, smooth=smooth)
    m.set_stroke(color=list(colors), width=width, opacity=opacity)
    return m


def glow(mob, color=None, width=None, layers=((7.0, 0.05), (3.5, 0.12), (1.0, 1.0))):
    """同一個形狀疊三層：粗而淡 → 細而亮。"""
    g = VGroup()
    for k, op in layers:
        c = mob.copy()
        for sm in c.family_members_with_points():
            w = width if width is not None else sm.get_stroke_width()
            sm.set_stroke(color=color, width=w * k, opacity=op * sm.get_stroke_opacity())
        g.add(c)
    return g


def _style(d, color, radius, opacity):
    """radius、color、opacity 都可以是逐點的陣列／串列。"""
    if np.ndim(radius):
        d.set_radii(np.asarray(radius, float))
    if isinstance(color, (list, tuple, np.ndarray)) or np.ndim(opacity):
        n = d.get_num_points()
        cols = list(color) if isinstance(color, (list, tuple, np.ndarray)) else [color] * n
        ops = list(np.broadcast_to(opacity, n))
        d.set_rgba_array_by_color(cols, ops)
    return d


def dots(pts, color=WHITE, radius=0.03, opacity=1.0, glow_factor=0.0):
    c0 = color if isinstance(color, str) else WHITE
    r0 = radius if not np.ndim(radius) else 0.03
    o0 = opacity if not np.ndim(opacity) else 1.0
    d = DotCloud(to3(pts), color=c0, radius=r0, opacity=o0, glow_factor=glow_factor)
    return _style(d, color, radius, opacity)


def glow_dots(pts, color=WHITE, radius=0.12, opacity=1.0, glow_factor=2.0):
    c0 = color if isinstance(color, str) else WHITE
    r0 = radius if not np.ndim(radius) else 0.12
    o0 = opacity if not np.ndim(opacity) else 1.0
    d = GlowDots(to3(pts), color=c0, radius=r0, opacity=o0, glow_factor=glow_factor)
    return _style(d, color, radius, opacity)


def fit(pts, width=None, height=None, center=ORIGIN):
    """把點縮放進 width×height 的框並置中（保持比例）。"""
    pts = np.asarray(pts, float)[:, :2]
    lo, hi = pts.min(0), pts.max(0)
    span = np.maximum(hi - lo, 1e-9)
    s = min((width or 1e9) / span[0], (height or 1e9) / span[1])
    out = (pts - (lo + hi) / 2) * s + np.asarray(center)[:2]
    return out


def contours(F, xs, ys, levels):
    """F[iy, ix] 在格點 (xs, ys) 上的等高線：回傳 [(level, (n,2) 陣列), ...]。"""
    import contourpy
    gen = contourpy.contour_generator(xs, ys, F)
    out = []
    for lv in levels:
        for seg in gen.lines(lv):
            if len(seg) > 2:
                out.append((lv, np.asarray(seg)))
    return out


def grid(f, xr, yr, n=400):
    xs = np.linspace(*xr, n)
    ys = np.linspace(*yr, int(n * (yr[1] - yr[0]) / (xr[1] - xr[0])) + 2)
    X, Y = np.meshgrid(xs, ys)
    return xs, ys, f(X, Y)


# ---------- 隨機漫步 ----------

STEPS4 = np.array([[1, 0], [-1, 0], [0, 1], [0, -1]])


def lattice_walk(n, seed):
    s = STEPS4[rng(seed).integers(0, 4, n)]
    return np.vstack([[0, 0], np.cumsum(s, 0)])


def gauss_walk(n, seed, dim=2):
    s = rng(seed).normal(size=(n, dim))
    return np.vstack([np.zeros(dim), np.cumsum(s, 0)])


def walk1d(n, seed, k=1):
    s = rng(seed).choice([-1, 1], size=(k, n))
    return np.hstack([np.zeros((k, 1)), np.cumsum(s, 1)])


D4 = [np.array(m) for m in ([[0, -1], [1, 0]], [[-1, 0], [0, -1]], [[0, 1], [-1, 0]],
                              [[1, 0], [0, -1]], [[-1, 0], [0, 1]], [[0, 1], [1, 0]], [[0, -1], [-1, 0]])]


def pivot_saw(n, seed, attempts=None):
    """方格上的自迴避漫步（樞軸法，從直棒開始）。"""
    r_ = rng(seed)
    p = np.c_[np.arange(n + 1), np.zeros(n + 1, int)]
    code = lambda q: q[:, 0] * 100003 + q[:, 1]
    for _ in range(attempts or 4 * n):
        k = int(r_.integers(1, n))
        g = D4[r_.integers(0, 7)]
        tail = p[k] + (p[k + 1:] - p[k]) @ g.T
        if not np.isin(code(tail), code(p[:k + 1])).any():
            p = np.vstack([p[:k + 1], tail])
    return p


# ---------- 桌布基底 ----------

class Plate(Scene):
    samples = 4
    pal = "NOCTURNE"
    plate = 0
    chapter = 1
    title = ""
    sub = ""          # 純文字說明（英文）
    sub_tex = ""      # 或者一行公式
    bg_center = (0.5, 0.5)
    bg_spread = 1.0
    caption_corner = DR
    panel = False     # 畫面很滿時在文字後面墊一塊柔邊底色

    @property
    def P(self):
        return PAL[self.pal]

    def C(self, i):
        return self.P["c"][i % len(self.P["c"])]

    def add(self, *mobs):
        """這台 wgpu 渲染器的毛病：DotCloud 緊接在 VMobject 後面畫會出錯，中間夾一張透明小圖就好。"""
        for m in mobs:
            if isinstance(m, DotCloud):
                super().add(separator())
            elif isinstance(m, Group) and any(isinstance(x, DotCloud) for x in m.get_family()):
                self.add(*m.submobjects)
                continue
            super().add(m)
            self.on_add(m)
        return self

    def on_add(self, m):
        """每加入一個物件（不含夾在中間的透明小圖）就呼叫一次；動畫版用它記下加入順序。"""

    def construct(self):
        self.h = {}       # art() 可以把之後要動的物件存在這裡（動畫版用）
        self.add(background(self.pal, self.bg_center, self.bg_spread, seed=self.plate))
        self.art()
        self.caption()

    def art(self):
        pass

    def caption(self):
        ink = self.P["ink"]
        tag = Text(f"PLATE {self.plate:02d}   ·   CHAPTER {ROMAN[self.chapter]}",
                   font=SANS, font_size=15)
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


# ---------- 輸出 ----------

def render_png(cls, res, out):
    """建構 Plate（或任何 Scene）並把最後畫面存成 RGB PNG。"""
    scene = cls(camera_config=dict(resolution=res),
                file_writer_config=dict(write_to_movie=False, save_last_frame=False))
    scene.setup()
    scene.construct()
    scene.update_frame(force_draw=True)
    w, h = scene.camera.get_pixel_shape()
    a = np.frombuffer(bytes(scene.camera.get_frame_bytes()), np.uint8).reshape(h, w, 4)
    Image.fromarray(a[..., :3].copy()).save(out)


def contact_sheet(files, out, cols=4, tw=480, th=270):
    """把多張圖縮小拼成一張，左上角標檔名前兩個字（編號）。"""
    from PIL import ImageDraw, ImageFont
    files = list(files)
    rows = (len(files) + cols - 1) // cols
    S = Image.new("RGB", (cols * tw, rows * th), "black")
    font = ImageFont.truetype("arial.ttf", 18)
    for i, f in enumerate(files):
        x, y = (i % cols) * tw, (i // cols) * th
        S.paste(Image.open(f).convert("RGB").resize((tw, th), Image.LANCZOS), (x, y))
        ImageDraw.Draw(S).text((x + 6, y + 4), Path(f).name[:2], fill="red", font=font)
    S.save(out)


# ---------- 動畫（約 5 秒：0–3 秒長出來，3–5 秒動態段；最後一格＝桌布） ----------

T_A, T_B = 3.0, 2.0
GROW = (0.25, 2.75)           # 「長出來」的時間窗
CAP = (2.1, 3.0)              # 文字淡入


def win(t0, t1, T, f=smooth):
    """把整段動畫的 alpha（0–1，對應 T 秒）換成只在 [t0, t1] 秒之間走完的 rate function。"""
    def rf(a):
        return f(float(np.clip((a * T - t0) / max(t1 - t0, 1e-6), 0, 1)))
    return rf


def reveal(m, rf, T):
    kw = dict(rate_func=rf, run_time=T)
    if hasattr(m, "_reveal"):
        return m._reveal(rf, T)
    if isinstance(m, ImageMobject):
        return FadeIn(m, **kw)
    if isinstance(m, DotCloud):
        return FadeIn(m, scale=4.0, **kw)          # 從 1/4 大小炸開
    subs = m.submobjects
    if len(subs) > 8 and all(hasattr(s, "key") for s in subs):   # 依 key 由小到大一顆顆亮起
        keys = np.array([s.key for s in subs], float)
        keys = (keys - keys.min()) / max(np.ptp(keys), 1e-9)
        return AnimationGroup(*[FadeIn(s, scale=2.0, rate_func=win(GROW[0] + 2.0 * k, GROW[0] + 2.0 * k + 0.5, T), run_time=T)
                                for s, k in zip(subs, keys)])
    if any(s.get_fill_opacity() > 0.02 for s in m.family_members_with_points()):
        return FadeIn(m, scale=1.25, **kw)         # 從 0.8 倍放大淡入
    return ShowCreation(m, lag_ratio=0, **kw)


def grow_points(m, keys, op, rf, T, w=0.03):
    """點雲依 keys（0–1）一顆顆亮起。"""
    keys = np.asarray(keys, float)

    def f(mob, a):
        mob.set_opacity(op * np.clip((rf(a) * (1 + w) - keys) / w, 0, 1))
    return UpdateFromAlphaFunc(m, f, run_time=T)


def along(P, u):
    i = u * (len(P) - 1)
    j = min(int(i), len(P) - 2)
    return P[j] + (P[j + 1] - P[j]) * (i - j)


def travelers(scene, polys, color, T, t0=0.0, dur=None, stagger=0.0, rf=None, radius=0.14, glow_factor=2.5):
    """幾顆光點各沿一條折線走（依點的編號均勻走，和 ShowCreation 同步）；還沒出發或已到的點隱藏。"""
    polys = [np.asarray(P, float)[:, :2] for P in polys]
    dur = dur or T - t0
    gd = glow_dots([P[0] for P in polys], color, radius=radius, glow_factor=glow_factor, opacity=0.0)
    scene.add(gd)

    def f(mob, a):
        pts, ops = [], []
        for i, P in enumerate(polys):
            u = rf(a) if rf else (a * T - t0 - i * stagger) / dur
            pts.append(along(P, float(np.clip(u, 0, 1))))
            ops.append(1.0 if 0 < u < 1 or (rf and u > 0) else 0.0)
        mob.set_points(to3(np.array(pts)))
        mob.set_opacity(np.array(ops))
    return UpdateFromAlphaFunc(gd, f, run_time=T)


def keep(m, fn, T):
    """記下 m 一家每個物件的點，之後每格用 fn(原點, alpha) 重設。"""
    fam = [s for s in m.get_family() if s.get_num_points()]
    orig = [s.get_points().copy() for s in fam]

    def f(mob, a):
        for s, p in zip(fam, orig):
            s.set_points(fn(p, a))
    return UpdateFromAlphaFunc(m, f, run_time=T)


def spin(items, angle=None, T=T_B):
    out = []
    for it in items:
        m, c = it[0], it[1]
        ang = it[2] if len(it) > 2 and angle is None else angle
        out.append(Rotate(m, ang, about_point=c, run_time=T, rate_func=smooth))
    return out


def mirror(items, axis, T=T_B):
    """每個 (物件, 中心) 沿平面上的 axis 方向翻 180° 再翻回來（3D 翻轉＝鏡射；youngart 20、30）。"""
    ax = np.array([*axis, 0.0])
    return [Rotate(m, PI, axis=ax, about_point=np.array([*c, 0.0]), run_time=T, rate_func=there_and_back)
            for m, c in items]


def twinkle(m, seed=0, amp=0.75, T=T_B):
    """點雲逐點閃爍；頭尾不變，接在「長出來」後面不會跳。"""
    base = m.get_opacities().copy()
    ph, fr = rng(seed).uniform(0, TAU, len(base)), rng(seed + 1).uniform(1.5, 4.0, len(base))

    def f(mob, a):
        mob.set_opacity(base * (1 - amp * np.sin(PI * a) * (0.5 - 0.5 * np.cos(fr * TAU * a + ph))))
    return UpdateFromAlphaFunc(m, f, run_time=T)


def pulse(m, seed=0, amp=0.45, T=T_B):
    """點雲逐點呼吸（半徑）。"""
    r0 = m.get_radii()[:, 0].copy()
    ph = rng(seed).uniform(0, TAU, len(r0))

    def f(mob, a):
        mob.set_radii(r0 * (1 + amp * np.sin(PI * a) * (0.5 - 0.5 * np.cos(2 * TAU * a + ph))))
    return UpdateFromAlphaFunc(m, f, run_time=T)


def keyed_reveal(m, keys, w=0.03):
    """把 m 的「長出來」換成依 key（0–1）一個個出現：keys 是陣列時對 m 的子物件，純量時對 m 本身。
    一個 updater、只改有變的子物件（2500 格也快），保留每個子物件原本的填色／筆畫透明度。用在 animate_plate 的 prep。"""
    subs = m.submobjects if np.ndim(keys) else [m]
    k = np.atleast_1d(np.asarray(keys, float))
    fo = np.array([x.get_fill_opacity() for x in subs])
    so = np.array([x.get_stroke_opacity() for x in subs])

    def make(rf, T):
        last = np.full(len(subs), -1.0)

        def f(mob, a):
            x = np.clip((rf(a) * (1 + w) - k) / w, 0, 1)
            for i in np.nonzero(np.abs(x - last) > 1e-3)[0]:
                subs[i].set_fill(opacity=fo[i] * x[i]).set_stroke(opacity=so[i] * x[i])
                last[i] = x[i]
        return UpdateFromAlphaFunc(m, f, run_time=T)
    m._reveal = make


def animate_plate(base, cfg=None):
    """把 Plate 子類別包成動畫版（最後一格＝桌布）。cfg：together、A(scene)、B(scene)、drift、
    prep(scene)（art() 之後立刻執行，可用 keyed_reveal 換掉某些物件的「長出來」）。"""
    cfg = cfg or {}

    class Anim(base):
        def on_add(self, m):
            ph = getattr(self, "_phase", None)
            if ph:
                self.rec[ph].append(m)

        def construct(self):
            self.h, self.rec = {}, {"bg": [], "art": [], "cap": []}
            self._phase = "bg"
            self.add(background(self.pal, self.bg_center, self.bg_spread, seed=self.plate))
            self._phase = "art"
            self.art()
            if "prep" in cfg:
                cfg["prep"](self)
            self._phase = "cap"
            self.caption()
            self._phase = None
            for m in self.rec["cap"]:
                m.fix_in_frame()
            for mob, keys, op in self.h.get("grow", []):
                mob._reveal = lambda rf, T, mob=mob, keys=keys, op=op: grow_points(mob, keys, op, rf, T)
            art = self.rec["art"]
            n = len(art)
            anims = [FadeIn(m, run_time=T_A, rate_func=win(0, 0.4, T_A)) for m in self.rec["bg"]]
            for i, m in enumerate(art):
                if cfg.get("together") and isinstance(m, DotCloud) and not hasattr(m, "_reveal"):
                    t0, t1 = (0.1, 0.6) if isinstance(m, GlowDots) and m.get_num_points() <= 3 else (2.0, 2.9)
                elif cfg.get("together") or n == 1:
                    t0, t1 = GROW
                else:
                    t0 = GROW[0] + 1.7 * i / (n - 1)
                    t1 = min(t0 + 1.0, T_A - 0.05)
                anims.append(reveal(m, win(t0, t1, T_A), T_A))
            anims += [FadeIn(m, shift=0.1 * UP, run_time=T_A, rate_func=win(*CAP, T_A)) for m in self.rec["cap"]]
            if "A" in cfg:
                anims += cfg["A"](self)
            if cfg.get("drift", True):
                self.frame.add_updater(lambda f, dt: f.scale(1 - 0.012 * dt))
            self.play(AnimationGroup(*anims, lag_ratio=0))
            B = cfg["B"](self) if "B" in cfg else []
            for m in art:  # 焦點光點（≤3 顆的 GlowDots）呼吸一下
                if isinstance(m, GlowDots) and m.get_num_points() <= 3:
                    r0 = m.get_radii().copy()
                    B.append(UpdateFromAlphaFunc(m, lambda mob, a, r0=r0: mob.set_radii(r0 * (1 + 0.35 * there_and_back(a))),
                                                 run_time=T_B))
            if B:
                self.play(AnimationGroup(*B, lag_ratio=0))
            else:
                self.wait(T_B)

    Anim.__name__ = Anim.__qualname__ = f"A{base.plate:02d}_{base.__name__}"
    return Anim



def render_video(cls, res, fps, outdir, stem, crf=16):
    """一個行程裡直接渲染成 mp4（不用 manimgl 指令）；呼叫前同樣要讓 manimlib 以「寫檔」初始化。"""
    Path(outdir).mkdir(parents=True, exist_ok=True)
    scene = cls(camera_config=dict(resolution=res, fps=fps),
                file_writer_config=dict(write_to_movie=True, save_last_frame=False, output_directory=str(outdir),
                                        file_name=stem, crf=crf, quiet=True))
    scene.run()
