"""
polygon_groups.py —— 群組（block）的畫法，以及「一整個分法」的組裝

一個群組就是一串點的編號：
    2 個點  -> 一條弦
    3 個以上 -> 凸多邊形（半透明填色）
所有分法都用「明確寫出的點編號」指定，不使用亂數，確保每次渲染結果相同。
"""
from manimlib import *

from combkit.draw.circle_points import make_circle_points, point_angle
from combkit.counting import blocks_cross, is_noncrossing  # noqa: F401  （給場景直接使用）
from combkit.draw.styles import (
    BLOCK_FILL_OPACITY,
    BLOCK_STROKE,
    CIRCLE_STROKE,
    CROSSING_COLOR,
    DOT_RADIUS,
    ROOT_COLOR,
    block_color,
)


def make_block(indices, coords, color=None, fill_opacity: float = BLOCK_FILL_OPACITY,
               stroke_width: float = BLOCK_STROKE) -> VMobject:
    """
    依點編號畫出一個群組。
        1 個點 -> 小圓點
        2 個點 -> 弦
        3 個以上 -> 多邊形
    顏色預設由點數決定（見 styles.block_color）。
    """
    indices = list(indices)
    if color is None:
        color = block_color(len(indices))
    pts = [coords[i] for i in indices]
    if len(pts) == 1:
        block = Dot(pts[0], radius=0.11, fill_color=color)
    elif len(pts) == 2:
        block = Line(pts[0], pts[1], stroke_color=color, stroke_width=stroke_width)
    else:
        block = Polygon(*pts, stroke_color=color, stroke_width=stroke_width,
                        fill_color=color, fill_opacity=fill_opacity)
    block.indices = tuple(indices)
    return block


def make_blocks(blocks, coords, **kwargs) -> VGroup:
    """一次畫出多個群組。blocks 是「點編號串列」的串列。"""
    return VGroup(*[make_block(b, coords, **kwargs) for b in blocks])


class Configuration(VGroup):
    """
    一個完整的分法：圓 + 群組 + 圓周點。
    屬性：
        .circle, .blocks, .dots : 對應的 mobject
        .coords                 : 各點座標（建立當下的值）
        .n_points, .block_indices
    子物件順序固定為 [circle, blocks, dots]，所以點永遠畫在最上層。
    """

    def __init__(self, n_points: int, blocks, radius: float = 2.5, center=ORIGIN,
                 dot_radius: float = None, stroke_width: float = None,
                 circle_stroke: float = None, fill_opacity: float = BLOCK_FILL_OPACITY,
                 root: int = None, colors=None, color_dots: bool = True):
        # 小圖自動用細一點的線、小一點的點
        scale = min(1.0, radius / 2.5)
        if dot_radius is None:
            dot_radius = max(0.035, DOT_RADIUS * (0.45 + 0.55 * scale))
        if stroke_width is None:
            stroke_width = max(2.2, BLOCK_STROKE * (0.35 + 0.65 * scale))
        if circle_stroke is None:
            circle_stroke = max(1.5, CIRCLE_STROKE * (0.5 + 0.5 * scale))

        circle, dots, coords = make_circle_points(
            n_points, radius=radius, center=center,
            dot_radius=dot_radius, circle_stroke=circle_stroke,
        )
        block_mobs = VGroup()
        for k, b in enumerate(blocks):
            color = colors[k] if colors else None
            block_mobs.add(make_block(b, coords, color=color,
                                      fill_opacity=fill_opacity, stroke_width=stroke_width))
        if color_dots:
            # 規則：同一群組的點與它的多邊形同色
            for block in block_mobs:
                for i in block.indices:
                    dots[i].set_fill(block.get_stroke_color())
        if root is not None:
            dots[root].set_fill(ROOT_COLOR)
            dots[root].scale(1.6)

        super().__init__(circle, block_mobs, dots)
        self.circle = circle
        self.blocks = block_mobs
        self.dots = dots
        self.coords = coords
        self.n_points = n_points
        self.block_indices = [tuple(b) for b in blocks]
        self.radius = radius


def make_configuration(n_points: int, blocks, **kwargs) -> Configuration:
    """建立一個分法（圓 + 群組 + 點）。blocks 例如 [(0, 1, 2), (3, 4)]。"""
    return Configuration(n_points, blocks, **kwargs)


def dot_color_anims(dots, indices, color) -> list:
    """把指定的點變成某個顏色（回傳動畫串列，放進 self.play）。"""
    return [dots[i].animate.set_fill(color) for i in indices]


def mark_crossing(block: VMobject) -> VMobject:
    """把一個群組塗成「交叉」的紅色。"""
    block.set_stroke(CROSSING_COLOR)
    if block.get_fill_opacity() > 0:
        block.set_fill(CROSSING_COLOR)
    return block


def highlight_block(block: VMobject, scale_factor: float = 1.08, color=None, **kwargs):
    """讓一個群組閃一下（放大再縮回）。回傳動畫。"""
    if color is None:
        color = block.get_stroke_color()
    return Indicate(block, scale_factor=scale_factor, color=color, **kwargs)


def circular_segment(i: int, j: int, n: int, radius: float, center=ORIGIN,
                     color=TEAL, opacity: float = 0.22, n_samples: int = 48) -> VMobject:
    """
    弓形區域：從第 i 個點沿圓周順時針走到第 j 個點的圓弧，加上弦 i-j 圍成的區域。
    （j <= i 時表示繞過起點，例如 i=9, j=0 on n=12 代表 9 -> 10 -> 11 -> 0。）
    這正是「根群組相鄰兩頂點之間」那一段點所在的區域。
    """
    center = np.array(center, dtype=float)
    steps = (j - i) % n
    if steps == 0:
        steps = n
    a0 = point_angle(i, n)
    a1 = a0 - TAU * steps / n
    angles = np.linspace(a0, a1, max(4, int(n_samples * steps / n) + 2))
    pts = [center + radius * np.array([np.cos(a), np.sin(a), 0.0]) for a in angles]
    region = Polygon(*pts, stroke_width=0, fill_color=color, fill_opacity=opacity)
    region.arc_mid = center + radius * np.array([
        np.cos((a0 + a1) / 2), np.sin((a0 + a1) / 2), 0.0
    ])
    region.steps = steps
    return region


def root_regions(root_block, n: int, radius: float, center=ORIGIN, colors=None,
                 opacity: float = 0.22) -> VGroup:
    """
    根群組（依順時針排序的點編號）切出的各個區域。
    回傳 VGroup，第 j 個元素是根群組第 j 與第 j+1 個頂點之間的弓形。
    """
    from combkit.draw.styles import REGION_COLORS
    colors = colors or REGION_COLORS
    root_block = list(root_block)
    regions = VGroup()
    for k in range(len(root_block)):
        a, b = root_block[k], root_block[(k + 1) % len(root_block)]
        regions.add(circular_segment(a, b, n, radius, center,
                                     color=colors[k % len(colors)], opacity=opacity))
    return regions
