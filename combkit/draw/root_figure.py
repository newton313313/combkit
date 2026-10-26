"""
root_figure.py —— Scene 05 之後反覆用到的「根分解」示意圖，以及群組小圖示

示意圖內容（全部手動指定，與講義的圖 5 相同的結構）：
    12 個圓周點；固定點是 0
    根群組：三角形 (0, 4, 9)
    R_1 = {1,2,3}   -> 三角形 (1,2,3)           r_1 = (0, 1)
    R_2 = {5,6,7,8} -> 兩條弦 (5,8), (6,7)      r_2 = (2, 0)
    R_3 = {10,11}   -> 一條弦 (10,11)           r_3 = (1, 0)
對應 v = (2, 3)、n = (3, 2)。
"""
from manimlib import *

from combkit.draw.polygon_groups import make_configuration, root_regions
from combkit.draw.styles import DOT_RADIUS, REGION_COLORS, block_color

ROOT_N_POINTS = 12
ROOT_BLOCK = (0, 4, 9)
ROOT_SUB_BLOCKS = [
    [(1, 2, 3)],
    [(5, 8), (6, 7)],
    [(10, 11)],
]
ROOT_ALT_R2 = [(5, 6), (7, 8)]       # R_2 的另一種分法
ROOT_R_VECTORS = [(0, 1), (2, 0), (1, 0)]   # 各段的群組清單 r_1, r_2, r_3
ROOT_N_VEC = (3, 2)
ROOT_V_VEC = (2, 3)


class RootFigure(VGroup):
    """
    屬性：
        .regions        VGroup，三個弓形區域
        .config         Configuration（圓 + 群組 + 點）
        .root           根群組（三角形）
        .subs           list[VGroup]，各段裡的子群組
        .region_labels  VGroup，R_1, R_2, R_3
    """

    def __init__(self, center=ORIGIN, radius: float = 1.9, label_font_size: float = 40,
                 with_subs: bool = True, label_buff: float = 0.42, region_opacity: float = 0.30):
        center = np.array(center, dtype=float)
        blocks = [ROOT_BLOCK]
        if with_subs:
            for group in ROOT_SUB_BLOCKS:
                blocks.extend(group)
        config = make_configuration(ROOT_N_POINTS, blocks, radius=radius, center=center, root=0)
        root = config.blocks[0]
        root.set_stroke(width=root.get_stroke_width() * 1.35)

        subs = []
        if with_subs:
            k = 1
            for group in ROOT_SUB_BLOCKS:
                subs.append(VGroup(*config.blocks[k:k + len(group)]))
                k += len(group)

        regions = root_regions(ROOT_BLOCK, ROOT_N_POINTS, radius, center, opacity=region_opacity)
        regions.set_z_index(-1)
        region_labels = VGroup(*[
            Tex(f"R_{j + 1}", font_size=label_font_size, fill_color=REGION_COLORS[j]).move_to(
                center + (radius + label_buff) * normalize(r.arc_mid - center))
            for j, r in enumerate(regions)
        ])

        super().__init__(regions, config, region_labels)
        self.regions = regions
        self.config = config
        self.root = root
        self.subs = subs
        self.region_labels = region_labels
        self.coords = config.coords
        self.dots = config.dots
        self.center_point = center
        self.radius = radius


def block_icon(size: int, radius: float = 0.28, color=None, stroke_width: float = 4,
               dot_radius: float = None) -> VGroup:
    """
    不放在圓上的「群組圖示」：正 size 邊形（size=2 時是一條線段）加上頂點。
    用來表示「一個 size 點群組」。
    """
    if color is None:
        color = block_color(size)
    if dot_radius is None:
        dot_radius = DOT_RADIUS * 0.75
    if size == 2:
        pts = [radius * LEFT + 0.35 * radius * DOWN, radius * RIGHT + 0.35 * radius * UP]
        shape = Line(pts[0], pts[1], stroke_color=color, stroke_width=stroke_width)
    else:
        pts = [radius * np.array([np.cos(a), np.sin(a), 0.0])
               for a in (PI / 2 - TAU * j / size for j in range(size))]
        shape = Polygon(*pts, stroke_color=color, stroke_width=stroke_width,
                        fill_color=color, fill_opacity=0.3)
    dots = VGroup(*[Dot(p, radius=dot_radius, fill_color=color) for p in pts])
    return VGroup(shape, dots)
