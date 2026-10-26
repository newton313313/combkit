"""
lattice_paths.py —— Raney 公式證明用的「路徑」圖形

一條路徑是一串「步」（見 counting.py 第 4 節）：
    ↑ 步       ：高度 +1，灰色
    第 i 型步  ：高度 1 - v_i（弦 -1、三角形 -2），顏色與第 i 型群組相同
第 j 個點畫在 (j, S_j)，每一步往右一格。
"""
from manimlib import *

from combkit.counting import LEAF, first_bad_index, path_heights, step_weight
from combkit.draw.styles import CROSSING_COLOR, block_color

UP_STEP_COLOR = GREY_A


def letter_color(letter: int, v_vec):
    return UP_STEP_COLOR if letter == LEAF else block_color(v_vec[letter - 1])


def weight_tex(w: int) -> str:
    return f"+{w}" if w > 0 else str(w)


def word_tiles(word, v_vec, size: float = 0.5, buff: float = 0.08, font_size: float = 26) -> VGroup:
    """把一串步畫成一排方塊，方塊上寫高度變化（+1、-1、-2）。"""
    tiles = VGroup()
    for letter in word:
        color = letter_color(letter, v_vec)
        square = RoundedRectangle(width=size, height=size, corner_radius=0.08,
                                  stroke_color=color, stroke_width=2.5,
                                  fill_color=color, fill_opacity=0.22)
        label = Tex(weight_tex(step_weight(letter, v_vec)), font_size=font_size, fill_color=color)
        label.move_to(square)
        tile = VGroup(square, label)
        tile.letter = letter
        tiles.add(tile)
    tiles.arrange(RIGHT, buff=buff)
    return tiles


class PathMob(VGroup):
    """
    一條路徑。origin 是「第 0 個點」的位置，start_height 是它的高度（只影響 .heights）。
    屬性：
        .steps    VGroup，每一步一條線段
        .dots     VGroup，每個點一個圓點（可關掉）
        .heights  各點的高度
    方法：
        .vertex(j)  第 j 個點目前的座標（物件移動後也正確）
        .y_of(h)    高度 h 目前在畫面上的 y 座標
    """

    def __init__(self, word, v_vec, unit: float = 0.4, origin=ORIGIN, start_height: int = 0,
                 stroke_width: float = 5, show_dots: bool = True, dot_radius: float = 0.05):
        origin = np.array(origin, dtype=float)
        hs = path_heights(word, v_vec, start=start_height)
        pts = [origin + unit * (j * RIGHT + (h - start_height) * UP) for j, h in enumerate(hs)]
        steps = VGroup(*[
            Line(pts[j], pts[j + 1], stroke_color=letter_color(word[j], v_vec), stroke_width=stroke_width)
            for j in range(len(word))
        ])
        dots = VGroup()
        if show_dots:
            dots = VGroup(*[Dot(p, radius=dot_radius, fill_color=WHITE) for p in pts])
            dots.set_z_index(2)
        super().__init__(steps, dots)
        self.steps = steps
        self.dots = dots
        self.word = tuple(word)
        self.heights = hs
        self.unit = unit

    def vertex(self, j: int) -> np.ndarray:
        if j < len(self.steps):
            return self.steps[j].get_start()
        return self.steps[-1].get_end()

    def y_of(self, h: float) -> float:
        return self.vertex(0)[1] + self.unit * (h - self.heights[0])

    def x_of(self, j: float) -> float:
        return self.vertex(0)[0] + self.unit * j


def level_line(path: PathMob, h: float, j0: float, j1: float, color=GREY_B,
               stroke_width: float = 2, dashed: bool = True) -> VMobject:
    """在路徑的座標系裡，從第 j0 格到第 j1 格、高度 h 的水平線。"""
    a = np.array([path.x_of(j0), path.y_of(h), 0.0])
    b = np.array([path.x_of(j1), path.y_of(h), 0.0])
    if dashed:
        return DashedLine(a, b, dash_length=0.08, stroke_color=color, stroke_width=stroke_width)
    return Line(a, b, stroke_color=color, stroke_width=stroke_width)


def path_thumbnail(word, v_vec, unit: float = 0.17, origin=ORIGIN, mark_bad: bool = True,
                   stroke_width: float = 3.2, dot_radius: float = 0.03) -> VGroup:
    """
    小圖：路徑 + 高度 0 的虛線；壞路徑在第一個高度 <= 0 的點畫紅點。
    屬性：.path、.good、.bad_dot（好路徑時為 None）
    """
    path = PathMob(word, v_vec, unit=unit, origin=origin, stroke_width=stroke_width,
                   dot_radius=dot_radius)
    base = level_line(path, 0, 0, len(word), color=GREY_C, stroke_width=1.5)
    group = VGroup(base, path)
    bad = first_bad_index(word, v_vec)
    group.path = path
    group.good = bad is None
    group.bad_dot = None
    if bad is not None and mark_bad:
        group.bad_dot = Dot(path.vertex(bad), radius=dot_radius * 2.4, fill_color=CROSSING_COLOR)
        group.bad_dot.set_z_index(3)
        group.add(group.bad_dot)
    return group
