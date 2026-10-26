"""
circle_points.py —— 圓、圓周上的點、弦

點的編號一律從 0 開始，從正上方出發、順時針排列。
（畫面上的標籤預設顯示 1, 2, 3, ...，與論文和講義一致。）
"""
from manimlib import *

from combkit.draw.styles import (
    BLOCK_STROKE,
    CIRCLE_COLOR,
    CIRCLE_STROKE,
    DOT_RADIUS,
    POINT_COLOR,
    TYPE1_COLOR,
)


def point_angle(i: int, n: int, start_angle: float = PI / 2, clockwise: bool = True) -> float:
    """第 i 個點（共 n 個）在圓上的角度。"""
    sign = -1 if clockwise else 1
    return start_angle + sign * TAU * i / n


def circle_coords(n: int, radius: float = 2.5, center=ORIGIN,
                  start_angle: float = PI / 2, clockwise: bool = True) -> np.ndarray:
    """回傳 n 個等距圓周點的座標，shape = (n, 3)。"""
    center = np.array(center, dtype=float)
    return np.array([
        center + radius * np.array([np.cos(a), np.sin(a), 0.0])
        for a in (point_angle(i, n, start_angle, clockwise) for i in range(n))
    ])


def make_circle_points(n: int, radius: float = 2.5, center=ORIGIN,
                       dot_radius: float = DOT_RADIUS, circle_color=CIRCLE_COLOR,
                       circle_stroke: float = CIRCLE_STROKE, point_color=POINT_COLOR):
    """
    建立一個圓與圓周上的 n 個點。
    回傳 (circle, dots, coords)：
        circle : Circle
        dots   : VGroup，dots[i] 是第 i 個點
        coords : np.ndarray，coords[i] 是第 i 個點的座標
    """
    center = np.array(center, dtype=float)
    circle = Circle(radius=radius, stroke_color=circle_color, stroke_width=circle_stroke)
    circle.move_to(center)
    coords = circle_coords(n, radius, center)
    dots = VGroup(*[Dot(p, radius=dot_radius, fill_color=point_color) for p in coords])
    dots.set_z_index(2)   # 點永遠畫在弦與多邊形的上層
    return circle, dots, coords


def chord_between(coords, i: int, j: int, color=TYPE1_COLOR, stroke_width: float = BLOCK_STROKE) -> Line:
    """連接第 i 與第 j 個點的弦。"""
    return Line(coords[i], coords[j], stroke_color=color, stroke_width=stroke_width)


def label_points(coords, center=ORIGIN, labels=None, buff: float = 0.32,
                 font_size: float = 26, color=GREY_A, start: int = 1) -> VGroup:
    """在每個點的外側放上標籤（預設 1, 2, 3, ...）。"""
    center = np.array(center, dtype=float)
    if labels is None:
        labels = [str(i + start) for i in range(len(coords))]
    result = VGroup()
    for p, text in zip(coords, labels):
        direction = normalize(np.array(p) - center)
        lab = Tex(text, font_size=font_size, fill_color=color)
        lab.move_to(np.array(p) + buff * direction)
        result.add(lab)
    return result


def chords_intersection(coords, chord_a, chord_b):
    """兩條弦（各為一對點編號）在圓內的交點；不相交則回傳 None。"""
    p1, p2 = np.array(coords[chord_a[0]]), np.array(coords[chord_a[1]])
    p3, p4 = np.array(coords[chord_b[0]]), np.array(coords[chord_b[1]])
    d1, d2 = p2 - p1, p4 - p3
    denom = d1[0] * d2[1] - d1[1] * d2[0]
    if abs(denom) < 1e-9:
        return None
    diff = p3 - p1
    t = (diff[0] * d2[1] - diff[1] * d2[0]) / denom
    s = (diff[0] * d1[1] - diff[1] * d1[0]) / denom
    if 1e-6 < t < 1 - 1e-6 and 1e-6 < s < 1 - 1e-6:
        return p1 + t * d1
    return None
