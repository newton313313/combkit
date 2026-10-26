"""
raney_draw.py —— 整數序列的 Raney 引理（動畫 R01–R06）用到的圖形

    IntPathMob(steps, unit, origin)   一條路徑：每步改變高度 x_j（任意整數）；.vertex(j)、.y_of(h)、.heights
    level_line(path, h, j0, j1)       在高度 h、第 j0 到 j1 點之間畫虛線
    step_tiles(seq)                   一排數字方塊（依步的值上色）
    CardMob(w)                        摸牌遊戲的牌：w ≥ 1 是特殊牌、0 是無用牌、None 是勝利牌
    draws_badge(k)                    「還能摸幾張」的小圓牌

顏色規則與 Chu 講義、lattice_paths.py 一致：步 x = 1 − v 對應 v 點群組，
所以 x = +1（↑ 步）灰、x = −1 藍（弦）、x = −2 橘（三角形）、x = −3 綠；x = 0 淺灰；x ≥ 2 黃（Chu 的路徑沒有這種步）。
"""
from manimlib import *

from combkit.draw.styles import block_color, zh

UP_COLOR = GREY_B
BIG_UP_COLOR = YELLOW_C
WIN_COLOR = GOLD_C
USELESS_COLOR = GREY_C


def step_color(x: int):
    if x == 1:
        return UP_COLOR
    if x >= 2:
        return BIG_UP_COLOR
    return block_color(1 - x)


class IntPathMob(VGroup):
    """
    一條路徑。origin 是第 0 點的位置（高度 start_height）。
    屬性：.steps（每步一條線段）、.dots、.heights（各點高度）
    方法：.vertex(j)（第 j 點目前的座標）、.y_of(h)（高度 h 目前在畫面上的 y）
    """

    def __init__(self, steps, unit: float = 0.4, origin=ORIGIN, start_height: int = 0,
                 stroke_width: float = 5, show_dots: bool = True, dot_radius: float = 0.05, colors=None):
        origin = np.array(origin, dtype=float)
        hs = [start_height]
        for x in steps:
            hs.append(hs[-1] + x)
        pts = [origin + unit * (j * RIGHT + (h - start_height) * UP) for j, h in enumerate(hs)]
        segs = VGroup(*[
            Line(pts[j], pts[j + 1], stroke_color=(colors[j] if colors else step_color(x)), stroke_width=stroke_width)
            for j, x in enumerate(steps)
        ])
        dots = VGroup()
        if show_dots:
            dots = VGroup(*[Dot(p, radius=dot_radius, fill_color=WHITE) for p in pts])
            dots.set_z_index(2)
        super().__init__(segs, dots)
        self.steps = segs
        self.dots = dots
        self.heights = hs
        self.unit = unit
        self.start_height = start_height
        self._anchor = Dot(origin, radius=0.001, fill_opacity=0, stroke_width=0)
        self.add(self._anchor)

    def vertex(self, j: int):
        """第 j 點目前的座標（整條路徑移動或縮放後也正確）。"""
        o = self._anchor.get_center()
        u = self.get_unit()
        return o + u * (j * RIGHT + (self.heights[j] - self.start_height) * UP)

    def get_unit(self) -> float:
        if len(self.steps) == 0:
            return self.unit
        seg = self.steps[0]
        return abs(seg.get_end()[0] - seg.get_start()[0])

    def y_of(self, h: float) -> float:
        return self._anchor.get_center()[1] + self.get_unit() * (h - self.start_height)


def level_line(path: IntPathMob, h: float, j0: float, j1: float, color=GREY_C, stroke_width: float = 2,
               dashed: bool = True) -> VMobject:
    x0 = path._anchor.get_center()[0] + path.get_unit() * j0
    x1 = path._anchor.get_center()[0] + path.get_unit() * j1
    y = path.y_of(h)
    cls = DashedLine if dashed else Line
    return cls(np.array([x0, y, 0]), np.array([x1, y, 0]), stroke_color=color, stroke_width=stroke_width)


def step_tiles(seq, size: float = 0.55, buff: float = 0.08, font_size: float = 30) -> VGroup:
    """一排方塊，每格寫一個步的值（Tex），邊框依值上色。"""
    tiles = VGroup()
    for x in seq:
        sq = Square(side_length=size, stroke_color=step_color(x), stroke_width=3,
                    fill_color=step_color(x), fill_opacity=0.18)
        lab = Tex(str(x), font_size=font_size).move_to(sq)
        tiles.add(VGroup(sq, lab))
    tiles.arrange(RIGHT, buff=buff)
    return tiles


class CardMob(VGroup):
    """
    摸牌遊戲的牌。w ≥ 1：特殊牌（寫 +w，顏色同 w 點群組）；w = 0：無用牌；w = None：勝利牌（星號）。
    label：牌的名字（例如 "S1"、"U1"），寫在牌的上緣。
    """

    def __init__(self, w, label: str = "", width: float = 0.8, height: float = 1.1, font_size: float = 30):
        if w is None:
            color, text = WIN_COLOR, r"\bigstar"
        elif w == 0:
            color, text = USELESS_COLOR, "0"
        else:
            color, text = block_color(w), "+%d" % w
        body = RoundedRectangle(width=width, height=height, corner_radius=0.1,
                                stroke_color=color, stroke_width=3, fill_color=color, fill_opacity=0.18)
        mark = Tex(text, font_size=font_size, fill_color=color if w is None else WHITE).move_to(body)
        parts = [body, mark]
        if label:
            lab = Tex(r"\mathrm{%s}" % label, font_size=font_size * 0.6, fill_color=GREY_A)
            lab.next_to(body.get_top(), DOWN, buff=0.06)
            mark.shift(0.08 * DOWN)
            parts.append(lab)
        super().__init__(*parts)
        self.w = w
        self.body = body


def draws_badge(k: int, radius: float = 0.24, font_size: float = 26) -> VGroup:
    """「還能摸 k 張」的小圓牌：k ≥ 1 綠、k = 0 紅。"""
    color = GREEN_C if k >= 1 else RED_C
    c = Circle(radius=radius, stroke_color=color, stroke_width=3, fill_color=color, fill_opacity=0.2)
    return VGroup(c, Tex(str(k), font_size=font_size).move_to(c))


# ---------------------------------------------------------------------------
# Raney (1960) 導讀（動畫 N01–N12）：字、樹、秩路徑
#   符號 j（j 元運算子）的顏色 = j 點群組的顏色（2 藍、3 黃、4 綠、5 紫），1 棕、0（引數）灰。
# ---------------------------------------------------------------------------
LEAF_COLOR = GREY_B
UNARY_COLOR = LIGHT_BROWN


def symbol_color(j: int):
    if j == 0:
        return LEAF_COLOR
    if j == 1:
        return UNARY_COLOR
    return block_color(j)


def word_tiles(word, size: float = 0.5, buff: float = 0.06, font_size: float = 28, sep_after=()) -> VGroup:
    """字（或字串）的方塊列；sep_after：在這些位置之後多空一點（用來把字列切開）。"""
    tiles = VGroup()
    for x in word:
        sq = Square(side_length=size, stroke_color=symbol_color(x), stroke_width=3,
                    fill_color=symbol_color(x), fill_opacity=0.2)
        tiles.add(VGroup(sq, Tex(str(x), font_size=font_size).move_to(sq)))
    tiles.arrange(RIGHT, buff=buff)
    for k in sorted(sep_after):
        for t in tiles[k + 1:]:
            t.shift(size * 0.5 * RIGHT)
    return tiles


class WordTreeMob(VGroup):
    """
    把字（前序記法）畫成平面樹：葉子依序等距，內部節點在第一個與最後一個小孩的中間；節點依元數上色。
    labels：與字等長的字串串列（運算子名稱、引數名稱），給了就寫在節點裡。
    屬性：.nodes（依前序，= 字的第 k 個符號）、.edges（.edges 的順序同建立順序）、.parent（list）
    方法：.node(k)（第 k 個符號的節點）
    """

    def __init__(self, word, x_gap: float = 0.55, y_gap: float = 0.7, root_pos=ORIGIN, labels=None,
                 node_size: float = 0.22, font_size: float = 24, edge_width: float = 3):
        word = tuple(word)
        parent = [None] * len(word)
        xs, depth = [0.0] * len(word), [0] * len(word)
        pos = [0]
        leaf = [0]

        def rec(p, d):
            k = pos[0]
            pos[0] += 1
            parent[k] = p
            depth[k] = d
            if word[k] == 0:
                xs[k] = leaf[0]
                leaf[0] += 1
                return k
            kids = [rec(k, d + 1) for _ in range(word[k])]
            xs[k] = (xs[kids[0]] + xs[kids[-1]]) / 2
            return k

        rec(None, 0)
        assert pos[0] == len(word), "不是字"
        root_pos = np.array(root_pos, dtype=float)
        pts = [root_pos + x_gap * (xs[k] - xs[0]) * RIGHT + y_gap * depth[k] * DOWN for k in range(len(word))]
        nodes = VGroup()
        for k, x in enumerate(word):
            col = symbol_color(x)
            if labels:
                if x == 0:
                    shape = Square(side_length=node_size * 1.9, stroke_color=col, stroke_width=2,
                                   fill_color=col, fill_opacity=0.15)
                else:
                    shape = Circle(radius=node_size, stroke_color=col, stroke_width=3, fill_color=col, fill_opacity=0.2)
                node = VGroup(shape, Tex(labels[k], font_size=font_size).move_to(shape))
            elif x == 0:
                node = Square(side_length=node_size * 0.9, stroke_color=col, stroke_width=2, fill_color=col,
                              fill_opacity=0.5)
            else:
                node = Circle(radius=node_size * 0.75, stroke_color=col, stroke_width=3, fill_color=col,
                              fill_opacity=0.35)
            node.move_to(pts[k])
            nodes.add(node)
        edges = VGroup(*[Line(pts[parent[k]], pts[k], stroke_color=GREY_C, stroke_width=edge_width)
                         for k in range(1, len(word))])
        for e, k in zip(edges, range(1, len(word))):
            e.put_start_and_end_on(nodes[parent[k]].get_center(), nodes[k].get_center())
        super().__init__(edges, nodes)
        nodes.set_z_index(1)
        self.nodes = nodes
        self.edges = edges
        self.parent = parent
        self.word = word

    def node(self, k: int):
        return self.nodes[k]


def forest_mob(words, gap: float = 0.6, **kwargs) -> VGroup:
    """字列（一串字）畫成森林，由左到右。"""
    return VGroup(*[WordTreeMob(w, **kwargs) for w in words]).arrange(RIGHT, buff=gap, aligned_edge=UP)


def rank_walk(word, unit: float = 0.35, origin=ORIGIN, **kwargs) -> IntPathMob:
    """Raney 的秩路徑：由左往右讀，符號 j 走 j − 1；線段依符號上色。"""
    return IntPathMob([x - 1 for x in word], unit=unit, origin=origin,
                      colors=[symbol_color(x) for x in word], **kwargs)


def chu_walk(word, unit: float = 0.35, origin=ORIGIN, **kwargs) -> IntPathMob:
    """Chu 講義的路徑：把字串倒過來讀，符號 j 走 1 − j；線段依符號上色。"""
    rev = list(reversed(word))
    return IntPathMob([1 - x for x in rev], unit=unit, origin=origin,
                      colors=[symbol_color(x) for x in rev], **kwargs)


def caption_zh(text: str, font_size: float = 26, color=GREY_A) -> Text:
    return zh(text, font_size=font_size, color=color)


def readings_row(seq, starts, center=ORIGIN, spacing: float = 2.5, unit: float = 0.22,
                 label_gap: float = 0.25, mark: str = "good", font_size: float = 26) -> VGroup:
    """
    把「從第 i 點開始讀」的路徑排成一列（i 取自 starts），各條路徑的第 0 點在同一高度，基線（高度 0）畫虛線。
    每格下方寫 i 與標記：mark = "good" 畫 ✓/✗（好讀法與否），"count" 寫正的點數，None 不畫。
    回傳 VGroup(cells)，每格是 VGroup(path, base, label, tag)，另有屬性 .paths、.tags。
    """
    from combkit.raney import positive_count, reading
    center = np.array(center, dtype=float)
    n = len(starts)
    M = len(seq)
    x0 = center[0] - spacing * (n - 1) / 2 - unit * M / 2
    paths = [IntPathMob(reading(seq, i), unit=unit, origin=np.array([x0 + spacing * k, center[1], 0]),
                        stroke_width=4, dot_radius=0.035) for k, i in enumerate(starts)]
    low = min(p.get_bottom()[1] for p in paths)
    cells, tags = VGroup(), VGroup()
    for p, i in zip(paths, starts):
        base = level_line(p, 0, 0, M)
        lab = Tex("i=%d" % i, font_size=font_size)
        lab.move_to(np.array([p.vertex(0)[0] + unit * M / 2, low - label_gap, 0]))
        r = reading(seq, i)
        if mark == "good":
            ok = positive_count(r) == M
            tag = Tex(r"\checkmark" if ok else r"\times", font_size=font_size * 1.2,
                      fill_color=GREEN if ok else RED)
        elif mark == "count":
            tag = Tex(str(positive_count(r)), font_size=font_size * 1.2, fill_color=YELLOW)
        else:
            tag = VGroup()
        tag.next_to(lab, DOWN, buff=0.12)
        cells.add(VGroup(p, base, lab, tag))
        tags.add(tag)
    cells.paths = paths
    cells.tags = tags
    return cells
