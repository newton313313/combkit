"""
範例 2：Young 表的列插入（Fulton《Young Tableaux》p.7 的例子）。
combkit.tableaux 算出新表與擠壓路徑，combkit.draw.mit_draw 畫 Young 表。

在 examples/ 資料夾裡執行：
    manimgl ex02_row_insertion.py RowInsertionExample
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from manimlib import *

from combkit.draw.mit_draw import YoungDiagramMob
from combkit.draw.scene_base import ChuScene
from combkit.tableaux import FULTON_INSERT_T, is_tableau, row_insert, shape

X = 2


class RowInsertionExample(ChuScene):
    def construct(self):
        # 1. 數學
        T = FULTON_INSERT_T
        T2, route, new_box = row_insert(T, X)
        assert is_tableau(T2) and new_box == (3, 2)          # 書上的新格子是第 4 列第 3 格（這裡 0 起算）

        # 2. 畫面
        self.show_header(f"列插入 T ← {X}")
        before = YoungDiagramMob(shape(T), cell_size=0.7, center=3.2 * LEFT, entries=T)
        after = YoungDiagramMob(shape(T2), cell_size=0.7, center=3.2 * RIGHT, entries=T2)
        arrow = Arrow(before.get_right(), after.get_left(), buff=0.4)
        label = Tex(rf"\leftarrow {X}").next_to(arrow, UP)
        self.play(FadeIn(before))
        self.say("程式算出新表與擠壓路徑，畫面只負責把它們標出來")
        self.play(GrowArrow(arrow), FadeIn(label), FadeIn(after))
        marks = VGroup(*[after.cell(i, j).copy().set_fill(YELLOW, 0.35) for i, j in route])
        self.play(FadeIn(marks))
        self.wait(2)
