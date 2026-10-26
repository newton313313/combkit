"""
範例 12：Schubert 胞腔與 Schubert 多項式（Fulton 第 9–10 章）——先用程式算出 w = 2 1 5 4 3 的 Schubert 多項式
（x1²x2x3 的係數是 2），再畫出：w 的列梯形矩陣（星號數 = 逆序數）、S_3 的 Bruhat 序、係數 2 的兩個 pipe dream。

在 examples/ 資料夾裡執行（這裡有 custom_config.yml）：
    manimgl ex12_schubert.py SchubertExample
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))     # 讓 combkit 可以 import（已 pip install 就不需要）

from manimlib import *

from combkit.draw.fulton_draw import assert_in_frame
from combkit.draw.scene_base import ChuScene
from combkit.draw.schubert_draw import BoxDiagramMob, CellMatrixMob, HasseMob, PipeDreamMob, poly_tex
from combkit.schubert import (bruhat_covers, bruhat_levels, flag_cell_pattern, length, pipe_dreams, pipe_weight,
                              schubert_poly)


class SchubertExample(ChuScene):
    def construct(self):
        # 1. 數學
        w = (2, 1, 5, 4, 3)
        S = schubert_poly(w)
        two = [D for D in pipe_dreams(w) if pipe_weight(D, 5) == (2, 1, 1, 0, 0)]
        assert S[(2, 1, 1, 0, 0)] == 2 == len(two) and length(w) == 4

        # 2. 畫面
        self.show_header("Schubert 胞腔與 pipe dream")
        mat = CellMatrixMob(flag_cell_pattern(w), cs=0.4).move_to(np.array([-4.6, 1.3, 0]))
        hasse = HasseMob(bruhat_levels(3), bruhat_covers(3), dx=1.0, dy=0.8, font_size=24).move_to(np.array([-4.4, -1.6, 0]))
        pds = VGroup(*(PipeDreamMob(D, 5, cs=0.42) for D in two)).arrange(RIGHT, buff=0.6).move_to(np.array([2.2, 0.9, 0]))
        box = BoxDiagramMob((2, 1), 2, 2, cs=0.35).next_to(pds, DOWN, buff=0.5)
        eq = poly_tex({(2, 1, 1, 0, 0): 2}, font_size=34).next_to(box, RIGHT, buff=0.6)
        assert_in_frame(mat, hasse, pds, box, eq)
        self.say("列梯形矩陣：星號個數 = 逆序數 ", "$\\ell(w)=4$")
        self.play(FadeIn(mat), run_time=1.0)
        self.play(LaggedStartMap(Indicate, mat.stars(), lag_ratio=0.2), run_time=1.2)
        self.say("$S_3$", " 的 Bruhat 序")
        self.play(FadeIn(hasse), run_time=1.0)
        self.say("兩個 pipe dream 都有權重 ", "$x_1^2x_2x_3$", "：係數就是 2")
        self.play(FadeIn(pds), FadeIn(box), Write(eq), run_time=1.5)
        self.wait(2.0)
