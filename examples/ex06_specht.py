"""
範例 6：Specht 模組（Fulton 第 7 章）——先用程式算出 λ=(2,1) 的 polytabloid v_T 與 (12) 的表示矩陣，再用 manim 畫出來。

在 examples/ 資料夾裡執行（這裡有 custom_config.yml）：
    manimgl ex06_specht.py SpechtExample
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))     # 讓 combkit 可以 import（已 pip install 就不需要）

from manimlib import *

from combkit.draw.fulton_draw import assert_in_frame, tableau_mob
from combkit.draw.rep_draw import MatrixMob, TabloidMob, signed_sum
from combkit.draw.scene_base import ChuScene
from combkit.rep_data import T_A, load, polytabloid, rep_matrix21


class SpechtExample(ChuScene):
    def construct(self):
        # 1. 數學：v_T = b_T{T}，以及 (12) 在 (v_a, v_b) 基底下的矩陣（第 j 欄 = 第 j 個基底的像）
        v = polytabloid(T_A)
        rho = rep_matrix21((2, 1, 3))
        assert rho == load("matrices21")["rho12"] and sorted(v.values()) == [-1, 1]

        # 2. 畫面
        self.show_header("Specht 模組與表示矩陣")
        t = tableau_mob(T_A, cs=0.6, center=np.array([-4.5, 1.2, 0.0]))
        terms = [(c, TabloidMob(k, cs=0.5)) for k, c in sorted(v.items(), key=lambda kv: -kv[1])]
        lc = signed_sum(terms).move_to(np.array([0.5, 1.2, 0.0]))
        eq = Tex(r"v_T=", font_size=40).next_to(lc, LEFT, buff=0.25)
        mat = MatrixMob(rho, center=np.array([1.0, -1.3, 0.0]))
        lab = Tex(r"\rho(12)=", font_size=40).next_to(mat, LEFT, buff=0.25)
        assert_in_frame(t, lc, eq, mat, lab)
        self.say("把 ", "$T$", " 的行做交錯和：", "$v_T=\\{T\\}-\\{(13)T\\}$")
        self.play(FadeIn(t), FadeIn(eq), LaggedStartMap(FadeIn, lc, lag_ratio=0.3), run_time=1.5)
        self.say("標準表的 ", "$v_T$", " 是基底；", "$(12)$", " 在這組基底下的矩陣")
        self.play(FadeIn(lab), FadeIn(mat), run_time=1.0)
        self.wait(2.5)
