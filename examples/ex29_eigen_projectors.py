"""
範例 29：本徵投影算符（Harter《Principles of Symmetry, Dynamics, and Spectroscopy》第 1 章 §1.2）——先用程式算：
M = [[4, 1], [3, 2]] 的久期方程 ε² − 6ε + 5 = 0（ε = 1、5）；Hamilton–Cayley 因子 𝔭₁ = M − 5、𝔭₅ = M − 1；
投影算符 𝒫₁ = ¼[[1, −1], [−3, 3]]、𝒫₅ = ¼[[3, 1], [3, 1]]（冪等、互相正交、和為 1）；譜分解 M = 𝒫₁ + 5𝒫₅ ⇒ M⁵⁰ = 𝒫₁ + 5⁵⁰𝒫₅（＝ 直接乘）；
預解式 (λ − M)⁻¹ 繞 1、繞 5 的圍道積分（64 點）數值地重現 𝒫₁、𝒫₅。最後畫出：兩個投影矩陣（行＝ket 上青色、列＝bra 上橘色）、
單位圓被 M 拉成的橢圓與兩條本徵方向、向量 (1, 0) = 𝒫₁v + 𝒫₅v 的分解，以及一個 60° 的 Stern–Gerlach 分析器（通過率 3/4）。

在 examples/ 資料夾裡執行（這裡有 custom_config.yml）：
    manimgl ex29_eigen_projectors.py EigenProjectorsExample
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))     # 讓 combkit 可以 import（已 pip install 就不需要）

import sympy as sp
from manimlib import *

from combkit import sym_mat as SM
from combkit.draw import sym_draw as SD
from combkit.draw.fulton_draw import assert_in_frame
from combkit.draw.scene_base import ChuScene
from combkit.draw.styles import zh


class EigenProjectorsExample(ChuScene):
    def construct(self):
        # 1. 數學
        M = SM.MB
        assert SM.secular(M) == sp.Poly(SM.EPS ** 2 - 6 * SM.EPS + 5, SM.EPS) and SM.distinct_eigs(M) == [1, 5]
        P = SM.projectors(M)
        assert P[1] == sp.Matrix([[1, -1], [-3, 3]]) / 4 and P[5] == sp.Matrix([[3, 1], [3, 1]]) / 4
        assert P[1] ** 2 == P[1] and P[1] * P[5] == sp.zeros(2) and P[1] + P[5] == sp.eye(2) and P[1] + 5 * P[5] == M
        assert SM.mat_func(M, lambda e: e ** 50) == M ** 50 == P[1] + 5 ** 50 * P[5]
        assert np.allclose(SM.contour(M, 1, 1), SM.to_np(P[1]), atol=1e-12) and np.allclose(SM.contour(M, 5, 1), SM.to_np(P[5]), atol=1e-12)
        v = sp.Matrix([1, 0])
        assert P[1] * v == sp.Matrix([sp.Rational(1, 4), -sp.Rational(3, 4)]) and P[5] * v == sp.Matrix([sp.Rational(3, 4)] * 2)
        assert SM.pass_fraction(60) == sp.Rational(3, 4) and [f"{x:.4f}" for x in SM.singular_values(M)] == ["5.3983", "0.9262"]

        # 2. 畫
        self.show_header("本徵投影算符")
        p1 = SD.labeled_mat(r"\mathcal P_1", P[1], font_size=30, colors={"col": {0: SD.KET_COLOR}, "row": {0: SD.BRA_COLOR}})
        p5 = SD.labeled_mat(r"\mathcal P_5", P[5], font_size=30, colors={"col": {0: SD.KET_COLOR}, "row": {0: SD.BRA_COLOR}})
        hc = Tex(r"\mathfrak p_1=\mathcal M-5,\quad \mathfrak p_5=\mathcal M-1", font_size=30)
        spec = Tex(r"\mathcal M^{50}=\mathcal P_1+5^{50}\mathcal P_5", font_size=32).set_color(SD.P2_COLOR)
        left = VGroup(hc, p1, p5, spec).arrange(DOWN, buff=0.35, aligned_edge=LEFT).to_corner(UL, buff=0.5).shift(0.55 * DOWN)
        plane = SD.vec_plane(unit=0.55, xr=3, yr=3).move_to(RIGHT * 1.6 + UP * 0.35)
        circ = SD.circle_image(plane, M, alpha=0.0, color=GREY_B, stroke_width=2)
        lines = VGroup(SD.eigen_line(plane, (1, 1), SD.P2_COLOR), SD.eigen_line(plane, (1, -3), SD.P1_COLOR))
        a1 = SD.vec_arrow(plane, (0.25, -0.75), SD.P1_COLOR, r"\mathcal P_1v", font_size=26)
        a5 = SD.vec_arrow(plane, (0.75, 0.75), SD.P2_COLOR, r"\mathcal P_5v", font_size=26)
        vv = SD.vec_arrow(plane, (1, 0), WHITE, "v", font_size=26)
        an = SD.analyzer(60, width=1.2, height=0.8, label=r"60^\circ", font_size=22).move_to(RIGHT * 5.0 + DOWN * 1.7)
        inb = SD.beam(SD.port(an, "in") + 0.9 * LEFT, SD.port(an, "in"), 1.0)
        outs = SD.split_beams(an, 0.75, length=0.9, spread=0.25)
        lab = zh("通過率 3/4", font_size=22).next_to(an, UP, buff=0.15)
        self.add(plane, circ, lines)
        self.play(FadeIn(left), run_time=0.8)
        self.play(UpdateFromAlphaFunc(circ, lambda m, a: m.become(SD.circle_image(plane, M / 2.4, alpha=a, color=WHITE))),
                  run_time=1.5)
        self.play(GrowArrow(vv[0]), FadeIn(vv[1]), run_time=0.5)
        self.play(FadeIn(a1), FadeIn(a5), run_time=0.8)
        self.play(FadeIn(VGroup(an, inb, outs, lab)), run_time=0.6)
        assert_in_frame(left, plane, an, lab)
        self.say("$\\mathcal M$", " 只把兩個分量各乘 1 與 5")
        self.wait(1.5)


if __name__ == "__main__":
    print("請用 manimgl 執行：manimgl ex29_eigen_projectors.py EigenProjectorsExample")
