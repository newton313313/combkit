"""
範例 25：Wilson–Fisher 不動點與 ε 展開（Rudnick & Gaspari《Elements of the Random Walk》第 12 章）——先用程式解一迴圈的遞迴方程：
dr/dℓ = 2r + 12u/(1+r)、du/dℓ = εu − 36u²/(1+r)²（書上 (12.37)(12.38) 的 d 併進 u；O(n)：12 → 4(n+2)、36 → 4(n+8)）。
ε = 1 的 Wilson–Fisher 點精確是 (−1/7, 1/49)（fsolve 一致），線性化矩陣 [[5/3, 14], [1/21, −1]]、特徵值 1.8968、−1.2301；
臨界線（二分法求 r_c(u0)）流進不動點，旁邊的流線往 r → ±∞ 跑掉。ν 的一階公式 ½ + (n+2)ε/(4(n+8))：ε = 1 時
n = 0（自迴避）0.5625、n = 1（Ising）0.5833——已知值（只引用）0.588、0.630。殼層的配對：5 個片段 → 3 個 type 1＋1 個 type 2
有 5!/(3!2!) = 10 種（列舉 945 種配對）。最後畫出 (r, u) 平面的流與 ν 的小表。

在 examples/ 資料夾裡執行（這裡有 custom_config.yml）：
    manimgl ex25_wilson_fisher.py WilsonFisherExample
"""
import sys
from fractions import Fraction
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))     # 讓 combkit 可以 import（已 pip install 就不需要）

from manimlib import *

from combkit import walk_eps as E
from combkit.draw import walk_draw as D
from combkit.draw.fulton_draw import assert_in_frame
from combkit.draw.scene_base import ChuScene
from combkit.draw.styles import zh


def flow_curve(ax, pts, color):
    """流線（一個 VMobject）＋末端的小箭頭。"""
    line = VMobject().set_points_as_corners([ax.c2p(r, u) for r, u in pts]).set_stroke(color, 3)
    b = ax.c2p(*pts[-1])
    a = ax.c2p(*pts[-3])
    a = b + 0.3 * (a - b) / max(np.linalg.norm(a - b), 1e-9)          # 固定長度 0.3 的箭頭
    tip = Arrow(a, b, buff=0, thickness=2.0, max_tip_length_to_length_ratio=0.8).set_color(color)
    return VGroup(line, tip)


class WilsonFisherExample(ChuScene):
    def construct(self):
        # 1. 數學
        fp = E.fixed_point(Fraction(1))
        assert fp == (Fraction(-1, 7), Fraction(1, 49)) and E.beta(fp[0], fp[1], Fraction(1)) == (0, 0)
        fn = E.fixed_point_numeric(1.0)
        assert abs(fn[0] + 1 / 7) < 1e-10 and abs(fn[1] - 1 / 49) < 1e-10
        assert np.allclose(E.jacobian(float(fp[0]), float(fp[1]), 1.0), [[5 / 3, 14], [1 / 21, -1]])
        assert [f"{x:.4f}" for x in E.eigen(1.0)] == ["1.8968", "-1.2301"] and E.eigen(1.0, 1, "gauss") == (2.0, 1.0)
        nus = [E.nu_first(1, n) for n in range(4)]
        assert nus == [Fraction(9, 16), Fraction(7, 12), Fraction(3, 5), Fraction(27, 44)]
        assert [f"{float(v):.4f}" for v in nus] == ["0.5625", "0.5833", "0.6000", "0.6136"]
        assert [E.AB(n) for n in (0, 1)] == [(8, 32), (12, 36)] and E.on_counts(1) == (6.0, 36.0)
        assert len(E.matchings(5)) == 945 and E.diagram_census(5)[(2, 1, 1, 1)] == E.pair_factor(3, 1) == 10
        starts = E.flow_starts()
        lines = [E.flow_line(r0, u0) for r0, u0 in starts]
        for k, pts in enumerate(lines):
            if k % 3 == 1:                                   # 臨界線流進 WF 點
                assert abs(pts[-1][0] + 1 / 7) < 2e-3 and abs(pts[-1][1] - 1 / 49) < 2e-4
            else:
                assert abs(abs(pts[-1][0]) - 0.45) < 1e-6

        # 2. 畫面
        self.show_header("Wilson–Fisher 不動點")
        ax = D.plot_axes((-0.45, 0.45, 0.15), (0, 0.05, 0.01), width=6.4, height=4.4, x_ticks=(-0.3, 0.0, 0.3),
                         y_ticks=(0.02, 0.04), x_label="r", y_label="u", font_size=18).move_to(np.array([-2.6, 0.15, 0]))
        cols = [GREY_B, YELLOW, GREY_B, GREY_B, ORANGE, GREY_B]
        flows = VGroup(*[flow_curve(ax, pts, c) for pts, c in zip(lines, cols)])
        g_dot = Dot(ax.c2p(0, 0), radius=0.08).set_color(WHITE)
        wf_dot = Dot(ax.c2p(float(fp[0]), float(fp[1])), radius=0.08).set_color(WHITE)
        wf_lab = Tex(r"\big(-\tfrac17,\tfrac1{49}\big)", font_size=24).next_to(wf_dot, UL, buff=0.08)
        g_lab = zh("Gaussian", font_size=16).next_to(g_dot, DR, buff=0.08)
        eqs = Tex(r"\frac{dr}{d\ell}=2r+\frac{12u}{1+r},\quad \frac{du}{d\ell}=\epsilon u-\frac{36u^2}{(1+r)^2}",
                  font_size=24).move_to(np.array([-2.6, 2.75, 0]))
        rows = [[r"n", r"\nu_1", r"\text{known}"], [r"0", r"0.5625", r"0.588"], [r"1", r"0.5833", r"0.630"],
                [r"2", r"0.6000", r"0.672"], [r"3", r"0.6136", r"0.711"]]
        tab = D.tex_table(rows, col_widths=(0.6, 1.2, 1.2), row_h=0.48, font_size=22).move_to(np.array([4.2, 0.6, 0]))
        tcap = Tex(r"\nu=\tfrac12+\tfrac{(n+2)\epsilon}{4(n+8)},\ \epsilon=1", font_size=22).next_to(tab, UP, buff=0.25)
        leg = zh("黃、橘：臨界線流進 WF 點；灰：離開。已知值只引用", font_size=18).move_to(np.array([0, -2.75, 0]))
        assert_in_frame(ax, flows, g_dot, wf_dot, wf_lab, g_lab, eqs, tab, tcap, leg)

        self.say("一迴圈的遞迴：", "$\\epsilon=4-d=1$")
        self.play(FadeIn(ax), FadeIn(eqs), FadeIn(g_dot), FadeIn(g_lab), run_time=1.5)
        self.play(LaggedStartMap(ShowCreation, flows, lag_ratio=0.3), run_time=4.0)
        self.say("臨界線流進 Wilson–Fisher 點")
        self.play(FadeIn(wf_dot), FadeIn(wf_lab), run_time=1.0)
        self.say("一階的 ", "$\\nu$", "：自迴避 0.5625、Ising 0.5833")
        self.play(FadeIn(tcap), FadeIn(tab), FadeIn(leg), run_time=1.5)
        self.wait(2.0)
