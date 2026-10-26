"""
範例 20：自迴避漫步與 Flory 指數（Rudnick & Gaspari《Elements of the Random Walk》第 7 章）——先用程式算出
二維正方格子自迴避漫步的精確條數 c_N（DFS 列舉＝文獻值＝暴力列舉全部 4^N 條）與 ⟨R²⟩、
1 − f 展開的逐階部分和（N = 4：256, 28, 124, 96, 100 = c_4，排容原理）、弧線圖的數目 (2k−1)!! 與 Dyson 的計數、
Flory 的 R ∝ m^{3/(d+2)}（數值最大化的斜率＝3/(d+2)），再畫出：⟨R²⟩ 對 N 的點（普通漫步是直線 N），
以及固定種子的 10 步自迴避漫步與普通漫步的足跡。

在 examples/ 資料夾裡執行（這裡有 custom_config.yml）：
    manimgl ex20_saw_flory.py SawFloryExample
"""
import math
import sys
import warnings
from fractions import Fraction
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))     # 讓 combkit 可以 import（已 pip install 就不需要）

from manimlib import *

from combkit import walk_saw as SW
from combkit.draw import lattice_draw as L
from combkit.draw import walk_draw as D
from combkit.draw.fulton_draw import assert_in_frame
from combkit.draw.scene_base import ChuScene
from combkit.draw.styles import zh


class SawFloryExample(ChuScene):
    def construct(self):
        # 1. 數學
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            c, sr = SW.saw_enum(14)
            assert list(c) == SW.SAW_KNOWN_C[:15]                                      # DFS = 文獻值
            assert SW.saw_brute(8) == (c[8], sr[8])                                      # = 暴力列舉 4^8 條
            assert SW.ie_partial(4) == [256, 28, 124, 96, 100] and SW.ie_partial(4)[-1] == c[4]
            assert SW.ie_terms(4) == SW.ie_terms_subsets(4)                              # 兩種排容算法
            full, irr, Dy = SW.dyson_counts(4)
            assert full == Dy == [1, 1, 3, 15, 105] and irr[1:] == [1, 2, 10, 74]
            assert [SW.flory_nu(d) for d in (1, 2, 3, 4)] == [1, Fraction(3, 4), Fraction(3, 5), Fraction(1, 2)]
            assert all(abs(SW.flory_slope(d) - float(SW.flory_nu(d))) < 1e-6 for d in (1, 2, 3, 4))
            assert f"{SW.two_nu_eff(14):.4f}" == "1.4556" and f"{float(SW.saw_R2(14)):.2f}" == "42.79"
            assert f"{SW.cmf_ratio(1000, 10 ** 6, 1):.4f}" == "0.3675"
            saw, rw = SW.sample_saw(10), SW.sample_rw(10)
            assert len(set(saw)) == 11 and sum(v * v for v in saw[-1]) == 18 and sum(v * v for v in rw[-1]) == 2
            R2 = [float(SW.saw_R2(n)) for n in range(15)]

        # 2. 畫面
        self.show_header("自迴避漫步：越走越遠，R 平方約正比於 N 的 3/2 次方")
        ax = D.plot_axes((0, 14, 2), (0, 45, 10), width=5.6, height=4.0, x_ticks=(0, 2, 4, 6, 8, 10, 12, 14),
                         y_ticks=(10, 20, 30, 40), x_label="N", y_label="", font_size=18)
        ax.move_to(np.array([-3.2, -0.2, 0]))
        rw_line = DashedLine(ax.c2p(0, 0), ax.c2p(14, 14)).set_stroke(GREY_B, 3)
        a = R2[14] / 14 ** 1.5
        ref = ax.get_graph(lambda n: a * n ** 1.5, x_range=(0.5, 14)).set_stroke(ORANGE, 2)
        dots = VGroup(*[Dot(ax.c2p(n, R2[n]), radius=0.05).set_color(YELLOW) for n in range(1, 15)])
        lab1 = zh("點：自迴避  虛線：普通漫步 N  橘：N 的 3/2 次方", font_size=18).next_to(ax, UP, buff=0.15)
        res = Tex(f"c_{{14}}={c[14]},\\ \\langle R^2\\rangle_{{14}}={R2[14]:.2f},\\ 2\\nu_{{\\rm eff}}={SW.two_nu_eff(14):.4f}",
                  font_size=22).next_to(lab1, UP, buff=0.12)

        cell = 0.36
        o1 = np.array([2.0, 0.9, 0])
        o2 = np.array([4.6, 0.9, 0])
        t1 = L.trail(saw, cell, o1 - cell * 0.0 * RIGHT, color=YELLOW, stroke_width=4)
        t2 = L.trail(rw, cell, o2 - cell * 2.0 * RIGHT, color=TEAL, stroke_width=4)
        c1 = zh("自迴避 10 步：R 平方 = 18", font_size=18).move_to(np.array([2.0, 3.0, 0]))
        c2 = zh("普通漫步 10 步：R 平方 = 2", font_size=18).move_to(np.array([4.9, 3.0, 0]))
        flory = Tex(r"\nu_F=\frac{3}{d+2}:\ 1,\ \frac34,\ \frac35,\ \frac12", font_size=26).move_to(np.array([3.4, -1.8, 0]))
        part = Tex(r"256\to28\to124\to96\to100=c_4", font_size=24).move_to(np.array([3.4, -2.6, 0]))
        assert_in_frame(ax, rw_line, ref, dots, lab1, res, t1, t2, c1, c2, flory, part)

        self.say("嚴格自迴避：DFS 列舉出 ", "$c_N$", "，和普通漫步比一比")
        self.play(FadeIn(ax), ShowCreation(rw_line), FadeIn(lab1), run_time=1.0)
        self.play(LaggedStartMap(FadeIn, dots, lag_ratio=0.08), ShowCreation(ref), FadeIn(res), run_time=1.8)
        self.say("被自己的足跡推著走：", "$\\langle R^2\\rangle\\propto N^{3/2}$")
        self.play(ShowCreation(t1), ShowCreation(t2), FadeIn(c1), FadeIn(c2), run_time=2.0)
        self.say("Flory：熵與排斥拔河；", "$1-f$", " 展開要全部加完才等於 ", "$c_4$")
        self.play(FadeIn(flory), FadeIn(part), run_time=1.0)
        self.wait(2.0)
