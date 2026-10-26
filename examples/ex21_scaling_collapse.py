"""
範例 21：標度入門與資料塌縮（Rudnick & Gaspari《Elements of the Random Walk》第 8 章）——先用程式算出
一維格子生成函數 g(z;q) = 1/(1 − 2z cos q) 在 z_c − z = 0.1, 0.05, 0.01 的值（g(z;0) = 5, 10, 50），
換成標度座標 x = q/√(z_c − z)、Y = (z_c − z)g 之後三條曲線疊成一條 2/(4 + x²)（(8.8)(8.9)）；
再核對長尾巴步伐的寬度 ∝ N（模擬＝三維 Cauchy 的理論值）、最陡下降 (8.20) 的更正版 e^{−x²/2(N+1)}＝精確二項式（書上少了 1−ν）、
重新標度後自迴避項的放大倍數 N^{(4−d)/2}（上臨界維度 4）與 Flory 的 ν = 3/(d+2)，最後畫出塌縮的動畫。

在 examples/ 資料夾裡執行（這裡有 custom_config.yml）：
    manimgl ex21_scaling_collapse.py ScalingCollapseExample
"""
import math
import sys
import warnings
from fractions import Fraction
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))     # 讓 combkit 可以 import（已 pip install 就不需要）

from manimlib import *

from combkit import walk_saw as SW
from combkit import walk_scale as WS
from combkit.draw import walk_draw as D
from combkit.draw.fulton_draw import assert_in_frame
from combkit.draw.scene_base import ChuScene
from combkit.draw.styles import zh

EPS = WS.COLLAPSE_EPS                     # z_c − z = 0.1, 0.05, 0.01
COLORS = (BLUE, GREEN, YELLOW)


class ScalingCollapseExample(ChuScene):
    def construct(self):
        # 1. 數學
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            assert [f"{WS.g_lattice(0.5 - e, 0.0):.0f}" for e in EPS] == ["5", "10", "50"]
            assert [f"{WS.collapse(e, 2):.4f}" for e in EPS] == ["0.2819", "0.2652", "0.2529"]
            assert f"{WS.collapse_limit(2):.4f}" == "0.2500" and abs(WS.collapse(1e-6, 2) - 0.25) < 1e-5
            for e in EPS:                                                                 # 恆等式：塌縮座標＝原曲線
                q = 0.3
                assert abs(WS.collapse(e, q / math.sqrt(e)) - e * WS.g_lattice(0.5 - e, q)) < 1e-12
            lev, uni = WS.median_width(1000, 1), WS.median_width(1000, None)
            assert f"{lev:.0f}" == "1822" and f"{uni:.2f}" == "28.06"
            assert abs(lev / WS.width_pred(1000, 1) - 1) < 0.04 and abs(uni / WS.width_pred(1000, None) - 1) < 0.03
            ex, fix, book = WS.ratio_exact(20, 8), WS.sd_ratio(20, 8, 0.5, 2, 0.5), WS.sd_ratio(20, 8, 0.5, 2, 0.5, book=True)
            assert (f"{ex:.4f}", f"{fix:.4f}", f"{book:.4f}") == ("0.2098", "0.2179", "0.0475")
            assert [WS.sa_amplitude(100, d) for d in range(1, 6)] == [1000, 100, 10, 1, Fraction(1, 10)]
            assert all(WS.nu_rescale(d) == SW.flory_nu(d) for d in (1, 2, 3, 4))
            assert WS.gamma_exp(WS.ETA_2D, WS.NU_2D) == Fraction(43, 32)
            xs_raw = [i / 25 for i in range(-25, 26)]                                    # q ∈ [−1, 1]
            raw = [[(q, WS.g_lattice(0.5 - e, q)) for q in xs_raw] for e in EPS]
            col = [[(q / math.sqrt(e), e * WS.g_lattice(0.5 - e, q)) for q in xs_raw if abs(q) <= 4 * math.sqrt(e)] for e in EPS]

        # 2. 畫面
        self.show_header("標度形式：三條曲線疊成一條")
        ax1 = D.plot_axes((-1, 1, 0.5), (0, 50, 10), width=4.6, height=3.4, x_ticks=(-1, 0, 1),
                          y_ticks=(10, 30, 50), x_label="q", y_label="", font_size=18)
        ax1.move_to(np.array([-3.6, 0.3, 0]))
        ax2 = D.plot_axes((-4, 4, 1), (0, 0.5, 0.1), width=4.6, height=3.4, x_ticks=(-4, -2, 0, 2, 4),
                          y_ticks=(0.25, 0.5), x_label="x", y_label="", font_size=18)
        ax2.move_to(np.array([2.6, 0.3, 0]))
        c1 = VGroup(*[VMobject().set_points_smoothly([ax1.c2p(a, b) for a, b in pts]).set_stroke(cl, 3)
                      for pts, cl in zip(raw, COLORS)])
        c2 = VGroup(*[VMobject().set_points_smoothly([ax2.c2p(a, b) for a, b in pts]).set_stroke(cl, 3)
                      for pts, cl in zip(col, COLORS)])
        lim = DashedVMobject(ax2.get_graph(lambda x: 2 / (4 + x * x), x_range=(-4, 4)).set_stroke(RED, 2), num_dashes=40)
        t1 = Tex(r"g(z;q)=\frac{1}{1-2z\cos q}", font_size=26).next_to(ax1, UP, buff=0.2)
        t2 = Tex(r"(z_c-z)\,g\ \ \text{vs}\ \ x=\frac{q}{\sqrt{z_c-z}}", font_size=26).next_to(ax2, UP, buff=0.2)
        leg = zh("藍、綠、黃：z_c − z = 0.1、0.05、0.01；紅虛線：2/(4+x²)", font_size=18).move_to(np.array([0, -1.95, 0]))
        nums = Tex(r"\text{width}_{1000}:1822\ \text{vs}\ 28.06,\quad "
                   r"\frac{C(20;8)}{C(20;0)}=0.2098\approx e^{-64/42}=0.2179,\quad N^{(4-d)/2}:1000,100,10,1,0.1",
                   font_size=20).move_to(np.array([0, -2.6, 0]))
        assert_in_frame(ax1, ax2, c1, c2, lim, t1, t2, leg, nums)

        self.say("一維格子：", "$z$", " 越靠近 ", "$z_c=\\tfrac12$", "，曲線越尖越高")
        self.play(FadeIn(ax1), FadeIn(t1), LaggedStartMap(ShowCreation, c1, lag_ratio=0.3), run_time=2.0)
        self.say("縱軸乘 ", "$z_c-z$", "、橫軸除以 ", "$\\sqrt{z_c-z}$")
        self.play(FadeIn(ax2), FadeIn(t2), run_time=0.8)
        self.play(*[TransformFromCopy(a, b) for a, b in zip(c1, c2)], run_time=2.0)
        self.play(ShowCreation(lim), FadeIn(leg), run_time=1.0)
        self.say("塌縮成一條：", "$g=|q|^{-2}f\\big(q(z_c-z)^{-1/2}\\big)$")
        self.play(FadeIn(nums), run_time=1.0)
        self.wait(2.0)
