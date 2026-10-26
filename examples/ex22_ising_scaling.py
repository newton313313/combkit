"""
範例 22：臨界現象與 n → 0（Rudnick & Gaspari《Elements of the Random Walk》第 9 章）——先用程式算出
平均場 m = tanh(m/T) 與二維 Ising（Onsager）的自發磁化（log–log 斜率 β = 0.5000 vs 0.1250）、一維 Ising 的 χ = Σ G(r) = e^{2K}、
三組精確指數的標度律（Rushbrooke、Fisher、hyperscaling、Widom；Fraction）、Landau 雙井的谷底 m* = √(−r/u)、
固定種子 Metropolis 的翻轉等待時間隨格子變長（133、608、2151 sweep），以及 2×3 梯子上 O(n) 關聯函數在 n → 0 時
只剩自迴避路徑（K + K³ + K⁵ = 0.65625；n = 1、2 用暴力加總核對閉式）；最後畫出 M_s(T) 與雙井。

在 examples/ 資料夾裡執行（這裡有 custom_config.yml）：
    manimgl ex22_ising_scaling.py IsingScalingExample
"""
import math
import sys
import warnings
from fractions import Fraction
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))     # 讓 combkit 可以 import（已 pip install 就不需要）

from manimlib import *

from combkit import walk_crit as WC
from combkit.draw import walk_draw as D
from combkit.draw.fulton_draw import assert_in_frame
from combkit.draw.scene_base import ChuScene
from combkit.draw.styles import zh

R_VALUES = (1.0, 0.0, -1.0)
R_COLORS = (GREEN, GREY_B, RED)


def _x_edge(r, fmax=1.4):
    """雙井曲線在 F = fmax 處的 |m|（只畫在座標框裡的那一段）。"""
    m2 = -r + math.sqrt(r * r + 4 * fmax)          # u m⁴/4 + r m²/2 = fmax，u = 1
    return math.sqrt(m2)


class IsingScalingExample(ChuScene):
    def construct(self):
        # 1. 數學
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            tab = WC.ms_table()
            assert [(f"{a:.4f}", f"{b:.4f}") for a, b in tab.values()] == [
                ("0.9575", "0.9980"), ("0.5254", "0.8959"), ("0.1725", "0.6859")]
            assert all(abs(WC.mf_m(x) - WC.mf_m_newton(x)) < 1e-12 for x in WC.MS_X)       # 兩種解法
            assert abs(WC.onsager_M(1.0) - WC.onsager_lowT(1.0)) < 1e-5                    # Onsager＝低溫展開
            b_mf, b_on = WC.beta_slope(WC.mf_m, 1.0), WC.beta_slope(WC.onsager_M, WC.TC_2D)
            assert (f"{b_mf:.4f}", f"{b_on:.4f}") == ("0.5000", "0.1250")
            assert f"{WC.ising1d_chi(1.0):.4f}" == f"{math.exp(2):.4f}" == "7.3891"
            assert abs(WC.ising1d_chi_partial(1.0, 200) - math.exp(2)) < 1e-12
            for e in WC.EXPONENTS.values():                                              # 四條標度律
                law = WC.laws(e)
                assert law["rushbrooke"] == 2 and law["fisher"] == law["widom"] == e["gamma"]
                assert law["hyper_alpha"] == e["alpha"] and law["hyper_delta"] == e["delta"]
            assert WC.widom_exponents(0, Fraction(15, 8)) == (Fraction(1, 8), Fraction(7, 4), 15)
            assert abs(WC.landau_min(-1.0) - 1) < 1e-12 and WC.landau_F(1.0, -1.0) == -0.25
            assert abs(WC.landau_min(-1.0) - WC.landau_closed(-1.0)[0]) < 1e-12
            taus = [WC.flip_tau(L) for L in WC.MC_L]
            assert [f"{t:.0f}" for t in taus] == ["133", "608", "2151"]
            for K in (0.3, 0.5):
                for pair in WC.LADDER_PAIRS:
                    assert abs(WC.ladder_G(1, K, pair) - WC.ladder_G_ising(K, pair)) < 1e-13
                    assert abs(WC.ladder_G(2, K, pair) - WC.ladder_G_xy(K, pair)) < 1e-12
                    assert abs(WC.ladder_G(0, K, pair) - WC.ladder_saw_sum(K, pair)) < 1e-15
            assert f"{WC.ladder_G(0, 0.5, (0, 3)):.5f}" == "0.65625" and WC.ladder_saws(0, 3) == [1, 3, 5]
            xs = [i / 40 * 0.95 for i in range(1, 41)] + [1 - 10 ** (-1.3 - 0.25 * k) for k in range(19)]
            mf_pts = [(x, WC.mf_m(x)) for x in xs] + [(1.0, 0.0), (1.2, 0.0)]
            on_pts = [(x, WC.onsager_M(x * WC.TC_2D)) for x in xs] + [(1.0, 0.0), (1.2, 0.0)]

        # 2. 畫面
        self.show_header("臨界現象：平均場 vs 二維 Ising")
        ax1 = D.plot_axes((0, 1.2, 0.2), (0, 1, 0.2), width=4.6, height=3.4, x_ticks=(0.5, 1.0),
                          y_ticks=(0.5, 1.0), x_label="T/T_c", y_label="", font_size=18)
        ax1.move_to(np.array([-3.6, 0.4, 0]))
        ax2 = D.plot_axes((-2, 2, 1), (-0.5, 1.5, 0.5), width=4.6, height=3.4, x_ticks=(-1, 0, 1),
                          y_ticks=(-0.5, 0.5, 1.0), x_label="m", y_label="", font_size=18)
        ax2.move_to(np.array([3.0, 0.4, 0]))
        c_mf = VMobject().set_points_as_corners([ax1.c2p(a, b) for a, b in mf_pts]).set_stroke(BLUE, 3)
        c_on = VMobject().set_points_as_corners([ax1.c2p(a, b) for a, b in on_pts]).set_stroke(YELLOW, 3)
        wells = VGroup(*[ax2.get_graph(lambda m, r=r: WC.landau_F(m, r), x_range=(-_x_edge(r), _x_edge(r), 0.02))
                         .set_stroke(cl, 3) for r, cl in zip(R_VALUES, R_COLORS)])
        balls = VGroup(*[Dot(ax2.c2p(s * WC.landau_min(-1.0), -0.25), radius=0.08).set_fill(RED_A)
                         for s in (-1, 1)])
        t1 = Tex(r"M_s\propto(T_c-T)^{\beta}", font_size=26).next_to(ax1, UP, buff=0.2)
        t2 = Tex(r"F=\frac r2m^2+\frac u4m^4", font_size=26).next_to(ax2, UP, buff=0.2)
        leg = zh("藍：平均場 β = 1/2；黃：Onsager β = 1/8；綠、灰、紅：r = 1、0、−1", font_size=18)
        leg.move_to(np.array([0, -1.95, 0]))
        nums = Tex(r"\tau_{\rm flip}(L=4,6,8)=133,\,608,\,2151;\quad "
                   r"\chi_{1D}=e^{2}=7.3891;\quad G_{03}(n\to0)=K+K^3+K^5=0.65625",
                   font_size=20).move_to(np.array([0, -2.6, 0]))
        assert_in_frame(ax1, ax2, c_mf, c_on, wells, balls, t1, t2, leg, nums)

        self.say("自發磁化：平均場與二維 Ising 的精確解")
        self.play(FadeIn(ax1), FadeIn(t1), ShowCreation(c_mf), ShowCreation(c_on), run_time=2.0)
        self.say("斜率 ", "$\\beta=0.5000$", " 與 ", "$0.1250$", "：不同的普適類")
        self.play(FadeIn(ax2), FadeIn(t2), LaggedStartMap(ShowCreation, wells, lag_ratio=0.3), run_time=2.0)
        self.say("$r<0$", "：兩個谷，", "$m^*=\\pm\\sqrt{-r/u}$", "——對稱破缺")
        self.play(FadeIn(balls), FadeIn(leg), run_time=1.0)
        self.say("$n\\to0$", "：只剩自迴避路徑")
        self.play(FadeIn(nums), run_time=1.0)
        self.wait(2.0)
