"""
範例 24：重整化的流（Rudnick & Gaspari《Elements of the Random Walk》第 11 章）——先用程式做一維的「抽格」：
Ising 環把奇數自旋加總掉，Σ e^{K S3(S2+S4)} = A e^{K' S2 S4}，A = 2√cosh 2K、K' = ½ ln cosh 2K；配分函數變成連乘，
4 個與 8 個自旋的環（K = 1）：暴力列舉＝轉移矩陣＝抽格連乘 = 121.2329、9160.4387（書上 (11.63) 的指數 N/2^{j+1} 會得 18.5498）。
K 一路流向 0：1、0.6625、0.3501、0.1137；換成 x = e^{−2K}：x' = 2x/(1+x²)，不動點 0（不穩定）與 1（穩定）。
一維 Gaussian 鏈積掉奇數格點再重新標度：r'' = r(r+4a)/a（a = 1），小 r 時每次 ×4 ⇒ Λ_r = 2、ν = 1/2；
蛛網 0.05 → 0.2025 → 0.8510 → 4.1282。最後畫出兩張蛛網圖。

在 examples/ 資料夾裡執行（這裡有 custom_config.yml）：
    manimgl ex24_decimation_flow.py DecimationFlowExample
"""
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))     # 讓 combkit 可以 import（已 pip install 就不需要）

from manimlib import *

from combkit import walk_rg as RG
from combkit.draw import walk_draw as D
from combkit.draw.fulton_draw import assert_in_frame
from combkit.draw.scene_base import ChuScene
from combkit.draw.styles import zh


def cobweb_group(ax, pts, color=YELLOW):
    """轉角點 [(x, y), …] → 一段一段的 Line（逐段播放）。"""
    return VGroup(*[Line(ax.c2p(*a), ax.c2p(*b)).set_stroke(color, 3) for a, b in zip(pts[:-1], pts[1:])])


def x_cobweb_points(x0, n):
    xs = RG.x_iterate(x0, n)
    pts = [(xs[0], 0.0)]
    for k in range(n):
        pts.append((xs[k], xs[k + 1]))
        pts.append((xs[k + 1], xs[k + 1]))
    return pts


class DecimationFlowExample(ChuScene):
    def construct(self):
        # 1. 數學
        assert (f"{RG.ising_A(1.0):.4f}", f"{RG.ising_Kp(1.0):.4f}") == ("3.8793", "0.6625")
        assert RG.ising_pair_check(1.0) < 1e-12
        for N in (4, 8):
            zb, zt, zd = RG.ring_Z_brute(N, 1.0), RG.ring_Z_transfer(N, 1.0), RG.ring_Z_decimation(N, 1.0)
            assert abs(zb - zt) < 1e-9 * zb and abs(zb - zd) < 1e-9 * zb
        assert [f"{RG.ring_Z_brute(N, 1.0):.4f}" for N in (4, 8)] == ["121.2329", "9160.4387"]
        assert f"{RG.ring_Z_decimation(4, 1.0, book=True):.4f}" == "18.5498"
        Ks = RG.ising_iterate(1.0, 3)
        assert [f"{k:.4f}" for k in Ks] == ["1.0000", "0.6625", "0.3501", "0.1137"]
        x0 = RG.x_of_K(1.0)
        assert all(abs(RG.x_map(RG.x_of_K(k)) - RG.x_of_K(RG.ising_Kp(k))) < 1e-12 for k in Ks)
        assert RG.x_slope(0.0) == 2.0 and RG.x_slope(1.0) == 0.0
        rs = RG.gauss_iterate(RG.COBWEB_R0, 3)
        assert [f"{r:.4f}" for r in rs] == ["0.0500", "0.2025", "0.8510", "4.1282"]
        assert [f"{v:.4f}" for v in RG.chain_lambdas()[:2]] == ["2.0000", "1.5000"]
        assert RG.ring_decimation_check(16, 0.5, 1.0, 0.2) < 1e-9                  # 環的行列式＝抽格（精確）
        xpts = x_cobweb_points(x0, 4)
        rpts = RG.cobweb_points(RG.COBWEB_R0, 3)

        # 2. 畫面
        self.show_header("抽格：參數的流")
        axl = D.plot_axes((0, 1, 0.25), (0, 1, 0.25), width=4.2, height=4.2, x_ticks=(0.5, 1.0), y_ticks=(0.5, 1.0),
                          x_label="x", y_label="", font_size=18).move_to(np.array([-3.4, 0.2, 0]))
        axr = D.plot_axes((0, 1, 0.25), (0, 5, 1), width=4.2, height=4.2, x_ticks=(0.5, 1.0), y_ticks=(1, 3, 5),
                          x_label="r", y_label="", font_size=18).move_to(np.array([3.4, 0.2, 0]))
        cl = VMobject().set_points_as_corners([axl.c2p(t, RG.x_map(t)) for t in np.linspace(0, 1, 40)]).set_stroke(BLUE, 3)
        dl = DashedLine(axl.c2p(0, 0), axl.c2p(1, 1)).set_stroke(GREY_B, 2)
        cr = VMobject().set_points_as_corners([axr.c2p(t, RG.rr_gauss(t)) for t in np.linspace(0, 0.95, 40)]).set_stroke(BLUE, 3)
        dr = DashedLine(axr.c2p(0, 0), axr.c2p(1, 1)).set_stroke(GREY_B, 2)
        wl, wr = cobweb_group(axl, xpts), cobweb_group(axr, rpts)
        tl = Tex(r"x'=\frac{2x}{1+x^2},\quad x=e^{-2K}", font_size=26).next_to(axl, UP, buff=0.2)
        tr = Tex(r"r''=r(r+4)", font_size=26).next_to(axr, UP, buff=0.2)
        leg = zh("左：一維 Ising 流向 x = 1（高溫）；右：一維 Gaussian 鏈，小 r 時每次乘 4", font_size=18).move_to(np.array([0, -2.35, 0]))
        nums = Tex(r"Z_4=121.2329,\ Z_8=9160.4387;\quad K_n=1,\,0.6625,\,0.3501,\,0.1137;\quad "
                   r"r=0.05,\,0.2025,\,0.8510,\,4.1282", font_size=20).move_to(np.array([0, -2.85, 0]))
        assert_in_frame(axl, axr, cl, dl, cr, dr, wl, wr, tl, tr, leg, nums)

        self.say("把奇數自旋加總掉：", "$K'=\\tfrac12\\ln\\cosh 2K$")
        self.play(FadeIn(axl), FadeIn(tl), ShowCreation(cl), ShowCreation(dl), run_time=1.5)
        self.play(LaggedStartMap(ShowCreation, wl, lag_ratio=0.5), run_time=3.0)
        self.say("Gaussian 鏈：", "$r''=r(r+4a)/a$", "，越跑越大")
        self.play(FadeIn(axr), FadeIn(tr), ShowCreation(cr), ShowCreation(dr), run_time=1.5)
        self.play(LaggedStartMap(ShowCreation, wr, lag_ratio=0.5), run_time=2.5)
        self.play(FadeIn(leg), FadeIn(nums), run_time=1.0)
        self.wait(2.0)
