"""
範例 23：關在盒子裡的自迴避漫步（Rudnick & Gaspari《Elements of the Random Walk》第 10 章）——先用程式在 L×L 環面（週期邊界）上
精確列舉自迴避漫步（位元遮罩 DFS；8 重對稱＝不用對稱的 DFS＝4^N 條暴力挑選）：4×4 環面的 c_N 在 N < 4 時和無限平面一樣（4、12、36），
N = 4 少了 4 條（直走一圈回到起點），走滿 16 格的有 5712 條（16 × 5712 = 91392＝四維 Gray 碼）；再由生成函數 Σ c_N z^N 算出
步數密度 ρ(z) = ⟨N⟩/V（3×3、4×4、5×5；比值＝−∂lnΓ/∂r 的差分）：漲落的峰 z* = 0.4443、0.4380、0.4272 隨格子變大往 z_c 移、越來越陡——
有限大小的「凝聚」相變。最後畫出 4×4 環面上一條走滿的漫步（越界的一步拆成兩個半段）與三條 ρ(z)。

在 examples/ 資料夾裡執行（這裡有 custom_config.yml）：
    manimgl ex23_torus_saw.py TorusSawExample
"""
import sys
import warnings
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))     # 讓 combkit 可以 import（已 pip install 就不需要）

from manimlib import *

from combkit import walk_cond as WK
from combkit import walk_crit as WC
from combkit.draw import walk_draw as D
from combkit.draw.fulton_draw import assert_in_frame
from combkit.draw.scene_base import ChuScene
from combkit.draw.styles import zh

CELL = 0.8
BOX_CENTER = np.array([-3.6, 0.3, 0])
L_COLORS = {3: BLUE, 4: GREEN, 5: RED}


def torus_point(x, y, side=4):
    origin = BOX_CENTER - (side - 1) / 2 * CELL * (RIGHT + UP)
    return origin + x * CELL * RIGHT + y * CELL * UP


def torus_picture(ex, side=4):
    box = Square(side_length=side * CELL).move_to(BOX_CENTER).set_stroke(GREY_B, 2)
    dots = VGroup(*[Dot(torus_point(x, y, side), radius=0.05).set_fill(GREY_C) for y in range(side) for x in range(side)])
    segs = VGroup(*[Line(torus_point(*a, side), torus_point(*b, side)).set_stroke(YELLOW if w else TEAL_B, 4)
                    for a, b, w in ex["segments"]])
    sx, sy = ex["start"]
    ex_, ey_ = ex["lift"][-1]
    ends = VGroup(Dot(torus_point(sx, sy, side), radius=0.1).set_fill(GREEN),
                  Dot(torus_point(ex_ % side, ey_ % side, side), radius=0.1).set_fill(RED))
    return box, dots, segs, ends


class TorusSawExample(ChuScene):
    def construct(self):
        # 1. 數學
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            c4 = WK.torus_counts(4)
            assert c4 == (1, 4, 12, 36, 96, 264, 624, 1536, 3072, 6240, 9840, 15912, 17472, 19200, 10656, 5712)
            assert c4 == WK.torus_counts_plain(4) and all(WK.torus_counts_brute(4, N) == c4[N] for N in range(7))
            assert [WK.torus_first_diff(L) for L in (3, 4, 5)] == [(3, 4), (4, 4), (5, 4)]        # N = L 時少 4 條
            assert list(c4[:4]) == list(WK.SAW_INF[:4]) and (c4[4], WK.SAW_INF[4]) == (96, 100)
            assert 16 * c4[15] == WK.GRAY_Q4 == 91392 and sum(c4) == 90677
            ex = WK.torus_example("full")
            assert ex["steps"] == "WSSEEENNNEEESSW" and ex["wraps"] == 3 and len(set(ex["sites"])) == 16
            rho = {L: [f"{WK.torus_rho(L, z):.4f}" for z in (0.2, WC.ZC_2D, 1.0)] for L in (3, 4, 5)}
            assert rho == {3: ["0.1391", "0.3374", "0.6823"], 4: ["0.0900", "0.3056", "0.7338"], 5: ["0.0603", "0.2769", "0.7565"]}
            assert all(abs(WK.torus_rho(L, z) - WK.torus_rho_fd(L, z)) < 1e-8 for L in (3, 4, 5) for z in (0.3, 0.6))
            peaks = [WK.torus_peak(L) for L in (3, 4, 5)]
            assert [f"{z:.4f}" for z, _ in peaks] == ["0.4443", "0.4380", "0.4272"]
            assert peaks[0][1] < peaks[1][1] < peaks[2][1]                                        # 峰隨格子長高
            curves = {L: WK.torus_curve(L, 60, 1.2) for L in (3, 4, 5)}

        # 2. 畫面
        self.show_header("環面上的自迴避漫步：走滿與凝聚")
        box, dots, segs, ends = torus_picture(ex)
        ax = D.plot_axes((0, 1.2, 0.2), (0, 1, 0.2), width=4.6, height=3.4, x_ticks=(0.4, 0.8, 1.2),
                         y_ticks=(0.5, 1.0), x_label="z", y_label="", font_size=18)
        ax.move_to(np.array([3.0, 0.3, 0]))
        lines = VGroup(*[VMobject().set_points_as_corners([ax.c2p(z, r) for z, r in zip(*curves[L])])
                         .set_stroke(L_COLORS[L], 3) for L in (3, 4, 5)])
        zc_line = DashedLine(ax.c2p(WC.ZC_2D, 0), ax.c2p(WC.ZC_2D, 1)).set_stroke(GREY_B, 2)
        t1 = Tex(r"4\times4:\ c_{15}=5712", font_size=26).next_to(box, UP, buff=0.25)
        t2 = Tex(r"\rho(z)=\langle N\rangle/V", font_size=26).next_to(ax, UP, buff=0.2)
        leg = zh("黃：越界的半段；藍、綠、紅：L = 3、4、5；虛線 z_c", font_size=18).move_to(np.array([0, -2.0, 0]))
        nums = Tex(r"c_N=1,4,12,36,96\ (\infty:\ 100);\quad 16\times5712=91392;\quad "
                   r"z^*=0.4443,\,0.4380,\,0.4272", font_size=20).move_to(np.array([0, -2.6, 0]))
        assert_in_frame(box, dots, segs, ends, ax, lines, zc_line, t1, t2, leg, nums)

        self.say("週期邊界：走出右邊，從左邊回來")
        self.play(FadeIn(box), FadeIn(dots), FadeIn(t1), run_time=1.0)
        self.play(LaggedStartMap(ShowCreation, segs, lag_ratio=0.15), FadeIn(ends), run_time=3.0)
        self.say("走滿 16 格的有 ", "$5712$", " 條：一盤義大利麵")
        self.play(FadeIn(ax), FadeIn(t2), ShowCreation(zc_line), run_time=1.0)
        self.play(LaggedStartMap(ShowCreation, lines, lag_ratio=0.3), run_time=2.5)
        self.say("格子越大，", "$\\rho(z)$", " 的轉折越陡")
        self.play(FadeIn(leg), FadeIn(nums), run_time=1.0)
        self.wait(2.0)
