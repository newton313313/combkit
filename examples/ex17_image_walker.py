"""
範例 17：鏡像走路人與吸收牆（Rudnick & Gaspari《Elements of the Random Walk》第 4 章）——先用程式算出
反射原理（暴力列舉 = 鏡像公式 C(N; x−x0) − C(N; x+x0) = 逐步遞迴）、x0 = 3 的存活機率與 x0·√(2/πN)、
兩牆之間 L = 6 的格子特徵值 2cos(mπ/6)、一維穩定態的「帳篷」與賭徒破產 3/5，再畫出：
一條碰過牆的路線與它的鏡像路線（Fig. 4.1）、存活機率表、帳篷長條圖與流進兩牆的比例。

在 examples/ 資料夾裡執行（這裡有 custom_config.yml）：
    manimgl ex17_image_walker.py ImageWalkerExample
"""
import sys
from fractions import Fraction
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))     # 讓 combkit 可以 import（已 pip install 就不需要）

from manimlib import *

from combkit import walk_bc as BC
from combkit.draw import walk_draw as D
from combkit.draw.fulton_draw import assert_in_frame
from combkit.draw.scene_base import ChuScene
from combkit.draw.styles import zh


class ImageWalkerExample(ChuScene):
    def construct(self):
        # 1. 數學
        for N in range(9):
            assert BC.survivors_brute(N, 1) == BC.survivors_image(N, 1) == BC.survivors_dp(N, 1)
        assert BC.survivors_brute(4, 1) == {1: 2, 3: 3, 5: 1}                       # 16 條裡活著 6 條
        path = [1, 0, 1, 2, 3]
        image = BC.reflect_first(path)
        assert image == [-1, 0, 1, 2, 3] and (path, image) in BC.reflection_pairs(4, 1)
        Ns = (10, 100, 1000)
        exact = [float(BC.survival_prob(N, 3)) for N in Ns]
        asym = [BC.survival_asym(N, 3) for N in Ns]
        assert [f"{v:.4f}" for v in exact] == ["0.6562", "0.2356", "0.0756"]
        assert [f"{v:.4f}" for v in asym] == ["0.7569", "0.2394", "0.0757"]
        assert [round(v, 4) for v in BC.box_eigs(6)] == [1.7321, 1.0, 0.0, -1.0, -1.7321]
        assert BC.box_dp(10, 6, 3) == {1: 81, 3: 162, 5: 81}
        w = BC.tent_lattice(5, 2)
        assert w == [0, Fraction(6, 5), Fraction(12, 5), Fraction(8, 5), Fraction(4, 5), 0]
        assert BC.wall_fluxes(5, 2) == (Fraction(3, 5), Fraction(2, 5)) == (BC.absorb_left(5, 2), 1 - BC.absorb_left(5, 2))

        # 2. 畫面
        self.show_header("鏡像走路人：碰過牆的路線 = 從鏡像點出發的路線")
        ux, uy = 0.55, 0.42
        o1 = np.array([-5.4, 0.5, 0])
        axes = D.walk_axes(4, 3, ux, uy, o1, font_size=18)
        wall = Line(D.path_point(0, 0, ux, uy, o1), D.path_point(4, 0, ux, uy, o1)).set_stroke(D.WRONG_COLOR, 5)
        p_real = D.walk_path(path, ux, uy, o1, color=D.WALK_COLOR)
        p_image = D.walk_path(image, ux, uy, o1, color=D.TAIL_COLOR, stroke_width=3)
        lab = Tex("x_0=1,\\ \\ -x_0=-1", font_size=24).next_to(axes, DOWN, buff=0.3)

        rows = [["N", "P_N", "x_0\\sqrt{2/\\pi N}"]] + [[str(N), f"{e:.4f}", f"{a:.4f}"] for N, e, a in zip(Ns, exact, asym)]
        table = D.tex_table(rows, [0.9, 1.3, 1.7], row_h=0.5, font_size=24)
        table.move_to(np.array([0.5, 1.2, 0]))
        t_lab = zh("x0 = 3 的存活機率", font_size=22).next_to(table, UP, buff=0.2)

        o3 = np.array([2.9, -1.6, 0])
        bars = D.bar_chart({x: float(w[x]) for x in range(1, 5)}, 0.6, 0.85, o3, labels=True,
                           label_fmt=lambda v: f"{v:.1f}", font_size=20)
        flux = Tex("\\text{left } \\tfrac35,\\ \\ \\text{right } \\tfrac25", font_size=26).next_to(bars, UP, buff=0.35)
        b_lab = zh("帳篷：L = 5、x0 = 2", font_size=22).next_to(flux, UP, buff=0.2)
        assert_in_frame(axes, wall, p_real, p_image, lab, table, t_lab, bars, flux, b_lab)

        self.say("碰過牆的路線：第一次碰牆前的那段對 0 鏡射")
        self.play(FadeIn(axes), ShowCreation(wall), run_time=0.8)
        self.play(ShowCreation(p_real), run_time=1.2)
        self.play(TransformFromCopy(p_real, p_image), FadeIn(lab), run_time=1.2)
        self.say("存活機率 ", "$\\approx x_0\\sqrt{2/\\pi N}$", "：一維終究會掉下去")
        self.play(FadeIn(t_lab), FadeIn(table), run_time=1.0)
        self.say("穩定態是帳篷形；流進兩牆的比例 3/5、2/5")
        self.play(FadeIn(b_lab), LaggedStartMap(FadeIn, bars, lag_ratio=0.2), FadeIn(flux), run_time=1.5)
        self.wait(2.0)
