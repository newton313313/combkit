"""
範例 19：隨機漫步的形狀（Rudnick & Gaspari《Elements of the Random Walk》第 6 章）——先用程式算出
迴轉半徑張量（定義式＝兩兩差＝Kramers 恆等式 (6.21)）、4 步「東東北西」的精確 T 與主迴轉半徑 3/5、1/5、
固定種子的二維 1000 步漫步的主軸（旋轉不變）、平均非球度 A_3 = 10/19（公式＝有限 N 的精確動差＝模擬）、
個別非球度的平均 A'_3 = 0.394274，再畫出：二維 1000 步的路線＋重心＋兩條主軸（半長 2R_i），以及三維 2 萬條 100 步的個別非球度直方圖。

在 examples/ 資料夾裡執行（這裡有 custom_config.yml）：
    manimgl ex19_walk_shape.py WalkShapeExample
"""
import math
import sys
from fractions import Fraction
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))     # 讓 combkit 可以 import（已 pip install 就不需要）

from manimlib import *

from combkit import walk_shape as SH
from combkit.draw import walk_draw as D
from combkit.draw.fulton_draw import assert_in_frame
from combkit.draw.scene_base import ChuScene
from combkit.draw.styles import zh


class WalkShapeExample(ChuScene):
    def construct(self):
        # 1. 數學
        se = SH.small_example()
        assert se["T"] == [[Fraction(14, 25), Fraction(3, 25)], [Fraction(3, 25), Fraction(6, 25)]]
        assert SH.tensor_kramers_frac(SH.SMALL_STEPS) == se["T"]                     # (6.21) = 定義式
        assert se["eig"] == (Fraction(3, 5), Fraction(1, 5)) and se["asph"] == Fraction(1, 4)
        pts = SH.fig64()
        fs = SH.fig64_stats()
        assert np.allclose(fs["Tkr"], fs["T"]) and np.allclose(fs["eig_rot"], fs["eig"])   # Kramers、旋轉不變
        assert [f"{v:.1f}" for v in fs["eig"]] == ["258.8", "31.6"] and f"{fs['angle']:.1f}" == "28.4"
        assert SH.A_d(3) == Fraction(10, 19) and SH.A_S21(3) == SH.A_d(3)
        st = SH.sim3d_stats()
        assert f"{st['A_ratio']:.4f}" == "0.5231" and f"{float(SH.A_exact(100, 3, 'lattice')):.4f}" == "0.5235"
        assert f"{SH.Aprime(3):.6f}" == "0.394274" and abs(SH.Aprime(3) - SH.Aprime_laplace(3)) < 1e-6
        hist = st["hist"]
        assert int(hist.sum()) == 20000 and len(hist) == 20

        # 2. 畫面
        self.show_header("隨機漫步的形狀：主迴轉半徑與非球度")
        c = fs["com"]
        unit = 0.095
        origin = np.array([-3.4, 0.2, 0])

        def P(x, y):
            return origin + unit * (x - c[0]) * RIGHT + unit * (y - c[1]) * UP

        path = VMobject().set_points_as_corners([P(x, y) for x, y in pts]).set_stroke(D.WALK_COLOR, 1.2, opacity=0.8)
        vals, vecs = SH.principal(fs["T"])
        axes_g = VGroup()
        for lam, col in zip(vals, (ORANGE, TEAL)):
            v = vecs[:, list(vals).index(lam)]
            h = 2 * math.sqrt(lam) * unit
            d = v[0] * RIGHT + v[1] * UP
            axes_g.add(Line(origin - h * d, origin + h * d).set_stroke(col, 5))
        dot = Dot(origin, radius=0.07).set_color(RED)
        lab = Tex(f"R_1^2={vals[0]:.1f},\\ R_2^2={vals[1]:.1f}", font_size=24).move_to(origin + 2.2 * UP)
        cap1 = zh("二維 1000 步：重心與兩條主軸", font_size=22).move_to(origin + 2.7 * UP)

        ax = D.plot_axes((0, 1, 0.2), (0, 2000, 500), width=5.4, height=3.6, x_ticks=(0, 0.2, 0.4, 0.6, 0.8, 1.0),
                         y_ticks=(500, 1000, 1500, 2000), x_label="A", y_label="", font_size=18)
        ax.move_to(np.array([3.3, -0.3, 0]))
        bars = VGroup()
        for i, h in enumerate(hist):
            if h == 0:
                continue
            r = Rectangle(width=0.05 * ax.x_axis.get_unit_size() * 0.92, height=h * ax.y_axis.get_unit_size())
            r.set_fill(D.BAR_COLOR, 0.85).set_stroke(width=0)
            r.move_to(ax.c2p(0.05 * i + 0.025, h / 2))
            bars.add(r)
        mark = DashedLine(ax.c2p(10 / 19, 0), ax.c2p(10 / 19, 2000)).set_stroke(ORANGE, 3)
        cap2 = zh("三維 2 萬條 100 步：個別非球度", font_size=22).next_to(ax, UP, buff=0.55)
        res = Tex(f"A_3=\\frac{{10}}{{19}}=0.5263,\\ \\text{{sim}}\\ {st['A_ratio']:.4f};\\ A'_3={SH.Aprime(3):.4f}",
                  font_size=22).next_to(ax, UP, buff=0.12)
        assert_in_frame(path, axes_g, dot, lab, cap1, ax, bars, mark, cap2, res)

        self.say("單獨一條漫步幾乎總是長條：主軸＝張量的特徵向量")
        self.play(ShowCreation(path), FadeIn(cap1), run_time=2.0)
        self.play(FadeIn(dot), *[ShowCreation(a) for a in axes_g], FadeIn(lab), run_time=1.2)
        self.say("三維平均非球度 ", "$10/19$", "，但個別形狀從圓球到直線都有")
        self.play(FadeIn(ax), FadeIn(cap2), run_time=0.8)
        self.play(LaggedStartMap(FadeIn, bars, lag_ratio=0.05), ShowCreation(mark), FadeIn(res), run_time=1.8)
        self.wait(2.0)
