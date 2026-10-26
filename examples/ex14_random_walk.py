"""
範例 14：一維隨機漫步（Rudnick & Gaspari《Elements of the Random Walk》第 1 章）——先用程式算出
十次擲幣 HTHHTHTTHH 的路線、N = 10 的終點分布（252 種回到原點）與 √N 定律，再畫出：硬幣列、走法圖、
40 條路線的扇形與 ±√N 包絡、二項分布長條與高斯曲線（含奇偶性的因子 2）。

在 examples/ 資料夾裡執行（這裡有 custom_config.yml）：
    manimgl ex14_random_walk.py RandomWalkExample
"""
import sys
from fractions import Fraction
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))     # 讓 combkit 可以 import（已 pip install 就不需要）

from manimlib import *

from combkit import walk1d as W
from combkit.draw import walk_draw as D
from combkit.draw.fulton_draw import assert_in_frame
from combkit.draw.scene_base import ChuScene


class RandomWalkExample(ChuScene):
    def construct(self):
        # 1. 數學
        assert W.MAIN_POS[-1] == 2 and W.MAIN_COINS.count("H") == 6
        assert W.walk_count(10, 0) == 252 and W.walk_prob(10, 0) == Fraction(63, 256)
        assert W.moment_exact(10, 2) == 10 == W.moment_via_cosh(10, 2)
        assert round(W.gauss_prob(10, 0), 4) == 0.2523
        walks = [W.positions(w) for w in W.ensemble(40, 10, seed=3)]

        # 2. 畫面
        self.show_header("一維隨機漫步：擲硬幣、√N、二項分布")
        coins = D.coin_row(W.MAIN_COINS, radius=0.24).move_to(np.array([-3.3, 2.7, 0]))
        o = np.array([-5.9, -0.3, 0])
        assert max(abs(x) for p in walks for x in p) <= 8
        axes = D.walk_axes(10, 8, ux=0.55, uy=0.28, origin=o, y_every=4)
        fan = D.fan(walks, ux=0.55, uy=0.28, origin=o)
        env = D.sqrt_envelope(10, 0.55, 0.28, o)
        path = D.walk_path(W.MAIN_POS, ux=0.55, uy=0.28, origin=o)
        o2 = np.array([3.6, -2.5, 0])
        probs = {d: float(W.walk_prob(10, d)) for d in range(-10, 11, 2)}
        bars = D.bar_chart(probs, ux=0.26, uy=14, origin=o2)
        gauss = D.curve(lambda u: W.gauss_prob(10, u), -10, 10, ux=0.26, uy=14, origin=o2)
        assert_in_frame(coins, axes, fan, env, path, bars, gauss)

        self.say("十次擲幣：正面向右、反面向左，最後停在 ", "$d=2$")
        self.play(LaggedStartMap(FadeIn, coins, lag_ratio=0.15), FadeIn(axes), run_time=1.5)
        self.play(ShowCreation(path[0]), FadeIn(path[1]), run_time=2.0)
        self.say("很多人一起走：寬度只長成 ", "$\\sqrt N$")
        self.play(FadeIn(fan), ShowCreation(env), run_time=1.5)
        self.say("終點分布：", "$C(10,0)=252$", "，柱高", "$\\approx 2\\times$", "高斯密度")
        self.play(FadeIn(bars[0]), run_time=1.0)
        self.play(ShowCreation(gauss), D.bar_of(bars, 0).animate.set_fill(YELLOW, 0.9), run_time=1.2)
        self.wait(2.0)
