"""
範例 15：格子上的隨機漫步與生成函數（Rudnick & Gaspari《Elements of the Random Walk》第 2 章）——先用程式算出
正方格子 N = 4 的走法數（遞迴 = 兩枚硬幣 = 布里淵區積分，回到原點 36）、主例 1, 4, 36, 400, 4900, 63504 與高斯近似
（含奇偶性的因子 2），再畫出：數字格與 Fig. 2.3 的四個鄰居、布里淵區上的 (χ/4)^N 熱圖、三種立方晶格的最近鄰。

在 examples/ 資料夾裡執行（這裡有 custom_config.yml）：
    manimgl ex15_lattice_gf.py LatticeGFExample
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))     # 讓 combkit 可以 import（已 pip install 就不需要）

from manimlib import *

from combkit import walk_gf as GF
from combkit.draw import lattice_draw as L
from combkit.draw.fulton_draw import assert_in_frame
from combkit.draw.scene_base import ChuScene


class LatticeGFExample(ChuScene):
    def construct(self):
        # 1. 數學
        grid = GF.count_grid(4)                                           # 遞迴 (2.7)
        assert grid[(0, 0)] == 36 == GF.count_sq(4, 0, 0) == round(GF.fourier_count_sq(4, 0, 0))
        assert sum(grid.values()) == 4 ** 4
        assert [GF.returns_sq(N) for N in range(0, 11, 2)] == [1, 4, 36, 400, 4900, 63504]
        assert round(GF.returns_sq(10) / GF.asym_returns_sq(10), 4) == 0.9513    # 高斯 × 2（奇偶性）
        assert [GF.CHI_SMALL_Q[k][0] for k in ("sc", "bcc", "fcc")] == [6, 8, 12]

        # 2. 畫面
        self.show_header("格子上的漫步：遞迴、傅立葉、普適性")
        o1 = np.array([-4.3, 0.5, 0])
        lines = L.lattice_lines(4, 0.5, o1)
        nums = L.count_grid(grid, 0.5, o1, font_size=20)
        center = L.grid_number(nums, 0, 0).set_color(YELLOW)
        star = L.neighbor_star(0.45, np.array([0.6, -2.3, 0]), labels=False)
        o2 = np.array([0.6, 0.5, 0])
        heat = L.heatmap(lambda x, y: (GF.chi_sq(x, y) / 4) ** 4, n=28, width=3.4, origin=o2)
        bz = L.bz_frame(3.4, o2)
        cubes = VGroup(*[L.lattice3d(k, 0.55, np.array([4.6, 2.0 - 1.75 * i, 0])) for i, k in enumerate(("sc", "bcc", "fcc"))])
        assert_in_frame(lines, nums, star, heat, bz, cubes)

        self.say("每一格＝四個鄰居的和：", "$C(4;\\vec 0)=36$")
        self.play(ShowCreation(lines), FadeIn(star), run_time=1.2)
        self.play(LaggedStartMap(FadeIn, nums, lag_ratio=0.02), run_time=1.5)
        self.play(Indicate(center), run_time=0.8)
        self.say("傅立葉：", "$c(N;\\vec q)=\\chi(\\vec q)^N$", "，四個角落也有峰")
        self.play(FadeIn(heat), FadeIn(bz), run_time=1.2)
        self.say("sc、bcc、fcc：", "$\\chi/Z\\approx1-q^2/6$", "，長路線都是同一個高斯")
        self.play(LaggedStartMap(FadeIn, cubes, lag_ratio=0.3), run_time=1.5)
        self.wait(2.0)
