"""
範例 16：回家與維度（Rudnick & Gaspari《Elements of the Random Walk》第 3 章）——先用程式算出
一維 N 步都不回家的走法數 1, 2, 2, 4, 6, 12, 20, 40, 70, 140, 252（生成函數 = 暴力列舉）、
一、二、三維「1000 步還沒回家」的機率（更新方程 = 逐步遞迴）、三種立方格子的 Watson 積分（Gamma 閉式 = Bessel 積分 = 級數）
與平均足跡，再畫出：二維 400 步的足跡、還沒回家的機率表、sc／bcc／fcc 與它們的逃脫機率。

在 examples/ 資料夾裡執行（這裡有 custom_config.yml）：
    manimgl ex16_recurrence.py RecurrenceExample
"""
import sys
from fractions import Fraction
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))     # 讓 combkit 可以 import（已 pip install 就不需要）

from manimlib import *

from combkit import walk_rec as R
from combkit.draw import lattice_draw as L
from combkit.draw import walk_draw as D
from combkit.draw.fulton_draw import assert_in_frame
from combkit.draw.scene_base import ChuScene


class RecurrenceExample(ChuScene):
    def construct(self):
        # 1. 數學
        NR = [R.never_return_1d(N) for N in range(11)]
        assert NR == [1, 2, 2, 4, 6, 12, 20, 40, 70, 140, 252] == R.never_return_1d_gf(10)
        assert all(R.never_return_1d(N) == R.never_return_1d_brute(N) for N in range(11))
        assert abs(R.no_return_probs("sq", 30)[30] - R.never_return_dp("sq", 30) / 4 ** 30) < 1e-12
        P = [R.no_return_probs(k, 1000)[1000] for k in ("1d", "sq", "sc")]
        assert [round(p, 4) for p in P] == [0.0252, 0.3218, 0.6685]
        assert abs(R.watson_sc_bessel() - R.G_C["sc"]) < 1e-10 and abs(R.watson_bcc_series() - R.G_C["bcc"]) < 1e-10
        esc = [R.escape(k) for k in ("sc", "bcc", "fcc")]
        assert [round(e, 4) for e in esc] == [0.6595, 0.7178, 0.7437]
        assert R.distinct_sites_1d(4) == R.distinct_sites_brute("1d", 4) == Fraction(27, 8)
        tr = R.trail("sq", R.TRAIL_N, R.TRAIL_SEED)
        sites = set(tr)
        assert len(sites) == R.TRAIL_SITES == 159

        # 2. 畫面
        self.show_header("回家與維度：一、二維一定回家，三維不一定")
        cell = 0.13
        xs, ys = [p[0] for p in tr], [p[1] for p in tr]
        cx, cy = (min(xs) + max(xs)) / 2, (min(ys) + max(ys)) / 2
        o1 = np.array([-4.2 - cx * cell, 0.1 - cy * cell, 0])
        cells = L.heat_cells({p: 1.0 for p in sites}, cell, o1, vmax=2.0)
        path = L.trail(tr, cell, o1, stroke_width=1.6)
        rows = [["d", "P_{1000}", "P_\\infty"], ["1", f"{P[0]:.4f}", "0"], ["2", f"{P[1]:.4f}", "0"],
                ["3", f"{P[2]:.4f}", f"{esc[0]:.4f}"]]
        table = D.tex_table(rows, [0.7, 1.5, 1.5], row_h=0.55, font_size=28).move_to(np.array([0.6, 0.6, 0]))
        cubes = VGroup()
        for i, (k, e) in enumerate(zip(("sc", "bcc", "fcc"), esc)):
            cube = L.lattice3d(k, 0.42, np.array([4.3, 2.1 - 1.75 * i, 0]))
            lab = Tex(f"{k}:\\ {e:.4f}", font_size=26).next_to(cube, RIGHT, buff=0.25)
            cubes.add(VGroup(cube, lab))
        assert_in_frame(cells, path, table, cubes)

        self.say("二維 400 步：踩過 159 個點；", "$P_N\\approx\\pi/\\ln N$", "，遲早回家")
        self.play(FadeIn(cells), run_time=0.8)
        self.play(ShowCreation(path), run_time=2.0)
        self.say("$N$", " 步都沒回家的機率：", "$P_\\infty=1/G(z_c;0)$")
        self.play(FadeIn(table), run_time=1.0)
        self.say("三種立方格子的逃脫機率（Watson 積分）")
        self.play(LaggedStartMap(FadeIn, cubes, lag_ratio=0.3), run_time=1.5)
        self.wait(2.0)
