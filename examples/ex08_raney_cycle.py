"""
範例 8：Raney (1960) 的循環引理（定理 2.1）——字串 σ = 0020030 的秩是 −2，
它的 7 種循環讀法 σg_k 中恰有 2 種是 2 字列（k = 2, 5）。先用 combkit.raney 算出來並和 Chu 講義的路徑核對，
再把每種讀法換成 Chu 講義的路徑（倒過來讀、符號 j 走 1−j）畫出來：好的讀法一直在 0 之上。

在 examples/ 資料夾裡執行（這裡有 custom_config.yml）：
    manimgl ex08_raney_cycle.py RaneyCycleExample
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))     # 讓 combkit 可以 import（已 pip install 就不需要）

from manimlib import *

from combkit.counting import good_rotations, raney_A
from combkit.draw.fulton_draw import assert_in_frame
from combkit.draw.lattice_paths import PathMob
from combkit.draw.scene_base import ChuScene
from combkit.raney import L_formula, cyclic_variant, list_shifts, rank, raney_proof_trace, to_chu_path

SIGMA = (0, 0, 2, 0, 0, 3, 0)
V = (2, 3)                       # Chu 的型：符號 2 = 弦（第 1 型）、符號 3 = 三角形（第 2 型）


class RaneyCycleExample(ChuScene):
    def construct(self):
        # 1. 數學：秩、好的 k、Raney 證明裡的 α 與 σ1，以及和 Chu 講義路徑的對照
        n = -rank(SIGMA)
        ks = list_shifts(SIGMA)
        trace = raney_proof_trace(SIGMA)
        assert n == 2 and ks == [2, 5] and trace["alpha"] == (3,) and trace["words"] == [(3, 0, 0, 0), (2, 0, 0)]
        m = len(SIGMA)
        assert good_rotations(to_chu_path(SIGMA, V), V) == sorted((m - k) % m for k in ks)
        assert L_formula(2, (0, 1, 1)) == raney_A(2, (1, 1), V) == 12

        # 2. 畫面：7 種讀法，各畫成 Chu 講義的路徑；好的讀法打勾
        self.show_header("Raney 的循環引理：秩 −2，恰有 2 種好讀法")
        self.say("$\\sigma=0020030$", "，從第 ", "$k$", " 項開始讀一圈，倒過來畫成路徑")
        cells = VGroup()
        for k in range(m):
            v = cyclic_variant(SIGMA, k)
            path = PathMob(to_chu_path(v, V), V, unit=0.2, stroke_width=4, dot_radius=0.035)
            base = DashedLine(path.vertex(0), path.vertex(0) + m * 0.2 * RIGHT, stroke_color=GREY_C, stroke_width=1.5)
            label = Tex(f"k={k}:\\ " + "".join(map(str, v)), font_size=24).next_to(path, DOWN, buff=0.15)
            good = k in ks
            mark = Tex(r"\checkmark" if good else r"\times", font_size=30,
                       fill_color=GREEN if good else RED).next_to(label, DOWN, buff=0.1)
            cells.add(VGroup(base, path, label, mark))
        cells.arrange_in_grid(n_rows=2, n_cols=4, h_buff=0.6, v_buff=0.45).move_to(0.6 * UP)
        result = Tex(r"7\cdot L(2;0,1,1)=2\cdot\frac{7!}{5!\,1!\,1!}=84\ \Rightarrow\ L(2;0,1,1)=12",
                     font_size=34).to_edge(DOWN, buff=1.1)
        assert_in_frame(cells, result)
        assert cells.get_bottom()[1] > result.get_top()[1] + 0.2
        self.play(LaggedStartMap(FadeIn, cells, lag_ratio=0.15), run_time=2.5)
        self.say("好的讀法是 ", "$k=2,5$", "：路徑出發後一直在 0 之上")
        self.play(*[Indicate(cells[k], color=YELLOW) for k in ks])
        self.say("定理 2.2：每個字串恰有 2 種好讀法，雙重計數")
        self.play(Write(result))
        self.wait(2)
