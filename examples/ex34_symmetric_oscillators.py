"""
範例 34：對稱群與耦合振盪子（symosc 專案）——先用程式算出三個擺的 D3 約化 A1 ⊕ E（頻率平方 1 與 2.5）、
帶相位耦合下 E 的裂距 √3 sin φ、Q8 的 8 個擺（E 出現兩次、四重簡併、Frobenius–Schur 指標 −1），再畫出：
三角形上的 E 模式在擺動、C3 的相量 (1, ω, ω²)、Q8 的 Cayley 圖與 E 的相量 v1。

在 examples/ 資料夾裡執行（這裡有 custom_config.yml）：
    manimgl ex34_symmetric_oscillators.py SymmetricOscillatorsExample
"""
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))     # 讓 combkit 可以 import（已 pip install 就不需要）

from manimlib import *

from combkit import symosc as S
from combkit.draw import osc_draw as O
from combkit.draw.fulton_draw import assert_in_frame
from combkit.draw.scene_base import ChuScene


class SymmetricOscillatorsExample(ChuScene):
    def construct(self):
        # 1. 數學
        assert S.decompose(S.d3_chars_of(S.d3_matrix), S.D3_CHARS, S.D3_SIZES) == {"A1": 1, "A2": 0, "E": 1}
        wa, we = S.d3_freqs(1.0, 0.5)
        assert (round(wa ** 2, 6), round(we ** 2, 6)) == (1.0, 2.5)
        lev = S.flux_levels(2.0, -0.5, 0.2)
        assert abs(lev[1] - lev[-1] - math.sqrt(3) * math.sin(0.2)) < 1e-12
        assert sorted(S.q8_levels(1.0, 0.05, 0.25, 0.5).values()) == [1.0, 2.2, 2.6, 3.2, 4.0]
        assert S.frobenius_schur(S.Q8_CHARS["E"], S.Q8_CLASSES, S.Q8_SIZES, lambda g: S.q8_mul(g, g)) == -1
        v1 = S.q8_E_basis()[0][1]

        # 2. 畫面
        self.show_header("對稱群與耦合振盪子：D3、C3、Q8")
        c1 = np.array([-4.3, 0.3, 0])
        tri = O.osc_triangle(1.3, c1)
        t = ValueTracker(0.0)
        m = O.masses_updater(lambda: [v * math.cos(we * t.get_value()) for v in (1, -0.5, -0.5)], 1.3, c1)
        c2 = np.array([-0.4, 0.3, 0])
        w = S.OMEGA
        tri2 = O.osc_triangle(1.3, c2, labels=False)
        ph = O.phasors([1, w, w * w], O.tri_points(1.3, c2), 0.5)
        pos = O.q8_positions(0.8, 1.7, np.array([4.2, 0.3, 0]))
        graph = O.q8_graph(pos, S.q8_mul, font_size=22, node_radius=0.2)
        ph8 = O.phasors(list(v1), [pos[h] for h in S.Q8_NAMES], 0.45, zero_dot=False)
        assert_in_frame(tri, tri2, ph, graph)

        self.say("D3：", "$E=(2,-1,-1)$", " 的模式，", "$\\omega_E^2=2.5$")
        self.play(FadeIn(tri[1]), FadeIn(tri[2]), FadeIn(m), run_time=0.8)
        self.play(t.animate.set_value(6.0), run_time=4.0, rate_func=linear)
        self.say("C3 的標籤：", "$(1,\\omega,\\omega^2)$", "，破缺時裂距 ", "$\\sqrt3\\sin\\varphi$")
        self.play(FadeIn(tri2), LaggedStartMap(GrowArrow, ph, lag_ratio=0.3), run_time=1.5)
        self.say("Q8：", "$E$", " 出現兩次、四重簡併；", "$v_1$", " 的相量")
        self.play(FadeIn(graph), run_time=1.0)
        self.play(FadeIn(ph8), run_time=1.0)
        self.wait(2.0)
