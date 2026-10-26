"""
範例 30：C_n 對稱的單擺環（Harter《Principles of Symmetry, Dynamics, and Spectroscopy》第 2 章 §2.6–2.7）——先用程式算：
C₃ 的投影算符「最小方程連乘」＝「特徵標加權平均」（𝒫^{0₃} = J/3，r𝒫^{k} = ε^k 𝒫^{k}）；n = 6 環（a = 1、b = 0）的 ω² = 0、1、3、4、3、1
（＝ 6×6 數值對角化）；n = 24、每隔一條彈簧變弱（ā = 3/2、a̲ = 1/2）裂成聲學帶 [0, 1] 與光學帶 [3, 4]（公式 ＝ 24×24 數值）。
最後畫出：6 擺環上的 k = 1 行波（六個相量鐘依序差 60°）、摺疊進 C₁₂ 區的 C₂₄ 頻譜 2 ∓ 2cos(k/2)（灰線）與對稱破缺後的兩條帶（點），
以及 C₃ 的特徵標表。

在 examples/ 資料夾裡執行（這裡有 custom_config.yml）：
    manimgl ex30_cn_ring.py CnRingExample
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))     # 讓 combkit 可以 import（已 pip install 就不需要）

import numpy as np
import sympy as sp
from manimlib import *

from combkit import sym_abel as SA
from combkit import sym_mat as SM
from combkit.draw import sym_draw as SD
from combkit.draw.fulton_draw import assert_in_frame
from combkit.draw.scene_base import ChuScene
from combkit.draw.styles import zh


class CnRingExample(ChuScene):
    def construct(self):
        # 1. 數學
        P = SA.cn_projectors_char(3)
        Q = SA.cn_projectors_minimal(3)
        R3 = SA.cyclic_rep(3)
        assert P[0] == sp.ones(3, 3) / 3
        for k in range(3):
            assert sp.simplify(sp.expand_complex(P[k] - Q[k])) == sp.zeros(3)
            assert sp.simplify(sp.expand_complex(R3 * P[k] - SA.eps(3, k) * P[k])) == sp.zeros(3)
        assert sp.simplify(sp.expand_complex(P[1] - SM.projectors(R3)[sp.nsimplify(SA.eps(3, 1))])) == sp.zeros(3)
        w6 = SA.ring_omega2(6, 1, 0)
        assert w6 == [0, 1, 3, 4, 3, 1]
        assert np.allclose(sorted(np.linalg.eigvalsh(np.array(SA.ring_matrix(6), dtype=float))), sorted(map(float, w6)))
        ks = [2 * np.pi * k / 12 for k in SA.brillouin(12)]
        bands = np.array([SA.broken_bands(k) for k in ks])
        assert np.allclose(np.sort(bands.ravel()), np.sort(np.linalg.eigvalsh(SA.broken_ring_matrix())), atol=1e-12)
        assert np.allclose([SA.broken_bands(k, 1, 1) for k in ks], [(2 - 2 * abs(np.cos(k / 2)), 2 + 2 * abs(np.cos(k / 2))) for k in ks])
        assert abs(bands[:, 0].max() - 1) < 1e-12 and abs(bands[:, 1].min() - 3) < 1e-12 and abs(bands[:, 1].max() - 4) < 1e-12

        # 2. 畫
        self.show_header("C_n 對稱的單擺環")
        ring = SD.PendulumRing(6, radius=1.5, unit=0.45, labels=True).move_to(LEFT * 4.3 + UP * 0.6)
        clocks = VGroup(*[SD.PhasorClock(0.34, label_tex=f"x_{p + 1}", font_size=20) for p in range(6)])
        clocks.arrange(RIGHT, buff=0.16).move_to(LEFT * 4.0 + DOWN * 2.0)
        ax = SD.plot_axes((-3.3, 3.3, 1), (0, 4.5, 1), width=5.2, height=3.0, x_label="k", y_label=r"\omega^2")
        ax.move_to(RIGHT * 2.9 + UP * 0.9)
        free = VGroup(*[SD.band_curve(ax, (lambda k, sg=sg: 2 + sg * 2 * np.cos(k / 2)), color=GREY_B, stroke_width=2)
                        for sg in (-1, 1)])                  # ā = a̲ = 1 時摺疊進 C₁₂ 區的 C₂₄ 頻譜：2 ∓ 2cos(k/2)
        ac = SD.band_dots(ax, ks, bands[:, 0], color=SD.WEAK_COLOR)
        op = SD.band_dots(ax, ks, bands[:, 1], color=SD.STRONG_COLOR)
        cap = zh("C24 → C12：聲學帶 [0, 1]、光學帶 [3, 4]", font_size=22).next_to(ax, DOWN, buff=0.3)
        tab = SD.tex_table([["", "r^0", "r^1", "r^2"], ["0_3", "1", "1", "1"],
                            ["1_3", "1", r"\varepsilon", r"\varepsilon^*"], ["2_3", "1", r"\varepsilon^*", r"\varepsilon"]],
                           [0.7, 0.6, 0.6, 0.6], row_h=0.36, font_size=22).move_to(RIGHT * 3.0 + DOWN * 2.05)
        for m in (ring, clocks, ax, cap, tab):
            assert_in_frame(m)
        self.add(ring, clocks)
        self.play(FadeIn(ax), ShowCreation(free), run_time=1.0)
        self.play(FadeIn(ac), FadeIn(op), FadeIn(cap), FadeIn(tab), run_time=1.0)
        v = SA.wave_vector(6, 1) * np.sqrt(6)                         # 分量 e^{−iπp/3}
        w = float(sp.sqrt(w6[1]))
        tr = ValueTracker(0.0)
        ring.add_updater(lambda m: m.set_displacements((np.exp(-1j * w * tr.get_value()) * v).real))
        clocks.add_updater(lambda m: [c.set_value(np.exp(-1j * w * tr.get_value()) * z) for c, z in zip(m, v)])
        self.say("k = 1 的行波：六個相量鐘依序差 60°，一起順時針轉")
        self.play(tr.animate.set_value(2 * np.pi), run_time=5, rate_func=linear)
        ring.clear_updaters()
        clocks.clear_updaters()
        self.wait(0.8)
