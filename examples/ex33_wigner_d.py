"""
範例 33：Wigner 的轉動矩陣 𝒟^j 兩種算法，與 Stern–Gerlach 的機率長條（Harter《Principles of Symmetry, Dynamics, and Spectroscopy》第 5 章 §5.4D–§5.5B）——
先用程式算：d^j_{mn}(β) 用 Wigner 的階乘公式 (5.4.45) 與矩陣指數 e^{−iβJ_y}（J 由升降算子 (5.4.23) 組出）各算一次，j = ½…4 都一致、
j = 20 誤差 < 10⁻⁹；第三種算法是 Schwinger 的兩個振子（把 (c a↑ + s a↓)^{j+n}(−s a↑ + c a↓)^{j−n} 展開）。
𝒟^{1/2} 轉 360° 是 −1（射線表示）；書上 (5.4.35) 的 𝒟² 有四格印錯：照印的矩陣在 β = 0.7 時 (𝒟𝒟ᵀ)₁₁ = 1.3174（不正交）。
Stern–Gerlach：j = 2、m = 1 的態轉 β 之後進未轉的分析器，五個通道的機率 |𝒟²_{m'1}(β)|²；β = 90° 時 0.25, 0.25, 0, 0.25, 0.25，
平均 ⟨J_z⟩ = cos β（(5.5.18)）。畫一個隨 β 傾斜的向量錐與五根機率長條。

在 examples/ 資料夾裡執行（這裡有 custom_config.yml）：
    manimgl ex33_wigner_d.py WignerDExample
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))     # 讓 combkit 可以 import（已 pip install 就不需要）

import numpy as np
from manimlib import *

from combkit import sym_rot as RT
from combkit.draw import sym_draw as SD
from combkit.draw.fulton_draw import assert_apart, assert_in_frame
from combkit.draw.scene_base import ChuScene
from combkit.draw.styles import zh

BETAS = np.linspace(0, 180, 61)


def check():
    """兩種（三種）算法一致、360° 的號、照印 (5.4.35) 不正交、Stern–Gerlach 的數字。回傳 (β 網格上的機率表, 卡片字串)。"""
    for j in (0.5, 1, 1.5, 2, 3, 4):
        for b in (23.0, 71.0, 140.0):
            assert np.allclose(RT.wigner_d(j, b), RT.D_expm(j, 0, b, 0).real, atol=1e-12)
            assert np.allclose(RT.wigner_d(j, b), RT.D_schwinger(j, b), atol=1e-12)
    assert np.allclose(RT.wigner_d(20, 60), RT.D_expm(20, 0, 60, 0).real, atol=1e-9)
    assert np.allclose(RT.wigner_D(0.5, 360, 0, 0), -np.eye(2)) and np.allclose(RT.wigner_D(1, 360, 0, 0), np.eye(3))
    assert f"{RT.d2_book_orthogonality(0.7):.4f}" == "1.3174"
    probs = np.array([RT.sg_probs(2, 1, b) for b in BETAS])
    assert np.allclose(probs.sum(axis=1), 1)
    ms = np.array([2, 1, 0, -1, -2])
    assert np.allclose(probs @ ms, np.cos(np.radians(BETAS)))                   # ⟨J_z⟩ = m cos β
    p90 = RT.sg_probs(2, 1, 90)
    card = " ".join(f"{x:.2f}" for x in p90)
    assert card == "0.25 0.25 0.00 0.25 0.25" and np.allclose(RT.sg_probs(2, 1, 180), [0, 0, 0, 1, 0])
    return probs, card


class WignerDExample(ChuScene):
    def construct(self):
        # 1. 數學
        probs, card = check()

        # 2. 畫
        self.show_header("Wigner 的 D 矩陣與 Stern–Gerlach")
        d1 = SD.labeled_mat(r"d^1(60^\circ)", [[f"{x:.4f}" for x in r] for r in RT.wigner_d(1, 60)], font_size=26)
        d1.move_to(LEFT * 3.4 + UP * 1.9)
        eq = Tex(r"\text{Wigner}=e^{-i\beta J_y}=\text{Schwinger},\quad \mathcal D^{1/2}(360^\circ)=-\mathbf 1", font_size=26)
        eq.move_to(LEFT * 3.0 + UP * 0.5)
        cone = SD.VectorCone(2, 1, unit=0.55, view=SD.View3D(-20, 75, center=LEFT * 4.2 + DOWN * 1.4))
        bars = SD.ProbBars(probs[0], labels=["2", "1", "0", "-1", "-2"], width=0.42, gap=0.16, height=2.4)
        bars.shift(RIGHT * 2.0 + DOWN * 2.2)
        lab = zh("轉過 β 的 |2, 1⟩ 進未轉的分析器", font_size=22).next_to(bars, UP, buff=2.7)
        beta_lab = DecimalNumber(0, num_decimal_places=0, font_size=28).next_to(bars, RIGHT, buff=0.5)
        beta_name = Tex(r"\beta=", font_size=28).next_to(beta_lab, LEFT, buff=0.08)
        for m in (d1, eq, cone, bars, lab):
            assert_in_frame(m)
        assert_apart(d1, bars)
        assert_apart(cone, bars)
        self.play(FadeIn(d1), FadeIn(eq), run_time=1.0)
        self.play(FadeIn(cone), FadeIn(bars), FadeIn(lab), FadeIn(beta_lab), FadeIn(beta_name), run_time=1.0)
        t = ValueTracker(0)

        def upd_bars(m):
            k = int(round(t.get_value() * (len(BETAS) - 1)))
            m.set_values(probs[k])

        cone.add_updater(lambda m: m.set_beta(180 * t.get_value()))
        bars.add_updater(upd_bars)
        beta_lab.add_updater(lambda m: m.set_value(180 * t.get_value()))
        self.play(t.animate.set_value(0.5), run_time=2.5, rate_func=linear)
        for m in (cone, bars, beta_lab):
            m.clear_updaters()
        self.say("$\\beta=90^\\circ$", "：機率 " + card.replace(" ", "、"))
        self.wait(0.5)
        for m, f in ((cone, lambda m: m.set_beta(180 * t.get_value())), (bars, upd_bars), (beta_lab, lambda m: m.set_value(180 * t.get_value()))):
            m.add_updater(f)
        self.play(t.animate.set_value(1.0), run_time=2.5, rate_func=linear)
        for m in (cone, bars, beta_lab):
            m.clear_updaters()
        self.say("$\\beta=180^\\circ$", "：全部進 ", "$m'=-1$", " 通道")
        self.wait(1.0)
