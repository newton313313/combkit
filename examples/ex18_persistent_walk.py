"""
範例 18：有慣性的漫步與電報方程（Rudnick & Gaspari《Elements of the Random Walk》第 5 章）——先用程式算出
「和上一步同方向」的機率 Q = 3/4 與 1/4 走 20 步的精確分布（左右兩個狀態的逐步遞迴＝暴力列舉＝轉移矩陣）、
精確變異數 vs 大 N 的 NQ/(1−Q)、持續長度 Q/(1−Q)，以及電報方程（v = γ = 1）更正後的 (5.61)：連續部分＋兩個尖峰、總數守恆，再畫出：
兩個醉漢的分布長條＋(5.27) 的高斯（格子上 ×2）、t = 1、3、8 的 Σ(x, t) 與尖峰重量 e^{−t}/2。

在 examples/ 資料夾裡執行（這裡有 custom_config.yml）：
    manimgl ex18_persistent_walk.py PersistentWalkExample
"""
import sys
from fractions import Fraction
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))     # 讓 combkit 可以 import（已 pip install 就不需要）

from manimlib import *

from combkit import walk_var as V
from combkit.draw import walk_draw as D
from combkit.draw.fulton_draw import assert_in_frame
from combkit.draw.scene_base import ChuScene
from combkit.draw.styles import zh

N = 20
QS = (Fraction(3, 4), Fraction(1, 4))
TIMES = (1, 3, 8)


class PersistentWalkExample(ChuScene):
    def construct(self):
        # 1. 數學
        dists = {}
        for Q in QS:
            P = V.persistent_pmf(N, Q)                                              # Fraction，精確
            assert sum(P.values()) == 1 and P == V.persistent_pmf(N, Q)
            assert V.persistent_var(N, Q) == V.persistent_var_formula(N, Q)
            dists[Q] = {x: float(w) for x, w in P.items()}
        for n in range(9):                                                          # 逐步遞迴 = 暴力列舉
            assert V.persistent_pmf(n, QS[0]) == V._persistent_brute(n, QS[0])
        assert abs(V.total_fourier(N, 0.3, 0.75) - V.total_eig(N, 0.3, 0.75)) < 1e-13   # = 轉移矩陣的特徵值
        var = [float(V.persistent_var(N, Q)) for Q in QS]
        asym = [float(V.persistent_var_asym(N, Q)) for Q in QS]
        assert [f"{v:.2f}" for v in var] == ["56.00", "7.11"] and [f"{v:.2f}" for v in asym] == ["60.00", "6.67"]
        assert V.persistence_length(QS[0]) == 3 and V.persistence_length(QS[1]) == Fraction(1, 3)
        spikes = [V.tele_spike(t) for t in TIMES]
        assert [f"{s:.4g}" for s in spikes] == ["0.1839", "0.02489", "0.0001677"]
        assert all(abs(V.tele_mass(t) - 1) < 1e-9 for t in TIMES)                   # 連續部分 + 尖峰 = 1
        assert [f"{float(V.tele_cont(0, t)):.4f}" for t in TIMES] == ["0.3368", "0.2199", "0.1388"]

        # 2. 畫面
        self.show_header("有慣性的漫步與電報方程")
        ux, uy = 0.12, 6.0
        origins = [np.array([-3.6, 0.4, 0]), np.array([-3.6, -2.6, 0])]
        colors = [D.BAR_COLOR, D.TAIL_COLOR]
        charts, curves, labels = VGroup(), VGroup(), VGroup()
        for Q, o, c, v_, a_ in zip(QS, origins, colors, var, asym):
            charts.add(D.bar_chart(dists[Q], ux, uy, o, color=c))
            curves.add(D.curve(lambda x, Q=Q: 2 * float(V.gauss_persistent(x, N, float(Q))), -20, 20, ux, uy, o))
            lp = V.persistence_length(Q)
            lp_tex = str(lp.numerator) if lp.denominator == 1 else f"{lp.numerator}/{lp.denominator}"
            lab = Tex(f"Q={Q.numerator}/{Q.denominator}:\\ \\mathrm{{var}}={v_:.2f}\\ (NQ/(1-Q)={a_:.2f}),\\ l_p={lp_tex}",
                      font_size=22)
            labels.add(lab.move_to(o + 2.2 * UP))

        ax = D.plot_axes((0, 9, 1), (0, 0.4, 0.1), width=5.2, height=3.6, x_ticks=(0, 3, 6, 9), y_ticks=(0.1, 0.2, 0.3),
                         x_label="x", y_label="\\Sigma", font_size=18)
        ax.move_to(np.array([3.6, 0.0, 0]))
        tele, arrows = VGroup(), VGroup()
        for t, s, col in zip(TIMES, spikes, (YELLOW, ORANGE, RED)):
            g = ax.get_graph(lambda x, t=t: float(V.tele_cont(min(x, t * (1 - 1e-9)), t)), x_range=(0, t * (1 - 1e-6), t / 120))
            tele.add(g.set_stroke(col, 3))
            top, side = {1: (0.36, UP), 3: (0.16, UP), 8: (0.10, LEFT)}[t]
            arr = Arrow(ax.c2p(t, 0), ax.c2p(t, top), buff=0).set_color(col)
            arrows.add(VGroup(arr, Tex(f"t={t}:\\ {s:.4g}", font_size=20).set_color(col)
                              .next_to(arr.get_end(), side, buff=0.08)))
        t_lab = zh("電報方程 v = γ = 1：分布＋尖峰重量", font_size=22).next_to(ax, UP, buff=0.25)
        assert_in_frame(charts, curves, labels, ax, tele, arrows, t_lab)

        self.say("記得上一步：", "$Q=3/4$", " 越走越寬、", "$Q=1/4$", " 原地打轉")
        self.play(*[LaggedStartMap(FadeIn, ch[0], lag_ratio=0.02) for ch in charts], FadeIn(labels), run_time=1.5)
        self.play(*[ShowCreation(c) for c in curves], run_time=1.2)
        self.say("電報方程：兩個尖峰以 ", "$e^{-\\gamma t}$", " 變小，後面變成高斯")
        self.play(FadeIn(ax), FadeIn(t_lab), run_time=0.8)
        self.play(*[ShowCreation(g) for g in tele], LaggedStartMap(FadeIn, arrows, lag_ratio=0.3), run_time=2.0)
        self.wait(2.0)
