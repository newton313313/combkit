"""
範例 36：不對稱陀螺的轉動能量（RE）曲面與量子化路徑（Harter《Principles of Symmetry, Dynamics, and Spectroscopy》第 7 章 §7.4A、
(7.4.3)–(7.4.5)、Fig. 7.4.1–7.4.2）——先用程式算：J = 10、A, B, C = 0.2, 0.4, 0.6 cm⁻¹ 的 21 個能階（2J+1 維對角化 ＝ (5.5.54) 的三對角矩陣），
分隔線 BJ(J+1) = 44 正好是一個本徵值；分隔線下面的能階成對（C2(x) 叢 (A1 B1)、(A2 B2)…）、上面也成對（C2(z) 叢 (B2 A1)、(B1 A2)…），
最低一對只差 8.6 × 10⁻⁷ cm⁻¹（26 kHz 的穿隧）。再在 RE 曲面 E(n̂) = J(J+1)(A n_x² + B n_y² + C n_z²) 上，於每個本徵能量畫等高線（自己的 marching squares），
這就是 Fig. 7.4.1(b) 的「量子化路徑」。畫：投影的 RE 曲面（半徑隨能量）、21 組路徑、鏡頭繞一圈，右邊是最低與最高兩對的能階。

在 examples/ 資料夾裡執行（這裡有 custom_config.yml）：
    manimgl ex36_re_surface.py RESurfaceExample
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))     # 讓 combkit 可以 import（已 pip install 就不需要）

import numpy as np
from manimlib import *

from combkit import sym_coup as SC
from combkit import sym_rot as RT
from combkit.draw import sym_draw as SD
from combkit.draw.fulton_draw import assert_apart, assert_in_frame
from combkit.draw.scene_base import ChuScene


def check():
    """21 個能階 ＝ 三對角公式；分隔線 44；成對的 D2 標記；最低一對的分裂；每個能量都有等高線。回傳 (能階, 路徑)。"""
    p = SC.FIG742
    lv = SC.rigid_levels(p["J"], p["A"], p["B"], p["C"])
    assert np.allclose(sorted(e for e, _ in lv), np.linalg.eigvalsh(RT.asym_matrix_formula(p["J"], p["A"], p["B"], p["C"])))
    f7 = SC.fig742()
    assert abs(f7["separatrix"] - 44) < 1e-12 and any(abs(e - 44) < 1e-9 for e, _ in lv)
    assert [(a[1], b[1]) for a, b in f7["low"][:2]] == [("A_1", "B_1"), ("A_2", "B_2")]
    assert [(a[1], b[1]) for a, b in f7["high"][:2]] == [("B_2", "A_1"), ("B_1", "A_2")]
    assert abs(f7["low_split"][0] - 8.6e-7) < 1e-8
    paths = SC.rigid_paths()
    assert len(paths) == 21 and all(len(ps) >= 1 for _, _, ps in paths)
    for e, _, ps in paths:
        for q in ps:
            assert np.abs(SC.re_rigid(q) - e).max() < 0.15                                      # 重取樣後仍在等高線上
    return lv, paths, f7


class RESurfaceExample(ChuScene):
    def construct(self):
        # 1. 數學
        lv, paths, f7 = check()

        # 2. 畫
        self.show_header("J = 10 不對稱陀螺的 RE 曲面")
        radius = lambda e: 1.2 + 1.0 * (e - 22.0) / 44.0
        cols = color_gradient([BLUE_B, TEAL, YELLOW, RED_B], len(paths))
        mesh = SD.RESurface.mesh_lines(lambda u: radius(SC.re_rigid(u)), 3, 4, 40)
        res = SD.RESurface([ps for _, _, ps in paths], [radius(e) for e, _, _ in paths], view=SD.View3D(-30, 65, center=LEFT * 2.6 + DOWN * 0.2),
                           colors=cols, mesh=mesh)
        low = [(round(a[0], 4), f"{a[1]}\\,{b[1]}") for a, b in f7["low"][:2]]
        high = [(round(b[0], 4), f"{a[1]}\\,{b[1]}") for a, b in f7["high"][:2]]
        tex = VGroup(*[Tex(rf"{e:.4f}\ \ ({lab})", font_size=26).set_color(c) for (e, lab), c in
                       zip(high + low[::-1], [RED_B, RED_B, BLUE_B, BLUE_B])]).arrange(DOWN, buff=0.28, aligned_edge=LEFT)
        sep = Tex(r"BJ(J+1)=44", font_size=26).set_color(GREY_A)
        side = VGroup(tex[:2], sep, tex[2:]).arrange(DOWN, buff=0.35, aligned_edge=LEFT).move_to(RIGHT * 4.3 + UP * 0.3)
        for m in (side,):
            assert_in_frame(m)
        self.play(FadeIn(res.mesh), FadeIn(res.axes), run_time=1.0)
        self.say("能量 ", "$E=J(J+1)(An_x^2+Bn_y^2+Cn_z^2)$", " 的徑向圖")
        self.play(*[ShowCreation(c) for c in res.curves], run_time=2.0)
        self.say("每個本徵能量一組等高線：21 條量子化路徑")
        self.play(FadeIn(side), run_time=1.0)
        self.say("分隔線下繞 ", "$\\pm x$", "、上繞 ", "$\\pm z$", "：能階兩兩成對")
        self.play(UpdateFromAlphaFunc(res, lambda m, a: m.set_view(-30 + 80 * a, 65)), run_time=3.0)
        self.say("最低一對只差 ", "$8.6\\times10^{-7}\\,\\mathrm{cm}^{-1}$", "：26 kHz 的穿隧")
        self.wait(1.0)
