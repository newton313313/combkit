"""
範例 32：由轉動矩陣產生八面體群 O，算特徵標表與一張相關表（Harter《Principles of Symmetry, Dynamics, and Spectroscopy》第 4 章 §4.1–§4.2）——
先用程式算：只給兩個 3×3 生成元（繞 x 軸、y 軸各轉 90°），`FiniteGroup.from_generators` 封閉出 24 個轉動；共軛類的大小 1, 8, 3, 6, 6；
繞稜的 180° 那一類的類算子 c_i 的最小方程 c⁵ − 40c³ + 144c = (c+6)(c−6)(c+2)(c−2)c（五個根 ⇒ 五個不可約表示）；
類代數算出的特徵標表 ＝ 書上的 (4.1.11)（差一個重新命名）＝ `combkit.sym_oct.O_CHARS`；四條體對角線的置換讓 O ≅ S₄；
再用 `sym_oct.correlation("D4")` 讀出 O ↓ D4 的相關表 (4.2.42a)，並用 Frobenius 互反（誘導表示的跡）核對它的每一行。
最後畫兩張表與一個轉動中的立方體（四條體對角線上色）。

在 examples/ 資料夾裡執行（這裡有 custom_config.yml）：
    manimgl ex32_octahedral_group.py OctahedralGroupExample
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))     # 讓 combkit 可以 import（已 pip install 就不需要）

import numpy as np
import sympy as sp
from manimlib import *

from combkit import sym_nonabel as SN
from combkit import sym_oct as SO
from combkit.draw import sym_draw as SD
from combkit.draw.fulton_draw import assert_apart, assert_in_frame
from combkit.draw.scene_base import ChuScene
from combkit.draw.styles import zh

CLASS_TEX = {(1, 1): r"\mathbf 1", (8, 3): r"8C_3", (3, 2): r"3C_2", (6, 4): r"6C_4", (6, 2): r"6C_2'"}


def build():
    """生成元 → (群, 依書上類序 1, r, R², R, i 排好的特徵標表, 表頭)；並核對類大小、最小方程、表 ＝ (4.1.11)、O ≅ S4。"""
    G = SN.FiniteGroup.from_generators({"R": SN.rot3((1, 0, 0), 90), "S": SN.rot3((0, 1, 0), 90)})
    assert G.order == 24
    sig = [(len(C), G.order_of(C[0])) for C in G.classes]
    order = [sig.index(k) for k in [(1, 1), (8, 3), (3, 2), (6, 4), (6, 2)]]       # 書上的類序
    assert sorted(len(C) for C in G.classes) == [1, 3, 6, 6, 8]
    c = sp.Symbol("c")
    ki = order[4]
    assert G.min_poly_class(ki).as_expr() == c ** 5 - 40 * c ** 3 + 144 * c        # (4.1.8)
    ells, T = G.char_table_class_algebra()
    Tb = sp.Matrix([[T[a, k] for k in order] for a in range(T.rows)])
    perm = SN.reorder(Tb, SO.O_CHARS_BOOK)                                            # 列換成書上的 A1 A2 E T1 T2
    Tb = sp.Matrix([list(Tb.row(a)) for a in perm])
    assert Tb == SO.O_CHARS and sum(int(e) ** 2 for e in ells) == 24
    # O ≅ S4：每個轉動對四條體對角線的置換都不同
    perms = set()
    for g in G.names:
        M = np.array(G.mats[g].evalf(), dtype=float)
        p = []
        for d in SO.DIAGONALS.values():
            v = M @ np.array(d, float)
            p.append(next(k for k, e in SO.DIAGONALS.items() if np.allclose(v, e) or np.allclose(v, -np.array(e))))
        perms.add(tuple(p))
    assert len(perms) == 24
    heads = [CLASS_TEX[sig[k]] for k in order]
    return G, Tb, heads


def correlation_d4():
    """O ↓ D4（(4.2.42a)）＝ 相關表的特徵標公式 ＝ Frobenius：每個 D4 irrep 誘導到 O 的分解。"""
    corr = SO.correlation("D4")
    assert corr == sp.Matrix(SO.CORR_BOOK["D4"])
    for lab, (f_mat, f_form, col) in SO.frobenius("D4").items():
        assert f_mat == f_form == col
    return corr


class OctahedralGroupExample(ChuScene):
    def construct(self):
        # 1. 數學
        G, Tb, heads = build()
        corr = correlation_d4()

        # 2. 畫
        self.show_header("八面體群 O 的特徵標表與相關表")
        names = [r"A_1", r"A_2", r"E", r"T_1", r"T_2"]
        t1 = SD.sym_table(heads, names, Tb.tolist(), cell_w=0.85, row_h=0.46, font_size=22, first_w=0.7)
        sub = SO.SUB["D4"][2]
        rows = [[("1" if x else r"\cdot") for x in r] for r in corr.tolist()]
        t2 = SD.sym_table(sub, names, rows, cell_w=0.62, row_h=0.46, font_size=22, first_w=0.7, corner=r"\downarrow D_4")
        t1.move_to(LEFT * 0.6 + UP * 0.9)
        t2.move_to(RIGHT * 4.6 + UP * 0.9)
        cube = SD.CubeView(size=0.85, view=SD.View3D(-25, 70, center=LEFT * 5.2 + UP * 0.6))
        lab1 = zh("生成元：繞 x、y 各轉 90°，封閉出 24 個轉動", font_size=22).next_to(t1, UP, buff=0.3)
        lab2 = zh("O ↓ D4（4.2.42a）", font_size=22).next_to(t2, UP, buff=0.3)
        eq = Tex(r"c_i^5-40c_i^3+144c_i=(c_i+6)(c_i-6)(c_i+2)(c_i-2)c_i=0", font_size=28).move_to(DOWN * 1.6)
        for m in (t1, t2, cube, lab1, lab2, eq):
            assert_in_frame(m)
        assert_apart(t1, t2)
        assert_apart(cube, t1)
        self.play(FadeIn(cube), FadeIn(lab1), FadeIn(t1), run_time=1.2)
        ax, deg, _ = SO.o_turn("r1")
        t = ValueTracker(0)
        cube.add_updater(lambda m: m.set_rotation(SD.rot_matrix(ax, t.get_value() * deg)))
        self.play(t.animate.set_value(1), run_time=1.5)
        cube.clear_updaters()
        self.play(FadeIn(eq), run_time=0.8)
        self.play(FadeIn(lab2), FadeIn(t2), run_time=1.0)
        hl = VGroup(*[SD.table_cell(t2, 4, j, len(sub) + 1) for j in range(len(sub) + 1)])
        self.say("$T_1\\downarrow D_4=A_2\\oplus E$", "：Frobenius 讀行就是誘導表示")
        self.play(hl.animate.set_color(SD.KET_COLOR), run_time=0.8)
        self.wait(1.0)
