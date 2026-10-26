"""
範例 31：由生成元算特徵標表（Harter《Principles of Symmetry, Dynamics, and Spectroscopy》第 3 章 §3.2、§3.5）——先用程式算：
只給兩個生成元的 2×2 矩陣（C3v：轉 120° 與一面鏡子；D4：轉 90° 與一面鏡子），`FiniteGroup.from_generators` 做封閉（6、8 個元素），
算共軛類、類代數的結構常數、類算子的最小方程（C3v 的 c³ − 9c；D4 要兩個 c³ − 4c 才分得完）、全對易冪等元，
再由 ℓ = (°G p₁)^{1/2}、χ* = p°G/ℓ 讀出特徵標表（＝不可約表示矩陣的跡，＝書上 (3.5.8)、(3.6.3)，差一個重新命名）；
最後畫出兩張表，並把 E 列的維度 2 與 Σℓ² = °G 標出來。

在 examples/ 資料夾裡執行（這裡有 custom_config.yml）：
    manimgl ex31_character_table.py CharacterTableExample
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))     # 讓 combkit 可以 import（已 pip install 就不需要）

import sympy as sp
from manimlib import *

from combkit import sym_nonabel as SN
from combkit.draw import sym_draw as SD
from combkit.draw.fulton_draw import assert_apart, assert_in_frame
from combkit.draw.scene_base import ChuScene
from combkit.draw.styles import zh


def canon(G, T):
    """特徵標表的標準形：行依 (類大小, 元素的階) 排、列排序——比較兩個群的表（差一個重新命名）。"""
    sig = [(len(C), G.order_of(C[0])) for C in G.classes]
    cols = sorted(range(len(sig)), key=lambda k: sig[k])
    return sorted(tuple(T[a, c] for c in cols) for a in range(T.rows))


def build(gens, book_group):
    """生成元 → (群, ℓ, 特徵標表, 類的 TeX 表頭)；並核對：類代數 ＝ 跡 ＝ 書上的表、Σℓ² = °G、列正交。"""
    G = SN.FiniteGroup.from_generators(gens)
    ells, T = G.char_table_class_algebra()
    assert sum(int(l_) ** 2 for l_ in ells) == G.order
    assert G.char_table_traces(G.irreps(sub_gen="s")) == T                        # 跡（不可約表示由基本算子求得）
    assert canon(G, T) == canon(book_group, book_group.char_table_class_algebra()[1])   # ＝ 書上的表
    sz = G.class_sizes()
    orth = sp.Matrix(T.rows, T.rows, lambda a, b: sum(sz[k] * T[a, k] * sp.conjugate(T[b, k]) for k in range(T.cols)))
    assert orth == G.order * sp.eye(T.rows)
    heads = [r"\{" + ",".join(G.tex[g] for g in C) + r"\}" for C in G.classes]
    return G, ells, T, heads


class CharacterTableExample(ChuScene):
    def construct(self):
        # 1. 數學
        G3, ell3, T3, h3 = build({"r": SN.rot2(120), "s": SN.refl2(0)}, SN.C3V)
        G4, ell4, T4, h4 = build({"R": SN.rot2(90), "s": SN.refl2(0)}, SN.D4)
        assert G3.order == 6 and G4.order == 8 and ell3 == [1, 1, 2] and ell4 == [1, 1, 1, 1, 2]
        c = sp.Symbol("c")
        k3 = next(i for i, C in enumerate(G3.classes) if len(C) == 3)
        assert G3.min_poly_class(k3).as_expr() == c ** 3 - 9 * c                     # (3.2.12)–(3.2.13)
        polys4 = [G4.min_poly_class(i).as_expr() for i, C in enumerate(G4.classes) if len(C) == 2]
        assert polys4 == [c ** 3 - 4 * c] * 3 and all(len(G4.class_idempotents_of(i)) == 3 for i in range(1, 5) if len(G4.classes[i]) == 2)

        # 2. 畫
        self.show_header("由生成元算特徵標表")
        names3 = [r"A_1", r"A_2", r"E"]
        names4 = [r"\chi^1", r"\chi^2", r"\chi^3", r"\chi^4", r"E"]
        t3 = SD.sym_table(h3, names3, T3.tolist(), cell_w=1.25, row_h=0.5, font_size=24, first_w=0.8)
        t4 = SD.sym_table(h4, names4, T4.tolist(), cell_w=1.05, row_h=0.46, font_size=20, first_w=0.8)
        t3.move_to(LEFT * 3.6 + UP * 0.6)
        t4.move_to(RIGHT * 2.9 + UP * 0.6)
        lab3 = zh("C3v：生成元 r（120°）與 s（鏡子），6 個元素", font_size=22).next_to(t3, UP, buff=0.3)
        lab4 = zh("D4：生成元 R（90°）與 s，8 個元素", font_size=22).next_to(t4, UP, buff=0.3)
        eq3 = Tex(r"c^3-9c=0,\quad 1^2+1^2+2^2=6", font_size=26).next_to(t3, DOWN, buff=0.35)
        eq4 = Tex(r"c^3-4c=0\ (\times2),\quad 1+1+1+1+2^2=8", font_size=26).next_to(t4, DOWN, buff=0.35)
        for m in (t3, t4, lab3, lab4, eq3, eq4):
            assert_in_frame(m)
        assert_apart(t3, t4)
        self.play(FadeIn(lab3), FadeIn(t3), run_time=1.0)
        self.play(FadeIn(lab4), FadeIn(t4), run_time=1.0)
        self.play(FadeIn(eq3), FadeIn(eq4), run_time=0.8)
        rows = VGroup(*[SD.table_cell(t3, 3, j, len(h3) + 1) for j in range(len(h3) + 1)],
                      *[SD.table_cell(t4, 5, j, len(h4) + 1) for j in range(len(h4) + 1)])
        self.say("只用類代數：", "$\\ell=({}^\\circ G\\,p_1)^{1/2}$", "，E 的維度是 2")
        self.play(rows.animate.set_color(SD.KET_COLOR), run_time=0.8)
        self.wait(1.0)
