"""
範例 35：投影算子算出的 Clebsch–Gordan 係數，把 T1u ⊗ T1u 約化成 A1g ⊕ Eg ⊕ T1g ⊕ T2g（Harter《Principles of Symmetry, Dynamics, and Spectroscopy》
第 6 章 §6.1B、(6.1.14)–(6.1.19)）——先用程式算：八面體群 O_h 的特徵標給出頻率 1, 1, 1, 1；投影公式 (6.1.18b)（每一份取 N_{nkl} ≠ 0 的第一組，
C_{kl,n} = √N > 0）算出 9 × 9 的耦合係數表，與書上 (6.1.19c) 逐格相同；C†(𝒟⊗𝒟)(R1)C 是四個方塊 1 ⊕ 𝒟^{Eg} ⊕ 𝒟^{T1g} ⊕ 𝒟^{T2g}（48 個元素都對）。
畫：左邊 9 × 9 的直積矩陣 𝒟^{T1u}(R1) ⊗ 𝒟^{T1u}(R1) 的非零格子，右邊乘上 C 之後的方塊對角，下面是耦合係數表的前三行。

在 examples/ 資料夾裡執行（這裡有 custom_config.yml）：
    manimgl ex35_clebsch_gordan.py ClebschGordanExample
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))     # 讓 combkit 可以 import（已 pip install 就不需要）

import numpy as np
from manimlib import *

from combkit import sym_prod as SP
from combkit.draw import sym_draw as SD
from combkit.draw.fulton_draw import assert_apart, assert_in_frame
from combkit.draw.scene_base import ChuScene


def tex_entry(x) -> str:
    """耦合係數 → TeX（±1/√2、1/√3、±1/√6、−2/√6；0 留空）。"""
    if abs(x) < 1e-12:
        return r"\cdot"
    table = {1 / np.sqrt(2): r"\tfrac{1}{\sqrt2}", 1 / np.sqrt(3): r"\tfrac{1}{\sqrt3}", 1 / np.sqrt(6): r"\tfrac{1}{\sqrt6}",
             2 / np.sqrt(6): r"\tfrac{2}{\sqrt6}"}
    for v, t in table.items():
        if abs(abs(x) - v) < 1e-12:
            return ("-" if x < 0 else "") + t
    raise ValueError(x)


def check():
    """頻率、投影法 ＝ 書上的表、約化方程、么正。回傳 (C, cols, 直積矩陣, 方塊)。"""
    S = SP.ohset("nonkosher")
    assert SP.diag_frequencies("T_{1u}", "T_{1u}") == {"A_{1g}": 1, "E_g": 1, "T_{1g}": 1, "T_{2g}": 1}
    C, cols = SP.cg_book("619c")
    assert SP.compare_table(C, cols, SP.CG_619C_BOOK)["equal"]
    assert SP.cg_reduces(S, "T_{1u}", "T_{1u}", C, cols)
    assert np.allclose(C.conj().T @ C, np.eye(9))
    K = np.kron(S.D["T_{1u}"]["R1"], S.D["T_{1u}"]["R1"]).real
    B = (C.conj().T @ K @ C).real
    assert np.allclose(B[1:3, 1:3], [[-0.5, -np.sqrt(3) / 2], [-np.sqrt(3) / 2, 0.5]])          # 𝒟^{Eg}(R1)
    assert np.allclose(B[6:, 6:], [[-1, 0, 0], [0, 0, 1], [0, -1, 0]])                          # 非 kosher 的 𝒟^{T2g}(R1)
    assert [SP.cg_norms(S, "T_{1u}", "T_{1u}", "T_{1g}")[k][1] for k in range(6)] == [0, 0, 0, 0, 0, 0.5]
    return C.real, cols, K, B


class ClebschGordanExample(ChuScene):
    def construct(self):
        # 1. 數學
        C, cols, K, B = check()

        # 2. 畫
        self.show_header("T1u 乘 T1u 的 Clebsch–Gordan 係數")
        left = SD.perm_grid((np.abs(K) > 1e-9).astype(int), cell=0.26, on_color=SD.KET_COLOR)
        left.move_to(LEFT * 4.3 + UP * 1.2)
        right = SD.perm_grid((np.abs(B) > 1e-9).astype(int), cell=0.26, on_color=SD.P2_COLOR)
        right.move_to(RIGHT * 0.2 + UP * 1.2)
        lab_l = Tex(r"\mathcal D^{T_{1u}}(R_1)\otimes\mathcal D^{T_{1u}}(R_1)", font_size=26).next_to(left, UP, buff=0.2)
        lab_r = Tex(r"\mathcal C^\dagger(\mathcal D\otimes\mathcal D)\,\mathcal C=1\oplus E_g\oplus T_{1g}\oplus T_{2g}", font_size=26).next_to(right, UP, buff=0.2)
        arrow = Arrow(left.get_right(), right.get_left(), buff=0.3).set_color(GREY_B)
        heads = [r"A_{1g}", r"E_g\,1", r"E_g\,2"]
        rows = [[tex_entry(C[r, k]) for k in range(3)] for r in range(9)]
        tab = SD.sym_table(heads, [f"{i}{j}" for i in (1, 2, 3) for j in (1, 2, 3)], rows, cell_w=0.9, row_h=0.32, font_size=20)
        tab.move_to(RIGHT * 4.9 + DOWN * 0.4)
        for m in (left, right, lab_l, lab_r, tab):
            assert_in_frame(m)
        assert_apart(left, right)
        assert_apart(right, tab)
        self.play(FadeIn(left), FadeIn(lab_l), run_time=1.0)
        self.say("直積矩陣：9 × 9，每列只有一個非零格")
        self.play(GrowArrow(arrow), FadeIn(right), FadeIn(lab_r), run_time=1.2)
        self.say("乘上耦合係數 ", "$\\mathcal C$", "：方塊大小 1、2、3、3")
        self.play(FadeIn(tab), run_time=1.0)
        self.say("前三行：", "$A_{1g}$", " 是 ", "$(11+22+33)/\\sqrt3$")
        self.wait(1.0)
