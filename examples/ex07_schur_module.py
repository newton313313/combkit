"""
範例 7：Schur 模 E^λ 的權重（Fulton 第 8 章）——先用程式列出 λ=(2,1)、m=3 的 8 個半標準表並按權重（content）分堆，
再用 manim 排成權重圖；x1x2x3 那一堆有 2 個，全部加起來就是 Schur 多項式 s_(2,1)(x1,x2,x3)。

在 examples/ 資料夾裡執行（這裡有 custom_config.yml）：
    manimgl ex07_schur_module.py SchurModuleExample
"""
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))     # 讓 combkit 可以 import（已 pip install 就不需要）

from manimlib import *

from combkit.draw.fulton_draw import assert_in_frame, tableau_mob
from combkit.draw.scene_base import ChuScene
from combkit.tableaux import content, hook_content_count, schur_poly, ssyt


class SchurModuleExample(ChuScene):
    def construct(self):
        # 1. 數學：8 個半標準表 = E^(2,1) 的基底（m = 3），按權重分 7 堆
        lam, m = (2, 1), 3
        piles = defaultdict(list)
        for T in ssyt(lam, m):
            piles[content(T, m)].append(T)
        assert len(ssyt(lam, m)) == hook_content_count(lam, m) == 8 and len(piles) == 7
        assert {w: len(ts) for w, ts in piles.items()} == {w: c for w, c in schur_poly(lam, m).items()}

        # 2. 畫面：權重 (a,b,c) 投影到平面（a+b+c=3 的三角座標）
        self.show_header("Schur 模的權重")
        groups = VGroup()
        for w, ts in piles.items():
            pos = np.array([1.6 * (w[0] - w[1]), 0.7 * (w[0] + w[1] - 2 * w[2]) + 0.2, 0.0])
            tabs = VGroup(*(tableau_mob(T, cs=0.32) for T in ts)).arrange(RIGHT, buff=0.15).move_to(pos)
            mono = "".join(f"x_{i + 1}^{{{e}}}" if e > 1 else f"x_{i + 1}" for i, e in enumerate(w) if e)
            lab = Tex(("2" if len(ts) == 2 else "") + mono, font_size=24).next_to(tabs, DOWN, buff=0.1)
            groups.add(VGroup(tabs, lab))
        assert_in_frame(groups)
        self.say("對角矩陣 ", "$\\mathrm{diag}(x_1,x_2,x_3)$", " 把 ", "$e_T$", " 乘上 ", "$x^T$", "：8 個表分成 7 堆")
        self.play(LaggedStartMap(FadeIn, groups, lag_ratio=0.15), run_time=2.0)
        self.say("中間那堆有 2 個；全部加起來就是 ", "$s_{(2,1)}(x_1,x_2,x_3)$")
        self.wait(2.5)
