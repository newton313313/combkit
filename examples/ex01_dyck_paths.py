"""
範例 1：先用程式算出 n=3 的所有 Dyck 路徑，再用 manim 畫出來。

在 examples/ 資料夾裡執行（這裡有 custom_config.yml）：
    manimgl ex01_dyck_paths.py DyckPathsExample
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))     # 讓 combkit 可以 import（已 pip install 就不需要）

from manimlib import *

from combkit.counting import binary_trees, catalan, dyck_words, tree_to_dyck
from combkit.draw.catalan_draw import dyck_item
from combkit.draw.scene_base import ChuScene

N = 3


class DyckPathsExample(ChuScene):
    def construct(self):
        # 1. 數學：由樹產生 Dyck 字，並和暴力列舉比對
        trees = binary_trees(N)
        words = [tree_to_dyck(t) for t in trees]
        assert len(words) == catalan(N) and sorted(words) == dyck_words(N)

        # 2. 畫面：同一棵樹交給畫圖工具，依根／左／右上色
        self.show_header(f"半長 {N} 的 Dyck 路徑")
        paths = VGroup(*[dyck_item(t, scale=2.0, colored=True) for t in trees]).arrange(RIGHT, buff=0.6,
                                                                                           aligned_edge=DOWN)
        labels = VGroup(*[Text(w, font_size=22).next_to(p, DOWN, buff=0.3) for w, p in zip(words, paths)])
        self.say("程式先列舉出所有路徑，manim 再把它們畫出來：共 ", f"$C_{N}={catalan(N)}$", " 條")
        self.play(LaggedStartMap(ShowCreation, paths, lag_ratio=0.2), run_time=2)
        self.play(FadeIn(labels))
        self.wait(2)
