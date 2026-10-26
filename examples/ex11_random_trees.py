"""
範例 11：用循環引理產生均勻隨機的樹——組成 (1; 0, 2, 1)（兩個二元、一個三元）的字有 21 個；
把 8 個符號隨機打亂，再轉到唯一的好讀法，就得到其中之一（機率都是 1/21）。
先用 combkit.raney_examples 抽樣並 assert，再用 combkit.draw.raney_draw 畫成樹與運算式。

在 examples/ 資料夾裡執行（這裡有 custom_config.yml）：
    manimgl ex11_random_trees.py RandomTreesExample
不畫圖、只要例題集：python -m combkit.raney_examples --md --seed 7
"""
import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))     # 讓 combkit 可以 import（已 pip install 就不需要）

from manimlib import *

from combkit.draw.fulton_draw import assert_in_frame
from combkit.draw.raney_draw import WordTreeMob, word_tiles
from combkit.draw.scene_base import ChuScene
from combkit.raney import L_formula, is_word
from combkit.raney_examples import random_word, word_to_expression

A = (0, 2, 1)


class RandomTreesExample(ChuScene):
    def construct(self):
        # 1. 數學：抽 4 次（固定 seed），每次都得到字；共有 21 個可能
        rng = random.Random(12)                    # 這個 seed 抽出 4 棵不同的樹
        draws = [random_word(A, rng, trace=True) for _ in range(4)]
        assert L_formula(1, A) == 21 and all(is_word(w) for w, _, _ in draws) and len({w for w, _, _ in draws}) == 4

        # 2. 畫面：打亂的字串 → 字 → 樹 → 運算式
        self.show_header("均勻隨機的樹：打亂 → 轉到好的起點")
        cols = VGroup()
        for w, s, k in draws:
            col = VGroup(word_tiles(s, size=0.26, font_size=16, buff=0.02), Tex(r"\downarrow", font_size=26),
                         word_tiles(w, size=0.26, font_size=16, buff=0.02),
                         WordTreeMob(w, x_gap=0.36, y_gap=0.45),
                         Tex(word_to_expression(w)[1], font_size=24)).arrange(DOWN, buff=0.22)
            cols.add(col)
        cols.arrange(RIGHT, buff=0.55, aligned_edge=UP).move_to(0.3 * UP)
        assert_in_frame(cols)
        self.say("組成 ", "$(1;0,2,1)$", "：21 個字裡均勻地挑一個")
        for c in cols:
            self.play(FadeIn(c[:3]), run_time=0.8)
            self.play(FadeIn(c[3:]), run_time=0.8)
        self.wait(2)
