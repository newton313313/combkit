"""
範例 5：Littlewood–Richardson 規則（Fulton 第 5 章）——先用程式找出 (3,2,1)/(2,1) 上內容 (2,1) 的 LR 斜表，
確認個數 = s(2,1)·s(2,1) 裡 s(3,2,1) 的係數 2，再用 manim 畫出來。

在 examples/ 資料夾裡執行（這裡有 custom_config.yml）：
    manimgl ex05_lr_rule.py LRRuleExample
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))     # 讓 combkit 可以 import（已 pip install 就不需要）

from manimlib import *

from combkit.draw.fulton_draw import SkewMob, assert_in_frame, word_tex
from combkit.draw.scene_base import ChuScene
from combkit.tableaux import is_reverse_lattice, lr_product, lr_tableaux, skew_tableaux, tableau_word

NU, LAM, MU = (3, 2, 1), (2, 1), (2, 1)


class LRRuleExample(ChuScene):
    def construct(self):
        # 1. 數學：所有內容 (2,1) 的斜表、其中字是反格字的（LR 斜表），以及 s(2,1)² 的展開
        cands = [S for S in skew_tableaux(NU, LAM, 2) if sorted(v for r in S for v in r if v) == [1, 1, 2]]
        good = lr_tableaux(NU, LAM, MU)
        assert len(cands) == 3 and len(good) == 2 == lr_product(LAM, MU)[NU]
        assert all(is_reverse_lattice(tableau_word(S)) == (S in good) for S in cands)

        # 2. 畫面：三個候選斜表，字是反格字的打勾
        self.show_header("LR 規則：數斜表就是算係數")
        self.say("形狀 ", "$(3,2,1)/(2,1)$", "、內容 ", "$(2,1)$", " 的斜表有 3 個；字倒著讀，1 不能落後 2")
        mobs = VGroup(*[SkewMob(S, cs=0.7) for S in cands]).arrange(RIGHT, buff=1.4).move_to(0.6 * UP)
        words = VGroup(*[word_tex(tableau_word(S), font_size=34).next_to(m, DOWN, buff=0.35)
                         for S, m in zip(cands, mobs)])
        marks = VGroup(*[Tex(r"\checkmark" if S in good else r"\times", font_size=40,
                             fill_color=GREEN if S in good else RED).next_to(w, DOWN, buff=0.25)
                         for S, w in zip(cands, words)])
        result = Tex(r"c^{(3,2,1)}_{(2,1),(2,1)}=2", font_size=44).move_to(2.3 * DOWN)
        assert_in_frame(mobs, words, marks, result)
        self.play(LaggedStartMap(FadeIn, mobs, lag_ratio=0.3), run_time=1.5)
        self.play(FadeIn(words), run_time=0.8)
        self.play(LaggedStartMap(FadeIn, marks, lag_ratio=0.3), run_time=1.0)
        self.say("兩個是 LR 斜表：", "$s_{(2,1)}^2$", " 裡 ", "$s_{(3,2,1)}$", " 的係數正是 2")
        self.play(Write(result))
        self.wait(2)
