"""
範例 4：R-S-K（Fulton 第 4 章）——先用程式算出排列對應的一對標準表 (P, Q)，確認反過來能還原，再用 manim 畫出來。

在 examples/ 資料夾裡執行（這裡有 custom_config.yml）：
    manimgl ex04_rsk.py RSKExample
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))     # 讓 combkit 可以 import（已 pip install 就不需要）

from manimlib import *

from combkit.draw.fulton_draw import WordMob, assert_in_frame, tableau_mob
from combkit.draw.scene_base import ChuScene
from combkit.tableaux import inverse_rsk, is_standard, rsk_steps

W = (4, 1, 3, 5, 2)


class RSKExample(ChuScene):
    def construct(self):
        # 1. 數學：逐步列插入，並確認 (P, Q) 都是標準表、反向擠出能還原排列
        steps = rsk_steps(W)
        P, Q = steps[-1][2], steps[-1][3]
        assert is_standard(P) and is_standard(Q) and inverse_rsk(P, Q) == (tuple(range(1, 6)), W)

        # 2. 畫面：上面是排列，左下長出 P、右下在新格子記下步數得到 Q
        self.show_header("R-S-K：排列 ↔ 一對標準表")
        word = WordMob(W, size=0.6, center=2.3 * UP)
        names = VGroup(Tex("P", font_size=44).move_to(np.array([-3.4, 1.35, 0.0])),
                       Tex("Q", font_size=44).move_to(np.array([2.1, 1.35, 0.0])))
        self.say("把排列的數依序列插入 ", "$P$", "；同時在 ", "$Q$", " 的新格子寫下這是第幾步")
        self.play(FadeIn(word), FadeIn(names))
        p_mob = q_mob = None
        for k, (_, _, Pk, Qk, box) in enumerate(steps):
            new_p = tableau_mob(Pk, cs=0.6, top_left=(-4.0, 0.9, 0))
            new_q = tableau_mob(Qk, cs=0.6, top_left=(1.5, 0.9, 0))
            new_q.fill_cells([box], YELLOW, 0.35)
            assert_in_frame(new_p, new_q)
            anims = [word.tiles[k].square.animate.set_fill(YELLOW, opacity=0.4)]
            if p_mob is None:
                anims += [FadeIn(new_p), FadeIn(new_q)]
            else:
                anims += [ReplacementTransform(p_mob, new_p), ReplacementTransform(q_mob, new_q)]
            self.play(*anims, run_time=0.9)
            p_mob, q_mob = new_p, new_q
        self.say("反過來：從 ", "$Q$", " 最大的數那一格對 ", "$P$", " 反向擠出，就能一步一步還原排列")
        self.wait(2.5)
