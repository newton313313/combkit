"""
範例 13：evacuation 與鑰匙（Fulton 附錄 A）——先用程式算出書上 p.184 的 T 與它的對偶表 T*（evacuation），
以及 p.210 的表 1122/233/4 的左右鑰匙 K_-、K_+，再用 manim 並排畫出來。

在 examples/ 資料夾裡執行（這裡有 custom_config.yml）：
    manimgl ex13_evacuation_keys.py EvacuationKeysExample
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))     # 讓 combkit 可以 import（已 pip install 就不需要）

from manimlib import *

from combkit.draw.fulton_draw import assert_in_frame, tableau_mob
from combkit.draw.scene_base import ChuScene
from combkit.fulton_appendix import evacuation, keys


class EvacuationKeysExample(ChuScene):
    def construct(self):
        # 1. 數學（書上的例子）
        T = ((1, 2, 2, 3), (2, 3, 5, 5), (4, 4, 6), (5, 6))
        Ts = evacuation(T, 6)
        assert Ts == ((1, 1, 2, 3), (2, 2, 3, 5), (4, 4, 5), (5, 6)) and evacuation(Ts, 6) == T
        U = ((1, 1, 2, 2), (2, 3, 3), (4,))
        Km, Kp = keys(U)
        assert Km == ((1, 1, 1, 1), (2, 2, 2), (4,)) and Kp == ((2, 2, 2, 2), (3, 3, 3), (4,))

        # 2. 畫面
        self.show_header("evacuation 與鑰匙")
        top = VGroup(tableau_mob(T, cs=0.42), Tex(r"\mapsto", font_size=40), tableau_mob(Ts, cs=0.42))
        top.arrange(RIGHT, buff=0.5).move_to(np.array([0, 1.6, 0]))
        bot = VGroup(tableau_mob(Km, cs=0.42), tableau_mob(U, cs=0.42), tableau_mob(Kp, cs=0.42))
        bot.arrange(RIGHT, buff=1.0).move_to(np.array([0, -1.2, 0]))
        labels = VGroup(*(Tex(s, font_size=30).next_to(m, DOWN, buff=0.2) for s, m in zip(("K_-(T)", "T", "K_+(T)"), bot)))
        assert_in_frame(top, bot, labels)
        self.say("evacuation：拿走左上角、滑動、空格填 ", "$x^*=7-x$")
        self.play(FadeIn(top), run_time=1.5)
        self.say("左右鑰匙：把行重新排隊，最左、最右那一行只跟長度有關")
        self.play(FadeIn(bot), FadeIn(labels), run_time=1.5)
        self.wait(2.0)
