"""
範例 10：整數序列的 Raney 引理——x = (2, −1, −3, 1, 2) 的總和是 1，5 種循環讀法中恰有一種「出發後高度都 > 0」，
起點是前 5 點的最後一個最低點（第 3 點）；加強版：正的點數 3,1,2,5,4 剛好是 1..5 各一次。
先用 combkit.raney 算出來並 assert，再用 combkit.draw.raney_draw 畫。

在 examples/ 資料夾裡執行（這裡有 custom_config.yml）：
    manimgl ex10_raney_lemma.py RaneyLemmaExample
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))     # 讓 combkit 可以 import（已 pip install 就不需要）

from manimlib import *

from combkit.draw.fulton_draw import assert_in_frame
from combkit.draw.raney_draw import IntPathMob, level_line, readings_row
from combkit.draw.scene_base import ChuScene
from combkit.raney import NOTE015_SEQ, good_starts, heights, last_minimum, positive_count, reading

SEQ = NOTE015_SEQ


class RaneyLemmaExample(ChuScene):
    def construct(self):
        # 1. 數學：好讀法只有一種，起點是最後一個最低點；加強版的正點數是 1..M 的排列
        M = len(SEQ)
        assert sum(SEQ) == 1 and good_starts(SEQ) == [last_minimum(SEQ) % M] == [3]
        counts = [positive_count(reading(SEQ, i)) for i in range(M)]
        assert counts == [3, 1, 2, 5, 4] and sorted(counts) == list(range(1, M + 1))

        # 2. 畫面：路徑、5 種讀法（✓/✗），再換成正點數
        self.show_header("Raney 引理：總和 1 的整數序列")
        path = IntPathMob(SEQ, unit=0.55, origin=np.array([-1.4, 1.2, 0]))
        base = level_line(path, 0, 0, M)
        low = Dot(path.vertex(3), radius=0.09, fill_color=MAROON_B)
        assert heights(SEQ)[3] == min(heights(SEQ)[1:])
        self.say("$x=(2,-1,-3,1,2)$", "：最後一個最低點在第 3 點")
        self.play(FadeIn(base), ShowCreation(path), run_time=1.5)
        self.play(FadeIn(low, scale=1.5))
        self.wait(1)
        row = readings_row(SEQ, range(M), center=np.array([0, -1.0, 0]), spacing=2.5, unit=0.22)
        assert_in_frame(path, row)
        self.say("從第 ", "$i$", " 點開始讀：只有 ", "$i=3$", " 出發後一直在 0 之上")
        self.play(LaggedStartMap(FadeIn, row, lag_ratio=0.2), run_time=2)
        self.wait(1.5)
        nums = readings_row(SEQ, range(M), center=np.array([0, -1.0, 0]), spacing=2.5, unit=0.22, mark="count").tags
        self.say("加強版：正的點數 ", "$3,1,2,5,4$", " 是 ", "$1$", " 到 ", "$5$", " 各一次")
        self.play(*[Transform(a, b) for a, b in zip(row.tags, nums)])
        self.wait(2)
