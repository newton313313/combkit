"""
範例 28：Young 表的藝術桌布（youngart 的做法）——YPlate 只要寫 art()；右上角是 PLATE NN · CHAPTER／APPENDIX。
這裡畫兩張：(3,3,3) 的 42 個標準表畫成 3×3 方框裡 Young 格的 42 條路（f^{(3,3,3)} = 42 = Gr(3,6) 的次數），
與 900 個數的隨機排列的 Q 表依「長出來的時間」上色成梯田（young_extras.rsk_rows）。輸出 ex28_*.png（960×540）。

在 examples/ 資料夾裡執行（這裡有 custom_config.yml）：
    python ex28_young_wallpaper.py
"""
import bisect
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))     # 讓 combkit 可以 import（已 pip install 就不需要）
ANIM = "--anim" in sys.argv                                        # 加 --anim 另出 5 秒動畫（1.26.0 起）
sys.argv = ["manimgl", "-w"]                                      # manimlib 以「寫檔、不開視窗」初始化

from combkit.schubert import rectangle_count                      # noqa: E402
from combkit.tableaux import standard_tableaux                    # noqa: E402
from combkit.draw.young_art import *                              # noqa: E402

RES = (960, 540)
W = tuple(int(x) + 1 for x in rng(28).permutation(900))


class FortyTwoPaths(YPlate):
    plate, chapter, pal = 43, 9, "ORCHID"
    title = "Forty-Two Paths"
    sub_tex = r"\deg Gr(3,6)=f^{(3,3,3)}=42"
    caption_corner = DL

    def art(self):
        parts = [(a, b, c) for a in range(4) for b in range(a + 1) for c in range(b + 1)]
        pos = {}
        for r in range(10):
            ps = sorted((p for p in parts if sum(p) == r), key=lambda p: p[0] - p[2] + 0.5 * p[1])
            for a, p in enumerate(ps):
                pos[p] = np.array([-5.6 + r * 1.32, 1.25 * (a - (len(ps) - 1) / 2) + 0.3])
        T = standard_tableaux((3, 3, 3))
        cols = ramp([self.C(2), self.C(0), self.C(1), self.C(3)], len(T))
        g = VGroup()
        for k, t in enumerate(T):
            where = {v: i for i, r in enumerate(t) for v in r}
            lam, seq = [0, 0, 0], [(0, 0, 0)]
            for v in range(1, 10):
                lam[where[v]] += 1
                seq.append(tuple(lam))
            pts = np.array([pos[p] for p in seq])
            pts[:, 1] += 0.028 * (k - 20.5) * np.sin(PI * np.arange(10) / 9)   # 兩端收攏、中間散開
            g.add(path(pts, cols[k], 1.2, 0.5, smooth=True))
        self.add(glow(g, layers=((5, 0.04), (2.2, 0.1), (1, 0.85))))
        self.add(glow_dots([pos[p] for p in parts], self.C(3), radius=0.1, glow_factor=2.2))


class GrowthRings(YPlate):
    plate, chapter, pal = 60, 12, "ORCHID"
    title = "Young Tableaux"
    sub = "the recording tableau Q of a random permutation of 900"

    def art(self):
        _, Q = rsk_rows(W)
        lam = [len(r) for r in Q]
        cs = 6.9 / max(lam[0], len(lam))
        o = (-0.3 - cs * lam[0] / 2, 0.15 + cs * len(lam) / 2)
        bands = ramp([self.C(3), self.C(2), self.C(0), self.C(1), "#5B2A86"], 14)
        cells = [(i, j) for i, r in enumerate(Q) for j in range(len(r))]
        sq = cell_squares(cells, cs, o, [bands[min(13, int(14 * np.sqrt((Q[i][j] - 1) / len(W))))] for i, j in cells],
                          0.95, gap=0.14)
        self.add(sq)
        self.h["keyed"] = (sq, [(Q[i][j] - 1) / len(W) for i, j in cells])   # 動畫版：依長出來的時間一格格出現
        for k in range(1, 14):                                       # 等高線：Q ≤ t 的格子是一個分割
            mu = [m for m in (bisect.bisect_right(r, 1 + len(W) * (k / 14) ** 2) for r in Q) if m]
            self.add(path(outline_pts(mu, cs, o)[1:-1], self.C(3), 0.8, 0.35))


if __name__ == "__main__":
    assert len(standard_tableaux((3, 3, 3))) == 42 == rectangle_count(3, 6)
    P, Q = rsk_rows(W)
    assert sorted(x for r in Q for x in r) == list(range(1, 901)) and [len(r) for r in P] == [len(r) for r in Q]
    assert len(P[0]) == len(patience(W)[1])                          # 第一列長 = 最長遞增子序列（Schensted）
    outs = []
    for cls in (FortyTwoPaths, GrowthRings):
        out = Path(f"ex28_{cls.__name__}.png")
        render_png(cls, RES, out)
        outs.append(out)
        print("wrote", out)
    contact_sheet(outs, Path("ex28_sheet.png"), cols=2)
    print("wrote ex28_sheet.png")
    if ANIM:   # 0–3 秒：格子依 R-S-K 的時間長出（keyed_reveal）；3–5 秒：鏡頭推近
        A = animate_plate(GrowthRings, dict(together=True, prep=lambda s: keyed_reveal(*s.h["keyed"])))
        render_video(A, (640, 360), 15, Path("."), "ex28_GrowthRings")
        print("wrote ex28_GrowthRings.mp4")
