"""
範例 26：藝術桌布（walkart 的做法）——Plate 基底類別只要寫 art()，背景漸層、右上角編號與右下角英文標題自動加上；
render_png() 直接從相機讀畫面存 PNG（不用 manimgl -s）。這裡畫兩張：400 條一維漫步的扇形（±√N 包絡）
與第 4 代 Koch 雪花（D = ln4/ln3 ≈ 1.2619，用 walk_rg.koch_points）。輸出 ex26_*.png（960×540；改 RES 可出 4K）。

在 examples/ 資料夾裡執行（這裡有 custom_config.yml）：
    python ex26_art_wallpaper.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))     # 讓 combkit 可以 import（已 pip install 就不需要）
sys.argv = ["manimgl", "-w"]                                       # manimlib 以「寫檔、不開視窗」初始化

from combkit import walk_rg as RG                                 # noqa: E402
from combkit.draw.art_kit import *                                # noqa: E402

RES = (960, 540)


class FanPlate(Plate):
    plate, chapter, pal = 1, 1, "NOCTURNE"
    title = "Four Hundred Walkers"
    sub_tex = r"\langle x_N^2\rangle = N"
    caption_corner = DL

    def art(self):
        N, K = 600, 400
        W = walk1d(N, 26, K)
        t = np.linspace(-6.6, 5.0, N + 1)
        sy = 3.2 / (3 * np.sqrt(N))
        cols = ramp([self.C(2), self.C(1), self.C(0)], 41)
        for w in W:
            k = int(np.clip((w[-1] / (2.5 * np.sqrt(N)) + 1) / 2 * 40, 0, 40))
            self.add(path(np.c_[t, w * sy], cols[k], 0.9, 0.22))
        for s in (1, -1):
            self.add(glow(path(np.c_[t, s * np.sqrt(np.arange(N + 1)) * sy], self.C(0), 1.6)))
        self.add(glow_dots([[-6.6, 0]], self.C(0), radius=0.35))   # 點雲前會自動夾透明小圖（渲染器的毛病）


class KochPlate(Plate):
    plate, chapter, pal = 2, 11, "GLACIER"
    title = "Koch Snowflake"
    sub_tex = r"D=\frac{\ln 4}{\ln 3}\approx 1.2619"

    def art(self):
        a = np.deg2rad([90, -30, 210])
        V = np.c_[np.cos(a), np.sin(a)] * 3.2 + [-1.0, 0.1]
        for g, col in zip(range(5), ramp([self.C(2), self.C(1), self.C(0)], 5)):
            P = np.vstack([RG.koch_points(g, V[i], V[(i + 1) % 3])[:-1] for i in range(3)])
            P = (P - [-1.0, 0.1]) * 0.84 ** (4 - g) + [-1.0, 0.1]
            self.add(path(np.vstack([P, P[:1]]), col, 1.0 + 0.4 * g, 0.4 + 0.15 * g))


if __name__ == "__main__":
    assert f"{np.log(4) / np.log(3):.4f}" == "1.2619"
    outs = []
    for cls in (FanPlate, KochPlate):
        out = Path(f"ex26_{cls.__name__}.png")
        render_png(cls, RES, out)
        outs.append(out)
        print("wrote", out)
    contact_sheet(outs, Path("ex26_sheet.png"), cols=2)
    print("wrote ex26_sheet.png")
