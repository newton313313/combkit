"""
範例 27：把藝術桌布變成約 5 秒的小動畫（walkart 動畫版的做法）——animate_plate(Plate 子類別, cfg) 照原本的 art()
建好物件（最後一格＝桌布），再依種類排出「長出來」的動畫：線條沿路徑畫出、點雲炸開、填色放大淡入；鏡頭緩慢推近、
文字固定在畫面上。cfg 裡 together=True 讓所有線一起畫，A(scene) 加 0–3 秒的動畫（這裡是 6 顆和畫線同步的領頭光點），
B(scene) 是 3–5 秒的動態段。render_video() 在同一個行程裡直接輸出 mp4。

在 examples/ 資料夾裡執行（這裡有 custom_config.yml）：
    python ex27_art_animation.py          → ex27_FanAnim.mp4（640×360、15 fps，約 5 秒）
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))     # 讓 combkit 可以 import（已 pip install 就不需要）
sys.argv = ["manimgl", "-w"]                                       # manimlib 以「寫檔、不開視窗」初始化

from combkit.draw.art_kit import *                                # noqa: E402


class FanPlate(Plate):
    plate, chapter, pal = 1, 1, "NOCTURNE"
    title = "Four Hundred Walkers"
    sub_tex = r"\langle x_N^2\rangle = N"
    caption_corner = DL

    def art(self):
        N, K = 600, 400
        W = walk1d(N, 27, K)
        t = np.linspace(-6.6, 5.0, N + 1)
        sy = 3.2 / (3 * np.sqrt(N))
        cols = ramp([self.C(2), self.C(1), self.C(0)], 41)
        g = VGroup()
        for w in W:
            k = int(np.clip((w[-1] / (2.5 * np.sqrt(N)) + 1) / 2 * 40, 0, 40))
            g.add(path(np.c_[t, w * sy], cols[k], 0.9, 0.22))
        self.add(g)
        self.h["walks"] = [np.c_[t, w * sy] for w in W[:6]]      # 給領頭光點用
        self.add(glow_dots([[-6.6, 0]], self.C(0), radius=0.35))


CFG = dict(together=True,
           A=lambda s: [travelers(s, s.h["walks"], s.C(0), T_A, rf=win(*GROW, T_A), radius=0.12)])

if __name__ == "__main__":
    cls = animate_plate(FanPlate, CFG)
    render_video(cls, (640, 360), 15, Path("."), "ex27_FanAnim")
    print("wrote ex27_FanAnim.mp4")
