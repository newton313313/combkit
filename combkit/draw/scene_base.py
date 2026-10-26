"""
scene_base.py —— 所有場景共用的基底類別

提供兩個小工具：
    self.show_header("場景標題")   在左上角放一個小標題
    self.say("一句說明")           在畫面底部顯示（或更換）字幕
"""
from manimlib import *

from combkit.draw.styles import MUTED_COLOR, mixed, zh


class ChuScene(Scene):
    header_font_size = 26
    caption_font_size = 32
    caption_y = -3.45

    def setup(self):
        super().setup()
        self.caption_mob = None
        self.header_mob = None

    # ------------------------------------------------------------------ 標題
    def show_header(self, text: str, animate: bool = True):
        header = zh(text, font_size=self.header_font_size, color=MUTED_COLOR)
        header.to_corner(UL, buff=0.35)
        header.fix_in_frame()
        if animate:
            self.play(FadeIn(header, shift=0.2 * RIGHT), run_time=0.6)
        else:
            self.add(header)
        self.header_mob = header
        return header

    # ------------------------------------------------------------------ 字幕
    def make_caption(self, *parts, color=WHITE, t2c=None):
        if len(parts) == 1 and isinstance(parts[0], Mobject):
            cap = parts[0]
        else:
            cap = mixed(*parts, font_size=self.caption_font_size, color=color, t2c=t2c)
        if cap.get_width() > FRAME_WIDTH - 1.0:
            cap.set_width(FRAME_WIDTH - 1.0)
        cap.move_to(self.caption_y * UP)
        cap.fix_in_frame()   # 鏡頭移動時字幕留在原位
        return cap

    def say(self, *parts, wait: float = 0.0, run_time: float = 0.6, color=WHITE, t2c=None):
        """
        在底部顯示字幕；已有字幕時會換掉舊的。
        parts 的寫法與 styles.mixed 相同：以 '$' 開頭的是數學式，其餘是中文。
        """
        cap = self.make_caption(*parts, color=color, t2c=t2c)
        if self.caption_mob is None:
            self.play(FadeIn(cap, shift=0.2 * UP), run_time=run_time)
        else:
            self.play(
                FadeOut(self.caption_mob, shift=0.2 * UP),
                FadeIn(cap, shift=0.2 * UP),
                run_time=run_time,
            )
        self.caption_mob = cap
        if wait > 0:
            self.wait(wait)
        return cap

    def clear_caption(self, run_time: float = 0.4):
        if self.caption_mob is not None:
            self.play(FadeOut(self.caption_mob), run_time=run_time)
            self.caption_mob = None

    # ------------------------------------------------------------------ 收尾
    def clear_stage(self, run_time: float = 0.7):
        """淡出標題以外的所有東西（字幕也會淡出）。"""
        mobs = [m for m in self.mobjects if m is not self.frame and m is not self.header_mob]
        if mobs:
            self.play(*[FadeOut(m) for m in mobs], run_time=run_time)
        self.caption_mob = None

    def fade_all(self, run_time: float = 0.8):
        """把畫面上所有東西淡出。"""
        mobs = [m for m in self.mobjects if m is not self.frame]
        if mobs:
            self.play(*[FadeOut(m) for m in mobs], run_time=run_time)
        self.caption_mob = None
        self.header_mob = None
