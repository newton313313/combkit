"""
styles.py —— 全專案共用的視覺規則

核心規則（所有場景都遵守）：
  * 圓周上的點：白色小圓點
  * 同一個群組（block）的多邊形／弦：同一種顏色
  * 群組顏色由「點數」決定：2 點=藍、3 點=黃、4 點=綠、5 點=紫
    因此在 v=(2,3) 的例子裡：第 1 型（弦）= 藍 = x_1，第 2 型（三角形）= 黃 = x_2
  * 交叉：永遠紅色
  * 固定點（root point）：永遠橘色、比較大
"""
from manimlib import *

# ---------------------------------------------------------------------------
# 字型
# ---------------------------------------------------------------------------
CJK_FONT = "Microsoft JhengHei"      # 微軟正黑體（Windows 內建）

# ---------------------------------------------------------------------------
# 顏色
# ---------------------------------------------------------------------------
POINT_COLOR = WHITE
CIRCLE_COLOR = GREY_C
TEXT_COLOR = WHITE
MUTED_COLOR = GREY_B

TYPE1_COLOR = BLUE        # 2 點群組（弦）          <-> x_1
TYPE2_COLOR = YELLOW      # 3 點群組（三角形）      <-> x_2
TYPE3_COLOR = GREEN       # 4 點群組（四邊形）
TYPE4_COLOR = PURPLE_B    # 5 點群組（五邊形）
GENERIC_COLOR = GREEN     # 公式裡「第 k 型」的代表色

CROSSING_COLOR = RED
HIGHLIGHT_COLOR = ORANGE  # 固定點
ROOT_COLOR = HIGHLIGHT_COLOR
OK_COLOR = "#7CE38B"
D_COLOR = TEAL_B          # 生成函數 D

REGION_COLORS = [TEAL, PINK, PURPLE_B]   # 根群組切出的區域 R_1, R_2, R_3

SIZE_TO_COLOR = {
    1: GREY_A,
    2: TYPE1_COLOR,
    3: TYPE2_COLOR,
    4: TYPE3_COLOR,
    5: TYPE4_COLOR,
}


def block_color(size: int):
    """依群組的點數回傳顏色。"""
    return SIZE_TO_COLOR.get(size, GREY_A)


# ---------------------------------------------------------------------------
# 尺寸
# ---------------------------------------------------------------------------
DOT_RADIUS = 0.075
BLOCK_STROKE = 5.0
BLOCK_FILL_OPACITY = 0.30
CIRCLE_STROKE = 2.5


# ---------------------------------------------------------------------------
# 文字
# ---------------------------------------------------------------------------
def zh(text: str, font_size: float = 32, color=TEXT_COLOR, **kwargs) -> Text:
    """中文（或中英混合）文字。"""
    return Text(text, font=CJK_FONT, font_size=font_size, fill_color=color, **kwargs)


def mixed(*parts, font_size: float = 32, color=TEXT_COLOR, buff: float = 0.12, t2c=None) -> VGroup:
    """
    把中文與數學式排成一行。
    以 '$' 開頭的字串當成 LaTeX 數學式，其餘當成中文。
        mixed("要 ", "$n_i$", " 個 ", "$v_i$", " 點群組")
    """
    t2c = t2c or {}
    mobs = []
    for part in parts:
        if isinstance(part, Mobject):
            mobs.append(part)
        elif part.startswith("$"):
            mobs.append(Tex(part.strip("$"), font_size=font_size * 1.3, fill_color=color, t2c=t2c))
        else:
            mobs.append(zh(part, font_size=font_size, color=color))
    return VGroup(*mobs).arrange(RIGHT, buff=buff)


def boxed(mob: Mobject, color=YELLOW, buff: float = 0.25, stroke_width: float = 3) -> VGroup:
    """在物件外面加一個框，回傳 VGroup(框, 物件)。"""
    rect = SurroundingRectangle(mob, buff=buff, color=color, stroke_width=stroke_width)
    return VGroup(rect, mob)
