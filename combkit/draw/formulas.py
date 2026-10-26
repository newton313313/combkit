"""
formulas.py —— 常用公式的 Tex 物件

顏色規則：x_1 與 n_1, v_1 用第 1 型的顏色（藍），x_2 與 n_2, v_2 用第 2 型的顏色（黃），
讓觀眾一眼看出變數對應哪一種多邊形。
"""
from math import factorial

from manimlib import *

from combkit.counting import closed_form, inner
from combkit.draw.styles import (
    D_COLOR,
    GENERIC_COLOR,
    TYPE1_COLOR,
    TYPE2_COLOR,
    block_color,
)

# 生成函數相關公式共用的著色表
GF_T2C = {
    "x_1": TYPE1_COLOR,
    "x_2": TYPE2_COLOR,
    "x_k": GENERIC_COLOR,
    "v_1": TYPE1_COLOR,
    "v_2": TYPE2_COLOR,
    "v_k": GENERIC_COLOR,
}


def catalan_formula(font_size: float = 48) -> Tex:
    """C_n = 1/(n+1) * binom(2n, n)"""
    return Tex(R"C_n=\frac{1}{n+1}\binom{2n}{n}", font_size=font_size)


def vec_str(values) -> str:
    return "(" + ",".join(str(x) for x in values) + ")"


def create_vector_formula(n_values, v_values, font_size: float = 44) -> Tex:
    """例如 \\vec n=(1,1), \\vec v=(2,3)。"""
    return Tex(
        R"\vec n=" + vec_str(n_values) + R",\qquad \vec v=" + vec_str(v_values),
        font_size=font_size,
    )


def create_total_points_formula(n_values, v_values, font_size: float = 44) -> Tex:
    """N = n_1 x v_1 + ... = 總點數（帶入數字）。"""
    terms = "+".join(fR"{n}\times{v}" for n, v in zip(n_values, v_values))
    return Tex(fR"N={terms}={inner(n_values, v_values)}", font_size=font_size)


def create_multinomial_formula(n_values, v_values, font_size: float = 44, steps: bool = True) -> Tex:
    """
    把數字代入論文式 (2)：
        C_v(n) = 1/(N+1) * binom(N+1; n_1, ..., n_k, L) = ... = 數值
    """
    N = inner(n_values, v_values)
    L = 1 + sum(n * (v - 1) for n, v in zip(n_values, v_values))
    lower = ",".join(str(n) for n in n_values) + f",{L}"
    facts = R"\,".join(f"{n}!" for n in n_values) + fR"\,{L}!"
    value = closed_form(n_values, v_values)
    tex = (
        fR"C_{{{vec_str(v_values)}}}{vec_str(n_values)}"
        fR"=\frac{{1}}{{{N + 1}}}\binom{{{N + 1}}}{{{lower}}}"
    )
    if steps:
        tex += fR"=\frac{{1}}{{{N + 1}}}\cdot\frac{{{N + 1}!}}{{{facts}}}"
    tex += f"={value}"
    assert factorial(N + 1) // (N + 1) // factorial(L) // _prod_fact(n_values) == value
    return Tex(tex, font_size=font_size)


def _prod_fact(values) -> int:
    out = 1
    for x in values:
        out *= factorial(x)
    return out


def general_formula(font_size: float = 44) -> Tex:
    """論文式 (2) 的一般形式。"""
    return Tex(
        R"C_{\vec v}(\vec n)=\frac{1}{N+1}\binom{N+1}{n_1,\,n_2,\,\ldots,\,n_k,\;L}",
        font_size=font_size,
    )


def gf_definition(font_size: float = 44) -> Tex:
    """論文式 (3)：生成函數 D 的定義。"""
    return Tex(
        R"D(\vec X)=\sum_{\vec n} d(\vec n)\,x_1^{n_1}x_2^{n_2}\cdots x_k^{n_k}",
        font_size=font_size,
        # 整個因子 x_i^{n_i} 一起上色（只給 x_i 上色會讓指數的位置跑掉）
        t2c={"D": D_COLOR, R"x_1^{n_1}": TYPE1_COLOR, R"x_2^{n_2}": TYPE2_COLOR,
             R"x_k^{n_k}": GENERIC_COLOR},
    )


def gf_equation(font_size: float = 48) -> Tex:
    """
    論文式 (6)：D = 1 + x_1 D^{v_1} + ... + x_k D^{v_k}。
    """
    return Tex(
        R"D=1+x_1D^{v_1}+x_2D^{v_2}+\cdots+x_kD^{v_k}",
        font_size=font_size, t2c={"D": D_COLOR, **GF_T2C},
    )


def recurrence_formula(font_size: float = 40) -> Tex:
    """論文式 (5)：遞迴關係。"""
    return Tex(
        R"d(\vec n)=\sum_{i=1}^{k}\ \ \sum\limits_{\vec r_1+\cdots+\vec r_{v_i}=\vec n-\vec e_i}\ "
        R"d(\vec r_1)\,d(\vec r_2)\cdots d(\vec r_{v_i})",
        font_size=font_size,
    )


def size_colored(text: str, size: int, font_size: float = 40) -> Tex:
    """用「size 點群組」的顏色畫一段 Tex。"""
    return Tex(text, font_size=font_size, fill_color=block_color(size))
