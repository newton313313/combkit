"""
raney_examples.py —— 用程式「製造」Raney (1960) 的例子（純 Python，不需要 manim）

核心想法：循環引理（Raney 定理 2.1）本身就是一個「均勻隨機樹產生器」——
把指定的符號隨機排成一列，再從 n 個好的起點中隨機挑一個轉過去，得到的 n 字列在所有同組成的字列中是均勻分布的
（每個字列恰好被 m 個（字串, 起點）的配對讀到，見 Raney 導讀第 6.1 節的雙重計數）。

    python utils/raney_examples.py                       自我檢查（均勻性用卡方式的計數檢查）
    python utils/raney_examples.py --md [--seed 7] [--count 3] [--out 檔名]
                                                         產生一份 Markdown 例題集（附詳解），預設印到螢幕

例題集的每一題都由程式算出答案並 assert，數字不會錯；換 --seed 就是一份新的例題集。
動畫 N06（scenes/sceneN06_random_examples.py）用這裡的產生器示範。
"""
from __future__ import annotations

import random
import sys
from fractions import Fraction
from math import factorial
from pathlib import Path

if __name__ == "__main__":                       # 直接執行時讓上一層（codes/）可以 import
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from combkit.raney import (
    L_formula, a0_of, a_from_chu, card_game_brute, card_game_formula, cyclic_variant, gould_G,
    head_ranks, heights, is_list, is_word, lagrange_coeff, last_minimum, length_of, list_shifts, lists_with,
    multinomial, partition_to_word, positive_count, raney_proof_trace, rank, reading, solve_w, word_to_partition, word_to_tree,
)

# ---------------------------------------------------------------------------
# 1. 隨機字串、隨機字列（均勻）、隨機樹
# ---------------------------------------------------------------------------


def multiset(n: int, a) -> list:
    """組成 (n; a) 的所有符號：a_0 個 0、a_1 個 1、……"""
    out = [0] * a0_of(n, a)
    for i, ai in enumerate(a, start=1):
        out += [i] * ai
    return out


def random_string(n: int, a, rng: random.Random) -> tuple:
    s = multiset(n, a)
    rng.shuffle(s)
    return tuple(s)


def random_list(n: int, a, rng: random.Random, trace: bool = False):
    """均勻隨機的 n 字列：隨機字串 → 隨機挑一個好的起點轉過去。trace=True 時一併回傳 (字串, 起點)。"""
    s = random_string(n, a, rng)
    k = rng.choice(list_shifts(s))
    out = cyclic_variant(s, k)
    return (out, s, k) if trace else out


def random_word(a, rng: random.Random, trace: bool = False):
    return random_list(1, a, rng, trace)


def random_composition(rng: random.Random, max_size: int = 7, max_arity: int = 3) -> tuple:
    """隨機的 (a_1, …, a_K)：Σ i·a_i ≤ max_size，至少一個非零。"""
    while True:
        a = tuple(rng.randint(0, max_size // i) for i in range(1, max_arity + 1))
        size = sum(i * x for i, x in enumerate(a, 1))
        if 0 < size <= max_size:
            return a


def random_sum_one(M: int, rng: random.Random, lo: int = -4, hi: int = 3) -> tuple:
    """總和為 1 的隨機整數序列（Raney 引理用）。"""
    while True:
        a = [rng.randint(lo, hi) for _ in range(M - 1)]
        last = 1 - sum(a)
        if lo - 2 <= last <= hi + 2:
            return tuple(a + [last])


def random_noncrossing(n_vec, v_vec, rng: random.Random) -> list:
    """均勻隨機的不交叉分法（Chu 的 (n⃗, v⃗)）：均勻隨機的字 → 「固定一個點」的一一對應。"""
    N = sum(n * v for n, v in zip(n_vec, v_vec))
    return word_to_partition(random_word(a_from_chu(n_vec, v_vec), rng), N)


# ---------------------------------------------------------------------------
# 2. 字 → 運算式
# ---------------------------------------------------------------------------

OPS = {1: "fgh", 2: "FGHPQ", 3: "TUVW", 4: "ABCD"}
ARGS = "abcdefghijklmnopqrstuvwxyz"


def word_to_expression(word) -> tuple:
    """把字寫成運算式：(前序寫法, 一般寫法)。一元用 f,g,h、二元 F,G,H,…、三元 T,U,V,…；引數 a,b,c,…"""
    counters = {}
    args = iter(ARGS)

    def name(k):
        counters[k] = counters.get(k, 0) + 1
        pool = OPS.get(k, "XYZ")
        return pool[(counters[k] - 1) % len(pool)]

    def rec(tree):
        k, kids = tree
        if k == 0:
            x = next(args)
            return x, x
        nm = name(k)
        parts = [rec(c) for c in kids]
        return nm + "".join(p for p, _ in parts), nm + "(" + ",".join(q for _, q in parts) + ")"

    return rec(word_to_tree(word))


def plain(s) -> str:
    return "".join(str(x) for x in s)


# ---------------------------------------------------------------------------
# 3. 例題集（Markdown）
# ---------------------------------------------------------------------------


def _ranks_line(word) -> str:
    return ", ".join(str(h) for h in head_ranks(word))


def ex_expression(rng) -> str:
    a = random_composition(rng, 7, 3)
    w = random_word(a, rng)
    pre, inf = word_to_expression(w)
    assert is_word(w) and len(w) == length_of(1, a)
    return (f"**題。** 運算式 `{inf}` 的前序寫法與型是什麼？是字嗎？\n\n"
            f"**解。** 前序寫法 `{pre}`，型 `{plain(w)}`。開頭的秩依序是 {_ranks_line(w)}："
            f"最後一項才第一次到 −1，所以是字（Rosenbloom 判別法）。組成 (1; {', '.join(map(str, a))})，"
            f"引數 a₀ = {a0_of(1, a)} 個，長度 m = {len(w)}。\n")


def ex_cycle(rng) -> str:
    while True:                                   # 挑「本來不是字列、至少 6 個符號」的字串，例子才有意思
        n = rng.choice([1, 2, 3])
        a = random_composition(rng, 6, 3)
        s = random_string(n, a, rng)
        if len(s) >= 6 and 0 not in list_shifts(s):
            break
    tr = raney_proof_trace(s)
    ks = list_shifts(s)
    assert len(ks) == n and sorted(tr["ks"]) == ks
    words = " | ".join(plain(w) for w in tr["words"])
    return (f"**題。** 字串 σ = `{plain(s)}` 的秩是多少？哪些 k 使 σg_k 是字列？照 Raney 的證明走一遍。\n\n"
            f"**解。** r(σ) = {rank(s)}，所以 n = {n}。好的 k = {ks}。"
            f"秩最大（同秩取最長）的一段 α = `{plain(tr['alpha']) or '（空）'}`（秩 {rank(tr['alpha'])}），"
            f"從它開始讀 σ₁ = σg_{tr['k1']} = `{plain(tr['sigma1'])}` = `{words}`；"
            f"輪換這 {n} 個字得到全部 {n} 個字列。\n")


def ex_count(rng) -> str:
    n = rng.choice([1, 1, 2, 3])
    a = random_composition(rng, 6, 3)
    a0, m = a0_of(n, a), length_of(n, a)
    M = multinomial(a0, *a)
    L = L_formula(n, a)
    assert L == len(lists_with(n, a))
    ex = ", ".join(plain(w) for w in lists_with(n, a)[:4])
    return (f"**題。** 組成 (n; a) = ({n}; {', '.join(map(str, a))}) 的 {n} 字列有幾個？\n\n"
            f"**解。** a₀ = {a0}、m = {m}、M = {M}。由定理 2.2，L = n·M/m = {n}·{M}/{m} = **{L}**"
            f"（程式窮舉也是 {L} 個，例如 {ex}…）。\n")


def ex_chu(rng) -> str:
    while True:                                   # 至少 6 個點
        v_vec = rng.choice([(2, 3), (2,), (3,), (1, 2), (2, 4)])
        n_vec = tuple(rng.randint(1, 3) for _ in v_vec)
        N = sum(n * v for n, v in zip(n_vec, v_vec))
        if 6 <= N <= 10:
            break
    blocks = random_noncrossing(n_vec, v_vec, rng)
    w = partition_to_word(blocks, N)
    shown = ", ".join("{" + ",".join(str(p + 1) for p in b) + "}" for b in blocks)
    L = L_formula(1, a_from_chu(n_vec, v_vec))
    return (f"**題。** 圓周 {N} 點（從正上方順時針編 1..{N}）的不交叉分法 {shown}"
            f"（n⃗ = {n_vec}、v⃗ = {v_vec}）對應哪一個 Raney 的字？\n\n"
            f"**解。** 點 1 所在的群組當根，依序寫出它切出的各段：`{plain(w)}`。"
            f"這種分法共有 C_v⃗(n⃗) = L(1; a) = {L} 種，這一個是用均勻隨機產生器挑出來的。\n")


def ex_lagrange(rng) -> str:
    coeffs = [rng.randint(1, 3)] + [rng.randint(0, 2) for _ in range(3)]
    N = 6
    phi = [Fraction(c) for c in coeffs] + [Fraction(0)] * (N + 1 - len(coeffs))
    w = solve_w(phi, N)
    lag = [lagrange_coeff(phi, n) for n in range(1, N + 1)]
    assert w[1:] == lag
    poly = " + ".join((str(c) if i == 0 else (("" if c == 1 else str(c)) + ("x" if i == 1 else f"x^{i}")))
                      for i, c in enumerate(coeffs) if c)
    return (f"**題。** 解 w = x·φ(w)，φ(x) = {poly}。前 {N} 項係數？\n\n"
            f"**解。** 由 Lagrange 反演 [xⁿ]w = (1/n)[xⁿ⁻¹]φⁿ（= 所有 n 個節點的平面樹的權重和）："
            f"{', '.join(str(x) for x in lag)}（逐次代入解方程得到相同的數）。\n")


def ex_convolution(rng) -> str:
    beta = rng.randint(1, 3)
    a1, a2 = Fraction(rng.randint(1, 5), rng.randint(1, 3)), Fraction(rng.randint(1, 5), rng.randint(1, 3))
    n = rng.randint(2, 4)
    lhs = gould_G(a1 + a2, n, beta, 1)
    terms = [gould_G(a1, k, beta, 1) * gould_G(a2, n - k, beta, 1) for k in range(n + 1)]
    assert lhs == sum(terms)
    return (f"**題。** 驗證 Hagen–Rothe（Raney (7.8)，γ = 1）：α₁ = {a1}、α₂ = {a2}、β = {beta}、n = {n}。\n\n"
            f"**解。** G(α₁+α₂, n) = {lhs}；右邊 Σ G(α₁,k)G(α₂,n−k) = {' + '.join(str(t) for t in terms)} = {sum(terms)}。"
            f"取 α 為正整數 t、β = v 時就是 Chu 的 A_t((n),(v))。\n")


def ex_sum_one(rng) -> str:
    M = rng.randint(4, 7)
    x = random_sum_one(M, rng)
    xm = last_minimum(x)
    good = [i for i in range(M) if positive_count(reading(x, i)) == M]
    assert good == [xm % M]
    counts = [positive_count(reading(x, i)) for i in range(M)]
    return (f"**題。** x = {x}（總和 1）。哪一種循環讀法出發後高度都 > 0？各讀法有幾個正的點？\n\n"
            f"**解。** 高度 {heights(x)[1:]}，最後一個最低點在第 {xm} 點，所以從第 {xm % M} 點開始讀。"
            f"從第 0..{M - 1} 點讀的正點數是 {counts}（1..{M} 各一次）。\n")


def ex_cards(rng) -> str:
    n = rng.randint(1, 3)
    m = rng.randint(n + 1, 5)
    while True:
        w = tuple(sorted(rng.randint(1, m) for _ in range(n)))
        if sum(w) == m:
            break
    got = len(card_game_brute(w, m))
    assert got == card_game_formula(w, m)
    return (f"**題。** 摸牌遊戲（洛谷 P6672）：w = {w}、m = {m}，牌都有標號。幾種排法能摸到勝利牌？\n\n"
            f"**解。** m!/(m−n+1) = {factorial(m)}/{m - n + 1} = **{got}**（窮舉也是 {got}）。"
            f"合格的牌序加上勝利牌就是一個 Raney 的字。\n")


SECTIONS = [
    ("運算式、前序記法與字", ex_expression),
    ("循環引理（定理 2.1）與 Raney 的證明", ex_cycle),
    ("數字列（定理 2.2）", ex_count),
    ("不交叉分法 ↔ 字", ex_chu),
    ("Lagrange 反演", ex_lagrange),
    ("卷積（Gould、Hagen–Rothe）", ex_convolution),
    ("整數序列的 Raney 引理", ex_sum_one),
    ("摸牌遊戲", ex_cards),
]


def examples_markdown(seed: int = 7, count: int = 3) -> str:
    rng = random.Random(seed)
    out = [f"# Raney (1960) 例題集（seed = {seed}）\n",
           "> 由 `codes/utils/raney_examples.py` 產生，每一題的答案都由程式算出並核對；換 `--seed` 就是新的一份。"
           "記號同《Raney (1960) 論文導讀》。\n"]
    for title, fn in SECTIONS:
        out.append(f"\n## {title}\n")
        for k in range(count):
            out.append(f"\n### {title} {k + 1}\n\n" + fn(rng))
    return "\n".join(out)


# ---------------------------------------------------------------------------
# 自我檢查
# ---------------------------------------------------------------------------


def _self_test():
    rng = random.Random(2026)
    # 均勻性：每個字列被抽到的次數差不多（期望值 T/L，容許 ±35%）
    for n, a in [(1, (0, 3)), (1, (1, 1, 1)), (2, (0, 1, 1)), (1, (0, 1, 1))]:
        allw = lists_with(n, a)
        T = 400 * len(allw)
        cnt = {w: 0 for w in allw}
        for _ in range(T):
            w = random_list(n, a, rng)
            assert is_list(w, n)
            cnt[w] += 1
        exp = T / len(allw)
        assert all(0.65 * exp < c < 1.35 * exp for c in cnt.values()), (n, a, cnt)
    # 隨機不交叉分法真的不交叉，而且組成對
    import combkit.counting as ct
    for _ in range(50):
        p = random_noncrossing((2, 1), (2, 3), rng)
        assert ct.is_noncrossing(p) and sorted(len(b) for b in p) == [2, 2, 3]
    # 運算式
    assert word_to_expression((2, 3, 0, 2, 0, 0, 0, 1, 0)) == ("FTaGbcdfe", "F(T(a,G(b,c),d),f(e))")
    assert word_to_expression((2, 2, 0, 0, 0)) == ("FGabc", "F(G(a,b),c)")
    for _ in range(30):
        x = random_sum_one(rng.randint(2, 7), rng)
        assert sum(x) == 1
    md = examples_markdown(7, 2)
    assert md.count("**解。**") == 2 * len(SECTIONS)
    print("Raney 例題產生器：均勻性、隨機不交叉分法、運算式、例題集全部通過")


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--md", action="store_true", help="產生 Markdown 例題集")
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--count", type=int, default=3)
    ap.add_argument("--out")
    args = ap.parse_args()
    if args.md:
        text = examples_markdown(args.seed, args.count)
        if args.out:
            Path(args.out).write_text(text, encoding="utf-8")
            print(f"寫出 {args.out}")
        else:
            sys.stdout.reconfigure(encoding="utf-8")
            print(text)
    else:
        _self_test()
