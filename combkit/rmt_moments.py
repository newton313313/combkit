"""
rmt_moments.py —— 隨機矩陣的動差與卡塔蘭數（純 Python，不需要 manim）

高斯變數數「所有配對」，大隨機矩陣的特徵值數「不交叉配對」：
    E X^{2k} = (2k-1)!!            （Wick：每一種配對貢獻 1）
    (1/N) E Tr J^{2k} -> C_k       （只有不交叉配對／樹形閉路徑存活）
本模組把有限 N 的動差「精確地」算成 1/N 的多項式：對 Tr J^{2k} 的 2k 個因子做 Wick 配對，
每一對有兩種接法（GOE）或一種接法（GUE），數獨立指標的類別數。

約定（N×N）：
    GOE  J 對稱實矩陣，E[J_ab J_cd] = (δ_ac δ_bd + δ_ad δ_bc)/N   （非對角方差 1/N、對角 2/N）
    GUE  H Hermite 矩陣，E[H_ab H_cd] = δ_ad δ_bc / N               （|H_ab|^2 的期望 1/N）
多項式表示成 {N 的冪次: 係數}，例如 GOE k=2 是 {1: 2, 0: 5, -1: 5}，即 E Tr J^4 = 2N + 5 + 5/N。

來源：spglass 專案（C:\\mypython\\spglass）先修講義的練習解答（notes/prerequisite_solutions.pdf，練習 3(a) 與附錄）。
"""
from __future__ import annotations

import itertools
from fractions import Fraction
from math import prod

from combkit.counting import (binary_to_plane, binary_trees, catalan, dyck_to_matching, matching_is_noncrossing,
                              perfect_matchings, plane_tree_dfs, plane_tree_walk)


# ------------------------------------------------------------------ 1. 高斯：所有配對
def double_factorial_odd(k: int) -> int:
    """(2k-1)!! = 1·3·5···(2k-1) = 2k 個點的配對數 = 標準高斯的 E X^{2k}。"""
    return prod(range(1, 2 * k, 2))


def pairing_counts(k: int) -> tuple[int, int]:
    """(全部配對數, 不交叉配對數) = ((2k-1)!!, C_k)，用暴力列舉。"""
    allm = perfect_matchings(2 * k)
    return len(allm), sum(1 for m in allm if matching_is_noncrossing(m))


# ------------------------------------------------------------------ 2. 有限 N 的精確動差
def _classes(n: int, unions) -> int:
    par = list(range(n))

    def f(x):
        while par[x] != x:
            par[x] = par[par[x]]
            x = par[x]
        return x
    for a, b in unions:
        par[f(a)] = f(b)
    return len({f(x) for x in range(n)})


def _moment_poly(k: int, twists) -> dict[int, int]:
    """E Tr M^{2k}：第 t 個因子是 M_{i_t, i_{t+1}}（指標 mod 2k）。twists 是每一對允許的接法。"""
    n = 2 * k
    out: dict[int, int] = {}
    for pm in perfect_matchings(n):
        for tw in itertools.product(twists, repeat=k):
            unions = []
            for (s, t), w in zip(pm, tw):
                a, b, c, d = s, (s + 1) % n, t, (t + 1) % n
                unions += [(a, c), (b, d)] if w == "straight" else [(a, d), (b, c)]
            p = _classes(n, unions) - k                # N^{類別數} × (1/N)^k
            out[p] = out.get(p, 0) + 1
    return dict(sorted(out.items(), reverse=True))


def goe_moment_poly(k: int) -> dict[int, int]:
    """GOE 的 E Tr J^{2k}，表示成 {N 的冪次: 係數}。k=2：2N + 5 + 5/N；k=3：5N + 22 + 52/N + 41/N^2。"""
    return _moment_poly(k, ("straight", "twisted"))


def gue_moment_poly(k: int) -> dict[int, int]:
    """GUE 的 E Tr H^{2k}（Harer–Zagier 展開）：N^{1-2g} 的係數是黏成虧格 g 曲面的方法數。k=2：2N + 1/N。"""
    return _moment_poly(k, ("twisted",))


def normalized_moment(poly: dict[int, int], N) -> Fraction:
    """(1/N) E Tr M^{2k} 在給定 N 的精確值（Fraction）。"""
    N = Fraction(N)
    return sum(c * N ** (p - 1) for p, c in poly.items())


def _goe_brute(k: int, N: int) -> Fraction:
    """獨立核對：對所有指標組直接把乘積寫成獨立高斯元素的單項式，用一維高斯動差取期望。"""
    total = Fraction(0)
    for idx in itertools.product(range(N), repeat=2 * k):
        cnt: dict[tuple[int, int], int] = {}
        for t in range(2 * k):
            a, b = idx[t], idx[(t + 1) % (2 * k)]
            e = (min(a, b), max(a, b))
            cnt[e] = cnt.get(e, 0) + 1
        if any(m % 2 for m in cnt.values()):
            continue
        term = Fraction(1)
        for (a, b), m in cnt.items():
            var = Fraction(2 if a == b else 1, N)
            term *= double_factorial_odd(m // 2) * var ** (m // 2)
        total += term
    return total / N


# ------------------------------------------------------------------ 3. 存活的路徑：樹形閉路徑
def tree_walks(k: int) -> list[list[int]]:
    """長度 2k 的「樹形閉路徑」（深度優先走一棵 k 條邊的平面樹再回到根），頂點按第一次到訪編號。共 C_k 條。"""
    return [plane_tree_walk(binary_to_plane(t)) for t in binary_trees(k)]


def walk_to_dyck(walk: list[int]) -> str:
    """樹形閉路徑 → Dyck 字：走到新頂點記 U，回到已到訪的頂點記 D。"""
    seen, out = {walk[0]}, []
    for v in walk[1:]:
        out.append("D" if v in seen else "U")
        seen.add(v)
    return "".join(out)


def walk_pairing(walk: list[int]) -> list[tuple[int, int]]:
    """樹形閉路徑裡同一條邊的兩次經過配成一對（第 t 步是邊 walk[t]–walk[t+1]）。"""
    first: dict[frozenset, int] = {}
    pairs = []
    for t in range(len(walk) - 1):
        e = frozenset((walk[t], walk[t + 1]))
        if e in first:
            pairs.append((first.pop(e), t))
        else:
            first[e] = t
    return sorted(pairs)


if __name__ == "__main__":
    # 1. 全部配對 vs 不交叉配對
    for k in range(1, 6):
        a, nc = pairing_counts(k)
        assert a == double_factorial_odd(k) and nc == catalan(k), k
    print("配對數：(2k-1)!! 與 C_k 一致（k=1..5）")

    # 2. 精確多項式
    assert goe_moment_poly(2) == {1: 2, 0: 5, -1: 5}
    assert goe_moment_poly(3) == {1: 5, 0: 22, -1: 52, -2: 41}
    assert gue_moment_poly(2) == {1: 2, -1: 1}
    assert gue_moment_poly(3) == {1: 5, -1: 10}
    for k in range(1, 5):
        g, u = goe_moment_poly(k), gue_moment_poly(k)
        assert g[1] == u[1] == catalan(k)                         # 領頭項 = 卡塔蘭數
        assert sum(g.values()) == double_factorial_odd(k) * 2 ** k   # 全部配對×接法：N=1 時 J ~ N(0,2)
        assert sum(u.values()) == double_factorial_odd(k)            # N=1 時 H ~ N(0,1)
        assert all((1 - p) % 2 == 0 for p in u)                      # GUE 只有 N^{1-2g}
    for k, N in [(2, 2), (2, 3), (3, 2), (3, 3), (2, 4)]:
        assert normalized_moment(goe_moment_poly(k), N) == _goe_brute(k, N), (k, N)
    print("GOE/GUE 精確動差：與逐指標暴力計算一致；領頭項都是 C_k")

    # 3. 樹形閉路徑 ↔ Dyck 路徑 ↔ 不交叉配對
    for k in range(1, 7):
        ws = tree_walks(k)
        assert len(ws) == catalan(k)
        for w, t in zip(ws, binary_trees(k)):
            d = walk_to_dyck(w)
            assert d == plane_tree_dfs(binary_to_plane(t))
            assert sorted(walk_pairing(w)) == sorted(dyck_to_matching(d))
            assert matching_is_noncrossing(walk_pairing(w))
            assert len(set(w)) == k + 1                               # 樹：k 條邊、k+1 個頂點
    print("樹形閉路徑：C_k 條，與 Dyck 路徑、不交叉配對一一對應（k=1..6）")
    print("rmt_moments：隨機矩陣動差與卡塔蘭數全部驗證通過")
