"""
raney.py —— G. N. Raney, *Functional composition patterns and power series reversion*,
Trans. AMS 94 (1960) 441–451 的組合與冪級數計算。純 Python，不需要 manim。

    python utils/raney.py        （在 codes/ 裡執行）自我檢查：每個主張都用暴力列舉再算一次

講義 notes/raney1960_導讀.tex 的數字與圖都出自這裡（經 notes/raney/make_figs.py）。

記號照 Raney：
    字串 σ   ：非負整數的有限序列，寫成 tuple，例如 (2, 0, 0)；第 x 項是「x 元運算子」，0 是「引數」。
    秩 r(σ)  ：Σ (σ(x) − 1)，也就是 Σ_i (i − 1) p_i(σ)。
    字       ：前序（波蘭）記法寫出來的一棵平面樹；r = −1，而且每個真開頭的秩 ≥ 0（Rosenbloom）。
    n 字列   ：n 個字接起來；r = −n，而且每個真開頭的秩 > −n。
    組成     ：(n; a_1, a_2, …)，n = 字數，a_i = 符號 i 出現的次數；a 寫成 tuple (a_1, …, a_K)。
    L(n; a)  ：組成為 (n; a) 的 n 字列個數（Raney 定理 2.2）。

和 Chu 講義（notes/chu1987_導讀.tex 第 6–7 節、counting.py 第 4 節）的對應：
    第 i 型群組（v_i 點）↔ 符號 v_i；n_i ↔ a_{v_i}；t ↔ n；N + t ↔ 字串長度 m；
    Chu 講義的「路徑」= 把 Raney 字串倒過來讀，符號 j 換成高度變化 1 − j（見 to_chu_path）。
"""
from __future__ import annotations

from fractions import Fraction
from itertools import product
from math import comb, factorial

# ---------------------------------------------------------------------------
# 1. 字串、秩、字（Raney 第 2 節）
# ---------------------------------------------------------------------------


def rank(s) -> int:
    """r(σ) = Σ (σ(x) − 1)。"""
    return sum(x - 1 for x in s)


def head_ranks(s) -> list:
    """[r(空), r(σ 的前 1 項), …, r(σ)]，也就是 Raney 的「秩路徑」。"""
    out = [0]
    for x in s:
        out.append(out[-1] + x - 1)
    return out


def is_list(s, n: int) -> bool:
    """s 是否為 n 字列（Rosenbloom 判別法的推廣）：r = −n，真開頭的秩都 > −n。"""
    hs = head_ranks(s)
    return hs[-1] == -n and all(h > -n for h in hs[:-1])


def is_word(s) -> bool:
    return is_list(s, 1)


def is_list_recursive(s, n: int) -> bool:
    """直接照定義：空字串是 0 字列；字 = 0，或 k α_1 … α_k（α_j 是字）；n 字列 = 字接 (n−1) 字列。"""
    s = tuple(s)
    if n == 0:
        return len(s) == 0
    for cut in range(1, len(s) + 1):
        if _is_word_recursive(s[:cut]) and is_list_recursive(s[cut:], n - 1):
            return True
    return False


def _is_word_recursive(s) -> bool:
    if not s:
        return False
    if s == (0,):
        return True
    k = s[0]
    return k > 0 and is_list_recursive(s[1:], k)


def split_words(s) -> list:
    """把字列切成字：每次取最短的、秩為 −1 的開頭（Raney 定理 2.1 證明裡的做法）。"""
    s = tuple(s)
    out, start, h = [], 0, 0
    for j, x in enumerate(s):
        h += x - 1
        if h == -1:
            out.append(s[start:j + 1])
            start, h = j + 1, 0
    assert start == len(s), "不是字列"
    return out


def word_to_tree(s):
    """字 → 樹：(arity, (子樹, …))；葉子是 (0, ())。"""
    tree, rest = _parse(tuple(s))
    assert not rest, "不是字"
    return tree


def _parse(s):
    k, rest = s[0], s[1:]
    kids = []
    for _ in range(k):
        kid, rest = _parse(rest)
        kids.append(kid)
    return (k, tuple(kids)), rest


def tree_to_word(tree) -> tuple:
    k, kids = tree
    out = (k,)
    for kid in kids:
        out += tree_to_word(kid)
    return out


def composition(s, K: int | None = None):
    """(n; a_1, …, a_K)：n = −r(σ)，a_i = 符號 i 的個數。"""
    K = max([x for x in s] + [1]) if K is None else K
    return -rank(s), tuple(sum(1 for x in s if x == i) for i in range(1, K + 1))


def a0_of(n: int, a) -> int:
    """引數（符號 0）的個數 a_0 = n + Σ (i − 1) a_i。"""
    return n + sum((i - 1) * ai for i, ai in enumerate(a, start=1))


def length_of(n: int, a) -> int:
    """字串長度 m = n + Σ i a_i。"""
    return n + sum(i * ai for i, ai in enumerate(a, start=1))


def multinomial(*parts) -> int:
    """M(a_0, a_1, …) = (Σ a_i)! / Π a_i!。"""
    if any(p < 0 for p in parts):
        return 0
    out = factorial(sum(parts))
    for p in parts:
        out //= factorial(p)
    return out


def strings_with(n: int, a) -> list:
    """所有組成為 (n; a) 的字串（字典序）。"""
    a0 = a0_of(n, a)
    if a0 < 0:
        return []
    letters = {0: a0}
    letters.update({i: ai for i, ai in enumerate(a, start=1) if ai})
    out = []

    def rec(prefix):
        if not any(letters.values()):
            out.append(tuple(prefix))
            return
        for x in sorted(letters):
            if letters[x]:
                letters[x] -= 1
                prefix.append(x)
                rec(prefix)
                prefix.pop()
                letters[x] += 1

    rec([])
    return out


def lists_with(n: int, a) -> list:
    return [s for s in strings_with(n, a) if is_list(s, n)]


def L_formula(n: int, a) -> int:
    """Raney 定理 2.2：m·L(n; a) = n·M(a_0, a_1, …)。"""
    a0 = a0_of(n, a)
    if a0 < 0:
        return 0
    m = length_of(n, a)
    if m == 0:
        return 1
    if n == 0:
        return 0
    value, rem = divmod(n * multinomial(a0, *a), m)
    assert rem == 0
    return value


def L_brute(n: int, a) -> int:
    return len(lists_with(n, a))


# ---------------------------------------------------------------------------
# 2. 循環變形與循環引理（Raney 定理 2.1）
# ---------------------------------------------------------------------------


def cyclic_variant(s, k: int) -> tuple:
    """σ g_k：(σ g_k)(x) = σ(x + k mod m)，也就是從第 k 項開始讀一圈。"""
    s = tuple(s)
    return s[k:] + s[:k]


def list_shifts(s) -> list:
    """所有使 σ g_k 成為 n 字列的 k（n = −r(σ)）。"""
    n = -rank(s)
    return [k for k in range(len(s)) if is_list(cyclic_variant(s, k), n)]


def raney_alpha(s):
    """Raney 證明定理 2.1 的第一步：在所有循環變形的所有開頭裡，取秩最大、同秩中最長的 α。
    回傳 (k, α)：α 是 σ g_k 的開頭；若有多個，取 k 最小的。"""
    best = None
    for k in range(len(s)):
        v = cyclic_variant(s, k)
        for j in range(len(v) + 1):
            key = (rank(v[:j]), j)
            if best is None or key > best[0]:
                best = (key, k, v[:j])
    return best[1], best[2]


def raney_proof_trace(s) -> dict:
    """照 Raney 的證明走一遍：α、σ_1、σ_1 = α_1 … α_n、n 個循環 σ_i，以及它們對應的 k。"""
    n = -rank(s)
    k1, alpha = raney_alpha(s)
    sigma1 = cyclic_variant(s, k1)
    words = split_words(sigma1)
    assert len(words) == n
    rotations, ks, pos = [], [], 0
    for i in range(n):
        rotations.append(tuple(x for w in words[i:] + words[:i] for x in w))
        ks.append((k1 + pos) % len(s))
        pos += len(words[i])
    return {"n": n, "k1": k1, "alpha": alpha, "sigma1": sigma1, "words": words,
            "rotations": rotations, "ks": ks}


# ---------------------------------------------------------------------------
# 3. 定理 2.3（Raney 公式）與定理 2.4（看第一個符號）
# ---------------------------------------------------------------------------


def split_head(s, m: int):
    """(m+n) 字列唯一的「m 字列開頭」：回傳 (前 m 個字, 其餘)。"""
    words = split_words(s)
    head = tuple(x for w in words[:m] for x in w)
    return head, tuple(s[len(head):])


def sub_sequences(a) -> list:
    return list(product(*[range(x + 1) for x in a]))


def thm23_rhs(m: int, n: int, a) -> int:
    return sum(L_formula(m, b) * L_formula(n, tuple(x - y for x, y in zip(a, b)))
               for b in sub_sequences(a))


def thm24_rhs(k: int, a) -> int:
    out = L_formula(k, a)
    for i, ai in enumerate(a, start=1):
        if ai > 0:
            out += L_formula(k + i, tuple(x - (j == i) for j, x in enumerate(a, start=1)))
    return out


# ---------------------------------------------------------------------------
# 4. 和 Chu 講義的路徑互換
# ---------------------------------------------------------------------------


def to_chu_path(s, v_vec) -> tuple:
    """Raney 字串 → Chu 講義（counting.py 第 4 節）的路徑字母：倒過來讀，0 → LEAF(=0)，符號 v_i → i。"""
    index = {v: i + 1 for i, v in enumerate(v_vec)}
    return tuple(0 if x == 0 else index[x] for x in reversed(s))


def from_chu_path(path, v_vec) -> tuple:
    return tuple(0 if x == 0 else v_vec[x - 1] for x in reversed(path))


def chu_heights(s) -> list:
    """Chu 講義路徑的高度（倒過來讀、每步 1 − j）。等於「後綴的秩」取負號。"""
    out = [0]
    for x in reversed(s):
        out.append(out[-1] + 1 - x)
    return out


def partition_to_word(blocks, N: int) -> tuple:
    """Chu 的不交叉分法 → Raney 的字（「固定一個點」的分解，Chu 講義第 4 節）。
    點是 0..N−1（講義上的 1..N）。點 0 所在的 v 點群組把其餘的點切成 v 段（可以是空的），
    根 = 符號 v，接著依順時針順序寫出 v 段各自的字；空段是一片葉子 0。"""
    blocks = [tuple(sorted(b)) for b in blocks]

    def rec(points):
        if not points:
            return (0,)
        first = points[0]
        block = next(b for b in blocks if first in b)
        cuts = list(block) + [points[-1] + 1]
        out = (len(block),)
        for a, b in zip(cuts, cuts[1:]):
            out += rec([p for p in points if a < p < b])
        return out

    return rec(list(range(N)))


def word_to_partition(s, N: int) -> list:
    """partition_to_word 的反函數。"""
    blocks = []

    def rec(points, s):
        if not points:
            assert s[0] == 0
            return s[1:]
        v, s = s[0], s[1:]
        block = [points[0]]
        rest = points[1:]
        for j in range(v):
            # 第 j 段：從 rest 開頭往後，長度由子字串決定
            sub_len = _segment_points(s)
            seg, rest = rest[:sub_len], rest[sub_len:]
            s = rec(seg, s)
            if j < v - 1:
                block.append(rest[0])
                rest = rest[1:]
        blocks.append(tuple(block))
        return s

    rest = rec(list(range(N)), tuple(s))
    assert rest == ()
    return sorted(blocks)


def _segment_points(s) -> int:
    """s 開頭那個字（一段）用掉幾個點 = 字裡各符號之和（每個 v 點群組佔 v 個點，葉子不佔點）。"""
    return sum(s[:_word_len(s)])


def _word_len(s) -> int:
    h = 0
    for j, x in enumerate(s):
        h += x - 1
        if h == -1:
            return j + 1
    raise ValueError("沒有完整的字")


def gessel_suffix_code(trees) -> tuple:
    """Gessel §4.5 的 suffix code：c(T) = c(T_1)…c(T_j) j（先寫子樹、最後寫根的小孩數）。"""
    out = ()
    for tr in trees:
        k, kids = tr
        out += gessel_suffix_code(kids) + (k,)
    return out


def prefix_code(trees) -> tuple:
    """Raney 的前序碼（n 字列）。"""
    out = ()
    for tr in trees:
        out += tree_to_word(tr)
    return out


def a_from_chu(n_vec, v_vec) -> tuple:
    """Chu 的 (n⃗, v⃗) → Raney 的 (a_1, …, a_K)，K = max v_i。"""
    K = max(v_vec)
    a = [0] * K
    for n, v in zip(n_vec, v_vec):
        a[v - 1] += n
    return tuple(a)


# ---------------------------------------------------------------------------
# 5. 截斷冪級數（係數用 Fraction）
# ---------------------------------------------------------------------------


def ps_mul(f, g, N: int) -> list:
    out = [Fraction(0)] * (N + 1)
    for i, fi in enumerate(f[:N + 1]):
        if fi:
            for j, gj in enumerate(g[:N + 1 - i]):
                out[i + j] += fi * gj
    return out


def ps_pow(f, m: int, N: int) -> list:
    out = [Fraction(1)] + [Fraction(0)] * N
    for _ in range(m):
        out = ps_mul(out, f, N)
    return out


def ps_compose(phi, w, N: int) -> list:
    """φ(w)，w 的常數項為 0。"""
    out = [Fraction(0)] * (N + 1)
    power = [Fraction(1)] + [Fraction(0)] * N
    for c in phi[:N + 1]:
        if c:
            out = [o + c * p for o, p in zip(out, power)]
        power = ps_mul(power, w, N)
    return out


def solve_w(phi, N: int) -> list:
    """Raney 式 (5.3)：w = x·Σ a_n w^n 的唯一解（逐次代入）。"""
    w = [Fraction(0)] * (N + 1)
    for _ in range(N + 1):
        phw = ps_compose(phi, w, N)
        w = [Fraction(0)] + phw[:N]
    return w


def lagrange_coeff(phi, n: int, t: int = 1) -> Fraction:
    """Lagrange（–Bürmann）：[x^n] w^t = (t/n) [x^{n−t}] φ^n（n ≥ 1）。t = 1 就是 Raney 式 (5.5)。"""
    if n - t < 0:
        return Fraction(0)
    return Fraction(t, n) * ps_pow(phi, n, n - t)[n - t]


def exp_series(c, N: int) -> list:
    """e^{c x} 的係數。"""
    return [Fraction(c) ** k / factorial(k) for k in range(N + 1)]


# ---------------------------------------------------------------------------
# 6. 多變數：L*（Raney 第 4 節）
# ---------------------------------------------------------------------------


def lstar_coeffs(K: int, max_size: int) -> dict:
    """L*(p) = L(1; p)，p = (p_1, …, p_K)，只取 Σ i p_i ≤ max_size。"""
    out = {}
    for p in product(*[range(max_size // i + 1) for i in range(1, K + 1)]):
        if sum(i * x for i, x in enumerate(p, start=1)) <= max_size:
            out[p] = L_formula(1, p)
    return out


def mv_mul(f: dict, g: dict, max_size: int) -> dict:
    out = {}
    for p, x in f.items():
        for q, y in g.items():
            r = tuple(a + b for a, b in zip(p, q))
            if sum(i * c for i, c in enumerate(r, start=1)) <= max_size:
                out[r] = out.get(r, 0) + x * y
    return out


def mv_pow(f: dict, n: int, K: int, max_size: int) -> dict:
    out = {tuple([0] * K): 1}
    for _ in range(n):
        out = mv_mul(out, f, max_size)
    return out


# ---------------------------------------------------------------------------
# 7. 卷積公式（Raney 第 7 節）
# ---------------------------------------------------------------------------


def a_table(phi, M: int, N: int) -> dict:
    """a(m, n) = [x^n] φ^m，0 ≤ m ≤ M、0 ≤ n ≤ N。"""
    out = {}
    for m in range(M + 1):
        pw = ps_pow(phi, m, N)
        for n in range(N + 1):
            out[(m, n)] = pw[n]
    return out


def b_from_a(phi, m: int, n: int) -> Fraction:
    """Raney 式 (7.2)：b(m, n) = m/(m+n)·a(m+n, n)，b(0, 0) = 1。"""
    if m == 0 and n == 0:
        return Fraction(1)
    return Fraction(m, m + n) * ps_pow(phi, m + n, n)[n]


def c_from_b(phi, alpha, k: int) -> Fraction:
    """Raney 式 (7.5)：c(α, k) = Σ_m α^m b(m, k−m)/m!。"""
    alpha = Fraction(alpha)
    return sum(alpha ** m * b_from_a(phi, m, k - m) / factorial(m) for m in range(k + 1))


def gould_G(alpha, n: int, beta, gamma) -> Fraction:
    """Raney 式 (7.7)：G(α, 0; β, γ) = 1，G(α, n; β, γ) = α/n!·Π_{m=1}^{n−1} (α + βn − γm)。"""
    alpha, beta, gamma = Fraction(alpha), Fraction(beta), Fraction(gamma)
    if n == 0:
        return Fraction(1)
    out = alpha / factorial(n)
    for m in range(1, n):
        out *= alpha + beta * n - gamma * m
    return out


def gould_phi(beta, gamma, N: int) -> list:
    """Raney 式 (7.9)：φ(x) = γ x e^{βx} / (e^{γx} − 1)（γ = 0 時是 e^{βx}）。"""
    beta, gamma = Fraction(beta), Fraction(gamma)
    if gamma == 0:
        return exp_series(beta, N)
    # (e^{γx} − 1)/(γx) = Σ γ^k x^k/(k+1)!，取倒數
    den = [gamma ** k / factorial(k + 1) for k in range(N + 1)]
    inv = [Fraction(0)] * (N + 1)
    inv[0] = 1 / den[0]
    for k in range(1, N + 1):
        inv[k] = -sum(den[j] * inv[k - j] for j in range(1, k + 1)) / den[0]
    return ps_mul(exp_series(beta, N), inv, N)


# ---------------------------------------------------------------------------
# 8. 第 6 節與 Problem 4983
# ---------------------------------------------------------------------------


def trinomial_coeff(n: int, k: int) -> Fraction:
    """Raney 式 (6.4)：w = 1 + x w^n 的 x^k 係數 = 1/(1+(n−1)k)·C(nk, k)。"""
    return Fraction(comb(n * k, k), 1 + (n - 1) * k)


def kelly_P(n: int, r: int) -> int:
    """Problem 4983（Klamkin–Shepp，Kelly 解）：n 個因子、每次 r+1 個相乘的非結合乘積個數。"""
    if (n - 1) % r:
        return 0
    k = (n - 1) // r
    value, rem = divmod(comb((r + 1) * k, k), n)
    assert rem == 0
    return value


def kelly_P_brute(n: int, r: int) -> int:
    """直接數：因子 a_1 … a_n，每一對括號恰好包 r+1 個項。"""
    from functools import lru_cache

    @lru_cache(None)
    def P(n):
        if n == 1:
            return 1
        total = 0
        for parts in _compositions_pos(n, r + 1):
            prod = 1
            for x in parts:
                prod *= P(x)
            total += prod
        return total

    return P(n)


def _compositions_pos(n: int, k: int):
    if k == 1:
        if n >= 1:
            yield (n,)
        return
    for first in range(1, n - k + 2):
        for rest in _compositions_pos(n - first, k - 1):
            yield (first,) + rest


def polygon_dissections_brute(sides: int, cell: int) -> int:
    """凸 sides 邊形用不交叉的對角線切成全是 cell 邊形的方法數（暴力：選對角線集合）。"""
    vs = range(sides)
    diags = [(i, j) for i in vs for j in vs if i < j and (j - i) % sides not in (1, sides - 1)]

    def cross(d, e):
        (a, b), (c, d2) = d, e
        if len({a, b, c, d2}) < 4:
            return False
        return (a < c < b) != (a < d2 < b)

    count = 0
    need = (sides - 2) // (cell - 2) - 1  # 切成 (sides−2)/(cell−2) 塊，要這麼多條對角線
    if (sides - 2) % (cell - 2):
        return 0
    from itertools import combinations
    for ds in combinations(diags, need):
        if any(cross(d, e) for d, e in combinations(ds, 2)):
            continue
        if all(len(face) == cell for face in _faces(sides, ds)):
            count += 1
    return count


def _faces(sides: int, diags) -> list:
    faces = [list(range(sides))]
    for (a, b) in diags:
        for f in faces:
            if a in f and b in f:
                ia, ib = f.index(a), f.index(b)
                if ia > ib:
                    ia, ib = ib, ia
                faces.remove(f)
                faces.append(f[ia:ib + 1])
                faces.append(f[ib:] + f[:ia + 1])
                break
    return faces


# ---------------------------------------------------------------------------
# 9. 整數序列的 Raney 引理（notes/ran_note015_改寫.tex；動畫 R01–R06）
#
#    序列 x_1..x_M 是一條路徑的步，S_i = 第 i 點的高度（S_0 = 0）。
#    「從第 i 點開始讀」= 依序讀 x_{i+1}, …, x_{i+M}（循環），i = 0..M−1。
#    好讀法：出發後每點的（相對）高度都 > 0。
# ---------------------------------------------------------------------------

NOTE015_SEQ = (2, -1, -3, 1, 2)          # 講義的主例：和為 1


def heights(seq) -> list:
    """[S_0, S_1, …, S_M]，S_0 = 0。"""
    out = [0]
    for x in seq:
        out.append(out[-1] + x)
    return out


def reading(seq, i: int) -> tuple:
    """從第 i 點開始讀一圈：x_{i+1}, …, x_{i+M}。"""
    seq = tuple(seq)
    return seq[i:] + seq[:i]


def positive_count(seq) -> int:
    """出發後高度 > 0 的點數。"""
    return sum(1 for h in heights(seq)[1:] if h > 0)


def is_good_reading(seq) -> bool:
    return positive_count(seq) == len(seq)


def good_starts(seq) -> list:
    """所有好讀法的起點 i。"""
    return [i for i in range(len(seq)) if is_good_reading(reading(seq, i))]


def last_minimum(seq) -> int:
    """S_1..S_M 中最後一個最小值的位置 x（1 起算）。從第 x 點開始讀就是好讀法（和為 1 時）。"""
    S = heights(seq)
    m = min(S[1:])
    return max(i for i in range(1, len(seq) + 1) if S[i] == m)


def positive_count_order(seq) -> list:
    """加強版的排序法：1..M 依 (S_k, −k) 由小到大。排第 r 位的起點有 M − r + 1 個正的點。"""
    S = heights(seq)
    return sorted(range(1, len(seq) + 1), key=lambda k: (S[k], -k))


def periodic_heights(seq, periods: int) -> list:
    """重複 periods 個週期的高度 [S_0, …, S_{periods·M}]。"""
    return heights(tuple(seq) * periods)


def stack_ops_to_steps(ops: str) -> tuple:
    """進棧 'P' → +1、出棧 'O' → −1。"""
    return tuple(1 if c == "P" else -1 for c in ops)


def stack_sequences(n: int) -> list:
    """所有合法的進出棧操作（n 次進、n 次出、出棧時棧不空），字典序。"""
    out = []

    def rec(s, pushed, popped):
        if len(s) == 2 * n:
            out.append(s)
            return
        if pushed < n:
            rec(s + "P", pushed + 1, popped)
        if popped < pushed:
            rec(s + "O", pushed, popped + 1)

    rec("", 0, 0)
    return out


def pop_order(ops: str) -> tuple:
    """操作序列 → 出棧順序（進棧順序 1, 2, …）。"""
    stack, out, nxt = [], [], 1
    for c in ops:
        if c == "P":
            stack.append(nxt)
            nxt += 1
        else:
            out.append(stack.pop())
    return tuple(out)


def draws_left(cards) -> list:
    """摸牌遊戲：一開始能摸 1 張，摸到 w 元的牌變化 w − 1。回傳 [1, 摸完第 1 張後, …]。
    （= Raney 的「還欠幾個引數」= 開頭的秩 + 1。）"""
    out = [1]
    for w in cards:
        out.append(out[-1] + w - 1)
    return out


def deck_ok(cards) -> bool:
    """前 m 張牌（勝利牌在最底下，不列入）能否摸到勝利牌：每張之後都還能再摸。"""
    return all(d >= 1 for d in draws_left(cards)[1:])


def card_game_brute(w, m: int) -> list:
    """所有合格的牌序（牌有標號）：特殊牌 S1..Sn（w_i），無用牌 U1..U(m−n)。回傳標號 tuple 的串列。"""
    from itertools import permutations
    n = len(w)
    labels = [f"S{i + 1}" for i in range(n)] + [f"U{i + 1}" for i in range(m - n)]
    value = dict(zip(labels, list(w) + [0] * (m - n)))
    return [p for p in permutations(labels) if deck_ok([value[c] for c in p])]


def card_game_formula(w, m: int) -> int:
    """m!/(m − n + 1)，n = w_i ≥ 1 的張數。"""
    n = sum(1 for x in w if x >= 1)
    value, rem = divmod(factorial(m), m - n + 1)
    assert rem == 0
    return value


# ---------------------------------------------------------------------------
# 自我檢查
# ---------------------------------------------------------------------------


def _self_test():
    # Rosenbloom 判別法 = 遞迴定義
    for m in range(0, 8):
        for s in product(range(4), repeat=m):
            r = rank(s)
            for n in range(0, 4):
                assert is_list(s, n) == is_list_recursive(s, n), (s, n)
            if r == -1:
                assert is_word(s) == _is_word_recursive(s)
    # 字 ↔ 樹
    for s in lists_with(1, (1, 2, 1)):
        assert tree_to_word(word_to_tree(s)) == s
    # 原文的三個例子（第 441 頁）
    assert (2, 2, 0, 0, 0) == tree_to_word((2, ((2, ((0, ()), (0, ()))), (0, ()))))
    for s in [(2, 2, 0, 0, 0), (2, 0, 2, 0, 0), (2, 3, 0, 2, 0, 0, 0, 1, 0)]:
        assert is_word(s)

    # 定理 2.1（循環引理）與 Raney 證明的每一步
    for n, a in [(1, (0, 2)), (2, (0, 1, 0)), (3, (1, 1, 1)), (2, (2, 0, 1)), (1, (2, 2)), (4, (0, 3))]:
        for s in strings_with(n, a):
            ks = list_shifts(s)
            assert len(ks) == n, (s, ks)
            tr = raney_proof_trace(s)
            assert sorted(tr["ks"]) == ks
            assert all(is_list(x, n) for x in tr["rotations"])
            assert all(is_word(w) for w in tr["words"])

    # 定理 2.2 = 暴力
    tested = 0
    for n in range(0, 4):
        for a in product(range(4), range(3), range(2)):
            if length_of(n, a) <= 10:
                assert L_formula(n, a) == L_brute(n, a), (n, a)
                tested += 1
    assert tested > 40
    # Cayley：L(1; 0, k) = C_k
    for k in range(8):
        assert L_formula(1, (0, k)) == comb(2 * k, k) // (k + 1)

    # 定理 2.3、2.4
    for a in product(range(3), range(3), range(2)):
        for m in range(0, 3):
            for n in range(0, 3):
                assert L_formula(m + n, a) == thm23_rhs(m, n, a)
        for k in range(0, 3):
            assert L_formula(k + 1, a) == thm24_rhs(k, a)
    # 定理 2.3 的切法是一一對應
    for s in lists_with(3, (1, 1)):
        h, rest = split_head(s, 1)
        assert is_list(h, 1) and is_list(rest, 2)

    # 和 Chu 講義的路徑：倒過來讀的字串 = 好路徑
    import combkit.counting as ct  # noqa: E402
    for t, n_vec, v_vec in [(1, (1, 1), (2, 3)), (2, (1, 1), (2, 3)), (1, (2, 1), (2, 3)), (3, (1, 2), (2, 3))]:
        a = a_from_chu(n_vec, v_vec)
        lists = lists_with(t, a)
        paths = ct.good_paths(t, n_vec, v_vec)
        assert sorted(to_chu_path(s, v_vec) for s in lists) == sorted(paths)
        assert len(lists) == ct.raney_A(t, n_vec, v_vec) == L_formula(t, a)
        for s in lists:
            assert from_chu_path(to_chu_path(s, v_vec), v_vec) == s
            assert chu_heights(s) == ct.path_heights(to_chu_path(s, v_vec), v_vec)
    # Chu 的不交叉分法 ↔ Raney 的字：一一對應，而且保持組成
    for n_vec, v_vec in [((1, 1), (2, 3)), ((2, 1), (2, 3)), ((3,), (2,)), ((3,), (3,)), ((1, 1, 1), (2, 3, 4)),
                         ((2, 2), (1, 2))]:
        N = ct.inner(n_vec, v_vec)
        parts = ct.enumerate_noncrossing(n_vec, v_vec)
        words = [partition_to_word(p, N) for p in parts]
        assert len(set(words)) == len(parts) == L_formula(1, a_from_chu(n_vec, v_vec))
        assert sorted(words) == lists_with(1, a_from_chu(n_vec, v_vec))
        for p, w in zip(parts, words):
            assert word_to_partition(w, N) == sorted(tuple(sorted(b)) for b in p)
    # Gessel 圖 3 的森林：suffix code 0 0 2 0 1
    fig3 = [(2, ((0, ()), (0, ()))), (1, ((0, ()),))]
    assert gessel_suffix_code(fig3) == (0, 0, 2, 0, 1) and prefix_code(fig3) == (2, 0, 0, 1, 0)
    # 一元運算子 ↔ Chu 的 1 點群組（單點區塊）：公式照樣成立
    for n_vec, v_vec in [((1, 1), (1, 2)), ((2, 1), (1, 2)), ((1, 1, 1), (1, 2, 3)), ((2, 2), (1, 2))]:
        assert ct.brute_force(n_vec, v_vec) == L_formula(1, a_from_chu(n_vec, v_vec)) == ct.closed_form(n_vec, v_vec)

    # 第 4 節：L* 的方程與冪
    K, S = 3, 7
    Ls = lstar_coeffs(K, S)
    rhs = {tuple([0] * K): 1}
    for i in range(1, K + 1):
        pw = mv_pow(Ls, i, K, S)
        for p, x in pw.items():
            q = tuple(c + (j == i) for j, c in enumerate(p, start=1))
            if sum(j * c for j, c in enumerate(q, start=1)) <= S:
                rhs[q] = rhs.get(q, 0) + x
    for p in Ls:
        assert Ls[p] == rhs.get(p, 0), p
    for n in range(4):
        pw = mv_pow(Ls, n, K, S)
        for p in Ls:
            assert pw.get(p, 0) == L_formula(n, p)

    # 第 5 節：Lagrange 反演
    N = 9
    examples = {
        "1/(1-x)": [Fraction(1)] * (N + 1),
        "(1+x)^2": [Fraction(1), Fraction(2), Fraction(1)] + [Fraction(0)] * (N - 2),
        "1+x+x^2": [Fraction(1), Fraction(1), Fraction(1)] + [Fraction(0)] * (N - 2),
        "e^x": exp_series(1, N),
        "2+3x+x^3": [Fraction(2), Fraction(3), Fraction(0), Fraction(1)] + [Fraction(0)] * (N - 3),
    }
    for name, phi in examples.items():
        w = solve_w(phi, N)
        for n in range(1, N + 1):
            assert w[n] == lagrange_coeff(phi, n), (name, n)
            for t in range(1, 4):
                assert ps_pow(w, t, N)[n] == lagrange_coeff(phi, n, t), (name, n, t)
    # (5.6)：(1/n)[x^{n−1}] z^n = Σ L(1; p) a^p（取 z = 2 + 3x + x^3 代入檢查）
    for n in range(1, 8):
        total = Fraction(0)
        phi = examples["2+3x+x^3"]
        for p in product(range(n), range(n), range(n)):
            if sum(i * c for i, c in enumerate(p, start=1)) == n - 1:
                p0 = 1 + sum((i - 1) * c for i, c in enumerate(p, start=1))
                total += L_formula(1, p) * phi[0] ** p0 * phi[1] ** p[0] * phi[2] ** p[1] * phi[3] ** p[2]
        assert total == lagrange_coeff(phi, n)

    # 第 6 節
    for n in range(2, 6):
        w = [Fraction(1)] + [Fraction(0)] * 8
        for _ in range(10):
            w = [Fraction(1)] + ps_mul([Fraction(0), Fraction(1)], ps_pow(w, n, 8), 8)[1:]
        for k in range(9):
            assert w[k] == trinomial_coeff(n, k) == L_formula(1, tuple(k if i == n else 0 for i in range(1, n + 1)))
    for r in range(1, 4):
        for n in range(1, 12):
            assert kelly_P(n, r) == kelly_P_brute(n, r)
    for r in (1, 2):
        for k in range(1, 4):
            n = r * k + 1
            assert polygon_dissections_brute(n + 1, r + 2) == kelly_P(n, r)

    # 第 7 節：定理 7.1、式 (7.6)、Gould (7.8)
    for phi in (examples["1/(1-x)"], examples["(1+x)^2"], examples["2+3x+x^3"]):
        for m1 in range(3):
            for m2 in range(3):
                for n in range(5):
                    lhs = b_from_a(phi, m1 + m2, n)
                    rhs = sum(b_from_a(phi, m1, n1) * b_from_a(phi, m2, n - n1) for n1 in range(n + 1))
                    assert lhs == rhs
    for beta, gamma in [(1, 0), (2, 1), (Fraction(3, 2), Fraction(1, 2)), (1, 1)]:
        phi = gould_phi(beta, gamma, 7)
        for alpha in (1, Fraction(2, 3), -2):
            for n in range(7):
                assert c_from_b(phi, alpha, n) == gould_G(alpha, n, beta, gamma), (beta, gamma, alpha, n)
        for a1, a2 in [(1, 2), (Fraction(1, 3), Fraction(5, 2))]:
            for n in range(7):
                assert gould_G(a1 + a2, n, beta, gamma) == sum(
                    gould_G(a1, k, beta, gamma) * gould_G(a2, n - k, beta, gamma) for k in range(n + 1))
    # γ = 1、α = t、β = v：就是 Chu 講義的 A_t(n)（只有一種群組）
    for v in (2, 3):
        for t in (1, 2, 3):
            for n in range(6):
                assert gould_G(t, n, v, 1) == ct.raney_A(t, (n,), (v,))
    print("Raney (1960)：字、循環引理、L(n;a)、L*、Lagrange 反演、第 6–7 節全部驗證通過")

    # 第 9 節：整數序列的 Raney 引理、加強版、兩題例題（notes/ran_note015_改寫.tex）
    import random
    rng = random.Random(15)
    for M in range(1, 8):
        for _ in range(1500):
            a = [rng.randint(-5, 5) for _ in range(M - 1)]
            a.append(1 - sum(a))
            assert len(good_starts(a)) == 1 and good_starts(a) == [last_minimum(a) % M]
            assert sorted(positive_count(reading(a, i)) for i in range(M)) == list(range(1, M + 1))
            for r, k in enumerate(positive_count_order(a), 1):
                assert positive_count(reading(a, k % M)) == M - r + 1
    a = NOTE015_SEQ
    assert heights(a) == [0, 2, 1, -2, -1, 1] and last_minimum(a) == 3 and good_starts(a) == [3]
    assert [positive_count(reading(a, i)) for i in range(5)] == [3, 1, 2, 5, 4]
    assert positive_count_order(a) == [3, 4, 5, 2, 1]
    assert periodic_heights(a, 3)[5::5] == [1, 2, 3]
    assert [positive_count(reading((Fraction(1, 2), Fraction(1, 2)), i)) for i in range(2)] == [2, 2]
    assert [positive_count(reading((2, 0), i)) for i in range(2)] == [2, 1]
    for n in range(0, 7):
        ops = stack_sequences(n)
        assert len(ops) == comb(2 * n, n) // (n + 1) == len({pop_order(o) for o in ops})
        for o in ops:
            assert is_good_reading((1,) + stack_ops_to_steps(o))
    for w, m in [((3,), 3), ((1, 2), 3), ((2, 2), 4), ((1, 1, 2), 4), ((2, 3), 5)]:
        decks = card_game_brute(w, m)
        assert len(decks) == card_game_formula(w, m)
        vals = dict(zip([f"S{i + 1}" for i in range(len(w))], w))
        for d in decks:
            word = tuple(vals.get(c, 0) for c in d) + (0,)
            assert is_word(word) and [h + 1 for h in head_ranks(word)][:-1] == draws_left(word[:-1])
    assert card_game_formula((0, 3), 3) == len(card_game_brute((0, 3), 3)) == 2
    print("Raney 引理（整數序列）、加強版、出棧、摸牌全部驗證通過")


if __name__ == "__main__":
    import sys
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
    _self_test()
