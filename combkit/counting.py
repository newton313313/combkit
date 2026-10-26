"""
counting.py —— 純 Python 的計數工具（不依賴 manim）

用三種彼此獨立的方法計算「圓周上不交叉分組」的數目，互相驗證：
  1. closed_form      : 論文式 (2) 的廣義卡塔蘭數公式
  2. brute_force      : 列舉所有分組方式，再過濾掉交叉的
  3. recurrence       : 論文式 (5) 的遞迴（固定一點的分解）
第 5 節是經典卡塔蘭物件（場景 00a–00m），第 6 節是 MIT 18.212 L1–4 的補充（場景 M01–M14）。

直接執行本檔會印出驗證表：
    python utils/counting.py
"""
from __future__ import annotations

from fractions import Fraction
from functools import lru_cache
from itertools import combinations, permutations, product
from math import factorial


# ---------------------------------------------------------------------------
# 1. 公式
# ---------------------------------------------------------------------------
def catalan(n: int) -> int:
    """經典卡塔蘭數 C_n = binom(2n, n) / (n + 1)。"""
    return factorial(2 * n) // (factorial(n) * factorial(n + 1))


def inner(n_vec, v_vec) -> int:
    """內積 <n, v> = sum n_i v_i，也就是圓周上的總點數 N。"""
    return sum(n * v for n, v in zip(n_vec, v_vec))


def closed_form(n_vec, v_vec) -> int:
    """
    論文式 (2)：
        C_v(n) = 1/(N+1) * (N+1)! / (n_1! ... n_k! * L!),
        N = <n, v>,  L = 1 + sum n_i (v_i - 1)
    """
    N = inner(n_vec, v_vec)
    L = 1 + sum(n * (v - 1) for n, v in zip(n_vec, v_vec))
    denom = factorial(L)
    for n in n_vec:
        denom *= factorial(n)
    value, rem = divmod(factorial(N + 1), denom * (N + 1))
    assert rem == 0
    return value


# ---------------------------------------------------------------------------
# 2. 暴力列舉
# ---------------------------------------------------------------------------
def blocks_cross(a, b) -> bool:
    """
    兩個群組 a, b（點的編號集合）是否交叉。
    交叉 <=> 沿圓周存在 a1 < b1 < a2 < b2 交錯出現（a1,a2 屬於 a；b1,b2 屬於 b）。
    """
    for a1, a2 in combinations(sorted(a), 2):
        inside = [x for x in b if a1 < x < a2]
        outside = [x for x in b if x < a1 or x > a2]
        if inside and outside:
            return True
    return False


def is_noncrossing(blocks) -> bool:
    return not any(blocks_cross(a, b) for a, b in combinations(blocks, 2))


def _partitions_with_sizes(points: tuple, sizes: tuple):
    """把 points 分成若干組，組的大小（不計順序）恰為 sizes 這個多重集合。"""
    if not points:
        if not sizes:
            yield []
        return
    first, rest = points[0], points[1:]
    for size in sorted(set(sizes)):
        remaining_sizes = list(sizes)
        remaining_sizes.remove(size)
        for others in combinations(rest, size - 1):
            block = (first,) + others
            left = tuple(p for p in rest if p not in others)
            for tail in _partitions_with_sizes(left, tuple(remaining_sizes)):
                yield [block] + tail


def enumerate_noncrossing(n_vec, v_vec):
    """列出所有不交叉分組。回傳 list，每個元素是 list[tuple]（各群組的點編號，0 起算）。"""
    N = inner(n_vec, v_vec)
    sizes = tuple(v for n, v in zip(n_vec, v_vec) for _ in range(n))
    return [
        part for part in _partitions_with_sizes(tuple(range(N)), sizes)
        if is_noncrossing(part)
    ]


def brute_force(n_vec, v_vec) -> int:
    return len(enumerate_noncrossing(n_vec, v_vec))


# ---------------------------------------------------------------------------
# 3. 遞迴（論文式 (5)）
# ---------------------------------------------------------------------------
def _compositions(target: tuple, parts: int):
    """把向量 target 拆成 parts 個非負整數向量的有序和。"""
    if parts == 1:
        yield (target,)
        return
    ranges = [range(t + 1) for t in target]
    for head in product(*ranges):
        rest = tuple(t - h for t, h in zip(target, head))
        for tail in _compositions(rest, parts - 1):
            yield (head,) + tail


def recurrence(n_vec, v_vec) -> int:
    """d(n) = sum_i sum_{r_1+...+r_{v_i} = n - e_i} d(r_1) ... d(r_{v_i})，d(0) = 1。"""
    v_vec = tuple(v_vec)

    @lru_cache(maxsize=None)
    def d(n: tuple) -> int:
        if all(x == 0 for x in n):
            return 1
        total = 0
        for i, v in enumerate(v_vec):
            if n[i] == 0:
                continue
            reduced = tuple(x - (1 if j == i else 0) for j, x in enumerate(n))
            for rs in _compositions(reduced, v):
                term = 1
                for r in rs:
                    term *= d(r)
                total += term
        return total

    return d(tuple(n_vec))


# ---------------------------------------------------------------------------
# 4. Raney 公式（論文式 (9)）的證明用到的「路徑」
#
#    一條路徑是一串「步」：
#        LEAF（↑ 步）：高度 +1
#        第 i 型步   ：高度 1 - v_i（弦 -1、三角形 -2）
#    用 n_i 個第 i 型步，加上 ell = t + sum n_i (v_i - 1) 個 ↑ 步，
#    總步數 N + t，終點高度恰好是 t。
#    「好路徑」：出發之後每一點的高度都 > 0。
# ---------------------------------------------------------------------------
LEAF = 0


def step_weight(letter: int, v_vec) -> int:
    return 1 if letter == LEAF else 1 - v_vec[letter - 1]


def path_heights(word, v_vec, start: int = 0) -> list:
    """回傳 [S_0, S_1, ..., S_M]，S_0 = start。"""
    hs = [start]
    for letter in word:
        hs.append(hs[-1] + step_weight(letter, v_vec))
    return hs


def is_good_path(word, v_vec) -> bool:
    return all(h > 0 for h in path_heights(word, v_vec)[1:])


def first_bad_index(word, v_vec):
    """第一個高度 <= 0 的點（不含起點）；好路徑回傳 None。"""
    for j, h in enumerate(path_heights(word, v_vec)):
        if j > 0 and h <= 0:
            return j
    return None


def multinomial_vec(M: int, n_vec) -> int:
    """binom(M; n) = M! / (n_1! ... n_k! (M - |n|)!)；不合法時回傳 0。"""
    rest = M - sum(n_vec)
    if rest < 0 or any(n < 0 for n in n_vec):
        return 0
    out = factorial(M) // factorial(rest)
    for n in n_vec:
        out //= factorial(n)
    return out


def raney_A(t: int, n_vec, v_vec) -> int:
    """A_t(n) = t / (N + t) * binom(N + t; n)。A_1 就是廣義卡塔蘭數。"""
    if any(n < 0 for n in n_vec):
        return 0
    M = inner(n_vec, v_vec) + t
    value, rem = divmod(t * multinomial_vec(M, n_vec), M)
    assert rem == 0
    return value


def word_content(word, k: int) -> tuple:
    """路徑用了幾個第 1, ..., k 型步。"""
    return tuple(sum(1 for x in word if x == i + 1) for i in range(k))


def path_words(t: int, n_vec, v_vec) -> list:
    """所有排法（依字典序排列，不含重複）。"""
    ell = t + sum(n * (v - 1) for n, v in zip(n_vec, v_vec))
    letters = [LEAF] * ell
    for i, n in enumerate(n_vec):
        letters += [i + 1] * n
    out = set()

    def rec(prefix, counts):
        if sum(counts.values()) == 0:
            out.add(tuple(prefix))
            return
        for letter in sorted(counts):
            if counts[letter]:
                counts[letter] -= 1
                prefix.append(letter)
                rec(prefix, counts)
                prefix.pop()
                counts[letter] += 1

    counts = {}
    for x in letters:
        counts[x] = counts.get(x, 0) + 1
    rec([], counts)
    return sorted(out)


def good_paths(t: int, n_vec, v_vec) -> list:
    return [w for w in path_words(t, n_vec, v_vec) if is_good_path(w, v_vec)]


def rotate(word, j: int) -> tuple:
    """從第 j 步開始讀（循環）。"""
    word = tuple(word)
    return word[j:] + word[:j]


def good_rotations(word, v_vec) -> list:
    """M 種讀法中，哪幾種是好路徑（回傳起點 j）。"""
    return [j for j in range(len(word)) if is_good_path(rotate(word, j), v_vec)]


def cut_index(word, s: int, v_vec) -> int:
    """最後一個高度 <= s 的點；好路徑在這裡切成「終點 s」與「終點 t」兩段。"""
    hs = path_heights(word, v_vec)
    return max(j for j, h in enumerate(hs) if h <= s)


def sub_vectors(n_vec):
    """所有 0 <= r <= n（逐分量）的向量 r，依字典序。"""
    return list(product(*[range(n + 1) for n in n_vec]))


# ---------------------------------------------------------------------------
# 5. 經典卡塔蘭物件（卡塔蘭基礎模組，場景 00a–00m 用）
#
#    所有物件都由同一種「樹」產生，因此彼此之間有固定的一一對應：
#        樹 t：None 表示空（葉子），(L, R) 表示一個內部節點，左右子樹為 L、R。
#        n 個內部節點的滿二元樹恰有 C_n 棵。
#    排列順序固定（左子樹由小到大），不使用亂數。
#
#    下面每一種物件都另外用「暴力列舉」獨立算一次，在自我檢查裡比對。
# ---------------------------------------------------------------------------
@lru_cache(maxsize=None)
def binary_trees(n: int) -> tuple:
    """n 個內部節點的所有滿二元樹，依「左子樹大小」由小到大排列。"""
    if n == 0:
        return (None,)
    out = []
    for k in range(n):
        for left in binary_trees(k):
            for right in binary_trees(n - 1 - k):
                out.append((left, right))
    return tuple(out)


def tree_size(t) -> int:
    """內部節點的個數。"""
    return 0 if t is None else 1 + tree_size(t[0]) + tree_size(t[1])


# ----------------------------------------------------------------- 括號
def tree_to_paren_tokens(t, factors="abcdefghij"):
    """
    完全括號化（最外層不加括號），用 · 表示每一次乘法。
    回傳 (tokens, parts)：
        tokens 例如 ['(', 'a', '·', 'b', ')', '·', '(', 'c', '·', 'd', ')']
        parts  每個 token 屬於 'root'（最後一次乘法）、'left' 或 'right'
    乘號依中序排列，正好對應樹的內部節點（中序）；因子依序對應葉子。
    """
    letters = iter(factors)

    def rec(node, wrap):
        if node is None:
            return [next(letters)]
        inner = rec(node[0], True) + ["·"] + rec(node[1], True)
        return ["("] + inner + [")"] if wrap else inner

    if t is None:
        return [next(letters)], ["root"]
    left = rec(t[0], True)
    right = rec(t[1], True)
    tokens = left + ["·"] + right
    parts = ["left"] * len(left) + ["root"] + ["right"] * len(right)
    return tokens, parts


def parenthesizations(n: int, factors="abcdefghij") -> list:
    """n+1 個因子的所有完全括號化（字串，最外層不加括號）。"""
    return ["".join(tree_to_paren_tokens(t, factors)[0]) for t in binary_trees(n)]


def parenthesizations_bruteforce(seq: str) -> set:
    """獨立的列舉：把因子串切成左右兩段、各自遞迴（字串形式與上面相同）。"""
    if len(seq) == 1:
        return {seq}
    out = set()
    for cut in range(1, len(seq)):
        for a in parenthesizations_bruteforce(seq[:cut]):
            for b in parenthesizations_bruteforce(seq[cut:]):
                wa = a if len(a) == 1 else f"({a})"
                wb = b if len(b) == 1 else f"({b})"
                out.add(f"{wa}·{wb}")
    return out


# ----------------------------------------------------------------- 三角剖分
def tree_to_triangulation(t, i: int, j: int) -> list:
    """
    凸多邊形頂點 i, i+1, ..., j（底邊為 i–j）的三角剖分，回傳對角線串列。
    根的三角形是 (i, k, j)，k = i + size(L) + 1；左邊多邊形 i..k、右邊 k..j。
    """
    if t is None:
        assert j == i + 1
        return []
    k = i + tree_size(t[0]) + 1
    diags = []
    if k - i >= 2:
        diags.append((i, k))
    if j - k >= 2:
        diags.append((k, j))
    return diags + tree_to_triangulation(t[0], i, k) + tree_to_triangulation(t[1], k, j)


def triangulation_apex(t, m: int) -> int:
    """凸 m 邊形（頂點 0..m-1，底邊 0–(m-1)）中，含底邊的那個三角形的第三個頂點。"""
    return tree_size(t[0]) + 1


def triangulations(m: int) -> list:
    """凸 m 邊形的所有三角剖分（對角線的 tuple），依「含底邊三角形的頂點」由左到右排列。"""
    return [tuple(sorted(tree_to_triangulation(t, 0, m - 1))) for t in binary_trees(m - 2)]


def diagonals_cross(a, b) -> bool:
    (p, q), (r, s) = sorted(a), sorted(b)
    return p < r < q < s or r < p < s < q


def triangulations_bruteforce(m: int) -> set:
    """獨立的列舉：從所有對角線裡挑 m-3 條互不相交的。"""
    diags = [(a, b) for a in range(m) for b in range(a + 2, m) if not (a == 0 and b == m - 1)]
    out = set()
    for combo in combinations(diags, m - 3):
        if all(not diagonals_cross(x, y) for x, y in combinations(combo, 2)):
            out.add(tuple(sorted(combo)))
    return out


# ----------------------------------------------------------------- Dyck 路徑與格子路徑
def tree_to_dyck(t) -> str:
    """Dyck 路徑（U 上坡、D 下坡）：t = (L, R) 對應 U·dyck(L)·D·dyck(R)。"""
    return "" if t is None else "U" + tree_to_dyck(t[0]) + "D" + tree_to_dyck(t[1])


def dyck_parts(t) -> list:
    """Dyck 路徑每一步屬於 'root'（第一個上坡與它對應的下坡）、'left' 或 'right'。"""
    if t is None:
        return []
    a, b = len(tree_to_dyck(t[0])), len(tree_to_dyck(t[1]))
    return ["root"] + ["left"] * a + ["root"] + ["right"] * b


def is_dyck(word: str) -> bool:
    h = 0
    for c in word:
        h += 1 if c == "U" else -1
        if h < 0:
            return False
    return h == 0


def dyck_words(n: int) -> list:
    """獨立的列舉：所有 n 個 U、n 個 D、從不低於 0 的字（字典序）。"""
    out = []
    for ups in combinations(range(2 * n), n):
        word = "".join("U" if i in ups else "D" for i in range(2 * n))
        if is_dyck(word):
            out.append(word)
    return sorted(out)


def dyck_to_lattice(word: str) -> str:
    """Dyck 上坡 -> 格子路徑往右 R；下坡 -> 往上 U。（不越過 y=x 的路徑。）"""
    return "".join("R" if c == "U" else "U" for c in word)


def lattice_words(n: int) -> list:
    """(0,0) 到 (n,n) 的所有單調格子路徑（n 個 R、n 個 U），字典序。"""
    out = []
    for rs in combinations(range(2 * n), n):
        out.append("".join("R" if i in rs else "U" for i in range(2 * n)))
    return sorted(out)


def lattice_points(word: str, start=(0, 0)) -> list:
    x, y = start
    pts = [(x, y)]
    for c in word:
        if c == "R":
            x += 1
        else:
            y += 1
        pts.append((x, y))
    return pts


def lattice_ok(word: str) -> bool:
    """從不走到對角線 y = x 的上方（可以碰到對角線）。"""
    return all(y <= x for x, y in lattice_points(word))


def first_bad_step(word: str):
    """第一次走到 y = x + 1 的那一步之後的位置（點的編號）；合法路徑回傳 None。"""
    for j, (x, y) in enumerate(lattice_points(word)):
        if y > x:
            return j
    return None


def reflect_tail(word: str) -> str:
    """
    反射原理：第一次碰到 y = x + 1 之後的部分，沿著 y = x + 1 翻轉（R 與 U 互換）。
    不合法的路徑（終點 (n,n)）因此一一對應到終點 (n-1, n+1) 的路徑。
    """
    j = first_bad_step(word)
    assert j is not None
    tail = "".join("U" if c == "R" else "R" for c in word[j:])
    return word[:j] + tail


# ----------------------------------------------------------------- 括號字串、投票、配對
def dyck_to_brackets(word: str) -> str:
    return "".join("(" if c == "U" else ")" for c in word)


def dyck_to_ballot(word: str) -> str:
    return "".join("A" if c == "U" else "B" for c in word)


def ballot_ok(word: str) -> bool:
    """卡塔蘭特例：開票過程中 A 的票數始終 >= B 的票數，而且最後同票。"""
    return is_dyck("".join("U" if c == "A" else "D" for c in word))


def bertrand_strict_count(a: int, b: int) -> int:
    """暴力列舉：A 得 a 票、B 得 b 票，開票過程中 A 始終「嚴格」領先的順序數。"""
    count = 0
    for pos in combinations(range(a + b), a):
        lead, ok = 0, True
        for i in range(a + b):
            lead += 1 if i in pos else -1
            if lead <= 0:
                ok = False
                break
        count += ok
    return count


def ballot_more_than_k_times(a: int, b: int, k: int) -> int:
    """暴力列舉：開票過程中 A 的票數始終「大於」B 的 k 倍的順序數（推廣的投票問題）。"""
    count = 0
    for pos in combinations(range(a + b), a):
        na = nb = 0
        ok = True
        for i in range(a + b):
            if i in pos:
                na += 1
            else:
                nb += 1
            if na <= k * nb:
                ok = False
                break
        count += ok
    return count


def dyck_to_matching(word: str) -> list:
    """每個上坡與它對應的下坡配成一對：圓周 2n 點的不交叉完美配對。"""
    stack, pairs = [], []
    for i, c in enumerate(word):
        if c == "U":
            stack.append(i)
        else:
            pairs.append((stack.pop(), i))
    return sorted(pairs)


def matching_parts(t) -> list:
    """配對裡每條弦屬於 'root'（含第 0 點的弦）、'left'（弦內側）或 'right'（外側）。"""
    word, parts = tree_to_dyck(t), dyck_parts(t)
    return [parts[i] for i, _ in dyck_to_matching(word)]


def series_square_shift(coeffs: list) -> list:
    """回傳 1 + x * C(x)^2 的前 len(coeffs) 個係數。"""
    m = len(coeffs)
    out = [0] * m
    out[0] = 1
    for n in range(1, m):
        out[n] = sum(coeffs[k] * coeffs[n - 1 - k] for k in range(n))
    return out


# ---------------------------------------------------------------------------
# 6. MIT 18.212 L1–4 補充（場景 M01–M14 用）
#
#    約定：
#      * 點、格子、列、行的編號一律 0 起算（畫面要顯示 1 起算時自己 +1）。
#      * 路徑字：'U' 上坡、'D' 下坡（與第 5 節相同）。
#      * 平面樹：一個頂點 = 它的小孩子樹組成的 tuple，葉子是 ()；n 條邊（n+1 個頂點）。
#      * 排列：1..n 的 tuple（也接受 "41532" 這種字串）。
#      * 整數分割 λ：由大到小的 tuple；Young 圖用英式（第 0 列在上、靠左）；
#        Young 表 T：列的 tuple，例如 ((1, 2, 4), (3, 5, 6))。
#      * 機率一律是精確分數 Fraction。
#    每一族物件在自我檢查裡都用兩種獨立的方法計數。
# ---------------------------------------------------------------------------

# ----------------------------------------------------------------- 生成函數（M01、M02）
def series_mul(a, b, terms: int | None = None) -> list:
    """兩個冪級數（係數串列）的 Cauchy 乘積前 terms 項：[x^n](ab) = sum_{k=0}^{n} a_k b_{n-k}。"""
    if terms is None:
        terms = min(len(a), len(b))
    return [sum(a[k] * b[n - k] for k in range(n + 1)) for n in range(terms)]


def binom_general(r, m: int) -> Fraction:
    """廣義二項係數 binom(r, m) = r (r-1) ... (r-m+1) / m!（r 可以是分數，例如 1/2）。"""
    out = Fraction(1)
    for i in range(m):
        out *= Fraction(r) - i
    return out / factorial(m)


def sqrt_1m4x_coeffs(terms: int) -> list:
    """sqrt(1-4x) = sum_m binom(1/2, m) (-4)^m x^m 的前 terms 個係數：1, -2, -2, -4, -10, ...（都是整數）。"""
    out = []
    for m in range(terms):
        c = binom_general(Fraction(1, 2), m) * (-4) ** m
        assert c.denominator == 1
        out.append(int(c))
    return out


def catalan_from_sqrt(terms: int) -> list:
    """C(x) = (1 - sqrt(1-4x)) / (2x)：C_n = -[x^{n+1}] sqrt(1-4x) / 2。"""
    s = sqrt_1m4x_coeffs(terms + 1)
    assert all(c % 2 == 0 for c in s[1:])
    return [-s[n + 1] // 2 for n in range(terms)]


# ----------------------------------------------------------------- 醉漢走路（M03）
def drunkard_fall_words(n: int) -> list:
    """暴力列舉：長度 2n+1 的 U/D 字，從 0 出發、在最後一步「第一次」走到 -1（掉下懸崖）。字典序。"""
    out = []
    for word in product("DU", repeat=2 * n + 1):
        h, first = 0, None
        for j, c in enumerate(word):
            h += 1 if c == "U" else -1
            if h < 0:
                first = j
                break
        if first == 2 * n:
            out.append("".join(word))
    return sorted(out)


def drunkard_fall_prob(n: int) -> Fraction:
    """第 2n+1 步第一次掉下去的機率 = C_n / 2^(2n+1)（Dyck 路徑接一步往下）。"""
    return Fraction(catalan(n), 2 ** (2 * n + 1))


def drunkard_partial_sums(N: int) -> list:
    """[S_1, ..., S_N]，S_K = sum_{n<K} C_n / 2^(2n+1) = 1/2, 5/8, 11/16, ...（精確分數，趨近 1）。"""
    out, s = [], Fraction(0)
    for n in range(N):
        s += drunkard_fall_prob(n)
        out.append(s)
    return out


# ----------------------------------------------------------------- 古典循環引理（M04）
def cycle_words(n: int) -> list:
    """
    n+1 個 +1（字母 0 = LEAF）與 n 個 -1（字母 1）的所有排法，長度 2n+1，共 binom(2n+1, n) 個。
    就是 path_words(1, (n,), (2,))；高度用 path_heights(w, (2,))、好讀法用 good_rotations(w, (2,))。
    """
    return path_words(1, (n,), (2,))


def cycle_classes(n: int) -> list:
    """
    依循環移位把 cycle_words(n) 分組：[(好的字, [rotate(好的字, j) for j in 0..2n]), ...]，
    依好的字的字典序。每組 2n+1 個不同的字、恰有一個好的（出發後前綴和都 > 0）。
    """
    v = (2,)
    out = {}
    for w in cycle_words(n):
        good = good_rotations(w, v)
        assert len(good) == 1
        g = rotate(w, good[0])
        out.setdefault(g, [rotate(g, j) for j in range(2 * n + 1)])
    return sorted(out.items())


def cycle_good_to_dyck(word) -> str:
    """好的字 = 一個 +1 接一條 Dyck 路徑：去掉第一個字母，0 -> 'U'、1 -> 'D'。"""
    assert word[0] == LEAF
    return "".join("U" if x == LEAF else "D" for x in word[1:])


# ----------------------------------------------------------------- 平面樹（M05、M07、M08）
def binary_to_plane(t) -> tuple:
    """
    滿二元樹（n 個內部節點）-> 平面樹（n+1 個頂點）：左小孩 = 第一個小孩、右小孩 = 下一個兄弟。
    二元樹的內部節點依序變成平面樹的非根頂點，平面樹的根是新加的。
    例：((None, (None, None)), None) -> (((), ()),)，也就是 根 -> A -> 兩片葉子。
    """
    def forest(node) -> list:
        out = []
        while node is not None:
            out.append(tuple(forest(node[0])))
            node = node[1]
        return out

    return tuple(forest(t))


def plane_to_binary(pt):
    """binary_to_plane 的反函數。"""
    if not pt:
        return None
    return (plane_to_binary(pt[0]), plane_to_binary(pt[1:]))


def plane_tree_edges(pt) -> int:
    """邊數（= 頂點數 - 1）。"""
    return sum(1 + plane_tree_edges(c) for c in pt)


def plane_trees(n: int) -> tuple:
    """n 條邊的所有平面樹，順序與 binary_trees(n) 相同（所以和 00k 的欄位一致）。"""
    return tuple(binary_to_plane(t) for t in binary_trees(n))


@lru_cache(maxsize=None)
def plane_trees_bruteforce(n: int) -> frozenset:
    """獨立的列舉：根有 k 個小孩，剩下 n-k 條邊任意分給 k 棵小孩子樹。"""
    if n == 0:
        return frozenset({()})
    out = set()
    for k in range(1, n + 1):
        for parts in _compositions((n - k,), k):
            for kids in product(*[plane_trees_bruteforce(p[0]) for p in parts]):
                out.add(tuple(kids))
    return frozenset(out)


def plane_tree_dfs(pt) -> str:
    """深度優先走訪：走向小孩記 'U'、回到父親記 'D'。plane_tree_dfs(binary_to_plane(t)) == tree_to_dyck(t)。"""
    return "".join("U" + plane_tree_dfs(c) + "D" for c in pt)


def plane_tree_preorder(pt) -> list:
    """頂點依前序（DFS 第一次到達的順序）編號 0, 1, ...（0 是根）；回傳 [(父親編號或 None, 深度), ...]。"""
    out = []

    def rec(node, parent, depth):
        idx = len(out)
        out.append((parent, depth))
        for c in node:
            rec(c, idx, depth + 1)

    rec(pt, None, 0)
    return out


def plane_tree_walk(pt) -> list:
    """沿輪廓走一圈依序經過的頂點（前序編號），長度 2n+1；第 k 步從 walk[k] 走到 walk[k+1]。"""
    walk, counter = [], [0]

    def rec(node):
        me = counter[0]
        walk.append(me)
        for c in node:
            counter[0] += 1
            rec(c)
            walk.append(me)

    rec(pt)
    return walk


# ----------------------------------------------------------------- 三角剖分的對偶樹（M07）
def triangulation_dual(t) -> list:
    """
    tree_to_triangulation(t, 0, m-1)（m = n+2，底邊 0–(m-1)）的對偶樹。每個三角形一筆，依「中序」排列
    （與 catalan_draw.BinaryTreeMob.nodes 相同順序；頂點為 k 的三角形是第 k-1 筆）：
        {"tri": (a, k, b), "base": (a, b), "left": 小孩, "right": 小孩, "parent": 編號或 None}
    base = 與父親共用的對角線（根三角形是底邊 (0, m-1)）；
    小孩 = ("tri", 三角形編號) 或 ("leaf", e)：葉子 e = 多邊形的邊 (e, e+1)（中序 = 由左到右 = BinaryTreeMob.leaves）。
    """
    n = tree_size(t)
    nodes = [None] * n

    def rec(node, a, b, parent):
        if node is None:
            assert b == a + 1
            return ("leaf", a)
        k = a + tree_size(node[0]) + 1
        idx = k - 1
        left = rec(node[0], a, k, idx)
        right = rec(node[1], k, b, idx)
        nodes[idx] = {"tri": (a, k, b), "base": (a, b), "left": left, "right": right, "parent": parent}
        return ("tri", idx)

    rec(t, 0, n + 1, None)
    return nodes


def triangulation_to_tree(m: int, diagonals):
    """
    獨立的反方向：只看對角線重建對偶二元樹。從底邊 (0, m-1) 出發，找含這條邊的三角形 (a, k, b)，
    往 (a, k)、(k, b) 兩邊遞迴；碰到多邊形的邊 (e, e+1) 就是葉子。
    """
    sides = {tuple(sorted(d)) for d in diagonals} | {(i, i + 1) for i in range(m - 1)}

    def rec(a, b):
        if b == a + 1:
            return None
        apex = [k for k in range(a + 1, b) if (a, k) in sides and (k, b) in sides]
        assert len(apex) == 1, (a, b, apex)
        return (rec(a, apex[0]), rec(apex[0], b))

    return rec(0, m - 1)


# ----------------------------------------------------------------- 直線上的完美配對與弧圖（M06、M08）
def perfect_matchings(n_points: int) -> list:
    """n_points 個點（0 起算）的所有完美配對，共 (n_points-1)!! 個；每個是排序好的 ((i, j), ...)，i < j。"""
    def rec(pts):
        if not pts:
            yield ()
            return
        a = pts[0]
        for idx in range(1, len(pts)):
            rest = pts[1:idx] + pts[idx + 1:]
            for m in rec(rest):
                yield ((a, pts[idx]),) + m

    return [tuple(sorted(m)) for m in rec(tuple(range(n_points)))]


def arcs_cross(p, q) -> bool:
    """直線上兩條弧交叉 <=> 端點交錯 i < j < k < l：(i, k)、(j, l)。"""
    (a, b), (c, d) = sorted(p), sorted(q)
    return a < c < b < d or c < a < d < b


def arcs_nest(p, q) -> bool:
    """直線上兩條弧巢狀 <=> 一條在另一條裡面 i < j < k < l：(i, l)、(j, k)。"""
    (a, b), (c, d) = sorted(p), sorted(q)
    return a < c < d < b or c < a < b < d


def matching_crossings(pairs) -> list:
    """交叉的弧：[(第 x 條, 第 y 條), ...]（pairs 裡的位置）。"""
    return [(x, y) for x, y in combinations(range(len(pairs)), 2) if arcs_cross(pairs[x], pairs[y])]


def matching_nestings(pairs) -> list:
    """巢狀的弧：[(外面那條, 裡面那條), ...]（pairs 裡的位置）。"""
    out = []
    for x, y in combinations(range(len(pairs)), 2):
        if arcs_nest(pairs[x], pairs[y]):
            out.append((x, y) if min(pairs[x]) < min(pairs[y]) else (y, x))
    return out


def matching_is_noncrossing(pairs) -> bool:
    return not matching_crossings(pairs)


def matching_is_nonnesting(pairs) -> bool:
    return not matching_nestings(pairs)


def dyck_to_nonnesting(word: str) -> list:
    """Dyck 路徑 -> 不巢狀配對：上坡依序排隊，每個下坡配給最早、還沒配的上坡（先進先出）。"""
    queue, pairs = [], []
    for i, c in enumerate(word):
        if c == "U":
            queue.append(i)
        else:
            pairs.append((queue.pop(0), i))
    return sorted(pairs)


def matching_to_dyck(pairs) -> str:
    """配對 -> 字：每條弧的左端記 'U'、右端記 'D'（不交叉或不巢狀配對都適用）。"""
    lefts = {min(p) for p in pairs}
    return "".join("U" if i in lefts else "D" for i in range(2 * len(pairs)))


# ----------------------------------------------------------------- 堆疊／佇列排序（M09）
def stack_sort(perm) -> tuple:
    """
    用一個堆疊排序，目標輸出 1, 2, ..., n（mit_plan.md 第 4 節 M09 的規則），每一步：
        堆疊頂端正好是下一個要輸出的數 -> 彈出到輸出      ('pop', x)
        否則把下一個輸入推入堆疊                            ('push', x)
        輸入用完而頂端不對                                  -> 失敗
    回傳 (ok, ops)；失敗時 ops 到卡住為止。可排序 <=> 避開 231。
    """
    inp, stack, ops, need = list(_as_tuple(perm)), [], [], 1
    while need <= len(perm):
        if stack and stack[-1] == need:
            ops.append(("pop", stack.pop()))
            need += 1
        elif inp:
            x = inp.pop(0)
            stack.append(x)
            ops.append(("push", x))
        else:
            return False, ops
    return True, ops


def queue_bypass_sort(perm) -> tuple:
    """
    用「佇列＋插隊線」排序（mit_plan.md 第 1 節第 5 點），每一步：
        佇列頭正好是下一個要輸出的數 -> 出列到輸出            ('pop', x)
        否則看下一個輸入：正好是下一個要輸出的 -> 走插隊線      ('bypass', x)
                          不是 -> 排到佇列尾                   ('push', x)
        輸入用完而佇列頭不對                                    -> 失敗
    回傳 (ok, ops)。可排序 <=> 避開 321。
    """
    inp, queue, ops, need = list(_as_tuple(perm)), [], [], 1
    while need <= len(perm):
        if queue and queue[0] == need:
            ops.append(("pop", queue.pop(0)))
            need += 1
        elif inp:
            x = inp.pop(0)
            if x == need:
                ops.append(("bypass", x))
                need += 1
            else:
                queue.append(x)
                ops.append(("push", x))
        else:
            return False, ops
    return True, ops


def sort_states(perm, ops) -> list:
    """
    重播 ops：回傳每一步之後的 (輸入, 容器, 輸出)，第 0 筆是開始時。
    容器：堆疊由底到頂、佇列由頭到尾（'pop' 拿走的一定是堆疊頂端／佇列頭）。
    """
    inp, box, out = list(_as_tuple(perm)), [], []
    states = [(tuple(inp), tuple(box), tuple(out))]
    for op, x in ops:
        if op == "push":
            assert inp[0] == x
            box.append(inp.pop(0))
        elif op == "pop":
            assert x in (box[0], box[-1])
            box.remove(x)
            out.append(x)
        elif op == "bypass":
            assert inp[0] == x
            out.append(inp.pop(0))
        else:
            raise ValueError(op)
        states.append((tuple(inp), tuple(box), tuple(out)))
    return states


# ----------------------------------------------------------------- 模式迴避（M10）
PATTERNS_3 = [tuple(p) for p in permutations((1, 2, 3))]     # 123, 132, 213, 231, 312, 321


def _as_tuple(w) -> tuple:
    return tuple(int(c) for c in w) if isinstance(w, str) else tuple(w)


def standardize(seq) -> tuple:
    """相對大小：(4, 5, 3) -> (2, 3, 1)。"""
    seq = _as_tuple(seq)
    ranks = sorted(seq)
    return tuple(ranks.index(x) + 1 for x in seq)


def contains_pattern(w, p):
    """
    w 含模式 p（有子序列和 p 的相對大小相同）時，回傳第一個出現的位置（0 起算的 index tuple，
    依 combinations 的順序），否則回傳 None。例：contains_pattern("41532", "231") == (0, 2, 3)（4, 5, 3）。
    """
    w, p = _as_tuple(w), _as_tuple(p)
    checks = [(a, b, p[a] < p[b]) for a, b in combinations(range(len(p)), 2)]
    for idx in combinations(range(len(w)), len(p)):
        if all((w[idx[a]] < w[idx[b]]) == lt for a, b, lt in checks):
            return idx
    return None


def avoiders(n: int, p) -> list:
    """1..n 的排列中避開 p 的（依字典序）。"""
    return [w for w in permutations(range(1, n + 1)) if contains_pattern(w, p) is None]


# ----------------------------------------------------------------- 整數分割、合成、Young 圖（M11–M14）
def integer_partitions(n: int) -> list:
    """n 的所有整數分割（由大到小的 tuple），反字典序：(5,), (4, 1), (3, 2), (3, 1, 1), ...。"""
    def rec(rem, mx):
        if rem == 0:
            yield ()
            return
        for first in range(min(rem, mx), 0, -1):
            for rest in rec(rem - first, first):
                yield (first,) + rest

    return list(rec(n, n))


def partition_count(n: int) -> int:
    """獨立的算法：Euler 五角數定理的遞迴 p(m) = sum_k (-1)^(k+1) [p(m - k(3k-1)/2) + p(m - k(3k+1)/2)]。"""
    p = [1] + [0] * n
    for m in range(1, n + 1):
        k, s = 1, 0
        while k * (3 * k - 1) // 2 <= m:
            sign = 1 if k % 2 else -1
            s += sign * p[m - k * (3 * k - 1) // 2]
            if k * (3 * k + 1) // 2 <= m:
                s += sign * p[m - k * (3 * k + 1) // 2]
            k += 1
        p[m] = s
    return p[n]


def integer_compositions(n: int) -> list:
    """n 的所有合成（有順序的正整數和），共 2^(n-1) 個（n >= 1）。字典序。"""
    if n == 0:
        return [()]
    out = []
    for cuts in product((0, 1), repeat=n - 1):      # 在 n 個 1 之間的 n-1 個空隙要不要切
        parts, cur = [], 1
        for c in cuts:
            if c:
                parts.append(cur)
                cur = 1
            else:
                cur += 1
        parts.append(cur)
        out.append(tuple(parts))
    return sorted(out)


def conjugate(lam) -> tuple:
    """共軛分割 λ'（Young 圖轉置）：λ'_j = 第 j 行的長度。"""
    lam = tuple(lam)
    return tuple(sum(1 for r in lam if r > j) for j in range(lam[0])) if lam else ()


def young_cells(lam) -> list:
    """所有格子 (i, j)（第 i 列、第 j 行，0 起算），逐列由左到右。"""
    return [(i, j) for i, r in enumerate(lam) for j in range(r)]


def hook_cells(lam, i: int, j: int) -> list:
    """鉤上的格子：[(i, j) 自己, 同列右邊（arm）..., 同行下面（leg）...]。"""
    conj = conjugate(lam)
    assert 0 <= i < len(lam) and 0 <= j < lam[i]
    return ([(i, j)] + [(i, jj) for jj in range(j + 1, lam[i])]
            + [(ii, j) for ii in range(i + 1, conj[j])])


def hook_length(lam, i: int, j: int) -> int:
    """鉤長 h(i, j) = arm + leg + 1（0 起算）。"""
    conj = conjugate(lam)
    assert 0 <= i < len(lam) and 0 <= j < lam[i]
    return (lam[i] - j - 1) + (conj[j] - i - 1) + 1


def hook_lengths(lam) -> list:
    """每一列的鉤長，例如 (4, 3, 1) -> [(6, 4, 3, 1), (4, 2, 1), (1,)]。"""
    return [tuple(hook_length(lam, i, j) for j in range(r)) for i, r in enumerate(lam)]


def hook_product(lam) -> int:
    """H(λ) = 所有格子鉤長的乘積。"""
    out = 1
    for i, j in young_cells(lam):
        out *= hook_length(lam, i, j)
    return out


def hook_formula(lam) -> int:
    """鉤長公式 f^λ = n! / H(λ)。"""
    value, rem = divmod(factorial(sum(lam)), hook_product(lam))
    assert rem == 0
    return value


def corners(lam) -> list:
    """角落（拿掉後仍是 Young 圖的格子）：[(i, λ_i - 1), ...]，由上到下。"""
    lam = tuple(lam)
    return [(i, r - 1) for i, r in enumerate(lam) if i == len(lam) - 1 or lam[i + 1] < r]


def remove_corner(lam, c) -> tuple:
    """λ - c：拿掉角落 c（列長變 0 就去掉那一列）。"""
    i, j = c
    assert c in corners(lam)
    out = list(lam)
    out[i] -= 1
    return tuple(r for r in out if r > 0)


def syt_count(lam) -> int:
    """角落遞迴：最大的數 n 一定在角落 -> f^λ = sum_{角落 c} f^{λ-c}，f^() = 1。（λ 可以是 list 或 tuple）"""
    return _syt_count(tuple(lam))


@lru_cache(maxsize=None)
def _syt_count(lam: tuple) -> int:
    if not lam:
        return 1
    return sum(_syt_count(remove_corner(lam, c)) for c in corners(lam))


def syt(lam) -> list:
    """
    形狀 λ 的所有標準 Young 表（每列 →、每行 ↓ 都遞增），依「列字」的字典序。
    例：syt((3, 2)) 第一個是 ((1, 2, 3), (4, 5))。
    """
    lam = tuple(lam)
    words = []

    def rec(counts, word):
        if len(word) == sum(lam):
            words.append(tuple(word))
            return
        for r in range(len(lam)):
            if counts[r] < lam[r] and (r == 0 or counts[r - 1] > counts[r]):
                counts[r] += 1
                word.append(r)
                rec(counts, word)
                word.pop()
                counts[r] -= 1

    rec([0] * len(lam), [])
    return [rowword_to_syt(w) for w in words]


def syt_bruteforce(lam) -> set:
    """獨立的列舉：把 1..n 的每一種排列逐列填入，留下列、行都遞增的。"""
    lam = tuple(lam)
    cells = young_cells(lam)
    out = set()
    for perm in permutations(range(1, len(cells) + 1)):
        val = dict(zip(cells, perm))
        if all(val[(i, j)] < val[(i, j + 1)] for i, j in cells if (i, j + 1) in val) and \
           all(val[(i, j)] < val[(i + 1, j)] for i, j in cells if (i + 1, j) in val):
            out.add(tuple(tuple(val[(i, j)] for j in range(r)) for i, r in enumerate(lam)))
    return out


def syt_to_rowword(T) -> tuple:
    """列字：第 k 個字母 = 數字 k+1 在第幾列（0 起算）。"""
    where = {x: r for r, row in enumerate(T) for x in row}
    return tuple(where[x] for x in range(1, len(where) + 1))


def rowword_to_syt(word) -> tuple:
    """syt_to_rowword 的反函數：依序把 1, 2, ... 放到列字指定那一列的尾端。"""
    rows = [[] for _ in range(max(word) + 1)] if word else []
    for k, r in enumerate(word):
        rows[r].append(k + 1)
    return tuple(tuple(row) for row in rows)


def rowword_counts(word, k: int | None = None) -> list:
    """
    累計次數表（k 維格子路徑）：第 t 筆 = 讀完前 t 個字母後各列的個數 (x_0, ..., x_{k-1})，第 0 筆全為 0。
    是 Young 表 <=> 每一筆都有 x_0 >= x_1 >= ... >= x_{k-1}（見 is_lattice_rowword）。
    """
    if k is None:
        k = max(word) + 1 if word else 0
    x = [0] * k
    out = [tuple(x)]
    for r in word:
        x[r] += 1
        out.append(tuple(x))
    return out


def is_lattice_rowword(word) -> bool:
    return all(all(a >= b for a, b in zip(x, x[1:])) for x in rowword_counts(word))


def syt_to_dyck(T) -> str:
    """形狀 (n, n)：數字 i 在第一列 -> 第 i 步 'U'、在第二列 -> 'D'。"""
    assert len(T) == 2 and len(T[0]) == len(T[1])
    return "".join("U" if r == 0 else "D" for r in syt_to_rowword(T))


def dyck_to_syt(word: str) -> tuple:
    """syt_to_dyck 的反函數：'UUDUDD' -> ((1, 2, 4), (3, 5, 6))。"""
    return rowword_to_syt(tuple(0 if c == "U" else 1 for c in word))


# ----------------------------------------------------------------- 鉤上漫步（M13、M14）
def corner_prob_formula(lam, c) -> Fraction:
    """(1/n) H(λ) / H(λ - c)。對所有角落加總 = 1 <=> 鉤長公式（對 n 歸納）。"""
    return Fraction(hook_product(lam), sum(lam) * hook_product(remove_corner(lam, c)))


def hook_walks(lam, start) -> list:
    """
    從格子 start 出發的所有鉤上漫步：每步跳到目前格子鉤上的另一格（機率 1/(h-1)），停在角落（h = 1）。
    回傳 [(格子 tuple, 權重), ...]，權重 = prod 1/(h(a_i) - 1)（不含選起點的 1/n）。
    """
    out = []

    def rec(path, w):
        i, j = path[-1]
        h = hook_length(lam, i, j)
        if h == 1:
            out.append((tuple(path), w))
            return
        for c in hook_cells(lam, i, j)[1:]:
            rec(path + [c], w / (h - 1))

    rec([tuple(start)], Fraction(1))
    return out


def hook_walk_probs(lam) -> dict:
    """{角落: 停在這個角落的機率}（精確分數）：起點均勻 1/n，列舉所有漫步加總。等於 corner_prob_formula。"""
    n = sum(lam)
    probs = {c: Fraction(0) for c in corners(lam)}
    for s in young_cells(lam):
        for path, w in hook_walks(lam, s):
            probs[path[-1]] += w / n
    return probs


def walk_rows_cols(path) -> tuple:
    """漫步的投影：(經過的列集合 A, 經過的行集合 B)，都是排序好的 tuple。"""
    return tuple(sorted({i for i, _ in path})), tuple(sorted({j for _, j in path}))


def hook_identity(lam, cell, c) -> tuple:
    """
    關鍵恆等式（c = (α, β) 是角落、cell = (a, b) 在 c 的左上方）：
        h(a, b) - 1 = (h(a, β) - 1) + (h(α, b) - 1)
    回傳三個數 (h(a,b)-1, h(a,β)-1, h(α,b)-1)。
    """
    (a, b), (al, be) = cell, c
    assert c in corners(lam) and a <= al and b <= be
    out = (hook_length(lam, a, b) - 1, hook_length(lam, a, be) - 1, hook_length(lam, al, b) - 1)
    assert out[0] == out[1] + out[2]
    return out


def hook_lemma_rhs(lam, c, A, B) -> Fraction:
    """
    引理的右邊：從 (min A, min B) 出發、經過的列集合 A、行集合 B、停在角落 c = (α, β) 的所有漫步權重和
        = prod_{i ∈ A - {α}} 1/(h(i, β) - 1) · prod_{j ∈ B - {β}} 1/(h(α, j) - 1)。
    """
    al, be = c
    out = Fraction(1)
    for i in A:
        if i != al:
            out /= hook_length(lam, i, be) - 1
    for j in B:
        if j != be:
            out /= hook_length(lam, al, j) - 1
    return out


def corner_prob_product(lam, c) -> Fraction:
    """M14 第 5 段：P(c) = (1/n) prod_{i<α} (1 + 1/(h(i,β)-1)) · prod_{j<β} (1 + 1/(h(α,j)-1))。"""
    al, be = c
    out = Fraction(1, sum(lam))
    for i in range(al):
        out *= 1 + Fraction(1, hook_length(lam, i, be) - 1)
    for j in range(be):
        out *= 1 + Fraction(1, hook_length(lam, al, j) - 1)
    return out


# ----------------------------------------------------------------- 貫穿的例子（mit_plan.md 第 3 節）
MIT_EXAMPLE_TREE = ((None, (None, None)), None)     # = binary_trees(3)[3]，00k 第 4 欄


# ---------------------------------------------------------------------------
# 自我檢查
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    cases = [
        ((1,), (2,)), ((2,), (2,)), ((3,), (2,)), ((4,), (2,)),   # 經典 Catalan: 1, 2, 5, 14
        ((2,), (3,)), ((3,), (3,)), ((2,), (4,)),                  # 高階 (Fuss-)Catalan
        ((1, 1), (2, 3)),                                          # 論文 Fig. 1  -> 5
        ((2, 1), (2, 3)), ((1, 2), (2, 3)), ((2, 2), (2, 3)),
        ((3, 2), (2, 3)),                                          # Scene05 的例子 (12 點)
        ((1, 1, 1), (2, 3, 4)), ((1, 1), (1, 3)), ((2, 1), (3, 4)),
    ]
    print(f"{'n':>12} {'v':>12} {'N':>3} {'公式':>7} {'暴力':>7} {'遞迴':>7}  結果")
    ok = True
    for n_vec, v_vec in cases:
        N = inner(n_vec, v_vec)
        f = closed_form(n_vec, v_vec)
        r = recurrence(n_vec, v_vec)
        b = brute_force(n_vec, v_vec) if N <= 12 else None
        good = (f == r) and (b is None or b == f)
        ok &= good
        print(f"{str(n_vec):>12} {str(v_vec):>12} {N:>3} {f:>7} {str(b):>7} {r:>7}  {'OK' if good else 'MISMATCH'}")
    assert [catalan(i) for i in range(6)] == [1, 1, 2, 5, 14, 42]
    assert all(closed_form((m,), (2,)) == catalan(m) for m in range(8))
    print("全部一致" if ok else "有不一致！")

    # ----------------------------------------------- Raney 公式證明的每一步
    print("\nRaney 公式的證明（暴力列舉檢查每一步）")
    raney_cases = [((1, 1), (2, 3)), ((2, 1), (2, 3)), ((1, 2), (2, 3)),
                   ((2,), (2,)), ((3,), (2,)), ((2,), (3,)), ((1, 1, 1), (2, 3, 4)),
                   ((1, 1), (1, 3))]
    for n_vec, v_vec in raney_cases:
        k = len(n_vec)
        for t in (1, 2, 3):
            words = path_words(t, n_vec, v_vec)
            N = inner(n_vec, v_vec)
            # (a) 排法總數 = binom(N+t; n)
            assert len(words) == multinomial_vec(N + t, n_vec)
            # (b) 循環引理：每條路徑恰有 t 種讀法是好路徑
            assert all(len(good_rotations(w, v_vec)) == t for w in words)
            # (c) 好路徑的數目 = A_t(n)
            assert len(good_paths(t, n_vec, v_vec)) == raney_A(t, n_vec, v_vec)
        # (d) A_1 = 廣義卡塔蘭數
        assert raney_A(1, n_vec, v_vec) == closed_form(n_vec, v_vec)
        for s, t in ((1, 1), (1, 2), (2, 1)):
            # (e) 切開是一一對應：好(s+t) = 好(s) 接 好(t)
            glued = set()
            for r in sub_vectors(n_vec):
                rest = tuple(a - b for a, b in zip(n_vec, r))
                for a in good_paths(s, r, v_vec):
                    for b in good_paths(t, rest, v_vec):
                        glued.add(a + b)
            target = good_paths(s + t, n_vec, v_vec)
            assert glued == set(target)
            assert all(cut_index(a + b, s, v_vec) == len(a)
                       for r in sub_vectors(n_vec)
                       for a in good_paths(s, r, v_vec)
                       for b in good_paths(t, tuple(x - y for x, y in zip(n_vec, r)), v_vec))
            # (f) Raney 公式本身
            conv = sum(raney_A(s, r, v_vec) * raney_A(t, tuple(a - b for a, b in zip(n_vec, r)), v_vec)
                       for r in sub_vectors(n_vec))
            assert raney_A(s + t, n_vec, v_vec) == conv
        print(f"  n={n_vec}, v={v_vec}: A_1..A_3 = "
              f"{[raney_A(t, n_vec, v_vec) for t in (1, 2, 3)]}  OK")
    print("Raney 公式證明的每一步都通過")

    # ----------------------------------------------- 卡塔蘭基礎模組（場景 00a–00m）
    print("\n卡塔蘭基礎模組（每種物件用兩種方法列舉，並檢查一一對應）")
    from math import comb
    assert [catalan(i) for i in range(11)] == [1, 1, 2, 5, 14, 42, 132, 429, 1430, 4862, 16796]
    for n in range(0, 7):
        cn = catalan(n)
        trees = binary_trees(n)
        assert len(trees) == len(set(trees)) == cn
        assert all(tree_size(t) == n for t in trees)
        if n >= 1:
            # Segner 遞迴、二項式差
            assert cn == sum(catalan(k) * catalan(n - 1 - k) for k in range(n))
            assert comb(2 * n, n) - comb(2 * n, n - 1) == cn
            # 凸 (n+2) 邊形的三角剖分
            tri = triangulations(n + 2)
            assert len(set(tri)) == cn and set(tri) == triangulations_bruteforce(n + 2)
            # n+1 個因子的完全括號化
            par = parenthesizations(n)
            assert len(set(par)) == cn and set(par) == parenthesizations_bruteforce("abcdefghij"[: n + 1])
        # Dyck 路徑（由樹產生 vs. 暴力列舉）
        dy = [tree_to_dyck(t) for t in trees]
        assert len(set(dy)) == cn and sorted(dy) == dyck_words(n)
        assert all(len(dyck_parts(t)) == 2 * n for t in trees)
        # 格子路徑：全部 binom(2n,n)、合法 C_n、不合法經反射一一對應到 binom(2n,n-1)
        lw = lattice_words(n)
        good = [w for w in lw if lattice_ok(w)]
        bad = [w for w in lw if not lattice_ok(w)]
        assert len(lw) == comb(2 * n, n) and len(good) == cn
        assert sorted(dyck_to_lattice(d) for d in dy) == sorted(good)
        if n >= 1:
            targets = set()
            for rs in combinations(range(2 * n), n - 1):
                targets.add("".join("R" if i in rs else "U" for i in range(2 * n)))
            images = [reflect_tail(w) for w in bad]
            assert len(set(images)) == len(bad) == comb(2 * n, n - 1)
            assert set(images) == targets
        # 平衡括號字串、投票（卡塔蘭特例）
        assert len({dyck_to_brackets(d) for d in dy}) == cn
        ballots = []
        for pos in combinations(range(2 * n), n):
            w = "".join("A" if i in pos else "B" for i in range(2 * n))
            if ballot_ok(w):
                ballots.append(w)
        assert sorted(ballots) == sorted(dyck_to_ballot(d) for d in dy)
        # 圓周 2n 點的不交叉完美配對（由樹產生 vs. 第 2 節的暴力列舉）
        mt = {tuple(dyck_to_matching(d)) for d in dy}
        assert len(mt) == cn
        if n >= 1:
            brute = {tuple(sorted(tuple(sorted(b)) for b in part))
                     for part in enumerate_noncrossing((n,), (2,))}
            assert mt == brute
            assert all(len(matching_parts(t)) == n for t in trees)
        print(f"  n={n}: C_n={cn:>4}   三角剖分、括號、二元樹、格子路徑、Dyck、投票、配對  OK")
    # 一般的 Bertrand 投票定理（嚴格領先）：(a-b)/(a+b) * binom(a+b, a)
    for a in range(1, 9):
        for b in range(0, a):
            assert bertrand_strict_count(a, b) * (a + b) == (a - b) * comb(a + b, a)
    # 來源論文寫的 k 倍推廣：a > k b 時，A 始終多於 B 的 k 倍的機率是 (a - k b)/(a + b)
    for k in (1, 2, 3):
        for a in range(1, 11):
            for b in range(0, 6):
                if a > k * b:
                    assert ballot_more_than_k_times(a, b, k) * (a + b) == (a - k * b) * comb(a + b, a)
    # 生成函數：C = 1 + x C^2
    cs = [catalan(k) for k in range(15)]
    assert series_square_shift(cs) == cs
    # 場景裡用到的具體例子
    assert parenthesizations(3) == ["a·(b·(c·d))", "a·((b·c)·d)", "(a·b)·(c·d)", "(a·(b·c))·d", "((a·b)·c)·d"]
    assert [tree_to_dyck(t) for t in binary_trees(3)] == ["UDUDUD", "UDUUDD", "UUDDUD", "UUDUDD", "UUUDDD"]
    assert [triangulation_apex(t, 6) for t in binary_trees(4)] == [1] * 5 + [2] * 2 + [3] * 2 + [4] * 5
    assert len(lattice_words(3)) == 20 and lattice_ok("RURRUU") and not lattice_ok("RUURRU")
    assert reflect_tail("RUURRU") == "RUUUUR"
    print("卡塔蘭基礎模組：所有物件的數目與一一對應都通過")

    # ----------------------------------------------- MIT 18.212 L1–4 補充（場景 M01–M14）
    print("\nMIT 18.212 L1–4 補充（每一族物件用兩種方法計數，並檢查一一對應）")
    from math import prod
    # 生成函數：Cauchy 乘積、sqrt(1-4x) 的二項級數
    cs = [catalan(k) for k in range(15)]
    sq = series_mul(cs, cs)
    assert sq[2] == 5 and all(sq[n] == cs[n + 1] for n in range(14))          # [x^n]C^2 = C_{n+1}
    assert sqrt_1m4x_coeffs(5) == [1, -2, -2, -4, -10]
    assert all(binom_general(Fraction(1, 2), m) * (-4) ** m == Fraction(-2, m) * comb(2 * m - 2, m - 1)
               for m in range(1, 15))
    s = sqrt_1m4x_coeffs(15)
    assert series_mul(s, s) == [1, -4] + [0] * 13                              # (sqrt(1-4x))^2 = 1-4x
    assert catalan_from_sqrt(14) == cs[:14]
    assert all(catalan(n) * (n + 1) == comb(2 * n, n) for n in range(15))
    # 醉漢走路：第 2n+1 步第一次掉下去的路徑 = Dyck 路徑 + 'D'
    for n in range(0, 7):
        assert drunkard_fall_words(n) == sorted(d + "D" for d in dyck_words(n))
        assert len(drunkard_fall_words(n)) == catalan(n)
    ps = drunkard_partial_sums(30)
    assert ps[:3] == [Fraction(1, 2), Fraction(5, 8), Fraction(11, 16)]
    assert all(1 - ps[K - 1] == Fraction(comb(2 * K, K), 4 ** K) for K in range(1, 31))   # 還沒掉下去的機率 -> 0
    assert all(p < 1 for p in ps)
    # 古典循環引理（v = (2,)）
    for n in range(1, 7):
        words = cycle_words(n)
        classes = cycle_classes(n)
        assert len(words) == comb(2 * n + 1, n) and len(classes) == catalan(n)
        assert sorted(w for _, rots in classes for w in rots) == words
        assert all(len(set(rots)) == 2 * n + 1 for _, rots in classes)
        assert sorted(cycle_good_to_dyck(g) for g, _ in classes) == dyck_words(n)
        assert comb(2 * n + 1, n) == (2 * n + 1) * catalan(n)
    assert len(cycle_words(3)) == 35 and len(cycle_classes(3)) == 5

    for n in range(0, 8):
        cn = catalan(n)
        trees = binary_trees(n)
        # 平面樹
        pts = plane_trees(n)
        assert len(set(pts)) == cn and set(pts) == plane_trees_bruteforce(n)
        assert all(plane_tree_edges(p) == n and plane_to_binary(p) == t for p, t in zip(pts, trees))
        assert all(plane_tree_dfs(p) == tree_to_dyck(t) for p, t in zip(pts, trees))
        for p in pts:
            walk, pre = plane_tree_walk(p), plane_tree_preorder(p)
            assert len(walk) == 2 * n + 1 and len(pre) == n + 1
            assert "".join("U" if pre[b][1] > pre[a][1] else "D" for a, b in zip(walk, walk[1:])) == plane_tree_dfs(p)
        # 三角剖分的對偶樹
        if n >= 1:
            for t in trees:
                m = n + 2
                diags = tree_to_triangulation(t, 0, m - 1)
                assert triangulation_to_tree(m, diags) == t
                dual = triangulation_dual(t)
                assert len(dual) == n and [d["parent"] for d in dual].count(None) == 1
                assert dual[tree_size(t[0])]["base"] == (0, m - 1)
                leaves = [c[1] for d in dual for c in (d["left"], d["right"]) if c[0] == "leaf"]
                assert sorted(leaves) == list(range(m - 1))
                shared = {tuple(sorted(d["base"])) for d in dual if d["parent"] is not None}
                assert shared == set(diags)
        if n <= 6:
            # 直線上 2n 點的完美配對：不交叉、不巢狀
            allm = perfect_matchings(2 * n)
            assert len(allm) == prod(range(1, 2 * n, 2))                           # (2n-1)!!
            nc = {m for m in allm if matching_is_noncrossing(m)}
            nn = {m for m in allm if matching_is_nonnesting(m)}
            dy = dyck_words(n)
            assert len(nc) == len(nn) == cn
            assert nc == {tuple(dyck_to_matching(d)) for d in dy}
            assert nn == {tuple(dyck_to_nonnesting(d)) for d in dy}
            assert all(matching_to_dyck(dyck_to_matching(d)) == d == matching_to_dyck(dyck_to_nonnesting(d))
                       for d in dy)
            if n >= 1:
                assert nc == {tuple(sorted(tuple(sorted(b)) for b in part))
                              for part in enumerate_noncrossing((n,), (2,))}
        # 堆疊／佇列排序、模式迴避
        perms = list(permutations(range(1, n + 1)))
        stack_ok = {p for p in perms if stack_sort(p)[0]}
        queue_ok = {p for p in perms if queue_bypass_sort(p)[0]}
        assert stack_ok == set(avoiders(n, (2, 3, 1)))
        assert queue_ok == set(avoiders(n, (3, 2, 1)))
        assert len(stack_ok) == len(queue_ok) == cn
        for p in perms:
            for sorter in (stack_sort, queue_bypass_sort):
                ok, ops = sorter(p)
                final = sort_states(p, ops)[-1]
                assert (final[2] == tuple(range(1, n + 1))) == ok
        for pat in PATTERNS_3:
            av = set(avoiders(n, pat))
            assert len(av) == cn
            assert all(standardize([p[i] for i in contains_pattern(p, pat)]) == pat
                       for p in perms if p not in av)
        print(f"  n={n}: C_n={cn:>4}   平面樹、對偶樹、弧圖、堆疊／佇列排序、6 種模式迴避  OK")

    # 整數分割、合成、Young 表、鉤長公式、鉤上漫步
    for n in range(0, 8):
        parts = integer_partitions(n)
        assert len(parts) == len(set(parts)) == partition_count(n)
        assert all(sum(p) == n and list(p) == sorted(p, reverse=True) for p in parts)
        comps = integer_compositions(n)
        assert len(comps) == (2 ** (n - 1) if n else 1)
        assert set(comps) == {q for p in parts for q in permutations(p)}
        assert sum(hook_formula(p) ** 2 for p in parts) == factorial(n)              # sum (f^λ)^2 = n!
        for lam in parts:
            if n == 0:
                continue
            assert conjugate(conjugate(lam)) == lam and sum(conjugate(lam)) == n
            tabs = syt(lam)
            assert len(set(tabs)) == len(tabs) == syt_count(lam) == hook_formula(lam)
            assert set(tabs) == syt_bruteforce(lam)
            assert all(rowword_to_syt(syt_to_rowword(T)) == T and is_lattice_rowword(syt_to_rowword(T))
                       for T in tabs)
            # 角落遞迴、化約恆等式、拿掉角落只改 c 左方與上方格子的鉤長（各減 1）
            cs_ = corners(lam)
            assert syt_count(lam) == sum(hook_formula(remove_corner(lam, c)) for c in cs_)
            assert sum(corner_prob_formula(lam, c) for c in cs_) == 1
            for c in cs_:
                sub = remove_corner(lam, c)
                changed = {(i, j) for i, j in young_cells(sub) if hook_length(sub, i, j) != hook_length(lam, i, j)}
                assert changed == {(i, c[1]) for i in range(c[0])} | {(c[0], j) for j in range(c[1])}
                assert all(hook_length(lam, i, j) - hook_length(sub, i, j) == 1 for i, j in changed)
            # 鉤上漫步：列舉所有漫步的精確機率 = (1/n) H(λ)/H(λ-c) = M14 的乘積式
            probs = hook_walk_probs(lam)
            assert sum(probs.values()) == 1
            assert all(probs[c] == corner_prob_formula(lam, c) == corner_prob_product(lam, c) for c in cs_)
            # 關鍵恆等式與引理（依起點、列集合 A、行集合 B 分組加總）
            if n <= 6:
                for start in young_cells(lam):
                    groups = {}
                    for path, w in hook_walks(lam, start):
                        A, B = walk_rows_cols(path)
                        assert path[-1] == (A[-1], B[-1]) and path[-1] in cs_
                        groups[(path[-1], A, B)] = groups.get((path[-1], A, B), 0) + w
                    for (c, A, B), total in groups.items():
                        assert total == hook_lemma_rhs(lam, c, A, B)
                for c in cs_:
                    for a in range(c[0] + 1):
                        for b in range(c[1] + 1):
                            hook_identity(lam, (a, b), c)
        # SYT(n,n) <-> Dyck，f^(n,n) = C_n，H((n,n)) = (n+1)! n!
        if 1 <= n <= 7:
            tabs = syt((n, n))
            assert len(tabs) == catalan(n) and hook_product((n, n)) == factorial(n + 1) * factorial(n)
            assert sorted(syt_to_dyck(T) for T in tabs) == dyck_words(n)
            assert all(dyck_to_syt(syt_to_dyck(T)) == T for T in tabs)
        print(f"  n={n}: p(n)={len(parts):>3}  合成 {len(comps):>3}  "
              f"整數分割、Young 表、鉤長公式、角落遞迴、鉤上漫步  OK")

    # 場景裡用到的具體數字
    assert integer_partitions(5) == [(5,), (4, 1), (3, 2), (3, 1, 1), (2, 2, 1), (2, 1, 1, 1), (1, 1, 1, 1, 1)]
    assert len(integer_compositions(5)) == 16 and (2, 3) in integer_compositions(5)
    assert syt((3, 2)) == [((1, 2, 3), (4, 5)), ((1, 2, 4), (3, 5)), ((1, 2, 5), (3, 4)),
                           ((1, 3, 4), (2, 5)), ((1, 3, 5), (2, 4))]
    assert hook_lengths((4, 3, 1)) == [(6, 4, 3, 1), (4, 2, 1), (1,)]
    assert hook_product((4, 3, 1)) == 576 and hook_formula((4, 3, 1)) == 70 == len(syt((4, 3, 1)))
    assert hook_product((3, 2)) == 24 and hook_formula((3, 2)) == 5
    assert hook_product((3, 3)) == 144 and hook_formula((3, 3)) == 5
    assert corners((3, 2)) == [(0, 2), (1, 1)]
    assert [remove_corner((3, 2), c) for c in corners((3, 2))] == [(2, 2), (3, 1)]
    assert syt_count((2, 2)) == 2 and syt_count((3, 1)) == 3
    assert hook_product((2, 2)) == 12 and hook_product((3, 1)) == 8
    assert hook_walk_probs((3, 2)) == {(0, 2): Fraction(2, 5), (1, 1): Fraction(3, 5)}
    assert contains_pattern("41532", "231") == (0, 2, 3)
    assert [len(avoiders(n, (2, 3, 1))) for n in (3, 4)] == [5, 14]
    assert not stack_sort((2, 3, 1))[0] and queue_bypass_sort((2, 3, 1))[0]
    assert not queue_bypass_sort((3, 2, 1))[0]
    assert len(perfect_matchings(6)) == 15
    assert matching_is_noncrossing(((0, 1), (2, 3), (4, 5))) and matching_is_nonnesting(((0, 1), (2, 3), (4, 5)))
    assert drunkard_fall_words(1) == ["UDD"] and drunkard_fall_words(2) == ["UDUDD", "UUDDD"]

    # 貫穿的例子（mit_plan.md 第 3 節）
    EX = binary_trees(3)[3]
    assert EX == MIT_EXAMPLE_TREE == ((None, (None, None)), None)
    assert parenthesizations(3)[3] == "(a·(b·c))·d"
    word = tree_to_dyck(EX)
    assert word == "UUDUDD" and dyck_to_ballot(word) == "AABABB" and dyck_to_brackets(word) == "(()())"
    diags = sorted(tree_to_triangulation(EX, 0, 4))
    assert diags == [(0, 3), (1, 3)] and [{a + 1, b + 1} for a, b in diags] == [{1, 4}, {2, 4}]
    assert triangulation_to_tree(5, diags) == EX
    assert [d["tri"] for d in triangulation_dual(EX)] == [(0, 1, 3), (1, 2, 3), (0, 3, 4)]
    assert binary_to_plane(EX) == (((), ()),) == plane_trees(3)[3]
    assert plane_tree_dfs(binary_to_plane(EX)) == "UUDUDD" and plane_tree_walk(binary_to_plane(EX)) == [0, 1, 2, 1, 3, 1, 0]
    assert [(a + 1, b + 1) for a, b in dyck_to_matching(word)] == [(1, 6), (2, 3), (4, 5)]
    assert dyck_to_syt(word) == ((1, 2, 4), (3, 5, 6))
    print("MIT 補充模組：所有計數與一一對應都通過")
