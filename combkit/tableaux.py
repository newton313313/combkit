"""
tableaux.py —— 純 Python 的 Young 表組合演算（不依賴 manim）

內容對應 W. Fulton《Young Tableaux》的 Notation 與第 1 章（§1.1 列插入、§1.2 jeu de taquin），
給場景 F01–F08 使用；第 7 節是第 2–4 章（字與 Knuth 等價、遞增子序列、R-S-K、矩陣球），給 F09–F19 使用。
所有慣例都照書：

資料表示
  分割 λ    : 由大到小的正整數 tuple，例如 (6, 4, 4, 2)；空分割是 ()。
  格子 box  : (i, j) = 第 i 列、第 j 行，0 起算。英式畫法：第 0 列在最上面、第 0 行在最左邊。
  共軛 λ̃    : conjugate(λ)，第 j 項 = Young 圖第 j 行的長度（書 p.2：(6,4,4,2) 的共軛是 (4,4,3,3,1,1)）。
  表 T      : 列的 tuple，每一列是正整數的 tuple，例如 ((1, 2, 2, 3), (2, 3, 5, 5), (4, 4, 6), (5, 6))；
              空表是 ()。（半標準）表 = 每列 → 弱遞增、每行 ↓ 嚴格遞增；標準表 = 1..n 各用一次。
  斜表 S    : 也是列的 tuple，但第 i 列的前 μ_i 格放 None，代表被挖掉的格子（λ/μ 裡的 μ），
              第 i 列的長度 = λ_i。例：書 p.12 的 (5,5,4,3,2)/(4,4,1) 斜表是
                  ((None, None, None, None, 2),
                   (None, None, None, None, 6),
                   (None, 2, 4, 4),
                   (2, 3, 6),
                   (5, 5))
              λ 與 μ 都能從 S 讀回來（skew_shape），所以同一個斜圖形有兩種寫法時
              （書 p.13：(3,1,1)/(2,1) 與 (3,2,1)/(2,2)）內角也分得清楚。可以有整列都是 None 的列，
              但不能有空列。μ = () 的斜表就是一般的表，表示法完全相同，所以所有函數兩種都能吃。
              需要字典 {(i, j): 數字}（只含有數字的格子）時用 skew_cells(S)，反過來用 skew_from_cells。
  多項式    : dict {指數 tuple (a_1, ..., a_m): 係數}，代表 Σ 係數 · x_1^a_1 ... x_m^a_m。

書上的定義（照抄成程式）
  列插入 T ← x（§1.1）：x 進第 0 列；若該列沒有 > x 的數就接在列尾；否則把「最左邊、嚴格大於 x」的數
      擠出來，換它插入下一列，依此類推。回傳新表、擠壓路徑（每列一格，由上到下）、新格子。
  表的字 w(T)：由「最下面一列」讀到「最上面一列」，每列由左到右。
  乘積 T·U：把 w(U) 的字母依序列插入 T。
  滑動（§1.2）：選一個內角（μ 的角）當洞，把洞右邊與下面兩個鄰居中「較小的」滑進洞；
      一樣大時滑「下面」那個；只有一個鄰居就滑那個；沒有鄰居（洞到了 λ 的外角）就把洞拿掉。
  反向滑動：洞往上或往左，跟兩個鄰居中「較大的」交換，一樣大時選「上面」；洞成為內角時停。
  T*U：一個「T 的行數 × U 的列數」的空白矩形當 μ，T 放在矩形下方、U 放在矩形右邊
      （T 在左下、U 在右上，書 p.15 的圖）。書上的 Claim 2、3：Rect 與內角的選法無關，Rect(T*U) = T·U。
  Schur 多項式 s_λ(x_1..x_m) = Σ_T x^T，T 跑遍形狀 λ、數字在 [m] = {1..m} 的表；
      Kostka 數 K_{λμ} = 形狀 λ、內容（各數字出現次數）為 μ 的表的個數。

直接執行本檔會用暴力列舉驗證第 1 章的所有敘述：
    python utils/tableaux.py
"""
from __future__ import annotations

from bisect import bisect_left, bisect_right
from collections import Counter
from itertools import combinations, combinations_with_replacement, permutations
from itertools import product as _cartesian     # 名稱 product 留給表的乘積 T·U
from typing import Callable, Dict, List, Optional, Sequence, Tuple

Box = Tuple[int, int]
Partition = Tuple[int, ...]
Tableau = Tuple[Tuple[int, ...], ...]
SkewTableau = Tuple[Tuple[Optional[int], ...], ...]
Poly = Dict[Tuple[int, ...], int]


def _rows(T) -> tuple:
    """把 list 的 list 之類的輸入整理成 tuple 的 tuple。"""
    return tuple(tuple(row) for row in T)


# ---------------------------------------------------------------------------
# 1. 分割與 Young 圖（Notation）
# ---------------------------------------------------------------------------
def partitions(n: int, max_parts: int | None = None) -> List[Partition]:
    """n 的所有分割（由大到小的 tuple），反字典序；max_parts 限制列數。例：partitions(3) = [(3,), (2, 1), (1, 1, 1)]。"""
    def rec(rem, mx, parts_left):
        if rem == 0:
            yield ()
            return
        if parts_left == 0:
            return
        for first in range(min(rem, mx), 0, -1):
            for rest in rec(rem - first, first, parts_left - 1):
                yield (first,) + rest

    return list(rec(n, n, n if max_parts is None else max_parts))


def is_partition(lam) -> bool:
    """λ 是否為由大到小（弱遞減）的正整數序列。"""
    lam = tuple(lam)
    return all(isinstance(a, int) and a > 0 for a in lam) and all(a >= b for a, b in zip(lam, lam[1:]))


def conjugate(lam) -> Partition:
    """共軛分割 λ̃（沿主對角線翻轉 Young 圖）：λ̃_j = 第 j 行的長度。"""
    lam = tuple(lam)
    return tuple(sum(1 for r in lam if r > j) for j in range(lam[0])) if lam else ()


def young_boxes(lam) -> List[Box]:
    """Young 圖的所有格子 (i, j)，逐列由左到右、由上到下。"""
    return [(i, j) for i, r in enumerate(lam) for j in range(r)]


def corners(lam) -> List[Box]:
    """角（右邊與下面都沒有格子的格子）[(i, λ_i - 1), ...]，由上到下。對 λ 就是外角，對 μ 就是內角。"""
    lam = tuple(lam)
    return [(i, r - 1) for i, r in enumerate(lam) if i == len(lam) - 1 or lam[i + 1] < r]


def addable_boxes(lam) -> List[Box]:
    """加上去之後仍是 Young 圖的格子，由上到下（最後一個是新的一列 (len(λ), 0)）。"""
    lam = tuple(lam)
    return [(i, lam[i] if i < len(lam) else 0) for i in range(len(lam) + 1)
            if i == 0 or (lam[i] if i < len(lam) else 0) < lam[i - 1]]


def contains(lam, mu) -> bool:
    """μ ⊂ λ（μ 的 Young 圖包含在 λ 的裡面），也就是對所有 i 有 μ_i <= λ_i。"""
    lam, mu = tuple(lam), tuple(mu)
    return len(mu) <= len(lam) and all(m <= l for m, l in zip(mu, lam))


def subpartitions(lam) -> List[Partition]:
    """所有 μ ⊂ λ（含 λ 自己與 ()）。"""
    lam = tuple(lam)

    def rec(i, prev):
        if i < len(lam):
            for v in range(min(lam[i], prev), 0, -1):
                for rest in rec(i + 1, v):
                    yield (v,) + rest
        yield ()

    return list(rec(0, lam[0] if lam else 0))


def skew_boxes(lam, mu) -> List[Box]:
    """斜圖形 λ/μ 的格子，逐列由左到右。"""
    mu = tuple(mu) + (0,) * (len(lam) - len(mu))
    return [(i, j) for i, r in enumerate(lam) for j in range(mu[i], r)]


def is_horizontal_strip(lam, mu) -> bool:
    """λ/μ 是水平帶（μ ⊂ λ，且沒有兩格在同一行）<=> λ_{i+1} <= μ_i。"""
    lam = tuple(lam)
    mu = tuple(mu) + (0,) * (len(lam) - len(mu))
    return contains(lam, mu) and all(lam[i + 1] <= mu[i] for i in range(len(lam) - 1))


def is_vertical_strip(lam, mu) -> bool:
    """λ/μ 是垂直帶（μ ⊂ λ，且沒有兩格在同一列）<=> λ_i - μ_i <= 1。"""
    lam = tuple(lam)
    mu = tuple(mu) + (0,) * (len(lam) - len(mu))
    return contains(lam, mu) and all(l - m <= 1 for l, m in zip(lam, mu))


def dominates(lam, mu) -> bool:
    """支配序 λ ⊵ μ：|λ| = |μ| 且每個前綴和 λ_1+...+λ_k >= μ_1+...+μ_k。"""
    lam, mu = tuple(lam), tuple(mu)
    if sum(lam) != sum(mu):
        return False
    k = max(len(lam), len(mu))
    lam, mu = lam + (0,) * (k - len(lam)), mu + (0,) * (k - len(mu))
    a = b = 0
    for x, y in zip(lam, mu):
        a, b = a + x, b + y
        if a < b:
            return False
    return True


# ---------------------------------------------------------------------------
# 2. 表（Notation）
# ---------------------------------------------------------------------------
def shape(T) -> Partition:
    """表（或斜表）的外形 λ = 各列長度。"""
    return tuple(len(row) for row in T)


def num_boxes(S) -> int:
    """有數字的格子數（斜表不算 None）。"""
    return sum(1 for row in S for v in row if v is not None)


def is_tableau(T) -> bool:
    """（半標準）Young 表：外形是分割、數字是正整數、每列弱遞增、每行嚴格遞增。"""
    T = _rows(T)
    if not is_partition(shape(T)):
        return False
    if any(not isinstance(v, int) or isinstance(v, bool) or v < 1 for row in T for v in row):
        return False
    rows_ok = all(a <= b for row in T for a, b in zip(row, row[1:]))
    cols_ok = all(T[i][j] < T[i + 1][j] for i in range(len(T) - 1) for j in range(len(T[i + 1])))
    return rows_ok and cols_ok


def is_standard(T) -> bool:
    """標準表：是表，且數字恰好是 1..n 各一次。"""
    return is_tableau(T) and sorted(v for row in T for v in row) == list(range(1, num_boxes(T) + 1))


def transpose(T) -> tuple:
    """轉置 T^τ（共軛圖上的編號）：T^τ 的 (j, i) = T 的 (i, j)。標準表轉置後仍是標準表，一般的表不一定。"""
    T = _rows(T)
    return tuple(tuple(T[i][j] for i in range(c)) for j, c in enumerate(conjugate(shape(T))))


def content(T, m: int | None = None) -> Tuple[int, ...]:
    """內容（= 單項式 x^T 的指數）：(1 出現幾次, 2 出現幾次, ..., m 出現幾次)。m 預設為最大的數。"""
    cnt = Counter(v for row in T for v in row if v is not None)
    if m is None:
        m = max(cnt) if cnt else 0
    assert all(v <= m for v in cnt), "表裡有大於 m 的數"
    return tuple(cnt[k] for k in range(1, m + 1))


def ssyt(lam, m: int | None = None, content: Sequence[int] | None = None) -> List[Tableau]:
    """
    形狀 λ、數字在 [m] 的所有表，依逐列讀法的字典序。
    給 content=(c_1, ..., c_m) 時只列出數字 k 恰好出現 c_k 次的（此時 m 預設為 len(content)）。
    例：ssyt((2, 1), 2) = [((1, 1), (2,)), ((1, 2), (2,))]。
    """
    lam = tuple(lam)
    need = list(content) if content is not None else None
    if m is None:
        assert need is not None, "要給 m 或 content"
        m = len(need)
    if need is not None and sum(need) != sum(lam):
        return []
    boxes = young_boxes(lam)
    conj = conjugate(lam)
    grid: Dict[Box, int] = {}
    out: List[Tableau] = []

    def rec(k):
        if k == len(boxes):
            out.append(tuple(tuple(grid[(i, j)] for j in range(r)) for i, r in enumerate(lam)))
            return
        i, j = boxes[k]
        lo = max(grid[(i, j - 1)] if j else 1, grid[(i - 1, j)] + 1 if i else 1)
        hi = m - (conj[j] - 1 - i)          # 這一行下面還有 conj[j]-1-i 格要嚴格遞增
        for v in range(lo, hi + 1):
            if need is not None:
                if v > len(need) or need[v - 1] == 0:
                    continue
                need[v - 1] -= 1
            grid[(i, j)] = v
            rec(k + 1)
            if need is not None:
                need[v - 1] += 1
        grid.pop((i, j), None)

    rec(0)
    return out


def standard_tableaux(lam) -> List[Tableau]:
    """形狀 λ 的所有標準表（= 內容為 (1, 1, ..., 1) 的表）。"""
    n = sum(lam)
    return ssyt(lam, content=(1,) * n)


def tableaux_of_size(n: int, m: int) -> List[Tableau]:
    """恰有 n 格、數字在 [m] 的所有表（所有形狀）。"""
    return [T for lam in partitions(n, max_parts=m) for T in ssyt(lam, m)]


# ---------------------------------------------------------------------------
# 3. 列插入（§1.1）
# ---------------------------------------------------------------------------
def bump_steps(T, x: int) -> List[Tuple[Box, int, Optional[int]]]:
    """
    列插入 T ← x 的每一步：[(格子, 放進去的數, 被擠出來的數 或 None), ...]，由上到下每列一筆。
    最後一筆的被擠出數是 None（那個數停在列尾或新的一列），它的格子就是新格子。
    例：bump_steps(((1, 2, 2, 3), (2, 3, 5, 5), (4, 4, 6), (5, 6)), 2)
        = [((0, 3), 2, 3), ((1, 2), 3, 5), ((2, 2), 5, 6), ((3, 2), 6, None)]（書 p.7–8）
    """
    rows = [list(row) for row in T]
    steps: List[Tuple[Box, int, Optional[int]]] = []
    r = 0
    while True:
        if r == len(rows):                  # 擠出最底下：自成新的一列
            steps.append(((r, 0), x, None))
            return steps
        row = rows[r]
        c = bisect_right(row, x)            # 最左邊、嚴格大於 x 的位置
        if c == len(row):                   # 沒有比 x 大的：接在列尾
            steps.append(((r, c), x, None))
            return steps
        y = row[c]
        row[c] = x
        steps.append(((r, c), x, y))
        x, r = y, r + 1


def row_insert(T, x: int) -> Tuple[Tableau, List[Box], Box]:
    """
    列插入 T ← x。回傳 (新表, 擠壓路徑, 新格子)。
    擠壓路徑 = 書上的 bumping route：每列一格、由上到下，是有數被擠出（或最後落下）的格子；
    最後一格就是新格子。
    """
    rows = [list(row) for row in T]
    steps = bump_steps(T, x)
    for (r, c), v, _ in steps:
        if r == len(rows):
            rows.append([])
        if c == len(rows[r]):
            rows[r].append(v)
        else:
            rows[r][c] = v
    route = [box for box, _, _ in steps]
    return _rows(rows), route, route[-1]


def reverse_bump(T, box: Box) -> Tuple[Tableau, int, List[Box]]:
    """
    反向擠壓（書 p.8）：從外角 box 的數 y 開始，在上一列找「最右邊、嚴格小於 y」的數，換上去並把它往上擠，
    一直到擠出第 0 列。回傳 (原來的表, 被擠出的數 x, 路徑)，路徑與 row_insert 一樣由上到下。
    也就是 row_insert(T0, x) == (T, 路徑, box)。
    """
    T = _rows(T)
    assert box in corners(shape(T)), "反向擠壓要從外角開始"
    rows = [list(row) for row in T]
    r, c = box
    y = rows[r].pop()
    route = [box]
    for rr in range(r - 1, -1, -1):
        row = rows[rr]
        k = bisect_left(row, y) - 1         # 最右邊、嚴格小於 y 的位置
        row[k], y = y, row[k]
        route.append((rr, k))
    route.reverse()
    return tuple(tuple(row) for row in rows if row), y, route


def tableau_word(S) -> Tuple[int, ...]:
    """表（或斜表）的字 w(T)：由最下面一列讀到最上面一列，每列由左到右（斜表略過 None）。"""
    return tuple(v for row in reversed(S) for v in row if v is not None)


def insert_word_steps(word, T=()) -> List[Tuple[int, Tableau, List[Box], Box]]:
    """把字母依序列插入 T，每一步記下 (字母, 插入後的表, 擠壓路徑, 新格子)。"""
    T = _rows(T)
    out = []
    for x in word:
        T, route, new = row_insert(T, x)
        out.append((x, T, route, new))
    return out


def insert_word(word, T=()) -> Tableau:
    """((T ← w_1) ← w_2) ← ... ← w_k；T 預設為空表，此時就是 P(w)。"""
    T = _rows(T)
    for x in word:
        T = row_insert(T, x)[0]
    return T


def product(T, U) -> Tableau:
    """乘積 T·U（§1.1）：把 U 的字 w(U)（由下而上、由左到右）依序列插入 T。"""
    return insert_word(tableau_word(U), T)


# ---------------------------------------------------------------------------
# 4. 斜表與 jeu de taquin（§1.2）
# ---------------------------------------------------------------------------
def skew_shape(S) -> Tuple[Partition, Partition]:
    """斜表的形狀 (λ, μ)：λ_i = 第 i 列長度，μ_i = 第 i 列開頭 None 的個數。"""
    lam = tuple(len(row) for row in S)
    mu = []
    for row in S:
        k = 0
        while k < len(row) and row[k] is None:
            k += 1
        mu.append(k)
    while mu and mu[-1] == 0:
        mu.pop()
    return lam, tuple(mu)


def skew_cells(S) -> Dict[Box, int]:
    """斜表裡有數字的格子 {(i, j): 數字}。"""
    return {(i, j): v for i, row in enumerate(S) for j, v in enumerate(row) if v is not None}


def skew_from_cells(lam, mu, cells: Dict[Box, int]) -> SkewTableau:
    """skew_cells 的反函數：給 λ、μ 與 λ/μ 每一格的數字，組回斜表。"""
    mu = tuple(mu) + (0,) * (len(lam) - len(mu))
    return tuple(tuple(None if j < mu[i] else cells[(i, j)] for j in range(r))
                 for i, r in enumerate(lam) if r > 0)


def is_skew_tableau(S) -> bool:
    """斜表：每列是「若干 None + 正整數」，λ、μ 都是分割且 μ ⊂ λ，數字每列弱遞增、每行嚴格遞增。"""
    S = _rows(S)
    if any(len(row) == 0 for row in S):
        return False
    lam, mu = skew_shape(S)
    if not (is_partition(lam) and is_partition(mu) and contains(lam, mu)):
        return False
    mu_full = mu + (0,) * (len(lam) - len(mu))
    for i, row in enumerate(S):
        vals = row[mu_full[i]:]
        if any(not isinstance(v, int) or isinstance(v, bool) or v < 1 for v in vals):
            return False
        if any(a > b for a, b in zip(vals, vals[1:])):
            return False
    cells = skew_cells(S)
    return all(v < cells[(i + 1, j)] for (i, j), v in cells.items() if (i + 1, j) in cells)


def inside_corners(S) -> List[Box]:
    """內角：μ 裡面、下面與右邊都不在 μ 的格子，由上到下。沒有內角 <=> S 已是一般的表。"""
    return corners(skew_shape(S)[1])


def outside_corners(S) -> List[Box]:
    """外角：λ 裡面、下面與右邊都不在 λ 的格子，由上到下（可能同時是內角，見書 p.13）。"""
    return corners(skew_shape(S)[0])


def forward_slide(S, corner: Box) -> Tuple[SkewTableau, List[Box]]:
    """
    從內角 corner 做一次滑動（挖洞）。回傳 (新斜表, 洞的路徑)。
    路徑由 corner 開始、到被拿掉的外角為止；第 k 步是「路徑[k+1] 上的數移到 路徑[k]」，
    移動的數就是原來 S 在 路徑[k+1] 的數（每格最多移一次）。
    規則：洞的右鄰與下鄰取較小者滑進來；一樣大取下面；只有一個就取那個；都沒有就停。
    """
    S = _rows(S)
    lam, mu = skew_shape(S)
    assert corner in corners(mu), "滑動要從內角開始"
    cells = skew_cells(S)
    hole = corner
    path = [hole]
    while True:
        i, j = hole
        right, below = cells.get((i, j + 1)), cells.get((i + 1, j))
        if right is None and below is None:
            break
        nxt = (i + 1, j) if right is None or (below is not None and below <= right) else (i, j + 1)
        cells[hole] = cells.pop(nxt)
        hole = nxt
        path.append(hole)
    assert hole in corners(lam)
    new_lam, new_mu = list(lam), list(mu)
    new_lam[hole[0]] -= 1
    new_mu[corner[0]] -= 1
    return skew_from_cells(new_lam, tuple(v for v in new_mu if v > 0), cells), path


def reverse_slide(S, box: Box) -> Tuple[SkewTableau, List[Box]]:
    """
    反向滑動（書 p.14）：在 λ 外加上格子 box 當洞，洞往上或往左跟較大的鄰居交換（一樣大選上面），
    直到洞成為內角（上面與左邊都沒有數字）為止，洞併入 μ。回傳 (新斜表, 洞的路徑)。
    若 forward_slide(S, c) == (S2, path)，則 reverse_slide(S2, path[-1]) == (S, path[::-1])。
    """
    S = _rows(S)
    lam, mu = skew_shape(S)
    assert box in addable_boxes(lam), "反向滑動要從 λ 外可加上的格子開始"
    cells = skew_cells(S)
    hole = box
    path = [hole]
    while True:
        i, j = hole
        up = cells.get((i - 1, j)) if i > 0 else None
        left = cells.get((i, j - 1)) if j > 0 else None
        if up is None and left is None:
            break
        nxt = (i - 1, j) if left is None or (up is not None and up >= left) else (i, j - 1)
        cells[hole] = cells.pop(nxt)
        hole = nxt
        path.append(hole)
    new_lam = list(lam) + [0]
    new_lam[box[0]] += 1
    new_mu = list(mu) + [0] * (len(new_lam) - len(mu))
    assert new_mu[hole[0]] == hole[1]
    new_mu[hole[0]] += 1
    return skew_from_cells(tuple(v for v in new_lam if v > 0), tuple(v for v in new_mu if v > 0), cells), path


def jdt_steps(S, choose: Callable[[List[Box]], Box] | None = None) -> List[Tuple[Box, List[Box], SkewTableau]]:
    """
    jeu de taquin：一直選內角滑動，直到沒有內角。回傳每一步 (內角, 洞的路徑, 滑動後的斜表)。
    choose(內角清單) 決定選哪個內角；預設選最下面那個（清單最後一個）。最後一步的斜表就是 Rect(S)。
    """
    S = _rows(S)
    steps = []
    while True:
        cs = inside_corners(S)
        if not cs:
            return steps
        c = choose(cs) if choose is not None else cs[-1]
        S, path = forward_slide(S, c)
        steps.append((c, path, S))


def rectify(S, choose: Callable[[List[Box]], Box] | None = None) -> Tableau:
    """Rect(S)：jeu de taquin 做到底得到的一般表（書的 Claim 2：與內角的選法無關）。"""
    steps = jdt_steps(S, choose)
    return steps[-1][2] if steps else _rows(S)


def star(*tabs) -> SkewTableau:
    """
    T*U（§1.2）：μ 是「T 的行數 × U 的列數」的空白矩形，T 放在矩形下面、U 放在矩形右邊（T 左下、U 右上）。
    也可以給更多個：star(T, U, V) 是書上練習 2 的 T*U*V（由左下到右上依序放 T、U、V），
    它的字是 w(T) w(U) w(V)，Rect 就是 T·U·V。
    例：star(((1, 2),), ((3,),)) = ((None, None, 3), (1, 2))。
    """
    tabs = [_rows(T) for T in tabs]
    cols = [len(T[0]) if T else 0 for T in tabs]
    offset = sum(cols)
    out = []
    for T, c in zip(reversed(tabs), reversed(cols)):   # 最後一個放最上面
        offset -= c
        out.extend((None,) * offset + row for row in T)
    return tuple(out)


def skew_tableaux(lam, mu, m: int) -> List[SkewTableau]:
    """形狀 λ/μ、數字在 [m] 的所有斜表，依逐列讀法的字典序。"""
    lam, mu = tuple(lam), tuple(mu)
    assert contains(lam, mu)
    boxes = skew_boxes(lam, mu)
    grid: Dict[Box, int] = {}
    out: List[SkewTableau] = []

    def rec(k):
        if k == len(boxes):
            out.append(skew_from_cells(lam, mu, grid))
            return
        i, j = boxes[k]
        lo = max(grid.get((i, j - 1), 1), grid.get((i - 1, j), 0) + 1)
        for v in range(lo, m + 1):
            grid[(i, j)] = v
            rec(k + 1)
        grid.pop((i, j), None)

    rec(0)
    return out


# ---------------------------------------------------------------------------
# 5. Schur 多項式與 Kostka 數（Notation p.3）
# ---------------------------------------------------------------------------
def schur_poly(lam, m: int) -> Poly:
    """s_λ(x_1, ..., x_m) = Σ_T x^T（T 跑遍形狀 λ、數字在 [m] 的表）。回傳 {指數: 係數}；λ 列數 > m 時是 {}（零）。"""
    out: Poly = {}
    for T in ssyt(lam, m):
        e = content(T, m)
        out[e] = out.get(e, 0) + 1
    return out


def complete_h(n: int, m: int) -> Poly:
    """完全齊次對稱多項式 h_n(x_1..x_m)：所有 n 次單項式，係數都是 1。"""
    out: Poly = {}
    for combo in combinations_with_replacement(range(m), n):
        e = [0] * m
        for k in combo:
            e[k] += 1
        out[tuple(e)] = 1
    return out


def elementary_e(n: int, m: int) -> Poly:
    """基本對稱多項式 e_n(x_1..x_m) = Σ_{i_1 < ... < i_n} x_{i_1} ... x_{i_n}。"""
    return {tuple(1 if k in combo else 0 for k in range(m)): 1 for combo in combinations(range(m), n)}


def kostka(lam, mu) -> int:
    """Kostka 數 K_{λμ}：形狀 λ、內容為 μ（數字 k 出現 μ_k 次）的表的個數。μ 可以是任意（弱）合成。"""
    return len(ssyt(lam, content=tuple(mu)))


# ---------------------------------------------------------------------------
# 6. 書上的例子（場景直接拿來用）
# ---------------------------------------------------------------------------
FULTON_TABLEAU = ((1, 2, 2, 3, 3, 5), (2, 3, 5, 5), (4, 4, 6, 6), (5, 6))        # p.2 左圖，x^T = x1 x2^3 x3^3 x4^2 x5^4 x6^3
FULTON_STANDARD = ((1, 3, 7, 12, 13, 15), (2, 5, 10, 14), (4, 8, 11, 16), (6, 9))  # p.2 右圖
FULTON_INSERT_T = ((1, 2, 2, 3), (2, 3, 5, 5), (4, 4, 6), (5, 6))               # p.7：插入 2；p.12、p.15 的 T
FULTON_PRODUCT_U = ((1, 3), (2,))                                                # p.12、p.15 的 U
FULTON_SKEW = ((None, None, None, None, 2),                                      # p.12–13：(5,5,4,3,2)/(4,4,1)
               (None, None, None, None, 6),
               (None, 2, 4, 4),
               (2, 3, 6),
               (5, 5))


# ---------------------------------------------------------------------------
# 7. 第 2–4 章：字、Knuth 等價、遞增子序列、RSK、矩陣球（場景 F09–F19）
# ---------------------------------------------------------------------------
Word = Tuple[int, ...]


def word_runs(w) -> List[Word]:
    """把字在「後一個字母嚴格變小」的地方切開，得到極大的弱遞增片段（依原來順序）。書 p.17。"""
    runs: List[list] = []
    for a in w:
        if runs and runs[-1][-1] <= a:
            runs[-1].append(a)
        else:
            runs.append([a])
    return [tuple(r) for r in runs]


def word_to_tableau(w) -> Optional[Tableau]:
    """w 若是某個表的字，傳回那個表（片段由下往上當各列），否則 None。書 p.17。"""
    T = tuple(reversed(word_runs(w)))
    return T if is_tableau(T) else None


def skew_from_word(w) -> SkewTableau:
    """字為 w 的一個斜表：每個弱遞增片段當一列，每列放在下一列的整個右上方。書 p.17–18。"""
    runs = word_runs(w)
    starts = [sum(len(r) for r in runs[:t]) for t in range(len(runs))]
    return tuple((None,) * s + r for s, r in reversed(list(zip(starts, runs))))


def column_word(S) -> Word:
    """行字 w_col(S)：由左到右逐行，每行由下往上讀（斜表略過 None）。書 p.27。"""
    S = _rows(S)
    width = max((len(r) for r in S), default=0)
    return tuple(S[i][j] for j in range(width) for i in reversed(range(len(S)))
                 if j < len(S[i]) and S[i][j] is not None)


KNUTH_RULES = ("K'", "K'⁻¹", "K''", "K''⁻¹")


def knuth_moves(w) -> List[Tuple[int, str, Word]]:
    """
    一步基本 Knuth 變換能到的字：[(i, 規則, 新字)]，i = 被改的三個字母的起點（0 起算）。書 p.19：
        K'    y z x → y x z  （x < y ≤ z）      K'⁻¹  反方向
        K''   x z y → z x y  （x ≤ y < z）      K''⁻¹ 反方向
    """
    w = tuple(w)
    out = []
    for i in range(len(w) - 2):
        a, b, c = w[i:i + 3]
        if c < a <= b:
            out.append((i, "K'", w[:i] + (a, c, b) + w[i + 3:]))
        if b < a <= c:
            out.append((i, "K'⁻¹", w[:i] + (a, c, b) + w[i + 3:]))
        if a <= c < b:
            out.append((i, "K''", w[:i] + (b, a, c) + w[i + 3:]))
        if b <= c < a:
            out.append((i, "K''⁻¹", w[:i] + (b, a, c) + w[i + 3:]))
    return out


def knuth_class(w) -> set:
    """與 w Knuth 等價的所有字（廣度優先搜尋；只適合短字）。"""
    w = tuple(w)
    seen, todo = {w}, [w]
    while todo:
        for _, _, v in knuth_moves(todo.pop()):
            if v not in seen:
                seen.add(v)
                todo.append(v)
    return seen


def knuth_path(w1, w2) -> Optional[List[Tuple[int, str, Word]]]:
    """從 w1 到 w2 的最短基本變換序列 [(i, 規則, 變換後的字)]；不等價時傳回 None。"""
    w1, w2 = tuple(w1), tuple(w2)
    prev = {w1: None}
    queue = [w1]
    for cur in queue:
        if cur == w2:
            path = []
            while prev[cur] is not None:
                parent, i, rule = prev[cur]
                path.append((i, rule, cur))
                cur = parent
            return path[::-1]
        for i, rule, v in knuth_moves(cur):
            if v not in prev:
                prev[v] = (cur, i, rule)
                queue.append(v)
    return None


def insertion_knuth_chain(T, x) -> List[Tuple[int, str, Word]]:
    """
    Knuth 對列插入的拆解（書 p.18–19）：從 w(T)·x 一步一步變到 w(T ← x)，每步一個基本變換 [(i, 規則, 新字)]。
    在每一列：x 先用 K' 一格一格往左移過 v（比被擠出的 x' 還右邊的數），x' 再用 K'' 往左移過 u，
    移到下一列（在字裡是它左邊那一段）的尾端，換 x' 繼續。
    """
    T = _rows(T)
    rows = [list(r) for r in T]
    w = list(tableau_word(T)) + [x]
    out = []

    def step(i, rule, new3):
        w[i:i + 3] = new3
        out.append((i, rule, tuple(w)))

    cur, p = x, len(w) - 1                      # p = 正在移動的字母的位置
    for r in range(len(rows)):
        row = rows[r]
        j = bisect_right(row, cur)
        if j == len(row):                       # 接在列尾：字已經是 w(T ← x)
            row.append(cur)
            break
        big, u, v = row[j], row[:j], row[j + 1:]
        start = p - len(row)                    # 這一列在字裡的起點
        for _ in v:                             # K'：cur 往左移過 v
            a, b, c = w[p - 2:p + 1]
            assert c < a <= b
            step(p - 2, "K'", [a, c, b])
            p -= 1
        for _ in u:                             # K''：big 往左移過 u
            a, b, c = w[p - 2:p + 1]
            assert a <= c < b
            step(p - 2, "K''", [b, a, c])
            p -= 1
        row[j] = cur
        cur, p = big, start                     # big 在這一列的最左邊 = 下一列的尾端之後
        assert w[p] == big
    else:
        rows.append([cur])                      # 自成最下面一列：已經在字的最前面
    assert tuple(w) == tableau_word(row_insert(T, x)[0])
    return out


def max_increasing(w, k: int = 1, decreasing: bool = False) -> Tuple[int, List[List[int]]]:
    """
    L(w, k)（書 p.30）：從 w 取出 k 個互不相交的弱遞增子序列，長度總和的最大值，以及一組達到它的子序列
    （每個是位置 list，0 起算）。decreasing=True 改成嚴格遞減子序列（第 3 章 Exercise 1）。
    動態規劃：狀態 = k 個子序列目前的最後一個字母（排序後），所以結果是精確的最大值。
    """
    empty = float("inf") if decreasing else 0
    best = {(empty,) * k: (0, ((),) * k)}
    for pos, a in enumerate(w):
        new = dict(best)
        for state, (cnt, seqs) in best.items():
            for s in range(k):
                if s and state[s] == state[s - 1]:
                    continue
                if (state[s] > a) if decreasing else (state[s] <= a):
                    last, sq = list(state), list(seqs)
                    last[s], sq[s] = a, seqs[s] + (pos,)
                    order = sorted(range(k), key=lambda t: (last[t], len(sq[t])))
                    key = tuple(last[t] for t in order)
                    if key not in new or new[key][0] < cnt + 1:
                        new[key] = (cnt + 1, tuple(sq[t] for t in order))
        best = new
    cnt, seqs = max(best.values(), key=lambda v: v[0])
    return cnt, sorted((list(s) for s in seqs if s), key=lambda s: (-len(s), s))


def greene_numbers(w, decreasing: bool = False) -> Partition:
    """L(w,1), L(w,2)-L(w,1), ...：依 Greene 定理就是 P(w) 的列長（decreasing=True 時是行長）。"""
    out, prev, k = [], 0, 1
    while prev < len(w):
        cur = max_increasing(w, k, decreasing)[0]
        out.append(cur - prev)
        prev, k = cur, k + 1
    return tuple(out)


def remove_largest(w, p: int = 1) -> Word:
    """去掉 p 個最大的字母（相同時右邊的算比較大，先去掉）。書 p.33 Lemma 3。"""
    order = sorted(range(len(w)), key=lambda i: (w[i], i))
    drop = set(order[len(w) - p:]) if p else set()
    return tuple(a for i, a in enumerate(w) if i not in drop)


def tableau_from_greene(w) -> Tableau:
    """
    第 3 章的唯一性證明當演算法（書 p.34）：只用 L(w,k) 決定形狀，一次去掉一個最大字母，
    看形狀少了哪一格，那一格就放那個字母。結果等於 P(w)。
    """
    w = tuple(w)
    if not w:
        return ()
    lam = greene_numbers(w)
    w0 = remove_largest(w)
    T0 = tableau_from_greene(w0)
    lam0 = shape(T0)
    i = next(i for i in range(len(lam)) if i >= len(lam0) or lam[i] > lam0[i])
    big = max(w)
    rows = [list(r) for r in T0] + ([[]] if i == len(lam0) else [])
    rows[i].append(big)
    return tuple(tuple(r) for r in rows)


def is_lex_array(top, bottom) -> bool:
    """雙列陣列是字典序：上列弱遞增，上列相同時下列弱遞增（書 p.39 (1)(2)）。"""
    pairs_ = list(zip(top, bottom))
    return len(top) == len(bottom) and all(a <= b for a, b in zip(pairs_, pairs_[1:]))


def rsk_steps(top, bottom=None) -> List[Tuple[int, int, Tableau, Tableau, Box]]:
    """
    R-S-K（書 p.36–39）：陣列 (top 在上、bottom 在下，須是字典序) → 每一步 (u, v, P_k, Q_k, 新格子)：
    v 列插入 P，新格子在 Q 裡放 u。bottom 省略時 top 當成字 w，Q 是記錄表（放 1..r）。
    """
    if bottom is None:
        top, bottom = tuple(range(1, len(top) + 1)), tuple(top)
    assert is_lex_array(top, bottom), "陣列要先排成字典序"
    P: Tableau = ()
    Q: list = []
    out = []
    for u, v in zip(top, bottom):
        P, _, (i, j) = row_insert(P, v)
        if i == len(Q):
            Q.append([])
        Q[i].append(u)
        assert len(Q[i]) == j + 1
        out.append((u, v, P, tuple(tuple(r) for r in Q), (i, j)))
    return out


def rsk(top, bottom=None) -> Tuple[Tableau, Tableau]:
    """R-S-K 對應的 (P, Q)；用法同 rsk_steps。"""
    steps = rsk_steps(top, bottom)
    return (steps[-1][2], steps[-1][3]) if steps else ((), ())


def inverse_rsk(P, Q) -> Tuple[Word, Word]:
    """由同形狀的 (P, Q) 還原字典序陣列 (top, bottom)：每次取 Q 最大的數（相同取最右邊）那一格，對 P 反向擠出。"""
    P, Q = _rows(P), [list(r) for r in _rows(Q)]
    assert shape(P) == shape(tuple(map(tuple, Q)))
    top, bottom = [], []
    while Q:
        big = max(v for r in Q for v in r)
        i, j = max(((i, j) for i, r in enumerate(Q) for j, v in enumerate(r) if v == big), key=lambda b: b[1])
        P, x, _ = reverse_bump(P, (i, j))
        Q[i].pop()
        if not Q[i]:
            Q.pop()
        top.append(big)
        bottom.append(x)
    return tuple(reversed(top)), tuple(reversed(bottom))


def matrix_to_array(A) -> Tuple[Word, Word]:
    """矩陣 A（a(i,j) = (i 在上、j 在下) 出現的次數，i、j 從 1 起算）→ 字典序陣列。書 p.41。"""
    top, bottom = [], []
    for i, row in enumerate(A, 1):
        for j, a in enumerate(row, 1):
            top += [i] * a
            bottom += [j] * a
    return tuple(top), tuple(bottom)


def array_to_matrix(top, bottom, m: int | None = None, n: int | None = None) -> Tuple[Tuple[int, ...], ...]:
    m = m or max(top, default=0)
    n = n or max(bottom, default=0)
    A = [[0] * n for _ in range(m)]
    for u, v in zip(top, bottom):
        A[u - 1][v - 1] += 1
    return tuple(map(tuple, A))


def rsk_matrix(A) -> Tuple[Tableau, Tableau]:
    """矩陣的 R-S-K：P 的數來自行號、Q 的數來自列號。"""
    return rsk(*matrix_to_array(A))


def transpose_matrix(A):
    return tuple(zip(*A)) if A else ()


def ball_numbering(A) -> Dict[Box, List[int]]:
    """
    矩陣球的編號（書 p.42）：(i, j) 放 a(i,j) 個球（0 起算）；球的號碼 = 比它西北方所有球的號碼都大的最小整數
    （同一格的球由西北往東南排，號碼連續）。回傳 {(i, j): [這一格各球的號碼]}。
    """
    nums: Dict[Box, List[int]] = {}
    for i, row in enumerate(A):
        for j, a in enumerate(row):
            if a:
                nw = [x for (r, c), xs in nums.items() if r <= i and c <= j for x in xs]
                base = max(nw, default=0) + 1
                nums[(i, j)] = list(range(base, base + a))
    return nums


def derived_matrix(A):
    """
    導出矩陣 A♭（書 p.43–44）：同號的 ℓ 個球由西南排到東北，除了最後一個，每個球的右邊、
    下一個球的正下方放一個新球（即 (這個球的列, 下一個球的行)）。
    """
    m, n = len(A), len(A[0]) if A else 0
    B = [[0] * n for _ in range(m)]
    by_num: Dict[int, List[Box]] = {}
    for pos, xs in ball_numbering(A).items():
        for x in xs:
            by_num.setdefault(x, []).append(pos)
    for balls in by_num.values():
        balls.sort(key=lambda b: b[1])          # 西南 → 東北
        for (r1, c1), (r2, c2) in zip(balls, balls[1:]):
            assert r1 > r2 and c1 < c2
            B[r1][c2] += 1
    return tuple(map(tuple, B))


def matrix_ball(A) -> Tuple[Tableau, Tableau, List[dict]]:
    """
    矩陣球構造（書 §4.2）：每一輪把球編號，P 的這一列 = 各號碼最左邊的行號、Q 的這一列 = 各號碼最上面的列號
    （1 起算），再換成導出矩陣，直到沒有球。回傳 (P, Q, 各輪 {"A", "balls", "P_row", "Q_row"})。
    """
    stages, P, Q = [], [], []
    while any(any(r) for r in A):
        balls = ball_numbering(A)
        top = max(x for xs in balls.values() for x in xs)
        p_row = tuple(min(c for (r, c), xs in balls.items() if k in xs) + 1 for k in range(1, top + 1))
        q_row = tuple(min(r for (r, c), xs in balls.items() if k in xs) + 1 for k in range(1, top + 1))
        stages.append({"A": tuple(map(tuple, A)), "balls": balls, "P_row": p_row, "Q_row": q_row})
        P.append(p_row)
        Q.append(q_row)
        A = derived_matrix(A)
    return tuple(P), tuple(Q), stages


def pieri_shapes(lam, p: int, vertical: bool = False) -> List[Partition]:
    """λ 加 p 格、任兩格不同行（水平條）或不同列（vertical=True，垂直條）得到的所有 μ。書 p.24 (4)(5)。"""
    lam = tuple(lam)
    out = []
    for mu in partitions(sum(lam) + p):
        if contains(mu, lam) and (is_vertical_strip(mu, lam) if vertical else is_horizontal_strip(mu, lam)):
            out.append(mu)
    return out


def poly_mul(f: Poly, g: Poly) -> Poly:
    out: Poly = {}
    for a, x in f.items():
        for b, y in g.items():
            e = tuple(s + t for s, t in zip(a, b))
            out[e] = out.get(e, 0) + x * y
    return {e: c for e, c in out.items() if c}


def poly_sum(polys) -> Poly:
    out: Poly = {}
    for f in polys:
        for e, c in f.items():
            out[e] = out.get(e, 0) + c
    return {e: c for e, c in out.items() if c}


def lex_greater(lam, mu) -> bool:
    """字典序 λ > μ：第一個不同的位置 λ_i > μ_i。書 p.26。"""
    k = max(len(lam), len(mu))
    a, b = tuple(lam) + (0,) * (k - len(lam)), tuple(mu) + (0,) * (k - len(mu))
    return a > b


def kostka_matrix(n: int) -> Tuple[List[Partition], List[List[int]]]:
    """n 的所有分割（字典序由小到大，同書 p.26）與 Kostka 矩陣 K[λ][μ]（下三角、對角線為 1）。"""
    parts = partitions(n)
    parts.sort(key=lambda p: p + (0,) * (n - len(p)))
    return parts, [[kostka(lam, mu) for mu in parts] for lam in parts]


def hook_content_count(lam, m: int) -> int:
    """Stanley 的鉤–內容公式（書 p.55）：數字在 [m] 的 λ 形表的個數 d_λ(m) = Π (m + j − i) / h(i,j)。"""
    from fractions import Fraction
    lam = tuple(lam)
    lt = conjugate(lam)
    val = Fraction(1)
    for i, j in young_boxes(lam):
        val *= Fraction(m + j - i, (lam[i] - j) + (lt[j] - i) - 1)
    assert val.denominator == 1
    return int(val)


def count_involutions(n: int) -> int:
    """S_n 裡對合（σ² = 1）的個數 = Σ_k n! / ((n−2k)! 2^k k!)（書 p.52 Exercise 6）。"""
    from math import factorial
    return sum(factorial(n) // (factorial(n - 2 * k) * 2 ** k * factorial(k)) for k in range(n // 2 + 1))


# 第 2–4 章書上的例子
CH2_WORD = (5, 6, 4, 4, 6, 2, 3, 5, 5, 1, 2, 2, 3)                # p.17：= w(FULTON_INSERT_T)
CH2_SLIDE_BEFORE = ((1, None, 4), (2, 3, 5))                       # p.20：u x y／v z 的數字版（u,v,x,y,z = 1,2,3,4,5）
CH2_SLIDE_AFTER = ((1, 3, 4), (2, None, 5))                        #   字 v x z u y = 23514 → v z u x y = 25134
CH2_SLIDE_CHAIN = ((2, 3, 5, 1, 4), (2, 3, 1, 5, 4), (2, 1, 3, 5, 4), (2, 1, 5, 3, 4), (2, 5, 1, 3, 4))
CH3_WORD = (1, 3, 4, 2, 3, 4, 1, 2, 2, 3, 3, 2)                    # p.30
CH3_WORD2 = (3, 4, 4, 2, 3, 3, 1, 1, 2, 2, 2, 3)                   # p.31：= w(P(CH3_WORD))
CH3_TABLEAU = ((1, 1, 2, 2, 2, 3), (2, 3, 3), (3, 4, 4))           # p.31
CH4_WORD = (5, 4, 8, 2, 3, 4, 1, 7, 5, 3, 1)                       # p.37
CH4_ARRAY = ((1, 1, 1, 2, 2, 3, 3, 3, 3), (1, 2, 2, 1, 2, 1, 1, 1, 2))   # p.40 Exercise 1
CH4_MATRIX = ((1, 2), (1, 1), (3, 1))                              # p.41、p.43
CH4_P = ((1, 1, 1, 1, 1, 2), (2, 2, 2))                            # p.40、p.43–44
CH4_Q = ((1, 1, 1, 2, 3, 3), (2, 3, 3))
CH4_HOOK_SHAPE = (6, 5, 5, 3)                                      # p.53–56：f^λ = 6,651,216、d_λ(5) = 3,360


# ---------------------------------------------------------------------------
# 8. 第 5–6 章：Littlewood–Richardson 規則、對稱多項式（場景 F20–F27）
# ---------------------------------------------------------------------------
def _skew_content(S) -> Tuple[int, ...]:
    c = Counter(v for r in S for v in r if v is not None)
    return tuple(c[i] for i in range(1, max(c, default=0) + 1))


def factorizations(V, lam, mu) -> List[Tuple[Tableau, Tableau]]:
    """書 p.61 的 𝒯(λ, μ, V)：所有 [T, U]，T 形狀 λ、U 形狀 μ、T·U = V（暴力列舉）。"""
    V = _rows(V)
    m = max(v for r in V for v in r)
    cv = Counter(v for r in V for v in r)
    out = []
    for T in ssyt(lam, m):
        ct = Counter(v for r in T for v in r)
        if any(ct[k] > cv[k] for k in ct):
            continue
        rest = tuple(cv[k] - ct[k] for k in range(1, m + 1))
        for U in ssyt(mu, content=rest):
            if product(T, U) == V:
                out.append((T, U))
    return out


def rectifying_skews(nu, lam, U0) -> List[SkewTableau]:
    """書 p.60 的 𝒮(ν/λ, U0)：形狀 ν/λ、rectification 等於 U0 的斜表。"""
    U0 = _rows(U0)
    m = max(v for r in U0 for v in r)
    want = _skew_content(U0)
    return [S for S in skew_tableaux(nu, lam, m) if _skew_content(S) == want and rectify(S) == U0]


def is_reverse_lattice(w) -> bool:
    """反格字（Yamanouchi 字，書 p.63）：從字尾往前讀到任何位置，1 的個數 ≥ 2 的個數 ≥ 3 的個數 ≥ …。"""
    cnt = Counter()
    for a in reversed(w):
        cnt[a] += 1
        if a > 1 and cnt[a] > cnt[a - 1]:
            return False
    return True


def lattice_tableau(mu) -> Tableau:
    """U(μ)：第 i 列全是 i 的表（LR 斜表 rectify 之後的樣子，書 p.65）。"""
    return tuple((i + 1,) * r for i, r in enumerate(mu))


def lr_tableaux(nu, lam, mu=None) -> List[SkewTableau]:
    """ν/λ 上的 Littlewood–Richardson 斜表：字是反格字（μ 給定時只取內容為 μ 的）。"""
    out = []
    for S in skew_tableaux(nu, lam, len(nu)):
        if is_reverse_lattice(tableau_word(S)) and (mu is None or _skew_content(S) == tuple(mu)):
            out.append(S)
    return out


def lr_coefficient(lam, mu, nu) -> int:
    """Littlewood–Richardson 數 c^ν_{λμ} = ν/λ 上內容為 μ 的 LR 斜表個數（書 p.64 命題 3）。"""
    lam, mu, nu = tuple(lam), tuple(mu), tuple(nu)
    if sum(nu) != sum(lam) + sum(mu) or not contains(nu, lam) or not contains(nu, mu):
        return 0
    return len(lr_tableaux(nu, lam, mu))


def lr_product(lam, mu) -> Dict[Partition, int]:
    """s_λ s_μ = Σ c^ν_{λμ} s_ν（書 p.66 (4)）：{ν: c}。"""
    out = {}
    for nu in partitions(sum(lam) + sum(mu)):
        c = lr_coefficient(lam, mu, nu)
        if c:
            out[nu] = c
    return out


def skew_decomposition(nu, lam) -> Dict[Partition, int]:
    """s_{ν/λ} = Σ_μ c^ν_{λμ} s_μ（書 p.65）：{μ: c}。"""
    cnt = Counter(_skew_content(S) for S in lr_tableaux(nu, lam))
    return dict(sorted(cnt.items(), key=lambda kv: kv[0], reverse=True))


def reverse_numbering(nu, lam) -> Dict[Box, int]:
    """ν/λ 的反向編號（書 p.68）：由上到下逐列，每列由右往左編 1, 2, …。"""
    out, k = {}, 0
    lam = tuple(lam) + (0,) * (len(nu) - len(lam))
    for i, r in enumerate(nu):
        for j in range(r - 1, lam[i] - 1, -1):
            k += 1
            out[(i, j)] = k
    return out


def remmel_whitney(nu, lam, mu) -> List[Tableau]:
    """Remmel–Whitney（書 p.68–69 命題 4）：μ 上滿足 (i)(ii) 的標準表；個數 = c^ν_{λμ}。"""
    num = reverse_numbering(nu, lam)
    pos = {k: b for b, k in num.items()}
    pairs_row = [(num[(i, j + 1)], num[(i, j)]) for (i, j) in num if (i, j + 1) in num]     # (k-1, k) 同一列
    pairs_col = [(num[(i - 1, j)], num[(i, j)]) for (i, j) in num if (i - 1, j) in num]     # j 在 k 正上方
    out = []
    for U in standard_tableaux(mu):
        where = {v: (i, j) for i, r in enumerate(U) for j, v in enumerate(r)}
        ok = all(where[k][0] <= where[a][0] and where[k][1] > where[a][1] for a, k in pairs_row)
        ok = ok and all(where[k][0] > where[a][0] and where[k][1] <= where[a][1] for a, k in pairs_col)
        if ok:
            out.append(U)
    assert pos
    return out


# ---- 第 6 章：對稱多項式（m 個變數；多項式 = {指數 tuple: 係數}）
def poly_scale(f: Poly, c) -> Poly:
    return {e: c * v for e, v in f.items() if c * v}


def poly_prod(polys, m: int) -> Poly:
    out: Poly = {(0,) * m: 1}
    for f in polys:
        out = poly_mul(out, f)
    return out


def monomial_sym(lam, m: int) -> Poly:
    """m_λ：x^λ 的所有不同排列之和（λ 的列數 > m 時是 0）。"""
    lam = tuple(lam)
    if len(lam) > m:
        return {}
    return {e: 1 for e in set(permutations(lam + (0,) * (m - len(lam))))}


def power_sum(r: int, m: int) -> Poly:
    return {tuple(r if k == i else 0 for k in range(m)): 1 for i in range(m)}


def h_lambda(lam, m):
    return poly_prod([complete_h(p, m) for p in lam], m)


def e_lambda(lam, m):
    return poly_prod([elementary_e(p, m) for p in lam], m)


def p_lambda(lam, m):
    return poly_prod([power_sum(p, m) for p in lam], m)


def z_lambda(lam) -> int:
    """z(λ) = Π r^{m_r} m_r!（書 p.74 (3)）；n!/z(λ) 是 S_n 中循環型為 λ 的排列個數。"""
    from math import factorial
    out = 1
    for r, k in Counter(lam).items():
        out *= r ** k * factorial(k)
    return out


def schur_expand(f: Poly, m: int) -> Dict[Partition, int]:
    """把對稱多項式寫成 Σ c_λ s_λ（每次扣掉字典序最大的單項式；需要 m ≥ 次數才完整）。"""
    f = dict(f)
    out = {}
    while f:
        top = max(f)
        lam = tuple(a for a in top if a)
        assert list(top) == sorted(top, reverse=True), "不是對稱多項式"
        c = f[top]
        out[lam] = c
        f = poly_sum([f, poly_scale(schur_poly(lam, m), -c)])
    return out


def _sign(p) -> int:
    s, seen = 1, set()
    for i in range(len(p)):
        if i not in seen:
            j, L = i, 0
            while j not in seen:
                seen.add(j)
                j, L = p[j], L + 1
            s *= (-1) ** (L - 1)
    return s


def poly_det(M, m: int) -> Poly:
    """多項式矩陣的行列式（Leibniz 展開；只適合小矩陣）。"""
    k = len(M)
    terms = []
    for p in permutations(range(k)):
        terms.append(poly_scale(poly_prod([M[i][p[i]] for i in range(k)], m), _sign(p)))
    return poly_sum(terms)


def jacobi_trudi(lam, m: int, dual: bool = False) -> Poly:
    """(5) s_λ = det(h_{λ_i+j-i})；dual=True 時是 (6) s_λ = det(e_{λ̃_i+j-i})（書 p.75）。"""
    lam = conjugate(lam) if dual else tuple(lam)
    base = elementary_e if dual else complete_h
    k = len(lam)

    def entry(p):
        return {} if p < 0 else ({(0,) * m: 1} if p == 0 else base(p, m))
    return poly_det([[entry(lam[i] + j - i) for j in range(k)] for i in range(k)], m)


def alternant(exps, m: int) -> Poly:
    """a(ℓ_1..ℓ_m) = det(x_j^{ℓ_i})；書 p.75 (7)：s_λ = a(λ+δ)/a(δ)，δ = (m-1, …, 1, 0)。"""
    out: Poly = {}
    for p in permutations(range(m)):
        e = [0] * m
        for i in range(m):
            e[p[i]] = exps[i]
        out[tuple(e)] = out.get(tuple(e), 0) + _sign(p)
    return {e: c for e, c in out.items() if c}


def principal_specialization(lam, m: int) -> Dict[int, int]:
    """s_λ(1, q, …, q^{m-1}) 當成 q 的多項式 {次數: 係數}（書 p.76 練習 5）。"""
    out: Dict[int, int] = {}
    for e, c in schur_poly(lam, m).items():
        d = sum(i * a for i, a in enumerate(e))
        out[d] = out.get(d, 0) + c
    return dict(sorted(out.items()))


def character_table(n: int) -> Tuple[List[Partition], List[List[int]]]:
    """χ^λ_μ（書 p.78 (11)：p_μ = Σ_λ χ^λ_μ s_λ）。回傳 (n 的分割, 表[λ][μ])，都用 partitions(n) 的順序。"""
    parts = partitions(n)
    table = [[0] * len(parts) for _ in parts]
    for b, mu in enumerate(parts):
        exp = schur_expand(p_lambda(mu, n), n)
        for a, lam in enumerate(parts):
            table[a][b] = exp.get(lam, 0)
    return parts, table


def schur_inner(f: Poly, g: Poly, m: int) -> int:
    """對稱函數的內積：以 s_λ 為正交規範基底（書 p.78；m 要 ≥ 次數）。"""
    a, b = schur_expand(f, m), schur_expand(g, m)
    return sum(c * b.get(lam, 0) for lam, c in a.items())


def omega(f: Poly, m: int) -> Poly:
    """對合 ω：s_λ ↦ s_λ̃（書 p.78）。"""
    return poly_sum(poly_scale(schur_poly(conjugate(lam), m), c) for lam, c in schur_expand(f, m).items())


# 第 5–6 章書上的例子
CH5_RL_GOOD = (2, 1, 3, 2, 1, 2, 1)                               # p.63：反格字
CH5_RL_BAD = (1, 2, 3, 2, 1, 2, 1)                                #   不是：最後 6 個字母裡 2 比 1 多
CH5_NU, CH5_LAM = (5, 4, 3, 2), (3, 3, 1)                         # p.63–65、p.68–70
CH5_DECOMP = {(5, 2): 1, (5, 1, 1): 1, (4, 1, 1, 1): 1, (2, 2, 2, 1): 1, (4, 3): 2, (3, 2, 1, 1): 2,
              (3, 3, 1): 2, (3, 2, 2): 2, (4, 2, 1): 3}
CH5_S21_SQUARED = {(4, 2): 1, (4, 1, 1): 1, (3, 3): 1, (3, 2, 1): 2, (3, 1, 1, 1): 1, (2, 2, 2): 1,
                   (2, 2, 1, 1): 1}                                # p.66–67
CH5_ZERO = ((2, 2), (3, 2), (3, 3, 3))                             # p.62–63：ν 穿出 λ*μ，所以 c = 0


# ---------------------------------------------------------------------------
# 自我檢查
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    import random
    from math import factorial

    try:
        from combkit.counting import hook_product, syt_count
    except ImportError:                     # 直接執行 python utils/tableaux.py 時，utils/ 本身在 sys.path 上
        from counting import hook_product, syt_count

    def pairs(seq):
        return list(zip(seq, seq[1:]))

    # ----------------------------------------------- Notation：分割、表、書上的例子
    print("Notation（分割、共軛、表、標準表）")
    assert partitions(4) == [(4,), (3, 1), (2, 2), (2, 1, 1), (1, 1, 1, 1)]
    assert [len(partitions(n)) for n in range(9)] == [1, 1, 2, 3, 5, 7, 11, 15, 22]
    assert conjugate((6, 4, 4, 2)) == (4, 4, 3, 3, 1, 1)
    assert all(conjugate(conjugate(lam)) == lam for n in range(9) for lam in partitions(n))
    assert is_tableau(FULTON_TABLEAU) and shape(FULTON_TABLEAU) == (6, 4, 4, 2)
    assert content(FULTON_TABLEAU) == (1, 3, 3, 2, 4, 3)
    assert is_standard(FULTON_STANDARD) and shape(FULTON_STANDARD) == (6, 4, 4, 2)
    assert is_standard(transpose(FULTON_STANDARD)) and not is_tableau(transpose(FULTON_TABLEAU))
    assert not is_tableau(((1, 2), (1,))) and not is_tableau(((2, 1),)) and not is_tableau(((1,), (2, 3)))
    # ssyt 的列舉 = 暴力（所有填法再篩選）
    for n in range(0, 5):
        for lam in partitions(n):
            for m in range(1, 4):
                brute = set()
                for vals in _cartesian(range(1, m + 1), repeat=n):
                    it = iter(vals)
                    T = tuple(tuple(next(it) for _ in range(r)) for r in lam)
                    if is_tableau(T):
                        brute.add(T)
                got = ssyt(lam, m)
                assert len(got) == len(set(got)) and set(got) == brute
    # 標準表的個數 = 鉤長公式 n!/H(λ)（H 取自 counting.py）= 角落遞迴
    n_syt = 0
    for n in range(0, 9):
        for lam in partitions(n):
            tabs = standard_tableaux(lam)
            assert all(is_standard(T) and shape(T) == lam for T in tabs) and len(set(tabs)) == len(tabs)
            assert len(tabs) == factorial(n) // hook_product(lam) == syt_count(lam)
            assert all(is_standard(transpose(T)) for T in tabs)
            n_syt += len(tabs)
    print(f"  標準表 {n_syt} 個（n <= 8），個數全部等於鉤長公式  OK")

    # 小的表全部列出來（之後的暴力檢查都用這些）
    small6 = [T for n in range(7) for T in tableaux_of_size(n, 4)]      # <= 6 格，數字 1..4
    small4 = [T for n in range(5) for T in tableaux_of_size(n, 3)]      # <= 4 格，數字 1..3
    small2 = [T for n in range(3) for T in tableaux_of_size(n, 3)]      # <= 2 格，數字 1..3
    pool = [T for n in range(8) for T in tableaux_of_size(n, 5)]        # <= 7 格，數字 1..5（隨機抽樣用）
    assert all(is_tableau(T) for T in small6 + pool)
    rng = random.Random(1987)

    # ----------------------------------------------- §1.1 列插入
    print("\n§1.1 列插入")
    T2, route, new = row_insert(FULTON_INSERT_T, 2)                     # 書 p.7–8 的例子
    assert T2 == ((1, 2, 2, 2), (2, 3, 3, 5), (4, 4, 5), (5, 6, 6))
    assert route == [(0, 3), (1, 2), (2, 2), (3, 2)] and new == (3, 2)
    assert bump_steps(FULTON_INSERT_T, 2) == [((0, 3), 2, 3), ((1, 2), 3, 5), ((2, 2), 5, 6), ((3, 2), 6, None)]
    assert reverse_bump(T2, (3, 2)) == (FULTON_INSERT_T, 2, route)
    assert row_insert((), 5) == (((5,),), [(0, 0)], (0, 0))
    assert row_insert(((1, 2),), 1) == (((1, 1), (2,)), [(0, 1), (1, 0)], (1, 0))
    count = 0
    for T in small6:
        for x in range(1, 6):
            U, route, new = row_insert(T, x)
            assert is_tableau(U) and num_boxes(U) == num_boxes(T) + 1
            assert sorted(v for row in U for v in row) == sorted([x] + [v for row in T for v in row])
            assert new == route[-1] and new in corners(shape(U))
            assert set(young_boxes(shape(U))) == set(young_boxes(shape(T))) | {new} != set(young_boxes(shape(T)))
            assert [b[0] for b in route] == list(range(len(route)))         # 每列一格、由上到下
            assert all(b2[1] <= b1[1] for b1, b2 in pairs(route))           # 路徑只會往左或垂直往下
            assert reverse_bump(U, new) == (T, x, route)                    # 由新格子還原
            count += 1
    print(f"  {count} 組 (T, x)：結果都是表、路徑每列一格且不往右、可由新格子反向擠壓還原  OK")

    # Row Bumping Lemma
    count = 0
    for T in small6:
        for x, x2 in _cartesian(range(1, 6), repeat=2):
            T1, R, B = row_insert(T, x)
            _, R2, B2 = row_insert(T1, x2)
            if x <= x2:     # (1) R 嚴格在 R' 左邊；B 嚴格在 B' 左邊、弱在 B' 下面
                assert len(R2) <= len(R) and all(R[k][1] < R2[k][1] for k in range(len(R2)))
                assert B[1] < B2[1] and B[0] >= B2[0]
            else:           # (2) R' 弱在 R 左邊；B' 弱在 B 左邊、嚴格在 B 下面
                assert len(R2) > len(R) and all(R2[k][1] <= R[k][1] for k in range(len(R)))
                assert B2[1] <= B[1] and B2[0] > B[0]
            count += 1
    print(f"  Row Bumping Lemma：{count} 組 (T, x, x') 全部符合  OK")

    # 命題：弱遞增序列 -> 新格子成水平帶（由左到右）；嚴格遞減 -> 垂直帶（由上到下）
    count = 0
    for T in small4:
        for p in range(1, 4):
            for seq in combinations_with_replacement(range(1, 5), p):          # x_1 <= ... <= x_p
                steps = insert_word_steps(seq, T)
                boxes = [s[3] for s in steps]
                assert all(b2[1] > b1[1] and b2[0] <= b1[0] for b1, b2 in pairs(boxes))
                assert is_horizontal_strip(shape(steps[-1][1]), shape(T))
                count += 1
            for seq in combinations(range(4, 0, -1), p):                        # x_1 > ... > x_p
                steps = insert_word_steps(seq, T)
                boxes = [s[3] for s in steps]
                assert all(b2[0] > b1[0] and b2[1] <= b1[1] for b1, b2 in pairs(boxes))
                assert is_vertical_strip(shape(steps[-1][1]), shape(T))
                count += 1
    # 反方向：U 與水平帶（垂直帶）μ/λ -> 唯一的 T 與 x_1 <= ... <= x_p（x_1 > ... > x_p）
    for U in small6:
        mu = shape(U)
        for lam in subpartitions(mu):
            for horizontal in (True, False):
                if not (is_horizontal_strip(mu, lam) if horizontal else is_vertical_strip(mu, lam)):
                    continue
                strip = skew_boxes(mu, lam)
                order = sorted(strip, key=lambda b: b[1]) if horizontal else sorted(strip)
                V, xs = U, []
                for b in reversed(order):           # 水平帶從最右邊、垂直帶從最下面開始反向擠壓
                    V, x, _ = reverse_bump(V, b)
                    xs.append(x)
                xs.reverse()
                assert shape(V) == lam and insert_word(xs, V) == U
                assert all((a <= b) if horizontal else (a > b) for a, b in pairs(xs))
                count += 1
    print(f"  水平帶／垂直帶命題（正反兩個方向）：{count} 個情形  OK")

    # ----------------------------------------------- 字與乘積
    print("\n字與乘積 T·U")
    assert tableau_word(FULTON_INSERT_T) == (5, 6, 4, 4, 6, 2, 3, 5, 5, 1, 2, 2, 3)
    assert tableau_word(FULTON_PRODUCT_U) == (2, 1, 3)
    steps = insert_word_steps((2, 1, 3), FULTON_INSERT_T)                # 書 p.12 的三步
    assert steps[0][1] == ((1, 2, 2, 2), (2, 3, 3, 5), (4, 4, 5), (5, 6, 6))
    assert steps[1][1] == ((1, 1, 2, 2), (2, 2, 3, 5), (3, 4, 5), (4, 6, 6), (5,))
    TU = ((1, 1, 2, 2, 3), (2, 2, 3, 5), (3, 4, 5), (4, 6, 6), (5,))
    assert steps[2][1] == TU == product(FULTON_INSERT_T, FULTON_PRODUCT_U)
    assert all(insert_word(tableau_word(T)) == T for T in small6)          # P(w(T)) = T
    assert all(product((), T) == T == product(T, ()) for T in small6)      # 空表是單位元
    count = 0
    for T in small4:
        for U in small4:
            P = product(T, U)
            assert is_tableau(P) and num_boxes(P) == num_boxes(T) + num_boxes(U)
            assert content(P, 5) == tuple(a + b for a, b in zip(content(T, 5), content(U, 5)))
            if num_boxes(U) == 1:
                assert P == row_insert(T, U[0][0])[0]
            count += 1
    assoc = 0
    for T in small2:
        for U in small2:
            TU_ = product(T, U)
            for V in small2:
                assert product(TU_, V) == product(T, product(U, V))
                assoc += 1
    for _ in range(1000):
        T, U, V = rng.choice(pool), rng.choice(pool), rng.choice(pool)
        assert product(product(T, U), V) == product(T, product(U, V))
        assoc += 1
    print(f"  {count} 組 (T, U) 的乘積都是表、格數與內容相加；結合律 {assoc} 組  OK")

    # ----------------------------------------------- §1.2 斜表與 jeu de taquin
    print("\n§1.2 斜表與 jeu de taquin")
    S = FULTON_SKEW                                                       # 書 p.12–13 的例子
    assert is_skew_tableau(S) and skew_shape(S) == ((5, 5, 4, 3, 2), (4, 4, 1))
    assert inside_corners(S) == [(1, 3), (2, 0)]
    assert outside_corners(S) == [(1, 4), (2, 3), (3, 2), (4, 1)]
    assert tableau_word(S) == (5, 5, 2, 3, 6, 2, 4, 4, 6, 2)
    S2, path = forward_slide(S, (2, 0))
    assert path == [(2, 0), (3, 0), (3, 1), (4, 1)]
    assert S2 == ((None, None, None, None, 2), (None, None, None, None, 6), (2, 2, 4, 4), (3, 5, 6), (5,))
    assert reverse_slide(S2, (4, 1)) == (S, path[::-1])
    # p.13：同一個斜圖形的兩種寫法，內角不同；(3,2,1)/(2,2) 的內角 (1,1) 同時是外角
    A = ((None, None, 1), (None,), (2,))                                  # (3,1,1)/(2,1)
    B = ((None, None, 1), (None, None), (2,))                             # (3,2,1)/(2,2)
    assert skew_cells(A) == skew_cells(B) and skew_shape(B) == ((3, 2, 1), (2, 2))
    assert inside_corners(A) == [(0, 1), (1, 0)] and inside_corners(B) == [(1, 1)]
    assert (1, 1) in outside_corners(B) and forward_slide(B, (1, 1)) == (A, [(1, 1)])
    assert rectify(A) == rectify(B) == ((1,), (2,))
    # T*U 的形狀（p.15 的圖與練習 1）
    TstarU = star(FULTON_INSERT_T, FULTON_PRODUCT_U)
    assert TstarU == ((None, None, None, None, 1, 3), (None, None, None, None, 2),
                      (1, 2, 2, 3), (2, 3, 5, 5), (4, 4, 6), (5, 6))
    assert skew_shape(TstarU) == ((6, 5, 4, 4, 3, 2), (4, 4))
    assert rectify(TstarU) == TU                                          # 練習 1 的答案 = p.12 的乘積
    assert star((), FULTON_PRODUCT_U) == FULTON_PRODUCT_U and star(FULTON_PRODUCT_U, ()) == FULTON_PRODUCT_U
    assert star(((1, 2),), ((3,),)) == ((None, None, 3), (1, 2))

    # 所有小斜表：每次滑動都得到斜表、可反向還原；Rect 與內角的選法無關（Claim 2）；Rect(S) = P(w(S))
    memo: Dict[tuple, frozenset] = {}

    def all_rects(S):
        if S not in memo:
            cs = inside_corners(S)
            memo[S] = frozenset([S]) if not cs else frozenset().union(*(all_rects(forward_slide(S, c)[0]) for c in cs))
        return memo[S]

    n_skew = n_slides = 0
    for n in range(1, 9):
        for lam in partitions(n):
            for mu in subpartitions(lam):
                if not mu or not 1 <= n - sum(mu) <= 5:
                    continue
                for S in skew_tableaux(lam, mu, 3):
                    assert is_skew_tableau(S) and skew_shape(S) == (lam, mu)
                    for c in inside_corners(S):
                        S2, path = forward_slide(S, c)
                        assert is_skew_tableau(S2) and num_boxes(S2) == num_boxes(S)
                        assert all((b[0] - a[0], b[1] - a[1]) in ((0, 1), (1, 0)) for a, b in pairs(path))
                        assert path[-1] in outside_corners(S)
                        assert reverse_slide(S2, path[-1]) == (S, path[::-1])
                        n_slides += 1
                    rects = all_rects(S)
                    R = rectify(S)
                    assert len(rects) == 1 and R in rects and is_tableau(R)
                    assert R == rectify(S, choose=lambda cs: cs[0])
                    assert R == insert_word(tableau_word(S))
                    n_skew += 1
    print(f"  {n_skew} 個小斜表、{n_slides} 次滑動：都是斜表、可反向滑回；所有內角順序的 Rect 都相同  OK")

    # Rect(T*U) = T·U（Claim 3），以及練習 2：Rect(T*U*V) = T·U·V
    count = 0
    for T in small4:
        for U in small4:
            S = star(T, U)
            assert is_skew_tableau(S) and tableau_word(S) == tableau_word(T) + tableau_word(U)
            assert rectify(S) == product(T, U) == rectify(S, choose=lambda cs: cs[0])
            count += 1
    for _ in range(1000):
        T, U, V = rng.choice(pool), rng.choice(pool), rng.choice(pool)
        assert rectify(star(T, U)) == product(T, U)
        assert rectify(star(T, U, V)) == product(product(T, U), V) == rectify(star(T, U, V), choose=rng.choice)
        count += 1
    print(f"  Rect(T*U) = T·U：{count} 組（含 T*U*V）  OK")

    # ----------------------------------------------- Schur 多項式與 Kostka 數
    print("\nSchur 多項式與 Kostka 數")
    s = schur_poly((2, 1), 3)
    assert s == {(2, 1, 0): 1, (2, 0, 1): 1, (1, 2, 0): 1, (0, 2, 1): 1, (1, 0, 2): 1, (0, 1, 2): 1, (1, 1, 1): 2}
    assert schur_poly((1, 1, 1, 1), 3) == {} and schur_poly((), 2) == {(0, 0): 1}
    count = 0
    for n in range(0, 6):
        for m in range(1, 6 if n <= 4 else 5):
            for lam in partitions(n):
                s = schur_poly(lam, m)
                for perm in permutations(range(m)):                   # 對稱
                    assert {tuple(e[p] for p in perm): c for e, c in s.items()} == s
                for e, c in s.items():                                # 係數 = Kostka 數
                    assert sum(e) == n and c == kostka(lam, tuple(sorted((a for a in e if a), reverse=True)))
                count += 1
            assert schur_poly((n,) if n else (), m) == complete_h(n, m)
            assert schur_poly((1,) * n, m) == elementary_e(n, m)
    for n in range(1, 7):
        for lam in partitions(n):
            assert kostka(lam, lam) == 1
            assert kostka(lam, (1,) * n) == syt_count(lam)
            for mu in partitions(n):
                assert (kostka(lam, mu) > 0) == dominates(lam, mu)
    assert kostka((3, 2), (2, 2, 1)) == 2 and kostka((3, 2), (1, 1, 1, 1, 1)) == 5
    print(f"  {count} 個 s_λ(x_1..x_m)：對稱、s_(n) = h_n、s_(1^n) = e_n、係數 = Kostka 數；"
          f"K_λλ = 1、K_λμ > 0 <=> λ 支配 μ  OK")

    print("Fulton 第 1 章：插入、乘積、jeu de taquin 全部驗證通過")

    # ================================================= 第 2–4 章
    def bf_increasing(w, k, decreasing=False):
        best = 0
        for lab in _cartesian(range(k + 1), repeat=len(w)):
            ok = True
            for s in range(1, k + 1):
                seq = [a for a, l in zip(w, lab) if l == s]
                ok = ok and all((a > b) if decreasing else (a <= b) for a, b in pairs(seq))
            if ok:
                best = max(best, sum(1 for l in lab if l))
        return best

    print("第 2 章（字、Knuth 等價、主定理、Pieri、Kostka、行字）")
    T = FULTON_INSERT_T
    assert tableau_word(T) == CH2_WORD and word_to_tableau(CH2_WORD) == T
    assert word_runs(CH2_WORD) == [(5, 6), (4, 4, 6), (2, 3, 5, 5), (1, 2, 2, 3)]
    assert word_to_tableau((1, 2, 3, 1)) is None and word_to_tableau((3, 4, 1, 4)) is None
    S = skew_from_word(CH3_WORD)
    assert is_skew_tableau(S) and tableau_word(S) == CH3_WORD
    chain = insertion_knuth_chain(T, 2)
    assert chain[-1][2] == tableau_word(row_insert(T, 2)[0]) == (5, 6, 6, 4, 4, 5, 2, 3, 3, 5, 1, 2, 2, 2)
    for _ in range(300):
        T1 = insert_word([random.randint(1, 5) for _ in range(random.randint(0, 9))])
        x = random.randint(1, 5)
        prev = tableau_word(T1) + (x,)
        for i, rule, wd in insertion_knuth_chain(T1, x):
            assert (i, rule, wd) in knuth_moves(prev)
            prev = wd
        assert prev == tableau_word(row_insert(T1, x)[0])
    for w in _cartesian(range(1, 4), repeat=4):          # 主定理：每個等價類恰有一個表的字
        cls = knuth_class(w)
        assert [v for v in cls if word_to_tableau(v) is not None] == [tableau_word(insert_word(w))]
        assert all(insert_word(v) == insert_word(w) for v in cls)
    for lam in [(2, 1), (2, 2), (3, 1), (3, 2), (2, 1, 1)]:   # 第 4 章 Exercise 16：等價類大小 = f^λ
        for T1 in standard_tableaux(lam)[:2] + ssyt(lam, 3)[:2]:
            assert len(knuth_class(tableau_word(T1))) == len(standard_tableaux(lam))
    assert CH2_SLIDE_CHAIN[0] == tableau_word(CH2_SLIDE_BEFORE) and CH2_SLIDE_CHAIN[-1] == tableau_word(CH2_SLIDE_AFTER)
    for a, b in pairs(CH2_SLIDE_CHAIN):
        assert b in [v for _, _, v in knuth_moves(a)]
    for S1 in skew_tableaux((3, 3, 2), (2, 1), 3)[:40]:     # Proposition 2：滑動保持 Knuth 等價
        for _, _, S2 in jdt_steps(S1):
            assert knuth_path(tableau_word(S1), tableau_word(S2)) is not None
    assert rectify(FULTON_SKEW) == insert_word(tableau_word(FULTON_SKEW))          # Corollary 1
    TU = product(FULTON_INSERT_T, FULTON_PRODUCT_U)                                 # Corollary 2：三種乘積
    assert TU == rectify(star(FULTON_INSERT_T, FULTON_PRODUCT_U))
    assert TU == insert_word(tableau_word(FULTON_INSERT_T) + tableau_word(FULTON_PRODUCT_U))
    assert pieri_shapes((2, 1), 2) == [(4, 1), (3, 2), (3, 1, 1), (2, 2, 1)]
    assert pieri_shapes((2, 1), 2, vertical=True) == [(3, 2), (3, 1, 1), (2, 2, 1), (2, 1, 1, 1)]
    for lam in [(1,), (2, 1), (2, 2), (3, 1)]:
        for p in (1, 2):
            assert poly_mul(schur_poly(lam, 3), complete_h(p, 3)) == poly_sum(schur_poly(mu, 3) for mu in pieri_shapes(lam, p))
            assert poly_mul(schur_poly(lam, 3), elementary_e(p, 3)) == poly_sum(
                schur_poly(mu, 3) for mu in pieri_shapes(lam, p, vertical=True))
    parts, K = kostka_matrix(4)
    assert parts == [(1, 1, 1, 1), (2, 1, 1), (2, 2), (3, 1), (4,)]
    for n in (4, 5, 6):
        parts, K = kostka_matrix(n)
        for a, lam in enumerate(parts):
            assert K[a][a] == 1 and all(K[a][b] == 0 for b in range(a + 1, len(parts)))
            for b, mu in enumerate(parts):
                assert (K[a][b] != 0) == dominates(lam, mu)          # Exercise 2
    assert dominates((3, 1), (2, 2)) and not contains((3, 1), (2, 2))
    assert lex_greater((4, 1, 1), (3, 3)) and not dominates((4, 1, 1), (3, 3)) and not dominates((3, 3), (4, 1, 1))
    assert not dominates((2, 2, 2), (3, 1, 1, 1)) and not dominates((3, 1, 1, 1), (2, 2, 2))
    for mu in [(2, 1), (2, 1, 1), (2, 2), (3, 1)]:             # (8) h_μ = Σ K_λμ s_λ
        h = complete_h(mu[0], 3)
        for part in mu[1:]:
            h = poly_mul(h, complete_h(part, 3))
        assert h == poly_sum(f for lam in partitions(sum(mu)) for f in [schur_poly(lam, 3)] * kostka(lam, mu))
    for lam in [(3, 2, 1), (4, 2)]:                           # K_λμ 與 μ 的順序無關
        for mu in set(permutations((3, 2, 1))):
            assert kostka(lam, mu) == kostka(lam, (3, 2, 1))
    for T1 in [FULTON_TABLEAU, FULTON_INSERT_T, CH3_TABLEAU]:   # (10) 行字 ≡ 列字；Lemma 1 切開
        assert insert_word(column_word(T1)) == T1
        for c in range(1, len(T1[0])):
            left = tuple(r[:c] for r in T1 if r[:c])
            right = tuple(r[c:] for r in T1 if r[c:])
            assert insert_word(tableau_word(left) + tableau_word(right)) == T1 == rectify(star(left, right))

    print("第 3 章（遞增子序列 L(w,k)、唯一性的證明）")
    assert [max_increasing(CH3_WORD, k)[0] for k in (1, 2, 3, 4)] == [6, 9, 12, 12]
    assert [max_increasing(CH3_WORD2, k)[0] for k in (1, 2, 3)] == [6, 9, 12]
    assert insert_word(CH3_WORD) == CH3_TABLEAU and tableau_word(CH3_TABLEAU) == CH3_WORD2
    for k in (1, 2, 3):
        cnt, seqs = max_increasing(CH3_WORD, k)
        assert len(seqs) <= k and sum(map(len, seqs)) == cnt and len({p for s in seqs for p in s}) == cnt
        assert all(CH3_WORD[a] <= CH3_WORD[b] for s in seqs for a, b in pairs(s))
    assert remove_largest((1, 3, 4, 2, 3, 1, 3), 3) == (1, 3, 2, 1)              # 書 p.33 的例子
    assert tableau_from_greene(CH3_WORD) == CH3_TABLEAU
    for _ in range(150):
        w = tuple(random.randint(1, 3) for _ in range(random.randint(1, 7)))
        for k in (1, 2):
            assert max_increasing(w, k)[0] == bf_increasing(w, k)
            assert max_increasing(w, k, decreasing=True)[0] == bf_increasing(w, k, decreasing=True)
        P1 = insert_word(w)
        assert greene_numbers(w) == shape(P1) and greene_numbers(w, decreasing=True) == conjugate(shape(P1))
        assert tableau_from_greene(w) == P1
        for _, _, v in knuth_moves(w):                       # Lemma 2、Lemma 3
            assert greene_numbers(v) == greene_numbers(w)
            for p in range(len(w) + 1):
                assert insert_word(remove_largest(v, p)) == insert_word(remove_largest(w, p))
    for w in permutations(range(1, 6)):                       # Erdős–Szekeres：長度 5 > 2²
        assert max_increasing(w)[0] > 2 or max_increasing(w, decreasing=True)[0] > 2

    print("第 4 章（R-S-K、對稱定理、矩陣球、應用）")
    P1, Q1 = rsk(CH4_WORD)
    assert is_tableau(P1) and is_standard(Q1) and shape(P1) == shape(Q1) == shape(insert_word(CH4_WORD))
    assert inverse_rsk(P1, Q1) == (tuple(range(1, 12)), CH4_WORD)
    assert matrix_to_array(CH4_MATRIX) == CH4_ARRAY and array_to_matrix(*CH4_ARRAY) == CH4_MATRIX
    assert rsk(*CH4_ARRAY) == (CH4_P, CH4_Q) == rsk_matrix(CH4_MATRIX)
    assert rsk_matrix(transpose_matrix(CH4_MATRIX)) == (CH4_Q, CH4_P)               # Exercise 2
    Pb, Qb, stages = matrix_ball(CH4_MATRIX)
    assert (Pb, Qb) == (CH4_P, CH4_Q) and len(stages) == 2
    assert (stages[0]["P_row"], stages[0]["Q_row"]) == ((1, 1, 1, 1, 1, 2), (1, 1, 1, 2, 3, 3))
    assert (stages[1]["P_row"], stages[1]["Q_row"]) == ((2, 2, 2), (2, 3, 3))
    for _ in range(300):
        m, n = random.randint(1, 4), random.randint(1, 4)
        A = tuple(tuple(random.choice((0, 0, 1, 1, 2)) for _ in range(n)) for _ in range(m))
        if not any(any(r) for r in A):
            continue
        P1, Q1 = rsk_matrix(A)
        assert matrix_ball(A)[:2] == (P1, Q1) and rsk_matrix(transpose_matrix(A)) == (Q1, P1)
        assert array_to_matrix(*inverse_rsk(P1, Q1), m, n) == A
        k = min(m, n)
        S = tuple(tuple(A[i][j] + A[j][i] if i != j else A[i][i] for j in range(k)) for i in range(k))
        if any(any(r) for r in S):                               # 對稱矩陣 ↔ P = Q；Exercise 4
            P1, Q1 = rsk_matrix(S)
            assert P1 == Q1
            tr = sum(S[i][i] for i in range(k)) + sum(derived_matrix(S)[i][i] for i in range(k))
            assert tr == len(P1[0])
            assert sum(S[i][i] for i in range(k)) == sum(c % 2 for c in conjugate(shape(P1)))
    for w in permutations(range(1, 6)):                          # 排列：P(w⁻¹) = Q(w)；Exercise 17
        inv = tuple(sorted(range(1, 6), key=lambda i: w[i - 1]))
        P1, Q1 = rsk(w)
        assert rsk(inv) == (Q1, P1)
        row_of = {v: i for i, r in enumerate(Q1) for v in r}
        assert all((w[i] > w[i + 1]) == (row_of[i + 2] > row_of[i + 1]) for i in range(4))
    for n in range(1, 7):                                        # (4)、(6)
        assert sum(len(standard_tableaux(l)) ** 2 for l in partitions(n)) == factorial(n)
        assert sum(len(standard_tableaux(l)) for l in partitions(n)) == count_involutions(n)
    for n, m in [(3, 2), (4, 3)]:                                # (5)
        assert sum(len(ssyt(l, m)) * len(standard_tableaux(l)) for l in partitions(n)) == m ** n
    assert sum(len(ssyt(l, 2)) ** 2 for l in partitions(3)) == 20     # (7)：2×2 矩陣、總和 3 的個數 = C(6,3)
    for lam, m in [((2, 1), 3), ((3, 2), 3), ((2, 2, 1), 4)]:
        assert hook_content_count(lam, m) == len(ssyt(lam, m))
    assert hook_content_count(CH4_HOOK_SHAPE, 5) == 3360
    assert factorial(19) // hook_product(CH4_HOOK_SHAPE) == 6651216 == syt_count(CH4_HOOK_SHAPE)
    print("Fulton 第 2–4 章：Knuth 等價、L(w,k)、R-S-K、矩陣球全部驗證通過")

    # ================================================= 第 5–6 章
    from math import prod as _prod

    def idet(M):
        return sum(_sign(p) * _prod(M[i][p[i]] for i in range(len(M))) for p in permutations(range(len(M))))

    print("第 5 章（Littlewood–Richardson 規則）")
    assert is_reverse_lattice(CH5_RL_GOOD) and not is_reverse_lattice(CH5_RL_BAD)
    assert len(lr_tableaux(CH5_NU, CH5_LAM)) == 15 == sum(CH5_DECOMP.values())
    assert skew_decomposition(CH5_NU, CH5_LAM) == CH5_DECOMP
    for mu, c in CH5_DECOMP.items():
        assert len(remmel_whitney(CH5_NU, CH5_LAM, mu)) == c                       # 命題 4
        assert all(rectify(S) == lattice_tableau(mu) for S in lr_tableaux(CH5_NU, CH5_LAM, mu))   # Lemma 1
    assert lr_product((2, 1), (2, 1)) == CH5_S21_SQUARED and lr_coefficient(*CH5_ZERO) == 0
    for lam, mu in [((2, 1), (2, 1)), ((2, 1), (1, 1)), ((2,), (2, 1)), ((1, 1), (2, 1)), ((2, 2), (1,))]:
        prod_ = lr_product(lam, mu)
        n_ = sum(lam) + sum(mu)
        assert schur_expand(poly_mul(schur_poly(lam, n_), schur_poly(mu, n_)), n_) == prod_
        assert prod_ == lr_product(mu, lam)                                           # c^ν_λμ = c^ν_μλ
        for nu, c in prod_.items():
            assert lr_coefficient(conjugate(lam), conjugate(mu), conjugate(nu)) == c   # 共軛
            for V in [standard_tableaux(nu)[0], standard_tableaux(nu)[-1], ssyt(nu, len(nu) + 1)[0]]:
                assert len(factorizations(V, lam, mu)) == c                         # 與 V 無關（推論 1）
            for U0 in [lattice_tableau(mu), standard_tableaux(mu)[-1]]:
                assert len(rectifying_skews(nu, lam, U0)) == c                       # 與 U0 無關
        lm = tuple(a + b for a, b in zip(lam + (0,) * 3, mu + (0,) * 3) if a + b)
        assert lr_coefficient(lam, mu, lm) == 1                                       # 練習 2
    for n_ in range(1, 6):                                                          # 練習 8：卡塔蘭數
        words = [w for w in set(permutations((1,) * n_ + (2,) * n_)) if is_reverse_lattice(w)]
        assert len(words) == factorial(2 * n_) // (factorial(n_ + 1) * factorial(n_))

    print("第 6 章（對稱多項式的基底、Jacobi–Trudi、對稱函數環）")
    m = 4
    parts4 = partitions(4)
    for fam in (monomial_sym, h_lambda, e_lambda):                                  # 命題 1：整係數基底
        M = [[schur_expand(fam(l, m), m).get(mu, 0) for mu in parts4] for l in parts4]
        assert abs(idet(M)) == 1
    assert idet([[schur_expand(p_lambda(l, m), m).get(mu, 0) for mu in parts4] for l in parts4]) != 0
    for k in range(1, 5):
        assert poly_sum(poly_scale(poly_mul(elementary_e(i, m) if i else {(0,) * m: 1},
                                            complete_h(k - i, m) if k - i else {(0,) * m: 1}), (-1) ** i)
                        for i in range(k + 1)) == {}
        assert poly_scale(elementary_e(k, m), k) == poly_sum(                       # Newton 恆等式（練習 1）
            poly_scale(poly_mul(elementary_e(k - i, m) if k - i else {(0,) * m: 1}, power_sum(i, m)), (-1) ** (i - 1))
            for i in range(1, k + 1))
        assert poly_scale(complete_h(k, m), k) == poly_sum(
            poly_mul(complete_h(k - i, m) if k - i else {(0,) * m: 1}, power_sum(i, m)) for i in range(1, k + 1))
        assert poly_scale(complete_h(k, m), factorial(k)) == poly_sum(               # Lemma 1：h_n = Σ p_λ / z(λ)
            poly_scale(p_lambda(l, m), factorial(k) // z_lambda(l)) for l in partitions(k))
        assert sum(factorial(k) // z_lambda(l) for l in partitions(k)) == factorial(k)
    for mu in parts4:
        assert schur_expand(h_lambda(mu, m), m) == {l: kostka(l, mu) for l in parts4 if kostka(l, mu)}
    for lam in [(2, 1), (3, 1), (2, 2), (2, 1, 1), (3, 2, 1), (4,), (1, 1, 1)]:      # (5)(6)(7)
        assert jacobi_trudi(lam, m) == schur_poly(lam, m) == jacobi_trudi(lam, m, dual=True)
    delta = (2, 1, 0)
    for lam in [(1,), (2, 1), (3, 1), (2, 2), (3, 2, 1), (2, 1, 1)]:
        lp = tuple(lam) + (0,) * (3 - len(lam))
        assert poly_mul(schur_poly(lam, 3), alternant(delta, 3)) == alternant(tuple(a + d for a, d in zip(lp, delta)), 3)
    for lam, mm in [((2, 1), 3), ((3, 1), 3), ((2, 2), 4), ((3, 2, 1), 4), (CH4_HOOK_SHAPE, 5)]:   # 練習 5、6
        ps = principal_specialization(lam, mm)
        assert sum(ps.values()) == hook_content_count(lam, mm)
        lp = tuple(lam) + (0,) * (mm - len(lam))
        lhs, rhs = {0: 1}, {sum(i * a for i, a in enumerate(lp)): 1}

        def umul(a, b):
            out_ = Counter()
            for i, x in a.items():
                for j, y in b.items():
                    out_[i + j] += x * y
            return {d: c for d, c in out_.items() if c}
        lhs = ps
        for i in range(mm):
            for j in range(i + 1, mm):
                lhs = umul(lhs, {0: 1, j - i: -1})
                rhs = umul(rhs, {0: 1, lp[i] - lp[j] + j - i: -1})
        assert lhs == rhs
    parts3, X = character_table(3)
    assert parts3 == [(3,), (2, 1), (1, 1, 1)] and X == [[1, 1, 1], [-1, 0, 2], [1, -1, 1]]
    for n_ in (3, 4, 5):                                                             # 特徵值的正交關係
        parts, X = character_table(n_)
        for a in range(len(parts)):
            for b in range(len(parts)):
                assert sum(X[k][a] * X[k][b] for k in range(len(parts))) == (z_lambda(parts[a]) if a == b else 0)
    for l in parts4:                                                                # 命題 3、推論 1（ω）
        for mu in parts4:
            assert schur_inner(h_lambda(l, m), monomial_sym(mu, m), m) == (1 if l == mu else 0)
            assert schur_inner(p_lambda(l, m), p_lambda(mu, m), m) == (z_lambda(l) if l == mu else 0)
        assert omega(h_lambda(l, m), m) == e_lambda(l, m)
        assert omega(p_lambda(l, m), m) == poly_scale(p_lambda(l, m), (-1) ** (sum(l) - len(l)))
    print("Fulton 第 5–6 章：Littlewood–Richardson 規則、對稱函數全部驗證通過")
