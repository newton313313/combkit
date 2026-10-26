"""
fulton_appendix.py —— Fulton《Young Tableaux》附錄 A（組合變形）與附錄 B（代數簇的拓撲）的計算（純 Python，場景 F56–F65）

表用 tableaux.py 的寫法（列的 tuple，格子 (列, 行) 0 起算）；字是 tuple。對偶字母 A* 用 [N] 上的 a ↦ N+1−a 表示。

附錄 A
    dual_word(w, N)            w* = 反過來讀並把每個字母換成對偶（§A.1）
    evacuation(T, N)           對偶表 T*（evacuation）；evacuation_steps(T, N) 逐步記錄（每步拿掉左上角、滑動、把 x* 放進空出的格）
    column_insert(T, x)        行插入（column bumping，§A.2）：回傳 (新表, 被擠的格子路徑, 新格子)
    column_insert_word(w)      依序行插入 w 的字母（從空表開始）
    lex_sort(u, v)、antilex_sort(u, v)    兩列陣列排成字典序／反字典序（上列遞減、相同時下列遞增）
    matrix_to_pairs(A)         矩陣 → 成對的 (u, v)（a_ij 個 (i, j)，1 起算）
    rsk_array(u, v)            R-S-K（先排成字典序）
    burge(A)                   Burge 對應（§A.4.1，作法 (1b)）：反字典序陣列，行插入 v_r, …, v_1、依序放 u_r, …, u_1 → (R, S)
    corners(A)                 四個角（§A.4.2 定理 1）：{'nw': (P,Q), 'se': (P*,Q*), 'sw': (R,S*), 'ne': (R*,S)}，由轉換過的陣列重算
    dual_rsk(A)                0-1 矩陣的對應（§A.4.3，Knuth）：依序行插入 v_1..v_r、在新格子的共軛位置放 u → {P, Q}，形狀互為共軛
    zero_one_matrices(rows, cols)   列和、行和給定的所有 0-1 矩陣（暴力）；gale_ryser(lam, mu) 存在 ⇔ μ ⊴ λ~（練習 20）
    compact_skew(col_lengths)  給定各行長度的「緊湊」斜形（相鄰兩行：左長則頂端對齊、右長則底端對齊）
    frank_tableau(T, col_lengths)   形狀為 compact_skew(col_lengths)、rectify 成 T 的唯一斜表（§A.5，暴力搜尋）
    keys(T)                    左右鑰匙 (K_-(T), K_+(T))（§A.5）
附錄 B
    betti_from_cells(counts)   胞腔個數（依維數）→ Betti 數 b_0, b_1, …（奇數為 0）與 Euler 示性數（Lemma 6）
    chern_projective(n)        c(T P^n) = (1+h)^{n+1} 截到 h^n 的係數；c_n = n+1 = χ(P^n)
    plucker_degree(r, m)       Gr(r, m) 在 Plücker 嵌入下的次數 = σ_1^{r(m-r)} = 矩形的標準表個數

直接執行：python utils/fulton_appendix.py（在 codes/ 內）跑自我測試（含書上的例子）。
"""
from __future__ import annotations

from itertools import combinations, permutations, product
from math import comb

try:
    from combkit.tableaux import (conjugate, dominates, kostka, partitions, rectify, row_insert, rsk, shape,
                                standard_tableaux, transpose)
except ImportError:                                   # 在 utils/ 裡直接執行時
    from tableaux import (conjugate, dominates, kostka, partitions, rectify, row_insert, rsk, shape,  # type: ignore
                          standard_tableaux, transpose)


# ---------------------------------------------------------------- §A.1 對偶與 evacuation
def dual_word(w, N):
    return tuple(N + 1 - x for x in reversed(w))


def _dual_tab(T, N):
    return tuple(tuple(N + 1 - x for x in r) for r in T)


def evacuation_steps(T, N):
    """[(ΔT, 空出的格子, 放進去的 x*), …]：每步拿掉左上角 x、把剩下的斜表滑成一般表，空出的格子放 N+1−x。"""
    T = tuple(tuple(r) for r in T)
    lam = shape(T)
    out_cells = {}
    steps = []
    cur = T
    while cur:
        x = cur[0][0]
        skew = ((None,) + cur[0][1:],) + cur[1:]
        nxt = rectify(skew) if sum(map(len, cur)) > 1 else ()
        old = {(i, j) for i, r in enumerate(cur) for j in range(len(r))}
        new = {(i, j) for i, r in enumerate(nxt) for j in range(len(r))}
        (b,) = old - new
        out_cells[b] = N + 1 - x
        steps.append((nxt, b, N + 1 - x))
        cur = nxt
    assert len(out_cells) == sum(lam)
    return steps


def evacuation(T, N):
    cells = {b: v for _, b, v in evacuation_steps(T, N)}
    lam = shape(T)
    return tuple(tuple(cells[(i, j)] for j in range(lam[i])) for i in range(len(lam)))


# ---------------------------------------------------------------- §A.2 行插入
def _strict_row_insert(T, x):
    """列嚴格遞增的陣列裡的列插入：每列擠掉第一個 ≥ x 的數。"""
    T = [list(r) for r in T]
    path = []
    for i, r in enumerate(T):
        j = next((k for k, y in enumerate(r) if y >= x), None)
        if j is None:
            r.append(x)
            path.append((i, len(r) - 1))
            return tuple(map(tuple, T)), path, (i, len(r) - 1)
        path.append((i, j))
        r[j], x = x, r[j]
    T.append([x])
    path.append((len(T) - 1, 0))
    return tuple(map(tuple, T)), path, (len(T) - 1, 0)


def column_insert(T, x):
    """行插入 x → T：在第一行擠掉最上面（最小）的 ≥ x 的數，被擠的數去下一行；沒有就放在該行底部。"""
    Tt = transpose(T) if T else ()
    Tt2, path, nb = _strict_row_insert(Tt, x)
    return transpose(Tt2), [(j, i) for i, j in path], (nb[1], nb[0])


def column_insert_word(w, T=()):
    for x in w:
        T = column_insert(T, x)[0]
    return T


# ---------------------------------------------------------------- §A.4 陣列與矩陣
def matrix_to_pairs(A):
    return [(i + 1, j + 1) for i, row in enumerate(A) for j, a in enumerate(row) for _ in range(a)]


def lex_sort(u, v):
    pr = sorted(zip(u, v))
    return tuple(p[0] for p in pr), tuple(p[1] for p in pr)


def antilex_sort(u, v):
    pr = sorted(zip(u, v), key=lambda p: (-p[0], p[1]))
    return tuple(p[0] for p in pr), tuple(p[1] for p in pr)


def rsk_array(u, v):
    u, v = lex_sort(u, v)
    return rsk(u, v)


def _place(Q, box, val):
    Q = [list(r) for r in Q]
    i, j = box
    if i == len(Q):
        Q.append([])
    assert len(Q[i]) == j
    Q[i].append(val)
    return tuple(map(tuple, Q))


def burge(A):
    """作法 (1b)：反字典序陣列 (u_1..u_r; v_1..v_r)，依序行插入 v_r, v_{r-1}, …, v_1，並依序把 u_r, …, u_1 放進新格子。"""
    u, v = antilex_sort(*zip(*matrix_to_pairs(A)))
    R, S = (), ()
    for k in range(len(v) - 1, -1, -1):
        R, _, nb = column_insert(R, v[k])
        S = _place(S, nb, u[k])
    return R, S


def corners(A):
    """§A.4.2 定理 1 的四個（弱序）配對，由轉換過的陣列用 R-S-K 重算；自我測試核對它們等於 (P,Q)、(P*,Q*)、(R,S*)、(R*,S)。"""
    pairs = matrix_to_pairs(A)
    m, n = len(A), len(A[0])
    u, v = zip(*pairs)
    du = tuple(m + 1 - x for x in u)
    dv = tuple(n + 1 - x for x in v)
    return {"nw": rsk_array(u, v), "se": rsk_array(du, dv), "sw": rsk_array(du, v), "ne": rsk_array(u, dv)}


def dual_rsk(A):
    """0-1 矩陣（§A.4.3）：字典序陣列，依序行插入 v_1, …, v_r，並把 u_i 放在新格子的共軛（轉置）位置。回傳 (P, Q)。"""
    u, v = lex_sort(*zip(*matrix_to_pairs(A)))
    P, Qt = (), ()
    for a, b in zip(u, v):
        P, _, nb = column_insert(P, b)
        Qt = _place(Qt, (nb[1], nb[0]), a)
    return P, Qt


def zero_one_matrices(rows, cols):
    m, n = len(rows), len(cols)
    out = []
    for choice in product(*(combinations(range(n), r) for r in rows)):
        A = [[1 if j in c else 0 for j in range(n)] for c in choice]
        if all(sum(A[i][j] for i in range(m)) == cols[j] for j in range(n)):
            out.append(A)
    return out


def gale_ryser(lam, mu) -> bool:
    """練習 20：列和 λ、行和 μ 的 0-1 矩陣存在 ⇔ μ 被 λ 的共軛支配（μ ⊴ λ~）。"""
    lt = conjugate(lam)
    return sum(lam) == sum(mu) and dominates(lt, tuple(sorted(mu, reverse=True)))


# ---------------------------------------------------------------- §A.5 鑰匙
def compact_skew(col_lengths):
    """回傳每一行的 (頂端列, 長度)：相鄰兩行左邊較長則頂端對齊，右邊較長則底端對齊；最後整體平移使最上面是第 0 列。"""
    tops = [0]
    for a, b in zip(col_lengths, col_lengths[1:]):
        t = tops[-1]
        tops.append(t if a >= b else t + a - b)
    m = min(tops)
    return [(t - m, L) for t, L in zip(tops, col_lengths)]


def _skew_from_columns(cols_geo, fill):
    rows = max(t + L for t, L in cols_geo)
    grid = [[None] * len(cols_geo) for _ in range(rows)]
    for j, ((t, L), col) in enumerate(zip(cols_geo, fill)):
        for k in range(L):
            grid[t + k][j] = col[k]
    # 斜表的寫法：每列左邊的空格是 None；右邊多餘的 None 去掉
    out = []
    for r in grid:
        while r and r[-1] is None:
            r = r[:-1]
        out.append(tuple(r))
    return tuple(out)


def frank_tableau(T, col_lengths):
    """形狀 compact_skew(col_lengths)、rectify 成 T 的斜表（唯一）；回傳各行的內容（由上而下）。"""
    T = tuple(tuple(r) for r in T)
    assert sorted(col_lengths) == sorted(len(c) for c in _columns(T))
    entries = sorted(x for r in T for x in r)
    geo = compact_skew(col_lengths)
    found = []

    def rec(j, remaining, cols):
        if j == len(geo):
            S = _skew_from_columns(geo, cols)
            if _is_skew_tableau(S) and rectify(S) == T:
                found.append(tuple(cols))
            return
        L = geo[j][1]
        for c in sorted(set(combinations(remaining, L))):
            rest = list(remaining)
            for x in c:
                rest.remove(x)
            rec(j + 1, tuple(rest), cols + [c])

    rec(0, tuple(entries), [])
    assert len(found) == 1, (col_lengths, found)
    return found[0]


def _is_skew_tableau(S):
    for i, r in enumerate(S):
        vals = [x for x in r if x is not None]
        if any(a > b for a, b in zip(vals, vals[1:])):
            return False
        for j, x in enumerate(r):
            if x is not None and i + 1 < len(S) and j < len(S[i + 1]) and S[i + 1][j] is not None and S[i + 1][j] <= x:
                return False
    # 行要連續：skew 的空格只能在每列的左邊
    for r in S:
        seen = False
        for x in r:
            if x is not None:
                seen = True
            elif seen:
                return False
    return True


def _columns(T):
    return [tuple(T[i][j] for i in range(len(T)) if j < len(T[i])) for j in range(len(T[0]))]


def keys(T):
    """左鑰匙 K_-(T)：長度 c 的行換成「最左行長度為 c 的 frank 斜表」的最左行；右鑰匙 K_+(T) 同理取最右行。"""
    cols = _columns(T)
    lengths = [len(c) for c in cols]
    left, right = {}, {}
    for c in set(lengths):
        others = list(lengths)
        others.remove(c)
        lf = frank_tableau(T, [c] + sorted(others, reverse=True))
        rt = frank_tableau(T, sorted(others, reverse=True) + [c])
        left[c], right[c] = lf[0], rt[-1]
    build = lambda d: transpose(tuple(tuple(sorted(d[L])) for L in lengths))
    return build(left), build(right)


# ---------------------------------------------------------------- 附錄 B
def betti_from_cells(counts):
    """counts[k] = 複維數 k 的胞腔個數 ⇒ b_{2k} = counts[k]、奇數 Betti 數為 0；χ = Σ counts。"""
    b = []
    for c in counts:
        b += [c, 0]
    return b[:-1], sum(counts)


def chern_projective(n):
    return [comb(n + 1, k) for k in range(n + 1)]


def plucker_degree(r, m) -> int:
    return len(standard_tableaux((m - r,) * r))


# ---------------------------------------------------------------- 自我測試
def check():
    from itertools import permutations as _perms
    # §A.1：書 p.184 的例子（m = 6，a* = 7 − a）：ΔT、第一個 1* 放在第 4 列第 2 格；T* = 6*6*5*4*/5*5*4*2*/3*3*2*/2*1*
    T = ((1, 2, 2, 3), (2, 3, 5, 5), (4, 4, 6), (5, 6))
    steps = evacuation_steps(T, 6)
    assert steps[0][0] == ((2, 2, 2, 3), (3, 4, 5, 5), (4, 6, 6), (5,)) and steps[0][1] == (3, 1) and steps[0][2] == 6
    Ts = evacuation(T, 6)
    assert Ts == ((1, 1, 2, 3), (2, 2, 3, 5), (4, 4, 5), (5, 6)) and evacuation(Ts, 6) == T
    # Duality theorem (3)：P(w*) = P(w)*（對所有 [4] 上長度 5 的字中抽樣）
    try:
        from combkit.tableaux import insert_word
    except ImportError:
        from tableaux import insert_word  # type: ignore
    for w in product(range(1, 4), repeat=5):
        assert insert_word(dual_word(w, 3)) == evacuation(insert_word(w), 3)
    # §A.2：行插入依序插入 v_r, …, v_1 得到同一個 P；相異字母時 P(w^rev) = P(w)^T
    for w in list(product(range(1, 4), repeat=5))[::7]:
        assert column_insert_word(tuple(reversed(w))) == insert_word(w)
    for w in _perms(range(1, 6)):
        assert insert_word(tuple(reversed(w))) == transpose(insert_word(w))
    assert column_insert(((1, 2), (3,)), 2)[0] == ((1, 2), (2, 3))      # 第一行擠掉 3，3 去第二行的底部
    # §A.4.1–4.2：書 p.199 的矩陣
    A = [[1, 2, 0], [0, 1, 1], [1, 1, 1]]
    P, Q = rsk_array(*zip(*matrix_to_pairs(A)))
    R, S = burge(A)
    assert shape(R) == shape(S)
    c = corners(A)
    assert c["nw"] == (P, Q) and c["se"] == (evacuation(P, 3), evacuation(Q, 3))
    assert c["sw"] == (R, evacuation(S, 3)) and c["ne"] == (evacuation(R, 3), S)
    # 再用所有 2×2、元素 0..2 的矩陣核對定理 1
    for a in product(range(3), repeat=4):
        B = [list(a[:2]), list(a[2:])]
        if not any(a):
            continue
        P2, Q2 = rsk_array(*zip(*matrix_to_pairs(B)))
        R2, S2 = burge(B)
        cc = corners(B)
        assert cc["se"] == (evacuation(P2, 2), evacuation(Q2, 2)) and cc["sw"] == (R2, evacuation(S2, 2))
        assert cc["ne"] == (evacuation(R2, 2), S2)
    # §A.4.3：0-1 矩陣 ↔ 共軛形狀的表對（3×3 全部 511 個非零矩陣，一一對應）
    seen = set()
    for bits in product((0, 1), repeat=9):
        if not any(bits):
            continue
        M = [list(bits[0:3]), list(bits[3:6]), list(bits[6:9])]
        P3, Q3 = dual_rsk(M)
        assert shape(P3) == conjugate(shape(Q3))
        assert sorted(x for r in P3 for x in r) == sorted(j for _, j in matrix_to_pairs(M))
        assert sorted(x for r in Q3 for x in r) == sorted(i for i, _ in matrix_to_pairs(M))
        seen.add((P3, Q3))
    assert len(seen) == 511
    # 練習 19、20：個數 = Σ_ν K_{ν~,λ} K_{ν,μ}；存在 ⇔ μ ⊴ λ~
    for n in (3, 4, 5):
        for lam in partitions(n):
            for mu in partitions(n):
                cnt = len(zero_one_matrices(lam, mu))
                assert cnt == sum(kostka(conjugate(nu), lam) * kostka(nu, mu) for nu in partitions(n))
                assert (cnt > 0) == gale_ryser(lam, mu)
    # §A.5：書 p.210 的例子
    T5 = ((1, 1, 2, 2), (2, 3, 3), (4,))
    Km, Kp = keys(T5)
    assert Km == ((1, 1, 1, 1), (2, 2, 2), (4,)) and Kp == ((2, 2, 2, 2), (3, 3, 3), (4,))
    # 附錄 B
    assert betti_from_cells([1, 1, 2, 1, 1]) == ([1, 0, 1, 0, 2, 0, 1, 0, 1], 6)
    assert chern_projective(2) == [1, 3, 3] and chern_projective(3)[3] == 4
    assert plucker_degree(2, 4) == 2 and plucker_degree(2, 5) == 5 and plucker_degree(3, 6) == 42
    return True


if __name__ == "__main__":
    check()
    print("Fulton 附錄 A–B：evacuation、行插入、Burge 與四個角、0-1 矩陣、鑰匙、胞腔與 Chern 類全部驗證通過")
