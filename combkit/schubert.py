"""
schubert.py —— Fulton《Young Tableaux》第 9–10 章（場景 F46–F55）的計算（純 Python，不依賴 manim）

置換用一行記法的 tuple（值 1..m）；w·s_i = 交換第 i、i+1 個位置（act_right）。多項式沿用 tableaux.Poly
（{指數 tuple: 係數}，tuple 長度 = 變數個數）。

    length(w)、inversions(w)、descents(w)、act_right(w, i)、w0(m)、perms(m)、lehmer_code(w)
    rank_matrix(w)            r_w(p,q) = #{i ≤ p : w(i) ≤ q}（Fulton §10.5）
    bruhat_leq(u, v)          u ≤ v ⇔ r_u ≥ r_v；bruhat_leq_tableau（練習 8：前綴排序後逐項比較）；bruhat_chain（Lemma 11）
    bruhat_levels(m)          依長度分層；bruhat_covers(m) 覆蓋關係（Hasse 圖的邊）
    flag_cell_pattern(w)      Schubert 胞腔 X°_w 的列梯形矩陣（'1' '*' '0'），星號數 = ℓ(w)（§10.2）
    grass_cell_pattern(lam, r, n)   Gr 的胞腔 Ω°_λ：第 i 列的 1 在第 n+i−λ_i 行，星號 n−λ_i 個（§9.4）
    gaussian_binomial(m, r)、q_factorial(m)   係數 list（q^0, q^1, …）；count_subspaces_F2、count_flags_F2 暴力數 F_2 上的點
    plucker(A)、plucker_from_coords(x, d, m)   最大子式；由 Plücker 座標還原矩陣（§9.1 式 (2)、練習 2）
    flag_ring_dim(a, b)       S(3;1,2) 在雙次數 (a,b) 的維數 = C(a+2,2)C(b+2,2) − C(a+1,2)C(b+1,2)
    box_product(lam, mu, r, n)   σ_λσ_μ = Σ c σ_ν（ν 在 r×n 方框內，LR 係數）；rectangle_count(r, m) 練習 17
    swap_vars(f, i)、divided_difference(f, i)   ∂_i f = (f − s_i f)/(x_i − x_{i+1})（§10.3）
    schubert_poly(w)          𝔖_w：由 𝔖_{w0} = x1^{m-1}…x_{m-1} 用 ∂_i 往下（§10.4）；schubert_path(m, steps) 逐步記錄
    pipe_dreams(w)            縮減 pipe dream（Billey–Jockusch–Stanley）：Σ x^D = 𝔖_w
    grassmannian_partition(w)、grassmannian_perm(lam, r, m)   命題 8：𝔖_w = s_λ(x_1..x_r)
    staircase_monomials(m)    H*(Fl(m)) = Z[x]/(e_1..e_m) 的基底 x^I，i_j ≤ m−j（命題 3）
    poly_str(f)               多項式的 LaTeX 字串（場景用）

直接執行：python utils/schubert.py（在 codes/ 內）跑自我測試。
"""
from __future__ import annotations

from functools import lru_cache
from itertools import combinations, permutations, product
from math import comb, factorial, prod

try:
    from combkit.tableaux import (conjugate, hook_content_count, lr_coefficient, poly_mul, poly_sum, schur_expand,
                                schur_poly, standard_tableaux)
except ImportError:                                   # 在 utils/ 裡直接執行時
    from tableaux import (conjugate, hook_content_count, lr_coefficient, poly_mul, poly_sum, schur_expand,  # type: ignore
                          schur_poly, standard_tableaux)


# ---------------------------------------------------------------- 置換
def length(w) -> int:
    return sum(w[i] > w[j] for i in range(len(w)) for j in range(i + 1, len(w)))


def inversions(w):
    return [(i + 1, j + 1) for i in range(len(w)) for j in range(i + 1, len(w)) if w[i] > w[j]]


def descents(w):
    return [i + 1 for i in range(len(w) - 1) if w[i] > w[i + 1]]


def act_right(w, i):
    """w·s_i：交換第 i、i+1 個位置的值（1 起算）。"""
    w = list(w)
    w[i - 1], w[i] = w[i], w[i - 1]
    return tuple(w)


def w0(m):
    return tuple(range(m, 0, -1))


def perms(m):
    return [tuple(p) for p in permutations(range(1, m + 1))]


def lehmer_code(w):
    return tuple(sum(w[j] < w[i] for j in range(i + 1, len(w))) for i in range(len(w)))


def perm_str(w) -> str:
    return "".join(map(str, w))


# ---------------------------------------------------------------- Bruhat 序（§10.5）
def rank_matrix(w):
    m = len(w)
    return [[sum(1 for i in range(p) if w[i] <= q) for q in range(1, m + 1)] for p in range(1, m + 1)]


def bruhat_leq(u, v) -> bool:
    ru, rv = rank_matrix(u), rank_matrix(v)
    return all(a >= b for x, y in zip(ru, rv) for a, b in zip(x, y))


def bruhat_leq_tableau(u, v) -> bool:
    """練習 8：每個 p，把 u(1..p)、v(1..p) 各自排序後，u 的每一項 ≤ v 的對應項。"""
    return all(all(a <= b for a, b in zip(sorted(u[:p]), sorted(v[:p]))) for p in range(1, len(u)))


def bruhat_chain(u, v):
    """Lemma 11 的演算法：從 v 一步步往下（每步交換一對、長度減 1）；到得了 u 就回傳整條鏈，否則 None。"""
    chain = [tuple(v)]
    while chain[-1] != tuple(u):
        x = chain[-1]
        j = next(i for i in range(len(u)) if u[i] != x[i])
        if x[j] < u[j]:
            return None
        ks = [k for k in range(j + 1, len(u)) if x[j] > x[k] >= u[j]]
        if not ks:
            return None
        k = ks[0]
        y = list(x)
        y[j], y[k] = y[k], y[j]
        chain.append(tuple(y))
    return chain


def bruhat_levels(m):
    lv = {}
    for w in perms(m):
        lv.setdefault(length(w), []).append(w)
    return [sorted(lv[k]) for k in sorted(lv)]


def bruhat_covers(m):
    lv = bruhat_levels(m)
    return [(u, v) for k in range(len(lv) - 1) for u in lv[k] for v in lv[k + 1] if bruhat_leq(u, v)]


# ---------------------------------------------------------------- 胞腔與點數
def flag_cell_pattern(w):
    """第 p 列在第 w(p) 行放 1，1 的右邊是 0，每個 1 的下方是 0，其餘是 *（Fulton p.157）。"""
    m = len(w)
    used, out = set(), []
    for p in range(m):
        row = []
        for j in range(1, m + 1):
            row.append("1" if j == w[p] else "*" if j < w[p] and j not in used else "0")
        used.add(w[p])
        out.append(row)
    return out


def grass_cell_pattern(lam, r, n):
    """Gr_r(C^{r+n}) 的胞腔 Ω°_λ（Fulton p.147）：第 i 列的 1 在第 n+i−λ_i 行。"""
    lam = tuple(lam) + (0,) * (r - len(lam))
    piv = [n + i + 1 - lam[i] for i in range(r)]
    return [["1" if j == piv[i] else "*" if j < piv[i] and j not in piv[:i] else "0" for j in range(1, r + n + 1)]
            for i in range(r)]


def _qpoly_mul(a, b):
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return out


def q_factorial(m):
    out = [1]
    for k in range(1, m + 1):
        out = _qpoly_mul(out, [1] * k)
    return out


def gaussian_binomial(m, r):
    """[m choose r]_q 的係數 = 各維數的 Schubert 胞腔個數（方框 r×(m−r) 裡 |λ| = k 的 λ 個數）。"""
    box = [lam for lam in product(range(m - r + 1), repeat=r) if list(lam) == sorted(lam, reverse=True)]
    out = [0] * (r * (m - r) + 1)
    for lam in box:
        out[sum(lam)] += 1
    return out


def _span(vecs):
    s = {0}
    for v in vecs:
        s |= {x ^ v for x in s}
    return frozenset(s)


def count_subspaces_F2(m, r) -> int:
    """F_2^m 裡 r 維子空間的個數（暴力：所有 r 個向量張成的空間去重）。"""
    return len({S for vs in combinations(range(1, 2 ** m), r) if len(S := _span(vs)) == 2 ** r})


def count_flags_F2(m) -> int:
    """F_2^m 裡完全旗 0 ⊂ V_1 ⊂ … ⊂ V_{m-1} ⊂ F_2^m 的個數（暴力）。"""
    layers = [{_span(vs) for vs in combinations(range(1, 2 ** m), k) if len(_span(vs)) == 2 ** k}
              for k in range(1, m)]
    cnt = {S: 1 for S in layers[0]}
    for k in range(1, m - 1):
        cnt = {T: sum(c for S, c in cnt.items() if S <= T) for T in layers[k]}
    return sum(cnt.values())


# ---------------------------------------------------------------- Plücker（§9.1）
def _det(M):
    if len(M) == 1:
        return M[0][0]
    return sum((-1) ** j * M[0][j] * _det([r[:j] + r[j + 1:] for r in M[1:]]) for j in range(len(M)))


def plucker(A):
    d, m = len(A), len(A[0])
    return {I: _det([[A[s][i - 1] for i in I] for s in range(d)]) for I in combinations(range(1, m + 1), d)}


def _coord(x, J):
    """可交錯下標的座標：重複下標 → 0，排序時每次對換變號。"""
    if len(set(J)) < len(J):
        return 0
    J = list(J)
    sgn = 1
    for i in range(len(J)):
        for j in range(len(J) - 1 - i):
            if J[j] > J[j + 1]:
                J[j], J[j + 1] = J[j + 1], J[j]
                sgn = -sgn
    return sgn * x[tuple(J)]


def plucker_from_coords(x, d, m):
    """Fulton p.133 式 (2)：取 x_I ≠ 0 的 I（正規化成 1），a_{s,t} = x_{I 的第 s 個換成 t}。"""
    from fractions import Fraction
    I = next(I for I in sorted(x) if x[I])
    c = Fraction(1, x[I])
    return [[c * _coord(x, I[:s] + (t,) + I[s + 1:]) for t in range(1, m + 1)] for s in range(d)]


def flag_ring_dim(a, b) -> int:
    """S(3;1,2) = Sym(E ⊕ ∧²E)/(關聯關係) 在雙次數 (a,b) 的維數；關係是非零因子，所以扣掉 (a−1,b−1) 的維數。"""
    return comb(a + 2, 2) * comb(b + 2, 2) - (comb(a + 1, 2) * comb(b + 1, 2) if a and b else 0)


# ---------------------------------------------------------------- Schubert 計算（§9.4）
def _partitions_in_box(k, r, n):
    return [lam for lam in product(range(n, -1, -1), repeat=r)
            if sum(lam) == k and all(lam[i] >= lam[i + 1] for i in range(r - 1))]


def box_product(lam, mu, r, n):
    lam = tuple(x for x in lam if x)
    mu = tuple(x for x in mu if x)
    out = {}
    for nu in _partitions_in_box(sum(lam) + sum(mu), r, n):
        nu = tuple(x for x in nu if x)
        c = lr_coefficient(lam, mu, nu)
        if c:
            out[nu] = c
    return out


def rectangle_count(r, m) -> int:
    """練習 17：與 r(m−r) 個一般的 (m−r) 維子空間都相交的 r 維子空間個數。"""
    num = factorial(r * (m - r)) * prod(factorial(i) for i in range(1, r))
    den = prod(factorial(m - i) for i in range(1, r + 1))
    assert num % den == 0
    return num // den


# ---------------------------------------------------------------- 差商與 Schubert 多項式（§10.3–10.4）
def swap_vars(f, i):
    out = {}
    for e, c in f.items():
        e = list(e)
        e[i - 1], e[i] = e[i], e[i - 1]
        out[tuple(e)] = out.get(tuple(e), 0) + c
    return {e: c for e, c in out.items() if c}


def divided_difference(f, i):
    """∂_i：x_i^a x_{i+1}^b ↦ Σ_{k=0}^{a−b−1} x_i^{a−1−k} x_{i+1}^{b+k}（a > b），a < b 時取負號，a = b 時為 0。"""
    out = {}
    for e, c in f.items():
        a, b = e[i - 1], e[i]
        if a == b:
            continue
        hi, lo, s = (a, b, 1) if a > b else (b, a, -1)
        for k in range(hi - lo):
            g = list(e)
            g[i - 1], g[i] = hi - 1 - k, lo + k
            out[tuple(g)] = out.get(tuple(g), 0) + s * c
    return {e: c for e, c in out.items() if c}


@lru_cache(maxsize=None)
def _schubert(w):
    m = len(w)
    if w == w0(m):
        return tuple(sorted({tuple(m - 1 - j for j in range(m)): 1}.items()))
    i = next(i for i in range(1, m) if w[i - 1] < w[i])
    return tuple(sorted(divided_difference(dict(_schubert(act_right(w, i))), i).items()))


def schubert_poly(w):
    """𝔖_w，用 len(w) 個變數（Lemma 10 與命題 5：𝔖_{w·s_i} = ∂_i 𝔖_w 當 w(i) > w(i+1)）。"""
    return dict(_schubert(tuple(w)))


def schubert_path(m, steps):
    """從 w0 出發依序做 ∂_i（i 在 steps 裡）；每步都要求目前 w(i) > w(i+1)。回傳 [(w, 𝔖_w), …]。"""
    w = w0(m)
    out = [(w, schubert_poly(w))]
    for i in steps:
        assert w[i - 1] > w[i]
        f = divided_difference(out[-1][1], i)
        w = act_right(w, i)
        assert f == schubert_poly(w)
        out.append((w, f))
    return out


def pipe_dreams(w):
    """縮減 pipe dream：階梯 {(i,j): i+j ≤ m} 中的十字集合 D，逐列由上而下、每列由右而左讀，
    十字 (i,j) 記 s_{i+j−1}，乘起來是 w 的縮減字（長度 ℓ(w)）。回傳 D 的 list（每個 D 是排好序的格子 tuple）。"""
    m = len(w)
    cells = [(i, j) for i in range(1, m) for j in range(1, m) if i + j <= m]
    out = []
    for D in combinations(cells, length(w)):
        x = tuple(range(1, m + 1))
        ok = True
        for (i, j) in sorted(D, key=lambda c: (c[0], -c[1])):
            k = i + j - 1
            if x[k - 1] > x[k]:                        # 不是縮減字
                ok = False
                break
            x = act_right(x, k)
        if ok and x == tuple(w):
            out.append(tuple(sorted(D)))
    return out


def pipe_exits(D, m):
    """每條 pipe 從第 i 列的左邊進來、往右或往上走（十字直走；彎管：從左進→往上出，從下進→往右出），
    回傳 tuple：第 i 條 pipe 從第幾行的頂端出去。縮減 pipe dream 時這就是 w（自我測試核對）。"""
    D = set(D)
    out = []
    for i in range(1, m + 1):
        r, c, d = i, 1, "R"                            # 目前所在格子、前進方向（R = 從左邊進來、U = 從下面進來）
        while True:
            cross = (r, c) in D and r + c <= m
            if (d == "R") != cross:                    # 從左進的彎管、或從下進的十字：往上
                if r == 1:
                    out.append(c)
                    break
                r, d = r - 1, "U"
            else:                                      # 從左進的十字、或從下進的彎管：往右
                c, d = c + 1, "R"
    return tuple(out)


def pipe_weight(D, m):
    e = [0] * m
    for i, _ in D:
        e[i - 1] += 1
    return tuple(e)


def grassmannian_partition(w):
    """w 只在 r 有下降 ⇒ λ = (w(r)−r, …, w(1)−1)（命題 8）；回傳 (λ, r)。"""
    d = descents(w)
    assert len(d) == 1
    r = d[0]
    return tuple(w[r - 1 - k] - (r - k) for k in range(r)), r


def grassmannian_perm(lam, r, m):
    lam = tuple(lam) + (0,) * (r - len(lam))
    first = [lam[r - 1 - i] + i + 1 for i in range(r)]
    return tuple(first + sorted(set(range(1, m + 1)) - set(first)))


def staircase_monomials(m):
    return [e + (0,) for e in product(*(range(m - j + 1) for j in range(1, m)))]


def poly_str(f, var="x") -> str:
    """{(2,1,0): 1, (1,1,1): 2} → 'x_1^{2}x_2+2x_1x_2x_3'（字典序由大到小）。"""
    terms = []
    for e, c in sorted(f.items(), reverse=True):
        mono = "".join(f"{var}_{i + 1}" + (f"^{{{k}}}" if k > 1 else "") for i, k in enumerate(e) if k)
        coef = "" if (abs(c) == 1 and mono) else str(abs(c))
        terms.append(("-" if c < 0 else "+") + coef + (mono or ""))
    s = "".join(terms) or "0"
    return s[1:] if s.startswith("+") else s


# ---------------------------------------------------------------- 自我測試
def check():
    # §9.1 Plücker：F45 的矩陣、練習 2
    A = [[1, 2, 0, 1], [0, 1, 3, 2]]
    x = plucker(A)
    assert [x[I] for I in sorted(x)] == [1, 3, 2, 6, 3, -3]
    rel = lambda x: x[(1, 2)] * x[(3, 4)] - x[(1, 3)] * x[(2, 4)] + x[(1, 4)] * x[(2, 3)]
    assert rel(x) == 0
    ex2 = {(1, 2): 1, (1, 3): 2, (1, 4): 1, (2, 3): 1, (2, 4): 2, (3, 4): 3}
    assert rel(ex2) == 1 * 3 - 2 * 2 + 1 * 1 == 0
    B = plucker_from_coords(ex2, 2, 4)
    assert B == [[1, 0, -1, -2], [0, 1, 2, 1]] and plucker(B) == ex2
    assert all(sum(B[s][t] * v[t] for t in range(4)) == 0 for s in range(2) for v in ((1, -2, 1, 0), (2, -1, 0, 1)))
    # 旗流形的座標環：S(3;1,2) 的雙次數 (a,b) 部分 = E^{(a+b,b)}
    for a in range(4):
        for b in range(4):
            assert flag_ring_dim(a, b) == hook_content_count((a + b, b), 3)
    assert flag_ring_dim(1, 1) == 8 and flag_ring_dim(2, 1) == 15
    # §9.4 Gr 的胞腔
    assert gaussian_binomial(4, 2) == [1, 1, 2, 1, 1] and sum(gaussian_binomial(4, 2)) == 6
    assert sum(c * 2 ** k for k, c in enumerate(gaussian_binomial(4, 2))) == count_subspaces_F2(4, 2) == 35
    assert count_subspaces_F2(3, 1) == 7 and count_subspaces_F2(4, 1) == 15
    pat = grass_cell_pattern((5, 3, 2, 2, 1), 5, 7)
    assert [r.count("*") for r in pat] == [2, 4, 5, 5, 6] and sum(r.count("*") for r in pat) == 35 - 13
    for lam in _partitions_in_box(0, 2, 2) + _partitions_in_box(1, 2, 2) + _partitions_in_box(2, 2, 2) \
            + _partitions_in_box(3, 2, 2) + _partitions_in_box(4, 2, 2):
        assert sum(r.count("*") for r in grass_cell_pattern(lam, 2, 2)) == 4 - sum(lam)
    # Schubert 計算：四條線、五條、LR 在方框裡
    one = (1,)
    s = {(): 1}
    for _ in range(4):
        s = poly_sum([{nu: c * k for nu, k in box_product(lam, one, 2, 2).items()} for lam, c in s.items()])
    assert s == {(2, 2): 2}
    assert rectangle_count(2, 4) == 2 == len(standard_tableaux((2, 2)))
    assert rectangle_count(2, 5) == 5 == len(standard_tableaux((3, 3)))
    assert rectangle_count(3, 6) == len(standard_tableaux((3, 3, 3))) == 42
    assert box_product((1,), (1,), 2, 2) == {(2,): 1, (1, 1): 1}
    assert box_product((2, 1), (2, 1), 3, 3) == {(3, 3): 1, (3, 2, 1): 2, (2, 2, 2): 1}
    # 對偶：λ 與 方框補集
    assert box_product((2, 1), (1,), 2, 2) == {(2, 2): 1}
    # §10.1–10.2 旗流形
    assert len(perms(3)) == 6 and q_factorial(3) == [1, 2, 2, 1] and q_factorial(4) == [1, 3, 5, 6, 5, 3, 1]
    assert [len(l) for l in bruhat_levels(3)] == [1, 2, 2, 1] and [len(l) for l in bruhat_levels(4)] == q_factorial(4)
    assert count_flags_F2(3) == 21 == sum(c * 2 ** k for k, c in enumerate(q_factorial(3)))
    assert count_flags_F2(4) == 315 == sum(c * 2 ** k for k, c in enumerate(q_factorial(4)))
    w = (4, 2, 6, 1, 3, 5)
    assert length(w) == 7 and sum(r.count("*") for r in flag_cell_pattern(w)) == 7
    assert all(sum(r.count("*") for r in flag_cell_pattern(v)) == length(v) for v in perms(4))
    assert len(staircase_monomials(3)) == 6 and len(staircase_monomials(4)) == 24
    # §10.3 差商
    X = lambda *e: {tuple(e): 1}
    assert divided_difference(X(2, 0, 0), 1) == {(1, 0, 0): 1, (0, 1, 0): 1}
    assert divided_difference(X(1, 1, 0), 1) == {} and divided_difference(X(0, 1, 0), 1) == {(0, 0, 0): -1}
    f = poly_sum([X(3, 1, 0), X(0, 2, 2), {(1, 0, 1): 5}])
    g = poly_sum([X(1, 2, 0), {(0, 0, 3): -2}])
    dd = divided_difference
    assert dd(dd(f, 1), 1) == {}
    assert dd(dd(dd(f, 1), 2), 1) == dd(dd(dd(f, 2), 1), 2)
    assert dd(poly_mul(f, g), 1) == poly_sum([poly_mul(dd(f, 1), g), poly_mul(swap_vars(f, 1), dd(g, 1))])
    # §10.4 Schubert 多項式
    S3 = {perm_str(v): schubert_poly(v) for v in perms(3)}
    assert S3 == {"123": {(0, 0, 0): 1}, "213": {(1, 0, 0): 1}, "132": {(1, 0, 0): 1, (0, 1, 0): 1},
                  "231": {(1, 1, 0): 1}, "312": {(2, 0, 0): 1}, "321": {(2, 1, 0): 1}}
    path = schubert_path(5, [1, 4, 3, 2, 3])                  # Fulton p.172：54321 → … → 41352
    assert [perm_str(v) for v, _ in path] == ["54321", "45321", "45312", "45132", "41532", "41352"]
    s41352 = schubert_poly((4, 1, 3, 5, 2))
    assert s41352 == {(3, 1, 1, 0, 0): 1, (3, 1, 0, 1, 0): 1, (3, 0, 1, 1, 0): 1}      # x1^3(x2x3+x2x4+x3x4)
    s21543 = schubert_poly((2, 1, 5, 4, 3))
    assert s21543[(2, 1, 1, 0, 0)] == 2 and all(c > 0 for c in s21543.values())
    # pipe dreams = Schubert 多項式（S4 全部、21543）
    for v in perms(4) + [(2, 1, 5, 4, 3), (4, 1, 3, 5, 2)]:
        cnt = {}
        for D in pipe_dreams(v):
            e = pipe_weight(D, len(v))
            cnt[e] = cnt.get(e, 0) + 1
        assert cnt == schubert_poly(v), v
    assert sum(pipe_weight(D, 5) == (2, 1, 1, 0, 0) for D in pipe_dreams((2, 1, 5, 4, 3))) == 2
    for v in perms(4) + [(2, 1, 5, 4, 3)]:
        assert all(pipe_exits(D, len(v)) == v for D in pipe_dreams(v)), v
    # §10.5 Bruhat：秩矩陣判準 = 前綴判準；Lemma 11 的鏈（Fulton p.175）
    for u in perms(4):
        for v in perms(4):
            assert bruhat_leq(u, v) == bruhat_leq_tableau(u, v) == (bruhat_chain(u, v) is not None)
    ch = bruhat_chain((4, 2, 8, 3, 6, 1, 7, 9, 5), (6, 7, 9, 2, 5, 1, 8, 3, 4))
    assert len(ch) == 9 and all(length(a) == length(b) + 1 for a, b in zip(ch, ch[1:]))
    assert bruhat_chain((4, 2, 8, 6, 7, 1, 9, 5, 3), (6, 7, 9, 2, 5, 1, 8, 3, 4)) is None
    assert len(bruhat_covers(3)) == 8
    # §10.6 Grassmann 置換 = Schur 多項式（S5 全部）
    for v in perms(5):
        if len(descents(v)) == 1:
            lam, r = grassmannian_partition(v)
            assert grassmannian_perm(lam, r, 5) == v
            sw = schubert_poly(v)
            sl = schur_poly(tuple(x for x in lam if x) or (0,), r) if any(lam) else {(0,) * r: 1}
            assert sw == {e + (0,) * (5 - r): c for e, c in sl.items()}, v
    assert grassmannian_partition((2, 4, 1, 3)) == ((2, 1), 2) and grassmannian_partition((1, 3, 4, 2)) == ((1, 1, 0), 3)
    assert schubert_poly((1, 3, 2, 4)) == {(1, 0, 0, 0): 1, (0, 1, 0, 0): 1}
    # 四條線再算一次：s_1(x1,x2)^4 = s_4 + 3 s_31 + 2 s_22（兩個變數）
    s1 = schur_poly((1,), 2)
    p4 = poly_mul(poly_mul(s1, s1), poly_mul(s1, s1))
    assert schur_expand(p4, 2) == {(4,): 1, (3, 1): 3, (2, 2): 2}
    assert schubert_poly((3, 4, 1, 2)) == {(2, 2, 0, 0): 1}
    return True


if __name__ == "__main__":
    check()
    print("Fulton 第 9–10 章：Plücker、Schubert 胞腔與計算、Bruhat 序、差商、Schubert 多項式全部驗證通過")
