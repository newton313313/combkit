"""
young_extras.py —— youngart（Fulton《Young Tableaux》全書 60 張 4K 藝術桌布）用到、可以重複使用的計算（Python＋numpy）

表用 tableaux.py 的寫法（列的 tuple，格子 (列, 行) 0 起算）；排列 w 是 1..n 的 tuple。

組合
    patience(w)                patience sorting：(每個元素放進的堆號, 一組最長遞增子序列的位置)
    slide_path(T)              拿掉左上角後做一次 jeu de taquin：(洞經過的格子, 新表)；evacuation 每一步就是這個
    bump_routes(w)             列插入 w 的每條擠壓路徑，與倒序做行插入的每條路徑；兩者得到同一個 P（附錄 A.2）
    shadow_lines(w)            Viennot 影子線：點 (i, w(i)) 依西北方最長鏈編號，同號的點連成階梯；條數 = 最長遞增子序列長
    four_corners(w)            從四個角各做一次 shadow_lines（翻轉列／行後計算再翻回來）：{"nw","ne","sw","se": [折線…]}
    rsk_rows(w)                排列的 R-S-K（bisect 版，n 上萬也很快）：回傳 (P, Q) 的列串列
幾何（numpy）
    spectral_xy(n, edges)      圖的譜嵌入（正規化拉普拉斯第 2、3 個特徵向量）
    pipe_route(D, m, i, cs)    pipe dream 第 i 條水管經過的點（十字直走，彎管畫成四分之一圓；原點是左上角、y 往上）
    ruling(fam, th, t)         單葉雙曲面 x²+y²−z² = 1 的兩族直線：fam = "A" 或 "B"、角度 th、高度 t
    line_meets_hyperboloid(p0, d)   直線 p0 + s·d 與雙曲面的交點，以及過交點的 B 族直線的角度（= 與 A 族三條線都相交的「橫截線」）
    ellipse_curve(a, b, phi, c)、line_curve(p0, ang)、cubic_curve(c, s, k)   平面曲線：(參數式, 隱函數, 參數範圍)
    meet(A, B)                 沿 A 的參數找 B 的隱函數變號再二分，回傳交點（Bézout 桌布：每對曲線 d·e 個實交點）

直接執行：python -m combkit.young_extras（在 combkit 的上層資料夾）跑自我測試。
"""
from __future__ import annotations

import bisect

import numpy as np

try:
    from combkit.fulton_appendix import column_insert, evacuation_steps
    from combkit.schubert import pipe_dreams, pipe_exits
    from combkit.tableaux import row_insert, rsk, shape
except ImportError:                                   # 在 combkit/ 裡直接執行時
    from fulton_appendix import column_insert, evacuation_steps  # type: ignore
    from schubert import pipe_dreams, pipe_exits  # type: ignore
    from tableaux import row_insert, rsk, shape  # type: ignore

TAU = 2 * np.pi


# ---------------------------------------------------------------- 組合
def patience(w):
    tops, tops_idx, pile, prev = [], [], [], []
    for k, x in enumerate(w):
        p = bisect.bisect_left(tops, x)
        if p == len(tops):
            tops.append(x)
            tops_idx.append(k)
        else:
            tops[p] = x
            tops_idx[p] = k
        pile.append(p)
        prev.append(tops_idx[p - 1] if p > 0 else -1)
    lis, k = [], tops_idx[-1]
    while k >= 0:
        lis.append(k)
        k = prev[k]
    return pile, lis[::-1]


def slide_path(T):
    """洞從 (0,0) 出發往右或往下（取較小者，相等取下面），走到外角為止。"""
    R = [list(r) for r in T]
    i = j = 0
    cells = [(0, 0)]
    while True:
        rt = R[i][j + 1] if j + 1 < len(R[i]) else None
        dn = R[i + 1][j] if i + 1 < len(R) and j < len(R[i + 1]) else None
        if rt is None and dn is None:
            break
        if dn is not None and (rt is None or dn <= rt):
            R[i][j], i = dn, i + 1
        else:
            R[i][j], j = rt, j + 1
        cells.append((i, j))
    del R[i][j]
    return cells, tuple(tuple(r) for r in R if r)


def bump_routes(w):
    P, rows = (), []
    for x in w:
        P, route, _ = row_insert(P, x)
        rows.append(route)
    Pc, cols = (), []
    for x in reversed(w):
        Pc, route, _ = column_insert(Pc, x)
        cols.append(route)
    return P, Pc, rows, cols


def shadow_lines(w):
    """回傳折線串列，每條是 (列, 行) 座標的陣列；邊界用 0.5 與 n+0.5。"""
    n = len(w)
    num = []
    for i in range(n):
        num.append(1 + max([num[j] for j in range(i) if w[j] < w[i]], default=0))
    lines = []
    for k in range(1, max(num) + 1):
        B = sorted((i + 1, w[i]) for i in range(n) if num[i] == k)
        pts = [(B[0][0], n + 0.5)]
        for a, (r, c) in enumerate(B):
            pts += [(r, c), (B[a + 1][0] if a + 1 < len(B) else n + 0.5, c)]
        lines.append(np.array(pts, float))
    return lines


def four_corners(w):
    n = len(w)
    out = {}
    for name, fr, fc in (("nw", 0, 0), ("ne", 0, 1), ("sw", 1, 0), ("se", 1, 1)):
        v = list(w[::-1]) if fr else list(w)
        v = [n + 1 - x for x in v] if fc else v
        out[name] = [np.c_[n + 1 - P[:, 0] if fr else P[:, 0], n + 1 - P[:, 1] if fc else P[:, 1]] for P in shadow_lines(v)]
    return out


def rsk_rows(w):
    P, Q = [], []
    for t, x in enumerate(w, 1):
        r = 0
        while True:
            if r == len(P):
                P.append([x])
                Q.append([t])
                break
            j = bisect.bisect_left(P[r], x)
            if j == len(P[r]):
                P[r].append(x)
                Q[r].append(t)
                break
            P[r][j], x = x, P[r][j]
            r += 1
    return P, Q


# ---------------------------------------------------------------- 幾何
def spectral_xy(n, edges):
    A = np.zeros((n, n))
    for a, b in edges:
        A[a, b] = A[b, a] = 1
    d = A.sum(1)
    L = np.eye(n) - A / np.sqrt(np.outer(d, d))
    _, V = np.linalg.eigh(L)
    return V[:, 1:3] / np.sqrt(d)[:, None]


def pipe_route(D, m, i, cs):
    D = set(D)
    pts = [np.array([0.0, -(i - 0.5) * cs])]
    r, c, d = i, 1, "R"
    while True:
        x0, y0 = (c - 1) * cs, -(r - 1) * cs
        cross = (r, c) in D and r + c <= m
        if d == "R" and not cross:          # 從左進、往上出：圓心在左上角
            for a in np.linspace(-np.pi / 2, 0, 7)[1:]:
                pts.append(np.array([x0, y0]) + cs / 2 * np.array([np.cos(a), np.sin(a)]))
            d = "U"
        elif d == "U" and not cross:        # 從下進、往右出：圓心在右下角
            for a in np.linspace(np.pi, np.pi / 2, 7)[1:]:
                pts.append(np.array([x0 + cs, y0 - cs]) + cs / 2 * np.array([np.cos(a), np.sin(a)]))
            d = "R"
        else:
            pts.append(np.array([x0 + cs, y0 - cs / 2]) if d == "R" else np.array([x0 + cs / 2, y0]))
        if d == "U":
            if r == 1:
                return np.array(pts)
            r -= 1
        else:
            c += 1


def ruling(fam, th, t):
    s = 1 if fam == "A" else -1
    return np.array([np.cos(th) - s * t * np.sin(th), np.sin(th) + s * t * np.cos(th), t])


def line_meets_hyperboloid(p0, d):
    p0, d = np.asarray(p0, float), np.asarray(d, float)
    a = d[0] ** 2 + d[1] ** 2 - d[2] ** 2
    b = 2 * (p0[0] * d[0] + p0[1] * d[1] - p0[2] * d[2])
    c = p0[0] ** 2 + p0[1] ** 2 - p0[2] ** 2 - 1
    disc = b * b - 4 * a * c
    if disc < 0:
        return [], []
    pts = [p0 + s * d for s in ((-b - np.sqrt(disc)) / (2 * a), (-b + np.sqrt(disc)) / (2 * a))]
    th = [np.arctan2((y + z * x) / (1 + z * z), (x - z * y) / (1 + z * z)) for x, y, z in pts]
    return pts, th


def ellipse_curve(a, b, phi, c=(0.0, 0.0)):
    c = np.array(c, float)
    R = np.array([[np.cos(phi), -np.sin(phi)], [np.sin(phi), np.cos(phi)]])
    return (lambda t: c + np.c_[a * np.cos(t), b * np.sin(t)] @ R.T,
            lambda p: (((p - c) @ R)[:, 0] / a) ** 2 + (((p - c) @ R)[:, 1] / b) ** 2 - 1, (0, TAU))


def line_curve(p0, ang):
    p0, d = np.array(p0, float), np.array([np.cos(ang), np.sin(ang)])
    return lambda t: p0 + np.outer(t, d), lambda p: (p - p0) @ np.array([-d[1], d[0]]), (-12, 12)


def cubic_curve(c=(0.0, 0.0), s=1.6, k=0.9):
    f = lambda x: c[1] + k * (((x - c[0]) / s) ** 3 - 3 * (x - c[0]) / s)
    return lambda t: np.c_[t, f(t)], lambda p: p[:, 1] - f(p[:, 0]), (-9, 9)


def meet(A, B, N=60000):
    par, F = A[0], B[1]
    t = np.linspace(*A[2], N)
    g = F(par(t))
    out = []
    for i in np.nonzero(np.sign(g[:-1]) * np.sign(g[1:]) < 0)[0]:
        lo, hi, s0 = t[i], t[i + 1], np.sign(g[i])
        for _ in range(50):
            mid = (lo + hi) / 2
            if np.sign(F(par(np.array([mid])))[0]) == s0:
                lo = mid
            else:
                hi = mid
        out.append(par(np.array([lo]))[0])
    return out


# ---------------------------------------------------------------- 自我測試
def check():
    from itertools import combinations
    rg = np.random.default_rng(0)
    # patience sorting 的堆數 = 最長遞增子序列長 = P 的第一列長
    for _ in range(40):
        w = tuple(int(x) + 1 for x in rg.permutation(12))
        pile, lis = patience(w)
        assert max(pile) + 1 == len(lis) == len(rsk(w)[0][0]) and all(w[a] < w[b] for a, b in zip(lis, lis[1:]))
    # evacuation（書 p.184）：每一步的滑動路徑走到 evacuation_steps 記下的格子
    cur = ((1, 2, 2, 3), (2, 3, 5, 5), (4, 4, 6), (5, 6))
    for nxt, b, _ in evacuation_steps(cur, 6):
        cells, nx = slide_path(cur)
        assert nx == nxt and cells[-1] == b
        cur = nxt
    # 列插入與倒序行插入得到同一個 P；四個角的影子線條數 = 第一列長／第一行長（Greene）
    for _ in range(10):
        w = tuple(int(x) + 1 for x in rg.permutation(15))
        P, Pc, rows, cols = bump_routes(w)
        lam = shape(rsk(w)[0])
        assert P == Pc == rsk(w)[0] and len(rows) == len(cols) == 15
        fc = four_corners(w)
        assert len(fc["nw"]) == len(fc["se"]) == lam[0] and len(fc["ne"]) == len(fc["sw"]) == len(lam)
        P2, Q2 = rsk_rows(w)
        assert (tuple(map(tuple, P2)), tuple(map(tuple, Q2))) == rsk(w)
    # pipe dream：水管的終點行 = pipe_exits
    for D in pipe_dreams((2, 1, 5, 4, 3)):
        ex = pipe_exits(D, 5)
        for i in range(1, 6):
            assert abs(pipe_route(D, 5, i, 1.0)[-1][0] - (ex[i - 1] - 0.5)) < 1e-9
    # 雙曲面：兩族直線在曲面上；過交點的 B 族直線與每條 A 族直線相交
    for th in np.linspace(0, TAU, 7):
        for fam in "AB":
            x, y, z = ruling(fam, th, 0.7)
            assert abs(x * x + y * y - z * z - 1) < 1e-12
    pts, th = line_meets_hyperboloid((-2.2, 0.35, -0.75), (1.0, 0.12, 0.42))
    assert len(pts) == 2
    for q, t0 in zip(pts, th):
        assert np.allclose(ruling("B", t0, q[2]), q)
        for tA in (0.4, 2.3, 4.3):
            M = np.c_[ruling("A", tA, 1) - ruling("A", tA, 0), -(ruling("B", t0, 1) - ruling("B", t0, 0))]
            sol, *_ = np.linalg.lstsq(M, ruling("B", t0, 0) - ruling("A", tA, 0), rcond=None)
            assert np.linalg.norm(ruling("A", tA, sol[0]) - ruling("B", t0, sol[1])) < 1e-9
    # Bézout：直線、兩個橢圓、三次曲線兩兩的實交點數 = 次數乘積
    c = (-0.6, 0.1)
    cv = [(line_curve((-0.6, 0.45), 0.33), 1), (ellipse_curve(3.2, 1.9, 0.0, c), 2), (ellipse_curve(4.0, 1.85, 0.5, c), 2),
          (cubic_curve(c), 3)]
    for (A, da), (B, db) in combinations(cv, 2):
        assert len(meet(A, B)) == da * db
    # 譜嵌入：路徑圖的第 2 特徵向量是單調的
    xy = spectral_xy(8, [(k, k + 1) for k in range(7)])
    assert np.all(np.diff(xy[:, 0]) > 0) or np.all(np.diff(xy[:, 0]) < 0)
    return True


if __name__ == "__main__":
    check()
    print("young_extras（youngart 桌布）：patience、滑動路徑、擠壓路徑、影子線、快速 R-S-K、pipe 水管、雙曲面直線、Bézout 交點全部驗證通過")
