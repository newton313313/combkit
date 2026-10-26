# 由 combkit/tools/sync_from_walk017.py 從 walk017 專案的 codes/utils/walk_rec.py 複製；請改 walk017 那份再同步。
"""
walk_rec.py —— 《Elements of the Random Walk》第 3 章 Generating functions II（回家、維度、走過幾個點）的數學
（純 Python，不 import manim；只依賴同資料夾的 walk_gf，以便原樣複製進 combkit）

記號照書：從原點出發；u_N = C(N; 0)/Z^N 是「第 N 步剛好在原點」的機率（Z = 每步的方向數）；
U(s) = Σ u_N s^N = G(z; 0)（s = z/z_c，z_c = 1/Z）；「到第 N 步都沒回過原點」的機率
P_N = [s^N] 1/((1−s)U(s))（(3.16)、(3.25)、(3.33)）；到訪次數的生成函數 (3.1)–(3.10)；
新點的機率 Δ_n = P_n（(3.43)–(3.48)）；平均足跡 δ_n = Σ_{n'≤n} Δ_{n'}（(3.42)）。
腳印算第 0..N 步（共 N+1 個；書上 (3.3)(3.5) 的乘積應從 i = 0 起）。書上式號與筆誤見 ref/md/ch03_formulas.md。
kind（格子）：'1d'（Z = 2）、'sq'（正方，4）、'sc'（簡單立方，6）、'bcc'（8）、'fcc'（12）；前四種是二部格子（奇數步回不了家）。

§3.2–3.3 到訪次數的生成函數 (3.1)–(3.10)，Fig. 3.1
    walks_1d(N, x=0)                全部 2^N 條一維走法（位置 list，長 N+1）
    visit_counts_1d(N, w, x=0)      {M: 走法數}：到訪 w 恰好 M 次（終點不限；暴力列舉）
    visit_counts_gf_1d(N, w, x=0)   同上，用 (3.9)–(3.10)：M = 0 是 E − G_xw·E/G，M ≥ 1 是 G_xw·E·(G−1)^{M−1}/G^{M+1}
                                    （G = G(z;0)、G_xw = G(z;w−x)、E = Σ_y G(z;w,y) = 1/(1−2z)；第二種算法）
    visit_poly(counts)              (3.1)：Σ_M C(M)(1−t)^M 展開成 t 的多項式（整數 list，t^0 起）
    counts_from_poly(p)             (3.2)：C(M) = (−1)^M (1/M!) d^M/dt^M |_{t=1}
    virial_terms_1d(N, w, x=0)      Fig. 3.1／(3.8) 的第 k 項（終點不限，不帶 (−t)^k 的號）= 有序 k 次到訪的總數 Σ_M C(M)·C(M,k)
§3.4 一維 (3.11)–(3.24)，Ex. 3.2–3.4
    never_return_1d(N)              N 步都不回原點的走法數 2·C(N−1, ⌊(N−1)/2⌋)：1, 2, 2, 4, 6, 12, 20, 40, 70, 140, 252
    never_return_1d_brute(N)        暴力列舉；never_return_1d_gf(nmax) 是 (3.13)(3.22) 的 [z^N] g/G（g = 1/(1−2z)、G = 1/√(1−4z²)）
    escape_1d(N)                    機率（Fraction）= C(2k,k)/4^k，k = ⌊N/2⌋；escape_1d_asym(N) = √(2/(πN))（更正後的 (3.24)）；
                                    escape_1d_book(N) = 1/√(πN)（書上印的，少了 z_c^{−1/2} = √2）
    escape_biased(p)                偏硬幣：1/G(z=1; 0) = √(1−4pq) = |p−q|；escape_biased_series(p, N) 有限 N；
                                    escape_biased_dp(p, N) 逐步遞迴（第二種算法）
§3.5 各維度的「還沒回家」(3.25)–(3.30)
    returns_count(kind, N)          C(N; 0)：'sc' 用 Σ_m C(N,m)·C(m,m/2)²·C(N−m,(N−m)/2)（xy 平面走 m 步）；
                                    'fcc' 用 Σ N!/(a!b!c!)·CT(a+c)CT(a+b)CT(b+c)（CT(m) = C(m, m/2)）
    return_probs(kind, nmax)        u_0..u_nmax（float；'sc' 用對數權重的浮點版，n ≤ 60 與整數版差 < 1e−12）
    first_return_probs(kind, nmax)  f_0..f_nmax：第一次回家在第 n 步（更新方程 u_n = Σ_{k≤n} f_k u_{n−k}）
    no_return_probs(kind, nmax)     P_0..P_nmax（float tuple）：P_n = 1 − Σ_{k≤n} f_k
    return_count_dist(kind, N, mmax)   N 步內剛好回家 m 次的機率（m = 0..mmax；有限 N 的 f_m）
    never_return_dp(kind, N)        N 步都不回原點的走法數（逐步遞迴、每步把原點清零；整數，第二種算法）
    no_return_sq_book(N) = π/ln N（書上的首項）；no_return_sq_refined(N) = π/(ln N + ln 8 + γ)（含常數項，U(s) ≈ (ln 8 − ln(1−s))/π）
    log10_steps_for(p)              π/ln N = p 時 N 的 log10（p = 0.01 → 136.4）
§3.6 維度的分水嶺 (3.29)–(3.33)，Ex. 3.5–3.6
    K_d(d)、K_d_gauss(d)            (3.29) 2π^{d/2}/Γ(d/2)；Ex. 3.5 用 ∫e^{−r²}d^dr = π^{d/2} 的數值積分核對
    tail_integral(d, eps)           ∫_eps^1 q^{d−1}/q² dq（閉式）；tail_integral_numeric(d, eps) 數值積分
    G0(kind, x)                     G(z; 0)，z = x·z_c：'1d' 1/√(1−x²)、'sq' (2/π)K(x)、'sc' ∫_0^∞ e^{−t}I_0(xt/3)³dt
    G0_series(kind, x, nmax)        Σ u_n x^n（第二種算法，x < 1）；i0e(y) = e^{−y}I_0(y)
    B1_SC = 9/π                     sc 的 G(z;0) ≈ G_c − B1·√(z_c − z)（(3.32)，c1 = 6、c2 = 1）
    sc_correction_coeff()           Ex. 3.6：P_N ≈ 1/G_c + A/√N，A = B1·√z_c/(G_c²√π) = 0.28696
§3.7 Watson 積分 (3.34)–(3.41)，Ex. 3.7–3.8
    watson_gamma(kind)              Gamma 閉式：sc √6/(32π³)Γ(1/24)Γ(5/24)Γ(7/24)Γ(11/24)、bcc Γ(1/4)⁴/(4π³)、fcc 9Γ(1/3)⁶/(2^{14/3}π⁴)
    watson_sc_bessel()              ∫_0^∞ e^{−t} I_0(t/3)³ dt（數值積分＋尾巴的漸近展開）
    watson_bcc_series(K=200)        Σ (C(2k,k)/4^k)³（＋Hurwitz ζ 的尾巴；hurwitz_zeta(s, a)）
    watson_fcc_series(nmax=60)      Σ u_n（精確的 fcc 回家次數）＋高斯尾巴（粗核對，誤差 < 2e−3）
    G_C                             {'sc': 1.516386…, 'bcc': 1.393204…, 'fcc': 1.344661…}；escape(kind) = 1/G_c
    f_m(kind, m)                    剛好回家 m 次的機率 R^m(1−R)（R = 1 − 1/G_c）；Σ f_m = 1、Σ m·f_m = G_c − 1
§3.8 走過幾個點 (3.42)–(3.58)，Ex. 3.9–3.12
    new_site_probs_brute(kind, n)   Δ_0..Δ_n（暴力列舉，Fraction）；distinct_sites_brute(kind, n) 平均足跡（Fraction）
    distinct_sites_1d(n)            一維精確 δ_n（Fraction）：(2n+1)c（n 偶）、2n·c（n 奇），c = C(2k,k)/4^k（第三種算法）
    distinct_sites(kind, nmax)      δ_0..δ_nmax = Σ Δ（float tuple）
    footprint_1d_asym(n) = √(8n/π)；footprint_2d_book(n) = πn/ln n；footprint_2d_mw(n) = πn/ln(8n)；三維 n/G_c
模擬（固定種子；格點編成整數，不畫圖）
    simulate_first_returns(kind, M, N, seed)   M 個走路人各走 N 步：第一次回家的步數（沒回家為 None）
    simulate_return_counts(kind, M, N, seed)   N 步內回家幾次（list）
    simulate_distinct(kind, M, n, seed)        n 步的足跡大小（list）
    trail(kind, n, seed)                       一條足跡（格點 tuple 的 list，長 n+1）
    SIM_*、TRAIL_*                             check() 鎖定的模擬結果（場景直接引用；見檔尾）

check() 驗證全部（書上的數字、筆誤的更正、plan_ch03.md 場景卡上的每個數字）；直接執行本檔會跑 check()。
"""
from __future__ import annotations

import math
import random
from fractions import Fraction
from functools import lru_cache
from itertools import product
from operator import mul

if __package__:
    from . import walk_gf as GF
else:                                   # python utils/walk_rec.py
    import walk_gf as GF

EULER_GAMMA = 0.5772156649015329
Z = {"1d": 2, "sq": 4, "sc": 6, "bcc": 8, "fcc": 12}
BIPARTITE = {"1d", "sq", "sc", "bcc"}
STEPS = GF.STEPS


def CT(m: int) -> int:
    """一維 m 步回到原點的走法數 C(m, m/2)（m 奇數為 0）。"""
    return math.comb(m, m // 2) if m % 2 == 0 else 0


# =========================================================================== 冪級數（list，z^0 起）


def _ser_mul(a, b, n: int) -> list:
    out = [0] * (n + 1)
    for i, ai in enumerate(a[: n + 1]):
        if ai:
            for j, bj in enumerate(b[: n + 1 - i]):
                out[i + j] += ai * bj
    return out


def _ser_inv(a, n: int) -> list:
    """1/a；a[0] = 1 時整數係數保持整數。"""
    assert a[0] == 1
    out = [0] * (n + 1)
    out[0] = 1
    for k in range(1, n + 1):
        out[k] = -sum(a[j] * out[k - j] for j in range(1, min(k, len(a) - 1) + 1))
    return out


def _simpson(f, a: float, b: float, n: int = 2000) -> float:
    n += n % 2
    h = (b - a) / n
    s = f(a) + f(b)
    s += 4 * sum(f(a + (2 * i - 1) * h) for i in range(1, n // 2 + 1))
    s += 2 * sum(f(a + 2 * i * h) for i in range(1, n // 2))
    return s * h / 3


# =========================================================================== §3.2–3.3 到訪次數


def walks_1d(N: int, x: int = 0) -> list[list[int]]:
    out = []
    for steps in product((1, -1), repeat=N):
        pos = [x]
        for s in steps:
            pos.append(pos[-1] + s)
        out.append(pos)
    return out


def visit_counts_1d(N: int, w: int, x: int = 0) -> dict:
    cnt = {}
    for pos in walks_1d(N, x):
        M = pos.count(w)
        cnt[M] = cnt.get(M, 0) + 1
    return dict(sorted(cnt.items()))


def _G1d(d: int, n: int) -> list[int]:
    """G(z; d) 的係數 C1(k, d)，k = 0..n。"""
    return [GF.count1(k, d) for k in range(n + 1)]


def visit_counts_gf_1d(N: int, w: int, x: int = 0) -> dict:
    G = _G1d(0, N)
    Ginv = _ser_inv(G, N)
    E = [2 ** k for k in range(N + 1)]
    A = _ser_mul(_G1d(w - x, N), E, N)                 # G(z; x, w)·Σ_y G(z; w, y)
    out = {}
    m0 = [e - c for e, c in zip(E, _ser_mul(A, Ginv, N))]
    if m0[N]:
        out[0] = m0[N]
    Gm1 = [g - (k == 0) for k, g in enumerate(G)]      # G_1 = G − 1
    term = _ser_mul(A, _ser_mul(Ginv, Ginv, N), N)     # M = 1：A/G²；之後每次乘 (G−1)/G
    M = 1
    while any(term):
        if term[N]:
            out[M] = term[N]
        term = _ser_mul(_ser_mul(term, Gm1, N), Ginv, N)
        M += 1
    return dict(sorted(out.items()))


def visit_poly(counts: dict) -> list[int]:
    p = [0] * (max(counts) + 1)
    for M, c in counts.items():
        for k in range(M + 1):
            p[k] += c * math.comb(M, k) * (-1) ** k
    return p


def counts_from_poly(p) -> dict:
    out = {}
    for M in range(len(p)):
        d = sum(p[k] * math.perm(k, M) for k in range(M, len(p)))       # d^M/dt^M Σ p_k t^k 在 t = 1
        v = Fraction((-1) ** M * d, math.factorial(M))
        assert v.denominator == 1
        if v:
            out[M] = int(v)
    return out


def virial_terms_1d(N: int, w: int, x: int = 0) -> list[int]:
    E = [2 ** k for k in range(N + 1)]
    B = [0] + _G1d(0, N)[1:]                            # N2..Nk ≥ 1（至少走一步才再到訪）
    terms = [2 ** N]
    cur = _ser_mul(_G1d(w - x, N), E, N)
    while any(cur):
        terms.append(cur[N])
        cur = _ser_mul(cur, B, N)
    while len(terms) > 1 and terms[-1] == 0:
        terms.pop()
    return terms


# =========================================================================== §3.4 一維


def never_return_1d(N: int) -> int:
    return 1 if N == 0 else 2 * math.comb(N - 1, (N - 1) // 2)


def never_return_1d_brute(N: int) -> int:
    return sum(1 for pos in walks_1d(N) if 0 not in pos[1:])


def never_return_1d_gf(nmax: int) -> list[int]:
    return _ser_mul([2 ** k for k in range(nmax + 1)], _ser_inv(_G1d(0, nmax), nmax), nmax)


def escape_1d(N: int) -> Fraction:
    return Fraction(never_return_1d(N), 2 ** N)


def escape_1d_asym(N: float) -> float:
    return math.sqrt(2 / (math.pi * N))


def escape_1d_book(N: float) -> float:
    return 1 / math.sqrt(math.pi * N)


def escape_biased(p: float) -> float:
    return math.sqrt(1 - 4 * p * (1 - p))


def escape_biased_series(p: float, N: int) -> float:
    """[z^N] (1/(1−z))/G_p(z; 0)，G_p = Σ C(2k,k)(pq)^k z^{2k}（機率的生成函數，z_c = 1）。"""
    pq = p * (1 - p)
    G = [0.0] * (N + 1)
    c = 1.0
    for k in range(N // 2 + 1):
        G[2 * k] = c
        c *= (2 * k + 1) * (2 * k + 2) / (k + 1) ** 2 * pq
    inv = [1.0] + [0.0] * N
    for k in range(1, N + 1):
        inv[k] = -sum(G[j] * inv[k - j] for j in range(2, k + 1, 2))
    return sum(inv)


def escape_biased_dp(p: float, N: int) -> float:
    cur = {0: 1.0}
    for _ in range(N):
        nxt = {}
        for x, w in cur.items():
            nxt[x + 1] = nxt.get(x + 1, 0.0) + p * w
            nxt[x - 1] = nxt.get(x - 1, 0.0) + (1 - p) * w
        nxt.pop(0, None)
        cur = nxt
    return sum(cur.values())


# =========================================================================== §3.5 還沒回家


@lru_cache(maxsize=None)
def returns_count(kind: str, N: int) -> int:
    if kind == "1d":
        return CT(N)
    if kind == "sq":
        return CT(N) ** 2
    if kind == "bcc":
        return CT(N) ** 3
    if kind == "sc":
        if N % 2:
            return 0
        return sum(math.comb(N, m) * CT(m) ** 2 * CT(N - m) for m in range(0, N + 1, 2))
    if kind == "fcc":
        f = math.factorial
        tot = 0
        for a in range(N + 1):
            for b in range(N + 1 - a):
                c = N - a - b
                t = CT(a + c) * CT(a + b) * CT(b + c)
                if t:
                    tot += f(N) // (f(a) * f(b) * f(c)) * t
        return tot
    raise ValueError(kind)


@lru_cache(maxsize=None)
def _return_probs(kind: str, nmax: int) -> tuple:
    u = [0.0] * (nmax + 1)
    if kind in ("1d", "sq", "bcc"):
        dim = {"1d": 1, "sq": 2, "bcc": 3}[kind]
        c = 1.0
        for k in range(nmax // 2 + 1):
            u[2 * k] = c ** dim
            c *= (2 * k + 1) / (2 * k + 2)
    elif kind == "sc":          # Σ_m C(n,m)(2/3)^m(1/3)^{n−m}·u_sq(m)·u_1d(n−m)（權重用對數，只加 |m − 2n/3| ≤ 12σ+6 的項）
        u2, u1 = _return_probs("sq", nmax), _return_probs("1d", nmax)
        lg = [math.lgamma(k + 1) for k in range(nmax + 1)]
        a, b = math.log(2 / 3), math.log(1 / 3)
        for n in range(0, nmax + 1, 2):
            h = 12 * math.sqrt(2 * n / 9) + 6
            lo, hi = max(0, int(2 * n / 3 - h)), min(n, int(2 * n / 3 + h))
            u[n] = sum(math.exp(lg[n] - lg[m] - lg[n - m] + m * a + (n - m) * b) * u2[m] * u1[n - m]
                       for m in range(lo - lo % 2, hi + 1, 2))
    else:
        for n in range(nmax + 1):
            u[n] = returns_count(kind, n) / Z[kind] ** n
    return tuple(u)


def return_probs(kind: str, nmax: int) -> list[float]:
    return list(_return_probs(kind, nmax))


@lru_cache(maxsize=None)
def first_return_probs(kind: str, nmax: int) -> tuple:
    """f_0..f_nmax：第一次回家在第 n 步的機率（更新方程 u_n = Σ_{k≤n} f_k u_{n−k}，f_0 = 0）。"""
    u = _return_probs(kind, nmax)
    step = 2 if kind in BIPARTITE else 1
    v = u[::step]
    K = len(v) - 1
    f = [0.0] * (K + 1)
    for k in range(1, K + 1):
        f[k] = v[k] - sum(map(mul, f[1:k], v[k - 1:0:-1]))
    out = [0.0] * (nmax + 1)
    for k in range(K + 1):
        out[k * step] = f[k]
    return tuple(out)


@lru_cache(maxsize=None)
def no_return_probs(kind: str, nmax: int) -> tuple:
    out, acc = [], 1.0
    for fn in first_return_probs(kind, nmax):
        acc -= fn
        out.append(acc)
    return tuple(out)


def return_count_dist(kind: str, N: int, mmax: int) -> list[float]:
    """N 步內剛好回家 m 次的機率（m = 0..mmax）：P(≥ m) = Σ_{t≤N} f^{*m}(t)。"""
    f = first_return_probs(kind, N)
    conv, ge = [1.0] + [0.0] * N, []
    for _ in range(mmax + 2):
        ge.append(sum(conv))
        conv = [sum(map(mul, conv[:n], f[n:0:-1])) for n in range(N + 1)]
    return [ge[m] - ge[m + 1] for m in range(mmax + 1)]


def never_return_dp(kind: str, N: int) -> int:
    steps = STEPS[kind]
    origin = (0,) * len(steps[0])
    cur = {origin: 1}
    for _ in range(N):
        nxt = {}
        for p, c in cur.items():
            for s in steps:
                q = tuple(a + b for a, b in zip(p, s))
                nxt[q] = nxt.get(q, 0) + c
        nxt.pop(origin, None)
        cur = nxt
    return sum(cur.values())


def no_return_sq_book(N: float) -> float:
    return math.pi / math.log(N)


def no_return_sq_refined(N: float) -> float:
    return math.pi / (math.log(N) + math.log(8) + EULER_GAMMA)


def log10_steps_for(p: float) -> float:
    return math.pi / p / math.log(10)


# =========================================================================== §3.6 維度


def K_d(d: int) -> float:
    return 2 * math.pi ** (d / 2) / math.gamma(d / 2)


def K_d_gauss(d: int) -> float:
    return math.pi ** (d / 2) / _simpson(lambda r: r ** (d - 1) * math.exp(-r * r), 0.0, 12.0, 6000)


def tail_integral(d: int, eps: float) -> float:
    if d == 2:
        return -math.log(eps)
    return (1 - eps ** (d - 2)) / (d - 2)


def tail_integral_numeric(d: int, eps: float) -> float:
    return _simpson(lambda v: math.exp((d - 2) * v), math.log(eps), 0.0, 4000)      # q = e^v


_I0_ASYM = [1.0, 1 / 8, 9 / 128, 225 / 3072, 11025 / 98304, 893025 / 3932160, 108056025 / 188743680]


def i0e(y: float) -> float:
    """e^{−y} I_0(y)（y ≥ 0）：小 y 用級數、大 y 用漸近展開。"""
    if y < 40:
        s, t, k, q = 1.0, 1.0, 0, y * y / 4
        while True:
            k += 1
            t *= q / (k * k)
            s += t
            if t < 1e-17 * s:
                return s * math.exp(-y)
    return sum(c / y ** k for k, c in enumerate(_I0_ASYM)) / math.sqrt(2 * math.pi * y)


def _cube_coeffs(n: int = 5) -> list[float]:
    a = _I0_ASYM[: n + 1]
    return _ser_mul(_ser_mul(a, a, n), a, n)


def G0(kind: str, x: float) -> float:
    if kind == "1d":
        return 1 / math.sqrt(1 - x * x)
    if kind == "sq":
        return 2 / math.pi * GF.ellipk(x)
    if kind != "sc":
        raise ValueError(kind)
    # ∫_0^∞ e^{−t(1−x)} [e^{−xt/3} I_0(xt/3)]³ dt；t = e^v − 1
    T = 3000.0 if x == 1 else 40 / (1 - x)
    V = math.log1p(T)
    val = _simpson(lambda v: math.exp(v) * math.exp(-(math.exp(v) - 1) * (1 - x)) * i0e(x * (math.exp(v) - 1) / 3) ** 3,
                   0.0, V, 6000)
    if x == 1:                                        # 尾巴：3(2π)^{−3/2} Σ b_k Y^{−k−1/2}/(k+1/2)，Y = T/3
        Y = T / 3
        b = _cube_coeffs()
        val += 3 * (2 * math.pi) ** -1.5 * sum(bk * Y ** (-k - 0.5) / (k + 0.5) for k, bk in enumerate(b))
    return val


def G0_series(kind: str, x: float, nmax: int = 4000) -> float:
    u = return_probs(kind, nmax)
    return sum(un * x ** n for n, un in enumerate(u))


B1_SC = 9 / math.pi


def sc_correction_coeff() -> float:
    return B1_SC * math.sqrt(1 / 6) / (G_C["sc"] ** 2 * math.sqrt(math.pi))


# =========================================================================== §3.7 Watson 積分


def watson_gamma(kind: str) -> float:
    g, pi = math.gamma, math.pi
    if kind == "sc":
        return math.sqrt(6) / (32 * pi ** 3) * g(1 / 24) * g(5 / 24) * g(7 / 24) * g(11 / 24)
    if kind == "bcc":
        return g(0.25) ** 4 / (4 * pi ** 3)
    if kind == "fcc":
        return 9 * g(1 / 3) ** 6 / (2 ** (14 / 3) * pi ** 4)
    raise ValueError(kind)


def watson_sc_bessel() -> float:
    return G0("sc", 1.0)


def hurwitz_zeta(s: float, a: float, M: int = 12) -> float:
    """ζ(s, a) = Σ_{j≥0} (a+j)^{−s}（Euler–Maclaurin）。"""
    tot = sum((a + j) ** -s for j in range(M))
    b = a + M
    tot += b ** (1 - s) / (s - 1) + 0.5 * b ** -s
    fact = s * b ** (-s - 1)
    for B2k, k in ((1 / 6, 1), (-1 / 30, 2), (1 / 42, 3), (-1 / 30, 4)):
        tot += B2k / math.factorial(2 * k) * fact
        fact *= (s + 2 * k - 1) * (s + 2 * k) / (b * b)
    return tot


def watson_bcc_series(K: int = 200) -> float:
    """Σ_{k≤K} c_k³ + 尾巴；c_k = C(2k,k)/4^k ≈ (πk)^{−1/2}(1 − 1/8k + 1/128k² + 5/1024k³)，
    c_k³ ≈ (πk)^{−3/2}(1 − 3/8k + 9/128k² + 7/1024k³)。"""
    s, c = 0.0, 1.0
    for k in range(K + 1):
        s += c ** 3
        c *= (2 * k + 1) / (2 * k + 2)
    a = K + 1
    tail = (hurwitz_zeta(1.5, a) - 3 / 8 * hurwitz_zeta(2.5, a) + 9 / 128 * hurwitz_zeta(3.5, a)
            + 7 / 1024 * hurwitz_zeta(4.5, a))
    return s + tail * math.pi ** -1.5


def watson_fcc_series(nmax: int = 60) -> float:
    """Σ_{n≤nmax} u_n + (1/√2)(3/2π)^{3/2} ζ(3/2, nmax+1)（高斯尾巴；fcc 單胞體積 2、每座標變異數 2/3）。"""
    s = sum(return_probs("fcc", nmax))
    return s + 2 ** -0.5 * (3 / (2 * math.pi)) ** 1.5 * hurwitz_zeta(1.5, nmax + 1)


G_C = {k: watson_gamma(k) for k in ("sc", "bcc", "fcc")}


def escape(kind: str) -> float:
    return 1 / G_C[kind]


def f_m(kind: str, m: int) -> float:
    R = 1 - 1 / G_C[kind]
    return R ** m * (1 - R)


# =========================================================================== §3.8 走過幾個點


def _all_walks(kind: str, n: int):
    steps = STEPS[kind]
    origin = (0,) * len(steps[0])
    for seq in product(steps, repeat=n):
        pos = [origin]
        for s in seq:
            pos.append(tuple(a + b for a, b in zip(pos[-1], s)))
        yield pos


def new_site_probs_brute(kind: str, n: int) -> list[Fraction]:
    tot = [0] * (n + 1)
    for pos in _all_walks(kind, n):
        seen = set()
        for k, p in enumerate(pos):
            if p not in seen:
                seen.add(p)
                tot[k] += 1
    return [Fraction(t, Z[kind] ** n) for t in tot]


def distinct_sites_brute(kind: str, n: int) -> Fraction:
    return Fraction(sum(len(set(pos)) for pos in _all_walks(kind, n)), Z[kind] ** n)


def distinct_sites_1d(n: int) -> Fraction:
    k = n // 2
    c = Fraction(math.comb(2 * k, k), 4 ** k)
    return (2 * n + 1) * c if n % 2 == 0 else 2 * n * c


def distinct_sites(kind: str, nmax: int) -> tuple:
    out, acc = [], 0.0
    for P in no_return_probs(kind, nmax):
        acc += P
        out.append(acc)
    return tuple(out)


def footprint_1d_asym(n: float) -> float:
    return math.sqrt(8 * n / math.pi)


def footprint_2d_book(n: float) -> float:
    return math.pi * n / math.log(n)


def footprint_2d_mw(n: float) -> float:
    return math.pi * n / math.log(8 * n)


# =========================================================================== 模擬（固定種子）


def _codes(kind: str, N: int) -> list[int]:
    B = 2 * N + 3
    return [sum(c * B ** i for i, c in enumerate(s)) for s in STEPS[kind]]


def simulate_first_returns(kind: str, M: int, N: int, seed: int) -> list:
    rng = random.Random(seed)
    codes = _codes(kind, N)
    Zk = len(codes)
    out = []
    for _ in range(M):
        p, t = 0, None
        for n in range(1, N + 1):
            p += codes[rng.randrange(Zk)]
            if p == 0:
                t = n
                break
        out.append(t)
    return out


def simulate_return_counts(kind: str, M: int, N: int, seed: int) -> list[int]:
    rng = random.Random(seed)
    codes = _codes(kind, N)
    Zk = len(codes)
    out = []
    for _ in range(M):
        p = r = 0
        for _ in range(N):
            p += codes[rng.randrange(Zk)]
            r += p == 0
        out.append(r)
    return out


def simulate_distinct(kind: str, M: int, n: int, seed: int) -> list[int]:
    rng = random.Random(seed)
    codes = _codes(kind, n)
    Zk = len(codes)
    out = []
    for _ in range(M):
        p, seen = 0, {0}
        for _ in range(n):
            p += codes[rng.randrange(Zk)]
            seen.add(p)
        out.append(len(seen))
    return out


def trail(kind: str, n: int, seed: int) -> list[tuple]:
    rng = random.Random(seed)
    steps = STEPS[kind]
    pos = [(0,) * len(steps[0])]
    for _ in range(n):
        s = steps[rng.randrange(len(steps))]
        pos.append(tuple(a + b for a, b in zip(pos[-1], s)))
    return pos


# =========================================================================== check


def _r(x, k):
    return [round(v, k) for v in x]


def check():
    F = Fraction
    # ---- §3.2–3.3：生成函數 (3.9)(3.10)＝暴力列舉；(3.1)(3.2)；Fig. 3.1／(3.8)
    for N in range(0, 9):
        for w in (-2, 0, 1, 3):
            c = visit_counts_1d(N, w)
            assert c == visit_counts_gf_1d(N, w), (N, w)
            p = visit_poly(c)
            assert counts_from_poly(p) == c
            assert virial_terms_1d(N, w) == [abs(a) for a in p], (N, w, p)
    assert visit_counts_1d(4, 1) == {0: 6, 1: 6, 2: 4}                                       # W302 主例
    assert visit_poly({0: 6, 1: 6, 2: 4}) == [16, -14, 4] and virial_terms_1d(4, 1) == [16, 14, 4]
    assert GF.count1(1, 1) * 2 ** 3 + GF.count1(3, 1) * 2 == 14 == 8 + 6
    assert GF.count1(1, 1) * GF.count1(2, 0) * 2 == 4
    assert visit_counts_1d(4, 0) == {1: 6, 2: 6, 3: 4}                                        # w = 起點：第 0 步也算
    # ---- §3.4 一維：三種算法；(3.23)(3.24) 的更正；偏硬幣
    NR = [1, 2, 2, 4, 6, 12, 20, 40, 70, 140, 252]
    assert [never_return_1d(N) for N in range(11)] == NR == never_return_1d_gf(10)
    assert all(never_return_1d(N) == never_return_1d_brute(N) == never_return_dp("1d", N) for N in range(15))
    assert all(never_return_1d(2 * k) == math.comb(2 * k, k) for k in range(30))               # = 第 2k 步在家的走法數
    assert escape_1d(10) == F(63, 256) == F(252, 1024)
    assert [round(float(escape_1d(N)), 4) for N in (10, 100)] == [0.2461, 0.0796]
    assert round(float(escape_1d(1000)), 6) == 0.025225
    assert _r([escape_1d_asym(N) for N in (10, 100, 1000)], 5) == [0.25231, 0.07979, 0.02523]
    assert _r([escape_1d_book(N) for N in (10, 100, 1000)], 5) == [0.17841, 0.05642, 0.01784]
    assert abs(GF.coeff_power_asym(10, 0.5, 0.5) / 2 ** 10 - escape_1d_asym(10)) < 1e-14        # z_c^{−1/2} = √2 不能少
    assert round(float(escape_1d(1000)) / escape_1d_asym(1000), 5) == 0.99975
    assert abs(escape_biased(0.6) - 0.2) < 1e-12 and escape_biased(1.0) == 1.0              # Ex. 3.2、3.3
    assert round(1 / escape_biased(0.6), 9) == 5 and 2 ** 4 - never_return_1d(4) == 10       # G(1;0) = 5；W301：16 條裡 10 條回過家
    assert sum(len(set(p)) for p in walks_1d(4)) == 54                                        # W307：16 條的足跡總和 = 16·27/8
    PB = [escape_biased_series(0.6, N) for N in (1, 2, 4, 10, 20, 50, 100)]
    assert all(abs(a - escape_biased_dp(0.6, N)) < 1e-12 for a, N in zip(PB, (1, 2, 4, 10, 20, 50, 100)))
    assert _r(PB, 4) == [1.0, 0.52, 0.4048, 0.294, 0.2427, 0.2094, 0.2016]
    assert abs(escape_biased_series(0.6, 400) - 0.2) < 1e-6 and round(1 - escape_biased(0.6), 12) == 0.8   # Ex. 3.4
    # ---- §3.5：回家次數兩種算法；更新方程＝逐步遞迴；二維的表
    assert all(returns_count("sc", n) == GF.returns_sc(n) for n in range(13))
    assert all(returns_count("fcc", n) == GF.returns_fcc(n) for n in range(7))
    assert [returns_count("fcc", n) for n in range(5)] == [1, 0, 12, 48, 540]
    usc = return_probs("sc", 60)
    assert all(abs(usc[n] - returns_count("sc", n) / 6 ** n) <= 1e-12 * usc[n] for n in range(61))
    assert [returns_count("sq", n) for n in range(0, 11, 2)] == [1, 4, 36, 400, 4900, 63504]
    for kind, Ns in (("sq", (2, 4, 10, 30)), ("sc", (2, 4, 10, 16)), ("fcc", (1, 2, 3, 6)), ("bcc", (2, 4, 8))):
        P = no_return_probs(kind, max(Ns))
        assert all(abs(P[N] - never_return_dp(kind, N) / Z[kind] ** N) < 1e-12 for N in Ns), kind
    fs = first_return_probs("sq", 10)                                                         # W304：第一次回家
    assert _r([fs[n] for n in (2, 4, 6, 8, 10)], 4) == [0.25, 0.0781, 0.043, 0.0286, 0.021] and round(sum(fs), 4) == 0.4208
    P2 = no_return_probs("sq", 10000)
    assert _r([P2[N] for N in (10, 100, 1000, 10000)], 4) == [0.5792, 0.4168, 0.3218, 0.2613]
    assert _r([no_return_sq_book(N) for N in (10, 100, 1000, 10000)], 4) == [1.3644, 0.6822, 0.4548, 0.3411]
    assert _r([no_return_sq_refined(N) for N in (10, 100, 1000, 10000)], 4) == [0.6335, 0.4326, 0.3285, 0.2647]
    assert round(log10_steps_for(0.01), 1) == 136.4
    P1, P3 = no_return_probs("1d", 1000), no_return_probs("sc", 1000)
    assert all(abs(P1[N] - float(escape_1d(N))) < 1e-12 for N in range(0, 1001, 37))
    assert _r([P3[N] for N in (10, 100, 1000)], 4) == [0.7472, 0.6881, 0.6685]
    # 模擬（固定種子）：1000 人走 1000 步
    for kind in ("1d", "sq", "sc"):
        fr = simulate_first_returns(kind, SIM_M, SIM_N, SIM_SEED[kind])
        assert sum(t is not None for t in fr) == SIM_RETURNED[kind]
        assert sum(t is not None and t <= 100 for t in fr) == SIM_RETURNED_100[kind]
        assert abs(SIM_RETURNED[kind] / SIM_M - (1 - no_return_probs(kind, 1000)[1000])) < 0.03
    assert _r([1 - no_return_probs(k, 1000)[1000] for k in ("1d", "sq", "sc")], 4) == [0.9748, 0.6782, 0.3315]
    assert _r([1 - no_return_probs(k, 1000)[100] for k in ("1d", "sq", "sc")], 4) == [0.9204, 0.5832, 0.3119]
    # ---- §3.6 維度
    assert _r([K_d(d) for d in (1, 2, 3, 4)], 3) == [2.0, 6.283, 12.566, 19.739]
    assert all(abs(K_d(d) - K_d_gauss(d)) < 1e-9 for d in (1, 2, 3, 4)) and abs(K_d(4) - 2 * math.pi ** 2) < 1e-12
    T = {d: [tail_integral(d, e) for e in (0.1, 0.01, 0.001)] for d in (1, 2, 3)}
    assert _r(T[1], 6) == [9, 99, 999] and _r(T[2], 3) == [2.303, 4.605, 6.908] and _r(T[3], 6) == [0.9, 0.99, 0.999]
    assert all(abs(tail_integral(d, e) - tail_integral_numeric(d, e)) < 1e-9 for d in (1, 2, 3) for e in (0.1, 0.01, 0.001))
    G0T = {k: [G0(k, x) for x in (0.9, 0.99, 0.999)] for k in ("1d", "sq", "sc")}
    assert _r(G0T["1d"], 3) == [2.294, 7.089, 22.366] and _r(G0T["sq"], 3) == [1.452, 2.137, 2.862]
    assert _r(G0T["sc"], 3) == [1.233, 1.409, 1.48]
    for k in ("sq", "sc"):
        assert all(abs(G0(k, x) - G0_series(k, x, 3000)) < 1e-10 for x in (0.9, 0.99)), k
    est = (G_C["sc"] - G0("sc", 0.99999)) / math.sqrt(0.00001 / 6)                         # (3.32)：B1 = 9/π（c2 = 1）
    assert abs(est - B1_SC) < 0.01 and round(B1_SC, 4) == 2.8648
    A = sc_correction_coeff()
    assert round(A, 5) == 0.28696 and abs((P3[1000] - 1 / G_C["sc"]) * math.sqrt(1000) - A) < 2e-4   # Ex. 3.6
    assert _r([1 / G_C["sc"] + A / math.sqrt(N) for N in (100, 1000)], 6) == [0.688159, 0.668537]
    # ---- §3.7 Watson：Gamma 閉式＝Bessel 積分（sc）＝級數（bcc）；fcc 粗核對
    assert abs(watson_sc_bessel() - G_C["sc"]) < 1e-10 and abs(watson_bcc_series() - G_C["bcc"]) < 1e-10
    assert abs(watson_fcc_series() - G_C["fcc"]) < 1e-3
    assert [round(G_C[k], 6) for k in ("sc", "bcc", "fcc")] == [1.516386, 1.393204, 1.344661]
    assert math.floor(G_C["bcc"] * 1e6) / 1e6 == 1.393203                                     # 書上印的截斷值
    assert [round(escape(k), 6) for k in ("sc", "bcc", "fcc")] == [0.659463, 0.71777, 0.743682]
    assert round(1 - escape("sc"), 6) == 0.340537 and round(1 - escape("sc"), 4) == 0.3405
    assert _r([f_m("sc", m) for m in range(5)], 4) == [0.6595, 0.2246, 0.0765, 0.026, 0.0089]
    assert abs(sum(f_m("sc", m) for m in range(400)) - 1) < 1e-12                            # Ex. 3.8
    assert abs(sum(m * f_m("sc", m) for m in range(400)) - (G_C["sc"] - 1)) < 1e-12
    assert round(G_C["sc"] - 1, 6) == 0.516386
    assert _r(return_count_dist("sc", 1000, 4), 4) == [0.6685, 0.2217, 0.0735, 0.0243, 0.008]
    rc = simulate_return_counts("sc", SIM_M, SIM_N, SIM_COUNTS_SEED)
    assert [rc.count(m) for m in range(4)] + [sum(r >= 4 for r in rc)] == SIM_COUNTS_SC
    assert round(sum(rc) / SIM_M, 3) == 0.509 and round(sum(return_probs("sc", 1000)) - 1, 4) == 0.4955
    # ---- §3.8 足跡：暴力列舉＝生成函數＝閉式
    D1 = new_site_probs_brute("1d", 4)
    assert D1 == [1, 1, F(1, 2), F(1, 2), F(3, 8)]
    assert [distinct_sites_brute("1d", n) for n in range(5)] == [1, 2, F(5, 2), 3, F(27, 8)]
    assert all(distinct_sites_1d(n) == distinct_sites_brute("1d", n) for n in range(13))
    assert all(abs(float(a) - b) < 1e-12 for a, b in zip(new_site_probs_brute("1d", 10), no_return_probs("1d", 10)))
    assert new_site_probs_brute("sq", 5) == [1, 1, F(3, 4), F(3, 4), F(43, 64), F(43, 64)]
    assert [distinct_sites_brute("sq", n) for n in range(6)] == [1, 2, F(11, 4), F(7, 2), F(267, 64), F(155, 32)]
    assert all(abs(float(distinct_sites_brute("sq", n)) - distinct_sites("sq", 5)[n]) < 1e-12 for n in range(6))
    assert abs(float(distinct_sites_1d(100)) - distinct_sites("1d", 100)[100]) < 1e-9
    assert round(float(distinct_sites_1d(100)), 3) == 15.997 and round(footprint_1d_asym(100), 3) == 15.958
    D2 = distinct_sites("sq", 10000)
    assert [round(D2[n], 1) for n in (10, 100, 1000, 10000)] == [7.9, 49.5, 363.6, 2874.6]
    assert [round(footprint_2d_book(n), 1) for n in (10, 100, 1000, 10000)] == [13.6, 68.2, 454.8, 3410.9]
    assert [round(footprint_2d_mw(n), 1) for n in (10, 100, 1000, 10000)] == [7.2, 47.0, 349.6, 2782.7]
    assert round(D2[400], 2) == 162.57 and len(set(trail("sq", TRAIL_N, TRAIL_SEED))) == TRAIL_SITES
    assert trail("sq", TRAIL_N, TRAIL_SEED)[-1] == TRAIL_END
    D3 = distinct_sites("sc", 1000)
    assert [round(D3[n] / n, 4) for n in (10, 100, 1000)] == [0.9051, 0.7229, 0.6782]
    assert round(1000 / G_C["sc"] + 2 * A * math.sqrt(1000), 1) == 677.6 and round(D3[1000], 1) == 678.2   # Ex. 3.11
    assert round(2 * A, 4) == 0.5739 and round(1000 / G_C["sc"], 1) == 659.5 and round(2 * A * math.sqrt(1000), 1) == 18.1
    assert round(escape("bcc"), 4) == 0.7178 and round(escape("sc"), 4) == 0.6595                      # Ex. 3.9、3.12
    assert round(sum(escape_biased_series(0.6, N) for N in range(101)), 2) == 24.95                     # Ex. 3.10
    for kind, n, tot in SIM_DISTINCT:
        assert sum(simulate_distinct(kind, 200, n, SIM_DISTINCT_SEED)) == tot
        assert abs(tot / 200 / distinct_sites(kind, n)[n] - 1) < 0.03
    return True


# 由 check() 鎖定的模擬結果（固定種子；改了取樣方法就要重新核對 plan_ch03.md）
SIM_M, SIM_N = 1000, 1000                               # 1000 人各走 1000 步
SIM_SEED = {"1d": 301, "sq": 302, "sc": 303}
SIM_RETURNED = {"1d": 977, "sq": 676, "sc": 325}        # 1000 步內回過家（理論 974.8／678.2／331.5 人）
SIM_RETURNED_100 = {"1d": 930, "sq": 579, "sc": 300}    # 100 步內（理論 920.4／583.2／311.9 人）
SIM_COUNTS_SEED = 306
SIM_COUNTS_SC = [672, 222, 60, 28, 18]                  # sc：1000 步內回家 0、1、2、3、≥4 次的人數
TRAIL_N, TRAIL_SEED, TRAIL_SITES, TRAIL_END = 400, 301, 159, (-8, 14)   # 二維足跡（平均 δ_400 = 162.57）
SIM_DISTINCT_SEED = 308                                 # 200 人的足跡總和（平均 15.665／365.68／679.705）
SIM_DISTINCT = [("1d", 100, 3133), ("sq", 1000, 73136), ("sc", 1000, 135941)]


if __name__ == "__main__":
    import sys
    sys.stdout.reconfigure(encoding="utf-8")
    check()
    print("walk_rec（第 3 章）：到訪次數、一維／二維／三維的回家、Watson 積分、足跡全部驗證通過")
