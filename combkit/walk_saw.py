# 由 combkit/tools/sync_from_walk017.py 從 walk017 專案的 codes/utils/walk_saw.py 複製；請改 walk017 那份再同步。
"""walk_saw.py —— 《Elements of the Random Walk》第 7 章 Path integrals and self-avoidance 的數學
（Python＋numpy／scipy（integrate、special、optimize），不 import manim；不依賴其他模組，以便原樣複製進 combkit）

約定（書上式號與筆誤見 ref/md/ch07_formulas.md）：
  格子：一維（每步 ±1）、二維正方格子（每步 E、N、W、S 之一，格距 a = 1）；C(N; r) = 從原點走 N 步到 r 的走法數。
  一條 N 步的漫步走過 N+1 個點 w_0 = 0, w_1, …, w_N；「嚴格」自迴避 = N+1 個點兩兩不同（踩到就停）；c_N = 二維 N 步自迴避漫步的條數，
  ⟨R²⟩ = 終點到起點的平方距離對 c_N 條平均（普通漫步 ⟨R²⟩ = N）。
  1 − f 展開：e^{−v} = 1 − f（(7.17)）；嚴格自迴避時「兩點重合」f = 1、否則 0，所以 Π_{i>j}(1 − f_ij) 展開的第 k 階
  T_k = Σ_walks C(n_w, k)（n_w = 這條漫步的重合點對數），部分和 P_K = Σ_{k≤K} (−1)^k T_k，全部加完 = c_N（排容原理）。
  弧線圖（chord diagram）：一條線上 2k 個點配成 k 條弧（點編號 1..2k，沿漫步的順序）；不可約（1PI）＝找不到一個切口把弧分成左右兩群。
  平均場 (7.57)–(7.60)：格點數 𝒩（程式裡叫 Ns）、步數 m、x = z/z_c；ratio = 𝒩 z_c^m C_MF(m) = 自迴避讓條數變成普通漫步的幾倍
  （書上把格點數也寫成 N；本模組 N 只當步數）。Flory：F(R) = −R²/m − u m²/R^d（(7.63)），最大值 R* ∝ m^{3/(d+2)}（(7.64)）。
  分數用 fractions.Fraction；回傳快取（lru_cache）陣列或 list 的函數不要就地修改結果。

§7.1 路徑和 (7.1)–(7.6)，Ex. 7.1
    count1(N, x)                    一維 C(N; x) = C(N, (N+x)/2)（奇偶不合或 |x| > N 為 0）
    count2(N, x, y)                 二維 C(N; (x,y)) = C1(N, x+y)·C1(N, x−y)（兩枚硬幣）
    count_brute(N, end, dim)        暴力列舉 (2·dim)^N 條，數終點 = end 的條數（第二種算法，N ≤ 8）
    split_terms(N, M, end, dim)     (7.3)：[(中間點 r, C(M; r), C(N−M; end−r))]（乘積非零者；dim=1 時 r 是整數，dim=2 時是 (x, y)）
    split_sum(N, M, end, dim)       Σ 上式的乘積（＝ C(N; end)，與 M 無關）
    multi_split(N, Ms, end, dim)    (7.5)：切成 len(Ms) 段（Σ Ms = N），逐段摺積、對所有中間點求和
    EX_SPLIT                        卡上的例子 dict：一維 N=4 回原點 1·1+2·2+1·1 = 6；二維 N=4 回原點 Σ_r C(2;r)² = 16+4·1+4·4 = 36；
                                    (7.5) 切 2+2+2：一維 20、二維 400；M = N−1 退化成 (2.7)：2·3 = 6、4·9 = 36
    gauss1(M, x, d=3, a=1, book=False, book711=False)
                                    (7.6) 一段的高斯「一個直角分量」：正確 (d/(2πMa²))^{1/2} e^{−d x²/(2Ma²)}（變異數 Ma²/d）；
                                    book=True 照書上的前因子 (2πMa²)^{−1/2}（少了 d）；book711=True 照 (7.11) 指數的分母 3Ma²（變異數 Ma²/2）
    gauss_conv(M1, M2, x, **kw)     ∫ gauss1(M1, x−r) gauss1(M2, r) dr（數值積分）
    conv_ratio(M1, M2, x=0.7, book=False)   一個分量的 gauss_conv / gauss1(M1+M2)：正確 = 1；書上的前因子 = 3^{−1/2}（三維 3^{−3/2} = 0.1925）
    lattice_var(M, d=3)             簡單立方格子 M 步後一個分量的變異數（精確，由一維邊際分布加總）＝ M/d
    EX71                            Ex. 7.1 的例子：M1 = 3、M2 = 6（三維）→ 分量變異數 1、2、3；峰高 0.3989、0.2821、0.2303
§7.2–7.2.1 自迴避漫步（二維正方格子）
    saw_enum(Nmax=14)               DFS 列舉（用對稱：第一步向東、第一次轉彎向北）→ (c, sumR2)：N = 0..Nmax 的 tuple（int；約 0.2 秒，快取）
    saw_c(N)、saw_R2(N)             c_N 與 ⟨R²⟩（Fraction）；saw_frac(N) = c_N/4^N（Fraction）
    SAW_KNOWN_C                     文獻值 c_0..c_16（OEIS A001411，核對用）
    all_walks_xy(N)                 全部 4^N 條的點座標 (x, y)：各為 (4^N, N+1) int16（步伐 E,N,W,S 的字典序；N ≤ 10，快取）
    n_coinc(N)                      每條漫步的重合點對數 n_w（(4^N,) int；快取）
    saw_brute(N)                    由全部 4^N 條數 n_w = 0 的條數與 ΣR²（第二種算法，N ≤ 10）
    rw_R2_brute(N)                  普通漫步 ⟨R²⟩ 的暴力平均（Fraction，＝ N）
    saw_list(N)                     N 步自迴避漫步的點列 list（N ≤ 5；DFS、步伐順序 E,N,W,S）；N=2：12 條、N=3：36 條（第一步向東的 9 條 ×4）
    sample_saw(N, seed)、sample_rw(N, seed)   隨機抽一條（自迴避：對普通漫步拒絕取樣＝均勻抽樣）→ 點列；預設種子見 SAMPLE_SEEDS
    two_nu_eff(N)                   局部斜率 ln(⟨R²⟩_N/⟨R²⟩_{N−2}) / ln(N/(N−2))（→ 2ν = 3/2）；mu_eff(N) = c_N/c_{N−1}（→ μ ≈ 2.638，只引用）
§7.2.2 1 − f 展開 (7.17)–(7.19)，Fig. 7.1–7.4，Ex. 7.2
    ie_terms(N)                     [T_0, T_1, …]（由每條漫步的 n_w：T_k = Σ C(n_w, k)，int）
    ie_terms_subsets(N)             同上，但「對每個重合點對的子集 S 數同時滿足 S 的漫步」（第二種算法，N ≤ 6）
    order1_formula(N)               T_1 = Σ_{j<i} 4^{N−(i−j)}·C(i−j; 0)，C(n; 0) = C(n, n/2)²（第三種算法，任意 N）
    ie_partial(N)                   部分和 [P_0, P_1, …]（P_K = Σ_{k≤K} (−1)^k T_k；最後一個 = c_N）
    weighted_count(N, f)            溫和自迴避 Σ_w (1 − f)^{n_w}（f = 1 → c_N、f = 0 → 4^N；f 可為 Fraction）
    order2_types(N)                 T_2 依兩個點對的相對位置分類：{'sequential', 'nested', 'crossed', 'shared'}（shared = 共用一個點）
    matchings(k)                    2k 個點（1..2k）的全部配對（tuple of (a, b)，a < b）；chord_count(k) = (2k−1)!!
    is_irreducible(m)               不可約（1PI）；irreducible_count(k)：1, 2, 10, 74（k = 1..4）
    reversal_classes(k)             把頭尾倒過來視為同一種之後的種類數（k = 2：3、k = 3：11）
    SECOND_ORDER                    書上 Fig. 7.3：(a) 巢狀 ((1,4),(2,3))、(b) 並排 ((1,2),(3,4))、(c) 交叉 ((1,3),(2,4))
    dyson_counts(kmax)              (全部 (2k−1)!!, 不可約 I_k, 由 1/(1 − I(x)) 展開的係數)：1/(1−I) 的係數 = (2k−1)!!（計數版的 Dyson 方程）
§7.2.3–7.2.4 生成函數、Dyson 方程 (7.20)–(7.37)，Fig. 7.5–7.10，Ex. 7.3
    convolve(a, b)                  (7.20) 數列的摺積 Σ_{M} a[M] b[N−M]（list）
    EX_CONV                         H1 = H2 = 2^K（一維 K 步的走法數）→ 摺積 (N+1)2^N：1, 4, 12, 32, 80 ↔ h1h2 = 1/(1−2z)²
    gf_value(coeffs, z)             Σ c_N z^N；gf_check(z)：摺積的級數 vs h1(z)h2(z)（z = 0.1：1.5625）
    fourier_conv_error(n, seed)     (7.24) 的離散版：循環摺積的 DFT = 兩個 DFT 的乘積（回傳最大誤差）
    dyson_partial(g0, sigma, n)     g0 Σ_{k<n} (σ g0)^k 的部分和（Fraction list）；dyson_closed(g0, sigma) = 1/(1/g0 − σ)（(7.37)）
    dyson_book35(g0, sigma) = g0/(1 + σ g0)（書上 (7.35) 的號）；dyson_book33(g0, sigma) = g0 + g0²σ/(1 − σ)（照 (7.33) 的 g0Σ^n g0）
    DYSON_TOY = (2, 1/4)            部分和 2, 3, 7/2, 15/4, 31/8 → 4；書上 (7.35) 得 4/3、(7.33) 得 10/3
    g0_1d(z, q) = 1/(1 − 2z cos q)；sigma1_1d(z) = ∫dq/2π g0 = 1/√(1−4z²)（閉式）；sigma1_1d_quad(z)（數值積分）
    g_dyson_1d(z, q, f, n=None)     接觸交互作用 −f 的一階自能：1/(g0^{-1} + f·sigma1)；n 給整數時改用幾何級數前 n 項
§7.2.5 冪次數算 (7.38)–(7.42)
    power_exp(l, d) = 1 − l(4−d)/2（Fraction）；power_exp_count(l, d) = −(2l−1) + dl/2（數線條與積分，第二種算法）
    power_table(dmax=5, lmax=4)     {d: [l = 0..lmax 的冪次]}；d = 3：1, 1/2, 0, −1/2, −1；d = 4：全是 1；d = 5：1, 3/2, 2, 5/2, 3
    loop_slope(case)                數值積分量出的冪次：'l1d1' → −1/2、'l1d3'（截斷 Λ=1 的奇異部分）→ 1/2、'l2d1'（巢狀）→ −2
§7.2.6 平均場 (7.43)–(7.61)，Ex. 7.4–7.5
    G_MF(x, Ns, f)                  (7.57) 第二行 ∫_0^∞ exp[Ns(−fλ² − λ(1−x))] dλ（閉式 erfcx）；G_MF_quad（數值積分）
    cmf_ratio(m, Ns, f, method)     (7.59) 更正版 × Ns z_c^m = E[e^{−fT²/Ns}]，T ~ Gamma(m+1)；method 'quad'、'grid'（梯形）、
                                    'pcf'（拋物柱函數，m 小）、'laplace'（鞍點近似）
    cmf_asym(m, Ns, f) = e^{−f m²/Ns}（(7.60)）；cmf_book59(m, Ns, f)：照書上 (7.59) 指數裡的 Nλ（f = 0 時得 Ns^{−(m+1)} 而不是 1）
    MF_CASES                        Ex. 7.4 的三組 (m, Ns, f)：(10, 100, 1)、(100, 10⁴, 1)、(1000, 10⁶, 1)（fm²/Ns 都是 1）
    ex75(Ns, x=1/2, f=1)            Ns(1−x)G_MF → 1（z < z_c：條數 z_c^{−m}/Ns，自迴避沒有效果）；ex75_big(Ns, x=3/2, f=1) = ln G_MF / Ns → (x−1)²/4f
    gmf_series(x, Ns, f, mmax)      Σ_m cmf_ratio(m)/Ns · x^m（(7.58) 更正版，核對 G_MF）
    eq755(z0, z, zc)                ∫_0^∞ e^{−t(z0−z)/zc} dt（數值）＝ zc/(z0−z)（書上說是 1/(z0−z)）
    G1_lattice(z, x)                一維格子的生成函數 λ^{|x|}/√(1−4z²)；eq747_residual(z, x, book)：(7.47) 兩邊之差（更正版 = 0）
    eq748_eigen(q, book)            ∇²ψ/ψ 的數值（二階差分）：e^{iqx} 得 −q²；書上的 e^{−qx} 得 +q²
§7.2.7 Flory (7.62)–(7.64)
    flory_F(R, m, d, u=1)；flory_Rstar(m, d, u=1) = (d u m³/2)^{1/(d+2)}（解析）；flory_Rstar_num（數值最大化）
    flory_nu(d) = 3/(d+2)（Fraction）；flory_slope(d)：log R* 對 log m 的數值斜率（m = 10³ → 10⁶）
    NU_QUOTED                       已知（只引用）：d=1：1、d=2：3/4（精確）、d=3：0.588（數值）、d≥4：1/2
summary()                         場景卡與講義用的全部數字（dict，字串已照卡上的小數位格式化）

check() 驗證全部（每個結果至少兩種獨立算法、筆誤的更正、plan_ch07.md 場景卡與講義的每個數字）；直接執行本檔會跑 check()（約 1 秒）。
`python walk_saw.py --summary` 印出 summary() 的全部數字。
最慢的呼叫：saw_enum(14) 約 0.2 秒、all_walks_xy(10)＋n_coinc(10) 約 0.3 秒（都有快取），其餘都很快。
"""
import itertools
import math
from fractions import Fraction
from functools import lru_cache

import numpy as np
from scipy import integrate, optimize, special

F = Fraction

# ---------------------------------------------------------------- §7.1 路徑和

STEPS2 = ((1, 0), (0, 1), (-1, 0), (0, -1))          # E, N, W, S


def count1(N, x):
    """一維 N 步從 0 到 x 的走法數。"""
    if N < 0 or abs(x) > N or (N + x) % 2:
        return 0
    return math.comb(N, (N + x) // 2)


def count2(N, x, y):
    """二維正方格子 N 步從原點到 (x, y) 的走法數（u = x+y、v = x−y 兩個獨立的一維）。"""
    return count1(N, x + y) * count1(N, x - y)


def _count(N, r, dim):
    return count1(N, r) if dim == 1 else count2(N, r[0], r[1])


def count_brute(N, end, dim):
    """暴力列舉 (2·dim)^N 條路線，數終點 = end 的條數。"""
    if dim == 1:
        return sum(1 for s in itertools.product((1, -1), repeat=N) if sum(s) == end)
    n = 0
    for s in itertools.product(STEPS2, repeat=N):
        if sum(a for a, _ in s) == end[0] and sum(b for _, b in s) == end[1]:
            n += 1
    return n


def _sites(M, dim):
    if dim == 1:
        return [r for r in range(-M, M + 1)]
    return [(x, y) for x in range(-M, M + 1) for y in range(-M, M + 1) if abs(x) + abs(y) <= M]


def split_terms(N, M, end, dim):
    """(7.3) 的每一項：(中間點 r, C(M; 0→r), C(N−M; r→end))，只列乘積非零者。"""
    out = []
    for r in _sites(M, dim):
        a = _count(M, r, dim)
        rest = (end - r) if dim == 1 else (end[0] - r[0], end[1] - r[1])
        b = _count(N - M, rest, dim)
        if a and b:
            out.append((r, a, b))
    return out


def split_sum(N, M, end, dim):
    return sum(a * b for _, a, b in split_terms(N, M, end, dim))


def multi_split(N, Ms, end, dim):
    """(7.5)：依序走 Ms[0], Ms[1], … 步，對每個中間點求和（逐段摺積的分布）。"""
    assert sum(Ms) == N
    dist = {0 if dim == 1 else (0, 0): 1}
    for M in Ms:
        new = {}
        for r0, c0 in dist.items():
            for s in _sites(M, dim):
                c = _count(M, s, dim)
                if not c:
                    continue
                r = r0 + s if dim == 1 else (r0[0] + s[0], r0[1] + s[1])
                new[r] = new.get(r, 0) + c0 * c
        dist = new
    return dist.get(end, 0)


EX_SPLIT = {
    "1d_terms": [(-2, 1, 1), (0, 2, 2), (2, 1, 1)],          # 1·1 + 2·2 + 1·1
    "1d_total": 6,
    "2d_center": 16, "2d_axis": 4, "2d_diag": 16,              # 4² + 4·1² + 4·2²
    "2d_total": 36,
    "split222_1d": 20, "split222_2d": 400,
    "eq27_1d": (2, 3), "eq27_2d": (4, 9),                      # M = N−1：鄰居數 × C(3; 鄰居)
}


def gauss1(M, x, d=3, a=1.0, book=False, book711=False):
    """(7.6) 一段高斯的一個直角分量。"""
    denom = (3.0 if book711 else 2.0) * M * a * a
    pref = (2 * math.pi * M * a * a) ** -0.5 if book else (d / (2 * math.pi * M * a * a)) ** 0.5
    return pref * np.exp(-d * np.asarray(x, dtype=float) ** 2 / denom)


def gauss_conv(M1, M2, x, **kw):
    val, _ = integrate.quad(lambda r: float(gauss1(M1, x - r, **kw)) * float(gauss1(M2, r, **kw)), -np.inf, np.inf)
    return val


def conv_ratio(M1, M2, x=0.7, book=False):
    return gauss_conv(M1, M2, x, book=book) / float(gauss1(M1 + M2, x, book=book))


def lattice_var(M, d=3):
    """簡單立方格子 M 步後 x 分量的變異數：x 分量每步以 1/d 的機率 ±1。"""
    tot = F(0)
    for k in range(M + 1):                                   # k 步沿 x 軸
        pk = F(math.comb(M, k)) * F(1, d) ** k * F(d - 1, d) ** (M - k)
        tot += pk * k                                          # 一維 k 步的 ⟨x²⟩ = k
    return tot


EX71 = {"M1": 3, "M2": 6, "var": (1, 2, 3)}

# ---------------------------------------------------------------- §7.2 自迴避漫步的列舉

SAW_KNOWN_C = [1, 4, 12, 36, 100, 284, 780, 2172, 5916, 16268, 44100, 120292, 324932, 881500,
               2374444, 6416596, 17245332]


@lru_cache(maxsize=None)
def saw_enum(Nmax=14):
    """二維自迴避漫步：c_N 與 Σ R²（N = 0..Nmax）。對稱：第一步向東（×4），第一次轉彎向北（×2）。"""
    off = Nmax + 2
    W = 2 * off + 1
    vis = bytearray(W * W)
    cnt = [0] * (Nmax + 1)
    s2 = [0] * (Nmax + 1)
    moves = STEPS2

    def rec(x, y, n):
        cnt[n] += 1
        s2[n] += x * x + y * y
        if n == Nmax:
            return
        for dx, dy in moves:
            nx = x + dx
            ny = y + dy
            k = (nx + off) * W + ny + off
            if not vis[k]:
                vis[k] = 1
                rec(nx, ny, n + 1)
                vis[k] = 0

    vis[off * W + off] = 1
    for k in range(1, Nmax + 1):                             # 先直走 k 步向東
        vis[(k + off) * W + off] = 1
        if k < Nmax:
            idx = (k + off) * W + 1 + off                      # 轉向北
            vis[idx] = 1
            rec(k, 1, k + 1)
            vis[idx] = 0
    c = [1] + [4 * (1 + 2 * cnt[n]) for n in range(1, Nmax + 1)]
    sr = [0] + [4 * (n * n + 2 * s2[n]) for n in range(1, Nmax + 1)]
    return tuple(c), tuple(sr)


def saw_c(N):
    return saw_enum(max(14, N))[0][N]


def saw_R2(N):
    c, sr = saw_enum(max(14, N))
    return F(sr[N], c[N]) if N else F(0)


def saw_frac(N):
    return F(saw_c(N), 4 ** N)


@lru_cache(maxsize=None)
def all_walks_xy(N):
    assert N <= 10
    idx = np.arange(4 ** N, dtype=np.int64)
    dx = np.array([1, 0, -1, 0], dtype=np.int16)
    dy = np.array([0, 1, 0, -1], dtype=np.int16)
    x = np.zeros((4 ** N, N + 1), dtype=np.int16)
    y = np.zeros((4 ** N, N + 1), dtype=np.int16)
    for s in range(N):
        d = (idx >> (2 * (N - 1 - s))) & 3
        x[:, s + 1] = x[:, s] + dx[d]
        y[:, s + 1] = y[:, s] + dy[d]
    return x, y


@lru_cache(maxsize=None)
def _keys(N):
    x, y = all_walks_xy(N)
    return (x.astype(np.int32) + N) * (2 * N + 1) + (y.astype(np.int32) + N)


def _pairs(N):
    """可能重合的點對 (j, i)：i − j 是 ≥ 2 的偶數（正方格子是二部圖）。"""
    return [(j, i) for j in range(N + 1) for i in range(j + 2, N + 1, 2)]


@lru_cache(maxsize=None)
def n_coinc(N):
    key = _keys(N)
    n = np.zeros(4 ** N, dtype=np.int32)
    for j, i in _pairs(N):
        n += key[:, j] == key[:, i]
    return n


def saw_brute(N):
    n = n_coinc(N)
    x, y = all_walks_xy(N)
    ok = n == 0
    r2 = x[:, -1].astype(np.int64) ** 2 + y[:, -1].astype(np.int64) ** 2
    return int(ok.sum()), int(r2[ok].sum())


def rw_R2_brute(N):
    x, y = all_walks_xy(N)
    r2 = x[:, -1].astype(np.int64) ** 2 + y[:, -1].astype(np.int64) ** 2
    return F(int(r2.sum()), 4 ** N)


def saw_list(N):
    assert N <= 5
    out = []

    def rec(path):
        if len(path) == N + 1:
            out.append(tuple(path))
            return
        x, y = path[-1]
        for dx, dy in STEPS2:
            p = (x + dx, y + dy)
            if p not in path:
                path.append(p)
                rec(path)
                path.pop()

    rec([(0, 0)])
    return out


SAMPLE_SEEDS = {"saw": 702, "rw": 702}


def sample_rw(N, seed=702):
    rng = np.random.default_rng(seed)
    d = rng.integers(0, 4, N)
    pts = [(0, 0)]
    for k in d:
        pts.append((pts[-1][0] + STEPS2[k][0], pts[-1][1] + STEPS2[k][1]))
    return pts


def sample_saw(N, seed=702):
    rng = np.random.default_rng(seed)
    while True:
        d = rng.integers(0, 4, N)
        pts = [(0, 0)]
        for k in d:
            pts.append((pts[-1][0] + STEPS2[k][0], pts[-1][1] + STEPS2[k][1]))
        if len(set(pts)) == N + 1:
            return pts


def two_nu_eff(N):
    return math.log(saw_R2(N) / saw_R2(N - 2)) / math.log(N / (N - 2))


def mu_eff(N):
    return saw_c(N) / saw_c(N - 1)


# ---------------------------------------------------------------- §7.2.2 1 − f 展開

def ie_terms(N):
    cnt = np.bincount(n_coinc(N))
    kmax = len(cnt) - 1
    return [sum(int(cnt[j]) * math.comb(j, k) for j in range(k, kmax + 1)) for k in range(kmax + 1)]


def ie_terms_subsets(N):
    assert N <= 6
    key = _keys(N)
    P = _pairs(N)
    eq = [key[:, j] == key[:, i] for j, i in P]
    T = [0] * (len(P) + 1)

    def rec(start, mask, k):
        T[k] += int(mask.sum())
        for p in range(start, len(P)):
            m2 = mask & eq[p]
            if m2.any():
                rec(p + 1, m2, k + 1)

    rec(0, np.ones(4 ** N, dtype=bool), 0)
    while len(T) > 1 and T[-1] == 0:
        T.pop()
    return T


def order1_formula(N):
    return sum(4 ** (N - (i - j)) * math.comb(i - j, (i - j) // 2) ** 2 for j, i in _pairs(N))


def ie_partial(N, terms=None):
    T = terms if terms is not None else ie_terms(N)
    out, s = [], 0
    for k, t in enumerate(T):
        s += (-1) ** k * t
        out.append(s)
    return out


def weighted_count(N, f):
    cnt = np.bincount(n_coinc(N))
    return sum(int(c) * (1 - f) ** j for j, c in enumerate(cnt))


def _ptype(p, q):
    (a, b), (c, e) = sorted([p, q])
    if len({a, b, c, e}) < 4:
        return "shared"
    if b < c:
        return "sequential"
    if e < b:
        return "nested"
    return "crossed"


def order2_types(N):
    key = _keys(N)
    P = _pairs(N)
    eq = [key[:, j] == key[:, i] for j, i in P]
    out = {"sequential": 0, "nested": 0, "crossed": 0, "shared": 0}
    for s in range(len(P)):
        for t in range(s + 1, len(P)):
            out[_ptype(P[s], P[t])] += int((eq[s] & eq[t]).sum())
    return out


def matchings(k):
    def rec(pts):
        if not pts:
            yield ()
            return
        a = pts[0]
        for idx in range(1, len(pts)):
            b = pts[idx]
            rest = pts[1:idx] + pts[idx + 1:]
            for m in rec(rest):
                yield ((a, b),) + m
    return [tuple(sorted(m)) for m in rec(tuple(range(1, 2 * k + 1)))]


def chord_count(k):
    return math.prod(range(1, 2 * k, 2))


def is_irreducible(m):
    n = 2 * len(m)
    for c in range(1, n):
        if not any(a <= c < b for a, b in m):
            return False
    return True


def irreducible_count(k):
    return sum(1 for m in matchings(k) if is_irreducible(m))


def reversal_classes(k):
    n = 2 * k + 1
    seen = set()
    for m in matchings(k):
        r = tuple(sorted((n - b, n - a) for a, b in m))
        seen.add(min(m, r))
    return len(seen)


SECOND_ORDER = {"a": ("nested", ((1, 4), (2, 3))), "b": ("sequential", ((1, 2), (3, 4))),
                "c": ("crossed", ((1, 3), (2, 4)))}


def dyson_counts(kmax=4):
    full = [chord_count(k) for k in range(kmax + 1)]
    irr = [0] + [irreducible_count(k) for k in range(1, kmax + 1)]
    D = [1] + [0] * kmax                                      # D = 1 + I·D
    for n in range(1, kmax + 1):
        D[n] = sum(irr[j] * D[n - j] for j in range(1, n + 1))
    return full, irr, D


# ---------------------------------------------------------------- §7.2.3–7.2.4 生成函數與 Dyson

def convolve(a, b):
    n = min(len(a), len(b))
    return [sum(a[M] * b[N - M] for M in range(N + 1)) for N in range(n)]


EX_CONV = {"H": [1, 2, 4, 8, 16], "conv": [1, 4, 12, 32, 80]}


def gf_value(coeffs, z):
    return sum(c * z ** n for n, c in enumerate(coeffs))


def gf_check(z=0.1, n=80):
    H = [2 ** k for k in range(n)]
    return gf_value(convolve(H, H), z), (1 / (1 - 2 * z)) ** 2


def fourier_conv_error(n=16, seed=7):
    rng = np.random.default_rng(seed)
    a, b = rng.normal(size=n), rng.normal(size=n)
    circ = np.array([sum(a[m] * b[(k - m) % n] for m in range(n)) for k in range(n)])
    return float(np.max(np.abs(np.fft.fft(circ) - np.fft.fft(a) * np.fft.fft(b))))


def dyson_partial(g0, sigma, n):
    g0, sigma = F(g0), F(sigma)
    out, s = [], F(0)
    for k in range(n):
        s += g0 * (sigma * g0) ** k
        out.append(s)
    return out


def dyson_closed(g0, sigma):
    return 1 / (1 / F(g0) - F(sigma))


def dyson_book35(g0, sigma):
    return F(g0) / (1 + F(sigma) * F(g0))


def dyson_book33(g0, sigma):
    g0, sigma = F(g0), F(sigma)
    return g0 + g0 * g0 * sigma / (1 - sigma)


DYSON_TOY = (2, F(1, 4))


def g0_1d(z, q):
    return 1.0 / (1.0 - 2.0 * z * np.cos(q))


def sigma1_1d(z):
    return 1.0 / math.sqrt(1.0 - 4.0 * z * z)


def sigma1_1d_quad(z):
    val, _ = integrate.quad(lambda q: g0_1d(z, q), -math.pi, math.pi)
    return val / (2 * math.pi)


def g_dyson_1d(z, q, f, n=None):
    g0 = g0_1d(z, q)
    s = -f * sigma1_1d(z)
    if n is None:
        return 1.0 / (1.0 / g0 - s)
    return sum(g0 * (s * g0) ** k for k in range(n))


# ---------------------------------------------------------------- §7.2.5 冪次

def power_exp(l, d):
    return 1 - F(l * (4 - d), 2)


def power_exp_count(l, d):
    return -(2 * l - 1) + F(d * l, 2)


def power_table(dmax=5, lmax=4):
    return {d: [power_exp(l, d) for l in range(lmax + 1)] for d in range(1, dmax + 1)}


def _loop_val(case, eps):
    if case == "l1d1":
        return integrate.quad(lambda q: 1 / (eps + q * q), -np.inf, np.inf)[0]
    if case == "l1d3":
        f = lambda q: 4 * math.pi * q * q / (eps + q * q)
        f0 = 4 * math.pi                                       # ε = 0 的被積函數
        return -integrate.quad(lambda q: f0 - f(q), 0, 1, points=[math.sqrt(eps)])[0]
    if case == "l2d1":
        a = integrate.quad(lambda q: 1 / (eps + q * q) ** 2, -np.inf, np.inf)[0]
        b = integrate.quad(lambda q: 1 / (eps + q * q), -np.inf, np.inf)[0]
        return a * b
    raise ValueError(case)


def loop_slope(case, e1=1e-6, e2=1e-8):
    v1, v2 = abs(_loop_val(case, e1)), abs(_loop_val(case, e2))
    return math.log(v1 / v2) / math.log(e1 / e2)


# ---------------------------------------------------------------- §7.2.6 平均場

def G_MF(x, Ns, f):
    A, B = Ns * f, Ns * (1 - x)
    return 0.5 * math.sqrt(math.pi / A) * float(special.erfcx(B / (2 * math.sqrt(A))))


def G_MF_quad(x, Ns, f):
    lam_star = max(0.0, (x - 1) / (2 * f))
    g = lambda lam: math.exp(Ns * (-f * lam * lam - lam * (1 - x)) - Ns * (-f * lam_star ** 2 - lam_star * (1 - x)))
    val = integrate.quad(g, 0, np.inf, points=None, limit=200)[0] if lam_star == 0 else \
        integrate.quad(g, 0, lam_star)[0] + integrate.quad(g, lam_star, np.inf)[0]
    return val * math.exp(Ns * (-f * lam_star ** 2 - lam_star * (1 - x)))


def _phi(t, m, Ns, f):
    return -f * t * t / Ns - t + m * np.log(t) - special.gammaln(m + 1)


def cmf_ratio(m, Ns, f, method="quad"):
    """Ns z_c^m C_MF(m) = (1/m!) ∫_0^∞ t^m e^{−t − f t²/Ns} dt（t = Ns λ）。"""
    a = f / Ns
    if m == 0:
        return integrate.quad(lambda t: math.exp(-t - a * t * t), 0, np.inf)[0]
    ts = (-1 + math.sqrt(1 + 8 * a * m)) / (4 * a) if a > 0 else float(m)   # 鞍點：m/t − 1 − 2at = 0
    sig = math.sqrt(m + 1.0)
    lo, hi = max(1e-12, ts - 40 * sig), ts + 40 * sig
    p0 = float(_phi(ts, m, Ns, f))
    if method == "quad":
        val = integrate.quad(lambda t: math.exp(float(_phi(t, m, Ns, f)) - p0), lo, hi, points=[ts], limit=200)[0]
        return val * math.exp(p0)
    if method == "grid":
        t = np.linspace(lo, hi, 40001)
        y = np.exp(_phi(t, m, Ns, f) - p0)
        return float(np.sum((y[1:] + y[:-1]) * np.diff(t)) / 2) * math.exp(p0)
    if method == "pcf":
        v = -(m + 1)
        return (2 * a) ** (v / 2) * math.exp(1 / (8 * a)) * float(special.pbdv(v, 1 / math.sqrt(2 * a))[0])
    if method == "laplace":
        p2 = -2 * a - m / ts ** 2
        return math.exp(p0) * math.sqrt(2 * math.pi / -p2)
    raise ValueError(method)


def cmf_asym(m, Ns, f):
    return math.exp(-f * m * m / Ns)


def cmf_book59(m, Ns, f):
    """照書上 (7.59)：指數 −Ns(fλ² + Ns λ) + m ln λ，乘 Ns^{m+1}/m!（對數空間數值積分）。"""
    lo = lambda lam: (m + 1) * math.log(Ns) - math.lgamma(m + 1) - Ns * (f * lam * lam + Ns * lam) + m * math.log(lam)
    ls = m / Ns ** 2
    p0 = lo(ls)
    val = integrate.quad(lambda lam: math.exp(lo(lam) - p0), 1e-300, 50 * (m + 1) / Ns ** 2, points=[ls], limit=200)[0]
    return val * math.exp(p0)


MF_CASES = ((10, 100, 1), (100, 10 ** 4, 1), (1000, 10 ** 6, 1))


def ex75(Ns, x=0.5, f=1.0):
    return Ns * (1 - x) * G_MF(x, Ns, f)


def ex75_big(Ns, x=1.5, f=1.0):
    A, B = Ns * f, Ns * (1 - x)
    # ln G = ln(½√(π/A)) + ln erfcx(B/2√A)；B < 0 時 erfcx 很大，用 erfc 的對數
    t = B / (2 * math.sqrt(A))
    lg = t * t + math.log(float(special.erfc(t)))
    return (math.log(0.5 * math.sqrt(math.pi / A)) + lg) / Ns


def gmf_series(x, Ns, f, mmax=400):
    return sum(cmf_ratio(m, Ns, f) / Ns * x ** m for m in range(mmax))


def eq755(z0, z, zc):
    return integrate.quad(lambda t: math.exp(-t * (z0 - z) / zc), 0, np.inf)[0]


def G1_lattice(z, x):
    s = math.sqrt(1 - 4 * z * z)
    lam = (1 - s) / (2 * z)
    return lam ** abs(x) / s


def eq747_residual(z, x, book=False):
    """一維、a = 1、z_c = 1/2、δ 換成 δ_{x,0}：更正版 (1−z/z_c)G − zΔG − δ；書上 zΔG − (z/z_c − 1)G − δ。"""
    G = lambda u: G1_lattice(z, u)
    lap = G(x + 1) + G(x - 1) - 2 * G(x)
    delta = 1.0 if x == 0 else 0.0
    if book:
        return z * lap - (2 * z - 1) * G(x) - delta
    return (1 - 2 * z) * G(x) - z * lap - delta


def eq748_eigen(q, book=False, x=0.3, h=1e-3):
    psi = (lambda u: np.exp(-q * u)) if book else (lambda u: np.exp(1j * q * u))
    lap = (psi(x + h) + psi(x - h) - 2 * psi(x)) / h ** 2
    return complex(lap / psi(x)).real


# ---------------------------------------------------------------- §7.2.7 Flory

def flory_F(R, m, d, u=1.0, c=1.0):
    """(7.63) 的指數；c = 熵那一項的係數（書上 1；高斯的正確係數 d/2 不改變指數）。"""
    return -c * R * R / m - u * m * m / R ** d


def flory_Rstar(m, d, u=1.0, c=1.0):
    return (d * u * m ** 3 / (2 * c)) ** (1 / (d + 2))


def flory_Rstar_num(m, d, u=1.0, c=1.0):
    r0 = flory_Rstar(m, d, u, c) * 1.7                         # 從偏離解析解的地方開始找
    res = optimize.minimize_scalar(lambda s: -flory_F(math.exp(s), m, d, u, c), bounds=(math.log(r0) - 5, math.log(r0) + 5),
                                   method="bounded", options={"xatol": 1e-12})
    return math.exp(res.x)


def flory_nu(d):
    return F(3, d + 2)


def flory_slope(d, m1=1e3, m2=1e6, c=1.0):
    return math.log(flory_Rstar_num(m2, d, c=c) / flory_Rstar_num(m1, d, c=c)) / math.log(m2 / m1)


NU_QUOTED = {1: 1.0, 2: 0.75, 3: 0.588, 4: 0.5}


# ---------------------------------------------------------------- 卡片與講義的數字

def _kernel_var(M, **kw):
    m0 = integrate.quad(lambda x: float(gauss1(M, x, **kw)), -np.inf, np.inf)[0]
    m2 = integrate.quad(lambda x: x * x * float(gauss1(M, x, **kw)), -np.inf, np.inf)[0]
    return m0, m2 / m0


def summary():
    s = {}
    s["W701_1d"] = "1·1+2·2+1·1 = 6"
    s["W701_2d"] = "16 + 4·1 + 4·4 = 36"
    s["W701_split222"] = (multi_split(6, [2, 2, 2], 0, 1), multi_split(6, [2, 2, 2], (0, 0), 2))
    s["W701_ex71_peaks"] = tuple(f"{float(gauss1(M, 0)):.4f}" for M in (3, 6, 9))
    s["W701_book_ratio"] = (f"{conv_ratio(3, 6, book=True):.4f}", f"{conv_ratio(3, 6, book=True) ** 3:.4f}")
    s["W701_var711"] = (f"{_kernel_var(6)[1]:.4f}", f"{_kernel_var(6, book711=True)[1]:.4f}")
    Ns_ = (1, 2, 3, 4, 6, 8, 10, 12, 14)
    s["W702_c"] = tuple(saw_c(n) for n in Ns_)
    s["W702_frac"] = tuple(f"{float(saw_frac(n)):.4f}" for n in Ns_)
    s["W702_R2"] = tuple(f"{float(saw_R2(n)):.4f}" for n in Ns_)
    s["W702_samples_R2"] = (sum(v * v for v in sample_saw(10)[-1]), sum(v * v for v in sample_rw(10)[-1]))
    s["W703_T4"] = tuple(ie_terms(4))
    s["W703_P4"] = tuple(ie_partial(4))
    s["W703_P6"] = tuple(ie_partial(6))
    s["W703_P8"] = tuple(ie_partial(8))
    s["W703_types4"] = order2_types(4)
    s["W703_chords"] = tuple(chord_count(k) for k in range(1, 5))
    s["W704_conv"] = tuple(convolve(EX_CONV["H"], EX_CONV["H"]))
    s["W704_gf"] = f"{gf_check(0.1)[1]:.4f}"
    s["W704_dyson"] = tuple(str(v) for v in dyson_partial(*DYSON_TOY, 5)) + (str(dyson_closed(*DYSON_TOY)),)
    s["W704_book"] = (str(dyson_book35(*DYSON_TOY)), str(dyson_book33(*DYSON_TOY)))
    s["W704_irr"] = tuple(dyson_counts(4)[1][1:])
    s["W705_power"] = {d: tuple(str(v) for v in row) for d, row in power_table().items()}
    s["W705_loops"] = tuple(f"{loop_slope(c):.2f}" for c in ("l1d1", "l1d3", "l2d1"))
    s["W705_sigma1"] = tuple(f"{sigma1_1d(z):.4f}" for z in (0.4, 0.49, 0.499))
    s["W705_mf"] = tuple(f"{cmf_ratio(*c):.4f}" for c in MF_CASES) + (f"{math.exp(-1):.4f}",)
    s["W705_ex75"] = tuple(f"{ex75(n):.4f}" for n in (10, 100, 1000, 10000))
    s["W705_ex75_big"] = tuple(f"{ex75_big(n):.4f}" for n in (100, 1000, 10000)) + ("0.0625",)
    s["W706_Rstar"] = tuple(f"{flory_Rstar(1000, d):.2f}" for d in (1, 2, 3, 4))
    s["W706_slopes"] = tuple(f"{flory_slope(d):.4f}" for d in (1, 2, 3, 4))
    s["W706_2nu"] = tuple(f"{two_nu_eff(n):.4f}" for n in (6, 8, 10, 12, 14))
    s["W706_2nu_odd"] = tuple(f"{two_nu_eff(n):.4f}" for n in (7, 9, 11, 13))
    s["W706_mu"] = f"{mu_eff(14):.4f}"
    return s


# ---------------------------------------------------------------- 自我測試

def check():
    # ---- §7.1 路徑和：切段＝公式＝暴力列舉
    for N in range(1, 9):
        for end in range(-N, N + 1):
            direct = count1(N, end)
            assert all(split_sum(N, M, end, 1) == direct for M in range(N + 1))
        assert count_brute(N, 0, 1) == count1(N, 0)
    for N in range(1, 7):
        for end in ((0, 0), (1, 1), (2, 0), (3, 1)):
            direct = count2(N, *end)
            assert all(split_sum(N, M, end, 2) == direct for M in range(N + 1))
        assert count_brute(N, (0, 0), 2) == count2(N, 0, 0) == (math.comb(N, N // 2) ** 2 if N % 2 == 0 else 0)
    assert split_terms(4, 2, 0, 1) == EX_SPLIT["1d_terms"] and split_sum(4, 2, 0, 1) == 6 == count_brute(4, 0, 1)
    t2 = split_terms(4, 2, (0, 0), 2)
    cen = sum(a * b for r, a, b in t2 if r == (0, 0))
    axis = sum(a * b for r, a, b in t2 if r != (0, 0) and 0 in r)
    diag = sum(a * b for r, a, b in t2 if 0 not in r)
    assert (cen, axis, diag) == (16, 4, 16) and cen + axis + diag == 36 == count_brute(4, (0, 0), 2) == math.comb(4, 2) ** 2
    assert all(a * b == 4 for r, a, b in t2 if 0 not in r) and all(a * b == 1 for r, a, b in t2 if r != (0, 0) and 0 in r)
    assert multi_split(6, [2, 2, 2], 0, 1) == 20 == math.comb(6, 3) == count_brute(6, 0, 1)
    assert multi_split(6, [2, 2, 2], (0, 0), 2) == 400 == math.comb(6, 3) ** 2 == count_brute(6, (0, 0), 2)
    assert multi_split(8, [1, 3, 4], (2, 0), 2) == count2(8, 2, 0) and multi_split(7, [2, 2, 3], 1, 1) == count1(7, 1)
    assert split_terms(4, 3, 0, 1) == [(-1, 3, 1), (1, 3, 1)]                                  # M = N−1 → (2.7)
    assert [(a, b) for _, a, b in split_terms(4, 3, (0, 0), 2)] == [(9, 1)] * 4
    # ---- Ex. 7.1：高斯的摺積還是高斯（變異數相加）；書上的前因子少了 3
    for M1, M2 in ((3, 6), (1, 2), (5, 7)):
        assert abs(conv_ratio(M1, M2) - 1) < 1e-9 and abs(conv_ratio(M1, M2, x=2.1) - 1) < 1e-9
        assert abs(conv_ratio(M1, M2, book=True) - 3 ** -0.5) < 1e-9
    assert f"{conv_ratio(3, 6, book=True):.4f}" == "0.5774" and f"{conv_ratio(3, 6, book=True) ** 3:.4f}" == "0.1925"
    for M, v in zip((3, 6, 9), EX71["var"]):
        m0, var = _kernel_var(M)
        assert abs(m0 - 1) < 1e-9 and abs(var - v) < 1e-9 and lattice_var(M) == v == F(M, 3)
        assert abs(float(gauss1(M, 0)) - 1 / math.sqrt(2 * math.pi * v)) < 1e-12
    assert tuple(f"{float(gauss1(M, 0)):.4f}" for M in (3, 6, 9)) == ("0.3989", "0.2821", "0.2303")
    m0b, varb = _kernel_var(6, book=True)
    assert abs(m0b - 3 ** -0.5) < 1e-9                                                         # 書上的前因子：面積 3^{−1/2}
    assert f"{_kernel_var(6, book711=True)[1]:.4f}" == "3.0000" and lattice_var(6) == 2       # (7.11)：變異數 M/2 ≠ M/3
    # ---- §7.2 列舉：DFS＝文獻＝暴力
    c, sr = saw_enum(14)
    assert list(c) == SAW_KNOWN_C[:15]
    for N in range(0, 11):
        assert saw_brute(N) == (c[N], sr[N])
    for N in (2, 4, 8):
        assert rw_R2_brute(N) == N
    assert len(saw_list(2)) == 12 and len(saw_list(3)) == 36 and len(saw_list(4)) == 100
    assert len([w for w in saw_list(3) if w[1] == (1, 0)]) == 9 and all(len(set(w)) == 4 for w in saw_list(3))
    assert sum(sum(x * x + y * y for x, y in [w[-1]]) for w in saw_list(4)) == sr[4]
    sw, rw = sample_saw(10), sample_rw(10)
    assert len(set(sw)) == 11 and sw[-1] == (-3, 3) and rw[-1] == (1, -1) and len(rw) - len(set(rw)) == 3   # 踩回 3 次
    assert f"{float(saw_R2(14)):.2f}" == "42.79" and _pairs(4) == [(0, 2), (0, 4), (1, 3), (2, 4)]
    Ns_ = (1, 2, 3, 4, 6, 8, 10, 12, 14)
    assert tuple(saw_c(n) for n in Ns_) == (4, 12, 36, 100, 780, 5916, 44100, 324932, 2374444)
    assert tuple(f"{float(saw_frac(n)):.4f}" for n in Ns_) == ("1.0000", "0.7500", "0.5625", "0.3906", "0.1904", "0.0903",
                                                              "0.0421", "0.0194", "0.0088")
    assert tuple(f"{float(saw_R2(n)):.4f}" for n in Ns_) == ("1.0000", "2.6667", "4.5556", "7.0400", "12.5744", "19.0128",
                                                            "26.2425", "34.1870", "42.7864")
    assert saw_R2(2) == F(8, 3) and saw_R2(3) == F(41, 9) and saw_R2(4) == F(176, 25)
    assert tuple(f"{two_nu_eff(n):.4f}" for n in (6, 8, 10, 12, 14)) == ("1.4306", "1.4372", "1.4442", "1.4505", "1.4556")
    assert tuple(f"{two_nu_eff(n):.4f}" for n in (7, 9, 11, 13)) == ("1.4459", "1.4528", "1.4562", "1.4593")
    assert f"{mu_eff(14):.4f}" == "2.6936"
    # ---- §7.2.2 1 − f 展開：逐條漫步＝逐個子集＝一階公式；全部加完＝c_N
    assert ie_terms(4) == [256, 228, 96, 28, 4] == ie_terms_subsets(4)
    assert ie_partial(4) == [256, 28, 124, 96, 100] and ie_partial(4)[-1] == c[4]
    assert ie_terms(3) == ie_terms_subsets(3) and ie_terms(5) == ie_terms_subsets(5) and ie_terms(6) == ie_terms_subsets(6)
    assert 3 * 4 ** 2 * 4 + 36 == 228 == order1_formula(4)
    for N in range(2, 9):
        assert order1_formula(N) == ie_terms(N)[1]
        assert ie_partial(N)[-1] == c[N] == weighted_count(N, 1) and weighted_count(N, 0) == 4 ** N
        assert sum(order2_types(N).values()) == (ie_terms(N) + [0, 0])[2]
    assert ie_partial(6) == [4096, -3152, 4116, -1492, 1996, 268, 928, 748, 784, 780]
    assert ie_partial(8)[:4] == [65536, -119332, 190940, -222036] and ie_partial(8)[-1] == 5916 and len(ie_partial(8)) == 17
    assert order2_types(4) == {"sequential": 0, "nested": 16, "crossed": 32, "shared": 48}
    assert weighted_count(4, F(1, 2)) == F(651, 4)
    # 弧線圖
    for k in range(1, 6):
        ms = matchings(k)
        assert len(ms) == chord_count(k) == math.factorial(2 * k) // (2 ** k * math.factorial(k)) and len(set(ms)) == len(ms)
    assert [chord_count(k) for k in range(1, 5)] == [1, 3, 15, 105]
    assert [irreducible_count(k) for k in range(1, 5)] == [1, 2, 10, 74]
    assert reversal_classes(2) == 3 and reversal_classes(3) == 11
    assert sorted(m for _, m in SECOND_ORDER.values()) == sorted(matchings(2))
    assert [is_irreducible(SECOND_ORDER[k][1]) for k in "abc"] == [True, False, True]
    assert {_ptype(*SECOND_ORDER[k][1]) for k in "abc"} == {"nested", "sequential", "crossed"}
    assert all(_ptype(*SECOND_ORDER[k][1]) == SECOND_ORDER[k][0] for k in "abc")
    assert chord_count(3) - irreducible_count(3) == 5 == 1 * 1 * 1 + 1 * 2 + 2 * 1                  # Ex. 7.2：可約 5 種
    full, irr, D = dyson_counts(4)
    assert D == full == [1, 1, 3, 15, 105] and irr == [0, 1, 2, 10, 74]
    # ---- §7.2.3–7.2.4 生成函數、Dyson
    H = EX_CONV["H"]
    assert convolve(H, H) == EX_CONV["conv"] == [(N + 1) * 2 ** N for N in range(5)]
    a, b = gf_check(0.1)
    assert abs(a - b) < 1e-12 and f"{b:.4f}" == "1.5625"
    assert fourier_conv_error() < 1e-12
    g0, sg = DYSON_TOY
    assert dyson_partial(g0, sg, 5) == [2, 3, F(7, 2), F(15, 4), F(31, 8)] and dyson_closed(g0, sg) == 4
    assert abs(float(dyson_partial(g0, sg, 60)[-1]) - 4) < 1e-12 and 4 == 1 / (F(1, 2) - F(1, 4))
    assert dyson_book35(g0, sg) == F(4, 3) and dyson_book33(g0, sg) == F(10, 3)
    s33 = g0 + sum(g0 * sg ** n * g0 for n in range(1, 80))                                   # 照 (7.33) 的 g0 Σ^n g0
    assert abs(float(s33) - 10 / 3) < 1e-12
    for z in (0.1, 0.3, 0.4, 0.49):
        assert abs(sigma1_1d(z) - sigma1_1d_quad(z)) < 1e-9
        conv = abs(0.05 * sigma1_1d(z) * g0_1d(z, 0.5)) < 1                                    # |Σ g0| < 1 才收斂
        assert conv == (z < 0.45)
        if conv:
            assert abs(g_dyson_1d(z, 0.5, 0.05) - g_dyson_1d(z, 0.5, 0.05, n=400)) < 1e-9
    assert tuple(f"{sigma1_1d(z):.4f}" for z in (0.4, 0.49, 0.499)) == ("1.6667", "5.0252", "15.8193")
    sl = math.log(sigma1_1d(0.4999) / sigma1_1d(0.49999)) / math.log(0.0001 / 0.00001)
    assert abs(sl - power_exp(1, 1)) < 1e-3                                                     # 一維一階：−1/2
    # ---- §7.2.5 冪次
    for d in range(1, 6):
        for l in range(0, 5):
            assert power_exp(l, d) == power_exp_count(l, d)
    pt = power_table()
    assert pt[3] == [1, F(1, 2), 0, F(-1, 2), -1] and pt[4] == [1] * 5 and pt[5] == [1, F(3, 2), 2, F(5, 2), 3]
    assert pt[1] == [1, F(-1, 2), -2, F(-7, 2), -5] and pt[2] == [1, 0, -1, -2, -3]
    assert abs(loop_slope("l1d1") - power_exp(1, 1)) < 1e-6 and abs(loop_slope("l2d1") - power_exp(2, 1)) < 1e-6
    assert abs(loop_slope("l1d3") - power_exp(1, 3)) < 1e-3
    assert tuple(f"{loop_slope(cc):.2f}" for cc in ("l1d1", "l1d3", "l2d1")) == ("-0.50", "0.50", "-2.00")
    # ---- §7.2.6 平均場
    for m, Ns, f in MF_CASES:
        q, g, lp = (cmf_ratio(m, Ns, f, mt) for mt in ("quad", "grid", "laplace"))
        assert abs(q / g - 1) < 1e-9 and abs(lp / q - 1) < 2e-3
    assert abs(cmf_ratio(10, 100, 1, "pcf") / cmf_ratio(10, 100, 1) - 1) < 1e-9
    assert tuple(f"{cmf_ratio(*cc):.4f}" for cc in MF_CASES) == ("0.3385", "0.3643", "0.3675") and f"{math.exp(-1):.4f}" == "0.3679"
    assert all(abs(cmf_asym(*cc) - math.exp(-1)) < 1e-15 for cc in MF_CASES)
    devs = [abs(cmf_ratio(*cc) / cmf_asym(*cc) - 1) for cc in MF_CASES]
    assert devs[0] > devs[1] > devs[2] and devs[2] < 2e-3                                       # Ex. 7.4：越大越準
    for m in (0, 1, 5, 20):
        assert abs(cmf_ratio(m, 100, 0) - 1) < 1e-9                                            # f = 0：普通漫步 z_c^{−m}/𝒩
    assert abs(cmf_book59(10, 100, 0) / 100.0 ** -11 - 1) < 1e-6                                # 書上 (7.59)：f = 0 也不對
    assert tuple(f"{ex75(n):.4f}" for n in (10, 100, 1000, 10000)) == ("0.6899", "0.9341", "0.9922", "0.9992")
    for x in (0.5, 1.5, 0.9):
        for Ns in (20, 200):
            assert abs(G_MF(x, Ns, 1.0) / G_MF_quad(x, Ns, 1.0) - 1) < 1e-8
    assert abs(gmf_series(0.5, 20, 1.0) / G_MF(0.5, 20, 1.0) - 1) < 1e-9                         # (7.58) 更正版
    for Ns in (100, 1000):                                                                      # Ex. 7.5 的展開 1/B − 2A/B³
        A_, B_ = Ns * 1.0, Ns * 0.5
        assert abs(G_MF(0.5, Ns, 1.0) - (1 / B_ - 2 * A_ / B_ ** 3)) < 15 * A_ ** 2 / B_ ** 5
    assert tuple(f"{ex75_big(n):.4f}" for n in (100, 1000, 10000)) == ("0.0452", "0.0596", "0.0621")
    assert abs(math.log(G_MF(1.5, 1000, 1.0)) / 1000 - ex75_big(1000)) < 1e-9
    assert abs(eq755(0.25, 0.2, 0.25) - 5) < 1e-9 and abs(0.25 / (0.25 - 0.2) - 5) < 1e-12 and abs(1 / (0.25 - 0.2) - 20) < 1e-9
    for z in (0.2, 0.4):
        for x in (0, 1, 2, 3):
            assert abs(eq747_residual(z, x)) < 1e-12
    assert f"{eq747_residual(0.4, 1, book=True):.4f}" == "0.3333"
    assert abs(eq748_eigen(2.0) + 4) < 1e-4 and abs(eq748_eigen(2.0, book=True) - 4) < 1e-4
    # ---- §7.2.7 Flory
    for d in (1, 2, 3, 4):
        for m in (1e2, 1e4, 1e6):
            assert abs(flory_Rstar_num(m, d) / flory_Rstar(m, d) - 1) < 1e-6
        assert abs(flory_slope(d) - float(flory_nu(d))) < 1e-6
        assert abs(flory_slope(d, c=d / 2) - float(flory_nu(d))) < 1e-6                         # 正確的高斯係數 d/2：指數不變
        R = flory_Rstar(1000, d)
        assert flory_F(R, 1000, d) > max(flory_F(R * 0.99, 1000, d), flory_F(R * 1.01, 1000, d))
    assert [flory_nu(d) for d in (1, 2, 3, 4)] == [1, F(3, 4), F(3, 5), F(1, 2)]
    assert tuple(f"{flory_slope(d):.4f}" for d in (1, 2, 3, 4)) == ("1.0000", "0.7500", "0.6000", "0.5000")
    assert tuple(f"{flory_Rstar(1000, d):.2f}" for d in (1, 2, 3, 4)) == ("793.70", "177.83", "68.43", "35.50")
    R2s = flory_Rstar(1000, 2)
    assert f"{flory_F(R2s, 1000, 2):.2f}" == "-63.25" and abs(R2s ** 2 / 1000 - 1000 ** 2 / R2s ** 2) < 1e-9   # d=2：兩項相等
    assert f"{R2s ** 2 / 1000:.2f}" == "31.62"
    assert float(flory_nu(2)) == NU_QUOTED[2] and abs(float(flory_nu(3)) - NU_QUOTED[3]) < 0.02
    return True


if __name__ == "__main__":
    import sys
    import time
    import warnings
    sys.stdout.reconfigure(encoding="utf-8")
    warnings.simplefilter("ignore")
    if "--summary" in sys.argv:
        for k, v in summary().items():
            print(k, v)
        sys.exit(0)
    t0 = time.time()
    check()
    print(f"walk_saw（第 7 章）：路徑和的切段、高斯摺積、自迴避漫步的精確列舉、1−f 展開（排容）、弧線圖與 Dyson 方程、冪次數算、"
          f"平均場、Flory 指數全部驗證通過（{time.time() - t0:.1f} 秒）")
