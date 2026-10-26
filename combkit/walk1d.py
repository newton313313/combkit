# 由 combkit/tools/sync_from_walk017.py 從 walk017 專案的 codes/utils/walk1d.py 複製；請改 walk017 那份再同步。
"""
walk1d.py —— 《Elements of the Random Walk》第 1 章的數學（純 Python，不 import manim）

一維簡單隨機漫步：每步長 l = 1，擲硬幣正面（H）向右 +1、反面（T）向左 −1；N 步後的位置 d（向右為正）。
畫面上、講義裡出現的每個數字都從這裡算，場景與 make_figs.py 再 assert。書上式號見 ref/md/ch01_formulas.md。

計數與機率（精確，用 int／Fraction）
    walk_count(N, d)            C(N,d)，(1.6)；d 與 N 不同奇偶或 |d| > N 時為 0
    walk_prob(N, d)             P(N,d) = C(N,d)/2^N（Fraction），(1.7)
    distribution(N)             {d: C(N,d)}，d = −N, −N+2, …, N
    biased_prob(N, n, p)        n 次正面的機率 C(N,n) p^n (1−p)^(N−n)（Exercise 1.1；p 可用 Fraction）
    all_walks(N)                全部 2^N 條走法（±1 的 tuple），依字典序（+1 在前）
    sequences_with_heads(N, n)  恰好 n 個正面的序列（字串，如 "HHT"），依字典序（H 在前）
    coins_to_steps("HTH")       → [1, −1, 1]；positions(steps) → 部分和 [0, x1, x1+x2, …]
    MAIN_COINS = "HTHHTHTTHH"   全系列主例：十次擲幣（6 正 4 反），MAIN_STEPS、MAIN_POS（終點 d = 2）

模擬（固定亂數種子，可重現）
    ensemble(M, N, seed)        M 人各走 N 步 → list[list[int]]（每人的步伐）
    finals(walks)               每人的終點；mean(xs)、mean_square(xs)（Fraction）
    histogram(xs)               {d: 人數}

動差（平均都是對 2^N 條走法平均）
    moment_exact(N, k)          平均 d^k（Fraction），直接對分布加總，(1.15)
    moment_formula(N, k)        封閉式：k 奇數 0；2: N；4: 3N²−2N；6: 15N³−30N²+16N；8: 105N⁴−420N³+588N²−272N
    moment_via_cosh(N, k)       (1.24)：(−2)^k · d^k/dy^k cosh(y/2)^N |_{y=0}，用冪級數相乘（與加總無關的第二種算法）
    mean_heads(N)               平均正面數 N/2，(1.18)
    pair_grid_counts(N)         (1.12) 的 N×N 格子：(對角格數 N, 非對角格數 N(N−1))
    fourth_moment_pairings(N)   Exercise 1.3：Δ_iΔ_jΔ_kΔ_l 平均不為 0 的指標組數 N + 3N(N−1) = 3N²−2N

Stirling（§1.3）
    ln_factorial(n)；stirling(n) = sqrt(2πn) n^n e^−n；ln_stirling(n)；stirling_ratio(n) = n!/stirling(n)
    gamma_integrand(w, N) = w^N e^−w；gamma_gauss(w, N) = exp(N ln N − N − (w−N)²/(2N))（Fig. 1.3 的虛線）
    BOOK_LN10 = (15.1044, 15.0961)  (1.35)(1.36)

高斯近似
    gauss_density(u, N)         (1.41)：exp(−u²/2N)/sqrt(2πN)（對 u = d/l 的機率密度）
    gauss_prob(N, d)            離散 P(N,d) 的近似 = 2·gauss_density（奇偶性的因子 2，見 (1.40) 的註）
    ln_binom_approx(N, m)       (1.39) 正確版：N ln2 − 2m²/N − ½ ln(πN/2)（≈ ln C(N, N/2+m)）

遞迴與擴散（§1.4）
    pascal_rows(Nmax)           由 (1.42) 一列列算出 [{d: C(N,d)}，N = 0..Nmax]
    diffusion_step(P)           P_{N+1}(y) = ½[P_N(y−1) + P_N(y+1)]（dict，值可為 Fraction）
    second_difference(f, y)     f(y−1) + f(y+1) − 2f(y)
    heat_kernel(y, N, l=1)      (1.48) 正確版：exp(−y²/(2Nl²))/sqrt(2πNl²)

傅立葉與生成函數（§1.5）
    chi(q)                      2 cos q（l = 1）
    chi_power_coeffs(N)         (e^{iq}+e^{−iq})^N 的展開係數 {d: C(N,d)}
    fourier_count(N, d)         (1.51)：(1/2π)∫_{−π}^{π} e^{−iqd}(2cos q)^N dq，用 M > N+|d| 點均分求積（精確）
    saddle_q0(N, d)             (1.54)：2^N/sqrt(2πN) · e^{−d²/2N}（只取 q = 0 的鞍點）
    saddle_both(N, d)           再加上 q = ±π 的鞍點：saddle_q0 · (1 + (−1)^{N+d})
    series_coeffs(f, n)         冪級數係數（f 為已知係數的生成函數名稱，見下）
    gz0_coeffs(n)               g(z,0) = 1/(1−2z) 的係數 [2^N]，(1.65)
    d2_gen_coeffs(n)            −∂²g/∂q²|_{q=0} = 2z/(1−2z)² 的係數 [N·2^N]（Σ_d d² C(N,d)，Exercise 1.6）
    d4_gen_coeffs(n)            ∂⁴g/∂q⁴|_{q=0} 的係數 [2^N(3N²−2N)]（Exercise 1.7）
    cauchy_coeff(f, N, r, M)    (1.58)：在半徑 r 的圓周上取 M 點，數值算 z^N 的係數（f 是複變函數）
    walk_saddle_z(N, q)         (1.60) 的精確鞍點 (N+1)/((N+2)(2−q²))；書上近似 1/(2−q²)

最陡下降（§1.6，f(z) = (zc − z)^−α）
    sd_exact(N, a, zc=1)        精確係數 Γ(N+a)/(Γ(a) N!) · zc^(−N−a)
    sd_approx(N, a, zc=1)       (1:S-15) 正確版：(1/sqrt(2π)) a^(−a+1/2) e^a (N+1)^(a−1) zc^(−N−a)
    sd_book(N, a, zc=1)         書上印的版本（前因子 1/(2π)，錯）
    sd_saddle_exact(N, a, zc=1) z* = (N+1)zc/(N+1+a)；sd_saddle_book(N, a, zc) = zc − a·zc/(N+1)
    sd_p(z, N, a, zc=1)         p(z) = −a ln(zc−z) − (N+1) ln z（z 可為複數）；sd_p2(z, …) = p''(z)
    sd_full(N, a, zc=1)         用精確鞍點的鞍點積分 e^{p(z*)}/sqrt(2π p''(z*))

check() 驗證全部（含書上的數字）；直接執行本檔會跑 check()。
"""
from __future__ import annotations

import cmath
import math
import random
from fractions import Fraction
from itertools import product

# --------------------------------------------------------------------------- 主例
MAIN_COINS = "HTHHTHTTHH"


def coins_to_steps(coins: str) -> list[int]:
    return [1 if c == "H" else -1 for c in coins]


def positions(steps) -> list[int]:
    out, x = [0], 0
    for s in steps:
        x += s
        out.append(x)
    return out


MAIN_STEPS = coins_to_steps(MAIN_COINS)
MAIN_POS = positions(MAIN_STEPS)


# --------------------------------------------------------------------------- 計數與機率
def walk_count(N: int, d: int) -> int:
    if abs(d) > N or (N + d) % 2:
        return 0
    return math.comb(N, (N + d) // 2)


def walk_prob(N: int, d: int) -> Fraction:
    return Fraction(walk_count(N, d), 2 ** N)


def distribution(N: int) -> dict[int, int]:
    return {d: walk_count(N, d) for d in range(-N, N + 1, 2)}


def biased_prob(N: int, n: int, p) -> Fraction:
    p = Fraction(p)
    return math.comb(N, n) * p ** n * (1 - p) ** (N - n)


def all_walks(N: int) -> list[tuple[int, ...]]:
    return list(product((1, -1), repeat=N))


def sequences_with_heads(N: int, n: int) -> list[str]:
    return ["".join(s) for s in product("HT", repeat=N) if s.count("H") == n]


# --------------------------------------------------------------------------- 模擬
def ensemble(M: int, N: int, seed: int) -> list[list[int]]:
    rng = random.Random(seed)
    return [[1 if rng.random() < 0.5 else -1 for _ in range(N)] for _ in range(M)]


def finals(walks) -> list[int]:
    return [sum(w) for w in walks]


def mean(xs) -> Fraction:
    return Fraction(sum(xs), len(xs))


def mean_square(xs) -> Fraction:
    return Fraction(sum(x * x for x in xs), len(xs))


def histogram(xs) -> dict[int, int]:
    h: dict[int, int] = {}
    for x in xs:
        h[x] = h.get(x, 0) + 1
    return dict(sorted(h.items()))


# --------------------------------------------------------------------------- 動差
def moment_exact(N: int, k: int) -> Fraction:
    return Fraction(sum(c * d ** k for d, c in distribution(N).items()), 2 ** N)


def moment_formula(N: int, k: int) -> int:
    if k % 2:
        return 0
    table = {0: [1], 2: [0, 1], 4: [0, -2, 3], 6: [0, 16, -30, 15], 8: [0, -272, 588, -420, 105]}
    return sum(c * N ** j for j, c in enumerate(table[k]))


def _series_mul(a, b, n):
    out = [Fraction(0)] * (n + 1)
    for i, x in enumerate(a[:n + 1]):
        if x:
            for j, y in enumerate(b[:n + 1 - i]):
                out[i + j] += x * y
    return out


def moment_via_cosh(N: int, k: int) -> Fraction:
    """(1.24)：(−2)^k · k! · [y^k] cosh(y/2)^N，冪級數一次一次相乘（不經過分布加總）。"""
    ch = [Fraction(0)] * (k + 1)
    for j in range(0, k + 1, 2):
        ch[j] = Fraction(1, 2 ** j * math.factorial(j))       # cosh(y/2) = Σ (y/2)^{2j}/(2j)!
    acc = [Fraction(1)] + [Fraction(0)] * k
    for _ in range(N):
        acc = _series_mul(acc, ch, k)
    return (-2) ** k * math.factorial(k) * acc[k]


def mean_heads(N: int) -> Fraction:
    return Fraction(sum(math.comb(N, m) * m for m in range(N + 1)), 2 ** N)


def pair_grid_counts(N: int) -> tuple[int, int]:
    return N, N * (N - 1)


def fourth_moment_pairings(N: int) -> int:
    """Δ_iΔ_jΔ_kΔ_l 的平均不為 0：四個指標全相同（N 組）或兩兩成對（3·N(N−1) 組）。"""
    return sum(1 for t in product(range(N), repeat=4)
               if all(t.count(i) % 2 == 0 for i in set(t)))


# --------------------------------------------------------------------------- Stirling
BOOK_LN10 = (15.1044, 15.0961)


def ln_factorial(n: int) -> float:
    return math.lgamma(n + 1)


def stirling(n: float) -> float:
    return math.sqrt(2 * math.pi * n) * n ** n * math.exp(-n)


def ln_stirling(n: float) -> float:
    return 0.5 * math.log(2 * math.pi * n) + n * math.log(n) - n


def stirling_ratio(n: int) -> float:
    return math.exp(ln_factorial(n) - ln_stirling(n))


def gamma_integrand(w: float, N: int) -> float:
    return math.exp(N * math.log(w) - w) if w > 0 else 0.0


def gamma_gauss(w: float, N: int) -> float:
    return math.exp(N * math.log(N) - N - (w - N) ** 2 / (2 * N))


# --------------------------------------------------------------------------- 高斯近似
def gauss_density(u: float, N: int) -> float:
    return math.exp(-u * u / (2 * N)) / math.sqrt(2 * math.pi * N)


def gauss_prob(N: int, d: int) -> float:
    return 2 * gauss_density(d, N)


def ln_binom_approx(N: int, m: float) -> float:
    return N * math.log(2) - 2 * m * m / N - 0.5 * math.log(math.pi * N / 2)


# --------------------------------------------------------------------------- 遞迴與擴散
def pascal_rows(Nmax: int) -> list[dict[int, int]]:
    rows = [{0: 1}]
    for N in range(1, Nmax + 1):
        prev = rows[-1]
        rows.append({y: prev.get(y - 1, 0) + prev.get(y + 1, 0) for y in range(-N, N + 1, 2)})
    return rows


def diffusion_step(P: dict) -> dict:
    """一步：每個人各以 ½ 往左右走。P 的值用 Fraction（或 int）就是精確的。"""
    ys = sorted({y + s for y in P for s in (-1, 1)})
    return {y: Fraction(P.get(y - 1, 0) + P.get(y + 1, 0)) / 2 for y in ys}


def second_difference(f, y):
    return f(y - 1) + f(y + 1) - 2 * f(y)


def heat_kernel(y: float, N: float, l: float = 1.0) -> float:
    return math.exp(-y * y / (2 * N * l * l)) / math.sqrt(2 * math.pi * N * l * l)


# --------------------------------------------------------------------------- 傅立葉與生成函數
def chi(q: float) -> float:
    return 2 * math.cos(q)


def chi_power_coeffs(N: int) -> dict[int, int]:
    """把 (e^{iq}+e^{−iq})^N 一次一次乘開，記錄 e^{iqd} 的係數。"""
    poly = {0: 1}
    for _ in range(N):
        new: dict[int, int] = {}
        for d, c in poly.items():
            for s in (1, -1):
                new[d + s] = new.get(d + s, 0) + c
        poly = new
    return dict(sorted(poly.items()))


def fourier_count(N: int, d: int, M: int | None = None) -> float:
    M = M or 2 * (N + abs(d)) + 2
    tot = 0j
    for k in range(M):
        q = -math.pi + 2 * math.pi * k / M
        tot += cmath.exp(-1j * q * d) * (2 * math.cos(q)) ** N
    return (tot / M).real


def saddle_q0(N: int, d: float) -> float:
    return 2 ** N / math.sqrt(2 * math.pi * N) * math.exp(-d * d / (2 * N))


def saddle_both(N: int, d: int) -> float:
    return saddle_q0(N, d) * (1 + (-1) ** (N + d))


def gz0_coeffs(n: int) -> list[int]:
    return [2 ** N for N in range(n + 1)]


def d2_gen_coeffs(n: int) -> list[int]:
    return [N * 2 ** N for N in range(n + 1)]


def d4_gen_coeffs(n: int) -> list[int]:
    return [2 ** N * (3 * N * N - 2 * N) for N in range(n + 1)]


def cauchy_coeff(f, N: int, r: float, M: int = 256) -> complex:
    tot = 0j
    for k in range(M):
        z = r * cmath.exp(2j * math.pi * k / M)
        tot += f(z) / z ** N
    return tot / M


def walk_saddle_z(N: int, q: float) -> tuple[float, float]:
    return (N + 1) / ((N + 2) * (2 - q * q)), 1 / (2 - q * q)


# --------------------------------------------------------------------------- 最陡下降
def sd_exact(N: int, a: float, zc: float = 1.0) -> float:
    return math.exp(math.lgamma(N + a) - math.lgamma(a) - math.lgamma(N + 1)) * zc ** (-N - a)


def sd_approx(N: int, a: float, zc: float = 1.0) -> float:
    return a ** (-a + 0.5) * math.exp(a) * (N + 1) ** (a - 1) * zc ** (-N - a) / math.sqrt(2 * math.pi)


def sd_book(N: int, a: float, zc: float = 1.0) -> float:
    return a ** (-a + 0.5) * math.exp(a) * (N + 1) ** (a - 1) * zc ** (-N - a) / (2 * math.pi)


def sd_saddle_exact(N: int, a: float, zc: float = 1.0) -> float:
    return (N + 1) * zc / (N + 1 + a)


def sd_saddle_book(N: int, a: float, zc: float = 1.0) -> float:
    return zc - a * zc / (N + 1)


def sd_p(z, N: int, a: float, zc: float = 1.0):
    if isinstance(z, complex):
        return -a * cmath.log(zc - z) - (N + 1) * cmath.log(z)
    return -a * math.log(zc - z) - (N + 1) * math.log(z)


def sd_p2(z, N: int, a: float, zc: float = 1.0):
    return a / (zc - z) ** 2 + (N + 1) / z ** 2


def sd_full(N: int, a: float, zc: float = 1.0) -> float:
    z = sd_saddle_exact(N, a, zc)
    return math.exp(sd_p(z, N, a, zc)) / math.sqrt(2 * math.pi * sd_p2(z, N, a, zc))


# --------------------------------------------------------------------------- 自我測試
def check():
    # ---- 主例與 §1.1
    assert MAIN_COINS.count("H") == 6 and MAIN_POS[-1] == 2 and len(MAIN_POS) == 11
    assert MAIN_POS == [0, 1, 0, 1, 2, 1, 2, 1, 0, 1, 2]
    assert sequences_with_heads(3, 2) == ["HHT", "HTH", "THH"]            # Fig. 1.2
    assert math.factorial(3) // math.factorial(2) == 3
    assert distribution(4) == {-4: 1, -2: 4, 0: 6, 2: 4, 4: 1}
    assert list(distribution(10).values()) == [1, 10, 45, 120, 210, 252, 210, 120, 45, 10, 1]
    assert walk_count(10, 0) == 252 and walk_prob(10, 0) == Fraction(63, 256)
    assert walk_count(10, 2) == 210 and walk_count(10, 3) == 0 and walk_count(10, 12) == 0
    for N in range(0, 9):                                                  # 公式 = 暴力列舉
        brute = histogram(finals(all_walks(N)))
        assert brute == {d: c for d, c in distribution(N).items() if c}, N
        assert sum(distribution(N).values()) == 2 ** N
    assert biased_prob(3, 2, Fraction(2, 3)) == Fraction(4, 9)
    assert sum(biased_prob(5, n, Fraction(1, 3)) for n in range(6)) == 1
    assert biased_prob(10, 5, Fraction(1, 2)) == walk_prob(10, 0)

    # ---- §1.2 動差
    for N in range(0, 13):
        for k in range(0, 9):
            e = moment_exact(N, k)
            assert e == moment_formula(N, k), (N, k)
            assert moment_via_cosh(N, k) == e, (N, k)
        assert mean_heads(N) == Fraction(N, 2)
    assert moment_exact(10, 2) == 10 and moment_exact(10, 4) == 280
    assert moment_exact(100, 4) == 29800 and moment_exact(100, 2) == 100
    assert Fraction(moment_formula(10, 4), 10 ** 2) == Fraction(14, 5)    # 2.8
    for N in (1, 2, 3, 5):
        assert fourth_moment_pairings(N) == N + 3 * N * (N - 1) == moment_formula(N, 4)
    assert pair_grid_counts(4) == (4, 12)
    for N in (3, 6):                                                       # (1.12)：交叉項平均 0
        ws = all_walks(N)
        assert all(sum(w[j] * w[k] for w in ws) == 0 for j in range(N) for k in range(N) if j != k)

    # ---- 模擬（固定種子；場景要用的數字）
    f20 = finals(ensemble(20, 10, seed=1))                                 # W103：20 人、10 步
    assert f20 == [2, 0, 6, 4, -4, -2, -2, 0, -4, -6, 0, -2, 2, 6, 2, 0, 0, 0, -6, 6]
    assert mean(f20) == Fraction(1, 10) and mean_square(f20) == Fraction(63, 5)
    assert histogram(f20) == {-6: 2, -4: 2, -2: 3, 0: 6, 2: 3, 4: 1, 6: 3}
    big = finals(ensemble(2000, 100, seed=7))                              # W104：2000 人、100 步
    assert round(float(mean_square(big)), 3) == 98.922 and round(float(mean(big)), 3) == 0.121

    # ---- §1.3 Stirling
    assert round(ln_factorial(10), 4) == BOOK_LN10[0]
    assert round(ln_stirling(10), 4) == BOOK_LN10[1]
    rel10 = (ln_factorial(10) - ln_stirling(10)) / ln_factorial(10)
    rel100 = (ln_factorial(100) - ln_stirling(100)) / ln_factorial(100)
    assert 5e-4 < rel10 < 6e-4 and 2e-6 < rel100 < 3e-6                   # five parts in 10^4；two in 10^6
    assert math.factorial(10) == 3628800 and round(stirling(10)) == 3598696
    assert [round(stirling_ratio(n), 4) for n in (1, 2, 5, 10, 100)] == [1.0844, 1.0422, 1.0168, 1.0084, 1.0008]
    assert round(gamma_integrand(10, 10)) == 453999                        # Fig. 1.3 的峰高
    assert abs(gamma_gauss(10, 10) - gamma_integrand(10, 10)) < 1e-6
    # 被積函數的面積 = 10!；高斯的面積 = Stirling
    h = 0.01
    area = sum(gamma_integrand(k * h, 10) for k in range(1, 6000)) * h
    garea = sum(gamma_gauss(k * h, 10) for k in range(-3000, 6000)) * h
    assert abs(area / 3628800 - 1) < 1e-6 and abs(garea / stirling(10) - 1) < 1e-6
    # Exercise 1.4：ln N! − ln (N−1)! = ln N
    assert abs(ln_factorial(10) - ln_factorial(9) - math.log(10)) < 1e-12

    # ---- 高斯近似（含奇偶因子 2）
    assert round(gauss_density(0, 10), 5) == 0.12616
    assert round(gauss_prob(10, 0), 4) == 0.2523 and round(float(walk_prob(10, 0)), 4) == 0.2461
    assert round(gauss_prob(10, 4), 4) == 0.1134 and round(float(walk_prob(10, 4)), 4) == 0.1172
    for N in (100, 400):                                                   # (1.39) 正確版，誤差 ~ 1/(4N)
        assert abs(math.log(walk_count(N, 0)) - ln_binom_approx(N, 0)) < 3e-3
        assert abs(math.log(walk_count(N, 20)) - ln_binom_approx(N, 10)) < 5e-3
    assert abs(sum(gauss_prob(1000, d) for d in range(-1000, 1001, 2)) - 1) < 1e-9

    # ---- §1.4 遞迴與擴散
    rows = pascal_rows(10)
    assert all(rows[N] == distribution(N) for N in range(11))
    assert rows[6][0] == 20 == rows[5][-1] + rows[5][1] and rows[5][-1] == 10
    assert math.comb(6, 3) == math.comb(5, 2) + math.comb(5, 3)          # (1.43)
    P = {0: Fraction(1)}
    for N in range(1, 11):
        P = diffusion_step(P)
        assert P == {d: walk_prob(N, d) for d in range(-N, N + 1, 2)}
    assert sum(P.values()) == 1
    # 二階差分 ≈ 二階導數（對熱核）
    f = lambda y: heat_kernel(y, 50)                                      # noqa: E731
    d2 = (heat_kernel(1e-3, 50) - 2 * heat_kernel(0, 50) + heat_kernel(-1e-3, 50)) / 1e-6
    assert abs(second_difference(f, 0) - d2) < 2e-4
    # 擴散方程：∂P/∂N = ½ ∂²P/∂y²（熱核，數值微分）
    for y in (0.0, 3.0, 7.0):
        dPdN = (heat_kernel(y, 50 + 1e-4) - heat_kernel(y, 50 - 1e-4)) / 2e-4
        d2P = (heat_kernel(y + 1e-3, 50) - 2 * heat_kernel(y, 50) + heat_kernel(y - 1e-3, 50)) / 1e-6
        assert abs(dPdN - 0.5 * d2P) < 1e-6

    # ---- §1.5 傅立葉與生成函數
    assert chi_power_coeffs(3) == {-3: 1, -1: 3, 1: 3, 3: 1}
    for N in range(0, 13):
        assert chi_power_coeffs(N) == {d: c for d, c in distribution(N).items()}
        for d in range(-N - 2, N + 3):
            assert abs(fourier_count(N, d) - walk_count(N, d)) < 1e-6, (N, d)
    assert round(saddle_q0(10, 2), 1) == 105.8 and round(saddle_both(10, 2), 1) == 211.5
    assert saddle_both(10, 3) == 0
    assert round(saddle_both(10, 0), 1) == 258.4 and walk_count(10, 0) == 252
    assert gz0_coeffs(10)[10] == 1024
    for N, c in enumerate(d2_gen_coeffs(12)):
        assert c == sum(d * d * walk_count(N, d) for d in range(-N, N + 1))
    for N, c in enumerate(d4_gen_coeffs(12)):
        assert c == sum(d ** 4 * walk_count(N, d) for d in range(-N, N + 1))
    assert d2_gen_coeffs(10)[10] == 10240 and d4_gen_coeffs(10)[10] == 286720
    c10 = cauchy_coeff(lambda z: 1 / (1 - 2 * z), 10, 0.25)
    assert abs(c10 - 1024) < 1e-6
    c10q = cauchy_coeff(lambda z: 1 / (1 - z * chi(0.3)), 10, 0.25)
    assert abs(c10q - chi(0.3) ** 10) < 1e-6
    zs, zb = walk_saddle_z(10, 0.0)
    assert abs(zs - 11 / 24) < 1e-12 and zb == 0.5
    # (1.61)–(1.62)：鞍點值的指數部分 (2−q²)^{N+1} 與 2^N e^{−Nq²/2}
    q = 0.1                                                                # 誤差是下一項 (N+1)q⁴/8
    assert abs(math.log((2 - q * q) ** 11) - (11 * math.log(2) - 11 * q * q / 2)) < 1.1 * 11 * q ** 4 / 8

    # ---- §1.6 最陡下降
    assert abs(sd_exact(10, 2) - 11) < 1e-9 and abs(sd_exact(10, 1) - 1) < 1e-12
    assert abs(sd_exact(10, 0.5) - math.comb(20, 10) / 4 ** 10) < 1e-12
    assert round(sd_approx(10, 2) / 11, 4) == 1.0422 and round(sd_book(10, 2) / 11, 4) == 0.4158
    assert round(sd_approx(10, 2), 2) == 11.46 and round(sd_book(10, 2), 2) == 4.57
    assert round(sd_saddle_exact(10, 2), 4) == 0.8462 and round(sd_saddle_book(10, 2), 4) == 0.8182
    assert round(sd_full(10, 2), 2) == 10.59
    # 鞍點：沿實軸是最低點、沿垂直方向是最高點
    zs = sd_saddle_exact(10, 2)
    re_line = [sd_p(zs + t, 10, 2) for t in (-0.05, 0.0, 0.05)]
    im_line = [sd_p(complex(zs, t), 10, 2).real for t in (-0.05, 0.0, 0.05)]
    assert re_line[1] < re_line[0] and re_line[1] < re_line[2]
    assert im_line[1] > im_line[0] and im_line[1] > im_line[2]
    # (1:S-15) 等於「把 1/Γ(α) 換成 Stirling 近似」
    for a in (0.5, 2.0, 3.0):
        st = a ** (-a + 0.5) * math.exp(a) / math.sqrt(2 * math.pi)
        assert abs(sd_approx(50, a) - st * 51 ** (a - 1)) < 1e-9
        assert abs(1 / math.gamma(a) / st - 1) < 0.17
    assert sd_approx(2000, 3) / sd_exact(2000, 3) < 1.03                 # 大 N：只剩 Γ(α) 的 Stirling 誤差

    # ---- plan_ch01.md 場景卡上的其他數字
    assert [biased_prob(3, n, Fraction(2, 3)) for n in range(4)] == [Fraction(k, 27) for k in (1, 6, 12, 8)]
    assert sum(sum(w) ** 2 for w in all_walks(4)) == 64                   # W104
    assert [moment_formula(N, 4) for N in range(1, 7)] == [1, 8, 21, 40, 65, 96]
    assert [Fraction(moment_formula(N, 4), N * N) for N in range(1, 11)] == [
        Fraction(a, b) for a, b in [(1, 1), (2, 1), (7, 3), (5, 2), (13, 5), (8, 3), (19, 7), (11, 4), (25, 9), (14, 5)]]
    assert round(math.log(walk_count(100, 0)), 4) == 66.7838 and round(ln_binom_approx(100, 0), 4) == 66.7863
    assert all(second_difference(lambda t: t * t, y) == 2 for y in range(-5, 6))
    assert moment_exact(30, 2) == 30 and round(float(walk_prob(30, 0)), 4) == 0.1445
    assert round(2 * heat_kernel(0, 30), 4) == 0.1457
    assert round(sd_exact(10, 0.5), 4) == 0.1762 and round(sd_approx(10, 0.5), 4) == 0.1983
    assert round(math.sqrt(float(mean_square(big))), 3) == 9.946
    assert round(100 * float(walk_prob(10, 0)), 1) == 24.6
    return True


if __name__ == "__main__":
    check()
    import sys
    sys.stdout.reconfigure(encoding="utf-8")
    print("walk1d（第 1 章）：計數、動差、Stirling、高斯近似、遞迴、傅立葉、生成函數、最陡下降全部驗證通過")
