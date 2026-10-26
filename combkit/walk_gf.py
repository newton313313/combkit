# 由 combkit/tools/sync_from_walk017.py 從 walk017 專案的 codes/utils/walk_gf.py 複製；請改 walk017 那份再同步。
"""
walk_gf.py —— 《Elements of the Random Walk》第 2 章 Generating functions I 的數學（純 Python，不 import manim）

格子上的隨機漫步：最近鄰格距 a = 1；C(N; r) 是從原點走 N 步到格點 r 的走法數（每步走到一個最近鄰）。
畫面上、講義裡出現的每個數字都從這裡算，場景與 make_figs.py 再 assert。書上式號與筆誤見 ref/md/ch02_formulas.md。
本模組不依賴 walk1d（一維的二項式計數用 math.comb 重寫），以便原樣複製進 combkit。

§2.1.1 生成函數入門
    exp_coeffs(n)                    e^x 的係數 [1/k!]（Fraction）
    legendre_gf(nmax)                由 (2.1) 1/sqrt(1−2xz+z²) 展開得 P_0..P_nmax（多項式＝係數 list，Fraction，x^0 起）
    legendre_bonnet(nmax)            Bonnet 遞迴 (n+1)P_{n+1} = (2n+1)xP_n − nP_{n−1}（第二種算法）
    half_binom_series(n)             (1−t)^{−1/2} 的係數 [C(2k,k)/4^k]：1, 1/2, 3/8, 5/16, …
    poly_eval(p, x)、poly_mul(p, q)、poly_int11(p)   多項式求值、相乘、∫_{−1}^{1}（精確）
    ex22_integral(z)                 Ex. 2.2：∫_{−1}^{1} (1−2xz+z²)^{−1/2} dx（閉式 (1+z−|1−z|)/z，|z|<1 時 = 2）
    ex23_closed(y, z)                Ex. 2.3：(1/sqrt(yz)) ln((1+sqrt(yz))/(1−sqrt(yz)))；係數 ex23_coeffs(n) = [2/(2k+1)]
    bessel_series(n, rho)、bessel_integral(n, rho)   Ex. 2.4：J_n 的級數與積分表示（互相核對）
    ideal_gas_Z(N, x)                Ex. 2.5：Z(N) = x^N/N!（x = 𝒞T^{3/2}V）；grand_Q(z, x) = e^{zx}；grand_probs(x, z, Nmax)
    survival_prob(n, z) = z^n(1−z)；survival_mean(z) = z/(1−z)
    mortal_walkers(M, z, seed)       M 個「每步以機率 z 存活」的走路人在正方格子上 → [(步數, 終點)]
    mortal_trail(z, seed, nmax)      一個走路人的足跡（格點 list）
    ellipk(k)                        第一類完全橢圓積分 K(k)（AGM）
    mortal_origin_exact(z)           停在原點的機率 (1−z)·G(z/4; 0) = (1−z)(2/π)K(z)；mortal_origin_series 用級數核對

§2.1.2 走法數、遞迴 (2.7)
    SQ_STEPS（E N W S）、SC_STEPS（6）、BCC_STEPS（(±1,±1,±1) 共 8）、FCC_STEPS（(±1,±1,0) 排列共 12）——整數座標
    count1(N, d)                     一維 C(N,d)
    count_grid(N, steps)             由 (2.7) 一步步算出 {格點: C(N;格點)}；count_grids(Nmax, steps) → 每一步的表
    count_sq(N, x, y)                正方格子「兩枚硬幣」公式 C1(N,x+y)·C1(N,x−y)
    returns_sq(N)                    C(N;0) = C(N,N/2)²：1, 4, 36, 400, 4900, 63504（主例）
    returns_sc(N)                    簡單立方：C(2n,n) Σ (n!/(j!k!l!))²：1, 6, 90, 1860, 44730, 1172556
    returns_bcc(N) = C(N,N/2)³；returns_fcc(N)（由 count_grid）
    UV_STEP                          兩枚硬幣 → 一步：(H,H)=E、(H,T)=N、(T,T)=W、(T,H)=S（u = x+y、v = x−y）
    trail_sq(n, seed)                正方格子上 n 步的隨機足跡（Fig. 2.1–2.2 的樣子）

§2.1.2 傅立葉 (2.10)–(2.13)
    chi_1d(q)、chi_sq(qx, qy)、chi_sc(qx, qy, qz)、chi_bcc(…)、chi_fcc(…)   結構函數（最近鄰距離 a = 1）
    chi_from_steps(q, steps, scale)  Σ e^{iq·w}（直接加總，核對閉式）
    CHI_SMALL_Q                      χ ≈ Z − c2 q²：{1D: (2, 1), 方: (4, 1), sc: (6, 1), bcc: (8, 4/3), fcc: (12, 2)}
    g_fourier(z, chi)                (2.13) 1/(1 − zχ)
    plane_wave_check(q)              cos(q(x−1)) + cos(q(x+1)) = 2cos q · cos(qx)（鄰居的和＝乘上 χ）

§2.1.3 取係數與高斯極限 (2.14)–(2.20)
    fourier_count_sq(N, x, y)        (2.17)：(1/(2π)²)∫∫ e^{−iq·r} χ^N（均分求積，精確）
    fourier_count_sc(N, r)           三維版
    chi2_terms()                     χ² 展開的 16 項（步對），常數項 4 個
    gauss_sc(N, r2)                  (2.20)：6^N (3/(2πN))^{3/2} e^{−3r²/(2N)}（書上版，沒有奇偶性的因子 2）
    gauss_sq(N, r2)                  正方格子：4^N/(πN) e^{−r²/N}
    gauss_lattice(N, kind, r2, parity)   一般：Z^N Ω (d/(2πN))^{d/2} e^{−d r²/(2N)} × 奇偶因子（二部圖：1+(−1)^{N+Σr}；fcc：1）
    CELL_VOLUME                      最近鄰距離 1 時的單胞體積：1D/方/sc 1、bcc 4/(3√3)、fcc 1/√2

§2.1.4 普適性與離開格子 (2.21)–(2.31)
    eq222(N, r, a, b)                (2.22) 右邊（指數用正確的 −3r²/(2a²N)）；eq221_first_order(…) 數值積分核對
    coin_quartic(N, d)               一維擲幣：(精確 P, 高斯×2, 高斯×2×(2.22) 的修正)（a² = 3、b = 1/12）
    p_gauss_book(r, r0)、p_gauss(r, r0)   Ex. 2.7：書上印的（∫d³r = 4π）與正確歸一化的高斯步長
    radial_integral(p, r0)           4π∫ r² p(r) dr（數值）
    p_cauchy(r, r0)                  (2.26) r0/(π²(r²+r0²)²)；chi_cauchy(q, r0) = e^{−q r0}（(2.30)）
    chi_cauchy_numeric(q, r0)        (2.27) 第二行的數值積分
    P_cauchy_N(r, N, r0)             (2.31)：N r0/(π²(r²+N²r0²)²)；P_cauchy_N_numeric 由 (2.31) 第一行數值積分
    cauchy_radial_cdf(R, r0)、cauchy_radial_median(r0)   一步的 |r| 分布與中位數
    steps_gauss_3d(n, r0, rng)、steps_cauchy_3d(n, r0, rng)   取樣（Cauchy：r0·Z/|W|）
    trail_3d(steps)；median_distance(kind, N, M, seed)

§2.1.5 擴散 (2.32)–(2.38)
    prob_grid_sq(N)                  P(N; r) = C(N; r)/4^N（Fraction）；diffuse_sq(P) 一步 (2.32)
    variance(P)                      Σ r² P
    ex28_sides(N, r)                 Ex. 2.8：∂P/∂N 與 (1/6)∇²P（a = 1，解析式）
    g_cont(z, q, d)、g_cont_book(z, q, d)   (2.38) 的傅立葉版：正確 1/(1−2dz+zq²)；書上的號 1/(1−2dz−zq²)
    G1_exact(z, x)、G1_continuum(z, x)      一維 G(z; x) 的精確值 λ^{|x|}/sqrt(1−4z²) 與連續極限 e^{−κ|x|}/(2zκ)

§2.2 補充 1：高斯積分
    gauss_moment(n, a)、gauss_moment_numeric(n, a)   (2:S1-6) 與數值積分
    gauss_gen(a, b)                  (2:S1-5) e^{b²/4a} sqrt(π/a)
    square_complete(a, b, x, book)   (2:S1-4) 右邊（book=True 用書上的 −b²/2a）

§2.3 補充 2：格子上的傅立葉
    SC_BASIS、BCC_BASIS、FCC_BASIS   (2:S2-9)–(2:S2-20)（最近鄰距離 1）；reciprocal(a1, a2, a3) → b1, b2, b3（(2:S2-7)）
    dot、cross、triple
    roots_sum(m, n)                  Σ_{k<n} e^{2πimk/n}（(2:S2-33)/(2:S2-34)）

§2.4 補充 3：冪級數係數的漸近
    series_coeffs_power(alpha, zc, nmax)   (zc − z)^{−α} 的係數（精確遞推）
    coeff_power_exact(N, alpha, zc)  (2:S3-8)；coeff_power_asym(N, alpha, zc)  (2:S3-10)
    coeff_two_poles(N, a, b, zc, zd)  (2:S3-4)
    coeff_log_exact(N, alpha, zc)    −(zc−z)^{−α} ln(zc−z) 的係數（級數相乘）
    coeff_log_asym(N, alpha, zc, book=False, terms=3)   (2:S3-11)：正確版第一項 −ln zc（書上 +ln zc）
    digamma(x)
    asym_returns_sq(N)               2·4^N/(πN)（偶數 N）；asym_returns_1d(N) = 2^N sqrt(2/(πN))

check() 驗證全部（含書上的數字、筆誤的更正與 plan_ch02.md 場景卡上的每個數字）；直接執行本檔會跑 check()。
"""
from __future__ import annotations

import cmath
import math
import random
from fractions import Fraction

# =========================================================================== §2.1.1 生成函數入門


def exp_coeffs(n: int) -> list[Fraction]:
    return [Fraction(1, math.factorial(k)) for k in range(n + 1)]


def half_binom_series(n: int) -> list[Fraction]:
    """(1 − t)^{−1/2} = Σ C(2k,k)/4^k t^k。"""
    return [Fraction(math.comb(2 * k, k), 4 ** k) for k in range(n + 1)]


def legendre_gf(nmax: int) -> list[list[Fraction]]:
    """把 (1 − t)^{−1/2}、t = 2xz − z² 展開，收集 z^n 的係數（x 的多項式）。
    t^k = Σ_j C(k,j)(2xz)^{k−j}(−z²)^j → z^{k+j} x^{k−j}；z^n x^m 來自 k = (n+m)/2、j = (n−m)/2。"""
    out = []
    for n in range(nmax + 1):
        p = [Fraction(0)] * (n + 1)
        for m in range(n % 2, n + 1, 2):
            k, j = (n + m) // 2, (n - m) // 2
            p[m] = Fraction(math.comb(2 * k, k), 4 ** k) * math.comb(k, j) * 2 ** m * (-1) ** j
        out.append(p)
    return out


def poly_mul(p, q):
    out = [Fraction(0)] * (len(p) + len(q) - 1)
    for i, a in enumerate(p):
        for j, b in enumerate(q):
            out[i + j] += a * b
    return out


def poly_eval(p, x):
    s = 0
    for c in reversed(p):
        s = s * x + c
    return s


def poly_int11(p) -> Fraction:
    """∫_{−1}^{1} p(x) dx（精確）。"""
    return sum(Fraction(2 * c, k + 1) for k, c in enumerate(p) if k % 2 == 0)


def legendre_bonnet(nmax: int) -> list[list[Fraction]]:
    P = [[Fraction(1)], [Fraction(0), Fraction(1)]]
    for n in range(1, nmax):
        xp = [Fraction(0)] + P[n]
        nxt = [Fraction(0)] * (n + 2)
        for k, c in enumerate(xp):
            nxt[k] += Fraction(2 * n + 1, n + 1) * c
        for k, c in enumerate(P[n - 1]):
            nxt[k] -= Fraction(n, n + 1) * c
        P.append(nxt)
    return P[: nmax + 1]


def ex22_integral(z: float) -> float:
    """∫_{−1}^{1} dx/sqrt(1 − 2xz + z²)，閉式：[−sqrt(1−2xz+z²)/z]_{−1}^{1} = (1 + z − |1 − z|)/z。"""
    return (1 + z - abs(1 - z)) / z


def ex23_closed(y: float, z: float) -> float:
    s = math.sqrt(y * z)
    return math.log((1 + s) / (1 - s)) / s


def ex23_coeffs(n: int) -> list[Fraction]:
    """ex23_closed = Σ_k 2/(2k+1) (yz)^k → ∫P_n P_m dx = 2δ_{nm}/(2n+1)。"""
    return [Fraction(2, 2 * k + 1) for k in range(n + 1)]


def _simpson(f, a: float, b: float, n: int = 2000) -> float:
    if n % 2:
        n += 1
    h = (b - a) / n
    s = f(a) + f(b)
    for k in range(1, n):
        s += (4 if k % 2 else 2) * f(a + k * h)
    return s * h / 3


def bessel_series(n: int, rho: float, terms: int = 40) -> float:
    n0 = abs(n)
    s = sum((-1) ** k * (rho / 2) ** (2 * k + n0) / (math.factorial(k) * math.factorial(k + n0)) for k in range(terms))
    return s * ((-1) ** n0 if n < 0 else 1)


def bessel_integral(n: int, rho: float) -> float:
    """Ex. 2.4：J_n(ρ) = (1/2π)∫_{−π}^{π} e^{i(ρ sinθ − nθ)} dθ（取實部；虛部為 0）。"""
    return _simpson(lambda t: math.cos(rho * math.sin(t) - n * t), -math.pi, math.pi, 400) / (2 * math.pi)


def ideal_gas_Z(N: int, x) -> Fraction:
    return Fraction(x) ** N / math.factorial(N)


def grand_Q(z: float, x: float) -> float:
    return math.exp(z * x)


def grand_probs(x: float, z: float, Nmax: int) -> list[float]:
    """巨正則系綜裡粒子數 N 的機率 Z(N) z^N / Q（理想氣體 → Poisson(zx)）。"""
    Q = grand_Q(z, x)
    return [float(ideal_gas_Z(N, Fraction(x))) * z ** N / Q for N in range(Nmax + 1)]


def survival_prob(n: int, z):
    return z ** n * (1 - z)


def survival_mean(z):
    return z / (1 - z)


SQ_STEPS = [(1, 0), (0, 1), (-1, 0), (0, -1)]          # E N W S
SQ_NAMES = "ENWS"


def mortal_walkers(M: int, z: float, seed: int) -> list[tuple[int, tuple[int, int]]]:
    rng = random.Random(seed)
    out = []
    for _ in range(M):
        x = y = n = 0
        while rng.random() < z:
            dx, dy = SQ_STEPS[rng.randrange(4)]
            x, y, n = x + dx, y + dy, n + 1
        out.append((n, (x, y)))
    return out


def mortal_trail(z: float, seed: int, nmax: int = 10 ** 6) -> list[tuple[int, int]]:
    rng = random.Random(seed)
    pts = [(0, 0)]
    while rng.random() < z and len(pts) <= nmax:
        dx, dy = SQ_STEPS[rng.randrange(4)]
        pts.append((pts[-1][0] + dx, pts[-1][1] + dy))
    return pts


def ellipk(k: float) -> float:
    """K(k) = π / (2·AGM(1, sqrt(1−k²)))。"""
    a, g = 1.0, math.sqrt(1 - k * k)
    for _ in range(60):
        a, g = (a + g) / 2, math.sqrt(a * g)
    return math.pi / (2 * a)


def mortal_origin_exact(z: float) -> float:
    """Σ_N z^N(1−z) C(N;0)/4^N = (1−z)·Σ_n C(2n,n)² (z/4)^{2n} = (1−z)(2/π)K(z)。"""
    return (1 - z) * 2 / math.pi * ellipk(z)


def mortal_origin_series(z: float, nmax: int = 4000) -> float:
    s, term = 0.0, 1.0                     # term = C(2n,n)² (z/4)^{2n}，用比值遞推
    for n in range(nmax):
        s += term
        term *= ((2 * n + 1) * (2 * n + 2) / ((n + 1) ** 2)) ** 2 * (z / 4) ** 2
    return (1 - z) * s


# =========================================================================== §2.1.2 走法數與遞迴

SC_STEPS = [(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)]
BCC_STEPS = [(i, j, k) for i in (1, -1) for j in (1, -1) for k in (1, -1)]
FCC_STEPS = [s for s in {(i, j, 0) for i in (1, -1) for j in (1, -1)}
             | {(i, 0, k) for i in (1, -1) for k in (1, -1)}
             | {(0, j, k) for j in (1, -1) for k in (1, -1)}]
FCC_STEPS.sort()
STEPS = {"1d": [(1,), (-1,)], "sq": SQ_STEPS, "sc": SC_STEPS, "bcc": BCC_STEPS, "fcc": FCC_STEPS}
UV_STEP = {("H", "H"): (1, 0), ("H", "T"): (0, 1), ("T", "T"): (-1, 0), ("T", "H"): (0, -1)}


def count1(N: int, d: int) -> int:
    if abs(d) > N or (N + d) % 2:
        return 0
    return math.comb(N, (N + d) // 2)


def _add(p, s):
    return tuple(a + b for a, b in zip(p, s))


def count_grid(N: int, steps=SQ_STEPS) -> dict:
    """(2.7)：C(N; y) = Σ_i C(N−1; y − s_i)（最後一步來自鄰居），C(0; y) = δ_{y,0}。"""
    dim = len(steps[0])
    cur = {(0,) * dim: 1}
    for _ in range(N):
        nxt = {}
        for p, c in cur.items():
            for s in steps:
                q = _add(p, s)
                nxt[q] = nxt.get(q, 0) + c
        cur = nxt
    return cur


def count_grids(Nmax: int, steps=SQ_STEPS) -> list[dict]:
    return [count_grid(N, steps) for N in range(Nmax + 1)]


def count_sq(N: int, x: int, y: int) -> int:
    return count1(N, x + y) * count1(N, x - y)


def returns_sq(N: int) -> int:
    return count1(N, 0) ** 2


def returns_sc(N: int) -> int:
    if N % 2:
        return 0
    n = N // 2
    s = sum((math.factorial(n) // (math.factorial(j) * math.factorial(k) * math.factorial(n - j - k))) ** 2
            for j in range(n + 1) for k in range(n + 1 - j))
    return math.comb(2 * n, n) * s


def returns_bcc(N: int) -> int:
    return count1(N, 0) ** 3


def returns_fcc(N: int) -> int:
    return count_grid(N, FCC_STEPS).get((0, 0, 0), 0)


def uv_walk(ucoins: str, vcoins: str) -> list[tuple[int, int]]:
    pts = [(0, 0)]
    for cu, cv in zip(ucoins, vcoins):
        pts.append(_add(pts[-1], UV_STEP[(cu, cv)]))
    return pts


def trail_sq(n: int, seed: int) -> list[tuple[int, int]]:
    rng = random.Random(seed)
    pts = [(0, 0)]
    for _ in range(n):
        pts.append(_add(pts[-1], SQ_STEPS[rng.randrange(4)]))
    return pts


# =========================================================================== §2.1.2 傅立葉

R3, R2 = math.sqrt(3), math.sqrt(2)


def chi_1d(q: float) -> float:
    return 2 * math.cos(q)


def chi_sq(qx: float, qy: float) -> float:
    return 2 * (math.cos(qx) + math.cos(qy))


def chi_sc(qx: float, qy: float, qz: float) -> float:
    return 2 * (math.cos(qx) + math.cos(qy) + math.cos(qz))


def chi_bcc(qx: float, qy: float, qz: float) -> float:
    return 8 * math.cos(qx / R3) * math.cos(qy / R3) * math.cos(qz / R3)


def chi_fcc(qx: float, qy: float, qz: float) -> float:
    cx, cy, cz = math.cos(qx / R2), math.cos(qy / R2), math.cos(qz / R2)
    return 4 * (cx * cy + cy * cz + cz * cx)


def chi_from_steps(q, steps, scale: float = 1.0) -> float:
    """Σ_i e^{iq·w_i}，w_i = scale·steps[i]（對稱的鄰居 → 實數）。"""
    s = sum(cmath.exp(1j * scale * sum(a * b for a, b in zip(q, w))) for w in steps)
    return s.real


CHI_SMALL_Q = {"1d": (2, Fraction(1)), "sq": (4, Fraction(1)), "sc": (6, Fraction(1)),
               "bcc": (8, Fraction(4, 3)), "fcc": (12, Fraction(2))}
DIM = {"1d": 1, "sq": 2, "sc": 3, "bcc": 3, "fcc": 3}
CELL_VOLUME = {"1d": 1.0, "sq": 1.0, "sc": 1.0, "bcc": 4 / (3 * R3), "fcc": 1 / R2}
BIPARTITE = {"1d": True, "sq": True, "sc": True, "bcc": True, "fcc": False}
STEP_SCALE = {"1d": 1.0, "sq": 1.0, "sc": 1.0, "bcc": 1 / R3, "fcc": 1 / R2}   # 整數座標 → 最近鄰距離 1


def g_fourier(z: float, chi: float) -> float:
    return 1 / (1 - z * chi)


def plane_wave_check(q: float, xs=range(-3, 4)) -> float:
    """max |cos(q(x−1)) + cos(q(x+1)) − 2cos q·cos(qx)|。"""
    return max(abs(math.cos(q * (x - 1)) + math.cos(q * (x + 1)) - 2 * math.cos(q) * math.cos(q * x)) for x in xs)


# =========================================================================== §2.1.3 取係數與高斯極限


def fourier_count_sq(N: int, x: int, y: int, M: int | None = None) -> float:
    """(2.17)（二維、a = 1）：均分 M×M 點求積；M > N + |x| + |y| 時精確（三角多項式）。"""
    M = M or N + abs(x) + abs(y) + 1
    s = 0.0
    for j in range(M):
        qx = -math.pi + 2 * math.pi * j / M
        for k in range(M):
            qy = -math.pi + 2 * math.pi * k / M
            s += (cmath.exp(-1j * (qx * x + qy * y)) * chi_sq(qx, qy) ** N).real
    return s / M ** 2


def fourier_count_sc(N: int, r=(0, 0, 0), M: int | None = None) -> float:
    M = M or N + sum(abs(t) for t in r) + 1
    qs = [-math.pi + 2 * math.pi * j / M for j in range(M)]
    cs = [math.cos(q) for q in qs]
    s = 0.0
    for a, qx in enumerate(qs):
        for b, qy in enumerate(qs):
            for c, qz in enumerate(qs):
                s += math.cos(qx * r[0] + qy * r[1] + qz * r[2]) * (2 * (cs[a] + cs[b] + cs[c])) ** N
    return s / M ** 3


def chi2_terms() -> list[tuple[str, str, bool]]:
    """χ² = (Σ_i e^{iq·w_i})² 的 16 項：(第一步, 第二步, 是否為常數項＝回到原點)。"""
    return [(SQ_NAMES[i], SQ_NAMES[j], _add(SQ_STEPS[i], SQ_STEPS[j]) == (0, 0))
            for i in range(4) for j in range(4)]


def gauss_sc(N: float, r2: float = 0.0) -> float:
    """(2.20)，a = 1（書上版）。"""
    return 6.0 ** N * (3 / (2 * math.pi * N)) ** 1.5 * math.exp(-3 * r2 / (2 * N))


def gauss_sq(N: float, r2: float = 0.0) -> float:
    return 4.0 ** N / (math.pi * N) * math.exp(-r2 / N)


def gauss_lattice(N: int, kind: str, r2: float = 0.0, parity: int = 0) -> float:
    """Z^N Ω (d/(2πN))^{d/2} e^{−d r²/(2N)} × 奇偶因子。r2 用最近鄰距離 1 的單位；
    parity = N + Σ(整數座標) 的奇偶（二部圖的格子：偶 → 2、奇 → 0；fcc 不是二部圖 → 1）。"""
    Z = CHI_SMALL_Q[kind][0]
    d = DIM[kind]
    base = Z ** N * CELL_VOLUME[kind] * (d / (2 * math.pi * N)) ** (d / 2) * math.exp(-d * r2 / (2 * N))
    if BIPARTITE[kind]:
        return base * (2 if parity % 2 == 0 else 0)
    return base


def asym_returns_sq(N: int) -> float:
    return 2 * 4.0 ** N / (math.pi * N) if N % 2 == 0 else 0.0


def asym_returns_1d(N: int) -> float:
    return 2.0 ** N * math.sqrt(2 / (math.pi * N)) if N % 2 == 0 else 0.0


# =========================================================================== §2.1.4 普適性與離開格子


def eq222(N: float, r: float, a: float, b: float) -> float:
    """(2.22) 右邊（指數 −3r²/(2a²N)；書上印成 2aN）。"""
    u = r * r / (N * a * a)
    return math.sqrt(6 * math.pi / (N * a * a)) * math.exp(-1.5 * u) * (1 - 27 * b / (N * a ** 4) * (1 - 6 * u + 3 * u * u))


def eq221_first_order(N: float, r: float, a: float, b: float) -> float:
    """∫ e^{−iqr − a²Nq²/6}(1 − bNq⁴) dq 的數值積分（核對 (2.22)）。"""
    al = a * a * N / 6
    L = 12 / math.sqrt(al)
    return _simpson(lambda q: math.cos(q * r) * math.exp(-al * q * q) * (1 - b * N * q ** 4), -L, L, 4000)


def coin_quartic(N: int, d: int) -> tuple[float, float, float]:
    """一維擲幣：cos^N q = e^{−Nq²/2 − Nq⁴/12 − …}，即 (2.21) 的 a² = 3、b = 1/12；兩個鞍點 → 乘 2/(2π)。"""
    exact = count1(N, d) / 2 ** N
    gauss = 2 * math.exp(-d * d / (2 * N)) / math.sqrt(2 * math.pi * N)
    corr = 2 / (2 * math.pi) * eq222(N, d, math.sqrt(3), 1 / 12)
    return exact, gauss, corr


def p_gauss_book(r: float, r0: float = 1.0) -> float:
    return math.sqrt(2 / math.pi) * r0 ** -3 * math.exp(-r * r / (2 * r0 * r0))


def p_gauss(r: float, r0: float = 1.0) -> float:
    return (2 * math.pi * r0 * r0) ** -1.5 * math.exp(-r * r / (2 * r0 * r0))


def radial_integral(p, r0: float = 1.0, R: float = 60.0, n: int = 20000) -> float:
    return 4 * math.pi * _simpson(lambda r: r * r * p(r, r0), 0.0, R * r0, n)


def p_cauchy(r: float, r0: float = 1.0) -> float:
    return r0 / (math.pi ** 2 * (r * r + r0 * r0) ** 2)


def chi_cauchy(q: float, r0: float = 1.0) -> float:
    return math.exp(-q * r0)


def chi_cauchy_numeric(q: float, r0: float = 1.0, R: float = 400.0, n: int = 200000) -> float:
    """(2.27) 第二行：(4/π)∫_0^∞ sin(qr) r r0 / (q (r²+r0²)²) dr。"""
    f = lambda r: math.sin(q * r) * r * r0 / (q * (r * r + r0 * r0) ** 2)
    return 4 / math.pi * _simpson(f, 0.0, R, n)


def P_cauchy_N(r: float, N: int, r0: float = 1.0) -> float:
    return N * r0 / (math.pi ** 2 * (r * r + N * N * r0 * r0) ** 2)


def P_cauchy_N_numeric(r: float, N: int, r0: float = 1.0) -> float:
    """(2.31) 第一行：(4π/(2π)³)∫_0^∞ e^{−Nqr0} q sin(qr)/r dq。"""
    L = 60 / (N * r0)
    return 4 * math.pi / (2 * math.pi) ** 3 * _simpson(lambda q: math.exp(-N * q * r0) * q * math.sin(q * r) / r, 0.0, L, 6000)


def cauchy_radial_cdf(R: float, r0: float = 1.0) -> float:
    u = R / r0
    return 2 / math.pi * (math.atan(u) - u / (1 + u * u))


def cauchy_radial_median(r0: float = 1.0) -> float:
    lo, hi = 0.0, 100.0
    for _ in range(200):
        mid = (lo + hi) / 2
        lo, hi = (mid, hi) if cauchy_radial_cdf(mid) < 0.5 else (lo, mid)
    return lo * r0


CHI3_MEDIAN = 1.5381722544550522                     # 自由度 3 的 χ 分布中位數（核對見 check）


def steps_gauss_3d(n: int, r0: float, rng: random.Random) -> list[tuple[float, float, float]]:
    return [(rng.gauss(0, r0), rng.gauss(0, r0), rng.gauss(0, r0)) for _ in range(n)]


def steps_cauchy_3d(n: int, r0: float, rng: random.Random) -> list[tuple[float, float, float]]:
    """三維各向同性 Cauchy（多變量 t，自由度 1）：r0·Z/|W|，Z ~ N(0, I₃)、W ~ N(0,1)；密度正是 (2.26)。"""
    out = []
    for _ in range(n):
        w = abs(rng.gauss(0, 1))
        out.append(tuple(r0 * rng.gauss(0, 1) / w for _ in range(3)))
    return out


def trail_3d(steps) -> list[tuple[float, float, float]]:
    pts = [(0.0, 0.0, 0.0)]
    for s in steps:
        pts.append(_add(pts[-1], s))
    return pts


def median_distance(kind: str, N: int, M: int, seed: int, r0: float = 1.0) -> float:
    rng = random.Random(seed)
    draw = steps_gauss_3d if kind == "gauss" else steps_cauchy_3d
    ds = sorted(math.dist((0, 0, 0), trail_3d(draw(N, r0, rng))[-1]) for _ in range(M))
    return (ds[M // 2 - 1] + ds[M // 2]) / 2 if M % 2 == 0 else ds[M // 2]


# =========================================================================== §2.1.5 擴散


def prob_grid_sq(N: int) -> dict:
    return {p: Fraction(c, 4 ** N) for p, c in count_grid(N, SQ_STEPS).items()}


def diffuse_sq(P: dict) -> dict:
    """(2.32)，d = 2：P(N; y) = ¼ Σ_nn P(N−1; y − w_i)。"""
    out = {}
    for p, v in P.items():
        for s in SQ_STEPS:
            q = _add(p, s)
            out[q] = out.get(q, 0) + v / 4
    return out


def variance(P: dict):
    return sum(v * sum(t * t for t in p) for p, v in P.items())


def ex28_sides(N: float, r: float) -> tuple[float, float]:
    """P = (3/(2πN))^{3/2} e^{−3r²/(2N)}（a = 1）：∂P/∂N 與 (1/6)∇²P（徑向拉普拉斯 P'' + 2P'/r）。"""
    P = (3 / (2 * math.pi * N)) ** 1.5 * math.exp(-3 * r * r / (2 * N))
    dN = (-1.5 / N + 1.5 * r * r / (N * N)) * P
    lap = (9 * r * r / (N * N) - 9 / N) * P              # (∂²_r + (2/r)∂_r) P
    return dN, lap / 6


def g_cont(z: float, q: float, d: int = 1) -> float:
    return 1 / (1 - 2 * d * z + z * q * q)


def g_cont_book(z: float, q: float, d: int = 1) -> float:
    return 1 / (1 - 2 * d * z - z * q * q)


def G1_exact(z: float, x: int) -> float:
    s = math.sqrt(1 - 4 * z * z)
    lam = (1 - s) / (2 * z)
    return lam ** abs(x) / s


def G1_fourier(z: float, x: int, n: int = 4000) -> float:
    return _simpson(lambda q: math.cos(q * x) / (1 - 2 * z * math.cos(q)), -math.pi, math.pi, n) / (2 * math.pi)


def G1_continuum(z: float, x: float) -> float:
    kappa = math.sqrt((1 - 2 * z) / z)
    return math.exp(-kappa * abs(x)) / (2 * z * kappa)


# =========================================================================== §2.2 高斯積分


def gauss_moment(n: int, a: float = 1.0) -> float:
    if n % 2:
        return 0.0
    return math.sqrt(math.pi / a) * math.factorial(n) / (math.factorial(n // 2) * (4 * a) ** (n // 2))


def gauss_moment_frac(n: int, a=1) -> Fraction:
    """(2:S1-6) 去掉 sqrt(π/a) 的有理部分 n!/((n/2)!(4a)^{n/2})。"""
    return Fraction(math.factorial(n), math.factorial(n // 2) * (4 * Fraction(a)) ** (n // 2)) if n % 2 == 0 else Fraction(0)


def gauss_moment_numeric(n: int, a: float = 1.0) -> float:
    L = 12 / math.sqrt(a)
    return _simpson(lambda x: x ** n * math.exp(-a * x * x), -L, L, 4000)


def gauss_gen(a: float, b: float) -> float:
    return math.exp(b * b / (4 * a)) * math.sqrt(math.pi / a)


def square_complete(a: float, b: float, x: float, book: bool = False) -> float:
    return a * (x - b / (2 * a)) ** 2 - b * b / ((2 if book else 4) * a)


# =========================================================================== §2.3 格子上的傅立葉


def dot(u, v):
    return sum(a * b for a, b in zip(u, v))


def cross(u, v):
    return (u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0])


def triple(a1, a2, a3):
    return dot(a1, cross(a2, a3))


def reciprocal(a1, a2, a3):
    om = triple(a1, a2, a3)
    k = 2 * math.pi / om
    return tuple(tuple(k * t for t in cross(u, v)) for u, v in ((a2, a3), (a3, a1), (a1, a2)))


SC_BASIS = ((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0))
BCC_BASIS = tuple(tuple(t / R3 for t in v) for v in ((1, 1, 1), (1, 1, -1), (1, -1, 1)))
FCC_BASIS = tuple(tuple(t / R2 for t in v) for v in ((1, 1, 0), (1, 0, 1), (0, 1, 1)))


def roots_sum(m: int, n: int) -> complex:
    return sum(cmath.exp(2j * math.pi * m * k / n) for k in range(n))


# =========================================================================== §2.4 冪級數係數的漸近


def series_coeffs_power(alpha: float, zc: float, nmax: int) -> list[float]:
    """(zc − z)^{−α} = zc^{−α} Σ_n [α(α+1)…(α+n−1)/n!] (z/zc)^n。"""
    out, c = [], zc ** -alpha
    for n in range(nmax + 1):
        out.append(c)
        c *= (alpha + n) / ((n + 1) * zc)
    return out


def coeff_power_exact(N: int, alpha: float, zc: float) -> float:
    return math.exp(math.lgamma(N + alpha) - math.lgamma(alpha) - math.lgamma(N + 1)) * zc ** (-N - alpha)


def coeff_power_asym(N: int, alpha: float, zc: float) -> float:
    return N ** (alpha - 1) * zc ** (-N - alpha) / math.gamma(alpha)


def coeff_two_poles(N: int, a: float, b: float, zc: float, zd: float) -> float:
    return a * zc ** -(N + 1) + b * zd ** -(N + 1)


def coeff_log_exact(N: int, alpha: float, zc: float) -> float:
    """−(zc − z)^{−α} ln(zc − z)；ln(zc − z) = ln zc − Σ_{k≥1} (z/zc)^k / k。"""
    A = series_coeffs_power(alpha, zc, N)
    L = [math.log(zc)] + [-(zc ** -k) / k for k in range(1, N + 1)]
    return -sum(A[j] * L[N - j] for j in range(N + 1))


def digamma(x: float) -> float:
    s = 0.0
    while x < 20:
        s -= 1 / x
        x += 1
    f = 1 / (x * x)
    return s + math.log(x) - 0.5 / x - f * (1 / 12 - f * (1 / 120 - f / 252))


def coeff_log_asym(N: int, alpha: float, zc: float, book: bool = False, terms: int = 3) -> float:
    """(2:S3-11)：d/dα [zc^{−N−α} N^{α−1}/Γ(α)] = [−ln zc − ψ(α) + ln N]·zc^{−N−α}N^{α−1}/Γ(α)；
    書上第一項印成 +ln zc。terms = 1 只留 ln N（第二行）。"""
    base = zc ** (-N - alpha) * N ** (alpha - 1) / math.gamma(alpha)
    if terms == 1:
        return math.log(N) * base
    return ((1 if book else -1) * math.log(zc) - digamma(alpha) + math.log(N)) * base


# =========================================================================== check


def check():
    # ---- §2.1.1 Legendre、e^x、Ex. 2.1–2.5
    F = Fraction
    assert exp_coeffs(4) == [1, 1, F(1, 2), F(1, 6), F(1, 24)]
    assert half_binom_series(3) == [1, F(1, 2), F(3, 8), F(5, 16)]
    L = legendre_gf(6)
    assert L == legendre_bonnet(6)
    assert L[0] == [1] and L[1] == [0, 1] and L[2] == [F(-1, 2), 0, F(3, 2)]                  # Ex. 2.1
    assert L[3] == [0, F(-3, 2), 0, F(5, 2)]
    assert L[5] == [0, F(15, 8), 0, F(-70, 8), 0, F(63, 8)]
    assert all(poly_eval(p, 1) == 1 for p in L)                                                 # x = 1：1/(1−z)
    assert [poly_int11(p) for p in L[:4]] == [2, 0, 0, 0]                                       # Ex. 2.2
    assert [poly_int11(poly_mul(L[n], L[m])) for n in range(4) for m in range(4)] == \
        [F(2, 2 * n + 1) if n == m else 0 for n in range(4) for m in range(4)]                 # Ex. 2.3
    assert ex23_coeffs(3) == [2, F(2, 3), F(2, 5), F(2, 7)]
    for z in (0.1, 0.5, 0.9):
        assert abs(ex22_integral(z) - 2) < 1e-12
        num = _simpson(lambda x: 1 / math.sqrt(1 - 2 * x * z + z * z), -1, 1, 2000)
        assert abs(num - 2) < 1e-6
    y, z = 0.3, 0.5
    num = _simpson(lambda x: 1 / math.sqrt((1 - 2 * x * y + y * y) * (1 - 2 * x * z + z * z)), -1, 1, 2000)
    assert abs(num - ex23_closed(y, z)) < 1e-8
    assert abs(ex23_closed(y, z) - sum(float(c) * (y * z) ** k for k, c in enumerate(ex23_coeffs(40)))) < 1e-12
    for n in (0, 1, 2, -3):                                                                     # Ex. 2.4
        assert abs(bessel_series(n, 2.5) - bessel_integral(n, 2.5)) < 1e-9
    assert round(bessel_series(0, 2.5), 4) == -0.0484
    assert [ideal_gas_Z(N, 2) for N in range(6)] == [1, 2, 2, F(4, 3), F(2, 3), F(4, 15)]          # Ex. 2.5
    gp = grand_probs(2, 1.0, 5)
    assert [round(p, 4) for p in gp] == [0.1353, 0.2707, 0.2707, 0.1804, 0.0902, 0.0361]
    assert round(grand_Q(1.0, 2.0), 3) == 7.389
    assert abs(sum(N * p for N, p in enumerate(grand_probs(2, 1.0, 60))) - 2) < 1e-12              # ⟨N⟩ = zx
    assert [survival_prob(n, F(9, 10)) for n in range(3)] == [F(1, 10), F(9, 100), F(81, 1000)]
    assert survival_mean(F(9, 10)) == 9
    mw = mortal_walkers(2000, 0.9, seed=2)
    assert round(sum(n for n, _ in mw) / 2000, 2) == MORTAL_MEAN
    assert sum(1 for _, p in mw if p == (0, 0)) == MORTAL_AT_ORIGIN
    assert round(mortal_origin_exact(0.9), 4) == 0.1452
    assert abs(mortal_origin_exact(0.9) - mortal_origin_series(0.9)) < 1e-12
    assert len(mortal_trail(0.9, seed=MORTAL_TRAIL_SEED)) - 1 == MORTAL_TRAIL_LEN

    # ---- §2.1.2 走法數（主例）與遞迴
    grids = count_grids(4)
    assert [g.get((0, 0), 0) for g in grids] == [1, 0, 4, 0, 36]
    assert [sum(g.values()) for g in grids] == [1, 4, 16, 64, 256]
    g4 = grids[4]
    assert (g4[(1, 1)], g4[(2, 0)], g4[(2, 2)], g4[(3, 1)], g4[(4, 0)]) == (24, 16, 6, 4, 1)
    assert grids[2] == {(0, 0): 4, (1, 1): 2, (1, -1): 2, (-1, 1): 2, (-1, -1): 2,
                        (2, 0): 1, (-2, 0): 1, (0, 2): 1, (0, -2): 1}
    for N in range(9):                                                    # 兩枚硬幣公式 = 遞迴
        g = count_grid(N)
        assert all(count_sq(N, x, y) == g.get((x, y), 0) for x in range(-N - 1, N + 2) for y in range(-N - 1, N + 2))
    assert [returns_sq(N) for N in range(0, 11, 2)] == [1, 4, 36, 400, 4900, 63504]
    assert returns_sq(10) == 252 ** 2 and count_sq(10, 1, 1) == 52920
    assert [returns_sc(N) for N in range(0, 11, 2)] == [1, 6, 90, 1860, 44730, 1172556]
    assert all(returns_sc(N) == count_grid(N, SC_STEPS).get((0, 0, 0), 0) for N in range(9))
    assert all(returns_bcc(N) == count_grid(N, BCC_STEPS).get((0, 0, 0), 0) for N in range(9))
    assert [returns_bcc(N) for N in range(0, 11, 2)] == [1, 8, 216, 8000, 343000, 16003008]
    assert [returns_fcc(N) for N in range(7)] == [1, 0, 12, 48, 540, 4320, 42240]
    assert len(FCC_STEPS) == 12 and len(BCC_STEPS) == 8
    assert uv_walk("HTHHTHTTHH", "HHTTHTHTTH")[-1] == (1, 1)

    # ---- §2.1.2 傅立葉、結構函數、Ex. 2.6
    for q in ((0.3, -0.7, 1.1), (2.0, 0.5, -2.9)):
        assert abs(chi_sc(*q) - chi_from_steps(q, SC_STEPS)) < 1e-12
        assert abs(chi_bcc(*q) - chi_from_steps(q, BCC_STEPS, 1 / R3)) < 1e-12
        assert abs(chi_fcc(*q) - chi_from_steps(q, FCC_STEPS, 1 / R2)) < 1e-12
        assert abs(chi_sq(q[0], q[1]) - chi_from_steps(q[:2], SQ_STEPS)) < 1e-12
    assert (chi_1d(0), chi_sq(0, 0), chi_sc(0, 0, 0), chi_bcc(0, 0, 0), chi_fcc(0, 0, 0)) == (2, 4, 6, 8, 12)
    assert chi_sq(math.pi, math.pi) == -4 and chi_sc(math.pi, math.pi, math.pi) == -6
    h = 1e-3                                                               # χ ≈ Z − c2 q²（任意方向）
    for kind, f in (("sc", chi_sc), ("bcc", chi_bcc), ("fcc", chi_fcc)):
        Z, c2 = CHI_SMALL_Q[kind]
        for u in ((1, 0, 0), (1 / R3, 1 / R3, 1 / R3), (0.6, 0.8, 0)):
            assert abs((Z - f(*(h * t for t in u))) / h ** 2 - float(c2)) < 1e-4, kind
        assert abs(float(c2) / Z - 1 / 6) < 1e-15                         # 正規化後都是 1 − q²/6
    assert abs((4 - chi_sq(h, 0)) / h ** 2 - 1) < 1e-4
    assert plane_wave_check(math.pi / 3) < 1e-12 and abs(chi_1d(math.pi / 3) - 1) < 1e-12
    assert [round(g_fourier(0.1, 4) * 0 + 4 ** N) for N in range(5)] == [1, 4, 16, 64, 256]
    assert abs(g_fourier(0.1, 4) - sum(4 ** N * 0.1 ** N for N in range(200))) < 1e-12

    # ---- §2.1.3 取係數：傅立葉逆變換 = 數格子
    t = chi2_terms()
    assert len(t) == 16 and sum(c for *_, c in t) == 4
    assert [c for a, b, c in t if c] == [True] * 4 and {a + b for a, b, c in t if c} == {"EW", "NS", "WE", "SN"}
    for N in range(0, 11, 2):
        assert round(fourier_count_sq(N, 0, 0)) == returns_sq(N)
        assert abs(fourier_count_sq(N, 0, 0) - returns_sq(N)) < 1e-6 * max(1, returns_sq(N))
    assert round(fourier_count_sq(4, 2, 0)) == 16 and round(fourier_count_sq(10, 1, 1)) == 52920
    assert round(fourier_count_sq(5, 1, 0)) == count_sq(5, 1, 0)
    for N in range(0, 11, 2):
        assert round(fourier_count_sc(N)) == returns_sc(N)

    # ---- §2.1.3 高斯極限 (2.18)–(2.20)；奇偶性的因子 2
    for q in (0.05, 0.1):
        assert abs(math.log((math.cos(q) + 2) / 3) - (-q * q / 6)) < q ** 4
    assert (round(gauss_sc(10)), round(2 * gauss_sc(10))) == (630848, 1261696)
    assert round(returns_sc(20) / (2 * gauss_sc(20)), 4) == 0.9635
    assert abs(gauss_lattice(10, "sc") - 2 * gauss_sc(10)) < 1e-6
    assert round(returns_sc(10) / (2 * gauss_sc(10)), 3) == 0.929
    assert round(returns_sc(10) / gauss_sc(10), 2) == 1.86
    assert round(gauss_sq(10)) == 33377 and round(asym_returns_sq(10)) == 66754
    assert [round(returns_sq(N) / asym_returns_sq(N), 4) for N in (10, 20, 50, 100)] == SQ_RATIOS
    assert abs(gauss_lattice(10, "sq") - asym_returns_sq(10)) < 1e-6
    assert gauss_lattice(9, "sq", 1.0, parity=9 + 1) > 0 and gauss_lattice(9, "sq", 0.0, parity=9) == 0
    assert round(math.cos(math.pi) ** 10) == 1                              # 角落 (π,π,π)：(−1)^10 = 1

    # ---- §2.1.4 普適性：三種立方晶格、(2.22) 的 q⁴ 修正、Cauchy
    for N in (2, 4, 6, 8, 10):
        P = {"sc": returns_sc(N) / 6 ** N, "bcc": returns_bcc(N) / 8 ** N, "fcc": returns_fcc(N) / 12 ** N}
        for kind, p in P.items():
            assert round(p, 5) == UNIV_TABLE[kind][N], (kind, N)
    lim = (3 / (2 * math.pi)) ** 1.5
    assert round(lim, 4) == 0.3299
    for kind in ("sc", "bcc", "fcc"):                                      # N^{3/2}P/(gΩ) → (3/2π)^{3/2}
        p = (returns_sc if kind == "sc" else returns_bcc if kind == "bcc" else returns_fcc)(10) / CHI_SMALL_Q[kind][0] ** 10
        g = 2 if BIPARTITE[kind] else 1
        assert round(10 ** 1.5 * p / (g * CELL_VOLUME[kind]), 3) == UNIV_SCALED10[kind]
    for N, r in ((10, 0.0), (10, 2.0), (30, 3.0)):
        a, b = math.sqrt(3), 1 / 12
        assert abs(eq222(N, r, a, b) - eq221_first_order(N, r, a, b)) < 1e-9
    ex, ga, co = coin_quartic(10, 0)
    assert (round(ex, 6), round(ga, 6), round(co, 6)) == (0.246094, 0.252313, 0.246005)
    ex, ga, co = coin_quartic(10, 2)
    assert (round(ex, 6), round(ga, 6), round(co, 6)) == (0.205078, 0.206577, 0.205268)
    ex, ga, co = coin_quartic(20, 0)
    assert (round(ex, 6), round(ga, 6), round(co, 6)) == (0.176197, 0.178412, 0.176182)
    book_exp = math.exp(-3 * 4 / (2 * math.sqrt(3) * 10))                  # 書上 2aN：r = 2 時的錯誤指數
    assert round(book_exp, 4) == 0.7072 and round(math.exp(-3 * 4 / (2 * 3 * 10)), 4) == 0.8187
    assert abs(radial_integral(p_gauss_book) - 4 * math.pi) < 1e-9          # Ex. 2.7 筆誤：積分是 4π
    assert abs(radial_integral(p_gauss) - 1) < 1e-9
    assert abs(radial_integral(p_cauchy, 1.0, 4000.0, 400000) - 1) < 1e-3   # (2.26) 歸一化（尾巴 ~1/R）
    for q in (0.5, 1.0, 2.0):
        assert abs(chi_cauchy_numeric(q) - chi_cauchy(q)) < 2e-5
    assert round(chi_cauchy(1.0), 6) == 0.367879
    for N, r in ((1, 0.5), (3, 2.0), (10, 7.0)):
        assert abs(P_cauchy_N_numeric(r, N) - P_cauchy_N(r, N)) < 1e-9
    um = cauchy_radial_median()
    assert round(um, 4) == CAUCHY_MEDIAN and abs(cauchy_radial_cdf(um) - 0.5) < 1e-12
    chi3_cdf = lambda x: math.erf(x / R2) - math.sqrt(2 / math.pi) * x * math.exp(-x * x / 2)
    assert abs(chi3_cdf(CHI3_MEDIAN) - 0.5) < 1e-12
    for kind, N, val in MEDIANS:
        assert round(median_distance(kind, N, 400, seed=5), 2) == val, (kind, N)

    # ---- §2.1.5 擴散 (2.32)–(2.38)、Ex. 2.8
    P = {(0, 0): Fraction(1)}
    for N in range(1, 21):
        P = diffuse_sq(P)
        assert variance(P) == N
    assert P == prob_grid_sq(20)
    assert round(float(P[(0, 0)]), 6) == 0.031045 and round(2 / (math.pi * 20), 6) == 0.031831
    assert P.get((1, 0), 0) == 0                                          # 棋盤格：奇偶性
    for N, r in ((5.0, 1.0), (10.0, 3.0), (40.0, 0.5)):
        a, b = ex28_sides(N, r)
        assert abs(a - b) < 1e-15
    nm = _simpson(lambda r: 4 * math.pi * r * r * (3 / (2 * math.pi * 10)) ** 1.5 * math.exp(-3 * r * r / 20), 0, 40, 4000)
    assert abs(nm - 1) < 1e-10                                             # Ex. 2.8 歸一化
    z, q = 0.4, 0.3
    assert (round(g_fourier(z, chi_1d(q)), 5), round(g_cont(z, q), 5), round(g_cont_book(z, q), 5)) == \
        (4.24213, 4.23729, 6.09756)
    for z in (0.3, 0.45):
        for x in (0, 1, 3):
            assert abs(G1_exact(z, x) - G1_fourier(z, x)) < 1e-9
    assert (round(G1_exact(0.45, 0), 4), round(G1_continuum(0.45, 0), 4)) == (2.2942, 2.357)

    # ---- §2.2 高斯積分
    for a in (1.0, 2.5):
        for n in range(0, 9):
            assert abs(gauss_moment(n, a) - gauss_moment_numeric(n, a)) < 1e-9
        for b in (0.5, 2.0):
            num = _simpson(lambda x: math.exp(-a * x * x + b * x), -15, 15, 6000)
            assert abs(num - gauss_gen(a, b)) < 1e-9
            for x in (-1.0, 0.0, 0.7):
                assert abs(square_complete(a, b, x) - (a * x * x - b * x)) < 1e-12
                assert abs(square_complete(a, b, x, book=True) - (a * x * x - b * x) + b * b / (4 * a)) < 1e-12  # 書上差 b²/4a
    assert [gauss_moment_frac(n) for n in (0, 2, 4, 6)] == [1, F(1, 2), F(3, 4), F(15, 8)]
    assert [round(gauss_moment(n), 4) for n in (0, 2, 4, 6)] == [1.7725, 0.8862, 1.3293, 3.3234]
    assert round(gauss_gen(1, 2), 4) == 4.818
    assert gauss_moment_frac(4) / gauss_moment_frac(2) ** 2 == 3                # 峰度 3（第 1 章 W105）

    # ---- §2.3 格子上的傅立葉
    for basis, om in ((SC_BASIS, 1.0), (BCC_BASIS, 4 / (3 * R3)), (FCC_BASIS, 1 / R2)):
        assert abs(abs(triple(*basis)) - om) < 1e-12
        assert all(abs(math.sqrt(dot(v, v)) - 1) < 1e-12 for v in basis)            # 都是最近鄰向量
        bs = reciprocal(*basis)
        for i in range(3):
            for j in range(3):
                assert abs(dot(basis[i], bs[j]) - (2 * math.pi if i == j else 0)) < 1e-12
        assert abs(abs(triple(*bs)) - (2 * math.pi) ** 3 / abs(triple(*basis))) < 1e-9
    assert round(CELL_VOLUME["bcc"], 4) == 0.7698 and round(CELL_VOLUME["fcc"], 4) == 0.7071
    assert abs(roots_sum(0, 8) - 8) < 1e-12 and all(abs(roots_sum(m, 8)) < 1e-12 for m in range(1, 8))
    assert all(abs(cmath.exp(1j * (math.pi / 4 + 2 * math.pi) * n) - cmath.exp(1j * math.pi / 4 * n)) < 1e-12
               for n in range(-5, 6))                                               # k 與 k + 2π 分不出來

    # ---- §2.4 冪級數係數的漸近
    assert coeff_two_poles(10, 1, 0, 0.5, 1) == 2048
    assert [round(coeff_two_poles(N, 1, 1, 0.5, 1.0)) for N in range(6)] == [3, 5, 9, 17, 33, 65]
    for N in (5, 20):
        assert abs(1 / (2 ** (N + 1) + 1) - 1 / (1 + 2 ** (N + 1))) < 1e-15
    for N, al, zc in ((10, 0.5, 0.25), (30, 2.5, 0.5), (7, 3.0, 1.0)):
        assert abs(coeff_power_exact(N, al, zc) / series_coeffs_power(al, zc, N)[N] - 1) < 1e-12
    assert all(abs(series_coeffs_power(0.5, 0.25, 12)[N] - 2 * math.comb(2 * N, N)) < 1e-6 * math.comb(2 * N, N)
               for N in range(13))                                                   # (1/4 − z)^{−1/2} = 2(1−4z)^{−1/2}
    assert math.comb(20, 10) == 184756 and round(4 ** 10 / math.sqrt(10 * math.pi), 2) == 187078.97
    assert round(math.comb(20, 10) / (4 ** 10 / math.sqrt(10 * math.pi)), 4) == 0.9876
    assert round(coeff_power_exact(1000, 0.5, 1.0) / coeff_power_asym(1000, 0.5, 1.0), 4) == 0.9999
    ex_log = coeff_log_exact(10, 1.0, 0.5)
    harm10 = sum(1 / k for k in range(1, 11))
    assert abs(ex_log - 2 ** 11 * (math.log(2) + harm10)) < 1e-9
    assert (round(ex_log, 1), round(coeff_log_asym(10, 1.0, 0.5), 1), round(coeff_log_asym(10, 1.0, 0.5, book=True), 1),
            round(coeff_log_asym(10, 1.0, 0.5, terms=1), 1)) == (7418.1, 7317.4, 4478.3, 4715.7)
    for N in (200, 1000):                                                           # 正確版的相對誤差 ~ 1/(2N ln…)
        assert abs(coeff_log_asym(N, 1.5, 0.5) / coeff_log_exact(N, 1.5, 0.5) - 1) < 2.0 / N
    assert abs(digamma(1.0) + 0.5772156649015329) < 1e-12
    for N, y in ((400, 0.5), (400, 1.5)):                                           # (2:S3-13) 被積分函數用 N^{−y−1}
        c = series_coeffs_power(-y, 0.5, N)[N]
        assert abs(c / (0.5 ** (-N + y) * N ** (-y - 1) / math.gamma(-y)) - 1) < 0.01
    assert [round(returns_sq(N) / asym_returns_sq(N), 4) for N in (10, 100)] == [0.9513, 0.995]
    assert round(math.comb(10, 5) / asym_returns_1d(10), 4) == 0.9754
    return True


# 由 check() 鎖定的模擬結果（固定種子；改了取樣方法就要重新核對 plan_ch02.md）
MORTAL_MEAN = 9.34                      # 理論 z/(1−z) = 9
MORTAL_AT_ORIGIN = 296                  # 296/2000 = 0.148；理論 0.1452
MORTAL_TRAIL_SEED = 25                  # 24 步，停在 (3, −3)
MORTAL_TRAIL_LEN = 24
SQ_RATIOS = [0.9513, 0.9753, 0.9901, 0.995]
UNIV_TABLE = {"sc": {2: 0.16667, 4: 0.06944, 6: 0.03987, 8: 0.02663, 10: 0.01939},
              "bcc": {2: 0.125, 4: 0.05273, 6: 0.03052, 8: 0.02044, 10: 0.0149},
              "fcc": {2: 0.08333, 4: 0.02604, 6: 0.01415, 8: 0.00939, 10: 0.00684}}
UNIV_SCALED10 = {"sc": 0.307, "bcc": 0.306, "fcc": 0.306}      # 極限 (3/2π)^{3/2} = 0.3299
CAUCHY_MEDIAN = 2.2644                  # 一步 Cauchy 的 |r| 中位數（r0 = 1）
MEDIANS = [("gauss", 1, 1.45), ("gauss", 10, 4.92), ("gauss", 100, 15.79),     # 理論 1.538√N
           ("cauchy", 1, 2.5), ("cauchy", 10, 23.28), ("cauchy", 100, 236.06)]  # 理論 2.264N


if __name__ == "__main__":
    import sys
    sys.stdout.reconfigure(encoding="utf-8")
    check()
    print("walk_gf（第 2 章）：Legendre、統計力學、格子計數、傅立葉、高斯極限、普適性、Cauchy、擴散、補充 1–3 全部驗證通過")
