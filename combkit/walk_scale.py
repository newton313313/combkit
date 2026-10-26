# 由 combkit/tools/sync_from_walk017.py 從 walk017 專案的 codes/utils/walk_scale.py 複製；請改 walk017 那份再同步。
"""
walk_scale.py —— 《Elements of the Random Walk》第 8 章 Properties of the random walk: introduction to scaling 的數學
（Python＋numpy／scipy（integrate、special、optimize），不 import manim；用到同資料夾的 walk_gf（Cauchy、冪級數係數、回原點計數）
與 walk_saw（二維自迴避的 c_N、Flory）——在 walk017 是 utils 套件、在 combkit 是 combkit 套件，都用相對 import；直接執行本檔時退回絕對 import）

約定（書上式號與筆誤見 ref/md/ch08_formulas.md）：
  長尾巴步伐（Pareto）：步長 r ≥ 1、P(步長 > r) = r^{−y}（徑向密度 y r^{−1−y}），方向在 d 維均勻；空間密度 p(r⃗) ∝ r^{−α}，α = d + y。
  書上 (8.3) 是 d = 3、α = 5 − δ，所以 y = 2 − δ；主例 y = 1（δ = 1、α = 4、χ(q) = 1 − (π/4)q + …）。y 越小尾巴越肥；y ≥ 2 回到高斯。
  一般 d 維格子：χ(q⃗) = a − b|q⃗|² + …（書上 §8.1.2 的 a 是 χ(0)，不是格距）；簡單立方（含一維、正方）a = 2d、b = 1、z_c = 1/(2d)；
  q 沿一個座標軸時 χ = 2cos q + 2(d − 1)。標度：τ = z_c − z；ν、η 是 (8.10) 的指數（普通漫步 ν = 1/2、η = 0）。
  最陡下降 (8.15)–(8.20)：A(x) = ln F(x) ≈ −Bx（x ≫ 1）；簡單立方 B = 2d（由 (8.8c) 的實空間形式）。
  重新標度 (8.21)–(8.25)：s′ = l₂s、r⃗′ = l₁r⃗；動能項乘 l₂/l₁²、自迴避項乘 l₁^d/l₂²。
  分數用 fractions.Fraction；回傳快取（lru_cache）陣列的函數不要就地修改結果。

§8.1.1 非高斯 (8.1)–(8.7)，Ex. 8.1
    pareto_alpha(y, d=3) = d + y；delta_of(y) = 2 − y（書上 d = 3 的 δ）
    moment2_trunc(y, R)             (8.2) 截到 R 的二階動差 ⟨r²⟩_R = y(R^{2−y} − 1)/(2−y)（y = 1：R − 1 = 9, 99, 999）；moment2_trunc_quad 數值積分
    chi_pareto(q, y, d=3)           χ(q) 的數值積分（d = 3：∫ y r^{−1−y} sin(qr)/(qr) dr；d = 1：cos(qr)）
    B_pareto(y, d=3)                (8.4) 的 B（閉式）：d = 3：y Γ(−1−y) cos(πy/2)（y = 1：π/4）；d = 1：Γ(1−y) cos(πy/2)（y = 1：π/2）
    B_pareto_quad(y, d=3)           同上的數值積分（第二種算法）
    chi_small(q, y, d=3) = 1 − B q^y；chi_ratio(q, y, d=3) = (1 − χ)/q^y（→ B）；chi_slope(y, d=3, q=1e-3)：ln(1−χ) 對 ln q 的斜率（→ y）
    sample_steps(rng, M, y, d=3)    M 個步伐 (M, d)；y = None 是長度 1 的普通步伐
    median_width(N, y, M=2000, seed=801)   M 人各走 N 步（三維）的 |R_N| 中位數（快取）
    width_pred(N, y)                理論值：y = 1：N·B·m₁（m₁ = 三維 Cauchy 的徑向中位數 2.2644，walk_gf）＝ 1.7785 N；
                                    y = None：√(N/3)·1.5382（Maxwell 中位數）＝ 0.8881 √N
    width_slope(y)                  N = 10, 100, 1000 的 log–log 斜率（y = 1 → 1、普通 → 1/2）
    trail(N, y, seed=801)           一條三維足跡 (N+1, 3)（畫圖取前兩欄）；TRAIL_SEED = 801
    levy_density(r, s, y)           三維：χ^N ≈ e^{−s q^y}（s = N B/A）的逆 Fourier（1/(2π²r))∫ q sin(qr) e^{−s q^y} dq（QAWF）
    levy_tail(r, s, y)              (8.5)／Ex. 8.1 由換變數得的尾巴 s Γ(2+y) sin(πy/2)/(2π²) r^{−3−y}；y = 1 時 = walk_gf.P_cauchy_N 的尾巴
    levy_slope(r, s, y)             ln levy_density 對 ln r 的局部斜率（→ −3 − y = −5 + δ）
    tail_exponent(delta, d=3)       Fraction：−(d + 2 − δ)（d = 3：−5 + δ）
    subst_integral(r, y, eps)       (8.6) ∫_0^∞ e^{ikr} k^y e^{−εk} dk（數值，複數）；subst_closed(r, y, eps) = Γ(y+1)/(ε − ir)^{y+1}
    subst_book(r, y, eps)           照書上 w = k/r（k = rw）：r^{y+1}∫ e^{i r² w} w^y e^{−εrw} dw（＝同一個積分，r 沒有被提出來）
§8.1.2 生成函數的標度形式 (8.8)–(8.14)
    zc(d) = Fraction(1, 2d)；AB(d) = (a, b) = (2d, 1)
    g_lattice(z, q, d=1)            1/(1 − zχ(q))（q 沿一軸）；g_88b(z, q, d)、g_88c(z, q, d)：(8.8b)(8.8c)
    collapse(eps, x, d=1)           z = z_c − eps、q = x√eps 時的 (z_c − z)·g（塌縮的縱軸）；collapse_limit(x, d) = (1/a)/(1 + b x²/a²)
    COLLAPSE_EPS = (0.1, 0.05, 0.01)、COLLAPSE_X = (0, 1, 2, 4)；一維 g(z;0) = 5, 10, 50；x = 2：0.2819、0.2652、0.2529 → 0.25
    g_general(tau, q, nu, eta)      (8.10) 的一個具體 f：f(x) = x^{2−η}/(1 + x²)^{1−η/2}（η = 0、ν = 1/2 就是 1/(q² + τ)）
    gamma_numeric(nu, eta)          ln g_general(τ, q → 0) 對 ln τ 的斜率（→ −γ，Ex. 8.3）
    G_lattice_sc(z, n)              三維簡單立方的精確格子 G(z; (n,0,0)) = ∫_0^∞ e^{−t} I_n(2zt) I_0(2zt)² dt
    xi_88c(x) = 1/√(6(1 − x))（x = z/z_c）；xi_lattice(x) = 1/arccosh((1 − 4z)/(2z))（格子上的精確衰減長度）
    G_yukawa(x, r)                  (8.8c) 回到實空間：(3/2π) e^{−r/ξ}/r（d = 3、a = 6、b = 1）；F3(u) = (3/2π) e^{−6u}（(8.13)，u = r√(z_c − z)）
    XI_CASES = (0.96, 0.99)         ξ = 2.0412、4.0825（格子 2.0205、4.0722）；r = 8 的格子/連續比 G_ratio(x, 8)
§8.1.2 最陡下降 (8.15)–(8.20)，Ex. 8.2–8.5
    ratio_exact(N, x) = C(N, (N+x)/2)/C(N, N/2)（一維精確，float）
    sd_K(nu, B, zc, book=False)     (8.20) 指數的係數：正確 (1−ν)((νz_c)^ν B)^{1/(1−ν)}；book=True 照書上（少了 1−ν）
    sd_ratio(N, r, nu, B, zc, book=False) = exp[−K (N+1)^{−ν/(1−ν)} r^{1/(1−ν)}]（一維 ν = 1/2、B = 2、z_c = 1/2：e^{−x²/2(N+1)}；書上 e^{−x²/(N+1)}）
    sd_delta(N, r, nu, B, zc)       (8.19)；sd_delta_num：直接解 dφ/dδ = 0，φ(δ) = −B r δ^ν − (N+1) ln(z_c − δ)（第二種算法）
    sd_phi_num(N, r, nu, B, zc)     φ(δ*) + (N+1) ln z_c（→ −K(N+1)^{−ν/(1−ν)} r^{1/(1−ν)}，N 大時）
    sd_K_dim(d) = d/2               簡單立方 ν = 1/2、B = 2d、z_c = 1/2d 時的 K：指數 −d r²/(2N)＝(2.20)
    SD_N = (20, 100)、SD_X = {20: (4, 8), 100: (8, 12)}  卡上的例子
    fisher(nu) = 1/(1 − ν)（Fraction）；FISHER_NU = (1/2, 3/5, 3/4)：2、5/2、4；fisher_num(nu)：由 sd_ratio 量 log 斜率
    gamma_exp(eta, nu) = (2 − η)ν（Fraction）；ETA_2D = 5/24、NU_2D = 3/4 → 43/32（已知值，只引用）
    growth(gamma, zc) = (x, p) = (1/z_c, γ − 1)（Ex. 8.4）；p_numeric(gamma, zc, N=1000)：(z_c − z)^{−γ} 係數的局部斜率（log_coeff_zN(N, gamma, zc) = ln[z_c^N·係數]）
    ratio_method(N)                 二維自迴避 c_N 的比值法：(r_N, μ̂_N, γ̂_N)；r_N = c_N/c_{N−1}、μ̂_N = (N r_N − (N−2) r_{N−2})/2、
                                    γ̂_N = 1 + (N/2)(c_N/(c_{N−2} μ²) − 1)（μ = MU_QUOTED）；c_15、c_16 用 walk_saw.SAW_KNOWN_C（文獻）
    RATIO_N = (10, 12, 14, 16)      μ̂ = 2.6592, 2.6529, 2.6483, 2.6454 → 2.638；γ̂ = 1.3552, 1.3519, 1.3496, 1.3483 → 43/32 = 1.34375
    MU_QUOTED = 2.638158、GAMMA_QUOTED = 43/32（只引用）
    theta(gamma, nu) = (γ − 1)/ν（Ex. 8.5：小 r 時條數 ∝ r^θ；二維自迴避 11/24）；ex85_s(d, nu) = d − 1/ν（F(x) ≈ F₀ + F₁x^s）
    return_slope(d, N)              普通漫步回原點條數 C(N;0)/(2d)^N 的局部 log 斜率（→ −dν = −d/2；Ex. 8.5 的 N^{−dν}）
§8.1.3 重新標度與 Flory (8.21)–(8.25)，Ex. 8.6–8.8
    kin_factor(l1, l2) = l2/l1²；int_factor(l1, l2, d) = l1^d/l2²；int_factor_book86(l1, l2, d) = l1^d/l2（Ex. 8.6 書上）
    rescale_check(l1, l2, d)        離散高斯鏈＋寬 w 的高斯「δ」：數值量出 (K/K′, V/V′)（＝ 上面兩個因子）
    sa_amplitude(N, d)              l₁ = √l₂、l₂ = 1/N 時自迴避項的放大倍數 N^{(4−d)/2}（Fraction 或 float）；N = 100、d = 1..5：1000, 100, 10, 1, 0.1
    UPPER_CRITICAL_D = 4
    flory_l1(l2, d) = l2^{3/(d+2)}；flory_l1_num：數值解 l2/l1² = l1^d/l2²；nu_rescale(d) = 3/(d+2)（＝ walk_saw.flory_nu）
    ex88_factor(l1, l2, q) = l1^q/l2²；ex88_amplitude(N, q) = N^{(4−q)/2}；ex88_relevant(d, q) = min(d, q) < 4；ex88_nu(q) = 3/(q+2)
summary()                         場景卡與講義用的全部數字（dict，字串已照卡上的小數位格式化）

check() 驗證全部（每個結果至少兩種獨立算法、筆誤的更正、plan_ch08.md 場景卡與講義的每個數字）；直接執行本檔會跑 check()（約 1–2 秒）。
`python walk_scale.py --summary` 印出 summary() 的全部數字。
最慢的呼叫：median_width(1000, …) 約 0.1 秒（快取）、G_lattice_sc 每次約 0.01 秒，其餘都很快。scipy 的 IntegrationWarning 不影響結果。
"""
from __future__ import annotations

import math
from fractions import Fraction
from functools import lru_cache

import numpy as np
from scipy import integrate, optimize, special

try:                                    # walk017：utils 套件；combkit：combkit 套件
    from . import walk_gf as GF
    from . import walk_saw as SW
except ImportError:                     # 直接執行本檔
    import walk_gf as GF
    import walk_saw as SW


# =========================================================================== §8.1.1 非高斯


def pareto_alpha(y, d=3):
    return d + y


def delta_of(y):
    return 2 - y


def moment2_trunc(y, R):
    """(8.2)：⟨r²⟩_R = ∫_1^R y r^{−1−y} r² dr。"""
    if y == 2:
        return 2 * math.log(R)
    return y * (R ** (2 - y) - 1) / (2 - y)


def moment2_trunc_quad(y, R):
    return integrate.quad(lambda r: y * r ** (1 - y), 1, R, limit=200)[0]


def chi_pareto(q, y, d=3):
    """χ(q) = ⟨e^{iq⃗·r⃗}⟩：d = 3 平均 sin(qr)/(qr)、d = 1 平均 cos(qr)（QAWF 振盪積分）。"""
    if d == 3:
        return integrate.quad(lambda r: y * r ** (-2 - y) / q, 1, np.inf, weight="sin", wvar=q, limit=400)[0]
    if d == 1:
        return integrate.quad(lambda r: y * r ** (-1 - y), 1, np.inf, weight="cos", wvar=q, limit=400)[0]
    raise ValueError("d = 1 或 3")


def B_pareto(y, d=3):
    if d == 3:
        return math.pi / 4 if y == 1 else y * special.gamma(-1 - y) * math.cos(math.pi * y / 2)
    if d == 1:
        return math.pi / 2 if y == 1 else special.gamma(1 - y) * math.cos(math.pi * y / 2)
    raise ValueError("d = 1 或 3")


def B_pareto_quad(y, d=3):
    """B = y ∫_0^∞ u^{−1−y}(1 − sin u/u) du（d = 3）或 y ∫_0^∞ u^{−1−y}(1 − cos u) du（d = 1）。"""
    if d == 3:
        k = lambda u: (1 - math.sin(u) / u) if u > 1e-3 else u * u / 6 - u ** 4 / 120
    else:
        k = lambda u: 1 - math.cos(u) if u > 1e-3 else u * u / 2 - u ** 4 / 24
    f = lambda u: y * u ** (-1 - y) * k(u)
    a = integrate.quad(f, 0, 50, limit=400)[0]
    # 尾端：u > 50 時 k = 1 − 振盪項；∫_50^∞ y u^{−1−y} du = 50^{−y}，振盪項用 QAWF
    if d == 3:
        osc = integrate.quad(lambda u: y * u ** (-2 - y), 50, np.inf, weight="sin", wvar=1.0, limit=400)[0]
    else:
        osc = integrate.quad(lambda u: y * u ** (-1 - y), 50, np.inf, weight="cos", wvar=1.0, limit=400)[0]
    return a + 50.0 ** (-y) - osc


def chi_small(q, y, d=3):
    return 1 - B_pareto(y, d) * q ** y


def chi_ratio(q, y, d=3):
    return (1 - chi_pareto(q, y, d)) / q ** y


def chi_slope(y, d=3, q=1e-3):
    return math.log((1 - chi_pareto(2 * q, y, d)) / (1 - chi_pareto(q, y, d))) / math.log(2)


def sample_steps(rng, M, y, d=3):
    v = rng.normal(size=(M, d))
    v /= np.linalg.norm(v, axis=1)[:, None]
    if y is None:
        return v
    return v * (rng.random(M) ** (-1.0 / y))[:, None]


WIDTH_N = (10, 100, 1000)
WIDTH_M = 2000
WIDTH_SEED = 801
TRAIL_SEED = 801


@lru_cache(maxsize=None)
def median_width(N, y, M=WIDTH_M, seed=WIDTH_SEED):
    rng = np.random.default_rng(seed)
    pos = np.zeros((M, 3))
    for _ in range(N):
        pos += sample_steps(rng, M, y)
    return float(np.median(np.linalg.norm(pos, axis=1)))


def width_pred(N, y):
    if y is None:
        return math.sqrt(N / 3) * GF.CHI3_MEDIAN
    if y == 1:
        return N * B_pareto(1) * GF.cauchy_radial_median(1.0)
    raise ValueError("只做 y = 1 與普通步伐")


def width_slope(y):
    xs = np.log(WIDTH_N)
    ys = np.log([median_width(N, y) for N in WIDTH_N])
    return float(np.polyfit(xs, ys, 1)[0])


@lru_cache(maxsize=None)
def trail(N, y, seed=TRAIL_SEED):
    rng = np.random.default_rng(seed)
    return np.vstack([np.zeros(3), np.cumsum(sample_steps(rng, N, y), axis=0)])


def levy_density(r, s, y):
    v = integrate.quad(lambda q: q * math.exp(-s * q ** y), 0, np.inf, weight="sin", wvar=r, limit=400)[0]
    return v / (2 * math.pi ** 2 * r)


def levy_tail(r, s, y):
    return s * special.gamma(2 + y) * math.sin(math.pi * y / 2) / (2 * math.pi ** 2) * r ** (-3 - y)


def levy_slope(r, s, y, h=1.05):
    return math.log(levy_density(r * h, s, y) / levy_density(r / h, s, y)) / math.log(h * h)


def tail_exponent(delta, d=3):
    return -(d + 2 - Fraction(delta))


def subst_integral(r, y, eps):
    re = integrate.quad(lambda k: k ** y * math.exp(-eps * k), 0, np.inf, weight="cos", wvar=r, limit=400)[0]
    im = integrate.quad(lambda k: k ** y * math.exp(-eps * k), 0, np.inf, weight="sin", wvar=r, limit=400)[0]
    return complex(re, im)


def subst_closed(r, y, eps):
    return special.gamma(y + 1) / complex(eps, -r) ** (y + 1)


def subst_book(r, y, eps):
    """照書上 w = k/r（k = rw、dk = r dw）：r^{y+1}∫ e^{i r² w} w^y e^{−ε r w} dw——r 仍留在指數裡。"""
    re = integrate.quad(lambda w: w ** y * math.exp(-eps * r * w), 0, np.inf, weight="cos", wvar=r * r, limit=400)[0]
    im = integrate.quad(lambda w: w ** y * math.exp(-eps * r * w), 0, np.inf, weight="sin", wvar=r * r, limit=400)[0]
    return r ** (y + 1) * complex(re, im)


SUBST_Y = 0.5
SUBST_EPS = 1e-3


# =========================================================================== §8.1.2 標度形式


def zc(d):
    return Fraction(1, 2 * d)


def AB(d):
    return 2 * d, 1


def chi_axis(q, d):
    return 2 * math.cos(q) + 2 * (d - 1)


def g_lattice(z, q, d=1):
    return 1 / (1 - z * chi_axis(q, d))


def g_88b(z, q, d=1):
    a, b = AB(d)
    return (1 / a) / (float(zc(d)) - z + z * b / a * q * q)


def g_88c(z, q, d=1):
    a, b = AB(d)
    return (1 / a) / (float(zc(d)) - z + b / a ** 2 * q * q)


COLLAPSE_EPS = (0.1, 0.05, 0.01)
COLLAPSE_X = (0, 1, 2, 4)


def collapse(eps, x, d=1):
    return eps * g_lattice(float(zc(d)) - eps, x * math.sqrt(eps), d)


def collapse_limit(x, d=1):
    a, b = AB(d)
    return (1 / a) / (1 + b * x * x / a ** 2)


def g_general(tau, q, nu, eta):
    x = q * tau ** (-nu)
    return q ** (-2 + eta) * x ** (2 - eta) / (1 + x * x) ** (1 - eta / 2)


def gamma_numeric(nu, eta, q=1e-12):
    t1, t2 = 1e-3, 1e-4
    return -math.log(g_general(t2, q, nu, eta) / g_general(t1, q, nu, eta)) / math.log(t2 / t1)


@lru_cache(maxsize=None)
def G_lattice_sc(z, n):
    """∫_0^∞ e^{−t} I_n(2zt) I_0(2zt)² dt（1/(1 − zχ) = ∫ e^{−t(1−zχ)} dt，每個方向 ∫ e^{2zt cos q} e^{iqn} dq/2π = I_n(2zt)）。"""
    f = lambda t: special.ive(n, 2 * z * t) * special.ive(0, 2 * z * t) ** 2 * math.exp(-t * (1 - 6 * z))
    return integrate.quad(f, 0, 50, limit=200)[0] + integrate.quad(f, 50, np.inf, limit=400)[0]


XI_CASES = (0.96, 0.99)


def xi_88c(x):
    return 1 / math.sqrt(6 * (1 - x))


def xi_lattice(x):
    z = x / 6
    return 1 / math.acosh((1 - 4 * z) / (2 * z))


def F3(u):
    return 3 / (2 * math.pi) * math.exp(-6 * u)


def G_yukawa(x, r):
    return 3 / (2 * math.pi) * math.exp(-r / xi_88c(x)) / r


def G_ratio(x, n):
    return G_lattice_sc(x / 6, n) / G_yukawa(x, n)


def G_yukawa_ft(x, r):
    """(8.8c) 的數值逆 Fourier：(1/(2π²r))∫ q sin(qr) (a/b)/(q² + m²) dq（第二種算法；m = 1/ξ）。"""
    m2 = 1 / xi_88c(x) ** 2
    return 6 * integrate.quad(lambda q: q / (q * q + m2), 0, np.inf, weight="sin", wvar=r, limit=400)[0] / (2 * math.pi ** 2 * r)


# =========================================================================== §8.1.2 最陡下降


def ratio_exact(N, x):
    return math.comb(N, (N + x) // 2) / math.comb(N, N // 2)


def sd_K(nu, B, zc_, book=False):
    nu, zc_ = float(nu), float(zc_)
    K = ((nu * zc_) ** nu * B) ** (1 / (1 - nu))
    return K if book else (1 - nu) * K


def sd_ratio(N, r, nu, B, zc_, book=False):
    nu = float(nu)
    return math.exp(-sd_K(nu, B, zc_, book) * (N + 1) ** (-nu / (1 - nu)) * r ** (1 / (1 - nu)))


def sd_delta(N, r, nu, B, zc_):
    nu, zc_ = float(nu), float(zc_)
    return ((N + 1) / (nu * B * zc_ * r)) ** (1 / (nu - 1))


def sd_delta_num(N, r, nu, B, zc_):
    nu, zc_ = float(nu), float(zc_)
    dphi = lambda dl: -nu * B * r * dl ** (nu - 1) + (N + 1) / (zc_ - dl)
    return optimize.brentq(dphi, 1e-14, zc_ * (1 - 1e-12), xtol=1e-16)


def sd_phi_num(N, r, nu, B, zc_):
    nu, zc_ = float(nu), float(zc_)
    dl = sd_delta_num(N, r, nu, B, zc_)
    return -B * r * dl ** nu - (N + 1) * math.log(zc_ - dl) + (N + 1) * math.log(zc_)


def sd_K_dim(d):
    return sd_K(Fraction(1, 2), 2 * d, zc(d))


SD_N = (20, 100)
SD_X = {20: (4, 8), 100: (8, 12)}


def fisher(nu):
    return 1 / (1 - Fraction(nu))


FISHER_NU = (Fraction(1, 2), Fraction(3, 5), Fraction(3, 4))


def fisher_num(nu, N=10 ** 6, r=1000.0):
    return math.log(math.log(sd_ratio(N, 2 * r, nu, 2, 0.5)) / math.log(sd_ratio(N, r, nu, 2, 0.5))) / math.log(2)


def gamma_exp(eta, nu):
    return (2 - Fraction(eta)) * Fraction(nu)


ETA_2D = Fraction(5, 24)
NU_2D = Fraction(3, 4)
MU_QUOTED = 2.638158
GAMMA_QUOTED = Fraction(43, 32)


def growth(gamma, zc_):
    return 1 / zc_, gamma - 1


def log_coeff_zN(N, gamma, zc_):
    """ln[z_c^N × (z^N 在 (z_c − z)^{−γ} 的係數)] = ln Γ(N+γ) − ln Γ(γ) − ln N! − γ ln z_c。"""
    g = float(gamma)
    return math.lgamma(N + g) - math.lgamma(g) - math.lgamma(N + 1) - g * math.log(zc_)


def p_numeric(gamma, zc_, N=1000):
    return (log_coeff_zN(2 * N, gamma, zc_) - log_coeff_zN(N, gamma, zc_)) / math.log(2)


def _cN(N):
    return SW.saw_c(N) if N <= 14 else SW.SAW_KNOWN_C[N]


RATIO_N = (10, 12, 14, 16)


def ratio_method(N, mu=MU_QUOTED):
    rN = _cN(N) / _cN(N - 1)
    rN2 = _cN(N - 2) / _cN(N - 3)
    mu_hat = (N * rN - (N - 2) * rN2) / 2
    gam_hat = 1 + N / 2 * (_cN(N) / (_cN(N - 2) * mu * mu) - 1)
    return rN, mu_hat, gam_hat


def theta(gamma, nu):
    return (Fraction(gamma) - 1) / Fraction(nu)


def ex85_s(d, nu):
    return d - 1 / Fraction(nu)


def _returns(d, N):
    if d == 1:
        return math.comb(N, N // 2)
    if d == 2:
        return GF.returns_sq(N)
    return GF.returns_sc(N)


def return_slope(d, N):
    a = _returns(d, N) / (2 * d) ** N
    b = _returns(d, 2 * N) / (2 * d) ** (2 * N)
    return math.log(b / a) / math.log(2)


# =========================================================================== §8.1.3 重新標度


def kin_factor(l1, l2):
    return l2 / l1 ** 2


def int_factor(l1, l2, d):
    return l1 ** d / l2 ** 2


def int_factor_book86(l1, l2, d):
    return l1 ** d / l2


def rescale_check(l1, l2, d, n=200, w=0.5, seed=804):
    """離散高斯鏈 r_i（s_i = iΔs）；K = Σ|Δr|²/Δs、V = Σ_{i≠j} Δs² g_w(r_i − r_j)（g_w 是寬 w 的歸一高斯＝δ 的近似）。
    新變數 r′ = l1 r、s′ = l2 s、w′ = l1 w 再算一次 → (K/K′, V/V′)。"""
    rng = np.random.default_rng(seed)
    ds = 1.0
    r = np.cumsum(rng.normal(size=(n, d)), axis=0)

    def KV(r, ds, w):
        K = float(np.sum(np.diff(r, axis=0) ** 2) / ds)
        D2 = np.sum((r[:, None, :] - r[None, :, :]) ** 2, axis=2)
        g = np.exp(-D2 / (2 * w * w)) / (2 * math.pi * w * w) ** (d / 2)
        np.fill_diagonal(g, 0.0)
        return K, float(ds * ds * g.sum())

    K, V = KV(r, ds, w)
    K2, V2 = KV(l1 * r, l2 * ds, l1 * w)
    return K / K2, V / V2


def sa_amplitude(N, d):
    s = math.isqrt(N)
    if s * s == N:                                       # N 是平方數：精確的 Fraction
        return Fraction(s) ** (4 - d)
    return N ** ((4 - d) / 2)


SA_N = 100
UPPER_CRITICAL_D = 4


def flory_l1(l2, d):
    return l2 ** (3 / (d + 2))


def flory_l1_num(l2, d):
    return optimize.brentq(lambda l1: math.log(l2 / l1 ** 2) - math.log(l1 ** d / l2 ** 2), 1e-12, 1.0)


def nu_rescale(d):
    return Fraction(3, d + 2)


def ex88_factor(l1, l2, q):
    return l1 ** q / l2 ** 2


def ex88_amplitude(N, q):
    return N ** ((4 - q) / 2)


def ex88_relevant(d, q):
    return min(d, q) < 4


def ex88_nu(q):
    return Fraction(3, q + 2)


# =========================================================================== 卡上的數字


def summary():
    S = {}
    S["moment2_y1"] = [moment2_trunc(1, R) for R in (10, 100, 1000)]                 # 9, 99, 999
    S["B3_y1"] = f"{B_pareto(1):.4f}"                                                 # π/4
    S["chi_ratio_y1"] = [f"{chi_ratio(q, 1):.4f}" for q in (0.1, 0.01, 0.001)]
    S["chi_y1_q01"] = f"{chi_pareto(0.01, 1):.5f}", f"{chi_small(0.01, 1):.5f}"
    S["chi_slope"] = {y: f"{chi_slope(y):.2f}" for y in (1, 1.5)}
    S["width_levy"] = [f"{median_width(N, 1):.0f}" for N in WIDTH_N]
    S["width_levy_perN"] = [f"{median_width(N, 1) / N:.2f}" for N in WIDTH_N]
    S["width_unit"] = [f"{median_width(N, None):.2f}" for N in WIDTH_N]
    S["width_unit_persqrt"] = [f"{median_width(N, None) / math.sqrt(N):.2f}" for N in WIDTH_N]
    S["width_pred"] = f"{width_pred(1, 1):.4f}", f"{width_pred(1, None):.4f}"
    S["width_slope"] = f"{width_slope(1):.2f}", f"{width_slope(None):.2f}"
    S["tail_ratio_y1"] = [f"{levy_density(r, 1, 1) / levy_tail(r, 1, 1):.4f}" for r in (10, 30, 100)]
    S["tail_ratio_y15"] = [f"{levy_density(r, 1, 1.5) / levy_tail(r, 1, 1.5):.4f}" for r in (10, 30, 100)]
    S["tail_slope"] = f"{levy_slope(100, 1, 1):.2f}", f"{levy_slope(100, 1, 1.5):.2f}"
    S["subst_gamma"] = f"{special.gamma(SUBST_Y + 1):.4f}"
    S["subst_ratio"] = f"{abs(subst_integral(2, SUBST_Y, SUBST_EPS)) / abs(subst_integral(1, SUBST_Y, SUBST_EPS)):.4f}"
    S["g0_1d"] = [f"{g_lattice(0.5 - e, 0.0):.0f}" for e in COLLAPSE_EPS]
    S["collapse"] = {e: [f"{collapse(e, x):.4f}" for x in COLLAPSE_X] for e in COLLAPSE_EPS}
    S["collapse_limit"] = [f"{collapse_limit(x):.4f}" for x in COLLAPSE_X]
    S["xi"] = {x: (f"{xi_88c(x):.4f}", f"{xi_lattice(x):.4f}") for x in XI_CASES}
    S["G_ratio8"] = {x: f"{G_ratio(x, 8):.4f}" for x in XI_CASES}
    S["sd"] = {N: [(x, f"{ratio_exact(N, x):.4f}", f"{sd_ratio(N, x, 0.5, 2, 0.5):.4f}",
                    f"{sd_ratio(N, x, 0.5, 2, 0.5, book=True):.4f}") for x in SD_X[N]] for N in SD_N}
    S["sd_K"] = sd_K(0.5, 2, 0.5), sd_K(0.5, 2, 0.5, book=True)                       # 1/2、1
    S["sd_K_dim3"] = sd_K_dim(3)                                                         # 3/2
    S["fisher"] = [fisher(nu) for nu in FISHER_NU] + [f"{1 / (1 - SW.NU_QUOTED[3]):.2f}"]
    S["gamma_2d"] = gamma_exp(ETA_2D, NU_2D)
    S["growth_2d"] = f"{MU_QUOTED:.3f}", GAMMA_QUOTED - 1
    S["ratio"] = {N: (f"{ratio_method(N)[1]:.4f}", f"{ratio_method(N)[2]:.4f}") for N in RATIO_N}
    S["theta_2d"] = theta(GAMMA_QUOTED, NU_2D)
    S["sa_amp"] = [sa_amplitude(SA_N, d) for d in range(1, 6)]
    S["flory_d3"] = f"{flory_l1(0.01, 3):.4f}", f"{1 / flory_l1(0.01, 3):.2f}"
    S["nu_rescale"] = [nu_rescale(d) for d in range(1, 5)]
    S["ex88_nu"] = [ex88_nu(q) for q in (1, 2, 3, 4)]
    return S


# =========================================================================== 自我測試


def _fmt(x, n):
    return f"{float(x):.{n}f}"


def check():
    import warnings
    warnings.simplefilter("ignore")
    F = Fraction

    # ---------------- §8.1.1：(8.2) 二階動差發散
    assert [moment2_trunc(1, R) for R in (10, 100, 1000)] == [9, 99, 999]
    for y in (0.5, 1, 1.5):
        for R in (10, 1000):
            assert abs(moment2_trunc(y, R) - moment2_trunc_quad(y, R)) < 1e-8 * moment2_trunc(y, R)
    assert pareto_alpha(1) == 4 and delta_of(1) == 1 and pareto_alpha(F(3, 2)) == F(9, 2)
    # (8.4)：B 的閉式＝數值積分；χ 的數值→ 1 − B q^y
    for y in (0.5, 1, 1.5):
        for d in (1, 3):
            assert abs(B_pareto(y, d) - B_pareto_quad(y, d)) < 1e-6, (y, d)
    assert _fmt(B_pareto(1), 4) == "0.7854" and abs(B_pareto(1, 1) - math.pi / 2) < 1e-15
    assert abs(B_pareto(1.0000001) - math.pi / 4) < 1e-6 and abs(B_pareto(0.9999999, 1) - math.pi / 2) < 1e-6
    assert [_fmt(chi_ratio(q, 1), 4) for q in (0.1, 0.01, 0.001)] == ["0.7687", "0.7837", "0.7852"]
    assert _fmt(chi_pareto(0.01, 1), 5) == "0.99216" and _fmt(chi_small(0.01, 1), 5) == "0.99215"
    assert abs(chi_ratio(1e-3, 0.5) - B_pareto(0.5)) < 1e-5
    assert abs(chi_ratio(1e-3, 1, d=1) - math.pi / 2) < 2e-3
    assert _fmt(chi_slope(1), 2) == "1.00" and _fmt(chi_slope(1.5), 2) == "1.49"
    # 一步的 χ(q) 也可由「蒙地卡羅平均 sin(qr)/(qr)」核對（第二種算法）
    rng = np.random.default_rng(1)
    rr = rng.random(400000) ** -1.0
    assert abs(np.mean(np.sin(0.1 * rr) / (0.1 * rr)) - chi_pareto(0.1, 1)) < 2e-3
    # ---------------- 模擬：寬度 ∝ N^{1/y}（固定種子）
    lev = [median_width(N, 1) for N in WIDTH_N]
    uni = [median_width(N, None) for N in WIDTH_N]
    assert [f"{v:.0f}" for v in lev] == ["18", "172", "1822"]
    assert [_fmt(v / N, 2) for v, N in zip(lev, WIDTH_N)] == ["1.77", "1.72", "1.82"]
    assert [_fmt(v, 2) for v in uni] == ["2.86", "8.98", "28.06"]
    assert [_fmt(v / math.sqrt(N), 2) for v, N in zip(uni, WIDTH_N)] == ["0.90", "0.90", "0.89"]
    assert _fmt(width_pred(1, 1), 4) == "1.7785" and _fmt(width_pred(1, None), 4) == "0.8881"
    for v, N in zip(lev, WIDTH_N):
        assert abs(v / width_pred(N, 1) - 1) < 0.04
    for v, N in zip(uni, WIDTH_N):
        assert abs(v / width_pred(N, None) - 1) < 0.03
    assert _fmt(width_slope(1), 2) == "1.01" and _fmt(width_slope(None), 2) == "0.50"
    # Maxwell 中位數（自由度 3 的 χ 分布）與三維 Cauchy 徑向中位數：兩種算法
    assert abs(special.gammainc(1.5, GF.CHI3_MEDIAN ** 2 / 2) - 0.5) < 1e-12
    m1 = GF.cauchy_radial_median(1.0)
    assert abs(GF.cauchy_radial_cdf(m1) - 0.5) < 1e-9 and abs(m1 - 2.2644) < 1e-4
    t = trail(1000, 1)
    assert t.shape == (1001, 3) and np.all(t[0] == 0)
    # ---------------- (8.5)、Ex. 8.1：尾巴 r^{−5+δ}：逆 Fourier（數值）vs 換變數的漸近式 vs 第 2 章的 Cauchy 閉式
    for r in (3, 10, 30, 100):
        assert abs(levy_density(r, 1, 1) / GF.P_cauchy_N(r, 1, 1.0) - 1) < 1e-6
        assert abs(levy_tail(r, 1, 1) - 1 / (math.pi ** 2 * r ** 4)) < 1e-15
    assert [_fmt(levy_density(r, 1, 1) / levy_tail(r, 1, 1), 4) for r in (10, 30, 100)] == ["0.9803", "0.9978", "0.9998"]
    assert [_fmt(levy_density(r, 1, 1.5) / levy_tail(r, 1, 1.5), 4) for r in (10, 30, 100)] == ["1.1755", "1.0316", "1.0051"]
    assert _fmt(levy_slope(100, 1, 1), 2) == "-4.00" and _fmt(levy_slope(100, 1, 1.5), 2) == "-4.51"
    assert tail_exponent(1) == -4 and tail_exponent(F(1, 2)) == F(-9, 2)
    assert tail_exponent(delta_of(F(3, 2))) == -3 - F(3, 2)
    # N 步：s = N B（χ^N ≈ e^{−N B q^y}），尾巴 ∝ N（固定 r、N 變 2 倍 → 尾巴 2 倍）
    assert abs(levy_density(200, 2, 1.5) / levy_density(200, 1, 1.5) - 2) < 0.02
    # (8.6)(8.7) 換變數 w = kr：|I(2r)|/|I(r)| = 2^{−(y+1)}；數值＝閉式；書上的 w = k/r 是同一個積分但 r 提不出來
    for r in (1, 2, 4):
        assert abs(subst_integral(r, SUBST_Y, 0.05) - subst_closed(r, SUBST_Y, 0.05)) < 1e-6
        assert abs(subst_book(r, SUBST_Y, 0.05) - subst_closed(r, SUBST_Y, 0.05)) < 1e-6
        assert abs(abs(subst_closed(r, SUBST_Y, 1e-9)) * r ** 1.5 - special.gamma(1.5)) < 1e-9
    assert _fmt(special.gamma(1.5), 4) == "0.8862"
    assert _fmt(abs(subst_integral(2, SUBST_Y, SUBST_EPS)) / abs(subst_integral(1, SUBST_Y, SUBST_EPS)), 4) == "0.3536"
    assert _fmt(2 ** -1.5, 4) == "0.3536"

    # ---------------- §8.1.2：(8.8) 與標度塌縮
    assert zc(1) == F(1, 2) and zc(3) == F(1, 6) and AB(3) == (6, 1)
    for d in (1, 2, 3):
        a, b = AB(d)
        qs = np.linspace(-1e-3, 1e-3, 5)
        assert np.allclose([chi_axis(q, d) for q in qs], [a - b * q * q for q in qs], atol=1e-12)
        z = float(zc(d)) * 0.999
        assert abs(g_88b(z, 1e-3, d) / g_lattice(z, 1e-3, d) - 1) < 1e-4                # (8.8a)→(8.8b) 是恆等式（到 q²）
        assert abs(g_88c(z, 1e-2, d) / g_lattice(z, 1e-2, d) - 1) < 2e-3
    assert [f"{g_lattice(0.5 - e, 0.0):.0f}" for e in COLLAPSE_EPS] == ["5", "10", "50"]
    want = {0.1: ["0.5000", "0.4172", "0.2819", "0.1317"], 0.05: ["0.5000", "0.4085", "0.2652", "0.1145"],
            0.01: ["0.5000", "0.4017", "0.2529", "0.1027"]}
    for e in COLLAPSE_EPS:
        assert [_fmt(collapse(e, x), 4) for x in COLLAPSE_X] == want[e]
    assert [_fmt(collapse_limit(x), 4) for x in COLLAPSE_X] == ["0.5000", "0.4000", "0.2500", "0.1000"]
    assert abs(collapse(1e-6, 2) - 0.25) < 1e-5
    # (8.9)：|q|^{−2} f(q τ^{−1/2})，f(x) = (1/a)/(b/a² + x^{−2})
    for (tau, q) in ((0.01, 0.3), (0.001, 0.05)):
        a, b = AB(1)
        x = q / math.sqrt(tau)
        assert abs(q ** -2 * (1 / a) / (b / a ** 2 + x ** -2) - g_88c(0.5 - tau, q)) < 1e-9
    # Ex. 8.3：γ = (2 − η)ν：數值斜率＝公式
    for nu, eta in ((0.5, 0.0), (0.75, 5 / 24), (0.588, 0.031)):
        assert abs(gamma_numeric(nu, eta) - (2 - eta) * nu) < 1e-6
    assert abs(g_general(0.01, 0.3, 0.5, 0.0) - 1 / (0.09 + 0.01)) < 1e-12
    # (8.11)–(8.14)：三維實空間 G ≈ (3/2π) e^{−r/ξ}/r，ξ ∝ (z_c − z)^{−1/2}
    assert [_fmt(xi_88c(x), 4) for x in XI_CASES] == ["2.0412", "4.0825"]
    assert [_fmt(xi_lattice(x), 4) for x in XI_CASES] == ["2.0205", "4.0722"]
    assert abs(xi_88c(0.99) / xi_88c(0.96) - 2) < 1e-12
    for x in XI_CASES:
        for r in (2.0, 8.0):
            assert abs(G_yukawa_ft(x, r) / G_yukawa(x, r) - 1) < 1e-6
    # 格子的精確 G 用另一種算法核對：回原點的級數 Σ z^N C(N; 0)（N ≤ 60 精確、之後用漸近式；x = 0.96）
    zz = 0.96 / 6
    ser = sum(GF.returns_sc(N) * zz ** N for N in range(0, 61, 2))
    ser += sum(math.exp(N * math.log(zz) + _log_returns_sc_asym(N)) for N in range(62, 4000, 2))
    assert abs(ser / G_lattice_sc(zz, 0) - 1) < 2e-4
    assert {x: _fmt(G_ratio(x, 8), 4) for x in XI_CASES} == {0.96: "1.0231", 0.99: "1.0178"}
    assert abs(G_ratio(0.99, 16) - 1) < 0.02
    # (8.13)：rG = F(r√(z_c − z))，F(u) = (3/2π)e^{−6u}
    for x in XI_CASES:
        u = 5 * math.sqrt((1 - x) / 6)
        assert abs(5 * G_yukawa(x, 5) - F3(u)) < 1e-12
    assert _fmt(F3(0), 4) == "0.4775"

    # ---------------- 最陡下降 (8.15)–(8.20)
    assert sd_K(0.5, 2, 0.5) == 0.5 and sd_K(0.5, 2, 0.5, book=True) == 1.0
    assert abs(sd_K_dim(3) - 1.5) < 1e-12 and abs(sd_K_dim(1) - 0.5) < 1e-12 and abs(sd_K_dim(2) - 1) < 1e-12
    want_sd = {20: [(4, "0.6818", "0.6832", "0.4668"), (8, "0.2098", "0.2179", "0.0475")],
               100: [(8, "0.7282", "0.7285", "0.5306"), (12, "0.4894", "0.4902", "0.2403")]}
    for N in SD_N:
        got = [(x, _fmt(ratio_exact(N, x), 4), _fmt(sd_ratio(N, x, 0.5, 2, 0.5), 4), _fmt(sd_ratio(N, x, 0.5, 2, 0.5, book=True), 4))
               for x in SD_X[N]]
        assert got == want_sd[N], got
    assert abs(math.log(ratio_exact(2000, 60)) / (-60 ** 2 / (2 * 2001)) - 1) < 0.01          # 精確 → e^{−x²/2N}
    # (8.19)：閉式 δ ＝ 數值解（φ′ = 0），(8.20) 指數 ＝ φ(δ*)（N 大時）
    for nu in (0.5, 0.6, 0.75):
        N, r = 10 ** 8, 3 * (10 ** 8) ** nu
        dl, dn = sd_delta(N, r, nu, 2, 0.5), sd_delta_num(N, r, nu, 2, 0.5)
        assert abs(dl / dn - 1) < 1e-3, (nu, dl, dn)
        ex = -sd_K(nu, 2, 0.5) * (N + 1) ** (-nu / (1 - nu)) * r ** (1 / (1 - nu))
        assert abs(sd_phi_num(N, r, nu, 2, 0.5) / ex - 1) < 1e-3
        bk = -sd_K(nu, 2, 0.5, book=True) * (N + 1) ** (-nu / (1 - nu)) * r ** (1 / (1 - nu))
        assert abs(bk / ex - 1 / (1 - nu)) < 1e-12                                       # 書上多出 1/(1−ν)
    # 書上 (8.18) 的號：ν r δ^{ν−1}A′ − (N+1)/z = 0、A′ = −B → 兩項都負，無解
    assert all(0.5 * 1 * dl ** -0.5 * (-2) - 101 / (0.5 - dl) < 0 for dl in (1e-6, 1e-3, 0.1, 0.4))
    # Fisher 指數 1/(1−ν)
    assert [fisher(nu) for nu in FISHER_NU] == [2, F(5, 2), 4]
    for nu in FISHER_NU:
        assert abs(fisher_num(nu) - float(fisher(nu))) < 1e-9
    assert f"{1 / (1 - SW.NU_QUOTED[3]):.2f}" == "2.43"
    # Ex. 8.3–8.4：γ = (2−η)ν；條數 ∝ x^N N^p，x = 1/z_c、p = γ − 1
    assert gamma_exp(0, F(1, 2)) == 1 and gamma_exp(ETA_2D, NU_2D) == GAMMA_QUOTED == F(43, 32)
    assert growth(F(43, 32), F(1, 4)) == (4, F(11, 32)) and _fmt(F(11, 32), 5) == "0.34375"
    for g in (1, F(43, 32), F(7, 6)):
        assert abs(p_numeric(g, 0.25) - float(g - 1)) < 2e-3
    for g in (1, F(43, 32)):
        assert abs(math.exp(log_coeff_zN(40, g, 0.25)) / (GF.series_coeffs_power(float(g), 0.25, 40)[40] * 0.25 ** 40) - 1) < 1e-9
    assert abs(p_numeric(1, 0.25)) < 1e-12                                             # 普通漫步：(2d)^N 剛好
    # 比值法（二維自迴避）
    assert all(SW.saw_c(N) == SW.SAW_KNOWN_C[N] for N in range(15))
    want_r = {10: ("2.6592", "1.3552"), 12: ("2.6529", "1.3519"), 14: ("2.6483", "1.3496"), 16: ("2.6454", "1.3483")}
    for N in RATIO_N:
        _, mu_h, ga_h = ratio_method(N)
        assert (_fmt(mu_h, 4), _fmt(ga_h, 4)) == want_r[N], (N, mu_h, ga_h)
    mus = [ratio_method(N)[1] for N in RATIO_N]
    gas = [ratio_method(N)[2] for N in RATIO_N]
    assert all(m > MU_QUOTED for m in mus) and mus[-1] < mus[1] and abs(mus[-1] - MU_QUOTED) < 0.01
    assert all(abs(g - 43 / 32) < 0.015 for g in gas) and gas[-1] < gas[0]
    assert abs(SW.mu_eff(14) - SW.saw_c(14) / SW.saw_c(13)) < 1e-12
    # Ex. 8.5：小 r 的條數 ∝ z_c^{−N} N^{−dν} r^θ，θ = (γ − 1)/ν；F(x) ≈ F₀ + F₁ x^{d − 1/ν}
    assert theta(1, F(1, 2)) == 0 and theta(GAMMA_QUOTED, NU_2D) == F(11, 24)
    assert ex85_s(3, F(1, 2)) == 1 and ex85_s(2, NU_2D) == F(2, 3)
    for d in (1, 2, 3):
        assert abs(return_slope(d, 200) + d / 2) < 0.01, (d, return_slope(d, 200))
    # (z_c − z)^{ν s} 的係數 ∝ N^{−νs−1} = N^{−dν}：用冪級數係數核對
    for d, nu in ((3, F(1, 2)), (2, F(3, 4))):
        s = ex85_s(d, nu)
        a = -float(nu * s)                                       # (1 − z)^{−a} 的係數 Γ(N+a)/(Γ(a) N!)
        lc = lambda N: math.lgamma(N + a) - math.lgamma(N + 1)
        assert abs((lc(2000) - lc(1000)) / math.log(2) + float(d * nu)) < 2e-3

    # ---------------- §8.1.3 重新標度
    for d in (1, 2, 3, 4):
        for l1, l2 in ((0.5, 0.3), (2.0, 5.0)):
            K, V = rescale_check(l1, l2, d)
            assert abs(K - kin_factor(l1, l2)) < 1e-9 and abs(V / int_factor(l1, l2, d) - 1) < 1e-9
    # Ex. 8.6：a′ = a l1/√l2 吸收動能項；l1 = √l2 → 自迴避項乘 l2^{(d−4)/2}（書上的 l1^d/l2 會得 l2^{(d−2)/2}）
    for d in range(1, 6):
        l2 = 0.37
        assert abs(int_factor(math.sqrt(l2), l2, d) - l2 ** ((d - 4) / 2)) < 1e-12
        assert abs(int_factor_book86(math.sqrt(l2), l2, d) - l2 ** ((d - 2) / 2)) < 1e-12
        assert abs(kin_factor(math.sqrt(l2), l2) - 1) < 1e-12
    assert [sa_amplitude(SA_N, d) for d in range(1, 6)] == [1000, 100, 10, 1, F(1, 10)]
    assert all(abs(float(sa_amplitude(SA_N, d)) - int_factor(math.sqrt(1 / SA_N), 1 / SA_N, d)) < 1e-9 for d in range(1, 6))
    assert [d for d in range(1, 9) if sa_amplitude(SA_N, d) == 1] == [UPPER_CRITICAL_D]
    # Ex. 8.7：l2 = 1/N → R′ = l1 R = R/√N 是一步的長度 → R ∝ √N
    assert abs(math.sqrt(1 / SA_N) * 10 - 1) < 1e-12
    # (8.24)(8.25)：l1 = l2^{3/(d+2)} → ν = 3/(d+2)（＝ walk_saw.flory_nu）
    for d in (1, 2, 3, 4):
        for l2 in (0.01, 0.3):
            assert abs(flory_l1(l2, d) / flory_l1_num(l2, d) - 1) < 1e-9
            l1 = flory_l1(l2, d)
            assert abs(kin_factor(l1, l2) / int_factor(l1, l2, d) - 1) < 1e-9
        assert nu_rescale(d) == SW.flory_nu(d)
    assert [nu_rescale(d) for d in range(1, 5)] == [1, F(3, 4), F(3, 5), F(1, 2)]
    assert _fmt(flory_l1(0.01, 3), 4) == "0.0631" and _fmt(1 / flory_l1(0.01, 3), 2) == "15.85"
    # Ex. 8.8：|r|^{−q}：因子 l1^q/l2² → l2^{(q−4)/2}；與 d 無關
    for q in (1, 2, 3, 5):
        l2 = 0.2
        assert abs(ex88_factor(math.sqrt(l2), l2, q) - l2 ** ((q - 4) / 2)) < 1e-12
        assert abs(ex88_amplitude(100, q) - 100 ** ((4 - q) / 2)) < 1e-9
    assert ex88_relevant(3, 1) and ex88_relevant(3, 6) and not ex88_relevant(5, 5) and ex88_relevant(6, 3)
    assert [ex88_nu(q) for q in (1, 2, 3, 4)] == [1, F(3, 4), F(3, 5), F(1, 2)]
    return True


def _log_returns_sc_asym(N):
    """C(N;0) 的漸近 2·6^N (3/(2πN))^{3/2}（偶數 N；check 的級數尾巴用）。"""
    return math.log(2) + N * math.log(6) + 1.5 * math.log(3 / (2 * math.pi * N))


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
    print(f"walk_scale（第 8 章）：長尾巴步伐的 χ(q) 與寬度 N^(1/y)、尾巴 r^(−5+δ) 與換變數、生成函數的標度塌縮與實空間 ξ、"
          f"最陡下降與 Fisher 指數、γ = (2−η)ν 與比值法、重新標度與上臨界維度 4、Flory 指數全部驗證通過（{time.time() - t0:.1f} 秒）")
