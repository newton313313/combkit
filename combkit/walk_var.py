# 由 combkit/tools/sync_from_walk017.py 從 walk017 專案的 codes/utils/walk_var.py 複製；請改 walk017 那份再同步。
"""
walk_var.py —— 《Elements of the Random Walk》第 5 章 Variations on the random walk 的數學
（Python＋numpy／scipy.special，不 import manim；只依賴同資料夾的 walk_bc（→ walk_gf、walk_rec），以便原樣複製進 combkit）

約定（書上式號與筆誤見 ref/md/ch05_formulas.md）：
  一維格子步長 l = 1；偏向漫步向右的機率 p；有慣性的漫步「和上一步同方向」的機率 Q，起點 x = 0、第一步左右各半（L = R = ½）。
  傅立葉照書上 Σ_x P(x) e^{iqx}。電報方程（一維）：速率 v、反轉率 γ（右行者每單位時間以 γ 變左行者）；t = 0 時全部在 0、左右各半。
  二維、三維「跑一跑、轉個彎」（run-and-tumble）：速率 v、以速率 γ 換成完全隨機的新方向（(5.69) Γ = γ − γ|0⟩⟨0|）。
  密度一律用標準歸一化（∫ 密度 d^d r = 1）；書上二維 (5.84)–(5.89) 是 2πρ（且有 2π 的筆誤）、三維 (5.98)–(5.103) 是 2π²ρ
  （所以書上說「r²Σ 積分＝π/2」）。徑向機率密度 P_d(r, t)：二維 2πrρ、三維 4πr²ρ（∫ P dr ＋ 尖峰 e^{−γt} = 1）。
  連續時間漫步：kind = 'poisson'（q(t) = p0 常數）或 'gamma2'（等待時間密度 t e^{−t}，q(t) = t/(1+t)，平均等 2）。
  函數多半接受 numpy 陣列（標「向量」者）。

§5.1 偏向的漫步 (5.1)–(5.10)，Fig. 5.1，Ex. 5.1
    biased_pmf(N, p)                {x: P}，x = 2n − N（(5.1)；p 是 Fraction 時精確）；biased_dp(N, p) 逐步遞迴 (5.2)（第二種）
    moments(P)                      {x: P} → (平均, 變異數)；biased_mean(N, p) = (2p−1)N、biased_var(N, p) = 4Np(1−p)
    biased_gauss(x, N, p)           (5.7) 高斯密度（向量）；格子機率 ≈ 2 × 它（奇偶性）
    taylor_residual(p, book=False)  (5.3) 用平滑測試函數核對一階項：更正 l(1−2p)∂P/∂y 殘差 ~1e−6；書上 l(2p−1) 殘差大
    diffusion_D(p, kind)            Ex. 5.1：'naive' 直接寫擴散方程 l²/2（變異數 N）；'careful' l²/2 − v²/2 = 2p(1−p)（v = 2p−1）
    width_ratio(p) = √(4p(1−p))     精確寬度 ÷ 擴散方程的寬度
    shift_pmf(N, b)                 (5.9)(5.10) 一維版「先被風吹 b、再擲公正硬幣」：P0(N; x − Nb)（list[(x, P)]）；shift_dp 逐步（第二種）
    shift_fourier(N, b, x)          (5.10)：(1/2π)∫(e^{−iqb}cos q)^N e^{iqx} dq（b 整數，均分求積）＝ C(N, (x−Nb+N)/2)/2^N
§5.2 有慣性的漫步 (5.11)–(5.38)，Ex. 5.2
    persistent_dp(N, Q)             {x: (L, R)}（(5.11)(5.12)；Q 是 Fraction 時精確）；persistent_pmf(N, Q) = {x: L+R}
    persistent_np(N, Q)             同上（numpy，大 N）→ (xs, L, R)
    persistent_var(N, Q)            精確變異數（DP）；persistent_var_formula(N, Q) = N(1+r)/(1−r) − 2r(1−r^N)/(1−r)²，r = 2Q−1
    persistent_var_asym(N, Q)       (5.27)：NQ/(1−Q)；persistence_length(Q) = Q/(1−Q)（(5.30) 的 n，單位 l）
    T_matrix(q, Q)、T_total(N, q, Q) (5.15)(5.16)：(1,1)·T^N·(½,½)
    eigs(q, Q)                      (5.17)(5.18)（complex）；eig_residual(lam, q, Q, book)：λ² − 2Q cos q λ + (2Q−1)（book：+(1−2Q)，(5.37)）
    lambda1_small(q, Q)             (5.23)：1 − Qq²/(2(1−Q))
    total_fourier(N, q, Q)          Σ_x P e^{iqx}（DP）；total_eig(N, q, Q) = c_N − (2Q−1)cos q·c_{N−1}，c_k = (λ1^{k+1}−λ2^{k+1})/(λ1−λ2)
    recursion2_residual(N, q, Q)    (5.35) R(N+2) − 2Q cos q R(N+1) − (1−2Q)R(N)（DP 的 R 分量）
    gauss_persistent(x, N, Q)       (5.27) 的高斯密度（向量）
§5.2.2 電報方程 (5.39)–(5.64)，Fig. 5.2–5.3，Ex. 5.3–5.6
    tele_cont(x, t, v, gam)         (5.61) 更正後的連續部分 (γe^{−γt}/2v)[I0(γs) + t I1(γs)/s]，s = √(t²−x²/v²)（向量；|x| ≥ vt 為 0）
    tele_book(x, t, v, gam)         照書上 (5.60)(5.61)（少 1/v）；tele_spike(t, gam) = e^{−γt}/2（每個尖峰）
    tele_mass(t, v, gam, book)      ∫Σ_cont dx + e^{−γt}（更正 = 1）
    tele_var(t, v, gam)             (v²/γ)[t − (1−e^{−2γt})/2γ]；tele_var_numeric 由 Σ 積分
    tele_fourier(q, t, v, gam)      (5.52) 反轉的閉式 e^{−γt}[cosh ωt + γ sinh ωt/ω]，ω² = γ²−q²v²；tele_fourier_numeric 由 Σ 積分；
                                    tele_laplace(q, lam, v, gam) = (λ+2γ)/(λ²+2γλ+q²v²)
    tele_D(v, gam) = v²/2γ、tele_lp(v, gam) = v/γ、tele_gauss(x, t, v, gam)   大 t 的高斯（Ex. 5.5–5.6）
    tele_lattice(t, dt, v, gam)     Ex. 5.3：Q = 1−γδt、l = vδt 的有慣性漫步 → (xs, P)；格點密度 ≈ P/(2l)
§5.2.3–5.2.5 二維、三維 (5.65)–(5.106)，Fig. 5.4–5.6，Ex. 5.7
    angle_avg_2d(lam, q, v, gam)    (5.77) 數值角度平均；angle_avg_2d_closed = 1/√((λ+γ)²+q²v²)
    angle_avg_3d(…)、angle_avg_3d_closed = arctan(vq/(λ+γ))/(vq)   (5.93)
    sm_residual(n, lam, q)          Ex. 5.7：角度離散成 n 格，(5.75) × (5.73) − I 的最大誤差
    rt2_cont(r, t, v, gam)          二維連續部分（更正）(γ/2πv) e^{−γt+γs/v}/s，s = √(v²t²−r²)（向量）
    rt2_pole、rt2_cut               ＝ (γ/πv)e^{−γt}cosh(γs/v)/s 與 −(γ/2πv)e^{−γt−γs/v}/s（更正後的 (5.84)、(5.89) ÷ 2π）
    rt2_book_pole、rt2_book_cut     照書上（2πρ 的歸一化）；rt2_book_spike = 2π（書上 δ(r−vt)/r 的係數 ÷ e^{−γt}；應為 1）
    rt2_radial(r, t) = 2πrρ；rt2_mass(t)；rt2_hankel(q, t) vs rt2_laplace(q, lam)（Talbot）
    hankel_pair(r)                  (5.88)：∫q J0(qr)[∫r'J0(qr')e^{−r'²/2}dr']dq（數值）= e^{−r²/2}（沒有 2π）
    rt3_pole_book(r, t)             (5.101)（書上歸一化 2π²ρ，式子本身對）；rt3_cut_book(r, t) (5.103) 第一項；rt3_once_book(r, t) 第三項
    rt3_5106(r, t, book)            (5.106)：書上 e^{−γut cot u}；更正 e^{+γut cot u}（＝ −(5.101)）
    rt3_density(r, t)、rt3_radial(r, t) = 4πr²ρ、rt3_r2sigma_book = (π/2)·rt3_radial（Fig. 5.6 的縱軸）
    rt3_once_radial(r, t)           只轉一次彎的人：γe^{−γt} r ln((vt+r)/(vt−r))/(v²t)（總量 γt e^{−γt}）
    rt3_mass(t)；rt3_fourier(q, t) vs rt3_laplace(q, lam)（Talbot）
    rt_mc(dim, t, M, seed, v, gam)  固定種子的粒子模擬（numpy）：一維回傳位置、二三維回傳距離；mc_fraction(vals, a, b)
    radial_fraction(dim, a, b, t)   公式：a ≤ r < b 的比例（一維用 |x|，不含尖峰）
    talbot(F, t, M=32)              數值 Laplace 反轉（Abate–Valkó 固定 Talbot 圍道）
§5.3 連續時間漫步 (5.107)–(5.133)，Ex. 5.8–5.9
    hazard(t, kind)、survival(t, kind)、wait_density(t, kind)   q(t)、a(t) = e^{−∫q}（(5.108)(5.115)）、a·q；survival_numeric 數值積分
    AB(lam, kind)、AB_numeric(lam, kind)    (5.117)(5.118) 的閉式與數值積分；A = 1 − λB（(5.119)）
    steps_prob(m, t, kind)          走了 m 步的機率：Poisson (p0t)^m e^{−p0t}/m!；gamma2 e^{−t}(t^{2m}/(2m)! + t^{2m+1}/(2m+1)!)
    steps_prob_laplace(m, t, kind)  Talbot 反轉 A^m B（(5.121)）；steps_prob_conv(m, t, kind) 數值卷積 (5.110)(5.111)
    mean_steps(t, kind)             Poisson p0t；gamma2 t/2 − 1/4 + e^{−2t}/4；mean_wait(kind) = B(0)
    ctrw_lattice(x, t, kind)        一維格子（χ = cos q）：Σ_m P_m(t)·C(m,(m+x)/2)/2^m；ctrw_bessel(x, t, p0) = e^{−p0t}I_|x|(p0t)（Ex. 5.8）
    ctrw_var(t, kind)、ctrw_var_numeric   精確變異數＝平均步數；ctrw_gauss(x, t, kind) (5.133)÷(2π)：√(B0/4πCt)e^{−B0x²/4Ct}，C = ½
    ctrw_transform(lam, q, kind)    (5.128) B/(1−χA)；ex58(lam, q, p0) = 1/(λ + p0(1−cos q))
    mortal_total(t, p0, p1) = e^{−p1t}（Ex. 5.9）；mortal_series 逐步數相加；mortal_laplace Talbot 反轉 B/(1−A)

check() 驗證全部（每個結果至少兩種獨立算法、書上的數字、筆誤的更正、plan_ch05.md 場景卡與講義的每個數字）；直接執行本檔會跑 check()（約 1 秒）。
最慢的呼叫：rt3_radial 每點約 0.1 毫秒（一條 120 點的曲線約 0.01 秒）、rt3_mass 約 0.03 秒、rt_mc(3, 3.0, 100000) 約 0.05 秒，其餘都很快。
"""
import cmath
import math
from fractions import Fraction

import numpy as np
from scipy import special as _sp

if __package__:
    from . import walk_bc as BC
else:                                   # python utils/walk_var.py
    import walk_bc as BC

PI = math.pi


def _simps(y, h):
    """均勻格點的 Simpson（len(y) 為奇數）。"""
    y = np.asarray(y)
    return h / 3 * (y[0] + y[-1] + 4 * y[1:-1:2].sum() + 2 * y[2:-1:2].sum())


def talbot(F, t: float, M: int = 32) -> float:
    """數值 Laplace 反轉 f(t)（Abate–Valkó 固定 Talbot；F 在圍道內只要奇點都在左邊、λ 大時解析）。"""
    r = 2 * M / (5 * t)
    s = 0.5 * (F(complex(r)) * cmath.exp(r * t)).real
    for k in range(1, M):
        th = k * PI / M
        cot = math.cos(th) / math.sin(th)
        lam = complex(r * th * cot, r * th)
        sig = th + (th * cot - 1) * cot
        s += (cmath.exp(lam * t) * F(lam) * complex(1, sig)).real
    return r / M * s


# =========================================================================== §5.1 偏向的漫步

def biased_pmf(N: int, p) -> dict:
    """(5.1)：{x = 2n − N: C(N,n) p^n (1−p)^{N−n}}。"""
    return {2 * n - N: math.comb(N, n) * p ** n * (1 - p) ** (N - n) for n in range(N + 1)}


def biased_dp(N: int, p) -> dict:
    """(5.2) 逐步遞迴：P(N+1; y) = p P(N; y−1) + (1−p) P(N; y+1)。"""
    cur = {0: p * 0 + 1}
    for _ in range(N):
        nxt = {}
        for x, w in cur.items():
            nxt[x + 1] = nxt.get(x + 1, 0) + p * w
            nxt[x - 1] = nxt.get(x - 1, 0) + (1 - p) * w
        cur = nxt
    return cur


def moments(P: dict):
    m = sum(x * w for x, w in P.items())
    return m, sum(x * x * w for x, w in P.items()) - m * m


def biased_mean(N, p, l=1):
    return (2 * p - 1) * N * l


def biased_var(N, p, l=1):
    return 4 * N * p * (1 - p) * l * l


def biased_gauss(x, N, p, l=1.0):
    """(5.7)：exp[−(x − (2p−1)Nl)²/(8Nl²p(1−p))]/√(8πNp(1−p)l²)。"""
    x = np.asarray(x, float)
    s2 = 8 * N * p * (1 - p) * l * l
    return np.exp(-(x - (2 * p - 1) * N * l) ** 2 / s2) / math.sqrt(PI * s2)


def taylor_residual(p: float, book: bool = False, l: float = 1.0, s: float = 30.0, y: float = 7.0) -> float:
    """(5.3)：|pP(y−l) + (1−p)P(y+l) − P(y) − [c ∂P + (l²/2)∂²P]|，P = e^{−y²/2s²}；c = l(1−2p)（更正）或 l(2p−1)（書上）。"""
    P = lambda u: math.exp(-u * u / (2 * s * s))
    d1 = -y / s ** 2 * P(y)
    d2 = (y * y / s ** 4 - 1 / s ** 2) * P(y)
    c = l * (2 * p - 1) if book else l * (1 - 2 * p)
    return abs(p * P(y - l) + (1 - p) * P(y + l) - P(y) - (c * d1 + l * l / 2 * d2))


def diffusion_D(p, l=1.0, kind="careful"):
    """Ex. 5.1：'naive' l²/2（直接把 P(N+1) − P(N) 換成 ∂P/∂N）；'careful' 加上 −½∂²P/∂N² ≈ −(v²/2)∂²P/∂y² → 2p(1−p)l²。"""
    v = (2 * p - 1) * l
    return l * l / 2 if kind == "naive" else l * l / 2 - v * v / 2


def width_ratio(p) -> float:
    return math.sqrt(4 * p * (1 - p))


def shift_pmf(N: int, b: float) -> list:
    """(5.10) 一維：每步先走 b、再 ±1（公正）→ [(x, P)]，x = k + Nb，P = C(N,(k+N)/2)/2^N。"""
    return [(k + N * b, math.comb(N, (k + N) // 2) / 2 ** N) for k in range(-N, N + 1, 2)]


def shift_dp(N: int, b: float) -> list:
    cur = {0.0: 1.0}
    for _ in range(N):
        nxt = {}
        for x, w in cur.items():
            for d in (-1, 1):
                y = round(x + b + d, 9)
                nxt[y] = nxt.get(y, 0) + w / 2
        cur = nxt
    return sorted(cur.items())


def shift_fourier(N: int, b: int, x: int) -> float:
    """(5.10)：(1/2π)∫_{−π}^{π} χ(q)^N e^{iqx} dq，χ = e^{−iqb} cos q（(5.9)；b 整數，均分求積精確）。"""
    K = 2 * (N * (abs(b) + 1) + abs(x)) + 2
    s = 0
    for j in range(K):
        q = -PI + 2 * PI * j / K
        s += (cmath.exp(-1j * q * b) * math.cos(q)) ** N * cmath.exp(1j * q * x)
    return (s / K).real


# =========================================================================== §5.2 有慣性的漫步

def persistent_dp(N: int, Q) -> dict:
    """(5.11)(5.12)：{x: (L, R)}，L/R = 最後一步向左／向右；N = 0 時 L = R = ½ 在 0。"""
    half = Fraction(1, 2) if isinstance(Q, Fraction) else 0.5
    cur = {0: (half, half)}
    for _ in range(N):
        nxt = {}
        for x, (Lw, Rw) in cur.items():
            a = nxt.get(x - 1, (0, 0))
            nxt[x - 1] = (a[0] + Q * Lw + (1 - Q) * Rw, a[1])
            b = nxt.get(x + 1, (0, 0))
            nxt[x + 1] = (b[0], b[1] + (1 - Q) * Lw + Q * Rw)
        cur = nxt
    return cur


def persistent_pmf(N: int, Q) -> dict:
    return {x: Lw + Rw for x, (Lw, Rw) in persistent_dp(N, Q).items()}


def persistent_np(N: int, Q: float):
    """numpy 版：(xs, L, R)，xs = −N..N。"""
    L = np.zeros(2 * N + 3)
    R = np.zeros(2 * N + 3)
    L[N + 1] = R[N + 1] = 0.5
    for _ in range(N):
        L2 = np.zeros_like(L)
        R2 = np.zeros_like(R)
        L2[:-1] = Q * L[1:] + (1 - Q) * R[1:]
        R2[1:] = (1 - Q) * L[:-1] + Q * R[:-1]
        L, R = L2, R2
    return np.arange(-N - 1, N + 2), L, R


def persistent_var(N: int, Q) -> float:
    if isinstance(Q, Fraction) or N <= 60:
        return moments(persistent_pmf(N, Q))[1]
    xs, L, R = persistent_np(N, Q)
    P = L + R
    return float((xs ** 2 * P).sum() - ((xs * P).sum()) ** 2)


def persistent_var_formula(N: int, Q):
    r = 2 * Q - 1
    if r == 1:
        return N * N
    return N * (1 + r) / (1 - r) - 2 * r * (1 - r ** N) / (1 - r) ** 2


def persistent_var_asym(N, Q):
    return N * Q / (1 - Q)


def persistence_length(Q):
    return Q / (1 - Q)


def T_matrix(q: float, Q: float):
    em, ep = cmath.exp(-1j * q), cmath.exp(1j * q)
    return [[Q * em, (1 - Q) * em], [(1 - Q) * ep, Q * ep]]


def T_total(N: int, q: float, Q: float) -> complex:
    T = T_matrix(q, Q)
    a, b = 0.5, 0.5
    for _ in range(N):
        a, b = T[0][0] * a + T[0][1] * b, T[1][0] * a + T[1][1] * b
    return a + b


def eigs(q: float, Q: float):
    c = Q * math.cos(q)
    d = cmath.sqrt((1 - Q) ** 2 - Q * Q * math.sin(q) ** 2)
    return c + d, c - d


def eig_residual(lam, q, Q, book=False):
    return abs(lam * lam - 2 * Q * math.cos(q) * lam + ((1 - 2 * Q) if book else (2 * Q - 1)))


def lambda1_small(q, Q):
    return 1 - Q * q * q / (2 * (1 - Q))


def total_fourier(N: int, q: float, Q: float, part: str = "LR") -> complex:
    xs, L, R = persistent_np(N, Q)
    P = (L if "L" in part else 0) + (R if "R" in part else 0)
    return complex((P * np.exp(1j * q * xs)).sum())


def total_eig(N: int, q: float, Q: float) -> complex:
    l1, l2 = eigs(q, Q)
    c = lambda k: (l1 ** (k + 1) - l2 ** (k + 1)) / (l1 - l2)
    return 1 if N == 0 else c(N) - (2 * Q - 1) * math.cos(q) * c(N - 1)


def recursion2_residual(N: int, q: float, Q: float, book: bool = False) -> float:
    RF = [total_fourier(n, q, Q, "R") for n in (N, N + 1, N + 2)]
    return abs(RF[2] - 2 * Q * math.cos(q) * RF[1] - (1 - 2 * Q) * RF[0])


def gauss_persistent(x, N, Q):
    x = np.asarray(x, float)
    s2 = N * Q / (1 - Q)
    return np.exp(-x * x / (2 * s2)) / math.sqrt(2 * PI * s2)


# =========================================================================== §5.2.2 電報方程

def tele_cont(x, t, v=1.0, gam=1.0):
    """(5.61) 更正後的連續部分（向量）：(γe^{−γt}/2v)[I0(γs) + t I1(γs)/s]，s = √(t² − x²/v²)。"""
    x = np.asarray(x, float)
    s = np.sqrt(np.maximum(t * t - (x / v) ** 2, 0.0))
    y = gam * s
    ratio = np.where(s > 1e-12, _sp.i1e(y) / np.maximum(s, 1e-300), gam / 2)
    val = gam / (2 * v) * np.exp(y - gam * t) * (_sp.i0e(y) + t * ratio)
    return np.where(np.abs(x) < v * t, val, 0.0)


def tele_book(x, t, v=1.0, gam=1.0):
    """照書上 (5.60)(5.61)：前面是 ½ 而不是 1/(2v)。"""
    return v * tele_cont(x, t, v, gam)


def tele_spike(t, gam=1.0):
    return math.exp(-gam * t) / 2


def tele_mass(t, v=1.0, gam=1.0, book=False, n=2000):
    xs = np.linspace(-v * t, v * t, n + 1)
    f = (tele_book if book else tele_cont)(xs, t, v, gam)
    f[0] = f[-1] = (tele_book if book else tele_cont)(v * t * (1 - 1e-12), t, v, gam)
    return float(_simps(f, xs[1] - xs[0])) + math.exp(-gam * t)


def tele_var(t, v=1.0, gam=1.0):
    return v * v / gam * (t - (1 - math.exp(-2 * gam * t)) / (2 * gam))


def tele_var_numeric(t, v=1.0, gam=1.0, n=2000):
    xs = np.linspace(-v * t, v * t, n + 1)
    f = tele_cont(xs, t, v, gam)
    f[0] = f[-1] = float(tele_cont(v * t * (1 - 1e-12), t, v, gam))
    return float(_simps(xs ** 2 * f, xs[1] - xs[0])) + math.exp(-gam * t) * (v * t) ** 2


def tele_fourier(q, t, v=1.0, gam=1.0):
    """σ(q, t) = e^{−γt}[cosh ωt + γ sinh(ωt)/ω]，ω² = γ² − q²v²（(5.52) 的 Laplace 反轉）。"""
    w2 = gam * gam - (q * v) ** 2
    if abs(w2) < 1e-14:
        return math.exp(-gam * t) * (1 + gam * t)
    if w2 > 0:
        w = math.sqrt(w2)
        return math.exp(-gam * t) * (math.cosh(w * t) + gam * math.sinh(w * t) / w)
    w = math.sqrt(-w2)
    return math.exp(-gam * t) * (math.cos(w * t) + gam * math.sin(w * t) / w)


def tele_fourier_numeric(q, t, v=1.0, gam=1.0, n=4000):
    xs = np.linspace(-v * t, v * t, n + 1)
    f = tele_cont(xs, t, v, gam)
    f[0] = f[-1] = float(tele_cont(v * t * (1 - 1e-12), t, v, gam))
    return float(_simps(f * np.cos(q * xs), xs[1] - xs[0])) + math.exp(-gam * t) * math.cos(q * v * t)


def tele_laplace(q, lam, v=1.0, gam=1.0):
    return (lam + 2 * gam) / (lam * lam + 2 * gam * lam + (q * v) ** 2)


def tele_D(v=1.0, gam=1.0):
    return v * v / (2 * gam)


def tele_lp(v=1.0, gam=1.0):
    return v / gam


def tele_gauss(x, t, v=1.0, gam=1.0):
    D = tele_D(v, gam)
    x = np.asarray(x, float)
    return np.exp(-x * x / (4 * D * t)) / math.sqrt(4 * PI * D * t)


def tele_lattice(t, dt, v=1.0, gam=1.0):
    """Ex. 5.3 的格子版：步長 l = vδt、Q = 1 − γδt、N = t/δt 步 → (xs, P)（xs 是實際座標）。"""
    N = int(round(t / dt))
    xs, L, R = persistent_np(N, 1 - gam * dt)
    return xs * v * dt, L + R


# =========================================================================== §5.2.3–5.2.5 二維、三維

def _sqrt_seg(mu, a):
    """√(μ² + a²)，分支切割放在 μ ∈ i[−a, a]（λ 大時 ≈ μ）。"""
    mu = complex(mu)
    return mu * cmath.sqrt(1 + (a / mu) ** 2) if mu != 0 else complex(a)


def angle_avg_2d(lam, q, v=1.0, gam=1.0, n=512):
    th = 2 * PI * np.arange(n) / n
    return complex(np.mean(1 / (lam + gam + 1j * v * q * np.cos(th))))


def angle_avg_2d_closed(lam, q, v=1.0, gam=1.0):
    return 1 / _sqrt_seg(lam + gam, q * v)


def angle_avg_3d(lam, q, v=1.0, gam=1.0, n=2000):
    u = np.linspace(-1, 1, n + 1)
    return complex(_simps(1 / (lam + gam + 1j * v * q * u), u[1] - u[0]) / 2)


def angle_avg_3d_closed(lam, q, v=1.0, gam=1.0):
    mu = complex(lam + gam)
    return 1 / mu if q == 0 else cmath.atan(v * q / mu) / (v * q)


def sm_residual(n=64, lam=0.7, q=1.3, v=1.0, gam=1.0) -> float:
    """Ex. 5.7：A = diag(λ+γ+ivq cos θ_k)，|0⟩ = (1,…,1)/√n；(5.75) 乘 (5.73) 與單位矩陣的最大差。"""
    th = 2 * PI * np.arange(n) / n
    A = np.diag(lam + gam + 1j * v * q * np.cos(th))
    z = np.ones((n, 1)) / math.sqrt(n)
    P0 = z @ z.T
    Ai = np.diag(1 / np.diag(A))
    inv = Ai + gam * Ai @ P0 @ Ai / (1 - gam * (z.T @ Ai @ z)[0, 0])
    return float(np.abs((A - gam * P0) @ inv - np.eye(n)).max())


def _s(r, t, v):
    return np.sqrt(np.maximum((v * t) ** 2 - np.asarray(r, float) ** 2, 1e-300))


def rt2_cont(r, t, v=1.0, gam=1.0):
    """二維密度的連續部分（更正）：(γ/2πv) e^{−γt + γs/v}/s，s = √(v²t² − r²)（向量；r ≥ vt 為 0）。"""
    r = np.asarray(r, float)
    s = _s(r, t, v)
    return np.where(r < v * t, gam / (2 * PI * v) * np.exp(-gam * t + gam * s / v) / s, 0.0)


def rt2_pole(r, t, v=1.0, gam=1.0):
    s = _s(r, t, v)
    return gam / (PI * v) * np.exp(-gam * t) * np.cosh(gam * s / v) / s


def rt2_cut(r, t, v=1.0, gam=1.0):
    s = _s(r, t, v)
    return -gam / (2 * PI * v) * np.exp(-gam * t - gam * s / v) / s


def rt2_book_pole(r, t, v=1.0, gam=1.0):
    """照書上 (5.84)：2γ√(π/2) e^{−γt} cosh(γ√(v²t²−r²))/(v√(v²t²−r²))（書上的歸一化是 2πρ）。"""
    s = _s(r, t, v)
    return 2 * gam * math.sqrt(PI / 2) * np.exp(-gam * t) * np.cosh(gam * s) / (v * s)


def rt2_book_cut(r, t, v=1.0, gam=1.0):
    """照書上 (5.89) 的連續部分：−2π γ e^{−γt} e^{−γs/v}/(v s)。"""
    s = _s(r, t, v)
    return -2 * PI * gam * np.exp(-gam * t - gam * s / v) / (v * s)


rt2_book_spike = 2 * PI         # 書上 (5.89)：2π e^{−γt} δ(r−vt)/r；在 2πρ 的歸一化下應為 1


def rt2_radial(r, t, v=1.0, gam=1.0):
    return 2 * PI * np.asarray(r, float) * rt2_cont(r, t, v, gam)


def _phi_grid(t, v, n):
    ph = np.linspace(0, PI / 2, n + 1)
    ph[-1] = PI / 2 - 1e-7                                  # 端點：r 略小於 vt（被積函數在 φ 變數裡有限）
    return ph, v * t * np.sin(ph), v * t * np.cos(ph)


def rt2_mass(t, v=1.0, gam=1.0, n=400, book=False):
    """∫ 2πrρ dr ＋ 尖峰（r = vt sin φ）。book=True：用書上 (5.84)+(5.89)（2πρ 的歸一化：∫ r(…) dr）。"""
    ph, r, jac = _phi_grid(t, v, n)
    if book:
        f = r * (rt2_book_pole(r, t, v, gam) + rt2_book_cut(r, t, v, gam)) * jac
        return float(_simps(f, ph[1] - ph[0])) + rt2_book_spike * math.exp(-gam * t)
    f = rt2_radial(r, t, v, gam) * jac
    return float(_simps(f, ph[1] - ph[0])) + math.exp(-gam * t)


def rt2_hankel(q, t, v=1.0, gam=1.0, n=400):
    ph, r, jac = _phi_grid(t, v, n)
    f = rt2_radial(r, t, v, gam) * _sp.j0(q * r) * jac
    return float(_simps(f, ph[1] - ph[0])) + math.exp(-gam * t) * float(_sp.j0(q * v * t))


def rt2_laplace(q, lam, v=1.0, gam=1.0):
    """(5.78)：1/(√((λ+γ)²+q²v²) − γ)。"""
    return 1 / (_sqrt_seg(lam + gam, q * v) - gam)


def hankel_pair(r, n=3000, qmax=14.0):
    """(5.88) 的核對：f(r') = e^{−r'²/2} → F(q) = ∫ r'J0(qr')f dr'（數值）→ ∫ q J0(qr) F(q) dq（數值）；應回到 f(r)。"""
    rr = np.linspace(0, qmax, n + 1)
    F = np.array([_simps(rr * _sp.j0(q * rr) * np.exp(-rr ** 2 / 2), rr[1] - rr[0]) for q in rr[::10]])
    qq = rr[::10]
    return float(_simps(qq * _sp.j0(qq * r) * F, qq[1] - qq[0]))


_XP = np.linspace(1e-9, PI / 2, 801)


def rt3_pole_book(r, t, v=1.0, gam=1.0):
    """(5.101)（書上歸一化 2π²ρ）：(γ/v)² ∫_0^{π/2} x³/(r sin²x) e^{−γt} e^{γxt cot x} sin(γxr/v) dx。"""
    x = _XP
    f = x ** 3 / np.sin(x) ** 2 * np.exp(gam * t * (x / np.tan(x) - 1)) * np.sin(gam * x * r / v)
    return (gam / v) ** 2 * float(_simps(f, x[1] - x[0])) / r


def rt3_5106(r, t, v=1.0, gam=1.0, book=True):
    """(5.106)：−e^{−γt}∫du/sin²u (γ/v)² u³ sin(γur/v)/r e^{∓γut cot u}（書上 −；更正 +，即 −(5.101)）。"""
    u = _XP
    sg = -1 if book else 1
    f = u ** 3 / np.sin(u) ** 2 * np.exp(-gam * t + sg * gam * u * t / np.tan(u)) * np.sin(gam * u * r / v)
    return -(gam / v) ** 2 * float(_simps(f, u[1] - u[0])) / r


def rt3_cut_book(r, t, v=1.0, gam=1.0, n=600):
    """(5.103) 第一項（x = tanh y）：e^{−γt}(1/r)(γ/v)² Re∫ sech²y (π/2 − iy)³ exp[−i(r − vt tanh y)(γ/v)(π/2 − iy)] dy。"""
    Y = 18.0 if r >= v * t else min(18.0, math.atanh(r / (v * t)))
    y = np.linspace(0, Y, n + 1)
    B = PI / 2 - 1j * y
    f = (B ** 3 * np.exp(-1j * (r - v * t * np.tanh(y)) * gam / v * B) / np.cosh(y) ** 2).real
    return math.exp(-gam * t) * (gam / v) ** 2 * float(_simps(f, y[1] - y[0])) / r


def rt3_once_book(r, t, v=1.0, gam=1.0):
    """(5.103) 第三項：e^{−γt}(γ/v)(1/(rvt))(π/2) ln((1+r/vt)/(1−r/vt))（只轉一次彎的人；書上歸一化）。"""
    u = r / (v * t)
    return math.exp(-gam * t) * gam / v / (r * v * t) * PI / 2 * math.log((1 + u) / (1 - u)) if u < 1 else 0.0


def rt3_density(r, t, v=1.0, gam=1.0):
    """三維密度的連續部分（標準歸一化）＝ [(5.101) + (5.103)]/(2π²)；r ≥ vt 為 0。"""
    if r >= v * t:
        return 0.0
    r = max(r, 1e-6)
    return (rt3_pole_book(r, t, v, gam) + rt3_cut_book(r, t, v, gam) + rt3_once_book(r, t, v, gam)) / (2 * PI * PI)


def rt3_radial(r, t, v=1.0, gam=1.0):
    """徑向機率密度 4πr²ρ（∫_0^{vt} dr ＝ 1 − e^{−γt}）；陣列逐點算。"""
    if np.ndim(r):
        return np.array([rt3_radial(float(x), t, v, gam) for x in np.ravel(r)]).reshape(np.shape(r))
    return 4 * PI * r * r * rt3_density(r, t, v, gam)


def rt3_r2sigma_book(r, t, v=1.0, gam=1.0):
    """書上 Fig. 5.6 的縱軸 r²Σ（Σ = 2π²ρ）＝ (π/2)·4πr²ρ。"""
    return PI / 2 * rt3_radial(r, t, v, gam)


def rt3_once_radial(r, t, v=1.0, gam=1.0):
    u = r / (v * t)
    return gam * math.exp(-gam * t) * r * math.log((1 + u) / (1 - u)) / (v * v * t) if u < 1 else 0.0


def rt3_mass(t, v=1.0, gam=1.0, n=300):
    ph, r, jac = _phi_grid(t, v, n)
    f = np.array([rt3_radial(float(x), t, v, gam) for x in r]) * jac
    return float(_simps(f, ph[1] - ph[0])) + math.exp(-gam * t)


def rt3_fourier(q, t, v=1.0, gam=1.0, n=300):
    ph, r, jac = _phi_grid(t, v, n)
    f = np.array([rt3_radial(float(x), t, v, gam) for x in r]) * np.sinc(q * r / PI) * jac
    return float(_simps(f, ph[1] - ph[0])) + math.exp(-gam * t) * float(np.sinc(q * v * t / PI))


def rt3_laplace(q, lam, v=1.0, gam=1.0):
    """(5.91)–(5.93)：g/(1 − γg)，g = arctan(vq/(λ+γ))/(vq)。"""
    g = angle_avg_3d_closed(lam, q, v, gam)
    return g / (1 - gam * g)


def rt_mc(dim: int, t: float, M: int = 100000, seed: int = 5, v: float = 1.0, gam: float = 1.0):
    """固定種子的粒子模擬：一維＝電報過程（以速率 γ 反轉）；二三維＝以速率 γ 換成隨機新方向。回傳一維位置或二三維距離。"""
    rng = np.random.default_rng(seed)

    def new_dir(m):
        if dim == 1:
            return rng.choice([-1.0, 1.0], size=(m, 1))
        d = rng.normal(size=(m, dim))
        return d / np.linalg.norm(d, axis=1, keepdims=True)

    pos = np.zeros((M, dim))
    d = new_dir(M)
    tl = np.full(M, float(t))
    turned = np.zeros(M, bool)
    while True:
        act = tl > 0
        if not act.any():
            break
        idx = np.nonzero(act)[0]
        tau = rng.exponential(1 / gam, size=idx.size)
        step = np.minimum(tau, tl[idx])
        pos[idx] += d[idx] * (v * step)[:, None]
        tumble = tau < tl[idx]
        tl[idx] -= step
        k = idx[tumble]
        turned[k] = True
        if dim == 1:
            d[k] = -d[k]
        else:
            d[k] = new_dir(k.size)
    out = pos[:, 0] if dim == 1 else np.linalg.norm(pos, axis=1)
    out[~turned] = np.sign(out[~turned]) * v * t * (1 + 1e-12)     # 沒轉過的（尖峰）：放在 vt 外側一點點，區間 [a, b ≤ vt) 不算它
    return out


def mc_fraction(vals, a, b, absval=True):
    x = np.abs(vals) if absval else vals
    return float(np.mean((x >= a) & (x < b)))


def radial_fraction(dim: int, a: float, b: float, t: float, v=1.0, gam=1.0, n=400) -> float:
    """公式：a ≤ r < b 的比例（連續部分；一維是 |x|）。"""
    b = min(b, v * t)
    if b <= a:
        return 0.0
    ph = np.linspace(math.asin(a / (v * t)), math.asin(b / (v * t)), n + 1)
    if ph[-1] >= PI / 2 - 1e-7:
        ph[-1] = PI / 2 - 1e-7
    r = v * t * np.sin(ph)
    jac = v * t * np.cos(ph)
    if dim == 1:
        f = 2 * tele_cont(r, t, v, gam)
    elif dim == 2:
        f = rt2_radial(r, t, v, gam)
    else:
        f = np.array([rt3_radial(float(x), t, v, gam) for x in r])
    return float(_simps(f * jac, ph[1] - ph[0]))


# =========================================================================== §5.3 連續時間漫步

P0 = 1.0


def hazard(t, kind="poisson", p0=P0):
    return p0 if kind == "poisson" else t / (1 + t)


def survival(t, kind="poisson", p0=P0):
    return math.exp(-p0 * t) if kind == "poisson" else (1 + t) * math.exp(-t)


def survival_numeric(t, kind="poisson", p0=P0, n=400):
    ts = np.linspace(0, t, n + 1)
    qs = np.array([hazard(x, kind, p0) for x in ts])
    return math.exp(-float(_simps(qs, ts[1] - ts[0])))


def wait_density(t, kind="poisson", p0=P0):
    return p0 * math.exp(-p0 * t) if kind == "poisson" else t * math.exp(-t)


def mean_wait(kind="poisson", p0=P0):
    """B(0) = ∫a dt = 平均等待時間。"""
    return 1 / p0 if kind == "poisson" else 2.0


def AB(lam, kind="poisson", p0=P0):
    if kind == "poisson":
        return p0 / (p0 + lam), 1 / (p0 + lam)
    return 1 / (1 + lam) ** 2, (2 + lam) / (1 + lam) ** 2


def AB_numeric(lam, kind="poisson", p0=P0, T=60.0, n=6000):
    ts = np.linspace(0, T, n + 1)
    a = np.array([survival(x, kind, p0) for x in ts])
    q = np.array([hazard(x, kind, p0) for x in ts])
    e = np.exp(-lam * ts)
    return float(_simps(a * q * e, ts[1] - ts[0])), float(_simps(a * e, ts[1] - ts[0]))


def steps_prob(m: int, t: float, kind="poisson", p0=P0) -> float:
    pois = lambda k, mu: math.exp(k * math.log(mu) - mu - math.lgamma(k + 1)) if mu > 0 else float(k == 0)
    if kind == "poisson":
        return pois(m, p0 * t)
    return pois(2 * m, t) + pois(2 * m + 1, t)


def steps_prob_laplace(m: int, t: float, kind="poisson", p0=P0) -> float:
    def F(lam):
        A, B = AB(lam, kind, p0)
        return A ** m * B
    return talbot(F, t)


def steps_prob_conv(m: int, t: float, kind="poisson", p0=P0, n=3000) -> float:
    """(5.110)(5.111)：f_0 = a，f_{k+1}(τ) = ∫_0^τ a(s)q(s) f_k(τ−s) ds（梯形法數值卷積）。"""
    ts = np.linspace(0, t, n + 1)
    h = ts[1] - ts[0]
    w = np.array([wait_density(x, kind, p0) for x in ts])
    f = np.array([survival(x, kind, p0) for x in ts])
    for _ in range(m):
        c = np.convolve(w, f)[: n + 1] * h
        c -= 0.5 * h * (w[0] * f + w * f[0])           # 梯形法的端點修正
        f = c
    return float(f[-1])


def mean_steps(t, kind="poisson", p0=P0):
    return p0 * t if kind == "poisson" else t / 2 - 0.25 + math.exp(-2 * t) / 4


def ctrw_lattice(x: int, t: float, kind="poisson", p0=P0, mmax: int = 120) -> float:
    s = 0.0
    for m in range(abs(x), mmax + 1):
        if (m + x) % 2 == 0:
            s += steps_prob(m, t, kind, p0) * math.comb(m, (m + x) // 2) / 2 ** m
    return s


def ctrw_bessel(x: int, t: float, p0=P0) -> float:
    """Ex. 5.8：1/(λ + p0(1−cos q)) 的反轉 = e^{−p0t} I_|x|(p0 t)（walk_bc.ine_list）。"""
    return BC.ine_list(abs(x), p0 * t)[abs(x)]


def ctrw_var(t, kind="poisson", p0=P0):
    return mean_steps(t, kind, p0)


def ctrw_var_numeric(t, kind="poisson", p0=P0, X=40):
    return sum(x * x * ctrw_lattice(x, t, kind, p0) for x in range(-X, X + 1))


def ctrw_gauss(x, t, kind="poisson", p0=P0, C=0.5):
    """(5.131)–(5.133) ÷ (2π)（一維）：√(B0/4πCt) e^{−B0x²/4Ct}，B0 = B(0)。"""
    B0 = mean_wait(kind, p0)
    x = np.asarray(x, float)
    return np.sqrt(B0 / (4 * PI * C * t)) * np.exp(-B0 * x * x / (4 * C * t))


def ctrw_transform(lam, q, kind="poisson", p0=P0):
    A, B = AB(lam, kind, p0)
    return B / (1 - math.cos(q) * A)


def ex58(lam, q, p0=P0):
    return 1 / (lam + p0 * (1 - math.cos(q)))


def mortal_total(t, p0=1.0, p1=0.1):
    return math.exp(-p1 * t)


def mortal_series(t, p0=1.0, p1=0.1, mmax=200):
    """Σ_m (還活著且走了 m 步) = Σ_m (p0t)^m e^{−(p0+p1)t}/m!（Ex. 5.9，第二種）。"""
    s, term = 0.0, math.exp(-(p0 + p1) * t)
    for m in range(mmax):
        s += term
        term *= p0 * t / (m + 1)
    return s


def mortal_laplace(t, p0=1.0, p1=0.1):
    """Talbot 反轉 Σ A^m B = B/(1 − A)，A = p0/(p0+p1+λ)、B = 1/(p0+p1+λ)（第三種）。"""
    def F(lam):
        A, B = p0 / (p0 + p1 + lam), 1 / (p0 + p1 + lam)
        return B / (1 - A)
    return talbot(F, t)


# =========================================================================== 自我測試

def _persistent_brute(N: int, Q: Fraction) -> dict:
    """暴力列舉 2^N 條路線（第一步左右各半，之後同方向的機率 Q）。"""
    out = {}
    for code in range(2 ** N):
        steps = [1 if (code >> k) & 1 else -1 for k in range(N)]
        w = Fraction(1, 2) if N else Fraction(1)
        for a, b in zip(steps, steps[1:]):
            w *= Q if a == b else 1 - Q
        x = sum(steps)
        out[x] = out.get(x, 0) + w
    return out


def check() -> bool:
    f4 = lambda x: f"{x:.4f}"
    F = Fraction
    # ---------------------------------------------------------------- §5.1 偏向的漫步
    p = F(3, 5)
    for N in range(13):
        P = biased_pmf(N, p)
        assert P == biased_dp(N, p) and sum(P.values()) == 1
        assert moments(P) == (biased_mean(N, p), biased_var(N, p))
    P10 = biased_pmf(10, p)
    assert moments(P10) == (2, F(48, 5)) and biased_var(10, F(1, 2)) == 10
    assert [f4(float(P10[x])) for x in (0, 2, 4)] == ["0.2007", "0.2508", "0.2150"]
    assert [f4(2 * float(biased_gauss(x, 10, 0.6))) for x in (0, 2, 4)] == ["0.2091", "0.2575", "0.2091"]
    peaks = []
    for N, mu, var in ((10, 2, "9.6"), (20, 4, "19.2"), (40, 8, "38.4")):
        P = biased_pmf(N, 0.6)
        xm = max(P, key=P.get)
        assert xm == mu == round(biased_mean(N, 0.6)) and f"{biased_var(N, 0.6):.1f}" == var
        peaks.append((f4(P[xm]), f4(2 * float(biased_gauss(xm, N, 0.6)))))
    assert peaks == [("0.2508", "0.2575"), ("0.1797", "0.1821"), ("0.1279", "0.1288")]
    assert taylor_residual(0.6) < 1e-5 < 1e-3 < taylor_residual(0.6, book=True)          # (5.3) 的號
    for pp in (0.5, 0.6, 0.9):                                                              # Ex. 5.1
        assert abs(diffusion_D(pp) - biased_var(10, pp) / 20) < 1e-14 and diffusion_D(pp, kind="naive") == 0.5
    assert f"{diffusion_D(0.6):.2f}" == "0.48" and f"{math.sqrt(10):.3f}" == "3.162" and f"{math.sqrt(9.6):.3f}" == "3.098"
    assert f4(width_ratio(0.6)) == "0.9798" and f"{width_ratio(0.9):.1f}" == "0.6"
    sp, sd = shift_pmf(10, 0.2), shift_dp(10, 0.2)                                          # (5.9)(5.10)
    assert len(sp) == len(sd) and all(abs(a[0] - b[0]) < 1e-9 and abs(a[1] - b[1]) < 1e-15 for a, b in zip(sp, sd))
    m_ = sum(x * w for x, w in sp)
    assert abs(m_ - 2) < 1e-12 and abs(sum(x * x * w for x, w in sp) - m_ ** 2 - 10) < 1e-9
    for N, b, x in ((6, 1, 8), (7, 1, 5), (5, 2, 9), (6, -1, -4)):
        k = x - N * b
        exact = math.comb(N, (k + N) // 2) / 2 ** N if (k + N) % 2 == 0 and abs(k) <= N else 0.0
        assert abs(shift_fourier(N, b, x) - exact) < 1e-12
    # ---------------------------------------------------------------- §5.2 有慣性的漫步
    Q = F(3, 4)
    table = {0: {0: (F(1, 2), F(1, 2))},
             1: {-1: (F(1, 2), 0), 1: (0, F(1, 2))},
             2: {-2: (F(3, 8), 0), 0: (F(1, 8), F(1, 8)), 2: (0, F(3, 8))},
             3: {-3: (F(9, 32), 0), -1: (F(1, 8), F(3, 32)), 1: (F(3, 32), F(1, 8)), 3: (0, F(9, 32))},
             4: {-4: (F(27, 128), 0), -2: (F(15, 128), F(9, 128)), 0: (F(13, 128), F(13, 128)),
                 2: (F(9, 128), F(15, 128)), 4: (0, F(27, 128))}}
    for N, tab in table.items():
        assert persistent_dp(N, Q) == tab
    for Qv in (F(1, 4), F(1, 2), F(3, 4), F(1, 10)):
        for N in range(11):
            assert persistent_pmf(N, Qv) == _persistent_brute(N, Qv)
            assert persistent_var(N, Qv) == persistent_var_formula(N, Qv)
    assert [persistent_var(N, F(3, 4)) for N in (1, 2, 3, 4)] == [1, 3, F(11, 2), F(33, 4)]
    assert [persistent_var(N, F(1, 2)) for N in (1, 2, 3, 4)] == [1, 2, 3, 4]
    assert [persistent_var(N, F(1, 4)) for N in (1, 2, 3, 4)] == [1, 1, F(3, 2), F(7, 4)]
    rows = {}
    for Qv in (0.75, 0.25):
        for N in (10, 100, 1000):
            ex = persistent_var(N, Qv)
            assert abs(ex - persistent_var_formula(N, Qv)) < 1e-8 * N
            rows[(Qv, N)] = (ex, persistent_var_asym(N, Qv))
    assert [f"{rows[(0.75, N)][0]:.2f}" for N in (10, 100, 1000)] == ["26.00", "296.00", "2996.00"]
    assert [round(rows[(0.75, N)][1], 9) for N in (10, 100, 1000)] == [30, 300, 3000]
    assert [f"{rows[(0.25, N)][0]:.3f}" for N in (10, 100, 1000)] == ["3.777", "33.778", "333.778"]
    assert [f"{rows[(0.25, N)][1]:.3f}" for N in (10, 100, 1000)] == ["3.333", "33.333", "333.333"]
    assert persistent_var_asym(1000, F(3, 4)) - persistent_var_formula(1000, F(3, 4)) - 4 < F(1, 10 ** 100)   # 差 → 4
    assert abs(persistent_var_asym(1000, F(1, 4)) - persistent_var_formula(1000, F(1, 4)) + F(4, 9)) < F(1, 10 ** 100)  # 差 → −4/9
    assert persistent_var_formula(10, 1) == 100 and persistent_var(10, 1.0) == 100              # Q = 1：直線
    assert persistent_var_formula(10, 0) == 0 and persistent_var_formula(11, 0) == 1             # Q = 0：來回
    assert persistent_var(10, 0.0) == 0 and persistent_var(11, 0.0) == 1
    assert persistence_length(F(3, 4)) == 3 and persistence_length(F(1, 4)) == F(1, 3) and persistence_length(F(1, 2)) == 1
    xs, Lp, Rp = persistent_np(100, 0.75)
    i0 = int(np.nonzero(xs == 0)[0][0])
    assert f"{Lp[i0] + Rp[i0]:.5f}" == "0.04587" and f"{2 * float(gauss_persistent(0, 100, 0.75)):.5f}" == "0.04607"
    for q in (0.1, 0.3, 0.5, 1.2):                                                          # (5.17)(5.18)(5.37)
        l1, l2 = eigs(q, 0.75)
        num = np.linalg.eigvals(np.array(T_matrix(q, 0.75)))                                # 第二種：數值特徵值
        assert min(abs(num - l1)) < 1e-12 and min(abs(num - l2)) < 1e-12 and abs(l1 - l2) > 0.1
        assert eig_residual(l1, q, 0.75) < 1e-14 and eig_residual(l2, q, 0.75) < 1e-14
        assert abs(eig_residual(l1, q, 0.75, book=True) - 1.0) < 1e-12                       # 書上 (5.37)：差 2|1−2Q| = 1
    l1, l2 = eigs(0.1, 0.75)
    assert (f4(l1.real), f4(l2.real), f"{lambda1_small(0.1, 0.75):.3f}") == ("0.9848", "0.5077", "0.985")
    assert abs(l1 - lambda1_small(0.1, 0.75)) < 1e-3 and abs(eigs(0.01, 0.75)[0] - lambda1_small(0.01, 0.75)) < 1e-6
    z1, z2 = eigs(0.5, 0.75)
    assert (f4(z1.real), f4(z1.imag), z2 == z1.conjugate()) == ("0.6582", "0.2584", True)
    for N in (1, 5, 20):
        a, b, c = total_fourier(N, 0.3, 0.75), T_total(N, 0.3, 0.75), total_eig(N, 0.3, 0.75)
        assert abs(a - b) < 1e-13 and abs(a - c) < 1e-13
    assert f4(total_fourier(5, 0.3, 0.75).real) == "0.5699" and f4(total_fourier(20, 0.3, 0.75).real) == "0.0388"
    assert all(recursion2_residual(N, q, Qv) < 1e-14 for N in (0, 3, 7) for q in (0.2, 0.9) for Qv in (0.3, 0.75))
    # ---------------------------------------------------------------- §5.2.2 電報方程
    ie = BC.ine_list(1, 2.5)
    assert abs(_sp.i0e(2.5) - ie[0]) < 1e-13 and abs(_sp.i1e(2.5) - ie[1]) < 1e-13
    for t in (1, 3, 5, 8):
        assert abs(tele_mass(t) - 1) < 1e-9 and abs(tele_mass(t, v=2.0, gam=0.7) - 1) < 1e-9
        assert abs(tele_var(t) - tele_var_numeric(t)) < 1e-9
        assert abs(tele_var(t, 2.0, 0.7) - tele_var_numeric(t, 2.0, 0.7)) < 1e-8
    assert f4(tele_mass(3, v=2.0, book=True)) == "1.9502"                                 # 書上少 1/v：連續部分多乘了 v
    assert float(tele_cont(4.0, 3, v=2.0)) > 0.06                                          # p. 107「x < t」應為 x < vt
    assert float(tele_cont(3.01, 3)) == 0 and float(tele_cont(6.01, 3, v=2.0)) == 0         # 因果性
    assert [f4(float(tele_cont(0, t))) for t in (1, 3, 5, 8)] == ["0.3368", "0.2199", "0.1738", "0.1388"]
    front = [float(tele_cont(t * (1 - 1e-12), t)) for t in (1, 3, 5, 8)]
    assert [f"{x:.4g}" for x in front] == ["0.2759", "0.06223", "0.01179", "0.0008387"]
    assert all(abs(fx - math.exp(-t) / 2 * (1 + t / 2)) < 1e-9 for fx, t in zip(front, (1, 3, 5, 8)))
    assert [f"{tele_spike(t):.4g}" for t in (1, 3, 5, 8)] == ["0.1839", "0.02489", "0.003369", "0.0001677"]
    assert [f4(1 - math.exp(-t)) for t in (1, 3, 5, 8)] == ["0.6321", "0.9502", "0.9933", "0.9997"]
    assert [f4(tele_var(t)) for t in (1, 3, 5, 8)] == ["0.5677", "2.5012", "4.5000", "7.5000"]
    assert [f4(float(tele_gauss(0, t))) for t in (1, 3, 5, 8)] == ["0.3989", "0.2303", "0.1784", "0.1410"]
    assert abs(float(tele_cont(0, 60)) / float(tele_gauss(0, 60)) - 1) < 0.01                 # Ex. 5.5
    assert tele_D() == 0.5 and tele_lp() == 1.0 and abs(tele_var(40) - (2 * tele_D() * 40 - 0.5)) < 1e-12   # Ex. 5.6
    for q in (0.5, 1.0, 2.0):
        a = tele_fourier(q, 3)
        assert abs(a - tele_fourier_numeric(q, 3)) < 1e-9
        assert abs(a - talbot(lambda lam: tele_laplace(q, lam), 3)) < 1e-8
    h, x0, t0 = 1e-3, 0.7, 2.0                                                                 # Ex. 5.4：Σ_tt + 2γΣ_t = v²Σ_xx
    S = lambda x, t: float(tele_cont(x, t))
    res = ((S(x0, t0 + h) - 2 * S(x0, t0) + S(x0, t0 - h)) + h * (S(x0, t0 + h) - S(x0, t0 - h))
           - (S(x0 + h, t0) - 2 * S(x0, t0) + S(x0 - h, t0))) / h ** 2
    assert abs(res) < 1e-6
    lat = []
    for dt in (0.01, 0.005):                                                                   # Ex. 5.3：格子極限
        xl, Pl = tele_lattice(1.0, dt)
        lat.append(Pl[int(np.argmin(np.abs(xl)))] / (2 * dt))
    rich = 2 * lat[1] - lat[0]
    assert (f4(lat[0]), f4(lat[1]), f4(rich)) == ("0.3392", "0.3380", "0.3368") and abs(rich - float(tele_cont(0, 1))) < 1e-4
    xl, Pl = tele_lattice(3.0, 0.005)
    assert abs(Pl[-2] - tele_spike(3)) < 1e-4                                                  # 格子的尖峰 (1−γδt)^N/2
    x1 = rt_mc(1, 3.0, 100000, 5)
    mc1 = (mc_fraction(x1, 0, 1.5), float(np.mean(np.abs(x1) > 3)))
    assert (f4(mc1[0]), f4(radial_fraction(1, 0, 1.5, 3)), f4(mc1[1])) == ("0.6060", "0.6051", "0.0488")
    assert abs(mc1[0] - radial_fraction(1, 0, 1.5, 3)) < 0.005 and abs(mc1[1] - math.exp(-3)) < 0.003
    assert abs(radial_fraction(1, 0, 3, 3) - (1 - math.exp(-3))) < 1e-6
    # ---------------------------------------------------------------- §5.2.3–5.2.5 二維、三維
    for lam in (1.0, 0.3 + 0.2j, 2 - 1j):
        for q in (0.4, 1.0, 2.5):
            assert abs(angle_avg_2d(lam, q) - angle_avg_2d_closed(lam, q)) < 1e-12
            assert abs(angle_avg_3d(lam, q) - angle_avg_3d_closed(lam, q)) < 1e-9
    assert (f4(angle_avg_2d_closed(1, 1).real), f4(angle_avg_3d_closed(1, 1).real)) == ("0.4472", "0.4636")
    assert sm_residual() < 1e-12 and sm_residual(40, 0.2 - 0.5j, 3.0, 1.5, 0.8) < 1e-12           # Ex. 5.7
    for x in (0.2, 0.5, 0.9):                                                                  # p. 115：切割兩邊的 arctan
        on = cmath.atan(1 / complex(1e-13, x))
        assert abs(on - (PI / 2 - 0.5j * math.log((1 + x) / (1 - x)))) < 1e-9
        assert abs(on - (PI / 2 - 1j * math.log((1 + x) * (1 - x)))) > 0.1
        assert abs(cmath.atan(1 / complex(-1e-13, x)) - (-PI / 2 - 0.5j * math.log((1 + x) / (1 - x)))) < 1e-9
    for t in (1, 3):
        assert abs(rt2_mass(t) - 1) < 1e-5
        rr = np.linspace(0.1, 0.95 * t, 7)
        assert np.allclose(rt2_pole(rr, t) + rt2_cut(rr, t), rt2_cont(rr, t), rtol=1e-12)
        assert np.allclose(rt2_book_pole(rr, t) / (2 * PI * rt2_pole(rr, t)), math.sqrt(PI / 2))   # (5.84) 多了 √(π/2)
        assert np.allclose(rt2_book_cut(rr, t) / (2 * PI * rt2_cut(rr, t)), 2 * PI)                # (5.89) 多了 2π
    assert (f4(rt2_mass(1, book=True)), f4(rt2_mass(3, book=True))) == ("1.9340", "1.2658")
    assert abs(rt2_book_pole(1.0, 2.0, v=2.0) / (2 * PI * rt2_pole(1.0, 2.0, v=2.0)) - math.sqrt(PI / 2)) > 0.5  # cosh 少了 /v
    assert abs(hankel_pair(1.0) - math.exp(-0.5)) < 1e-6 and abs(hankel_pair(0.4) - math.exp(-0.08)) < 1e-6   # (5.88) 沒有 2π
    for q in (0.5, 1.0, 2.0):
        assert abs(rt2_hankel(q, 3) - talbot(lambda lam: rt2_laplace(q, lam), 3)) < 2e-6
    assert [f4(float(rt2_radial(r, 3))) for r in (1, 2, 2.5, 2.9)] == ["0.2978", "0.4167", "0.3941", "0.4052"]
    for t in (1, 3, 5, 8):
        assert abs(rt3_mass(t) - 1) < 1e-5
    for q in (0.5, 1.0, 2.0):
        assert abs(rt3_fourier(q, 3) - talbot(lambda lam: rt3_laplace(q, lam), 3)) < 1e-5
    for r in (3.5, 5.0):                                                                       # 因果性：r > vt 時兩項抵銷
        pole, cut = rt3_pole_book(r, 3), rt3_cut_book(r, 3)
        assert abs(pole + cut) < 1e-6 and abs(pole) > 5e-3
        assert abs(rt3_5106(r, 3, book=False) - cut) < 1e-6 and abs(rt3_5106(r, 3, book=True) - cut) > 5e-3   # (5.106) 的號
    assert (f4(rt3_pole_book(3.5, 3)), f4(rt3_cut_book(3.5, 3))) == ("0.0113", "-0.0113")
    for t in (1, 3):                                                                           # (5.103) 第三項＝只轉一次彎（直接算）
        for r in (0.3 * t, 0.8 * t):
            assert abs(4 * PI * r * r * rt3_once_book(r, t) / (2 * PI * PI) - rt3_once_radial(r, t)) < 1e-12
    ph, rr, jac = _phi_grid(3, 1.0, 400)
    once = float(_simps(np.array([rt3_once_radial(float(x), 3) for x in rr]) * jac, ph[1] - ph[0]))
    assert abs(once - 3 * math.exp(-3)) < 1e-5
    assert (f4(math.exp(-3)), f4(3 * math.exp(-3)), f4(1 - 4 * math.exp(-3))) == ("0.0498", "0.1494", "0.8009")
    assert [f4(float(rt3_radial(r, 3))) for r in (1, 1.5, 2, 2.5)] == ["0.2765", "0.4558", "0.5120", "0.4221"]
    fig56 = {1: ["0.1472", "0.5334", "1.0275", "1.3699"], 3: ["0.1729", "0.4558", "0.4828", "0.3594"],
             5: ["0.1869", "0.3630", "0.2169", "0.0926"], 8: ["0.1948", "0.2411", "0.0622", "0.0118"]}
    for t, row in fig56.items():
        assert [f4(float(rt3_radial(u * t, t))) for u in (0.25, 0.5, 0.75, 0.9)] == row
        assert abs(rt3_r2sigma_book(0.5 * t, t) - PI / 2 * rt3_radial(0.5 * t, t)) < 1e-15
    fr = {}
    for d in (2, 3):
        vals = rt_mc(d, 3.0, 100000, 5)
        fr[d] = (mc_fraction(vals, 0, 1.5), radial_fraction(d, 0, 1.5, 3), mc_fraction(vals, 2.5, 3), radial_fraction(d, 2.5, 3, 3),
                 float(np.mean(vals > 3)))
        assert abs(fr[d][0] - fr[d][1]) < 0.005 and abs(fr[d][2] - fr[d][3]) < 0.005 and abs(fr[d][4] - math.exp(-3)) < 0.003
    assert [f4(x) for x in fr[2][:2]] == ["0.3335", "0.3310"] and [f4(x) for x in fr[3][:2]] == ["0.2880", "0.2881"]
    assert (f4(fr[2][4]), f4(fr[3][4])) == ("0.0493", "0.0497")
    # ---------------------------------------------------------------- §5.3 連續時間漫步
    for kind in ("poisson", "gamma2"):
        for t in (0.5, 3.0):
            assert abs(survival(t, kind) - survival_numeric(t, kind)) < 1e-9
        for lam in (0.3, 0.7, 2.0):
            A, B = AB(lam, kind)
            An, Bn = AB_numeric(lam, kind)
            assert abs(A - An) < 1e-8 and abs(B - Bn) < 1e-8 and abs(A - (1 - lam * B)) < 1e-14       # (5.119)
        for m in range(5):
            a = steps_prob(m, 3.0, kind)
            assert abs(a - steps_prob_laplace(m, 3.0, kind)) < 1e-9 and abs(a - steps_prob_conv(m, 3.0, kind)) < 2e-6
        assert abs(sum(steps_prob(m, 3.0, kind) for m in range(60)) - 1) < 1e-14                    # (5.122)
        assert abs(talbot(lambda lam: AB(lam, kind)[1] / (1 - AB(lam, kind)[0]), 3.0) - 1) < 1e-9
        assert abs(mean_steps(3.0, kind) - sum(m * steps_prob(m, 3.0, kind) for m in range(60))) < 1e-12
        for t in (10.0, 40.0):
            assert abs(ctrw_var(t, kind) - ctrw_var_numeric(t, kind, X=60 if t < 20 else 90)) < 1e-8
    assert (f4(survival(3, "poisson")), f4(survival(3, "gamma2"))) == ("0.0498", "0.1991")
    assert (mean_wait("poisson"), mean_wait("gamma2")) == (AB(0, "poisson")[1], AB(0, "gamma2")[1]) == (1, 2)   # B(0)＝平均等待
    assert [f4(x) for x in AB(0.7, "poisson")] == ["0.5882", "0.5882"] and [f4(x) for x in AB(0.7, "gamma2")] == ["0.3460", "0.9343"]
    assert [f4(steps_prob(m, 3, "poisson")) for m in range(7)] == ["0.0498", "0.1494", "0.2240", "0.2240", "0.1680", "0.1008", "0.0504"]
    assert [f4(steps_prob(m, 3, "gamma2")) for m in range(5)] == ["0.1991", "0.4481", "0.2689", "0.0720", "0.0108"]
    assert (f4(mean_steps(3, "poisson")), f4(mean_steps(3, "gamma2"))) == ("3.0000", "1.2506")
    for x in range(-4, 5):                                                                     # Ex. 5.8 的格子
        assert abs(ctrw_lattice(x, 3.0) - ctrw_bessel(x, 3.0)) < 1e-13
    assert [f4(ctrw_bessel(x, 3.0)) for x in (0, 1, 2)] == ["0.2430", "0.1968", "0.1118"]
    for lam, q in ((0.4, 0.7), (1.3, 2.0)):
        assert abs(ctrw_transform(lam, q) - ex58(lam, q)) < 1e-14
    assert (f4(ctrw_lattice(0, 10, "poisson")), f4(float(ctrw_gauss(0, 10, "poisson")))) == ("0.1278", "0.1262")
    assert (f4(ctrw_lattice(0, 10, "gamma2")), f4(float(ctrw_gauss(0, 10, "gamma2")))) == ("0.1805", "0.1784")
    assert abs(ctrw_lattice(0, 200, "gamma2", mmax=400) / float(ctrw_gauss(0, 200, "gamma2")) - 1) < 0.01
    assert [f"{ctrw_var(t, 'poisson'):.2f}" for t in (10, 40)] == ["10.00", "40.00"]
    assert [f"{ctrw_var(t, 'gamma2'):.2f}" for t in (10, 40)] == ["4.75", "19.75"]
    assert (f"{math.sqrt(10):.3f}", f"{math.sqrt(40):.3f}") == ("3.162", "6.325")                 # p. 124：寬度 ∝ √t，不是 ∝ t
    assert abs(math.sqrt(ctrw_var(40, "poisson") / ctrw_var(10, "poisson")) - 2) < 1e-12
    xs_ = np.linspace(-40, 40, 4001)
    assert abs(_simps(ctrw_gauss(xs_, 10, "gamma2"), xs_[1] - xs_[0]) - 1) < 1e-10              # (5.133) ÷ (2π) 才是機率
    B0, C, t = 2.0, 0.5, 10.0
    assert abs(_simps(np.sqrt(PI * B0 / (C * t)) * np.exp(-B0 * xs_ ** 2 / (4 * C * t)), xs_[1] - xs_[0]) - 2 * PI) < 1e-9
    for t in (3.0, 10.0):                                                                      # Ex. 5.9
        assert abs(mortal_total(t) - mortal_series(t)) < 1e-13 and abs(mortal_total(t) - mortal_laplace(t)) < 1e-9
    assert f4(mortal_total(10)) == "0.3679" and f"{math.log(2) / 0.1:.3f}" == "6.931"
    return True


if __name__ == "__main__":
    import sys
    sys.stdout.reconfigure(encoding="utf-8")
    check()
    print("walk_var（第 5 章）：偏向漫步、有慣性的漫步、電報方程、二三維的跑與轉、連續時間漫步全部驗證通過")
