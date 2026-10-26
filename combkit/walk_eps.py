# 由 combkit/tools/sync_from_walk017.py 從 walk017 專案的 codes/utils/walk_eps.py 複製；請改 walk017 那份再同步。
"""
walk_eps.py —— 《Elements of the Random Walk》第 12 章 More on the renormalization group 的數學（全書最後一章）
（Python＋numpy／scipy（integrate、optimize、special），不 import manim；用到同資料夾的 walk_saw（c_N、count2）、walk_crit（ν 已知值、
z_c、Flory 差距）與 walk_rg（ln4/ln3、複本公式）——在 walk017 是 utils 套件、在 combkit 是 combkit 套件，都用相對 import；直接執行本檔時退回絕對 import）

約定（書上式號、筆誤與〔注〕見 ref/md/ch12_formulas.md）：
  一分量（Ising 型）ℋ = Σ_q ½(r + q²)|S(q)|² + (u/N)Σ S⁴δ（(12.3)）；布里淵區取半徑 Λ = e^{−ℓ} 的球；殼層的 q² 換成 1（r + 1）。
  遞迴方程用 d 吸收進 u 的寫法（本系列的 u ＝ 書上 (12.37) 的 u·d；書上 (12.66)(12.67) 本來就是這個寫法）：
    dr/dℓ = 2r + A u/(1+r)、du/dℓ = εu − B u²/(1+r)²，A = 4(n+2)、B = 4(n+8)（Ising n = 1：A = 12、B = 36），ε = 4 − d。
  分數用 Fraction（Fraction 不能直接 :.4f，要 float()）；回傳快取（lru_cache）的函數不要就地修改結果。

§12.1–12.2.4 前半：「怎樣正確地積掉自由度」(12.5)–(12.10)、殼層 (12.20)–(12.25)、重新標度 (12.29)–(12.36)          [W1201]
    TOY = (r, u) = (0.2, 0.01)；toy_H(s, S', r, u) = r s²/2 + u s⁴ + (r+1)S'²/2 + 6u s²S'²（s 在內區、S' 在殼層，只留 6u 頂點）
    toy_A(r) = √(2π/(r+1)) = 2.2882（(12.7) 的 A）；toy_Heff(s, r, u) = r s²/2 + u s⁴ + ½ln(1 + 12u s²/(r+1))（閉式）＝ toy_Heff_numeric（對 S' 數值積分）：
                                   ℋ'(1) = 0.1577
    toy_coeffs(r, u) = (r', u', c6) = (r + 12u/(1+r), u − 36u²/(1+r)², (12u/(1+r))³/6) = 0.3000、0.0075、0.000167（新的 s⁶ 項，書上丟掉）
                                   ＝ toy_coeffs_numeric（擬合）；正是 (12.28) 的 12、36（玩具裡沒有 dδℓ）
    toy_Z(r, u) → (二重積分, A·∫e^{−ℋ'}) = 9.0517 = 2.2882 × 3.9558（(12.8)(12.9) Z = A(K')Z'(K')）
    shell_fraction(d, δℓ) = 1 − e^{−dδℓ}、shell_linear = dδℓ、shell_fraction_lattice(d, δℓ, L=200)（數格點，d = 2、3；ball_count(d, R)）：
                                   SHELL_CASES (3, 0.1) → 0.2592／0.3000、(3, 0.01) → 0.0296／0.0300、(4, 0.01) → 0.0392／0.0400；格點 d = 3 差 < 2e−4
    shell_dq(ℓ, δℓ) = (e^{−ℓ} − e^{−ℓ−δℓ}, e^{−ℓ}δℓ)：ℓ = 0.5、δℓ = 0.01 → 0.006035、0.006065（Λ = e^{−0.5} = 0.6065）
    rescale_factors(d, δℓ, book=False) → (q² 項, r 項, u 項) 的乘數：d = 3、δℓ = 0.1 → 1.0000、1.2214（e^{2δℓ}）、1.1052（e^{(4−d)δℓ}）；
                                   book=True 照 (12.33) 的 q → q'e^{−dδℓ}：q² 項 0.6703（不再是 1）；(12.32)：e^{dδℓ} = 1.0305 vs 書上 1 + δℓ = 1.0100（d = 3、δℓ = 0.01）
§12.2.4 後半、§12.7：配對、組合因子、連結叢集 (12.12)–(12.15)、(12:S-1)–(12:S-7)，Ex. 12.1、12.9–12.11                   [W1202]
    每個四次片段（Fig. 12.6(b)）有 2 條殼層腳；配對後片段連成環：type 1 = 1 個片段自己接、type 2 = 2 個片段雙線、type 3 = 三角環…
    matchings(n)（全部 (2n−1)!! 個配對：n = 3、4、5、6 → 15、105、945、10395）；ring_sizes(m, n)；sizes_of(n1, n2, n3)
    type_census(n) = {環的大小: 配對數}：n = 3：(1,1,1) 1、(2,1) 6、(3) 8（共 15）；n = 4：1、12、32、12、48（共 105）；n = 5：(2,1,1,1) 20
    diagram_census(n) = {環的大小: 不同的片段圖數}：n = 3：1、3、1；n = 5：(2,1,1,1) 10（＝ 5!/(3!2!)；書上第二次寫 5!/(2!1!) = 60）
    pair_factor(n1, n2) = n!/(n1! n2! 2^{n2})（(12:S-4)）：(3,1) 10、(2,1) 6、(0,2) 3、(4,1) 15、(2,2) 45、(0,3) 15（＝列舉，n ≤ 6 全部）
    ring_aut(k) = f_k（列舉排列）：1、2、6、8、10（f₃ = 6 = 3!；k ≥ 4 是 2k 不是 k!）；general_factor({k: n_k}) = n!/Π(n_k! f_k^{n_k})（(12:S-6)）
    ring_weight(k) = 2^{k−1}(k−1)! = 1、2、8、48、384（k 個片段連成一個環的腳配對數）
    egf_check(x) → Σ(2n−1)!! xⁿ/n! = exp(Σ w_k x^k/k!) = (1−2x)^{−1/2}：x = 0.1 → 1.1180（三種都一樣；「指數化只留連通圖」的最簡形）
    LC = (a, σ²) = (0.1, 1)：X = −a s²（s 是 Gaussian，σ² = 1/(1+r)，r = 0）
    lc_average(a, σ², method)：'closed'＝'quad'＝'series'（Σ⟨Xⁿ⟩/n!，所有圖）＝'linked'（exp Σκ_k/k!，只有連通圖）= 0.9129；ln = −0.0912
    lc_order(a, σ², n) → (全部, 連通, 不連通)：n = 2 → 0.0150 = 0.0100 + 0.0050；lc_moments、lc_cumulants_rings（環）＝ lc_cumulants_moments（遞迴）：
                                   κ_k = −0.1000、0.0200、−0.0080、0.0048
    LC_MAT、lc_multi(A, method)    三個殼層變數：det(1+2A)^{−1/2} ＝ 三維 Gauss–Hermite ＝ 環的和 = 0.8087（Ex. 12.1 的數值核對）
    on_counts(n) → (2n+4, 4(n+8))（Ex. 12.5–12.7，einsum 縮併；n = 1..4 → (6,36)、(8,40)、(10,44)、(12,48)，線性外插 n = 0 → 4、32）
§12.2.4–12.4：一迴圈遞迴 (12.16)–(12.38)、(12.61)–(12.72)，Fig. 12.13，Ex. 12.2–12.4、12.8                            [W1203–W1205]
    AB(n) = (4(n+2), 4(n+8))：n = 0..3 → (8,32)、(12,36)、(16,40)、(20,44)
    beta(r, u, ε, n=1)（本系列的形）；beta_book(r, u, ε)（書上 (12.37)＋更正的 (12.38)：d 留在式子裡）
    fixed_point(ε, n=1)（精確，ε 可用 Fraction）r* = −Aε/(2B+Aε)、u* = 4Bε/(2B+Aε)²；Ising r* = −ε/(6+ε)、u* = ε/(6+ε)²：
                                   ε = 1 → −1/7 = −0.1429、1/49 = 0.0204；ε = ½ → −1/13 = −0.0769、2/169 = 0.0118；ε = 0.1 → −1/61 = −0.0164、10/3721 = 0.0027
                                   ＝ fixed_point_numeric（fsolve）；fixed_point_lowest(ε) = (−Aε/2B, ε/B)：−0.1667、0.0278；−0.0833、0.0139；−0.0167、0.0028
                                   O(n)、ε = 1：n = 0 (−1/9, 2/81)、2 (−1/6, 5/288)、3 (−5/27, 11/729)
    jacobian(r, u, ε, n)＝ jacobian_fd；ε = 1 的 WF 點精確 [[5/3, 14], [1/21, −1]]；eigen(ε, n, which='wf'|'gauss')：
                                   WF 1.8968、−1.2301（= (1±√22)/3）；Gaussian (2, ε)（u 方向相關 ⇒ ε > 0 時 Gaussian 點不穩定）
    lambda_r_closed(ε, n)＝ numpy 特徵值；nu_eigen(ε, n) = 1/λ_r ＝ nu_from_flow(ε, n)（(12.53)–(12.58)：量 ℓ*，第三種算法）
    nu_first(ε, n) = ½ + (n+2)ε/(4(n+8))（Fraction；(12.59)(12.68)）；nu_unexpanded(ε, n) = 1/(2 − Aε/B)（Ising = 3/(6−ε) = Flory 3/(d+2)）
                                   Ising（一階／特徵值／不展開）：ε = 0.1：0.5083／0.5076／0.5085；ε = ½：0.5417／0.5256／0.5455；ε = 1：0.5833／0.5272／0.6000
                                   （特徵值與一階差 O(ε²)：一迴圈方程只準到 ε 的一次）
                                   ε = 1、n = 0..3：一階 0.5625、0.5833、0.6000、0.6136（9/16、7/12、3/5、27/44）；特徵值 0.5205、0.5272、0.5325、0.5368；
                                   不展開 0.5714、0.6000、0.6250、0.6471；NU_3D_CITED（只引用）0.588、0.630、0.672、0.711
    NU_EPS = (1/10, 1/2, 1)
    flow_line(r0, u0, ε=1, n=1, lmax=10, npts=40, box)（(r, u) 平面一條流線，(k, 2) 陣列）；FLOW_U0 = (0.002, 0.04)；
                                   flow_starts(ε, n)：u0 各取 r_c(u0) + (−0.02, 0, +0.02) 共 6 條（中間那條是臨界線，流進 WF 點；其餘跑出 |r| = 0.45）
    flow_end(r0, u0, ε, n) = +1／−1／0；r_critical(u0, ε, n)（二分法，快取）vs r_critical_linear = −A u0/(2−ε)（Ex. 12.4）：
                                   u0 = 10⁻⁴、ε = 0.1：−0.000631 vs −0.000632 vs 一階 −6u0 = −0.000600（書上 r = t − 6ud）；u0 = 10⁻³、ε = 1：−0.0111 vs −0.0120
    u_analytic(ℓ, u0, ε, B=36)（(12.41)）＝ u_numeric（solve_ivp）；U_CASES：u0 = 0.005：ℓ = 1、2、5 → 0.0104、0.0172、0.0270；
                                   u0 = 0.06 → 0.0346、0.0300、0.0279 → u* = 1/36 = 0.0278；u_from_1240(u0, u1, ε)：(12.40) 的積分 = ℓ（Ex. 12.2）
                                   Ex. 12.3：u = 0 的斜率 +ε（不穩定）、u* 的斜率 −ε；u0 = 10⁻⁶：u(2)/u0 = 7.3874 ≈ e² = 7.3891
    int_1248(u0, u1, ε, d, book) ＝ int_1249（−⅓ln(ε − 36ud)|）：u0 = 0.001 → 0.005、ε = 1、d = 3：0.2207；照書上 (12.48) 81.7804
    t_closed(ℓ, t0, u0, ε, n)（(12.51)，指數 A/B）＝ t_numeric：ℓ = 3、t0 = 10⁻⁴、u0 = 0.01 → 0.0203；大 ℓ 斜率 2 − ε/3
    flory_nu(d) = 3/(d+2)；flory_series(k) = ½(1/6)^k：½、1/12、1/72、1/432（3/(6−ε) 的展開；一階 ε = 1 → 7/12 = 0.5833，係數 1/12 ≠ n = 0 的 1/16）
    saw_size(N, ν) = N^ν：N = 1000 → 31.62（½）、48.70（9/16）、58.08（0.588）
§12.5–12.6：圖解法 (12.81)–(12.89)                                                                                  [W1206]
    g_numeric(q0, g0, ε, B=36, η=0)（(12.86) 從 q0 = 1 積下來）＝ u_analytic(−ln q0)：q0 = e^{−5} = 0.0067 → 0.0270；g* = ε/B = 0.0278
    pc_exponent(d, book=False) = d − 3（u²r^{d−3}）；book=True 照書上 2d − 3；pc_rhs(d, book) = d − 3（r·(ur^{−(4−d)/2})²）、書上 r^{−1}·(…)² → d − 5
    watermelon(d, r)（Fig. 12.27 外動量 0，實空間 ∫G³，d < 3）；pc_numeric(d) = ln(I(4)/I(1))/ln 4：d = 2 → −1.0000、d = 2.5 → −0.5000
全書回顧（W1206 最後一步）
    BOOK_RECAP = [(章, 主角, 關鍵結果的 Tex)] 12 列；recap_numbers() ＝ RECAP_EXPECT：
      1：252/1024 = 0.2461 vs 2/√(20π) = 0.2523；2：C(N/2)² = 1、4、36、400、4900、63504（＝ walk_saw.count2）；3：sc 回家 0.3405（Watson：Gamma 閉式＝積分）；
      4：x0 = 2、L = 5 → 3/5（線性方程）；5：Q = ½、N = 1000：2996（閉式＝雙重和）；6：10/19 = 0.5263；7：c₁₄ = 2374444；8：μ = 1/z_c = 2.638；
      9：Flory 差 0.0204（2%）；10：ν' = 1/2；11：1.2619（＝ walk_rg.replica_dim(4, 3)）；12：0.5833、0.5625
    BOOK_SCENES = {章: 幕數}（12、12、8、8、7、7、6、4、6、5、6、6，共 87）；scene_range(ch) → 'W101–W112'、…、'W1201–W1206'
summary()                         場景卡與講義用的全部數字（dict）

check() 驗證全部（每個結果至少兩種獨立算法、筆誤的數值、plan_ch12.md 場景卡與講義的每個數字）；直接執行本檔會跑 check()（約 1–2 秒）。
`python walk_eps.py --summary` 印出 summary()。最慢的呼叫：r_critical 每個約 0.2 秒（快取）、toy_Z 約 0.1 秒；其餘都很快。
"""
import math
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from itertools import combinations, permutations

import numpy as np
from scipy import integrate, optimize, special

try:                                    # walk017：utils 套件；combkit：combkit 套件
    from . import walk_saw as SW
    from . import walk_crit as WC
    from . import walk_rg as RG
except ImportError:                     # 直接執行本檔
    import walk_saw as SW
    import walk_crit as WC
    import walk_rg as RG

F = Fraction


def _f4(x):
    return f"{float(x):.4f}"


def _close(a, b, tol=1e-9):
    a, b = float(a), float(b)
    return abs(a - b) <= tol * max(1.0, abs(b))


# ======================================================================
# §12.1–12.2.4 前半：兩變數玩具、殼層體積、重新標度  [W1201]
# ======================================================================
TOY = (0.2, 0.01)                       # (r, u)：卡上的例子


def toy_H(s, Sp, r, u):
    """一個「內區」變數 s、一個「殼層」變數 S'：ℋ = r s²/2 + u s⁴ + (r+1) S'²/2 + 6u s² S'²（只留 (12.11) 的 6u 頂點）。"""
    return 0.5 * r * s * s + u * s ** 4 + 0.5 * (r + 1) * Sp * Sp + 6 * u * s * s * Sp * Sp


def toy_A(r):
    """(12.7) 的 A：S' 的 Gaussian 正規化 √(2π/(r+1))。"""
    return math.sqrt(2 * math.pi / (r + 1))


def toy_Heff(s, r, u):
    """閉式：ℋ'(s) = −ln[∫dS' e^{−ℋ} / A] = r s²/2 + u s⁴ + ½ ln(1 + 12u s²/(r+1))。"""
    return 0.5 * r * s * s + u * s ** 4 + 0.5 * math.log1p(12 * u * s * s / (r + 1))


def toy_Heff_numeric(s, r, u):
    """同上，對 S' 做數值積分（第二種算法）。"""
    val, _ = integrate.quad(lambda x: math.exp(-toy_H(s, x, r, u)), -np.inf, np.inf, epsabs=1e-14, epsrel=1e-13)
    return -math.log(val / toy_A(r))


def toy_coeffs(r, u):
    """ℋ' 的 Taylor 係數 → (r', u', c6)：r' = r + 12u/(1+r)、u' = u − 36u²/(1+r)²、s⁶ 的新項 c6 = (12u/(1+r))³/6。"""
    x = 12 * u / (1 + r)
    return r + x, u - 36 * u * u / (1 + r) ** 2, x ** 3 / 6


def toy_coeffs_numeric(r, u):
    """數值 ℋ'(s) 對 s² 做多項式擬合（第二種算法）→ (r', u')。"""
    xs = np.linspace(0.0, 0.4, 41)
    ys = np.array([toy_Heff_numeric(math.sqrt(x), r, u) for x in xs])
    c = np.polynomial.polynomial.polyfit(xs, ys, 7)
    return 2 * c[1], c[2]


def toy_Z(r, u):
    """(12.8)(12.9)：Z = ∫∫ e^{−ℋ} ds dS'（二重積分）與 A·∫ e^{−ℋ'(s)} ds（先積 S'）。"""
    z2, _ = integrate.dblquad(lambda Sp, s: math.exp(-toy_H(s, Sp, r, u)), -12, 12, -12, 12, epsabs=1e-12, epsrel=1e-11)
    z1, _ = integrate.quad(lambda s: math.exp(-toy_Heff(s, r, u)), -np.inf, np.inf, epsabs=1e-13, epsrel=1e-12)
    return z2, toy_A(r) * z1


def shell_fraction(d, dl):
    """殼層 e^{−ℓ−δℓ} < |q| < e^{−ℓ} 佔布里淵區的比例 1 − e^{−dδℓ}（(12.20)–(12.25)）。"""
    return 1 - math.exp(-d * dl)


def shell_linear(d, dl):
    return d * dl


def ball_count(d, R):
    """|x| ≤ R 的整數格點數（d = 2、3；d = 3 對 z 用解析計數）。"""
    m = int(math.floor(R))
    xs = np.arange(-m, m + 1)
    if d == 2:
        rest = R * R - xs * xs
        return int(np.sum(2 * np.floor(np.sqrt(np.maximum(rest, 0))) + 1))
    X, Y = np.meshgrid(xs, xs, indexing="ij")
    rest = R * R - X * X - Y * Y
    ok = rest >= 0
    return int(np.sum(2 * np.floor(np.sqrt(rest[ok])) + 1))


def shell_fraction_lattice(d, dl, L=200):
    """布里淵區裡「正好 N 個 q」：用半徑 L 的球內格點數數殼層的比例（第二種算法）。"""
    return 1 - ball_count(d, L * math.exp(-dl)) / ball_count(d, L)


SHELL_CASES = ((3, 0.1), (3, 0.01), (4, 0.01))   # 卡上 (d, δℓ)：0.2592／0.3000、0.0296／0.0300、0.0392／0.0400


def shell_dq(ell, dl):
    """(12.24)：dq = e^{−ℓ} − e^{−ℓ−δℓ} ≈ e^{−ℓ}δℓ → (精確, 一階)。"""
    return math.exp(-ell) - math.exp(-ell - dl), math.exp(-ell) * dl


def rescale_factors(d, dl, book=False):
    """(12.32)–(12.34) 之後各項的乘數：q² 項、r 項、u 項。正確：q → q'e^{−δℓ}、S → S''e^{δℓ}、N → N'e^{dδℓ}
    ⇒ (1, e^{2δℓ}, e^{(4−d)δℓ})；book=True 照書上 (12.33) q → q'e^{−dδℓ}：q² 項變成 e^{−2(d−1)δℓ}。"""
    qf = math.exp(-d * dl) if book else math.exp(-dl)
    sf = math.exp(dl)
    return qf * qf * sf * sf, sf * sf, sf ** 4 * math.exp(-d * dl)


# ======================================================================
# §12.2.4 後半、§12.7：配對、組合因子、連結叢集  [W1202]
# ======================================================================
def _matchings(legs):
    if not legs:
        yield ()
        return
    a = legs[0]
    for i in range(1, len(legs)):
        rest = legs[1:i] + legs[i + 1:]
        for m in _matchings(rest):
            yield ((a, legs[i]),) + m


@lru_cache(maxsize=None)
def matchings(n):
    """n 個四次片段（每個 2 條殼層腳：片段 i 的腳是 2i、2i+1）的全部配對：(2n−1)!! 個。"""
    return tuple(_matchings(tuple(range(2 * n))))


def ring_sizes(m, n):
    """一個配對把片段連成的環（每個片段剛好 2 條腳 ⇒ 每個連通塊都是環）的大小，由大到小；
    1 = type 1（自己的兩腳相接）、2 = type 2（兩個片段雙線相連）、3 = type 3（三角環）…"""
    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    for a, b in m:
        ra, rb = find(a // 2), find(b // 2)
        if ra != rb:
            parent[ra] = rb
    return tuple(sorted(Counter(find(i) for i in range(n)).values(), reverse=True))


@lru_cache(maxsize=None)
def type_census(n):
    """{環的大小 tuple: 配對數}（腳的層次）；總和 = (2n−1)!!。"""
    return dict(Counter(ring_sizes(m, n) for m in matchings(n)))


def _graph_of(m):
    return tuple(sorted(tuple(sorted((a // 2, b // 2))) for a, b in m))


@lru_cache(maxsize=None)
def diagram_census(n):
    """{環的大小: 不同的「片段圖」數}（只看哪些片段相連、不管用哪一條腳；第二種數法）。"""
    seen = {}
    for m in matchings(n):
        seen.setdefault(ring_sizes(m, n), set()).add(_graph_of(m))
    return {k: len(v) for k, v in seen.items()}


def sizes_of(n1, n2=0, n3=0):
    return tuple([3] * n3 + [2] * n2 + [1] * n1)


def pair_factor(n1, n2):
    """(12:S-4)：n!/(n1! n2! (2!)^{n2})，n = n1 + 2n2。"""
    n = n1 + 2 * n2
    return math.factorial(n) // (math.factorial(n1) * math.factorial(n2) * 2 ** n2)


def ring_aut(k):
    """f_k：k 個頂點排成環的「換頂點但配對不變」的方法數（列舉 k! 個排列）：1、2、6、8、10（k ≥ 3 是 2k）。"""
    edges = sorted(tuple(sorted((i, (i + 1) % k))) for i in range(k)) if k > 1 else [(0, 0)]
    if k == 2:
        edges = [(0, 1), (0, 1)]
    cnt = 0
    for p in permutations(range(k)):
        e2 = sorted(tuple(sorted((p[a], p[b]))) for a, b in edges)
        cnt += e2 == edges
    return cnt


def general_factor(counts):
    """(12:S-6)：n!/Π(n_k! f_k^{n_k})，counts = {k: n_k}。"""
    n = sum(k * c for k, c in counts.items())
    den = 1
    for k, c in counts.items():
        den *= math.factorial(c) * ring_aut(k) ** c
    return F(math.factorial(n), den)


def ring_weight(k):
    """一組 k 個（標了號的）片段連成一個環的腳配對數 w_k = 2^{k−1}(k−1)!：1、2、8、48。"""
    return 2 ** (k - 1) * math.factorial(k - 1)


def egf_check(x, nmax=40):
    """連結叢集的最簡形：Σ_n (2n−1)!! xⁿ/n! = exp(Σ_k w_k x^k/k!) = (1−2x)^{−1/2} → (左, 右, 閉式)。"""
    left = sum(special.factorial2(2 * n - 1, exact=True) * x ** n / math.factorial(n) for n in range(nmax))
    right = math.exp(sum(ring_weight(k) * x ** k / math.factorial(k) for k in range(1, nmax)))
    return left, right, (1 - 2 * x) ** -0.5


LC = (0.1, 1.0)                         # X = −a s²，s 是變異數 σ² = 1/(1+r) 的 Gaussian（r = 0）


def lc_moments(a, sig2, nmax):
    """⟨Xⁿ⟩ = (−a)ⁿ (2n−1)!! σ^{2n}（全部的圖，含不連通的），n = 0..nmax。"""
    return [(-a) ** n * special.factorial2(2 * n - 1, exact=True) * sig2 ** n if n else 1.0 for n in range(nmax + 1)]


def lc_cumulants_rings(a, sig2, kmax):
    """連通的圖＝k 個片段連成一個環：κ_k = (−aσ²)^k · 2^{k−1}(k−1)!（k = 1..kmax）。"""
    return [(-a * sig2) ** k * ring_weight(k) for k in range(1, kmax + 1)]


def lc_cumulants_moments(a, sig2, kmax):
    """由動差遞迴求累積量（第二種算法）：κ_n = m_n − Σ_{j<n} C(n−1, j−1) κ_j m_{n−j}。"""
    m = lc_moments(a, sig2, kmax)
    k = [0.0] * (kmax + 1)
    for n in range(1, kmax + 1):
        k[n] = m[n] - sum(math.comb(n - 1, j - 1) * k[j] * m[n - j] for j in range(1, n))
    return k[1:]


def lc_average(a, sig2, method="closed"):
    """⟨e^X⟩：'closed' (1+2aσ²)^{−1/2}、'quad'（Gauss–Hermite）、'series'（Σ⟨Xⁿ⟩/n!）、'linked'（exp Σ κ_k/k!）。"""
    if method == "closed":
        return (1 + 2 * a * sig2) ** -0.5
    if method == "quad":
        x, w = np.polynomial.hermite.hermgauss(60)
        return float(np.sum(w * np.exp(-a * 2 * sig2 * x * x)) / math.sqrt(math.pi))
    if method == "series":
        return sum(mn / math.factorial(n) for n, mn in enumerate(lc_moments(a, sig2, 40)))
    return math.exp(sum(kk / math.factorial(k + 1) for k, kk in enumerate(lc_cumulants_rings(a, sig2, 40))))


def lc_order(a, sig2, n):
    """第 n 階：(全部 ⟨Xⁿ⟩/n!, 連通 κ_n/n!, 不連通 = 差)。"""
    allp = lc_moments(a, sig2, n)[n] / math.factorial(n)
    conn = lc_cumulants_rings(a, sig2, n)[n - 1] / math.factorial(n)
    return allp, conn, allp - conn


LC_MAT = np.array([[0.10, 0.02, 0.0], [0.02, 0.08, 0.01], [0.0, 0.01, 0.05]])


def lc_multi(A=LC_MAT, method="det"):
    """三個殼層變數（σ² = 1）：⟨e^{−sᵀAs}⟩：'det' det(1+2A)^{−1/2}、'quad'（三維 Gauss–Hermite）、'linked'（環的和 Σ(−1)^k 2^{k−1} tr A^k / k）。"""
    A = np.asarray(A, float)
    if method == "det":
        return float(np.linalg.det(np.eye(len(A)) + 2 * A) ** -0.5)
    if method == "quad":
        x, w = np.polynomial.hermite.hermgauss(24)
        s = math.sqrt(2) * x
        S = np.stack(np.meshgrid(s, s, s, indexing="ij"), -1)
        W = w[:, None, None] * w[None, :, None] * w[None, None, :]
        q = np.einsum("...i,ij,...j->...", S, A, S)
        return float(np.sum(W * np.exp(-q)) / math.pi ** 1.5)
    tot, P = 0.0, np.eye(len(A))
    for k in range(1, 60):
        P = P @ A
        tot += (-1) ** k * 2 ** (k - 1) * np.trace(P) / k
    return math.exp(tot)


def on_counts(n):
    """O(n) 頂點 (S·S)(S·S) 的縮併（Ex. 12.5–12.7，numpy einsum）→ (2n+4, 4(n+8))：
    一個頂點的兩條腳在殼層（6 種選法）縮成 δ 的係數；兩個頂點各 2 條腳互接（6×6×2 = 72 種）的四腳張量對稱化後是 (S·S)² 的幾倍、再乘 ½（二階的 1/2!）。"""
    dlt = np.eye(n)
    V = np.einsum("ab,cd->abcd", dlt, dlt)
    L1 = "abcd"
    two = np.zeros((n, n))
    for i, j in combinations(range(4), 2):
        s = list(L1)
        s[j] = s[i]
        out = "".join(L1[k] for k in range(4) if k not in (i, j))
        two += np.einsum("".join(s) + "->" + out, V)
    c2 = two[0, 0]
    assert np.allclose(two, c2 * dlt)
    L2 = "efgh"
    T = np.zeros((n,) * 4)
    for i, j in combinations(range(4), 2):
        for k, l in combinations(range(4), 2):
            for kk, ll in ((k, l), (l, k)):
                s1, s2 = list(L1), list(L2)
                s2[kk], s2[ll] = s1[i], s1[j]
                out = "".join(L1[m] for m in range(4) if m not in (i, j)) + "".join(L2[m] for m in range(4) if m not in (k, l))
                T += np.einsum("".join(s1) + "," + "".join(s2) + "->" + out, V, V)
    Ts = sum(np.transpose(T, p) for p in permutations(range(4))) / 24
    Sd = sum(np.transpose(V, p) for p in permutations(range(4))) / 24
    c4 = np.sum(Ts * Sd) / np.sum(Sd * Sd)
    assert np.allclose(Ts, c4 * Sd)
    return float(c2), float(c4 / 2)


# ======================================================================
# §12.2.4–12.4：遞迴方程、不動點、線性化、ν  [W1203–W1205]
# ======================================================================
def AB(n=1):
    """(12.66)(12.67) 的係數 A = 4(n+2)、B = 4(n+8)；Ising（n = 1）12、36。"""
    return 4 * (n + 2), 4 * (n + 8)


def beta(r, u, eps, n=1):
    """dr/dℓ = 2r + A u/(1+r)、du/dℓ = εu − B u²/(1+r)²（本系列：書上 (12.37) 的 d 吸收進 u）。"""
    A, B = AB(n)
    return 2 * r + A * u / (1 + r), eps * u - B * u * u / (1 + r) ** 2


def beta_book(r, u, eps):
    """書上 (12.37)＋更正後的 (12.38)：d = 4 − ε 留在式子裡：2r + 12ud/(1+r)、εu − 36d u²/(1+r)²。"""
    d = 4 - eps
    return 2 * r + 12 * u * d / (1 + r), eps * u - 36 * d * u * u / (1 + r) ** 2


def fixed_point(eps, n=1):
    """Wilson–Fisher 不動點（精確，eps 可用 Fraction）：r* = −Aε/(2B+Aε)、u* = 4Bε/(2B+Aε)²；
    Ising：r* = −ε/(6+ε)、u* = ε/(6+ε)²（ε = 1：−1/7、1/49）。"""
    A, B = AB(n)
    den = 2 * B + A * eps
    return -A * eps / den, 4 * B * eps / den ** 2


def fixed_point_numeric(eps, n=1):
    A, B = AB(n)
    sol = optimize.fsolve(lambda v: beta(v[0], v[1], eps, n), [-A * eps / (2 * B), eps / B], xtol=1e-14)
    return float(sol[0]), float(sol[1])


def fixed_point_lowest(eps, n=1):
    """一階：r* ≈ −Aε/(2B)、u* ≈ ε/B（Ising −ε/6、ε/36；書上 u* → ε/36d）。"""
    A, B = AB(n)
    return -F(A, 2 * B) * eps, F(1, B) * eps


def jacobian(r, u, eps, n=1):
    A, B = AB(n)
    return np.array([[2 - A * u / (1 + r) ** 2, A / (1 + r)],
                     [2 * B * u * u / (1 + r) ** 3, eps - 2 * B * u / (1 + r) ** 2]], float)


def jacobian_fd(r, u, eps, n=1, h=1e-6):
    J = np.zeros((2, 2))
    for j, (dr, du) in enumerate(((h, 0), (0, h))):
        p = np.array(beta(r + dr, u + du, eps, n))
        m = np.array(beta(r - dr, u - du, eps, n))
        J[:, j] = (p - m) / (2 * h)
    return J


def eigen(eps, n=1, which="wf"):
    """線性化矩陣的特徵值（由大到小）：'gauss' 在 (0, 0) 是 (2, ε)；'wf' 在 Wilson–Fisher 點。"""
    r, u = (0.0, 0.0) if which == "gauss" else fixed_point(float(eps), n)
    return tuple(sorted(np.linalg.eigvals(jacobian(r, u, float(eps), n)).real, reverse=True))


def lambda_r_closed(eps, n=1):
    """WF 點的相關特徵值（閉式）：矩陣 [[2 − Aε/B, ·], [·, −ε]]、非對角乘積 2Aε²/B。Ising ε = 1：(1+√22)/3 = 1.8968。"""
    A, B = AB(n)
    a = 2 - A * eps / B
    return ((a - eps) + math.sqrt((a + eps) ** 2 + 8 * A * eps * eps / B)) / 2


def nu_first(eps, n=1):
    """(12.59)(12.68)：ν = ½ + (n+2)ε/(4(n+8))（Fraction）。"""
    return F(1, 2) + F(n + 2, 4 * (n + 8)) * eps


def nu_unexpanded(eps, n=1):
    """書上 (12.52)–(12.54) 不展開：ν = 1/(2 − Aε/B)（Ising 1/(2 − ε/3) = 3/(6−ε)，剛好是 Flory 3/(d+2)）。"""
    A, B = AB(n)
    return 1 / (2 - F(A, B) * eps)


def nu_eigen(eps, n=1):
    """完整的一迴圈線性化：ν = 1/λ_r（數值特徵值）。"""
    return 1 / eigen(eps, n)[0]


def nu_from_flow(eps, n=1, deltas=(1e-5, 1e-7), Dist=0.05):
    """(12.53)–(12.58) 的數值版：從 WF 點沿相關方向偏 δ 出發，量到 |r − r*| = Dist 的 ℓ*；ν = −Δℓ*/Δln δ（第三種算法）。"""
    r0, u0 = fixed_point(float(eps), n)
    w, v = np.linalg.eig(jacobian(r0, u0, float(eps), n))
    vec = v[:, int(np.argmax(w.real))].real
    vec = vec / vec[0]
    ells = []
    for dl in deltas:
        ev = lambda l, y: abs(y[0] - r0) - Dist
        ev.terminal = True
        sol = integrate.solve_ivp(lambda l, y: beta(y[0], y[1], eps, n), (0, 200), [r0 + dl * vec[0], u0 + dl * vec[1]],
                                  events=ev, rtol=1e-11, atol=1e-14)
        ells.append(sol.t_events[0][0])
    return -(ells[1] - ells[0]) / (math.log(deltas[1]) - math.log(deltas[0]))


NU_3D_CITED = {0: 0.588, 1: 0.630, 2: 0.672, 3: 0.711}     # 三維的已知 ν（只引用）：自迴避、Ising、XY、Heisenberg
NU_EPS = (F(1, 10), F(1, 2), F(1))                          # 卡上的 ε


def flow_line(r0, u0, eps=1.0, n=1, lmax=10.0, npts=40, box=(-0.45, 0.45, -0.001, 0.06)):
    """(r, u) 平面的一條流線（solve_ivp），跑出 box = (rmin, rmax, umin, umax) 就停 → (k, 2) 陣列（≤ npts 點）。"""
    def ev(l, y):
        return min(y[0] - box[0], box[1] - y[0], y[1] - box[2], box[3] - y[1])
    ev.terminal = True
    sol = integrate.solve_ivp(lambda l, y: beta(y[0], y[1], eps, n), (0, lmax), [r0, u0], events=ev,
                              dense_output=True, rtol=1e-9, atol=1e-12)
    ts = np.linspace(0, sol.t[-1], npts)
    return sol.sol(ts).T


def flow_end(r0, u0, eps=1.0, n=1, lmax=40.0):
    """往上（+1）或往下（−1）跑出去，或停在 WF 附近（0）。"""
    def up(l, y):
        return y[0] - 0.5
    def dn(l, y):
        return y[0] + 0.5
    up.terminal = dn.terminal = True
    sol = integrate.solve_ivp(lambda l, y: beta(y[0], y[1], eps, n), (0, lmax), [r0, u0], events=(up, dn), rtol=1e-9, atol=1e-13)
    if len(sol.t_events[0]):
        return 1
    if len(sol.t_events[1]):
        return -1
    return 0


@lru_cache(maxsize=None)
def r_critical(u0, eps=1.0, n=1):
    """Ex. 12.4：臨界的 r(0)（二分法：再大就流向 +∞、再小就流向 −）。"""
    lin = r_critical_linear(u0, eps, n)
    lo, hi = max(-0.45, 1.5 * lin - 1e-3), min(0.05, 0.5 * lin + 1e-3) if u0 < 0.01 else 0.05
    assert flow_end(lo, u0, eps, n) < 0 < flow_end(hi, u0, eps, n)
    for _ in range(36):
        mid = (lo + hi) / 2
        if flow_end(mid, u0, eps, n) > 0:
            hi = mid
        else:
            lo = mid
    return (lo + hi) / 2


def r_critical_linear(u0, eps=1.0, n=1):
    """Gaussian 點的 ε 特徵向量 (A, ε − 2)：r_c ≈ −A u0/(2 − ε)；ε → 0 為 −(A/2)u0（Ising −6u0，書上 r = t − 6ud）。"""
    A, B = AB(n)
    return -A * u0 / (2 - eps)


FLOW_U0 = (0.002, 0.04)                 # 卡上流線的起點 u0（ε = 1、Ising）；r0 = r_c(u0) + (−0.02, 0, 0.02)


def flow_starts(eps=1.0, n=1):
    return [(r_critical(u0, eps, n) + dlt, u0) for u0 in FLOW_U0 for dlt in (-0.02, 0.0, 0.02)]


def u_analytic(ell, u0, eps, B=36):
    """(12.41)（d 吸收進 u）：u(ℓ) = εu0 e^{εℓ}/(ε + B u0(e^{εℓ} − 1))（r = 0）。"""
    e = math.exp(eps * ell)
    return eps * u0 * e / (eps + B * u0 * (e - 1))


def u_numeric(ell, u0, eps, B=36):
    sol = integrate.solve_ivp(lambda l, y: [eps * y[0] - B * y[0] ** 2], (0, ell), [u0], rtol=1e-12, atol=1e-15)
    return float(sol.y[0, -1])


def u_from_1240(u0, u1, eps, B=36):
    """(12.40) 左邊的積分 ∫_{u0}^{u1} du/(εu − Bu²)（數值積分；應等於 ℓ）。"""
    return integrate.quad(lambda x: 1 / (eps * x - B * x * x), u0, u1, epsabs=1e-13, epsrel=1e-12)[0]


U_CASES = ((0.005, 1.0), (0.06, 1.0))   # (u0, ε)：卡上的 u(ℓ)，ℓ = 1、2、5 → u* = 1/36


def int_1248(u0, u1, eps, d, book=False):
    """(12.47)–(12.49)：12d∫u dℓ = 12d∫ u du/(εu − 36du²) = −⅓ ln(ε − 36ud)|；book=True 照 (12.48) 少了分子的 u'。"""
    if book:
        f = lambda x: 12 * d / (eps * x - 36 * d * x * x)
    else:
        f = lambda x: 12 * d * x / (eps * x - 36 * d * x * x)
    return integrate.quad(f, u0, u1, epsabs=1e-13, epsrel=1e-12)[0]


def int_1249(u0, u1, eps, d):
    return -(math.log(eps - 36 * u1 * d) - math.log(eps - 36 * u0 * d)) / 3


def t_closed(ell, t0, u0, eps, n=1):
    """(12.51)（O(n)：指數 A/B）：t = t0 e^{2ℓ}[ε/(ε + B u0(e^{εℓ} − 1))]^{A/B}；大 ℓ ∝ e^{(2 − Aε/B)ℓ}。"""
    A, B = AB(n)
    return t0 * math.exp(2 * ell) * (eps / (eps + B * u0 * math.expm1(eps * ell))) ** (A / B)


def t_numeric(ell, t0, u0, eps, n=1):
    """(12.45) dt/dℓ = 2t − A u(ℓ) t 與 du/dℓ = εu − Bu² 一起數值積分。"""
    A, B = AB(n)
    sol = integrate.solve_ivp(lambda l, y: [2 * y[0] - A * y[1] * y[0], eps * y[1] - B * y[1] ** 2], (0, ell), [t0, u0],
                              rtol=1e-12, atol=1e-16)
    return float(sol.y[0, -1])


def flory_nu(d):
    return F(3, 1) / (d + 2)


def flory_series(k):
    """3/(6 − ε) = ½ Σ (ε/6)^k 的第 k 個係數：½、1/12、1/72、1/432。"""
    return F(1, 2) * F(1, 6) ** k


# ======================================================================
# §12.5–12.6：圖解法的 g0(q0) 方程與冪次數算  [W1206]
# ======================================================================
def g_numeric(q0, g0, eps, B=36, eta=0.0):
    """(12.86) q0 dg/dq0 = (−4 + d − 2η)g + B g²（−4 + d = −ε），從 q0 = 1 積到 q0。"""
    sol = integrate.solve_ivp(lambda q, y: [((-eps - 2 * eta) * y[0] + B * y[0] ** 2) / q], (1.0, q0), [g0],
                              rtol=1e-12, atol=1e-15)
    return float(sol.y[0, -1])


def pc_exponent(d, book=False):
    """(12.89) 的冪次：u² r^{d−3}（＝ r·(u r^{−(4−d)/2})²）；book=True 照書上 u² r^{2d−3}。"""
    d = F(d)
    return 2 * d - 3 if book else d - 3


def pc_rhs(d, book=False):
    """右邊的組合：r^{+1}·(u r^{−(4−d)/2})² 的 r 冪次 d − 3；書上 r^{−1}·(…)² 是 d − 5。"""
    d = F(d)
    return (-1 if book else 1) - (4 - d)


def _G_x(x, r, d):
    nu = d / 2 - 1
    k = math.sqrt(r)
    return (2 * math.pi) ** (-d / 2) * (k / x) ** nu * special.kv(nu, k * x)


def watermelon(d, r):
    """Fig. 12.27 的二迴圈（外動量 0）＝ ∫ dᵈx G(x)³（實空間；d < 3 時收斂）。"""
    Sd = 2 * math.pi ** (d / 2) / math.gamma(d / 2)
    f = lambda x: Sd * x ** (d - 1) * _G_x(x, r, d) ** 3
    a = integrate.quad(f, 0, 1 / math.sqrt(r), limit=200, epsabs=0, epsrel=1e-11)[0]
    b = integrate.quad(f, 1 / math.sqrt(r), np.inf, limit=200, epsabs=0, epsrel=1e-11)[0]
    return a + b


def pc_numeric(d, r1=1.0, r2=4.0):
    """數值冪次：ln(I(r2)/I(r1))/ln(r2/r1) → d − 3（d = 2：−1、d = 2.5：−0.5）。"""
    return math.log(watermelon(d, r2) / watermelon(d, r1)) / math.log(r2 / r1)


# ======================================================================
# 全書回顧（W1206 最後一步）
# ======================================================================
def _watson_sc():
    """sc 的 u(1)（Watson 積分）：Gamma 函數閉式與 ∫e^{−t} I0(t/3)³ dt（兩種算法）。"""
    g = special.gamma
    closed = math.sqrt(6) / (32 * math.pi ** 3) * g(1 / 24) * g(5 / 24) * g(7 / 24) * g(11 / 24)
    num = integrate.quad(lambda t: special.ive(0, t / 3) ** 3, 0, np.inf, limit=400)[0]     # e^{−t} I0(t/3)³ = ive³
    return closed, num


def _gambler(L=5, x0=2):
    """一維 0..L 兩端吸收，從 x0 出發先到 0 的機率（線性方程）。"""
    M = np.zeros((L + 1, L + 1))
    b = np.zeros(L + 1)
    M[0, 0] = M[L, L] = 1
    b[0] = 1
    for i in range(1, L):
        M[i, i], M[i, i - 1], M[i, i + 1] = 1, -0.5, -0.5
    return float(np.linalg.solve(M, b)[x0])


def _persist_var(N=1000, Q=0.5):
    """持續漫步 ⟨x²⟩：閉式 N(1+Q)/(1−Q) − 2Q(1−Q^N)/(1−Q)² 與 Σ_{ij} Q^{|i−j|}（兩種算法）。"""
    closed = N * (1 + Q) / (1 - Q) - 2 * Q * (1 - Q ** N) / (1 - Q) ** 2
    k = np.arange(1, N)
    direct = N + 2 * float(np.sum((N - k) * Q ** k))
    return closed, direct


def recap_numbers():
    """全書回顧表每一列的數字（字串）；check() 用兩種算法核對。"""
    wc, wn = _watson_sc()
    pc, pd = _persist_var()
    return {
        1: (_f4(math.comb(10, 5) / 2 ** 10), _f4(2 / math.sqrt(2 * math.pi * 10))),
        2: tuple(str(math.comb(N, N // 2) ** 2) for N in (0, 2, 4, 6, 8, 10)),
        3: (_f4(1 - 1 / wc),),
        4: (str(F(5 - 2, 5)),),
        5: (f"{pc:.0f}",),
        6: (str(F(4 + 2 * 3, 4 + 5 * 3)), _f4(F(10, 19))),
        7: (str(SW.SAW_KNOWN_C[14]),),
        8: (f"{1 / WC.ZC_2D:.3f}",),
        9: (f"{WC.flory_gap_3d():.4f}",),
        10: (str(F(1, 2)),),
        11: (_f4(RG.LN4_LN3),),
        12: (_f4(nu_first(1, 1)), _f4(nu_first(1, 0))),
    }


BOOK_RECAP = [   # (章, 主角, 關鍵結果)；數字見 recap_numbers()
    (1, "一維漫步", r"\langle x^2\rangle=N;\ P_{10}(0)=\tfrac{252}{1024}=0.2461\ \text{vs}\ 0.2523"),
    (2, "生成函數與格子", r"C(N;\vec 0)=\binom{N}{N/2}^2:\ 1,4,36,400,4900,63504"),
    (3, "回家問題", r"\text{P\'olya: } d=1,2\ \text{return};\ \text{sc } 0.3405"),
    (4, "牆與邊界", r"\text{images};\ x_0=2,\ L=5:\ 3/5"),
    (5, "偏向與持續", r"l_p=\tfrac{Q}{1-Q}l;\ Q=\tfrac12,\ N=1000:\ 2996"),
    (6, "漫步的形狀", r"A_d=\tfrac{4+2d}{4+5d};\ A_3=\tfrac{10}{19}=0.5263"),
    (7, "自迴避", r"c_{14}=2374444;\ \nu_F=\tfrac{3}{d+2}"),
    (8, "長尾與標度", r"d_c=4;\ \mu=2.638"),
    (9, "臨界現象", r"n\to0;\ \nu_F=0.600\ \text{vs}\ 0.588\ (2\%)"),
    (10, "受限的高分子", r"\rho=|r|/2u;\ \nu'=\tfrac12"),
    (11, "重整化群", r"\tfrac{\ln4}{\ln3}=1.2619;\ \nu=1/\Lambda_r"),
    (12, "ε 展開", r"\nu=\tfrac12+\tfrac{(n+2)\epsilon}{4(n+8)}:\ 0.5833,\ 0.5625"),
]
BOOK_SCENES = {1: 12, 2: 12, 3: 8, 4: 8, 5: 7, 6: 7, 7: 6, 8: 4, 9: 6, 10: 5, 11: 6, 12: 6}   # 各章幕數（共 87）


def scene_range(ch):
    """第 ch 章的場景代號範圍，例如 'W101–W112'、'W1201–W1206'。"""
    return f"W{ch}01–W{ch}{BOOK_SCENES[ch]:02d}"


RECAP_EXPECT = {1: ("0.2461", "0.2523"), 2: ("1", "4", "36", "400", "4900", "63504"), 3: ("0.3405",), 4: ("3/5",),
                5: ("2996",), 6: ("10/19", "0.5263"), 7: ("2374444",), 8: ("2.638",), 9: ("0.0204",), 10: ("1/2",),
                11: ("1.2619",), 12: ("0.5833", "0.5625")}


def saw_size(N, nu):
    """R ∝ N^ν 的例子：N = 1000 → √N = 31.62、N^{0.5625} = 48.70、N^{0.588} = 58.08。"""
    return N ** float(nu)


# ======================================================================
# 摘要與自我測試
# ======================================================================
def _f6(x):
    return f"{float(x):.6f}"


def summary():
    """場景卡與講義用的全部數字（字串已照卡上的小數位格式化）。"""
    r, u = TOY
    rp, up, c6 = toy_coeffs(r, u)
    z2, z1 = toy_Z(r, u)
    out = {"toy": (_f4(rp), _f4(up), _f6(c6), _f4(toy_Heff(1.0, r, u)), _f4(z1), _f4(toy_A(r)), _f4(z1 / toy_A(r)))}
    out["shell"] = {c: (_f4(shell_fraction(*c)), _f4(shell_linear(*c))) for c in SHELL_CASES}
    out["dq"] = tuple(_f6(x) for x in shell_dq(0.5, 0.01))
    out["rescale"] = tuple(_f4(x) for x in rescale_factors(3, 0.1)) + (_f4(rescale_factors(3, 0.1, book=True)[0]),)
    out["census"] = {n: type_census(n) for n in (3, 4, 5)}
    out["diagrams"] = {n: diagram_census(n) for n in (3, 4, 5, 6)}
    out["f_k"] = tuple(ring_aut(k) for k in range(1, 6))
    out["w_k"] = tuple(ring_weight(k) for k in range(1, 6))
    out["egf"] = _f4(egf_check(0.1)[2])
    a, s2 = LC
    out["lc"] = (_f4(lc_average(a, s2)), _f4(math.log(lc_average(a, s2))), tuple(_f4(x) for x in lc_order(a, s2, 2)),
                 tuple(_f4(x) for x in lc_cumulants_rings(a, s2, 4)), _f4(lc_multi()))
    out["on"] = {n: (2 * n + 4, 4 * (n + 8)) for n in range(4)}
    out["fp"] = {float(e): tuple(_f4(x) for x in fixed_point(e) + fixed_point_lowest(e)) for e in NU_EPS}
    out["fp_on"] = {n: fixed_point(F(1), n) for n in range(4)}
    out["eig"] = tuple(_f4(x) for x in eigen(1.0) + eigen(1.0, 1, "gauss"))
    out["nu_ising"] = {float(e): (_f4(nu_first(e)), _f4(nu_eigen(float(e))), _f4(nu_unexpanded(e))) for e in NU_EPS}
    out["nu_on"] = {n: (_f4(nu_first(1, n)), _f4(nu_eigen(1.0, n)), _f4(nu_unexpanded(1, n)), NU_3D_CITED[n]) for n in range(4)}
    out["rc"] = (_f6(r_critical(1e-4, 0.1)), _f6(r_critical_linear(1e-4, 0.1)), _f6(-6e-4),
                 _f4(r_critical(1e-3, 1.0)), _f4(r_critical_linear(1e-3, 1.0)))
    out["u"] = {u0: tuple(_f4(u_analytic(l, u0, e)) for l in (1, 2, 5)) for u0, e in U_CASES}
    out["u_tiny"] = (_f4(u_analytic(2, 1e-6, 1.0) / 1e-6), _f4(math.exp(2)))
    out["1248"] = (_f4(int_1249(0.001, 0.005, 1.0, 3)), _f4(int_1248(0.001, 0.005, 1.0, 3, book=True)))
    out["t"] = _f4(t_closed(3, 1e-4, 0.01, 1.0))
    out["flory"] = (flory_nu(3), tuple(flory_series(k) for k in range(4)))
    out["saw_size"] = tuple(f"{saw_size(1000, x):.2f}" for x in (F(1, 2), nu_first(1, 0), 0.588))
    out["g"] = (_f4(g_numeric(math.exp(-5), 0.005, 1.0)), _f4(F(1, 36)))
    out["pc"] = {d: (pc_exponent(d), pc_exponent(d, True), pc_rhs(d, True)) for d in (2, F(5, 2), 3, 4)}
    out["pc_num"] = (_f4(pc_numeric(2)), _f4(pc_numeric(2.5)))
    out["recap"] = recap_numbers()
    return out


def check():
    import warnings
    warnings.simplefilter("ignore")
    # ---- W1201：兩變數玩具（閉式＝數值積分＝擬合）、殼層（指數＝線性＝格點計數）、重新標度
    r, u = TOY
    rp, up, c6 = toy_coeffs(r, u)
    rn, un = toy_coeffs_numeric(r, u)
    assert (_f4(rp), _f4(up), _f6(c6)) == ("0.3000", "0.0075", "0.000167")
    assert _close(rp, rn, 1e-7) and _close(up, un, 1e-6)
    for s in (0.3, 1.0, 2.0):
        assert _close(toy_Heff(s, r, u), toy_Heff_numeric(s, r, u), 1e-10)
    assert _f4(toy_Heff(1.0, r, u)) == "0.1577"
    z2, z1 = toy_Z(r, u)
    assert _close(z2, z1, 1e-8) and (_f4(z1), _f4(toy_A(r)), _f4(z1 / toy_A(r))) == ("9.0517", "2.2882", "3.9558")
    exp = {(3, 0.1): ("0.2592", "0.3000"), (3, 0.01): ("0.0296", "0.0300"), (4, 0.01): ("0.0392", "0.0400")}
    for (d, dl), s in exp.items():
        f = shell_fraction(d, dl)
        radial = 1 - integrate.quad(lambda q: d * q ** (d - 1), 0, math.exp(-dl))[0]
        assert (_f4(f), _f4(shell_linear(d, dl))) == s and _close(f, radial, 1e-12)
        if d == 3:
            assert abs(shell_fraction_lattice(3, dl) - f) < 2e-4
    assert abs(shell_fraction_lattice(2, 0.1, 400) - shell_fraction(2, 0.1)) < 2e-4
    dq = shell_dq(0.5, 0.01)
    assert tuple(_f6(x) for x in dq) == ("0.006035", "0.006065") and _f4(math.exp(-0.5)) == "0.6065"
    qf, rf, uf = rescale_factors(3, 0.1)
    assert _close(qf, 1) and (_f4(rf), _f4(uf), _f4(rescale_factors(3, 0.1, book=True)[0])) == ("1.2214", "1.1052", "0.6703")
    for d in (2, 3, 4):                                  # 一階：(1, 1 + 2δℓ, 1 + (4−d)δℓ)
        q1, r1, u1 = rescale_factors(d, 1e-4)
        assert abs(q1 - 1) < 1e-12 and abs(r1 - (1 + 2e-4)) < 1e-7 and abs(u1 - (1 + (4 - d) * 1e-4)) < 1e-7
    assert _f4(rescale_factors(3, 0.01, book=False)[1]) == "1.0202"
    assert (_f4(math.exp(3 * 0.01)), _f4(1 + 0.01)) == ("1.0305", "1.0100")       # (12.32) 的筆誤：e^{dδℓ} vs 1 + δℓ
    # ---- W1202：配對列舉＝組合公式、連結叢集
    for n in range(1, 7):
        dfac = special.factorial2(2 * n - 1, exact=True)
        assert len(matchings(n)) == dfac == sum(type_census(n).values())
        tc, dc = type_census(n), diagram_census(n)
        for t, cnt in tc.items():
            cts = Counter(t)
            legs = 1
            for k, c in cts.items():
                legs *= (1 if k == 1 else 2 if k == 2 else 2 ** k) ** c
            assert cnt == dc[t] * legs and dc[t] == general_factor(dict(cts))
        for n2 in range(n // 2 + 1):
            assert pair_factor(n - 2 * n2, n2) == dc[sizes_of(n - 2 * n2, n2)]
    for k in range(1, 7):
        assert type_census(k)[(k,)] == ring_weight(k)
    assert len(matchings(5)) == 945 and type_census(5)[(2, 1, 1, 1)] == 20 and diagram_census(5)[(2, 1, 1, 1)] == 10
    assert pair_factor(3, 1) == F(math.factorial(5), math.factorial(3) * math.factorial(2)) == 10
    assert F(math.factorial(5), math.factorial(2) * math.factorial(1)) == 60                     # p. 320 的筆誤
    assert type_census(3) == {(1, 1, 1): 1, (2, 1): 6, (3,): 8} and diagram_census(3) == {(1, 1, 1): 1, (2, 1): 3, (3,): 1}
    assert type_census(4) == {(1, 1, 1, 1): 1, (2, 1, 1): 12, (3, 1): 32, (2, 2): 12, (4,): 48}
    assert (pair_factor(2, 1), pair_factor(0, 2), pair_factor(4, 1), pair_factor(2, 2), pair_factor(0, 3)) == (6, 3, 15, 45, 15)
    assert tuple(ring_aut(k) for k in range(1, 6)) == (1, 2, 6, 8, 10) and ring_aut(3) == math.factorial(3) != 0
    assert ring_aut(4) == 8 != math.factorial(4)
    assert tuple(ring_weight(k) for k in range(1, 6)) == (1, 2, 8, 48, 384)
    left, right, closed = egf_check(0.1)
    assert _close(left, closed, 1e-12) and _close(right, closed, 1e-12) and _f4(closed) == "1.1180"
    a, s2 = LC
    vals = [lc_average(a, s2, m) for m in ("closed", "quad", "series", "linked")]
    assert all(_close(v, vals[0], 1e-12) for v in vals) and _f4(vals[0]) == "0.9129"
    assert _f4(math.log(vals[0])) == "-0.0912" and _close(math.log(vals[0]), -0.5 * math.log(1.2), 1e-12)
    assert tuple(_f4(x) for x in lc_order(a, s2, 2)) == ("0.0150", "0.0100", "0.0050")
    kr, km = lc_cumulants_rings(a, s2, 8), lc_cumulants_moments(a, s2, 8)
    assert all(_close(x, y, 1e-10) for x, y in zip(kr, km))
    assert tuple(_f4(x) for x in kr[:4]) == ("-0.1000", "0.0200", "-0.0080", "0.0048")
    al, co, di = lc_order(a, s2, 3)
    assert _close(di, kr[0] * kr[1] / 2 + kr[0] ** 3 / 6, 1e-12)                                  # 不連通＝連通的乘積
    mv = [lc_multi(method=m) for m in ("det", "quad", "linked")]
    assert all(_close(v, mv[0], 1e-10) for v in mv) and _f4(mv[0]) == "0.8087"
    oc = {n: on_counts(n) for n in (1, 2, 3, 4)}
    for n, (c2, c4) in oc.items():
        assert _close(c2, 2 * n + 4, 1e-12) and _close(c4, 4 * (n + 8), 1e-12)
        assert (2 * (2 * n + 4), round(c4)) == AB(n)
    p2 = np.polyfit([1, 2, 3, 4], [oc[n][0] for n in (1, 2, 3, 4)], 1)
    p4 = np.polyfit([1, 2, 3, 4], [oc[n][1] for n in (1, 2, 3, 4)], 1)
    assert _close(np.polyval(p2, 0), 4, 1e-9) and _close(np.polyval(p4, 0), 32, 1e-9) and AB(0) == (8, 32)
    # ---- W1203：不動點（Fraction 精確＝fsolve）、線性化（閉式＝numpy＝差分）
    assert AB(1) == (12, 36) and tuple(AB(n) for n in range(4)) == ((8, 32), (12, 36), (16, 40), (20, 44))
    fps = {F(1): ((F(-1, 7), F(1, 49)), ("-0.1429", "0.0204", "-0.1667", "0.0278")),
           F(1, 2): ((F(-1, 13), F(2, 169)), ("-0.0769", "0.0118", "-0.0833", "0.0139")),
           F(1, 10): ((F(-1, 61), F(10, 3721)), ("-0.0164", "0.0027", "-0.0167", "0.0028"))}
    for e, (fr, s) in fps.items():
        fp = fixed_point(e)
        assert fp == fr and beta(fp[0], fp[1], e) == (0, 0) and fp == (-e / (6 + e), e / (6 + e) ** 2)
        assert tuple(_f4(x) for x in fp + fixed_point_lowest(e)) == s
        fn = fixed_point_numeric(float(e))
        assert _close(fn[0], fp[0], 1e-10) and _close(fn[1], fp[1], 1e-10)
        rb, ub = fixed_point(e)                                 # 書上的 d 形：u_book = u/d → 一階 ε/36d
        db = 4 - e
        assert all(abs(x) < 1e-15 for x in beta_book(float(rb), float(ub / db), float(e)))
    assert {n: fixed_point(F(1), n) for n in range(4)} == {0: (F(-1, 9), F(2, 81)), 1: (F(-1, 7), F(1, 49)),
                                                          2: (F(-1, 6), F(5, 288)), 3: (F(-5, 27), F(11, 729))}
    for n in range(4):
        for e in (0.1, 0.5, 1.0):
            rs, us = fixed_point(e, n)
            A, B = AB(n)
            assert _close(us / (1 + rs) ** 2, e / B, 1e-12)                    # ⇒ ∂ṙ/∂r = 2 − Aε/B 精確
            J, Jf = jacobian(rs, us, e, n), jacobian_fd(rs, us, e, n)
            assert np.allclose(J, Jf, atol=1e-7)
            assert _close(eigen(e, n)[0], lambda_r_closed(e, n), 1e-12)
    r1, u1 = fixed_point(F(1))
    Jx = [[2 - F(12, 36), 12 / (1 + r1)], [2 * 36 * u1 ** 2 / (1 + r1) ** 3, 1 - 2 * 36 * u1 / (1 + r1) ** 2]]
    assert Jx == [[F(5, 3), F(14)], [F(1, 21), F(-1)]] and np.allclose(jacobian(float(r1), float(u1), 1.0), np.array(Jx, float))
    assert tuple(_f4(x) for x in eigen(1.0)) == ("1.8968", "-1.2301") and _close(eigen(1.0)[0], (1 + math.sqrt(22)) / 3, 1e-12)
    assert eigen(1.0, 1, "gauss") == (2.0, 1.0) and eigen(0.5, 1, "gauss") == (2.0, 0.5)
    starts = flow_starts()
    for (r0, u0), k in zip(starts, range(6)):
        p = flow_line(r0, u0)
        assert len(p) == 40
        if k % 3 == 1:                                          # 臨界線流進 WF 點
            assert abs(p[-1][0] - float(r1)) < 2e-3 and abs(p[-1][1] - float(u1)) < 2e-4
        else:
            assert abs(abs(p[-1][0]) - 0.45) < 1e-6
    # ---- W1204：u(ℓ)、t(ℓ)、ν、Ex. 12.2–12.4
    exp_u = {0.005: ("0.0104", "0.0172", "0.0270"), 0.06: ("0.0346", "0.0300", "0.0279")}
    for u0, e in U_CASES:
        vals = tuple(u_analytic(l, u0, e) for l in (1, 2, 5))
        assert tuple(_f4(x) for x in vals) == exp_u[u0]
        assert all(_close(v, u_numeric(l, u0, e), 1e-9) for v, l in zip(vals, (1, 2, 5)))
        assert _close(u_from_1240(u0, vals[1], e), 2.0, 1e-9)                                        # Ex. 12.2
        assert abs(u_analytic(40, u0, e) - 1 / 36) < 1e-12
    assert _f4(F(1, 36)) == "0.0278"
    h = 1e-7
    bu = lambda x: 1.0 * x - 36 * x * x
    assert _close((bu(h) - bu(-h)) / (2 * h), 1.0, 1e-9) and _close((bu(1 / 36 + h) - bu(1 / 36 - h)) / (2 * h), -1.0, 1e-7)
    assert (_f4(u_analytic(2, 1e-6, 1.0) / 1e-6), _f4(math.exp(2))) == ("7.3874", "7.3891")       # Ex. 12.3
    i_ok, i_closed, i_book = int_1248(0.001, 0.005, 1.0, 3), int_1249(0.001, 0.005, 1.0, 3), int_1248(0.001, 0.005, 1.0, 3, True)
    assert _close(i_ok, i_closed, 1e-10) and (_f4(i_ok), _f4(i_book)) == ("0.2207", "81.7804")
    for n in (0, 1):
        assert _close(t_closed(3, 1e-4, 0.01, 1.0, n), t_numeric(3, 1e-4, 0.01, 1.0, n), 1e-8)
    assert _f4(t_closed(3, 1e-4, 0.01, 1.0)) == "0.0203"
    slope = math.log(t_closed(31, 1e-20, 0.01, 1.0) / t_closed(30, 1e-20, 0.01, 1.0))
    assert _close(slope, 2 - 1 / 3, 1e-9)
    exp_nu = {F(1, 10): ("0.5083", "0.5076", "0.5085"), F(1, 2): ("0.5417", "0.5256", "0.5455"), F(1): ("0.5833", "0.5272", "0.6000")}
    for e, s in exp_nu.items():
        assert (_f4(nu_first(e)), _f4(nu_eigen(float(e))), _f4(nu_unexpanded(e))) == s
        assert nu_unexpanded(e) == flory_nu(4 - e)                     # n = 1 不展開 = Flory 3/(d+2)
    assert nu_first(1) == F(7, 12) and nu_first(1, 0) == F(9, 16) and _close(nu_from_flow(1.0), nu_eigen(1.0), 1e-5)
    d1, d2 = nu_eigen(0.01) - float(nu_first(F(1, 100))), nu_eigen(0.02) - float(nu_first(F(1, 50)))
    assert abs(d2 / d1 - 4) < 0.1                                      # 差是 O(ε²)
    u1d = float(nu_unexpanded(F(1, 100))) - float(nu_first(F(1, 100)))
    assert 0 < u1d < 1e-3
    rc = r_critical(1e-4, 0.1)
    assert (_f6(rc), _f6(r_critical_linear(1e-4, 0.1)), _f6(-6e-4)) == ("-0.000631", "-0.000632", "-0.000600")
    assert flow_end(rc + 1e-7, 1e-4, 0.1) == 1 and flow_end(rc - 1e-7, 1e-4, 0.1) == -1
    rc1 = r_critical(1e-3, 1.0)
    assert (_f4(rc1), _f4(r_critical_linear(1e-3, 1.0))) == ("-0.0111", "-0.0120")
    # ---- W1205：O(n) 表、Flory 的 ε 展開
    exp_on = {0: ("0.5625", "0.5205", "0.5714"), 1: ("0.5833", "0.5272", "0.6000"),
              2: ("0.6000", "0.5325", "0.6250"), 3: ("0.6136", "0.5368", "0.6471")}
    for n, s in exp_on.items():
        assert (_f4(nu_first(1, n)), _f4(nu_eigen(1.0, n)), _f4(nu_unexpanded(1, n))) == s
        assert _close(nu_eigen(1.0, n), 1 / lambda_r_closed(1.0, n), 1e-12)
    assert tuple(nu_first(1, n) for n in range(4)) == (F(9, 16), F(7, 12), F(3, 5), F(27, 44))
    assert NU_3D_CITED[0] == WC.NU_KNOWN[3] == 0.588
    assert flory_nu(3) == F(3, 5) and tuple(flory_series(k) for k in range(4)) == (F(1, 2), F(1, 12), F(1, 72), F(1, 432))
    for k in range(1, 4):                                              # 數值 Taylor 係數
        e = 1e-3
        partial = sum(float(flory_series(j)) * e ** j for j in range(k))
        assert _close((3 / (6 - e) - partial) / e ** k, float(flory_series(k)), 2e-2)
    assert sum(flory_series(k) for k in range(2)) == F(7, 12) == nu_first(1, 1) and F(1, 12) != F(1, 16)
    assert tuple(f"{saw_size(1000, x):.2f}" for x in (F(1, 2), nu_first(1, 0), 0.588)) == ("31.62", "48.70", "58.08")
    # ---- W1206：(12.86) ＝ (12.39)、冪次數算
    for l in (1.0, 5.0):
        assert _close(g_numeric(math.exp(-l), 0.005, 1.0), u_analytic(l, 0.005, 1.0), 1e-8)
    assert _f4(g_numeric(math.exp(-5), 0.005, 1.0)) == "0.0270" and _f4(math.exp(-5)) == "0.0067"
    for d in (2, F(5, 2), 3, 4):
        assert pc_exponent(d) == pc_rhs(d) == F(d) - 3 and pc_exponent(d, True) == 2 * F(d) - 3 and pc_rhs(d, True) == F(d) - 5
        assert pc_exponent(d) == 3 * F(d) / 2 - F(d) / 2 - 3                     # 3 個 dᵈQ、1 個 δ、3 個 G
    assert (_f4(pc_numeric(2)), _f4(pc_numeric(2.5))) == ("-1.0000", "-0.5000")
    # ---- 全書回顧
    rn = recap_numbers()
    assert rn == RECAP_EXPECT
    wc, wn = _watson_sc()
    assert _close(wc, wn, 1e-8)
    assert _close(_gambler(), 0.6, 1e-12)
    pcl, pdi = _persist_var()
    assert _close(pcl, pdi, 1e-10) and _close(pcl, 2996, 1e-12)
    assert all(SW.count2(N, 0, 0) == math.comb(N, N // 2) ** 2 for N in (0, 2, 4, 6, 8, 10))
    assert _close(RG.replica_dim(4, 3), RG.LN4_LN3, 1e-12)
    assert len(BOOK_RECAP) == 12 and [row[0] for row in BOOK_RECAP] == list(range(1, 13))
    assert sum(BOOK_SCENES.values()) == 87 and scene_range(1) == "W101–W112" and scene_range(12) == "W1201–W1206"
    assert scene_range(8) == "W801–W804" and scene_range(10) == "W1001–W1005"
    for ch, vals in RECAP_EXPECT.items():
        assert all(v in BOOK_RECAP[ch - 1][2] or v.replace("/", "}{") in BOOK_RECAP[ch - 1][2] for v in vals if ch not in (4, 6, 9, 10)), ch


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
    print(f"walk_eps（第 12 章）：兩變數玩具（閉式＝數值積分＝擬合）、殼層體積（指數＝線性＝格點計數）、重新標度、配對列舉＝組合公式（(12:S-4)(12:S-6)、f₃ = 6）、"
          f"連結叢集（四種算法、三變數 det＝求積＝環的和）、O(n) 縮併 2n+4 與 4(n+8)、不動點（Fraction＝fsolve）、線性化（閉式＝numpy＝差分＝流的時間）、"
          f"u(ℓ) 與 t(ℓ)（解析＝數值）、ν 的一階公式與 Flory 展開、(12.86)＝(12.39)、冪次數算（Fraction＝數值積分）、全書回顧的數字全部驗證通過（{time.time() - t0:.1f} 秒）")
