# 由 combkit/tools/sync_from_walk017.py 從 walk017 專案的 codes/utils/walk_rg.py 複製；請改 walk017 那份再同步。
"""
walk_rg.py —— 《Elements of the Random Walk》第 11 章 Scaling, fractals, and renormalization 的數學
（Python＋numpy／scipy（integrate、linalg），不 import manim；用到同資料夾的 walk_saw（四方向步伐、全部 4^N 條漫步、2ν_eff）與
walk_crit（ν 的已知值、指數表、Landau 的根、Widom 指數、Ising 列舉）——在 walk017 是 utils 套件、在 combkit 是 combkit 套件，
都用相對 import；直接執行本檔時退回絕對 import）

約定（書上式號、筆誤與〔注〕見 ref/md/ch11_formulas.md）：
  Landau 用第 9 章的寫法 F/V = r m²/2 + u m⁴/4 − h m（書上 (11.9) 是 V[r m² + u m⁴ − h m]，常數不影響標度），u = 1。
  Gaussian 模型（連續）只有 q = 0 的模態受外場：m = h/r（r > 0）。一維 Gaussian 鏈 ℋ = Σ[r/2 S_i² + a/2 (S_i − S_{i+1})²]，
  書上的 a 取 A_BOOK = 1（Fig. 11.6 的座標其實是 r/a）。一維 Ising：K = J/k_BT、H = h/k_BT，週期環。
  抽格 b = 2；連續重整化 b = e^ℓ；特徵值 Λ 用 e^{Λℓ}（離散時 Λ = ln(特徵值)/ln b）。指數用 Fraction。
  回傳快取（lru_cache）陣列的函數不要就地修改結果（已設成唯讀）。

§11.1.1 碎形 (11.1)–(11.6)，Fig. 11.1–11.5，Ex. 11.1–11.3                                           [W1101]
    LN4_LN3 = 1.2619、LN2_LN3 = 0.6309
    koch_points(gen, p0=(0,0), p1=(1,0))   第 gen 代頂點 (4^gen + 1, 2)（突起在前進方向左邊；gen ≤ 5 給畫圖用就夠）
    koch_counts(gen) = (4^gen, (4/3)^gen, 3^−gen)（Fraction）：總長 1, 1.3333, 1.7778, 2.3704, 3.1605, 4.2140（gen 0..5）→ ∞
    cantor_intervals(gen)          [(a, b)]（Fraction，2^gen 段）；cantor_counts(gen) = (2^gen, (2/3)^gen, 3^−gen)：
                                   總長 1, 0.6667, 0.4444, 0.2963, 0.1975, 0.1317 → 0
    koch_mass(R)、cantor_mass(R)   以起點為圓心、半徑 R 內的份量（全部 = 1）；mass_table(kind) → [(k, R = 3^−k, n)]（k = 0..4）：
                                   Koch n = 1, 1/4, 1/16, …（精確：第 2–4 段都在 R = 1/3 外，30° 錐）、Cantor n = 1, 1/2, 1/4, …
    mass_radius_dim(kind, Rs=None) ln n 對 ln R 的斜率：R = 3^−k → 1.2619、0.6309（精確）；R = 2^−j（j = 0..8）的 Koch → 1.2571（對數週期抖動）
    box_counts(kind, ks, base)、box_dim(kind)   盒子計數（第二種算法；格線位移 BOX_OFFSET）：
                                   Koch ε = 3^−k（k = 2..5）盒子數 22、83、325、1291 → 斜率 1.2362；Cantor ε = 2^−k（k = 2..10）→ 0.6254
    REPLICAS、replica_dim(ℛ, ℳ) = ln ℛ/ln ℳ（(11.6)）：line (3,3) 1.0000、square (9,3) 2.0000、koch (4,3) 1.2619、
                                   cantor (2,3) 0.6309、sierpinski (3,2) 1.5850
    rw_paths(Nmax, M, seed=1101)、rw_rg2(N)（M = 400 條二維普通漫步前 N 步的 ⟨R_g²⟩）、rg2_exact(N) = N(N+2)/(6(N+1))（Fraction；
                                   N ≤ 5 與 4^N 全部列舉相同）；rw_df() = ln N 對 ln R_g 的斜率（RW_NS = 64, 256, 1024, 4096）：
                                   固定種子 2.0258、公式 2.0069 → 2（普通漫步的碎形維度）；rw_trail(N=1500)：畫圖用一條
    saw_df(d) = 1/ν：d = 2 → 4/3、d = 3 → 1.70（ν = 0.588 只引用）、d = 4 → 2；saw_df_eff(N) = 2/(2ν_eff)：
                                   N = 10、12、14 → 1.3848、1.3788、1.3740（→ 4/3；二維精確列舉，walk_saw）
§11.1.2 被動的重整化：Landau (11.9)–(11.17) 與 Gaussian (11.18)–(11.25)，Ex. 11.4                       [W1102]
    landau_m(r, h)、landau_m_bisect（二分法）＝ walk_crit.landau_min：m(−0.25, 0) = 0.5000 = √(−r)；m(0, 0.008) = 0.2000、m(0, 0.125) = 0.5000（h^{1/3}）
    landau_scaled(x) = f(x)：y + y³ = x（y = m/r^{1/2}、x = h/r^{3/2}，r > 0）：f(1) = 0.6823、f(2) = 1、f(10) = 2；f(x) → x^{1/3}
    COLLAPSE_R = (0.25, 0.5, 1)；landau_curve(r, n=40)（h ∈ [0, 2]）、collapse_curve(r)：三條疊成一條；
                                   COLLAPSE_POINTS (0.25, 0.25)、(0.5, 0.7071)、(1, 2)：m = 0.5000、0.7071、1.0000，都是 x = 2 → y = 1
    mf_rescale_check()             m(r, h) = b^{−d/4} m(r b^{d/2}, h b^{3d/4})（(11.11)）的誤差
    book1116_limit(r, h) → (f, 書上 h^{1/3} r, 更正 h^{1/3} r^{−1/2})：r = 1e−4、h = 1e−3 → 9.9667、0.0000、10.0000
    gauss_m(r, h) = h/r；gauss_rescale_check(d)：(11.22) 與 (11.24) 的誤差（任何 b、d）
    gauss_exponents(d)（Ex. 11.4）：β = (d−2)/4、γ = 1、δ = (d+2)/(d−2)、α = 2 − d/2、ν = 1/2：d = 3 → 1/4、1、5（α = 1/2）；d = 4 → 1/2、1、3（＝平均場）
    scaling_combos(d) = (3/2, (d+2)/4, 1/2, (d−2)/4)：d = 4 時 Landau 與 Gaussian 的組合相同
§11.2 一維 Gaussian 鏈的抽格 (11.26)–(11.49)，Fig. 11.6–11.7，Ex. 11.5–11.6                          [W1103]
    square_check(r=0.5, a=1, S2=0.3, S4=−0.7, h=0) → (數值積分, 完成平方)：1.224819、1.224819；square_book_form（(11.30) 照書，對的）
    decimate_gauss(r, a, h) = (r(r+4a)/(r+2a), a²/(r+2a), h(r+4a)/(r+2a))（(11.32)(11.33)，未重新標度）
    spin_factor(r, a) = √((r+2a)/a)（舊自旋 = 因子 × 新自旋；小 r → √2）
    rr_gauss(r, a=1) = r(r+4a)/a（(11.36)）；hh_gauss(r, h, a=1) = h(r+4a)/√(a(r+2a))（(11.38)；hh_gauss(0, 1) = 2.8284 = 2^{3/2}）
    gauss_rescaled(r, a, h)        decimate_gauss × 自旋因子（第二種算法）＝ (rr_gauss, a, hh_gauss)
    ring_lnZ(N, r, a, h)           Gaussian 環的 ln Z（行列式）；ring_decimation_check(N, r, a, h, book=False)：
                                   ln Z_N ＝ (N/4)ln(2π/(r+2a)) + (N/2)h²/(2(r+2a)) + ln Z_{N/2}(r', a', h') 精確成立；
                                   書上 (11.31) 的 (2π/(r+2a))^{N/2} 多了 (N/4)ln(2π/(r+2a))（N = 8、r = 0.5：1.8432）
    gauss_m_invariance(r, h)       h/r ＝ 自旋因子 × h''/r''（精確）
    book1147_m(r, h, n) → (書上 2^{−n/2} f, 更正 2^{+n/2} f, h/r)：r = 0.01、h = 0.001、n = 1 → 0.0500、0.1000、0.1000
    gauss_iterate(r0, n)、COBWEB_R0 = 0.05：0.0500、0.2025、0.8510、4.1282、33.5553；cobweb_points(r0, n=3)：蛛網折線（起點 (r0, 0)，6 點）
    rr_curve(n=40)                 Fig. 11.6 的 (r, r'')，r ∈ [0, 1.7417]（r'' 到 10）
    small_r_table(n=5, r=h=1e−4)   [(k, r^(k)/(4^k r), h^(k)/(2^{3k/2} h), 自旋總因子/2^{k/2})]：k = 1..5：
                                   r：1.0000、1.0001、1.0005、1.0021、1.0086；h：都是 1.0000；自旋：1.0000、1.0001、1.0005、1.0021、1.0085
    gauss_fixed_points() = (0, −3a)（−3a 不物理）；chain_lambdas() = (2.0000, 1.5000, 0.5000)：Λ_r = 2、Λ_h = 3/2 = (d+2)/2、自旋 (2−d)/2（d = 1）
§11.2 例題：一維 Ising 的抽格 (11.50)–(11.63)，Ex. 11.7–11.9                                        [W1104]
    ising_A(K) = 2√cosh 2K、ising_Kp(K) = ½ ln cosh 2K：K = 1 → 3.8793、0.6625；ising_pair_check(K)（四種 (S2, S4) 的誤差）
    book1153(K) → (書上左邊 cosh 2K, 書上右邊 2cosh 2K, 正確左邊 2cosh 2K)：K = 1 → 3.7622、7.5244、7.5244
    ex118(K, H) → (A, K1', K2', H' = H + 2K2')（Ex. 11.8）：K = 1、H = 0.1 → 3.8897、0.6602、0.0482、0.1964；ex118_check
    ising_iterate(K0=1, n) = 1.0000、0.6625、0.3501、0.1137、0.0128、0.0002（→ 0）
    ring_Z_brute(N, K, H=0)（2^N 列舉）＝ ring_Z_transfer(N, K) = (2cosh K)^N + (2sinh K)^N ＝ ring_Z_decimation(N, K)（(11.63) 連乘到剩
                                   1 個自旋 2e^{K_{k+1}}）：K = 1：N = 4 → 121.2329（＝第 9 章 W902 的 4 自旋環）、N = 8 → 9160.4387；
                                   book=True 照書上指數 N/2^{j+1}：18.5498、143.2705
    x_of_K(K) = e^{−2K}、x_map(x) = 2x/(1+x²)、x_iterate(x0, n)：x0 = e^{−2}：0.1353、0.2658、0.4965、0.7966、0.9747、0.9997（→ 1）
    x_fixed_points() = (0, 1)；x_slope(x) = 2(1−x²)/(1+x²)²：x = 0 → 2（不穩定，T = 0）、x = 1 → 0（穩定，T = ∞）
    xi_exact(K) = −1/ln tanh K；xi_ratio(K) = ξ/e^{2K}：K = 1、2、3：ξ = 3.6719、27.2960、201.7140；e^{2K}/2 = 3.6945、27.2991、201.7144；
                                   比 0.4969、0.4999、0.5000；tanh_rg_check(K)：tanh K' = tanh² K ⇒ ξ(K') = ξ(K)/2（精確）
    free_energy_ratio(K, N=None) = ln Z/(NK)：K = 1、2、5 → 1.1269、1.0091、1.0000（F → −NJ）；N = 16 的環（抽格連乘）1.1277；
                                   free_energy_series(K) = Σ_j 2^{−j} ln A(K_j)（第二種算法）
§11.2.1–11.2.2 動量殼與連續重整化 (11.64)–(11.87)，Fig. 11.8                                        [W1105]
    alpha_q(q, a=1) = a(1 − cos q)（≈ aq²/2）；alpha_q_book = −2a(1 − cos q)；alpha_table(N=8) → [(k, q = 2πk/8, α, aq²/2, 書上)]：
                                   α = 0, 0.2929, 1, 1.7071, 2；aq²/2 = 0, 0.3084, 1.2337, 2.7758, 4.9348；書上 0, −0.5858, −2, −3.4142, −4
    fourier_check(N=8, r=0.5, seed=1105) → (相對誤差, 實空間 3.0394, 照書上 α −4.3508)
    shell_fraction(d, b=2) = 1 − b^{−d}：d = 2 → 3/4、d = 3 → 7/8
    shell_rg(r, h, d, b=2) = (b²r, b^{(d+2)/2}h)（(11.71)(11.72)；d = 3：h × 5.6569）；shell_rg_steps（q→q/b、N→b^dN、S→bS 一步步，第二種算法）
    flow(r0, h0, d, ℓ) = (r0e^{2ℓ}, h0e^{(d+2)ℓ/2})（(11.82)(11.83)）＝ flow_numeric（solve_ivp）；ℓ = ln 2 ＝ shell_rg(b = 2)
                                   r0 = h0 = 0.01、d = 3：ℓ = 1 → 0.0739、0.1218；ℓ = 2 → 0.5460、1.4841；flow_line(r0, h0, d, ℓmax, n=40)
    m24(r, h, d, F)（(11.24)，預設 F(x) = x/(1+x²)）；pde_residual(d)：(11.87) 的殘差（中央差分，< 1e−7）；char_check(d)：(11.86) 與 ℓ 無關
§11.3 不動點、線性化與臨界指數 (11.88)–(11.114)，Fig. 11.9–11.10                                    [W1106]
    TOY_A = [[1, 2], [2, −2]]（示意）；toy_eigen() = [(2, (2, 1)), (−3, (1, −2))]（相關、無關）；toy_flow(c1, c2, ℓ)＝toy_flow_expm；
                                   TOY_STARTS、toy_flow_line(c1, c2, ℓmax=1, n=30)（c1 = 0：留在臨界線上流進不動點）
    rg_exponents(Λr, Λh, d)（(11.105)–(11.108)＋ν = 1/Λr）＝ rg_exponents_widom（walk_crit.widom_exponents，Δ = Λh/Λr）；lambdas_from(ν, β, d)
    RG_CASES：gauss3 (2, 5/2, 3) → α 1/2、β 1/4、γ 1、δ 5、ν 1/2；gauss4 (2, 3, 4) → 0、1/2、1、3、1/2（＝WC.EXPONENTS['mf']）；
              ising2d (1, 15/8, 2) → 0、1/8、7/4、15、1；saw2d (4/3, 91/48, 2) → 1/2、5/64、43/32、91/5、3/4（都＝walk_crit.EXPONENTS）
    gauss_lambdas(d) = (2, (d+2)/2)
    xi_from_flow(r, Λr) = (ℓ*, ξ)（K = D = 1）＝ xi_from_flow_numeric：r = 0.01：Λr = 2 → 2.3026、10.0000；Λr = 1 → 4.6052、100.0000
summary()                         場景卡與講義用的全部數字（dict，字串已照卡上的小數位格式化）

check() 驗證全部（每個結果至少兩種獨立算法、筆誤的更正、plan_ch11.md 場景卡與講義的每個數字）；直接執行本檔會跑 check()（約 1 秒）。
`python walk_rg.py --summary` 印出 summary() 的全部數字。最慢的呼叫：box_counts('koch') 約 0.3 秒、rw_rg2 第一次約 0.2 秒（快取），其餘都很快。
"""
import math
from fractions import Fraction
from functools import lru_cache

import numpy as np
from scipy import integrate, linalg

try:                                    # walk017：utils 套件；combkit：combkit 套件
    from . import walk_saw as SW
    from . import walk_crit as WC
except ImportError:                     # 直接執行本檔
    import walk_saw as SW
    import walk_crit as WC

F = Fraction


def _ro(a):
    a = np.asarray(a)
    a.flags.writeable = False
    return a


def _f4(x):
    return f"{float(x):.4f}"


# ======================================================================
# §11.1.1 碎形：Koch 曲線、Cantor 集、維度  [W1101]
# ======================================================================
LN4_LN3 = math.log(4) / math.log(3)          # 1.2619
LN2_LN3 = math.log(2) / math.log(3)          # 0.6309
_ROT60 = np.array([[0.5, -math.sqrt(3) / 2], [math.sqrt(3) / 2, 0.5]])


@lru_cache(maxsize=None)
def _koch_unit(gen):
    pts = np.array([[0.0, 0.0], [1.0, 0.0]])
    for _ in range(gen):
        a, b = pts[:-1], pts[1:]
        d = (b - a) / 3
        p1 = a + d
        p3 = a + 2 * d
        p2 = p1 + d @ _ROT60.T
        new = np.empty((4 * len(a) + 1, 2))
        new[0:-1:4] = a
        new[1::4] = p1
        new[2::4] = p2
        new[3::4] = p3
        new[-1] = pts[-1]
        pts = new
    return _ro(pts)


def koch_points(gen, p0=(0.0, 0.0), p1=(1.0, 0.0)):
    """第 gen 代 Koch 曲線的頂點 (4^gen + 1, 2)，從 p0 到 p1（突起在前進方向的左邊）。"""
    u = _koch_unit(gen)
    p0 = np.asarray(p0, float)
    v = np.asarray(p1, float) - p0
    R = np.array([[v[0], -v[1]], [v[1], v[0]]])
    return p0 + u @ R.T


def koch_counts(gen):
    """(段數 4^gen, 總長 (4/3)^gen, 每段長 3^−gen)（Fraction）。"""
    return 4 ** gen, F(4, 3) ** gen, F(1, 3) ** gen


def cantor_intervals(gen):
    """第 gen 代 Cantor 集的區間 [(a, b)]（Fraction，2^gen 個）。"""
    iv = [(F(0), F(1))]
    for _ in range(gen):
        nxt = []
        for a, b in iv:
            d = (b - a) / 3
            nxt += [(a, a + d), (b - d, b)]
        iv = nxt
    return iv


def cantor_counts(gen):
    """(段數 2^gen, 總長 (2/3)^gen, 每段長 3^−gen)（Fraction）。"""
    return 2 ** gen, F(2, 3) ** gen, F(1, 3) ** gen


@lru_cache(maxsize=None)
def _cantor_mid(gen):
    a = np.array([0.0])
    w = 1.0
    for _ in range(gen):
        w /= 3
        a = np.concatenate([a, a + 2 * w])
    return _ro(a + w / 2)                     # 各區間中點


MASS_GEN = 8


def koch_mass(R, gen=MASS_GEN):
    """以起點 (0, 0) 為圓心、半徑 R 的圓內 Koch 曲線的份量（第 gen 代，以段中點判斷；全長 = 1）。"""
    p = _koch_unit(gen)
    mid = (p[:-1] + p[1:]) / 2
    return float(np.mean(np.hypot(mid[:, 0], mid[:, 1]) < R))


def cantor_mass(R, gen=12):
    """以左端點 0 為中心、半徑 R 的區間內 Cantor 集的份量（第 gen 代；全部 = 1）。"""
    return float(np.mean(_cantor_mid(gen) < R))


MASS_K = (0, 1, 2, 3, 4)                     # R = 3^−k


def mass_table(kind="koch"):
    """[(k, R = 3^−k, n(R))]：Koch n = 4^−k、Cantor n = 2^−k（(11.1)：n(R/3) = n(R)/4）。"""
    fn = koch_mass if kind == "koch" else cantor_mass
    return [(k, 3.0 ** -k, fn(3.0 ** -k * (1 + 1e-9))) for k in MASS_K]


def _slope(x, y):
    return float(np.polyfit(np.log(x), np.log(y), 1)[0])


def mass_radius_dim(kind="koch", Rs=None):
    """ln n(R) 對 ln R 的最小平方斜率。Rs 預設 3^−k（k = 0..4，正好 ln4/ln3、ln2/ln3）；
    給其他 R（例如 2^−j）會看到對數週期的抖動。"""
    fn = koch_mass if kind == "koch" else cantor_mass
    Rs = [3.0 ** -k * (1 + 1e-9) for k in MASS_K] if Rs is None else Rs
    return _slope(Rs, [fn(R) for R in Rs])


BOX_OFFSET = (0.0123, 0.0371)                # 盒子格線的位移（避免頂點剛好落在格線上）
BOX_K = (2, 3, 4, 5)                         # Koch：ε = 3^−k
BOX_K_CANTOR = (2, 3, 4, 5, 6, 7, 8, 9, 10)   # Cantor：ε = 2^−k（base = 2，不和 1/3 對齊）


def _sample_koch(gen=9, per=2):
    p = _koch_unit(gen)
    t = (np.arange(per) + 0.5) / per
    a, b = p[:-1], p[1:]
    return (a[:, None, :] + t[None, :, None] * (b - a)[:, None, :]).reshape(-1, 2)


def box_counts(kind="koch", ks=BOX_K, base=3):
    """[(ε, 盒子數)]：邊長 ε = base^−k 的方格（Cantor 用線段）裡有曲線的個數（格線位移 BOX_OFFSET）。"""
    if kind == "koch":
        pts = _sample_koch()
        out = []
        for k in ks:
            e = float(base) ** -k
            idx = np.floor((pts - np.array(BOX_OFFSET)) / e).astype(np.int64)
            out.append((e, len(np.unique(idx[:, 0] * 10 ** 7 + idx[:, 1]))))
        return out
    pts = _cantor_mid(14)
    out = []
    for k in ks:
        e = float(base) ** -k
        out.append((e, len(np.unique(np.floor((pts - BOX_OFFSET[0]) / e)))))
    return out


def box_dim(kind="koch", ks=None, base=None):
    """盒子計數的維度：ln N(ε) 對 ln(1/ε) 的斜率（第二種算法；Koch 1.2619、Cantor 0.6309 附近）。"""
    if ks is None:
        ks, base = (BOX_K, 3) if kind == "koch" else (BOX_K_CANTOR, 2)
    bc = box_counts(kind, ks, base)
    return _slope([1 / e for e, _ in bc], [n for _, n in bc])


REPLICAS = {"line": (3, 3), "square": (9, 3), "koch": (4, 3), "cantor": (2, 3), "sierpinski": (3, 2)}


def replica_dim(nrep, mag):
    """(11.6)：p = ln ℛ / ln ℳ（ℛ 個複本，每個放大 ℳ 倍變回原物）。"""
    return math.log(nrep) / math.log(mag)


RW_SEED = 1101
RW_NS = (64, 256, 1024, 4096)
RW_M = 400


def rw_paths(Nmax=RW_NS[-1], M=RW_M, seed=RW_SEED):
    """M 條二維正方格子普通漫步（每步 E N W S 各 1/4），形狀 (M, Nmax+1, 2)。"""
    rng = np.random.default_rng(seed)
    st = np.array(SW.STEPS2)[rng.integers(0, 4, size=(M, Nmax))]
    out = np.zeros((M, Nmax + 1, 2), np.int64)
    out[:, 1:] = np.cumsum(st, axis=1)
    return out


@lru_cache(maxsize=None)
def rw_rg2(N, M=RW_M, seed=RW_SEED):
    """前 N 步（N+1 個點）的迴轉半徑平方 R_g² 的 M 條平均（固定種子）。"""
    P = rw_paths(max(RW_NS), M, seed)[:, :N + 1].astype(float)
    c = P.mean(axis=1, keepdims=True)
    return float(((P - c) ** 2).sum(axis=2).mean())


def rg2_exact(N):
    """⟨R_g²⟩ = N(N+2)/(6(N+1))（每步長 1，任何維度；Σ_{i<j}|i−j| / (N+1)²）。"""
    return F(N * (N + 2), 6 * (N + 1))


def rw_df(Ns=RW_NS, exact=False):
    """普通漫步的質量－半徑：ln N 對 ln R_g 的斜率（→ 2）。exact=True 用 rg2_exact。"""
    R = [math.sqrt(float(rg2_exact(N))) if exact else math.sqrt(rw_rg2(N)) for N in Ns]
    return _slope(R, Ns)


def rw_trail(N=1500, seed=RW_SEED):
    """畫圖用的一條普通漫步（N+1 個格點座標）。"""
    return rw_paths(N, 1, seed)[0]


def saw_df(d):
    """自迴避漫步的碎形維度 1/ν（d = 1、2、4 是 Fraction；三維 ν = 0.588 只引用）。"""
    nu = WC.NU_KNOWN[d]
    return 1 / nu if isinstance(nu, F) else 1.0 / nu


SAW_NEFF = (10, 12, 14)


def saw_df_eff(N):
    """有限 N 的有效維度 2/(2ν_eff)（walk_saw.two_nu_eff，二維精確列舉）。"""
    return 2.0 / SW.two_nu_eff(N)


# ======================================================================
# §11.1.2 「被動」的重整化：Landau 與 Gaussian 模型  [W1102]
# ======================================================================
def landau_m(r, h, u=1.0):
    """F/V = r m²/2 + u m⁴/4 − h m 的平衡磁化（h ≥ 0；u m³ + r m = h 的實根中 F 最小者）。"""
    roots = np.roots([u, 0.0, r, -h])
    re = [float(z.real) for z in roots if abs(z.imag) < 1e-9]
    Fv = [r * m * m / 2 + u * m ** 4 / 4 - h * m for m in re]
    return max(m for m, f in zip(re, Fv) if f - min(Fv) < 1e-14)


def landau_m_bisect(r, h, u=1.0):
    """同上，用二分法（第二種算法；h > 0 時根唯一的正根，或 r < 0、h = 0 時的正根）。"""
    if h == 0 and r >= 0:
        return 0.0
    lo, hi = 1e-300 if h == 0 else 0.0, 10.0 + abs(r) + abs(h)
    g = lambda m: u * m ** 3 + r * m - h
    if h == 0:
        lo = math.sqrt(-r / u) * 0.5
    for _ in range(200):
        mid = (lo + hi) / 2
        if g(mid) > 0:
            hi = mid
        else:
            lo = mid
    return (lo + hi) / 2


def landau_scaled(x, sign=1):
    """標度函數 f(x)：y = m/|r|^{1/2}、x = h/|r|^{3/2}（u = 1）滿足 sign·y + y³ = x（(11.16)）。"""
    return landau_m(float(sign), float(x))


COLLAPSE_R = (0.25, 0.5, 1.0)
COLLAPSE_HMAX = 2.0
COLLAPSE_POINTS = ((0.25, 0.25), (0.5, 2 ** -0.5), (1.0, 2.0))   # 都是 x = 2 → y = 1


def landau_curve(r, n=40, hmax=COLLAPSE_HMAX):
    """畫圖用 (h, m(r, h))，n 點。"""
    hs = np.linspace(0, hmax, n)
    return [(float(h), landau_m(r, float(h))) for h in hs]


def collapse_curve(r, n=40, hmax=COLLAPSE_HMAX):
    """同上換成 (h/r^{3/2}, m/r^{1/2})（r > 0）。"""
    return [(h / r ** 1.5, m / r ** 0.5) for h, m in landau_curve(r, n, hmax)]


def mf_rescale_check(r=0.3, h=0.2, d=3, bs=(0.5, 2.0, 3.0)):
    """(11.10)–(11.11)：m(r, h) = b^{−d/4} m(r b^{d/2}, h b^{3d/4})（u 不變）的最大誤差。"""
    m0 = landau_m(r, h)
    return max(abs(b ** (-d / 4) * landau_m(r * b ** (d / 2), h * b ** (3 * d / 4)) - m0) for b in bs)


def book1116_limit(r, h):
    """p. 263：書上說 f(h/r^{3/2}) → h^{1/3} r；應為 → h^{1/3} r^{−1/2}（即 f(x) ≈ x^{1/3}）。
    回傳 (f(h/r^{3/2}), 書上 h^{1/3} r, 更正 h^{1/3} r^{−1/2})。"""
    return landau_scaled(h / r ** 1.5), h ** (1 / 3) * r, h ** (1 / 3) * r ** -0.5


def gauss_m(r, h):
    """Gaussian 模型（只有 q = 0 的模態受外場）：m = h/r（r > 0）。"""
    return h / r


def gauss_rescale_check(d, r=0.3, h=0.2, bs=(0.5, 2.0, 3.0)):
    """(11.22)：m = b^{(2−d)/2} m(r b², h b^{(d+2)/2}) 的最大誤差（m = h/r）；(11.24) 的形狀 r^{(d−2)/4}F(h/r^{(d+2)/4}) 也核對。"""
    m0 = gauss_m(r, h)
    e1 = max(abs(b ** ((2 - d) / 2) * gauss_m(r * b * b, h * b ** ((d + 2) / 2)) - m0) for b in bs)
    e2 = abs(r ** ((d - 2) / 4) * (h / r ** ((d + 2) / 4)) - m0)          # F(x) = x
    return max(e1, e2)


def gauss_exponents(d):
    """Ex. 11.4：由 (11.24) 讀出 β = (d−2)/4、γ = 1、δ = (d+2)/(d−2)、α = 2 − d/2、ν = 1/2（Fraction）。"""
    d = F(d)
    return dict(beta=(d - 2) / 4, gamma=F(1), delta=(d + 2) / (d - 2), alpha=2 - d / 2, nu=F(1, 2))


def scaling_combos(d):
    """(Landau 的 h/r^{3/2} 指數, Gaussian 的 h/r^{(d+2)/4} 指數, Landau 前因子 1/2, Gaussian 前因子 (d−2)/4)：d = 4 時兩兩相同。"""
    d = F(d)
    return F(3, 2), (d + 2) / 4, F(1, 2), (d - 2) / 4


# ======================================================================
# §11.2 一維 Gaussian 鏈的抽格  [W1103]
# ======================================================================
A_BOOK = 1.0                                  # 書上的 a（Fig. 11.6 的座標是 r/a）


def square_check(r=0.5, a=1.0, S2=0.3, S4=-0.7, h=0.0):
    """(11.28)→(11.30)：∫dS3 exp(−r/2 S3² − a/2(S3−S2)² − a/2(S3−S4)² + h S3) 的 (數值積分, 完成平方的閉式)。"""
    f = lambda s: math.exp(-r / 2 * s * s - a / 2 * (s - S2) ** 2 - a / 2 * (s - S4) ** 2 + h * s)
    q = integrate.quad(f, -np.inf, np.inf, epsabs=1e-13, epsrel=1e-12)[0]
    A = r + 2 * a
    B = a * (S2 + S4) + h
    closed = math.sqrt(2 * math.pi / A) * math.exp(B * B / (2 * A) - a / 2 * (S2 * S2 + S4 * S4))
    return q, closed


def square_book_form(r=0.5, a=1.0, S2=0.3, S4=-0.7):
    """(11.30) 照書上的寫法：√(2π/(r+2a)) exp[−a²(S2−S4)²/(2(r+2a)) − ar(S2²+S4²)/(2(r+2a))]（對的）。"""
    A = r + 2 * a
    return math.sqrt(2 * math.pi / A) * math.exp(-a * a * (S2 - S4) ** 2 / (2 * A) - a * r * (S2 * S2 + S4 * S4) / (2 * A))


def decimate_gauss(r, a=A_BOOK, h=0.0):
    """抽掉奇數格點、還沒重新標度：(r', a', h')＝(r(r+4a)/(r+2a), a²/(r+2a), h(r+4a)/(r+2a))（(11.32)(11.33)）。"""
    A = r + 2 * a
    return r * (r + 4 * a) / A, a * a / A, h * (r + 4 * a) / A


def spin_factor(r, a=A_BOOK):
    """(11.34) 的自旋因子 √((r+2a)/a)（舊自旋 = 因子 × 新自旋；小 r → √2）。"""
    return math.sqrt((r + 2 * a) / a)


def rr_gauss(r, a=A_BOOK):
    """(11.36)：r'' = r(r+4a)/a（a'' = a）。"""
    return r * (r + 4 * a) / a


def hh_gauss(r, h, a=A_BOOK):
    """(11.38)：h'' = h(r+4a)/√(a(r+2a))（小 r → 2^{3/2} h）。"""
    return h * (r + 4 * a) / math.sqrt(a * (r + 2 * a))


def gauss_rescaled(r, a=A_BOOK, h=0.0):
    """decimate_gauss 再乘自旋因子的平方 / 一次方（第二種算法）→ (r'', a'', h'')。"""
    rp, ap, hp = decimate_gauss(r, a, h)
    c = spin_factor(r, a)
    return rp * c * c, ap * c * c, hp * c


def ring_matrix(N, r, a=A_BOOK):
    M = np.zeros((N, N))
    for i in range(N):
        M[i, i] += r + 2 * a
        M[i, (i + 1) % N] -= a
        M[(i + 1) % N, i] -= a
    return M


def ring_lnZ(N, r, a=A_BOOK, h=0.0):
    """N 個自旋的 Gaussian 環（週期，N ≥ 3）：ln Z = (N/2)ln 2π − ½ ln det M + ½ h²·1ᵀM⁻¹1。"""
    M = ring_matrix(N, r, a)
    one = np.ones(N)
    return N / 2 * math.log(2 * math.pi) - 0.5 * np.linalg.slogdet(M)[1] + 0.5 * h * h * one @ np.linalg.solve(M, one)


def ring_decimation_check(N=16, r=0.5, a=A_BOOK, h=0.2, book=False):
    """ln Z_N(r, a, h) ＝ (N/2)[½ ln(2π/(r+2a)) + h²/(2(r+2a))] + ln Z_{N/2}(r', a', h') 的誤差（精確，不必小 r）。
    N/2 個奇數格點各給 √(2π/(r+2a))，所以前因子是 (2π/(r+2a))^{N/4}；book=True 照書上 (11.31) 的 ^{N/2}
    （誤差 = (N/4) ln(2π/(r+2a))）。"""
    rp, ap, hp = decimate_gauss(r, a, h)
    A = r + 2 * a
    pw = 1.0 if book else 0.5
    rhs = N / 2 * (pw * math.log(2 * math.pi / A) + h * h / (2 * A)) + ring_lnZ(N // 2, rp, ap, hp)
    return abs(ring_lnZ(N, r, a, h) - rhs)


def gauss_m_invariance(r, h, a=A_BOOK):
    """m = h/r 與「自旋因子 × h''/r''」（一次抽格＋重新標度）相同：回傳差。"""
    return abs(gauss_m(r, h) - spin_factor(r, a) * gauss_m(rr_gauss(r, a), hh_gauss(r, h, a)))


def book1147_m(r, h, n):
    """照書上 (11.47) m = 2^{−n/2} f(r 2^{2n}, h 2^{3n/2})（f 用精確的 h/r）；正確是 2^{+n/2}。回傳 (書上, 更正, 精確 h/r)。"""
    f = gauss_m(r * 4 ** n, h * 2 ** (1.5 * n))
    return 2 ** (-n / 2) * f, 2 ** (n / 2) * f, gauss_m(r, h)


def gauss_iterate(r0, n, a=A_BOOK):
    """[r^(0), …, r^(n)]（(11.39)）。"""
    out = [r0]
    for _ in range(n):
        out.append(rr_gauss(out[-1], a))
    return out


COBWEB_R0 = 0.05
COBWEB_N = 3


def cobweb_points(r0=COBWEB_R0, n=COBWEB_N, a=A_BOOK):
    """Fig. 11.7 的蛛網折線：(r0, 0) → (r0, r1) → (r1, r1) → (r1, r2) → …（2n 個轉角＋起點）。"""
    rs = gauss_iterate(r0, n, a)
    pts = [(rs[0], 0.0)]
    for k in range(n):
        pts.append((rs[k], rs[k + 1]))
        if k < n - 1:
            pts.append((rs[k + 1], rs[k + 1]))
    return pts


def rr_curve(n=40, xmax=None, a=A_BOOK, ymax=10.0):
    """Fig. 11.6 的曲線 (r, r'')，r 從 0 到 r''= ymax 的地方（a = 1：r = −2 + √14 = 1.7417）。"""
    xmax = (-2 * a + math.sqrt(4 * a * a + a * ymax)) if xmax is None else xmax
    return [(float(x), rr_gauss(float(x), a)) for x in np.linspace(0, xmax, n)]


SMALL_R = 1e-4
SMALL_H = 1e-4


def small_r_table(n=5, r=SMALL_R, h=SMALL_H, a=A_BOOK):
    """[(k, r^(k)/(4^k r), h^(k)/(2^{3k/2} h), 自旋總因子/2^{k/2})]（(11.41)–(11.45)；小 r 時都 ≈ 1）。"""
    out = []
    rr, hh, sf = r, h, 1.0
    for k in range(1, n + 1):
        sf *= spin_factor(rr, a)
        hh = hh_gauss(rr, hh, a)
        rr = rr_gauss(rr, a)
        out.append((k, rr / (4 ** k * r), hh / (2 ** (1.5 * k) * h), sf / 2 ** (k / 2)))
    return out


def gauss_fixed_points(a=A_BOOK):
    """r'' = r 的根：0（不穩定，斜率 4 = 2^{Λ_r}）與 −3a（r + 2a < 0，積分發散，不物理）。"""
    return 0.0, -3.0 * a


def chain_lambdas():
    """抽格 b = 2 的特徵值 → Λ = ln(特徵值)/ln b：r 的 4 → 2、h 的 2^{3/2} → 3/2、自旋 √2 → 1/2 = (2−d)/2（d = 1）。"""
    eps = 1e-7
    er = (rr_gauss(eps) - rr_gauss(-eps)) / (2 * eps)
    eh = hh_gauss(0.0, 1.0)
    es = spin_factor(0.0)
    return tuple(math.log(e) / math.log(2) for e in (er, eh, es))


# ======================================================================
# §11.2 例題：一維 Ising 的抽格  [W1104]
# ======================================================================
def ising_A(K):
    """(11.54)：A = 2√cosh 2K。"""
    return 2 * math.sqrt(math.cosh(2 * K))


def ising_Kp(K):
    """(11.55)：K' = ½ ln cosh 2K。"""
    return 0.5 * math.log(math.cosh(2 * K))


def ising_pair_check(K):
    """(11.52)–(11.57)：Σ_{S3} e^{K S3(S2+S4)} 與 A e^{K' S2 S4} 在四種 (S2, S4) 的最大差。"""
    A, Kp = ising_A(K), ising_Kp(K)
    err = 0.0
    for s2 in (1, -1):
        for s4 in (1, -1):
            lhs = sum(math.exp(K * s3 * (s2 + s4)) for s3 in (1, -1))
            err = max(err, abs(lhs - A * math.exp(Kp * s2 * s4)))
    return err


def book1153(K):
    """(11.53) 照書上：左邊 cosh K(S2+S4)（同號時 cosh 2K），右邊 2cosh 2K → (左邊, 書上的右邊, 正確左邊 2cosh K(S2+S4))。"""
    return math.cosh(2 * K), 2 * math.cosh(2 * K), 2 * math.cosh(2 * K)


def ex118(K, H):
    """Ex. 11.8：Σ_{S3} e^{S3(K(S2+S4)+H)} = A e^{K1' S2S4 + K2'(S2+S4)}
    → (A, K1', K2', H' = H + 2K2')（H' 是偶數格點的新外場：自己的 H 加兩邊各 K2'）。"""
    cp, cm, c0 = math.cosh(2 * K + H), math.cosh(2 * K - H), math.cosh(H)
    K2 = 0.25 * math.log(cp / cm)
    K1 = 0.25 * math.log(cp * cm / c0 ** 2)
    A = 2 * (cp * cm * c0 ** 2) ** 0.25
    return A, K1, K2, H + 2 * K2


def ex118_check(K, H):
    """Ex. 11.8 的四種 (S2, S4) 直接加總 vs 公式的最大差。"""
    A, K1, K2, _ = ex118(K, H)
    err = 0.0
    for s2 in (1, -1):
        for s4 in (1, -1):
            lhs = sum(math.exp(s3 * (K * (s2 + s4) + H)) for s3 in (1, -1))
            err = max(err, abs(lhs - A * math.exp(K1 * s2 * s4 + K2 * (s2 + s4))))
    return err


ISING_K0 = 1.0
ISING_FIELD = 0.1


def ising_iterate(K0=ISING_K0, n=5):
    """(11.62)：K_{n+1} = ½ ln cosh 2K_n 的 [K_1 = K0, K_2, …]（→ 0）。"""
    out = [K0]
    for _ in range(n):
        out.append(ising_Kp(out[-1]))
    return out


def ring_Z_brute(N, K, H=0.0):
    """N 個 ±1 自旋的環：逐一列舉 2^N 個組態（N ≤ 16）。"""
    s = 1 - 2 * ((np.arange(2 ** N)[:, None] >> np.arange(N)[None, :]) & 1)
    bond = (s * np.roll(s, -1, axis=1)).sum(1)
    return float(np.exp(K * bond + H * s.sum(1)).sum())


def ring_Z_transfer(N, K):
    """轉移矩陣：Z = (2cosh K)^N + (2sinh K)^N。"""
    return (2 * math.cosh(K)) ** N + (2 * math.sinh(K)) ** N


def ring_Z_decimation(N, K, book=False):
    """(11.63)：N = 2^k 的環一直抽格到剩 1 個自旋：Z = Π_{j=1}^{k} A(K_j)^{N/2^j} · 2e^{K_{k+1}}
    （最後 1 個自旋自己和自己的鍵 K S1S1 = K）。book=True 照書上指數 N/2^{j+1}。"""
    k = int(round(math.log2(N)))
    assert 2 ** k == N
    Ks = ising_iterate(K, k)
    lnZ = sum((N / 2 ** (j + 1 if book else j)) * math.log(ising_A(Ks[j - 1])) for j in range(1, k + 1))
    return math.exp(lnZ + math.log(2) + Ks[k])


RING_NS = (4, 8, 16)


def x_of_K(K):
    """x = e^{−2K}（Ex. 11.9）。"""
    return math.exp(-2 * K)


def x_map(x):
    """Ex. 11.9：x' = 2x/(1+x²)（由 e^{2K'} = cosh 2K）。"""
    return 2 * x / (1 + x * x)


def x_iterate(x0, n):
    out = [x0]
    for _ in range(n):
        out.append(x_map(out[-1]))
    return out


def x_fixed_points():
    """x' = x 的根：x = 0（T = 0，不穩定，斜率 2）與 x = 1（T = ∞，穩定，斜率 0）（x = −1 不物理）。"""
    return 0.0, 1.0


def x_slope(x):
    """dx'/dx = 2(1−x²)/(1+x²)²：x = 0 → 2、x = 1 → 0。"""
    return 2 * (1 - x * x) / (1 + x * x) ** 2


def xi_exact(K):
    """一維 Ising 的關聯長度 ξ = −1/ln tanh K（＝ walk_crit.ising1d_xi）。"""
    return -1 / math.log(math.tanh(K))


XI_K = (1.0, 2.0, 3.0)


def xi_ratio(K):
    """ξ / e^{2K}（→ ½；Ex. 11.9(b)：ξ ∝ 1/x = e^{2J/k_BT}）。"""
    return xi_exact(K) / math.exp(2 * K)


def tanh_rg_check(K):
    """tanh K' = tanh² K（精確）⇒ ξ(K') = ξ(K)/2：回傳 (|tanhK' − tanh²K|, ξ(K)/ξ(K') − 2)。"""
    Kp = ising_Kp(K)
    return abs(math.tanh(Kp) - math.tanh(K) ** 2), xi_exact(K) / xi_exact(Kp) - 2


FREE_K = (1.0, 2.0, 5.0)


def free_energy_ratio(K, N=None):
    """−F/(NJ) = ln Z/(NK)：N = None 用無限長 ln(2cosh K)/K；N = 2^k 用抽格的連乘（→ 1：F → −NJ）。"""
    if N is None:
        return math.log(2 * math.cosh(K)) / K
    k = int(round(math.log2(N)))
    Ks = ising_iterate(K, k)
    lnZ = sum((N / 2 ** j) * math.log(ising_A(Ks[j - 1])) for j in range(1, k + 1)) + math.log(2) + Ks[k]
    return lnZ / (N * K)


def free_energy_series(K, jmax=60):
    """無限長鏈：ln Z/N = Σ_j 2^{−j} ln A(K_j)（(11.63) 取 N → ∞；第二種算法）。"""
    Ks = ising_iterate(K, jmax)
    return sum(2.0 ** -j * math.log(ising_A(Ks[j - 1])) for j in range(1, jmax + 1))


# ======================================================================
# §11.2.1–11.2.2 動量殼與連續重整化  [W1105]
# ======================================================================
def alpha_q(q, a=A_BOOK, c=1.0):
    """一維 Gaussian 鏈的 Fourier 係數：(a/2)Σ(S_i − S_{i+1})² = Σ_q a(1 − cos qc)|S(q)|²；小 q：a c² q²/2。"""
    return a * (1 - math.cos(q * c))


def alpha_q_book(q, a=A_BOOK, c=1.0):
    """p. 273 書上 α(q) = −2a(1 − cos qc)（負號與 2 倍都不對）。"""
    return -2 * a * (1 - math.cos(q * c))


ALPHA_N = 8


def alpha_table(N=ALPHA_N, a=A_BOOK):
    """[(k, q = 2πk/N, a(1−cos q), aq²/2, 書上)]，k = 0..N/2。"""
    return [(k, 2 * math.pi * k / N, alpha_q(2 * math.pi * k / N, a), a * (2 * math.pi * k / N) ** 2 / 2,
             alpha_q_book(2 * math.pi * k / N, a)) for k in range(N // 2 + 1)]


def fourier_check(N=ALPHA_N, r=0.5, a=A_BOOK, seed=1105):
    """隨機實數自旋：Σ_i [r/2 S_i² + a/2(S_i−S_{i+1})²] ＝ Σ_q [r/2 + a(1−cos q)]|S(q)|²（S(q) = N^{−1/2}Σ S_i e^{−iqi}）的相對誤差；
    以及書上的 α 得到的值。回傳 (誤差, 實空間值, 書上 α 的值)。"""
    S = np.random.default_rng(seed).normal(size=N)
    real = float(np.sum(r / 2 * S * S + a / 2 * (S - np.roll(S, -1)) ** 2))
    Sq = np.fft.fft(S) / math.sqrt(N)
    q = 2 * np.pi * np.arange(N) / N
    mom = float(np.sum((r / 2 + a * (1 - np.cos(q))) * np.abs(Sq) ** 2))
    book = float(np.sum((r / 2 - 2 * a * (1 - np.cos(q))) * np.abs(Sq) ** 2))
    return abs(mom - real) / real, real, book


def shell_fraction(d, b=2):
    """布里淵區每邊縮成 1/b：拿掉的模態比例 1 − b^{−d}（d = 2：3/4）。"""
    return 1 - F(1, b ** d)


def shell_rg(r, h, d, b=2):
    """(11.71)(11.72) 與 b 的推廣：(b² r, b^{(d+2)/2} h)。"""
    return b * b * r, b ** ((d + 2) / 2) * h


def shell_rg_steps(r, h, d, b=2):
    """一步一步：q → q/b（Jq² → Jq²/b²）、N → b^d N、S → bS（保持 J）→ (r', h')（第二種算法）。"""
    Jcoef = 1.0 / b ** 2                      # 換 q 之後的 J 係數
    s = 1.0 / math.sqrt(Jcoef)                # 自旋因子 b
    return r * s * s, h * math.sqrt(b ** d) * s


FLOW_R0 = 0.01
FLOW_H0 = 0.01
FLOW_ELL = (1.0, 2.0)


def flow(r0, h0, d, ell):
    """(11.82)(11.83) 的解：r(ℓ) = r0 e^{2ℓ}、h(ℓ) = h0 e^{(d+2)ℓ/2}。"""
    return r0 * math.exp(2 * ell), h0 * math.exp((d + 2) * ell / 2)


def flow_numeric(r0, h0, d, ell):
    """同上，solve_ivp（RK45）數值積分（第二種算法）。"""
    sol = integrate.solve_ivp(lambda t, y: [2 * y[0], (d + 2) / 2 * y[1]], (0, ell), [r0, h0],
                              rtol=1e-11, atol=1e-14)
    return float(sol.y[0, -1]), float(sol.y[1, -1])


def flow_line(r0, h0, d, ell_max, n=40):
    """畫圖用流線 [(r(ℓ), h(ℓ))]，ℓ ∈ [0, ell_max]，n 點（h ∝ |r|^{(d+2)/4}）。"""
    return [flow(r0, h0, d, float(t)) for t in np.linspace(0, ell_max, n)]


def m24(r, h, d, Ffn=None):
    """(11.24)：m = r^{(d−2)/4} F(h/r^{(d+2)/4})（預設 F(x) = x/(1+x²)，任意可微的 F 都行）。"""
    Ffn = Ffn or (lambda x: x / (1 + x * x))
    return r ** ((d - 2) / 4) * Ffn(h / r ** ((d + 2) / 4))


def pde_residual(d, r=0.4, h=0.3, eps=1e-5, Ffn=None):
    """(11.87)：[−(d−2)/2 + 2r ∂_r + (d+2)/2 h ∂_h] m 用中央差分（m = (11.24)）；回傳殘差 ÷ |m|。"""
    m = lambda rr, hh: m24(rr, hh, d, Ffn)
    dr = (m(r + eps, h) - m(r - eps, h)) / (2 * eps)
    dh = (m(r, h + eps) - m(r, h - eps)) / (2 * eps)
    res = -(d - 2) / 2 * m(r, h) + 2 * r * dr + (d + 2) / 2 * h * dh
    return abs(res) / abs(m(r, h))


def char_check(d, r=0.4, h=0.3, ells=(0.5, 1.0, 2.0), Ffn=None):
    """特徵線：e^{(2−d)ℓ/2} m(r(ℓ), h(ℓ)) 與 ℓ 無關（(11.86)）：最大相對差。"""
    m0 = m24(r, h, d, Ffn)
    return max(abs(math.exp((2 - d) * t / 2) * m24(*flow(r, h, d, t), d, Ffn) - m0) / abs(m0) for t in ells)


# ======================================================================
# §11.3 不動點、線性化與臨界指數  [W1106]
# ======================================================================
TOY_A = ((F(1), F(2)), (F(2), F(-2)))         # 示意用的 2×2 線性化矩陣（對稱）


def toy_eigen(A=TOY_A):
    """2×2 Fraction 矩陣的特徵值與特徵向量（特徵多項式 λ² − tr λ + det，判別式要是完全平方）：
    TOY_A → [(2, (2, 1)), (−3, (1, −2))]：相關（Λ > 0）與無關（Λ < 0）。"""
    (a, b), (c, d) = A
    tr, det = a + d, a * d - b * c
    disc = tr * tr - 4 * det
    sq = F(math.isqrt(disc.numerator), math.isqrt(disc.denominator))
    assert sq * sq == disc
    out = []
    for lam in ((tr + sq) / 2, (tr - sq) / 2):
        v = (b, lam - a) if b != 0 else (lam - d, c)
        g = math.gcd(int(v[0]), int(v[1])) if all(F(x).denominator == 1 for x in v) else 1
        v = (F(v[0]) / g, F(v[1]) / g)
        if v[0] < 0:
            v = (-v[0], -v[1])
        out.append((lam, v))
    return out


def toy_flow(c1, c2, ell):
    """(11.95)：δκ(ℓ) = c1 e^{2ℓ}(2, 1) + c2 e^{−3ℓ}(1, −2)。"""
    (l1, v1), (l2, v2) = toy_eigen()
    return tuple(float(c1 * math.exp(float(l1) * ell) * v1[i] + c2 * math.exp(float(l2) * ell) * v2[i]) for i in range(2))


def toy_flow_expm(c1, c2, ell):
    """同上，用矩陣指數 e^{Aℓ} δκ(0)（第二種算法）。"""
    A = np.array([[float(x) for x in row] for row in TOY_A])
    (l1, v1), (l2, v2) = toy_eigen()
    k0 = c1 * np.array(v1, float) + c2 * np.array(v2, float)
    return tuple(float(x) for x in linalg.expm(A * ell) @ k0)


TOY_STARTS = ((0.0, 1.0), (0.05, 1.0), (-0.05, 1.0), (0.2, 0.5), (-0.2, 0.5))   # (c1, c2)


def toy_flow_line(c1, c2, ell_max=1.0, n=30):
    """畫圖用流線（n 點）；c1 = 0 的線留在臨界線（無關方向）上流進不動點。"""
    return [toy_flow(c1, c2, float(t)) for t in np.linspace(0, ell_max, n)]


def rg_exponents(Lr, Lh, d):
    """(11.105)–(11.108)＋ν = 1/Λ_r → dict(alpha, beta, gamma, delta, nu)（Fraction）。"""
    Lr, Lh, d = F(Lr), F(Lh), F(d)
    return dict(alpha=2 - d / Lr, beta=(d - Lh) / Lr, gamma=-(d - 2 * Lh) / Lr, delta=Lh / (d - Lh), nu=1 / Lr)


def rg_exponents_widom(Lr, Lh, d):
    """第二種算法：α = 2 − d/Λ_r、能隙 Δ = Λ_h/Λ_r 代入 walk_crit.widom_exponents（第 9 章 Ex. 9.1–9.2）。"""
    alpha = 2 - F(d) / F(Lr)
    beta, gamma, delta = WC.widom_exponents(alpha, F(Lh) / F(Lr))
    return dict(alpha=alpha, beta=beta, gamma=gamma, delta=delta, nu=1 / F(Lr))


def lambdas_from(nu, beta, d):
    """反過來：Λ_r = 1/ν、Λ_h = d − β/ν（第三種核對）。"""
    return 1 / F(nu), F(d) - F(beta) / F(nu)


RG_CASES = {
    "gauss3": (F(2), F(5, 2), 3),
    "gauss4": (F(2), F(3), 4),
    "ising2d": (F(1), F(15, 8), 2),
    "saw2d": (F(4, 3), F(91, 48), 2),
}


def gauss_lambdas(d):
    """(11.109)(11.110)：Λ_r = 2、Λ_h = (d+2)/2。"""
    return F(2), F(d + 2, 2)


XI_R = 0.01


def xi_from_flow(r, Lr, D=1.0, Kc=1.0):
    """(11.99)–(11.103)：K r e^{Λ_r ℓ*} = D ⇒ ξ = e^{ℓ*} = (K r/D)^{−1/Λ_r} → (ℓ*, ξ)。"""
    ell = math.log(D / (Kc * r)) / float(Lr)
    return ell, math.exp(ell)


def xi_from_flow_numeric(r, Lr, D=1.0, Kc=1.0):
    """同上，數值積分 dδ/dℓ = Λ_r δ（δ(0) = K r）直到 δ = D（事件）（第二種算法）。"""
    ev = lambda t, y: y[0] - D
    ev.terminal = True
    sol = integrate.solve_ivp(lambda t, y: [float(Lr) * y[0]], (0, 100), [Kc * r], events=ev, rtol=1e-11, atol=1e-14)
    ell = float(sol.t_events[0][0])
    return ell, math.exp(ell)


# ======================================================================
# summary / check
# ======================================================================
def summary():
    S = {}
    # W1101
    S["ln4_ln3"] = _f4(LN4_LN3)
    S["ln2_ln3"] = _f4(LN2_LN3)
    S["koch_counts"] = [(g, *map(str, koch_counts(g))) for g in range(6)]
    S["koch_len"] = [_f4(koch_counts(g)[1]) for g in range(6)]
    S["cantor_counts"] = [(g, *map(str, cantor_counts(g))) for g in range(6)]
    S["cantor_len"] = [_f4(cantor_counts(g)[1]) for g in range(6)]
    S["koch_mass"] = [(k, _f4(R), _f4(n)) for k, R, n in mass_table("koch")]
    S["cantor_mass"] = [(k, _f4(R), _f4(n)) for k, R, n in mass_table("cantor")]
    S["mass_dim"] = (_f4(mass_radius_dim("koch")), _f4(mass_radius_dim("cantor")))
    S["box_koch"] = box_counts("koch")
    S["box_cantor"] = box_counts("cantor", BOX_K_CANTOR, 2)
    S["box_dim"] = (_f4(box_dim("koch")), _f4(box_dim("cantor")))
    S["replica"] = {k: _f4(replica_dim(*v)) for k, v in REPLICAS.items()}
    S["rw_rg2"] = [(N, f"{rw_rg2(N):.2f}", f"{float(rg2_exact(N)):.2f}") for N in RW_NS]
    S["rw_df"] = (_f4(rw_df()), _f4(rw_df(exact=True)))
    S["saw_df"] = (str(saw_df(2)), f"{saw_df(3):.2f}", str(saw_df(4)))
    S["saw_df_eff"] = [(N, _f4(saw_df_eff(N))) for N in SAW_NEFF]
    # W1102
    S["landau_f"] = [(x, _f4(landau_scaled(x))) for x in (1, 2, 10)]
    S["landau_pts"] = [(r, _f4(h), _f4(landau_m(r, h)), _f4(landau_m(r, h) / r ** 0.5)) for r, h in COLLAPSE_POINTS]
    S["landau_m_r0"] = _f4(landau_m(-0.25, 0.0))
    S["landau_m_h"] = [(h, _f4(landau_m(0.0, h))) for h in (0.008, 0.125)]
    S["book1116"] = tuple(_f4(v) for v in book1116_limit(1e-4, 1e-3))
    S["gauss_exp3"] = {k: str(v) for k, v in gauss_exponents(3).items()}
    S["gauss_exp4"] = {k: str(v) for k, v in gauss_exponents(4).items()}
    # W1103
    q, c = square_check()
    S["square"] = (f"{q:.6f}", f"{c:.6f}")
    S["rr_curve_end"] = _f4(rr_curve()[-1][0])
    S["cobweb"] = [_f4(x) for x in gauss_iterate(COBWEB_R0, 4)]
    S["small_r"] = [(k, *(f"{v:.4f}" for v in t)) for k, *t in small_r_table()]
    S["chain_lambdas"] = tuple(_f4(v) for v in chain_lambdas())
    S["hh_small"] = _f4(hh_gauss(0.0, 1.0))
    S["book1147"] = tuple(_f4(v) for v in book1147_m(0.01, 0.001, 1))
    S["gauss_fixed"] = gauss_fixed_points()
    # W1104
    S["ising_A_Kp"] = (_f4(ising_A(1.0)), _f4(ising_Kp(1.0)))
    S["ising_K"] = [_f4(k) for k in ising_iterate(1.0, 5)]
    S["book1153"] = tuple(_f4(v) for v in book1153(1.0))
    S["ex118"] = tuple(_f4(v) for v in ex118(ISING_K0, ISING_FIELD))
    S["ringZ"] = [(N, _f4(ring_Z_brute(N, 1.0)), _f4(ring_Z_transfer(N, 1.0)), _f4(ring_Z_decimation(N, 1.0)),
                   _f4(ring_Z_decimation(N, 1.0, book=True))) for N in RING_NS]
    S["x_iter"] = [_f4(x) for x in x_iterate(x_of_K(1.0), 5)]
    S["xi"] = [(K, _f4(xi_exact(K)), _f4(math.exp(2 * K) / 2), _f4(xi_ratio(K))) for K in XI_K]
    S["free"] = [(K, _f4(free_energy_ratio(K))) for K in FREE_K]
    S["free_ring16"] = _f4(free_energy_ratio(1.0, 16))
    # W1105
    S["alpha"] = [(k, _f4(q), _f4(al), _f4(sm), _f4(bk)) for k, q, al, sm, bk in alpha_table()]
    fc = fourier_check()
    S["fourier"] = (f"{fc[0]:.1e}", _f4(fc[1]), _f4(fc[2]))
    S["shell_frac"] = (str(shell_fraction(2)), str(shell_fraction(3)))
    S["shell_h3"] = _f4(shell_rg(1.0, 1.0, 3)[1])
    S["flow"] = [(t, _f4(flow(FLOW_R0, FLOW_H0, 3, t)[0]), _f4(flow(FLOW_R0, FLOW_H0, 3, t)[1])) for t in FLOW_ELL]
    S["pde"] = [(d, f"{pde_residual(d):.1e}") for d in (1, 3, 4)]
    # W1106
    S["toy"] = [(str(l), tuple(map(str, v))) for l, v in toy_eigen()]
    S["rg"] = {k: {e: str(v) for e, v in rg_exponents(*c).items()} for k, c in RG_CASES.items()}
    S["xi_flow"] = [(str(Lr), _f4(xi_from_flow(XI_R, Lr)[0]), _f4(xi_from_flow(XI_R, Lr)[1])) for Lr in (F(2), F(1))]
    return S


def _close(a, b, tol=1e-10):
    return abs(a - b) <= tol * max(1.0, abs(a), abs(b))


def check():
    S = summary()
    # ---- W1101：Koch／Cantor 的段數與長度（頂點陣列 ＝ Fraction 公式）
    for g in range(7):
        n, Ltot, seg = koch_counts(g)
        p = koch_points(g)
        lens = np.hypot(*np.diff(p, axis=0).T)
        assert len(p) - 1 == n == 4 ** g and np.allclose(lens, float(seg)) and _close(lens.sum(), float(Ltot))
        assert np.allclose(p[0], (0, 0)) and np.allclose(p[-1], (1, 0))
        n, Ltot, seg = cantor_counts(g)
        iv = cantor_intervals(g)
        assert len(iv) == n and sum(b - a for a, b in iv) == Ltot and all(b - a == seg for a, b in iv)
    q = koch_points(3, (1.0, 2.0), (1.0, 5.0))                      # 任意端點：長度 ×3、起終點
    assert np.allclose(q[0], (1, 2)) and np.allclose(q[-1], (1, 5)) and _close(np.hypot(*np.diff(q, axis=0).T).sum(), 3 * 64 / 27)
    assert [s for s in S["koch_len"]] == ['1.0000', '1.3333', '1.7778', '2.3704', '3.1605', '4.2140']
    assert [s for s in S["cantor_len"]] == ['1.0000', '0.6667', '0.4444', '0.2963', '0.1975', '0.1317']
    # 質量－半徑（以端點為圓心）：n(3^−k) = 4^−k、2^−k 精確（Koch 的 30° 錐：第 2–4 段都在 R = 1/3 之外）
    for k, R, n in mass_table("koch"):
        assert _close(n, 4.0 ** -k, 1e-12)
    for k, R, n in mass_table("cantor"):
        assert _close(n, 2.0 ** -k, 1e-12)
    assert _close(mass_radius_dim("koch"), LN4_LN3, 1e-9) and _close(mass_radius_dim("cantor"), LN2_LN3, 1e-9)
    assert _f4(LN4_LN3) == "1.2619" and _f4(LN2_LN3) == "0.6309"
    assert _f4(mass_radius_dim("koch", [2.0 ** -j for j in range(9)])) == "1.2571"     # 非 3 的冪：對數週期的抖動
    # 盒子計數（第二種算法）與複本公式（第三種）
    assert [n for _, n in box_counts("koch", BOX_K, 3)] == [22, 83, 325, 1291] and _f4(box_dim("koch")) == "1.2362"
    assert _f4(box_dim("cantor")) == "0.6254"
    assert abs(box_dim("koch") - LN4_LN3) < 0.03 and abs(box_dim("cantor") - LN2_LN3) < 0.01
    rep = {k: _f4(replica_dim(*v)) for k, v in REPLICAS.items()}
    assert rep == {"line": "1.0000", "square": "2.0000", "koch": "1.2619", "cantor": "0.6309", "sierpinski": "1.5850"}
    assert _close(replica_dim(4, 3), LN4_LN3) and _close(replica_dim(2, 3), LN2_LN3)
    # 普通漫步：⟨R_g²⟩ 公式 ＝ 4^N 全部列舉（N ≤ 5）；固定種子的斜率
    for N in range(1, 6):
        x, y = (np.asarray(v, float) for v in SW.all_walks_xy(N))
        rg2 = ((x - x.mean(1, keepdims=True)) ** 2 + (y - y.mean(1, keepdims=True)) ** 2).mean(1).mean()
        assert _close(float(rg2), float(rg2_exact(N)), 1e-12)
    assert rg2_exact(4) == F(4, 5)
    for N in RW_NS:
        assert abs(rw_rg2(N) / float(rg2_exact(N)) - 1) < 0.06
    assert _f4(rw_df()) == "2.0258" and _f4(rw_df(exact=True)) == "2.0069"
    assert saw_df(2) == F(4, 3) and f"{saw_df(3):.2f}" == "1.70" and saw_df(4) == 2 and saw_df(1) == 1
    assert [_f4(saw_df_eff(N)) for N in SAW_NEFF] == ["1.3848", "1.3788", "1.3740"]
    assert _close(2 / saw_df_eff(14), SW.two_nu_eff(14)) and f"{SW.two_nu_eff(14):.4f}" == "1.4556"
    tr = rw_trail()
    assert tr.shape == (1501, 2) and np.all(np.abs(np.diff(tr, axis=0)).sum(1) == 1)

    # ---- W1102：Landau（三次方程 ＝ 二分法 ＝ walk_crit.landau_min）與塌縮
    for r, h in ((0.25, 0.25), (0.5, 0.3), (-0.25, 0.0), (-0.3, 0.1), (0.0, 0.125), (1.0, 2.0)):
        m1 = landau_m(r, h)
        assert _close(m1, landau_m_bisect(r, h), 1e-10) and _close(m1, WC.landau_min(r, 1.0, h), 1e-10)
    assert _close(landau_scaled(2), 1.0) and _close(landau_scaled(10), 2.0)
    cub = (0.5 + math.sqrt(0.25 + 1 / 27)) ** (1 / 3) - (math.sqrt(0.25 + 1 / 27) - 0.5) ** (1 / 3)   # y³ + y = 1（Cardano）
    assert _close(landau_scaled(1), cub) and _f4(cub) == "0.6823"
    for r in COLLAPSE_R:
        for (h, m), (x, y) in zip(landau_curve(r), collapse_curve(r)):
            assert _close(m / r ** 0.5, landau_scaled(h / r ** 1.5), 1e-9) and _close(y, landau_scaled(x), 1e-9)
    for r, h in COLLAPSE_POINTS:
        assert _close(h / r ** 1.5, 2.0) and _close(landau_m(r, h), r ** 0.5)
    assert [_f4(landau_m(r, h)) for r, h in COLLAPSE_POINTS] == ["0.5000", "0.7071", "1.0000"]
    assert _f4(landau_m(-0.25, 0.0)) == "0.5000" and _f4(landau_m(0.0, 0.008)) == "0.2000" and _f4(landau_m(0.0, 0.125)) == "0.5000"
    assert mf_rescale_check() < 1e-10 and mf_rescale_check(-0.2, 0.05, 4) < 1e-10
    fx, bk, ok = book1116_limit(1e-4, 1e-3)
    assert _f4(fx) == "9.9667" and _f4(ok) == "10.0000" and bk == 1e-3 ** (1 / 3) * 1e-4 and _f4(bk) == "0.0000"
    for d in (1, 2, 3, 4, 5):
        assert gauss_rescale_check(d) < 1e-12
    assert gauss_exponents(3) == dict(beta=F(1, 4), gamma=F(1), delta=F(5), alpha=F(1, 2), nu=F(1, 2))
    assert gauss_exponents(4) == dict(beta=F(1, 2), gamma=F(1), delta=F(3), alpha=F(0), nu=F(1, 2))
    a3, g3, p3, q3 = scaling_combos(3)
    a4, g4, p4, q4 = scaling_combos(4)
    assert (a3, g3, p3, q3) == (F(3, 2), F(5, 4), F(1, 2), F(1, 4)) and a4 == g4 and p4 == q4
    for d in (3, 4, 5):                                             # Ex. 11.4 ＝ (11.111)–(11.114)
        e = rg_exponents(*gauss_lambdas(d), d)
        assert all(e[k] == gauss_exponents(d)[k] for k in e)

    # ---- W1103：一維 Gaussian 鏈
    for (r, a, S2, S4, h) in ((0.5, 1.0, 0.3, -0.7, 0.0), (0.1, 2.0, 1.2, 0.4, 0.0), (0.5, 1.0, 0.3, -0.7, 0.4)):
        qv, cv = square_check(r, a, S2, S4, h)
        assert _close(qv, cv, 1e-9)
        if h == 0:
            assert _close(cv, square_book_form(r, a, S2, S4), 1e-12)
    assert S["square"] == ("1.224819", "1.224819")
    for r in (0.05, 0.5, 3.0):
        for a in (1.0, 2.5):
            r2, a2, h2 = gauss_rescaled(r, a, 0.3)
            assert _close(r2, rr_gauss(r, a)) and _close(a2, a) and _close(h2, hh_gauss(r, 0.3, a))
            assert gauss_m_invariance(r, 0.3, a) < 1e-12
            for N in (8, 16, 32):
                assert ring_decimation_check(N, r, a, 0.2) < 1e-9
                assert _close(ring_decimation_check(N, r, a, 0.2, book=True), abs(N / 4 * math.log(2 * math.pi / (r + 2 * a))), 1e-8)
    assert _f4(ring_decimation_check(8, 0.5, 1.0, 0.0, book=True)) == "1.8432"
    assert _f4(rr_curve()[-1][0]) == "1.7417" and _close(rr_curve()[-1][0], -2 + math.sqrt(14)) and _close(rr_curve()[-1][1], 10.0)
    cw = gauss_iterate(COBWEB_R0, 4)
    assert [_f4(x) for x in cw] == ["0.0500", "0.2025", "0.8510", "4.1282", "33.5553"]
    cp = cobweb_points()
    assert len(cp) == 2 * COBWEB_N and cp[0] == (COBWEB_R0, 0.0) and _close(cp[-1][1], cw[COBWEB_N])
    st = small_r_table()
    assert [tuple(f"{v:.4f}" for v in t[1:]) for t in st] == [("1.0000", "1.0000", "1.0000"), ("1.0001", "1.0000", "1.0001"),
                                                            ("1.0005", "1.0000", "1.0005"), ("1.0021", "1.0000", "1.0021"),
                                                            ("1.0086", "1.0000", "1.0085")]
    assert [_f4(v) for v in chain_lambdas()] == ["2.0000", "1.5000", "0.5000"]
    assert gauss_lambdas(1) == (F(2), F(3, 2)) and _f4(hh_gauss(0.0, 1.0)) == "2.8284" and _close(hh_gauss(0.0, 1.0), 2 ** 1.5)
    for x in gauss_fixed_points():
        assert _close(rr_gauss(x), x)
    bk, ok, ex = book1147_m(0.01, 0.001, 1)
    assert (_f4(bk), _f4(ok), _f4(ex)) == ("0.0500", "0.1000", "0.1000")
    for n in (1, 2, 3):
        bk, ok, ex = book1147_m(0.01, 0.001, n)
        assert _close(ok, ex) and _close(bk, ex * 2.0 ** -n)

    # ---- W1104：一維 Ising
    for K in (0.3, 1.0, 2.0):
        assert ising_pair_check(K) < 1e-12 and ex118_check(K, 0.1) < 1e-12 and ex118_check(K, -0.7) < 1e-12
        A0, K10, K20, H0 = ex118(K, 0.0)
        assert _close(A0, ising_A(K)) and _close(K10, ising_Kp(K)) and abs(K20) < 1e-15 and abs(H0) < 1e-15
        assert _close(x_map(x_of_K(K)), x_of_K(ising_Kp(K)))
        e1, e2 = tanh_rg_check(K)
        assert e1 < 1e-14 and abs(e2) < 1e-9
        assert _close(xi_exact(K), WC.ising1d_xi(K))
    assert (_f4(ising_A(1.0)), _f4(ising_Kp(1.0))) == ("3.8793", "0.6625")
    assert [_f4(k) for k in ising_iterate(1.0, 5)] == ["1.0000", "0.6625", "0.3501", "0.1137", "0.0128", "0.0002"]
    assert tuple(_f4(v) for v in book1153(1.0)) == ("3.7622", "7.5244", "7.5244")
    assert tuple(_f4(v) for v in ex118(ISING_K0, ISING_FIELD)) == ("3.8897", "0.6602", "0.0482", "0.1964")
    for K, H in ((1.0, 0.1), (0.5, -0.3)):                          # 加外場的環：Z_8(K, H) = A⁴ Z_4(K1', H')
        A, K1, K2, Hp = ex118(K, H)
        assert _close(ring_Z_brute(8, K, H), A ** 4 * ring_Z_brute(4, K1, Hp), 1e-11)
    for N in RING_NS:
        for K in (0.4, 1.0):
            zb = ring_Z_brute(N, K)
            assert _close(zb, ring_Z_transfer(N, K), 1e-11) and _close(zb, ring_Z_decimation(N, K), 1e-11)
    assert _close(ring_Z_brute(4, 1.0), WC.ising_exact((4, WC.ring_edges(4)), 1.0)["Z"], 1e-12)
    assert [(_f4(ring_Z_brute(N, 1.0)), _f4(ring_Z_decimation(N, 1.0, book=True))) for N in (4, 8)] == \
        [("121.2329", "18.5498"), ("9160.4387", "143.2705")]
    assert [_f4(x) for x in x_iterate(x_of_K(1.0), 5)] == ["0.1353", "0.2658", "0.4965", "0.7966", "0.9747", "0.9997"]
    for x in x_fixed_points():
        assert _close(x_map(x), x)
        eps = 1e-6
        assert _close((x_map(x + eps) - x_map(x - eps)) / (2 * eps), x_slope(x), 1e-6) or abs(x_slope(x)) < 1e-9
    assert x_slope(0.0) == 2.0 and x_slope(1.0) == 0.0
    assert [(_f4(xi_exact(K)), _f4(math.exp(2 * K) / 2), _f4(xi_ratio(K))) for K in XI_K] == \
        [("3.6719", "3.6945", "0.4969"), ("27.2960", "27.2991", "0.4999"), ("201.7140", "201.7144", "0.5000")]
    assert [_f4(free_energy_ratio(K)) for K in FREE_K] == ["1.1269", "1.0091", "1.0000"]
    for K in FREE_K:
        assert _close(free_energy_series(K), K * free_energy_ratio(K), 1e-12)
    assert _close(free_energy_ratio(1.0, 16), math.log(ring_Z_brute(16, 1.0)) / 16, 1e-12) and _f4(free_energy_ratio(1.0, 16)) == "1.1277"

    # ---- W1105：動量殼與連續流
    at = alpha_table()
    assert [(_f4(al), _f4(sm), _f4(bk)) for _, _, al, sm, bk in at] == [
        ("0.0000", "0.0000", "-0.0000"), ("0.2929", "0.3084", "-0.5858"), ("1.0000", "1.2337", "-2.0000"),
        ("1.7071", "2.7758", "-3.4142"), ("2.0000", "4.9348", "-4.0000")]
    assert [_f4(q) for _, q, *_ in at] == ["0.0000", "0.7854", "1.5708", "2.3562", "3.1416"]
    err, real, book = fourier_check()
    assert err < 1e-12 and (_f4(real), _f4(book)) == ("3.0394", "-4.3508")
    for seed in (1, 2, 3):
        assert fourier_check(16, 0.3, 1.7, seed)[0] < 1e-12
    assert shell_fraction(2) == F(3, 4) and shell_fraction(3) == F(7, 8) and shell_fraction(1) == F(1, 2)
    for d in (1, 2, 3, 4):
        for b in (2, 3):
            assert all(_close(x, y) for x, y in zip(shell_rg(0.3, 0.2, d, b), shell_rg_steps(0.3, 0.2, d, b)))
        rl, hl = flow(0.3, 0.2, d, math.log(2))                     # 連續 ℓ = ln 2 ＝ 一次 b = 2
        assert _close(rl, shell_rg(0.3, 0.2, d)[0]) and _close(hl, shell_rg(0.3, 0.2, d)[1])
        for t in (0.5, 1.0, 2.0):
            assert all(_close(x, y, 1e-8) for x, y in zip(flow(FLOW_R0, FLOW_H0, d, t), flow_numeric(FLOW_R0, FLOW_H0, d, t)))
    assert _f4(shell_rg(1.0, 1.0, 3)[1]) == "5.6569"
    assert [(_f4(flow(FLOW_R0, FLOW_H0, 3, t)[0]), _f4(flow(FLOW_R0, FLOW_H0, 3, t)[1])) for t in FLOW_ELL] == \
        [("0.0739", "0.1218"), ("0.5460", "1.4841")]
    fl = flow_line(FLOW_R0, FLOW_H0, 3, 2.0)
    assert all(_close(h / FLOW_H0, (r / FLOW_R0) ** 1.25) for r, h in fl)                 # h ∝ r^{(d+2)/4}
    for d in (1, 3, 4):
        assert pde_residual(d) < 1e-7 and pde_residual(d, Ffn=lambda x: x) < 1e-7 and char_check(d) < 1e-12
        assert pde_residual(d, Ffn=lambda x: math.atan(x)) < 1e-7
    wrong = lambda rr, hh: rr ** 0.5 * (hh / rr ** 1.5) / (1 + (hh / rr ** 1.5) ** 2)   # Landau 的組合放進 d = 3 的 (11.87)：不滿足
    eps = 1e-5
    res = -0.5 * wrong(0.4, 0.3) + 2 * 0.4 * (wrong(0.4 + eps, 0.3) - wrong(0.4 - eps, 0.3)) / (2 * eps) + \
        2.5 * 0.3 * (wrong(0.4, 0.3 + eps) - wrong(0.4, 0.3 - eps)) / (2 * eps)
    assert abs(res) / abs(wrong(0.4, 0.3)) > 0.1

    # ---- W1106：線性化、特徵值、指數
    te = toy_eigen()
    assert te == [(F(2), (F(2), F(1))), (F(-3), (F(1), F(-2)))]
    A = np.array([[float(x) for x in row] for row in TOY_A])
    w, V = np.linalg.eigh(A)
    assert np.allclose(sorted(w), [-3, 2])
    for lam, v in te:
        assert all(sum(TOY_A[i][j] * v[j] for j in range(2)) == lam * v[i] for i in range(2))
    for c1, c2 in TOY_STARTS:
        for t in (0.3, 1.0):
            assert all(_close(x, y, 1e-10) for x, y in zip(toy_flow(c1, c2, t), toy_flow_expm(c1, c2, t)))
    assert max(abs(x) for x in toy_flow(0.0, 1.0, 5.0)) < 1e-6                          # 臨界線上：流進不動點
    for name, (Lr, Lh, d) in RG_CASES.items():
        e1, e2 = rg_exponents(Lr, Lh, d), rg_exponents_widom(Lr, Lh, d)
        assert e1 == e2 and e1["alpha"] + 2 * e1["beta"] + e1["gamma"] == 2
        assert lambdas_from(e1["nu"], e1["beta"], d) == (Lr, Lh)
    for name, key in (("ising2d", "ising2d"), ("saw2d", "saw2d"), ("gauss4", "mf")):
        e1 = rg_exponents(*RG_CASES[name])
        assert all(e1[k] == WC.EXPONENTS[key][k] for k in ("alpha", "beta", "gamma", "delta", "nu"))
    assert rg_exponents(*RG_CASES["gauss3"]) == dict(alpha=F(1, 2), beta=F(1, 4), gamma=F(1), delta=F(5), nu=F(1, 2))
    assert RG_CASES["gauss3"][:2] == gauss_lambdas(3) and RG_CASES["gauss4"][:2] == gauss_lambdas(4)
    for Lr, s in ((F(2), ("2.3026", "10.0000")), (F(1), ("4.6052", "100.0000"))):
        e, x = xi_from_flow(XI_R, Lr)
        en, xn = xi_from_flow_numeric(XI_R, Lr)
        assert (_f4(e), _f4(x)) == s and _close(e, en, 1e-8) and _close(x, XI_R ** (-1 / float(Lr)))



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
    print(f"walk_rg（第 11 章）：Koch 曲線與 Cantor 集的質量－半徑＝盒子計數＝複本公式、普通漫步與自迴避的碎形維度、Landau 與 Gaussian 的標度塌縮、"
          f"一維 Gaussian 鏈的抽格（完成平方＝數值積分、環的行列式）、一維 Ising 的抽格（暴力＝轉移矩陣＝連乘）、x 的遞迴與 ξ、"
          f"動量殼與連續流（解析＝數值）、(11.87) 的驗證、線性化與指數公式（Fraction＝Widom）全部驗證通過（{time.time() - t0:.1f} 秒）")
