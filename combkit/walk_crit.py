# 由 combkit/tools/sync_from_walk017.py 從 walk017 專案的 codes/utils/walk_crit.py 複製；請改 walk017 那份再同步。
"""
walk_crit.py —— 《Elements of the Random Walk》第 9 章 Scaling of walks and critical phenomena 的數學
（Python＋numpy／scipy（integrate、optimize、special、linalg），不 import manim；用到同資料夾的 walk_saw（二維自迴避的 c_N、⟨R²⟩_N、
全部 4^N 條漫步的重合點對數、配對）與 walk_scale（μ 的引用值）——在 walk017 是 utils 套件、在 combkit 是 combkit 套件，都用相對 import；
直接執行本檔時退回絕對 import）

約定（書上式號與筆誤見 ref/md/ch09_formulas.md）：
  J = k_B = 1。Ising 自旋 s = ±1，K = J/T，外場寫成無因次的 H = h/k_BT：權重 e^{K Σ s_i s_j + H Σ s_i}；m = M/自旋數。
  平均場 m = tanh(m/T)（T 以 T_c 為單位，T_c = 1）；二維 Ising（正方格子）T_c = 2/ln(1+√2) = 2.2692，Onsager M_s = (1 − sinh^{−4}(2/T))^{1/8}。
  一維 Ising 鏈：t = tanh K，G(r) = t^{|r|}，ξ = −1/ln t，χ = Σ_r G(r) = (1+t)/(1−t) = e^{2K}。
  臨界指數 α β γ δ ν η 用 Fraction；d 是空間維度。Landau：F = r m²/2 + u m⁴/4 − h m（書上 (9.35) 是 u = 1/3，即 m⁴/12）。
  自迴避（二維正方格子）：c_N、⟨R²⟩_N 來自 walk_saw（N ≤ 14 精確）；z 是每步的逸度，Σ c_N z^N 的奇點 z_c = 1/μ = 0.3791
  （書上 (9.37) 寫成 Σ Γ_N (z/z_0)^N，所以奇點在 z_c = z_0/μ）。
  O(n) 小圖：自旋長度 |S|² = n（每個分量 ⟨S_α²⟩ = 1）、權重 Π_bonds (1 + K S_i·S_j)；n → 0 時和 e^{K S·S} 一樣只剩自迴避路徑。
  回傳快取（lru_cache）陣列的函數不要就地修改結果（已設成唯讀）。

§9.1–9.2 ⟨R²⟩ ∝ N^{2ν} 與自發磁化 (9.1)–(9.3)，Fig. 9.1                                    [W901]
    flory_p(d) = 6/(d+2)（(9.2)，Fraction）＝ 2·walk_saw.flory_nu(d)；NU_KNOWN（ν 已知值；三維 0.588 只引用，2ν = 1.176）
    flory_gap_3d() = 0.600/0.588 − 1 = 0.0204（書上說 within a percent or so；實際 2%）
    mf_m(T)、mf_m_newton(T)        平均場 m = tanh(m/T) 的正根（brentq／牛頓，兩種算法；T ≥ 1 → 0）；mf_m_asym(T) = √(3(1−T))
    TC_2D = 2.269185；onsager_M(T)  Onsager 的二維 Ising 自發磁化；onsager_lowT(T) = 1 − 2u² − 8u³（u = e^{−4/T}，低溫團簇展開）
    beta_slope(fn, Tc, eps=1e-4)   ln M 對 ln(1 − T/T_c) 的局部斜率：平均場 0.5000、Onsager 0.1250
    MS_X = (0.5, 0.9, 0.99)        卡上的 T/T_c；ms_table() → {x: (平均場, Onsager)}：0.9575/0.9980、0.5254/0.8959、0.1725/0.6859
    BETA_CITED                     三維的 β（只引用）：Ising 0.326、n = 0 0.302（書上 Table 9.1 寫 0.308）
§9.2.1–9.2.2 O(n) 模型、對稱與對稱破缺 (9.4)–(9.15)，Fig. 9.2                                [W902]
    ring_edges(Ns)、torus_edges(L)；PAIR、RING4（4 個自旋圍成一圈＝2×2 方塊）、TORUS4（4×4 週期）＝ (自旋數, 鍵)
    on_energy(S, edges) = −Σ S_i·S_j（(9.4)）；rotation_spins(n=2, Ns=8, seed=902) → (S, S')：8 個平面單位向量一圈（大致朝右）與全部轉 120° 之後
    rotation_check() → (E, E', M⃗, M⃗')：E = −6.3851 不變、|M⃗| = 6.7943 不變、方向 −13.1° → 106.9°（angle_deg(v)）；n = 3、4 也核對
    mc_snapshot(L=8, T, sweeps=400)  畫箭頭用的一張自旋圖 (L, L)（±1，從 +1 開始）：SNAP_T = (1.5, 3.5) → 朝上 64、33 個（共 64）
    ising_exact(graph, K, H=0)     逐一列舉 2^Ns 個組態 → dict：Z、m = ⟨M⟩/Ns、absm = ⟨|M|⟩/Ns、M2 = ⟨M²⟩、e（每個自旋的能量）、
                                   PM = {M: 機率}、corr0 = [⟨s_0 s_j⟩]；ring4_closed(K) = (Z, absm, P(M=4), P(M=2), P(M=0))（閉式）
    TINY_K = (1.0, 0.2)            卡上的兩個溫度（T = 1、5）：4 自旋 Z = 121.2329、17.3497；⟨|m|⟩ = 0.9337、0.4871；
                                   P(M = ±4) = 0.4504、0.1283；P(±2) = 0.0330、0.2306；P(0) = 0.0333、0.2823；⟨M⟩ = 0；2 自旋 K = 1：Z = 6.1723、⟨|m|⟩ = 0.8808
    FIELD_H = (0.1, 0.01)          外場：K = 1 時 ⟨m⟩：4 自旋 0.3506、0.0367；4×4 週期 0.9210、0.1584（格子越大，同一個小 H 越接近飽和；
                                   無限大格子 T = 1 的 M_s = onsager_M(1.0) = 0.9993）
    mc_run(L, T, sweeps, replicas, seed=902)   二維 Ising L×L 週期、棋盤式 Metropolis（numpy 向量化、全部從 +1 開始）→ (m, e)：
                                   各為 (sweeps, replicas)；快取、唯讀
    flip_count(ms, m0=0.5)         m 在 > m0 與 < −m0 之間來回的次數（遲滯，忽略 0 附近的抖動）
    flip_tau(L, T=2.0)             平均翻轉等待時間（sweep）＝ 4000 × 50 個複本 ÷ 翻轉次數：L = 4、6、8 → 133、608、2151（各約 0.3–0.5 秒，快取）
    mc_trace(L, T, n_pts=400)      畫圖用的 m(t)（(sweep, m) 各 n_pts 點）：MC_TRACES = ((8, 3.0), (4, 2.0), (8, 2.0))；trace_replica(L, T)：
                                   T < T_c 取第一個翻轉次數 = 複本平均的複本（典型的一條）：第 25、11 個；高溫取第 0 個
                                   trace_flips(L, T)：(4, 2.0) 30 次、(8, 2.0) 2 次；(8, 3.0) 在 0 附近亂跳（平均 0.0229、⟨|m|⟩ = 0.3470）
    mc_check_L4()                  L = 4、T = 2 的 Metropolis 平均 vs 精確列舉（⟨|m|⟩、e）：驗證取樣器
§9.2.2–9.2.3 臨界指數、標度律、Widom 形式 (9.13)–(9.27)，Table 9.1，Ex. 9.1–9.2               [W903]
    ising1d_G(r, K)、ising1d_xi(K)、ising1d_chi(K)、ising1d_chi_partial(K, R)（Σ_{|r|≤R}）、ising1d_transfer_G(r, K)（轉移矩陣，第二種算法）
    ising1d_ring_closed(N, K) → (G 列表, χ_N)：環的閉式 G = (t^r + t^{N−r})/(1 + t^N)；與 ising_exact(環) 一致
    K1D = 1.0                      卡上：t = 0.7616，G(r) = 1、0.7616、0.5800、0.4417、0.3364（r = 0..4），ξ = 3.6719，χ = e² = 7.3891；
                                   部分和 R = 5、10、20：5.7520、6.9696、7.3615；16 自旋的環 χ_16 = 7.2021
    ising1d_xi2_moment(K) = Σr²G/ΣG = 2t/(1−t)²（Table 9.1 的比值是長度的平方）；÷ ξ² → 1.9877、1.9998、2.0000（K = 1、2、3）
    FISHER1D_K = (1, 2, 3)         χ/ξ = 2.0123、2.0002、2.0000 → 2（χ ∝ ξ^{2−η}、一維 η = 1）
    EXPONENTS = {"mf", "ising2d", "saw2d"}   各為 dict(alpha, beta, gamma, delta, nu, eta, d)（Fraction；平均場的 hyperscaling 用 d = 4）
    BOOK_3D、CITED_3D              書上 Table 9.1 的三維 n = 0 欄（α = 0.12、β = 0.308、γ = 1.19、δ = 4.865、ν = 0.602、η = 0.023）；文獻值（只引用）
    laws(e) → dict(rushbrooke = α+2β+γ, fisher = (2−η)ν, hyper_alpha = 2 − dν, hyper_delta = (d+2−η)/(d−2+η), widom = β(δ−1))
    from_nu_eta(nu, eta, d) → (α, β, γ, δ)（由 ν、η 推出，第二種算法）
    book3d_alpha() → (由 (9.16) 2 − 2β − γ, 由 (9.18) 2 − dν, (9.16) 左邊)：0.194、0.194、1.926（書上 α = 0.12）
    widom_exponents(alpha, Delta) → (β, γ, δ)：β = 2 − α − Δ、γ = Δ − β、δ = Δ/β（Ex. 9.1–9.2；Fraction）
    widom_G(r, h, alpha, Delta)    一個具體的齊次函數 |r|^{2−α} f(h/|r|^Δ)（r < 0；r = 0 取極限）；widom_slopes(alpha, Delta, method)：
                                   數值斜率 (β, γ, 1/δ)（'exact' 用 f 的導數、'fd' 用有限差分）；α = 0、Δ = 15/8 → 0.1250、1.7500、0.0667
§9.3 Ginzburg–Landau–Wilson (9.28)–(9.35)，Ex. 9.3–9.4                                      [W904]
    fourier_check(N=8, seed=904)   (9.30)–(9.32) 的離散版：Σ S_i² = N Σ_q|S(q)|²、Σ(∇S)² = N Σ_q 4sin²(q/2)|S(q)|²、Σ S_i⁴ = N Σ_{q1+…+q4≡0} SSSS
                                   → 三個相對誤差（< 1e-12）
    landau_F(m, r, u=1, h=0)；landau_min(r, u=1, h=0)（三次方程的實根中 F 最小者）；r = −1：m* = 1、F_min = −0.25；r = 1：m* = 0
    LANDAU_R = (1.0, 0.0, -1.0)    雙井的三條曲線；LANDAU_H = 0.1 時 r = −1 的 m* = 1.0467
    ex94_m(r, h)（u = 1/3，書上 (9.35)）；ex94_closed(r) = (√(−3r), −3r²/4)；ex94_scaling(x, r) = m/|r|^{1/2}（h = x|r|^{3/2}）：
                                   x = 1：2.1038（r = −0.01、−0.0001 一樣）；ex94_slopes() → (β, Δ, α, 1/δ, γ) 數值：0.5000、1.5000、0.0000、0.3333、1.0000
    zmf_avg(V, r=-1, u=1, h=0, what)   (9.34) 的有限體積平均場積分（n = 1）：what = 'm'、'absm'、'f'（−ln Z/V）；method 'quad'／'grid'
    ZMF_V = (10, 100, 1000)        h = 0.01（ZMF_H）：⟨m⟩ = 0.0869、0.7576、1.0042（極小值 m* = 1.0050）；h = 0：⟨m⟩ = 0、⟨|m|⟩ = 0.8888、0.9922、0.9992；
                                   −ln Z/V = −0.2719、−0.2397、−0.2478（→ F_min = −0.25）
§9.4 由磁鐵讀出 ⟨R²⟩ (9.36)–(9.41)                                                          [W905]
    ZC_2D = 1/μ = 0.379052；saw_gf(z, Nmax=14) → (χ, Σ c_N⟨R²⟩_N z^N, ξ²)；saw_gf_realspace(z, Nmax=8)：由終點的 G(z; r) 算（第二種算法）
    SAW_Z = (0.1, 0.15, 0.2)       ξ² = 0.6347、1.2799、2.4593（N ≤ 14 的部分和）；saw_table() → [(N, c_N, c_N⟨R²⟩_N)]，N = 1..6：
                                   (1, 4, 4)、(2, 12, 32)、(3, 36, 164)、(4, 100, 704)、(5, 284, 2716)、(6, 780, 9808)
    coef_log(N, s, zc)、coef_log_rec（遞迴，第二種算法）、coef_exact = e^{coef_log}、coef_asym(N, s, zc) = N^{s−1} z_c^{−N−s}/Γ(s)；
    coef_ratio(N, s) = exact/asym   S_GAMMA = 43/32、S_GAMMA2NU = 91/32，COEF_N = (10, 100, 1000)：s = 43/32：1.0227、1.0023、1.0002；
                                   s = 91/32：1.2760、1.0264、1.0026
    coef_ratio_book(N, s, delta)   照書上 N^{δ−1}：δ = 91/5（二維）時 N = 10、100 → 1.4e-17、1.9e-34（完全不對）
    r2_prefactor(gamma, nu, zc, K1=1, book=False)   (9.41) K₁Γ(γ)/(z_c^{2ν}Γ(γ+2ν))；book=True 照書上 z_c^γ
    r2_ratio(N, gamma, nu, zc, rec=False)   兩個係數的比 ÷ (前因子 · N^{2ν})（→ 1）：N = 10、100、1000：1.2476、1.0240、1.0024
    zc_pow() = (z_c^{−2ν}, z_c^{−γ}) = (μ^{3/2}, μ^{43/32}) = 4.2850、3.6823（正確 vs 書上）；前因子 2.1962 vs 書上 1.8873（K₁ = 1）
    SAW_N = (10, 12, 14)           2ν_eff = 1.4442、1.4505、1.4556 → 1.5（walk_saw.two_nu_eff）；⟨R²⟩_N/N^{3/2} = 0.8299、0.8224、0.8168
§9.5–9.6 O(n) → 自迴避、高斯積分、Wick (9.42)–(9.77)，(9:S-1)–(9:S-22)，Ex. 9.5–9.6        [W906]
    hs_lhs(s, u=1) = e^{−us²/4}；hs_rhs(s, u=1)：(1/√(πu))∫dφ e^{−φ²/u − iφs} 的數值積分（(9.46)）；s = 1、2：0.7788、0.3679
    hs_book72(u=1, method)         照書上 (9.72) 從 0 積到 ∞：0.3894 − 0.2395i（閉式 ½e^{−u/4}(1 − i erfi(√u/2))；正確 0.7788）
    visit_weight(k, u=1) = e^{−uk²/4}；intersection_factor(u=1) = e^{−u}/(e^{−u/4})² = e^{−u/2} = 0.6065；pair_sum(k) = k²
    weak_saw_count(N, u, sign=-1)  Σ_w e^{−(u/2)n_w}（n_w = 重合點對數，walk_saw）：WEAK_N = 8，WEAK_U = 0、1、4、∞ → 65536、23695.11、7988.15、5916；
                                   sign=+1 照書上 (9.75) 少了負號：u = 1 → 6.2e+05（越相交越多）
    heat_kernel(x, eps)、heat_kernel_fourier(x, eps)（數值逆 Fourier）、heat_kernel_book(x, eps, d)；kernel_mass_book(eps, d) = (ε/d)^{d/2}
                                   （ε = 0.1：d = 1 → 0.3162、d = 3 → 0.00609；正確的核積分是 1、變異數 ε）
    GAUSS_T = ((2,−1),(−1,2))、GAUSS_H = (1, 0)   (9:S-11)：gauss_I() = 2π/√3·e^{1/3} = 5.0627；gauss_I_quad()（dblquad）
    gauss_moment(T, a, b)          ⟨x₁^a x₂^b⟩（Gauss–Hermite，精確）：⟨x₁x₂⟩ = 1/3、⟨x₁²⟩ = 2/3、⟨x₁²x₂²⟩ = 2/3、⟨x₁⁴⟩ = 4/3
    wick4(C, r, s, t, u)、wick4_book(...)（照 (9:S-19) 第三行）：T = ((2,−1),(−1,3)) 的 ⟨x₁²x₂²⟩ = 0.32（書上的寫法得 0.44）
    wick_count(m) = (2m)!/(2^m m!) = 1、3、15、105（m = 1..4）；wick_count_enum（walk_saw.matchings 的個數）
    trotter_err(eps, order=1|2)    ‖e^{(X+Y)ε} − e^{Xε}e^{Yε}(e^{−ε²[X,Y]/2})‖，X = −iφ（對角，種子 906）、Y = ½∇²（6 點環）：TROTTER_EPS = 0.1、0.05、0.025：
                                   一階 0.005784、0.001502、0.0003827（比 3.85、3.92 → 4：∝ ε²）；乘上交換子項 0.0009465、0.000124、1.589e-05（比 7.63、7.81 → 8：∝ ε³）
    trotter_lead(eps) = (ε²/2)‖[X,Y]‖   一階誤差 ÷ 領頭項 = 0.9268、0.9627、0.9811 → 1
    trotter_global(M, t=1)         M 步的總誤差（∝ 1/M）：TROTTER_M = 10、20、40：0.02366、0.0117、0.005821
    ex95_check(r=1)                Ex. 9.5：V = (r/2 − ½∇²)I + iφ，∫_0^∞ e^{−tV}dt 是 V 的左、右逆（三個誤差 < 1e-10）
    LADDER（2×3 梯子，6 點 7 鍵；格點 (x, y) → x + 3y：0–3 是左邊的橫檔、0–5 是對角）  LADDER_K = 0.5、LADDER_PAIRS = ((0, 3), (0, 5))
                                   ladder_G(n, K, pair)：閉式 G_{03} = (K+K³+K⁵+nK⁵)/(1+2nK⁴+nK⁶)、G_{05} = (3K³+K⁵)/(1+2nK⁴+nK⁶)；
                                   ladder_G_ising(K, pair)（n = 1 暴力 64 組態）、ladder_G_xy(K, pair)（n = 2，角度格點，精確）、ladder_saws(a, b)（DFS 長度表）
                                   K = 0.5：G_{03} = 0.6027（n=1）、0.5610（n=2）、0.65625（n→0 ＝ 自迴避 K+K³+K⁵）；G_{05} = 0.3562、0.3171、0.40625（3K³+K⁵）
summary()                         場景卡與講義用的全部數字（dict，字串已照卡上的小數位格式化）

check() 驗證全部（每個結果至少兩種獨立算法、筆誤的更正、plan_ch09.md 場景卡與講義的每個數字）；直接執行本檔會跑 check()（約 4 秒）。
`python walk_crit.py --summary` 印出 summary() 的全部數字。
最慢的呼叫：flip_tau(L)／mc_trace(L, T) 每個 (L, T) 第一次約 0.3–0.5 秒（4000 sweep × 50 複本，快取）、mc_snapshot 約 0.03 秒、
weak_saw_count(8, u) 約 0.01 秒、ladder_G_xy 約 5 毫秒，其餘都很快。scipy 的 IntegrationWarning 不影響結果。
"""
from __future__ import annotations

import math
from fractions import Fraction
from functools import lru_cache

import numpy as np
from scipy import integrate, linalg, optimize, special

try:                                    # walk017：utils 套件；combkit：combkit 套件
    from . import walk_saw as SW
    from . import walk_scale as WS
except ImportError:                     # 直接執行本檔
    import walk_saw as SW
    import walk_scale as WS

F = Fraction


def _ro(a):
    a = np.asarray(a)
    a.flags.writeable = False
    return a


# ---------------------------------------------------------------- §9.1–9.2  [W901]
def flory_p(d):
    """(9.2) p = 6/(d+2)（Fraction）＝ 2ν_Flory。"""
    return F(6, d + 2)


NU_KNOWN = {1: F(1), 2: F(3, 4), 3: 0.588, 4: F(1, 2)}      # 三維只引用（與 walk_saw.NU_QUOTED 一致）


def flory_gap_3d(nu3=0.588):
    """三維 Flory 3/5 比已知 ν 大多少（相對）。"""
    return 0.6 / nu3 - 1


def mf_m(T):
    """平均場 m = tanh(m/T) 的正根（T_c = 1；brentq）。"""
    if T >= 1:
        return 0.0
    return optimize.brentq(lambda m: m - math.tanh(m / T), 1e-300, 1.0, xtol=1e-16, maxiter=500)


def mf_m_newton(T, n=500):
    """同上，用牛頓法從 m = 1 往下（g(m) = m − tanh(m/T) 在 m > 0 是凸的，單調收斂；第二種算法）。"""
    if T >= 1:
        return 0.0
    m = 1.0
    for _ in range(n):
        th = math.tanh(m / T)
        step = (m - th) / (1 - (1 - th * th) / T)
        m -= step
        if abs(step) < 1e-16:
            break
    return m


def mf_m_asym(T):
    """T → 1⁻ 的領頭項 √(3(1−T))（β = 1/2）。"""
    return math.sqrt(3 * (1 - T)) if T < 1 else 0.0


TC_2D = 2 / math.log(1 + math.sqrt(2))


def onsager_M(T):
    """二維 Ising（正方格子，J = 1）的自發磁化 (1 − sinh^{−4}(2/T))^{1/8}。"""
    if T >= TC_2D:
        return 0.0
    return (1 - math.sinh(2 / T) ** -4) ** 0.125


def onsager_lowT(T):
    """低溫團簇展開：翻一個自旋（能量 8，磁化少 2）＋相鄰兩個（能量 12、兩種方向、少 4）：1 − 2u² − 8u³，u = e^{−4/T}。"""
    u = math.exp(-4 / T)
    return 1 - 2 * u ** 2 - 8 * u ** 3


def beta_slope(fn, Tc, eps=1e-4):
    """ln M 對 ln(1 − T/T_c) 的局部斜率（取 1 − T/T_c = eps 與 eps/10）。"""
    a, b = fn(Tc * (1 - eps)), fn(Tc * (1 - eps / 10))
    return (math.log(a) - math.log(b)) / math.log(10)


MS_X = (0.5, 0.9, 0.99)
BETA_CITED = {"3d_ising": 0.326, "3d_saw": 0.302, "book_3d_saw": 0.308}


def ms_table():
    return {x: (mf_m(x), onsager_M(x * TC_2D)) for x in MS_X}


# ---------------------------------------------------------------- §9.2.1–9.2.2  [W902]
def ring_edges(Ns):
    return tuple((i, (i + 1) % Ns) for i in range(Ns))


def torus_edges(L):
    """L×L 週期正方格子（L ≥ 3）的 2L² 條鍵；格點 (x, y) 的編號 x + L·y。"""
    e = []
    for y in range(L):
        for x in range(L):
            i = x + L * y
            e.append((i, (x + 1) % L + L * y))
            e.append((i, x + L * ((y + 1) % L)))
    return tuple(e)


PAIR = (2, ((0, 1),))
RING4 = (4, ring_edges(4))
TORUS4 = (16, torus_edges(4))


def on_energy(S, edges):
    """(9.4) H = −Σ_{鍵} S_i·S_j（J = 1）。"""
    S = np.asarray(S, float)
    return -float(sum(np.dot(S[i], S[j]) for i, j in edges))


ROT_N, ROT_NS, ROT_SPREAD, ROT_ANGLE = 2, 8, 0.6, 2 * math.pi / 3


def rotation_spins(n=ROT_N, Ns=ROT_NS, seed=902, spread=ROT_SPREAD):
    """Ns 個單位 n 維向量排成一圈（大致朝 x 方向、亂度 spread）與「全部一起轉」之後的樣子 → (S, S')，各 (Ns, n)。
    n = 2 轉 ROT_ANGLE = 120°（好畫）；n ≥ 3 用隨機的旋轉矩陣（QR）。"""
    rng = np.random.default_rng(seed)
    S = np.zeros((Ns, n))
    S[:, 0] = 1.0
    S += spread * rng.normal(size=(Ns, n))
    S /= np.linalg.norm(S, axis=1, keepdims=True)
    if n == 2:
        c, s_ = math.cos(ROT_ANGLE), math.sin(ROT_ANGLE)
        Q = np.array([[c, -s_], [s_, c]])
    else:
        Q, R = np.linalg.qr(rng.normal(size=(n, n)))
        Q = Q * np.sign(np.diag(R))
        if np.linalg.det(Q) < 0:
            Q[:, 0] *= -1
    return _ro(S), _ro(S @ Q.T)


def rotation_check(n=ROT_N, Ns=ROT_NS, seed=902, spread=ROT_SPREAD):
    """(E, E', M⃗, M⃗')：能量 (9.4) 不變、|M⃗| 不變、方向轉了。"""
    S, S2 = rotation_spins(n, Ns, seed, spread)
    ed = ring_edges(Ns)
    return on_energy(S, ed), on_energy(S2, ed), S.sum(0), S2.sum(0)


def angle_deg(v):
    """二維向量的方向角（度，−180..180）。"""
    return math.degrees(math.atan2(v[1], v[0]))


@lru_cache(maxsize=None)
def _spins(Ns):
    a = np.arange(2 ** Ns)[:, None] >> np.arange(Ns)
    return _ro((1 - 2 * (a & 1)).astype(np.int16))


def ising_exact(graph, K, H=0.0):
    """逐一列舉 2^Ns 個 ±1 組態（Ns ≤ 16），權重 e^{K Σ_鍵 s s + H Σ s}。"""
    Ns, edges = graph
    s = _spins(Ns).astype(np.int64)
    bond = np.zeros(len(s), np.int64)
    for i, j in edges:
        bond += s[:, i] * s[:, j]
    M = s.sum(1)
    ex = K * bond + H * M
    c = ex.max()
    w = np.exp(ex - c)
    Zr = w.sum()
    PM = {int(v): float(w[M == v].sum() / Zr) for v in np.unique(M)}
    return dict(Z=float(Zr * math.exp(c)), m=float((w * M).sum() / Zr / Ns),
                absm=float((w * np.abs(M)).sum() / Zr / Ns), M2=float((w * M * M).sum() / Zr),
                e=float(-(w * bond).sum() / Zr / Ns), PM=PM,
                corr0=[float((w * s[:, 0] * s[:, j]).sum() / Zr) for j in range(Ns)])


def ring4_closed(K):
    """4 個自旋一圈（H = 0）的閉式：(Z, ⟨|m|⟩, P(M=4), P(M=2), P(M=0))。"""
    a, b = math.exp(4 * K), math.exp(-4 * K)
    Z = 2 * a + 12 + 2 * b
    return Z, (2 * a + 4) / Z, a / Z, 4 / Z, (4 + 2 * b) / Z


def pair_closed(K):
    """2 個自旋（一條鍵）：(Z, ⟨|m|⟩) = (2e^K + 2e^{−K}, e^K/(e^K + e^{−K}))。"""
    return 2 * math.exp(K) + 2 * math.exp(-K), math.exp(K) / (math.exp(K) + math.exp(-K))


TINY_K = (1.0, 0.2)
FIELD_H = (0.1, 0.01)

MC_SEED = 902
MC_SWEEPS = 4000
MC_REPLICAS = 50
MC_M0 = 0.5
MC_T = 2.0
MC_L = (4, 6, 8)
MC_TRACES = ((8, 3.0), (4, 2.0), (8, 2.0))


def _metropolis(L, T, sweeps, replicas, seed):
    rng = np.random.default_rng(seed)
    s = np.ones((replicas, L, L), dtype=np.int32)
    ii, jj = np.indices((L, L))
    masks = [((ii + jj) % 2 == c)[None] for c in (0, 1)]
    acc = np.minimum(1.0, np.exp(-np.arange(-8, 9, 4) / T))         # ΔE = −8, −4, 0, 4, 8
    ms = np.empty((sweeps, replicas))
    es = np.empty((sweeps, replicas))
    for t in range(sweeps):
        for mk in masks:
            nb = np.roll(s, 1, 1) + np.roll(s, -1, 1) + np.roll(s, 1, 2) + np.roll(s, -1, 2)
            flip = (rng.random(s.shape) < acc[(2 * s * nb + 8) // 4]) & mk
            s[flip] *= -1
        ms[t] = s.mean(axis=(1, 2))
        es[t] = -(s * (np.roll(s, 1, 1) + np.roll(s, 1, 2))).sum(axis=(1, 2)) / (L * L)
    return ms, es, s


@lru_cache(maxsize=None)
def mc_run(L, T, sweeps=MC_SWEEPS, replicas=MC_REPLICAS, seed=MC_SEED):
    """二維 Ising L×L 週期、棋盤式 Metropolis（全部從 +1 開始）→ (m, e)：各 (sweeps, replicas)，每個 sweep 之後記一次。"""
    ms, es, _ = _metropolis(L, T, sweeps, replicas, seed)
    return _ro(ms), _ro(es)


@lru_cache(maxsize=None)
def mc_snapshot(L=8, T=MC_T, sweeps=400, seed=MC_SEED):
    """單一系統（從 +1 開始）跑 sweeps 次之後的自旋 (L, L)（±1；畫箭頭用）。"""
    return _ro(_metropolis(L, T, sweeps, 1, seed)[2][0])


SNAP_T = (1.5, 3.5)


def flip_count(ms, m0=MC_M0):
    """m 在「> m0」與「< −m0」兩個狀態之間換邊的次數（從 + 開始；0 附近的抖動不算）。"""
    ms = np.asarray(ms)
    cols = ms.T if ms.ndim == 2 else ms[None]
    n = 0
    for col in cols:
        sg = np.concatenate(([1.0], np.sign(col[np.abs(col) > m0])))
        n += int(np.count_nonzero(sg[1:] != sg[:-1]))
    return n


def flip_tau(L, T=MC_T, sweeps=MC_SWEEPS, replicas=MC_REPLICAS):
    """平均翻轉等待時間（sweep）＝ 總 sweep 數 × 複本數 ÷ 翻轉次數。"""
    ms, _ = mc_run(L, T, sweeps, replicas)
    return sweeps * replicas / flip_count(ms)


def trace_replica(L, T):
    """畫哪一個複本：T < T_c 時取第一個「翻轉次數 = 50 個複本的平均（四捨五入）」的複本（典型的一條），否則第 0 個。"""
    ms, _ = mc_run(L, T)
    if T >= TC_2D:
        return 0
    n = [flip_count(ms[:, r]) for r in range(ms.shape[1])]
    return n.index(round(sum(n) / len(n)))


def mc_trace(L, T, n_pts=400):
    """畫圖用：trace_replica 那一條的 m(t)，均勻取 n_pts 點 → (sweep 陣列, m 陣列)。"""
    ms, _ = mc_run(L, T)
    idx = np.linspace(0, len(ms) - 1, n_pts).round().astype(int)
    return idx + 1, ms[idx, trace_replica(L, T)]


def trace_flips(L, T):
    return flip_count(mc_run(L, T)[0][:, trace_replica(L, T)])


def mc_check_L4(T=MC_T):
    """L = 4 的 Metropolis 平均 vs 精確列舉 → ((⟨|m|⟩_MC, 精確), (e_MC, 精確))。"""
    ms, es = mc_run(4, T)
    ex = ising_exact(TORUS4, 1 / T)
    return (float(np.abs(ms[200:]).mean()), ex["absm"]), (float(es[200:].mean()), ex["e"])


# ---------------------------------------------------------------- §9.2.2–9.2.3  [W903]
def ising1d_G(r, K):
    return math.tanh(K) ** abs(r)


def ising1d_xi(K):
    return -1 / math.log(math.tanh(K))


def ising1d_chi(K):
    t = math.tanh(K)
    return (1 + t) / (1 - t)


def ising1d_chi_partial(K, R):
    return sum(ising1d_G(r, K) for r in range(-R, R + 1))


def ising1d_transfer_G(r, K):
    """轉移矩陣 [[e^K, e^{−K}], [e^{−K}, e^K]] 的兩個本徵值之比的 r 次方（無窮長鏈）。"""
    lam = np.linalg.eigvalsh(np.array([[math.exp(K), math.exp(-K)], [math.exp(-K), math.exp(K)]]))
    return float((lam[0] / lam[1]) ** abs(r))


def ising1d_ring_closed(N, K):
    """N 個自旋的環：G(r) = (t^r + t^{N−r})/(1 + t^N)（r = 0..N−1）、χ_N = Σ_r G(r)。"""
    t = math.tanh(K)
    G = [(t ** r + t ** (N - r)) / (1 + t ** N) for r in range(N)]
    return G, sum(G)


def ising1d_xi2_moment(K):
    """Σ r² G / Σ G（r 從 −∞ 到 ∞）＝ 2t/(1−t)²（Table 9.1 的比值是長度的平方）。"""
    t = math.tanh(K)
    return 2 * t / (1 - t) ** 2


def ising1d_xi2_moment_sum(K, R=4000):
    num = sum(r * r * ising1d_G(r, K) for r in range(-R, R + 1))
    return num / ising1d_chi_partial(K, R)


K1D = 1.0
FISHER1D_K = (1.0, 2.0, 3.0)

EXPONENTS = {
    "mf": dict(alpha=F(0), beta=F(1, 2), gamma=F(1), delta=F(3), nu=F(1, 2), eta=F(0), d=4),
    "ising2d": dict(alpha=F(0), beta=F(1, 8), gamma=F(7, 4), delta=F(15), nu=F(1), eta=F(1, 4), d=2),
    "saw2d": dict(alpha=F(1, 2), beta=F(5, 64), gamma=F(43, 32), delta=F(91, 5), nu=F(3, 4), eta=F(5, 24), d=2),
}
BOOK_3D = dict(alpha=0.12, beta=0.308, gamma=1.19, delta=4.865, nu=0.602, eta=0.023, d=3)
CITED_3D = dict(nu=0.5876, gamma=1.1570, eta=0.031, alpha=0.237, beta=0.302)   # 三維 n = 0 的文獻值（只引用）


def laws(e):
    a, b, g, de, nu, eta, d = (e[k] for k in ("alpha", "beta", "gamma", "delta", "nu", "eta", "d"))
    return dict(rushbrooke=a + 2 * b + g, fisher=(2 - eta) * nu, hyper_alpha=2 - d * nu,
                hyper_delta=(d + 2 - eta) / (d - 2 + eta), widom=b * (de - 1))


def from_nu_eta(nu, eta, d):
    """由 ν、η、d 推出 (α, β, γ, δ)：hyperscaling、Fisher、β = ν(d−2+η)/2。"""
    return 2 - d * nu, nu * (d - 2 + eta) / 2, (2 - eta) * nu, (d + 2 - eta) / (d - 2 + eta)


def book3d_alpha():
    e = BOOK_3D
    return 2 - 2 * e["beta"] - e["gamma"], 2 - 3 * e["nu"], e["alpha"] + 2 * e["beta"] + e["gamma"]


def widom_exponents(alpha, Delta):
    """Ex. 9.1–9.2：β = 2 − α − Δ、γ = Δ − β、δ = Δ/β。"""
    alpha, Delta = F(alpha), F(Delta)
    beta = 2 - alpha - Delta
    return beta, Delta - beta, Delta / beta


def _wf(x, q):
    """f(x) = −(1+x²)^{q/2} − x(1+x²)^{(q−1)/2} 與它的一、二階導數。"""
    a = 1 + x * x
    f = -a ** (q / 2) - x * a ** ((q - 1) / 2)
    f1 = -q * x * a ** (q / 2 - 1) - a ** ((q - 1) / 2) - (q - 1) * x * x * a ** ((q - 3) / 2)
    f2 = (-q * a ** (q / 2 - 1) - q * (q - 2) * x * x * a ** (q / 2 - 2)
          - 3 * (q - 1) * x * a ** ((q - 3) / 2) - (q - 1) * (q - 3) * x ** 3 * a ** ((q - 5) / 2))
    return f, f1, f2


def widom_G(r, h, alpha, Delta):
    """一個具體的 Widom 形式 (9.23)：G = |r|^{2−α} f(h/|r|^Δ)（r < 0；f 見 _wf，q = (2−α)/Δ）；r = 0 取極限 −(|h|^q + h|h|^{q−1})。"""
    alpha, Delta = float(alpha), float(Delta)
    q = (2 - alpha) / Delta
    if r == 0:
        return -(abs(h) ** q + h * abs(h) ** (q - 1))
    ar = abs(r)
    return ar ** (2 - alpha) * _wf(h / ar ** Delta, q)[0]


def widom_slopes(alpha, Delta, method="exact"):
    """數值斜率 (β, γ, 1/δ)：M(r, 0) ∝ |r|^β、χ(r, 0) ∝ |r|^{−γ}、M(0, h) ∝ h^{1/δ}（r = −1e-2 與 −1e-4、h = 1e-3 與 1e-5）。
    method='exact' 用 f 的解析導數；'fd' 用有限差分（第二種算法）。"""
    a, D = float(alpha), float(Delta)
    q = (2 - a) / D

    def M(r, h):
        if method == "exact":
            if r == 0:
                return 2 * q * abs(h) ** (q - 1)
            return -abs(r) ** (2 - a - D) * _wf(h / abs(r) ** D, q)[1]
        dh = 1e-6 * (abs(r) ** D if r else abs(h))
        return -(widom_G(r, h + dh, a, D) - widom_G(r, h - dh, a, D)) / (2 * dh)

    def chi(r):
        if method == "exact":
            return -abs(r) ** (2 - a - 2 * D) * _wf(0.0, q)[2]
        dh = 1e-3 * abs(r) ** D
        return -(widom_G(r, dh, a, D) - 2 * widom_G(r, 0.0, a, D) + widom_G(r, -dh, a, D)) / dh ** 2

    ln10 = 2 * math.log(10)
    b = (math.log(M(-1e-2, 0.0)) - math.log(M(-1e-4, 0.0))) / ln10
    g = -(math.log(chi(-1e-2)) - math.log(chi(-1e-4))) / ln10
    di = (math.log(M(0, 1e-3)) - math.log(M(0, 1e-5))) / ln10
    return b, g, di


WIDOM_CASE = (F(0), F(15, 8))      # 二維 Ising：α = 0、Δ = β + γ = 15/8


# ---------------------------------------------------------------- §9.3  [W904]
def fourier_check(N=8, seed=904):
    """(9.30)–(9.32) 的離散版（一維 N 點環，S(q) = (1/N)Σ S_i e^{iqr_i}）→ 三個相對誤差。"""
    rng = np.random.default_rng(seed)
    S = rng.normal(size=N)
    Sq = np.fft.ifft(S)                                   # (1/N) Σ_i S_i e^{+2πi k i/N}
    q = 2 * np.pi * np.arange(N) / N
    e2 = abs(np.sum(S ** 2) - N * np.sum(np.abs(Sq) ** 2).real) / np.sum(S ** 2)
    g = np.sum((np.roll(S, -1) - S) ** 2)
    eg = abs(g - N * np.sum(4 * np.sin(q / 2) ** 2 * np.abs(Sq) ** 2)) / g
    P = np.array([sum(Sq[k1] * Sq[(k - k1) % N] for k1 in range(N)) for k in range(N)])
    r4 = N * sum(P[k] * P[(-k) % N] for k in range(N))
    e4 = abs(np.sum(S ** 4) - r4) / np.sum(S ** 4)
    return float(e2), float(eg), float(e4)


def landau_F(m, r, u=1.0, h=0.0):
    return r * m * m / 2 + u * m ** 4 / 4 - h * m


def landau_min(r, u=1.0, h=0.0):
    """∂F/∂m = u m³ + r m − h = 0 的實根中 F 最小者（h = 0、r < 0 時取正根）。"""
    roots = np.roots([u, 0.0, r, -h])
    re = sorted(float(z.real) for z in roots if abs(z.imag) < 1e-9)
    Fs = [landau_F(m, r, u, h) for m in re]
    best = min(Fs)
    cands = [m for m, f in zip(re, Fs) if f - best < 1e-14]
    return max(cands)


def landau_closed(r, u=1.0):
    """h = 0 的閉式：(m*, F_min) = (√(−r/u), −r²/4u)（r < 0），r ≥ 0 為 (0, 0)。"""
    return (math.sqrt(-r / u), -r * r / (4 * u)) if r < 0 else (0.0, 0.0)


LANDAU_R = (1.0, 0.0, -1.0)
LANDAU_H = 0.1
EX94_U = 1 / 3                      # 書上 (9.35) 的 m⁴/12


def ex94_m(r, h):
    return landau_min(r, EX94_U, h)


def ex94_closed(r):
    """(9.35) h = 0、r < 0：(m*, F_min) = (√(−3r), −3r²/4)。"""
    return math.sqrt(-3 * r), -0.75 * r * r


def ex94_scaling(x, r):
    """m/|r|^{1/2}，h = x|r|^{3/2}（r < 0）→ ℳ(x)（與 r 無關）。"""
    return ex94_m(r, x * abs(r) ** 1.5) / abs(r) ** 0.5


def ex94_slopes():
    """數值：(β, Δ, α, 1/δ, γ)。β 由 m*(r)、Δ 由「固定 m/|r|^β 所需的 h」、α 由 F_min(r)、1/δ 由 m(0, h)、γ 由 χ = 1/F''(m*)。"""
    r1, r2 = -1e-2, -1e-4
    l = math.log(r1 / r2)
    beta = math.log(ex94_m(r1, 0) / ex94_m(r2, 0)) / l
    # Δ：令 m/|r|^{1/2} = 2 的 h（解 h = F'(m) 的反函數：h = r m + m³/3）
    hh = [r * (2 * abs(r) ** 0.5) + (2 * abs(r) ** 0.5) ** 3 / 3 for r in (r1, r2)]
    Delta = math.log(hh[0] / hh[1]) / l
    Fm = [landau_F(ex94_m(r, 0), r, EX94_U) for r in (r1, r2)]
    alpha = 2 - math.log(Fm[0] / Fm[1]) / l
    di = math.log(ex94_m(0.0, 1e-3) / ex94_m(0.0, 1e-6)) / math.log(1e3)
    chi = [1 / (r + ex94_m(r, 0) ** 2) for r in (1e-2, 1e-4)]       # r > 0：χ = 1/r
    gamma = -math.log(chi[0] / chi[1]) / math.log(100)
    return beta, Delta, alpha, di, gamma


ZMF_V = (10, 100, 1000)
ZMF_H = 0.01


def zmf_avg(V, r=-1.0, u=1.0, h=0.0, what="m", method="quad"):
    """(9.34) n = 1 的有限體積積分 ∫ dm e^{−V F(m)}：what = 'm'（⟨m⟩）、'absm'（⟨|m|⟩）、'f'（−ln Z/V）。"""
    Fmin = landau_F(landau_min(r, u, h), r, u, h)

    def w(m):
        return math.exp(-V * (landau_F(m, r, u, h) - Fmin))

    if method == "quad":
        pts = [-1.0, 0.0, 1.0]
        Z = integrate.quad(w, -3, 3, points=pts, limit=400, epsabs=0, epsrel=1e-12)[0]
        if what == "f":
            return Fmin - math.log(Z) / V
        g = (lambda m: m * w(m)) if what == "m" else (lambda m: abs(m) * w(m))
        return integrate.quad(g, -3, 3, points=pts, limit=400, epsabs=1e-300, epsrel=1e-12)[0] / Z
    x = np.linspace(-3, 3, 600001)
    wx = np.exp(-V * (landau_F(x, r, u, h) - Fmin))
    Z = integrate.trapezoid(wx, x)
    if what == "f":
        return Fmin - math.log(Z) / V
    return integrate.trapezoid((x if what == "m" else np.abs(x)) * wx, x) / Z


# ---------------------------------------------------------------- §9.4  [W905]
ZC_2D = 1 / WS.MU_QUOTED
SAW_Z = (0.1, 0.15, 0.2)
S_GAMMA = F(43, 32)
S_GAMMA2NU = F(43, 32) + F(3, 2)
COEF_N = (10, 100, 1000)
SAW_N = (10, 12, 14)


def saw_gf(z, Nmax=14):
    """(9.37)：χ = Σ c_N z^N、分子 Σ c_N⟨R²⟩_N z^N、ξ² = 分子/χ（walk_saw 的精確列舉）。"""
    c, s = SW.saw_enum(Nmax)
    chi = sum(c[N] * z ** N for N in range(Nmax + 1))
    num = sum(s[N] * z ** N for N in range(Nmax + 1))
    return chi, num, num / chi


def saw_gf_realspace(z, Nmax=8):
    """由終點：G(z; r⃗) = Σ_N (N 步到 r⃗ 的自迴避條數) z^N，再算 Σ r²G / Σ G（全部 4^N 條暴力篩選；第二種算法）。"""
    G = {(0, 0): 1.0}
    for N in range(1, Nmax + 1):
        xs, ys = SW.all_walks_xy(N)
        ok = SW.n_coinc(N) == 0
        ex, ey = xs[ok, -1].astype(int), ys[ok, -1].astype(int)
        keys, cnt = np.unique(np.stack([ex, ey], 1), axis=0, return_counts=True)
        for (a, b), n in zip(keys, cnt):
            G[(int(a), int(b))] = G.get((int(a), int(b)), 0.0) + int(n) * z ** N
    chi = sum(G.values())
    return chi, sum((a * a + b * b) * g for (a, b), g in G.items()) / chi


def saw_table(Nmax=6):
    c, s = SW.saw_enum(14)
    return [(N, c[N], s[N]) for N in range(1, Nmax + 1)]


def coef_log(N, s, zc):
    """(z_c − z)^{−s} 的 z^N 係數的對數：ln[z_c^{−N−s} Γ(N+s)/(Γ(s) N!)]（lgamma）。"""
    s = float(s)
    return float(special.gammaln(N + s) - special.gammaln(s) - special.gammaln(N + 1)) - (N + s) * math.log(zc)


def coef_log_rec(N, s, zc):
    """同上，由遞迴 a_0 = z_c^{−s}、a_k = a_{k−1}(k−1+s)/(k z_c) 的對數和（第二種算法）。"""
    s = float(s)
    return -s * math.log(zc) + math.fsum(math.log((k - 1 + s) / (k * zc)) for k in range(1, N + 1))


def coef_exact(N, s, zc):
    return math.exp(coef_log(N, s, zc))


def coef_rec(N, s, zc):
    return math.exp(coef_log_rec(N, s, zc))


def coef_asym(N, s, zc):
    s = float(s)
    return N ** (s - 1) * zc ** (-N - s) / math.gamma(s)


def coef_ratio(N, s):
    """exact ÷ N^{s−1} z_c^{−N−s}/Γ(s)（z_c 消掉）。"""
    s = float(s)
    return math.exp(coef_log(N, s, 1.0) - (s - 1) * math.log(N) + math.lgamma(s))


def coef_ratio_book(N, s, delta):
    """照書上 N^{δ−1}（而不是 N^{s−1}）的比值。"""
    s = float(s)
    return math.exp(coef_log(N, s, 1.0) - (float(delta) - 1) * math.log(N) + math.lgamma(s))


def r2_prefactor(gamma, nu, zc, K1=1.0, book=False):
    """(9.41) 的前因子 K₁Γ(γ)/(z_c^{2ν} Γ(γ+2ν))；book=True 照書上 z_c^γ。"""
    g, n2 = float(gamma), 2 * float(nu)
    return K1 * math.gamma(g) / (zc ** (g if book else n2) * math.gamma(g + n2))


def r2_ratio(N, gamma, nu, zc, rec=False):
    """[z^N]K₁K₂(z_c−z)^{−γ−2ν} ÷ [z^N]K₂(z_c−z)^{−γ} ÷ (前因子·N^{2ν})（K₁ = K₂ = 1；→ 1）；rec=True 用遞迴的係數。"""
    g, n2 = float(gamma), 2 * float(nu)
    cl = coef_log_rec if rec else coef_log
    return math.exp(cl(N, g + n2, zc) - cl(N, g, zc) - math.log(r2_prefactor(gamma, nu, zc)) - n2 * math.log(N))


def zc_pow():
    """(z_c^{−2ν}, z_c^{−γ}) = (μ^{3/2}, μ^{43/32})（二維）。"""
    return ZC_2D ** -1.5, ZC_2D ** -float(S_GAMMA)


def saw_amp(N):
    c, s = SW.saw_enum(14)
    return s[N] / c[N] / N ** 1.5


# ---------------------------------------------------------------- §9.5–9.6  [W906]
def hs_lhs(s, u=1.0):
    return math.exp(-u * s * s / 4)


def hs_rhs(s, u=1.0):
    """(9.46) 右邊：(1/√(πu))∫_{−∞}^{∞} dφ e^{−φ²/u} cos(φs)（虛部是奇函數，積分為 0）。"""
    v = integrate.quad(lambda p: math.exp(-p * p / u) * math.cos(p * s), -np.inf, np.inf, epsabs=1e-14)[0]
    return v / math.sqrt(math.pi * u)


def hs_book72(u=1.0, method="quad"):
    """照書上 (9.72) 從 0 積到 ∞：(1/√(πu))∫_0^∞ e^{−φ²/u − iφ} dφ（複數）。閉式 ½e^{−u/4}(1 − i·erfi(√u/2))。"""
    if method == "closed":
        return 0.5 * math.exp(-u / 4) * complex(1, -float(special.erfi(math.sqrt(u) / 2)))
    re = integrate.quad(lambda p: math.exp(-p * p / u) * math.cos(p), 0, np.inf, epsabs=1e-14)[0]
    im = -integrate.quad(lambda p: math.exp(-p * p / u) * math.sin(p), 0, np.inf, epsabs=1e-14)[0]
    return complex(re, im) / math.sqrt(math.pi * u)


def visit_weight(k, u=1.0):
    """一個格點被走過 k 次：(1/√(πu))∫dφ e^{−φ²/u − ikφ} = e^{−uk²/4}（(9.72) k = 1、(9.73) k = 2）。"""
    return math.exp(-u * k * k / 4)


def intersection_factor(u=1.0):
    return visit_weight(2, u) / visit_weight(1, u) ** 2


def pair_sum(k):
    """同一格點走過 k 次時 Σ_{l,l'} δ（有序對、含 l = l'）＝ k²＝ k ＋ 2·C(k,2)。"""
    return k * k


def weak_saw_count(N, u, sign=-1):
    """Σ_w e^{sign·(u/2)·n_w}（N 步、全部 4^N 條；u = inf 時只剩 n_w = 0 的自迴避）。"""
    nc = SW.n_coinc(N)
    if math.isinf(u):
        return float(np.count_nonzero(nc == 0)) if sign < 0 else math.inf
    return float(np.exp(sign * 0.5 * u * nc).sum())


WEAK_N = 8
WEAK_U = (0.0, 1.0, 4.0, math.inf)


def heat_kernel(x, eps):
    """e^{½∇²ε} 的核（一維）：(2πε)^{−1/2} e^{−x²/2ε}。"""
    return math.exp(-x * x / (2 * eps)) / math.sqrt(2 * math.pi * eps)


def heat_kernel_fourier(x, eps):
    """(1/2π)∫dk e^{ikx} e^{−εk²/2}（數值；第二種算法）。"""
    return integrate.quad(lambda k: math.cos(k * x) * math.exp(-eps * k * k / 2), -np.inf, np.inf)[0] / (2 * math.pi)


def heat_kernel_book(x, eps, d=1):
    """照書上 (9.70)：(2π)^{−d/2} exp(−(d/2)x²/ε)（一個方向的因子）。"""
    return (2 * math.pi) ** (-d / 2) * math.exp(-d * x * x / (2 * eps))


def kernel_mass_book(eps, d=1):
    """書上的核對 d 維空間積分 = (ε/d)^{d/2}（數值，d = 1 的一維積分取 d 次方）。"""
    one = integrate.quad(lambda x: math.exp(-d * x * x / (2 * eps)), -np.inf, np.inf)[0] / math.sqrt(2 * math.pi)
    return one ** d


GAUSS_T = ((2.0, -1.0), (-1.0, 2.0))
GAUSS_H = (1.0, 0.0)


def gauss_I(T=GAUSS_T, h=GAUSS_H):
    """(9:S-11)：(2π)^{N/2}/√det T · exp(½ hᵀT^{−1}h)。"""
    T, h = np.array(T, float), np.array(h, float)
    return float((2 * math.pi) ** (len(h) / 2) / math.sqrt(np.linalg.det(T)) * math.exp(0.5 * h @ np.linalg.solve(T, h)))


def gauss_I_quad(T=GAUSS_T, h=GAUSS_H):
    """同上，2×2 直接數值積分（dblquad）。"""
    (a, b), (c, d) = T
    h1, h2 = h
    f = lambda y, x: math.exp(-0.5 * (a * x * x + (b + c) * x * y + d * y * y) + h1 * x + h2 * y)
    return integrate.dblquad(f, -12, 12, -12, 12, epsabs=1e-12, epsrel=1e-12)[0]


def gauss_moment(T, a, b, n=20):
    """ρ ∝ exp(−½xᵀTx) 下的 ⟨x₁^a x₂^b⟩（2×2；Gauss–Hermite 張量求積，對多項式精確）。"""
    C = np.linalg.inv(np.array(T, float))
    Lc = np.linalg.cholesky(C)
    z, w = np.polynomial.hermite_e.hermegauss(n)
    w = w / math.sqrt(2 * math.pi)
    Z1, Z2 = np.meshgrid(z, z, indexing="ij")
    W = np.outer(w, w)
    x1 = Lc[0, 0] * Z1
    x2 = Lc[1, 0] * Z1 + Lc[1, 1] * Z2
    return float((W * x1 ** a * x2 ** b).sum())


def wick4(C, r, s, t, u):
    """(9:S-20)：⟨x_r x_s x_t x_u⟩ = C_rs C_tu + C_rt C_su + C_ru C_ts。"""
    return C[r][s] * C[t][u] + C[r][t] * C[s][u] + C[r][u] * C[t][s]


def wick4_book(C, r, s, t, u):
    """照書上 (9:S-19) 第三行括號裡的 ⟨x_r x_s⟩（應為 ⟨x_t x_u⟩）再微分得到的結果。"""
    return C[r][s] * C[r][s] + C[t][r] * C[u][s] + C[t][s] * C[u][r]


def wick_count(m):
    return math.factorial(2 * m) // (2 ** m * math.factorial(m))


def wick_count_enum(m):
    return len(SW.matchings(m))


TROTTER_EPS = (0.1, 0.05, 0.025)
TROTTER_M = (10, 20, 40)


@lru_cache(maxsize=None)
def trotter_mats(N=6, seed=906):
    """X = −iφ（對角，φ ~ N(0,1) 固定種子）、Y = ½∇²（N 點環的離散 Laplacian）。"""
    rng = np.random.default_rng(seed)
    phi = rng.normal(size=N)
    lap = -2 * np.eye(N) + np.roll(np.eye(N), 1, 1) + np.roll(np.eye(N), -1, 1)
    return _ro(-1j * np.diag(phi)), _ro(0.5 * lap)


def trotter_err(eps, order=1):
    """‖e^{(X+Y)ε} − e^{Xε}e^{Yε}‖（order=1）；order=2 再乘 e^{−(ε²/2)[X,Y]}（Ex. 9.6 的下一項）。"""
    X, Y = trotter_mats()
    ex = linalg.expm((X + Y) * eps)
    ap = linalg.expm(X * eps) @ linalg.expm(Y * eps)
    if order == 2:
        ap = ap @ linalg.expm(-(eps ** 2 / 2) * (X @ Y - Y @ X))
    return float(np.linalg.norm(ex - ap, 2))


def trotter_lead(eps):
    """一階誤差的領頭項 (ε²/2)‖[X,Y]‖。"""
    X, Y = trotter_mats()
    return float(eps ** 2 / 2 * np.linalg.norm(X @ Y - Y @ X, 2))


def trotter_global(M, t=1.0):
    X, Y = trotter_mats()
    eps = t / M
    step = linalg.expm(X * eps) @ linalg.expm(Y * eps)
    return float(np.linalg.norm(linalg.expm((X + Y) * t) - np.linalg.matrix_power(step, M), 2))


def ex95_check(r=1.0):
    """Ex. 9.5：V = (r/2 − ½∇²)I + iφ（6 點環、φ 同 trotter_mats）；V⁻¹ = ∫_0^∞ e^{−tV}dt（quad_vec）
    → (‖V·V⁻¹ − I‖, ‖V⁻¹·V − I‖, ‖V⁻¹ − inv(V)‖)。"""
    X, Y = trotter_mats()
    Id = np.eye(len(X))
    V = (r / 2) * Id - Y - X
    Vi = integrate.quad_vec(lambda t: linalg.expm(-t * V), 0, 80, epsabs=1e-13, epsrel=1e-12)[0]
    return (float(np.linalg.norm(V @ Vi - Id, 2)), float(np.linalg.norm(Vi @ V - Id, 2)),
            float(np.linalg.norm(Vi - np.linalg.inv(V), 2)))


LADDER = (6, ((0, 1), (1, 2), (3, 4), (4, 5), (0, 3), (1, 4), (2, 5)))     # 格點 (x, y) → x + 3y
LADDER_K = 0.5
LADDER_PAIRS = ((0, 3), (0, 5))


def ladder_G(n, K, pair):
    """閉式（|S|² = n、權重 Π(1 + K S·S)）：只有「端點度數 1、其他度數 0 或 2」的鍵集合存活；每個閉圈乘 n。"""
    Z = 1 + 2 * n * K ** 4 + n * K ** 6
    num = {(0, 3): K + K ** 3 + K ** 5 + n * K ** 5, (0, 5): 3 * K ** 3 + K ** 5}[tuple(pair)]
    return num / Z


def ladder_G_ising(K, pair):
    """n = 1：64 個 ±1 組態暴力加總，權重 Π(1 + K s_i s_j)。"""
    s = _spins(6).astype(float)
    W = np.ones(len(s))
    for i, j in LADDER[1]:
        W = W * (1 + K * s[:, i] * s[:, j])
    a, b = pair
    return float((W * s[:, a] * s[:, b]).sum() / W.sum())


def ladder_G_xy(K, pair, m=8):
    """n = 2：S = √2(cos θ, sin θ)，6 個角度各取 m 個等距點（被積函數是低次三角多項式，梯形法精確）。"""
    th = 2 * np.pi * np.arange(m) / m
    grids = np.meshgrid(*([th] * 6), indexing="ij", sparse=True)
    W = 1.0
    for i, j in LADDER[1]:
        W = W * (1 + 2 * K * np.cos(grids[i] - grids[j]))
    W = np.broadcast_to(W, (m,) * 6)
    a, b = pair
    obs = 2 * np.cos(grids[a]) * np.cos(grids[b])
    return float((W * obs).sum() / W.sum())


def ladder_saws(a, b):
    """梯子上從 a 到 b 的自迴避路徑長度（DFS）。"""
    nb = {i: [] for i in range(6)}
    for i, j in LADDER[1]:
        nb[i].append(j)
        nb[j].append(i)
    out = []

    def dfs(v, seen):
        if v == b:
            out.append(len(seen) - 1)
            return
        for w in nb[v]:
            if w not in seen:
                dfs(w, seen | {w})

    dfs(a, {a})
    return sorted(out)


def ladder_saw_sum(K, pair):
    return sum(K ** L for L in ladder_saws(*pair))


# ---------------------------------------------------------------- summary / check
def _f(x, nd=4):
    return f"{float(x):.{nd}f}"


def summary():
    S = {}
    # W901
    S["flory_p"] = [str(flory_p(d)) for d in (1, 2, 3, 4)]
    S["flory_nu"] = [str(SW.flory_nu(d)) for d in (1, 2, 3, 4)]
    S["flory_gap_3d"] = _f(flory_gap_3d(), 4)
    S["ms_table"] = {x: (_f(a), _f(b)) for x, (a, b) in ms_table().items()}
    S["beta_slopes"] = (_f(beta_slope(mf_m, 1.0)), _f(beta_slope(onsager_M, TC_2D)))
    S["TC_2D"] = _f(TC_2D)
    S["two_nu_eff"] = [_f(SW.two_nu_eff(N)) for N in SAW_N]
    # W902
    E0, E1, M0, M1 = rotation_check()
    S["rotation"] = (_f(E0), _f(E1), _f(np.linalg.norm(M0)), _f(np.linalg.norm(M1)), f"{angle_deg(M0):.1f}", f"{angle_deg(M1):.1f}")
    S["snapshot_up"] = {T: int((mc_snapshot(8, T) > 0).sum()) for T in SNAP_T}
    S["onsager_T1"] = _f(onsager_M(1.0))
    S["nu3_2"] = f"{2 * NU_KNOWN[3]:.3f}"
    S["ring4"] = {K: tuple(_f(v) for v in ring4_closed(K)) for K in TINY_K}
    S["pair"] = tuple(_f(v) for v in pair_closed(1.0))
    S["field_ring4"] = [_f(ising_exact(RING4, 1.0, H)["m"]) for H in FIELD_H]
    S["field_torus4"] = [_f(ising_exact(TORUS4, 1.0, H)["m"]) for H in FIELD_H]
    S["flip_tau"] = [f"{flip_tau(L):.0f}" for L in MC_L]
    S["trace_flips"] = {LT: trace_flips(*LT) for LT in MC_TRACES}
    S["trace_replica"] = {LT: trace_replica(*LT) for LT in MC_TRACES}
    hi = mc_run(8, 3.0)[0][:, 0]
    S["hi_trace"] = (_f(hi.mean()), _f(np.abs(hi).mean()))
    (am, ae), (em, ee) = mc_check_L4()
    S["mc_L4"] = (_f(am), _f(ae), _f(em), _f(ee))
    S["MS_T2"] = _f(onsager_M(2.0))
    # W903
    S["G1d"] = [_f(ising1d_G(r, K1D)) for r in range(5)]
    S["xi1d"] = _f(ising1d_xi(K1D))
    S["chi1d"] = _f(ising1d_chi(K1D))
    S["chi1d_partial"] = [_f(ising1d_chi_partial(K1D, R)) for R in (5, 10, 20)]
    S["chi_ring16"] = _f(ising1d_ring_closed(16, K1D)[1])
    S["chi_over_xi"] = [_f(ising1d_chi(K) / ising1d_xi(K)) for K in FISHER1D_K]
    S["xi2_moment_over_xi2"] = [_f(ising1d_xi2_moment(K) / ising1d_xi(K) ** 2) for K in FISHER1D_K]
    S["laws"] = {k: {kk: str(v) for kk, v in laws(e).items()} for k, e in EXPONENTS.items()}
    S["book3d"] = tuple(_f(v, 3) for v in book3d_alpha())
    S["book3d_laws"] = {k: _f(v, 3) for k, v in laws(BOOK_3D).items()}
    S["widom"] = tuple(str(v) for v in widom_exponents(*WIDOM_CASE))
    S["widom_slopes"] = tuple(_f(v) for v in widom_slopes(*WIDOM_CASE))
    # W904
    S["fourier_err"] = tuple(f"{v:.1e}" for v in fourier_check())
    S["landau"] = {r: (_f(landau_min(r)), _f(landau_F(landau_min(r), r))) for r in LANDAU_R}
    S["landau_h"] = _f(landau_min(-1.0, 1.0, LANDAU_H))
    S["ex94_closed"] = tuple(_f(v) for v in ex94_closed(-1.0))
    S["ex94_scaling"] = (_f(ex94_scaling(1.0, -1e-2)), _f(ex94_scaling(1.0, -1e-4)))
    S["ex94_slopes"] = tuple(_f(v) for v in ex94_slopes())
    S["zmf_m"] = [_f(zmf_avg(V, h=ZMF_H)) for V in ZMF_V]
    S["zmf_mstar"] = _f(landau_min(-1.0, 1.0, ZMF_H))
    S["zmf_absm"] = [_f(zmf_avg(V, what="absm")) for V in ZMF_V]
    S["zmf_f"] = [_f(zmf_avg(V, what="f")) for V in ZMF_V]
    # W905
    S["ZC_2D"] = _f(ZC_2D, 6)
    S["saw_xi2"] = [_f(saw_gf(z)[2]) for z in SAW_Z]
    S["saw_table"] = saw_table()
    S["coef_ratio"] = {str(s): [_f(coef_ratio(N, s)) for N in COEF_N] for s in (S_GAMMA, S_GAMMA2NU)}
    S["coef_ratio_book"] = [f"{coef_ratio_book(N, S_GAMMA, F(91, 5)):.1e}" for N in (10, 100)]
    S["r2_ratio"] = [_f(r2_ratio(N, S_GAMMA, F(3, 4), ZC_2D)) for N in COEF_N]
    S["zc_pow"] = tuple(_f(v) for v in zc_pow())
    S["r2_prefactor"] = (_f(r2_prefactor(S_GAMMA, F(3, 4), ZC_2D)), _f(r2_prefactor(S_GAMMA, F(3, 4), ZC_2D, book=True)))
    S["saw_amp"] = [_f(saw_amp(N)) for N in SAW_N]
    # W906
    S["hs"] = [_f(hs_rhs(s)) for s in (1, 2)]
    b72 = hs_book72()
    S["hs_book72"] = (_f(b72.real), _f(b72.imag))
    S["intersection"] = _f(intersection_factor())
    S["weak_saw"] = [f"{weak_saw_count(WEAK_N, u):.2f}" for u in WEAK_U]
    S["weak_saw_book"] = f"{weak_saw_count(WEAK_N, 1.0, sign=+1):.1e}"
    S["kernel_mass_book"] = (_f(kernel_mass_book(0.1, 1)), _f(kernel_mass_book(0.1, 3), 5))
    S["gauss_I"] = _f(gauss_I())
    S["gauss_mom"] = tuple(_f(gauss_moment(GAUSS_T, a, b)) for a, b in ((1, 1), (2, 0), (2, 2), (4, 0)))
    S["wick_count"] = [wick_count(m) for m in (1, 2, 3, 4)]
    T2 = ((2.0, -1.0), (-1.0, 3.0))
    C2 = np.linalg.inv(np.array(T2))
    S["wick_book"] = (_f(gauss_moment(T2, 2, 2)), _f(wick4(C2, 0, 0, 1, 1)), _f(wick4_book(C2, 0, 0, 1, 1)))
    S["trotter1"] = [f"{trotter_err(e, 1):.4g}" for e in TROTTER_EPS]
    S["trotter2"] = [f"{trotter_err(e, 2):.4g}" for e in TROTTER_EPS]
    S["trotter_global"] = [f"{trotter_global(M):.4g}" for M in TROTTER_M]
    e1 = [trotter_err(e, 1) for e in TROTTER_EPS]
    e2 = [trotter_err(e, 2) for e in TROTTER_EPS]
    S["trotter_ratios"] = ([f"{e1[i] / e1[i + 1]:.2f}" for i in range(2)], [f"{e2[i] / e2[i + 1]:.2f}" for i in range(2)])
    S["trotter_lead"] = [_f(e1[i] / trotter_lead(TROTTER_EPS[i])) for i in range(3)]
    S["ex95"] = tuple(f"{v:.1e}" for v in ex95_check())
    S["ladder"] = {p: [_f(ladder_G(1, LADDER_K, p)), _f(ladder_G(2, LADDER_K, p)), f"{ladder_G(0, LADDER_K, p):.5f}"]
                   for p in LADDER_PAIRS}
    S["ladder_saws"] = {p: ladder_saws(*p) for p in LADDER_PAIRS}
    return S


def check():
    import warnings
    warnings.simplefilter("ignore")
    # ------------------------------------------------ W901
    for d in (1, 2, 3, 4):
        assert flory_p(d) == 2 * SW.flory_nu(d) == F(6, d + 2)
    assert [str(flory_p(d)) for d in (1, 2, 3, 4)] == ["2", "3/2", "6/5", "1"]
    for d in (1, 2, 4):
        assert flory_p(d) == 2 * NU_KNOWN[d]
    assert NU_KNOWN[3] == SW.NU_QUOTED[3] and f"{flory_gap_3d():.4f}" == "0.0204"
    for T in (0.3, 0.5, 0.9, 0.99, 0.9999):
        assert abs(mf_m(T) - mf_m_newton(T)) < 1e-12
        assert abs(mf_m(T) - math.tanh(mf_m(T) / T)) < 1e-14
    assert abs(mf_m(1 - 1e-6) / mf_m_asym(1 - 1e-6) - 1) < 1e-5
    assert abs(math.sinh(2 / TC_2D) - 1) < 1e-14 and f"{TC_2D:.4f}" == "2.2692"
    u = math.exp(-4 / 1.0)
    assert abs(onsager_M(1.0) - onsager_lowT(1.0)) < 40 * u ** 4
    u = math.exp(-4 / 0.5)
    assert abs((onsager_M(0.5) - onsager_lowT(0.5)) / u ** 4 + 34) < 0.1      # 下一項 −34u⁴
    tab = ms_table()
    want = {0.5: ("0.9575", "0.9980"), 0.9: ("0.5254", "0.8959"), 0.99: ("0.1725", "0.6859")}
    for x in MS_X:
        assert (_f(tab[x][0]), _f(tab[x][1])) == want[x], (x, tab[x])
    assert _f(beta_slope(mf_m, 1.0)) == "0.5000" and _f(beta_slope(mf_m_newton, 1.0)) == "0.5000"
    assert _f(beta_slope(onsager_M, TC_2D)) == "0.1250"
    assert [_f(SW.two_nu_eff(N)) for N in SAW_N] == ["1.4442", "1.4505", "1.4556"]
    assert BETA_CITED == {"3d_ising": 0.326, "3d_saw": 0.302, "book_3d_saw": 0.308} == {**BETA_CITED, "book_3d_saw": BOOK_3D["beta"]}
    # ------------------------------------------------ W902
    for n in (2, 3, 4):
        E0, E1, M0, M1 = rotation_check(n=n)
        assert abs(E0 - E1) < 1e-12 and abs(np.linalg.norm(M0) - np.linalg.norm(M1)) < 1e-12
        S_, S2_ = rotation_spins(n=n)
        assert abs(on_energy(S_, ring_edges(ROT_NS)) + sum(float(np.dot(S_[i], S_[(i + 1) % ROT_NS])) for i in range(ROT_NS))) < 1e-12
    E0, E1, M0, M1 = rotation_check()
    assert _f(E0) == _f(E1) == "-6.3851" and _f(np.linalg.norm(M0)) == _f(np.linalg.norm(M1)) == "6.7943"
    assert (f"{angle_deg(M0):.1f}", f"{angle_deg(M1):.1f}") == ("-13.1", "106.9")
    assert abs((angle_deg(M1) - angle_deg(M0)) - 120) < 1e-9
    assert {T: int((mc_snapshot(8, T) > 0).sum()) for T in SNAP_T} == {1.5: 64, 3.5: 33}
    assert _f(onsager_M(1.0)) == "0.9993" and f"{2 * NU_KNOWN[3]:.3f}" == "1.176" and np.linalg.norm(M0 - M1) > 0.1
    for K in TINY_K:
        ex, cl = ising_exact(RING4, K), ring4_closed(K)
        assert abs(ex["Z"] - cl[0]) < 1e-9 * cl[0] and abs(ex["absm"] - cl[1]) < 1e-14 and abs(ex["m"]) < 1e-15
        assert abs(ex["PM"][4] - cl[2]) < 1e-14 and abs(ex["PM"][-4] - cl[2]) < 1e-14
        assert abs(ex["PM"][2] - cl[3]) < 1e-14 and abs(ex["PM"][0] - cl[4]) < 1e-14
    assert [tuple(_f(v) for v in ring4_closed(K)) for K in TINY_K] == [
        ("121.2329", "0.9337", "0.4504", "0.0330", "0.0333"), ("17.3497", "0.4871", "0.1283", "0.2306", "0.2823")]
    ex = ising_exact(PAIR, 1.0)
    assert (_f(ex["Z"]), _f(ex["absm"])) == tuple(_f(v) for v in pair_closed(1.0)) == ("6.1723", "0.8808")
    assert abs(ex["m"]) < 1e-15
    fr = [ising_exact(RING4, 1.0, H)["m"] for H in FIELD_H]
    ft = [ising_exact(TORUS4, 1.0, H)["m"] for H in FIELD_H]
    assert [_f(v) for v in fr] == ["0.3506", "0.0367"] and [_f(v) for v in ft] == ["0.9210", "0.1584"]
    for H in FIELD_H:            # 第二種算法：⟨M⟩ = ∂ ln Z/∂H（差分）
        for g in (RING4, TORUS4):
            dH = 1e-6
            d = (math.log(ising_exact(g, 1.0, H + dH)["Z"]) - math.log(ising_exact(g, 1.0, H - dH)["Z"])) / (2 * dH)
            assert abs(d / g[0] - ising_exact(g, 1.0, H)["m"]) < 1e-7
    # 線性區：H 小時 ⟨M⟩ ≈ ⟨M²⟩·H（漲落–響應）
    assert abs(ising_exact(RING4, 1.0, 1e-5)["m"] * 4 / (ising_exact(RING4, 1.0)["M2"] * 1e-5) - 1) < 1e-6
    taus = [flip_tau(L) for L in MC_L]
    assert [f"{t:.0f}" for t in taus] == ["133", "608", "2151"], taus
    assert taus[0] < taus[1] < taus[2] and taus[1] / taus[0] > 3 and taus[2] / taus[1] > 3
    (am, ae), (em, ee) = mc_check_L4()            # 取樣器 vs 精確列舉
    assert abs(am - ae) < 0.01 and abs(em - ee) < 0.01, (am, ae, em, ee)
    assert (_f(am), _f(ae), _f(em), _f(ee)) == ("0.9182", "0.9189", "-1.7528", "-1.7554")
    assert {LT: trace_flips(*LT) for LT in MC_TRACES[1:]} == {(4, 2.0): 30, (8, 2.0): 2}
    assert {LT: trace_replica(*LT) for LT in MC_TRACES} == {(8, 3.0): 0, (4, 2.0): 25, (8, 2.0): 11}
    hi = mc_run(8, 3.0)[0][:, 0]
    assert (_f(hi.mean()), _f(np.abs(hi).mean())) == ("0.0229", "0.3470")      # 高溫：在 0 附近亂跳
    assert len(mc_trace(8, 2.0)[0]) == 400 and mc_trace(8, 2.0)[0][-1] == MC_SWEEPS
    assert _f(onsager_M(2.0)) == "0.9113"
    # ------------------------------------------------ W903
    K = K1D
    assert [_f(ising1d_G(r, K)) for r in range(5)] == ["1.0000", "0.7616", "0.5800", "0.4417", "0.3364"]
    for r in range(8):
        assert abs(ising1d_G(r, K) - ising1d_transfer_G(r, K)) < 1e-13
    ring = ising_exact((12, ring_edges(12)), K)
    G12, chi12 = ising1d_ring_closed(12, K)
    assert max(abs(a - b) for a, b in zip(ring["corr0"], G12)) < 1e-12 and abs(ring["M2"] / 12 - chi12) < 1e-12
    assert _f(ising1d_xi(K)) == "3.6719" and _f(ising1d_chi(K)) == _f(math.exp(2 * K)) == "7.3891"
    assert [_f(ising1d_chi_partial(K, R)) for R in (5, 10, 20)] == ["5.7520", "6.9696", "7.3615"]
    assert abs(ising1d_chi_partial(K, 200) - ising1d_chi(K)) < 1e-12
    assert _f(ising1d_ring_closed(16, K)[1]) == "7.2021"
    assert abs(ising_exact((16, ring_edges(16)), K)["M2"] / 16 - ising1d_ring_closed(16, K)[1]) < 1e-11   # 2^16 個組態
    assert [_f(ising1d_chi(k) / ising1d_xi(k)) for k in FISHER1D_K] == ["2.0123", "2.0002", "2.0000"]
    for k in (1.0, 2.0):
        assert abs(ising1d_xi2_moment(k) - ising1d_xi2_moment_sum(k)) < 1e-8 * ising1d_xi2_moment(k)
    xr = [ising1d_xi2_moment(k) / ising1d_xi(k) ** 2 for k in FISHER1D_K]
    assert [_f(v) for v in xr] == ["1.9877", "1.9998", "2.0000"]             # → 2：Σr²G/ΣG ∝ ξ²（不是 ξ）
    assert ising1d_xi2_moment(3.0) / ising1d_xi(3.0) > 100                  # 若是 ξ 就該是常數
    for name, e in EXPONENTS.items():
        L = laws(e)
        assert L["rushbrooke"] == 2 and L["fisher"] == e["gamma"] and L["widom"] == e["gamma"]
        assert L["hyper_alpha"] == e["alpha"] and L["hyper_delta"] == e["delta"]
        assert from_nu_eta(e["nu"], e["eta"], e["d"]) == (e["alpha"], e["beta"], e["gamma"], e["delta"])
        b, g, de = widom_exponents(e["alpha"], e["beta"] + e["gamma"])
        assert (b, g, de) == (e["beta"], e["gamma"], e["delta"])
    assert 2 - 3 * EXPONENTS["mf"]["nu"] == F(1, 2)                         # 平均場在 d = 3 不滿足 hyperscaling
    sw = EXPONENTS["saw2d"]
    assert sw["gamma"] == WS.GAMMA_QUOTED and sw["eta"] == WS.ETA_2D and sw["nu"] == WS.NU_2D
    assert WS.gamma_exp(sw["eta"], sw["nu"]) == sw["gamma"]
    assert BOOK_3D == dict(alpha=0.12, beta=0.308, gamma=1.19, delta=4.865, nu=0.602, eta=0.023, d=3)
    a_r, a_h, rush = book3d_alpha()
    assert (_f(a_r, 3), _f(a_h, 3), _f(rush, 3)) == ("0.194", "0.194", "1.926")
    L3 = laws(BOOK_3D)
    assert abs(L3["fisher"] - BOOK_3D["gamma"]) < 1e-3 and abs(L3["hyper_delta"] - BOOK_3D["delta"]) < 1e-3
    assert abs(L3["widom"] - BOOK_3D["gamma"]) < 2e-3 and abs(from_nu_eta(0.602, 0.023, 3)[1] - 0.308) < 1e-3
    assert abs(BOOK_3D["alpha"] - a_h) > 0.07                                 # 書上 α = 0.12 與同表其他值不合
    assert widom_exponents(*WIDOM_CASE) == (F(1, 8), F(7, 4), F(15))
    for case in (WIDOM_CASE, (F(0), F(3, 2)), (F(1, 2), F(91, 64))):
        bx, gx, dx = widom_exponents(*case)
        for meth in ("exact", "fd"):
            b, g, di = widom_slopes(*case, method=meth)
            assert abs(b - bx) < 1e-6 and abs(g - gx) < 1e-5 and abs(di - 1 / dx) < 1e-6, (case, meth, b, g, di)
    assert tuple(_f(v) for v in widom_slopes(*WIDOM_CASE)) == ("0.1250", "1.7500", "0.0667")
    # ------------------------------------------------ W904
    assert max(fourier_check()) < 1e-12
    for r in LANDAU_R:
        m, Fm = landau_closed(r)
        assert abs(landau_min(r) - m) < 1e-12 and abs(landau_F(landau_min(r), r) - Fm) < 1e-12
    xs = np.linspace(-2, 2, 400001)          # 第二種算法：格點上直接找最小
    for r in LANDAU_R:
        assert abs(abs(xs[np.argmin(landau_F(xs, r))]) - landau_closed(r)[0]) < 2e-5
    assert abs(xs[np.argmin(landau_F(xs, -1.0, 1.0, LANDAU_H))] - landau_min(-1.0, 1.0, LANDAU_H)) < 2e-5
    assert _f(landau_min(-1.0, 1.0, LANDAU_H)) == "1.0467"
    assert (_f(landau_min(-1.0)), _f(landau_F(1.0, -1.0))) == ("1.0000", "-0.2500")
    for r in (-1.0, -0.01):
        assert abs(ex94_m(r, 0) - ex94_closed(r)[0]) < 1e-12
        assert abs(landau_F(ex94_m(r, 0), r, EX94_U) - ex94_closed(r)[1]) < 1e-14
    assert tuple(_f(v) for v in ex94_closed(-1.0)) == ("1.7321", "-0.7500")
    sc = (ex94_scaling(1.0, -1e-2), ex94_scaling(1.0, -1e-4))
    assert abs(sc[0] - sc[1]) < 1e-10 and _f(sc[0]) == "2.1038"
    assert abs(sc[0] ** 3 / 3 - sc[0] - 1) < 1e-12                       # ℳ³/3 − ℳ = x
    assert abs(ex94_m(0.0, 1e-3) - (3e-3) ** (1 / 3)) < 1e-12
    chi_lo, chi_hi = 1 / (-0.01 + ex94_m(-0.01, 0) ** 2), 1 / (0.01 + ex94_m(0.01, 0) ** 2)
    assert abs(chi_lo - 50) < 1e-9 and abs(chi_hi - 100) < 1e-9                    # χ = 1/(2|r|)、1/r：振幅比 2
    assert tuple(_f(v) for v in ex94_slopes()) == ("0.5000", "1.5000", "0.0000", "0.3333", "1.0000")
    assert widom_exponents(0, F(3, 2)) == (F(1, 2), F(1), F(3))          # Ex. 9.4 的 α、Δ → β、γ、δ
    for V in ZMF_V:
        for what, kw in (("m", dict(h=ZMF_H)), ("absm", {}), ("f", {})):
            a, b = zmf_avg(V, what=what, **kw), zmf_avg(V, what=what, method="grid", **kw)
            assert abs(a - b) < 1e-8, (V, what, a, b)
        assert abs(zmf_avg(V)) < 1e-12                                    # h = 0：⟨m⟩ = 0（對稱）
    assert [_f(zmf_avg(V, h=ZMF_H)) for V in ZMF_V] == ["0.0869", "0.7576", "1.0042"]
    assert _f(landau_min(-1.0, 1.0, ZMF_H)) == "1.0050"
    assert [_f(zmf_avg(V, what="absm")) for V in ZMF_V] == ["0.8888", "0.9922", "0.9992"]
    assert [_f(zmf_avg(V, what="f")) for V in ZMF_V] == ["-0.2719", "-0.2397", "-0.2478"]
    # 兩個井的近似：⟨m⟩_V ≈ m* tanh(V h m*)
    ms_ = landau_min(-1.0, 1.0, ZMF_H)
    assert abs(zmf_avg(1000, h=ZMF_H) - ms_ * math.tanh(1000 * ZMF_H * ms_)) < 1e-3
    # ------------------------------------------------ W905
    for z in SAW_Z:
        chi, num, xi2 = saw_gf(z, 8)
        chi2, xi2b = saw_gf_realspace(z, 8)
        assert abs(chi - chi2) < 1e-12 * chi and abs(xi2 - xi2b) < 1e-12 * xi2
    assert [_f(saw_gf(z)[2]) for z in SAW_Z] == ["0.6347", "1.2799", "2.4593"]
    assert abs(saw_gf(0.1, 12)[2] - saw_gf(0.1, 14)[2]) < 1e-5            # z = 0.1 已收斂
    assert saw_table() == [(1, 4, 4), (2, 12, 32), (3, 36, 164), (4, 100, 704), (5, 284, 2716), (6, 780, 9808)]
    assert f"{ZC_2D:.4f}" == "0.3791"
    for s in (S_GAMMA, S_GAMMA2NU, F(1, 2), F(3)):
        for N in (0, 1, 7, 50, 1000):
            assert abs(coef_log(N, s, ZC_2D) - coef_log_rec(N, s, ZC_2D)) < 1e-9
        assert abs(coef_exact(30, s, ZC_2D) / coef_asym(30, s, ZC_2D) - coef_ratio(30, s)) < 1e-12
    assert abs(coef_exact(5, 3, 1.0) - math.comb(7, 2)) < 1e-9                      # (1−z)^{−3}：C(N+2, 2)
    assert {str(s): [_f(coef_ratio(N, s)) for N in COEF_N] for s in (S_GAMMA, S_GAMMA2NU)} == {
        "43/32": ["1.0227", "1.0023", "1.0002"], "91/32": ["1.2760", "1.0264", "1.0026"]}
    assert coef_ratio_book(100, S_GAMMA, F(91, 5)) < 1e-20 and coef_ratio_book(10, S_GAMMA, F(3)) < 0.1
    assert [f"{coef_ratio_book(N, S_GAMMA, F(91, 5)):.1e}" for N in (10, 100)] == ["1.4e-17", "1.9e-34"]
    assert [_f(r2_ratio(N, S_GAMMA, F(3, 4), ZC_2D)) for N in COEF_N] == ["1.2476", "1.0240", "1.0024"]
    for N in COEF_N:
        assert abs(r2_ratio(N, S_GAMMA, F(3, 4), ZC_2D) - r2_ratio(N, S_GAMMA, F(3, 4), ZC_2D, rec=True)) < 1e-10
    pw = zc_pow()
    assert (_f(pw[0]), _f(pw[1])) == ("4.2850", "3.6823")
    assert (_f(r2_prefactor(S_GAMMA, F(3, 4), ZC_2D)), _f(r2_prefactor(S_GAMMA, F(3, 4), ZC_2D, book=True))) == ("2.1962", "1.8873")
    assert abs(r2_prefactor(S_GAMMA, F(3, 4), ZC_2D) / r2_prefactor(S_GAMMA, F(3, 4), ZC_2D, book=True)
               - pw[0] / pw[1]) < 1e-12
    # 書上的 z_c^γ：用 z_c = 1 與 z_c = 1/μ 算出來的比值不一樣（前因子必須是 z_c^{−2ν}）
    for zc in (1.0, ZC_2D, 0.25):
        assert abs(r2_ratio(10 ** 5, S_GAMMA, F(3, 4), zc) - 1) < 1e-4
        book = r2_ratio(10 ** 5, S_GAMMA, F(3, 4), zc) * r2_prefactor(S_GAMMA, F(3, 4), zc) / r2_prefactor(S_GAMMA, F(3, 4), zc, book=True)
        assert (abs(book - 1) < 1e-4) == (zc == 1.0)                         # 書上的 z_c^γ 只在 z_c = 1 時碰巧對
    # K₁：若 ξ = K₁(z_c−z)^{−ν}，則 ξ²χ 的係數是 K₁²K₂（K₁ = 2 時差 2 倍）
    assert abs(2 ** 2 * coef_exact(50, S_GAMMA2NU, ZC_2D) / (2 * coef_exact(50, S_GAMMA2NU, ZC_2D)) - 2) < 1e-12
    # VΓ_N/z_0^N：Σ Γ_N (z/z_0)^N 的 z^N 係數（z_0 = 1/4）
    c14 = SW.saw_enum(14)[0]
    zz = np.polynomial.polynomial.Polynomial([c14[N] * 4.0 ** N for N in range(15)])
    assert abs(zz.coef[10] - c14[10] / 0.25 ** 10) < 1e-3
    assert [_f(saw_amp(N)) for N in SAW_N] == ["0.8299", "0.8224", "0.8168"]
    # ------------------------------------------------ W906
    for s in (0.0, 0.5, 1.0, 2.0, 3.0):
        for u in (0.5, 1.0, 2.0):
            assert abs(hs_rhs(s, u) - hs_lhs(s, u)) < 1e-12
    assert [_f(hs_rhs(s)) for s in (1, 2)] == ["0.7788", "0.3679"]
    b1, b2 = hs_book72(), hs_book72(method="closed")
    assert abs(b1 - b2) < 1e-12 and (_f(b1.real), _f(b1.imag)) == ("0.3894", "-0.2395")
    for k in (1, 2, 3):
        v = integrate.quad(lambda p: math.exp(-p * p) * math.cos(k * p), -np.inf, np.inf)[0] / math.sqrt(math.pi)
        assert abs(v - visit_weight(k)) < 1e-12
        assert pair_sum(k) == k + 2 * math.comb(k, 2)
        assert abs(visit_weight(k) - visit_weight(1) ** k * intersection_factor() ** math.comb(k, 2)) < 1e-15
    assert _f(intersection_factor()) == _f(math.exp(-0.5)) == "0.6065"
    wk = [weak_saw_count(WEAK_N, u) for u in WEAK_U]
    for u in WEAK_U[:3]:
        assert abs(weak_saw_count(WEAK_N, u) - SW.weighted_count(WEAK_N, 1 - math.exp(-u / 2))) < 1e-6
    assert [f"{v:.2f}" for v in wk] == ["65536.00", "23695.11", "7988.15", "5916.00"]
    assert wk[0] == 4 ** WEAK_N and wk[-1] == SW.saw_c(WEAK_N) == SW.SAW_KNOWN_C[WEAK_N]
    assert abs(weak_saw_count(WEAK_N, 60.0) - SW.saw_c(WEAK_N)) < 1e-6
    assert f"{weak_saw_count(WEAK_N, 1.0, sign=+1):.1e}" == "6.2e+05"
    for x in (0.0, 0.3, 1.0):
        assert abs(heat_kernel(x, 0.1) - heat_kernel_fourier(x, 0.1)) < 1e-10
    assert abs(integrate.quad(lambda x: heat_kernel(x, 0.1), -np.inf, np.inf)[0] - 1) < 1e-12
    assert abs(integrate.quad(lambda x: x * x * heat_kernel(x, 0.1), -np.inf, np.inf)[0] - 0.1) < 1e-12
    for d in (1, 2, 3):
        assert abs(kernel_mass_book(0.1, d) - (0.1 / d) ** (d / 2)) < 1e-12
    assert (_f(kernel_mass_book(0.1, 1)), _f(kernel_mass_book(0.1, 3), 5)) == ("0.3162", "0.00609")
    gI = gauss_I()
    assert abs(gI - gauss_I_quad()) < 1e-8 and abs(gI - 2 * math.pi / math.sqrt(3) * math.exp(1 / 3)) < 1e-12
    assert _f(gI) == "5.0627"
    assert abs(gauss_I(((2.0, -1.0), (-1.0, 3.0)), (0.5, -1.0)) - gauss_I_quad(((2.0, -1.0), (-1.0, 3.0)), (0.5, -1.0))) < 1e-8
    C = np.linalg.inv(np.array(GAUSS_T))
    mom = [gauss_moment(GAUSS_T, a, b) for a, b in ((1, 1), (2, 0), (2, 2), (4, 0))]
    assert max(abs(m - v) for m, v in zip(mom, (1 / 3, 2 / 3, 2 / 3, 4 / 3))) < 1e-12
    assert abs(mom[2] - wick4(C, 0, 0, 1, 1)) < 1e-12 and abs(mom[3] - wick4(C, 0, 0, 0, 0)) < 1e-12
    assert abs(gauss_moment(GAUSS_T, 3, 1) - wick4(C, 0, 0, 0, 1)) < 1e-12
    T2 = ((2.0, -1.0), (-1.0, 3.0))
    C2 = np.linalg.inv(np.array(T2))
    m22 = gauss_moment(T2, 2, 2)
    assert abs(m22 - wick4(C2, 0, 0, 1, 1)) < 1e-12 and abs(m22 - 0.32) < 1e-12
    assert abs(wick4_book(C2, 0, 0, 1, 1) - 0.44) < 1e-12                       # 書上 (9:S-19) 第三行
    assert [wick_count(m) for m in (1, 2, 3, 4)] == [wick_count_enum(m) for m in (1, 2, 3, 4)] == [1, 3, 15, 105]
    assert max(ex95_check()) < 1e-10 and max(ex95_check(0.3)) < 1e-10                     # Ex. 9.5
    e1 = [trotter_err(e, 1) for e in TROTTER_EPS]
    e2 = [trotter_err(e, 2) for e in TROTTER_EPS]
    assert [f"{v:.4g}" for v in e1] == ["0.005784", "0.001502", "0.0003827"], e1
    assert [f"{v:.4g}" for v in e2] == ["0.0009465", "0.000124", "1.589e-05"], e2
    assert [f"{e1[i] / e1[i + 1]:.2f}" for i in range(2)] == ["3.85", "3.92"]
    assert [f"{e2[i] / e2[i + 1]:.2f}" for i in range(2)] == ["7.63", "7.81"]
    assert abs(e1[-1] / trotter_lead(TROTTER_EPS[-1]) - 1) < 0.05
    assert [_f(e1[i] / trotter_lead(TROTTER_EPS[i])) for i in range(3)] == ["0.9268", "0.9627", "0.9811"]
    eg = [trotter_global(M) for M in TROTTER_M]
    assert [f"{v:.4g}" for v in eg] == ["0.02366", "0.0117", "0.005821"] and 1.95 < eg[1] / eg[2] < 2.05
    for p in LADDER_PAIRS:
        for Kl in (0.3, LADDER_K, 0.9):
            assert abs(ladder_G(1, Kl, p) - ladder_G_ising(Kl, p)) < 1e-13
            assert abs(ladder_G(2, Kl, p) - ladder_G_xy(Kl, p)) < 1e-12
            assert abs(ladder_G(0, Kl, p) - ladder_saw_sum(Kl, p)) < 1e-15
    assert ladder_saws(0, 3) == [1, 3, 5] and ladder_saws(0, 5) == [3, 3, 3, 5]
    assert {p: [_f(ladder_G(1, LADDER_K, p)), _f(ladder_G(2, LADDER_K, p)), f"{ladder_G(0, LADDER_K, p):.5f}"]
            for p in LADDER_PAIRS} == {(0, 3): ["0.6027", "0.5610", "0.65625"], (0, 5): ["0.3562", "0.3171", "0.40625"]}
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
    print(f"walk_crit（第 9 章）：平均場與 Onsager 的 β、小 Ising 的精確列舉與 Metropolis 翻轉時間、一維 Ising 的 G、ξ、χ、"
          f"標度律與 Widom 形式、Landau 雙井與 (9.34)、(z_c−z)^(−s) 的係數與 ⟨R²⟩ ∝ N^(2ν)、輔助場與高斯積分、Wick、Trotter、"
          f"相交權重與 n → 0 的自迴避全部驗證通過（{time.time() - t0:.1f} 秒）")
