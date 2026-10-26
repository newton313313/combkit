# 由 combkit/tools/sync_from_walk017.py 從 walk017 專案的 codes/utils/walk_cond.py 複製；請改 walk017 那份再同步。
"""
walk_cond.py —— 《Elements of the Random Walk》第 10 章 Walks and the O(n) model: mean field theory and spin waves 的數學
（Python＋numpy／scipy（integrate、optimize、special、linalg），不 import manim；用到同資料夾的 walk_saw（無限平面的 c_N、全部 4^N 條漫步、
(7.57) 的 G_MF）、walk_crit（z_c = 1/μ、指數表）、walk_rec（sc 的 G_c、回家機率）、walk_gf（sc 的 C(N; 0)）——在 walk017 是 utils 套件、
在 combkit 是 combkit 套件，都用相對 import；直接執行本檔時退回絕對 import）

約定（書上式號、筆誤與〔注〕見 ref/md/ch10_formulas.md）：
  這一章的 r 定義成 e^{−r} = z/z_c（(10.4)、(10.49)：G = Σ_N C(N) e^{−rN}；第 9 章是 e^{−rN/2}，差 2 倍）；z 是每步的逸度。
  環面（torus）：L×L 正方格子、週期邊界，格點 s = x + L·y（0 ≤ x, y < L），V = L²；從格點 0 出發的自迴避漫步條數 c_N（N = 0..V−1）。
  畫環面上的漫步：格點 (x, y) 畫在格子中心，盒子是 [−½, L−½]²；越界的一步拆成兩個半段（走出去半格、從對邊進來半格）。
  平均場：f(λ) = e^{−N0(rλ/2 + uλ²/4)}（(10.15)），I(r) = ∫_0^∞ f dλ；Γ_MF = (N0/2)·I（更正；書上寫 N0/4），V = N0（格距 a = 1）。
  自旋波：λ1 = r/2 + q²/2 + u|M|²/2（橫向，n−1 重）、λ2 = r/2 + q²/2 + 3u|M|²/2（縱向，1 重）；平均場的鞍點 |M|² = |r|/u（r < 0）。
  相關洞：三維 L³ 週期格子、格子 Laplace q̂² = Σ 4 sin²(q_i/2)、去掉 q = 0；D(Δ) = (1/V)Σ_{q≠0} e^{iq·Δ}[1/(r+q̂²+u|M|²) − 1/(r+q̂²+3u|M|²)]；
  r < 0 的鞍點上方括號 = 2|r|/(q̂²(q̂²+2|r|))。sc（簡單立方）普通漫步每步 1/6；P_n(R) = n 步後在 R 的機率；「離起點 R 的步數」= Σ_{n≤N} P_n(R)。
  分數用 Fraction；回傳快取（lru_cache）陣列、tuple 或 dict 的函數不要就地修改結果。

§10.1 關在盒子裡：環面上的自迴避漫步 (10.1)，Fig. 10.1                                           [W1001]
    TORI = (3, 4, 5)；torus_neighbors(L)：每個格點的 (E, N, W, S) 鄰居
    torus_counts(L)               c_0..c_{V−1}（位元遮罩 DFS、8 重對稱、快取；L = 5 約 1 秒，其餘瞬間）
                                  3×3：1, 4, 12, 32, 80, 160, 256, 296, 168（共 1009）
                                  4×4：1, 4, 12, 36, 96, 264, 624, 1536, 3072, 6240, 9840, 15912, 17472, 19200, 10656, 5712（共 90677）
                                  5×5：…, c_24 = 236456（共 22094601）
    torus_counts_plain(L, Nmax=None, sym=1)   不用對稱（sym=1）或只固定第一步（sym=4）的 DFS（第二種算法）
    torus_counts_brute(L, N)      由 walk_saw.all_walks_xy 的全部 4^N 條，數「各點模 L 互不相同」的條數（第三種算法，N ≤ 8）
    SAW_INF = walk_saw.SAW_KNOWN_C（無限平面 c_0..c_16）；torus_vs_inf(L)：[(N, c_N(環面), c_N(無限))]（N ≤ min(V−1, 16)）
    torus_first_diff(L) = (L, 4)：N < L 完全一樣（4、12、36…）；N = L 時剛好少 4 條（直走一圈回到起點）
    GRAY_Q4 = 91392               4×4 環面＝四維超立方體：16 × c_15 = 91392（四維 Gray 碼的條數，OEIS A091299）
    torus_walks(L, N)             從格點 0 出發、長 N 的全部環面自迴避漫步（格點 tuple；L ≤ 4）
    torus_steps(path, L)、torus_lift(path, L, start)、torus_wraps(path, L, start)、torus_segments(path, L, start)
                                  步伐字串（E N W S）、攤開到無限平面的座標、越界次數、畫圖用線段 [((x0,y0), (x1,y1), 越界?)]
    torus_example(name, L=4, start=(1, 1), seed=1004)   卡上的三條（dict：N、sites、steps、lift、wraps、segments、start、L；快取，不要修改）：
                                  'short' N = 3 'WNE'（不越界，20 條候選）、'mid' N = 8 'WSSEENNE'（越界 2 次，1352 條候選）、
                                  'full' N = 15 'WSSEEENNNEEESSW'（越界 3 次，走滿 16 格；5712 條候選）；起點畫在 (1, 1)
    torus_mean_N(L, z)            ⟨N⟩（4×4：z = 0.2、z_c、1.0 → 1.4400、4.8899、11.7409；最多 15）
    book101_mean_N(L, z)          照書上 (10.1) 的 z^n（n → 0）：Σ_N c_N 與 z 無關（4×4 永遠是 90677）→ ⟨N⟩ = 0
§10.2 O(n) 模型的平均場 (10.8)–(10.15)，註 2–3                                                     [W1002]
    MF_U = 1、MF_N0 = 100、MF_R = 0.5（r < 0 的例子用 −0.5）
    mf_f(lam, r, u, N0)；mf_I(r, u, N0, method='erfcx'|'quad')；mf_lnI(r, u, N0)（大 N0 用對數）
    gamma_mf(r, u, N0, book=False) = (N0/2)·I（book=True 照書上 N0/4）；ln_gamma_mf(r, u, N0)
    gamma_n(n, r, u, N0, method='quad'|'pcf')   有限 n 的 N0⟨M1²⟩_n（徑向積分；'pcf' 用拋物柱函數 D_v，第二種算法）；
                                  gamma_n1_direct(r, u, N0)：n = 1 直接一維積分（第三種算法）
                                  MF_NS = (3, 2, 1, 0.1, 0.01)：1.7183、1.7634、1.8131、1.8624、1.8676 → n → 0：1.8682 ＝ (N0/2)I；書上 N0/4 得 0.9341
                                  u = 0（高斯）時任何 n 都是 1/r = 2（書上 1/(2r) = 1）：gamma_gauss_check()
                                  MF_NS_NEG = (0.01, 0.001, 1e−4, 1e−6)，r = −0.5：Γ_n/Γ_MF = 0.3299、0.8311、0.9801、0.9998（Γ_MF = 9179.7）：
                                  r < 0 要 n·e^{N0r²/4u} ≪ 1 才到極限（偏差 ∝ n）——必須先取 n → 0（p. 234 的極限順序）
    s_integral(lam, eps=0.1)      ∫dS e^{iSλ}/(S − iε)（註 2 的寫法，QAWF 數值）：λ = 1 → 5.6853i ＝ 2πi e^{−ε}；λ = −1 → 0
                                  （書上說「對 S 的積分是 iπ」：那是主值，λ < 0 時是 −iπ 而不是 0）
    mf_asym_I(r, u, N0)、gamma_mf_asym(r, u, N0)   r > 0：I ≈ 2/(N0 r)、Γ ≈ 1/r；r < 0：I ≈ 2√(π/(N0u)) e^{N0r²/4u}、Γ ≈ √(πN0/u) e^{N0r²/4u}
                                  MF_N0S = (100, 1000, 10000)：r = 0.5：Γ_MF = 1.8682、1.9844、1.9984 → 2；
                                  r = −0.5：精確 ÷ 漸近 = 0.9998、1.0000、1.0000；ln Γ_MF / N0 = 0.0912、0.0665、0.0630 → r²/4u = 0.0625
    fn3_check(r, u, N0)           註 3 的對照：(更正的 Γ_MF, N0·G_MF(1−r, N0, u), 同上 quad, 書上 (N0/4)I, N0·G_MF(1−2r, N0, 4u))；
                                  r = 0.5：1.8682、1.8682、1.8682、0.9341、0.9341（x = z/z_c = 1/2 就是第 7 章 Ex. 7.5 的 x）
§10.3 自旋波 (10.16)–(10.36)，Ex. 10.1–10.4                                                        [W1003]
    wall_M(x, r=−1, u=1) = √(|r|/u) tanh(x√(|r|/2))（(10.17) 在 M(0) = 0 的牆邊的解；Ex. 10.1 的例子）：
                                  WALL_X = (1, 2, 4) → 0.6089、0.8884、0.9930；癒合長度 √(2/|r|) = 1.4142
    wall_bvp(r, u, X=20, n=2001)  有限差分＋牛頓解邊值問題 → (x, M)；wall_errors() → (BVP 與 tanh 的最大差 2.8e−06, tanh 代回 (10.17) 的殘差 8.5e−08)
    ex101_grad_check(seed=1003)   離散化 ℋ 的數值梯度 ＝ 格距 × [−∇²M + rM + uM|M|² − h⃗]（(10.17)；回傳相對誤差 < 1e−7）
    ex102_coeffs(seed=1003)       |M+εσ| 的多項式係數（數值內插 vs 公式）：三階 u(M·σ)|σ|²、四階 u|σ|⁴/4（Ex. 10.2；誤差）
    SW_R = −1、SW_U = 1、SW_Q2 = 1、SW_M = (2/3, 2/3, 1/3)（|M|² = 1 = |r|/u）：sw_lambdas() = (1/2, 3/2)
    sw_X(M, r, q2, u)、sw_Xinv(M, r, q2, u)（(10.27)、(10.30)；Fraction 或 float 都可以）：
                                  X = (1/18)[[17,8,4],[8,17,4],[4,4,11]]、X⁻¹ = (1/27)[[38,−16,−8],[−16,38,−8],[−8,−8,50]]、X X⁻¹ = I（精確）
    sw_check(ns=(1..6), seed=1003)   隨機 M、r、q²、u：特徵值、重數 (n−1, 1)、(10.31) 左右逆、det、trace 的最大誤差（< 1e−12）
    mode_table(ns=(3, 2, 1, 0))   [(n, 橫向重數 n−1, 縱向 1, det X = λ1^{n−1}λ2, (det X)^{−1/2})]：
                                  det = 3/8、3/4、3/2、3；(det)^{−1/2} = 1.6330、1.1547、0.8165、0.5774（n → 0 ＝ √(λ1/λ2) = 1/√3；書上 √(λ2/λ1) = 1.7321）
    gauss_n2_check()              n = 2 的高斯積分（(10.32)）：數值 3.9891 ＝ π/√det·e^{+h·X⁻¹·h/4}；書上的 −1/4 得 3.2988（M = (3/5, 4/5)、h = (1/2, 0)）
§10.3 後半–§10.4 相關洞 (10.37)–(10.48)、(10.76)–(10.84)，Fig. 10.2–10.3，Ex. 10.5                  [W1004]
    HOLE_L = 64、HOLE_R = −1（u = 1、|M|² = 1、2|r| = 2）
    D_box(L, r, u=1)              三維週期盒子的 D（FFT；去掉 q = 0；快取、唯讀）；D_axis(L, r)：沿 x 軸 Δ = 0..L/2
                                  HOLE_D = (0, 1, 2, 4, 8, 16, 32) → 0.1083、0.0613、0.0359、0.0167、0.0065、0.0016、−0.0001（嚴格遞減）；
                                  1 − D = 0.8917、0.9387、0.9641、0.9833、0.9935、0.9984、1.0001（洞：頭尾靠近時少 11%）
    D_direct(L, r, delta)         直接對 q 求和（第二種算法）；D_inf(delta, r)：無限大格子（Bessel 積分，第三種算法）
    SC_EWALD = 2.837297；box_shift(L, r, delta) = −2.837297/(4πL) + Δ²/(6L³) + 1/(L³·2|r|)：D_box − D_inf 與它相差 < 5e−6
    D_book46(L, r)                照書上 (10.46)（分子少了 2|r|）＝ D_box/(2|r|)：Δ = 0 → 0.0541（應為 0.1083）
    coulomb_fit(L, r, dmin=4, dmax=16) → (c0, c)：−D(Δ) ≈ c0 − c/Δ：c0 = 0.0035、c = 0.0801（1/(4π) = 0.0796，比 1.0064）
    D_cont(rho, r)、D_cont_quad(rho, r)   Ex. 10.5 連續版 (1/(4πρ))(1 − e^{−√(2|r|)ρ})：ρ = 1, 2, 4, 8 → 0.0602、0.0374、0.0198、0.0099
    case1_check(L=32, r=1)        (10.44) r > 0：盒子裡的 G 與無限大格子相差（Δ ≤ 4 的最大相對差 1.2e−11）
    sc_prob(n, R)                 sc 普通漫步 P_n(R)（R 是整數＝(R, 0, 0) 或三元組；三個軸分步數＋兩枚硬幣，精確）；
                                  sc_visits_dp(N, R) → (總和, [P_0..P_N])（三維格點逐步遞迴，第二種算法，N ≤ 40）
    sc_visits(N, R)               Σ_{n≤N} P_n(R)；VIS_N = 100、VIS_R = (1, 2, 4, 8, 16)：0.4509、0.1931、0.0610、0.0103、0.000183
    erfc_form(R, N) = (3/(2πR)) erfc(R√(3/(2(N+1))))（(10.83) 更正版，sc 格子 z_c ↔ 6）：0.4121、0.1744、0.0586、0.0100、0.000174；
                                  精確 ÷ erfc 形 = 1.094、1.108、1.042、1.028、1.052
    erfc_form_book(R, N) = (3/(2πR)) erfc(3R²/(N+1))（照書上的引數）：R = 8 → 4.3e−04、R = 16 → 1.7e−28（完全不對）
    eq1080_check()、eq1083_check(R, a)   (10.80) 的係數（級數＝書上左邊×z_c^{−N}＝右邊）；(10.82)→(10.83) 的 Fourier 積分（QAWF＝erfc 閉式）
    book45_argmax(r, u)           照書上 (10.45) 的指數極大在 |r|/(2u) = 0.5；(10.15) 的 f 在 |r|/u = 1（r = −1）
    green_sc(R)                   sc 格子的 G(z_c; R)（Bessel 積分；R = 0 ＝ walk_rec.G_C['sc'] = 1.5164、R = e1 ＝ G_c − 1）：
                                  VIS_R → 0.5164、0.2573、0.1217、0.0599、0.0299；× 2πR/3 = 1.0815、1.0779、1.0198、1.0041、1.0009 → 1；
                                  green_sc_series(R, N=400)：精確部分和＋連續尾巴（第二種算法，差 < 2e−3）
§10.3.1 Redner–Reynolds：漫步當作巨正則系綜 (10.49)–(10.75)，Table 10.1                             [W1005]
    torus_moments(L, z) → (Γ, ⟨N⟩, Var N)（精確多項式）；torus_rho(L, z) = ⟨N⟩/V；torus_rho_fd(L, z)：−∂lnΓ/∂r 的差分（第二種算法）
    RHO_Z = (0.2, z_c, 0.6, 1.0)  torus_rho：3×3：0.1391、0.3374、0.5284、0.6823；4×4：0.0900、0.3056、0.5585、0.7338；
                                  5×5：0.0603、0.2769、0.5764、0.7565（z → ∞：(V−1)/V）
    torus_peak(L) → (z*, max Var N/V)：(0.4443, 0.4289)、(0.4380, 0.5846)、(0.4272, 0.7242)（∂²Ω/∂r² 的峰隨 V 長高、z* 往 z_c 移）
    torus_curve(L, n=60, zmax=1.2)   畫圖用 (z, ρ) 各 n 點
    mf_rho(r, u, N0)、mf_chi(r, u, N0)   平均場 ρ = ⟨N⟩/V = ⟨λ⟩/2、χ' = ∂²Ω/∂r² = (N0/4)Var λ（閉式）；mf_rho_fd、mf_chi_fd（ln Γ 的差分）、mf_rho_quad
                                  MF_RHO_N0 = 10⁶、u = 1：r = −0.5、−0.1、−0.01 → ρ = 0.2500、0.0500、0.0050（= |r|/2u）；r = 0.1、0.5 → 1.0e−05、2.0e−06（= 1/(Vr) → 0）
                                  斜率 d ln ρ / d ln|r| = 1.0000 → β' = 1；χ' = 0.5000（r = −0.1）、1.0e−04、4.0e−06（r = 0.1、0.5 → 0）→ γ' = α' = 0
    book55_rho(r, u)              照書上 (10.55) 的負號：ρ = r/(2u) < 0（r = −0.5 → −0.2500）；lnG53(r, N0, sign)：(10.53) 的 r < 0 式（sign=+1 才和精確值相符）
    torus_var_fd(L, z)            Var N/V ＝ ∂²lnΓ/∂r²/V 的差分（(10.59)，第二種算法）
    rr_exponents(nu, d) → dict(beta = νd − 1, gamma = 2 − νd, alpha = 2 − νd, nu = ν)（(10.73)、(10.67)–(10.72)）
    table101(alpha, nu) → dict(beta = 1 − α, gamma = α, alpha = α, nu = ν)（Table 10.1；第二種算法，α 取 walk_crit.EXPONENTS）
    RR_CASES = {'mf': (1/2, 4), 'd2': (3/4, 2), 'd3': (0.588, 3)}（三維 ν = 0.588 只引用）：
                                  β' = 1、1/2、0.764；γ' = α' = 0、1/2、0.236；ν' = 1/2、3/4、0.588
    rr_laws(e, d) → (α'+2β'+γ', 2−α', dν', 4ν')：全部 = 2、2 − α' = dν'；書上 (10.75) 的 d_cν'：d = 2 得 3 ≠ 3/2、d = 3 得 2.352 ≠ 1.764
    dc_mean_field() = (2 − α'_MF)/ν'_MF = 4（(10.66)）
summary()                         場景卡與講義用的全部數字（dict，字串已照卡上的小數位格式化）

check() 驗證全部（每個結果至少兩種獨立算法、筆誤的更正、plan_ch10.md 場景卡與講義的每個數字）；直接執行本檔會跑 check()（約 4 秒）。
`python walk_cond.py --summary` 印出 summary() 的全部數字。
最慢的呼叫：torus_counts(5) 約 1 秒（快取）、D_box(64) 約 0.1 秒（快取）、green_sc／D_inf 每點約 5 毫秒，其餘都很快。
scipy 的 IntegrationWarning 不影響結果。
"""
from __future__ import annotations

import math
from fractions import Fraction
from functools import lru_cache

import numpy as np
from scipy import integrate, linalg, optimize, special

try:                                    # walk017：utils 套件；combkit：combkit 套件
    from . import walk_saw as SW
    from . import walk_crit as WC
    from . import walk_rec as WR
    from . import walk_gf as GF
except ImportError:                     # 直接執行本檔
    import walk_saw as SW
    import walk_crit as WC
    import walk_rec as WR
    import walk_gf as GF

F = Fraction


def _ro(a):
    a = np.asarray(a)
    a.flags.writeable = False
    return a


def _f(x, nd=4):
    return f"{float(x):.{nd}f}"


def _e(x, nd=1):
    return f"{float(x):.{nd}e}"


# ---------------------------------------------------------------- §10.1  [W1001]
TORI = (3, 4, 5)
STEP_DIRS = ((1, 0), (0, 1), (-1, 0), (0, -1))      # E N W S
STEP_CHARS = "ENWS"
SAW_INF = tuple(SW.SAW_KNOWN_C)
GRAY_Q4 = 91392                                      # OEIS A091299：四維超立方體的有向 Hamilton 路徑（Gray 碼）


@lru_cache(maxsize=None)
def torus_neighbors(L):
    return tuple(tuple(((x + dx) % L) + L * ((y + dy) % L) for dx, dy in STEP_DIRS)
                 for y in range(L) for x in range(L))


@lru_cache(maxsize=None)
def torus_counts(L):
    """8 重對稱：第一步向東、第一次轉彎向北的 ×8；一路直走的 ×4。"""
    V = L * L
    nb = torus_neighbors(L)
    c = [0] * V
    c[0] = 1

    def dfs(s, mask, n):
        c[n] += 8
        for t in nb[s]:
            b = 1 << t
            if not mask & b:
                dfs(t, mask | b, n + 1)

    s, mask = 0, 1
    for k in range(1, L):
        s = nb[s][0]
        mask |= 1 << s
        c[k] += 4
        t = nb[s][1]
        dfs(t, mask | (1 << t), k + 1)
    return tuple(c)


def torus_counts_plain(L, Nmax=None, sym=1):
    """sym=1：完全不用對稱；sym=4：只固定第一步向東（×4）。Nmax 限制深度。"""
    V = L * L
    Nmax = V - 1 if Nmax is None else Nmax
    nb = torus_neighbors(L)
    c = [0] * (Nmax + 1)

    def dfs(s, mask, n):
        c[n] += sym
        if n == Nmax:
            return
        for t in nb[s]:
            b = 1 << t
            if not mask & b:
                dfs(t, mask | b, n + 1)

    if sym == 1:
        dfs(0, 1, 0)
    else:
        s1 = nb[0][0]
        dfs(s1, 1 | (1 << s1), 1)
        c[0] = 1
    return tuple(c)


def torus_counts_brute(L, N):
    x, y = SW.all_walks_xy(N)
    res = (x.astype(np.int64) % L) + L * (y.astype(np.int64) % L)
    res = np.sort(res, axis=1)
    return int(np.all(np.diff(res, axis=1) != 0, axis=1).sum())


def torus_vs_inf(L):
    c = torus_counts(L)
    return [(N, c[N], SAW_INF[N]) for N in range(min(len(c), len(SAW_INF)))]


def torus_first_diff(L):
    rows = torus_vs_inf(L)
    N = next(n for n, a, b in rows if a != b)
    return N, SAW_INF[N] - torus_counts(L)[N]


@lru_cache(maxsize=None)
def torus_walks(L, N):
    assert L <= 4 or N <= 10
    nb = torus_neighbors(L)
    out = []
    path = [0]

    def dfs(s, mask, n):
        if n == N:
            out.append(tuple(path))
            return
        for t in nb[s]:
            b = 1 << t
            if not mask & b:
                path.append(t)
                dfs(t, mask | b, n + 1)
                path.pop()

    dfs(0, 1, 0)
    return tuple(out)


def _step_of(a, b, L):
    dx = (b % L - a % L) % L
    dy = (b // L - a // L) % L
    return {(1, 0): 0, (0, 1): 1, (L - 1, 0): 2, (0, L - 1): 3}[(dx, dy)]


def torus_steps(path, L):
    return "".join(STEP_CHARS[_step_of(a, b, L)] for a, b in zip(path, path[1:]))


def torus_lift(path, L, start=(0, 0)):
    x, y = start
    pts = [(x, y)]
    for a, b in zip(path, path[1:]):
        dx, dy = STEP_DIRS[_step_of(a, b, L)]
        x, y = x + dx, y + dy
        pts.append((x, y))
    return pts


def torus_segments(path, L, start=(0, 0)):
    pts = torus_lift(path, L, start)
    segs = []
    for (x0, y0), (x1, y1) in zip(pts, pts[1:]):
        a = (x0 % L, y0 % L)
        b = (x1 % L, y1 % L)
        dx, dy = x1 - x0, y1 - y0
        if abs(b[0] - a[0]) + abs(b[1] - a[1]) == 1:
            segs.append(((float(a[0]), float(a[1])), (float(b[0]), float(b[1])), False))
        else:
            segs.append(((float(a[0]), float(a[1])), (a[0] + dx / 2, a[1] + dy / 2), True))
            segs.append(((b[0] - dx / 2, b[1] - dy / 2), (float(b[0]), float(b[1])), True))
    return segs


def torus_wraps(path, L, start=(0, 0)):
    return sum(1 for s in torus_segments(path, L, start) if s[2]) // 2


TORUS_SEED = 1004
TORUS_START = (1, 1)
TORUS_EXAMPLE_N = {"short": 3, "mid": 8, "full": 15}
TORUS_EXAMPLE_WRAPS = {"short": 0, "mid": 2, "full": None}


@lru_cache(maxsize=None)
def torus_example(name, L=4, start=TORUS_START, seed=TORUS_SEED):
    N = TORUS_EXAMPLE_N[name]
    want = TORUS_EXAMPLE_WRAPS[name]
    cand = [w for w in torus_walks(L, N) if want is None or torus_wraps(w, L, start) == want]
    rng = np.random.default_rng(seed)
    w = cand[int(rng.integers(len(cand)))]
    return dict(N=N, sites=w, steps=torus_steps(w, L), lift=tuple(torus_lift(w, L, start)),
                wraps=torus_wraps(w, L, start), segments=tuple(torus_segments(w, L, start)), start=start, L=L,
                n_candidates=len(cand))


def torus_moments(L, z):
    c = np.array(torus_counts(L), dtype=float)
    n = np.arange(len(c), dtype=float)
    z = np.asarray(z, dtype=float)
    w = c * np.power.outer(z, n) if z.ndim else c * z ** n
    G = w.sum(axis=-1)
    m1 = (w * n).sum(axis=-1) / G
    m2 = (w * n * n).sum(axis=-1) / G
    return G, m1, m2 - m1 * m1


def torus_mean_N(L, z):
    return float(torus_moments(L, z)[1])


def book101_mean_N(L, z, n_comp=0):
    """照書上 (10.1) 的 z^n：G = Σ_N c_N z^n 與 N 無關，⟨N⟩ = z ∂ln G/∂z 只剩 n = 0。"""
    c = torus_counts(L)
    G = lambda zz: sum(c) * zz ** n_comp
    h = 1e-6 * z
    return z * (math.log(G(z + h)) - math.log(G(z - h))) / (2 * h)


# ---------------------------------------------------------------- §10.2  [W1002]
MF_U = 1.0
MF_N0 = 100
MF_R = 0.5
MF_N0S = (100, 1000, 10000)
MF_NS = (3, 2, 1, 0.1, 0.01)
MF_NS_NEG = (0.01, 0.001, 0.0001, 1e-6)


def mf_f(lam, r, u=MF_U, N0=MF_N0):
    return math.exp(-N0 * (r * lam / 2 + u * lam * lam / 4))


def _mf_x(r, u, N0):
    return r * math.sqrt(N0 / (4 * u))


def _mf_scale(r, u, N0):
    """被積函數 f(λ) 的寬度：r > 0 時在 0 附近以 min(1/B, 1/√A) 衰減；r ≤ 0 時是 λ* = |r|/u 附近的高斯寬 1/√A。"""
    A, B = N0 * u / 4, N0 * r / 2
    return min(1 / B, 1 / math.sqrt(A)) if r > 0 else 1 / math.sqrt(A)


def mf_I(r, u=MF_U, N0=MF_N0, method="erfcx"):
    if method == "erfcx":
        return math.sqrt(math.pi / (N0 * u)) * float(special.erfcx(_mf_x(r, u, N0)))
    ls = max(0.0, -r / u)
    top = -N0 * (r * ls / 2 + u * ls * ls / 4)
    g = lambda lam: math.exp(-N0 * (r * lam / 2 + u * lam * lam / 4) - top)
    s = _mf_scale(r, u, N0)
    b = ls + 80 * s
    val = integrate.quad(g, 0, b, points=[ls] if ls > 0 else None, limit=400, epsabs=0, epsrel=1e-13)[0]
    return val * math.exp(top)


def mf_lnI(r, u=MF_U, N0=MF_N0):
    x = _mf_x(r, u, N0)
    base = 0.5 * math.log(math.pi / (N0 * u))
    if x >= 0:
        return base + math.log(float(special.erfcx(x)))
    return base + math.log(2) + x * x + math.log1p(-0.5 * float(special.erfcx(-x)) * math.exp(-x * x))


def gamma_mf(r, u=MF_U, N0=MF_N0, book=False):
    return (N0 / 4 if book else N0 / 2) * mf_I(r, u, N0)


def ln_gamma_mf(r, u=MF_U, N0=MF_N0):
    return math.log(N0 / 2) + mf_lnI(r, u, N0)


def gamma_n(n, r=MF_R, u=MF_U, N0=MF_N0, method="quad"):
    """有限 n：N0⟨M1²⟩ = (N0/n)·∫λ^{n/2} f / ∫λ^{n/2−1} f（徑向座標 d^nM ∝ λ^{n/2−1}dλ）。"""
    if method == "pcf":
        A, B = N0 * u / 4, N0 * r / 2
        xx = B / math.sqrt(2 * A)

        def mom(s):                     # ∫_0^∞ λ^{s−1} e^{−Aλ²−Bλ} dλ（拋物柱函數）
            return (2 * A) ** (-s / 2) * math.gamma(s) * math.exp(B * B / (8 * A)) * special.pbdv(-s, xx)[0]
        return N0 / n * mom(n / 2 + 1) / mom(n / 2)
    f = lambda lam: mf_f(lam, r, u, N0)
    ls = max(0.0, -r / u)
    s = _mf_scale(r, u, N0)
    b = ls + 80 * s
    kw = dict(limit=400, epsabs=0, epsrel=1e-12)
    pts = [ls] if ls > 0 else None
    num = integrate.quad(lambda lam: lam ** (n / 2) * f(lam), 0, b, points=pts, **kw)[0]
    c = min(1.0, b)
    den = (integrate.quad(lambda lam: lam ** (n / 2 - 1) * (f(lam) - 1.0), 0, c, **kw)[0] + c ** (n / 2) * 2 / n
           + (integrate.quad(lambda lam: lam ** (n / 2 - 1) * f(lam), c, b, points=[p for p in (pts or []) if c < p < b] or None, **kw)[0]
              if b > c else 0.0))
    return N0 / n * num / den


def gamma_n1_direct(r=MF_R, u=MF_U, N0=MF_N0):
    g = lambda M: math.exp(-N0 * (r * M * M / 2 + u * M ** 4 / 4))
    lim = 2 * math.sqrt(max(0.0, -r / u)) + 20 / (N0 * u) ** 0.25 + 10 / math.sqrt(N0 * max(r, 0.0) + 1e-300) * (r > 0)
    lim = min(lim, 50.0)
    num = integrate.quad(lambda M: M * M * g(M), -lim, lim, limit=400)[0]
    den = integrate.quad(g, -lim, lim, limit=400)[0]
    return N0 * num / den


def gamma_gauss_check(r=MF_R, N0=MF_N0, ns=(3, 1, 0.01)):
    """u → 0（高斯）：任何 n 都是 1/r；更正的 (10.15) 給 1/r、書上 1/(2r)。"""
    tiny = 1e-12
    return [gamma_n(n, r, tiny, N0) for n in ns], gamma_mf(r, tiny, N0), gamma_mf(r, tiny, N0, book=True)


def s_integral(lam, eps=0.1):
    a = abs(lam)
    i_sin = 2 * integrate.quad(lambda S: S / (S * S + eps * eps), 0, np.inf, weight="sin", wvar=a)[0]
    i_cos = 2 * integrate.quad(lambda S: eps / (S * S + eps * eps), 0, np.inf, weight="cos", wvar=a)[0]
    return 1j * (math.copysign(i_sin, lam) + i_cos)


def mf_asym_I(r, u=MF_U, N0=MF_N0):
    return 2 / (N0 * r) if r > 0 else 2 * math.sqrt(math.pi / (N0 * u)) * math.exp(N0 * r * r / (4 * u))


def gamma_mf_asym(r, u=MF_U, N0=MF_N0):
    return 1 / r if r > 0 else math.sqrt(math.pi * N0 / u) * math.exp(N0 * r * r / (4 * u))


def ln_gamma_mf_asym(r, u=MF_U, N0=MF_N0):
    return -math.log(r) if r > 0 else 0.5 * math.log(math.pi * N0 / u) + N0 * r * r / (4 * u)


def fn3_check(r=MF_R, u=MF_U, N0=MF_N0):
    return (gamma_mf(r, u, N0), N0 * SW.G_MF(1 - r, N0, u), N0 * SW.G_MF_quad(1 - r, N0, u),
            gamma_mf(r, u, N0, book=True), N0 * SW.G_MF(1 - 2 * r, N0, 4 * u))


# ---------------------------------------------------------------- §10.3  [W1003]
WALL_R, WALL_U = -1.0, 1.0
WALL_X = (1, 2, 4)
SW_R, SW_U, SW_Q2 = F(-1), F(1), F(1)
SW_M = (F(2, 3), F(2, 3), F(1, 3))
SEED_SW = 1003


def wall_M(x, r=WALL_R, u=WALL_U):
    return math.sqrt(-r / u) * math.tanh(x * math.sqrt(-r / 2))


def wall_bvp(r=WALL_R, u=WALL_U, X=20.0, n=2001):
    x = np.linspace(0.0, X, n)
    h = x[1] - x[0]
    M0 = math.sqrt(-r / u)
    M = M0 * np.minimum(1.0, x)
    for _ in range(30):
        Mi = M[1:-1]
        res = -(M[2:] - 2 * Mi + M[:-2]) / h ** 2 + r * Mi + u * Mi ** 3
        ab = np.zeros((3, n - 2))
        ab[0, 1:] = -1 / h ** 2
        ab[1, :] = 2 / h ** 2 + r + 3 * u * Mi ** 2
        ab[2, :-1] = -1 / h ** 2
        d = linalg.solve_banded((1, 1), ab, -res)
        M[1:-1] += d
        if np.max(np.abs(d)) < 1e-14:
            break
    return x, M


def wall_errors(r=WALL_R, u=WALL_U):
    x, M = wall_bvp(r, u)
    err_bvp = float(np.max(np.abs(M - np.array([wall_M(t, r, u) for t in x]))))
    hh = 1e-3
    xs = np.linspace(0.1, 8, 80)
    Mf = lambda t: wall_M(t, r, u)
    resid = max(abs(-(Mf(t + hh) - 2 * Mf(t) + Mf(t - hh)) / hh ** 2 + r * Mf(t) + u * Mf(t) ** 3) for t in xs)
    return err_bvp, float(resid)


def ex101_grad_check(seed=SEED_SW, n=2, sites=12, r=-0.8, u=1.3, hstep=0.5):
    """一維週期鏈、n 分量：ℋ = Σ hstep[(r/2)|M|² + (u/4)|M|⁴ − h·M + ½|ΔM/hstep|²]；數值梯度 vs (10.17) 的離散版。"""
    rng = np.random.default_rng(seed)
    M = rng.normal(size=(sites, n))
    hf = rng.normal(size=n)

    def H(Mv):
        m2 = (Mv * Mv).sum(1)
        grad = ((np.roll(Mv, -1, 0) - Mv) ** 2).sum(1) / hstep ** 2
        return hstep * float(np.sum(r / 2 * m2 + u / 4 * m2 * m2 - Mv @ hf + 0.5 * grad))

    lap = (np.roll(M, -1, 0) - 2 * M + np.roll(M, 1, 0)) / hstep ** 2
    eos = -lap + r * M + u * M * (M * M).sum(1, keepdims=True) - hf
    num = np.zeros_like(M)
    e = 1e-6
    for i in range(sites):
        for a in range(n):
            Mp, Mm = M.copy(), M.copy()
            Mp[i, a] += e
            Mm[i, a] -= e
            num[i, a] = (H(Mp) - H(Mm)) / (2 * e)
    return float(np.max(np.abs(num - hstep * eos)) / np.max(np.abs(hstep * eos)))


def ex102_coeffs(seed=SEED_SW, n=3, r=-0.7, u=1.3):
    """g(ε) = (r/2)|M+εσ|² + (u/4)|M+εσ|⁴ − h·(M+εσ) 的 ε^k 係數：數值（5 點內插）vs 公式；回傳最大誤差。"""
    rng = np.random.default_rng(seed)
    M, s, hf = rng.normal(size=n), rng.normal(size=n), rng.normal(size=n)
    g = lambda e: r / 2 * np.dot(M + e * s, M + e * s) + u / 4 * np.dot(M + e * s, M + e * s) ** 2 - np.dot(hf, M + e * s)
    es = np.array([-2.0, -1.0, 0.0, 1.0, 2.0])
    coef = np.linalg.solve(np.vander(es, 5, increasing=True), np.array([g(e) for e in es]))
    M2, s2, Ms = M @ M, s @ s, M @ s
    formula = [r / 2 * M2 + u / 4 * M2 * M2 - hf @ M, (r + u * M2) * Ms - hf @ s,
               r / 2 * s2 + u / 2 * M2 * s2 + u * Ms * Ms, u * Ms * s2, u / 4 * s2 * s2]
    return float(np.max(np.abs(coef - np.array(formula))))


def sw_lambdas(r=SW_R, q2=SW_Q2, u=SW_U, M2=None):
    M2 = sum(m * m for m in SW_M) if M2 is None else M2
    l1 = r / 2 + q2 / 2 + u * M2 / 2
    return l1, l1 + u * M2


def sw_X(M=SW_M, r=SW_R, q2=SW_Q2, u=SW_U):
    M2 = sum(m * m for m in M)
    l1 = r / 2 + q2 / 2 + u * M2 / 2
    n = len(M)
    return [[(l1 if i == j else 0) + u * M[i] * M[j] for j in range(n)] for i in range(n)]


def sw_Xinv(M=SW_M, r=SW_R, q2=SW_Q2, u=SW_U):
    M2 = sum(m * m for m in M)
    l1, l2 = sw_lambdas(r, q2, u, M2)
    n = len(M)
    return [[(1 / l1 if i == j else 0) - u * M[i] * M[j] / (l1 * l2) for j in range(n)] for i in range(n)]


def _matmul(A, B):
    return [[sum(A[i][k] * B[k][j] for k in range(len(B))) for j in range(len(B[0]))] for i in range(len(A))]


def sw_check(ns=(1, 2, 3, 4, 5, 6), seed=SEED_SW, trials=4):
    rng = np.random.default_rng(seed)
    err = dict(eig=0.0, inv=0.0, det=0.0, tr=0.0)
    mult_ok = True
    for n in ns:
        for _ in range(trials):
            M = rng.normal(size=n)
            r, q2, u = rng.uniform(-2, 2), rng.uniform(0, 3), rng.uniform(0.2, 2)
            M2 = float(M @ M)
            l1, l2 = sw_lambdas(r, q2, u, M2)
            if abs(l1) < 0.05:
                continue
            X = np.array(sw_X(tuple(M), r, q2, u))
            Xi = np.array(sw_Xinv(tuple(M), r, q2, u))
            ev = np.sort(np.linalg.eigvalsh(X))
            want = np.sort(np.array([l1] * (n - 1) + [l2]))
            err["eig"] = max(err["eig"], float(np.max(np.abs(ev - want))))
            mult_ok &= int(np.sum(np.abs(ev - l1) < 1e-9)) == n - 1 and int(np.sum(np.abs(ev - l2) < 1e-9)) == 1
            I = np.eye(n)
            err["inv"] = max(err["inv"], float(np.max(np.abs(Xi @ X - I))), float(np.max(np.abs(X @ Xi - I))),
                             float(np.max(np.abs(Xi - np.linalg.inv(X)))))
            err["det"] = max(err["det"], abs(np.linalg.det(X) / (l1 ** (n - 1) * l2) - 1))
            err["tr"] = max(err["tr"], abs(np.trace(X) - ((n - 1) * l1 + l2)))
    return err, mult_ok


def mode_table(ns=(3, 2, 1, 0)):
    l1, l2 = sw_lambdas()
    return [(n, n - 1, 1, l2 * l1 ** (n - 1), float(l2 * l1 ** (n - 1)) ** -0.5) for n in ns]


GAUSS_M = (0.6, 0.8)
GAUSS_H = (0.5, 0.0)


def gauss_n2_check(npts=1201, lim=12.0):
    M = GAUSS_M
    X = np.array(sw_X(M, -1.0, 1.0, 1.0))
    Xi = np.array(sw_Xinv(M, -1.0, 1.0, 1.0))
    h = np.array(GAUSS_H)
    t = np.linspace(-lim, lim, npts)
    dx = t[1] - t[0]
    S1, S2 = np.meshgrid(t, t, indexing="ij")
    E = -(X[0, 0] * S1 * S1 + 2 * X[0, 1] * S1 * S2 + X[1, 1] * S2 * S2) + h[0] * S1 + h[1] * S2
    num = float(np.exp(E).sum() * dx * dx)
    pref = math.pi / math.sqrt(np.linalg.det(X))
    q = float(h @ Xi @ h)
    return num, pref * math.exp(q / 4), pref * math.exp(-q / 4)


# ---------------------------------------------------------------- §10.3–10.4  [W1004]
HOLE_L = 64
HOLE_R = -1.0
HOLE_D = (0, 1, 2, 4, 8, 16, 32)
SC_EWALD = 2.837297                     # 簡單立方週期 Coulomb 和（去掉零模）的常數：G_box(0) − G_∞(0) ≈ −2.837297/(4πL)
CONT_RHO = (1, 2, 4, 8)


def _q2_lattice(L):
    q = 2 * np.pi * np.fft.fftfreq(L)
    s = 2 - 2 * np.cos(q)
    return s[:, None, None] + s[None, :, None] + s[None, None, :]


def _hole_weight(q2, r, u, M2):
    with np.errstate(divide="ignore", invalid="ignore"):
        w = 1 / (r + q2 + u * M2) - 1 / (r + q2 + 3 * u * M2)
    return w


@lru_cache(maxsize=None)
def D_box(L=HOLE_L, r=HOLE_R, u=1.0):
    M2 = -r / u
    w = _hole_weight(_q2_lattice(L), r, u, M2)
    w[0, 0, 0] = 0.0
    return _ro(np.fft.ifftn(w).real)


def D_axis(L=HOLE_L, r=HOLE_R):
    return D_box(L, r)[: L // 2 + 1, 0, 0].copy()


def D_direct(L, r, delta, u=1.0):
    q = 2 * np.pi * np.arange(L) / L
    s = 2 - 2 * np.cos(q)
    q2 = s[:, None, None] + s[None, :, None] + s[None, None, :]
    w = _hole_weight(q2, r, u, -r / u)
    w[0, 0, 0] = 0.0
    ph = (delta[0] * q)[:, None, None] + (delta[1] * q)[None, :, None] + (delta[2] * q)[None, None, :]
    return float((w * np.cos(ph)).sum() / L ** 3)


def D_book46(L=HOLE_L, r=HOLE_R):
    """照書上 (10.46)：Σψψ'/(q²(2|r|+q²))（少了分子的 2|r|）＝ D/(2|r|)。"""
    q2 = _q2_lattice(L)
    with np.errstate(divide="ignore"):
        w = 1 / (q2 * (2 * abs(r) + q2))
    w[0, 0, 0] = 0.0
    return np.fft.ifftn(w).real[: L // 2 + 1, 0, 0]


def _bessel_tail3(a, S, scale):
    """∏ ive(a_i, scale·t) 在 t > S 的積分（大 t 漸近：(2π·scale·t)^{−3/2}(1 − Σ(4a²−1)/(8·scale·t))）。"""
    c = sum((4 * x * x - 1) / 8 for x in a) / scale
    k = (2 * math.pi * scale) ** -1.5
    return k * (2 * S ** -0.5 - c * (2 / 3) * S ** -1.5)


def D_inf(delta, r=HOLE_R, u=1.0):
    a = [abs(int(v)) for v in delta]
    k2 = 2 * abs(r)
    f = lambda t: (1 - math.exp(-k2 * t)) * special.ive(a[0], 2 * t) * special.ive(a[1], 2 * t) * special.ive(a[2], 2 * t)
    S = 4000.0
    v = integrate.quad(f, 0, 40, limit=400, epsabs=1e-14)[0] + integrate.quad(f, 40, S, limit=400, epsabs=1e-14)[0]
    return v + _bessel_tail3(a, S, 2.0)


def box_shift(L, r, delta):
    d2 = sum(int(v) ** 2 for v in delta)
    return -SC_EWALD / (4 * math.pi * L) + d2 / (6 * L ** 3) + 1 / (L ** 3 * 2 * abs(r))


def coulomb_fit(L=HOLE_L, r=HOLE_R, dmin=4, dmax=16):
    ax = D_axis(L, r)
    d = np.arange(dmin, dmax + 1)
    A = np.vstack([np.ones(len(d)), -1.0 / d]).T
    c0, c = np.linalg.lstsq(A, -ax[d], rcond=None)[0]
    return float(c0), float(c)


def D_cont(rho, r=HOLE_R):
    k = math.sqrt(2 * abs(r))
    return (1 - math.exp(-k * rho)) / (4 * math.pi * rho)


def D_cont_quad(rho, r=HOLE_R):
    k2 = 2 * abs(r)
    v = integrate.quad(lambda q: k2 / (q * (q * q + k2)) if q > 0 else 0.0, 0, np.inf, weight="sin", wvar=rho)[0]
    return v / (2 * math.pi ** 2 * rho)


def _yukawa_inf(delta, r):
    a = [abs(int(v)) for v in delta]
    f = lambda t: math.exp(-r * t) * special.ive(a[0], 2 * t) * special.ive(a[1], 2 * t) * special.ive(a[2], 2 * t)
    return integrate.quad(f, 0, 60, limit=400, epsabs=1e-15)[0] + integrate.quad(f, 60, np.inf, limit=400)[0]


def case1_check(L=32, r=1.0, dmax=4):
    G = np.fft.ifftn(1 / (r + _q2_lattice(L))).real
    return max(abs(G[d, 0, 0] / _yukawa_inf((d, 0, 0), r) - 1) for d in range(dmax + 1))


VIS_N = 100
VIS_R = (1, 2, 4, 8, 16)


def _vec(R):
    return (int(R), 0, 0) if np.isscalar(R) else tuple(int(v) for v in R)


def _p1(m, k):
    m = np.asarray(m)
    k = abs(int(k))
    ok = (m >= k) & ((m + k) % 2 == 0)
    mm = np.where(ok, m, k)
    val = np.exp(special.gammaln(mm + 1) - special.gammaln((mm + k) // 2 + 1) - special.gammaln((mm - k) // 2 + 1) - mm * math.log(2))
    return np.where(ok, val, 0.0)


def sc_prob(n, R):
    x, y, z = _vec(R)
    a = np.arange(n + 1)
    b = n - a
    logw = special.gammaln(n + 1) - special.gammaln(a + 1) - special.gammaln(b + 1) + a * math.log(1 / 3) + b * math.log(2 / 3)
    return float(np.sum(np.exp(logw) * _p1(a, x) * _p1(b, y + z) * _p1(b, y - z)))


def sc_visits(N, R):
    return sum(sc_prob(n, R) for n in range(N + 1))


def sc_visits_dp(N, R):
    assert N <= 40
    x, y, z = _vec(R)
    W = 2 * N + 1
    P = np.zeros((W, W, W))
    P[N, N, N] = 1.0
    tot = P[N + x, N + y, N + z] if max(abs(x), abs(y), abs(z)) <= N else 0.0
    probs = [tot]
    for _ in range(N):
        P = (np.roll(P, 1, 0) + np.roll(P, -1, 0) + np.roll(P, 1, 1) + np.roll(P, -1, 1)
             + np.roll(P, 1, 2) + np.roll(P, -1, 2)) / 6
        probs.append(P[N + x, N + y, N + z])
    return float(sum(probs)), probs


def erfc_form(R, N):
    R = float(np.linalg.norm(_vec(R)))
    return 3 / (2 * math.pi * R) * math.erfc(R * math.sqrt(3 / (2 * (N + 1))))


def erfc_form_book(R, N):
    R = float(np.linalg.norm(_vec(R)))
    return 3 / (2 * math.pi * R) * math.erfc(3 * R * R / (N + 1))


def eq1080_check(k2=0.7, zc=6.0, N=12):
    """(10.80)：z^N 的係數（級數摺積）＝ z_c^{−N} ×（書上的左邊）＝ z_c^{−N}·[1 − (1+k²/z_c)^{−(N+1)}]/(z_c k²)。"""
    coef = sum((k2 + zc) ** -(n + 1) * zc ** -(N - n + 1) for n in range(N + 1))
    lhs = 1 / zc / (zc + k2) * sum((zc / (zc + k2)) ** n for n in range(N + 1))
    rhs = (1 - (1 + k2 / zc) ** -(N + 1)) / (zc * k2)
    return coef, zc ** -N * lhs, zc ** -N * rhs


def eq1083_check(R=3.0, a=2.0):
    """(10.82)→(10.83) 的 Fourier 積分：∫d³k/(2π)³ e^{ik·R}(1 − e^{−ak²})/k² ＝ erfc(R/(2√a))/(4πR)（QAWF 數值 vs 閉式）。"""
    v = integrate.quad(lambda k: (1 - math.exp(-a * k * k)) / k if k > 0 else 0.0, 0, np.inf, weight="sin", wvar=R)[0]
    return v / (2 * math.pi ** 2 * R), math.erfc(R / (2 * math.sqrt(a))) / (4 * math.pi * R)


def book45_argmax(r=HOLE_R, u=1.0):
    """(10.45) 的指數：照書上 −(rλ + uλ²) 的極大在 |r|/(2u)；(10.15) 的 f 在 |r|/u（數值極值）。"""
    b = optimize.minimize_scalar(lambda lam: r * lam + u * lam * lam, bounds=(0, 10), method="bounded", options=dict(xatol=1e-12)).x
    f = optimize.minimize_scalar(lambda lam: r * lam / 2 + u * lam * lam / 4, bounds=(0, 10), method="bounded", options=dict(xatol=1e-12)).x
    return float(b), float(f)


def green_sc(R):
    a = [abs(v) for v in _vec(R)]
    f = lambda s: 3 * special.ive(a[0], s) * special.ive(a[1], s) * special.ive(a[2], s)
    S = 2000.0
    v = integrate.quad(f, 0, 50, limit=400, epsabs=1e-14)[0] + integrate.quad(f, 50, S, limit=400, epsabs=1e-14)[0]
    return v + 3 * _bessel_tail3(a, S, 1.0)


def green_sc_series(R, N=400):
    """第二種算法：精確的 Σ_{n≤N} P_n(R) ＋ 連續極限的尾巴 (3/(2πR)) erf(R√(3/(2(N+1))))。"""
    Rn = float(np.linalg.norm(_vec(R)))
    return sc_visits(N, R) + 3 / (2 * math.pi * Rn) * math.erf(Rn * math.sqrt(3 / (2 * (N + 1))))


# ---------------------------------------------------------------- §10.3.1  [W1005]
RHO_Z = (0.2, WC.ZC_2D, 0.6, 1.0)
MEAN_N_Z = (0.2, WC.ZC_2D, 1.0)
MF_RHO_N0 = 10 ** 6
MF_RHO_R = (-0.5, -0.1, -0.01, 0.1, 0.5)
RR_CASES = {"mf": (F(1, 2), 4), "d2": (F(3, 4), 2), "d3": (0.588, 3)}


def torus_rho(L, z):
    return torus_moments(L, z)[1] / (L * L)


def torus_rho_fd(L, z, h=1e-5):
    """⟨N⟩ = −∂lnΓ/∂r，e^{−r} = z/z_c（中央差分）。"""
    zc = WC.ZC_2D
    r = math.log(zc / z)
    lnG = lambda rr: math.log(float(torus_moments(L, zc * math.exp(-rr))[0]))
    return -(lnG(r + h) - lnG(r - h)) / (2 * h) / (L * L)


def torus_var_fd(L, z, h=1e-4):
    zc = WC.ZC_2D
    r = math.log(zc / z)
    lnG = lambda rr: math.log(float(torus_moments(L, zc * math.exp(-rr))[0]))
    return (lnG(r + h) - 2 * lnG(r) + lnG(r - h)) / h ** 2 / (L * L)


@lru_cache(maxsize=None)
def torus_peak(L):
    res = optimize.minimize_scalar(lambda z: -torus_moments(L, z)[2] / (L * L), bounds=(0.25, 1.0), method="bounded",
                                   options=dict(xatol=1e-10))
    return float(res.x), float(-res.fun)


def torus_curve(L, n=60, zmax=1.2):
    zs = np.linspace(0.0, zmax, n)
    return zs, np.array(torus_rho(L, zs), dtype=float)


def _mf_lam_mean(r, u, N0):
    x = _mf_x(r, u, N0)
    ex = float(special.erfcx(x))
    return 2 / math.sqrt(N0 * u) * (1 / (math.sqrt(math.pi) * ex) - x)


def mf_rho(r, u=MF_U, N0=MF_RHO_N0):
    return _mf_lam_mean(r, u, N0) / 2


def mf_chi(r, u=MF_U, N0=MF_RHO_N0):
    A, B = N0 * u / 4, N0 * r / 2
    lam = _mf_lam_mean(r, u, N0)
    var = (1 - B * lam) / (2 * A) - lam * lam
    return N0 / 4 * var


def mf_rho_fd(r, u=MF_U, N0=MF_RHO_N0, h=None):
    h = 1e-4 * max(abs(r), 1e-3) if h is None else h
    return -(ln_gamma_mf(r + h, u, N0) - ln_gamma_mf(r - h, u, N0)) / (2 * h) / N0


def mf_chi_fd(r, u=MF_U, N0=MF_RHO_N0, h=None):
    h = 1e-3 * max(abs(r), 1e-3) if h is None else h
    return (ln_gamma_mf(r + h, u, N0) - 2 * ln_gamma_mf(r, u, N0) + ln_gamma_mf(r - h, u, N0)) / h ** 2 / N0


def mf_rho_quad(r, u=MF_U, N0=10 ** 4):
    f = lambda lam: math.exp(-N0 * (r * lam / 2 + u * lam * lam / 4) - N0 * (r * r / (4 * u) if r < 0 else 0.0))
    ls = max(0.0, -r / u)
    s = _mf_scale(r, u, N0)
    b = ls + 80 * s
    kw = dict(points=[ls] if ls > 0 else None, limit=400, epsabs=0, epsrel=1e-12)
    m0 = integrate.quad(f, 0, b, **kw)[0]
    m1 = integrate.quad(lambda lam: lam * f(lam), 0, b, **kw)[0]
    return m1 / m0 / 2


def book55_rho(r, u=MF_U):
    """書上 (10.55)：Ω = (1/V){ln V − V r²/(4u) + …} → ρ = −∂Ω/∂r = r/(2u)（r < 0 時是負的）。"""
    return r / (2 * u)


def lnG53(r, N0, u=MF_U, sign=+1):
    """(10.53) r < 0（照書上的 N0/4 正規化，Γ 多乘一個 N0）：ln[(1/2)(π/u)^{1/2} N0^{3/2} e^{±N0r²/4u}]。"""
    return math.log(0.5 * math.sqrt(math.pi / u)) + 1.5 * math.log(N0) + sign * N0 * r * r / (4 * u)


def rr_exponents(nu, d):
    nd = nu * d
    return dict(beta=nd - 1, gamma=2 - nd, alpha=2 - nd, nu=nu)


def table101(alpha, nu):
    return dict(beta=1 - alpha, gamma=alpha, alpha=alpha, nu=nu)


def rr_laws(e, d):
    return (e["alpha"] + 2 * e["beta"] + e["gamma"], 2 - e["alpha"], d * e["nu"], 4 * e["nu"])


def dc_mean_field():
    e = rr_exponents(*RR_CASES["mf"])
    return (2 - e["alpha"]) / e["nu"]


# ---------------------------------------------------------------- summary / check
def _rr_str(v):
    return str(v) if isinstance(v, Fraction) else f"{v:.3f}"


def summary():
    S = {}
    # W1001
    for L in TORI:
        S[f"torus{L}"] = list(torus_counts(L))
        S[f"torus{L}_total"] = sum(torus_counts(L))
    S["saw_inf"] = list(SAW_INF[:16])
    S["first_diff"] = {L: torus_first_diff(L) for L in TORI}
    S["gray_q4"] = (16 * torus_counts(4)[15], GRAY_Q4)
    for name in TORUS_EXAMPLE_N:
        ex = torus_example(name)
        S[f"example_{name}"] = (ex["N"], ex["steps"], ex["wraps"], ex["lift"][-1], ex["n_candidates"])
    S["mean_N_4"] = [_f(torus_mean_N(4, z)) for z in MEAN_N_Z]
    S["book101_mean_N"] = _f(book101_mean_N(4, 0.3))
    # W1002
    S["gamma_n"] = [_f(gamma_n(n)) for n in MF_NS]
    S["gamma_mf"] = (_f(gamma_mf(MF_R)), _f(gamma_mf(MF_R, book=True)))
    S["gamma_n_neg"] = [_f(gamma_n(n, -MF_R) / gamma_mf(-MF_R)) for n in MF_NS_NEG]
    S["gamma_mf_neg"] = f"{gamma_mf(-MF_R):.1f}"
    gg = gamma_gauss_check()
    S["gauss_u0"] = ([_f(v) for v in gg[0]], _f(gg[1]), _f(gg[2]))
    S["s_integral"] = (f"{s_integral(1.0).imag:.4f}", f"{abs(s_integral(-1.0)):.4f}", f"{2 * math.pi * math.exp(-0.1):.4f}")
    S["gamma_pos"] = [_f(gamma_mf(MF_R, MF_U, N0)) for N0 in MF_N0S]
    S["ratio_neg"] = [_f(math.exp(ln_gamma_mf(-MF_R, MF_U, N0) - ln_gamma_mf_asym(-MF_R, MF_U, N0))) for N0 in MF_N0S]
    S["lnG_neg"] = [_f(ln_gamma_mf(-MF_R, MF_U, N0) / N0) for N0 in MF_N0S]
    S["fn3"] = [_f(v) for v in fn3_check()]
    S["fn3_neg"] = [_f(v) for v in fn3_check(-MF_R)]
    # W1003
    S["wall"] = [_f(wall_M(x)) for x in WALL_X] + [_f(math.sqrt(2 / abs(WALL_R)))]
    S["wall_err"] = [_e(v) for v in wall_errors()]
    S["sw_lambdas"] = [str(v) for v in sw_lambdas()]
    S["sw_X18"] = [[str(18 * v) for v in row] for row in sw_X()]
    S["sw_Xinv27"] = [[str(27 * v) for v in row] for row in sw_Xinv()]
    S["mode_table"] = [(n, a, b, str(d), _f(s)) for n, a, b, d, s in mode_table()]
    S["book33"] = _f(math.sqrt(float(sw_lambdas()[1] / sw_lambdas()[0])))
    S["gauss_n2"] = [_f(v) for v in gauss_n2_check()]
    # W1004
    S["hole_D"] = [_f(D_axis()[d]) for d in HOLE_D]
    S["hole_1mD"] = [_f(1 - D_axis()[d]) for d in HOLE_D]
    S["book46_D0"] = _f(D_book46()[0])
    c0, c = coulomb_fit()
    S["coulomb"] = (_f(c0), _f(c), _f(1 / (4 * math.pi)), _f(4 * math.pi * c), _f(SC_EWALD / (4 * math.pi * HOLE_L)))
    S["D_cont"] = [_f(D_cont(p)) for p in CONT_RHO]
    S["D_inf"] = [_f(D_inf((p, 0, 0))) for p in CONT_RHO]
    S["case1"] = _e(case1_check())
    S["visits"] = [f"{sc_visits(VIS_N, R):.4f}" if R < 16 else f"{sc_visits(VIS_N, R):.6f}" for R in VIS_R]
    S["erfc_form"] = [f"{erfc_form(R, VIS_N):.4f}" if R < 16 else f"{erfc_form(R, VIS_N):.6f}" for R in VIS_R]
    S["visit_ratio"] = [f"{sc_visits(VIS_N, R) / erfc_form(R, VIS_N):.3f}" for R in VIS_R]
    S["erfc_book"] = [_e(erfc_form_book(R, VIS_N)) for R in (8, 16)]
    S["green"] = [_f(green_sc(R)) for R in VIS_R]
    S["green_ratio"] = [_f(green_sc(R) * 2 * math.pi * R / 3) for R in VIS_R]
    # W1005
    S["rho"] = {L: [_f(torus_rho(L, z)) for z in RHO_Z] for L in TORI}
    S["rho_limit"] = {L: str(F(L * L - 1, L * L)) for L in TORI}
    S["peak"] = {L: (_f(torus_peak(L)[0]), _f(torus_peak(L)[1])) for L in TORI}
    S["mf_rho"] = [_f(mf_rho(r)) if r < 0 else _e(mf_rho(r)) for r in MF_RHO_R]
    S["mf_slope"] = _f(math.log(mf_rho(-0.1) / mf_rho(-0.01)) / math.log(10))
    S["mf_chi"] = [_f(mf_chi(-0.1)), _e(mf_chi(0.1)), _e(mf_chi(0.5))]
    S["book55"] = _f(book55_rho(-0.5))
    S["rr"] = {k: {q: _rr_str(v) for q, v in rr_exponents(*a).items()} for k, a in RR_CASES.items()}
    S["rr_laws"] = {k: [_rr_str(v) for v in rr_laws(rr_exponents(*a), a[1])] for k, a in RR_CASES.items()}
    S["dc"] = str(dc_mean_field())
    return S


def check():
    # ---- W1001：環面列舉（8 重對稱 DFS ＝ 不用對稱的 DFS ＝ 4^N 暴力 ＝ 文獻）
    c3, c4, c5 = (torus_counts(L) for L in TORI)
    assert c3 == (1, 4, 12, 32, 80, 160, 256, 296, 168) and sum(c3) == 1009
    assert c4 == (1, 4, 12, 36, 96, 264, 624, 1536, 3072, 6240, 9840, 15912, 17472, 19200, 10656, 5712) and sum(c4) == 90677
    assert c5[-1] == 236456 and sum(c5) == 22094601 and len(c5) == 25
    assert torus_counts_plain(3) == c3 and torus_counts_plain(4) == c4 and torus_counts_plain(4, sym=4) == c4
    assert torus_counts_plain(5, Nmax=16, sym=4) == c5[:17]
    for L in TORI:
        for N in range(9):
            assert torus_counts_brute(L, N) == torus_counts(L)[N], (L, N)
        assert torus_first_diff(L) == (L, 4)
        assert all(a == b for n, a, b in torus_vs_inf(L) if n < L)
        assert all(a < b for n, a, b in torus_vs_inf(L) if n >= L)
    assert 16 * c4[15] == GRAY_Q4
    assert SW.saw_enum(15)[0][:16] == SAW_INF[:16] and SAW_INF[15] == 6416596      # DFS ＝ 文獻（OEIS A001411）
    assert (c3[3], c4[4], c5[5]) == (32, 96, 280)
    for name, N in TORUS_EXAMPLE_N.items():
        ex = torus_example(name)
        assert len(ex["sites"]) == N + 1 == len(set(ex["sites"])) and len(ex["steps"]) == N
        L = ex["L"]
        sx, sy = ex["start"]
        wrapped = [((x % L), (y % L)) for x, y in ex["lift"]]
        assert len(set(wrapped)) == N + 1
        assert wrapped == [(((s % L) + sx) % L, ((s // L) + sy) % L) for s in ex["sites"]]
        if TORUS_EXAMPLE_WRAPS[name] is not None:
            assert ex["wraps"] == TORUS_EXAMPLE_WRAPS[name]
        assert all(abs(a[0] - b[0]) + abs(a[1] - b[1]) in (0.5, 1.0) for a, b, _ in ex["segments"])
    assert len(torus_walks(4, 15)) == c4[15] and len(torus_walks(4, 8)) == c4[8]
    for z in MEAN_N_Z:
        assert abs(torus_mean_N(4, z) - 16 * torus_rho_fd(4, z)) < 1e-6
    S = summary()
    assert S["mean_N_4"] == ["1.4400", "4.8899", "11.7409"], S["mean_N_4"]
    assert [S[f"example_{k}"][:4] for k in TORUS_EXAMPLE_N] == [(3, "WNE", 0, (1, 2)), (8, "WSSEENNE", 2, (3, 1)),
                                                               (15, "WSSEEENNNEEESSW", 3, (5, 0))]
    assert S["book101_mean_N"] == "0.0000"

    # ---- W1002：平均場（erfcx ＝ quad；徑向 quad ＝ 拋物柱 ＝ 一維直接；n → 0 ＝ N0/2·I；漸近）
    for r in (MF_R, -MF_R, 1.0, -1.0):
        for N0 in (10, 100, 1000):
            assert abs(mf_I(r, MF_U, N0, "quad") / mf_I(r, MF_U, N0) - 1) < 1e-10
            assert abs(math.log(mf_I(r, MF_U, N0)) - mf_lnI(r, MF_U, N0)) < 1e-10
    for n in MF_NS + (0.001,):
        assert abs(gamma_n(n, method="pcf") / gamma_n(n) - 1) < 1e-8, n
    assert abs(gamma_n(1) / gamma_n1_direct() - 1) < 1e-9
    assert abs(gamma_n(1, -MF_R) / gamma_n1_direct(-MF_R) - 1) < 1e-8
    assert abs(gamma_n(1e-4) / gamma_mf(MF_R) - 1) < 1e-4 and abs(gamma_n(1e-7, -MF_R) / gamma_mf(-MF_R) - 1) < 1e-4
    for n in MF_NS_NEG:                   # r < 0：要 n·e^{N0r²/4u} ≪ 1 才到 n → 0 的極限（先取 n → 0）
        assert abs(gamma_n(n, -MF_R, method="pcf") / gamma_n(n, -MF_R) - 1) < 1e-8
    dev = [1 - gamma_n(n, -MF_R) / gamma_mf(-MF_R) for n in MF_NS_NEG]
    assert 0.0095 < dev[3] / dev[2] < 0.0108                     # 偏差 ∝ n
    assert S["gamma_n_neg"] == ["0.3299", "0.8311", "0.9801", "0.9998"] and S["gamma_mf_neg"] == "9179.7", S["gamma_n_neg"]
    assert S["gamma_n"] == ["1.7183", "1.7634", "1.8131", "1.8624", "1.8676"], S["gamma_n"]
    assert S["gamma_mf"] == ("1.8682", "0.9341"), S["gamma_mf"]
    assert S["gauss_u0"] == (["2.0000", "2.0000", "2.0000"], "2.0000", "1.0000"), S["gauss_u0"]
    assert S["s_integral"] == ("5.6853", "0.0000", "5.6853"), S["s_integral"]
    assert abs(s_integral(1.0, 1e-3).imag - 2 * math.pi) < 1e-2          # ε → 0：2πi（不是 iπ）
    assert S["gamma_pos"] == ["1.8682", "1.9844", "1.9984"], S["gamma_pos"]
    assert S["ratio_neg"] == ["0.9998", "1.0000", "1.0000"], S["ratio_neg"]
    assert S["lnG_neg"] == ["0.0912", "0.0665", "0.0630"], S["lnG_neg"]
    assert abs(gamma_mf(MF_R, MF_U, 10 ** 6) * MF_R - 1) < 1e-5
    assert S["fn3"] == ["1.8682"] * 3 + ["0.9341"] * 2, S["fn3"]
    a, b, c, d, e = fn3_check(-MF_R)
    assert abs(b / a - 1) < 1e-12 and abs(c / a - 1) < 1e-9 and abs(e / d - 1) < 1e-12 and abs(a / d - 2) < 1e-12
    assert abs(SW.ex75(MF_N0, 0.5, 1.0) - (1 - 0.5) * MF_N0 * SW.G_MF(0.5, MF_N0, 1.0)) < 1e-12

    # ---- W1003：自旋波（tanh ＝ BVP；梯度；多項式係數；特徵值／逆矩陣 Fraction ＝ numpy）
    assert S["wall"] == ["0.6089", "0.8884", "0.9930", "1.4142"], S["wall"]
    eb, er = wall_errors()
    assert eb < 1e-5 and er < 1e-5 and S["wall_err"] == ["2.8e-06", "8.5e-08"], (eb, er)
    assert ex101_grad_check() < 1e-7
    assert ex102_coeffs() < 1e-11
    l1, l2 = sw_lambdas()
    assert (l1, l2) == (F(1, 2), F(3, 2))
    X, Xi = sw_X(), sw_Xinv()
    assert [[18 * v for v in row] for row in X] == [[17, 8, 4], [8, 17, 4], [4, 4, 11]]
    assert [[27 * v for v in row] for row in Xi] == [[38, -16, -8], [-16, 38, -8], [-8, -8, 50]]
    I3 = [[F(int(i == j)) for j in range(3)] for i in range(3)]
    assert _matmul(X, Xi) == I3 and _matmul(Xi, X) == I3
    ev = np.sort(np.linalg.eigvalsh(np.array(X, dtype=float)))
    assert np.allclose(ev, [0.5, 0.5, 1.5], atol=1e-13)
    err, mult = sw_check()
    assert mult and max(err.values()) < 1e-11, err
    mt = mode_table()
    assert [(n, a, b, str(d)) for n, a, b, d, _ in mt] == [(3, 2, 1, "3/8"), (2, 1, 1, "3/4"), (1, 0, 1, "3/2"), (0, -1, 1, "3")]
    for n in (1, 2, 3):                    # 整數 n：numpy 的 det
        Mv = np.array([0.6, 0.8, 0.0][:n]) if n > 1 else np.array([1.0])
        Xn = np.array(sw_X(tuple(Mv), -1.0, 1.0, 1.0))
        assert abs(np.linalg.det(Xn) - float(l2 * l1 ** (n - 1))) < 1e-12
    assert [s for *_, s in S["mode_table"]] == ["1.6330", "1.1547", "0.8165", "0.5774"]
    assert abs(mt[-1][4] - math.sqrt(float(l1 / l2))) < 1e-15 and S["book33"] == "1.7321"
    g = gauss_n2_check()
    assert abs(g[0] / g[1] - 1) < 1e-10 and S["gauss_n2"] == ["3.9891", "3.9891", "3.2988"], S["gauss_n2"]

    # ---- W1004：相關洞（FFT ＝ 直接求和 ＝ 無限格子＋Ewald；連續 quad ＝ 閉式；Coulomb 擬合）
    ax = D_axis()
    for dlt in ((0, 0, 0), (3, 0, 0), (2, 1, 1), (8, 5, 0)):
        assert abs(D_direct(HOLE_L, HOLE_R, dlt) - D_box()[dlt]) < 1e-12
    for dlt in ((0, 0, 0), (1, 0, 0), (2, 0, 0), (4, 0, 0), (8, 0, 0), (3, 2, 1)):
        assert abs(D_box()[dlt] - D_inf(dlt) - box_shift(HOLE_L, HOLE_R, dlt)) < 5e-6, dlt
    assert np.all(np.diff(ax) < 0)                                     # D 隨距離遞減（書上說 increases）
    assert abs(D_box().sum()) < 1e-10                                    # 零模去掉：平均為 0
    assert S["hole_D"] == ["0.1083", "0.0613", "0.0359", "0.0167", "0.0065", "0.0016", "-0.0001"], S["hole_D"]
    assert S["hole_1mD"] == ["0.8917", "0.9387", "0.9641", "0.9833", "0.9935", "0.9984", "1.0001"], S["hole_1mD"]
    assert np.allclose(D_book46() * 2 * abs(HOLE_R), ax, atol=1e-14) and S["book46_D0"] == "0.0541"
    assert S["coulomb"] == ("0.0035", "0.0801", "0.0796", "1.0064", "0.0035"), S["coulomb"]
    for p in CONT_RHO:
        assert abs(D_cont_quad(p) - D_cont(p)) < 1e-9
    assert S["D_cont"] == ["0.0602", "0.0374", "0.0198", "0.0099"], S["D_cont"]
    assert S["D_inf"] == ["0.0648", "0.0394", "0.0202", "0.0100"], S["D_inf"]
    assert abs(D_cont(1e4) * 4 * math.pi * 1e4 - 1) < 1e-12
    assert case1_check() < 1e-10 and S["case1"] == "1.2e-11"
    for n in (0, 2, 4, 10, 20):
        assert abs(sc_prob(n, 0) - GF.returns_sc(n) / 6 ** n) < 1e-14
    up = WR.return_probs("sc", 30)
    assert max(abs(sc_prob(n, 0) - up[n]) for n in range(31)) < 1e-13
    for R in VIS_R + ((3, 2, 1),):
        tot, probs = sc_visits_dp(30, R)
        assert abs(tot - sc_visits(30, R)) < 1e-13
    assert abs(green_sc(0) - WR.G_C["sc"]) < 1e-9 and abs(green_sc(1) - (WR.G_C["sc"] - 1)) < 1e-9
    assert _f(WR.G_C["sc"]) == "1.5164" and _f(green_sc(0)) == "1.5164"
    for R in VIS_R:
        assert abs(green_sc_series(R) / green_sc(R) - 1) < 2e-3, R
    assert S["visits"] == ["0.4509", "0.1931", "0.0610", "0.0103", "0.000183"], S["visits"]
    assert S["erfc_form"] == ["0.4121", "0.1744", "0.0586", "0.0100", "0.000174"], S["erfc_form"]
    assert S["visit_ratio"] == ["1.094", "1.108", "1.042", "1.028", "1.052"], S["visit_ratio"]
    assert S["erfc_book"] == ["4.3e-04", "1.7e-28"], S["erfc_book"]
    c1, c2, c3 = eq1080_check()
    assert abs(c2 / c1 - 1) < 1e-12 and abs(c3 / c1 - 1) < 1e-12
    for R, aa in ((3.0, 2.0), (1.0, 50.0), (8.0, 101 / 6)):
        q, cl = eq1083_check(R, aa)
        assert abs(q - cl) < 1e-10, (R, aa, q, cl)
    bm, fm = book45_argmax()
    assert abs(bm - 0.5) < 1e-6 and abs(fm - 1.0) < 1e-6
    assert S["green"] == ["0.5164", "0.2573", "0.1217", "0.0599", "0.0299"], S["green"]
    assert S["green_ratio"] == ["1.0815", "1.0779", "1.0198", "1.0041", "1.0009"], S["green_ratio"]

    # ---- W1005：Redner–Reynolds（環面：比值 ＝ 差分；平均場：閉式 ＝ 差分 ＝ quad；指數兩條路）
    for L in TORI:
        for z in RHO_Z + (0.45, 1.5):
            assert abs(torus_rho(L, z) - torus_rho_fd(L, z)) < 1e-8
            assert abs(torus_moments(L, z)[2] / L ** 2 - torus_var_fd(L, z)) < 1e-5
        assert abs(torus_rho(L, 1e6) - (L * L - 1) / (L * L)) < 1e-5
    assert S["rho"] == {3: ["0.1391", "0.3374", "0.5284", "0.6823"], 4: ["0.0900", "0.3056", "0.5585", "0.7338"],
                        5: ["0.0603", "0.2769", "0.5764", "0.7565"]}, S["rho"]
    assert S["peak"] == {3: ("0.4443", "0.4289"), 4: ("0.4380", "0.5846"), 5: ("0.4272", "0.7242")}, S["peak"]
    for r in MF_RHO_R:
        assert abs(mf_rho_fd(r) / mf_rho(r) - 1) < 1e-6, r
    for r in (-0.5, -0.1, 0.1, 0.5):
        assert abs(mf_rho_quad(r) / mf_rho(r, MF_U, 10 ** 4) - 1) < 1e-9
        assert abs(mf_chi_fd(r, MF_U, 10 ** 4) / mf_chi(r, MF_U, 10 ** 4) - 1) < 1e-4, r
    assert S["mf_rho"] == ["0.2500", "0.0500", "0.0050", "1.0e-05", "2.0e-06"], S["mf_rho"]
    assert S["mf_slope"] == "1.0000" and S["mf_chi"] == ["0.5000", "1.0e-04", "4.0e-06"], (S["mf_slope"], S["mf_chi"])
    assert S["book55"] == "-0.2500"
    N0 = 1000
    exact53 = math.log(N0) + math.log(N0 / 4) + mf_lnI(-0.5, MF_U, N0)       # 書上的正規化：N0 × (N0/4)I
    assert abs(exact53 - lnG53(-0.5, N0, sign=+1)) < 1e-3 and lnG53(-0.5, N0, sign=-1) < exact53 - 100
    for k, (nu, d) in RR_CASES.items():
        e = rr_exponents(nu, d)
        alpha = WC.EXPONENTS["mf"]["alpha"] if k == "mf" else (WC.EXPONENTS["saw2d"]["alpha"] if k == "d2" else 2 - d * nu)
        if k != "d3":
            assert alpha == 2 - d * nu
        t = table101(alpha, nu)
        assert all(abs(e[q] - t[q]) < 1e-12 for q in e)
        lw = rr_laws(e, d)
        assert abs(lw[0] - 2) < 1e-12 and abs(lw[1] - lw[2]) < 1e-12
    assert rr_exponents(*RR_CASES["mf"]) == dict(beta=1, gamma=0, alpha=0, nu=F(1, 2))
    assert S["rr"] == {"mf": dict(beta="1", gamma="0", alpha="0", nu="1/2"), "d2": dict(beta="1/2", gamma="1/2", alpha="1/2", nu="3/4"),
                       "d3": dict(beta="0.764", gamma="0.236", alpha="0.236", nu="0.588")}, S["rr"]
    assert S["rr_laws"] == {"mf": ["2", "2", "2", "2"], "d2": ["2", "3/2", "3/2", "3"], "d3": ["2.000", "1.764", "1.764", "2.352"]}, S["rr_laws"]
    assert S["dc"] == "4"
    assert WC.EXPONENTS["saw2d"]["nu"] == RR_CASES["d2"][0] and WC.NU_KNOWN[3] == RR_CASES["d3"][0]
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
    print(f"walk_cond（第 10 章）：環面上的自迴避列舉（三種算法、與無限平面並排）、平均場積分與 n → 0、註 3 的對照、自旋波矩陣的特徵值與逆矩陣、"
          f"牆邊的 tanh、相關洞的 Coulomb 形與 Ewald 常數、erfc 形的步數、Redner–Reynolds 的 ρ 與指數全部驗證通過（{time.time() - t0:.1f} 秒）")
