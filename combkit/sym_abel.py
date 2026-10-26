# 由 combkit/tools/sync_from_harter_sym.py 從 harter_sym 專案的 codes/utils/sym_abel.py 複製；請改 harter_sym 那份再同步。
r"""
sym_abel.py —— Harter《Principles of Symmetry, Dynamics, and Spectroscopy》（Wiley 1993）第 2 章
Basic Theory and Applications of Symmetry Representations（Abelian Symmetry Groups）的數學
（純 Python＋numpy＋sympy＋scipy，不 import manim；「最小方程」的投影算符借用同資料夾的 sym_mat）。
書上的式號、筆誤與〔注〕見 ref/md/ch02_formulas.md（相對於 C:\mypython\manim\projects\harter_sym）。

約定
  * 精確值用 sympy（Rational、sqrt、I、單位根 ε = e^{2πi/n} 展成 cos + i sin）；浮點用 numpy。f-string 不能直接 :.4f 一個 sympy 數或
    Fraction，要 float()。卡片上的小數一律是 f"{float(x):.4f}"。ħ = 1。
  * 運動方程一律寫成 ẍ = −𝒜x（加速度矩陣 𝒜，本徵值 = ω²）；書上 (2.8.1a) 漏了負號（見筆誤）。
  * 兩擺 (2.2.1c)：a = k/mℓ²（書上印 2k/mℓ²，是筆誤；本徵值 α₋ = 2a + b = 2k/mℓ² + g/ℓ 才對），b = g/ℓ。
  * 單擺環 (2.7.1)：|r^p⟩ = r^p|1⟩，p = 0..n−1；𝓡(r) 把第 p 個基底送到第 p+1 個（(2.6.7)：𝓡(r)x = (x_n, x_1, …)）。
    |k_n⟩ 的分量 e^{−ik_n p}/√n（(2.7.3)），k_n = 2πk/n；ω²(k_n) = b + 2a(1 − cos k_n)（(2.7.4)）。
  * 群乘法表 table：dict {(g, h): gh}（列 g、行 h，先 h 後 g 作用）。

§2.1–2.2 群（H201–H202）
    ROT3 = {"1", "Rx", "Ry", "Rz"} 的 3×3 轉動矩陣（繞 x、y、z 轉 180°）；D2_NAMES；D2_TABLE_BOOK（(2.1.3) 照抄的 4 列）
    table_from_matrices(mats)、table_from_rows(names, rows)、table_rows(table, names)；(2.1.1) RzRy = Rx
    C6_NAMES/C6_BOOK（(2.2.10)：A 生成，A⁶ = 1）、NOT_GROUP_NAMES/NOT_GROUP_BOOK（(2.2.11)）
    axioms(table) → {closure, associative, identity, inverses, latin}；assoc_failures(table) → [(x, y, z)]；is_abelian(table)
        (2.2.10)：四條公理都成立、Abel；(2.2.11)：封閉、單位元、反元素、每列每行不重複都成立，**只有結合律不成立**：
        (ab)c = dc = a ≠ a(bc) = af = c；216 組 (x, y, z) 中 assoc_failures 有 100 組；它也可交換（對角線全是 1）
    cyclic_table(n)（"1", "r", "r^2", …）、C2_TABLE、C3_TABLE；element_order(table, g)、order_counts(table)
    GROUP_COUNTS（階 1–64 的群個數，OEIS A000001 引用值）：1–16 → 1,1,1,2,1,2,1,5,2,2,1,5,1,2,1,14；24 → 15、32 → 51、48 → 52、64 → 267
        （書上 Fig. 2.2.2 印 294 ＝ Miller 1930 的舊數；Hall–Senior 1964 更正為 267）；GROUP_COUNT_64_BOOK = 294
    abelian_count(n) = Π_p p(e_p)（分割數）：8 → 3、16 → 5、24 → 3、32 → 7、64 → 11、72 → 6（＝ ABELIAN_COUNTS 引用值）；abelian_types(n)
        1–16 的 Abel 群個數 1,1,1,2,1,1,1,3,2,1,1,2,1,1,1,5；1–64 中有 28 個階數只有 Abel 群（1, 2, 3, 4, 5, 7, 9, 11, …, 61）
§2.2–2.3 兩個耦合擺（H202–H204）：PEND2 = {a: 21/200, b: 1}（ω₊ = 1、ω₋ = 11/10）
    accel2(a, b) = [[a+b, −a], [−a, a+b]]；R2 = [[0, 1], [1, 0]]；commutator(A, B)（[𝒜, R] = 0 ⇔ 對角相等、非對角相等）
    c2_projectors() = {1: 𝒫⁺, −1: 𝒫⁻}（最小方程 R² − 1 = 0，用 sym_mat.projectors）＝ c2_projectors_char()：
        𝒫⁺ = ½[[1, 1], [1, 1]]、𝒫⁻ = ½[[1, −1], [−1, 1]]；本徵向量 (1, ±1)/√2；α₊ = b = 1、α₋ = 2a + b = 121/100
    normal2(a, b) → (α₊, α₋)；beat_times2(a, b) = (π/Δω, 2π/Δω) = (10π, 20π) = (31.4159, 62.8319)（Δω = 1/10）
    Modal(𝒜, x0, v0)：ẍ = −𝒜x 的解析解（可呼叫 x(t)、.velocity(t)，給 updater 用，很快）＝ ode_solve(𝒜, x0, v0, ts)（solve_ivp）
    DETUNED = (a, b₁, b₂) = (21/200, 1, 221/200)（沒有 C₂ 對稱）：tan 2θ = 2a/(b₂ − b₁) = 2，x₂ 最大振幅 sin 2θ = 2/√5 = 0.8944（不到 1）
§2.3(a)(b) 組態空間與相量（H204）：LISSA = {weak: 21/20, mid: 13/10, strong: 2}（ω₋/ω₊，ω₊ = 1、b = 1）
    a_for_ratio(ρ) = (ρ² − 1)/2：weak 41/800、mid 69/200 = 0.3450、strong 3/2
    lissajous(t, ρ) = ((cos t + cos ρt)/2, (cos t − cos ρt)/2)（從 (1, 0) 靜止放手）＝ ode 數值；closure(ρ) →
        weak：n₊ = 20、n₋ = 21、t_LJ = 40π = 125.6637、拍 1 次；mid：10、13、20π = 62.8319、3 拍（t_beat = 20π/3）；strong：1、2、2π = 6.2832、1 拍
    ellipse_axes(a, b, V0=½)：位能等高線的半軸 1/√α₊、1/√α₋（沿 (1,1)、(1,−1) 方向，轉 45°）：strong 1、0.5000；mid 1、0.7692；weak 1、0.9524
    quad_potential(a, b) = 2×2 的力常數矩陣 𝓕（= accel2，m = ℓ = 1）
    PHASOR_NU = (½, 1) Hz（(+)、(−) 模，Fig. 2.3.7）；phasors(t) → (z₁, z₂)；PHASOR_TIMES = (0, ¼, ½, ¾, 1) s：
        |z₁| = 1、0.9239、0.7071、0.3827、0；|z₂| = 0、0.3827、0.7071、0.9239、1；z₁ = −i cot(πt/2) z₂（x₁ 的相量順時針領先 90°）
§2.4–2.5 不可約表示、部分解題（H205）
    C2_CHARS（Table 2.4.1）；O_R3 = (2.4.5) 的 𝒪(R)；rep_projectors(reps, chars)、multiplicities(reps, chars)：f⁺ = 2、f⁻ = 1
        𝒪(𝒫⁺) = [[½,0,½],[0,1,0],[½,0,½]]、𝒪(𝒫⁻) = [[½,0,−½],[0,0,0],[−½,0,½]]；T3（(2.4.7)）：T3†𝒪(R)T3 = diag(1, 1, −1)
    PEND3 = {m: 1, M: 2, k: 1, l: 1, g: 1} → a = 1、A = ½、b = 1；accel3() = [[2, −1, 0], [−½, 2, −½], [0, −1, 2]]（(2.5.1c)）
        reduce3() = T3†𝒜T3 = [[a+b, −√2a, 0], [−√2A, 2A+b, 0], [0, 0, a+b]] = [[2, −√2, 0], [−√2/2, 2, 0], [0, 0, 2]]（(2.5.3)，(2,2) 更正為 2A+b）
        omega2_3() = 1、3（偶，2×2 還要解）、2（奇）；ω = 1.0000、1.7321、1.4142；向量 (1,1,1)、(1,−1,1)、(1,0,−1)
§2.6–2.7 C₃、C_n 環（H206–H207）
    eps(n, k) = e^{2πik/n}（精確）；cyclic_rep(n) = 𝓡(r)；cn_char_table(n)（(2.6.5)、Fig. 2.7.2）
    cn_projectors_char(n)（(2.6.6b)）＝ cn_projectors_minimal(n)（(2.3.2) 的連乘）＝ sym_mat.projectors(𝓡(r))（n = 3）：
        𝒫^{0₃} = ⅓J；𝒫^{1₃} = ⅓[[1, ε, ε*], [ε*, 1, ε], [ε, ε*, 1]]（ε = −½ + (√3/2)i；＝書上 (2.6.8) 的 𝒫²，ε₂ = ε、ε₃ = ε*）
        cn_projectors_np(n, method)（浮點版，大 n 用）
    RING3 = {a: 1, b: 1}；ring_matrix(n, a, b)（循環三對角）；ring_omega2(n, a, b)；T†𝓡(r)T = diag(1, ε, ε*)（(2.6.9)）
        C₃：ω² = 1、4、4（b、3a+b）→ ω = 1、2、2；駐波 (2,−1,−1)/√6、(0,−1,1)/√2；行波 |1₃⟩ = (1, ε*, ε)/√3
        wave_vector(n, k)、cos_wave(n, k)、sin_wave(n, k)（(2.7.3)、(2.7.5)）
    n = 6（a = 1、b = 0）：ω² = 0、1、3、4、3、1 → ω = 0、1、1.7321、2、1.7321、1
    dispersion(k, a, b) = √(2a(1 − cos k) + b) = 2√a sin(k/2)（b = 0；(2.7.6)）；小 k：√a k − √a k³/24（相速度 √a = 1）
    brillouin(n)：−n/2 < k ≤ n/2（n 偶）、|k| ≤ (n−1)/2（n 奇）——各 n 個（書上兩端都寫 <，是筆誤）；alias(n, k)
    ALIAS = (5, (2, 7, −3))：cos(2πkp/5)，p = 0..4 → 1.0000、−0.8090、0.3090、0.3090、−0.8090（三個 k 一樣）
§2.7A–B 對稱破缺與奔馳波（H208）：BROKEN = {n: 24, abar: 3/2, alow: ½, b: 0}
    broken_ring_matrix(n, ā, a̲, b)（彈簧 p→p+1：p 偶 a̲、p 奇 ā；(2.7.8)）、broken_block(k, ā, a̲, b)（(2.7.9) 加回 b）
    broken_bands(k, ā, a̲, b) = b + ā + a̲ ∓ √(ā² + a̲² + 2āa̲ cos k)（(2.7.10) 加回 b）＝ 24×24 數值對角化：
        聲學帶 ω² ∈ [0, 1]、光學帶 [3, 4]；k = π：2a̲ = 1、2ā = 3；C₂₄ 的 ω² = 2（二重）裂成 1、3
    broken_standing(n) = (soft, hard)：k = π 的駐波 (1,−1,−1,1,…)/√n 只拉軟彈簧 → ω² = 2a̲ = 1、(1,1,−1,−1,…)/√n 只拉硬彈簧 → 2ā = 3
    GALLOP = {a: (6/5, 3/5), b: (1, 9/10)}；swr(I, R) = (I − R)/(I + R) = 1/3、1/19（包絡極大 I+R、極小 I−R；Δ = 極小/極大）
        gallop_psi(x, t, I, R)、gallop_envelope(x, I, R)（(2.7.14)）、gallop_zeros(t, Δ)（Re ψ 的零點）
        zero_velocity(t, Δ) = cΔ/(Δ² sin²ωt + cos²ωt)：在 cΔ 與 c/Δ 之間 —— (a) 0.3333c–3c、(b) 0.0526c–19c（＝零點的數值微分）
        SWR_LIST = (1, 0.5, 0.1, 0, −0.1, −0.5, −1)（Fig. 2.7.9）
§2.8–2.11 其他 Abel 對稱（H209）
    D2_REP4（(2.8.2b)，由群表＝正則表示）＝ D2_REP4_BOOK（(2.8.2a) 照抄）；D2_CHARS（(2.8.7)：A1、B2、B1、A2 對 1、Rx、Ry、Rz）
    d2_projectors_product()（(2.8.5b)：𝒫^{x±}𝒫^{y±}）＝ d2_projectors_char()；𝒫¹ = ¼J（全部 ¼）
    D2_ABC = (A, a, b, c) = (11/4, −½, −1, −¼)（單擺 g/ℓ = 1 加彈簧 k_a = ½、k_b = 1、k_c = ¼）；d2_matrix()
        ω² = A+a+b+c = 1、A−a+b−c = 5/2、A+a−b−c = 7/2、A−a−b+c = 4 → ω = 1.0000、1.5811、1.8708、2.0000；向量 (1,±1,±1,±1)/2
    abelian_group(orders)、abelian_char_table(orders)、outer_char_table(T1, T2)（(2.10) 的外積＝Kronecker 積）、is_cyclic(orders)
        ABELIAN8：C8 階數分布 {1:1, 2:1, 4:2, 8:4}、C4×C2 {1:1, 2:3, 4:4}、C2×C2×C2 {1:1, 2:7}；C3×C2 ≅ C6、C2×C2 ≇ C4
        abelian_types(n)：8 → (8)、(4,2)、(2,2,2)；72 → 6 種
    POINT3（3×3：1、I、sxy、sxz、syz、Rx、Ry、Rz、S4）；I·σxy = Rz（(2.11.4)）、σxz·σyz = Rz（(2.11.5)）；is_group(names)
        C2H = (1, Rz, I, sxy)、C2V = (1, Rz, sxz, syz)、D2 都是群且 ≅ C2×C2（point_table(names)）；書上的 C2V_BOOK = (1, Rz, sxy, syz) 不封閉（sxy·syz = Ry）
        S4（90° 轉動反演）生成 4 階循環群 ≅ C4；32 個晶體點群中 16 個 Abel（引用）
§2.12A C₁₂「時鐘烷」與 Bloch 波（H210）：E = ħ²/2μr² = 1，n = 12；V(φ) = −V cos 12φ ⇒ G(±12) = −V/2，其餘 0
    rotor_matrix(V, mmax)（m 基底，(2.12.12)）、rotor_block(κ, V, N)（只含 m ≡ κ mod 12，(2.12.16)）、rotor_bands(V, nb)
    WEAK_V = 2（G(12) = −1）：pair_exact(5) = 37 ∓ √145 → ε(±5) = 24.9584、ε(±7) = 49.0416；pair_pert(5) = 25 − 1/24 = 24.9583、49.0417；
        全基底數值 24.9546、49.0384（多出 |17⟩、|−19⟩ 的修正 −1/264、−1/312）；m = ±6：ε(c₆) = 36 + G = 35、ε(s₆) = 37（一階），
        數值 34.9965、36.9965，能隙 2|G(12)| = 2（數值 2.0000）；V = 2 最低 8 個（不計簡併）：−0.0139、0.9857、3.9844、8.9815、15.9750、24.9546、34.9965、36.9965
    STRONG_V = 360（緊束縛）：tight_fit() → 最低帶 H = −208.6223、S = 0.0897（帶寬 4S = 0.3588）；12 個能階與 H − 2S cos k 的最大偏差 0.0004
    tb_matrix(n, H, S)（(2.12.23) 第二式）、tb_levels(n, H, S) = H − 2S cos(2πk/n)：C₃ → H−2S、H+S、H+S（書上 Fig. 2.12.4 印 H+2S）；C₄ → H−2S、H、H、H+2S
    mathieu_edges(Vp, nb)（縮放單位 ψ'' + (ε + V' cos 2t)ψ = 0，V' = V/36、ε' = ε/36）＝ scipy 的 a_r(q)、b_r(q)（q = V'/2）：
        帶 j = [a_j, b_{j+1}]；V' = 0 時能隙起於 1、4、9（Fig. 2.12.5）；rotor_edges(V, nb) = 36 × 上面（轉子單位）
    square_well_G(m, U)：n 方井位能的 Fourier 係數 G(0) = U/2、G((2N+1)n) = (−1)^{N+1}U/((2N+1)π)（書上印 (−1)^N 4π/(2N+1)）＝ FFT
    Kronig–Penney（Fig. 2.12.6）：ψ'' + (E − U(x))ψ = 0，井寬 = 位障寬 = π/2（週期 π）；kp_half_trace(E, U)（閉式）＝ kp_half_trace_ode（solve_ivp）
        kp_edges(U, Emax)、kp_gap_near(U, E)；能隙閉合（偶然簡併）KP_CROSSINGS (U, E) = (12, 16)、(20, 36)、(32, 36)、(28, 64)、(48, 64)、(60, 64)
        （U 差 ±1 能隙就打開 0.13–0.45）；seg_matrix(E, U, L)（(ψ, ψ') 傳遞矩陣）
§2.12B 兩能階與氨（H211）：NH3 = {H: 0, S: 1}
    two_level(pE, H, S) = H ∓ √(S² + (pE)²)（(2.12.34)；(2.12.33) 的矩陣 [[H − pE, −S], [−S, H + pE]] 數值對角化相同）
        pE = 0、1、2、3 → ±1.0000、±1.4142、±2.2361、±3.1623；小場 ≈ S + (pE)²/2S、大場 ≈ pE
    tunnel_populations(t) = (cos² St, sin² St)：週期 π/S = 3.1416、上→下 π/2 = 1.5708；分裂 2S = 2
    NH3_NU = 23.8701 GHz（(J,K) = (3,3) 反轉線，氨邁射）→ 2S = hν = 9.872e−05 eV、週期 1/ν = 41.89 ps（引用值）
    lz_formula(β) = e^{−πS²/β}（Landau–Zener，pE = βt，書外）＝ lz_numeric(β)：β = 4 → 0.4559（快：跳到另一支）；β = ¼ → 3.5e−06（慢：絕熱跟隨）
§2.12C–D 散射、穿越矩陣、WKB（H212）：WELL = {V: 25, a: 1}（ħ = 1、2m = 1；井內位能 0、井外 V；e = E − V > 0）
    wavenumbers(e) = (l, k) = (√(e+V), √e)（(2.12.37) 更正為 /ħ）
    s_matrix(e)（(2.12.38)，D 的負號已更正）＝ s_matrix_direct(e)（4×4 邊界條件數值解）；s_matrix(e, book=True) = −S
    phase_shifts(e) = (δ₊, δ₋)：tan(ka + δ₊) = (l/k) tan la、tan(ka + δ₋) = −(l/k) cot la；S 的本徵值 e^{2iδ±}、本徵向量 (1, ±1)/√2
        e = 5：δ₊ = −0.2919、δ₋ = −1.0675；|S₁₁|² = 0.5098（反射）、|S₂₁|² = 0.4902（穿透）
    transmission(e) = 1/(1 + V² sin²(2la)/(4e(e+V)))；resonances()：2la = nπ ⇒ e = 14.4784、36.6850、63.8264（T = 1，δ₊ − δ₋ ≡ π/2 mod π）
    crossing_well(e)（(2.12.51)）＝ crossing_step 兩個相乘（(2.12.50b)）；det C = 1；s_from_c(C)（S₁₁ = C₂₁/C₁₁，書上印 C₂₁/C₁₂）、c_from_s(S)
    bound_geometric(V, w)（習題 2.12.2：k = √V sin δ 與 k = (nπ − 2δ)/w 的交點）＝ bound_fd(V, w)：V = 25、w = 2 →
        k = 1.3064、2.5957、3.8375、4.9063；E = 1.7068、6.7379、14.7262、24.0717（井底起算；＝偶／奇超越方程）
    c_barrier(θ)、c_valley(α)（(2.12.54)(2.12.55)）；½ tr(C_b C_v) = √(1+θ²) cos α（(2.12.57)）；wkb_cluster(…)（(2.12.60)）
    wkb_vs_exact(Vp)（Mathieu 最低帶：WKB 帶寬 4S 與精確帶寬之比）：V' = 10 → 精確 0.0100、WKB 0.0089、比 0.8958（θ = 601.85）
習題：problem_233()、response_234(ω, …)、problem_261()、problem_262()、problem_271(n)、problem_273()、problem_274()、problem_275()、
    problem_281()、problem_282()、problem_283()、problem_2121()、problem_2122()、problem_2123() —— 數字見各函式說明與 summary()
    （2.12.3 的第二種算法：line_nodes／line_levels（節點計數＋二分）對 ring4_fd／line4_fd（有限差分））
check()：每個上卡片／講義的數字用兩種以上算法核對（約 3 秒）；summary() → {名稱: 字串}。
直接執行：python utils/sym_abel.py（在 codes/ 裡）印一句成功訊息；--summary 印出全部數字。
"""
import math
from collections import Counter
from itertools import product

import numpy as np
import sympy as sp
from scipy import integrate, linalg, optimize, special

try:
    from . import sym_mat as SM
except ImportError:                      # 直接執行 python utils/sym_abel.py
    import sym_mat as SM

I = sp.I
R_ = sp.Rational
HBAR = 1


def _f4(x) -> str:
    return f"{float(x):.4f}"


def _close(a, b, tol: float = 1e-10) -> bool:
    return np.allclose(np.asarray(a, dtype=complex), np.asarray(b, dtype=complex), atol=tol, rtol=0)


def _zero(M) -> bool:
    return sp.simplify(sp.expand_complex(M)) == sp.zeros(*M.shape)


def to_np(M) -> np.ndarray:
    return np.array([[complex(sp.N(x, 30)) for x in row] for row in sp.Matrix(M).tolist()], dtype=complex)


def _real(M) -> np.ndarray:
    if isinstance(M, sp.MatrixBase):
        return np.array(M.evalf(), dtype=float)
    return np.asarray(M, dtype=float)


# =============================================================================
# §2.1–2.2 群、乘法表、公理（H201–H202）
# =============================================================================
ROT3 = {"1": sp.diag(1, 1, 1), "Rx": sp.diag(1, -1, -1), "Ry": sp.diag(-1, 1, -1), "Rz": sp.diag(-1, -1, 1)}
D2_NAMES = ("1", "Rx", "Ry", "Rz")
D2_TABLE_BOOK = (("1", "Rx", "Ry", "Rz"),            # (2.1.3)：列 g、行 h 的格子 = gh
                 ("Rx", "1", "Rz", "Ry"),
                 ("Ry", "Rz", "1", "Rx"),
                 ("Rz", "Ry", "Rx", "1"))
C6_NAMES = ("1", "A", "B", "C", "D", "F")
C6_BOOK = ("1ABCDF", "ABCDF1", "BCDF1A", "CDF1AB", "DF1ABC", "F1ABCD")          # (2.2.10)
NOT_GROUP_NAMES = ("1", "a", "b", "c", "d", "f")
NOT_GROUP_BOOK = ("1abcdf", "a1dbfc", "bd1fca", "cbf1ad", "dfca1b", "fcadb1")  # (2.2.11)


def table_from_rows(names, rows) -> dict:
    """照抄的表 → {(g, h): gh}。rows[i][j] 是第 i 列（g = names[i]）、第 j 行（h = names[j]）的格子。"""
    return {(g, h): rows[i][j] for i, g in enumerate(names) for j, h in enumerate(names)}


def table_from_matrices(mats: dict) -> dict:
    """由表示矩陣算乘法表：gh = 與 mats[g]·mats[h] 相等的那個名字。"""
    out = {}
    for g, Mg in mats.items():
        for h, Mh in mats.items():
            P = Mg * Mh
            out[(g, h)] = next(k for k, Mk in mats.items() if Mk == P)
    return out


def table_rows(table: dict, names) -> list:
    return [[table[(g, h)] for h in names] for g in names]


def _elements(table: dict) -> list:
    seen = []
    for (g, _h) in table:
        if g not in seen:
            seen.append(g)
    return seen


def assoc_failures(table: dict) -> list:
    """所有 (x, y, z) 使 (xy)z ≠ x(yz)。"""
    els = _elements(table)
    return [(x, y, z) for x in els for y in els for z in els
            if table[(table[(x, y)], z)] != table[(x, table[(y, z)])]]


def axioms(table: dict) -> dict:
    """Table 2.2.1 的四條公理（外加「每列每行每個元素恰好一次」）。"""
    els = _elements(table)
    closure = all(v in els for v in table.values())
    ident = next((e for e in els if all(table[(e, g)] == g == table[(g, e)] for g in els)), None)
    inverses = ident is not None and all(any(table[(g, h)] == ident == table[(h, g)] for h in els) for g in els)
    latin = all(sorted(table[(g, h)] for h in els) == sorted(els) and sorted(table[(h, g)] for h in els) == sorted(els)
                for g in els)
    return {"closure": closure, "associative": not assoc_failures(table), "identity": ident,
            "inverses": inverses, "latin": latin}


def is_abelian(table: dict) -> bool:
    return all(table[(g, h)] == table[(h, g)] for (g, h) in table)


def cyclic_names(n: int) -> list:
    return ["1", "r"] + [f"r^{k}" for k in range(2, n)] if n > 1 else ["1"]


def cyclic_table(n: int) -> dict:
    """C_n = {1, r, …, r^{n−1}}：r^i r^j = r^{(i+j) mod n}（(2.2.12)、(2.6.2)）。"""
    nm = cyclic_names(n)
    return {(nm[i], nm[j]): nm[(i + j) % n] for i in range(n) for j in range(n)}


C2_TABLE = cyclic_table(2)
C3_TABLE = cyclic_table(3)


def _identity(table: dict):
    els = _elements(table)
    return next((e for e in els if all(table[(e, g)] == g == table[(g, e)] for g in els)), None)


def element_order(table: dict, g) -> int:
    e = _identity(table)
    k, x = 1, g
    while x != e:
        x = table[(x, g)]
        k += 1
    return k


def order_counts(table: dict) -> dict:
    return dict(sorted(Counter(element_order(table, g) for g in _elements(table)).items()))


# 階 1–64 的群個數（OEIS A000001；Fig. 2.2.2），與 Abel 群個數（OEIS A000688，用來核對 abelian_count）
GROUP_COUNTS = dict(enumerate(
    [1, 1, 1, 2, 1, 2, 1, 5, 2, 2, 1, 5, 1, 2, 1, 14, 1, 5, 1, 5, 2, 2, 1, 15, 2, 2, 5, 4, 1, 4, 1, 51,
     1, 2, 1, 14, 1, 2, 2, 14, 1, 6, 1, 4, 2, 2, 1, 52, 2, 5, 1, 5, 1, 15, 2, 13, 2, 2, 1, 13, 1, 2, 4, 267], start=1))
GROUP_COUNT_64_BOOK = 294          # 書上 Fig. 2.2.2 的數（G. A. Miller 1930）；Hall–Senior（1964）更正為 267
ABELIAN_COUNTS = {1: 1, 2: 1, 4: 2, 8: 3, 9: 2, 12: 2, 16: 5, 24: 3, 27: 3, 32: 7, 36: 4, 40: 3, 48: 5, 54: 3,
                  56: 3, 64: 11, 72: 6}


def abelian_count(n: int) -> int:
    """有限 Abel 群的個數 = Π_p p(e_p)（n = Π p^{e_p}，p(e) 是 e 的分割數）——§2.10 的定理。"""
    out = 1
    for _p, e in sp.factorint(n).items():
        out *= sp.npartitions(e)
    return int(out)


def abelian_types(n: int) -> list:
    """n 階 Abel 群的全部結構（質數冪循環群的外積，大到小）。"""
    per_prime = []
    for p, e in sorted(sp.factorint(n).items()):
        per_prime.append([tuple(p ** k for k in part) for part in _partitions(e)])
    out = []
    for combo in product(*per_prime):
        out.append(tuple(sorted(sum(combo, ()), reverse=True)))
    return sorted(out, key=lambda t: (len(t), [-x for x in t]))


def _partitions(e: int, mx: int | None = None) -> list:
    mx = e if mx is None else mx
    if e == 0:
        return [()]
    return [(k,) + rest for k in range(min(e, mx), 0, -1) for rest in _partitions(e - k, k)]


# =============================================================================
# §2.2–2.3 兩個耦合擺（H202–H204）
# =============================================================================
R2 = sp.Matrix([[0, 1], [1, 0]])
PEND2 = {"a": R_(21, 200), "b": 1}          # a = k/mℓ²、b = g/ℓ：ω₊ = 1、ω₋ = 11/10


def accel2(a=PEND2["a"], b=PEND2["b"]) -> sp.Matrix:
    """(2.2.1c)：[[a+b, −a], [−a, a+b]]。"""
    a, b = sp.nsimplify(a), sp.nsimplify(b)
    return sp.Matrix([[a + b, -a], [-a, a + b]])


def commutator(A, B):
    return A * B - B * A


def c2_projectors(R=R2) -> dict:
    """最小方程 R² − 1 = 0 的 Hamilton–Cayley 投影算符（第 1 章的同一招，sym_mat.projectors）：{1: 𝒫⁺, −1: 𝒫⁻}。"""
    return SM.projectors(R)


def c2_projectors_char(R=R2) -> dict:
    """由特徵標：𝒫^± = (1 ± R)/2（(2.3.2)_x）。"""
    E = sp.eye(R.shape[0])
    return {1: (E + R) / 2, -1: (E - R) / 2}


def normal2(a=PEND2["a"], b=PEND2["b"]) -> tuple:
    """(α₊, α₋) = (b, 2a + b)（p. 71）。"""
    a, b = sp.nsimplify(a), sp.nsimplify(b)
    return b, 2 * a + b


def beat_times2(a=PEND2["a"], b=PEND2["b"]) -> tuple:
    """完全轉移 π/Δω、拍的週期 2π/Δω（精確 sympy）。"""
    ap, am = normal2(a, b)
    dw = sp.sqrt(am) - sp.sqrt(ap)
    return sp.pi / dw, 2 * sp.pi / dw


class Modal:
    """ẍ = −𝒜x 的解析解：x(t) = Σ_j e_j (c_j cos ω_j t + s_j sin ω_j t/ω_j)（ω_j = 0 時 c_j + s_j t）。給 updater 用（純 numpy）。"""

    def __init__(self, A, x0, v0=None):
        A = _real(A)
        w2, V = np.linalg.eig(A)
        w2, V = w2.real, V.real
        self.V = V
        self.w = np.sqrt(np.clip(w2, 0, None))
        x0 = np.asarray(x0, dtype=float)
        v0 = np.zeros_like(x0) if v0 is None else np.asarray(v0, dtype=float)
        self.c = np.linalg.solve(V, x0)
        self.s = np.linalg.solve(V, v0)

    def __call__(self, t):
        t = np.asarray(t, dtype=float)
        w = self.w
        with np.errstate(invalid="ignore", divide="ignore"):
            sw = np.where(w > 1e-12, np.sin(np.multiply.outer(t, w)) / np.where(w > 1e-12, w, 1), np.multiply.outer(t, np.ones_like(w)))
        q = self.c * np.cos(np.multiply.outer(t, w)) + self.s * sw
        return q @ self.V.T

    def velocity(self, t):
        t = np.asarray(t, dtype=float)
        w = self.w
        q = -self.c * w * np.sin(np.multiply.outer(t, w)) + self.s * np.cos(np.multiply.outer(t, w))
        return q @ self.V.T


def ode_solve(A, x0, v0, ts) -> np.ndarray:
    """數值積分 ẍ = −𝒜x（solve_ivp, DOP853）→ 陣列 (len(ts), n)。"""
    A = _real(A)
    n = A.shape[0]
    y0 = np.concatenate([np.asarray(x0, float), np.zeros(n) if v0 is None else np.asarray(v0, float)])
    sol = integrate.solve_ivp(lambda t, y: np.concatenate([y[n:], -A @ y[:n]]), (ts[0], ts[-1]), y0,
                              t_eval=ts, method="DOP853", rtol=1e-11, atol=1e-12)
    return sol.y[:n].T


DETUNED = (R_(21, 200), 1, R_(221, 200))     # (a, b₁, b₂)：右擺擺長不同，沒有 C₂ 對稱


def detuned_matrix(a=DETUNED[0], b1=DETUNED[1], b2=DETUNED[2]) -> sp.Matrix:
    return sp.Matrix([[a + b1, -a], [-a, a + b2]])


def detuned_max_transfer(a=DETUNED[0], b1=DETUNED[1], b2=DETUNED[2]):
    """x(0) = (1, 0)：x₂(t) = (sin 2θ/2)(cos ω₁t − cos ω₂t)，最大振幅 sin 2θ，tan 2θ = 2a/(b₂ − b₁)。"""
    t2 = 2 * a / (b2 - b1)
    return sp.simplify(t2 / sp.sqrt(1 + t2 ** 2))


# =============================================================================
# §2.3(a)(b) Lissajous 與相量（H204）
# =============================================================================
LISSA = {"weak": R_(21, 20), "mid": R_(13, 10), "strong": sp.Integer(2)}     # ω₋/ω₊（ω₊ = 1，b = 1）


def a_for_ratio(rho, b=1):
    """ω₋/ω₊ = ρ ⇔ 2a + b = ρ²b ⇒ a = (ρ² − 1)b/2。"""
    return sp.nsimplify((sp.nsimplify(rho) ** 2 - 1) * b / 2)


def lissajous(t, rho):
    """從 x = (1, 0) 靜止放手：x₁ = (cos t + cos ρt)/2、x₂ = (cos t − cos ρt)/2（ω₊ = 1）。"""
    t = np.asarray(t, dtype=float)
    r = float(rho)
    return (np.cos(t) + np.cos(r * t)) / 2, (np.cos(t) - np.cos(r * t)) / 2


def lissajous_numeric(ts, rho) -> np.ndarray:
    return ode_solve(accel2(a_for_ratio(rho), 1), (1, 0), (0, 0), ts)


def closure(rho) -> dict:
    """(2.3.11)–(2.3.13)：ω₊/ω₋ = n₊/n₋（互質）⇒ t_LJ = 2πn₊/ω₊ = 2πn₋/ω₋ = (n₋ − n₊)·t_beat（拍數 n₋ − n₊）。"""
    rho = sp.nsimplify(rho)
    q = sp.Rational(1, 1) / rho                 # ω₊/ω₋
    n_plus, n_minus = int(sp.fraction(q)[0]), int(sp.fraction(q)[1])
    t_lj = 2 * sp.pi * n_plus                   # ω₊ = 1
    t_beat = 2 * sp.pi / (rho - 1)
    return {"n_plus": n_plus, "n_minus": n_minus, "t_LJ": t_lj, "t_beat": t_beat, "beats": n_minus - n_plus}


def quad_potential(a=PEND2["a"], b=PEND2["b"]) -> sp.Matrix:
    """(2.3.9) 的力常數矩陣 𝓕（m = ℓ = 1 時 = 加速度矩陣）；V = ½ xᵀ𝓕x。"""
    return accel2(a, b)


def ellipse_axes(a, b=1, V0=R_(1, 2)) -> tuple:
    """V = V0 的等高線：沿 (1,1) 的半軸 √(2V0/α₊)、沿 (1,−1) 的半軸 √(2V0/α₋)。"""
    ap, am = normal2(a, b)
    return sp.sqrt(2 * V0 / ap), sp.sqrt(2 * V0 / am)


PHASOR_NU = (R_(1, 2), sp.Integer(1))           # Fig. 2.3.7：(+) 模 ½ Hz、(−) 模 1 Hz
PHASOR_TIMES = (0, R_(1, 4), R_(1, 2), R_(3, 4), 1)


def phasors(t, nu=PHASOR_NU) -> tuple:
    """𝒞 = Ce^{−iωt}：z₁ = (e^{−iω₊t} + e^{−iω₋t})/2、z₂ = (e^{−iω₊t} − e^{−iω₋t})/2（t = 0：等量的 (+)、(−)，x = (1, 0)）。"""
    wp, wm = 2 * math.pi * float(nu[0]), 2 * math.pi * float(nu[1])
    t = float(t)
    a, b = np.exp(-1j * wp * t), np.exp(-1j * wm * t)
    return (a + b) / 2, (a - b) / 2


# =============================================================================
# §2.4–2.5 不可約表示、部分解題（H205）
# =============================================================================
C2_CHARS = {"+": {"1": 1, "R": 1}, "-": {"1": 1, "R": -1}}       # Table 2.4.1
O_R3 = sp.Matrix([[0, 0, 1], [0, 1, 0], [1, 0, 0]])                # (2.4.5)
O_REP3 = {"1": sp.eye(3), "R": O_R3}
T3 = sp.Matrix([[1 / sp.sqrt(2), 0, 1 / sp.sqrt(2)], [0, 1, 0], [1 / sp.sqrt(2), 0, -1 / sp.sqrt(2)]])   # (2.4.7) = (2.5.2)


def rep_projectors(reps: dict, chars: dict) -> dict:
    """𝒪(𝒫^α) = (1/°G) Σ_g D^α(g)* 𝒪(g)。"""
    n = len(reps)
    return {al: sp.simplify(sum((sp.conjugate(ch[g]) * reps[g] for g in reps), sp.zeros(*next(iter(reps.values())).shape)) / n)
            for al, ch in chars.items()}


def multiplicities(reps: dict, chars: dict) -> dict:
    """重複頻率 f^α = (1/°G) Σ_g D^α(g)* χ(g)，χ(g) = tr 𝒪(g)。"""
    n = len(reps)
    return {al: sp.simplify(sum(sp.conjugate(ch[g]) * reps[g].trace() for g in reps) / n) for al, ch in chars.items()}


PEND3 = {"m": 1, "M": 2, "k": 1, "l": 1, "g": 1}


def accel3(m=1, M=2, k=1, l=1, g=1) -> sp.Matrix:
    """(2.5.1c)：a = k/mℓ²、A = k/Mℓ²、b = g/ℓ。"""
    a, A, b = R_(k, m * l * l), R_(k, M * l * l), R_(g, l)
    return sp.Matrix([[a + b, -a, 0], [-A, 2 * A + b, -A], [0, -a, a + b]])


def reduce3(**kw) -> sp.Matrix:
    """(2.5.3)：T3†𝒜T3（部分約化：2×2 ⊕ 1×1）。"""
    return sp.simplify(T3.H * accel3(**kw) * T3)


def omega2_3(**kw) -> list:
    """三擺的 ω²（精確，由小到大）。"""
    return sorted(accel3(**kw).eigenvals(), key=float)


# =============================================================================
# §2.6–2.7 C_n：單位根、投影算符、特徵標表、單擺環（H206–H207）
# =============================================================================
def eps(n: int, k: int = 1):
    """e^{2πik/n}（精確，展成 cos + i sin）。"""
    return sp.nsimplify(sp.expand_complex(sp.exp(2 * sp.pi * I * sp.Rational(k % n, n)))) if n <= 4 else \
        sp.expand_complex(sp.exp(2 * sp.pi * I * sp.Rational(k % n, n)))


def cyclic_rep(n: int) -> sp.Matrix:
    """𝓡(r)：r|p⟩ = |p+1⟩（(2.6.7)：𝓡(r)x = (x_n, x_1, x_2, …)）。"""
    return sp.Matrix(n, n, lambda i, j: 1 if i == (j + 1) % n else 0)


def cn_char_table(n: int) -> sp.Matrix:
    """D^{k_n}(r^p) = e^{2πikp/n}：列 k = 0..n−1、行 p = 0..n−1（(2.6.5)、Fig. 2.7.2）。"""
    return sp.Matrix(n, n, lambda k, p: eps(n, k * p))


def cn_projectors_char(n: int) -> dict:
    """(2.6.6b)：𝒫^{k_n} = (1/n) Σ_m e^{−2πikm/n} r^m。"""
    R = cyclic_rep(n)
    return {k: sp.simplify(sp.expand_complex(sum((sp.conjugate(eps(n, k * m)) * R ** m for m in range(n)), sp.zeros(n)) / n))
            for k in range(n)}


def cn_projectors_minimal(n: int) -> dict:
    """(2.3.2)：𝒫^{ε_k} = Π_{l≠k}(r − ε_l)/Π_{l≠k}(ε_k − ε_l)，r^n − 1 = Π(r − ε_l)。"""
    R = cyclic_rep(n)
    E = sp.eye(n)
    out = {}
    for k in range(n):
        num, den = E, sp.Integer(1)
        for l in range(n):
            if l != k:
                num = num * (R - eps(n, l) * E)
                den = den * (eps(n, k) - eps(n, l))
        out[k] = sp.simplify(sp.expand_complex(num / den))
    return out


def cn_projectors_np(n: int, method: str = "char") -> dict:
    """浮點版（大 n 用）：method = "char"（(2.6.6b)）或 "minimal"（(2.3.2) 的連乘）。"""
    R = _real(cyclic_rep(n))
    z = np.exp(2j * np.pi * np.arange(n) / n)
    E = np.eye(n)
    out = {}
    for k in range(n):
        if method == "char":
            out[k] = sum(np.conj(z[k]) ** m * np.linalg.matrix_power(R, m) for m in range(n)) / n
        else:
            P = E.astype(complex)
            for l in range(n):
                if l != k:
                    P = P @ (R - z[l] * E) / (z[k] - z[l])
            out[k] = P
    return out


RING3 ={"a": 1, "b": 1}


def ring_matrix(n: int, a=1, b=0) -> sp.Matrix:
    """(2.6.1)、(2.7.1)：對角 2a + b、最近鄰 −a（環）。n = 2 時兩根彈簧並聯：[[2a+b, −2a], [−2a, 2a+b]]。"""
    a, b = sp.nsimplify(a), sp.nsimplify(b)
    M = sp.zeros(n)
    for p in range(n):
        M[p, p] += 2 * a + b
        M[p, (p + 1) % n] += -a
        M[p, (p - 1) % n] += -a
    return M


def ring_omega2(n: int, a=1, b=0) -> list:
    """(2.7.4)：ω²(k_n) = b + 2a(1 − cos 2πk/n)，k = 0..n−1（精確）。"""
    a, b = sp.nsimplify(a), sp.nsimplify(b)
    return [sp.expand(b + 2 * a * (1 - sp.cos(2 * sp.pi * sp.Rational(k, n)))) for k in range(n)]


def wave_vector(n: int, k: int) -> np.ndarray:
    """(2.7.3)：|k_n⟩ 的分量 e^{−ik_n p}/√n。"""
    p = np.arange(n)
    return np.exp(-2j * np.pi * k * p / n) / np.sqrt(n)


def cos_wave(n: int, k: int) -> np.ndarray:
    """(2.7.5a)：√(2/n) cos(k_n p)（k = 0 或 n/2 時再除以 √2 才歸一）。"""
    return np.sqrt(2 / n) * np.cos(2 * np.pi * k * np.arange(n) / n)


def sin_wave(n: int, k: int) -> np.ndarray:
    """(2.7.5b)：√(2/n) sin(k_n p)。"""
    return np.sqrt(2 / n) * np.sin(2 * np.pi * k * np.arange(n) / n)


C3_STANDING = (sp.Matrix([2, -1, -1]) / sp.sqrt(6), sp.Matrix([0, -1, 1]) / sp.sqrt(2))
C3_MOVING = sp.Matrix([1, sp.conjugate(eps(3)), eps(3)]) / sp.sqrt(3)       # |1₃⟩ = |e₂⟩


def dispersion(k, a=1, b=0):
    """(2.7.6)：ω(k) = √(2a(1 − cos k) + b)。"""
    k = np.asarray(k, dtype=float)
    return np.sqrt(2 * float(a) * (1 - np.cos(k)) + float(b))


def dispersion_series(k, a=1):
    return math.sqrt(a) * k - math.sqrt(a) * k ** 3 / 24


def brillouin(n: int) -> list:
    """第一 Brillouin 區的 n 個整數 k：−n/2 < k ≤ n/2（n 偶）、−(n−1)/2 ≤ k ≤ (n−1)/2（n 奇）。"""
    return list(range(-((n - 1) // 2), n // 2 + 1))


def alias(n: int, k: int) -> int:
    """把 k 換成第一 Brillouin 區裡等價的整數（(2.7.7)）。"""
    z = k % n
    return z - n if z > n // 2 else z


ALIAS = (5, (2, 7, -3))


def alias_values(n: int = 5, k: int = 2) -> np.ndarray:
    return np.cos(2 * np.pi * k * np.arange(n) / n)


# =============================================================================
# §2.7A–B 對稱破缺與奔馳波（H208）
# =============================================================================
BROKEN = {"n": 24, "abar": R_(3, 2), "alow": R_(1, 2), "b": 0}


def broken_ring_matrix(n=24, abar=R_(3, 2), alow=R_(1, 2), b=0) -> np.ndarray:
    """(2.7.8)：彈簧 p ↔ p+1 的常數 a̲（p 偶）或 ā（p 奇）。"""
    M = np.zeros((n, n))
    for p in range(n):
        s = float(alow) if p % 2 == 0 else float(abar)
        q = (p + 1) % n
        M[p, p] += s
        M[q, q] += s
        M[p, q] -= s
        M[q, p] -= s
    return M + float(b) * np.eye(n)


def broken_block(k, abar=R_(3, 2), alow=R_(1, 2), b=0) -> np.ndarray:
    """(2.7.9)（對角加回 b）：基底 |k⟩（偶數格點）、|k'⟩（奇數格點）。"""
    ab, al = float(abar), float(alow)
    off = -(ab * np.exp(1j * k) + al)
    return np.array([[ab + al + float(b), off], [np.conj(off), ab + al + float(b)]])


def broken_bands(k, abar=R_(3, 2), alow=R_(1, 2), b=0) -> tuple:
    """(2.7.10)（加回 b）：ω² = b + ā + a̲ ∓ √(ā² + a̲² + 2āa̲ cos k)。回傳（聲學, 光學）。"""
    ab, al, b = float(abar), float(alow), float(b)
    r = np.sqrt(ab * ab + al * al + 2 * ab * al * np.cos(k))
    return b + ab + al - r, b + ab + al + r


def broken_standing(n: int = 24) -> tuple:
    """k = π（C_{n/2} 的 6₁₂）的兩個駐波：(1, −1, −1, 1, …) 只拉軟彈簧（ω² = 2a̲ + b）、(1, 1, −1, −1, …) 只拉硬彈簧（2ā + b）；
    格點 0、1 之間是軟彈簧（(2.7.11b) 的 s、c，相差一個原點的選法）。"""
    p = np.arange(n)
    soft = np.sqrt(2) * np.cos(np.pi * p / 2 + np.pi / 4) / np.sqrt(n)
    hard = np.sqrt(2) * np.cos(np.pi * p / 2 - np.pi / 4) / np.sqrt(n)
    return soft, hard


GALLOP ={"a": (R_(6, 5), R_(3, 5)), "b": (sp.Integer(1), R_(9, 10))}   # Fig. 2.7.8：(I, R)
SWR_LIST = (1, 0.5, 0.1, 0, -0.1, -0.5, -1)                           # Fig. 2.7.9


def swr(I_, R_amp):
    """(2.7.13)：Δ = (I − R)/(I + R) = |ψ| 的極小/極大。"""
    return sp.nsimplify((I_ - R_amp) / (I_ + R_amp))


def gallop_psi(x, t, I_, R_amp, k=1.0, w=1.0):
    """(2.7.12)：ψ = (Ie^{ikx} + Re^{−ikx})e^{−iωt}。"""
    x = np.asarray(x, dtype=float)
    return (float(I_) * np.exp(1j * k * x) + float(R_amp) * np.exp(-1j * k * x)) * np.exp(-1j * w * t)


def gallop_envelope(x, I_, R_amp, k=1.0):
    """(2.7.14)：|ψ| = [I² + R² + 2IR cos 2kx]^{1/2}。"""
    I_, R_amp = float(I_), float(R_amp)
    return np.sqrt(I_ * I_ + R_amp * R_amp + 2 * I_ * R_amp * np.cos(2 * k * np.asarray(x, float)))


def gallop_zeros(t, delta, k=1.0, w=1.0, span=2 * math.pi) -> np.ndarray:
    """Re ψ = 0 的位置（I = (1+Δ)/2、R = (1−Δ)/2）：tan kx = −cot(ωt)/Δ；在 [0, span) 內。"""
    d = float(delta)
    wt = w * t
    if abs(d) < 1e-15:
        base = math.pi / 2
    else:
        base = math.atan2(-math.cos(wt), d * math.sin(wt))
    xs = (base + np.pi * np.arange(-4, 2 * int(span / math.pi) + 4)) / k
    return np.sort(xs[(xs >= 0) & (xs < span)])


def zero_velocity(t, delta, c=1.0, w=1.0):
    """零點的速度 v = cΔ/(Δ² sin²ωt + cos²ωt)：在 cΔ（ωt = 0）與 c/Δ（ωt = π/2）之間（習題 2.7.2）。"""
    d = float(delta)
    return c * d / (d * d * math.sin(w * t) ** 2 + math.cos(w * t) ** 2)


# =============================================================================
# §2.8–2.11 D₂、外積、8 階 Abel 群、點群（H209）
# =============================================================================
D2_TABLE = table_from_rows(D2_NAMES, D2_TABLE_BOOK)


def regular_rep(table: dict, names) -> dict:
    """正則表示：g|h⟩ = |gh⟩（基底順序 names）。"""
    n = len(names)
    return {g: sp.Matrix(n, n, lambda i, j: 1 if table[(g, names[j])] == names[i] else 0) for g in names}


D2_REP4 = regular_rep(D2_TABLE, D2_NAMES)                     # (2.8.2b)，基底 |1⟩ = |1⟩、|2⟩ = |Rx⟩、|3⟩ = |Ry⟩、|4⟩ = |Rz⟩
_D2_BOOK_ACTION = {"1": (1, 2, 3, 4), "Rx": (2, 1, 4, 3), "Ry": (3, 4, 1, 2), "Rz": (4, 3, 2, 1)}   # (2.8.2a)：g|j⟩ 的編號
D2_REP4_BOOK = {g: sp.Matrix(4, 4, lambda i, j: 1 if act[j] == i + 1 else 0) for g, act in _D2_BOOK_ACTION.items()}
D2_CHARS = {"A1": {"1": 1, "Rx": 1, "Ry": 1, "Rz": 1}, "B2": {"1": 1, "Rx": -1, "Ry": 1, "Rz": -1},
            "B1": {"1": 1, "Rx": 1, "Ry": -1, "Rz": -1}, "A2": {"1": 1, "Rx": -1, "Ry": -1, "Rz": 1}}   # (2.8.7)


def d2_projectors_product() -> dict:
    """(2.8.3)–(2.8.5b)：𝒫¹ = 𝒫^{x+}𝒫^{y+}、𝒫² = 𝒫^{x−}𝒫^{y+}、𝒫³ = 𝒫^{x+}𝒫^{y−}、𝒫⁴ = 𝒫^{x−}𝒫^{y−}（以 A1、B2、B1、A2 命名）。"""
    E = sp.eye(4)
    px = {1: (E + D2_REP4["Rx"]) / 2, -1: (E - D2_REP4["Rx"]) / 2}
    py = {1: (E + D2_REP4["Ry"]) / 2, -1: (E - D2_REP4["Ry"]) / 2}
    return {"A1": px[1] * py[1], "B2": px[-1] * py[1], "B1": px[1] * py[-1], "A2": px[-1] * py[-1]}


def d2_projectors_char() -> dict:
    return rep_projectors(D2_REP4, D2_CHARS)


D2_ABC = (R_(11, 4), R_(-1, 2), sp.Integer(-1), R_(-1, 4))   # (A, a, b, c)：g/ℓ = 1，k_a = ½、k_b = 1、k_c = ¼


def d2_matrix(A=D2_ABC[0], a=D2_ABC[1], b=D2_ABC[2], c=D2_ABC[3]) -> sp.Matrix:
    """(2.8.1a)（ẍ = −𝒜x）：[[A, a, b, c], [a, A, c, b], [b, c, A, a], [c, b, a, A]]。"""
    return sp.Matrix([[A, a, b, c], [a, A, c, b], [b, c, A, a], [c, b, a, A]])


def d2_omega2(A=D2_ABC[0], a=D2_ABC[1], b=D2_ABC[2], c=D2_ABC[3]) -> dict:
    """Fig. 2.8.2：ω²(α) = Σ_g D^α(g)·(第一列的係數)。"""
    coef = {"1": A, "Rx": a, "Ry": b, "Rz": c}
    return {al: sum(ch[g] * coef[g] for g in coef) for al, ch in D2_CHARS.items()}


def abelian_group(orders) -> tuple:
    """C_{n1} × C_{n2} × …：元素是指數組 (e1, e2, …)，乘法＝逐分量加（mod n_i）。回傳 (elements, table)。"""
    els = list(product(*[range(n) for n in orders]))
    tab = {(g, h): tuple((x + y) % n for x, y, n in zip(g, h, orders)) for g in els for h in els}
    return els, tab


def abelian_char_table(orders) -> tuple:
    """外積的特徵標表：χ_k(g) = Π_i ε_{n_i}^{k_i g_i}。回傳 (elements, irreps, sympy Matrix[列 irrep, 行 元素])。"""
    els = list(product(*[range(n) for n in orders]))
    M = sp.Matrix(len(els), len(els), lambda a, b: sp.prod([eps(n, k * g) for k, g, n in zip(els[a], els[b], orders)]))
    return els, els, M


def outer_char_table(T1: sp.Matrix, T2: sp.Matrix) -> sp.Matrix:
    """(2.10)：H × K 的特徵標表 = 兩張表的 Kronecker 積（列 (α, β)、行 (h, k)）。"""
    return sp.Matrix(sp.kronecker_product(T1, T2))


ABELIAN8 = {"C8": (8,), "C4xC2": (4, 2), "C2xC2xC2": (2, 2, 2)}


def is_cyclic(orders) -> bool:
    """有一個元素的階 = 群的階 ⇔ 循環（⇔ 各 n_i 兩兩互質）。"""
    els, tab = abelian_group(orders)
    N = len(els)
    return any(element_order(tab, g) == N for g in els)


POINT3 = {"1": sp.eye(3), "I": -sp.eye(3), "sxy": sp.diag(1, 1, -1), "sxz": sp.diag(1, -1, 1), "syz": sp.diag(-1, 1, 1),
          "Rx": ROT3["Rx"], "Ry": ROT3["Ry"], "Rz": ROT3["Rz"],
          "S4": sp.Matrix([[0, 1, 0], [-1, 0, 0], [0, 0, -1]])}      # S4 = I·Rz(90°)：90° 轉動反演（C_{2s}）
C2H = ("1", "Rz", "I", "sxy")
C2V = ("1", "Rz", "sxz", "syz")
C2V_BOOK = ("1", "Rz", "sxy", "syz")
CRYSTAL_POINT_GROUPS, ABELIAN_CRYSTAL_POINT_GROUPS = 32, 16


def is_group(names) -> bool:
    mats = [POINT3[x] for x in names]
    return all(any(A * B == C for C in mats) for A in mats for B in mats)


def generated_order(M: sp.Matrix) -> int:
    k, X = 1, M
    while X != sp.eye(M.shape[0]):
        X = X * M
        k += 1
    return k


def point_table(names) -> dict:
    """點群的乘法表（由 POINT3 的矩陣）。"""
    return table_from_matrices({x: POINT3[x] for x in names})


# =============================================================================
# §2.12A C₁₂ 轉子：自由、弱位能、緊束縛、Mathieu、Kronig–Penney（H210）
# =============================================================================
ROTOR_N = 12
WEAK_V, STRONG_V = 2, 360                   # V(φ) = −V cos 12φ（E = ħ²/2μr² = 1）


def free_levels(mmax: int = 12) -> dict:
    """(2.12.5)：ε_m = m²E。"""
    return {m: m * m for m in range(-mmax, mmax + 1)}


def rotor_matrix(V=WEAK_V, mmax: int = 60, n: int = ROTOR_N) -> tuple:
    """(2.12.12)：⟨m'|H|m⟩ = δ m² + G(m' − m)，G(±n) = −V/2。回傳 (ms, H)。"""
    ms = np.arange(-mmax, mmax + 1)
    H = np.diag(ms.astype(float) ** 2)
    for i in range(len(ms)):
        for j in range(len(ms)):
            if abs(ms[i] - ms[j]) == n:
                H[i, j] = -float(V) / 2
    return ms, H


def rotor_block(kappa: int, V=WEAK_V, N: int = 8, n: int = ROTOR_N) -> tuple:
    """(2.12.16)：只含 m = κ + nN'（|N'| ≤ N）的子矩陣（C_n 的 irrep κ）。回傳 (ms, H)。"""
    ms = np.array([kappa + n * j for j in range(-N, N + 1)])
    H = np.diag(ms.astype(float) ** 2)
    for i in range(len(ms) - 1):
        H[i, i + 1] = H[i + 1, i] = -float(V) / 2
    return ms, H


def rotor_bands(V=WEAK_V, nb: int = 3, n: int = ROTOR_N, N: int = 8) -> dict:
    """{κ: 最低 nb 個本徵值}，κ 走第一 Brillouin 區（−5..6）。"""
    return {k: np.linalg.eigvalsh(rotor_block(k, V, N, n)[1])[:nb] for k in brillouin(n)}


def pair_exact(m: int = 5, G=-1, n: int = ROTOR_N) -> tuple:
    """(2.12.17)：|m⟩、|m−n⟩ 的 2×2：(m² + (m−n)²)/2 ∓ √(((m−n)² − m²)²/4 + G²)。"""
    e1, e2 = m * m, (m - n) ** 2
    c, h = sp.Rational(e1 + e2, 2), sp.Rational(abs(e2 - e1), 2)
    r = sp.sqrt(h * h + sp.nsimplify(G) ** 2)
    return c - r, c + r


def pair_pert(m: int = 5, G=-1, n: int = ROTOR_N) -> tuple:
    """(2.12.17) 的近似：m² − G²/((m−n)² − m²)、(m−n)² + G²/(…)。"""
    e1, e2 = m * m, (m - n) ** 2
    d = sp.Integer(e2 - e1)
    G2 = sp.nsimplify(G) ** 2
    return e1 - G2 / d, e2 + G2 / d


def tb_matrix(n: int, H=sp.Symbol("H"), S=sp.Symbol("S")) -> sp.Matrix:
    """(2.12.23) 第二式：對角 H、最近鄰 −S（環）。"""
    M = sp.zeros(n)
    for p in range(n):
        M[p, p] = H
        M[p, (p + 1) % n] += -S
        M[p, (p - 1) % n] += -S
    return M


def tb_levels(n: int, H=sp.Symbol("H"), S=sp.Symbol("S")) -> list:
    """(2.12.25)：ε_k = H − 2S cos(2πk/n)，k = 0..n−1。"""
    return [sp.nsimplify(H - 2 * S * sp.cos(2 * sp.pi * sp.Rational(k, n))) for k in range(n)]


def tight_fit(V=STRONG_V, n: int = ROTOR_N) -> dict:
    """最低帶的 n 個能階（每個 irrep κ 一個）→ H = 平均、S = (ε(π) − ε(0))/4、與 H − 2S cos k 的最大偏差。"""
    lev = {k: rotor_bands(V, 1, n)[k][0] for k in brillouin(n)}
    e0, epi = lev[0], lev[n // 2]
    S = (epi - e0) / 4
    H = (epi + e0) / 2
    dev = max(abs(lev[k] - (H - 2 * S * math.cos(2 * math.pi * k / n))) for k in lev)
    return {"levels": lev, "H": H, "S": S, "width": epi - e0, "dev": dev, "mean": float(np.mean(list(lev.values())))}


def mathieu_edges(Vp: float, nb: int = 3, N: int = 30) -> list:
    """縮放的 Mathieu 方程 ψ'' + (ε + V' cos 2t)ψ = 0 的帶邊：A（週期 π，基底 e^{2ijt}）與 B（反週期，e^{i(2j+1)t}），
    排序後兩兩一帶：[(下緣, 上緣), …]（共 nb 帶）。"""
    def ev(qs):
        H = np.diag(np.array(qs, float) ** 2)
        for i in range(len(qs) - 1):
            H[i, i + 1] = H[i + 1, i] = -Vp / 2
        return np.linalg.eigvalsh(H)
    a_lv = ev([2 * j for j in range(-N, N + 1)])
    b_lv = ev([2 * j + 1 for j in range(-N - 1, N + 1)])
    allv = np.sort(np.concatenate([a_lv, b_lv]))[: 2 * nb]
    return [(allv[2 * j], allv[2 * j + 1]) for j in range(nb)]


def mathieu_edges_scipy(Vp: float, nb: int = 3) -> list:
    """同上，用 scipy.special.mathieu_a/b（q = V'/2，帶 j = [a_j, b_{j+1}]）。"""
    q = Vp / 2
    return [(special.mathieu_a(j, q), special.mathieu_b(j + 1, q)) for j in range(nb)]


def rotor_edges(V: float, nb: int = 3, n: int = ROTOR_N) -> list:
    """C_n 轉子（E = 1）的帶邊 = (n/2)²·Mathieu 帶邊（V' = V/(n/2)²）：n = 12 時乘 36（Fig. 2.12.5 的縱軸）。"""
    s = (n / 2) ** 2
    return [(s * lo, s * hi) for lo, hi in mathieu_edges(V / s, nb)]


def square_well_G(m: int, U=1, n: int = ROTOR_N):
    """n 方井位能 V(φ) = U(nφ)（|θ| < π/2 為 0、其餘 U）的 Fourier 係數 G(m) = (1/2π)∫e^{−imφ}V dφ（精確）。"""
    U = sp.nsimplify(U)
    if m == 0:
        return U / 2
    if m % n:
        return sp.Integer(0)
    j = abs(m) // n
    if j % 2 == 0:
        return sp.Integer(0)
    N = (j - 1) // 2
    return (-1) ** (N + 1) * U / (j * sp.pi)


def square_well_G_fft(m: int, U=1.0, n: int = ROTOR_N, M: int = 4096 * 3) -> complex:
    phi = (np.arange(M) + 0.5) * 2 * np.pi / M
    th = np.mod(n * phi + np.pi, 2 * np.pi) - np.pi
    V = np.where(np.abs(th) < np.pi / 2, 0.0, float(U))
    return complex(np.mean(np.exp(-1j * m * phi) * V))


KP_W = KP_B = math.pi / 2                 # 井寬、位障寬（週期 π）
KP_CROSSINGS = ((12, 16), (20, 36), (32, 36), (28, 64), (48, 64), (60, 64))   # (U, E)：E = (2n_w)²、U = E − (2n_p)²


def seg_matrix(E: float, U: float, L: float) -> np.ndarray:
    """常數位能區段的 (ψ, ψ') 傳遞矩陣（ψ'' = (U − E)ψ）。"""
    d = E - U
    if d > 1e-14:
        k = math.sqrt(d)
        return np.array([[math.cos(k * L), math.sin(k * L) / k], [-k * math.sin(k * L), math.cos(k * L)]])
    if d < -1e-14:
        q = math.sqrt(-d)
        return np.array([[math.cosh(q * L), math.sinh(q * L) / q], [q * math.sinh(q * L), math.cosh(q * L)]])
    return np.array([[1.0, L], [0.0, 1.0]])


def kp_half_trace(E: float, U: float, w: float = KP_W, b: float = KP_B) -> float:
    """½ tr(M_位障 M_井) = cos(qL)（Kronig–Penney 方程；|·| ≤ 1 是能帶）。"""
    M = seg_matrix(E, U, b) @ seg_matrix(E, 0.0, w)
    return 0.5 * float(np.trace(M))


def kp_half_trace_ode(E: float, U: float, w: float = KP_W, b: float = KP_B) -> float:
    """同上，用 solve_ivp 積分一個週期（單擺式的 monodromy）。"""
    def f(x, y):
        Ux = 0.0 if x < w else U
        return [y[1], (Ux - E) * y[0]]
    tr = 0.0
    for k, y0 in enumerate(([1.0, 0.0], [0.0, 1.0])):
        s1 = integrate.solve_ivp(f, (0, w), y0, rtol=1e-11, atol=1e-12, method="DOP853")
        s2 = integrate.solve_ivp(f, (w, w + b), s1.y[:, -1], rtol=1e-11, atol=1e-12, method="DOP853")
        tr += s2.y[k, -1]
    return 0.5 * tr


def kp_edges(U: float, Emax: float = 60.0, n: int = 6000) -> list:
    """½ tr = ±1 的所有根（E < Emax），由小到大。"""
    Es = np.linspace(1e-6, Emax, n)
    out = []
    for target in (1.0, -1.0):
        f = np.array([kp_half_trace(E, U) - target for E in Es])
        for i in np.flatnonzero(np.sign(f[:-1]) * np.sign(f[1:]) < 0):
            out.append(optimize.brentq(lambda E: kp_half_trace(E, U) - target, Es[i], Es[i + 1], xtol=1e-13))
    return sorted(out)


def kp_gap_near(U: float, E0: float, half: float = 1.5) -> float:
    """E0 附近最窄的「帶外」區間寬度（偶然簡併處 → 0）：在 [E0 − half, E0 + half] 內 |½tr| > 1 的最長連續段長度。"""
    Es = np.linspace(E0 - half, E0 + half, 2001)
    out = np.array([abs(kp_half_trace(E, U)) > 1 + 1e-12 for E in Es])
    best, run = 0, 0
    for o in out:
        run = run + 1 if o else 0
        best = max(best, run)
    return best * (Es[1] - Es[0])


# =============================================================================
# §2.12B 兩能階系統：氨（H211）
# =============================================================================
NH3 = {"H": 0, "S": 1}
NH3_NU = 23.8701e9                        # Hz，NH₃ (J,K) = (3,3) 反轉線（氨邁射）——引用值
H_PLANCK_EV = 4.135667696e-15             # eV·s
LZ_BETAS = (4.0, 0.25)


def two_level(pE, H=0, S=1) -> tuple:
    """(2.12.34)：ε = H ∓ √(S² + (pE)²)。"""
    r = sp.sqrt(sp.nsimplify(S) ** 2 + sp.nsimplify(pE) ** 2)
    return H - r, H + r


def two_level_matrix(pE, H=0, S=1) -> np.ndarray:
    """(2.12.33)：[[H − pE, −S], [−S, H + pE]]（基底 up、dn）。"""
    return np.array([[H - float(pE), -float(S)], [-float(S), H + float(pE)]])


def tunnel_populations(t, S=1) -> tuple:
    """從 |up⟩ 出發：P_up = cos²(St/ħ)、P_dn = sin²(St/ħ)（拍頻 2S/ħ）。"""
    t = np.asarray(t, dtype=float)
    return np.cos(float(S) * t) ** 2, np.sin(float(S) * t) ** 2


def tunnel_numeric(ts, H=0, S=1) -> np.ndarray:
    """|ψ(t)⟩ = expm(−i𝓗t)|up⟩ 的 P_dn。"""
    Hm = two_level_matrix(0, H, S)
    return np.array([abs((linalg.expm(-1j * Hm * t) @ np.array([1, 0]))[1]) ** 2 for t in ts])


def nh3_numbers() -> dict:
    return {"2S_eV": H_PLANCK_EV * NH3_NU, "period_ps": 1e12 / NH3_NU, "wavenumber_cm": NH3_NU / 2.99792458e10}


def lz_formula(beta: float, S: float = 1.0) -> float:
    """Landau–Zener（書外）：pE = βt，對角差 2βt ⇒ 跳到另一支的機率 e^{−πS²/β}。"""
    return math.exp(-math.pi * S * S / beta)


def lz_numeric(beta: float, S: float = 1.0, T: float | None = None) -> float:
    """數值積分 iψ' = 𝓗(t)ψ，pE = βt，t ∈ [−T, T]；從下面那支的本徵態出發，最後在上面那支的機率。"""
    if T is None:
        T = 80.0 / beta if beta > 1 else 60.0
    # 中點 Magnus 積分：每一步 U = exp(−i𝓗(t_mid)dt)，𝓗 = −Sσx − βtσz（2×2 閉式）
    dt = 0.25 / (1 + beta * T) / 4
    nst = int(2 * T / dt)
    dt = 2 * T / nst
    tm = -T + (np.arange(nst) + 0.5) * dt
    hx, hz = -S * np.ones(nst), -beta * tm
    hn = np.hypot(hx, hz)
    c, sn = np.cos(hn * dt), np.sin(hn * dt) / hn
    psi = np.linalg.eigh(two_level_matrix(-beta * T, 0, S))[1][:, 0].astype(complex)
    a0, a1 = psi[0], psi[1]
    for j in range(nst):
        u, x_, z_ = c[j], sn[j] * hx[j], sn[j] * hz[j]
        a0, a1 = (u - 1j * z_) * a0 - 1j * x_ * a1, -1j * x_ * a0 + (u + 1j * z_) * a1
    v = np.linalg.eigh(two_level_matrix(beta * T, 0, S))[1]
    return float(abs(np.vdot(v[:, 1], np.array([a0, a1]))) ** 2)


# =============================================================================
# §2.12C–D 散射、穿越矩陣、WKB（H212）
# =============================================================================
WELL = {"V": 25, "a": 1}                  # ħ = 1、2m = 1；井內 0、井外 V；e = E − V


def wavenumbers(e, V=WELL["V"]) -> tuple:
    """(2.12.37)（更正為 /ħ）：l = √(2mE)/ħ = √(e + V)（井內）、k = √(2m(E − V))/ħ = √e（井外）。"""
    return math.sqrt(e + V), math.sqrt(e)


def s_matrix(e, V=WELL["V"], a=WELL["a"], book: bool = False) -> np.ndarray:
    """(2.12.38)：S = [[i(l² − k²) sin 2la, 2lk], [2lk, i(l² − k²) sin 2la]]/D，
    D = e^{2ika}(2lk cos 2la − i(l² + k²) sin 2la)；book=True 用書上多一個負號的 D（得 −S）。"""
    l, k = wavenumbers(e, V)
    D = np.exp(2j * k * a) * (2 * l * k * np.cos(2 * l * a) - 1j * (l * l + k * k) * np.sin(2 * l * a))
    if book:
        D = -D
    s11 = 1j * (l * l - k * k) * np.sin(2 * l * a) / D
    return np.array([[s11, 2 * l * k / D], [2 * l * k / D, s11]])


def s_matrix_direct(e, V=WELL["V"], a=WELL["a"]) -> np.ndarray:
    """直接解 x = ±a 的四個邊界條件（ψ₁ = I₁e^{ikx} + O₁e^{−ikx}、ψ_w = Ae^{ilx} + Be^{−ilx}、ψ₂ = I₂e^{−ikx} + O₂e^{ikx}）。"""
    l, k = wavenumbers(e, V)
    S = np.zeros((2, 2), complex)
    for col, (I1, I2) in enumerate(((1, 0), (0, 1))):
        M = np.zeros((4, 4), complex)
        r = np.zeros(4, complex)
        x = -a
        M[0] = [np.exp(-1j * k * x), 0, -np.exp(1j * l * x), -np.exp(-1j * l * x)]
        r[0] = -I1 * np.exp(1j * k * x)
        M[1] = [-1j * k * np.exp(-1j * k * x), 0, -1j * l * np.exp(1j * l * x), 1j * l * np.exp(-1j * l * x)]
        r[1] = -I1 * 1j * k * np.exp(1j * k * x)
        x = a
        M[2] = [0, np.exp(1j * k * x), -np.exp(1j * l * x), -np.exp(-1j * l * x)]
        r[2] = -I2 * np.exp(-1j * k * x)
        M[3] = [0, 1j * k * np.exp(1j * k * x), -1j * l * np.exp(1j * l * x), 1j * l * np.exp(-1j * l * x)]
        r[3] = I2 * 1j * k * np.exp(-1j * k * x)
        s = np.linalg.solve(M, r)
        S[0, col], S[1, col] = s[0], s[1]
    return S


def phase_shifts(e, V=WELL["V"], a=WELL["a"]) -> tuple:
    """偶：井內 cos lx、井外 cos(kx + δ₊) ⇒ tan(ka + δ₊) = (l/k) tan la；奇：sin lx、−cos(kx + δ₋) ⇒ tan(ka + δ₋) = −(l/k) cot la。
    回傳 (δ₊, δ₋)，各在 (−π/2, π/2]（mod π）。"""
    l, k = wavenumbers(e, V)
    dp = math.atan(l / k * math.tan(l * a)) - k * a
    dm = math.atan(-l / k / math.tan(l * a)) - k * a
    wrap = lambda d: (d + math.pi / 2) % math.pi - math.pi / 2     # noqa: E731
    return wrap(dp), wrap(dm)


def transmission(e, V=WELL["V"], a=WELL["a"]):
    """|S₂₁|² = 1/(1 + V² sin²(2la)/(4e(e + V)))。"""
    e = np.asarray(e, dtype=float)
    l = np.sqrt(e + V)
    return 1 / (1 + V * V * np.sin(2 * l * a) ** 2 / (4 * e * (e + V)))


def resonances(V=WELL["V"], a=WELL["a"], emax: float = 70) -> list:
    """2la = nπ ⇒ e = (nπ/2a)² − V（T = 1，Ramsauer–Townsend）。"""
    out = []
    for nn in range(1, 100):
        e = (nn * math.pi / (2 * a)) ** 2 - V
        if 0 < e < emax:
            out.append(e)
    return out


def crossing_step(g1, g2, x) -> np.ndarray:
    """(2.12.50b)：(R₁, L₁) = C (R₂, L₂)，ψ_j = R_j e^{g_j x} + L_j e^{−g_j x}，g_j = [2m(V_j − E)/ħ²]^{1/2}（更正 ħ²）。"""
    g1, g2 = complex(g1), complex(g2)
    return np.array([[np.exp(-(g1 - g2) * x) * (g1 + g2) / (2 * g1), np.exp(-(g1 + g2) * x) * (g1 - g2) / (2 * g1)],
                     [np.exp((g1 + g2) * x) * (g1 - g2) / (2 * g1), np.exp((g1 - g2) * x) * (g1 + g2) / (2 * g1)]])


def crossing_well(e, V=WELL["V"], a=WELL["a"]) -> np.ndarray:
    """(2.12.51)：(I₁, O₁) = C (O₂, I₂)。"""
    l, k = wavenumbers(e, V)
    s, c = np.sin(2 * l * a), np.cos(2 * l * a)
    p, m = (l * l + k * k) / (2 * k * l), (l * l - k * k) / (2 * k * l)
    return np.array([[np.exp(2j * k * a) * (c - 1j * p * s), -1j * m * s],
                     [1j * m * s, np.exp(-2j * k * a) * (c + 1j * p * s)]])


def crossing_well_steps(e, V=WELL["V"], a=WELL["a"]) -> np.ndarray:
    """兩個 (2.12.50b) 的乘積：外 → 井（x = −a）、井 → 外（x = a）；g = ik（外）、il（井）。"""
    l, k = wavenumbers(e, V)
    return crossing_step(1j * k, 1j * l, -a) @ crossing_step(1j * l, 1j * k, a)


def s_from_c(C) -> np.ndarray:
    """(2.12.52a)（更正 S₁₁ = C₂₁/C₁₁）。"""
    d = np.linalg.det(C)
    return np.array([[C[1, 0] / C[0, 0], d / C[0, 0]], [1 / C[0, 0], -C[0, 1] / C[0, 0]]])


def c_from_s(S) -> np.ndarray:
    """(2.12.52b)。"""
    d = np.linalg.det(S)
    return np.array([[-np.conj(d) / np.conj(S[1, 0]), np.conj(S[0, 0]) / np.conj(S[1, 0])],
                     [S[0, 0] / S[1, 0], -d / S[1, 0]]])


def bound_geometric(V=25.0, w=2.0) -> list:
    """習題 2.12.2（ħ = 1、m = ½）：井內 sin(kx + δ)，sin δ = k/√V，k = (nπ − 2δ)/w。回傳 [(k, δ, E = k²)]。"""
    out = []
    for nn in range(1, 50):
        f = lambda d: math.sqrt(V) * math.sin(d) - (nn * math.pi - 2 * d) / w      # noqa: E731
        if f(math.pi / 2) < 0:
            break
        d = optimize.brentq(f, 0.0, math.pi / 2, xtol=1e-14)
        k = math.sqrt(V) * math.sin(d)
        out.append((k, d, k * k))
    return out


def bound_fd(V=25.0, w=2.0, pad: float = 10.0, N: int = 4000) -> np.ndarray:
    """有限差分（Richardson 外插）：−ψ'' + U ψ = Eψ，井 [0, w] 內 U = 0、外面 U = V；回傳 E < V 的本徵值。"""
    def ev(M):                                     # 格點在格子中央，井的邊界落在格子邊界上（二階精度）
        h = w / M
        npad = int(round(pad / h))
        x = (np.arange(-npad, M + npad) + 0.5) * h
        U = np.where((x > 0) & (x < w), 0.0, V)
        n_ = len(x)
        return linalg.eigh_tridiagonal(2 / h ** 2 + U, -np.ones(n_ - 1) / h ** 2, select="v", select_range=(-1, V - 1e-6))[0]
    e1, e2 = ev(N // 2), ev(N)
    n = min(len(e1), len(e2))
    return (4 * e2[:n] - e1[:n]) / 3


def c_barrier(theta) -> np.ndarray:
    """(2.12.54a)。"""
    r = math.sqrt(1 + theta * theta)
    return np.array([[r, 1j * theta], [-1j * theta, r]])


def c_valley(alpha) -> np.ndarray:
    """(2.12.55a)。"""
    return np.array([[np.exp(1j * alpha), 0], [0, np.exp(-1j * alpha)]])


def wkb_cluster(Ej, hw, theta, n, j=0) -> list:
    """(2.12.60)：E = E_j − (−1)^j (ħω/π)[1 + θ²]^{−1/2} cos(2πm/n)，m = 0..n−1。"""
    return [Ej - (-1) ** j * hw / math.pi / math.sqrt(1 + theta * theta) * math.cos(2 * math.pi * m / n) for m in range(n)]


def wkb_vs_exact(Vp: float = 10.0) -> dict:
    """Mathieu（縮放）位能 U(t) = −V' cos 2t（ħ = 1、2m = 1）的最低帶：
    精確帶寬 b₁ − a₀ vs WKB 4S，S = (ħω/2π)/√(1 + θ²)，θ = exp ∫|k|dt（在帶中心能量 E 的兩個轉折點之間，跨過 t = π/2 的位障）。"""
    a0, b1 = mathieu_edges_scipy(Vp, 1)[0]
    E = (a0 + b1) / 2
    tp = 0.5 * math.acos(-E / Vp)                       # U(tp) = E
    kap = lambda t: math.sqrt(max(-Vp * math.cos(2 * t) - E, 0.0))      # noqa: E731
    theta = math.exp(integrate.quad(kap, tp, math.pi - tp, limit=200)[0])
    # 古典週期：dt/dτ = 2√(E − U)（m = ½），週期 = 2∫dt/v
    vel = lambda t: 2 * math.sqrt(max(E + Vp * math.cos(2 * t), 1e-300))  # noqa: E731
    period = 2 * integrate.quad(lambda t: 1 / vel(t), -tp, tp, limit=200)[0]
    hw = 2 * math.pi / period
    S = hw / (2 * math.pi) / math.sqrt(1 + theta * theta)
    return {"exact": b1 - a0, "wkb": 4 * S, "ratio": 4 * S / (b1 - a0), "theta": theta, "hw": hw, "E": E}


# =============================================================================
# 習題
# =============================================================================
def problem_233(A=1.0, B=0.3, D=1.0) -> dict:
    """H = A/2(p₁² + q₁²) + D/2(p₂² + q₂²) + B(q₁q₂ + p₁p₂)（C = 0）：𝒜 = [[A² + B², B(A + D)], [B(A + D), D² + B²]]；
    C₂ 要 A = D；ω = A ± B = 1.3、0.7；拍週期 2π/0.6 = 10.4720 s、每拍 ω̄/Δω = 5/3 = 1.6667 次擺動；
    Poincaré 週期 20π = 62.8319 s（13 : 7）、含 6 拍。"""
    Aa = np.array([[A * A + B * B, B * (A + D)], [B * (A + D), D * D + B * B]])
    w = np.sort(np.sqrt(np.linalg.eigvalsh(Aa)))
    wl, wh = sp.nsimplify(A - B), sp.nsimplify(A + B)
    ratio = wh / wl                                   # 13/7
    p, q = sp.fraction(ratio)
    tp = 2 * sp.pi * q / wl                           # 7 個慢週期 = 13 個快週期
    return {"accel": Aa, "omega": w, "omega_exact": (wl, wh), "T_beat": 2 * sp.pi / (wh - wl),
            "wiggles": ((wh + wl) / 2) / (wh - wl), "T_poincare": sp.simplify(tp), "beats": int(p - q)}


def response_234(w, a=3.0, b=1.0, d=0.1, e=0.1, f=(1.0, 0.0)) -> np.ndarray:
    """習題 2.3.4 的穩態複振幅 x = (𝒂 − ω² + iω𝒅)⁻¹ f（x(t) = Re[x e^{iωt}]）。"""
    Am = np.array([[a, b], [b, a]], complex)
    Dm = np.array([[d, e], [e, d]], complex)
    return np.linalg.solve(Am - w * w * np.eye(2) + 1j * w * Dm, np.asarray(f, complex))


def response_234_modal(w, a=3.0, b=1.0, d=0.1, e=0.1, f=(1.0, 0.0)) -> np.ndarray:
    """同上，用 C₂ 模：x± = f±/(a ± b − ω² + iω(d ± e))，f± = (f₁ ± f₂)/√2。"""
    fp, fm = (f[0] + f[1]) / math.sqrt(2), (f[0] - f[1]) / math.sqrt(2)
    xp = fp / (a + b - w * w + 1j * w * (d + e))
    xm = fm / (a - b - w * w + 1j * w * (d - e)) if abs(fm) > 0 else 0.0
    return np.array([xp + xm, xp - xm]) / math.sqrt(2)


def problem_261() -> dict:
    """C₂：[[a, b], [c, d]] 與 R 對易 ⇔ a = d、b = c；C₃：與 r 對易 ⇔ 循環矩陣 [[A, B, C], [C, A, B], [B, C, A]]
    （E = I = A、F = G = B、D = H = C）；Hermitian ⇒ A 實、C = B*（B 可以是複數！）。"""
    a, b, c, d = sp.symbols("a b c d")
    H2 = sp.Matrix([[a, b], [c, d]])
    sol2 = sp.solve(list(commutator(H2, R2)), [c, d], dict=True)[0]
    s = sp.symbols("A B C D E F G H I_")
    H3 = sp.Matrix(3, 3, s)
    sol3 = sp.solve(list(commutator(H3, cyclic_rep(3))), s[3:], dict=True)[0]
    return {"C2": sol2, "C3": sol3, "C3_matrix": H3.subs(sol3)}


def problem_262(m1=10, m2=1, k1=20, k12=1) -> dict:
    """座標換尺度 y_i = √m_i x_i 後要求對角相等：(k₁ + k₁₂)/m₁ = (k₂ + k₁₂)/m₂ ⇒ k₂ = (k₁ + k₁₂)m₂/m₁ − k₁₂。
    k₁ = 20：k₂ = 1.1（ω² = 2.1 ∓ 1/√10 = 1.7838、2.4162）；k₁ = 8：k₂ = −0.1 < 0（不可能；需 k₁ ≥ k₁₂(m₁/m₂ − 1) = 9）。"""
    k2 = sp.nsimplify(R_(k1 + k12, 1) * m2 / m1 - k12)
    a11 = sp.nsimplify(R_(k1 + k12, m1))
    off = sp.nsimplify(k12 / sp.sqrt(m1 * m2))
    return {"k2": k2, "omega2": (a11 - off, a11 + off), "k1_min": k12 * (R_(m1, m2) - 1)}


def problem_271(n: int, a=1, b=0) -> list:
    """習題 2.7.1(e)：ω_k = √(b + 2a(1 − cos 2πk/n))，k = 0..n−1。n = 4：0、√2、√2、2；n = 5：0、1.1756、1.9021（各二重）；
    n = 8：0、0.7654、1.4142、1.8478（二重）、2。"""
    return [sp.sqrt(x) for x in ring_omega2(n, a, b)]


def problem_273(abar=R_(6, 5), alow=R_(4, 5)) -> dict:
    """C₈ → C₄（ā = 1.2、a̲ = 0.8、b = 0）：k = 0：ω² = 0、4；k = ±π/2：2 ∓ √2.08 = 0.5578、3.4422（二重）；k = π：2a̲ = 1.6、2ā = 2.4。"""
    out = {}
    for kk in range(-1, 3):
        out[kk] = broken_bands(2 * math.pi * kk / 4, abar, alow)
    num = np.sort(np.linalg.eigvalsh(broken_ring_matrix(8, abar, alow)))
    return {"bands": out, "numeric": num}


def problem_274(m: int = 4, a=1.0, b=0.0) -> dict:
    """兩牆之間 m 個擺 ＝ C_{2(m+1)} 環中「牆是節點」的那些正弦駐波：ω_j² = b + 4a sin²(jπ/(2(m+1)))，j = 1..m。
    m = 4（C₁₀）、a = 1、b = 0：ω = 0.6180、1.1756、1.6180、1.9021。"""
    A = np.diag([2 * a + b] * m) - a * (np.eye(m, k=1) + np.eye(m, k=-1))
    num = np.sqrt(np.sort(np.linalg.eigvalsh(A)))
    formula = np.array([math.sqrt(b + 4 * a * math.sin(j * math.pi / (2 * (m + 1))) ** 2) for j in range(1, m + 1)])
    return {"numeric": num, "formula": formula, "ring": 2 * (m + 1)}


def problem_275(a=-3.0, s=1.2, t=0.8) -> dict:
    """ẍ = 𝓜x，𝓜 = [[a, s, t], [t, a, s], [s, t, a]]（循環 ⇒ 有 C₃ 對稱，但沒有反射對稱）：
    λ_k = a + s e^{−ik} + t e^{ik}；k = 0：a + s + t（實）；k = ±2π/3：a − (s+t)/2 ∓ i(√3/2)(s − t)（s ≠ t 時是複數 ⇒ 一個模成長＝不穩定）。
    例 a = −3、s = 1.2、t = 0.8：λ₀ = −1（ω = 1，穩定）；λ = −4 ± 0.3464i ⇒ 成長率 0.0865。"""
    M = np.array([[a, s, t], [t, a, s], [s, t, a]])
    lam = np.linalg.eigvals(M)
    rates = [max(np.sqrt(complex(x)).real, (-np.sqrt(complex(x))).real) for x in lam]
    formula = [a + s * np.exp(-1j * k) + t * np.exp(1j * k) for k in (0, 2 * math.pi / 3, -2 * math.pi / 3)]
    return {"lambda": lam, "formula": np.array(formula), "growth": max(rates)}


def problem_281(j=1.0, k=2.0, alpha_deg=45.0, m=1.0) -> dict:
    """「Double Trouble」（模型假設：上下各一條繩子繞過定滑輪，兩根 j 串聯成 j/2）：
    V = (k/2)(x₁² + x₂²) + (j/2)[cos²α (x₁ + x₂)² + sin²α (y₁ − y₂)²]（x₁ 向右、x₂ 向左、y₁ 向上、y₂ 向下）。
    j = 1、k = 2、α = 45°：x 方向 ω² = k + 2j cos²α = 3（同向內縮）、k = 2；y 方向 2j sin²α = 1、0（繩子滑過滑輪，零頻）。
    對稱群 D₂ ≅ C₂ × C₂；(d) 加上 z 方向的鏡面 ⇒ C₂ × C₂ × C₂（D_{2h}，8 個 irrep）。"""
    c2, s2 = math.cos(math.radians(alpha_deg)) ** 2, math.sin(math.radians(alpha_deg)) ** 2
    K = np.zeros((4, 4))                  # 順序 x1, x2, y1, y2
    K[0, 0] = K[1, 1] = k + j * c2
    K[0, 1] = K[1, 0] = j * c2
    K[2, 2] = K[3, 3] = j * s2
    K[2, 3] = K[3, 2] = -j * s2
    w2 = np.sort(np.linalg.eigvalsh(K / m))
    return {"K": K, "omega2": w2}


def problem_282() -> dict:
    """n = 8：3（C8、C4×C2、C2×C2×C2）；9：2；12：2；16：5；24：3。"""
    return {n: abelian_types(n) for n in (8, 9, 12, 16, 24)}


def problem_283() -> dict:
    """12 質量環的方形駐波 s_p = −S(2π(p + ½)/12)（前 3 個在下、中間 6 個在上、後 3 個在下）：
    s_p = Σ_{k=1,3,5} c_k cos(k x_p)，c₁ = −(√6 + √2)/3 = −1.2879、c₃ = √2/3 = 0.4714、c₅ = −(√6 − √2)/3 = −0.3451
    （無窮方波：−4/π = −1.2732、4/3π = 0.4244、−4/5π = −0.2546）。"""
    x = [2 * sp.pi * (p + R_(1, 2)) / 12 for p in range(12)]
    s = [-1, -1, -1, 1, 1, 1, 1, 1, 1, -1, -1, -1]
    c = {k: sp.sqrtdenest(sp.nsimplify(sp.simplify(sum(sv * sp.cos(k * xv) for sv, xv in zip(s, x)) * R_(2, 12)))) for k in range(1, 7)}
    cont = {1: -4 / sp.pi, 3: 4 / (3 * sp.pi), 5: -4 / (5 * sp.pi)}
    return {"s": s, "c": c, "continuum": cont}


def problem_2121(Emax: float = 75.0) -> list:
    """Kronig–Penney（a = b = π/2）的偶然簡併：E = (2n_w)²、U = E − (2n_p)²（1 ≤ n_p < n_w）。"""
    out = []
    for nw in range(1, 10):
        E = (2 * nw) ** 2
        if E >= Emax:
            break
        for npl in range(1, nw):
            out.append((E - (2 * npl) ** 2, E))
    return sorted(out, key=lambda t: (t[1], t[0]))


def problem_2122(V: float = 25.0, w: float = 2.0) -> dict:
    return {"geometric": bound_geometric(V, w), "fd": bound_fd(V, w)}


def line_nodes(E, U, nw=4, w=KP_W, b=KP_B) -> int:
    """直線上 nw 個方井（外面無限長的 U）：從左邊衰減的解的節點數 ＝ 低於 E 的本徵值個數（Sturm 振盪定理）。"""
    q = math.sqrt(U - E)
    psi, dpsi = 1.0, q                            # 左邊 e^{qx}
    nodes = 0
    k = math.sqrt(E)
    for i in range(nw):
        phi = math.atan2(psi, dpsi / k)
        nodes += math.floor((phi + k * w) / math.pi) - math.floor(phi / math.pi)
        psi, dpsi = seg_matrix(E, 0.0, w) @ np.array([psi, dpsi])
        if i < nw - 1:
            p2, d2 = seg_matrix(E, U, b) @ np.array([psi, dpsi])
            nodes += int(np.sign(p2) != np.sign(psi))
            psi, dpsi = p2, d2
    A = psi + dpsi / q
    return nodes + int(np.sign(A) != np.sign(psi))


def line_levels(U: float = 60.0, nw: int = 4, count: int = 16) -> list:
    """二分法：第 i 個本徵值 = 使 line_nodes(E) 由 i 變 i+1 的 E（到 1e-12）。"""
    out = []
    for i in range(count):
        lo, hi = 1e-9, U - 1e-9
        if line_nodes(hi, U, nw) <= i:
            break
        for _ in range(60):
            mid = (lo + hi) / 2
            lo, hi = (mid, hi) if line_nodes(mid, U, nw) <= i else (lo, mid)
        out.append((lo + hi) / 2)
    return out


def problem_2123(U: float = 60.0, Emax: float = 80.0) -> dict:
    """(a) 圓上 4 個方井（C₄ 環，週期 π）：½tr = cos(2πm/4) = 1、0、−1；(b) 直線上 4 個方井（外面是無限長的 U）：節點計數＋二分。
    E < U 有 4 簇（各 4 個）：E ≈ 2.943、11.683、25.884、44.581；環是 1 + 2 + 1（H − 2S、H、H、H + 2S），
    直線 4 個都不簡併（H − 2S cos(jπ/5)，j = 1..4）；第 4 簇 S = 0.0118（環 44.557、44.581×2、44.604；直線 44.562、44.574、44.588、44.600）。"""
    Es = np.linspace(1e-6, min(U, Emax) - 1e-6, 8000)
    ring = []
    for target, mult in ((1.0, 1), (0.0, 2), (-1.0, 1)):
        f = np.array([kp_half_trace(E, U) - target for E in Es])
        for i in np.flatnonzero(np.sign(f[:-1]) * np.sign(f[1:]) < 0):
            ring += [optimize.brentq(lambda E: kp_half_trace(E, U) - target, Es[i], Es[i + 1], xtol=1e-13)] * mult
    return {"ring": sorted(ring), "line": line_levels(U)}


def ring4_fd(U: float = 60.0, N: int = 4000) -> np.ndarray:
    """問題 2.12.3(a) 的有限差分（週期邊界，Richardson）：最低 8 個本徵值。"""
    from scipy.sparse import diags
    from scipy.sparse.linalg import eigsh

    def ev(N_):
        L = 4 * math.pi
        h = L / N_
        x = (np.arange(N_) + 0.5) * h
        Ux = np.where(np.mod(x, math.pi) < KP_W, 0.0, U)
        Hm = diags([2 / h ** 2 + Ux, -np.ones(N_ - 1) / h ** 2, -np.ones(N_ - 1) / h ** 2, [-1 / h ** 2], [-1 / h ** 2]],
                   [0, 1, -1, N_ - 1, -(N_ - 1)], format="csc")
        return np.sort(eigsh(Hm, k=8, sigma=-1.0, which="LM")[0])
    e1, e2 = ev(N), ev(2 * N)
    return (4 * e2 - e1) / 3


def line4_fd(U: float = 60.0, pad: float = 2.5, N: int = 6000) -> np.ndarray:
    def ev(M):                                     # M 個格子 / 每個井
        h = KP_W / M
        npad = int(round(pad / h))
        x = (np.arange(-npad, 7 * M + npad) + 0.5) * h
        inside = (x > 0) & (x < 4 * math.pi - KP_B) & (np.mod(x, math.pi) < KP_W)
        Ux = np.where(inside, 0.0, U)
        n_ = len(x)
        return linalg.eigh_tridiagonal(2 / h ** 2 + Ux, -np.ones(n_ - 1) / h ** 2, select="i", select_range=(0, 7))[0]
    e1, e2 = ev(N // 8), ev(N // 4)
    return (4 * e2 - e1) / 3


# =============================================================================
# summary 與 check
# =============================================================================
def summary() -> dict:
    s = {}
    s["NOT_GROUP_ASSOC_FAILS"] = len(assoc_failures(table_from_rows(NOT_GROUP_NAMES, NOT_GROUP_BOOK)))
    s["group_counts_1_16"] = [GROUP_COUNTS[n] for n in range(1, 17)]
    s["abelian_1_16"] = [abelian_count(n) for n in range(1, 17)]
    s["beat_times2"] = [_f4(x) for x in beat_times2()]
    s["detuned_max"] = _f4(detuned_max_transfer())
    s["lissa_a"] = {k: _f4(a_for_ratio(v)) for k, v in LISSA.items()}
    s["closure"] = {k: (closure(v)["n_plus"], closure(v)["n_minus"], _f4(closure(v)["t_LJ"]), closure(v)["beats"],
                        _f4(closure(v)["t_beat"])) for k, v in LISSA.items()}
    s["ellipse"] = {k: [_f4(x) for x in ellipse_axes(a_for_ratio(v))] for k, v in LISSA.items()}
    s["phasors"] = [(_f4(abs(phasors(t)[0])), _f4(abs(phasors(t)[1]))) for t in PHASOR_TIMES]
    s["omega3"] = [_f4(sp.sqrt(x)) for x in omega2_3()]
    s["ring6"] = [_f4(sp.sqrt(x)) for x in ring_omega2(6, 1, 0)]
    s["alias"] = [_f4(x) for x in alias_values()]
    s["d2_omega"] = {k: _f4(sp.sqrt(v)) for k, v in d2_omega2().items()}
    s["pair5_exact"] = [_f4(x) for x in pair_exact(5)]
    s["pair5_pert"] = [_f4(x) for x in pair_pert(5)]
    b5 = np.linalg.eigvalsh(rotor_matrix(WEAK_V, 60)[1])
    ms_ = sorted(set(np.round(b5, 6)))
    s["rotor_weak_lowest"] = [f"{x:.4f}" for x in ms_[:16]]
    bands = rotor_bands(WEAK_V, 2)
    s["gap6"] = (f"{bands[6][0]:.4f}", f"{bands[6][1]:.4f}", f"{bands[6][1] - bands[6][0]:.4f}")
    s["k5"] = (f"{bands[5][0]:.4f}", f"{bands[5][1]:.4f}")
    tf = tight_fit()
    s["tight"] = (f"{tf['H']:.4f}", f"{tf['S']:.4f}", f"{tf['width']:.4f}", f"{tf['dev']:.2e}")
    s["mathieu0"] = [tuple(f"{x:.4f}" for x in e) for e in mathieu_edges(0.0, 3)]
    s["two_level"] = [f"{float(two_level(p)[1]):.4f}" for p in (0, 1, 2, 3)]
    s["nh3"] = {k: f"{v:.4e}" for k, v in nh3_numbers().items()}
    s["lz"] = [(f"{lz_formula(b):.4e}", f"{lz_numeric(b):.4e}") for b in LZ_BETAS]
    s["phase5"] = [_f4(x) for x in phase_shifts(5.0)]
    Sm = s_matrix(5.0)
    s["S5"] = (_f4(abs(Sm[0, 0]) ** 2), _f4(abs(Sm[1, 0]) ** 2))
    s["res"] = [_f4(x) for x in resonances()]
    s["bound"] = [(_f4(k), _f4(E)) for k, _d, E in bound_geometric()]
    s["wkb"] = {k: f"{v:.4f}" for k, v in wkb_vs_exact().items()}
    s["p233"] = {k: (str(v) if not isinstance(v, np.ndarray) else np.round(v, 4).tolist()) for k, v in problem_233().items()}
    s["p262"] = problem_262()
    s["p271"] = {n: [_f4(x) for x in problem_271(n)] for n in (4, 5, 8)}
    p3 = problem_273()
    s["p273"] = [f"{x:.4f}" for x in np.sqrt(np.clip(p3["numeric"], 0, None))]
    s["p274"] = [f"{x:.4f}" for x in problem_274()["formula"]]
    s["p275"] = f"{problem_275()['growth']:.4f}"
    s["p281"] = problem_281()["omega2"].round(6).tolist()
    s["p283"] = {k: _f4(v) for k, v in problem_283()["c"].items()}
    s["p2123"] = {k: [f"{x:.3f}" for x in v] for k, v in problem_2123().items()}
    return s


def check():
    """每個上卡片／講義的數字用兩種以上算法核對。"""
    f4 = _f4
    s2, s3 = sp.sqrt(2), sp.sqrt(3)
    e3 = eps(3)
    rng = np.random.default_rng(2)

    # ---- §2.1–2.2 群（H201）
    D2T = table_from_rows(D2_NAMES, D2_TABLE_BOOK)
    assert table_from_matrices(ROT3) == D2T == table_from_matrices(D2_REP4_BOOK)        # 轉動矩陣＝書上的表＝位置態的置換
    assert ROT3["Rz"] * ROT3["Ry"] == ROT3["Rx"] and D2T[("Rz", "Ry")] == "Rx"           # (2.1.1)
    ax = axioms(D2T)
    assert ax == {"closure": True, "associative": True, "identity": "1", "inverses": True, "latin": True} and is_abelian(D2T)
    assert order_counts(D2T) == {1: 1, 2: 3}
    C6 = table_from_rows(C6_NAMES, C6_BOOK)
    ren = dict(zip(C6_NAMES, cyclic_names(6)))
    assert {(ren[g], ren[h]): ren[v] for (g, h), v in C6.items()} == cyclic_table(6)
    assert all(axioms(C6)[k] for k in ("closure", "associative", "inverses", "latin")) and element_order(C6, "A") == 6
    NG = table_from_rows(NOT_GROUP_NAMES, NOT_GROUP_BOOK)
    axn = axioms(NG)
    assert axn == {"closure": True, "associative": False, "identity": "1", "inverses": True, "latin": True} and is_abelian(NG)
    assert NG[(NG[("a", "b")], "c")] == "a" and NG[("a", NG[("b", "c")])] == "c" and NG[("a", "b")] == "d" and NG[("b", "c")] == "f"
    assert ("a", "b", "c") in assoc_failures(NG) and len(assoc_failures(NG)) == 100
    for t in (C2_TABLE, C3_TABLE, cyclic_table(5)):
        assert all(axioms(t)[k] for k in ("closure", "associative", "inverses", "latin")) and is_abelian(t)
    assert [GROUP_COUNTS[n] for n in range(1, 17)] == [1, 1, 1, 2, 1, 2, 1, 5, 2, 2, 1, 5, 1, 2, 1, 14]
    assert (GROUP_COUNTS[24], GROUP_COUNTS[32], GROUP_COUNTS[48], GROUP_COUNTS[64]) == (15, 51, 52, 267) != (0, 0, 0, GROUP_COUNT_64_BOOK)
    assert all(abelian_count(n) == c for n, c in ABELIAN_COUNTS.items())                 # 分割數公式 ＝ 引用值
    assert all(abelian_count(n) == len(abelian_types(n)) <= GROUP_COUNTS.get(n, 10 ** 9) for n in range(1, 73))
    assert [abelian_count(n) for n in range(1, 17)] == [1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 2, 1, 1, 1, 5]
    assert sum(1 for n in range(1, 65) if GROUP_COUNTS[n] == abelian_count(n)) == 28      # 「所有群都可交換」的階（1, 2, 3, 4, 5, 7, 9, …, 61）

    # ---- §2.2–2.3 兩擺（H202–H203）
    a, b = PEND2["a"], PEND2["b"]
    A2 = accel2()
    assert A2 == sp.Matrix([[R_(221, 200), R_(-21, 200)], [R_(-21, 200), R_(221, 200)]])
    assert commutator(A2, R2) == sp.zeros(2) and R2.H * R2 == sp.eye(2) and R2 * A2 * R2.H == A2           # (2.2.6)、(2.2.8)、(2.2.9)
    assert commutator(sp.Matrix([[1, 2], [3, 4]]), R2) != sp.zeros(2)
    P = c2_projectors()
    Pc = c2_projectors_char()
    assert P == Pc and P[1] == sp.Matrix([[1, 1], [1, 1]]) / 2 and P[-1] == sp.Matrix([[1, -1], [-1, 1]]) / 2
    assert P[1] * P[-1] == sp.zeros(2) and P[1] ** 2 == P[1] and P[1] + P[-1] == sp.eye(2) and P[1] - P[-1] == R2   # (2.3.3)、(2.4.1)、(2.4.2)
    ep, em = sp.Matrix([1, 1]) / s2, sp.Matrix([1, -1]) / s2
    assert A2 * ep == b * ep and A2 * em == (2 * a + b) * em and (ep.T * A2 * em)[0] == 0              # (2.3.4)–(2.3.5)
    assert normal2() == (1, R_(121, 100)) and set(A2.eigenvals()) == {1, R_(121, 100)}
    assert sp.sqrt(normal2()[1]) == R_(11, 10)
    bt = beat_times2()
    assert bt == (10 * sp.pi, 20 * sp.pi) and [f4(x) for x in bt] == ["31.4159", "62.8319"]
    ts = np.linspace(0, 70, 701)
    md = Modal(A2, (1, 0))
    assert _close(md(ts), ode_solve(A2, (1, 0), (0, 0), ts), 1e-8)
    assert _close(md(float(bt[0])), [0, 1], 1e-9) and _close(md(float(bt[1])), [1, 0], 1e-9)
    x1 = md(ts)[:, 0]
    assert _close(x1, np.cos(0.05 * ts) * np.cos(1.05 * ts), 1e-12)                          # 拍：cos(Δω t/2)cos(ω̄t)
    assert _close(md(ts)[:, 1], np.sin(0.05 * ts) * np.sin(1.05 * ts), 1e-12)
    assert (f4(math.pi), f4(math.pi / 2)) == ("3.1416", "1.5708")
    assert cn_char_table(4)[1, :] == sp.Matrix([[1, I, -1, -I]])
    Ad = detuned_matrix()
    assert detuned_max_transfer() == 2 / sp.sqrt(5) and f4(detuned_max_transfer()) == "0.8944"
    wv, Vv = np.linalg.eigh(_real(Ad))
    th = 0.5 * math.atan(2 * float(DETUNED[0]) / float(DETUNED[2] - DETUNED[1]))
    tt = np.linspace(0, 4000, 400001)
    x2 = Modal(Ad, (1, 0))(tt)[:, 1]
    assert _close(x2, math.sin(2 * th) / 2 * (np.cos(np.sqrt(wv[0]) * tt) - np.cos(np.sqrt(wv[1]) * tt)), 1e-9)
    assert 0.89 < np.max(np.abs(x2)) <= 0.8945

    # ---- §2.3(a)(b) Lissajous、相量（H204）
    assert {k: a_for_ratio(v) for k, v in LISSA.items()} == {"weak": R_(41, 800), "mid": R_(69, 200), "strong": R_(3, 2)}
    want = {"weak": (20, 21, 40 * sp.pi, 1), "mid": (10, 13, 20 * sp.pi, 3), "strong": (1, 2, 2 * sp.pi, 1)}
    for key, rho in LISSA.items():
        c = closure(rho)
        assert (c["n_plus"], c["n_minus"], c["t_LJ"], c["beats"]) == want[key]
        assert sp.simplify(c["t_LJ"] - c["beats"] * c["t_beat"]) == 0                          # (2.3.13)
        tl = float(c["t_LJ"])
        tsl = np.linspace(0, tl, 2001)
        an = np.array(lissajous(tsl, rho)).T
        assert _close(an, lissajous_numeric(tsl, rho), 1e-8)
        mo = Modal(accel2(a_for_ratio(rho), 1), (1, 0))
        assert _close(mo(tl), [1, 0], 1e-9) and _close(mo.velocity(tl), [0, 0], 1e-9)          # 閉合
        assert np.min(np.hypot(an[100:-100, 0] - 1, an[100:-100, 1])) > 1e-3                 # 中途沒有回到起點
    assert [f4(closure(v)["t_LJ"]) for v in LISSA.values()] == ["125.6637", "62.8319", "6.2832"]
    assert f4(closure(LISSA["mid"])["t_beat"]) == "20.9440"
    ea = {k: ellipse_axes(a_for_ratio(v)) for k, v in LISSA.items()}
    assert ea["strong"] == (1, R_(1, 2)) and ea["mid"] == (1, R_(10, 13)) and ea["weak"] == (1, R_(20, 21))
    assert [f4(x) for x in (ea["mid"][1], ea["weak"][1])] == ["0.7692", "0.9524"]
    F = _real(quad_potential(a_for_ratio(2), 1))
    for th_ in np.linspace(0, 2 * np.pi, 7):
        u = np.array([math.cos(th_), math.sin(th_)]) @ np.array([[1, 1], [0.5, -0.5]]) / math.sqrt(2)
        assert abs(0.5 * u @ F @ u - 0.5) < 1e-12                                            # 半軸 1、½ 的橢圓上 V = ½
    zs = [phasors(t) for t in PHASOR_TIMES]
    assert [f4(abs(z1)) for z1, _ in zs] == ["1.0000", "0.9239", "0.7071", "0.3827", "0.0000"]
    assert [f4(abs(z2)) for _, z2 in zs] == ["0.0000", "0.3827", "0.7071", "0.9239", "1.0000"]
    for t in PHASOR_TIMES[1:4]:
        z1, z2 = phasors(t)
        assert abs(z1 - (-1j / math.tan(math.pi * float(t) / 2)) * z2) < 1e-12                # x₁ 領先 90°
    assert abs(zs[4][0]) < 1e-12 and abs(zs[4][1] + 1) < 1e-12
    mo = Modal(accel2(R_(3, 2) * sp.pi ** 2, sp.pi ** 2), (1, 0))                              # ω₊ = π、ω₋ = 2π
    for t in PHASOR_TIMES:
        z1, z2 = phasors(t)
        assert _close(mo(float(t)), [z1.real, z2.real], 1e-10)

    # ---- §2.4–2.5（H205）
    assert O_R3 ** 2 == sp.eye(3)
    OP = rep_projectors(O_REP3, C2_CHARS)
    assert OP["+"] == sp.Matrix([[R_(1, 2), 0, R_(1, 2)], [0, 1, 0], [R_(1, 2), 0, R_(1, 2)]])
    assert OP["-"] == sp.Matrix([[R_(1, 2), 0, R_(-1, 2)], [0, 0, 0], [R_(-1, 2), 0, R_(1, 2)]])
    assert c2_projectors(O_R3) == {1: OP["+"], -1: OP["-"]}                                    # 最小方程 ＝ 特徵標
    assert multiplicities(O_REP3, C2_CHARS) == {"+": 2, "-": 1} and OP["+"].trace() == 2
    assert sp.simplify(T3.H * T3) == sp.eye(3) and sp.simplify(T3.H * O_R3 * T3) == sp.diag(1, 1, -1)   # (2.4.8)
    A3 = accel3()
    assert A3 == sp.Matrix([[2, -1, 0], [R_(-1, 2), 2, R_(-1, 2)], [0, -1, 2]]) and commutator(A3, O_R3) == sp.zeros(3)
    red = reduce3()
    assert red == sp.Matrix([[2, -s2, 0], [-s2 / 2, 2, 0], [0, 0, 2]])
    aa, AA, bb = sp.symbols("a A b")
    Ag = sp.Matrix([[aa + bb, -aa, 0], [-AA, 2 * AA + bb, -AA], [0, -aa, aa + bb]])
    assert sp.simplify(T3.H * Ag * T3 - sp.Matrix([[aa + bb, -s2 * aa, 0], [-s2 * AA, 2 * AA + bb, 0], [0, 0, aa + bb]])) == sp.zeros(3)
    assert omega2_3() == [1, 2, 3] and set(red[:2, :2].eigenvals()) == {1, 3}
    for vv, lam in (((1, 1, 1), 1), ((1, -1, 1), 3), ((1, 0, -1), 2)):
        assert A3 * sp.Matrix(vv) == lam * sp.Matrix(vv)
    assert [f4(sp.sqrt(x)) for x in omega2_3()] == ["1.0000", "1.4142", "1.7321"]

    # ---- §2.6 C₃（H206）
    assert e3 == R_(-1, 2) + s3 * I / 2 and cyclic_rep(3) == sp.Matrix([[0, 0, 1], [1, 0, 0], [0, 1, 0]])
    Pch, Pmin = cn_projectors_char(3), cn_projectors_minimal(3)
    Psm = SM.projectors(cyclic_rep(3))
    for k in range(3):
        assert _zero(Pch[k] - Pmin[k]) and _zero(Pch[k] - Psm[sp.nsimplify(eps(3, k))])
    assert Pch[0] == sp.ones(3, 3) / 3
    assert _zero(Pch[1] - sp.Matrix([[1, e3, e3.conjugate()], [e3.conjugate(), 1, e3], [e3, e3.conjugate(), 1]]) / 3)   # (2.6.8)
    R3 = cyclic_rep(3)
    for k in range(3):
        assert _zero(Pch[k] * Pch[k] - Pch[k]) and _zero(R3 * Pch[k] - eps(3, k) * Pch[k])
        for j in range(3):
            if j != k:
                assert _zero(Pch[k] * Pch[j])
    assert _zero(sum(Pch.values(), sp.zeros(3)) - sp.eye(3))
    assert _zero(cn_char_table(3) - sp.Matrix([[1, 1, 1], [1, e3, e3.conjugate()], [1, e3.conjugate(), e3]]))   # (2.6.5)
    Tm = sp.Matrix.hstack(*[Pch[k][:, 0] * s3 for k in range(3)])
    assert _zero(Tm.H * Tm - sp.eye(3)) and _zero(Tm.H * R3 * Tm - sp.diag(1, e3, e3.conjugate()))       # (2.6.9)
    A33 = ring_matrix(3, 1, 1)
    assert A33 == sp.Matrix([[3, -1, -1], [-1, 3, -1], [-1, -1, 3]]) and _zero(Tm.H * A33 * Tm - sp.diag(1, 4, 4))
    assert ring_omega2(3, 1, 1) == [1, 4, 4] and _zero(Pch[1][:, 0] * s3 - C3_MOVING)
    x0b = np.ones(3) / math.sqrt(3) + math.sqrt(2) * _real(C3_STANDING[0]).ravel()           # |e₁⟩ + √2|e_Re⟩（p. 91 的拍）
    tb3 = np.linspace(0, 7, 15)
    assert _close(x0b, [math.sqrt(3), 0, 0], 1e-12)
    assert _close(Modal(A33, x0b)(tb3)[:, 0], (np.cos(tb3) + 2 * np.cos(2 * tb3)) / math.sqrt(3), 1e-10)
    assert _close(Modal(A33, x0b)(tb3)[:, 1], (np.cos(tb3) - np.cos(2 * tb3)) / math.sqrt(3), 1e-10)
    for v in C3_STANDING:
        assert A33 * v == 4 * v and sp.simplify((v.T * v)[0]) == 1
    e2n = to_np(C3_MOVING).ravel()
    for t in (0.3, 1.1):
        lhs = (np.exp(-2j * t) * e2n).real
        rhs = (math.cos(2 * t) * _real(C3_STANDING[0]).ravel() + math.sin(2 * t) * _real(C3_STANDING[1]).ravel()) / math.sqrt(2)
        assert _close(lhs, rhs, 1e-12)                                                       # Re[e^{−iωt}|e₂⟩]（p. 90）
    for n in (2, 4):
        pc, pm = cn_projectors_char(n), cn_projectors_minimal(n)
        assert all(_zero(pc[k] - pm[k]) for k in range(n))
    for n in (3, 5, 6, 12):
        pc, pm = cn_projectors_np(n, "char"), cn_projectors_np(n, "minimal")
        assert all(_close(pc[k], pm[k], 1e-9) for k in range(n)) and _close(sum(pc.values()), np.eye(n), 1e-12)
    assert all(_close(cn_projectors_np(3)[k], to_np(Pch[k]), 1e-12) for k in range(3))

    # ---- §2.7 C_n 環（H207）
    w6 = ring_omega2(6, 1, 0)
    assert w6 == [0, 1, 3, 4, 3, 1] and _close(sorted(np.linalg.eigvalsh(_real(ring_matrix(6)))), sorted(map(float, w6)), 1e-12)
    assert [f4(sp.sqrt(x)) for x in w6] == ["0.0000", "1.0000", "1.7321", "2.0000", "1.7321", "1.0000"]
    for n in (5, 8, 24):
        Am = _real(ring_matrix(n, 1, R_(1, 3)))
        w2n = [float(x) for x in ring_omega2(n, 1, R_(1, 3))]
        for k in range(n):
            v = wave_vector(n, k)
            assert _close(Am @ v, w2n[k] * v, 1e-12)
        assert _close(Am @ cos_wave(n, 2), w2n[2] * cos_wave(n, 2), 1e-12)
    kk = np.linspace(-np.pi, np.pi, 101)
    assert _close(dispersion(kk), 2 * np.abs(np.sin(kk / 2)), 1e-12) and abs(dispersion(0.1) - dispersion_series(0.1)) < 1e-7
    assert all(len(brillouin(n)) == n and sorted(alias(n, k) for k in range(n)) == brillouin(n) for n in range(1, 31))
    assert brillouin(24) == list(range(-11, 13)) and brillouin(5) == [-2, -1, 0, 1, 2]
    assert alias(5, 7) == alias(5, -3) == 2 and _close(alias_values(5, 7), alias_values(5, 2), 1e-12)
    assert _close(alias_values(5, -3), alias_values(5, 2), 1e-12)
    assert [f4(x) for x in alias_values()] == ["1.0000", "-0.8090", "0.3090", "0.3090", "-0.8090"]

    # ---- §2.7A–B 對稱破缺、奔馳波（H208）
    Mb = broken_ring_matrix()
    ev24 = np.sort(np.linalg.eigvalsh(Mb))
    ks = [2 * math.pi * k / 12 for k in brillouin(12)]
    form = np.sort(np.concatenate([broken_bands(k) for k in ks]))
    blk = np.sort(np.concatenate([np.linalg.eigvalsh(broken_block(k)) for k in ks]))
    assert _close(ev24, form, 1e-12) and _close(ev24, blk, 1e-12)
    assert _close(broken_bands(0.0), (0, 4), 1e-12) and _close(broken_bands(math.pi), (1, 3), 1e-12)
    assert abs(max(broken_bands(k)[0] for k in ks) - 1) < 1e-12 and abs(min(broken_bands(k)[1] for k in ks) - 3) < 1e-12
    assert _close(np.sort(np.linalg.eigvalsh(broken_ring_matrix(24, 1, 1))), sorted(map(float, ring_omega2(24, 1, 0))), 1e-12)
    assert float(ring_omega2(24, 1, 0)[6]) == 2.0
    soft, hard = broken_standing()
    assert _close(Mb @ soft, 1.0 * soft, 1e-12) and _close(Mb @ hard, 3.0 * hard, 1e-12) and abs(soft @ soft - 1) < 1e-12
    assert _close(broken_ring_matrix(24, R_(3, 2), R_(1, 2), R_(1, 5)), Mb + 0.2 * np.eye(24))
    assert _close(broken_bands(0.7, b=0.2), np.array(broken_bands(0.7)) + 0.2, 1e-12)        # (2.7.10) 加回 b
    assert swr(*GALLOP["a"]) == R_(1, 3) and swr(*GALLOP["b"]) == R_(1, 19)
    xs = np.linspace(0, np.pi, 2001)
    for key in ("a", "b"):
        I_, Ra = GALLOP[key]
        env = gallop_envelope(xs, I_, Ra)
        assert abs(env.max() - float(I_ + Ra)) < 1e-12 and abs(env.min() - float(I_ - Ra)) < 1e-6
        assert abs(env.min() / env.max() - float(swr(I_, Ra))) < 1e-6
        tt_ = np.linspace(0, 2 * np.pi, 37)
        assert all(_close(gallop_envelope(xs[::50], I_, Ra), np.abs(gallop_psi(xs[::50], t, I_, Ra)), 1e-12) for t in tt_)
    for d in (R_(1, 3), R_(1, 19), 0.5, -0.5):
        Ii, Rr = (1 + float(d)) / 2, (1 - float(d)) / 2
        for t in (0.1, 0.7, 1.3):
            zx = gallop_zeros(t, d)
            assert len(zx) == 2 and np.max(np.abs(gallop_psi(zx, t, Ii, Rr).real)) < 1e-12
            h = 1e-6
            x0 = zx[0]
            xp = gallop_zeros(t + h, d)
            x1_ = xp[np.argmin(np.abs(xp - x0))]
            assert abs((x1_ - x0) / h - zero_velocity(t, d)) < 1e-4 * max(1, abs(zero_velocity(t, d)))
        vs = [zero_velocity(t, d) for t in np.linspace(0, math.pi, 3601)]
        assert abs(min(vs, key=abs) - float(d)) < 1e-9 and abs(max(vs, key=abs) - 1 / float(d)) < 1e-6
    assert [f4(x) for x in (R_(1, 3), 3, R_(1, 19), 19)] == ["0.3333", "3.0000", "0.0526", "19.0000"]
    assert len(gallop_zeros(0.4, 0)) == 2 and len(gallop_zeros(0.4, 1)) == 2

    # ---- §2.8–2.11（H209）
    assert D2_REP4 == D2_REP4_BOOK and table_from_matrices(D2_REP4) == D2T
    Pp, Pq = d2_projectors_product(), d2_projectors_char()
    assert Pp == Pq and Pp["A1"] == sp.ones(4, 4) / 4
    for al, Pa in Pp.items():
        assert Pa * Pa == Pa and Pa.trace() == 1
        for g in D2_NAMES:
            assert D2_REP4[g] * Pa == D2_CHARS[al][g] * Pa                                  # 本徵算符（p. 107）
        for be, Pb in Pp.items():
            if be != al:
                assert Pa * Pb == sp.zeros(4)
    assert sum(Pp.values(), sp.zeros(4)) == sp.eye(4)
    Ad2 = d2_matrix()
    assert all(commutator(Ad2, D2_REP4[g]) == sp.zeros(4) for g in D2_NAMES)
    om = d2_omega2()
    assert om == {"A1": 1, "B2": R_(5, 2), "B1": R_(7, 2), "A2": 4} and set(Ad2.eigenvals()) == set(om.values())
    for al, Pa in Pp.items():
        v = Pa[:, 0] * 2
        assert Ad2 * v == om[al] * v
    assert [f4(sp.sqrt(v)) for v in om.values()] == ["1.0000", "1.5811", "1.8708", "2.0000"]
    C2T = sp.Matrix([[1, 1], [1, -1]])
    assert outer_char_table(C2T, C2T) == sp.Matrix([[D2_CHARS[al][g] for g in D2_NAMES] for al in D2_CHARS])   # (2.10)
    want8 = {"C8": {1: 1, 2: 1, 4: 2, 8: 4}, "C4xC2": {1: 1, 2: 3, 4: 4}, "C2xC2xC2": {1: 1, 2: 7}}
    for key, orders in ABELIAN8.items():
        els, tab = abelian_group(orders)
        assert order_counts(tab) == want8[key] and all(axioms(tab)[k] for k in ("associative", "inverses", "latin"))
        _e, _i, CT = abelian_char_table(orders)
        CTn = to_np(CT)
        assert _close(CTn @ CTn.conj().T, 8 * np.eye(8), 1e-12)                              # 列正交
    _e, _i, C42 = abelian_char_table((4, 2))
    assert C42 == outer_char_table(cn_char_table(4), cn_char_table(2))
    _e, _i, C222 = abelian_char_table((2, 2, 2))
    assert C222 == outer_char_table(outer_char_table(C2T, C2T), C2T)
    assert is_cyclic((3, 2)) and not is_cyclic((2, 2)) and is_cyclic((8,)) and not is_cyclic((4, 2)) and is_cyclic((9, 8))
    assert abelian_types(8) == [(8,), (4, 2), (2, 2, 2)] and len(abelian_types(72)) == 6
    assert POINT3["I"] * POINT3["sxy"] == POINT3["Rz"] and POINT3["sxz"] * POINT3["syz"] == POINT3["Rz"]   # (2.11.4)(2.11.5)
    for grp in (C2H, C2V, D2_NAMES):
        assert is_group(grp) and order_counts(point_table(grp)) == {1: 1, 2: 3} and is_abelian(point_table(grp))
    assert not is_group(C2V_BOOK) and POINT3["sxy"] * POINT3["syz"] == POINT3["Ry"]
    assert generated_order(POINT3["S4"]) == 4 and POINT3["S4"] ** 2 == POINT3["Rz"]

    # ---- §2.12A 轉子（H210）
    ms_, Hf = rotor_matrix(WEAK_V, 60)
    full = np.sort(np.linalg.eigvalsh(Hf))[:30]
    blocks = np.sort(np.concatenate([np.linalg.eigvalsh(rotor_block(k, WEAK_V, 6)[1]) for k in brillouin(12)]))[:30]
    assert _close(full, blocks, 1e-9)                                                          # 只耦合同一個 irrep
    assert pair_exact(5) == (37 - sp.sqrt(145), 37 + sp.sqrt(145)) and pair_pert(5) == (25 - R_(1, 24), 49 + R_(1, 24))
    assert [f4(x) for x in pair_exact(5)] == ["24.9584", "49.0416"] and [f4(x) for x in pair_pert(5)] == ["24.9583", "49.0417"]
    bw = rotor_bands(WEAK_V, 2)
    assert [f"{x:.4f}" for x in bw[5]] == ["24.9546", "49.0384"]
    assert abs(bw[5][0] - (float(pair_exact(5)[0]) - 1 / 264)) < 2e-4 and abs(bw[5][1] - (float(pair_exact(5)[1]) - 1 / 312)) < 2e-4
    assert [f"{x:.4f}" for x in bw[6]] == ["34.9965", "36.9965"] and f"{bw[6][1] - bw[6][0]:.4f}" == "2.0000"
    q = WEAK_V / 72
    assert abs(36 * special.mathieu_b(1, q) - bw[6][0]) < 1e-8 and abs(36 * special.mathieu_a(1, q) - bw[6][1]) < 1e-8
    assert abs(36 * special.mathieu_a(0, q) - bw[0][0]) < 1e-8
    tf = tight_fit()
    assert (f"{tf['H']:.4f}", f"{tf['S']:.4f}", f"{tf['width']:.4f}", f"{tf['dev']:.4f}") == ("-208.6223", "0.0897", "0.3588", "0.0004")
    assert all(len(v) == 3 for v in rotor_bands(0, 3).values()) and _close(rotor_bands(0, 1)[6], [36], 1e-12)
    me = mathieu_edges_scipy(STRONG_V / 36, 1)[0]
    assert abs(36 * (me[1] - me[0]) - tf["width"]) < 1e-8 and abs(tf["mean"] - tf["H"]) < 5e-4
    assert _close(rotor_edges(STRONG_V, 1)[0], (36 * me[0], 36 * me[1]), 1e-8)
    Hs, Ss = sp.symbols("H S")
    assert tb_levels(3) == [Hs - 2 * Ss, Hs + Ss, Hs + Ss] and tb_levels(4) == [Hs - 2 * Ss, Hs, Hs + 2 * Ss, Hs]
    for n in (3, 4, 12):
        assert sorted(tb_matrix(n, 0, 1).eigenvals(multiple=True), key=float) == sorted(tb_levels(n, 0, 1), key=float)
    for Vp in (0.5, 2.0, 5.0, 10.0, 20.0):
        assert _close(mathieu_edges(Vp, 3), mathieu_edges_scipy(Vp, 3), 1e-8)
    assert _close(mathieu_edges(0.0, 3), [(0, 1), (1, 4), (4, 9)], 1e-12)
    for m in (0, 12, 24, 36, 60, 13):
        assert abs(complex(square_well_G(m, 1)) - square_well_G_fft(m, 1.0)) < 2e-4
    assert square_well_G(12, 1) == -1 / sp.pi and square_well_G(36, 1) == 1 / (3 * sp.pi)
    for E, U in ((3.0, 20.0), (30.0, 20.0), (50.0, 12.0)):
        assert abs(kp_half_trace(E, U) - kp_half_trace_ode(E, U)) < 1e-8
    assert problem_2121() == list(KP_CROSSINGS)
    for U, E in KP_CROSSINGS:
        assert abs(abs(kp_half_trace(E, U)) - 1) < 1e-12 and kp_gap_near(U, E) < 2e-3
        assert kp_gap_near(U - 1, E) > 0.1 and kp_gap_near(U + 1, E) > 0.1                    # U 偏一點，能隙就打開
    assert abs(kp_half_trace(36, 20) + 1) < 1e-12
    ed = kp_edges(20.0, 60.0)
    assert len(ed) == 13 and all(abs(e - 36) > 5 for e in ed)                                 # 36 是切點（不是根）：帶隙閉合

    # ---- §2.12B 兩能階（H211）
    for pE in (0, 1, 2, 3, 0.37):
        lo, hi = two_level(pE)
        assert _close(np.linalg.eigvalsh(two_level_matrix(pE)), [float(lo), float(hi)], 1e-12)
    assert [f4(two_level(p)[1]) for p in (0, 1, 2, 3)] == ["1.0000", "1.4142", "2.2361", "3.1623"]
    assert abs(float(two_level(0.1)[1]) - (1 + 0.01 / 2)) < 2e-5 and abs(float(two_level(10)[1]) - 10) < 0.06
    tsn = np.linspace(0, 4, 9)
    assert _close(tunnel_populations(tsn)[1], tunnel_numeric(tsn), 1e-12)
    assert _close(tunnel_populations(math.pi)[0], 1, 1e-12) and _close(tunnel_populations(math.pi / 2)[1], 1, 1e-12)
    nh = nh3_numbers()
    assert (f"{nh['2S_eV']:.3e}", f"{nh['period_ps']:.2f}", f"{nh['wavenumber_cm']:.4f}") == ("9.872e-05", "41.89", "0.7962")
    assert abs(lz_numeric(4.0) - lz_formula(4.0)) < 2e-3 and f4(lz_formula(4.0)) == "0.4559"
    assert lz_numeric(0.25) < 1e-4 and f"{lz_formula(0.25):.1e}" == "3.5e-06"

    # ---- §2.12C–D 散射、穿越矩陣（H212）
    for e in (0.5, 5.0, 14.4784, 30.0):
        S = s_matrix(e)
        assert _close(S, s_matrix_direct(e), 1e-10) and _close(s_matrix(e, book=True), -S, 1e-12)
        assert _close(S.conj().T @ S, np.eye(2), 1e-12)
        dp, dm = phase_shifts(e)
        assert _close(S @ np.array([1, 1]), np.exp(2j * dp) * np.array([1, 1]), 1e-10)
        assert _close(S @ np.array([1, -1]), np.exp(2j * dm) * np.array([1, -1]), 1e-10)
        C = crossing_well(e)
        assert _close(C, crossing_well_steps(e), 1e-10) and abs(np.linalg.det(C) - 1) < 1e-12
        assert _close(s_from_c(C), S, 1e-10) and _close(c_from_s(S), C, 1e-10)
        assert abs(transmission(e) - abs(S[1, 0]) ** 2) < 1e-12 and abs(transmission(e) - 1 / abs(C[0, 0]) ** 2) < 1e-12
        assert abs(C[1, 0] / C[0, 1] - S[0, 0]) > 1e-3 or abs(S[0, 0]) < 1e-9                  # 書上的 C₂₁/C₁₂ 不對
    assert _close(s_matrix_direct(3.0, V=0.0), [[0, 1], [1, 0]], 1e-12) and _close(s_matrix(3.0, V=0.0, book=True), [[0, -1], [-1, 0]], 1e-12)
    assert [f4(x) for x in phase_shifts(5.0)] == ["-0.2919", "-1.0675"]
    S5 = s_matrix(5.0)
    assert (f4(abs(S5[0, 0]) ** 2), f4(abs(S5[1, 0]) ** 2)) == ("0.5098", "0.4902")
    res = resonances()
    assert [f4(x) for x in res] == ["14.4784", "36.6850", "63.8264"]
    for e in res:
        dp, dm = phase_shifts(e)
        assert abs(transmission(e) - 1) < 1e-12 and abs(math.cos(dp - dm)) < 1e-9                # Ramsauer–Townsend
    dp_, dm_, r_ = sp.symbols("delta_p delta_m r", real=True)
    Am_ = sp.Symbol("A_m")
    Ap_ = Am_ * sp.exp(I * (dp_ - dm_))
    ch1 = (Ap_ * sp.cos(r_ + dp_) + Am_ * sp.cos(r_ + dm_)) / s2
    ch2 = (Ap_ * sp.cos(r_ + dp_) - Am_ * sp.cos(r_ + dm_)) / s2
    for a1, a2 in ((0.7, -0.4), (1.1, 0.9)):                                                  # (2.12.48)：通道 2 只有 e^{ikr}
        vals = [complex((ch2 * sp.exp(-I * r_)).subs({r_: rv, dp_: a1, dm_: a2, Am_: 1})) for rv in (0.0, 0.8, 2.1)]
        assert max(abs(v - vals[0]) for v in vals) < 1e-12
    want49 = Am_ / s2 * (sp.cos(dp_ - dm_) * sp.exp(I * r_) * sp.exp(I * dp_) + sp.exp(-I * r_) * sp.exp(-I * dm_))
    for rv, a1, a2 in ((0.3, 0.7, -0.4), (1.1, -0.2, 0.9)):
        sub = {r_: rv, dp_: a1, dm_: a2, Am_: 1}
        assert abs(complex(ch1.subs(sub)) - complex(want49.subs(sub))) < 1e-12                  # (2.12.49) 的前因子是 A⁻/√2
    bg = bound_geometric()
    bf = bound_fd()
    assert len(bg) == len(bf) == 4 and _close([E for _k, _d, E in bg], bf, 2e-4)
    for kq, _d, _E in bg:
        kap = math.sqrt(25 - kq * kq)
        assert min(abs(kq * math.tan(kq) - kap), abs(-kq / math.tan(kq) - kap)) < 1e-9          # 偶／奇超越方程
    assert [(f4(k), f4(E)) for k, _d, E in bg] == [("1.3064", "1.7068"), ("2.5957", "6.7379"), ("3.8375", "14.7262"), ("4.9063", "24.0717")]
    for th_, al_ in ((3.0, 0.4), (50.0, 2.0)):
        Cc = c_barrier(th_) @ c_valley(al_)
        assert abs(np.trace(Cc) / 2 - math.sqrt(1 + th_ ** 2) * math.cos(al_)) < 1e-10           # (2.12.57)
        assert abs(Cc[1, 0] - (-1j * th_ * np.exp(1j * al_))) < 1e-10 and abs(np.linalg.det(Cc) - 1) < 1e-9   # (2.12.56) 的 (2,1)
    cl = wkb_cluster(10.0, 2.0, 30.0, 12, j=0)
    Sx = 2.0 / (2 * math.pi) / math.sqrt(1 + 900)
    assert _close(cl, [10.0 - 2 * Sx * math.cos(2 * math.pi * m / 12) for m in range(12)], 1e-12)  # (2.12.60) ＝ H − 2S cos k
    wk = wkb_vs_exact(10.0)
    assert abs(wk["exact"] - tf["width"] / 36) < 1e-9 and f4(wk["ratio"]) == "0.8958"
    assert (f4(wk["exact"]), f4(wk["wkb"]), f"{wk['theta']:.2f}") == ("0.0100", "0.0089", "601.85")

    # ---- 習題
    p233 = problem_233()
    q1, q2, p1, p2 = sp.symbols("q1 q2 p1 p2")
    Ah, Bh, Dh = sp.symbols("A B D")
    Hh = Ah / 2 * (p1 ** 2 + q1 ** 2) + Dh / 2 * (p2 ** 2 + q2 ** 2) + Bh * (q1 * q2 + p1 * p2)
    qd = [sp.diff(Hh, p1), sp.diff(Hh, p2)]
    pdot = {p1: -sp.diff(Hh, q1), p2: -sp.diff(Hh, q2)}
    qdd = [sp.expand(sum(sp.diff(x, pp) * pdot[pp] for pp in (p1, p2))) for x in qd]
    acc = sp.Matrix([[-sp.diff(x, q) for q in (q1, q2)] for x in qdd])
    assert sp.simplify(acc - sp.Matrix([[Ah ** 2 + Bh ** 2, Bh * (Ah + Dh)], [Bh * (Ah + Dh), Dh ** 2 + Bh ** 2]])) == sp.zeros(2)
    assert _close(p233["accel"], np.array(acc.subs({Ah: 1, Bh: 0.3, Dh: 1}).evalf(), dtype=float), 1e-12)
    assert _close(p233["omega"], [0.7, 1.3], 1e-12) and p233["omega_exact"] == (R_(7, 10), R_(13, 10))
    assert (f4(p233["T_beat"]), f4(p233["wiggles"]), f4(p233["T_poincare"]), p233["beats"]) == ("10.4720", "1.6667", "62.8319", 6)
    Hm = np.array([[1, 0.3, 0, 0], [0.3, 1, 0, 0], [0, 0, 1, 0.3], [0, 0, 0.3, 1]])        # (q1, q2, p1, p2) 的二次型
    J = np.block([[np.zeros((2, 2)), np.eye(2)], [-np.eye(2), np.zeros((2, 2))]])
    y20 = linalg.expm(J @ Hm * 20 * math.pi) @ np.array([1, 0, 0, 0])
    assert _close(y20, [1, 0, 0, 0], 1e-9)                                                    # Poincaré 重現
    for f_, ww in (((1, 1), 2.0), ((1, -1), 1.0), ((1, 0), 1.7), ((1, 0), 0.5)):
        assert _close(response_234(ww, f=f_), response_234_modal(ww, f=f_), 1e-12)
    assert abs(abs(response_234(2.0, f=(1, 1))[0]) - 2.5) < 1e-12 and abs(abs(response_234(1.0, f=(1, -1))[0]) - 1) < 1e-12
    assert set(np.linalg.eigvalsh(np.array([[3, 10], [10, 3]]))) == {-7.0, 13.0}              # 照印的 b = 10：一個模不穩定
    p261 = problem_261()
    a_, b_, c_, d_ = sp.symbols("a b c d")
    assert p261["C2"] == {c_: b_, d_: a_}
    H3c = p261["C3_matrix"]
    assert commutator(H3c, cyclic_rep(3)) == sp.zeros(3) and H3c[1, 1] == H3c[0, 0] and H3c[1, 2] == H3c[0, 1] and H3c[1, 0] == H3c[0, 2]
    Bc = sp.Rational(1, 2) + I
    Hc = sp.Matrix([[2, Bc, sp.conjugate(Bc)], [sp.conjugate(Bc), 2, Bc], [Bc, sp.conjugate(Bc), 2]])
    assert Hc == Hc.H and commutator(Hc, cyclic_rep(3)) == sp.zeros(3)
    assert len(set(np.round(np.linalg.eigvalsh(to_np(Hc)), 9))) == 3                         # 複數 B 拆開 ±k 的簡併
    p262 = problem_262()
    assert p262["k2"] == R_(11, 10) and problem_262(k1=8)["k2"] == R_(-1, 10) and p262["k1_min"] == 9
    Fm = np.array([[21.0, -1.0], [-1.0, 2.1]])
    assert _close(np.sort(np.linalg.eigvals(np.diag([0.1, 1.0]) @ Fm).real), [float(x) for x in p262["omega2"]], 1e-12)
    assert [f4(x) for x in p262["omega2"]] == ["1.7838", "2.4162"]
    assert [f4(x) for x in problem_271(4)] == ["0.0000", "1.4142", "2.0000", "1.4142"]
    assert [f4(x) for x in problem_271(5)] == ["0.0000", "1.1756", "1.9021", "1.9021", "1.1756"]
    assert [f4(x) for x in problem_271(8)] == ["0.0000", "0.7654", "1.4142", "1.8478", "2.0000", "1.8478", "1.4142", "0.7654"]
    p273 = problem_273()
    form273 = np.sort(np.concatenate([p273["bands"][k] for k in p273["bands"]]))
    assert _close(form273, p273["numeric"], 1e-12)
    assert [f"{x:.4f}" for x in np.sqrt(np.clip(p273["numeric"], 0, None))] == ["0.0000", "0.7468", "0.7468", "1.2649", "1.5492", "1.8553", "1.8553", "2.0000"]
    p274 = problem_274()
    assert _close(p274["numeric"], p274["formula"], 1e-12) and [f4(x) for x in p274["formula"]] == ["0.6180", "1.1756", "1.6180", "1.9021"]
    A10 = _real(ring_matrix(10, 1, 0))
    for j in range(1, 5):
        v = np.sin(np.pi * j * np.arange(10) / 5)
        assert _close(A10 @ v, p274["formula"][j - 1] ** 2 * v, 1e-12) and abs(v[0]) + abs(v[5]) < 1e-12   # 牆＝節點
    p275 = problem_275()
    assert _close(np.sort_complex(p275["lambda"]), np.sort_complex(p275["formula"]), 1e-12) and f4(p275["growth"]) == "0.0865"
    p281 = problem_281()
    assert _close(p281["omega2"], [0, 1, 2, 3], 1e-12)
    Rlr = np.array([[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, -1], [0, 0, -1, 0]])
    Rud = np.diag([1, 1, -1, -1])
    for Rm in (Rlr, Rud, Rlr @ Rud):
        assert _close(Rm @ p281["K"], p281["K"] @ Rm, 1e-12) and _close(Rm @ Rm, np.eye(4), 0)
    p282 = problem_282()
    assert {n: len(v) for n, v in p282.items()} == {8: 3, 9: 2, 12: 2, 16: 5, 24: 3}
    assert p282[16] == [(16,), (8, 2), (4, 4), (4, 2, 2), (2, 2, 2, 2)]
    p283 = problem_283()
    c283 = p283["c"]
    assert all(sp.simplify(c283[k] - v) == 0 for k, v in ((1, -(sp.sqrt(6) + s2) / 3), (3, s2 / 3), (5, -(sp.sqrt(6) - s2) / 3)))
    assert all(c283[k] == 0 for k in (2, 4, 6))
    assert [f4(c283[k]) for k in (1, 3, 5)] == ["-1.2879", "0.4714", "-0.3451"]
    assert [f4(v) for v in p283["continuum"].values()] == ["-1.2732", "0.4244", "-0.2546"]
    xp_ = [2 * sp.pi * (p + R_(1, 2)) / 12 for p in range(12)]
    assert all(sp.simplify(sum(c283[k] * sp.cos(k * xv) for k in (1, 3, 5)) - sv) == 0 for xv, sv in zip(xp_, p283["s"]))
    p2123 = problem_2123()
    ringl, linel = p2123["ring"], p2123["line"]
    assert len(ringl) == len(linel) == 16
    assert _close(ringl[:8], ring4_fd(), 2e-3) and _close(linel[:8], line4_fd(), 2e-3)
    Hc4, Sc4 = (ringl[12] + ringl[15]) / 2, (ringl[15] - ringl[12]) / 4
    assert _close(linel[12:], [Hc4 - 2 * Sc4 * math.cos(j * math.pi / 5) for j in range(1, 5)], 1e-3)  # 開鏈 H − 2S cos(jπ/5)
    assert abs(ringl[13] - ringl[14]) < 1e-9 and f"{Sc4:.4f}" == "0.0118"
    assert [f"{x:.3f}" for x in ringl[12:]] == ["44.557", "44.581", "44.581", "44.604"]
    assert [f"{x:.3f}" for x in linel[12:]] == ["44.562", "44.574", "44.588", "44.600"]
    assert [f"{np.mean(ringl[4 * c:4 * c + 4]):.3f}" for c in range(4)] == ["2.943", "11.683", "25.884", "44.581"]
    assert bg and rng is not None


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
    print(f"sym_abel（第 2 章）：D₂／C_n 群表與公理（含「不是群」的表）、群個數、C₂／C₃／C_n／D₂ 投影算符（最小方程＝特徵標）、"
          f"兩擺拍頻與 Lissajous（解析＝數值積分）、相量、三擺部分約化、單擺環色散與混疊、對稱破缺能帶、奔馳波 SWR、"
          f"8 階 Abel 群與點群、C₁₂ 轉子（2×2＝微擾＝全基底＝Mathieu）、緊束縛、Kronig–Penney、氨兩能階與 Landau–Zener、"
          f"方井 S 矩陣（公式＝直接解＝相移）、穿越矩陣、WKB、習題 2.3.3–2.12.3 全部核對通過（{time.time() - t0:.1f} 秒）。")
