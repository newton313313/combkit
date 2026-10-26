# 由 combkit/tools/sync_from_harter_sym.py 從 harter_sym 專案的 codes/utils/sym_oct.py 複製；請改 harter_sym 那份再同步。
r"""
sym_oct.py —— Harter《Principles of Symmetry, Dynamics, and Spectroscopy》（Wiley 1993）第 4 章
Theory and Applications of Higher Finite Symmetry and Induced Representations 的數學
（純 Python＋numpy＋sympy＋scipy，不 import manim；有限群工具全部借用同資料夾的 sym_nonabel（FiniteGroup、類代數、冪等元、
基本算子、相關表、外積、Hamilton 轉），不重寫）。書上的式號、筆誤與〔注〕見 ref/md/ch04_formulas.md（相對於 C:\mypython\manim\projects\harter_sym）。

約定
  * 群乘積 ab = 先 b 後 a；SN.FiniteGroup 的所有方法照用（O.mul、O.table、O.classes、O.class_mul、O.frequencies、O.correlation…）。
    精確值用 sympy；卡片上的小數一律 f"{float(x):.4f}"（f-string 不能直接 :.4f 一個 sympy 數）。運動方程 m ẍ = −F x（本徵值 ω²）；ħ = 1。
  * 元素名（O_NAMES，書末表 F.2.1 的順序；TeX 在 O.tex[名]）：1, r1..r4, r1^2..r4^2, R1^2, R2^2, R3^2, R1, R2, R3, R1^3, R2^3, R3^3, i1..i6。
    R_k ＝ 繞 x_k 軸逆時針 90°；r_k ＝ 逆時針 120°，軸 R_AXES = {1: (1,1,1), 2: (−1,−1,1), 3: (1,−1,−1), 4: (−1,1,−1)}（一個正四面體的四個頂點）；
    i_k ＝ 180°，軸 I_AXES = {1: (1,0,1), 2: (1,0,−1), 3: (1,1,0), 4: (1,−1,0), 5: (0,1,−1), 6: (0,1,1)}（i3、i4 在 xy 面，i1、i2 在 xz 面，
    i5、i6 在 yz 面）。由 (4.2.14)、(4.2.16)（全部 24 個）、(4.3.9)、(4.3.11) 與書上的乘積 r3R2 = R1、R1²i4 = R3、R1R2 = r1、r1i1 = R1³ 核對。
    O_h 的元素："I"、"I.r1"、"I.R1^2"…（TeX "I r_1"）。O_h 的 irrep 名 "A_{1g}"…"T_{2u}"（OH_IRREP_NAMES）。
  * C4 的約定照本章書上：|1_4⟩ = (−x̂ − iŷ)/√2，𝒟^{1_4}(R3) = −i（(4.2.35)、(4.2.43)、(4.3.32ₓ)；第 3 章的本系列約定是 +i，見 ch04_formulas〔注 3〕）。
    三方行進波仿 (4.2.36)：|1_3⟩ = (−v1 − iv2)/√2，𝒟^{1_3}(r1) = ε* = e^{−2πi/3}。相關表與這個號無關。
  * D2 = {1, R1², R2², R3²} 的標記照 (4.2.10)：B1 對 R1² 為 +1、B2 對 R2² 為 +1、A2 對 R3² 為 +1。

§4.1 群（H401–H403）
  O、OH（= O × C_i）、TD（R、i 換成 IR、Ii；O.iso_to(TD, O_TO_TD)）、TGRP（T，12 個：1、r、r²、R²）、TH（T × C_i）；O_TEX、O_CLASS_TEX = [1, c_r, c_{R²}, c_R, c_i]
  O_TURNS = {名: (軸, 角)}、o_turn(g) → (單位軸, 角, 半角)、o_quat(g)、element_of(3×3) → 名、
  nomogram_product(a, b) → SN.hamilton_product 的 dict ＋ name、sign（自旋 ½ 的號）：r3R2 = R1（+）、R1²i4 = R3（−，號依 i4 軸取向）、R1R2 = r1、
    r1i1 = R1³、r1R3 = i1（O_PRODUCTS_BOOK）；COSETS_439_BOOK（(4.3.9)）
  S4：DIAGONALS = {1: (1,−1,1), 2: (−1,1,1), 3: (1,1,−1), 4: (1,1,1)}、diag_perm(g) → {k: k'}、cycle_str(perm)（"[132]" = 1→3→2→1）、
    O_PERM_BOOK（(4.2.28) 的 24 個置換）：r1 = [132]、R1² = [13][24]、R3 = [1423]、i4 = [12]；類大小 O_CLASS_SIZES = [1, 8, 3, 6, 6]
  類代數 (4.1.6)：class_alg_rows()、O_CLASS_ALG_BOOK（(i, j) → 類係數；c_{R²}c_i = 2c_R + c_i、c_r² = 8 + 4c_r + 8c_{R²}、c_i² = 6 + 3c_r + 2c_{R²}）
  ci_powers() ＝ CI_POWERS：c_i² = 6 + 3c_r + 2c_{R²}、c_i³ = 16c_R + 20c_i、c_i⁴ = 120 + 108c_r + 104c_{R²}、c_i⁵ = 640c_R + 656c_i（書上 CI5_BOOK 錯）；
    O.min_poly_class(4) = c⁵ − 40c³ + 144c = (c+6)(c−6)(c+2)(c−2)c；CI_ROOTS = {6: A1, −6: A2, 2: T2, −2: T1, 0: E}
  o_class_idempotents() ＝ O_CLASS_IDEMP_BOOK（(4.1.10)）：ℙ^{T2} = (3 − c_{R²} − c_R + c_i)/8 …；o_central(α)（群代數元素）
  O_CHARS（(4.1.11)，列 O_IRREP_NAMES = A1 A2 E T1 T2，行 1 r R² R i）＝ o_char_table_class_algebra() ＝ 跡；O_IRREP_ALT（Γ1–Γ5）；O_DIMS；
    O_CLASS_EIG（°c_kχ^α_k/ℓ^α；c_i 行 6, −6, 0, −2, 2）
  OH_CHARS（(4.1.16)，10×10 ＝ C_i ⊗ O）、OH_CLASS_TEX；T_CHARS（(4.1.20)，列 A, ¹E, ²E, T；ε = EPS3）、T_IRREP_NAMES、TH_CHARS；GOLDEN = 1.6180
§4.2 不可約表示與子群鏈（H404–H406）
  irrep(α, g, basis) → Matrix：basis = "tetragonal"（O⊃D4⊃D2，𝒟^{T1} ＝ 3×3 轉動矩陣、𝒟^{T2} ＝ χ^{A2}·S𝒟^{T1}S（S_KOSHER = diag(1,−1,1)，(4.2.16)）、
    𝒟^E 用二次式基底 (2z²−x²−y², x²−y²)）、"trigonal"（O⊃D3⊃C2，T43 = 𝒯(4←3)）、"c4"（O⊃D4⊃C4，TMOV4 = (4.2.36) 的 𝒯）、"c3"（O⊃D3⊃C3）；
    α 可帶 g/u 給 O_h 的元素（𝒟^{αu}(Ig) = −𝒟^α(g)）。字典版：TETRA、TRIG、MOV4、MOV3（{α: {g: Matrix}}），BASES[basis]
  標記：TETRA_LABELS（T1：(E, B1), (E, B2), (A2, A2)；T2：…(B2, A2)；E：(A1, A1), (B1, A1)）、TRIG_LABELS、MOV4_LABELS（T1：1_4, 3_4, 0_4）、MOV3_LABELS
  書上的矩陣：T1_BOOK_REPS（(4.2.14)）、T2_BOOK（(4.2.16) 全部 24 個，signed_perm(code)）、E_BOOK_REPS（(4.2.19)）、E_COSET_CLASSES（24 個 𝒟^E 只有 6 種）、
    TRIG_T1_BOOK（(4.2.28) 類代表）、TRIG_T1_BOOK_TYPOS（r2²、r2 的 (2,1) 印 ±√3/3 應 ±√3/6；R2² 的 (2,3) 印 √2/2 應 √2/3）、
    MOV4_T1_BOOK（(4.2.35)）、D436_BOOK（照印 diag(i, −i, 1)，錯；𝒟(R3) = diag(−i, i, 1)）；TRIG_U_T1（投影算出 ＝ T43）、TRIG_U_T2
  子群鏈分裂：d4_idempotents()（"E1" = P^E p^{B1}、"E2"、"A1"、"A2"、"B1"、"B2"，(4.2.8)）、d3_idempotents()（"EA"、"EB"、"A2B"、"A1A"，(4.2.24)）、
    chain_split(α, "tetragonal"/"trigonal") → 非零的 ℙ^α·P（(4.2.9)、(4.2.18a)、(4.2.25)）、P410_BOOK（(4.2.10)，再除以 8）、
    p418(which, prefactor)（(4.2.18b,c)，正確 2/3，書印 3/4；單位元係數 1/12 對 3/32）、elementary_from_irrep(D, ℓ)（(ℓ/24)Σ𝒟_jk* g）
  子群（元素清單）：C4Z、D4Z、D2N（正規）、D2P（D2'）、D3T、C3T、C2T = {1, i4}、C2Z = {1, R3²}、C1、T_NAMES；O_h 內：C4V、C3V、C2V（(4.2.45)）、
    C2V_B（(4.3.65)）、D4H、D3D、D2H。SUB[key] = (群 "O"/"OH", 元素, 標記, 特徵標)；sub_chars(key)、sub_irreps(key)（矩陣）
  correlation(key) → 相關表（列 O 或 O_h 的 irrep，行 SUB[key][2]）＝ CORR_BOOK[key]：(4.2.42a) D4、(4.2.42b) C4、(4.2.44) D3、C3、
    (4.2.46a–c) C2v、C3v、C4v、(4.2.48a,b) D2'、D2；另有 D4h（(4.3.53)）、D3d（(4.3.86)）。conjugates(H)：D2' 有三個共軛、D2 正規
§4.3 誘導表示與六點模型（H407–H410）
  induced(G, sub, Dsub, leaders) → {g: Matrix}（(4.3.34)）、induced_char(G, sub, χ)（公式）、frobenius(key) → {a: (由矩陣, 由公式, 相關表的行)}（三者相同）
  SITES（|1⟩…|6⟩ = +z, −z, +x, −x, +y, −y）、SITE_LEADERS = [1, R1², r1, r2, r1², r2²]（(4.3.10)）、site_perm(g)、site_of(g)
  six_H(H, T, S)（(4.3.3)）、six_accel(k, t)（(4.3.2)：h = 2k + t、s = k/2）、six_induced(label)（0_4↑O：跡 (6, 0, 2, 2, 0) ⇒ A1 ⊕ E ⊕ T1；
    ℐ(i4) ＝ SIX_INDUCED_I4_BOOK；1_4↑O ＝ T1 ⊕ T2，ℐ(i4) ＝ SIX_INDUCED_1_4_I4_BOOK）、SIX_LEVELS（H+T+4S、H−T、H+T−2S）、
    SIX_VECTORS_BOOK ＝ six_projected(α, j)（(ℓ/24)Σ𝒟* g|1⟩）＝ six_coset_vector(α, j)（只對 6 個陪集代表）、six_coset_energy(α)（(4.3.46)；
    E：H + T − 2S，書 (4.3.46ₓ) 印 H − T + 2S）、six_numeric(H, T, S)（S = −1：−4, 0, 0, 0, 2, 2）
  d6_induced_examples()（(4.3.49)–(4.3.51)：0_2↑D6 = A1+E1+E2+B1、1_2↑D6 = A2+B2+E1+E2、0_6 = A1+A2、3_6 = B1+B2；正則 0_1↑O = (1,1,2,3,3)）
  外場（本系列取 SCAN = {H: 0, S: −1, T: 0, Qmax: 8, Dmax: 8, n: 41}；書上只說 S 為負）：
    V_Q：vq_matrix(Q, H, S)（(4.3.52)）、VQ_BASIS、VQ_LABELS、vq_reduced() ＝ VQ_REDUCED_BOOK（(4.3.56)）、VQ_SLOPES（A2u 1、A1g 1/3、Eg[A1g] 2/3、Eu 0、B1g 0）、
      vq_levels(Q) → [(E, 標記)]、vq_scan(Qs) → {"Q", "A1g-", "A1g+", "Eu", "A2u", "B1g"}（Fig. 4.3.5）、vq_avoided() = (2S, 4√2|S|)（S = −1：Q = −2、間距 5.6569）、
      VQ_HIGHQ、VQ_460_BOOK（(4.3.58)–(4.3.60)：極區 (1,1,0,0,0,0)/√2、赤道 (0,0,1,1,1,1)/2）
    V_D：vd_matrix(D, H, S)（(4.3.61)）、VD_BASIS（(4.3.63) 更正：T1u[A'] = (1,−1,0,0,0,0)/√2；照印 VD_T1U_BOOK）、VD_REDUCED（(4.3.66) 更正：H + 4S）、
      VD_SECOND_ORDER（D² 係數：A1g 1/12S、T1u 1/4S、Eg −1/3S；S = −1 ⇒ −1/12、−1/4、1/3）、vd_scan(Ds) → {"D", "A'1", "A'2", "A'3", "E", "B'"}（Fig. 4.3.7；
      D = 8：−8.6264、−1.7677、8.3941、0、2）、VD_469、VD_469_RESULT（(4.3.69)）
    磁場 C4（本系列模型：赤道環的穿隧加 Peierls 相位 β）：vb_matrix(β)、vb_levels(β) → {1_4: H+2S sin β, 3_4: H−2S sin β, 2_4: H−2S cos β, T_1(0_4): H−T,
      A_1(0_4), E(0_4)}（β = 0.3：−0.5910、0.5910、1.9107、0、−3.9407、2.0301；一階斜率 ±2S）
    (111) 電場：v111_matrix(d)（(4.3.74)）、V111_BASIS（(4.3.81)、(4.3.83)、(4.3.84)）、V111_A、V111_E（(4.3.80)、(4.3.82)）、
      v111_levels(d)（d = 2：A' −4.8284、0.8284；E −1.2361、3.2361 各二重）、T485_BOOK（(4.3.85)）、
      v111_tunnel() → ({A1g: H+2S+2S', A2u: H+2S−2S', Eg: H−S−S', Eu: H−S+S'}, 矩陣)（(4.3.87)–(4.3.88)）
§4.4 六氟化物 XF6（H411–H412）
  XF6_NAMES（21 個座標：A 軌道 A, R1²A, R1A, R2A, R1³A, R2³A；B 軌道 B, r1B, r2B, r1²B, r4²B, R2²B, R1B, R2B, R3B, R1³B, R2³B, i3B；C1, C2, C3，
    ＝ (4.4.22) 的欄序）、XF6_TEX、XF6_COORDS_LIST（(名, 原子 0–5 = SITES 上的 F、6 = X, 方向)）＝ xf6_coords()、A_LEADERS、B_LEADERS、F_SITES
  xf6_force(k, j, b)（幾何；彎曲 b：每個 90° 角 V = (b/4)Δθ²）＝ xf6_force_from_rows(k, j, b)（只用 (4.4.22) 三列＋對稱）、xf6_rows(k, j, b, book)、
    xf6_mass(m, M)、xf6_rep(g)（21×21 的 Γ(g)）、xf6_translation(axis)、xf6_isotypic(α)、xf6_classify()（A：A1g+Eg+T1u；B：T1g+T2g+T1u+T2u；C：T1u；
    全部 A1g+Eg+T1g+T2g+T2u+3T1u）、xf6_projected(α, j, k, orbit)（(4.4.13)–(4.4.17)）
  xf6_numeric(k, j, b, m, M)（21 個 ω²，6 個 0）＝ xf6_formulas(k, j, b, m, M)（A1g (j+4k)/m、Eg (j+k)/m、T2g 2(k+b)/m、T2u (k+b)/m、T1u±）、
    XF6_DEG、xf6_coset_omega2(k, j, b)（(4.4.23)–(4.4.29) 的陪集和；Q_AA = 2k+j、Q_AB = −√2k、Q_BB = k+b）
  xf6_t1u_matrix(k, j, b, m, M, book)（(4.4.33) 更正：右上 2(j−b)(M+6m)/3m）＝ xf6_t1u_reduced(...)（由 21 維 |c⟩、|c'⟩ 直接算）、
    xf6_secular(...) → {S, P, lp, lm}（(4.4.34)）、xf6_c_vectors(m, M) → (|c⟩, |c'⟩, |t⟩)
  XF6 = {UF6: m 19, M 240, b 0；SF6: M 32, b 0.2, k 0.18；SF6_34}、UF6_NU、SF6_NU、SF6_ISO_EXP
  xf6_scan(kjs, m, M, b) → {"kj", "A_{1g}", …}（ω 以 √(j/m) 為單位；UF6 k/j = 0：A1g = Eg = 1、T1u+ = 1.0763；SF6 k = 0：T1u± = 1.5802、0.6045）
  uf6_fit()（624/667 ⇒ k/j = 0.1958）、sf6_predict()（A1g = 775 定尺度：T1u+ 955、Eg 642、T2g 515、T2u 364、T1u− 525 cm⁻¹）
  SF6（k = 0.18j、b = 0.2j、m = 19、M = 32）：(4.4.38a) [[0.0670, 4.098], [0.000439, 0.1121]]、S = 0.1791、P = 0.0057、λ± = 0.1375j、0.0415j；
    xf6_isotope() = (1.0172, 1.0065)（書 1.0173、1.0064；實驗 965/948 = 1.0179、615/612 = 1.0049）；
    xf6_t1u_vectors() → {"+": eps 0.0682、eps_p 0.00117、center −0.1339、axial 0.1058、equatorial 0.0035；"-": 0.1135、−0.000706、0.0805、0.0909、−0.0794}；
    xf6_modes() → [(標記, ω, 21 維位移)]（A1g、Eg×2、T2g、T2u、T1u±、T1g 轉動、T1u 平移；給 SD.XF6Model.set_displacements）；ω(T2g)/ω(T2u) = √2 = 1.4142
§4.4C–D 量子振動（H413）
  ladder(nmax) → (a, a†)、ladder_commutator(nmax) = diag(1, …, 1, −nmax)、osc_levels(nmax) = n + ½、degeneracy(ℓ, n) = C(ℓ+n−1, n) ＝ degeneracy_count(ℓ, n)、
  PASCAL_464、harmonic_states(ℓ, n)（ℓ = 3、n = 2：6 組，含 (0,1,1)）、state_count(Emax)（SF6，Fig. 4.4.6 的頻率 SF6_NU_FIG：≤ 1000 cm⁻¹ 40 個、≤ 2000 965 個、
  ≤ 3000 11674 個）、turning_point(n, ℓ, ω)（(4.4.71)；n = 2、ℓ = 3：2.6458）
習題：problem_411（最少生成元 C3 1、C3v 2、D4 2、O 2、Oh 2）、problem_412（N(r1) = C3、N(R3) = C4）、problem_413（D2 四個（一個正規）、D3 四個、D4 三個，都不正規）、
  problem_421（(a)–(f) 相關表；PROBLEM_421_LABELS）、problem_423（籠子：CAGES、CAGE_GROUPS；p 軌域分裂的維度）、problem_431（0_3↑O 跡 8,2,0,0,0 ⇒ A1+A2+T1+T2；
  A1(D4)↑O 跡 3,0,3,1,1 ⇒ A1+E）、problem_432（E(C3v)↑Oh = Eg+T1g+T2g+Eu+T1u+T2u；Eg(D4h)↑Oh = T1g+T2g）、problem_433（立方八井）、
  problem_434（立方八面體 12 井）、problem_435（Y 群↓D5、C5、D3、C3；Y_CHARS、D5_CHARS）、problem_436(F, axis, kind)（立方八井加外場）、
  xy4_fit(mol)（習題 4.4.1，兩參數擬合）、problem_442（立方 X8 與立方烷的模分類）
check()：每個上卡片／講義的數字用兩種以上算法核對（約 11 秒）；summary() → {名稱: 值}。
直接執行：python utils/sym_oct.py（在 codes/ 裡）印一句成功訊息；--summary 印主要數字。
"""
import itertools
import math
from functools import lru_cache

import numpy as np
import sympy as sp
from scipy.optimize import brentq

try:
    from . import sym_nonabel as SN
except ImportError:                      # 直接執行 python utils/sym_oct.py
    import sym_nonabel as SN

R_ = sp.Rational
IU = sp.I
SQ2, SQ3, SQ6 = sp.sqrt(2), sp.sqrt(3), sp.sqrt(6)
EPS3 = SN.EPS3                           # ε = e^{2πi/3}
Hs, Ss, Ts, Qs, Ds, ds = sp.symbols("H S T Q D d")   # 六點模型的符號（H 對角、S 最近鄰、T 對面穿隧、Q／D／d 外場）


def _close(a, b, tol: float = 1e-9) -> bool:
    return np.allclose(np.asarray(a, dtype=complex), np.asarray(b, dtype=complex), atol=tol, rtol=0)


def _f4(x) -> str:
    return f"{float(x):.4f}"


def _num(M) -> np.ndarray:
    return np.array(sp.Matrix(M).evalf(), dtype=complex)


def _zero(M) -> bool:
    return all(sp.simplify(x) == 0 for x in sp.Matrix(M))


# ======================================================================== §4.1A 八面體群 O
O_NAMES = ["1", "r1", "r2", "r3", "r4", "r1^2", "r2^2", "r3^2", "r4^2", "R1^2", "R2^2", "R3^2",
           "R1", "R2", "R3", "R1^3", "R2^3", "R3^3", "i1", "i2", "i3", "i4", "i5", "i6"]      # 書末表 F.2.1 的順序
R_AXES = {1: (1, 1, 1), 2: (-1, -1, 1), 3: (1, -1, -1), 4: (-1, 1, -1)}   # r_k：繞這個方向逆時針 120°（一個正四面體的四個頂點）
X_AXES = {1: (1, 0, 0), 2: (0, 1, 0), 3: (0, 0, 1)}                       # R_k：繞 x_k 軸逆時針 90°
I_AXES = {1: (1, 0, 1), 2: (1, 0, -1), 3: (1, 1, 0), 4: (1, -1, 0), 5: (0, 1, -1), 6: (0, 1, 1)}   # i_k：180°


def _o_turns() -> dict:
    t = {"1": ((0, 0, 1), 0)}
    for k, a in R_AXES.items():
        t[f"r{k}"], t[f"r{k}^2"] = (a, 120), (a, 240)
    for k, a in X_AXES.items():
        t[f"R{k}"], t[f"R{k}^2"], t[f"R{k}^3"] = (a, 90), (a, 180), (a, 270)
    for k, a in I_AXES.items():
        t[f"i{k}"] = (a, 180)
    return {g: t[g] for g in O_NAMES}


O_TURNS = _o_turns()                     # {名: (軸, 逆時針角度)}；Hamilton 弧長 = 角度/2


def _rot_int(axis, deg) -> sp.Matrix:
    M = SN.quat_to_rot(SN.quat(axis, deg))
    Mi = np.rint(M).astype(int)
    assert np.allclose(M, Mi)
    return sp.Matrix(Mi.tolist())


def _o_tex(g: str) -> str:
    if g == "1":
        return r"\mathbf 1"
    k = g[1]
    return f"{g[0]}_{k}" + (g[2:] if len(g) > 2 else "")


O_TEX = {g: _o_tex(g) for g in O_NAMES}
O = SN.FiniteGroup({g: _rot_int(*O_TURNS[g]) for g in O_NAMES}, O_TEX, "O")
INV3 = -sp.eye(3)
OH = SN.direct_product(O, "I", INV3, "I", "Oh")          # (4.1.13)：O_h = O × C_i；名稱 "I"、"I.r1"…
O_CLASS_TEX = [r"\mathbf 1", "c_r", "c_{R^2}", "c_R", "c_i"]
O_CLASS_SIZES = [1, 8, 3, 6, 6]

# 由乘法圖（nomogram）讀出、書上明寫的乘積（ab = 先 b 後 a）
O_PRODUCTS_BOOK = [("r3", "R2", "R1"), ("R1^2", "i4", "R3"), ("R1", "R2", "r1"), ("r1", "i1", "R1^3"), ("r1", "R3", "i1")]
# (4.3.9)：C4 = {1, R3, R3², R3³} 的陪集
COSETS_439_BOOK = {"R1^2": ["R1^2", "i4", "R2^2", "i3"], "r1": ["r1", "i1", "r4", "R2"], "r2": ["r2", "i2", "r3", "R2^3"],
                   "r1^2": ["r1^2", "R1^3", "r3^2", "i6"], "r2^2": ["r2^2", "R1", "r4^2", "i5"]}

# S4：四條體對角線（Fig. 4.1.1(c) 的四個數字），(4.2.28) 與附錄表上的循環記號
DIAGONALS = {1: (1, -1, 1), 2: (-1, 1, 1), 3: (1, 1, -1), 4: (1, 1, 1)}
O_PERM_BOOK = {"1": "", "i4": "[12]", "R1^2": "[13][24]", "R3": "[1423]", "r1": "[132]", "i5": "[13]", "r4": "[234]",
               "i6": "[24]", "r1^2": "[123]", "i2": "[23]", "r2^2": "[142]", "R2^3": "[1342]", "R2^2": "[14][23]",
               "R3^3": "[1324]", "R3^2": "[12][34]", "i3": "[34]", "r2": "[124]", "R1": "[1234]", "r3": "[143]",
               "R1^3": "[1432]", "r3^2": "[134]", "i1": "[14]", "r4^2": "[243]", "R2": "[1243]"}


def diag_perm(g: str) -> dict:
    """g 對四條體對角線的置換 {k: k'}（g·d_k = ±d_k'）。"""
    M = O.num[g].real
    out = {}
    for k, d in DIAGONALS.items():
        v = M @ np.array(d, float)
        out[k] = next(kk for kk, dd in DIAGONALS.items() if np.allclose(v, dd) or np.allclose(v, -np.array(dd)))
    return out


def cycle_str(perm: dict) -> str:
    """{1: 3, 3: 2, 2: 1, 4: 4} → "[132]"（每個循環從最小的數字寫起，長度 1 的省略）。"""
    seen, out = set(), ""
    for k in sorted(perm):
        if k in seen or perm[k] == k:
            seen.add(k)
            continue
        cyc, x = [], k
        while x not in seen:
            seen.add(x)
            cyc.append(x)
            x = perm[x]
        out += "[" + "".join(map(str, cyc)) + "]"
    return out


def o_turn(g: str) -> tuple:
    """(單位軸 np, 角度, Hamilton 弧長 = 角度/2)。"""
    a, deg = O_TURNS[g]
    a = np.asarray(a, float)
    return a / np.linalg.norm(a), deg, deg / 2


def o_quat(g: str) -> np.ndarray:
    return SN.quat(*O_TURNS[g])


def element_of(Rm) -> str:
    """3×3 轉動（或 O_h 的矩陣）→ 名稱（O 或 OH）。"""
    G = OH
    return next(g for g in G.names if _close(G.num[g], np.asarray(Rm, complex)))


def nomogram_product(a: str, b: str) -> dict:
    """ab（先 b 後 a）的 Hamilton 弧相接（Fig. 4.1.4）：SN.hamilton_product 的結果＋ 'name'（乘積的名稱）＋ 'sign'（四元數的號，自旋 ½）。"""
    ax_b, deg_b = O_TURNS[b]
    ax_a, deg_a = O_TURNS[a]
    hp = SN.hamilton_product(ax_b, deg_b, ax_a, deg_a)
    name = element_of(SN.quat_to_rot(hp["q"]))
    qn = o_quat(name)
    hp["name"] = name
    hp["sign"] = 1 if np.allclose(hp["q"], qn) else -1
    assert np.allclose(hp["q"], hp["sign"] * qn)
    return hp


# ---------------------------------------------------------------- 類代數 (4.1.6)–(4.1.10)
O_CLASS_ALG_BOOK = {(1, 1): [8, 4, 8, 0, 0], (1, 2): [0, 3, 0, 0, 0], (1, 3): [0, 0, 0, 4, 4], (1, 4): [0, 0, 0, 4, 4],
                    (2, 2): [3, 0, 2, 0, 0], (2, 3): [0, 0, 0, 1, 2], (2, 4): [0, 0, 0, 2, 1],
                    (3, 3): [6, 3, 2, 0, 0], (3, 4): [0, 3, 4, 0, 0], (4, 4): [6, 3, 2, 0, 0]}   # 類 1, c_r, c_R², c_R, c_i
CI_POWERS = [[1, 0, 0, 0, 0], [0, 0, 0, 0, 1], [6, 3, 2, 0, 0], [0, 0, 0, 16, 20], [120, 108, 104, 0, 0], [0, 0, 0, 640, 656]]
CI5_BOOK = [2400, 2160, 2080, 320, 256]          # 書上 (4.1.7) 第二行（把 20c_i^4 代進去的錯誤結果）
CI_ROOTS = {6: "A_1", -6: "A_2", 0: "E", -2: "T_1", 2: "T_2"}      # (4.1.8)：c_i 的根 ↔ irrep
O_IRREP_NAMES = ["A_1", "A_2", "E", "T_1", "T_2"]
O_IRREP_ALT = ["\\Gamma_1", "\\Gamma_2", "\\Gamma_3", "\\Gamma_4", "\\Gamma_5"]
O_CHARS_BOOK = [[1, 1, 1, 1, 1], [1, 1, 1, -1, -1], [2, -1, 2, 0, 0], [3, 0, -1, 1, -1], [3, 0, -1, -1, 1]]   # (4.1.11)
O_CHARS = sp.Matrix(O_CHARS_BOOK)
O_DIMS = [1, 1, 2, 3, 3]
O_CLASS_IDEMP_BOOK = {"A_1": [R_(1, 24)] * 5, "A_2": [R_(1, 24)] * 3 + [-R_(1, 24)] * 2,
                      "E": [R_(1, 6), -R_(1, 12), R_(1, 6), 0, 0], "T_1": [R_(3, 8), 0, -R_(1, 8), R_(1, 8), -R_(1, 8)],
                      "T_2": [R_(3, 8), 0, -R_(1, 8), -R_(1, 8), R_(1, 8)]}       # (4.1.10)（類係數）
O_CLASS_EIG = [[1, 8, 3, 6, 6], [1, 8, 3, -6, -6], [1, -4, 3, 0, 0], [1, 0, -1, 2, -2], [1, 0, -1, -2, 2]]  # °c_kχ^α_k/ℓ^α


def class_alg_rows() -> list:
    """類代數表（列、行依 O_CLASS_TEX）：rows[i][j] = c_i c_j 的類係數向量。"""
    return [[list(O.struct[i, j]) for j in range(5)] for i in range(5)]


def ci_powers(n: int = 5) -> list:
    """c_i^0 … c_i^n 的類係數向量（(4.1.7) 的正確版）。"""
    ci = sp.Matrix([0, 0, 0, 0, 1])
    out = [sp.Matrix([1, 0, 0, 0, 0])]
    for _ in range(n):
        out.append(O.class_mul(out[-1], ci))
    return [list(v) for v in out]


def o_class_idempotents() -> dict:
    """{irrep: ℙ^α 的類係數向量}（由 c_i 的最小方程 (4.1.8)–(4.1.9) 算，依 CI_ROOTS 命名）。"""
    return {CI_ROOTS[int(lam)]: list(v) for lam, v in O.class_idempotents_of(4).items()}


def o_char_table_class_algebra() -> sp.Matrix:
    """只用類代數（§3.5A）算的特徵標表，列依書 (A1, A2, E, T1, T2)。"""
    return SN.book_order(O, O_CHARS_BOOK)[1]


# ---------------------------------------------------------------- O_h、T_d、T、T_h (4.1.13)–(4.1.20)
OH_IRREP_NAMES = ["A_{1g}", "A_{2g}", "E_g", "T_{1g}", "T_{2g}", "A_{1u}", "A_{2u}", "E_u", "T_{1u}", "T_{2u}"]
OH_CHARS = SN.outer_table(SN.C2_CHARS, O_CHARS)          # (4.1.16)：行 1, r, R², R, i, I, Ir, IR², IR, Ii
OH_CLASS_TEX = O_CLASS_TEX + ["I", "Ic_r", "Ic_{R^2}", "Ic_R", "Ic_i"]

TD_NAMES = O_NAMES[:12] + ["I." + g for g in O_NAMES[12:]]       # (§4.1C)：R_k、i_k 換成 IR_k、Ii_k
TD = SN.FiniteGroup({g: OH.mats[g] for g in TD_NAMES}, {g: OH.tex[g] for g in TD_NAMES}, "Td")
O_TO_TD = {g: (g if g in TD_NAMES else "I." + g) for g in O_NAMES}
T_NAMES = O_NAMES[:12]                                            # T = {1, r_k, r_k², R_k²}
TGRP = SN.FiniteGroup({g: O.mats[g] for g in T_NAMES}, {g: O.tex[g] for g in T_NAMES}, "T")
TH = SN.direct_product(TGRP, "I", INV3, "I", "Th")                # (4.1.19)
T_IRREP_NAMES = ["A", r"{}^1E", r"{}^2E", "T"]                    # 書上第 2、3 列都印 χ^ε（第 3 列應為共軛 ε*）
T_CHARS = sp.Matrix([[1, 1, 1, 1], [1, EPS3, sp.conjugate(EPS3), 1], [1, sp.conjugate(EPS3), EPS3, 1], [3, 0, 0, -1]])  # (4.1.20)
TH_IRREP_NAMES = [n + p for p in ("_g", "_u") for n in ("A", r"{}^1E", r"{}^2E", "T")]
TH_CHARS = SN.outer_table(SN.C2_CHARS, T_CHARS)
GOLDEN = (1 + sp.sqrt(5)) / 2                                     # Fig. 4.1.6：1.6180…


def oh_split(g: str) -> tuple:
    """O_h 的元素名 → (轉動部分名稱, 是否帶反演)。"""
    if g == "I":
        return "1", True
    if g.startswith("I."):
        return g[2:], True
    return g, False


# ======================================================================== §4.2 不可約表示與子群鏈
S_KOSHER = sp.diag(1, -1, 1)
E_BASIS_Q = [sp.diag(-1, -1, 2) / SQ6, sp.diag(1, -1, 0) / SQ2]   # E 的二次式基底：|E,1⟩ ~ 2z²−x²−y²（D4 的 A1）、|E,2⟩ ~ x²−y²（B1）


def chi_A2(g: str) -> int:
    return 1 if O.class_of(g) in (0, 1, 2) else -1


def _irrep_T1(g):
    return O.mats[g]


def _irrep_T2(g):
    return chi_A2(g) * S_KOSHER * O.mats[g] * S_KOSHER        # 「正統」（kosher）𝒟^{T2}（(4.2.16)；(4.2.15) 再把第 2 基底變號）


def _irrep_E(g):
    M = O.mats[g]
    return sp.Matrix(2, 2, lambda i, j: sp.nsimplify(sp.simplify((E_BASIS_Q[i] * M * E_BASIS_Q[j] * M.T).trace())))


def _irrep_1d(row):
    return lambda g: sp.Matrix([[O_CHARS[row, O.class_of(g)]]])


TETRA = {"A_1": {g: sp.Matrix([[1]]) for g in O_NAMES}, "A_2": {g: sp.Matrix([[chi_A2(g)]]) for g in O_NAMES},
         "E": {g: _irrep_E(g) for g in O_NAMES}, "T_1": {g: _irrep_T1(g) for g in O_NAMES},
         "T_2": {g: _irrep_T2(g) for g in O_NAMES}}                 # O ⊃ D4 ⊃ D2 駐波（四方）基底
TETRA_LABELS = {"A_1": [("A_1", "A_1")], "A_2": [("B_1", "A_1")], "E": [("A_1", "A_1"), ("B_1", "A_1")],
                "T_1": [("E", "B_1"), ("E", "B_2"), ("A_2", "A_2")], "T_2": [("E", "B_1"), ("E", "B_2"), ("B_2", "A_2")]}

# (4.2.16) 照抄（從頁面影像 p017 讀出全部 24 個）：帶號置換，第 i 列的非零元在 |c|-1 行、號 = sign(c)
T2_BOOK = {"1": (1, 2, 3), "R1^2": (1, -2, -3), "r1": (3, -1, -2), "r2": (-3, -1, 2), "r1^2": (-2, -3, 1), "r2^2": (-2, 3, -1),
           "R3^2": (-1, -2, 3), "R2^2": (-1, 2, -3), "r4": (3, 1, 2), "r3": (-3, 1, -2), "r3^2": (2, -3, -1), "r4^2": (2, 3, 1),
           "R3": (-2, 1, -3), "i4": (-2, -1, 3), "i1": (-3, 2, -1), "i2": (3, 2, 1), "R1^3": (-1, 3, -2), "R1": (-1, -3, 2),
           "R3^3": (2, -1, -3), "i3": (2, 1, 3), "R2": (-3, -2, 1), "R2^3": (3, -2, -1), "i6": (1, 3, 2), "i5": (1, -3, -2)}


def signed_perm(code) -> sp.Matrix:
    M = sp.zeros(3)
    for i, c in enumerate(code):
        M[i, abs(c) - 1] = 1 if c > 0 else -1
    return M


T1_BOOK_REPS = {"r1": (3, 1, 2), "R1^2": (1, -2, -3), "R1": (1, -3, 2), "R3": (-2, 1, 3), "i1": (3, -2, 1)}   # (4.2.14)
HALF, RT3_2 = R_(1, 2), SQ3 / 2
E_BOOK_REPS = {"1": [[1, 0], [0, 1]], "R1^2": [[1, 0], [0, 1]], "r1": [[-HALF, -RT3_2], [RT3_2, -HALF]],
               "R3": [[1, 0], [0, -1]], "R1": [[-HALF, -RT3_2], [-RT3_2, HALF]], "i1": [[-HALF, RT3_2], [RT3_2, HALF]]}  # (4.2.19)
E_COSET_CLASSES = [["1", "R1^2", "R2^2", "R3^2"], ["r1", "r2", "r3", "r4"], ["r1^2", "r2^2", "r3^2", "r4^2"],
                   ["R3", "R3^3", "i3", "i4"], ["i1", "i2", "R2", "R2^3"], ["R1", "R1^3", "i5", "i6"]]   # 24 個 𝒟^E 只有 6 種

# ---------------------------------------------------------------- 子群（O 與 O_h 內）與它們的特徵標（行依元素順序）
C4Z = ["1", "R3", "R3^2", "R3^3"]
D4Z = ["1", "R3^2", "R3", "R3^3", "R1^2", "R2^2", "i3", "i4"]
D2N = ["1", "R1^2", "R2^2", "R3^2"]                       # 正規子群 D2
D2P = ["1", "R3^2", "i3", "i4"]                           # D2'（不正規）
D3T = ["1", "r1", "r1^2", "i2", "i4", "i5"]
C3T = ["1", "r1", "r1^2"]
C2T = ["1", "i4"]
C2Z = ["1", "R3^2"]
C1 = ["1"]
C4V = ["1", "R3^2", "R3", "R3^3", "I.R1^2", "I.R2^2", "I.i3", "I.i4"]
C3V = ["1", "r1", "r1^2", "I.i2", "I.i4", "I.i5"]
C2V = ["1", "i4", "I.R3^2", "I.i3"]                       # (4.2.45)
C2V_B = ["1", "R3^2", "I.R2^2", "I.R1^2"]                 # (4.3.65)
D4H = D4Z + ["I" if g == "1" else "I." + g for g in D4Z]
D3D = D3T + ["I" if g == "1" else "I." + g for g in D3T]
D2H = D2N + ["I" if g == "1" else "I." + g for g in D2N]
_m = -IU
SUB = {   # 名稱: (群 O 或 OH, 元素, irrep 標記, 特徵標（列 irrep、行元素）)
    "C4": ("O", C4Z, ["0_4", "1_4", "2_4", "3_4"], [[1, 1, 1, 1], [1, -IU, -1, IU], [1, -1, 1, -1], [1, IU, -1, -IU]]),   # (4.2.43)
    "D4": ("O", D4Z, ["A_1", "A_2", "B_1", "B_2", "E"],
           [[1] * 8, [1, 1, 1, 1, -1, -1, -1, -1], [1, 1, -1, -1, 1, 1, -1, -1], [1, 1, -1, -1, -1, -1, 1, 1], [2, -2, 0, 0, 0, 0, 0, 0]]),
    "D2": ("O", D2N, ["A_1", "B_1", "B_2", "A_2"], [[1, 1, 1, 1], [1, 1, -1, -1], [1, -1, 1, -1], [1, -1, -1, 1]]),
    "D2'": ("O", D2P, ["A_1", "B_1", "A_2", "B_2"], [[1, 1, 1, 1], [1, -1, 1, -1], [1, 1, -1, -1], [1, -1, -1, 1]]),   # (4.2.47)
    "D3": ("O", D3T, ["A_1", "A_2", "E"], [[1] * 6, [1, 1, 1, -1, -1, -1], [2, -1, -1, 0, 0, 0]]),
    "C3": ("O", C3T, ["0_3", "1_3", "2_3"], [[1, 1, 1], [1, sp.conjugate(EPS3), EPS3], [1, EPS3, sp.conjugate(EPS3)]]),
    "C2": ("O", C2T, ["A", "B"], [[1, 1], [1, -1]]),
    "C2z": ("O", C2Z, ["0_2", "1_2"], [[1, 1], [1, -1]]),
    "C1": ("O", C1, ["0_1"], [[1]]),
    "T": ("O", T_NAMES, T_IRREP_NAMES, None),
    "C2v": ("OH", C2V, ["A'", "B'", "A''", "B''"], [[1, 1, 1, 1], [1, -1, 1, -1], [1, 1, -1, -1], [1, -1, -1, 1]]),
    "C3v": ("OH", C3V, ["A'", "A''", "E"], [[1] * 6, [1, 1, 1, -1, -1, -1], [2, -1, -1, 0, 0, 0]]),
    "C4v": ("OH", C4V, ["A'", "B'", "A''", "B''", "E"],
            [[1] * 8, [1, 1, -1, -1, 1, 1, -1, -1], [1, 1, 1, 1, -1, -1, -1, -1], [1, 1, -1, -1, -1, -1, 1, 1], [2, -2, 0, 0, 0, 0, 0, 0]]),
    "C2vB": ("OH", C2V_B, ["A'", "B'", "A''", "B''"], [[1, 1, 1, 1], [1, -1, 1, -1], [1, 1, -1, -1], [1, -1, -1, 1]]),
}
D4H_NAMES = ["A_{1g}", "A_{2g}", "B_{1g}", "B_{2g}", "E_g", "A_{1u}", "A_{2u}", "B_{1u}", "B_{2u}", "E_u"]
D3D_NAMES = ["A_{1g}", "A_{2g}", "E_g", "A_{1u}", "A_{2u}", "E_u"]


def _sub_chars(key: str) -> sp.Matrix:
    grp, el, labels, ch = SUB[key]
    if key == "T":
        return sp.Matrix([[T_CHARS[a, TGRP.class_of(h)] for h in el] for a in range(4)])
    return sp.Matrix(ch)


def _d4h_chars() -> sp.Matrix:
    base = _sub_chars("D4")
    return sp.Matrix([[base[a, k % 8] * (1 if (k < 8 or p == 0) else -1) for k in range(16)] for p in (0, 1) for a in range(5)])


def _d3d_chars() -> sp.Matrix:
    base = _sub_chars("D3")
    return sp.Matrix([[base[a, k % 6] * (1 if (k < 6 or p == 0) else -1) for k in range(12)] for p in (0, 1) for a in range(3)])


SUB["D4h"] = ("OH", D4H, D4H_NAMES, None)
SUB["D3d"] = ("OH", D3D, D3D_NAMES, None)


def sub_chars(key: str) -> sp.Matrix:
    """子群的特徵標表（列依 SUB[key][2] 的標記、行依 SUB[key][1] 的元素）。"""
    if key == "D4h":
        return _d4h_chars()
    if key == "D3d":
        return _d3d_chars()
    return _sub_chars(key)


def correlation(key: str) -> sp.Matrix:
    """相關表 𝒟^α↓H（列 O 或 O_h 的 irrep（書上順序）、行 H 的 irrep，SUB[key][2]）。"""
    grp, el, labels, _ = SUB[key]
    G, tab = (O, O_CHARS) if grp == "O" else (OH, OH_CHARS)
    return G.correlation(el, sub_chars(key), tab)


CORR_BOOK = {   # 書上的相關表（列 A1 A2 E T1 T2，或 O_h 的十個）
    "D4": [[1, 0, 0, 0, 0], [0, 0, 1, 0, 0], [1, 0, 1, 0, 0], [0, 1, 0, 0, 1], [0, 0, 0, 1, 1]],           # (4.2.42a)
    "C4": [[1, 0, 0, 0], [0, 0, 1, 0], [1, 0, 1, 0], [1, 1, 0, 1], [0, 1, 1, 1]],                          # (4.2.42b)
    "D3": [[1, 0, 0], [0, 1, 0], [0, 0, 1], [0, 1, 1], [1, 0, 1]],                                         # (4.2.44a)
    "C3": [[1, 0, 0], [1, 0, 0], [0, 1, 1], [1, 1, 1], [1, 1, 1]],                                         # (4.2.44b)
    "C2v": [[1, 0, 0, 0], [0, 1, 0, 0], [1, 1, 0, 0], [0, 1, 1, 1], [1, 0, 1, 1],
            [0, 0, 1, 0], [0, 0, 0, 1], [0, 0, 1, 1], [1, 1, 0, 1], [1, 1, 1, 0]],                          # (4.2.46a)
    "C3v": [[1, 0, 0], [0, 1, 0], [0, 0, 1], [0, 1, 1], [1, 0, 1], [0, 1, 0], [1, 0, 0], [0, 0, 1], [1, 0, 1], [0, 1, 1]],  # (4.2.46b)
    "C4v": [[1, 0, 0, 0, 0], [0, 1, 0, 0, 0], [1, 1, 0, 0, 0], [0, 0, 1, 0, 1], [0, 0, 0, 1, 1],
            [0, 0, 1, 0, 0], [0, 0, 0, 1, 0], [0, 0, 1, 1, 0], [1, 0, 0, 0, 1], [0, 1, 0, 0, 1]],          # (4.2.46c)
    "D2'": [[1, 0, 0, 0], [0, 0, 1, 0], [1, 0, 1, 0], [0, 1, 1, 1], [1, 1, 0, 1]],                         # (4.2.48a)
    "D2": [[1, 0, 0, 0], [1, 0, 0, 0], [2, 0, 0, 0], [0, 1, 1, 1], [0, 1, 1, 1]],                          # (4.2.48b)（行 A1 B1 B2 A2）
}


def conjugates(H: list) -> list:
    """H 的共軛子群 tHt⁻¹（相異的，依出現順序）。"""
    out = []
    for t in O.names:
        Hc = sorted({O.conj(t, h) for h in H}, key=O.idx.get)
        if Hc not in out:
            out.append(Hc)
    return out


# ---------------------------------------------------------------- 子群鏈的冪等元分裂 (4.2.6)–(4.2.18)
def ga_proj(G, el: list, chars: list, ell: int = 1) -> dict:
    """子群（或群）的投影算子 (ℓ/°H)Σ χ(h)* h（群代數元素）。"""
    return G.ga({h: ell * sp.conjugate(c) / len(el) for h, c in zip(el, chars)})


def o_central(alpha: str) -> dict:
    """ℙ^α（O 的群代數元素，(4.1.10)）。"""
    v = O_CLASS_IDEMP_BOOK[alpha]
    return O.from_class_vec(v)


def d4_idempotents() -> dict:
    """D4 ⊃ D2 的六個不可約冪等元 (4.2.8)：P^E_1 = P^E p^{B1}、P^E_2 = P^E p^{B2}、P^{A1}、P^{A2}、P^{B1}、P^{B2}。"""
    ch4, ch2 = sub_chars("D4"), sub_chars("D2")
    P = {lab: ga_proj(O, D4Z, list(ch4.row(a)), 2 if lab == "E" else 1) for a, lab in enumerate(SUB["D4"][2])}
    p = {lab: ga_proj(O, D2N, list(ch2.row(a))) for a, lab in enumerate(SUB["D2"][2])}
    return {"E1": O.ga_mul(P["E"], p["B_1"]), "E2": O.ga_mul(P["E"], p["B_2"]), "A1": P["A_1"], "A2": P["A_2"],
            "B1": P["B_1"], "B2": P["B_2"], "p": p, "P": P}


def d3_idempotents() -> dict:
    """D3 ⊃ C2 的四個不可約冪等元 (4.2.24)：P^E P^A、P^E P^B、P^{A2}P^B、P^{A1}P^A（C2 = {1, i4}）。"""
    ch3, ch2 = sub_chars("D3"), sub_chars("C2")
    P = {lab: ga_proj(O, D3T, list(ch3.row(a)), 2 if lab == "E" else 1) for a, lab in enumerate(SUB["D3"][2])}
    p = {lab: ga_proj(O, C2T, list(ch2.row(a))) for a, lab in enumerate(SUB["C2"][2])}
    return {"EA": O.ga_mul(P["E"], p["A"]), "EB": O.ga_mul(P["E"], p["B"]), "A2B": O.ga_mul(P["A_2"], p["B"]),
            "A1A": O.ga_mul(P["A_1"], p["A"])}


def chain_split(alpha: str, chain: str = "tetragonal") -> list:
    """ℙ^α 乘上子群鏈的不可約冪等元，留下不為零的（(4.2.9)、(4.2.18a)、(4.2.25)）。"""
    ids = d4_idempotents() if chain == "tetragonal" else d3_idempotents()
    keys = ["E1", "E2", "A1", "A2", "B1", "B2"] if chain == "tetragonal" else ["EA", "EB", "A2B", "A1A"]
    return O.split(o_central(alpha), [ids[k] for k in keys])


def elementary_from_irrep(D: dict, ell: int) -> dict:
    """P^α_{jk} = (ℓ/°G) Σ 𝒟_{jk}(g)* g（(3.4.19)、(4.4.1)），{(j, k): 群代數元素}。"""
    return {(j, k): O.ga({g: ell * sp.conjugate(D[g][j, k]) / 24 for g in O_NAMES}) for j in range(ell) for k in range(ell)}


P410_BOOK = {   # (4.2.10)
    0: {"1": 1, "R1^2": 1, "R2^2": -1, "R3^2": -1, "R1": 1, "R1^3": 1, "i5": -1, "i6": -1},
    1: {"1": 1, "R1^2": -1, "R2^2": 1, "R3^2": -1, "R2": 1, "R2^3": 1, "i1": -1, "i2": -1},
    2: {"1": 1, "R1^2": -1, "R2^2": -1, "R3^2": 1, "R3": 1, "R3^3": 1, "i3": -1, "i4": -1}}


def p418(which: str = "A1", prefactor=R_(2, 3)) -> dict:
    """(4.2.18b,c)：prefactor·(P^{a}(D4) − ½ r1 P^{a} − ½ r1² P^{a})（書上印 3/4，正確 2/3）。"""
    Pa = d4_idempotents()["P"]["A_1" if which == "A1" else "B_1"]
    X = O.ga_add(Pa, O.ga_mul(O.el("r1"), Pa), O.ga_mul(O.el("r1^2"), Pa), coeffs=[1, -R_(1, 2), -R_(1, 2)])
    return O.ga_scale(X, prefactor)


# ---------------------------------------------------------------- 三方（D3 ⊃ C2）基底、行進波基底
T43 = sp.Matrix([[1 / SQ2, 1 / SQ6, 1 / SQ3], [-1 / SQ2, 1 / SQ6, 1 / SQ3], [0, -2 / SQ6, 1 / SQ3]])   # (4.2.31)：行 v1 v2 v3
TMOV4 = sp.Matrix([[-1 / SQ2, 1 / SQ2, 0], [-IU / SQ2, -IU / SQ2, 0], [0, 0, 1]])                     # (4.2.36)：行 1_4 3_4 0_4
TMOV2 = TMOV4[:2, :2]                                                                                  # (x, y) → (1_4, 3_4)


def _trig_basis(Dg: dict) -> sp.Matrix:
    """在 3 維表示 Dg 裡找 D3 ⊃ C2 的基底：(E A, E B, A2 B)（T1）或 (A1 A, E A, E B)（T2）。
    E 夥伴的相位：𝒟(r1) 的 E 塊 ＝ rot(120°)（與 T1 的三方 E 塊相同），C2 = {1, i4}。每個向量第一個非零分量取正。"""
    def proj(chars, ell):
        return sum((ell * c * Dg[h] for h, c in zip(D3T, chars)), sp.zeros(3)) / 6

    def unit(M):
        for k in range(3):
            v = M[:, k]
            if not _zero(v):
                v = sp.simplify(v / sp.sqrt(sp.simplify((v.T * v)[0])))
                s = next(x for x in v if sp.simplify(x) != 0)
                return sp.simplify(v * sp.sign(s))
        return None

    pA = (sp.eye(3) + Dg["i4"]) / 2
    eA = unit(proj([2, -1, -1, 0, 0, 0], 1) * pA)
    eB = sp.simplify((Dg["r1"] * eA + eA / 2) / (SQ3 / 2))
    a1 = unit(proj([1] * 6, 1))
    a2 = unit(proj([1, 1, 1, -1, -1, -1], 1))
    if a1 is not None:
        return sp.Matrix.hstack(a1, eA, eB)
    return sp.Matrix.hstack(eA, eB, a2)


def _conj_rep(D: dict, U) -> dict:
    Ud = U.H
    return {g: (Ud * D[g] * U).applyfunc(lambda x: sp.radsimp(sp.expand(x))) for g in D}


TRIG_U_T1 = _trig_basis(TETRA["T_1"])        # ＝ T43（由投影算出，check 核對）
TRIG_U_T2 = _trig_basis(TETRA["T_2"])        # 行 |T2 A1 A⟩、|T2 E A⟩、|T2 E B⟩（(4.2.33)）
TRIG = {"A_1": TETRA["A_1"], "A_2": TETRA["A_2"], "E": TETRA["E"],
        "T_1": _conj_rep(TETRA["T_1"], T43), "T_2": _conj_rep(TETRA["T_2"], TRIG_U_T2)}       # O ⊃ D3 ⊃ C2
TRIG_LABELS = {"A_1": [("A_1", "A")], "A_2": [("A_2", "B")], "E": [("E", "A"), ("E", "B")],
               "T_1": [("E", "A"), ("E", "B"), ("A_2", "B")], "T_2": [("A_1", "A"), ("E", "A"), ("E", "B")]}
TRIG_T1_BOOK = {   # (4.2.28) 的類代表（照印）
    "i4": [[1, 0, 0], [0, -1, 0], [0, 0, -1]], "r1": [[-HALF, -RT3_2, 0], [RT3_2, -HALF, 0], [0, 0, 1]],
    "R1^2": [[0, SQ3 / 3, SQ6 / 3], [SQ3 / 3, -R_(2, 3), SQ2 / 3], [SQ6 / 3, SQ2 / 3, -R_(1, 3)]],
    "R3": [[0, -SQ3 / 3, -SQ6 / 3], [SQ3 / 3, R_(2, 3), -SQ2 / 3], [SQ6 / 3, -SQ2 / 3, R_(1, 3)]]}
TRIG_T1_BOOK_TYPOS = {("r2^2", 1, 0): (SQ3 / 3, SQ3 / 6), ("r2", 1, 0): (-SQ3 / 3, -SQ3 / 6), ("R2^2", 1, 2): (SQ2 / 2, SQ2 / 3)}  # (書, 應為)

MOV4 = {"A_1": TETRA["A_1"], "A_2": TETRA["A_2"], "E": TETRA["E"],
        "T_1": _conj_rep(TETRA["T_1"], TMOV4), "T_2": _conj_rep(TETRA["T_2"], TMOV4)}         # O ⊃ D4 ⊃ C4 行進波
MOV4_LABELS = {"A_1": ["0_4"], "A_2": ["2_4"], "E": ["0_4", "2_4"], "T_1": ["1_4", "3_4", "0_4"], "T_2": ["1_4", "3_4", "2_4"]}
MOV4_T1_BOOK = {   # (4.2.35) 的類代表
    "R1^2": [[0, -1, 0], [-1, 0, 0], [0, 0, -1]], "R3": [[-IU, 0, 0], [0, IU, 0], [0, 0, 1]],
    "i4": [[0, -IU, 0], [IU, 0, 0], [0, 0, -1]],
    "r1": [[-IU / 2, IU / 2, -1 / SQ2], [-IU / 2, IU / 2, 1 / SQ2], [-IU / SQ2, -IU / SQ2, 0]],
    "R1": [[HALF, -HALF, -IU / SQ2], [-HALF, HALF, -IU / SQ2], [-IU / SQ2, -IU / SQ2, 0]],
    "i1": [[-HALF, -HALF, -1 / SQ2], [-HALF, -HALF, 1 / SQ2], [-1 / SQ2, 1 / SQ2, 0]]}
D436_BOOK = sp.diag(IU, -IU, 1)                 # (4.2.36) 照印（錯）；正確 diag(−i, i, 1)

TMOV3 = TMOV4                                   # 三方行進波：|1_3⟩ = (−v1 − i v2)/√2、|2_3⟩ = (v1 − i v2)/√2、|0_3⟩ = v3（仿 (4.2.36)）
_P3_T2 = sp.Matrix([[0, 1, 0], [0, 0, 1], [1, 0, 0]])   # (A1, E_A, E_B) → (E_A, E_B, A1) 的重排
MOV3 = {"A_1": TETRA["A_1"], "A_2": TETRA["A_2"], "E": _conj_rep(TETRA["E"], TMOV2),
        "T_1": _conj_rep(TRIG["T_1"], TMOV3),
        "T_2": _conj_rep(TRIG["T_2"], sp.Matrix([[0, 0, 1], [-1 / SQ2, 1 / SQ2, 0], [-IU / SQ2, -IU / SQ2, 0]]))}  # O ⊃ D3 ⊃ C3
MOV3_LABELS = {"A_1": ["0_3"], "A_2": ["0_3"], "E": ["1_3", "2_3"], "T_1": ["1_3", "2_3", "0_3"], "T_2": ["1_3", "2_3", "0_3"]}
BASES = {"tetragonal": TETRA, "trigonal": TRIG, "c4": MOV4, "c3": MOV3}


def irrep(alpha: str, g: str, basis: str = "tetragonal") -> sp.Matrix:
    """O 的不可約表示矩陣 𝒟^α(g)（basis：tetragonal（O⊃D4⊃D2）、trigonal（O⊃D3⊃C2）、c4（O⊃D4⊃C4）、c3（O⊃D3⊃C3））；
    O_h 的元素（"I.r1"…）與 α 帶 g/u（"T_{1u}"）也可以：𝒟^{αu}(Ig) = −𝒟^α(g)（(4.1.18)）。"""
    rot, inv = oh_split(g)
    base, parity = alpha, "g"
    if alpha.endswith("g}") or alpha.endswith("u}"):
        base, parity = alpha[:-2] + "}", alpha[-2]
        base = base.replace("{", "").replace("}", "")
    elif alpha.endswith("_g") or alpha.endswith("_u"):
        base, parity = alpha[:-2], alpha[-1]
    M = BASES[basis][base][rot]
    return -M if (inv and parity == "u") else M


# ======================================================================== §4.3 誘導表示、Frobenius、六點模型
SITES = [(0, 0, 1), (0, 0, -1), (1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0)]   # |1⟩…|6⟩ = +z, −z, +x, −x, +y, −y
SITE_LEADERS = ["1", "R1^2", "r1", "r2", "r1^2", "r2^2"]                       # (4.3.10)：|g⟩ = g|1⟩


def site_perm(g: str, G=None) -> sp.Matrix:
    """g 在六個頂點上的置換矩陣（⟨i|g|j⟩ = 1 若 g·site_j = site_i；O 或 O_h 的元素）。"""
    G = G or (OH if g not in O.names else O)
    M = G.num[g].real
    P = sp.zeros(6)
    for j, s in enumerate(SITES):
        v = M @ np.array(s, float)
        i = next(i for i, t in enumerate(SITES) if np.allclose(v, t))
        P[i, j] = 1
    return P


def site_of(g: str) -> int:
    v = (OH.num[g].real @ np.array([0, 0, 1.0]))
    return next(i for i, t in enumerate(SITES) if np.allclose(v, t))


def induced(G, sub: list, Dsub: dict, leaders: list | None = None) -> dict:
    """誘導表示 (4.3.34)：ℐ_{ij}(g) = D^a(l_i⁻¹ g l_j)（若在 H 內）否則 0；Dsub = {h: d×d 矩陣}；leaders 預設取每個左陪集
    在 names 中最前的元素。回傳 {g: (n·d)×(n·d) Matrix}。"""
    if leaders is None:
        leaders = [C[0] for C in G.cosets(sub, "left")]
    d = next(iter(Dsub.values())).shape[0]
    n = len(leaders)
    Hs_ = set(sub)
    out = {}
    for g in G.names:
        M = sp.zeros(n * d)
        for i, li in enumerate(leaders):
            for j, lj in enumerate(leaders):
                x = G.mul(G.inv(li), g, lj)
                if x in Hs_:
                    M[i * d:(i + 1) * d, j * d:(j + 1) * d] = Dsub[x]
        out[g] = M
    return out


def induced_char(G, sub: list, chi: dict) -> list:
    """誘導特徵標（公式）：χ(g) = (1/°H) Σ_{x∈G, x⁻¹gx∈H} χ_a(x⁻¹gx)，每個類一個值（依 G.classes）。"""
    out = []
    Hs_ = set(sub)
    for C in G.classes:
        g = C[0]
        s = 0
        for x in G.names:
            y = G.conj(G.inv(x), g)
            if y in Hs_:
                s += chi[y]
        out.append(sp.nsimplify(sp.simplify(s / len(sub))))
    return out


def sub_irreps(key: str) -> dict:
    """子群的不可約表示矩陣 {標記: {h: Matrix}}（一維由特徵標；D4 的 E ＝ T1 的 xy 塊；D3 的 E ＝ 三方 T1 的 E 塊；
    T 的 T ＝ T1；C3v／C4v 的 E 同理）。"""
    grp, el, labels, _ = SUB[key]
    ch = sub_chars(key)
    out = {}
    for a, lab in enumerate(labels):
        if ch[a, 0] == 1:
            out[lab] = {h: sp.Matrix([[ch[a, p]]]) for p, h in enumerate(el)}
        elif key in ("D4", "C4v"):
            out[lab] = {h: sp.Matrix(OH.mats[h][:2, :2]) for h in el}
        elif key in ("D3", "C3v"):
            out[lab] = {h: TRIG["T_1"][oh_split(h)[0]][:2, :2] for h in el}
        elif key == "T":
            out[lab] = {h: O.mats[h] for h in el}
        else:
            raise ValueError(key)
    return out


def frobenius(key: str) -> dict:
    """Frobenius 互反 (4.3.29)：對 H 的每個 irrep a，用誘導矩陣的跡與誘導特徵標公式各算一次 f^α(a↑G)，
    與相關表的第 a 行比較。回傳 {a: (由矩陣, 由公式, 相關表的行)}（三者都是 list）。"""
    grp, el, labels, _ = SUB[key]
    G, tab = (O, O_CHARS) if grp == "O" else (OH, OH_CHARS)
    corr = correlation(key)
    irr = sub_irreps(key)
    out = {}
    for b, lab in enumerate(labels):
        Dm = irr[lab]
        ind = induced(G, el, Dm)
        tr = [sp.nsimplify(sp.simplify(ind[C[0]].trace())) for C in G.classes]
        f_mat = G.frequencies(tr, tab)
        f_form = G.frequencies(induced_char(G, el, {h: Dm[h].trace() for h in el}), tab)
        out[lab] = (f_mat, f_form, list(corr.col(b)))
    return out


# ---------------------------------------------------------------- 六點八面體模型 (4.3.2)–(4.3.28)
def six_H(H=Hs, T=Ts, S=Ss) -> sp.Matrix:
    """(4.3.3)：對角 H、對面（1-2、3-4、5-6）T、其餘 S。"""
    return sp.Matrix(6, 6, lambda i, j: H if i == j else (T if i // 2 == j // 2 else S))


def six_accel(k=sp.Symbol("k"), t=sp.Symbol("t")) -> sp.Matrix:
    """(4.3.2)：古典版 −⟨a⟩ = (h, t, s)，h = 2k + t、s = k/2。"""
    return six_H(2 * k + t, t, k / 2)


SIX_LEVELS = {"A_1": Hs + Ts + 4 * Ss, "T_1": Hs - Ts, "E": Hs + Ts - 2 * Ss}          # (4.3.28)
SIX_VECTORS_BOOK = {   # (4.3.20)、(4.3.22)–(4.3.27)
    ("A_1", "0_4"): sp.Matrix([1, 1, 1, 1, 1, 1]) / SQ6, ("E", "0_4"): sp.Matrix([2, 2, -1, -1, -1, -1]) / (2 * SQ3),
    ("E", "2_4"): sp.Matrix([0, 0, 1, 1, -1, -1]) / 2, ("T_1", "1_4"): sp.Matrix([0, 0, -1, 1, -IU, IU]) / 2,
    ("T_1", "3_4"): sp.Matrix([0, 0, 1, -1, -IU, IU]) / 2, ("T_1", "0_4"): sp.Matrix([1, -1, 0, 0, 0, 0]) / SQ2,
    ("T_1", "x"): sp.Matrix([0, 0, 1, -1, 0, 0]) / SQ2, ("T_1", "y"): sp.Matrix([0, 0, 0, 0, 1, -1]) / SQ2}
SIX_INDUCED_I4_BOOK = sp.Matrix([[0, 1, 0, 0, 0, 0], [1, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 1], [0, 0, 0, 0, 1, 0],
                                 [0, 0, 0, 1, 0, 0], [0, 0, 1, 0, 0, 0]])                         # (4.3.12)
SIX_INDUCED_1_4_I4_BOOK = sp.Matrix([[0, IU, 0, 0, 0, 0], [-IU, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, -IU], [0, 0, 0, 0, -IU, 0],
                                     [0, 0, 0, IU, 0, 0], [0, 0, IU, 0, 0, 0]])                   # (4.3.34ₓ)


def six_induced(label: str = "0_4") -> dict:
    """D^{label}(C4)↑O，陪集代表 SITE_LEADERS（(4.3.9)–(4.3.12)、(4.3.34ₓ)）。"""
    return induced(O, C4Z, sub_irreps("C4")[label], SITE_LEADERS)


def six_projected(alpha: str, j: int) -> sp.Matrix:
    """P^α_{j,0_4}|1⟩ 歸一化（(4.3.19)、(3.4.19)：(ℓ/24) Σ_g 𝒟_{j a}(g)* g|1⟩；a = 0_4；T1 用行進波 (4.2.35)、E 用四方）。"""
    D = MOV4[alpha]
    ell = D["1"].shape[0]
    a = MOV4_LABELS[alpha].index("0_4")
    v = sp.zeros(6, 1)
    for g in O_NAMES:
        v[site_of(g)] += ell * sp.conjugate(D[g][j, a]) / 24
    nrm = sp.sqrt(sp.simplify((v.H * v)[0]))
    return sp.simplify(v / nrm)


def six_coset_vector(alpha: str, j: int) -> sp.Matrix:
    """陪集分解 (4.3.44)：(ℓ°H/°G)^{1/2} Σ_l 𝒟_{ja}(l)* |l⟩（只對六個陪集代表求和）。"""
    D = MOV4[alpha]
    ell = D["1"].shape[0]
    a = MOV4_LABELS[alpha].index("0_4")
    v = sp.Matrix([sp.conjugate(D[l][j, a]) for l in SITE_LEADERS])
    return sp.simplify(sp.sqrt(R_(ell * 4, 24)) * v)


def six_coset_energy(alpha: str, H=Hs, T=Ts, S=Ss):
    """(4.3.46)：⟨ja|H|ja⟩ = Σ_l ⟨1|H|l⟩ 𝒟_{aa}(l)*（只用 (4.3.3) 的第一列）。"""
    D = MOV4[alpha]
    a = MOV4_LABELS[alpha].index("0_4")
    row = list(six_H(H, T, S).row(0))
    return sp.simplify(sum(row[site_of(l)] * sp.conjugate(D[l][a, a]) for l in SITE_LEADERS))


SIX_E_BOOK_446X = Hs - Ts + 2 * Ss            # (4.3.46ₓ) 照印（錯）


def six_numeric(H: float = 0.0, T: float = 0.0, S: float = -1.0) -> np.ndarray:
    return np.linalg.eigvalsh(np.array(six_H(H, T, S), dtype=float))


# ---------------------------------------------------------------- §4.3D 自發對稱破缺：其他誘導表示
def d6_induced_examples() -> dict:
    """(4.3.49)–(4.3.51) 與正則表示：{名: 由誘導特徵標公式算的頻率（D6 依 SN.D6_IRREP_NAMES；O 依書）}。"""
    _, T6, _ = SN.book_order(SN.D6, SN.D6_CHARS_BOOK)
    C6 = SN.D6.subgroup_of("h")
    out = {}
    for lab, chi in (("0_2", {"1": 1, "p3": 1}), ("1_2", {"1": 1, "p3": -1})):
        out[lab] = SN.D6.frequencies(induced_char(SN.D6, ["1", "p3"], chi), T6)
    for k in range(6):
        chi = {h: SN.SA.eps(6, k * p) for p, h in enumerate(C6)}
        out[f"{k}_6"] = SN.D6.frequencies(induced_char(SN.D6, C6, chi), T6)
    out["regular_O"] = O.frequencies(induced_char(O, ["1"], {"1": 1}), O_CHARS)
    return out


# ---------------------------------------------------------------- §4.3E 外加對稱破缺
def vq_matrix(Q=Qs, H=Hs, S=Ss, T=0) -> sp.Matrix:
    """(4.3.52)：⟨H + V_Q⟩（|1⟩、|2⟩ 加 Q；書上取 T = 0）。"""
    return six_H(H, T, S) + sp.diag(Q, Q, 0, 0, 0, 0)


VQ_BASIS = sp.Matrix.hstack(sp.Matrix([0, 0, 1, -1, 0, 0]) / SQ2, sp.Matrix([0, 0, 0, 0, 1, -1]) / SQ2,
                            sp.Matrix([1, -1, 0, 0, 0, 0]) / SQ2, sp.Matrix([1] * 6) / SQ6,
                            sp.Matrix([2, 2, -1, -1, -1, -1]) / (2 * SQ3), sp.Matrix([0, 0, 1, 1, -1, -1]) / 2)   # (4.3.54)–(4.3.55)
VQ_LABELS = ["E_u", "E_u", "A_{2u}", "A_{1g}", "A_{1g}", "B_{1g}"]
VQ_REDUCED_BOOK = sp.Matrix([[Hs, 0, 0, 0, 0, 0], [0, Hs, 0, 0, 0, 0], [0, 0, Hs + Qs, 0, 0, 0],
                             [0, 0, 0, Hs + 4 * Ss + Qs / 3, SQ2 * Qs / 3, 0], [0, 0, 0, SQ2 * Qs / 3, Hs - 2 * Ss + 2 * Qs / 3, 0],
                             [0, 0, 0, 0, 0, Hs - 2 * Ss]])                                       # (4.3.56)
VQ_SLOPES = {"A_{2u}": 1, "A_{1g}": R_(1, 3), "E_g[A_{1g}]": R_(2, 3), "E_u": 0, "B_{1g}": 0}
VQ_HIGHQ = sp.Matrix([[1 / SQ3, SQ2 / SQ3], [SQ2 / SQ3, -1 / SQ3]])                               # (4.3.58)–(4.3.60) 的變換
VQ_460_BOOK = sp.Matrix([[Hs + Qs, Ss * sp.sqrt(8)], [Ss * sp.sqrt(8), Hs + 2 * Ss]])             # (4.3.60)
SCAN = dict(H=0.0, S=-1.0, T=0.0, Qmax=8.0, Dmax=8.0, n=41)     # Fig. 4.3.5、4.3.7 的本系列數字（書上只說「固定的負 S」）


def vq_reduced(Q=Qs, H=Hs, S=Ss) -> sp.Matrix:
    return sp.simplify(VQ_BASIS.T * vq_matrix(Q, H, S) * VQ_BASIS)


def vq_levels(Q: float, H: float = SCAN["H"], S: float = SCAN["S"]) -> list:
    """[(能量, 標記)]（由小到大）：E_u 二重、A_{2u}、兩個 A_{1g}（2×2）、B_{1g}。"""
    m2 = np.array([[H + 4 * S + Q / 3, math.sqrt(2) * Q / 3], [math.sqrt(2) * Q / 3, H - 2 * S + 2 * Q / 3]])
    a, b = np.linalg.eigvalsh(m2)
    out = [(H, "E_u"), (H, "E_u"), (H + Q, "A_{2u}"), (a, "A_{1g}"), (b, "A_{1g}"), (H - 2 * S, "B_{1g}")]
    return sorted(out, key=lambda t: t[0])


def vq_scan(Qs_=None, H: float = SCAN["H"], S: float = SCAN["S"]) -> dict:
    """Fig. 4.3.5：{標記: array}（"A1g-"、"A1g+" 是 2×2 的低、高根；"Eu"、"A2u"、"B1g"）。"""
    Qs_ = np.linspace(-SCAN["Qmax"], SCAN["Qmax"], SCAN["n"]) if Qs_ is None else np.asarray(Qs_, float)
    out = {k: [] for k in ("Q", "A1g-", "A1g+", "Eu", "A2u", "B1g")}
    for Q in Qs_:
        m2 = np.array([[H + 4 * S + Q / 3, math.sqrt(2) * Q / 3], [math.sqrt(2) * Q / 3, H - 2 * S + 2 * Q / 3]])
        a, b = np.linalg.eigvalsh(m2)
        for k, v in (("Q", Q), ("A1g-", a), ("A1g+", b), ("Eu", H), ("A2u", H + Q), ("B1g", H - 2 * S)):
            out[k].append(v)
    return {k: np.array(v) for k, v in out.items()}


def vq_avoided(S=Ss) -> tuple:
    """兩個 A_{1g} 的迴避交叉：最小間距在 Q = 2S，間距 4√2|S|（由 (4.3.60)）。"""
    gap2 = lambda Q: (-6 * S + Q / 3) ** 2 + 8 * Q ** 2 / 9      # noqa: E731
    Qm = sp.solve(sp.diff(gap2(Qs), Qs), Qs)[0]
    return sp.simplify(Qm), sp.simplify(sp.sqrt(gap2(Qm)))


def vd_matrix(D=Ds, H=Hs, S=Ss, T=0) -> sp.Matrix:
    """(4.3.61)：|1⟩ 降 D、|2⟩ 升 D。"""
    return six_H(H, T, S) + sp.diag(-D, D, 0, 0, 0, 0)


VD_BASIS = sp.Matrix.hstack(sp.Matrix([0, 0, 1, -1, 0, 0]) / SQ2, sp.Matrix([0, 0, 0, 0, 1, -1]) / SQ2,
                            sp.Matrix([1, -1, 0, 0, 0, 0]) / SQ2, sp.Matrix([1] * 6) / SQ6,
                            sp.Matrix([2, 2, -1, -1, -1, -1]) / (2 * SQ3), sp.Matrix([0, 0, 1, 1, -1, -1]) / 2)   # (4.3.63)–(4.3.64)（更正）
VD_T1U_BOOK = sp.Matrix([1, 1, 0, 0, 0, 0]) / SQ2          # (4.3.63) 照印（錯）
VD_LABELS = ["E", "E", "A'", "A'", "A'", "B'"]
VD_REDUCED = sp.Matrix([[Hs, 0, 0, 0, 0, 0], [0, Hs, 0, 0, 0, 0],
                        [0, 0, Hs, -Ds / SQ3, -SQ2 * Ds / SQ3, 0], [0, 0, -Ds / SQ3, Hs + 4 * Ss, 0, 0],
                        [0, 0, -SQ2 * Ds / SQ3, 0, Hs - 2 * Ss, 0], [0, 0, 0, 0, 0, Hs - 2 * Ss]])   # (4.3.66) 更正（書上 (4,4) 印 H − 4S）
VD_469 = sp.Matrix([[1 / SQ2, 1 / SQ6, 1 / SQ3], [0, 2 / SQ6, -1 / SQ3], [-1 / SQ2, 1 / SQ6, 1 / SQ3]])   # (4.3.69) 左邊的變換
VD_469_RESULT = sp.Matrix([[Hs - Ds, 2 * Ss, 0], [2 * Ss, Hs + 2 * Ss, 2 * Ss], [0, 2 * Ss, Hs + Ds]])
VD_SECOND_ORDER = {"A_{1g}": 1 / (12 * Ss), "T_{1u}": 1 / (4 * Ss), "E_g": -1 / (3 * Ss)}   # D² 的係數（二階微擾）


def vd_reduced(D=Ds, H=Hs, S=Ss) -> sp.Matrix:
    return sp.simplify(VD_BASIS.T * vd_matrix(D, H, S) * VD_BASIS)


def vd_scan(Ds_=None, H: float = SCAN["H"], S: float = SCAN["S"]) -> dict:
    """Fig. 4.3.7：{"D", "A'1", "A'2", "A'3"（3×3 由小到大）, "E", "B'"}。"""
    Ds_ = np.linspace(-SCAN["Dmax"], SCAN["Dmax"], SCAN["n"]) if Ds_ is None else np.asarray(Ds_, float)
    out = {k: [] for k in ("D", "A'1", "A'2", "A'3", "E", "B'")}
    for Dv in Ds_:
        m3 = np.array([[H, -Dv / math.sqrt(3), -math.sqrt(2) * Dv / math.sqrt(3)], [-Dv / math.sqrt(3), H + 4 * S, 0],
                       [-math.sqrt(2) * Dv / math.sqrt(3), 0, H - 2 * S]])
        w = np.linalg.eigvalsh(m3)
        for k, v in (("D", Dv), ("A'1", w[0]), ("A'2", w[1]), ("A'3", w[2]), ("E", H), ("B'", H - 2 * S)):
            out[k].append(v)
    return {k: np.array(v) for k, v in out.items()}


def vb_matrix(beta: float, H: float = SCAN["H"], S: float = SCAN["S"], T: float = SCAN["T"]) -> np.ndarray:
    """磁場 C4（本系列的模型，書上只給相關表 (4.3.73)）：赤道環 +x→+y→−x→−y（逆時針看 +z）的穿隧加 Peierls 相位：
    ⟨下一個|H|這一個⟩ = S e^{iβ}；極點的穿隧不變。保有 C4、破壞 σ_v。"""
    M = np.array(six_H(H, T, S), dtype=complex)
    ring = [2, 4, 3, 5]                                  # +x, +y, −x, −y 的索引
    for p in range(4):
        a, b = ring[p], ring[(p + 1) % 4]
        M[b, a] = S * np.exp(1j * beta)
        M[a, b] = S * np.exp(-1j * beta)
    return M


def vb_levels(beta: float, H: float = SCAN["H"], S: float = SCAN["S"], T: float = SCAN["T"]) -> dict:
    """{C4 標記: 能量}：1_4 = H + 2S sin β、3_4 = H − 2S sin β、2_4 = H − 2S cos β、T1(0_4) = H − T、
    A1／E(0_4) = 2×2 [[H+T, 2√2S], [2√2S, H + 2S cos β]] 的兩根（β = 0 時 H+T+4S、H+T−2S）。"""
    m2 = np.array([[H + T, 2 * math.sqrt(2) * S], [2 * math.sqrt(2) * S, H + 2 * S * math.cos(beta)]])
    lo, hi = np.linalg.eigvalsh(m2)
    return {"1_4": H + 2 * S * math.sin(beta), "3_4": H - 2 * S * math.sin(beta), "2_4": H - 2 * S * math.cos(beta),
            "T_1(0_4)": H - T, "A_1(0_4)": lo, "E(0_4)": hi}


def v111_matrix(d=ds, H=Hs, S=Ss, T=0) -> sp.Matrix:
    """(4.3.74)：沿 (111) 的電場，|1⟩|3⟩|5⟩（+z +x +y）降 d、|2⟩|4⟩|6⟩ 升 d。"""
    return six_H(H, T, S) + sp.diag(-d, d, -d, d, -d, d)


V111_BASIS = {"A1g": sp.Matrix([1] * 6) / SQ6, "T1uA": sp.Matrix([1, -1, 1, -1, 1, -1]) / SQ6,
              "Eg1": sp.Matrix([2, 2, -1, -1, -1, -1]) / sp.sqrt(12), "T1uE1": sp.Matrix([2, -2, -1, 1, -1, 1]) / sp.sqrt(12),
              "Eg2": sp.Matrix([0, 0, 1, 1, -1, -1]) / 2, "T1uE2": sp.Matrix([0, 0, 1, -1, -1, 1]) / 2}   # (4.3.81)–(4.3.84)
V111_A = sp.Matrix([[Hs + 4 * Ss, -ds], [-ds, Hs]])          # (4.3.80)
V111_E = sp.Matrix([[Hs - 2 * Ss, -ds], [-ds, Hs]])          # (4.3.82)
T485_BOOK = sp.Matrix([[-1 / SQ6, 1 / SQ2, 1 / SQ3], [-1 / SQ6, -1 / SQ2, 1 / SQ3], [2 / SQ6, 0, 1 / SQ3]])   # (4.3.85)


def v111_levels(d: float, H: float = SCAN["H"], S: float = SCAN["S"]) -> list:
    a = np.linalg.eigvalsh(np.array([[H + 4 * S, -d], [-d, H]]))
    e = np.linalg.eigvalsh(np.array([[H - 2 * S, -d], [-d, H]]))
    return sorted([(a[0], "A'"), (a[1], "A'"), (e[0], "E"), (e[0], "E"), (e[1], "E"), (e[1], "E")], key=lambda t: t[0])


def v111_tunnel(Sp=sp.Symbol("S'"), S=Ss, H=Hs) -> dict:
    """(4.3.87)–(4.3.88)：順場—逆場的鄰居穿隧 S'、同側的鄰居 S（D3d）。{標記: 能量}（sympy）。"""
    down = {0, 2, 4}
    M = sp.Matrix(6, 6, lambda i, j: H if i == j else (0 if i // 2 == j // 2 else (S if ((i in down) == (j in down)) else Sp)))
    vec = {"A_{1g}": V111_BASIS["A1g"], "A_{2u}": V111_BASIS["T1uA"], "E_g": V111_BASIS["Eg1"], "E_u": V111_BASIS["T1uE1"]}
    return {k: sp.simplify((v.T * M * v)[0]) for k, v in vec.items()}, M


# ======================================================================== §4.4 六氟化物 XF6
F_SITES = SITES                                              # F 原子 0..5 在 +z, −z, +x, −x, +y, −y；X（中心）是第 6 個
A_LEADERS = ["1", "R1^2", "R1", "R2", "R1^3", "R2^3"]        # (4.4.22) 的欄序
B_LEADERS = ["1", "r1", "r2", "r1^2", "r4^2", "R2^2", "R1", "R2", "R3", "R1^3", "R2^3", "i3"]
XF6_NAMES = [("A" if l == "1" else l + "A") for l in A_LEADERS] + [("B" if l == "1" else l + "B") for l in B_LEADERS] + ["C1", "C2", "C3"]
XF6_TEX = [r"|A\rangle" if l == "1" else rf"|{O_TEX[l]}A\rangle" for l in A_LEADERS] + \
          [r"|B\rangle" if l == "1" else rf"|{O_TEX[l]}B\rangle" for l in B_LEADERS] + [r"|C_1\rangle", r"|C_2\rangle", r"|C_3\rangle"]
XF6_COORDS_LIST = [   # (名稱, 原子, 位移方向)：A = 徑向、B = 切向、C = 中心原子（Fig. 4.4.1）
    ("A", 0, (0, 0, 1)), ("R1^2A", 1, (0, 0, -1)), ("R1A", 5, (0, -1, 0)), ("R2A", 2, (1, 0, 0)), ("R1^3A", 4, (0, 1, 0)),
    ("R2^3A", 3, (-1, 0, 0)), ("B", 0, (1, 0, 0)), ("r1B", 2, (0, 1, 0)), ("r2B", 3, (0, 1, 0)), ("r1^2B", 4, (0, 0, 1)),
    ("r4^2B", 5, (0, 0, 1)), ("R2^2B", 1, (-1, 0, 0)), ("R1B", 5, (1, 0, 0)), ("R2B", 2, (0, 0, -1)), ("R3B", 0, (0, 1, 0)),
    ("R1^3B", 4, (1, 0, 0)), ("R2^3B", 3, (0, 0, 1)), ("i3B", 1, (0, 1, 0)), ("C1", 6, (1, 0, 0)), ("C2", 6, (0, 1, 0)),
    ("C3", 6, (0, 0, 1))]
XF6 = {"UF6": dict(m=19, M=240, b=0.0), "SF6": dict(m=19, M=32, b=0.2, k=0.18), "SF6_34": dict(m=19, M=34, b=0.2, k=0.18)}
UF6_NU = {"A_{1g}": 667, "T_{1u}+": 624, "E_g": 535, "T_{2g}": 201, "T_{1u}-": 184, "T_{2u}": 140}    # cm⁻¹（Fig. 4.4.5）
SF6_NU = {"T_{1u}+": 965, "A_{1g}": 775, "E_g": 644, "T_{1u}-": 617, "T_{2g}": 524, "T_{2u}": 363}    # cm⁻¹（Fig. 4.4.6）
SF6_ISO_EXP = {"+": (965, 948), "-": (615, 612)}                                                       # (4.4.36b)


def xf6_coords() -> list:
    """由群作用產生的 21 個座標 [(名稱, 原子, 向量)]：|lA⟩ = l·（+z 上的 ẑ）、|lB⟩ = l·（+z 上的 x̂）（＝ XF6_COORDS_LIST）。"""
    out = []
    for l in A_LEADERS:
        M = O.num[l].real
        v = M @ np.array([0, 0, 1.0])
        out.append(("A" if l == "1" else l + "A", next(i for i, s in enumerate(F_SITES) if np.allclose(v, s)), tuple(int(round(x)) for x in v)))
    for l in B_LEADERS:
        M = O.num[l].real
        site = M @ np.array([0, 0, 1.0])
        v = M @ np.array([1.0, 0, 0])
        out.append(("B" if l == "1" else l + "B", next(i for i, s in enumerate(F_SITES) if np.allclose(site, s)),
                    tuple(int(round(x)) for x in v)))
    out += [("C1", 6, (1, 0, 0)), ("C2", 6, (0, 1, 0)), ("C3", 6, (0, 0, 1))]
    return out


def _atom_pos(a: int) -> np.ndarray:
    return np.zeros(3) if a == 6 else np.array(F_SITES[a], float)


def _xf6_basis() -> np.ndarray:
    """21×21：座標 → 原子笛卡兒位移（7 個原子 × 3）。"""
    B = np.zeros((21, 21))
    for c, (_, a, v) in enumerate(XF6_COORDS_LIST):
        B[3 * a:3 * a + 3, c] = v
    return B


def xf6_force(k: float, j: float, b: float = 0.0) -> np.ndarray:
    """由彈簧幾何算 21×21 力矩陣（座標依 XF6_NAMES）：F–F 最近鄰彈簧 k（12 根）、X–F 徑向 j（6 根）、
    彎曲 b：每個 90° 的 F–X–F 角 V = (b/4)Δθ²（Δθ = −[e₁·(u₂ − u_X) + e₂·(u₁ − u_X)]，使 ⟨B|F|B⟩ 的彎曲部分 = b）。"""
    Fc = np.zeros((21, 21))

    def spring(a1, a2, kk):
        n = _atom_pos(a2) - _atom_pos(a1)
        n = n / np.linalg.norm(n)
        g = np.zeros(21)
        g[3 * a2:3 * a2 + 3] += n
        g[3 * a1:3 * a1 + 3] -= n
        return kk * np.outer(g, g)

    for a1 in range(6):
        Fc += spring(6, a1, j)
        for a2 in range(a1 + 1, 6):
            e1, e2 = _atom_pos(a1), _atom_pos(a2)
            if abs(e1 @ e2) < 1e-9:
                Fc += spring(a1, a2, k)
                g = np.zeros(21)
                g[3 * a2:3 * a2 + 3] -= e1
                g[3 * a1:3 * a1 + 3] -= e2
                g[18:21] += e1 + e2
                Fc += (b / 2) * np.outer(g, g)
    B = _xf6_basis()
    return B.T @ Fc @ B


def xf6_mass(m: float = 19, M: float = 240) -> np.ndarray:
    return np.diag([m] * 18 + [M] * 3)


@lru_cache(maxsize=None)
def xf6_rep(g: str) -> np.ndarray:
    """O_h 元素 g 在 21 個座標上的表示 Γ(g)（Γ e_c = Σ_c' (g v_c)·v_c' e_c'；有快取，不要改回傳的陣列）。"""
    Mg = OH.num[g].real
    G = np.zeros((21, 21))
    for c, (_, a, v) in enumerate(XF6_COORDS_LIST):
        pa = Mg @ _atom_pos(a)
        gv = Mg @ np.array(v, float)
        for c2, (_, a2, v2) in enumerate(XF6_COORDS_LIST):
            if np.allclose(_atom_pos(a2), pa):
                G[c2, c] = gv @ np.array(v2, float)
    return G


def xf6_rows(k, j, b, book: bool = False) -> tuple:
    """(4.4.22) 的三列（⟨A|、⟨B|、⟨C1|；只有書上印的部分，其餘 None）。book=True 照印：⟨B|F|C1⟩ = −b/2、⟨C1|F|C1⟩ = 2(j+b)。"""
    rowA = [2 * k + j, 0, k / 2, k / 2, k / 2, k / 2, 0, 0, 0, -k / 2, -k / 2, 0, 0, k / 2, 0, 0, -k / 2, 0, 0, 0, -j]
    rowB = [None] * 6 + [k + b, 0, 0, 0, 0, 0, 0, -(k + b) / 2, 0, 0, -(k + b) / 2, 0, (-b / 2 if book else -b), 0, 0]
    rowC = [None] * 18 + [(2 * (j + b) if book else 2 * j + 4 * b), 0, 0]
    return rowA, rowB, rowC


def xf6_force_from_rows(k, j, b, book: bool = False) -> np.ndarray:
    """只用 (4.4.22) 的三列＋對稱 F = Γ(g) F Γ(g)ᵀ ＋ Hermitian 重建整個 21×21（(4.4.2) 的想法）。"""
    rowA, rowB, rowC = xf6_rows(k, j, b, book)
    F = np.full((21, 21), np.nan)
    F[0, :] = rowA
    for c, l in enumerate(A_LEADERS):
        F[c, :] = np.array(rowA) @ xf6_rep(l).T
    F[:, 0:6] = F[0:6, :].T
    rb = np.array([F[6, c] if x is None else x for c, x in enumerate(rowB)], float)
    for c, l in enumerate(B_LEADERS):
        F[6 + c, :] = rb @ xf6_rep(l).T
    F[:, 6:18] = F[6:18, :].T
    rc = np.array([F[18, c] if x is None else x for c, x in enumerate(rowC)], float)
    for c, l in enumerate(["1", "R3", "R2^3"]):
        F[18 + c, :] = rc @ xf6_rep(l).T
    return F


def xf6_translation(axis: int = 0) -> np.ndarray:
    """沿 x_{axis+1} 的剛體平移（21 維）：(4.4.20) 的 |t⟩。"""
    e = np.eye(3)[axis]
    return np.array([np.array(v, float) @ e for _, _, v in XF6_COORDS_LIST])


def xf6_numeric(k: float, j: float = 1.0, b: float = 0.0, m: float = 19, M: float = 240) -> np.ndarray:
    """21 個 ω²（由小到大；6 個 0）。"""
    Mi = np.diag(1 / np.sqrt(np.diag(xf6_mass(m, M))))
    w = np.linalg.eigvalsh(Mi @ xf6_force(k, j, b) @ Mi)
    return np.where(np.abs(w) < 1e-10, 0.0, w)


def xf6_secular(k, j=1, b=0, m=19, M=240) -> dict:
    """(4.4.34)：S、P、λ± = ω²（更正後的 (4.4.33)）。sympy 或浮點都可以。"""
    S = (3 * k + j + b) / m + (2 * j + 4 * b) / M
    P = (k * j + 2 * k * b + j * b) * (M + 6 * m) / (m ** 2 * M)
    if isinstance(S, sp.Basic):
        r = sp.sqrt(S ** 2 - 4 * P)
    else:
        r = math.sqrt(S * S - 4 * P)
    return dict(S=S, P=P, lp=(S + r) / 2, lm=(S - r) / 2)


def xf6_formulas(k, j=1, b=0, m=19, M=240) -> dict:
    """ω²：A1g (j+4k)/m、Eg (j+k)/m、T2g 2(k+b)/m、T2u (k+b)/m、T1u± ((4.4.23)–(4.4.26)、(4.4.34))。"""
    sec = xf6_secular(k, j, b, m, M)
    return {"A_{1g}": (j + 4 * k) / m, "E_g": (j + k) / m, "T_{2g}": 2 * (k + b) / m, "T_{2u}": (k + b) / m,
            "T_{1u}+": sec["lp"], "T_{1u}-": sec["lm"]}


XF6_DEG = {"A_{1g}": 1, "E_g": 2, "T_{2g}": 3, "T_{2u}": 3, "T_{1u}+": 3, "T_{1u}-": 3}


def xf6_t1u_matrix(k, j=1, b=0, m=19, M=240, book: bool = False):
    """(4.4.33) 的 2×2 加速度矩陣（基底 |c⟩、|c'⟩）；book=True 照印右上 2(j+b)(M+6m)/3m。"""
    up = 2 * ((j + b) if book else (j - b)) * (M + 6 * m) / (3 * m)
    return [[(9 * k + 2 * j + b) / (3 * m), up], [(j - b) / (3 * m * M), (j + 2 * b) * (M + 6 * m) / (3 * m * M)]]


def xf6_c_vectors(m=19, M=240) -> tuple:
    """(4.4.18)–(4.4.20) 的 x 夥伴（21 維）：|c⟩（軸向 F 1、赤道 F −½、中心 0）、|c'⟩（F 全部 M、中心 −6m）、|t⟩。"""
    c, cp = np.zeros(21), np.zeros(21)
    for n, (_, a, v) in enumerate(XF6_COORDS_LIST):
        vx = v[0]
        if a == 6:
            cp[n] = -6 * m * vx
        elif a in (2, 3):
            c[n] = vx
            cp[n] = M * vx
        else:
            c[n] = -0.5 * vx
            cp[n] = M * vx
    return c, cp, xf6_translation(0)


def xf6_t1u_reduced(k, j=1.0, b=0.0, m=19, M=240) -> tuple:
    """用 21 維的 |c⟩、|c'⟩ 直接算 ⟨F⟩、⟨m⟩ 的 2×2（(4.4.29)–(4.4.32)）與 ⟨m⁻¹F⟩。"""
    c, cp, _ = xf6_c_vectors(m, M)
    F, Mm = xf6_force(k, j, b), xf6_mass(m, M)
    V = np.stack([c, cp], 1)
    Fr, Mr = V.T @ F @ V, V.T @ Mm @ V
    return Fr, Mr, np.linalg.solve(Mr, Fr)


def xf6_t1u_vectors(k=0.18, j=1.0, b=0.2, m=19, M=32) -> dict:
    """(4.4.37)–(4.4.41)、Fig. 4.4.7：{'+'/'-': dict(lam, eps, eps_p, center, axial, equatorial)}（m 歸一化 (4.4.40)）。"""
    A = np.array(xf6_t1u_matrix(k, j, b, m, M), float)
    w, V = np.linalg.eig(A)
    order = np.argsort(-w.real)
    N = np.diag([3 * m, 6 * m * M * (M + 6 * m)])
    out = {}
    for tag, i in zip("+-", order):
        v = V[:, i].real
        v = v / math.sqrt(v @ N @ v)
        if v[0] < 0:
            v = -v
        out[tag] = dict(lam=w[i].real, eps=v[0], eps_p=v[1], center=v[1] * (-6 * m), axial=v[0] + v[1] * M,
                        equatorial=-0.5 * v[0] + v[1] * M)
    return out


def xf6_isotope(k=0.18, j=1.0, b=0.2, m=19, M1=32, M2=34) -> tuple:
    """(4.4.36a)：ω±(M1)/ω±(M2)。"""
    s1, s2 = xf6_secular(k, j, b, m, M1), xf6_secular(k, j, b, m, M2)
    return math.sqrt(s1["lp"] / s2["lp"]), math.sqrt(s1["lm"] / s2["lm"])


def xf6_scan(kjs=None, m=19, M=240, b=0.0) -> dict:
    """Fig. 4.4.5（UF6：m = 19、M = 240、b = 0）／4.4.6（SF6：M = 32、b = 0.2）：ω 以 √(j/m) 為單位，k/j ∈ [0, 1]。"""
    kjs = np.linspace(0, 1, 41) if kjs is None else np.asarray(kjs, float)
    out = {key: [] for key in ("kj",) + tuple(XF6_DEG)}
    for kj in kjs:
        f = xf6_formulas(kj, 1.0, b, m, M)
        out["kj"].append(kj)
        for key in XF6_DEG:
            out[key].append(math.sqrt(max(f[key], 0.0) * m))
    return {key: np.array(v) for key, v in out.items()}


def uf6_fit(m=19, M=240) -> dict:
    """用 UF6 的 ω(T1u+)/ω(A1g) = 624/667 定 k/j（b = 0）；再用 A1g 定 j 的尺度，預測其他頻率（cm⁻¹）。"""
    target = UF6_NU["T_{1u}+"] / UF6_NU["A_{1g}"]

    def ratio(kj):
        f = xf6_formulas(kj, 1.0, 0.0, m, M)
        return math.sqrt(f["T_{1u}+"] / f["A_{1g}"]) - target
    kj = brentq(ratio, 1e-4, 1.0)
    f = xf6_formulas(kj, 1.0, 0.0, m, M)
    scale = UF6_NU["A_{1g}"] / math.sqrt(f["A_{1g}"])
    return dict(kj=kj, pred={key: scale * math.sqrt(v) for key, v in f.items()})


def sf6_predict(kj=0.18, b=0.2, m=19, M=32) -> dict:
    """k/j = 0.18、b = 0.2 的 SF6：ω（單位 √(j/m)）與「用 A1g = 775 cm⁻¹ 定尺度」的預測。"""
    f = xf6_formulas(kj, 1.0, b, m, M)
    w = {key: math.sqrt(v * m) for key, v in f.items()}
    scale = SF6_NU["A_{1g}"] / w["A_{1g}"]
    return dict(omega=w, pred={key: scale * v for key, v in w.items()})


def xf6_isotypic(alpha: str) -> np.ndarray:
    """O_h irrep α 的同型投影 P^α = (ℓ/48) Σ χ^α(g) Γ(g)（21×21）。"""
    a = OH_IRREP_NAMES.index(alpha)
    ell = int(OH_CHARS[a, 0])
    P = np.zeros((21, 21))
    for g in OH.names:
        P += ell * float(OH_CHARS[a, OH.class_of(g)]) * xf6_rep(g)
    return P / 48


def xf6_classify() -> dict:
    """(4.4.3)–(4.4.9)：{軌道: O_h irrep 的頻率}（用 Γ 的跡）。"""
    out = {}
    for orb, idx in (("A", range(0, 6)), ("B", range(6, 18)), ("C", range(18, 21)), ("all", range(21))):
        idx = list(idx)
        tr = [float(np.trace(xf6_rep(C[0])[np.ix_(idx, idx)])) for C in OH.classes]
        out[orb] = [int(round(float(x))) for x in OH.frequencies([sp.nsimplify(round(t)) for t in tr], OH_CHARS)]
    return out


@lru_cache(maxsize=None)
def _irrep_num(alpha: str, g: str, basis: str = "tetragonal") -> np.ndarray:
    return np.array(irrep(alpha, g, basis).evalf(), dtype=complex)


def xf6_projected(alpha: str, j: int, k: int, orbit: str = "A") -> np.ndarray:
    """P^α_{jk}|A⟩ 或 |B⟩（(4.4.10)–(4.4.17)，四方基底的 O_h irreps），歸一化的 21 維向量。"""
    e = np.zeros(21)
    e[0 if orbit == "A" else 6] = 1
    v = np.zeros(21)
    for g in OH.names:
        v += (_irrep_num(alpha, g)[j, k].conjugate() * (xf6_rep(g) @ e)).real
    n = np.linalg.norm(v)
    return v / n if n > 1e-12 else v


def xf6_modes(k=0.18, j=1.0, b=0.2, m=19, M=32) -> list:
    """動畫用的簡正模（每種取第一個夥伴，座標依 XF6_NAMES）：[(標記, ω, 21 維位移)]；
    A1g、Eg（兩個夥伴）、T2g、T2u 由投影 (4.4.15)–(4.4.16)，T1u± 由 2×2 的本徵向量 (4.4.37)，再加 T1g（轉動）與 T1u（平移）。"""
    f = xf6_formulas(k, j, b, m, M)
    c, cp, t = xf6_c_vectors(m, M)
    tv = xf6_t1u_vectors(k, j, b, m, M)
    out = [("A_{1g}", math.sqrt(f["A_{1g}"]), xf6_projected("A_{1g}", 0, 0, "A")),
           ("E_g", math.sqrt(f["E_g"]), xf6_projected("E_g", 0, 0, "A")),
           ("E_g", math.sqrt(f["E_g"]), xf6_projected("E_g", 1, 0, "A")),
           ("T_{2g}", math.sqrt(f["T_{2g}"]), xf6_projected("T_{2g}", 0, 1, "B")),
           ("T_{2u}", math.sqrt(f["T_{2u}"]), xf6_projected("T_{2u}", 0, 0, "B")),
           ("T_{1u}+", math.sqrt(f["T_{1u}+"]), tv["+"]["eps"] * c + tv["+"]["eps_p"] * cp),
           ("T_{1u}-", math.sqrt(f["T_{1u}-"]), tv["-"]["eps"] * c + tv["-"]["eps_p"] * cp),
           ("T_{1g}", 0.0, xf6_projected("T_{1g}", 0, 1, "B")),
           ("T_{1u}(t)", 0.0, t / np.linalg.norm(t))]
    return out


# ======================================================================== §4.4C–D 正則座標與量子振動
def ladder(nmax: int) -> tuple:
    """截斷的 a、a†（(nmax+1)×(nmax+1)，基底 |0⟩…|nmax⟩）：a|n⟩ = √n|n−1⟩（(4.4.62)）。"""
    a = np.diag(np.sqrt(np.arange(1, nmax + 1)), 1)
    return a, a.T.copy()


def ladder_commutator(nmax: int) -> np.ndarray:
    """[a, a†] 在截斷空間：diag(1, …, 1, −nmax)（最後一格是截斷造成的）。"""
    a, ad = ladder(nmax)
    return a @ ad - ad @ a


def osc_levels(nmax: int, omega: float = 1.0) -> np.ndarray:
    """H = ħω(a†a + ½) 的本徵值 (4.4.52)、(4.4.63)。"""
    a, ad = ladder(nmax)
    return np.linalg.eigvalsh(omega * (ad @ a + 0.5 * np.eye(nmax + 1)))


def degeneracy(ell: int, n: int) -> int:
    """(4.4.65)：d(n) = (ℓ+n−1)!/(n!(ℓ−1)!)。"""
    return math.comb(ell + n - 1, n)


def degeneracy_count(ell: int, n: int) -> int:
    """數 ℓ 個非負整數和為 n 的組數（(4.4.64) 的另一種算法）。"""
    return sum(1 for t in itertools.product(range(n + 1), repeat=ell) if sum(t) == n)


def harmonic_states(ell: int, n: int) -> list:
    """ℓ 重簡併模的第 n 諧波的佔據數組（例：ℓ = 3、n = 2 → 6 組）。"""
    return [t for t in itertools.product(range(n, -1, -1), repeat=ell) if sum(t) == n]


PASCAL_464 = [[1, 1, 1, 1, 1], [1, 2, 3, 4, 5], [1, 3, 6, 10, 15], [1, 4, 10, 20, 35], [1, 5, 15, 35, 70]]   # (4.4.64)：列 n、行 ℓ=1..5
SF6_NU_FIG = {"nu1": (775, 1), "nu2": (644, 2), "nu3": (965, 3), "nu4": (617, 3), "nu5": (524, 3), "nu6": (363, 3)}  # Fig. 4.4.6 的值、ℓ


def state_count(Emax: float, modes: dict | None = None) -> int:
    """能量 ≤ Emax（cm⁻¹，零點能不計）的振動態數（含簡併，(4.4.66)）。"""
    modes = modes or SF6_NU_FIG
    vals = list(modes.values())

    def rec(i, E):
        if i == len(vals):
            return 1
        nu, ell = vals[i]
        tot, n = 0, 0
        while E + n * nu <= Emax + 1e-9:
            tot += degeneracy(ell, n) * rec(i + 1, E + n * nu)
            n += 1
        return tot
    return rec(0, 0.0)


def turning_point(n: int, ell: int, omega: float, hbar: float = 1.0) -> float:
    """(4.4.71)：q(TP) = [(2n + ℓ)ħ/ω]^{1/2}。"""
    return math.sqrt((2 * n + ell) * hbar / omega)


# ======================================================================== 習題
def _closure(G, gens: list) -> set:
    S = {G.e}
    frontier = [G.e]
    while frontier:
        new = []
        for x in frontier:
            for g in gens:
                y = G.table[g, x]
                if y not in S:
                    S.add(y)
                    new.append(y)
        frontier = new
    return S


def min_generators(G) -> int:
    """最少幾個生成元（習題 4.1.1；暴力搜尋 1、2、3 個）。"""
    for r in (1, 2, 3):
        for combo in itertools.combinations(G.names, r):
            if len(_closure(G, list(combo))) == G.order:
                return r
    return 4


def problem_411() -> dict:
    C3 = SN.FiniteGroup({g: O.mats[g] for g in C3T}, label="C3")
    D4g = SN.FiniteGroup({g: O.mats[g] for g in D4Z}, label="D4")
    return {"C3": min_generators(C3), "C3v": min_generators(SN.C3V), "D4": min_generators(D4g), "O": min_generators(O),
            "Oh": min_generators(OH)}


def problem_412() -> dict:
    return {"r1": O.normalizer("r1"), "R3": O.normalizer("R3")}


def subgroups_of_order(G, order: int) -> list:
    """G 中所有 order 階的子群（由兩個元素生成的；O 的 D2、D3、D4 都是）。"""
    found = []
    for a, b in itertools.combinations(G.names, 2):
        Sg = _closure(G, [a, b])
        if len(Sg) == order:
            Sg = sorted(Sg, key=G.idx.get)
            if Sg not in found:
                found.append(Sg)
    for a in G.names:
        Sg = sorted(_closure(G, [a]), key=G.idx.get)
        if len(Sg) == order and Sg not in found:
            found.append(Sg)
    return found


def problem_413() -> dict:
    """D2（四階、三個 180° 轉）、D3、D4 子群與是否正規。"""
    def is_d2(H):
        return len(H) == 4 and all(O.order_of(h) <= 2 for h in H)

    def is_dn(H, n):
        return len(H) == 2 * n and not SN.FiniteGroup({h: O.mats[h] for h in H}).is_abelian()

    d2 = [H for H in subgroups_of_order(O, 4) if is_d2(H)]
    d3 = [H for H in subgroups_of_order(O, 6) if is_dn(H, 3)]
    d4 = [H for H in subgroups_of_order(O, 8) if is_dn(H, 4)]
    return {"D2": [(H, O.is_normal(H)) for H in d2], "D3": [(H, O.is_normal(H)) for H in d3],
            "D4": [(H, O.is_normal(H)) for H in d4]}


TD_IRREP_NAMES = O_IRREP_NAMES


def _corr(G, tab, el, chars) -> sp.Matrix:
    return G.correlation(el, sp.Matrix(chars), tab)


def problem_421() -> dict:
    """習題 4.2.1 的相關表（列：大群的 irrep；行：子群的 irrep；標記見 PROBLEM_421_LABELS）。"""
    td_tab = O_CHARS                                       # T_d ≅ O（R → IR、i → Ii）
    C3vT = ["1", "r1", "r1^2", "I.i2", "I.i4", "I.i5"]
    CvT = ["1", "I.i4"]
    D2dT = ["1", "R3^2", "I.R3", "I.R3^3", "R1^2", "R2^2", "I.i3", "I.i4"]
    C2vT = ["1", "R3^2", "I.i3", "I.i4"]
    c3v = [[1] * 6, [1, 1, 1, -1, -1, -1], [2, -1, -1, 0, 0, 0]]
    d2d = [[1] * 8, [1, 1, -1, -1, 1, 1, -1, -1], [1, 1, 1, 1, -1, -1, -1, -1], [1, 1, -1, -1, -1, -1, 1, 1], [2, -2, 0, 0, 0, 0, 0, 0]]
    c2v = [[1, 1, 1, 1], [1, 1, -1, -1], [1, -1, 1, -1], [1, -1, -1, 1]]
    out = {"a": (_corr(TD, td_tab, C3vT, c3v), _corr(TD, td_tab, CvT, [[1, 1], [1, -1]])),
           "b": (_corr(TD, td_tab, D2dT, d2d), _corr(TD, td_tab, C2vT, c2v)),
           "c": (_corr(TD, td_tab, C3vT, c3v), _corr(TD, td_tab, C3T, sub_chars("C3").tolist())),
           "d": (correlation("D4h"), OH.correlation(D2H, sp.Matrix([[x * (1 if (k < 4 or p == 0) else -1)
                                                                       for k, x in enumerate(r + r)]
                                                                      for p in (0, 1) for r in sub_chars("D2").tolist()]), OH_CHARS)),
           "e": (correlation("C4v"), OH.correlation(C4Z, sub_chars("C4"), OH_CHARS)),
           "f": (OH.correlation(T_NAMES + ["I" if g == "1" else "I." + g for g in T_NAMES],
                                sp.Matrix([[TH_CHARS[a, TH.class_of(h)] for h in TH.names] for a in range(8)]), OH_CHARS),
                 TH.correlation(D2N, sub_chars("D2"), TH_CHARS))}
    return out


PROBLEM_421_LABELS = {"a": (["A_1", "A_2", "E"], ["A'", "A''"]), "b": (["A_1", "B_1", "A_2", "B_2", "E"], ["A_1", "A_2", "B_1", "B_2"]),
                      "c": (["A_1", "A_2", "E"], ["0_3", "1_3", "2_3"]), "d": (D4H_NAMES, [n + p + "}" for p in ("g", "u") for n in ("A_{1", "B_{1", "B_{2", "A_{2")]),
                      "e": (["A'", "B'", "A''", "B''", "E"], ["0_4", "1_4", "2_4", "3_4"]), "f": (TH_IRREP_NAMES, ["A_1", "B_1", "B_2", "A_2"])}
CUBE_VERTS = [(x, y, z) for x in (1, -1) for y in (1, -1) for z in (1, -1)]


def cage_symmetry(extra: list, field=None) -> list:
    """立方體八頂點（都有電荷）再加 extra 上的額外電荷（頂點的 ±1 座標），field = 磁場（軸向量）：
    回傳 O_h 中保持此電荷分佈（與磁場）的元素。"""
    ex = [np.array(p, float) for p in extra]
    out = []
    for g in OH.names:
        Mg = OH.num[g].real
        ok = all(any(np.allclose(Mg @ p, q) for q in ex) for p in ex)
        if ok and field is not None:
            B = np.array(field, float)
            ok = np.allclose(np.linalg.det(Mg) * (Mg @ B), B)
        if ok:
            out.append(g)
    return out


def split_t1u(el: list) -> list:
    """p 軌域（O_h 的 T1u）在子群 el 下分成幾個能階：回傳不可約成分的維度（由子群自己的類代數特徵標表）。"""
    Hg = SN.FiniteGroup({h: OH.mats[h] for h in el})
    ells, tab = Hg.char_table_class_algebra()
    chi = [OH_CHARS[8, OH.class_of(C[0])] for C in Hg.classes]
    f = Hg.frequencies(chi, tab)
    dims = []
    for a, n in enumerate(f):
        dims += [int(ells[a])] * int(n)
    return sorted(dims)


CAGES = {"a": ([], None), "b": ([(1, 1, 1)], None), "c": ([(1, 1, 1), (-1, -1, -1)], None),
         "d": ([(1, 1, 1), (-1, -1, 1)], None), "e": ([(1, 1, 1), (-1, 1, 1), (1, -1, 1), (-1, -1, 1)], None),
         "f": ([(1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1)], None),
         "g": ([(1, 1, 1), (-1, -1, -1), (1, -1, -1), (-1, 1, 1)], None), "h": ([(1, 1, 1)], (1, 1, 1))}
CAGE_GROUPS = {"a": "O_h", "b": "C_{3v}", "c": "D_{3d}", "d": "C_{2v}", "e": "C_{4v}", "f": "T_d", "g": "D_{2h}", "h": "C_3"}


def problem_423() -> dict:
    """{小題: (對稱群的階, 群名, T1u 分裂後的維度)}。"""
    out = {}
    for key, (extra, field) in CAGES.items():
        el = cage_symmetry(extra if extra else [], field) if extra else OH.names
        out[key] = (len(el), CAGE_GROUPS[key], split_t1u(el))
    return out


def problem_431() -> dict:
    """主 C3 誘導 0_3↑O（8 維）與 A1(D4)↑O（3 維）：{名: (在 1, r1, R1², R1, i1 的跡, 頻率)}。"""
    out = {}
    for name, sub in (("C3", C3T), ("D4", D4Z)):
        ind = induced(O, sub, {h: sp.Matrix([[1]]) for h in sub})
        tr = [ind[g].trace() for g in ("1", "r1", "R1^2", "R1", "i1")]
        f = O.frequencies([ind[C[0]].trace() for C in O.classes], O_CHARS)
        out[name] = (tr, f)
    return out


def problem_432() -> dict:
    """E(C3v)↑O_h（16 維）與 E_g(D4h)↑O_h（6 維）：{名: (在 1, r1, IR1², IR1, Ii1 的跡, 頻率)}。"""
    out = {}
    e3 = sub_irreps("C3v")["E"]
    ind = induced(OH, C3V, e3)
    out["C3v"] = ([ind[g].trace() for g in ("1", "r1", "I.R1^2", "I.R1", "I.i1")],
                  OH.frequencies([ind[C[0]].trace() for C in OH.classes], OH_CHARS))
    e4 = {h: O.mats[oh_split(h)[0]][:2, :2] for h in D4H}                 # E_g：反演不變號
    ind = induced(OH, D4H, e4)
    out["D4h"] = ([ind[g].trace() for g in ("1", "r1", "I.R1^2", "I.R1", "I.i1")],
                  OH.frequencies([ind[C[0]].trace() for C in OH.classes], OH_CHARS))
    return out


def _perm_rep_points(points: list, g: str) -> sp.Matrix:
    Mg = OH.num[g].real
    n = len(points)
    P = sp.zeros(n)
    for j, p in enumerate(points):
        v = Mg @ np.array(p, float)
        i = next(i for i, q in enumerate(points) if np.allclose(v, q))
        P[i, j] = 1
    return P


def _levels_by_irrep(points: list, dists: dict) -> dict:
    """點集上的穿隧模型 H − Σ s·A_d（A_d：距離 d 的鄰接矩陣）：{O_h irrep: 能量}（每個 irrep 只出現一次時 ＝ tr(P^αH)/tr(P^α)）。"""
    n = len(points)
    P = np.array(points, float)
    Hm = sp.eye(n) * Hs
    for dist, coef in dists.items():
        A = sp.Matrix(n, n, lambda i, j: 1 if (i != j and abs(np.linalg.norm(P[i] - P[j]) - dist) < 1e-9) else 0)
        Hm -= coef * A
    reps = {g: _perm_rep_points(points, g) for g in OH.names}
    out = {}
    for a, name in enumerate(OH_IRREP_NAMES):
        ell = OH_CHARS[a, 0]
        Pa = sum((ell * OH_CHARS[a, OH.class_of(g)] * reps[g] for g in OH.names), sp.zeros(n)) / 48
        tr = Pa.trace()
        if tr != 0:
            out[name] = (int(tr / ell), sp.simplify((Pa * Hm).trace() / tr))
    return out, Hm


def problem_433() -> dict:
    """立方體八井（局部 C3v）：{O_h irrep: (重數, 能量)}，最近鄰 −S（稜長 2）、次近鄰 −T（面對角 2√2）。"""
    return _levels_by_irrep(CUBE_VERTS, {2.0: Ss, 2 * math.sqrt(2): Ts})[0]


CUBOCT = [p for p in itertools.product((1, 0, -1), repeat=3) if sum(abs(x) for x in p) == 2]


def problem_434() -> dict:
    """立方八面體 12 井（局部 C2v）：最近鄰 −S（距離 √2）、次近鄰 −T（距離 2）。"""
    return _levels_by_irrep(CUBOCT, {math.sqrt(2): Ss, 2.0: Ts})[0]


GP, GM = (1 + sp.sqrt(5)) / 2, (1 - sp.sqrt(5)) / 2
Y_CHARS = sp.Matrix([[1, 1, 1, 1, 1], [3, GP, GM, 0, -1], [3, GM, GP, 0, -1], [4, -1, -1, 1, 0], [5, 0, 0, -1, 1]])   # 習題 4.3.5
Y_CLASS_SIZES = [1, 12, 12, 20, 15]
Y_IRREP_NAMES = ["A", "T_1", "T_3", "G", "H"]
D5_CHARS = sp.Matrix([[1, 1, 1, 1], [1, 1, 1, -1], [2, -GM, -GP, 0], [2, -GP, -GM, 0]])   # 類 0°, 72°, 144°, 180°（大小 1, 2, 2, 5）


def problem_435() -> dict:
    """Y ↓ D5、C5、D3、C3 的相關表（列 A, T1, T3, G, H）。"""
    def corr(sub_classes, sub_chars, sizes):
        n = sum(sizes)
        return sp.Matrix(5, sub_chars.rows, lambda a, b: sp.nsimplify(sp.simplify(
            sum(sizes[p] * Y_CHARS[a, sub_classes[p]] * sp.conjugate(sub_chars[b, p]) for p in range(len(sizes))) / n)))
    eps5 = [sp.exp(2 * sp.pi * IU * k / 5) for k in range(5)]
    c5 = sp.Matrix(5, 5, lambda k, p: sp.expand_complex(eps5[(k * p) % 5]))
    c3 = sp.Matrix([[1, 1, 1], [1, EPS3, sp.conjugate(EPS3)], [1, sp.conjugate(EPS3), EPS3]])
    d3 = sp.Matrix([[1, 1, 1], [1, 1, -1], [2, -1, 0]])
    return {"D5": corr([0, 1, 2, 4], D5_CHARS, [1, 2, 2, 5]), "C5": corr([0, 1, 2, 2, 1], c5, [1] * 5),
            "D3": corr([0, 3, 4], d3, [1, 2, 3]), "C3": corr([0, 3, 3], c3, [1, 1, 1])}


def problem_436(F: float, axis=(0, 0, 1), kind: str = "dipole", S: float = 1.0, T: float = 0.0) -> np.ndarray:
    """立方八井加外場的能階（H = 0、穿隧 −S、−T）：dipole 位能 −F (n·r)/|r|、quadrupole −F (n·r)²/|r|²。"""
    n = np.array(axis, float)
    n /= np.linalg.norm(n)
    P = np.array(CUBE_VERTS, float)
    Hm = np.zeros((8, 8))
    for i in range(8):
        for j in range(8):
            d = np.linalg.norm(P[i] - P[j])
            Hm[i, j] = -S if abs(d - 2) < 1e-9 else (-T if abs(d - 2 * math.sqrt(2)) < 1e-9 else 0.0)
        x = P[i] @ n / np.linalg.norm(P[i])
        Hm[i, i] += -F * x if kind == "dipole" else -F * x * x
    return np.linalg.eigvalsh(Hm)


AMU = 1.66053907e-27
ATOM_MASS = {"H": 1.008, "D": 2.014, "C": 12.011, "F": 18.998, "Si": 28.086}
XY4_NU = {"SiF4": ("Si", "F", (801, 264, 1032, 382)), "CF4": ("C", "F", (908, 435, 1283, 632)),
          "CD4": ("C", "D", (2069, 1092, 2259, 996)), "CH4": ("C", "H", (2917, 1534, 3019, 1306))}   # ν1 ν2 ν3 ν4（cm⁻¹，習題 4.4.1）


def xy4_omega2(j: float, k: float, mX: float, mY: float) -> dict:
    """四面體 XY4 的中心力模型（X–Y 鍵 j、Y–Y 稜 k）：{ν1 (A1), ν2 (E), ν3 (T2 高), ν4 (T2 低): ω²}（數值 15×15）。"""
    Y = [np.array(v, float) / math.sqrt(3) for v in R_AXES.values()]
    pos = Y + [np.zeros(3)]
    Fm = np.zeros((15, 15))

    def spring(a1, a2, kk):
        nvec = pos[a2] - pos[a1]
        nvec /= np.linalg.norm(nvec)
        g = np.zeros(15)
        g[3 * a2:3 * a2 + 3] += nvec
        g[3 * a1:3 * a1 + 3] -= nvec
        return kk * np.outer(g, g)
    for a in range(4):
        Fm += spring(4, a, j)
        for b in range(a + 1, 4):
            Fm += spring(a, b, k)
    Mi = np.diag(1 / np.sqrt([mY] * 12 + [mX] * 3))
    w = np.sort(np.linalg.eigvalsh(Mi @ Fm @ Mi))
    w = w[np.abs(w) > 1e-9 * max(1.0, w.max())]
    vals, cnt = [], []
    for x in w:
        if vals and abs(x - vals[-1]) < 1e-7 * max(1.0, abs(x)):
            cnt[-1] += 1
        else:
            vals.append(x)
            cnt.append(1)
    a1 = next(v for v, c in zip(vals, cnt) if c == 1)
    e = next(v for v, c in zip(vals, cnt) if c == 2)
    t2 = sorted(v for v, c in zip(vals, cnt) if c == 3)
    return {"nu1": a1, "nu2": e, "nu3": t2[1], "nu4": t2[0]}


def xy4_fit(mol: str = "CF4") -> dict:
    """習題 4.4.1：兩個彈簧常數（j、k，單位 N/m）用最小平方擬合四個基頻（log ω）；回傳 j、k、預測與誤差。"""
    from scipy.optimize import least_squares
    X, Y, nus = XY4_NU[mol]
    mX, mY = ATOM_MASS[X] * AMU, ATOM_MASS[Y] * AMU
    target = np.array(nus, float)

    def pred(p):
        j, k = p
        w2 = xy4_omega2(j, k, mX, mY)
        return np.array([math.sqrt(max(w2[n], 1e-30)) / (2 * math.pi * SN.C_LIGHT) for n in ("nu1", "nu2", "nu3", "nu4")])
    res = least_squares(lambda p: np.log(pred(np.exp(p))) - np.log(target), x0=np.log([500.0, 50.0]))
    j, k = np.exp(res.x)
    pr = pred((j, k))
    return dict(j=j, k=k, pred=pr, obs=target, rel=pr / target - 1)


def problem_442() -> dict:
    """立方 X8（每個頂點局部 C3v、3 個位移）與立方烷 C8H8 的模分類：{名: O_h irrep 的頻率（扣掉平移 T1u 與轉動 T1g）}。"""
    pts = CUBE_VERTS
    tr = []
    for C in OH.classes:
        g = C[0]
        Mg = OH.num[g].real
        P = np.array(_perm_rep_points(pts, g), dtype=float)
        tr.append(sp.nsimplify(round(float(np.trace(P)) * float(np.trace(Mg)))))
    f = [int(x) for x in OH.frequencies(tr, OH_CHARS)]
    vib8 = list(f)
    vib8[OH_IRREP_NAMES.index("T_{1u}")] -= 1
    vib8[OH_IRREP_NAMES.index("T_{1g}")] -= 1
    cub = [2 * x for x in f]
    cub[OH_IRREP_NAMES.index("T_{1u}")] -= 1
    cub[OH_IRREP_NAMES.index("T_{1g}")] -= 1
    return {"X8_all": f, "X8_vib": vib8, "C8H8_vib": cub}




def xf6_coset_omega2(k, j=1.0, b=0.0, m=19) -> dict:
    """陪集和 (4.3.46)、(4.4.23)–(4.4.29)：只用 ⟨A|F|·⟩、⟨B|F|·⟩ 兩列與四方 irreps 的一個分量求 ω²，以及 Q_AA、Q_AB、Q_BB。"""
    F = xf6_force(k, j, b)
    rA, rB = F[0], F[6]
    D = {a: {g: np.array(TETRA[a][g].evalf(), float) for g in O_NAMES} for a in ("E", "T_1", "T_2")}
    A = {l: c for c, l in enumerate(A_LEADERS)}
    B = {l: 6 + c for c, l in enumerate(B_LEADERS)}
    out = {"A_{1g}": sum(rA[A[l]] for l in A_LEADERS) / m, "E_g": sum(rA[A[l]] * D["E"][l][0, 0] for l in A_LEADERS) / m,
           "T_{2u}": sum(rB[B[l]] * D["T_2"][l][0, 0] for l in B_LEADERS) / m,
           "T_{2g}": sum(rB[B[l]] * D["T_2"][l][1, 1] for l in B_LEADERS) / m}
    out["Q_AA"] = sum(rA[A[l]] * D["T_1"][l][2, 2] for l in A_LEADERS)
    out["Q_AB"] = sum(rA[B[l]] * D["T_1"][l][2, 0] for l in B_LEADERS) / math.sqrt(2)
    out["Q_AB_book"] = sum(rA[A[l]] * D["T_1"][l][2, 0] for l in A_LEADERS) / math.sqrt(2)     # 書上印 ⟨A|Q|lA⟩
    out["Q_BB"] = sum(rB[B[l]] * D["T_1"][l][0, 0] for l in B_LEADERS)
    return out



def summary() -> dict:
    sec = xf6_secular(0.18, 1, 0.2, 19, 32)
    iso = xf6_isotope()
    tv = xf6_t1u_vectors()
    return {
        "O 特徵標表 (A1 A2 E T1 T2 × 1 r R² R i)": O_CHARS.tolist(),
        "c_i 最小方程": str(O.min_poly_class(4).as_expr()),
        "六點能階": {k: str(v) for k, v in SIX_LEVELS.items()},
        "0_4↑O 的跡 (1 r R² R i)": [six_induced("0_4")[C[0]].trace() for C in O.classes],
        "SF6 T1u λ± (j=1)": (_f4(sec["lp"]), _f4(sec["lm"])),
        "同位素比 ω±(32)/ω±(34)": (_f4(iso[0]), _f4(iso[1])),
        "T1u 向量 (+)": {k: _f4(v) for k, v in tv["+"].items()},
        "T1u 向量 (−)": {k: _f4(v) for k, v in tv["-"].items()},
        "UF6 擬合": {k: (_f4(v) if not isinstance(v, dict) else {kk: f"{vv:.1f}" for kk, vv in v.items()}) for k, v in uf6_fit().items()},
        "SF6 預測": {kk: f"{vv:.1f}" for kk, vv in sf6_predict()["pred"].items()},
    }


# ======================================================================== 自我測試
def _numrep_ok(D: dict, G=None) -> bool:
    G = G or O
    N = {g: np.array(D[g].evalf(), dtype=complex) for g in G.names}
    return all(np.allclose(N[a] @ N[b], N[G.table[a, b]]) for a in G.names for b in G.names)


def _elementary_multi(Ps: list) -> dict:
    """基本算子（(3.3.11)）：P_{1j} = P_1 g_j P_j / √λ_j（g_j 是第一個使 P_1 g P_j ≠ 0 的元素，N N† = λ P_1），
    P_{j1} = P_{1j}†、P_{ij} = P_{i1}P_{1j}。單項式（每列一個非零元）的表示需要每個 j 各自的參考元素。"""
    ell = len(Ps)
    P1j = {0: Ps[0]}
    for j in range(1, ell):
        for g in O_NAMES:
            N = O.ga_mul(O.ga_mul(Ps[0], O.el(g)), Ps[j])
            if N:
                break
        NN = O.ga_mul(N, O.ga_dag(N))
        g0 = next(iter(Ps[0]))
        lam = sp.nsimplify(sp.simplify(NN[g0] / Ps[0][g0]))
        P1j[j] = O.ga_scale(N, 1 / sp.sqrt(lam))
    Pi1 = {j: O.ga_dag(P1j[j]) for j in range(ell)}
    Pi1[0] = Ps[0]
    return {(i, j): (Ps[0] if i == j == 0 else O.ga_mul(Pi1[i], P1j[j])) for i in range(ell) for j in range(ell)}


def _check_group():
    assert O.order == 24 and OH.order == 48 and len(set(O_NAMES)) == 24
    for a, b, c in O_PRODUCTS_BOOK:                                          # 乘法表 ＝ Hamilton 弧（四元數）
        assert O.mul(a, b) == c and nomogram_product(a, b)["name"] == c
    hp = nomogram_product("r3", "R2")
    assert abs(hp["deg"] - 90) < 1e-9 and np.allclose(np.abs(hp["axis"]), [1, 0, 0])
    for l, C in COSETS_439_BOOK.items():                                     # (4.3.9)
        assert [O.mul(l, h) for h in C4Z] == C
    perms = {g: cycle_str(diag_perm(g)) for g in O_NAMES}                    # S4：四條體對角線
    assert perms == O_PERM_BOOK and len({tuple(sorted(diag_perm(g).items())) for g in O_NAMES}) == 24
    assert all(diag_perm(O.mul(a, b)) == {k: diag_perm(a)[diag_perm(b)[k]] for k in range(1, 5)} for a in O_NAMES for b in O_NAMES)
    assert [len(C) for C in O.classes] == O_CLASS_SIZES and [O.order_of(C[0]) for C in O.classes] == [1, 3, 2, 4, 2]
    assert O.classes[1] == ["r1", "r2", "r3", "r4", "r1^2", "r2^2", "r3^2", "r4^2"]
    for g in O_NAMES:                                                        # 軸與角 ＝ 矩陣
        ax, deg, half = o_turn(g)
        if g != "1":
            a2, d2 = SN.axis_angle(O.num[g].real)
            assert abs(d2 - min(deg, 360 - deg)) < 1e-9 and half == deg / 2
    # 類代數 (4.1.6)：結構常數 ＝ 照抄的表；元素個數
    for (i, j), v in O_CLASS_ALG_BOOK.items():
        assert list(O.struct[i, j]) == v and list(O.struct[j, i]) == v
        assert sum(x * s for x, s in zip(v, O_CLASS_SIZES)) == O_CLASS_SIZES[i] * O_CLASS_SIZES[j]
    assert ci_powers() == CI_POWERS                                          # (4.1.7)（正確版）
    assert CI5_BOOK == [20 * x for x in CI_POWERS[4][:3]] + [320, 256]       # 書上把 20c_i^4 代進去
    assert CI_POWERS[5] == [0, 0, 0, 320 + 20 * 16, 256 + 20 * 20]           # 應為 320c_R + 20c_i^3 + 256c_i = 640c_R + 656c_i
    c = sp.Symbol("c")
    assert O.min_poly_class(4).as_expr() == c ** 5 - 40 * c ** 3 + 144 * c   # (4.1.8)
    assert sp.factor(c ** 5 - 40 * c ** 3 + 144 * c) == c * (c - 6) * (c - 2) * (c + 2) * (c + 6)
    ev = sorted(np.linalg.eigvals(np.array(O.class_matrix(4), float)).real)  # 第二種：類代數的正則表示的本徵值
    assert np.allclose(ev, [-6, -2, 0, 2, 6])
    assert [a - 40 * b + 144 * d for a, b, d in zip(CI_POWERS[5], CI_POWERS[3], CI_POWERS[1])] == [0] * 5
    # (4.1.9)：ℙ^{T2} = (c_i+2)c_i(c_i²−36)/(−256)
    num = [CI_POWERS[4][k] + 2 * CI_POWERS[3][k] - 36 * CI_POWERS[2][k] - 72 * CI_POWERS[1][k] for k in range(5)]
    assert [R_(x, -256) for x in num] == O_CLASS_IDEMP_BOOK["T_2"]
    assert o_class_idempotents() == O_CLASS_IDEMP_BOOK                       # (4.1.10)
    assert o_char_table_class_algebra() == O_CHARS                           # 特徵標：類代數 ＝ (4.1.11)
    assert O.char_table_traces([TETRA[a] for a in O_IRREP_NAMES]) == O_CHARS  # ＝ 跡
    for a, name in enumerate(O_IRREP_NAMES):                                 # ℙ^α 的係數 ∝ χ^α
        assert [x * 24 / O_DIMS[a] for x in O_CLASS_IDEMP_BOOK[name]] == list(O_CHARS.row(a))
        assert [O_CLASS_SIZES[k] * O_CHARS[a, k] / O_DIMS[a] for k in range(5)] == O_CLASS_EIG[a]
    assert [O_CLASS_EIG[a][4] for a in range(5)] == [6, -6, 0, -2, 2] and sum(d * d for d in O_DIMS) == 24
    # O_h、T_d、T、T_h
    assert [len(C) for C in OH.classes] == O_CLASS_SIZES * 2
    assert OH.char_table_traces([{g: irrep(n, g) for g in OH.names} for n in OH_IRREP_NAMES]) == OH_CHARS   # (4.1.16)
    assert OH_CHARS.row(8).tolist()[0] == [3, 0, -1, 1, -1, -3, 0, 1, -1, 1]
    assert O.iso_to(TD, O_TO_TD) and SN.book_order(TD, O_CHARS_BOOK)[1] == O_CHARS
    assert TGRP.order == 12 and [len(C) for C in TGRP.classes] == [1, 4, 4, 3]
    assert TGRP.class_of("r1") != TGRP.class_of("r1^2") and O.class_of("r1") == O.class_of("r1^2")
    Tt = TGRP.char_table_class_algebra()[1]
    assert sorted(SN.reorder(Tt, T_CHARS.tolist())) == [0, 1, 2, 3]
    for a in (1, 2):                                                         # 一維複數表示是同態
        chi = {g: T_CHARS[a, TGRP.class_of(g)] for g in T_NAMES}
        assert all(sp.simplify(chi[x] * chi[y] - chi[TGRP.table[x, y]]) == 0 for x in T_NAMES for y in T_NAMES)
    assert TH.order == 24 and TH_CHARS.shape == (8, 8)
    assert _f4(GOLDEN) == "1.6180"


def _check_irreps():
    for basis, D in BASES.items():
        for a in O_IRREP_NAMES:
            assert _numrep_ok(D[a]), (basis, a)
            assert [sp.nsimplify(sp.simplify(D[a][C[0]].trace())) for C in O.classes] == list(O_CHARS.row(O_IRREP_NAMES.index(a)))
            U = {g: np.array(D[a][g].evalf(), complex) for g in O_NAMES}
            assert all(np.allclose(U[g] @ U[g].conj().T, np.eye(U[g].shape[0])) for g in O_NAMES)
    for g, code in T1_BOOK_REPS.items():                                    # (4.2.14)
        assert TETRA["T_1"][g] == signed_perm(code)
    for g, code in T2_BOOK.items():                                          # (4.2.16) 全部 24 個
        assert TETRA["T_2"][g] == signed_perm(code)
    for g in O_NAMES:                                                        # (4.2.15) 版本只差第二基底變號
        alt = (1 if O.class_of(g) < 3 else -1) * O.mats[g]
        assert S_KOSHER * alt * S_KOSHER == TETRA["T_2"][g]
    for g, M in E_BOOK_REPS.items():                                         # (4.2.19)
        assert TETRA["E"][g] == sp.Matrix(M)
    for C in E_COSET_CLASSES:
        assert all(TETRA["E"][g] == TETRA["E"][C[0]] for g in C)
    for g in D4Z:                                                            # (4.2.13)：在 D4 塊對角、在 D2 對角
        for a in ("T_1", "T_2"):
            M = TETRA[a][g]
            assert M[0, 2] == M[1, 2] == M[2, 0] == M[2, 1] == 0
        assert TETRA["E"][g].is_diagonal()
    for g in D2N:
        assert all(TETRA[a][g].is_diagonal() for a in ("E", "T_1", "T_2"))
    # 子群鏈分裂 ＝ 由矩陣的 (ℓ/24)Σ𝒟_jj* g
    for a, ell in (("T_1", 3), ("T_2", 3), ("E", 2)):
        sp_ = chain_split(a)
        Pij = elementary_from_irrep(TETRA[a], ell)
        assert len(sp_) == ell and all(O.ga_eq(sp_[j], Pij[j, j]) for j in range(ell)), a
        assert all(O.ga_eq(O.ga_mul(Pij[j, j], Pij[j, j]), Pij[j, j]) for j in range(ell))
    sT1 = chain_split("T_1")
    for j in range(3):                                                       # (4.2.10)
        assert O.ga_eq(sT1[j], {g: R_(c, 8) for g, c in P410_BOOK[j].items()})
    sE = chain_split("E")                                                    # (4.2.18b,c)：前置因子 2/3
    assert O.ga_eq(sE[0], p418("A1")) and O.ga_eq(sE[1], p418("B1"))
    bad = p418("A1", R_(3, 4))
    assert not O.ga_eq(O.ga_mul(bad, bad), bad) and bad["1"] == R_(3, 32) and sE[0]["1"] == R_(1, 12)
    Pd, Rd = d4_idempotents(), O.split(o_central("T_1"), [d4_idempotents()[k] for k in ["E1", "E2", "A1", "A2", "B1", "B2"]])
    assert len(Rd) == 3 and not O.ga_mul(o_central("T_1"), Pd["A1"])        # (4.2.9)：三項
    # SN 的第二階段（基本算子 P_{1j} = P_1 g_j P_j/√λ，每個 j 各取一個參考元素）＝ 四方矩陣（差一個對角的 ±1）
    Dsn = O.irrep_from_elementary(_elementary_multi(sT1))
    sg = [1] + [next(sp.sign(Dsn[g][0, j] / TETRA["T_1"][g][0, j]) for g in O_NAMES if TETRA["T_1"][g][0, j] != 0) for j in (1, 2)]
    Ph = sp.diag(*sg)
    assert all(Dsn[g] == Ph * TETRA["T_1"][g] * Ph for g in O_NAMES)
    # 三方：投影出來的基底 ＝ (4.2.26)／(4.2.31)；類代表 ＝ (4.2.28)；三個印錯的元素
    assert sp.simplify(TRIG_U_T1 - T43) == sp.zeros(3)
    for g, M in TRIG_T1_BOOK.items():
        assert sp.simplify(TRIG["T_1"][g] - sp.Matrix(M)) == sp.zeros(3)
    for (g, i, j), (book, right) in TRIG_T1_BOOK_TYPOS.items():
        assert sp.simplify(TRIG["T_1"][g][i, j] - right) == 0 and sp.simplify(TRIG["T_1"][g][i, j] - book) != 0
    assert sp.simplify(TRIG["T_1"]["r2"] - TRIG["T_1"]["r2^2"].T) == sp.zeros(3)
    for h in D3T:
        for a in ("T_1", "T_2"):
            M = TRIG[a][h]
            blk = (2, (0, 1)) if a == "T_1" else (0, (1, 2))
            assert all(M[blk[0], q] == 0 and M[q, blk[0]] == 0 for q in blk[1])
        assert sp.simplify(TRIG["T_2"][h][1:, 1:] - TRIG["T_1"][h][:2, :2]) == sp.zeros(2)
    assert all(TRIG[a]["i4"].is_diagonal() for a in ("E", "T_1", "T_2"))
    for a, ell in (("T_1", 3), ("T_2", 3)):                                  # 三方鏈的分裂 ＝ 三方矩陣的投影
        Pij = elementary_from_irrep(TRIG[a], ell)
        sp3 = chain_split(a, "trigonal")
        assert len(sp3) == 3 and all(any(O.ga_eq(x, Pij[j, j]) for x in sp3) for j in range(3))
    # 行進波：(4.2.35)；(4.2.36) 的對角；C4、C3 對角化
    for g, M in MOV4_T1_BOOK.items():
        assert sp.simplify(MOV4["T_1"][g] - sp.Matrix(M)) == sp.zeros(3), g
    assert MOV4["T_1"]["R3"] == sp.diag(-IU, IU, 1) and MOV4["T_1"]["R3"] != D436_BOOK and MOV4["T_1"]["R3^3"] == D436_BOOK
    for g in C4Z:
        assert all(MOV4[a][g].is_diagonal() for a in ("T_1", "T_2", "E"))
    assert MOV4["T_2"]["R3"] == sp.diag(-IU, IU, -1)
    for g in C3T:
        assert all(sp.simplify(MOV3[a][g] - sp.diag(*[MOV3[a][g][q, q] for q in range(MOV3[a][g].rows)])) == sp.zeros(MOV3[a][g].rows)
                   for a in ("T_1", "T_2", "E"))
    e1 = sp.conjugate(EPS3)
    assert sp.simplify(MOV3["T_1"]["r1"] - sp.diag(e1, EPS3, 1)) == sp.zeros(3)
    assert sp.simplify(MOV3["T_2"]["r1"] - sp.diag(e1, EPS3, 1)) == sp.zeros(3)
    assert sp.simplify(MOV3["E"]["r1"] - sp.diag(e1, EPS3)) == sp.zeros(2)
    assert sp.simplify(MOV4["T_1"]["R3"] - TMOV4.H * TETRA["T_1"]["R3"] * TMOV4) == sp.zeros(3)
    # (4.2.20)–(4.2.22)、(4.2.29)–(4.2.30)：T1u 的四方與三方
    assert irrep("T_{1u}", "I.r1")[1, 0] == -1 and irrep("T_{1u}", "I.R1^2") == sp.diag(-1, 1, 1)
    assert irrep("T_{1u}", "I.i4") == sp.Matrix([[0, 1, 0], [1, 0, 0], [0, 0, 1]])
    assert irrep("T_{1u}", "I.i4", "trigonal") == sp.diag(-1, 1, 1)
    assert sp.simplify(irrep("T_{1u}", "I.R1^2", "trigonal") + TRIG_T1_BOOK_R12()) == sp.zeros(3)


def TRIG_T1_BOOK_R12():
    return sp.Matrix(TRIG_T1_BOOK["R1^2"])


def _check_corr_frob():
    for key, tab in CORR_BOOK.items():
        assert correlation(key) == sp.Matrix(tab), key
    d4h = correlation("D4h")                                                 # (4.3.53)
    assert list(d4h.row(0)) == [1] + [0] * 9 and list(d4h.row(8)) == [0, 0, 0, 0, 0, 0, 1, 0, 0, 1]
    assert list(d4h.row(2)) == [1, 0, 1, 0, 0, 0, 0, 0, 0, 0]
    d3d = correlation("D3d")                                                 # (4.3.86)
    assert list(d3d.row(0)) == [1, 0, 0, 0, 0, 0] and list(d3d.row(8)) == [0, 0, 0, 0, 1, 1] and list(d3d.row(2)) == [0, 0, 1, 0, 0, 0]
    c4 = correlation("C4")                                                   # (4.3.73)、(4.2.43) 的 T1 列
    assert list(c4.row(3)) == [1, 1, 0, 1] and [O_CHARS[3, O.class_of(h)] for h in C4Z] == [3, 1, -1, 1]
    assert O.is_normal(D2N) and not O.is_normal(D2P)
    assert sorted(conjugates(D2P)) == sorted([D2P, ["1", "R2^2", "i1", "i2"], ["1", "R1^2", "i5", "i6"]])
    assert conjugates(D2N) == [D2N]
    for key in ("C4", "D4", "D2", "D2'", "D3", "C3", "C2", "C2z", "C1", "T", "C2v", "C3v", "C4v"):
        G = O if SUB[key][0] == "O" else OH
        for lab, (fm, ff, fc) in frobenius(key).items():
            assert fm == ff == fc, (key, lab)
            d = sub_irreps(key)[lab][SUB[key][1][0]].shape[0]
            tab = O_CHARS if G is O else OH_CHARS
            assert sum(f * tab[a, 0] for a, f in enumerate(fm)) == G.order // len(SUB[key][1]) * d
    fr = frobenius("C4v")
    assert fr["A'"][0] == [1, 0, 1, 0, 0, 0, 0, 0, 1, 0] and fr["E"][0] == [0, 0, 0, 1, 1, 0, 0, 0, 1, 1]   # (4.3.47)、(4.3.48)
    ex = d6_induced_examples()                                               # (4.3.49)–(4.3.51)、正則表示
    nm = SN.D6_IRREP_NAMES
    pick = lambda f: sorted(n for n, x in zip(nm, f) for _ in range(int(x)))   # noqa: E731
    assert pick(ex["0_2"]) == sorted(["A_1", "E_1", "E_2", "B_1"]) and pick(ex["1_2"]) == sorted(["B_2", "E_2", "E_1", "A_2"])
    assert pick(ex["0_6"]) == ["A_1", "A_2"] and pick(ex["3_6"]) == ["B_1", "B_2"]
    assert pick(ex["1_6"]) == ["E_1"] and pick(ex["2_6"]) == ["E_2"]
    assert ex["regular_O"] == [1, 1, 2, 3, 3]


def _six_geometric_accel(k: float, t: float) -> np.ndarray:
    """(4.3.2b) 的第二種算法：六個質點沿軸滑動，最近鄰彈簧 k、對面彈簧 t 的力矩陣（V = ½xᵀFx）。"""
    P = [np.array(s, float) for s in SITES]
    Fm = np.zeros((6, 6))
    for a in range(6):
        for b in range(a + 1, 6):
            nvec = P[b] - P[a]
            nvec /= np.linalg.norm(nvec)
            kk = t if np.allclose(P[a], -P[b]) else k
            g = np.zeros(6)
            g[b] += nvec @ P[b]
            g[a] -= nvec @ P[a]
            Fm += kk * np.outer(g, g)
    return Fm


def _check_six():
    ind = six_induced("0_4")
    assert all(ind[g] == site_perm(g) for g in O_NAMES)                    # 誘導表示 ＝ 頂點的置換
    assert [ind[C[0]].trace() for C in O.classes] == [6, 0, 2, 2, 0]
    assert O.frequencies([6, 0, 2, 2, 0], O_CHARS) == [1, 0, 1, 1, 0]       # A1 ⊕ E ⊕ T1
    assert ind["i4"] == SIX_INDUCED_I4_BOOK                                  # (4.3.12)
    assert [site_of(g) for g in SITE_LEADERS] == list(range(6))
    ind1 = six_induced("1_4")
    assert ind1["i4"] == SIX_INDUCED_1_4_I4_BOOK                             # (4.3.34ₓ)
    assert O.frequencies([ind1[C[0]].trace() for C in O.classes], O_CHARS) == [0, 0, 0, 1, 1]   # (4.3.40)：T1 ⊕ T2
    P14 = sum((sp.conjugate(c) * MOV4["T_1"][h] for h, c in zip(C4Z, [1, -IU, -1, IU])), sp.zeros(3)) / 4
    assert P14 == sp.diag(1, 0, 0)                                           # (4.3.39)
    Hm = six_H()
    ev = Hm.eigenvals()
    assert {sp.simplify(k): v for k, v in ev.items()} == {sp.simplify(SIX_LEVELS["A_1"]): 1, sp.simplify(SIX_LEVELS["T_1"]): 3,
                                                          sp.simplify(SIX_LEVELS["E"]): 2}
    for (a, lab), vb in SIX_VECTORS_BOOK.items():
        if lab in ("x", "y"):
            continue
        j = MOV4_LABELS[a].index(lab)
        v1, v2 = six_projected(a, j), six_coset_vector(a, j)
        assert sp.simplify(v1 - vb) == sp.zeros(6, 1) and sp.simplify(v2 - vb) == sp.zeros(6, 1), (a, lab)
        assert sp.simplify(Hm * vb - SIX_LEVELS[a] * vb) == sp.zeros(6, 1)
        assert sp.simplify(six_coset_energy(a) - SIX_LEVELS[a]) == 0
    x = (-SIX_VECTORS_BOOK[("T_1", "1_4")] + SIX_VECTORS_BOOK[("T_1", "3_4")]) / SQ2      # (4.3.27)
    y = IU * (SIX_VECTORS_BOOK[("T_1", "1_4")] + SIX_VECTORS_BOOK[("T_1", "3_4")]) / SQ2
    assert sp.simplify(x - SIX_VECTORS_BOOK[("T_1", "x")]) == sp.zeros(6, 1) and sp.simplify(y - SIX_VECTORS_BOOK[("T_1", "y")]) == sp.zeros(6, 1)
    assert sp.simplify(SIX_E_BOOK_446X - SIX_LEVELS["E"]) != 0
    assert np.allclose(six_numeric(0, 0, -1), [-4, 0, 0, 0, 2, 2])
    assert sp.simplify(SIX_LEVELS["A_1"] + 3 * SIX_LEVELS["T_1"] + 2 * SIX_LEVELS["E"] - 6 * Hs) == 0
    for kk, tt in ((1.0, 0.3), (2.0, 0.0)):                                  # (4.3.2b)：h = 2k + t、s = k/2
        assert np.allclose(_six_geometric_accel(kk, tt), np.array(six_accel(kk, tt), float))
    assert {sp.simplify(k_): v for k_, v in six_accel().eigenvals().items()} == {
        4 * sp.Symbol("k") + 2 * sp.Symbol("t"): 1, 2 * sp.Symbol("k"): 3, sp.Symbol("k") + 2 * sp.Symbol("t"): 2}


def _check_pert():
    assert sp.simplify(vq_reduced() - VQ_REDUCED_BOOK) == sp.zeros(6)     # (4.3.56)
    blk = VQ_REDUCED_BOOK[3:5, 3:5]
    assert sp.simplify(VQ_HIGHQ.T * blk * VQ_HIGHQ - VQ_460_BOOK) == sp.zeros(2)    # (4.3.60)
    v58 = VQ_BASIS[:, 3:5] * VQ_HIGHQ[:, 0]
    v59 = VQ_BASIS[:, 3:5] * VQ_HIGHQ[:, 1]
    assert sp.simplify(v58 - sp.Matrix([1, 1, 0, 0, 0, 0]) / SQ2) == sp.zeros(6, 1)   # (4.3.58b)
    assert sp.simplify(v59 - sp.Matrix([0, 0, 1, 1, 1, 1]) / 2) == sp.zeros(6, 1)     # (4.3.59b)
    for Q in (-6.0, -2.0, 0.5, 3.0):
        M = np.array(vq_matrix(Q, 0, -1), float)
        assert np.allclose(np.linalg.eigvalsh(M), [e for e, _ in vq_levels(Q)])
    h = 1e-6                                                                 # 一階斜率
    sl = {"A_{2u}": 1, "E_u": 0, "B_{1g}": 0}
    sc = vq_scan([-h, h])
    assert abs((sc["A1g-"][1] - sc["A1g-"][0]) / (2 * h) - 1 / 3) < 1e-6 and abs((sc["A1g+"][1] - sc["A1g+"][0]) / (2 * h) - 2 / 3) < 1e-6
    assert abs((sc["A2u"][1] - sc["A2u"][0]) / (2 * h) - sl["A_{2u}"]) < 1e-9
    Qm, gap = vq_avoided()
    assert Qm == 2 * Ss and sp.simplify(gap.subs(Ss, -1) - 4 * SQ2) == 0 and sp.simplify(gap.subs(Ss, R_(1, 2)) - 2 * SQ2) == 0
    qs = np.linspace(-4, 0, 4001)
    g = vq_scan(qs)
    i = int(np.argmin(g["A1g+"] - g["A1g-"]))
    assert abs(qs[i] + 2) < 2e-3 and _f4((g["A1g+"] - g["A1g-"])[i]) == "5.6569"
    # V_D
    red = vd_reduced()
    assert sp.simplify(red - VD_REDUCED) == sp.zeros(6)                     # (4.3.66) 更正後
    assert red[3, 3] != Hs - 4 * Ss
    assert (VD_T1U_BOOK.T * VD_BASIS[:, 3])[0] != 0                          # 照印的 T1u 與 A1g 不正交
    assert sp.simplify((VD_T1U_BOOK.T * sp.diag(-Ds, Ds, 0, 0, 0, 0) * VD_BASIS[:, 3])[0]) == 0   # 照印算不出 −D/√3
    blk3 = VD_REDUCED[2:5, 2:5]
    assert sp.simplify(VD_469 * blk3 * VD_469.T - VD_469_RESULT) == sp.zeros(3)  # (4.3.69)
    B3 = VD_BASIS[:, 2:5]
    e1 = B3 * VD_469[0, :].T                                                 # (4.3.70)：|1⟩
    e2 = B3 * VD_469[2, :].T                                                 # (4.3.71) 更正：係數 (−1/√2, 1/√6, 1/√3)
    e2_book = B3 * sp.Matrix([1 / SQ2, 1 / SQ6, 1 / SQ3])
    assert sp.simplify(e1 - sp.Matrix([1, 0, 0, 0, 0, 0])) == sp.zeros(6, 1)
    assert sp.simplify(e2 - sp.Matrix([0, 1, 0, 0, 0, 0])) == sp.zeros(6, 1) and sp.simplify(e2_book - e1) == sp.zeros(6, 1)
    eq = VD_BASIS[:, 3:5] * sp.Matrix([2 / SQ6, -1 / SQ3])                   # (4.3.72)
    assert sp.simplify(eq - sp.Matrix([0, 0, 1, 1, 1, 1]) / 2) == sp.zeros(6, 1)
    assert sp.simplify((eq.T * vd_matrix() * eq)[0] - (Hs + 2 * Ss)) == 0
    for Dv in (-5.0, 0.7, 4.0):
        sc = vd_scan([Dv])
        lv = sorted([sc[k][0] for k in ("A'1", "A'2", "A'3", "E", "B'")] + [sc["E"][0]])
        assert np.allclose(np.linalg.eigvalsh(np.array(vd_matrix(Dv, 0, -1), float)), lv)
    h = 1e-3                                                                 # 二階係數（S = −1）
    s0, s1 = vd_scan([0.0]), vd_scan([h])
    c2 = {k: (s1[k][0] - s0[k][0]) / h ** 2 for k in ("A'1", "A'2", "A'3")}
    assert abs(c2["A'1"] + 1 / 12) < 1e-5 and abs(c2["A'2"] + 1 / 4) < 1e-5 and abs(c2["A'3"] - 1 / 3) < 1e-5
    assert {k: v.subs(Ss, -1) for k, v in VD_SECOND_ORDER.items()} == {"A_{1g}": R_(-1, 12), "T_{1u}": R_(-1, 4), "E_g": R_(1, 3)}
    # 磁場 C4：數值 ＝ 公式；保有 C4、破壞 σ_v；一階 Zeeman 斜率 ±2S
    for beta in (0.0, 0.3, 1.1):
        M = vb_matrix(beta)
        assert np.allclose(np.linalg.eigvalsh(M), sorted(vb_levels(beta).values()))
        R3 = np.array(site_perm("R3"), float)
        assert np.allclose(R3 @ M @ R3.T, M)
    sv = np.array(site_perm("I.R1^2"), float)
    assert not np.allclose(sv @ vb_matrix(0.3) @ sv.T, vb_matrix(0.3))
    lv0 = vb_levels(0.0)
    assert np.isclose(lv0["A_1(0_4)"], -4) and np.isclose(lv0["E(0_4)"], 2) and np.isclose(lv0["2_4"], 2) and np.isclose(lv0["1_4"], 0)
    h = 1e-6
    assert abs((vb_levels(h)["1_4"] - vb_levels(-h)["1_4"]) / (2 * h) + 2) < 1e-6
    # (111) 電場
    Mv = v111_matrix()
    B = V111_BASIS
    assert sp.simplify(sp.Matrix([[(B[a].T * Mv * B[b])[0] for b in ("A1g", "T1uA")] for a in ("A1g", "T1uA")]) - V111_A) == sp.zeros(2)
    for p, q in (("Eg1", "T1uE1"), ("Eg2", "T1uE2")):
        assert sp.simplify(sp.Matrix([[(B[a].T * Mv * B[b])[0] for b in (p, q)] for a in (p, q)]) - V111_E) == sp.zeros(2)
    assert sp.simplify((B["A1g"].T * Mv * B["Eg1"])[0]) == 0 and sp.simplify((B["T1uA"].T * Mv * B["T1uE1"])[0]) == 0
    for dv in (0.5, 2.0):
        assert np.allclose(np.linalg.eigvalsh(np.array(v111_matrix(dv, 0, -1), float)), [e for e, _ in v111_levels(dv)])
    t1 = [TETRA["T_1"]["1"][:, q] for q in range(3)]
    assert sp.simplify(T485_BOOK[:, 2] - (t1[0] + t1[1] + t1[2]) / SQ3) == sp.zeros(3, 1)            # (4.3.78)
    assert sp.simplify(T485_BOOK[:, 0] - (-t1[0] - t1[1] + 2 * t1[2]) / SQ6) == sp.zeros(3, 1)       # (4.3.79a)
    P_A = sum((irrep("T_{1u}", g) for g in C3V), sp.zeros(3)) / 6
    assert sp.simplify(P_A * sp.Matrix([0, 0, 1]) * SQ3 - T485_BOOK[:, 2]) == sp.zeros(3, 1)
    assert T485_BOOK != T43
    lev, Mt = v111_tunnel()
    Spr = sp.Symbol("S'")
    assert lev == {"A_{1g}": Hs + 2 * Ss + 2 * Spr, "A_{2u}": Hs + 2 * Ss - 2 * Spr, "E_g": Hs - Ss - Spr, "E_u": Hs - Ss + Spr}
    lv = {sp.simplify(k_): v for k_, v in Mt.eigenvals().items()}
    assert lv == {Hs + 2 * Ss + 2 * Spr: 1, Hs + 2 * Ss - 2 * Spr: 1, Hs - Ss - Spr: 2, Hs - Ss + Spr: 2}


def _check_xf6():
    assert xf6_coords() == XF6_COORDS_LIST and len(XF6_NAMES) == 21
    rng = np.random.default_rng(4)
    for k, j, b in ((0.18, 1.0, 0.2), tuple(rng.uniform(0.1, 2, 3))):
        F = xf6_force(k, j, b)
        assert np.allclose(F, F.T)
        rA, rB, rC = xf6_rows(k, j, b)
        assert np.allclose(F[0], rA) and np.allclose(F[6, 6:], rB[6:]) and np.allclose(F[18, 18:], rC[18:])
        assert np.allclose(xf6_force_from_rows(k, j, b), F)                  # 三列＋對稱 ＝ 幾何
        Fb = xf6_force_from_rows(k, j, b, book=True)
        for ax in range(3):
            t = xf6_translation(ax)
            assert np.allclose(F @ t, 0) and not np.allclose(Fb @ t, 0)      # 照印的 ⟨B|F|C1⟩、⟨C1|F|C1⟩ 平移不是零模
        for g in ("r1", "R3", "i4", "I.R1", "I"):
            G = xf6_rep(g)
            assert np.allclose(G @ F @ G.T, F) and np.allclose(G @ G.T, np.eye(21))
        for m, M in ((19, 240), (19, 32)):
            w = xf6_numeric(k, j, b, m, M)
            assert np.sum(np.abs(w) < 1e-9) == 6
            f = xf6_formulas(k, j, b, m, M)
            pred = sorted([0.0] * 6 + [v for key, v in f.items() for _ in range(XF6_DEG[key])])
            assert np.allclose(w, pred)
            Fr, Mr, A = xf6_t1u_reduced(k, j, b, m, M)                       # (4.4.29)–(4.4.33)
            assert np.allclose(Fr, [[9 * k + 2 * j + b, 2 * (j - b) * (M + 6 * m)], [2 * (j - b) * (M + 6 * m), (2 * j + 4 * b) * (M + 6 * m) ** 2]])
            assert np.allclose(Mr, np.diag([3 * m, 6 * m * M * (M + 6 * m)]))
            assert np.allclose(A, xf6_t1u_matrix(k, j, b, m, M))
            if abs(b) > 1e-12:
                assert not np.allclose(A, xf6_t1u_matrix(k, j, b, m, M, book=True))
            sec = xf6_secular(k, j, b, m, M)
            assert np.isclose(np.trace(A), sec["S"]) and np.isclose(np.linalg.det(A), sec["P"])
        cs = xf6_coset_omega2(k, j, b)
        f = xf6_formulas(k, j, b)
        for key in ("A_{1g}", "E_g", "T_{2u}", "T_{2g}"):
            assert np.isclose(cs[key], f[key])
        assert np.isclose(cs["Q_AA"], 2 * k + j) and np.isclose(cs["Q_AB"], -math.sqrt(2) * k) and np.isclose(cs["Q_BB"], k + b)
        assert np.isclose(cs["Q_AB_book"], 0)
        assert np.isclose(2 * cs["Q_AA"] - 2 * math.sqrt(2) * cs["Q_AB"] + cs["Q_BB"], 9 * k + 2 * j + b)
    k_, j_, b_ = sp.symbols("k j b", positive=True)
    m_, M_ = sp.symbols("m M", positive=True)
    Am = sp.Matrix(xf6_t1u_matrix(k_, j_, b_, m_, M_))
    sec = xf6_secular(k_, j_, b_, m_, M_)
    assert sp.simplify(Am.trace() - sec["S"]) == 0 and sp.simplify(Am.det() - sec["P"]) == 0       # (4.4.34c,d)
    assert sp.simplify(sp.Matrix(xf6_t1u_matrix(k_, j_, b_, m_, M_, True)).det() - sec["P"]) != 0
    cl = xf6_classify()                                                      # (4.4.3)–(4.4.9)
    assert cl["A"] == [1, 0, 1, 0, 0, 0, 0, 0, 1, 0] and cl["B"] == [0, 0, 0, 1, 1, 0, 0, 0, 1, 1] and cl["C"] == [0] * 8 + [1, 0]
    assert cl["all"] == [1, 0, 1, 1, 1, 0, 0, 0, 3, 1]
    for a in OH_IRREP_NAMES:                                                 # 同型投影與 F 對易
        P = xf6_isotypic(a)
        assert np.allclose(P @ P, P) and np.allclose(P @ xf6_force(0.18, 1, 0.2), xf6_force(0.18, 1, 0.2) @ P)
    # SF6 的數字（(4.4.36)–(4.4.41)、Fig. 4.4.7）
    A = np.array(xf6_t1u_matrix(0.18, 1, 0.2, 19, 32), float)
    assert [_f4(A[0, 0]), f"{A[0, 1]:.3f}", f"{A[1, 0]:.6f}", _f4(A[1, 1])] == ["0.0670", "4.098", "0.000439", "0.1121"]
    sec = xf6_secular(0.18, 1, 0.2, 19, 32)
    assert (_f4(sec["S"]), _f4(sec["P"]), _f4(sec["lp"]), _f4(sec["lm"])) == ("0.1791", "0.0057", "0.1375", "0.0415")
    iso = xf6_isotope()
    assert (_f4(iso[0]), _f4(iso[1])) == ("1.0172", "1.0065") and abs(iso[0] - 1.0173) < 2e-4 and abs(iso[1] - 1.0064) < 2e-4
    w32, w34 = xf6_numeric(0.18, 1, 0.2, 19, 32), xf6_numeric(0.18, 1, 0.2, 19, 34)     # 第二種：21×21
    t32 = [x for x in w32 if any(np.isclose(x, v) for v in (sec["lp"], sec["lm"]))]
    sec34 = xf6_secular(0.18, 1, 0.2, 19, 34)
    assert np.isclose(max(t32), sec["lp"]) and any(np.isclose(x, sec34["lp"]) for x in w34)
    assert (_f4(965 / 948), _f4(615 / 612)) == ("1.0179", "1.0049")
    tv = xf6_t1u_vectors()
    assert [_f4(tv["+"][k]) for k in ("eps", "center", "axial", "equatorial")] == ["0.0682", "-0.1339", "0.1058", "0.0035"]
    assert f"{tv['+']['eps_p']:.5f}" == "0.00117" and f"{tv['-']['eps_p']:.6f}" == "-0.000706"
    assert [_f4(tv["-"][k]) for k in ("eps", "center", "axial", "equatorial")] == ["0.1135", "0.0805", "0.0909", "-0.0794"]
    assert _f4(0.068 * (-0.5) + 0.0012 * 32) == "0.0044"                     # 書上 (4.4.41) 用取整的係數得 0.004
    modes = xf6_modes()
    Fm, Mm = xf6_force(0.18, 1, 0.2), xf6_mass(19, 32)
    for lab, om, v in modes:
        assert np.allclose(Fm @ v, om ** 2 * Mm @ v, atol=1e-9), lab
    assert _f4(math.sqrt(2)) == "1.4142" and (f"{201 / 140:.2f}", f"{524 / 363:.2f}") == ("1.44", "1.44")
    sc = xf6_scan(np.array([0.0, 1.0]))
    assert _f4(sc["A_{1g}"][0]) == "1.0000" and _f4(sc["E_g"][0]) == "1.0000" and _f4(sc["T_{1u}+"][0]) == "1.0763"
    assert sc["T_{2g}"][0] == 0 and sc["T_{1u}-"][0] == 0 and _f4(sc["A_{1g}"][1]) == "2.2361"
    sc6 = xf6_scan(np.array([0.0]), 19, 32, 0.2)
    assert _f4(sc6["T_{1u}+"][0]) == "1.5802" and _f4(sc6["T_{1u}-"][0]) == "0.6045"
    fit = uf6_fit()
    assert _f4(fit["kj"]) == "0.1958" and abs(fit["pred"]["T_{1u}+"] - 624) < 1e-6
    pr = sf6_predict()["pred"]
    assert [f"{pr[k]:.0f}" for k in ("T_{1u}+", "E_g", "T_{2g}", "T_{2u}", "T_{1u}-")] == ["955", "642", "515", "364", "525"]


def _check_quantum():
    assert np.allclose(ladder_commutator(6), np.diag([1.0] * 6 + [-6.0]))
    assert np.allclose(osc_levels(6), np.arange(7) + 0.5)
    a, ad = ladder(12)
    for n in range(1, 5):                                                    # (4.4.59)
        lhs = a @ np.linalg.matrix_power(ad, n) - np.linalg.matrix_power(ad, n) @ a
        rhs = n * np.linalg.matrix_power(ad, n - 1)
        assert np.allclose(lhs[:, :12 - n], rhs[:, :12 - n])
        lhs2 = a @ a @ np.linalg.matrix_power(ad, n) - np.linalg.matrix_power(ad, n) @ a @ a
        rhs2 = n * (n - 1) * (np.linalg.matrix_power(ad, n - 2) if n >= 2 else 0 * ad) + 2 * n * np.linalg.matrix_power(ad, n - 1) @ a
        book = rhs2 + np.linalg.matrix_power(ad, n) @ a @ a
        assert np.allclose(lhs2[:, :10 - n], rhs2[:, :10 - n]) and not np.allclose(lhs2[:, :10 - n], book[:, :10 - n])
    for n in range(5):
        assert [degeneracy(ell, n) for ell in range(1, 6)] == PASCAL_464[n] == [degeneracy_count(ell, n) for ell in range(1, 6)]
    hs = harmonic_states(3, 2)
    assert len(hs) == 6 and (0, 1, 1) in hs and (0, 0, 1) not in hs and len(harmonic_states(2, 2)) == 3
    assert state_count(1000) == 1 + 15 + 6 + 9 + 9                         # 基態、六個基頻（15 態）、2ν6、ν5+ν6、ν4+ν6
    assert _f4(turning_point(0, 1, 1.0)) == "1.0000" and _f4(turning_point(2, 3, 1.0)) == "2.6458"


def _check_problems():
    assert problem_411() == {"C3": 1, "C3v": 2, "D4": 2, "O": 2, "Oh": 2}
    p2 = problem_412()
    assert p2["r1"] == ["1", "r1", "r1^2"] and sorted(p2["R3"]) == sorted(C4Z)
    p3 = problem_413()
    assert len(p3["D2"]) == 4 and [n for _, n in p3["D2"]].count(True) == 1 and len(p3["D3"]) == 4 and len(p3["D4"]) == 3
    assert not any(n for _, n in p3["D3"] + p3["D4"])
    p21 = problem_421()
    assert p21["a"][0] == sp.Matrix([[1, 0, 0], [0, 1, 0], [0, 0, 1], [0, 1, 1], [1, 0, 1]])
    assert p21["a"][1] == sp.Matrix([[1, 0], [0, 1], [1, 1], [1, 2], [2, 1]])
    assert p21["b"][0] == sp.Matrix([[1, 0, 0, 0, 0], [0, 1, 0, 0, 0], [1, 1, 0, 0, 0], [0, 0, 1, 0, 1], [0, 0, 0, 1, 1]])
    assert p21["e"][1][:5, :] == correlation("C4") and p21["e"][1][5:, :] == correlation("C4")
    assert all(sum(r) >= 1 for key in p21 for M in p21[key] for r in M.tolist())
    p23 = problem_423()
    assert {k: (v[0], v[2]) for k, v in p23.items()} == {"a": (48, [3]), "b": (6, [1, 2]), "c": (12, [1, 2]), "d": (4, [1, 1, 1]),
                                                        "e": (8, [1, 2]), "f": (24, [3]), "g": (8, [1, 1, 1]), "h": (3, [1, 1, 1])}
    p31 = problem_431()
    assert p31["C3"] == ([8, 2, 0, 0, 0], [1, 1, 0, 1, 1]) and p31["D4"] == ([3, 0, 3, 1, 1], [1, 0, 1, 0, 0])
    p32 = problem_432()
    assert p32["C3v"][1] == [0, 0, 1, 1, 1, 0, 0, 1, 1, 1] and p32["D4h"][1] == [0, 0, 0, 1, 1, 0, 0, 0, 0, 0]
    p33 = problem_433()
    assert p33 == {"A_{1g}": (1, Hs - 3 * Ss - 3 * Ts), "T_{2g}": (1, Hs + Ss + Ts), "A_{2u}": (1, Hs + 3 * Ss - 3 * Ts),
                   "T_{1u}": (1, Hs - Ss + Ts)}
    p34 = problem_434()
    assert sum(v[0] * OH_CHARS[OH_IRREP_NAMES.index(k), 0] for k, v in p34.items()) == 12 and all(v[0] == 1 for v in p34.values())
    p35 = problem_435()
    assert p35["D5"] == sp.Matrix([[1, 0, 0, 0], [0, 1, 1, 0], [0, 1, 0, 1], [0, 0, 1, 1], [1, 0, 1, 1]])
    assert p35["D3"] == sp.Matrix([[1, 0, 0], [0, 1, 1], [0, 1, 1], [1, 1, 1], [1, 0, 2]])
    assert p35["C5"] == sp.Matrix([[1, 0, 0, 0, 0], [1, 1, 0, 0, 1], [1, 0, 1, 1, 0], [0, 1, 1, 1, 1], [1, 1, 1, 1, 1]])
    assert p35["C3"] == sp.Matrix([[1, 0, 0], [1, 1, 1], [1, 1, 1], [2, 1, 1], [1, 2, 2]])
    assert sum(Y_CLASS_SIZES) == 60 and sum(x * x for x in Y_CHARS.col(0)) == 60
    w0 = problem_436(0.0)
    wq = problem_436(0.7, kind="quadrupole")
    assert np.allclose(wq, w0 - 0.7 / 3)                                     # 四次軸上的四極場：八個頂點等價，只平移
    wd = problem_436(0.7)
    assert not np.allclose(np.diff(wd), np.diff(w0))
    wq3 = problem_436(0.7, axis=(1, 1, 1), kind="quadrupole")
    assert len(np.unique(np.round(wq3 - wq3.min(), 9))) > len(np.unique(np.round(w0, 9))) - 1
    p42 = problem_442()
    assert p42["X8_all"] == [1, 0, 1, 1, 2, 0, 1, 1, 2, 1] and sum(p42["C8H8_vib"]) == 18
    assert sum(n * OH_CHARS[a, 0] for a, n in enumerate(p42["C8H8_vib"])) == 42
    fit = xy4_fit("CF4")
    assert fit["j"] > 0 and abs(fit["rel"]).max() < 0.5


def _check_cards():
    """場景卡上其餘的數字（plan_ch04.md §6）。"""
    hp = nomogram_product("R1^2", "i4")
    assert hp["name"] == "R3" and hp["sign"] == -1 and abs(hp["deg"] - 90) < 1e-9 and nomogram_product("r3", "R2")["sign"] == 1
    assert [cycle_str(diag_perm(g)) for g in ("r1", "R1^2", "R3", "i4")] == ["[132]", "[13][24]", "[1423]", "[12]"]
    lv = vb_levels(0.3)
    assert [_f4(lv[k]) for k in ("1_4", "3_4", "2_4", "T_1(0_4)", "A_1(0_4)", "E(0_4)")] ==         ["-0.5910", "0.5910", "1.9107", "0.0000", "-3.9407", "2.0301"]
    sc = vd_scan([8.0])
    assert [_f4(sc[k][0]) for k in ("A'1", "A'2", "A'3", "E", "B'")] == ["-8.6264", "-1.7677", "8.3941", "0.0000", "2.0000"]
    assert [_f4(e) for e, _ in v111_levels(2.0)] == ["-4.8284", "-1.2361", "-1.2361", "0.8284", "3.2361", "3.2361"]
    assert [e for e, _ in vq_levels(0.0)] == [-4.0, 0, 0, 0, 2, 2]

    def coeffs(v):
        return {XF6_NAMES[i]: round(float(x), 6) for i, x in enumerate(v) if abs(x) > 1e-9}
    assert coeffs(xf6_projected("T_{1u}", 0, 0, "B")) == {"B": 0.5, "R2^2B": -0.5, "R1B": 0.5, "R1^3B": 0.5}        # (4.4.13a)
    assert coeffs(xf6_projected("T_{1u}", 1, 0, "B")) == {"r1B": 0.5, "r2B": 0.5, "R3B": 0.5, "i3B": 0.5}           # (4.4.13b)
    assert coeffs(xf6_projected("T_{1u}", 2, 0, "B")) == {"r1^2B": 0.5, "r4^2B": 0.5, "R2B": -0.5, "R2^3B": 0.5}    # (4.4.13c)
    assert coeffs(xf6_projected("T_{1g}", 0, 1, "B")) == {"r1^2B": 0.5, "r4^2B": -0.5, "R3B": -0.5, "i3B": 0.5}     # (4.4.14)
    r = round(1 / math.sqrt(2), 6)
    assert coeffs(xf6_projected("T_{1u}", 0, 2, "A")) == {"R2A": r, "R2^3A": -r}                                  # (4.4.17)
    a6, a12 = round(1 / math.sqrt(6), 6), round(1 / math.sqrt(12), 6)
    assert set(coeffs(xf6_projected("A_{1g}", 0, 0, "A")).values()) == {a6}
    assert coeffs(xf6_projected("E_g", 0, 0, "A")) == {"A": 2 * a12, "R1^2A": 2 * a12, "R1A": -a12, "R2A": -a12, "R1^3A": -a12, "R2^3A": -a12}
    c, cp, t = xf6_c_vectors(19, 32)
    mass = np.diag(xf6_mass(19, 32))
    assert np.isclose((mass * c) @ xf6_translation(0), 0) and np.isclose((mass * cp) @ xf6_translation(0), 0)   # (4.4.18)–(4.4.19) 動量 0
    assert state_count(2000) == 965 and state_count(3000) == 11674
    assert [O.mul(f"R{k}^2", "i1") for k in (1, 2, 3)] == ["R2", "i2", "R2^3"]                       # (4.1.2)
    assert len(_closure(O, ["R1", "R2"])) == 24                                                        # 習題 4.1.1
    assert [len(O.classes[O.class_of(g)]) for g in ("r1", "R3")] == [24 // len(O.normalizer(g)) for g in ("r1", "R3")] == [8, 6]
    fit = xy4_fit("CF4")                                                                               # 習題 4.4.1
    assert (f"{fit['j']:.1f}", f"{fit['k']:.1f}") == ("485.6", "222.2") and [f"{x:.0f}" for x in fit["pred"]] == ["1108", "446", "1207", "537"]
    rels = np.concatenate([np.abs(xy4_fit(mol)["rel"]) for mol in XY4_NU])
    assert f"{100 * rels.max():.0f}" == "29"
    ex = [math.sqrt(xf6_secular(0.18, 1, 0.2, 18.998, 31.972)[k] / xf6_secular(0.18, 1, 0.2, 18.998, 33.968)[k]) for k in ("lp", "lm")]
    k17 = xf6_isotope(k=0.17)
    assert _f4(ex[0]) == "1.0172" and (_f4(k17[0]), _f4(k17[1])) == ("1.0174", "1.0063")                # 〔注 8〕


def check():
    """每個上卡片／講義的數字用兩種以上獨立算法核對（約 10 秒）。"""
    _check_group()
    _check_irreps()
    _check_corr_frob()
    _check_six()
    _check_pert()
    _check_xf6()
    _check_quantum()
    _check_problems()
    _check_cards()


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
    print(f"sym_oct（第 4 章）：O／Oh／Td／T／Th 由轉動矩陣產生（乘法＝Hamilton 弧、S4 置換）、類代數與最小方程、特徵標表（類代數＝跡）、"
          f"四條子群鏈的分裂冪等元與不可約表示（四方／三方駐波、C4／C3 行進波）、相關表、共軛與正規子群、誘導表示與 Frobenius 互反、"
          f"六點穿隧模型與外場（V_Q、V_D、磁場、(111)）、XF6 的 21 座標簡正模＝對稱公式（含筆誤更正）與同位素位移、升降算子與簡併度、"
          f"習題 4.1.1–4.4.2 全部核對通過（{time.time() - t0:.1f} 秒）。")
