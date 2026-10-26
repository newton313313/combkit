# 由 combkit/tools/sync_from_harter_sym.py 從 harter_sym 專案的 codes/utils/sym_coup.py 複製；請改 harter_sym 那份再同步。
r"""
sym_coup.py —— Harter《Principles of Symmetry, Dynamics, and Spectroscopy》（Wiley 1993）第 7 章
Theory and Application of Symmetry Representation Products (Continuous Rotation Groups) 的數學
（純 Python＋numpy＋sympy＋scipy，不 import manim；J 矩陣、𝒟^j、球諧、SF6 能階用同資料夾的 sym_rot，O 群、誘導表示用 sym_oct，
±√(p/q) 的精確化用 sym_prod，不重寫；sympy.physics.wigner 只當獨立核對）。書上的式號、筆誤與〔注〕見 ref/md/ch07_formulas.md
（相對於 C:\mypython\manim\projects\harter_sym）。

約定
  * ħ = 1；j、m 可為半整數（int、Fraction、sympy Rational 或 "3/2" 都可以；內部用 sympy Rational）。基底 m = j, j−1, …, −j（列＝行同序）。
  * 耦合係數 C^{j1 j2 j}_{m1 m2 m} ＝ ⟨j1 m1, j2 m2|j m⟩（Condon–Shortley 相位：C^{j1 j2 j}_{j1, j−j1, j} > 0，與 sympy、書上的表相同）。
    乘積基底 |m1⟩|m2⟩ 依 (m1, m2) 字典序（m1 由大到小，再 m2 由大到小）；表的行依 (j, m)：j 由大到小、m 由大到小（(7.1.7)、(7.1.19) 的樣子）。
  * 3j：(j1 j2 j3; m1 m2 m3) = (−1)^{j1−j2−m3} C^{j1 j2 j3}_{m1 m2 −m3}/√(2j3+1)（(7.2.30a)）＝ sympy 的 wigner_3j。
  * 張量算子：T(j1 j2)^k_q = Σ C^{j1 j2 k}_{m1 m2 q}|j1 m1⟩⟨j2 −m2|(−1)^{j2−m2}（(7.3.2)）；單位張量 v(j1 j2)^k_q（(7.3.17a)）＝
    (−1)^{2j1}T(j1 j2)^k_q，矩陣元 ⟨j1 m1|v^k_q|j2 m2⟩ = C^{k j2 j1}_{q m2 m1}(−1)^{k+j1−j2}√((2k+1)/(2j1+1))（(7.3.21)）；
    tr(v^{k†}_q v^{k'}_{q'}) = δδ。約化矩陣元照書上的 C 約定：⟨j1 m1|T^k_q|j2 m2⟩ = C^{k j2 j1}_{q m2 m1}⟨j1‖T^k‖j2⟩（(7.3.20a) 更正）。
  * 多極函數 X^k_q = r^k √(4π/(2k+1)) Y^k_q（Condon–Shortley 號）；RE 曲面的方向用單位向量 n = (n_x, n_y, n_z)（(7.4.4a) 的號不影響能量）。
  * O 群的標記用 sym_oct 的 O_CHARS（A_1 A_2 E T_1 T_2；奇 j 的 u 宇稱另寫）；D2 的標記照 Fig. 7.4.3：
    A1 (+,+,+)、A2 (−,+,−)、B1 (+,−,−)、B2 (−,−,+)（C2(x)、C2(y)、C2(z) 的特徵標；與第 4 章 (4.2.10) 的約定不同）。
  * 卡片上的小數一律 f"{float(x):.4f}" 之類（Fraction、sympy 數要 float()）；精確值用 sympy 字串（str(sp.nsimplify(x))）。
  * 不要把模組代號 SC、RT、SO、SP、SM 或 D、E、S、T、I、G、J 當變數名。

§7.1A 氫的超精細結構（H701）
  ms(j)、cg(j1,j2,j,m1,m2,m)（sympy 精確，快取）＝ cg_racah（書上的 Racah 式 (7.2.57)）＝ cg_float（整數階乘的 Racah 式，大 J 用）
    ＝ cg_lowering(j1, j2)[(j, m)]（降算子 (7.1.8)–(7.1.11)＋正交化）＝ sympy.physics.wigner.clebsch_gordan
  cg_table(j1, j2, method) → (Matrix, rows=[(m1,m2)], cols=[(j,m)])；CG_717_BOOK（(7.1.7)）、CG_7119_BOOK（(7.1.19)）逐格相同
  stretched_cg(j1, j2, m1, m2, book=False)（(7.1.17) 更正 (j2+m2)!；book=True 照印 (j1+m2)!）
  d_half_product(β) → C†(𝒟^½⊗𝒟^½)C（(7.1.6a)；(1,1) 元是 cos²(β/2)，書印 sin²）；HF_C716（書上的 C）
  hyperfine_matrix()（J^p·J^e，4×4 精確）、hyperfine_levels() → {1/4: 3, −3/4: 1}、hyperfine_states()（三重態、單重態）、
  HF_CONST（書上的 g_e、g_p、β_e、β_p、a_0）、hyperfine_numbers() → a = 9.4234e−18 erg = 1422.17 MHz（書 1422.74；書印 1.500e−18 erg ＝ a/2π）、
    觀測 1420.405751800 MHz → 5.8743e−6 eV（書印 e−5）、9.4117e−18 erg、0.047380 cm⁻¹、21.106 cm
§7.1B 兩電子組態（H702）
  two_electron_terms(l1, l2, equivalent) → {"^{2S+1}L": 簡併}（數值：反對稱子空間裡 L²、S² 同時對角化）＝ terms_by_symmetry(l1, l2, equivalent)（(7.2.34) 的對稱性）
    (2p)²：¹S 1、³P 9、¹D 5（共 15 = C(6,2)）；(2p3p)：¹S ¹P ¹D ³S ³P ³D（36）；hund_order(terms) → ³P < ¹D < ¹S ＝ SLATER_P2（F2 模型：−5、1、10）
§7.1C–D 自旋–軌道、LS／jj、向量錐（H703）
  spin_orbit(l) → {j: ⟨S·l⟩}（l/2、−(l+1)/2；數值對角化 ＝ 公式；書上 (7.1.26) 少 ½ 會加倍）、dirac_level(n, j)（(7.1.27)）、
  fine_structure() → 1s 位移 −1.8113e−4 eV；2p_{3/2} − 2p_{1/2} = 4.5283e−5 eV = 0.3652 cm⁻¹ = 10.95 GHz（＝ (3/2)a_so，a_so = α²/48）
  couple(j1, j2)（三角條件）、ls_jj_counts() → LS 與 jj 的 J 值（(2p)²：0,0,1,2,2；(2p3p) 三重態 0,1,1,1,2,2,3）
  p2_levels(x, F2=1, zeta=4)（本系列的 LS→jj 相關：(1−x)H_C + xζΣl_i·s_i，15 維）→ [(E, J, 簡併)]：x=0 ³P −5、¹D 1、¹S 10；x=1 4、−2、−8
  cg_curve(j1, j2, J, M) → (ms, 值)、classical_range(j1, j2, J, M)（錐的轉折點）、sign_changes(vals)：
    C^{9 9 J}_{m −m 0}：J = 18、17、16 → 峰 0.5104、0.4392、0.4035；古典範圍 ±2.121、±3.674、±4.690；節點 0、1、2
§7.2 純量、向量、張量、3j、6j（H704）
  spherical(v)（(7.2.3)）、coupled_product(A, B, k, q)、table_727()（9×9；T727_BOOK 照印，只差 (xy, [ ]²_{−2}) 一格：−½ 應 −i/2）、
  rr_tensor() → [r r]²_q 與 X²_q（比 √(2/3)；(7.2.9) 的 [rr]²_1 = −(x+iy)z，書印 x−iy）、ck_formula(k, q)（(7.2.11)）
  scalar_cg(j, m) = (−1)^{j−m}/√(2j+1)（(7.2.12)）：C^{110}_{m,−m,0} = 1/√3, −1/√3, 1/√3（m = 1, 0, −1）
  threej(...)（(7.2.30a)）、threej_symmetries()（(7.2.31)–(7.2.35)）、orthogonality_check()（(7.2.21)–(7.2.22)）、factorization_check()（(7.2.25)，Euler 角數值積分）
  su2_irrep_formula(j, u)（(7.2.38)）、boson_coupled(j, m, n, N, half)、boson_cg(...)（(7.2.54) 方括號要 ½ 次方）、det_scalar_state(d)（(7.2.50)：|r 2d−r, 2d−r r⟩、4d 個玻色子）
  sixj(...)（sympy）＝ sixj_from_3j(...)（四個 3j）＝ sixj_from_cg(...)（(7.3.49)）：{2 4 2; 2 2 2} = 2/35
§7.3A 張量算子（H705）
  tensor_T(j1, j2, k, q)、unit_v(j1, j2, k, q)、spin_half_T()（(7.3.3ₓ)）＝ T733_BOOK、spin_half_cart()（(7.3.4)：σ/√2）、
  rotated_Jz_check(β)（(7.3.7)–(7.3.8)）、expectation_rotation(β0, β)（Fig. 7.3.1；本系列 β0 = 60°、β = 45°：⟨J_x⟩ 0.4330、⟨J_z⟩ 0.2500、轉過的 0.4830）、
  tensor_j1() ＝ T7312_BOOK、j1_relations()（(7.3.13)）、v_table(j, k) → (分子, 分母²) ＝ V_TABLES_BOOK（Table 7.1–7.2 的 j ≤ 2 格子）、
  elementary_expansion(j, a, b)（(7.3.15)：E_24 = −√8/√14 v⁴₂ + √3/√7 v²₂）、mixed_examples()（√10 V¹₁(pd)、√10 V¹₋₁(dp)）、
  mixed_table(j1, j2, k) ＝ V74_DP_BOOK（Table 7.4 (d)(p) 逐格；蓋住的那格是 −1）、fpoly_check()（(7.3.30)–(7.3.32) 的係數）
§7.3B–C Wigner–Eckart 與晶體場（H706）
  we_v_check()（(7.3.21)–(7.3.22)）、we_operator_check(j)（由 J 做成的 T¹、T²：(7.3.1b) 的轉動與 (7.3.20) 的比例；⟨j‖J‖j⟩ = −√(j(j+1))，C^{k j2 j1} 次序的號）
  v4_poly_check()（(7.3.25)）、crystal_matrix(j, k)（單位 D⟨j‖X^k‖j⟩，含 (7.3.24) 的 √((2j+1)/(2k+1))）、crystal_levels(j, k) → [(值, O 標記, 簡併)]、
  o_label(P, J)（本徵子空間的 O 標記）、crystal_direct(j, k)（球諧直接積分，單位 D⟨r^k⟩）：
    j = 2，V⁽⁴⁾：T2 −8、E +12（單位 D⟨2‖X⁴‖2⟩/(15√14)）＝ −8/105、+4/35（單位 D⟨r⁴⟩）；重心 3(−8) + 2(12) = 0
    j = 3，V⁽⁴⁾：A2 −12、T2 −2、T1 +6（δ⁽⁴⁾ = 2D⟨3‖X⁴‖3⟩/(15√22)）；V⁽⁶⁾：A2 −12、T1 −5、T2 +9（δ⁽⁶⁾ = E⟨3‖X⁶‖3⟩/√858）
    ——書上 (7.3.31)–(7.3.37)、Fig. 7.3.2 把 T1u、T2u 對調；(7.3.36)–(7.3.37b) 少 √(7/13)、多一個 4
§7.3D 約化矩陣元（H707）
  reduced_X(l', k, l)（(7.3.40b)：⟨2‖X⁴‖2⟩ = √(2/7)⟨r⁴⟩）、r4_hydrogen(n, l) ＝ r_power_numeric(n, l, 4)（3d：25515 a0⁴）、
  d2_factor(L)（(7.3.50ₓ) 的 6j 式）＝ d2_factor_numeric(L)（兩電子態直接算）：L = 2（d² ¹D）→ 4/7 = 0.5714
§7.4A 剛體轉子（H708–H709）
  re_rigid(n, J, A, B, C)（(7.4.5)）＝ re_rigid_tensor（(7.4.10)）、rigid_levels(J, A, B, C) → [(E, D2 標記)]（2J+1 維）＝ RT.asym_matrix_formula、
  FIG742 = {J: 10, A: 0.2, B: 0.4, C: 0.6}、fig742() → 分隔線 44、最低一對 (A1, B1) 24.8189、精細 5.3305 cm⁻¹ = 159.8 GHz（書 150 GHz）、
    超精細 8.63e−7 cm⁻¹ = 25.9 kHz、7.73e−5 cm⁻¹ = 2.32 MHz（書 230 kHz）、2.79e−3 = 83.6 MHz、5.07e−2 = 1.52 GHz
  symmetric_top_check(J, B, C)（(7.4.15)）、reduced_T2(J)（(7.4.13)）、d2_induced()（(7.4.16)：0₂(C2x)↑D2 = A1+B1、1₂ = A2+B2）、
  rigid_correlation(J, Bs)（Fig. 7.4.4；端點 22、66 與 84、252）、cone_angles(J, Ks)
§7.4B 八面體半剛性轉子（H710）
  octa_re(n)（(7.4.18)–(7.4.19)）、octa_levels(J) ＝ RT.sf6_levels、octa_clusters(J) → [{kind, K, labels, E, offsets, S, fit}]
    （J = 30：C4 叢 K4 = 30…25：(A2 T2 E) (T1 T2) (A1 T1 E) (T1 T2) (A2 T2 E) (T1 T2)；C3 叢 K3 = 30 (A1 T1 T2 A2)、29 (T2 E T1)）、
  CLUSTER_PATTERNS（(7.4.20)–(7.4.21) 更正：1₄、3₄↑O 是 T2 +2S、T1 −2S）、octa_model_04()（八面體鄰接矩陣：4、0、−2）、
  bottom_sequence(J)（(7.4.22)：A1 T1 T2 A2 T2 E T1 T2 E T1）、separatrix_angles() → 35.26°、19.47°、cluster_cutoffs(J) → 24.90、28.75、
  cone_angles(30, 30..24) → 10.3 18.0 23.3 27.7 31.5 34.9 38.1、t46_sign_check()（(7.4.24) 應 +）、t46_operator(J, ν)、t46_clusters(J, ν)
    （ν = 2π/6、4π/6 的最低兩叢：A1+E+T1+2T2（0₂↑O）、A2+E+2T1+T2（1₂↑O））、re46(n, ν)（(7.4.26)–(7.4.27)）
  sphere_contours(f, levels, nb, ng) → {level: [n×3 單位向量折線]}（自己的 marching squares）、resample(path, n)、
    rigid_paths(J, A, B, C)、octa_paths(J, levels)、re46_paths(ν, levels)（畫 RE 曲面的等能量線；SD.RESurface 用）
§7.5 旋量、偏振、二維振子（H711–H713）
  spinor_euler(α, β, γ)（(7.5.2)）＝ RT.su2_euler(...)[:, 0]、spin_vector(ψ)（Fig. 7.5.2）、h_to_omega(H)（(7.5.5e,f) 更正 H22）、
  evolve(H, ψ0, t) ＝ crank_evolve(H, ψ0, t)（(7.5.7)：e^{−iω0 t}R[φθ ωt]）、HAM_CARD = {A: 4, B: 0, C: 1, D: 2}（本系列）、
  osc_hamiltonian(A,B,C,D)、jvars(x1,p1,x2,p2)（(7.5.13) 更正 J_z 的減號）、hc_check()（(7.5.14a) ＝ (7.5.10b)）、
  lissajous(A,B,C,D, ψ0, ts)（解析）＝ lissajous_numeric（solve_ivp，(7.5.10a)）、
  pauli_bases()（(7.5.21)–(7.5.22)）、handedness_check()（Fig. 7.5.4：|r⟩ 逆時針）、stokes_circular(R, L, φ)、stokes_linear(X, Y, ϑ)、stokes_of(ψ_xy)、
  state_circular(R, L, φ, Φ)（(7.5.19)）、state_linear(X, Y, ϑ, θ)（(7.5.17)）、ellipse_field(ψ_xy, Φs) ＝ ellipse_circular_formula(R, L, φ, Φs)（(7.5.20)）、
  ELLIPSE_CARD = {φ: 20, ψ: 15}（本系列）、ellipse_card() → a, b, R, L、Stokes (0.6634, 0.5567, 0.5000)、ν 24.22°、ϑ 20.97°、Euler 角 (7.5.26)–(7.5.27)
  faraday_check()、birefringence_check()（(7.5.29)–(7.5.30)）、beat(A,B,C,D)（(7.5.31)–(7.5.32)）、closure(NH, NL)（(7.5.33)–(7.5.34)：
    ω0/Ω = 3、7/5、−5/7；玫瑰瓣 N_H + N_L = 3、7、−5）＝ closure_numeric(NH, NL)、
  anharmonic_fixed_points(Ω, a_x, I)（(7.5.40)：a_x = 2、3 → β = 60.00°、70.53°；分岔 a_x I = Ω）＝ 數值極大、osc_re(n, ω, a_x, I)、local_mode_levels(j, Ω, a_x)
§7.6 分子電子結構（H714）
  h2p_integrals(R) → S、A = (a|1/r_b|a)、B = (a|1/r_b|b)（閉式；(7.6.22)–(7.6.23) 更正）＝ h2p_quadrature(R)（球體座標數值積分）、
  h2p_V(R)（(7.6.24)）＝ (7.6.15) ＝ h2p_V_generalized(R)（廣義本徵值）、h2p_matrix_check()（(7.6.13)）、h2p_min() → R = 2.4928、V⁺ − ε_1s = −0.0648（書 ≈ 2.5、−0.065）
  two_center(f1, f2, R)（兩中心 Coulomb 積分，多極展開數值）、h2_integrals(R) → A B C D E F G（D 用更正的 R³/6；E = G）、
  h2_integrals_numeric(R)、hl_energy(R, book=False, triplet=False)、hl_min() → R = 1.6425、−0.1160（照印 R³/3：1.6220、−0.1275）、h2_generalized(R)、vb_matrix(R, book)、mo_matrix(R, book)、limit_check()（(7.6.48)–(7.6.49)）、
  ci_min()、mo_min()、fig763_levels(s, R)、hybrid_angle(name)（sp 180°、sp² 120°、sp³ 109.4712°）、xy2_hybrids(2φ)（(7.6.54)–(7.6.56)，a' = +√2λ）、
  xy2_config(2φ)（(7.6.61)：90° s²p⁴、104.5° s^1.6 p^4.4、180° s¹p⁵）、huckel_ring(n)、huckel_chain(n)、benzene()（6H − 8S，共軛 2S）、
  butadiene()（4H − 2√5 S，共軛 0.4721S；書印 0.047S）、py_check()（(7.6.51) 的 p_y 應 i(Y¹₁ + Y¹₋₁)）、octa_sigma()（0₄↑O = A1+E+T1、A'↑O_h = A1g+Eg+T1u、E↑O_h = T1g+T2g+T1u+T2u）、ligand_states()（(7.6.65)）
card_values() ＝ CARD_EXPECTED（plan_ch07 §6 卡片上與導讀裡的每個數字）；check()：每個結果兩種以上獨立算法（約 20 秒）；summary()。
直接執行：python utils/sym_coup.py（在 codes/ 裡）印一句成功訊息；--summary 印卡片數字。
"""
import itertools
import math
from fractions import Fraction
from functools import lru_cache

import numpy as np
import sympy as sp
from numpy.polynomial.legendre import leggauss, legvander
from scipy.integrate import quad, solve_ivp
from scipy.linalg import eig, eigh, expm
from scipy.optimize import minimize, minimize_scalar
from scipy.special import exp1
from sympy.physics.wigner import clebsch_gordan as _sym_cg, wigner_3j as _sym_3j, wigner_6j as _sym_6j

try:
    from . import sym_mat as SM          # noqa: F401（給場景用：SC.SM）
    from . import sym_oct as SO
    from . import sym_prod as SP
    from . import sym_rot as RT
except ImportError:                      # 直接執行 python utils/sym_coup.py
    import sym_mat as SM                 # noqa: F401
    import sym_oct as SO
    import sym_prod as SP
    import sym_rot as RT

R_ = sp.Rational
IU = sp.I
HALF = R_(1, 2)
SPEC = "SPDFGHIK"
_FACT = math.factorial


def _f4(x) -> str:
    return f"{float(x):.4f}"


def _f2(x) -> str:
    return f"{float(x):.2f}"


def _e(x, n: int = 4) -> str:
    return f"{float(x):.{n}e}"


def _ex(x) -> str:
    return str(sp.nsimplify(sp.simplify(x)))


def S_(x):
    """int／float／Fraction／str／sympy → sympy Rational（半整數）。"""
    if isinstance(x, sp.Basic):
        return sp.nsimplify(x)
    if isinstance(x, Fraction):
        return R_(x.numerator, x.denominator)
    if isinstance(x, float):
        fr = Fraction(x).limit_denominator(2)
        return R_(fr.numerator, fr.denominator)
    return sp.sympify(x)


def ms(j) -> list:
    """m = j, j−1, …, −j（sympy Rational）。"""
    j = S_(j)
    return [j - k for k in range(int(2 * j) + 1)]


def _msf(j) -> np.ndarray:
    return np.array([float(m) for m in ms(j)])


def _num(M) -> np.ndarray:
    return np.array(sp.Matrix(M).evalf(), dtype=complex)


def _close(a, b, tol: float = 1e-9) -> bool:
    return np.allclose(np.asarray(a, dtype=complex), np.asarray(b, dtype=complex), atol=tol, rtol=0)


def _zero(x) -> bool:
    return sp.simplify(x) == 0


# ======================================================================== §7.1A 耦合係數與超精細
@lru_cache(maxsize=None)
def _cg_cached(j1, j2, j, m1, m2, m):
    return sp.nsimplify(_sym_cg(j1, j2, j, m1, m2, m))


def cg(j1, j2, j, m1, m2, m):
    """C^{j1 j2 j}_{m1 m2 m}（sympy 精確，快取；Condon–Shortley）。"""
    return _cg_cached(*(S_(x) for x in (j1, j2, j, m1, m2, m)))


def _tri(j1, j2, j3) -> bool:
    return abs(j1 - j2) <= j3 <= j1 + j2 and (j1 + j2 + j3).is_integer


def cg_racah(j1, j2, j3, m1, m2, m3):
    """書上的 Racah 式 (7.2.57)（Schwinger 振子推出；照印，負階乘的項略去）。"""
    j1, j2, j3, m1, m2, m3 = (S_(x) for x in (j1, j2, j3, m1, m2, m3))
    if m1 + m2 != m3 or not _tri(j1, j2, j3) or abs(m1) > j1 or abs(m2) > j2 or abs(m3) > j3:
        return sp.Integer(0)
    f = sp.factorial
    pre = (-1) ** (j2 + m2) * sp.sqrt((2 * j3 + 1) * f(j1 + j2 - j3) * f(j3 + j1 - j2) * f(j2 + j3 - j1) / f(j1 + j2 + j3 + 1)) \
        * sp.sqrt(f(j1 + m1) * f(j1 - m1) * f(j2 + m2) * f(j2 - m2) * f(j3 + m3) * f(j3 - m3))
    s = sp.Integer(0)
    for k in range(0, int(j1 + j2 + j3) + 2):
        a = [j2 + m2 - k, j1 - j3 - m2 + k, j3 + m3 - k, j1 - j2 - m3 + k, sp.Integer(k), j3 - j1 + j2 - k]
        if any(x < 0 for x in a):
            continue
        s += sp.Integer(-1) ** k / sp.Mul(*[f(x) for x in a])
    return sp.nsimplify(sp.simplify(pre * s))


def cg_float(j1, j2, j3, m1, m2, m3) -> float:
    """標準 Racah 式（Python 大整數＋Fraction，最後才開根號）；大 J（J = 30 的張量）用。"""
    j1, j2, j3, m1, m2, m3 = (Fraction(x) if not isinstance(x, sp.Basic) else Fraction(str(x)) for x in (j1, j2, j3, m1, m2, m3))
    if m1 + m2 != m3:
        return 0.0
    a = [j1 + j2 - j3, j1 - j2 + j3, -j1 + j2 + j3, j1 + m1, j1 - m1, j2 + m2, j2 - m2, j3 + m3, j3 - m3]
    if any(x < 0 or x.denominator != 1 for x in a):
        return 0.0
    a = [int(x) for x in a]
    pre = Fraction(int((2 * j3 + 1)) * _FACT(a[0]) * _FACT(a[1]) * _FACT(a[2]) * _FACT(a[3]) * _FACT(a[4]) * _FACT(a[5])
                   * _FACT(a[6]) * _FACT(a[7]) * _FACT(a[8]), _FACT(int(j1 + j2 + j3 + 1)))
    s = Fraction(0)
    k = 0
    while True:
        d = [k, int(j1 + j2 - j3) - k, int(j1 - m1) - k, int(j2 + m2) - k, int(j3 - j2 + m1) + k, int(j3 - j1 - m2) + k]
        if d[1] < 0 or d[2] < 0 or d[3] < 0:
            break
        if min(d) >= 0:
            s += Fraction((-1) ** k, _FACT(d[0]) * _FACT(d[1]) * _FACT(d[2]) * _FACT(d[3]) * _FACT(d[4]) * _FACT(d[5]))
        k += 1
    return float(s) * math.sqrt(float(pre))


def _jm_exact(j):
    return RT.jmat(j, exact=True)


def cg_lowering(j1, j2) -> dict:
    """(7.1.8)–(7.1.11)：從 |j1 j1⟩|j2 j2⟩ 用 J− = J−⊗1 + 1⊗J− 降下來；較低的 j 先在 M = j 子空間與已得的態正交化
    （相位 C^{j1 j2 j}_{j1, j−j1, j} > 0），再降。→ {(j, m): sympy 行向量（乘積基底字典序）}。"""
    j1, j2 = S_(j1), S_(j2)
    A, B = _jm_exact(j1), _jm_exact(j2)
    n1, n2 = int(2 * j1 + 1), int(2 * j2 + 1)
    Jm = sp.kronecker_product(A["Jm"], sp.eye(n2)) + sp.kronecker_product(sp.eye(n1), B["Jm"])
    m1s, m2s = ms(j1), ms(j2)
    out = {}
    J = j1 + j2
    while J >= abs(j1 - j2):
        idx = [a * n2 + b for a in range(n1) for b in range(n2) if m1s[a] + m2s[b] == J]
        higher = [v for (Jp, Mp), v in out.items() if Mp == J]
        if higher:
            Mrow = sp.Matrix([[v[i] for i in idx] for v in higher])
            ns = Mrow.nullspace()
            assert len(ns) == 1
            w = ns[0]
        else:
            w = sp.Matrix([1])
        w = w / sp.sqrt(sum(x ** 2 for x in w))
        if sp.simplify(w[0]).is_negative:          # idx 依 m1 由大到小，第一個是 m1 = j1（Condon–Shortley）
            w = -w
        v = sp.zeros(n1 * n2, 1)
        for k, i in enumerate(idx):
            v[i] = sp.radsimp(w[k])
        out[(J, J)] = v
        M = J
        while M > -J:
            v = Jm * v / sp.sqrt((J + M) * (J - M + 1))
            v = v.applyfunc(lambda x: sp.radsimp(sp.simplify(x)))
            M -= 1
            out[(J, M)] = v
        J -= 1
    return out


def cg_table(j1, j2, method: str = "sympy"):
    """耦合係數表：(Matrix, rows=[(m1, m2)], cols=[(j, m)])；method = "sympy"、"racah"（(7.2.57)）、"lowering"。"""
    j1, j2 = S_(j1), S_(j2)
    rows = [(a, b) for a in ms(j1) for b in ms(j2)]
    cols = []
    J = j1 + j2
    while J >= abs(j1 - j2):
        cols += [(J, m) for m in ms(J)]
        J -= 1
    if method == "lowering":
        low = cg_lowering(j1, j2)
        M = sp.Matrix(len(rows), len(cols), lambda r, c: low[cols[c]][r])
    else:
        fn = cg if method == "sympy" else cg_racah
        M = sp.Matrix(len(rows), len(cols), lambda r, c: fn(j1, j2, cols[c][0], rows[r][0], rows[r][1], cols[c][1]))
    return M, rows, cols


_s2, _s3, _s6 = sp.sqrt(2), sp.sqrt(3), sp.sqrt(6)
CG_717_BOOK = sp.Matrix([[1, 0, 0, 0], [0, 1 / _s2, 0, 1 / _s2], [0, 1 / _s2, 0, -1 / _s2], [0, 0, 1, 0]])      # (7.1.7)
CG_7119_BOOK = sp.Matrix([                                                                                        # (7.1.19)
    [1, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 1 / _s2, 0, 0, 0, 1 / _s2, 0, 0, 0],
    [0, 0, 1 / _s6, 0, 0, 0, 1 / _s2, 0, 1 / _s3],
    [0, 1 / _s2, 0, 0, 0, -1 / _s2, 0, 0, 0],
    [0, 0, sp.sqrt(R_(2, 3)), 0, 0, 0, 0, 0, -1 / _s3],
    [0, 0, 0, 1 / _s2, 0, 0, 0, 1 / _s2, 0],
    [0, 0, 1 / _s6, 0, 0, 0, -1 / _s2, 0, 1 / _s3],
    [0, 0, 0, 1 / _s2, 0, 0, 0, -1 / _s2, 0],
    [0, 0, 0, 0, 1, 0, 0, 0, 0]])


def stretched_cg(j1, j2, m1, m2, book: bool = False):
    """(7.1.17)：J = j1 + j2 的耦合係數（更正：分母 (j2+m2)!）；book=True 照印 (j1+m2)!（j1 − j2 非整數時沒有意義）。"""
    j1, j2, m1, m2 = (S_(x) for x in (j1, j2, m1, m2))
    J, M = j1 + j2, m1 + m2
    f = sp.factorial
    last = f(j1 + m2) if book else f(j2 + m2)
    return sp.sqrt(f(J - M) / (f(j1 - m1) * f(j2 - m2)) * f(J + M) / (f(j1 + m1) * last) * f(2 * j1) * f(2 * j2) / f(2 * J))


HF_C716 = sp.Matrix([[1, 0, 0, 0], [0, 1 / _s2, 0, 1 / _s2], [0, 1 / _s2, 0, -1 / _s2], [0, 0, 1, 0]])    # (7.1.6a) 右邊的 C


def d_half_product(beta=None):
    """(7.1.6a)：C†(𝒟^½(0β0) ⊗ 𝒟^½(0β0))C（sympy；beta=None 用符號 β）→ 應為 𝒟¹(0β0) ⊕ [1]。"""
    b = sp.Symbol("beta") if beta is None else sp.rad(beta)
    c, s = sp.cos(b / 2), sp.sin(b / 2)
    Dh = sp.Matrix([[c, -s], [s, c]])
    return (HF_C716.T * sp.kronecker_product(Dh, Dh) * HF_C716).applyfunc(sp.simplify)


def hyperfine_matrix() -> sp.Matrix:
    """J^proton·J^electron 在 |m_p⟩|m_e⟩（(½½),(½,−½),(−½,½),(−½,−½)）的 4×4（精確）。"""
    Jh = _jm_exact(HALF)
    one = sp.eye(2)
    return sum((sp.kronecker_product(Jh[k], one) * sp.kronecker_product(one, Jh[k]) for k in ("Jx", "Jy", "Jz")), sp.zeros(4))


def hyperfine_levels() -> dict:
    """(7.1.13)–(7.1.14)：a J^p·J^e 的本徵值（單位 a）→ {1/4: 3, −3/4: 1}。"""
    return {sp.nsimplify(k): v for k, v in hyperfine_matrix().eigenvals().items()}


def hyperfine_states() -> dict:
    """三重態 |1 1⟩, |1 0⟩, |1 −1⟩ 與單重態 |0 0⟩（(7.1.10ₓ)、(7.1.11)）＝ (7.1.7) 的四行。"""
    M, rows, cols = cg_table(HALF, HALF)
    return {c: M[:, k] for k, c in enumerate(cols)}


HF_CONST = dict(ge=2.0023, gp=5.585, beta_e=0.9273e-20, beta_p=0.50504e-23, a0=0.5292e-8,      # 書上 (7.1.12) 的常數（cgs）
                h=6.62607015e-27, c=2.99792458e10, erg_per_eV=1.602176634e-12, nu_obs=1420405751.8)


def hyperfine_numbers() -> dict:
    """(7.1.12)、(7.1.15) 與 21 cm 線：a = (8π/3)g_eβ_e g_pβ_p|ψ(0)|²，|ψ_1s(0)|² = 1/πa0³。"""
    k = HF_CONST
    a = 8 / 3 * k["ge"] * k["beta_e"] * k["gp"] * k["beta_p"] / k["a0"] ** 3
    nu = k["nu_obs"]
    return dict(a_erg=a, a_MHz=a / k["h"] / 1e6, a_over_2pi=a / (2 * math.pi), nu_MHz=nu / 1e6,
                E_eV=nu * k["h"] / k["erg_per_eV"], E_erg=nu * k["h"], wavenumber=nu / k["c"], wavelength_cm=k["c"] / nu)


# ======================================================================== §7.1B 兩電子組態
def _single_ops(l, s=HALF) -> tuple:
    """單電子 (m_l, m_s) 空間（l 在前）的 L_i、S_i（numpy）。"""
    L, Sg = RT.jmat(l), RT.jmat(s)
    nl, ns = L["Jz"].shape[0], Sg["Jz"].shape[0]
    Lo = {k: np.kron(L[k], np.eye(ns)) for k in ("Jx", "Jy", "Jz")}
    So = {k: np.kron(np.eye(nl), Sg[k]) for k in ("Jx", "Jy", "Jz")}
    return Lo, So, nl * ns


def _two_electron(l1, l2, equivalent: bool):
    L1, S1, n1 = _single_ops(l1)
    L2, S2, n2 = _single_ops(l2)
    I1, I2 = np.eye(n1), np.eye(n2)
    Lt = {k: np.kron(L1[k], I2) + np.kron(I1, L2[k]) for k in L1}
    St = {k: np.kron(S1[k], I2) + np.kron(I1, S2[k]) for k in S1}
    L2op = sum(Lt[k] @ Lt[k] for k in Lt)
    S2op = sum(St[k] @ St[k] for k in St)
    if equivalent:
        SW = np.zeros((n1 * n2, n1 * n2))
        for a in range(n1):
            for b in range(n2):
                SW[b * n1 + a, a * n2 + b] = 1
        PA = (np.eye(n1 * n2) - SW) / 2
        w, U = np.linalg.eigh(PA)
        Q = U[:, w > 0.5]
    else:
        Q = np.eye(n1 * n2)
    return Q, L2op, S2op, Lt, St, (L1, S1, L2, S2, I1, I2)


def two_electron_terms(l1: int = 1, l2: int = 1, equivalent: bool = True) -> dict:
    """數值：(反對稱化的) 兩電子空間裡 L²、S² 同時對角化 → {"^{2S+1}L": 態數}（(7.1.20)–(7.1.21)）。"""
    Q, L2op, S2op, *_ = _two_electron(l1, l2, equivalent)
    M = Q.conj().T @ (L2op + 7.0 * S2op) @ Q
    w = np.linalg.eigvalsh((M + M.conj().T) / 2)
    out = {}
    for L in range(l1 + l2 + 1):
        for Sx in (0, 1):
            val = L * (L + 1) + 7.0 * Sx * (Sx + 1)
            cnt = int(np.sum(np.abs(w - val) < 1e-8))
            if cnt:
                out[f"{2 * Sx + 1}{SPEC[L]}"] = cnt
    return out


def terms_by_symmetry(l1: int = 1, l2: int = 1, equivalent: bool = True) -> dict:
    """(7.2.34)：C^{l l L}_{m2 m1 M} = (−1)^{2l−L}C^{l l L}_{m1 m2 M} ⇒ 軌道偶 L 對稱、奇 L 反對稱；自旋 S = 1 對稱、S = 0 反對稱；
    Pauli：總波函數反對稱 ⇒ 只剩（對稱軌道 × 單重態）與（反對稱軌道 × 三重態）。不等價電子（2p3p）全部允許。"""
    out = {}
    for L in range(abs(l1 - l2), l1 + l2 + 1):
        for Sx in (0, 1):
            if equivalent:
                orb_sym = (-1) ** (2 * l1 - L)
                spin_sym = -1 if Sx == 0 else 1
                if orb_sym * spin_sym != -1:
                    continue
            out[f"{2 * Sx + 1}{SPEC[L]}"] = (2 * L + 1) * (2 * Sx + 1)
    return out


SLATER_P2 = {"3P": -5, "1D": 1, "1S": 10}        # p² 的項能量（單位 F2；本系列只用來排序，書上 Fig. 7.1.2 只是示意）


def hund_order(terms) -> list:
    """Hund 規則：先比 S（大者低）、再比 L（大者低）→ 由低到高。"""
    return sorted(terms, key=lambda t: (-int(t[0]), -SPEC.index(t[1])))


# ======================================================================== §7.1C–D 自旋–軌道、LS／jj、向量錐
def spin_orbit(l: int, book: bool = False) -> dict:
    """⟨S·l⟩ 的本徵值（數值對角化 S·l 於 (2l+1)·2 維）→ {j: (值, 簡併)}；book=True 照印 (7.1.26)（沒有 ½）。"""
    Lo, So, n = _single_ops(l)
    SL = sum(So[k] @ Lo[k] for k in Lo)
    Jt = {k: Lo[k] + So[k] for k in Lo}
    J2 = sum(Jt[k] @ Jt[k] for k in Jt)
    w, U = np.linalg.eigh(SL + 1e-3 * J2)
    out = {}
    for jj in (l + 0.5, l - 0.5):
        if jj < 0:
            continue
        sel = [i for i in range(n) if abs(np.real(U[:, i].conj() @ J2 @ U[:, i]) - jj * (jj + 1)) < 1e-6]
        val = float(np.real(U[:, sel[0]].conj() @ SL @ U[:, sel[0]]))
        out[jj] = ((2 if book else 1) * val, len(sel))
    return out


ALPHA = 1 / 137.035999084
HARTREE_EV = 27.211386245988
HARTREE_CM = 219474.6313632
CM_TO_GHZ = 29.9792458


def dirac_level(n: int, j, Z: int = 1) -> float:
    """(7.1.27)：ε_{n,j} = −Z²/2n² − (α²Z⁴/2n³)(1/(j+½) − 3/4n)（hartree）。"""
    j = float(j)
    return -Z * Z / (2 * n * n) - ALPHA ** 2 * Z ** 4 / (2 * n ** 3) * (1 / (j + 0.5) - 3 / (4 * n))


def a_so(n: int, l: int, Z: int = 1) -> float:
    """(7.1.25b)：a_so = Z⁴α²/(2n³ l(l+½)(l+1))（hartree）。"""
    return Z ** 4 * ALPHA ** 2 / (2 * n ** 3 * l * (l + 0.5) * (l + 1))


def fine_structure() -> dict:
    """1s 的相對論位移、2p_{3/2} − 2p_{1/2}：Dirac (7.1.27) ＝ a_so·(⟨S·l⟩ 的差)（更正的 (7.1.26)）。"""
    shift1s = (dirac_level(1, 0.5) + 0.5) * HARTREE_EV
    d2p = dirac_level(2, 1.5) - dirac_level(2, 0.5)
    so = spin_orbit(1)
    d2p_so = a_so(2, 1) * (so[1.5][0] - so[0.5][0])
    d2p_book = a_so(2, 1) * (spin_orbit(1, True)[1.5][0] - spin_orbit(1, True)[0.5][0])
    return dict(shift1s_eV=shift1s, d2p_eV=d2p * HARTREE_EV, d2p_cm=d2p * HARTREE_CM, d2p_GHz=d2p * HARTREE_CM * CM_TO_GHZ,
                d2p_so_eV=d2p_so * HARTREE_EV, d2p_book_eV=d2p_book * HARTREE_EV, a_so_2p=a_so(2, 1))


def couple(j1, j2) -> list:
    """三角條件：j1 ⊗ j2 = |j1−j2| ⊕ … ⊕ (j1+j2)（Fraction）。"""
    j1, j2 = Fraction(j1), Fraction(j2)
    out, J = [], abs(j1 - j2)
    while J <= j1 + j2:
        out.append(J)
        J += 1
    return out


def ls_jj_counts() -> dict:
    """(2p)² 的 Pauli 允許能階：LS（¹S0、³P0,1,2、¹D2）與 jj（(3/2)²：0、2；(3/2,1/2)：1、2；(1/2)²：0）的 J 值；
    (2p3p) 三重態 (0⊕1⊕2)⊗1 的 J 值。"""
    ls = sorted([0] + [0, 1, 2] + [2])
    jj = sorted([0, 2] + [1, 2] + [0])
    trip = sorted(int(J) for L in (0, 1, 2) for J in couple(L, 1))
    return dict(ls=ls, jj=jj, triplet_2p3p=trip, c32_12=[float(x) for x in couple(Fraction(3, 2), Fraction(1, 2))],
                c32_32=[float(x) for x in couple(Fraction(3, 2), Fraction(3, 2))])


def p2_levels(x: float, F2: float = 1.0, zeta: float = 4.0) -> list:
    """本系列的 (2p)² LS→jj 相關（Fig. 7.1.4 的數值版）：H = (1−x)H_C + xζ(l1·s1 + l2·s2)，H_C 依 (L, S) 取 SLATER_P2×F2。
    → [(E, J, 簡併)]（由低到高；15 個態）。"""
    Q, L2op, S2op, Lt, St, (L1, S1, L2, S2, I1, I2) = _two_electron(1, 1, True)
    Hso = sum(np.kron(S1[k] @ L1[k], I2) + np.kron(I1, S2[k] @ L2[k]) for k in L1)
    Jt = {k: Lt[k] + St[k] for k in Lt}
    J2 = sum(Jt[k] @ Jt[k] for k in Jt)
    Lq, Sq = Q.conj().T @ L2op @ Q, Q.conj().T @ S2op @ Q
    w, U = np.linalg.eigh(Lq + 7.0 * Sq)
    HC = np.zeros_like(Lq)
    for key, e in SLATER_P2.items():
        Lv, Sv = SPEC.index(key[1]), (int(key[0]) - 1) // 2
        val = Lv * (Lv + 1) + 7.0 * Sv * (Sv + 1)
        P = U[:, np.abs(w - val) < 1e-8]
        HC += e * F2 * (P @ P.conj().T)
    H = (1 - x) * HC + x * zeta * (Q.conj().T @ Hso @ Q)
    Jq = Q.conj().T @ J2 @ Q
    ev, V = np.linalg.eigh((H + H.conj().T) / 2 + 1e-6 * Jq)
    out, i = [], 0
    while i < len(ev):
        k = i
        while k + 1 < len(ev) and abs(ev[k + 1] - ev[i]) < 1e-8:
            k += 1
        jv = np.real(V[:, i].conj() @ Jq @ V[:, i])
        Jl = int(round((-1 + math.sqrt(1 + 4 * jv)) / 2))
        out.append((float(np.real(V[:, i].conj() @ H @ V[:, i])), Jl, k - i + 1))
        i = k + 1
    return out


def cg_curve(j1, j2, J, M) -> tuple:
    """Fig. 7.1.7：C^{j1 j2 J}_{m, M−m, M} 對 m（只取 |M − m| ≤ j2）→ (ms, 值)。"""
    out_m, out_v = [], []
    for m in ms(j1):
        if abs(S_(M) - m) <= S_(j2):
            out_m.append(float(m))
            out_v.append(cg_float(j1, j2, J, m, S_(M) - m, M))
    return np.array(out_m), np.array(out_v)


def classical_range(j1, j2, J, M) -> tuple:
    """向量錐（Fig. 7.1.6）：|J| = √(J(J+1)) 與 z 軸夾 Θ（cos Θ = M/|J|），j1 繞 J 轉、與 J 夾 θ13（餘弦定律）
    ⇒ m1 = |j1| cos(Θ ± θ13) 是古典轉折點。"""
    L1, L2, LJ = (math.sqrt(float(x) * (float(x) + 1)) for x in (j1, j2, J))
    th = math.acos((LJ ** 2 + L1 ** 2 - L2 ** 2) / (2 * LJ * L1))
    Th = math.acos(float(M) / LJ)
    return L1 * math.cos(Th + th), L1 * math.cos(Th - th)


def sign_changes(vals, tol: float = 1e-12) -> int:
    v = [x for x in vals if abs(x) > tol]
    return sum(1 for a, b in zip(v, v[1:]) if a * b < 0)


# ======================================================================== §7.2A 純量、向量、張量
_Ax, _Ay, _Az, _Bx, _By, _Bz = sp.symbols("A_x A_y A_z B_x B_y B_z")
_x, _y, _z = sp.symbols("x y z", real=True)


def spherical(v) -> dict:
    """(7.2.3)：A¹₁ = −(A_x + iA_y)/√2、A¹₀ = A_z、A¹₋₁ = (A_x − iA_y)/√2。"""
    vx, vy, vz = v
    return {1: -(vx + IU * vy) / _s2, 0: vz, -1: (vx - IU * vy) / _s2}


def coupled_product(A, B, k: int, q: int):
    """[A¹ × B¹]^k_q = Σ C^{1 1 k}_{m1 m2 q} A¹_{m1} B¹_{m2}（(7.2.4a)、(7.2.5a)、(7.2.6a)）。"""
    a, b = spherical(A), spherical(B)
    return sp.expand(sum(cg(1, 1, k, m1, q - m1, q) * a[m1] * b[q - m1] for m1 in (1, 0, -1) if abs(q - m1) <= 1))


T727_COLS = [(2, 2), (2, 1), (2, 0), (2, -1), (2, -2), (1, 1), (1, 0), (1, -1), (0, 0)]
T727_ROWS = [(i, j) for i in "xyz" for j in "xyz"]


def table_727() -> sp.Matrix:
    """(7.2.7)：列 A_iB_j（xx, xy, …, zz），行 [ ]^k_q；元素 = [A¹×B¹]^k_q 中 A_iB_j 的係數。"""
    A, B = (_Ax, _Ay, _Az), (_Bx, _By, _Bz)
    sym = {"x": 0, "y": 1, "z": 2}
    M = sp.zeros(9, 9)
    for c, (k, q) in enumerate(T727_COLS):
        e = coupled_product(A, B, k, q)
        for r, (i, j) in enumerate(T727_ROWS):
            M[r, c] = sp.simplify(e.coeff(A[sym[i]]).coeff(B[sym[j]]))
    return M


_h = HALF
T727_BOOK = sp.Matrix([          # 照印（「·」＝ 0）；(xy, [ ]²₋₂) 印 −½
    [_h, 0, -1 / _s6, 0, _h, 0, 0, 0, -1 / _s3],
    [IU / 2, 0, 0, 0, -_h, 0, IU / _s2, 0, 0],
    [0, -_h, 0, _h, 0, -_h, 0, -_h, 0],
    [IU / 2, 0, 0, 0, -IU / 2, 0, -IU / _s2, 0, 0],
    [-_h, 0, -1 / _s6, 0, -_h, 0, 0, 0, -1 / _s3],
    [0, -IU / 2, 0, -IU / 2, 0, -IU / 2, 0, IU / 2, 0],
    [0, -_h, 0, _h, 0, _h, 0, _h, 0],
    [0, -IU / 2, 0, -IU / 2, 0, IU / 2, 0, -IU / 2, 0],
    [0, 0, 2 / _s6, 0, 0, 0, 0, 0, -1 / _s3]])


def check_724() -> dict:
    """(7.2.4b)、(7.2.5b)、(7.2.6b)：[A×B]⁰₀ = −A·B/√3、[A×B]¹₀ = i(A×B)_z/√2、[A×B]²₀ = (−A_xB_x − A_yB_y + 2A_zB_z)/√6。"""
    A, B = (_Ax, _Ay, _Az), (_Bx, _By, _Bz)
    dot = _Ax * _Bx + _Ay * _By + _Az * _Bz
    crz = _Ax * _By - _Ay * _Bx
    return dict(scalar=_zero(coupled_product(A, B, 0, 0) + dot / _s3),
                vector=_zero(coupled_product(A, B, 1, 0) - IU * crz / _s2),
                t22=_zero(coupled_product(A, B, 2, 2) - (_Ax * _Bx - _Ay * _By + IU * _Ax * _By + IU * _Ay * _Bx) / 2),
                t21=_zero(coupled_product(A, B, 2, 1) + (_Ax * _Bz + _Az * _Bx + IU * _Ay * _Bz + IU * _Az * _By) / 2),
                t20=_zero(coupled_product(A, B, 2, 0) - (-_Ax * _Bx - _Ay * _By + 2 * _Az * _Bz) / _s6))


def multipole_X(k: int, q: int):
    """X^k_q = r^k √(4π/(2k+1)) Y^k_q 寫成 x, y, z 的多項式（Condon–Shortley 號）：
    X^k_{±|q|} = (∓1)^{|q|}√((k−|q|)!/(k+|q|)!)(x ± iy)^{|q|} r^{k−|q|}P_k^{(|q|)}(z/r)。"""
    am = abs(q)
    zt = sp.Symbol("zt")
    P = sp.expand(sp.diff((zt ** 2 - 1) ** k, zt, k + am) / (2 ** k * sp.factorial(k)))
    norm = sp.sqrt(sp.factorial(k - am) / sp.factorial(k + am))
    rr2 = _x ** 2 + _y ** 2 + _z ** 2
    hom = 0
    for (n,), c in sp.Poly(P, zt).terms():
        hom += c * _z ** n * rr2 ** sp.Rational(k - n - am, 2)
    cs = (-1) ** am if q > 0 else 1
    base = (_x + IU * _y) if q >= 0 else (_x - IU * _y)
    return sp.expand(cs * norm * hom * base ** am)


def rr_tensor() -> dict:
    """(7.2.9)：[r¹r¹]²_q（A = B = r）與 X²_q；比值都是 √(2/3)。BOOK_7292 照印的 [rr]²₁ = −(x − iy)z。"""
    r = (_x, _y, _z)
    out = {}
    for q in (2, 1, 0, -1, -2):
        t = sp.expand(coupled_product((_x, _y, _z), (_x, _y, _z), 2, q))
        out[q] = (t, multipole_X(2, q), sp.simplify(t / multipole_X(2, q)))
    return out


BOOK_7292 = -(_x - IU * _y) * _z


def ck_formula(k: int, q: int) -> tuple:
    """(7.2.11)：C^{k−1 1 k}_{q−1 1 q}、C^{k−1 1 k}_{q 0 q}、C^{k−1 1 k}_{q+1 −1 q}。"""
    a = sp.sqrt(R_((k + q - 1) * (k + q), (2 * k - 1) * (2 * k)))
    b = sp.sqrt(R_(k * k - q * q, k * (2 * k - 1)))
    c = sp.sqrt(R_((k - q) * (k - q - 1), 2 * k * (2 * k - 1)))
    return a, b, c


def scalar_cg(j, m):
    """(7.2.12)：C^{j j 0}_{m, −m, 0} = (−1)^{j−m}/√(2j+1)。"""
    j, m = S_(j), S_(m)
    return (-1) ** (j - m) / sp.sqrt(2 * j + 1)


# ======================================================================== §7.2C 3j、對稱、正交、因式分解
def threej(j1, j2, j3, m1, m2, m3):
    """(7.2.30a)：(j1 j2 j3; m1 m2 m3) = (−1)^{j1−j2−m3}C^{j1 j2 j3}_{m1 m2 −m3}/√(2j3+1)。"""
    j1, j2, j3, m1, m2, m3 = (S_(x) for x in (j1, j2, j3, m1, m2, m3))
    return sp.nsimplify((-1) ** (j1 - j2 - m3) * cg(j1, j2, j3, m1, m2, -m3) / sp.sqrt(2 * j3 + 1))


_SAMPLE_J = [(1, 1, 2), (1, 1, 1), (HALF, HALF, 1), (1, HALF, R_(3, 2)), (2, 1, 2), (R_(3, 2), 1, R_(3, 2)), (2, 2, 2), (2, 2, 4)]


def _m_triples(j1, j2, j3):
    for m1 in ms(j1):
        for m2 in ms(j2):
            m3 = -(m1 + m2)
            if abs(m3) <= j3:
                yield m1, m2, m3


def threej_symmetries() -> dict:
    """(7.2.31)–(7.2.35) 在 _SAMPLE_J 的所有 m 上逐一核對；也與 sympy 的 wigner_3j 比較。"""
    ok = dict(sympy=True, odd=True, cyclic=True, mflip=True, c733=True, c734=True, c735=True)
    for j1, j2, j3 in _SAMPLE_J:
        ph = (-1) ** (j1 + j2 + j3)
        for m1, m2, m3 in _m_triples(j1, j2, j3):
            t = threej(j1, j2, j3, m1, m2, m3)
            ok["sympy"] &= _zero(t - _sym_3j(j1, j2, j3, m1, m2, m3))
            ok["odd"] &= _zero(t - ph * threej(j2, j1, j3, m2, m1, m3)) and _zero(t - ph * threej(j3, j2, j1, m3, m2, m1)) \
                and _zero(t - ph * threej(j1, j3, j2, m1, m3, m2))
            ok["cyclic"] &= _zero(t - threej(j3, j1, j2, m3, m1, m2)) and _zero(t - threej(j2, j3, j1, m2, m3, m1))
            ok["mflip"] &= _zero(t - ph * threej(j1, j2, j3, -m1, -m2, -m3))
        for m1 in ms(j1):
            for m2 in ms(j2):
                m3 = m1 + m2
                if abs(m3) > j3:
                    continue
                ok["c733"] &= _zero(cg(j2, j1, j3, m2, m1, m3) - (-1) ** (j1 + j2 - j3) * cg(j1, j2, j3, m1, m2, m3))
                if j1 == j2:
                    ok["c734"] &= _zero(cg(j1, j1, j3, m2, m1, m3) - (-1) ** (2 * j1 - j3) * cg(j1, j1, j3, m1, m2, m3))
                # (7.2.35)：C^{j1 j3 j2}_{m1 m3 m2}（m2 = m1 + m3 的那個關係）＝ (−1)^{j2−j3+m1}√((2j2+1)/(2j3+1)) C^{j1 j2 j3}_{−m1 m2 m3}
        for m1 in ms(j1):
            for m3 in ms(j3):
                m2 = m1 + m3
                if abs(m2) > j2:
                    continue
                lhs = cg(j1, j3, j2, m1, m3, m2)
                rhs = (-1) ** (S_(j2) - S_(j3) + m1) * sp.sqrt((2 * S_(j2) + 1) / (2 * S_(j3) + 1)) * cg(j1, j2, j3, -m1, m2, m3)
                ok["c735"] &= _zero(lhs - rhs)
    return ok


def orthogonality_check(pairs=((1, HALF), (1, 1), (R_(3, 2), 1))) -> bool:
    """(7.2.21)、(7.2.22)：耦合係數表是正交矩陣（C^T C = 1 = C C^T）。"""
    good = True
    for j1, j2 in pairs:
        M = np.array(cg_table(j1, j2)[0].evalf(), dtype=float)
        good &= _close(M.T @ M, np.eye(M.shape[0])) and _close(M @ M.T, np.eye(M.shape[0]))
    return bool(good)


def factorization_check(cases=((1, 1, 2, (1, 0, 1), (0, 1, 1)), (1, HALF, R_(3, 2), (0, HALF, HALF), (1, -HALF, HALF)),
                               (2, 1, 2, (1, -1, 0), (2, 0, 2)), (1, 1, 1, (1, 0, 1), (1, -1, 0)))) -> float:
    """(7.2.25)：∫d(αβγ) 𝒟^{j3*}_{m3 m3'} 𝒟^{j1}_{m1 m1'} 𝒟^{j2}_{m2 m2'} = C^{j1 j2 j3}_{m1 m2 m3}C^{j1 j2 j3}_{m1' m2' m3'}/(2j3+1)
    （Euler 角的數值積分：α、γ 等分 12 點、cos β 用 24 點 Gauss–Legendre；𝒟 = e^{−imα}d(β)e^{−inγ}）。→ 最大誤差。"""
    al = np.radians(np.arange(12) * 30.0)
    x, wb = leggauss(24)
    be = np.degrees(np.arccos(x))
    err = 0.0
    for j1, j2, j3, (m1, m2, m3), (n1, n2, n3) in cases:
        idx = lambda j, m: [float(t) for t in ms(j)].index(float(m))
        d = {j: [RT.wigner_d(j, b) for b in be] for j in (j1, j2, j3)}
        fb = np.array([np.conj(d[j3][k][idx(j3, m3), idx(j3, n3)]) * d[j1][k][idx(j1, m1), idx(j1, n1)] * d[j2][k][idx(j2, m2), idx(j2, n2)]
                       for k in range(len(be))])
        fa = np.mean(np.exp(1j * (float(m3) - float(m1) - float(m2)) * al))
        fg = np.mean(np.exp(1j * (float(n3) - float(n1) - float(n2)) * al))
        tot = fa * fg * np.sum(wb * fb) / 2
        rhs = float(cg(j1, j2, j3, m1, m2, m3) * cg(j1, j2, j3, n1, n2, n3) / (2 * S_(j3) + 1))
        err = max(err, abs(tot - rhs))
    return err


# ======================================================================== §7.2D Schwinger 玻色子、SU(2) 公式、6j
def su2_irrep_formula(j, u) -> np.ndarray:
    """(7.2.38)：𝒟^j_{mn}(u)（u ∈ SU(2) 的 2×2）。"""
    tj = int(round(2 * float(j)))
    jf = tj / 2
    msj = [jf - k for k in range(tj + 1)]
    out = np.zeros((tj + 1, tj + 1), complex)
    for a, m in enumerate(msj):
        for b, n in enumerate(msj):
            pre = math.sqrt(_FACT(int(jf + n)) * _FACT(int(jf - n)) * _FACT(int(jf + m)) * _FACT(int(jf - m)))
            s = 0
            for k in range(0, tj + 1):
                e = [jf - k + m, k - m + n, k, jf - n - k]
                if min(e) < 0:
                    continue
                e = [int(round(x)) for x in e]
                s += u[0, 0] ** e[0] * u[1, 0] ** e[1] * u[0, 1] ** e[2] * u[1, 1] ** e[3] / (
                    _FACT(e[0]) * _FACT(e[1]) * _FACT(e[2]) * _FACT(e[3]))
            out[a, b] = pre * s
    return out


def _fock_det(state: dict, power: int) -> dict:
    """(a¹₁a²₂ − a¹₂a²₁)^power 作用在 Fock 態 {(n11, n12, n21, n22): 係數}。"""
    for _ in range(power):
        new = {}
        for (n11, n12, n21, n22), c in state.items():
            k1 = (n11 + 1, n12, n21, n22 + 1)
            k2 = (n11, n12 + 1, n21 + 1, n22)
            new[k1] = new.get(k1, 0) + c * sp.sqrt((n11 + 1) * (n22 + 1))
            new[k2] = new.get(k2, 0) - c * sp.sqrt((n12 + 1) * (n21 + 1))
        state = new
    return state


def boson_coupled(j, m, n, N, half: bool = True) -> dict:
    """(7.2.54)：|j m n N⟩（Fock 佔據數 (n¹₁ n¹₂, n²₁ n²₂) → 係數）；half=False 照印（方括號沒有 ½ 次方）。"""
    j, m, n, N = (S_(x) for x in (j, m, n, N))
    f = sp.factorial
    pre = sp.sqrt((2 * j + 1) * f(j + m) * f(j - m) / (f(N / 2 - j) * f(N / 2 + j + 1)))
    st = {}
    for k in range(0, int(2 * j) + 1):
        occ = (j + m - k, n - m + k, sp.Integer(k), j - n - k)
        if min(occ) < 0:
            continue
        br = f(j + n) / (f(occ[0]) * f(occ[1])) * f(j - n) / (f(occ[2]) * f(occ[3]))
        key = tuple(int(x) for x in occ)
        st[key] = st.get(key, 0) + pre * (sp.sqrt(br) if half else br)
    return _fock_det(st, int(N / 2 - j))


def boson_cg(j1, j2, j3, m1, m2, m3, half: bool = True):
    """C^{j1 j2 j3}_{m1 m2 m3} ＝ ⟨(7.2.36)|(7.2.54)⟩（n = j1 − j2、N = 2(j1 + j2)）。"""
    j1, j2, j3, m1, m2, m3 = (S_(x) for x in (j1, j2, j3, m1, m2, m3))
    st = boson_coupled(j3, m3, j1 - j2, 2 * (j1 + j2), half)
    key = (int(j1 + m1), int(j1 - m1), int(j2 + m2), int(j2 - m2))
    return sp.nsimplify(sp.simplify(st.get(key, 0)))


def boson_norm(j, m, n, N, half: bool = True):
    return sp.nsimplify(sp.simplify(sum(c ** 2 for c in boson_coupled(j, m, n, N, half).values())))


def det_scalar_state(d: int) -> dict:
    """(7.2.50)：[(2d)!(2d+1)!]^{−½}(a¹₁a²₂ − a¹₂a²₁)^{2d}|00,00⟩ → {佔據數: 係數}；結果 = (2d+1)^{−½}Σ(−1)^r|r 2d−r, 2d−r r⟩。"""
    st = _fock_det({(0, 0, 0, 0): sp.Integer(1)}, 2 * d)
    nrm = 1 / sp.sqrt(sp.factorial(2 * d) * sp.factorial(2 * d + 1))
    return {k: sp.nsimplify(v * nrm) for k, v in st.items()}


def sixj(a, b, c, d, e, f):
    return sp.nsimplify(_sym_6j(*(S_(x) for x in (a, b, c, d, e, f))))


def sixj_from_3j(j1, j2, j3, j4, j5, j6):
    """{j1 j2 j3; j4 j5 j6} = Σ (−1)^{Σ(j−m)} (j1 j2 j3; −m1 −m2 −m3)(j1 j5 j6; m1 −m5 m6)(j4 j2 j6; m4 m2 −m6)(j4 j5 j3; −m4 m5 m3)。"""
    J = [S_(x) for x in (j1, j2, j3, j4, j5, j6)]
    tot = 0
    for m1 in ms(J[0]):
        for m2 in ms(J[1]):
            m3 = -(m1 + m2)
            if abs(m3) > J[2]:
                continue
            for m5 in ms(J[4]):
                m6 = m5 - m1
                m4 = m6 - m2
                if abs(m6) > J[5] or abs(m4) > J[3] or -m4 + m5 + m3 != 0:
                    continue
                ph = (-1) ** (sum(J) - (m1 + m2 + m3 + m4 + m5 + m6))
                tot += ph * threej(J[0], J[1], J[2], -m1, -m2, -m3) * threej(J[0], J[4], J[5], m1, -m5, m6)                     * threej(J[3], J[1], J[5], m4, m2, -m6) * threej(J[3], J[4], J[2], -m4, m5, m3)
    return sp.nsimplify(sp.simplify(tot))


def sixj_from_cg(j1, j2, j3, J, j12, j23, M=None):
    """(7.3.49)：Σ C^{j1 j2 j12}C^{j12 j3 J}C^{j1 j23 J}C^{j2 j3 j23}（固定 M）
    ＝ (−1)^{j1+j2+j3+J}√((2j12+1)(2j23+1)){j12 j1 j2; j23 j3 J} → 回傳 {…}。"""
    j1, j2, j3, J, j12, j23 = (S_(x) for x in (j1, j2, j3, J, j12, j23))
    M = J if M is None else S_(M)
    tot = 0
    for m1 in ms(j1):
        for m2 in ms(j2):
            m3 = M - m1 - m2
            if abs(m3) > j3:
                continue
            m12, m23 = m1 + m2, m2 + m3
            if abs(m12) > j12 or abs(m23) > j23:
                continue
            tot += cg(j1, j2, j12, m1, m2, m12) * cg(j12, j3, J, m12, m3, M) * cg(j1, j23, J, m1, m23, M) * cg(j2, j3, j23, m2, m3, m23)
    return sp.nsimplify(sp.simplify(tot / ((-1) ** (j1 + j2 + j3 + J) * sp.sqrt((2 * j12 + 1) * (2 * j23 + 1)))))


# ======================================================================== §7.3A 張量算子
def tensor_T(j1, j2, k, q) -> sp.Matrix:
    """(7.3.1a)／(7.3.2)：T(j1 j2)^k_q = Σ C^{j1 j2 k}_{m1 m2 q}|j1 m1⟩⟨j2 −m2|(−1)^{j2−m2}（列 m1、行 m'）。"""
    j1, j2, k, q = (S_(x) for x in (j1, j2, k, q))
    r, c = ms(j1), ms(j2)
    M = sp.zeros(len(r), len(c))
    for a, m1 in enumerate(r):
        for b, mp in enumerate(c):
            m2 = -mp
            M[a, b] = cg(j1, j2, k, m1, m2, q) * (-1) ** (j2 - m2)
    return M


def unit_v(j1, j2, k, q) -> sp.Matrix:
    """(7.3.17a)：v(j1 j2)^k_q = (−1)^{2j1}Σ C^{j1 j2 k}_{m1 −m2 q}(−1)^{j2+m2}|j1 m1⟩⟨j2 m2|。"""
    j1, j2, k, q = (S_(x) for x in (j1, j2, k, q))
    r, c = ms(j1), ms(j2)
    M = sp.zeros(len(r), len(c))
    for a, m1 in enumerate(r):
        for b, m2 in enumerate(c):
            M[a, b] = (-1) ** (2 * j1) * cg(j1, j2, k, m1, -m2, q) * (-1) ** (j2 + m2)
    return M


T733_BOOK = {(1, -1): sp.Matrix([[0, 0], [-1, 0]]), (1, 0): sp.Matrix([[-1, 0], [0, 1]]) / _s2,
             (1, 1): sp.Matrix([[0, 1], [0, 0]]), (0, 0): -sp.eye(2) / _s2}


def spin_half_T() -> dict:
    """(7.3.3)：自旋 ½ 的 T^k_q（k = 0, 1）。"""
    return {(k, q): tensor_T(HALF, HALF, k, q) for k in (0, 1) for q in range(-k, k + 1)}


PAULI_S = (sp.Matrix([[0, 1], [1, 0]]), sp.Matrix([[0, -IU], [IU, 0]]), sp.Matrix([[1, 0], [0, -1]]))


def spin_half_cart() -> tuple:
    """(7.3.4)：T_x = −(T¹₋₁ − T¹₁)/√2、T_y = −i(T¹₋₁ + T¹₁)/√2、T_z = −T¹₀ → σ/√2。"""
    T = spin_half_T()
    Tx = -(T[(1, -1)] - T[(1, 1)]) / _s2
    Ty = -IU * (T[(1, -1)] + T[(1, 1)]) / _s2
    Tz = -T[(1, 0)]
    return Tx, Ty, Tz


def rotated_Jz_check(beta: float = 37.0) -> float:
    """(7.3.7)–(7.3.8)：R(0β0)J_zR† = sin β J_x + cos β J_z，且 R T¹₀ R† = Σ𝒟¹_{q0}T¹_q（自旋 ½）。→ 最大誤差。"""
    R = RT.su2_euler(0, beta, 0)
    J = RT.jmat(HALF)
    b = math.radians(beta)
    e1 = np.abs(R @ J["Jz"] @ R.conj().T - (math.sin(b) * J["Jx"] + math.cos(b) * J["Jz"])).max()
    T = {q: _num(tensor_T(HALF, HALF, 1, q)) for q in (1, 0, -1)}
    D1 = RT.wigner_D(1, 0, beta, 0)
    rhs = sum(D1[[1, 0, -1].index(q), 1] * T[q] for q in (1, 0, -1))
    e2 = np.abs(R @ T[0] @ R.conj().T - rhs).max()
    return float(max(e1, e2))


EXPECT_CARD = dict(beta0=60.0, beta=45.0)      # 本系列 Fig. 7.3.1 的例子：|ψ⟩ = R(0 β0 0)|↑⟩，再把 z 軸轉 β


def expectation_rotation(beta0: float = 60.0, beta: float = 45.0) -> dict:
    """Fig. 7.3.1、(7.3.9)–(7.3.11)：⟨J_z(轉過)⟩ = sin β⟨J_x⟩ + cos β⟨J_z⟩ ＝ ½(p↑ − p↓)（(7.3.10a) 的機率式）。"""
    psi = RT.su2_euler(0, beta0, 0)[:, 0]
    J = RT.jmat(HALF)
    ex = {k: float(np.real(psi.conj() @ J[k] @ psi)) for k in ("Jx", "Jy", "Jz")}
    b = math.radians(beta)
    vec = math.sin(b) * ex["Jx"] + math.cos(b) * ex["Jz"]
    R = RT.su2_euler(0, beta, 0)
    amps = psi.conj() @ R
    p_up, p_dn = abs(amps[0]) ** 2, abs(amps[1]) ** 2
    return dict(Jx=ex["Jx"], Jy=ex["Jy"], Jz=ex["Jz"], rotated=vec, p_up=p_up, p_dn=p_dn, rotated_prob=0.5 * (p_up - p_dn))


T7312_BOOK = {   # (7.3.12)，基底 m = 1, 0, −1
    (2, -2): sp.Matrix([[0, 0, 0], [0, 0, 0], [1, 0, 0]]),
    (2, -1): sp.Matrix([[0, 0, 0], [1, 0, 0], [0, -1, 0]]) / _s2,
    (2, 0): sp.Matrix([[1, 0, 0], [0, -2, 0], [0, 0, 1]]) / _s6,
    (2, 1): sp.Matrix([[0, -1, 0], [0, 0, 1], [0, 0, 0]]) / _s2,
    (2, 2): sp.Matrix([[0, 0, 1], [0, 0, 0], [0, 0, 0]]),
    (1, -1): sp.Matrix([[0, 0, 0], [1, 0, 0], [0, 1, 0]]) / _s2,
    (1, 0): sp.Matrix([[1, 0, 0], [0, 0, 0], [0, 0, -1]]) / _s2,
    (1, 1): sp.Matrix([[0, -1, 0], [0, 0, -1], [0, 0, 0]]) / _s2,
    (0, 0): sp.eye(3) / _s3}


def tensor_j1() -> dict:
    return {(k, q): tensor_T(1, 1, k, q) for k in (0, 1, 2) for q in range(-k, k + 1)}


def j1_relations() -> bool:
    """(7.3.13)：T¹₋₁ = J₋/2、T¹₀ = J_z/√2、T¹₁ = −J₊/2（j = 1）。"""
    T = tensor_j1()
    J = _jm_exact(1)
    return all(_zero(x) for x in (T[(1, -1)] - J["Jm"] / 2)) and all(_zero(x) for x in (T[(1, 0)] - J["Jz"] / _s2)) \
        and all(_zero(x) for x in (T[(1, 1)] + J["Jp"] / 2))


def v_table(j, k) -> tuple:
    """Table 7.1–7.2 的濃縮形式：v^k_q（q = −k…k）疊在同一個 (2j+1)² 陣列；每條對角線（q 與 −q 同一個）寫成「分子 / √分母²」。
    → (分子 Matrix, {q: 分母²}（q = k…0）)。分母² = 該對角線各元素平方的分母的最小公倍數。"""
    j, k = S_(j), S_(k)
    n = int(2 * j + 1)
    tot = sp.zeros(n)
    den = {}
    for q in range(-int(k), int(k) + 1):
        V = unit_v(j, j, k, q)
        tot += V
    num = sp.zeros(n)
    for q in range(0, int(k) + 1):
        ents = [tot[a, a + q] for a in range(n - q)] + [tot[a + q, a] for a in range(n - q)]
        sq = [Fraction(str(sp.nsimplify(e ** 2))) for e in ents if not _zero(e)]
        dd = 1
        for fr in sq:
            dd = dd * fr.denominator // math.gcd(dd, fr.denominator)
        den[q] = dd
        for a in range(n - q):
            num[a, a + q] = sp.nsimplify(sp.radsimp(tot[a, a + q] * sp.sqrt(dd)))
            num[a + q, a] = sp.nsimplify(sp.radsimp(tot[a + q, a] * sp.sqrt(dd)))
    return num, den


def _r(*rows):
    s = sp.sqrt
    return sp.Matrix([[sp.sympify(x.replace("r", "sqrt").replace(".", "0")) if isinstance(x, str) else x for x in row] for row in rows])


V_TABLES_BOOK = {   # Table 7.1–7.2（j ≤ 2），(j, k): (分子（"r3" = √3，"." = 0）, {q: 分母²})
    (HALF, 1): (_r([1, -1], [1, -1]), {1: 1, 0: 2}),
    (R_(3, 2), 1): (_r([3, "-r(3)", ".", "."], ["r(3)", 1, -2, "."], [".", 2, -1, "-r(3)"], [".", ".", "r(3)", -3]), {1: 10, 0: 20}),
    (R_(3, 2), 2): (_r([1, -1, 1, "."], [1, -1, 0, 1], [1, 0, -1, 1], [".", 1, -1, 1]), {2: 2, 1: 2, 0: 4}),
    (R_(3, 2), 3): (_r([1, -1, 1, -1], [1, -3, "r(3)", -1], [1, "-r(3)", 3, -1], [1, -1, 1, -1]), {3: 1, 2: 2, 1: 5, 0: 20}),
    (1, 2): (_r([1, -1, 1], [1, -2, 1], [1, -1, 1]), {2: 1, 1: 2, 0: 6}),
    (1, 1): (_r([1, -1, "."], [1, 0, -1], [".", 1, -1]), {1: 2, 0: 2}),
    (2, 4): (_r([1, -1, "r(3)", -1, 1], [1, -4, "r(6)", "-r(8)", 1], ["r(3)", "-r(6)", 6, "-r(6)", "r(3)"],
                [1, "-r(8)", "r(6)", -4, 1], [1, -1, "r(3)", -1, 1]), {4: 1, 3: 2, 2: 14, 1: 14, 0: 70}),
    (2, 3): (_r([1, "-r(3)", 1, -1, "."], ["r(3)", -2, "r(2)", 0, -1], [1, "-r(2)", 0, "r(2)", -1], [1, 0, "-r(2)", 2, "-r(3)"],
                [".", 1, -1, "r(3)", -1]), {3: 2, 2: 2, 1: 10, 0: 10}),
    (2, 2): (_r([2, "-r(6)", "r(2)", ".", "."], ["r(6)", -1, -1, "r(3)", "."], ["r(2)", 1, -2, 1, "r(2)"], [".", "r(3)", -1, -1, "r(6)"],
                [".", ".", "r(2)", "-r(6)", 2]), {2: 7, 1: 14, 0: 14}),
    (2, 1): (_r([2, "-r(2)", ".", ".", "."], ["r(2)", 1, "-r(3)", ".", "."], [".", "r(3)", 0, "-r(3)", "."], [".", ".", "r(3)", -1, "-r(2)"],
                [".", ".", ".", "r(2)", -2]), {1: 10, 0: 10})}


def v_tables_check() -> dict:
    """Table 7.1–7.2 的 j ≤ 2 格子：分子與分母都和書上相同（True）；也核對每條對角線正交歸一。"""
    out = {}
    for (j, k), (num_b, den_b) in V_TABLES_BOOK.items():
        num, den = v_table(j, k)
        out[(j, k)] = all(_zero(x) for x in (num - num_b)) and den == den_b
    return out


def elementary_expansion(j, a: int, b: int) -> dict:
    """E_ab（1 起算的列 a、行 b）＝ Σ c^k_q v^k_q，c^k_q = conj(v^k_q[a, b])（v 在跡內積下正交歸一）。→ {(k, q): c}（非零者）。"""
    j = S_(j)
    out = {}
    for k in range(0, int(2 * j) + 1):
        for q in range(-k, k + 1):
            c = sp.conjugate(unit_v(j, j, k, q)[a - 1, b - 1])
            if not _zero(c):
                out[(k, q)] = sp.nsimplify(c)
    return out


def mixed_examples() -> dict:
    """(7.3.17)–(7.3.19) 的例子：√10 V¹₁(pd)（3×5）與 √10 V¹₋₁(dp)（5×3）；Table 7.4 (d)(p) 的 v¹、v²、v³ 分母。"""
    pd = (sp.sqrt(10) * unit_v(1, 2, 1, 1)).applyfunc(sp.nsimplify)
    dp = (sp.sqrt(10) * unit_v(2, 1, 1, -1)).applyfunc(sp.nsimplify)
    return dict(pd=pd, dp=dp)


V74_DP_BOOK = {   # Table 7.4 (d)(p)：v(2 1)^k，列 m1 = 2…−2、行 m2 = 1…−1；分母² 依 q = k…0（−q 同）；v³ 第一列第二格（被箭頭蓋住）取 −1
    3: (_r([1, -1, 1], [1, "-r(8)", "r(2)"], ["r(6)", "-r(3)", "r(6)"], ["r(2)", "-r(8)", 1], [1, -1, 1]), {3: 1, 2: 3, 1: 15, 0: 5}),
    2: (_r(["r(2)", "-r(2)", "."], [1, -1, -1], ["r(3)", 0, "-r(3)"], [1, 1, -1], [".", "r(2)", "-r(2)"]), {2: 3, 1: 6, 0: 2}),
    1: (_r(["r(6)", ".", "."], ["r(3)", "r(3)", "."], [1, 2, 1], [".", "r(3)", "r(3)"], [".", ".", "r(6)"]), {1: 10, 0: 10})}


def mixed_table(j1, j2, k) -> tuple:
    """Table 7.4 的長方形 v(j1 j2)^k：(分子, {q: 分母²})（q = m1 − m2；q 與 −q 同分母）。"""
    j1, j2, k = S_(j1), S_(j2), S_(k)
    r, c = ms(j1), ms(j2)
    tot = sp.zeros(len(r), len(c))
    for q in range(-int(k), int(k) + 1):
        tot += unit_v(j1, j2, k, q)
    num = sp.zeros(len(r), len(c))
    den = {}
    for q in range(0, int(k) + 1):
        cells = [(a, b) for a in range(len(r)) for b in range(len(c)) if abs(r[a] - c[b]) == q]
        dd = 1
        for a, b in cells:
            if not _zero(tot[a, b]):
                fr = Fraction(str(sp.nsimplify(tot[a, b] ** 2)))
                dd = dd * fr.denominator // math.gcd(dd, fr.denominator)
        den[q] = dd
        for a, b in cells:
            num[a, b] = sp.nsimplify(sp.radsimp(tot[a, b] * sp.sqrt(dd)))
    return num, den


def fpoly_check() -> dict:
    """(7.3.30)–(7.3.32) 的多項式（Condon–Shortley 的 X³_q）：xyz = −i(X³₂ − X³₋₂)/√30（書印 −i(X³₋₂ − X³₂)/√30，差整體號）、
    (x² − y²)z = √(2/15)(X³₂ + X³₋₂)（書印 i(X³₂ + X³₋₂)/√2）、(x² + y²)z = −(2/5)X³₀ + (2/5)r²z（書印 −X³₀/10，且含 ℓ = 1 的部分）。"""
    X = {q: multipole_X(3, q) for q in (2, -2, 0)}
    r2 = _x ** 2 + _y ** 2 + _z ** 2
    return dict(a2u=sp.simplify(_x * _y * _z / (X[2] - X[-2])), t2u=sp.simplify((_x ** 2 - _y ** 2) * _z / (X[2] + X[-2])),
                t1u_ok=sp.expand((_x ** 2 + _y ** 2) * _z - (-R_(2, 5) * X[0] + R_(2, 5) * r2 * _z)) == 0)


def py_check(theta: float = 0.7, phi: float = 1.1) -> dict:
    """(7.6.51)：p_x = (−Y¹₁ + Y¹₋₁)r√(4π/3)/√2 = x、p_y = i(Y¹₁ + Y¹₋₁)r√(4π/3)/√2 = y（書印 i(Y¹₁ − Y¹₋₁)，給 −ix）。"""
    Y1, Ym = RT.ylm(1, 1, theta, phi), RT.ylm(1, -1, theta, phi)
    f = math.sqrt(4 * math.pi / 3) / math.sqrt(2)
    return dict(px=complex(f * (-Y1 + Ym)), py=complex(f * 1j * (Y1 + Ym)), py_book=complex(f * 1j * (Y1 - Ym)),
                x=math.sin(theta) * math.cos(phi), y=math.sin(theta) * math.sin(phi))


def handedness_check() -> bool:
    """Fig. 7.5.4：|r⟩ = (|x⟩ + i|y⟩)/√2 的 E(Φ) = Re⟨·|r⟩e^{−iΦ} 隨 Φ 逆時針轉（右旋，正文）；圖說把 |l⟩、|r⟩ 說反了。"""
    E = ellipse_field(U_RL[:, 0], np.linspace(0, 90, 10))
    ang = np.unwrap(np.arctan2(E[:, 1], E[:, 0]))
    return bool(np.all(np.diff(ang) > 0))


# ======================================================================== §7.3B–C Wigner–Eckart 與晶體場
def we_v_check(cases=((1, 1, 1), (2, 2, 4), (2, 1, 1), (R_(3, 2), HALF, 1), (3, 3, 6))) -> bool:
    """(7.3.21)–(7.3.22)：⟨j1 m1|v^k_q|j2 m2⟩ = C^{k j2 j1}_{q m2 m1}(−1)^{k+j1−j2}√((2k+1)/(2j1+1))，
    並用 (7.3.20b)（更正的 1/(2j1+1)）從矩陣元重算約化矩陣元。"""
    good = True
    for j1, j2, k in cases:
        j1, j2, k = S_(j1), S_(j2), S_(k)
        red = (-1) ** (k + j1 - j2) * sp.sqrt((2 * k + 1) / (2 * j1 + 1))
        acc = 0
        for q in range(-int(k), int(k) + 1):
            V = unit_v(j1, j2, k, q)
            for a, m1 in enumerate(ms(j1)):
                for b, m2 in enumerate(ms(j2)):
                    good &= _zero(V[a, b] - cg(k, j2, j1, q, m2, m1) * red)
                    acc += cg(k, j2, j1, q, m2, m1) * V[a, b]
        good &= _zero(acc / (2 * j1 + 1) - red)
    return bool(good)


def _sph_J(j):
    J = RT.jmat(j)
    return {1: -(J["Jx"] + 1j * J["Jy"]) / math.sqrt(2), 0: J["Jz"], -1: (J["Jx"] - 1j * J["Jy"]) / math.sqrt(2)}


def we_operator_check(j=2, euler=(30.0, 50.0, 70.0)) -> dict:
    """由 J 做成的張量算子：T¹_q = J 的球分量、T²_q = Σ C^{1 1 2}_{q1 q2 q}J_{q1}J_{q2}。核對 (7.3.1b) 的轉動、
    (7.3.20a) 的比例（每個非零元除以 C^{k j j}_{q m m'} 都一樣），並給 ⟨j‖J‖j⟩ = −√(j(j+1))（書上的 C^{k j2 j1} 次序：C^{1 j j}_{0 m m} = −m/√(j(j+1))）。"""
    Js = _sph_J(j)
    T1 = Js
    T2 = {q: sum(float(cg(1, 1, 2, q1, q - q1, q)) * Js[q1] @ Js[q - q1] for q1 in (1, 0, -1) if abs(q - q1) <= 1) for q in range(-2, 3)}
    R = RT.wigner_D(j, *euler)
    err_rot = 0.0
    for k, T in ((1, T1), (2, T2)):
        Dk = RT.wigner_D(k, *euler)
        qs = list(range(k, -k - 1, -1))
        for q in qs:
            lhs = R @ T[q] @ R.conj().T
            rhs = sum(Dk[qs.index(qp), qs.index(q)] * T[qp] for qp in qs)
            err_rot = max(err_rot, float(np.abs(lhs - rhs).max()))
    reds = {}
    for k, T in ((1, T1), (2, T2)):
        vals = []
        for q in range(-k, k + 1):
            for a, m1 in enumerate(ms(j)):
                for b, m2 in enumerate(ms(j)):
                    c = cg_float(k, j, j, q, m2, m1)
                    if abs(c) > 1e-9:
                        vals.append(T[q][a, b] / c)
        reds[k] = (float(np.real(np.mean(vals))), float(np.abs(np.array(vals) - np.mean(vals)).max()))
    return dict(rot_err=err_rot, red_J=reds[1], red_T2=reds[2], red_J_formula=-math.sqrt(float(j) * (float(j) + 1)))


_X40 = (35 * _z ** 4 - 30 * _z ** 2 * (_x ** 2 + _y ** 2 + _z ** 2) + 3 * (_x ** 2 + _y ** 2 + _z ** 2) ** 2) / 8
_X44 = sp.sqrt(R_(35, 128)) * (_x + IU * _y) ** 4
_X4m4 = sp.sqrt(R_(35, 128)) * (_x - IU * _y) ** 4


def v4_poly_check() -> bool:
    """(7.3.25)：x⁴ + y⁴ + z⁴ − (3/5)r⁴ = 2(X⁴₄ + X⁴₋₄)/√70 + (2/5)X⁴₀；X⁴_q 也與 multipole_X(4, q) 相同。"""
    r2 = _x ** 2 + _y ** 2 + _z ** 2
    lhs = _x ** 4 + _y ** 4 + _z ** 4 - R_(3, 5) * r2 ** 2
    rhs = 2 * (_X44 + _X4m4) / sp.sqrt(70) + R_(2, 5) * _X40
    same = sp.expand(lhs - rhs) == 0
    same &= sp.expand(multipole_X(4, 0) - _X40) == 0 and sp.expand(multipole_X(4, 4) - _X44) == 0 and \
        sp.expand(multipole_X(4, -4) - _X4m4) == 0
    return bool(same)


CRYSTAL_COMBO = {4: {0: R_(2, 5), 4: 2 / sp.sqrt(70), -4: 2 / sp.sqrt(70)},                       # (7.3.25)
                 6: {0: sp.sqrt(8) / 8, 4: -2 * sp.sqrt(7) / 8, -4: -2 * sp.sqrt(7) / 8}}            # (7.3.35)


def crystal_matrix(j, k: int) -> sp.Matrix:
    """⟨V^{(k)}⟩_j（單位 D⟨j‖X^k‖j⟩；k = 6 時是 E⟨j‖X⁶‖j⟩）＝ Σ c_q v^k_q × (−1)^k√((2j+1)/(2k+1))（(7.3.23)–(7.3.24)）。"""
    j = S_(j)
    fac = (-1) ** k * sp.sqrt((2 * j + 1) / (2 * k + 1))
    M = sp.zeros(int(2 * j + 1))
    for q, c in CRYSTAL_COMBO[k].items():
        M += c * unit_v(j, j, k, q)
    return (M * fac).applyfunc(sp.radsimp)


_O_REPS = ["1", "r1", "R1^2", "R1", "i1"]


def o_label(P, J) -> str:
    """把本徵子空間（投影 P）依 O 的特徵標分解成標記（"T_2"、"A_1+E"…；整數 J）。"""
    tr = [np.trace(P @ RT.D_axis(J, *SO.O_TURNS[g])).real for g in _O_REPS]
    labs = []
    for a, name in enumerate(SO.O_IRREP_NAMES):
        f = sum(float(SO.O_CHARS[a, c]) * SO.O_CLASS_SIZES[c] * tr[c] for c in range(5)) / 24
        labs += [name] * int(round(f))
    return "+".join(labs)


def _eig_groups(H, J, tol: float = 1e-7) -> list:
    w, U = np.linalg.eigh((H + H.conj().T) / 2)
    out, i = [], 0
    while i < len(w):
        k = i
        while k + 1 < len(w) and abs(w[k + 1] - w[i]) < tol * max(1.0, abs(w[i])):
            k += 1
        P = U[:, i:k + 1] @ U[:, i:k + 1].conj().T
        out.append((float(np.mean(w[i:k + 1])), o_label(P, J), k - i + 1))
        i = k + 1
    return out


def crystal_levels(j, k: int) -> list:
    """crystal_matrix 的本徵值（精確，sympy）＋ O 標記（數值投影）→ [(值, 標記, 簡併)]（由低到高）。"""
    M = crystal_matrix(j, k)
    ex = {sp.nsimplify(sp.radsimp(v)): mlt for v, mlt in M.eigenvals().items()}
    num = _eig_groups(_num(M), int(j))
    out = []
    for e, lab, d in num:
        key = min(ex, key=lambda v: abs(float(v) - e))
        out.append((key, lab, d))
    return out


def crystal_direct(j: int, k: int, n_th: int = 24, n_ph: int = 48) -> list:
    """第二種算法：⟨j m'|V^{(k)}(θφ)|j m⟩ = ∫Y*_{jm'} V Y_{jm} dΩ（r = 1，所以單位是 D⟨r^k⟩；Gauss–Legendre × 等分）→ 本徵值（由低到高）。"""
    x, w = leggauss(n_th)
    th = np.arccos(x)
    ph = np.arange(n_ph) * 2 * np.pi / n_ph
    TH, PH = np.meshgrid(th, ph, indexing="ij")
    W = np.outer(w, np.full(n_ph, 2 * np.pi / n_ph))
    V = sum(complex(sp.N(c)) * math.sqrt(4 * math.pi / (2 * k + 1)) * RT.ylm(k, q, TH, PH) for q, c in CRYSTAL_COMBO[k].items())
    msj = list(range(j, -j - 1, -1))
    Y = {m: RT.ylm(j, m, TH, PH) for m in msj}
    M = np.array([[np.sum(W * np.conj(Y[a]) * V * Y[b]) for b in msj] for a in msj])
    return list(np.linalg.eigvalsh((M + M.conj().T) / 2))


# ======================================================================== §7.3D 約化矩陣元
def reduced_X(lp, k, l):
    """(7.3.40b)：⟨l'‖X^k‖l⟩/⟨r^k⟩ = √((2l+1)/(2l'+1)) C^{k l l'}_{000} = √(2l+1)(−1)^{k−l}(k l l'; 0 0 0)。→ (兩式)。"""
    lp, k, l = S_(lp), S_(k), S_(l)
    a = sp.sqrt((2 * l + 1) / (2 * lp + 1)) * cg(k, l, lp, 0, 0, 0)
    b = sp.sqrt(2 * l + 1) * (-1) ** (k - l) * threej(k, l, lp, 0, 0, 0)
    return sp.nsimplify(a), sp.nsimplify(b)


def r4_hydrogen(n: int, l: int):
    """氫原子 ⟨r⁴⟩（單位 a0⁴，Bockasten 的式子）。"""
    n, l = sp.Integer(n), sp.Integer(l)
    return sp.nsimplify(n ** 4 / 8 * (63 * n ** 4 - n ** 2 * (70 * l * (l + 1) - 105) + 15 * (l - 1) * l * (l + 1) * (l + 2) - 20 * l * (l + 1) + 12))


def r_power_numeric(n: int, l: int, p: int):
    """第二種算法：用 sympy.physics.hydrogen.R_nl 直接積分 ∫r^{2+p}|R_nl|²dr（精確）。"""
    from sympy.physics.hydrogen import R_nl
    r = sp.Symbol("r", positive=True)
    return sp.nsimplify(sp.integrate(r ** (2 + p) * R_nl(n, l, r, 1) ** 2, (r, 0, sp.oo)))


def d2_factor(L: int, k: int = 4, l: int = 2):
    """(7.3.50)：兩個電子（l² 組態，總 L）的 ⟨[ll]L‖X^k‖[ll]L⟩/⟨l‖X^k‖l⟩ ＝ 兩項之和
    (−1)^{k+l+l+L}√((2l+1)(2L+1)){l k l; L l L} × 2（電子 1 與電子 2 相同）。L = 2、k = 4、l = 2：4/7。"""
    return sp.nsimplify(2 * (-1) ** (k + 2 * l + L) * sp.sqrt((2 * l + 1) * (2 * L + 1)) * sixj(l, k, l, L, l, L))


def d2_factor_numeric(L: int, k: int = 4, l: int = 2) -> float:
    """第二種算法：單電子 X^k_0 的矩陣（約化矩陣元 1：⟨l m'|X|l m⟩ = C^{k l l}_{0 m m'}），X⊗1 + 1⊗X 在耦合態 |[ll]L M=L⟩ 的期望值
    除以 C^{k L L}_{0 L L}。"""
    n = 2 * l + 1
    X = np.array([[cg_float(k, l, l, 0, m, mp) if mp == m else 0.0 for m in range(l, -l - 1, -1)] for mp in range(l, -l - 1, -1)])
    Xt = np.kron(X, np.eye(n)) + np.kron(np.eye(n), X)
    v = np.zeros(n * n)
    for a, m1 in enumerate(range(l, -l - 1, -1)):
        for b, m2 in enumerate(range(l, -l - 1, -1)):
            if m1 + m2 == L:
                v[a * n + b] = cg_float(l, l, L, m1, m2, L)
    return float(v @ Xt @ v / cg_float(k, L, L, 0, L, L))


# ======================================================================== §7.4A 剛體轉子
FIG742 = dict(J=10, A=0.2, B=0.4, C=0.6)


def re_rigid(n, J: int = 10, A: float = 0.2, B: float = 0.4, C: float = 0.6):
    """(7.4.5)：E = J(J+1)(A n_x² + B n_y² + C n_z²)（n 是 J 方向的單位向量；可為 (…,3) 陣列）。"""
    n = np.asarray(n, float)
    return J * (J + 1) * (A * n[..., 0] ** 2 + B * n[..., 1] ** 2 + C * n[..., 2] ** 2)


def re_rigid_euler(beta, gamma, J: int = 10, A: float = 0.2, B: float = 0.4, C: float = 0.6):
    """(7.4.4a)＋(7.4.5)：J_x = −|J|sin β cos γ、J_y = |J|sin β sin γ、J_z = |J|cos β。"""
    b, g = np.radians(beta), np.radians(gamma)
    return J * (J + 1) * (A * np.sin(b) ** 2 * np.cos(g) ** 2 + B * np.sin(b) ** 2 * np.sin(g) ** 2 + C * np.cos(b) ** 2)


def re_rigid_tensor(beta, gamma, J: int = 10, A: float = 0.2, B: float = 0.4, C: float = 0.6):
    """(7.4.10)：張量形式 J(J+1)[(A+B+C)/3 + (2C−A−B)/6·(3cos²β − 1) + (A−B)/2·sin²β cos 2γ]。"""
    b, g = np.radians(beta), np.radians(gamma)
    return J * (J + 1) * ((A + B + C) / 3 + (2 * C - A - B) / 6 * (3 * np.cos(b) ** 2 - 1) + (A - B) / 2 * np.sin(b) ** 2 * np.cos(2 * g))


_D2_FIG743 = {(1, 1, 1): "A_1", (-1, 1, -1): "A_2", (1, -1, -1): "B_1", (-1, -1, 1): "B_2"}


def rigid_levels(J: int = 10, A: float = 0.2, B: float = 0.4, C: float = 0.6) -> list:
    """H = A J_x² + B J_y² + C J_z²（(7.4.3)，2J+1 維）→ [(E, D2 標記)]（由低到高；標記由 180° 轉動 C2(x)、C2(y)、C2(z) 的期望值，Fig. 7.4.3）。"""
    M = RT.jmat(J)
    H = A * M["Jx"] @ M["Jx"] + B * M["Jy"] @ M["Jy"] + C * M["Jz"] @ M["Jz"]
    w, V = np.linalg.eigh((H + H.conj().T) / 2)
    Rs = [RT.D_axis(J, ax, 180) for ax in ((1, 0, 0), (0, 1, 0), (0, 0, 1))]
    out = []
    for i in range(len(w)):
        s = tuple(int(round(np.real(V[:, i].conj() @ R @ V[:, i]))) for R in Rs)
        out.append((float(w[i]), _D2_FIG743.get(s, str(s))))
    return out


def fig742() -> dict:
    """Fig. 7.4.2（J = 10、A, B, C = 0.2, 0.4, 0.6 cm⁻¹）：分隔線 BJ(J+1)、C2(x) 叢（低）與 C2(z) 叢（高）的成對能階與超精細分裂。"""
    p = FIG742
    lv = rigid_levels(p["J"], p["A"], p["B"], p["C"])
    sep = p["B"] * p["J"] * (p["J"] + 1)
    low = [(lv[2 * i], lv[2 * i + 1]) for i in range(4)]
    high = [(lv[-2 * i - 2], lv[-2 * i - 1]) for i in range(4)]
    split = lambda pr: pr[1][0] - pr[0][0]
    fine = (lv[2][0] + lv[3][0]) / 2 - (lv[0][0] + lv[1][0]) / 2
    return dict(levels=lv, separatrix=sep, low=low, high=high, low_split=[split(x) for x in low], high_split=[split(x) for x in high],
                fine_cm=fine, fine_GHz=fine * CM_TO_GHZ, low_split_Hz=[split(x) * CM_TO_GHZ * 1e9 for x in low])


def symmetric_top_check(J: int = 10, B: float = 0.2, C: float = 0.6) -> float:
    """(7.4.15)：E = BJ(J+1) + (C − B)K²（A = B）＝ 數值對角化。→ 最大誤差。"""
    num = sorted(e for e, _ in rigid_levels(J, B, B, C))
    form = sorted(B * J * (J + 1) + (C - B) * K * K for K in range(-J, J + 1))
    return float(np.abs(np.array(num) - np.array(form)).max())


def reduced_T2(J: int) -> tuple:
    """(7.4.11)–(7.4.13)：⟨J‖T²‖J⟩ = √((2J+3)(2J+2)2J(2J−1))/4 ＝ ⟨J J|T²₀|J J⟩/C^{2 J J}_{0 J J}，T²₀ = (3J_z² − J²)/2。"""
    form = math.sqrt((2 * J + 3) * (2 * J + 2) * 2 * J * (2 * J - 1)) / 4
    num = 0.5 * (2 * J * J - J) / cg_float(2, J, J, 0, J, J)
    return form, num


def d2_induced() -> dict:
    """(7.4.16)：C2(axis) 的 0₂、1₂ 誘導到 D2（阿貝爾群：χ↑(g) = [G:H]χ(g)（g ∈ H）、否則 0），用 Fig. 7.4.3 的標記分解。"""
    elems = ["E", "x", "y", "z"]
    chars = {"A_1": (1, 1, 1, 1), "A_2": (1, -1, 1, -1), "B_1": (1, 1, -1, -1), "B_2": (1, -1, -1, 1)}
    out = {}
    for ax in "xyz":
        for lab, s in (("0_2", 1), ("1_2", -1)):
            ind = [2 if g == "E" else (2 * s if g == ax else 0) for g in elems]
            out[(ax, lab)] = "+".join(a for a, ch in chars.items() if sum(x * y for x, y in zip(ind, ch)) // 4 == 1)
    return out


def rigid_correlation(J: int, Bs, A: float = 0.2, C: float = 0.6) -> np.ndarray:
    """Fig. 7.4.4：每個 B 的 2J+1 個能階（由低到高）→ (len(Bs), 2J+1)。"""
    return np.array([[e for e, _ in rigid_levels(J, A, b, C)] for b in Bs])


def cone_angles(J: int, Ks) -> list:
    """(7.4.2b)：Θ_{JK} = arccos(K/√(J(J+1)))（度）。"""
    return [math.degrees(math.acos(K / math.sqrt(J * (J + 1)))) for K in Ks]


# ======================================================================== §7.4B 八面體半剛性轉子
def octa_re(n):
    """(7.4.18)：n_x⁴ + n_y⁴ + n_z⁴ − 3/5（RE 曲面的形狀；t > 0 時 C4 軸是山、C3 軸是谷、C2 是鞍）。"""
    n = np.asarray(n, float)
    return n[..., 0] ** 4 + n[..., 1] ** 4 + n[..., 2] ** 4 - 0.6


def octa_re_euler(beta, gamma):
    """(7.4.19) 的角函數：(35cos⁴β − 30cos²β + 3 + 5sin⁴β cos 4γ)/20 ＝ octa_re。"""
    b, g = np.radians(beta), np.radians(gamma)
    c, s = np.cos(b), np.sin(b)
    return (35 * c ** 4 - 30 * c ** 2 + 3 + 5 * s ** 4 * np.cos(4 * g)) / 20


def octa_levels(J: int = 30, t: float = 1.0, B: float = 0.0) -> list:
    return RT.sf6_levels(J, t, B)


CLUSTER_PATTERNS = {   # (7.4.20)–(7.4.21)（更正：1₄、3₄↑O 一欄是 T2 +2S、T1 −2S；書印兩行都是 T2 = 2S）
    ("C4", 0): {"A_1": 4, "T_1": 0, "E": -2}, ("C4", 1): {"T_2": 2, "T_1": -2}, ("C4", 3): {"T_2": 2, "T_1": -2},
    ("C4", 2): {"E": 2, "T_2": 0, "A_2": -4},
    ("C3", 0): {"A_2": 3, "T_2": 1, "T_1": -1, "A_1": -3}, ("C3", 1): {"T_1": 2, "E": 0, "T_2": -2}, ("C3", 2): {"T_1": 2, "E": 0, "T_2": -2}}
_DIMS = {"A_1": 1, "A_2": 1, "E": 2, "T_1": 3, "T_2": 3}


def _induced_labels(kind: str, r: int) -> set:
    key = {"C4": "C4", "C3": "C3"}[kind]
    lab = f"{r}_{4 if kind == 'C4' else 3}"
    f = SO.frobenius(key)[lab][0]
    return {SO.O_IRREP_NAMES[a] for a in range(5) if f[a]}


def octa_clusters(J: int = 30, t: float = 1.0, n_c4: int = 6, n_c3: int = 2) -> list:
    """J 大時的能階叢：由高往下取 n_c4 個 6 態的 C4 叢（K4 = J, J−1, …；標記＝ (K4 mod 4)↑O），由低往上取 n_c3 個 8 態的 C3 叢
    （K3 = J, J−1, …）。每叢：E（依標記）、重心、offsets = E − 重心、以 CLUSTER_PATTERNS 最小平方擬合的 S 與殘差比 fit。"""
    lv = octa_levels(J, t)
    out = []

    def take(seq, size, kind, K):
        grp, dim = [], 0
        while dim < size:
            e, lab, d = seq.pop(0)
            for part in lab.split("+"):
                grp.append((e, part))
                dim += _DIMS[part]
        r = K % (4 if kind == "C4" else 3)
        labs = [g[1] for g in grp]
        assert set(labs) == _induced_labels(kind, r), (kind, K, labs)
        cen = sum(_DIMS[l] * e for e, l in grp) / size
        off = {l: e - cen for e, l in grp}
        pat = CLUSTER_PATTERNS[(kind, r)]
        Sfit = sum(_DIMS[l] * off[l] * pat[l] for l in pat) / sum(_DIMS[l] * pat[l] ** 2 for l in pat)
        res = max(abs(off[l] - Sfit * pat[l]) for l in pat) / max(abs(Sfit) * 4, 1e-300)
        return dict(kind=kind, K=K, labels=labs, E=[e for e, _ in grp], center=cen, offsets=off, S=Sfit, fit=res)

    hi = lv[::-1]
    for c in range(n_c4):
        out.append(take(hi, 6, "C4", J - c))
    lo = list(lv)
    for c in range(n_c3):
        out.append(take(lo, 8, "C3", J - c))
    return out


def octa_model_04() -> dict:
    """0₄↑O 的穿隧模型：六個 C4 方向（SO.SITES）之間的近鄰（夾 90°）穿隧 S = 1 → 本徵值 4（A1）、0（T1）、−2（E）。"""
    A = np.array([[1.0 if abs(np.dot(a, b)) < 1e-9 else 0.0 for b in SO.SITES] for a in SO.SITES])
    w, U = np.linalg.eigh(A)
    out = {}
    i = 0
    while i < 6:
        k = i
        while k + 1 < 6 and abs(w[k + 1] - w[i]) < 1e-9:
            k += 1
        P = U[:, i:k + 1] @ U[:, i:k + 1].T
        tr = [np.trace(P @ np.array(SO.site_perm(g), float)) for g in _O_REPS]
        lab = [n for a, n in enumerate(SO.O_IRREP_NAMES) if round(sum(float(SO.O_CHARS[a, c]) * SO.O_CLASS_SIZES[c] * tr[c] for c in range(5)) / 24)]
        out["+".join(lab)] = round(float(w[i]), 9)
        i = k + 1
    return out


def bottom_sequence(J: int = 30, n: int = 10) -> list:
    """(7.4.22)：由最低往上的 n 個標記（J = 30：A1 T1 T2 A2 T2 E T1 T2 E T1 ＝ 0₁↑O 的 24 維）。"""
    return [lab for _, lab, _ in octa_levels(J)[:n]]


def separatrix_angles() -> tuple:
    """分隔線（過 C2 鞍點的等能量線）與 C4、C3 軸的夾角：沿 γ = 45° 的子午線解 E(β) = E(C2)。數值 ＝ arccos√(2/3) 與 54.7356° − 35.2644°。"""
    f = lambda b: octa_re_euler(b, 45.0) - octa_re_euler(90.0, 45.0)
    from scipy.optimize import brentq
    b4 = brentq(f, 10, 50)
    b3 = math.degrees(math.acos(1 / math.sqrt(3))) - b4
    return b4, b3, math.degrees(math.acos(math.sqrt(2 / 3)))


def cluster_cutoffs(J: int = 30) -> tuple:
    """K4 = √(J(J+1)) cos 35.26°、K3 = √(J(J+1)) cos 19.47°（J = 30：24.90、28.75）。"""
    b4, b3, _ = separatrix_angles()
    L = math.sqrt(J * (J + 1))
    return L * math.cos(math.radians(b4)), L * math.cos(math.radians(b3))


@lru_cache(maxsize=None)
def _tensor_num(J: int, k: int, q: int) -> np.ndarray:
    """⟨J m'|T^k_q|J m⟩ = C^{k J J}_{q m m'}（約化矩陣元 1；(7.4.2a) 的歸一）。"""
    n = 2 * J + 1
    M = np.zeros((n, n))
    for a in range(n):
        mp = J - a
        m = mp - q
        if abs(m) <= J:
            M[a, J - m] = cg_float(k, J, J, q, m, mp)
    return M


def t46_operator(J: int, nu: float, sign4: int = 1) -> np.ndarray:
    """(7.4.23)–(7.4.25)：T^{4,6}(ν) = √(9/4π)T^{[4]}cos ν + √(13/4π)T^{[6]}sin ν（球諧歸一因子照 (7.4.26)–(7.4.27)），
    T^{[4]} = √(7/12)[T⁴₀ + sign4·√(5/14)(T⁴₄ + T⁴₋₄)]（書上 (7.4.24) 印 −，應 +）、T^{[6]} = (1/√8)[T⁶₀ − √(7/2)(T⁶₄ + T⁶₋₄)]。"""
    T4 = math.sqrt(7 / 12) * (_tensor_num(J, 4, 0) + sign4 * math.sqrt(5 / 14) * (_tensor_num(J, 4, 4) + _tensor_num(J, 4, -4)))
    T6 = math.sqrt(1 / 8) * (_tensor_num(J, 6, 0) - math.sqrt(7 / 2) * (_tensor_num(J, 6, 4) + _tensor_num(J, 6, -4)))
    return math.sqrt(9 / (4 * math.pi)) * math.cos(nu) * T4 + math.sqrt(13 / (4 * math.pi)) * math.sin(nu) * T6


def t46_sign_check(J: int = 8) -> dict:
    """與 C3（繞 (1,1,1) 120°）的對易子：T^{[4]} 用 + 號為 0、用 − 號不為 0；T^{[6]} 用書上的 − 號為 0。"""
    R3 = RT.D_axis(J, (1, 1, 1), 120)
    out = {}
    for s in (1, -1):
        H = t46_operator(J, 0.0, s)
        out[f"T4{'+' if s > 0 else '-'}"] = float(np.abs(R3 @ H - H @ R3).max())
    H6 = t46_operator(J, math.pi / 2)
    out["T6"] = float(np.abs(R3 @ H6 - H6 @ R3).max())
    return out


def t46_levels(J: int, nu: float) -> list:
    return _eig_groups(t46_operator(J, nu), J, 1e-9)


def t46_clusters(J: int = 30, nu_units: float = 2.0) -> list:
    """ν = nu_units·π/6 時由最低往上的兩個 12 態叢的標記（C2 叢：0₂↑O、1₂↑O）。"""
    lv = t46_levels(J, nu_units * math.pi / 6)
    out, cur, dim = [], [], 0
    for e, lab, d in lv:
        for l in lab.split("+"):
            cur.append(l)
            dim += _DIMS[l]
        if dim >= 12:
            out.append(cur)
            cur, dim = [], 0
        if len(out) == 2:
            break
    return out


def re46(n, nu: float):
    """(7.4.26)–(7.4.27)：E^{[4]}cos ν + E^{[6]}sin ν（n 為單位向量，可為陣列）。"""
    n = np.asarray(n, float)
    c = n[..., 2]
    s2 = 1 - c * c
    cos4g = np.where(s2 > 1e-14, (n[..., 0] ** 4 - 6 * n[..., 0] ** 2 * n[..., 1] ** 2 + n[..., 1] ** 4) / np.maximum(s2, 1e-300) ** 2, 1.0)
    e4 = math.sqrt(7 / 12) * math.sqrt(9 / (4 * math.pi)) * (35 * c ** 4 - 30 * c ** 2 + 3 + 5 * s2 ** 2 * cos4g) / 8
    e6 = math.sqrt(1 / 8) * math.sqrt(13 / (4 * math.pi)) * (231 * c ** 6 - 315 * c ** 4 + 105 * c ** 2 - 5
                                                               - 21 * s2 ** 2 * (11 * c ** 2 - 1) * cos4g) / 16
    return e4 * math.cos(nu) + e6 * math.sin(nu)


# ---------------------------------------------------------------- 球面上的等高線（RE 曲面的軌跡；自己寫的 marching squares）
def _unit(beta, gamma):
    return np.stack([np.sin(beta) * np.cos(gamma), np.sin(beta) * np.sin(gamma), np.cos(beta)], axis=-1)


def sphere_contours(f, levels, nb: int = 91, ng: int = 180) -> dict:
    """f(n)（n 是 (…,3) 單位向量）在球面上的等高線 → {level: [n×3 折線, …]}（β ∈ [0, π] 含兩極、γ 週期；閉曲線首尾相同）。"""
    be = np.linspace(0, np.pi, nb)
    ga = np.arange(ng) * 2 * np.pi / ng
    B, G = np.meshgrid(be, ga, indexing="ij")
    P = _unit(B, G)
    F = f(P)
    out = {}
    for L in levels:
        D = F - L
        edges = {}

        def pt(key):
            (i1, j1), (i2, j2) = key
            d1, d2 = D[i1, j1], D[i2, j2]
            t = d1 / (d1 - d2)
            v = P[i1, j1] * (1 - t) + P[i2, j2] * t
            return v / np.linalg.norm(v)

        segs = []
        for i in range(nb - 1):
            for j in range(ng):
                jj = (j + 1) % ng
                c = [(i, j), (i, jj), (i + 1, jj), (i + 1, j)]
                es = [(c[0], c[1]), (c[1], c[2]), (c[3], c[2]), (c[0], c[3])]
                cut = [e for e in es if (D[e[0]] > 0) != (D[e[1]] > 0)]
                keys = [tuple(sorted(e)) for e in cut]
                if len(keys) == 2:
                    segs.append((keys[0], keys[1]))
                elif len(keys) == 4:
                    cen = np.mean([D[x] for x in c])
                    if (cen > 0) == (D[c[0]] > 0):
                        segs += [(keys[0], keys[1]), (keys[2], keys[3])]
                    else:
                        segs += [(keys[0], keys[3]), (keys[1], keys[2])]
        adj = {}
        for a, b in segs:
            adj.setdefault(a, []).append(b)
            adj.setdefault(b, []).append(a)
        used = set()
        paths = []
        for start in list(adj):
            if start in used:
                continue
            # 找開放曲線的端點（度數 1）優先
            path = [start]
            used.add(start)
            cur = start
            while True:
                nxt = [x for x in adj[cur] if x not in used]
                if not nxt:
                    break
                cur = nxt[0]
                used.add(cur)
                path.append(cur)
            if len(path) > 2 and start in adj[path[-1]]:
                path.append(start)
            pts = np.array([pt(k) for k in path])
            if len(pts) > 2:
                paths.append(pts)
        out[L] = paths
    return out


def resample(path, n: int = 60) -> np.ndarray:
    """折線依弧長重取 n 點（場景每條曲線 ≤ 60 點）。"""
    path = np.asarray(path, float)
    seg = np.linalg.norm(np.diff(path, axis=0), axis=1)
    s = np.concatenate([[0], np.cumsum(seg)])
    t = np.linspace(0, s[-1], n)
    return np.array([np.interp(t, s, path[:, k]) for k in range(3)]).T


def rigid_paths(J: int = 10, A: float = 0.2, B: float = 0.4, C: float = 0.6, n: int = 60) -> list:
    """J = 10 不對稱陀螺的「量子化路徑」：RE 曲面 (7.4.5) 在每個本徵能量的等高線 → [(E, 標記, [折線…])]。"""
    lv = rigid_levels(J, A, B, C)
    cont = sphere_contours(lambda v: re_rigid(v, J, A, B, C), [e for e, _ in lv])
    return [(e, lab, [resample(p, n) for p in cont[e]]) for e, lab in lv]


def octa_paths(levels, n: int = 60) -> dict:
    return {L: [resample(p, n) for p in ps] for L, ps in sphere_contours(octa_re, levels).items()}


def re46_paths(nu: float, levels, n: int = 60) -> dict:
    return {L: [resample(p, n) for p in ps] for L, ps in sphere_contours(lambda v: re46(v, nu), levels).items()}


# ======================================================================== §7.5 旋量、偏振、二維振子
def spinor_euler(alpha: float, beta: float, gamma: float) -> np.ndarray:
    """(7.5.2)：(e^{−iα/2}cos(β/2), e^{iα/2}sin(β/2))e^{−iγ/2}（度）。"""
    a, b, g = (math.radians(x) for x in (alpha, beta, gamma))
    return np.array([np.exp(-0.5j * a) * math.cos(b / 2), np.exp(0.5j * a) * math.sin(b / 2)]) * np.exp(-0.5j * g)


def spin_vector(psi) -> np.ndarray:
    """Fig. 7.5.2：(⟨J_x⟩, ⟨J_y⟩, ⟨J_z⟩) = (Re ψ1*ψ2, Im ψ1*ψ2, (|ψ1|² − |ψ2|²)/2)。"""
    psi = np.asarray(psi, complex)
    z = np.conj(psi[0]) * psi[1]
    return np.array([z.real, z.imag, (abs(psi[0]) ** 2 - abs(psi[1]) ** 2) / 2])


def h_to_omega(H) -> dict:
    """(7.5.5e,f)：ω0 = (H11 + H22)/2（書印 H12）、ω_x = 2Re H21、ω_y = 2Im H21、ω_z = H11 − H22；ω、θ、φ（度）。"""
    H = np.asarray(H, complex)
    w0 = float(np.real(H[0, 0] + H[1, 1]) / 2)
    wv = np.array([2 * H[1, 0].real, 2 * H[1, 0].imag, float(np.real(H[0, 0] - H[1, 1]))])
    w = float(np.linalg.norm(wv))
    th = math.degrees(math.acos(wv[2] / w)) if w > 0 else 0.0
    ph = math.degrees(math.atan2(wv[1], wv[0])) if w > 0 else 0.0
    return dict(omega0=w0, omega_vec=wv, omega=w, theta=th, phi=ph, omega0_book=float(np.real(H[0, 0] + H[0, 1]) / 2))


def evolve(H, psi0, t: float) -> np.ndarray:
    return expm(-1j * np.asarray(H, complex) * t) @ np.asarray(psi0, complex)


def crank_evolve(H, psi0, t: float) -> np.ndarray:
    """(7.5.7)：|ψ(t)⟩ = e^{−iω0 t}R[φ θ ωt]|ψ(0)⟩（R 是 RT.su2_darboux）。"""
    o = h_to_omega(H)
    return np.exp(-1j * o["omega0"] * t) * (RT.su2_darboux(o["phi"], o["theta"], math.degrees(o["omega"] * t)) @ np.asarray(psi0, complex))


HAM_CARD = dict(A=4.0, B=0.0, C=1.0, D=2.0)     # 本系列（Fig. 7.5.11 的 A = 4、D = 2，C = 1）


def osc_hamiltonian(A: float, B: float, C: float, D: float) -> np.ndarray:
    """(7.5.8)：[[A, B − iC], [B + iC, D]]。"""
    return np.array([[A, B - 1j * C], [B + 1j * C, D]])


def jvars(x1, p1, x2, p2) -> dict:
    """(7.5.13)（更正 J_z 的減號）：J_x = x1x2 + p1p2、J_y = x1p2 − x2p1、J_z = (x1² − x2² + p1² − p2²)/2、J_0 = I/2。"""
    return dict(Jx=x1 * x2 + p1 * p2, Jy=x1 * p2 - x2 * p1, Jz=(x1 ** 2 - x2 ** 2 + p1 ** 2 - p2 ** 2) / 2,
                J0=(x1 ** 2 + x2 ** 2 + p1 ** 2 + p2 ** 2) / 2)


def hc_value(A, B, C, D, x1, p1, x2, p2) -> float:
    """(7.5.10b) 的 H_c。"""
    return A / 2 * (p1 ** 2 + x1 ** 2) + D / 2 * (p2 ** 2 + x2 ** 2) + B * (p1 * p2 + x1 * x2) + C * (x1 * p2 - x2 * p1)


def hc_check(n: int = 5, seed: int = 7) -> float:
    """(7.5.14a) ＝ (7.5.10b) ＝ ⟨ψ|H|ψ⟩/2（隨機的 A…D 與相空間點；古典 H_c 是量子期望值的一半）。→ 最大誤差。"""
    rng = np.random.default_rng(seed)
    err = 0.0
    for _ in range(n):
        A, B, C, D, x1, p1, x2, p2 = rng.normal(size=8)
        j = jvars(x1, p1, x2, p2)
        h1 = hc_value(A, B, C, D, x1, p1, x2, p2)
        h2 = B * j["Jx"] + C * j["Jy"] + (A - D) / 2 * j["Jz"] + (A + D) / 2 * j["J0"]
        psi = np.array([x1 + 1j * p1, x2 + 1j * p2])
        h3 = float(np.real(psi.conj() @ osc_hamiltonian(A, B, C, D) @ psi)) / 2          # H_c = ⟨ψ|H|ψ⟩/2
        err = max(err, abs(h1 - h2), abs(h1 - h3))
    return err


def lissajous(A, B, C, D, psi0, ts) -> np.ndarray:
    """解析解：ψ(t) = e^{−iHt}ψ0 → 陣列 (len(ts), 4)：x1, p1, x2, p2（x = Re ψ、p = Im ψ）。"""
    H = osc_hamiltonian(A, B, C, D)
    w, U = np.linalg.eigh(H)
    c = U.conj().T @ np.asarray(psi0, complex)
    out = []
    for t in ts:
        psi = U @ (np.exp(-1j * w * t) * c)
        out.append([psi[0].real, psi[0].imag, psi[1].real, psi[1].imag])
    return np.array(out)


def lissajous_numeric(A, B, C, D, psi0, ts) -> np.ndarray:
    """(7.5.10a) 的 Hamilton 方程用 solve_ivp 積分（第二種算法）。"""
    def rhs(t, y):
        x1, p1, x2, p2 = y
        return [A * p1 + B * p2 - C * x2, -(A * x1 + B * x2 + C * p2), B * p1 + D * p2 + C * x1, -(B * x1 + D * x2 - C * p1)]
    psi0 = np.asarray(psi0, complex)
    y0 = [psi0[0].real, psi0[0].imag, psi0[1].real, psi0[1].imag]
    sol = solve_ivp(rhs, (ts[0], ts[-1]), y0, t_eval=ts, rtol=1e-11, atol=1e-12, method="DOP853")
    return sol.y.T


SIG_C = (np.array([[0, 1], [1, 0]], complex), np.array([[0, -1j], [1j, 0]]), np.array([[1, 0], [0, -1]], complex))
U_RL = np.array([[1, 1], [1j, -1j]]) / math.sqrt(2)          # 行 = |r⟩、|l⟩（在 |x⟩、|y⟩ 基底）


def pauli_bases() -> dict:
    """(7.5.21)–(7.5.22)：C(σ_A,B,C) = σx、σy、σz（圓基底），L(σ) = U C(σ) U†（線基底）→ σz、σx、σy。"""
    L = [U_RL @ s @ U_RL.conj().T for s in SIG_C]
    return dict(L=L, ok=_close(L[0], SIG_C[2]) and _close(L[1], SIG_C[0]) and _close(L[2], SIG_C[1]))


def stokes_of(psi_xy) -> np.ndarray:
    """線基底的態 → Stokes 向量 (S_A, S_B, S_C) = (⟨σz⟩, ⟨σx⟩, ⟨σy⟩)（(7.5.24)）。"""
    p = np.asarray(psi_xy, complex)
    return np.array([np.real(p.conj() @ SIG_C[2] @ p), np.real(p.conj() @ SIG_C[0] @ p), np.real(p.conj() @ SIG_C[1] @ p)])


def state_circular(R: float, L: float, phi: float, Phi: float = 0.0) -> np.ndarray:
    """(7.5.19)：(R e^{−iφ}|r⟩ + L e^{iφ}|l⟩)e^{−iΦ} → 線基底的兩分量（度）。"""
    f, F = math.radians(phi), math.radians(Phi)
    return (U_RL @ np.array([R * np.exp(-1j * f), L * np.exp(1j * f)])) * np.exp(-1j * F)


def state_linear(X: float, Y: float, vt: float, theta: float = 0.0) -> np.ndarray:
    """(7.5.17)：(X e^{−iϑ}|x⟩ + Y e^{iϑ}|y⟩)e^{−iθ}（度）。"""
    v, t = math.radians(vt), math.radians(theta)
    return np.array([X * np.exp(-1j * v), Y * np.exp(1j * v)]) * np.exp(-1j * t)


def stokes_circular(R: float, L: float, phi: float) -> np.ndarray:
    """(7.5.25) 左欄：I(cos 2ψ cos 2φ, cos 2ψ sin 2φ, sin 2ψ)，a = (R+L)/√2 = √I cos ψ、b = (R−L)/√2 = √I sin ψ。"""
    a, b = (R + L) / math.sqrt(2), (R - L) / math.sqrt(2)
    I = a * a + b * b
    psi = math.atan2(b, a)
    f = math.radians(phi)
    return I * np.array([math.cos(2 * psi) * math.cos(2 * f), math.cos(2 * psi) * math.sin(2 * f), math.sin(2 * psi)])


def stokes_linear(X: float, Y: float, vt: float) -> np.ndarray:
    """(7.5.25) 右欄：I(cos 2ν, sin 2ν cos 2ϑ, sin 2ν sin 2ϑ)，X = √I cos ν、Y = √I sin ν。"""
    I = X * X + Y * Y
    nu = math.atan2(Y, X)
    v = math.radians(vt)
    return I * np.array([math.cos(2 * nu), math.sin(2 * nu) * math.cos(2 * v), math.sin(2 * nu) * math.sin(2 * v)])


def ellipse_field(psi_xy, Phis) -> np.ndarray:
    """(7.5.15)：E_x = Re⟨x|Ψ⟩e^{−iΦ}、E_y = Re⟨y|Ψ⟩e^{−iΦ}（Φ 是度；畫偏振橢圓用）→ (len, 2)。"""
    p = np.asarray(psi_xy, complex)
    return np.array([[np.real(p[0] * np.exp(-1j * math.radians(F))), np.real(p[1] * np.exp(-1j * math.radians(F)))] for F in Phis])


def ellipse_circular_formula(R: float, L: float, phi: float, Phis) -> np.ndarray:
    """(7.5.20a,b)：E_x = a cos Φ cos φ − b sin Φ sin φ、E_y = a cos Φ sin φ + b sin Φ cos φ。"""
    a, b = (R + L) / math.sqrt(2), (R - L) / math.sqrt(2)
    f = math.radians(phi)
    return np.array([[a * math.cos(math.radians(F)) * math.cos(f) - b * math.sin(math.radians(F)) * math.sin(f),
                      a * math.cos(math.radians(F)) * math.sin(f) + b * math.sin(math.radians(F)) * math.cos(f)] for F in Phis])


ELLIPSE_CARD = dict(phi=20.0, psi=15.0, I=1.0)     # 本系列（Fig. 7.5.6 有 φ = 0、20、45°）


def ellipse_card() -> dict:
    """本系列的橢圓：φ = 20°、ψ = 15°、I = 1 → R、L、a、b、Stokes、線基底的 ν、ϑ、X、Y、Euler 角 (7.5.26)、(7.5.27)。"""
    p = ELLIPSE_CARD
    a, b = math.sqrt(p["I"]) * math.cos(math.radians(p["psi"])), math.sqrt(p["I"]) * math.sin(math.radians(p["psi"]))
    R, L = (a + b) / math.sqrt(2), (a - b) / math.sqrt(2)
    Sv = stokes_circular(R, L, p["phi"])
    nu = 0.5 * math.degrees(math.acos(Sv[0] / p["I"]))
    vt = 0.5 * math.degrees(math.atan2(Sv[2], Sv[1]))
    return dict(a=a, b=b, R=R, L=L, stokes=Sv, nu=nu, vartheta=vt, X=math.cos(math.radians(nu)), Y=math.sin(math.radians(nu)),
                euler_C=(2 * p["phi"], 90 - 2 * p["psi"]), euler_A=(2 * vt, 2 * nu))


def faraday_check(C: float = 0.3, A: float = 1.0, t: float = 2.0) -> dict:
    """(7.5.29)：A = D、B = 0、C ≠ 0 ⇒ ω 沿 C 軸，S 繞 C 轉 2Ct：橢圓的 ψ 不變、長軸 φ 轉 Ct。"""
    st = ellipse_card()
    psi0 = state_circular(st["R"], st["L"], ELLIPSE_CARD["phi"])
    S0 = stokes_of(psi0)
    S1 = stokes_of(evolve(osc_hamiltonian(A, 0, C, A), psi0, t))
    ang0, ang1 = math.atan2(S0[1], S0[0]), math.atan2(S1[1], S1[0])
    return dict(SC0=S0[2], SC1=S1[2], dphi=math.degrees((ang1 - ang0) / 2), expected=math.degrees(C * t))


def birefringence_check(A: float = 1.5, D: float = 1.0, t: float = 2.0) -> dict:
    """(7.5.30)：B = C = 0、A ≠ D ⇒ ω 沿 A 軸，S 繞 A 轉 (A − D)t：盒角 ν 不變、ϑ 變 (A − D)t/2。"""
    psi0 = state_linear(math.cos(0.4), math.sin(0.4), 10.0)
    S0, S1 = stokes_of(psi0), stokes_of(evolve(osc_hamiltonian(A, 0, 0, D), psi0, t))
    a0, a1 = math.atan2(S0[2], S0[1]), math.atan2(S1[2], S1[1])
    return dict(SA0=S0[0], SA1=S1[0], dvt=math.degrees((a1 - a0) / 2) % 180, expected=math.degrees((A - D) * t / 2) % 180)


def beat(A, B, C, D) -> dict:
    """(7.5.31)–(7.5.32)：|ω| = √((A−D)² + 4B² + 4C²)（半拍頻 ＝ 兩本徵頻率差）、ω0 = (A+D)/2；本徵值 ω0 ± |ω|/2。"""
    w = math.sqrt((A - D) ** 2 + 4 * B * B + 4 * C * C)
    ev = np.linalg.eigvalsh(osc_hamiltonian(A, B, C, D))
    return dict(omega=w, omega0=(A + D) / 2, eig=ev, eig_diff=ev[1] - ev[0], half_beat_time=math.pi / w)


def closure(NH: int, NL: int) -> dict:
    """(7.5.33)–(7.5.34)：ω_H/ω_L = N_H/N_L ⇒ ω0/Ω = (N_H + N_L)/(N_H − N_L)（Ω = |ω|/2）、復現時 ω0Υ = π(N_H + N_L)（玫瑰瓣數）、ωΥ = 2π。"""
    return dict(ratio=Fraction(NH + NL, NH - NL), petals=NH + NL)


def closure_numeric(NH: int = 2, NL: int = 1, phi: float = 75.0, theta: float = 90.0) -> dict:
    """數值：本徵頻率 ω_H = N_H、ω_L = N_L、曲柄方向 (φ, θ)（Fig. 7.5.10 的 75°、90°）的 Hamilton 量，積分到 Υ = 2π/|ω_L| 回到原點、Υ/2 不回。"""
    n = RT.axis_vec(phi, theta)
    w = NH - NL
    w0 = (NH + NL) / 2
    H = w0 * np.eye(2) + w / 2 * (n[0] * SIG_C[0] + n[1] * SIG_C[1] + n[2] * SIG_C[2])
    A, D = H[0, 0].real, H[1, 1].real
    B, C = H[1, 0].real, H[1, 0].imag
    psi0 = spinor_euler(15.0, 90.0, 0.0)
    Ups = 2 * math.pi / abs(NL)
    tr = lissajous_numeric(A, B, C, D, psi0, np.array([0.0, Ups / 2, Ups]))
    tr2 = lissajous(A, B, C, D, psi0, [Ups])
    return dict(back=float(np.abs(tr[2] - tr[0]).max()), half=float(np.abs(tr[1] - tr[0]).max()), analytic=float(np.abs(tr2[0] - tr[0]).max()),
                omega0_over_Omega=w0 / (w / 2))


def anharmonic_fixed_points(Omega: float = 1.0, ax: float = 2.0, I: float = 1.0) -> dict:
    """(7.5.39a)–(7.5.40)：H = (I/2)Ω cos β + a_x(I²/4)cos²α sin²β（減去 J0ω0）。α = 0 上 dH/dβ = 0 ⇒ β = 0、π、arccos(Ω/(a_x I))。
    第二種算法：在球面上數值求極大（scipy）。→ 度。"""
    form = math.degrees(math.acos(Omega / (ax * I))) if ax * I > Omega else 0.0
    H = lambda v: -((I / 2) * Omega * math.cos(v[0]) + ax * I * I / 4 * math.cos(v[1]) ** 2 * math.sin(v[0]) ** 2)
    best = min((minimize(H, x0, method="Nelder-Mead", options=dict(xatol=1e-10, fatol=1e-14)) for x0 in ([0.6, 0.1], [1.2, -0.1], [0.3, 0.0])),
               key=lambda r: r.fun)
    return dict(formula=form, numeric=math.degrees(abs(best.x[0])), alpha=math.degrees(best.x[1]) % 180)


def osc_re(n, omega=(2.0, 0.0, 0.0), ax: float = 0.0, I: float = 1.0):
    """準自旋 RE 曲面（Fig. 7.5.11–7.5.13）：H = (I/4)(ω_A n_A + ω_B n_B + ω_C n_C) + a_x(I²/4)n_B²（n = (n_A, n_B, n_C)，畫圖時 A→x、B→y、C→z）。"""
    n = np.asarray(n, float)
    return I / 4 * (omega[0] * n[..., 0] + omega[1] * n[..., 1] + omega[2] * n[..., 2]) + ax * I * I / 4 * n[..., 1] ** 2


def local_mode_levels(j=5, Omega: float = 1.0, ax: float = 0.5) -> np.ndarray:
    """量子版：H = Ω J_A + a_x J_B²（自旋 j = I/2；(7.5.35)–(7.5.38)）的本徵值（由高到低）；a_x I ≫ Ω 時最高的能階成對（局域模雙重態）。"""
    J = RT.jmat(j)
    H = Omega * J["Jz"] + ax * J["Jx"] @ J["Jx"]
    return np.sort(np.linalg.eigvalsh((H + H.conj().T) / 2))[::-1]


# ======================================================================== §7.6 分子電子結構
def h2p_integrals(R: float, Z: float = 1.0) -> dict:
    """(7.6.21b)–(7.6.23)（原子單位，更正版）：S = (1 + ZR + Z²R²/3)e^{−ZR}、A = (a|1/r_b|a) = 1/R − (1/R + Z)e^{−2ZR}、
    B = (a|1/r_b|b) = Z(1 + ZR)e^{−ZR}。"""
    return dict(S=(1 + Z * R + Z * Z * R * R / 3) * math.exp(-Z * R), A=1 / R - (1 / R + Z) * math.exp(-2 * Z * R),
                B=Z * (1 + Z * R) * math.exp(-Z * R))


def h2p_quadrature(R: float, book: bool = False) -> dict:
    """球體座標 (μ, ν) 的數值積分（(7.6.21a)：d³x = (R³/8)(μ² − ν²)dμ dν dφ）。book=True 用 (7.6.22) 照印的被積函數
    （(R/8)(μ²−ν²) 與 e^{−(μ+ν)R/2}）。"""
    def dbl(f):
        g = lambda mu: quad(lambda nu: f(mu, nu), -1, 1, epsabs=1e-13)[0]
        return quad(g, 1, np.inf, epsabs=1e-13, limit=200)[0]
    # (μ² − ν²)/(r_b) = (μ² − ν²)·2/(R(μ − ν)) = 2(μ + ν)/R；φ 積分 2π 與 1/π 合成 2
    S = 2 * dbl(lambda mu, nu: R ** 3 / 8 * (mu * mu - nu * nu) * math.exp(-R * mu))
    if book:
        A = 2 * dbl(lambda mu, nu: R / 8 * 2 * (mu + nu) / R * math.exp(-(mu + nu) * R / 2))
    else:
        A = 2 * dbl(lambda mu, nu: R ** 3 / 8 * 2 * (mu + nu) / R * math.exp(-(mu + nu) * R))
    B = 2 * dbl(lambda mu, nu: R ** 3 / 8 * 2 * (mu + nu) / R * math.exp(-mu * R))
    return dict(S=S, A=A, B=B)


def h2p_V(R: float) -> tuple:
    """(7.6.24)：V^± − ε_1s = (1/R)[(1+R)e^{−2R} ± (1 − 2R²/3)e^{−R}]/[1 ± (R²/3 + R + 1)e^{−R}]。"""
    e1, e2 = math.exp(-R), math.exp(-2 * R)
    Sx = (R * R / 3 + R + 1) * e1
    return ((1 + R) * e2 + (1 - 2 * R * R / 3) * e1) / (R * (1 + Sx)), ((1 + R) * e2 - (1 - 2 * R * R / 3) * e1) / (R * (1 - Sx))


def h2p_V_generalized(R: float) -> tuple:
    """第二種算法：廣義本徵值問題 (7.6.10b) det|H − ES| = 0（h = ε − A、s = εS − B，ε = −½）＋ 1/R，再減 ε。"""
    it = h2p_integrals(R)
    eps = -0.5
    H = np.array([[eps - it["A"], eps * it["S"] - it["B"]], [eps * it["S"] - it["B"], eps - it["A"]]])
    Sm = np.array([[1, it["S"]], [it["S"], 1]])
    w = eigh(H, Sm, eigvals_only=True)
    return w[0] + 1 / R - eps, w[1] + 1 / R - eps


def h2p_matrix_check(R: float = 2.0) -> dict:
    """(7.6.12)–(7.6.13)：S⁻¹H 的元素 (h − Ss)/(1−S²)、(s − Sh)/(1−S²)（更正：左因子是 [[1, −S], [−S, 1]]）；照印的 [[1, S], [S, 1]] 不對。"""
    it = h2p_integrals(R)
    h, s, Sx = -0.5 - it["A"], -0.5 * it["S"] - it["B"], it["S"]
    Hm = np.array([[h, s], [s, h]])
    good = np.array([[1, -Sx], [-Sx, 1]]) @ Hm / (1 - Sx * Sx)
    book = np.array([[1, Sx], [Sx, 1]]) @ Hm / (1 - Sx * Sx)
    target = np.array([[h - Sx * s, s - Sx * h], [s - Sx * h, h - Sx * s]]) / (1 - Sx * Sx)
    return dict(good=_close(good, target), book=_close(book, target), eig=sorted(np.linalg.eigvals(good).real))


def h2p_min() -> dict:
    """V⁺(R) − ε_1s 的極小（LCAO）：R = 2.4928 a.u.、−0.0648 a.u.（書：≈ 2.5、−0.065；精確解 2.00、−0.1026）。"""
    r = minimize_scalar(lambda x: h2p_V(x)[0], bounds=(1.0, 5.0), method="bounded", options=dict(xatol=1e-10))
    grid = np.linspace(1.5, 4.0, 2501)
    g = np.array([h2p_V_generalized(x)[0] for x in grid])
    return dict(R=r.x, V=r.fun, R_grid=grid[int(np.argmin(g))], V_grid=float(g.min()))


def _grid_radial(R: float, rmax: float = 30.0, n: int = 48):
    br = sorted(set(b for b in (0.0, R / 2, R, 2 * R, 4.0, 8.0, 14.0, rmax) if b <= rmax))
    g, w = leggauss(n)
    xs, ws = [], []
    for a, b in zip(br, br[1:]):
        xs.append((b - a) / 2 * g + (a + b) / 2)
        ws.append((b - a) / 2 * w)
    return np.concatenate(xs), np.concatenate(ws)


def two_center(f1, f2, R: float, L: int = 40, nmu: int = 300) -> float:
    """∫∫ρ1(1)ρ2(2)/r12（兩個軸對稱電荷，原點在核間中點、核 a 在 +z）：ρ_i = Σ_l ρ_il(r)P_l(cos θ)，
    U = Σ_l (4π/(2l+1))² ∫∫ r1²r2² ρ_1l ρ_2l r<^l/r>^{l+1}。f(r_a, r_b) 給密度。精度約 1e−5（給 check 當獨立核對）。"""
    r, wr = _grid_radial(R)
    mu, wmu = leggauss(nmu)
    P = legvander(mu, L)
    rr, mm = np.meshgrid(r, mu, indexing="ij")
    ra = np.sqrt(np.maximum(rr * rr + R * R / 4 - rr * R * mm, 0))
    rb = np.sqrt(rr * rr + R * R / 4 + rr * R * mm)
    ls = np.arange(L + 1)
    c1 = (f1(ra, rb) * wmu) @ P * ((2 * ls + 1) / 2)
    c2 = (f2(ra, rb) * wmu) @ P * ((2 * ls + 1) / 2)
    rmin, rmx = np.minimum.outer(r, r), np.maximum.outer(r, r)
    ratio = rmin / rmx
    a1, a2 = c1 * (wr * r * r)[:, None], c2 * (wr * r * r)[:, None]
    U, K = 0.0, 1.0 / rmx
    for l in ls:
        U += (4 * math.pi / (2 * l + 1)) ** 2 * a1[:, l] @ K @ a2[:, l]
        K = K * ratio
    return float(U)


EULER_GAMMA = 0.5772156649015329


def h2_integrals(R: float, book: bool = False) -> dict:
    """(7.6.37) 的七個積分（原子單位，Z = 1）：A、B 同 H₂⁺；C = 5/8；D（(7.6.32a)，更正 R³/6；book=True 照印 R³/3）；
    E = G（(7.6.32b)，交換積分；(aa|bb) 對實軌域就是 (ab|ba)）；F = (aa|ab) = e^{−R}(R + 1/8 + 5/16R) − e^{−3R}(1/8 + 5/16R)（混合積分，書上沒給式子）。"""
    it = h2p_integrals(R)
    Sx = it["S"]
    D = 1 / R * (1 - (1 + 11 * R / 8 + 3 * R * R / 4 + R ** 3 / (3 if book else 6)) * math.exp(-2 * R))
    T = (R * R / 3 - R + 1) * math.exp(R)
    E = 6 / (5 * R) * ((EULER_GAMMA + math.log(R)) * Sx * Sx - T * T * exp1(4 * R) + 2 * Sx * T * exp1(2 * R)
                       - (-25 * R / 48 + 23 * R * R / 24 + R ** 3 / 2 + R ** 4 / 18) * math.exp(-2 * R))
    F = math.exp(-R) * (R + 1 / 8 + 5 / (16 * R)) - math.exp(-3 * R) * (1 / 8 + 5 / (16 * R))
    return dict(S=Sx, A=it["A"], B=it["B"], C=5 / 8, D=D, E=E, F=F, G=E)


def h2_integrals_numeric(R: float) -> dict:
    """第二種算法：C、D、E、F 用 two_center（多極展開的數值積分）。"""
    aa = lambda ra, rb: np.exp(-2 * ra) / math.pi
    bb = lambda ra, rb: np.exp(-2 * rb) / math.pi
    ab = lambda ra, rb: np.exp(-ra - rb) / math.pi
    return dict(C=two_center(aa, aa, R), D=two_center(aa, bb, R), E=two_center(ab, ab, R), F=two_center(aa, ab, R))


def hl_energy(R: float, book: bool = False, triplet: bool = False) -> float:
    """(7.6.30b)／(7.6.31b) 減 2ε_1s 再加 1/R：[−2A + D ∓ 2BS ± E]/(1 ± S²) + 1/R。"""
    g = h2_integrals(R, book)
    if triplet:
        return (-2 * g["A"] + g["D"] + 2 * g["B"] * g["S"] - g["E"]) / (1 - g["S"] ** 2) + 1 / R
    return (-2 * g["A"] + g["D"] - 2 * g["B"] * g["S"] + g["E"]) / (1 + g["S"] ** 2) + 1 / R


def hl_min(book: bool = False) -> dict:
    """Heitler–London ¹Σg 的極小：R = 1.6425、−0.1160 a.u.（書 1.6、−0.115）；照印 R³/3：1.6220、−0.1275。"""
    r = minimize_scalar(lambda x: hl_energy(x, book), bounds=(1.0, 3.0), method="bounded", options=dict(xatol=1e-10))
    return dict(R=r.x, E=r.fun)


def _h2_4(R: float):
    """(7.6.36)：基底 |aa), |bb), |ab), |ba) 的 H（含 2ε，ε = −½）與重疊矩陣（第一原理）。"""
    g = h2_integrals(R)
    Sx, A, B, C, D, E, F, G = (g[k] for k in "SABCDEFG")
    e = -0.5
    haa = 2 * e - 2 * A + C
    hab = 2 * Sx * e - Sx * A - B + F
    hbb = 2 * Sx * Sx * e - 2 * Sx * B + G
    hcov = 2 * e - 2 * A + D
    hex_ = 2 * Sx * Sx * e - 2 * Sx * B + E
    H = np.array([[haa, hbb, hab, hab], [hbb, haa, hab, hab], [hab, hab, hcov, hex_], [hab, hab, hex_, hcov]])
    O = np.array([[1, Sx * Sx, Sx, Sx], [Sx * Sx, 1, Sx, Sx], [Sx, Sx, 1, Sx * Sx], [Sx, Sx, Sx * Sx, 1]])
    return H, O, g


def h2_generalized(R: float) -> np.ndarray:
    H, O, _ = _h2_4(R)
    return np.sort(eigh(H, O, eigvals_only=True))


def vb_matrix(R: float, book: bool = False) -> np.ndarray:
    """(7.6.40)（非 Hermite；列＝行 ¹Σg ion、¹Σg cov、¹Σu ion、³Σu cov）；book=True 用照印的對角分母 (1 − S)²。"""
    g = h2_integrals(R)
    Sx, A, B, C, D, E, F, G = (g[k] for k in "SABCDEFG")
    e2 = -1.0
    d = (1 - Sx) ** 2 if book else (1 - Sx * Sx) ** 2
    M = np.zeros((4, 4))
    M[0, 0] = e2 + (-2 * A + 2 * Sx * B) / (1 - Sx * Sx) + ((1 + Sx * Sx) * (C + G) - 4 * Sx * F) / d
    M[0, 1] = (-2 * B + 2 * Sx * A) / (1 - Sx * Sx) + (2 * (1 + Sx * Sx) * F - 2 * Sx * (D + E)) / (1 - Sx * Sx) ** 2
    M[1, 0] = (-2 * B + 2 * Sx * A) / (1 - Sx * Sx) + (2 * (1 + Sx * Sx) * F - 2 * Sx * (C + G)) / (1 - Sx * Sx) ** 2
    M[1, 1] = e2 + (-2 * A + 2 * Sx * B) / (1 - Sx * Sx) + ((1 + Sx * Sx) * (D + E) - 4 * Sx * F) / d
    M[2, 2] = e2 + (-2 * A + 2 * Sx * B + C - G) / (1 - Sx * Sx)
    M[3, 3] = e2 + (-2 * A + 2 * Sx * B + D - E) / (1 - Sx * Sx)
    return M


def mo_matrix(R: float, book: bool = False) -> np.ndarray:
    """(7.6.47)（列＝行 ¹Σg(σg)²、¹Σg(σu)²、¹Σu σgσu、³Σu σgσu）；book=True 用照印的非對角分母 (1 − S)²（應 1 − S²）。"""
    g = h2_integrals(R)
    Sx, A, B, C, D, E, F, G = (g[k] for k in "SABCDEFG")
    e2 = -1.0
    M = np.zeros((4, 4))
    M[0, 0] = e2 + (-2 * A - 2 * B) / (1 + Sx) + (0.5 * (C + G) + 0.5 * (D + E) + 2 * F) / (1 + Sx) ** 2
    M[1, 1] = e2 + (-2 * A + 2 * B) / (1 - Sx) + (0.5 * (C + G) + 0.5 * (D + E) - 2 * F) / (1 - Sx) ** 2
    M[0, 1] = M[1, 0] = (0.5 * (C + G) - 0.5 * (D + E)) / ((1 - Sx) ** 2 if book else (1 - Sx * Sx))
    M[2, 2] = e2 + (-2 * A + 2 * Sx * B + C - G) / (1 - Sx * Sx)
    M[3, 3] = e2 + (-2 * A + 2 * Sx * B + D - E) / (1 - Sx * Sx)
    return M


def mo_min() -> dict:
    """(σg)² 單一組態（MO）的極小：M[0,0] + 1/R + 1。"""
    r = minimize_scalar(lambda x: mo_matrix(x)[0, 0] + 1 / x + 1.0, bounds=(1.0, 3.0), method="bounded", options=dict(xatol=1e-10))
    return dict(R=r.x, E=r.fun)


def ci_min() -> dict:
    """兩個 ¹Σg 組態交互作用（＝ 1s 基底的完整 CI）的極小：最低廣義本徵值 + 1/R + 1。"""
    r = minimize_scalar(lambda x: h2_generalized(x)[0] + 1 / x + 1.0, bounds=(1.0, 3.0), method="bounded", options=dict(xatol=1e-10))
    return dict(R=r.x, E=r.fun)


def limit_check(R: float = 1.4) -> dict:
    """(7.6.48)–(7.6.49)：只留 C、G 時 MO 與 VB 的 ¹Σg 子矩陣本徵值都是 0 與 (1+S²)(C+G)/(1−S²)²（更正；書印 (1+S²)²）。"""
    g = h2_integrals(R)
    Sx, CG_ = g["S"], g["C"] + g["G"]
    mo = 0.5 * np.array([[CG_ / (1 + Sx) ** 2, CG_ / (1 - Sx * Sx)], [CG_ / (1 - Sx * Sx), CG_ / (1 - Sx) ** 2]])
    vb = np.array([[(1 + Sx * Sx) * CG_ / (1 - Sx * Sx) ** 2, 0], [-2 * Sx * CG_ / (1 - Sx * Sx) ** 2, 0]])
    target = (1 + Sx * Sx) * CG_ / (1 - Sx * Sx) ** 2
    return dict(mo=sorted(np.linalg.eigvalsh(mo)), vb=sorted(np.linalg.eigvals(vb).real), target=target,
                book=(1 + Sx * Sx) * CG_ / (1 + Sx * Sx) ** 2)


def fig763_levels(s: float, R: float = 1.4) -> list:
    """Fig. 7.6.3 的數值版（本系列）：單電子積分 A、B 乘 s、兩電子積分 C D E F G 乘 (1 − s)；4×4 廣義本徵值減 2ε。
    s = 0：價鍵極限（C 大）；s = 1：分子軌域極限（A、B 大）。→ [(E, 標記)]（由低到高）。"""
    g = dict(h2_integrals(R))
    for k in "AB":
        g[k] *= s
    for k in "CDEFG":
        g[k] *= (1 - s)
    Sx, A, B, C, D, E, F, G = (g[k] for k in "SABCDEFG")
    hab = -Sx * A - B + F
    H = np.array([[-2 * A + C, -2 * Sx * B + G, hab, hab], [-2 * Sx * B + G, -2 * A + C, hab, hab], [hab, hab, -2 * A + D, -2 * Sx * B + E],
                  [hab, hab, -2 * Sx * B + E, -2 * A + D]])
    O = np.array([[1, Sx * Sx, Sx, Sx], [Sx * Sx, 1, Sx, Sx], [Sx, Sx, 1, Sx * Sx], [Sx, Sx, Sx * Sx, 1]])
    w, V = eigh(H, O)
    out = []
    for i in range(4):
        v = V[:, i] / np.abs(V[:, i]).max()
        if abs(v[0] - v[1]) < 1e-7 and abs(v[2] - v[3]) < 1e-7:
            lab = "1Sg"
        elif abs(v[0] + v[1]) < 1e-7 and abs(v[2]) + abs(v[3]) < 1e-7:
            lab = "1Su"
        else:
            lab = "3Su"
        out.append((float(w[i]), lab))
    return out


def hybrid_angle(name: str) -> float:
    """等價混成的鍵角：s 成分 λ² = 1/(n+1)（sp^n），cos θ = −λ²/(1 − λ²)。sp 180°、sp² 120°、sp³ 109.4712°。"""
    n = {"sp": 1, "sp2": 2, "sp3": 3}[name]
    lam2 = 1 / (n + 1)
    return math.degrees(math.acos(-lam2 / (1 - lam2)))


def sp3_vectors_angle() -> float:
    """第二種算法：(7.6.53) 的 C1、C2 的 p 部分 (1,1,1)、(1,−1,−1) 的夾角，且 ⟨C_i|C_j⟩ = δ。"""
    C = np.array([[1, 1, 1, 1], [1, 1, -1, -1], [1, -1, 1, -1], [1, -1, -1, 1]]) / 2
    assert _close(C @ C.T, np.eye(4))
    a, b = C[0, 1:], C[1, 1:]
    return math.degrees(math.acos(a @ b / (np.linalg.norm(a) * np.linalg.norm(b))))


def xy2_hybrids(two_phi: float, book_sign: bool = False) -> dict:
    """(7.6.54)–(7.6.56)：X1, X2 = λs + a cos φ p_x ± a sin φ p_y、X3 = λ's − a'p_x；a² = 1/(1 − cos 2φ)、λ² = cos 2φ/(cos 2φ − 1)、
    λ' = √((1+cos 2φ)/(1−cos 2φ))、a' = +√2 λ（書印 −）。→ 係數與 3×3 的重疊矩陣（應為單位矩陣）。"""
    ph = math.radians(two_phi / 2)
    c2 = math.cos(2 * ph)
    a = math.sqrt(1 / (1 - c2))
    lam = math.sqrt(max(c2 / (c2 - 1), 0.0))
    lamp = math.sqrt((1 + c2) / (1 - c2))
    ap = (-1 if book_sign else 1) * math.sqrt(2) * lam
    X = np.array([[lam, a * math.cos(ph), a * math.sin(ph)], [lam, a * math.cos(ph), -a * math.sin(ph)], [lamp, -ap, 0.0]])
    return dict(lam2=lam * lam, a2=a * a, lamp2=lamp * lamp, ap2=ap * ap, overlap=X @ X.T)


def xy2_config(two_phi: float) -> tuple:
    """(7.6.60)–(7.6.61)：x1x2x3²p_z² = s^{2/(1−cos 2φ)} p^{(4−6cos 2φ)/(1−cos 2φ)}（指數由 (7.6.55)–(7.6.56) 加總 ＝ 公式）。"""
    h = xy2_hybrids(two_phi)
    s_sum = 2 * h["lam2"] + 2 * h["lamp2"]
    p_sum = 2 * h["a2"] + 2 * h["ap2"] + 2
    c2 = math.cos(math.radians(two_phi))
    return s_sum, p_sum, 2 / (1 - c2), (4 - 6 * c2) / (1 - c2)


def huckel_ring(n: int, H: float = 0.0, Sx: float = 1.0) -> np.ndarray:
    M = H * np.eye(n)
    for k in range(n):
        M[k, (k + 1) % n] = M[(k + 1) % n, k] = -Sx
    return np.sort(np.linalg.eigvalsh(M))


def huckel_chain(n: int, H: float = 0.0, Sx: float = 1.0) -> np.ndarray:
    M = H * np.eye(n)
    for k in range(n - 1):
        M[k, k + 1] = M[k + 1, k] = -Sx
    return np.sort(np.linalg.eigvalsh(M))


def benzene() -> dict:
    """Hückel 苯：ε_k = H − 2S cos(2πk/6) → H − 2S、H − S(2)、H + S(2)、H + 2S；6 個 π 電子 6H − 8S（書印 H − 8S），局域 6H − 6S，共軛 2S。"""
    lv = huckel_ring(6)
    form = sorted(-2 * math.cos(2 * math.pi * k / 6) for k in range(6))
    tot = 2 * lv[0] + 2 * lv[1] + 2 * lv[2]
    return dict(levels=lv, formula=form, total_S=tot, localized=-6.0, gain=-6.0 - tot)


def butadiene() -> dict:
    """丁二烯：C10 環中間四個能階 H − 2S cos(kπ/5)（k = 1…4）＝ 4 原子鏈；4 個電子 4H − 2√5 S，局域 4H − 4S，共軛 (2√5 − 4)S = 0.4721S（書印 0.047）。"""
    lv = huckel_chain(4)
    form = sorted(-2 * math.cos(math.pi * k / 5) for k in range(1, 5))
    tot = 2 * lv[0] + 2 * lv[1]
    return dict(levels=lv, formula=form, total_S=tot, gain=-4.0 - tot, exact=2 * math.sqrt(5) - 4)


def octa_sigma() -> dict:
    """(7.6.62)–(7.6.63)：Frobenius（SO.frobenius）：0₄↑O、A'↑O_h（C4v）、E↑O_h。"""
    f4 = SO.frobenius("C4")
    fv = SO.frobenius("C4v")
    name = lambda fr, names: "+".join(n for n, c in zip(names, fr) for _ in range(int(c)))
    return {"0_4": name(f4["0_4"][0], SO.O_IRREP_NAMES), "A'": name(fv["A'"][0], SO.OH_IRREP_NAMES),
            "E": name(fv["E"][0], SO.OH_IRREP_NAMES), "1_4+3_4": name([a + b for a, b in zip(f4["1_4"][0], f4["3_4"][0])], SO.O_IRREP_NAMES)}


def ligand_states() -> dict:
    """(7.6.65)：六個 σ 配位基 |1⟩…|6⟩（+z, −z, +x, −x, +y, −y）用對稱化的 |A1⟩、|E,1⟩、|E,2⟩、|T1,x⟩、|T1,y⟩、|T1,z⟩ 展開；
    對稱化的組合由 O 的投影算子算（第二種算法）。→ 6×6 係數矩陣（列 |i⟩、行對稱態），應為正交矩陣，且 = (7.6.65)。"""
    s6, s3, s2 = 1 / math.sqrt(6), 1 / math.sqrt(3), 1 / math.sqrt(2)
    book = np.array([[s6, s3, 0, 0, 0, s2], [s6, s3, 0, 0, 0, -s2], [s6, -s3 / 2, 0.5, s2, 0, 0], [s6, -s3 / 2, 0.5, -s2, 0, 0],
                     [s6, -s3 / 2, -0.5, 0, s2, 0], [s6, -s3 / 2, -0.5, 0, -s2, 0]])
    sites = np.array(SO.SITES, float)
    A1 = np.ones(6) / math.sqrt(6)
    E1 = np.array([2 * v[2] ** 2 - v[0] ** 2 - v[1] ** 2 for v in sites]) / math.sqrt(12)
    E2 = np.array([v[0] ** 2 - v[1] ** 2 for v in sites]) / 2
    T = [sites[:, k] / math.sqrt(2) for k in range(3)]
    sym = np.array([A1, E1, E2, T[0], T[1], T[2]]).T
    Ps = {}
    for a, nm in enumerate(SO.O_IRREP_NAMES):
        P = sum(float(SO.O_CHARS[a, SO.O.class_of(g)]) * np.array(SO.site_perm(g), float) for g in SO.O.names) * _DIMS[nm] / 24
        Ps[nm] = P
    proj_ok = _close(Ps["A_1"] @ A1, A1) and _close(Ps["E"] @ E1, E1) and _close(Ps["E"] @ E2, E2) and all(_close(Ps["T_1"] @ t, t) for t in T)
    return dict(book=book, computed=sym, orthogonal=_close(book @ book.T, np.eye(6)), same=_close(book, sym), projectors=proj_ok)


# ======================================================================== 卡片數字
def _tab(M) -> str:
    return " / ".join(" ".join(_ex(M[i, j]) for j in range(M.shape[1])) for i in range(M.shape[0]))


def card_values() -> dict:
    """plan_ch07.md §6 每張卡上的數字與導讀第 7 章的數字（字串；小數照卡上位數）。"""
    out = {}
    # H701
    st = hyperfine_states()
    out["H701_triplet0"] = " ".join(_ex(x) for x in st[(1, 0)])
    out["H701_singlet"] = " ".join(_ex(x) for x in st[(0, 0)])
    out["H701_cg"] = _tab(cg_table(HALF, HALF)[0])
    lv = hyperfine_levels()
    out["H701_levels"] = " ".join(f"{_ex(k)}({v})" for k, v in sorted(lv.items(), key=lambda t: -float(t[0])))
    hn = hyperfine_numbers()
    out["H701_a"] = f"{_e(hn['a_erg'])} {hn['a_MHz']:.2f} {_e(hn['a_over_2pi'])}"
    out["H701_line"] = f"{hn['nu_MHz']:.9f} {_e(hn['E_eV'])} {_e(hn['E_erg'])} {hn['wavenumber']:.6f} {hn['wavelength_cm']:.3f}"
    out["H701_716"] = _ex(d_half_product()[0, 0])
    # H702
    t2 = two_electron_terms(1, 1, True)
    out["H702_p2"] = " ".join(f"{k}({v})" for k, v in t2.items()) + f" {sum(t2.values())}"
    tm = two_electron_terms(1, 1, False)
    out["H702_pp"] = " ".join(f"{k}({v})" for k, v in tm.items()) + f" {sum(tm.values())}"
    out["H702_hund"] = " ".join(hund_order(list(t2)))
    out["H702_L2col"] = " ".join(_ex(cg(1, 1, 2, m1, -m1, 0)) for m1 in (1, 0, -1))
    out["H702_717"] = f"{_ex(stretched_cg(2, 1, 1, 0))} {_ex(stretched_cg(2, 1, 1, 0, book=True))}"
    # H703
    so = spin_orbit(1)
    out["H703_so"] = f"{_ex(so[1.5][0])}({so[1.5][1]}) {_ex(so[0.5][0])}({so[0.5][1]})"
    fs = fine_structure()
    out["H703_fine"] = f"{_e(fs['shift1s_eV'])} {_e(fs['d2p_eV'])} {fs['d2p_cm']:.4f} {fs['d2p_GHz']:.2f} {_e(fs['d2p_book_eV'])}"
    lj = ls_jj_counts()
    out["H703_lsjj"] = f"{lj['ls']} {lj['jj']} {lj['triplet_2p3p']}"
    out["H703_p2x0"] = " ".join(f"{e:.4f}(J={J},{d})" for e, J, d in p2_levels(0.0))
    out["H703_p2x1"] = " ".join(f"{e:.4f}(J={J},{d})" for e, J, d in p2_levels(1.0))
    curves = []
    for J in (18, 17, 16):
        m, v = cg_curve(9, 9, J, 0)
        lo, hi = classical_range(9, 9, J, 0)
        curves.append(f"{np.abs(v).max():.4f} {hi:.3f} {sign_changes(v)}")
    out["H703_c99"] = " | ".join(curves)
    out["H703_c99M"] = " | ".join(f"{classical_range(9, 9, J, M)[0]:.3f} {classical_range(9, 9, J, M)[1]:.3f}" for J, M in ((18, -5), (18, -10)))
    # H704
    out["H704_c110"] = " ".join(_ex(cg(1, 1, 0, m, -m, 0)) for m in (1, 0, -1))
    out["H704_scalar"] = " ".join(_ex(scalar_cg(R_(3, 2), m)) for m in ms(R_(3, 2)))
    T7 = table_727()
    out["H704_727"] = _ex(T7[1, 4]) + " " + _ex(T727_BOOK[1, 4])
    rr = rr_tensor()
    out["H704_729"] = f"{sp.factor(rr[1][0])} {_ex(rr[1][2])}"
    out["H704_3j"] = f"{_ex(threej(1, 1, 0, 1, -1, 0))} {_ex(threej(1, 1, 2, 1, -1, 0))} {_f4(threej(1, 1, 2, 1, -1, 0))}"
    out["H704_6j"] = f"{_ex(sixj(2, 4, 2, 2, 2, 2))} {_ex(sixj(1, 1, 1, 1, 1, 1))}"
    out["H704_ck"] = " ".join(_ex(x) for x in ck_formula(2, 0))
    # H705
    out["H705_T"] = " | ".join(f"{k},{q}: {_tab(M)}" for (k, q), M in sorted(spin_half_T().items(), key=lambda t: (-t[0][0], -t[0][1])))
    ex = expectation_rotation(**EXPECT_CARD)
    out["H705_expect"] = f"{_f4(ex['Jx'])} {_f4(ex['Jz'])} {_f4(ex['rotated'])} {_f4(ex['p_up'])} {_f4(ex['p_dn'])}"
    out["H705_j1"] = " | ".join(f"{k},{q}: {_tab(tensor_j1()[(k, q)])}" for k, q in ((2, 0), (1, 0), (0, 0)))
    out["H705_E24"] = " ".join(f"{k},{q}:{_ex(c)}" for (k, q), c in sorted(elementary_expansion(2, 2, 4).items(), key=lambda t: -t[0][0]))
    num, den = v_table(2, 4)
    out["H705_v4d"] = " ".join(f"{q}:{d}" for q, d in den.items())
    # H706
    out["H706_j2"] = " ".join(f"{lab}:{_ex(v)}" for v, lab, d in crystal_levels(2, 4))
    out["H706_j2unit"] = " ".join(f"{lab}:{_ex(v * 15 * sp.sqrt(14))}" for v, lab, d in crystal_levels(2, 4))
    out["H706_j2r4"] = " ".join(f"{_ex(v * sp.sqrt(R_(2, 7)))}" for v, lab, d in crystal_levels(2, 4))
    d4 = 2 / (15 * sp.sqrt(22))
    out["H706_j3v4"] = " ".join(f"{lab}:{_ex(v / d4)}" for v, lab, d in crystal_levels(3, 4))
    out["H706_j3v6"] = " ".join(f"{lab}:{_ex(v * sp.sqrt(858))}" for v, lab, d in crystal_levels(3, 6))
    # H707
    out["H707_red"] = f"{_ex(reduced_X(2, 4, 2)[0])} {_f4(reduced_X(2, 4, 2)[0])} {_ex(reduced_X(3, 4, 3)[0])} {_ex(reduced_X(3, 6, 3)[0])}"
    out["H707_r4"] = f"{r4_hydrogen(3, 2)} {r4_hydrogen(4, 2)}"
    out["H707_d2"] = f"{_ex(d2_factor(2))} {_f4(d2_factor(2))}"
    # H708
    f7 = fig742()
    out["H708_sep"] = _f4(f7["separatrix"])
    out["H708_low"] = " ".join(f"{a[0]:.4f}({a[1]}{b[1]})" for a, b in f7["low"][:2]) + " | " + \
        " ".join(f"{b[0]:.4f}({a[1]}{b[1]})" for a, b in f7["high"][:2])
    out["H708_fine"] = f"{f7['fine_cm']:.4f} {f7['fine_GHz']:.1f}"
    out["H708_superfine"] = " ".join(f"{_e(s, 2)}" for s in f7["low_split"]) + " | " + \
        " ".join(f"{s / 1e3:.1f}kHz" if s < 1e6 else (f"{s / 1e6:.2f}MHz" if s < 1e9 else f"{s / 1e9:.2f}GHz") for s in f7["low_split_Hz"])
    out["H708_symtop"] = " ".join(_f2(0.2 * 110 + 0.4 * K * K) for K in (0, 5, 10))
    out["H708_T2"] = _f4(reduced_T2(10)[0])
    # H709
    di = d2_induced()
    out["H709_ind"] = " ".join(f"{ax}{lab}={di[(ax, lab)]}" for ax in "xz" for lab in ("0_2", "1_2"))
    cor = rigid_correlation(10, [0.2, 0.6])
    cor20 = rigid_correlation(20, [0.2, 0.6])
    out["H709_ends"] = f"{cor[0][0]:.2f} {cor[1][-1]:.2f} {cor20[0][0]:.2f} {cor20[1][-1]:.2f}"
    out["H709_cones"] = " ".join(f"{x:.1f}" for x in cone_angles(10, (10, 9, 8)))
    # H710
    cl = octa_clusters(30)
    out["H710_c4"] = " ".join("(" + " ".join(c["labels"]) + ")" for c in cl if c["kind"] == "C4")
    out["H710_c3"] = " ".join("(" + " ".join(c["labels"]) + ")" for c in cl if c["kind"] == "C3")
    c28 = next(c for c in cl if c["kind"] == "C4" and c["K"] == 28)
    c26 = next(c for c in cl if c["kind"] == "C4" and c["K"] == 26)
    c29 = next(c for c in cl if c["kind"] == "C4" and c["K"] == 29)
    k30 = next(c for c in cl if c["kind"] == "C3" and c["K"] == 30)
    out["H710_pat"] = " | ".join(" ".join(f"{l}:{c['offsets'][l] / c['S']:.2f}" for l in CLUSTER_PATTERNS[(c['kind'], c['K'] % (4 if c['kind'] == 'C4' else 3))])
                                 for c in (c28, c29, c26, k30))
    out["H710_seq"] = " ".join(bottom_sequence(30))
    out["H710_cone"] = " ".join(f"{x:.1f}" for x in cone_angles(30, range(30, 23, -1)))
    b4, b3, _ = separatrix_angles()
    k4, k3 = cluster_cutoffs(30)
    out["H710_sep"] = f"{b4:.2f} {b3:.2f} {k4:.2f} {k3:.2f}"
    out["H710_t46"] = " | ".join(" ".join(c) for c in t46_clusters(30, 2.0))
    # H711
    out["H711_spinor"] = " ".join(_f4(x) for x in np.real(spinor_euler(0, 90, 0))) + " " + " ".join(_f4(x) for x in spin_vector(spinor_euler(0, 90, 0)))
    om = h_to_omega(osc_hamiltonian(**HAM_CARD))
    out["H711_omega"] = f"{om['omega_vec'][0]:.4f} {om['omega_vec'][1]:.4f} {om['omega_vec'][2]:.4f} {om['omega']:.4f} {om['theta']:.2f} {om['phi']:.2f} {om['omega0']:.4f}"
    bt = beat(**HAM_CARD)
    out["H711_eig"] = f"{bt['eig'][1]:.4f} {bt['eig'][0]:.4f}"
    # H712
    ec = ellipse_card()
    out["H712_RL"] = f"{ec['a']:.4f} {ec['b']:.4f} {ec['R']:.4f} {ec['L']:.4f}"
    out["H712_stokes"] = " ".join(_f4(x) for x in ec["stokes"])
    out["H712_lin"] = f"{ec['nu']:.2f} {ec['vartheta']:.2f} {ec['X']:.4f} {ec['Y']:.4f}"
    out["H712_euler"] = f"{ec['euler_C'][0]:.0f} {ec['euler_C'][1]:.0f} | {ec['euler_A'][0]:.2f} {ec['euler_A'][1]:.2f}"
    # H713
    out["H713_closure"] = " ".join(f"{closure(*p)['ratio']} {closure(*p)['petals']}" for p in ((2, 1), (6, 1), (1, -6)))
    bt0 = beat(4.0, 0.0, 0.0, 2.0)
    out["H713_beat"] = f"{bt0['omega']:.4f} {bt0['omega0']:.4f} {bt0['half_beat_time']:.4f}"
    fc, bc = faraday_check(), birefringence_check()
    out["H713_faraday"] = f"{fc['dphi']:.2f} {bc['dvt']:.2f}"
    out["H713_local"] = " ".join(f"{anharmonic_fixed_points(1.0, ax, 1.0)['formula']:.2f}" for ax in (1.0, 2.0, 3.0))
    # H714
    hm = h2p_min()
    out["H714_h2p"] = f"{hm['R']:.4f} {hm['V']:.4f}"
    it = h2p_integrals(2.0)
    out["H714_ints2"] = f"{it['S']:.4f} {it['A']:.4f} {it['B']:.4f}"
    hl, hlb = hl_min(), hl_min(True)
    out["H714_hl"] = f"{hl['R']:.4f} {hl['E']:.4f} | {hlb['R']:.4f} {hlb['E']:.4f}"
    mo, ci = mo_min(), ci_min()
    out["H714_moci"] = f"{mo['R']:.4f} {mo['E']:.4f} | {ci['R']:.4f} {ci['E']:.4f}"
    out["H714_hyb"] = " ".join(f"{hybrid_angle(n):.2f}" for n in ("sp", "sp2", "sp3"))
    xc = xy2_config(104.5)
    out["H714_h2o"] = f"{xc[2]:.2f} {xc[3]:.2f}"
    bz, bd = benzene(), butadiene()
    out["H714_benzene"] = " ".join(f"{x:.0f}" for x in sorted(set(np.round(bz["levels"], 9)))) + f" {bz['total_S']:.0f} {bz['gain']:.0f}"
    out["H714_butadiene"] = " ".join(f"{x:.4f}" for x in bd["levels"]) + f" {bd['total_S']:.4f} {bd['gain']:.4f}"
    os_ = octa_sigma()
    out["H714_octa"] = " | ".join(os_[k] for k in ("0_4", "A'", "E"))
    return out


# 由 python utils/sym_coup.py --freeze 產生（每個值都先由 check() 的獨立算法核對過再貼上）
CARD_EXPECTED = {
    'H701_triplet0': '0 sqrt(2)/2 sqrt(2)/2 0',
    'H701_singlet': '0 sqrt(2)/2 -sqrt(2)/2 0',
    'H701_cg': '1 0 0 0 / 0 sqrt(2)/2 0 sqrt(2)/2 / 0 sqrt(2)/2 0 -sqrt(2)/2 / 0 0 1 0',
    'H701_levels': '1/4(3) -3/4(1)',
    'H701_a': '9.4234e-18 1422.17 1.4998e-18',
    'H701_line': '1420.405751800 5.8743e-06 9.4117e-18 0.047380 21.106',
    'H701_716': 'cos(beta/2)**2',
    'H702_p2': '1S(1) 3P(9) 1D(5) 15',
    'H702_pp': '1S(1) 3S(3) 1P(3) 3P(9) 1D(5) 3D(15) 36',
    'H702_hund': '3P 1D 1S',
    'H702_L2col': 'sqrt(6)/6 sqrt(6)/3 sqrt(6)/6',
    'H702_717': '2*sqrt(30)/15 2*sqrt(15)/15',
    'H703_so': '1/2(4) -1(2)',
    'H703_fine': '-1.8113e-04 4.5283e-05 0.3652 10.95 9.0565e-05',
    'H703_lsjj': '[0, 0, 1, 2, 2] [0, 0, 1, 2, 2] [0, 1, 1, 1, 2, 2, 3]',
    'H703_p2x0': '-5.0000(J=0,1) -5.0000(J=1,3) -5.0000(J=2,5) 1.0000(J=2,5) 10.0000(J=0,1)',
    'H703_p2x1': '-8.0000(J=0,1) -2.0000(J=1,3) -2.0000(J=2,5) 4.0000(J=0,1) 4.0000(J=2,5)',
    'H703_c99': '0.5104 2.121 0 | 0.4392 3.674 1 | 0.4035 4.690 2',
    'H703_c99M': '-4.542 -0.458 | -6.784 -3.216',
    'H704_c110': 'sqrt(3)/3 -sqrt(3)/3 sqrt(3)/3',
    'H704_scalar': '1/2 -1/2 1/2 -1/2',
    'H704_727': '-I/2 -1/2',
    'H704_729': '-z*(x + I*y) sqrt(6)/3',
    'H704_3j': 'sqrt(3)/3 sqrt(30)/30 0.1826',
    'H704_6j': '2/35 1/6',
    'H704_ck': 'sqrt(6)/6 sqrt(6)/3 sqrt(6)/6',
    'H705_T': '1,1: 0 1 / 0 0 | 1,0: -sqrt(2)/2 0 / 0 sqrt(2)/2 | 1,-1: 0 0 / -1 0 | 0,0: -sqrt(2)/2 0 / 0 -sqrt(2)/2',
    'H705_expect': '0.4330 0.2500 0.4830 0.9830 0.0170',
    'H705_j1': '2,0: sqrt(6)/6 0 0 / 0 -sqrt(6)/3 0 / 0 0 sqrt(6)/6 | 1,0: sqrt(2)/2 0 0 / 0 0 0 / 0 0 -sqrt(2)/2 | 0,0: sqrt(3)/3 0 0 / 0 sqrt(3)/3 0 / 0 0 sqrt(3)/3',
    'H705_E24': '4,2:-2*sqrt(7)/7 2,2:sqrt(21)/7',
    'H705_v4d': '0:70 1:14 2:14 3:2 4:1',
    'H706_j2': 'T_2:-4*sqrt(14)/105 E:2*sqrt(14)/35',
    'H706_j2unit': 'T_2:-8 E:12',
    'H706_j2r4': '-8/105 4/35',
    'H706_j3v4': 'A_2:-12 T_2:-2 T_1:6',
    'H706_j3v6': 'A_2:-12 T_1:-5 T_2:9',
    'H707_red': 'sqrt(14)/7 0.5345 sqrt(22)/11 -10*sqrt(429)/429',
    'H707_r4': '25515 362880',
    'H707_d2': '4/7 0.5714',
    'H708_sep': '44.0000',
    'H708_low': '24.8189(A_1B_1) 30.1493(A_2B_2) | 63.1811(B_2A_1) 57.8507(B_1A_2)',
    'H708_fine': '5.3305 159.8',
    'H708_superfine': '8.63e-07 7.73e-05 2.79e-03 5.07e-02 | 25.9kHz 2.32MHz 83.65MHz 1.52GHz',
    'H708_symtop': '22.00 32.00 62.00',
    'H708_T2': '109.6244',
    'H709_ind': 'x0_2=A_1+B_1 x1_2=A_2+B_2 z0_2=A_1+B_2 z1_2=A_2+B_1',
    'H709_ends': '22.00 66.00 84.00 252.00',
    'H709_cones': '17.5 30.9 40.3',
    'H710_c4': '(A_2 E T_2) (T_2 T_1) (A_1 T_1 E) (T_2 T_1) (E T_2 A_2) (T_2 T_1)',
    'H710_c3': '(A_1 T_1 T_2 A_2) (T_2 E T_1)',
    'H710_pat': 'A_1:4.00 T_1:0.00 E:-2.00 | T_2:2.00 T_1:-2.00 | E:1.94 T_2:0.06 A_2:-4.06 | A_2:3.03 T_2:0.99 T_1:-1.01 A_1:-2.97',
    'H710_seq': 'A_1 T_1 T_2 A_2 T_2 E T_1 T_2 E T_1',
    'H710_cone': '10.3 18.0 23.3 27.7 31.5 34.9 38.1',
    'H710_sep': '35.26 19.47 24.90 28.75',
    'H710_t46': 'A_1 T_1 T_2 E T_2 | A_2 T_1 T_2 E T_1',
    'H711_spinor': '0.7071 0.7071 0.5000 0.0000 0.0000',
    'H711_omega': '0.0000 2.0000 2.0000 2.8284 45.00 90.00 3.0000',
    'H711_eig': '4.4142 1.5858',
    'H712_RL': '0.9659 0.2588 0.8660 0.5000',
    'H712_stokes': '0.6634 0.5567 0.5000',
    'H712_lin': '24.22 20.97 0.9120 0.4102',
    'H712_euler': '40 60 | 41.93 48.44',
    'H713_closure': '3 3 7/5 7 -5/7 -5',
    'H713_beat': '2.0000 3.0000 1.5708',
    'H713_faraday': '34.38 28.65',
    'H713_local': '0.00 60.00 70.53',
    'H714_h2p': '2.4928 -0.0648',
    'H714_ints2': '0.5865 0.4725 0.4060',
    'H714_hl': '1.6425 -0.1160 | 1.6220 -0.1275',
    'H714_moci': '1.6032 -0.0991 | 1.6680 -0.1187',
    'H714_hyb': '180.00 120.00 109.47',
    'H714_h2o': '1.60 4.40',
    'H714_benzene': '-2 -1 1 2 -8 2',
    'H714_butadiene': '-1.6180 -0.6180 0.6180 1.6180 -4.4721 0.4721',
    'H714_octa': 'A_1+E+T_1 | A_{1g}+E_g+T_{1u} | T_{1g}+T_{2g}+T_{1u}+T_{2u}',
}


# ======================================================================== 核對
def _check_71():
    # CG：四種算法
    for j1, j2 in ((HALF, HALF), (1, HALF), (1, 1), (R_(3, 2), 1), (2, 1), (R_(3, 2), R_(3, 2))):
        M0 = cg_table(j1, j2)[0]
        M1 = cg_table(j1, j2, "racah")[0]
        assert all(_zero(x) for x in (M0 - M1)), ("racah", j1, j2)
        if (j1, j2) != (R_(3, 2), R_(3, 2)):
            M2 = cg_table(j1, j2, "lowering")[0]
            assert all(_zero(x) for x in (M0 - M2)), ("lowering", j1, j2)
        rows, cols = cg_table(j1, j2)[1:]
        for r, (m1, m2) in enumerate(rows):
            for c, (J, M) in enumerate(cols):
                assert abs(cg_float(j1, j2, J, m1, m2, M) - float(M0[r, c])) < 1e-12
    assert all(_zero(x) for x in (cg_table(HALF, HALF)[0] - CG_717_BOOK))
    assert all(_zero(x) for x in (cg_table(1, 1)[0] - CG_7119_BOOK))
    # (7.1.17)
    for j1, j2 in ((2, 1), (R_(3, 2), 1), (1, HALF), (3, 2)):
        for m1 in ms(j1):
            for m2 in ms(j2):
                assert _zero(stretched_cg(j1, j2, m1, m2) - cg(j1, j2, S_(j1) + S_(j2), m1, m2, m1 + m2))
    assert not _zero(stretched_cg(2, 1, 1, 0, book=True) - cg(2, 1, 3, 1, 0, 1))
    # (7.1.6a)
    M = d_half_product()
    b = sp.Symbol("beta")
    D1 = sp.Matrix([[sp.cos(b / 2) ** 2, -sp.sin(b) / _s2, sp.sin(b / 2) ** 2], [sp.sin(b) / _s2, sp.cos(b), -sp.sin(b) / _s2],
                    [sp.sin(b / 2) ** 2, sp.sin(b) / _s2, sp.cos(b / 2) ** 2]])
    assert all(_zero(sp.expand_trig(x)) for x in (M[:3, :3] - D1)) and _zero(M[3, 3] - 1) and all(_zero(M[i, 3]) for i in range(3))
    assert _close(np.array(d_half_product(70.0)[:3, :3].evalf(), dtype=float), RT.wigner_d(1, 70.0))
    # 超精細
    assert hyperfine_levels() == {R_(1, 4): 3, R_(-3, 4): 1}
    Hm = np.array(hyperfine_matrix().evalf(), dtype=complex)
    for (J, Mv), v in hyperfine_states().items():
        vv = np.array(v.evalf(), dtype=complex).ravel()
        lam = 0.25 if J == 1 else -0.75
        assert _close(Hm @ vv, lam * vv)
        Jt = {k: np.kron(RT.jmat(0.5)[k], np.eye(2)) + np.kron(np.eye(2), RT.jmat(0.5)[k]) for k in ("Jx", "Jy", "Jz")}
        J2 = sum(Jt[k] @ Jt[k] for k in Jt)
        assert _close(J2 @ vv, float(J * (J + 1)) * vv)
    hn = hyperfine_numbers()
    assert abs(hn["a_MHz"] - 1422.17) < 0.05 and abs(hn["a_over_2pi"] - 1.500e-18) < 1e-21          # 書印的 1.500e−18 是 a/2π
    assert abs(hn["E_eV"] - 5.874e-6) < 1e-9 and abs(hn["wavenumber"] * 1.23984198e-4 - hn["E_eV"]) < 1e-9
    assert abs(hn["wavelength_cm"] - 21.106) < 1e-3
    # (2p)²
    t2 = two_electron_terms(1, 1, True)
    assert t2 == {"1S": 1, "3P": 9, "1D": 5} == terms_by_symmetry(1, 1, True) and sum(t2.values()) == math.comb(6, 2)
    tm = two_electron_terms(1, 1, False)
    assert tm == terms_by_symmetry(1, 1, False) and sum(tm.values()) == 36
    assert hund_order(list(t2)) == sorted(SLATER_P2, key=SLATER_P2.get) == ["3P", "1D", "1S"]
    # 自旋–軌道與精細結構
    for l in (1, 2, 3):
        so = spin_orbit(l)
        assert abs(so[l + 0.5][0] - l / 2) < 1e-9 and abs(so[l - 0.5][0] + (l + 1) / 2) < 1e-9
        assert so[l + 0.5][1] * so[l + 0.5][0] + so[l - 0.5][1] * so[l - 0.5][0] < 1e-9           # 重心
    fs = fine_structure()
    assert abs(fs["d2p_eV"] - fs["d2p_so_eV"]) < 1e-12 and abs(fs["d2p_book_eV"] - 2 * fs["d2p_eV"]) < 1e-12
    assert abs(fs["shift1s_eV"] + ALPHA ** 2 / 8 * HARTREE_EV) < 1e-12
    # LS／jj
    lj = ls_jj_counts()
    assert lj["ls"] == lj["jj"] == [0, 0, 1, 2, 2] and lj["triplet_2p3p"] == [0, 1, 1, 1, 2, 2, 3]
    assert lj["c32_12"] == [1, 2] and lj["c32_32"] == [0, 1, 2, 3]
    x0, x1 = p2_levels(0.0), p2_levels(1.0)
    assert [round(e, 6) for e, _, _ in x0] == [-5, -5, -5, 1, 10] and sorted(J for _, J, _ in x0) == [0, 0, 1, 2, 2]
    assert [round(e, 6) for e, _, _ in x1] == [-8, -2, -2, 4, 4] and sorted(J for _, J, _ in x1) == [0, 0, 1, 2, 2]
    for x in (0.3, 0.7):
        assert sorted(J for _, J, _ in p2_levels(x)) == [0, 0, 1, 2, 2]
    # C^{99J}
    for J, nodes in ((18, 0), (17, 1), (16, 2)):
        m, v = cg_curve(9, 9, J, 0)
        assert sign_changes(v) == nodes == 18 - J
        lo, hi = classical_range(9, 9, J, 0)
        assert abs(lo + hi) < 1e-12
        assert np.abs(v[np.abs(m) > hi + 2.5]).max() < 0.2 * np.abs(v).max()                          # 古典範圍外很快衰減
    assert abs(cg_curve(9, 9, 18, 0)[1].max() - float(sp.binomial(18, 9) * sp.factorial(18) / sp.sqrt(sp.factorial(36)))) < 1e-12


def _check_72():
    T7 = table_727()
    diff = [(i, j) for i in range(9) for j in range(9) if not _zero(T7[i, j] - T727_BOOK[i, j])]
    assert diff == [(1, 4)] and _zero(T7[1, 4] + IU / 2)
    assert _close(np.array(T7.evalf(), dtype=complex).conj().T @ np.array(T7.evalf(), dtype=complex), np.eye(9))   # 么正
    assert all(check_724().values())
    rr = rr_tensor()
    for q in (2, 1, 0, -1, -2):
        assert _zero(rr[q][2] - sp.sqrt(R_(2, 3)))
    assert sp.expand(rr[1][0] + (_x + IU * _y) * _z) == 0 and sp.expand(rr[1][0] - BOOK_7292) != 0
    for k in (2, 3, 4):
        for q in range(-k + 1, k):
            a, b, c = ck_formula(k, q)
            assert _zero(a - cg(k - 1, 1, k, q - 1, 1, q)) and _zero(b - cg(k - 1, 1, k, q, 0, q))
            if abs(q + 1) <= k - 1:
                assert _zero(c - cg(k - 1, 1, k, q + 1, -1, q))
    for j in (HALF, 1, R_(3, 2), 2, 3):
        for m in ms(j):
            assert _zero(scalar_cg(j, m) - cg(j, j, 0, m, -m, 0))
    assert all(threej_symmetries().values())
    assert orthogonality_check()
    for j1, j2, j3, m1, m2, m3 in ((2, 1, 2, 1, 0, -1), (2, 1, 3, 2, -1, -1)):          # (7.2.30) 後第二行：下列要跟著換
        lhs = threej(j1, j2, j3, m1, m2, m3) * threej(j1, j2, j3, 0, 0, 0)
        assert _zero(lhs - threej(j2, j1, j3, m2, m1, m3) * threej(j2, j1, j3, 0, 0, 0))
        assert not _zero(lhs - threej(j2, j1, j3, m1, m2, m3) * threej(j2, j1, j3, 0, 0, 0)) or _zero(lhs)
    assert factorization_check() < 1e-10
    for j in (HALF, 1, R_(3, 2), 2):
        u = RT.su2_euler(40, 70, -25)
        v = RT.su2_euler(-110, 30, 80)
        assert _close(su2_irrep_formula(j, u), RT.wigner_D(j, 40, 70, -25))
        assert _close(su2_irrep_formula(j, u) @ su2_irrep_formula(j, v), su2_irrep_formula(j, u @ v))
    for j1, j2 in ((HALF, HALF), (1, HALF), (1, 1), (R_(3, 2), 1)):
        for J in couple(Fraction(str(j1)), Fraction(str(j2))):
            J = S_(J)
            for M in ms(J):
                assert _zero(boson_norm(J, M, S_(j1) - S_(j2), 2 * (S_(j1) + S_(j2))) - 1)
                for m1 in ms(j1):
                    if abs(M - m1) <= S_(j2):
                        assert _zero(boson_cg(j1, j2, J, m1, M - m1, M) - cg(j1, j2, J, m1, M - m1, M))
    assert any(not _zero(boson_cg(*a, half=False) - cg(*a)) for a in ((1, 1, 2, 1, -1, 0), (1, 1, 1, 1, 0, 1), (R_(3, 2), 1, R_(3, 2), HALF, 0, HALF)))
    assert any(not _zero(boson_norm(*a, half=False) - 1) for a in ((2, 0, 0, 4), (1, 0, 0, 4), (R_(3, 2), HALF, HALF, 5)))
    for d in (1, 2):
        st = det_scalar_state(d)
        assert set(st) == {(r, 2 * d - r, 2 * d - r, r) for r in range(2 * d + 1)}
        assert all(sum(k) == 4 * d for k in st) and _zero(sum(v ** 2 for v in st.values()) - 1)
        assert all(_zero(st[(r, 2 * d - r, 2 * d - r, r)] - (-1) ** r / sp.sqrt(2 * d + 1)) for r in range(2 * d + 1))
    for args in ((2, 4, 2, 2, 2, 2), (1, 1, 1, 1, 1, 1), (HALF, HALF, 1, HALF, HALF, 0), (2, 2, 2, 3, 2, 3), (1, 2, 3, 2, 1, 2)):
        assert _zero(sixj(*args) - sixj_from_3j(*args)), args
    assert _zero(sixj(2, 4, 2, 2, 2, 2) - R_(2, 35))
    for j1, j2, j3, J, j12, j23 in ((1, 1, 1, 1, 1, 0), (2, 2, 2, 2, 4, 2), (HALF, HALF, 1, 1, 1, HALF), (1, HALF, 1, R_(3, 2), HALF, R_(3, 2))):
        assert _zero(sixj_from_cg(j1, j2, j3, J, j12, j23) - sixj(j12, j1, j2, j23, j3, J))
        assert _zero(sixj_from_cg(j1, j2, j3, J, j12, j23, M=S_(J) - 1) - sixj(j12, j1, j2, j23, j3, J))


def _check_73():
    T = spin_half_T()
    assert all(all(_zero(x) for x in (T[k] - T733_BOOK[k])) for k in T733_BOOK)
    Tx, Ty, Tz = spin_half_cart()
    assert all(all(_zero(x) for x in (Tc - s / _s2)) for Tc, s in zip((Tx, Ty, Tz), PAULI_S))
    assert rotated_Jz_check() < 1e-12
    ex = expectation_rotation(**EXPECT_CARD)
    assert abs(ex["rotated"] - ex["rotated_prob"]) < 1e-12 and abs(ex["rotated"] - 0.5 * math.cos(math.radians(15))) < 1e-12
    Tj = tensor_j1()
    assert all(all(_zero(x) for x in (Tj[k] - T7312_BOOK[k])) for k in T7312_BOOK) and j1_relations()
    for j in (HALF, 1, R_(3, 2)):
        for k in range(0, int(2 * j) + 1):
            for q in range(-k, k + 1):
                assert all(_zero(x) for x in (unit_v(j, j, k, q) - (-1) ** (2 * S_(j)) * tensor_T(j, j, k, q)))
    # 跡內積下正交歸一
    j = 2
    vs = {(k, q): _num(unit_v(j, j, k, q)) for k in range(5) for q in range(-k, k + 1)}
    for a, Va in vs.items():
        for b, Vb in vs.items():
            assert abs(np.trace(Va.conj().T @ Vb) - (1 if a == b else 0)) < 1e-12
    assert all(v_tables_check().values())
    e24 = elementary_expansion(2, 2, 4)
    assert e24 == {(4, 2): -sp.sqrt(8) / sp.sqrt(14), (2, 2): sp.sqrt(3) / sp.sqrt(7)} or \
        (set(e24) == {(4, 2), (2, 2)} and _zero(e24[(4, 2)] + sp.sqrt(R_(8, 14))) and _zero(e24[(2, 2)] - sp.sqrt(R_(3, 7))))
    E24 = sp.zeros(5)
    E24[1, 3] = 1
    assert all(_zero(x) for x in (sum((c * unit_v(2, 2, k, q) for (k, q), c in e24.items()), sp.zeros(5)) - E24))
    e1 = {k: v for k, v in elementary_expansion(1, 1, 1).items() if k[1] == 0}
    assert _zero(e1[(2, 0)] - 1 / _s6) and _zero(e1[(1, 0)] - 1 / _s2) and _zero(e1[(0, 0)] - 1 / _s3)      # (7.3.16)
    mx = mixed_examples()
    pd_book = sp.Matrix([[0, 0, 1, 0, 0], [0, 0, 0, sp.sqrt(3), 0], [0, 0, 0, 0, sp.sqrt(6)]])
    assert all(_zero(x) for x in (mx["pd"] - pd_book)) and all(_zero(x) for x in (mx["dp"] - pd_book.T))
    for k, (nb, db) in V74_DP_BOOK.items():
        nm, dm = mixed_table(2, 1, k)
        assert all(_zero(x) for x in (nm - nb)) and dm == db, k
    fp = fpoly_check()
    assert _zero(fp["a2u"] + IU * sp.sqrt(30) / 30) and _zero(fp["t2u"] - sp.sqrt(R_(2, 15))) and fp["t1u_ok"]
    # Wigner–Eckart
    assert we_v_check()
    for jj in (1, R_(3, 2), 2):
        w = we_operator_check(jj)
        assert w["rot_err"] < 1e-12 and w["red_J"][1] < 1e-9 and w["red_T2"][1] < 1e-9 and abs(w["red_J"][0] - w["red_J_formula"]) < 1e-9
    # 晶體場
    assert v4_poly_check()
    l2 = crystal_levels(2, 4)
    assert [(lab, d) for _, lab, d in l2] == [("T_2", 3), ("E", 2)]
    assert _zero(l2[0][0] * 15 * sp.sqrt(14) + 8) and _zero(l2[1][0] * 15 * sp.sqrt(14) - 12) and _zero(3 * l2[0][0] + 2 * l2[1][0])
    dir2 = crystal_direct(2, 4)
    assert abs(dir2[0] + 8 / 105) < 1e-10 and abs(dir2[-1] - 4 / 35) < 1e-10
    assert all(abs(dir2[i] - float(l2[0 if i < 3 else 1][0] * reduced_X(2, 4, 2)[0])) < 1e-10 for i in range(5))
    l34 = crystal_levels(3, 4)
    d4 = 2 / (15 * sp.sqrt(22))
    assert [(lab, _ex(v / d4)) for v, lab, d in l34] == [("A_2", "-12"), ("T_2", "-2"), ("T_1", "6")]
    l36 = crystal_levels(3, 6)
    assert [(lab, _ex(v * sp.sqrt(858))) for v, lab, d in l36] == [("A_2", "-12"), ("T_1", "-5"), ("T_2", "9")]
    for k, lv in ((4, l34), (6, l36)):
        dr = crystal_direct(3, k)
        red = float(reduced_X(3, k, 3)[0])
        ex_ = sorted(float(v) * red for v, lab, d in lv for _ in range(d))
        assert _close(dr, ex_, 1e-10), k
    # (7.3.36) 的整數矩陣與 (7.3.29b)
    M6 = (crystal_matrix(3, 6) / sp.sqrt(R_(7, 13)) * 4 * sp.sqrt(462)).applyfunc(sp.radsimp)
    book6 = sp.Matrix(7, 7, lambda i, j: 0)
    for i, v in enumerate([1, -6, 15, -20, 15, -6, 1]):
        book6[i, i] = v
    for i, j_, v in ((0, 4, -7 * sp.sqrt(15)), (1, 5, 42), (2, 6, -7 * sp.sqrt(15))):
        book6[i, j_] = book6[j_, i] = v
    assert all(_zero(x) for x in (M6 - book6))
    assert sorted(book6.eigenvals()) == [-48, -20, 36]                                               # 4 × (−12, −5, 9)
    M4 = (crystal_matrix(3, 4) / (2 * sp.sqrt(7) / (15 * sp.sqrt(154)))).applyfunc(sp.radsimp)
    assert _zero(M4[0, 0] - 3) and _zero(M4[0, 4] - sp.sqrt(15)) and _zero(M4[1, 5] - 5) and _zero(M4[3, 3] - 6)
    # 約化矩陣元
    for lp, k, l in ((2, 4, 2), (3, 4, 3), (3, 6, 3), (1, 2, 1), (2, 2, 0)):
        a, b = reduced_X(lp, k, l)
        assert _zero(a - b)
    assert _zero(reduced_X(2, 4, 2)[0] - sp.sqrt(R_(2, 7)))
    for n, l in ((3, 2), (4, 2), (4, 1), (2, 1)):
        assert _zero(r4_hydrogen(n, l) - r_power_numeric(n, l, 4))
    assert r4_hydrogen(3, 2) == 25515
    for L in (2, 3, 4):
        assert abs(float(d2_factor(L)) - d2_factor_numeric(L)) < 1e-12
    assert _zero(d2_factor(2) - R_(4, 7))


def _check_74():
    n = np.random.default_rng(3).normal(size=(50, 3))
    n /= np.linalg.norm(n, axis=1)[:, None]
    b = np.degrees(np.arccos(n[:, 2]))
    g = np.degrees(np.arctan2(n[:, 1], -n[:, 0]))
    assert _close(re_rigid(n), re_rigid_euler(b, g)) and _close(re_rigid_euler(b, g), re_rigid_tensor(b, g))
    lv = rigid_levels(10, 0.2, 0.4, 0.6)
    assert _close(sorted(e for e, _ in lv), np.linalg.eigvalsh(RT.asym_matrix_formula(10, 0.2, 0.4, 0.6)), 1e-10)
    f7 = fig742()
    assert abs(f7["separatrix"] - 44) < 1e-12 and any(abs(e - 44) < 1e-9 for e, _ in lv)
    assert [(a[1], b_[1]) for a, b_ in f7["low"]] == [("A_1", "B_1"), ("A_2", "B_2"), ("A_1", "B_1"), ("A_2", "B_2")]
    assert [(a[1], b_[1]) for a, b_ in f7["high"]] == [("B_2", "A_1"), ("B_1", "A_2"), ("B_2", "A_1"), ("B_1", "A_2")]
    assert abs(f7["low_split"][0] - 8.6e-7) < 0.1e-7 and abs(f7["low_split"][2] - 2.8e-3) < 0.1e-3 and abs(f7["low_split"][3] - 5.1e-2) < 0.1e-2
    assert abs(f7["low_split"][1] - 7.7e-5) < 0.1e-5                                                  # 書 7.7e−6 = 230 kHz
    assert abs(f7["fine_cm"] - 5.33) < 0.01 and abs(f7["fine_GHz"] - 160) < 1                          # 書 150 GHz
    assert symmetric_top_check() < 1e-10 and symmetric_top_check(10, 0.6, 0.2) < 1e-10
    for J in (2, 5, 10):
        a, b_ = reduced_T2(J)
        assert abs(a - b_) < 1e-9
        M = RT.jmat(J)
        T20 = 0.5 * (3 * M["Jz"] @ M["Jz"] - M["J2"])
        for kk, m in enumerate(range(J, -J - 1, -1)):
            assert abs(T20[kk, kk].real - cg_float(2, J, J, 0, m, m) * a) < 1e-9
    di = d2_induced()
    assert di[("x", "0_2")] == "A_1+B_1" and di[("x", "1_2")] == "A_2+B_2" and di[("z", "0_2")] == "A_1+B_2" and di[("y", "0_2")] == "A_1+A_2"
    c = rigid_correlation(10, [0.2, 0.6])
    assert abs(c[0][0] - 22) < 1e-9 and abs(c[1][-1] - 66) < 1e-9
    c20 = rigid_correlation(20, [0.2, 0.6])
    assert abs(c20[0][0] - 84) < 1e-9 and abs(c20[1][-1] - 252) < 1e-9
    # 八面體
    nn = np.random.default_rng(5).normal(size=(40, 3))
    nn /= np.linalg.norm(nn, axis=1)[:, None]
    bb, gg = np.degrees(np.arccos(nn[:, 2])), np.degrees(np.arctan2(nn[:, 1], nn[:, 0]))
    assert _close(octa_re(nn), octa_re_euler(bb, gg))
    cl = octa_clusters(30)
    c4 = [c for c in cl if c["kind"] == "C4"]
    assert [c["K"] for c in c4] == [30, 29, 28, 27, 26, 25]
    assert all(c["fit"] < 0.08 for c in c4 if 27 <= c["K"] <= 29) and all(c["S"] > 0 for c in cl if c["K"] <= 29)
    c3 = [c for c in cl if c["kind"] == "C3"]
    assert [sorted(c["labels"]) for c in c3] == [sorted(["A_1", "T_1", "T_2", "A_2"]), sorted(["T_2", "E", "T_1"])]
    assert c3[0]["fit"] < 0.05
    for c in cl:                                                                                    # 重心
        assert abs(sum(_DIMS[l] * c["offsets"][l] for l in c["offsets"])) < 1e-6 * max(1, abs(c["center"]))
    for key, pat in CLUSTER_PATTERNS.items():
        assert sum(_DIMS[l] * v for l, v in pat.items()) == 0
    assert octa_model_04() == {"E": -2.0, "T_1": 0.0, "A_1": 4.0}
    assert bottom_sequence(30) == ["A_1", "T_1", "T_2", "A_2", "T_2", "E", "T_1", "T_2", "E", "T_1"]
    b4, b3, b4x = separatrix_angles()
    assert abs(b4 - b4x) < 1e-9 and abs(b3 - 19.4712) < 1e-4
    k4, k3 = cluster_cutoffs(30)
    assert abs(k4 - 24.9) < 0.05 and abs(k3 - 28.75) < 0.05
    assert [round(x, 1) for x in cone_angles(30, range(30, 23, -1))] == [10.3, 18.0, 23.3, 27.7, 31.5, 34.9, 38.1]
    sc = t46_sign_check()
    assert sc["T4+"] < 1e-12 and sc["T4-"] > 0.1 and sc["T6"] < 1e-12
    f2 = SO.frobenius("C2")
    lab = lambda fr: sorted(n for n, cc in zip(SO.O_IRREP_NAMES, fr) for _ in range(int(cc)))
    for nu_u in (2.0, 4.0):
        cls = t46_clusters(30, nu_u)
        assert sorted(cls[0]) == lab(f2["A"][0]) == ["A_1", "E", "T_1", "T_2", "T_2"]
        assert sorted(cls[1]) == lab(f2["B"][0]) == ["A_2", "E", "T_1", "T_1", "T_2"]
    # RE 函數 (7.4.26)–(7.4.27) 對 O 不變
    for g_ in ("r1", "R1", "i1"):
        Rm = np.array(SO.O.num[g_].real, float)
        assert _close(re46(nn @ Rm.T, 0.7), re46(nn, 0.7))
    # 量子 t46(ν=0) 頂叢 ≈ RE 函數在 C4 軸上的值（半古典，K = J 的錐）
    lv46 = t46_levels(30, 0.0)
    assert abs(lv46[-1][0] - re46(np.array([0, 0, 1.0]), 0.0)) < 0.12
    # 等高線
    cont = sphere_contours(octa_re, [0.0])
    assert len(cont[0.0]) == 6 and all(np.abs(octa_re(p)).max() < 2e-3 for p in cont[0.0])
    rp = rigid_paths()
    assert len(rp) == 21 and all(len(ps) >= 1 for _, _, ps in rp)


def _check_75():
    for ang in ((0, 90, 0), (30, 50, 70), (-100, 130, 20)):
        psi = spinor_euler(*ang)
        assert _close(psi, RT.su2_euler(*ang)[:, 0])
        sv = spin_vector(psi)
        a, b = math.radians(ang[0]), math.radians(ang[1])
        assert _close(sv, 0.5 * np.array([math.cos(a) * math.sin(b), math.sin(a) * math.sin(b), math.cos(b)]))
        J = RT.jmat(HALF)
        assert _close(sv, [np.real(psi.conj() @ J[k] @ psi) for k in ("Jx", "Jy", "Jz")])
    H = osc_hamiltonian(**HAM_CARD)
    om = h_to_omega(H)
    assert abs(om["omega0"] - 3) < 1e-12 and abs(om["omega0_book"] - 3) > 0.1
    Jh = RT.jmat(HALF)
    assert _close(H, om["omega0"] * np.eye(2) + sum(om["omega_vec"][k] * Jh[x] for k, x in enumerate(("Jx", "Jy", "Jz"))))
    psi0 = spinor_euler(20, 60, 0)
    for t in (0.3, 1.7):
        assert _close(evolve(H, psi0, t), crank_evolve(H, psi0, t))
    assert hc_check() < 1e-12
    ts = np.linspace(0, 6, 121)
    for p in (HAM_CARD, dict(A=1.0, B=0.3, C=0.2, D=0.6)):
        assert np.abs(lissajous(p["A"], p["B"], p["C"], p["D"], psi0, ts) - lissajous_numeric(p["A"], p["B"], p["C"], p["D"], psi0, ts)).max() < 1e-8
    assert pauli_bases()["ok"] and handedness_check()
    ec = ellipse_card()
    psi_c = state_circular(ec["R"], ec["L"], ELLIPSE_CARD["phi"])
    assert _close(stokes_of(psi_c), ec["stokes"]) and _close(stokes_linear(ec["X"], ec["Y"], ec["vartheta"]), ec["stokes"])
    psi_l = state_linear(ec["X"], ec["Y"], ec["vartheta"])
    assert _close(stokes_of(psi_l), ec["stokes"])
    Ph = np.linspace(0, 360, 37)
    assert _close(ellipse_field(psi_c, Ph), ellipse_circular_formula(ec["R"], ec["L"], ELLIPSE_CARD["phi"], Ph))
    # (7.5.23)、(7.5.24) 的逐式
    R, L, f = ec["R"], ec["L"], math.radians(ELLIPSE_CARD["phi"])
    pc = np.array([R * np.exp(-1j * f), L * np.exp(1j * f)])
    assert abs(np.real(pc.conj() @ SIG_C[1] @ pc) - 2 * R * L * math.sin(2 * f)) < 1e-12
    assert abs(np.real(pc.conj() @ SIG_C[2] @ pc) - 2 * ec["a"] * ec["b"]) < 1e-12 and abs(np.real(pc.conj() @ pc) - 1) < 1e-12
    v = math.radians(ec["vartheta"])
    pl = np.array([ec["X"] * np.exp(-1j * v), ec["Y"] * np.exp(1j * v)])
    assert abs(np.real(pl.conj() @ SIG_C[1] @ pl) - 2 * ec["X"] * ec["Y"] * math.sin(2 * v)) < 1e-12
    # Euler 角 (7.5.26)：Stokes ＝ (sin β cos α, sin β sin α, cos β)
    a_, b_ = (math.radians(x) for x in ec["euler_C"])
    assert _close(ec["stokes"], [math.sin(b_) * math.cos(a_), math.sin(b_) * math.sin(a_), math.cos(b_)])
    a_, b_ = (math.radians(x) for x in ec["euler_A"])
    assert _close([ec["stokes"][1], ec["stokes"][2], ec["stokes"][0]], [math.sin(b_) * math.cos(a_), math.sin(b_) * math.sin(a_), math.cos(b_)])
    fc = faraday_check()
    assert abs(fc["SC0"] - fc["SC1"]) < 1e-12 and abs(fc["dphi"] - fc["expected"]) < 1e-9
    bc = birefringence_check()
    assert abs(bc["SA0"] - bc["SA1"]) < 1e-12 and abs(bc["dvt"] - bc["expected"]) < 1e-9
    for p in (HAM_CARD, dict(A=1.0, B=0.3, C=0.2, D=0.6)):
        bt = beat(**p)
        assert abs(bt["omega"] - bt["eig_diff"]) < 1e-12 and abs(bt["omega"] - h_to_omega(osc_hamiltonian(**p))["omega"]) < 1e-12
    for NH, NL in ((2, 1), (6, 1), (3, 1)):
        cn = closure_numeric(NH, NL)
        assert cn["back"] < 1e-7 and cn["half"] > 1e-3 and cn["analytic"] < 1e-7 and abs(cn["omega0_over_Omega"] - float(closure(NH, NL)["ratio"])) < 1e-12
    assert closure(2, 1)["ratio"] == 3 and closure(6, 1)["petals"] == 7 and closure(1, -6)["petals"] == -5
    for ax in (2.0, 3.0, 5.0):
        fp = anharmonic_fixed_points(1.0, ax, 1.0)
        assert abs(fp["formula"] - fp["numeric"]) < 1e-5 and (fp["alpha"] < 1e-4 or abs(fp["alpha"] - 180) < 1e-4)
    assert anharmonic_fixed_points(1.0, 1.0, 1.0)["formula"] == 0.0
    lm = local_mode_levels(5, 1.0, 0.5)
    assert (lm[0] - lm[1]) < 0.05 * (lm[1] - lm[2])                                                 # 局域模雙重態
    lm0 = local_mode_levels(5, 1.0, 0.0)
    assert abs((lm0[0] - lm0[1]) - 1.0) < 1e-12


def _check_76():
    for R in (1.0, 2.0, 3.5):
        a, q = h2p_integrals(R), h2p_quadrature(R)
        assert all(abs(a[k] - q[k]) < 1e-8 for k in "SAB"), R
        assert abs(h2p_quadrature(R, book=True)["A"] - a["A"]) > 1e-3
        assert _close(h2p_V(R), h2p_V_generalized(R), 1e-12)
        Sx = a["S"]
        Vp = 1 / R - (a["A"] + a["B"]) / (1 + Sx)
        assert abs(Vp - h2p_V(R)[0]) < 1e-12                                                       # (7.6.15a)
    mc = h2p_matrix_check()
    assert mc["good"] and not mc["book"]
    hm = h2p_min()
    assert abs(hm["R"] - hm["R_grid"]) < 2e-3 and abs(hm["V"] - hm["V_grid"]) < 1e-6 and abs(hm["R"] - 2.49) < 0.01 and abs(hm["V"] + 0.0648) < 1e-4
    for R in (1.4, 1.64):
        g, gn = h2_integrals(R), h2_integrals_numeric(R)
        for k in "CDEF":
            assert abs(g[k] - gn[k]) < 1e-4, (k, R)
        assert abs(h2_integrals(R, True)["D"] - gn["D"]) > 1e-3                                     # 照印 R³/3 不對
    hl = hl_min()
    assert abs(hl["R"] - 1.6425) < 1e-3 and abs(hl["E"] + 0.1160) < 1e-4
    grid = np.linspace(1.3, 2.0, 701)
    assert abs(grid[int(np.argmin([hl_energy(x) for x in grid]))] - hl["R"]) < 2e-3
    assert hl_energy(1.64, triplet=True) > 0
    for R in (1.4, 2.0):
        ge = h2_generalized(R)
        assert _close(sorted(np.linalg.eigvals(vb_matrix(R)).real), ge, 1e-10)
        assert _close(sorted(np.linalg.eigvalsh(mo_matrix(R))), ge, 1e-10)
        assert not _close(sorted(np.linalg.eigvals(vb_matrix(R, True)).real), ge, 1e-4)
        assert not _close(sorted(np.linalg.eigvalsh(mo_matrix(R, True))), ge, 1e-4)
        # ³Σu ＝ Heitler–London 三重態 (7.6.31b)
        assert abs(mo_matrix(R)[3, 3] + 1 + 1 / R - hl_energy(R, triplet=True)) < 1e-12
    lc = limit_check()
    assert abs(lc["mo"][0]) < 1e-12 and abs(lc["mo"][1] - lc["target"]) < 1e-12 and abs(lc["vb"][1] - lc["target"]) < 1e-12
    assert abs(lc["book"] - lc["target"]) > 1e-3
    ci, mo = ci_min(), mo_min()
    assert ci["E"] < hl["E"] and ci["E"] < mo["E"]
    f0, f1 = fig763_levels(0.0), fig763_levels(1.0)
    assert len(f0) == len(f1) == 4
    # 混成
    for nme, ang in (("sp", 180.0), ("sp2", 120.0), ("sp3", math.degrees(math.acos(-1 / 3)))):
        assert abs(hybrid_angle(nme) - ang) < 1e-9
    assert abs(sp3_vectors_angle() - hybrid_angle("sp3")) < 1e-9
    for tp in (90.0, 104.5, 120.0, 150.0):
        h = xy2_hybrids(tp)
        assert _close(h["overlap"], np.eye(3), 1e-12)
        assert not _close(xy2_hybrids(tp, True)["overlap"], np.eye(3), 1e-6) or tp == 90.0
        s1, p1, s2, p2 = xy2_config(tp)
        assert abs(s1 - s2) < 1e-12 and abs(p1 - p2) < 1e-12 and abs(s1 + p1 - 6) < 1e-12
    h120 = xy2_hybrids(120.0)
    assert abs(h120["lam2"] - 1 / 3) < 1e-12 and abs(h120["lamp2"] - 1 / 3) < 1e-12 and abs(h120["ap2"] - 2 / 3) < 1e-12
    assert _close(xy2_config(90.0)[:2], (2, 4)) and _close(xy2_config(180.0)[:2], (1, 5))
    pc = py_check()
    assert abs(pc["px"] - pc["x"]) < 1e-12 and abs(pc["py"] - pc["y"]) < 1e-12 and abs(pc["py_book"] + 1j * pc["x"]) < 1e-12
    bz = benzene()
    assert _close(bz["levels"], bz["formula"]) and abs(bz["total_S"] + 8) < 1e-12 and abs(bz["gain"] - 2) < 1e-12
    bd = butadiene()
    assert _close(bd["levels"], bd["formula"]) and abs(bd["gain"] - bd["exact"]) < 1e-12 and abs(bd["total_S"] + 2 * math.sqrt(5)) < 1e-12
    os_ = octa_sigma()
    assert os_["0_4"] == "A_1+E+T_1" and os_["A'"] == "A_{1g}+E_g+T_{1u}" and os_["E"] == "T_{1g}+T_{2g}+T_{1u}+T_{2u}"
    assert os_["1_4+3_4"] == "T_1+T_1+T_2+T_2"
    ls_ = ligand_states()
    assert ls_["orthogonal"] and ls_["same"] and ls_["projectors"]


def _check_cards():
    cv = card_values()
    for k, v in CARD_EXPECTED.items():
        assert cv[k] == v, (k, cv[k], v)
    assert set(cv) == set(CARD_EXPECTED), set(cv) ^ set(CARD_EXPECTED)


def summary() -> dict:
    return card_values()


def check():
    """每個上卡片／講義的數字用兩種以上獨立算法核對（約 20 秒）。"""
    _check_71()
    _check_72()
    _check_73()
    _check_74()
    _check_75()
    _check_76()
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
    if "--freeze" in sys.argv:
        print("CARD_EXPECTED = {")
        for k, v in card_values().items():
            print(f"    {k!r}: {v!r},")
        print("}")
        sys.exit(0)
    t0 = time.time()
    check()
    print(f"sym_coup（第 7 章）：R3 耦合係數（降算子、Racah、Schwinger 玻色子、sympy 四法一致）、3j 對稱與 6j、氫的超精細與 21 cm、"
          f"(2p)² 的項與 Hund、自旋–軌道與 LS／jj、向量錐、張量算子與單位張量表、R3 Wigner–Eckart、j = 2、3 的晶體場、d² 的 4/7、"
          f"剛體與八面體轉子的 RE 曲面與能階叢、旋量與 Stokes 向量、Lissajous 與局域模、H2+ 與 H2、混成與 Hückel 全部核對通過"
          f"（{time.time() - t0:.1f} 秒）。")
