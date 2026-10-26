# 由 combkit/tools/sync_from_harter_sym.py 從 harter_sym 專案的 codes/utils/sym_rot.py 複製；請改 harter_sym 那份再同步。
r"""
sym_rot.py —— Harter《Principles of Symmetry, Dynamics, and Spectroscopy》（Wiley 1993）第 5 章
Representations of Continuous Rotation Groups and Applications 的數學（純 Python＋numpy＋sympy＋scipy，不 import manim；
有限群工具借用同資料夾的 sym_mat（譜分解）、sym_nonabel（Hamilton 轉、FiniteGroup、D6、四元數群）、sym_oct（O、O_h、T_d、Y 的特徵標），不重寫）。
書上的式號、筆誤與〔注〕見 ref/md/ch05_formulas.md（相對於 C:\mypython\manim\projects\harter_sym）。

約定
  * 角度一律用「度」（Y、ylm 的 θ, φ 用弧度）；ħ = 1。精確值用 sympy；卡片上的小數一律 f"{float(x):.4f}"（sympy 數、Fraction 要 float()）。
  * Euler 角 (αβγ)：R(αβγ) = R(α00)R(0β0)R(00γ) = R_z(α)R_y(β)R_z(γ)（實驗室的 z、y、z 曲柄；最右邊先做，(5.3.3)）。
    本體曲柄依序 z̄(α)、ȳ(β)、z̄(γ) 得到同一個 R（順序反過來，(5.3.5)）。0 ≤ β ≤ 180 時 euler_from_matrix 唯一（β = 0、180 取 γ = 0）。
  * Darboux [φθω]：軸 n̂ = (cos φ sin θ, sin φ sin θ, cos θ)，逆時針轉 ω（(5.5.5c)）。
  * SU(2)：𝒟^{1/2}(αβγ) = e^{−iαJz}e^{−iβJy}e^{−iγJz}（(5.4.30)）＝ SN.su2(四元數)，四元數 q = (cos ω/2, sin ω/2 n̂)，𝒟^{1/2} = q₀1 − i q⃗·σ。
  * 𝒟^j_{m'm}(αβγ) = ⟨j m'|R|j m⟩，基底 m = j, j−1, …, −j（列＝行同序）；= e^{−im'α} d^j_{m'm}(β) e^{−imγ}（(5.4.45)，標準 Wigner 相位）。
  * 書上的生成元 𝒢_x, 𝒢_y, 𝒢_z（3×3 實反對稱，[𝒢_x, 𝒢_y] = +𝒢_z），G = J/iħ。本體角動量用 Jahn 的反轉算子 J_x̄ = J_x*（矩陣共軛），[J_x̄, J_ȳ] = −iJ_z̄。
  * 球諧 Y^ℓ_m 含 Condon–Shortley 號（(5.5.71)、(5.6.14)）；√(4π)Y^ℓ_m(β, α) = 𝒟^{ℓ*}_{m0}(αβ·)√(2ℓ+1)（(5.5.65)；書上寫 Y^j_m）。
  * 射線代數（§5.7）：每個元素取一個 SU(2) 代表（D6 的 h^k ＝ 繞 z 轉 60k°（0 < 60k < 360）、ρ 的軸方位角 D6S_AXES；O 的元素 ＝ SO.O_TURNS 的轉角），乘積 ab = ±c。
  * 群乘積 ab = 先 b 後 a；不要把模組代號 RT、SO、SN、SM 當變數名。

§5.1–5.2 O3、O2、D∞h（H501）
  INV3、SIGMA_XY、I_Z（= σ_xy I = 繞 z 的 180°，det +1）、is_orthogonal、n_params(n)、N_PARAMS = [0, 1, 3, 6, 10, 15]（O1…O6）
  o2_rot(φ)（(5.1.10a)）、o2_refl(φ)（(5.1.10b)，＝ 對 φ/2 的直線反射；φ = 60°：[[0.5000, 0.8660], [0.8660, −0.5000]]）、o2_standing(m)（(5.2.3)，sympy）
  O2_LABELS（A1 Σ⁺、A2 Σ⁻、E1 Π、E2 Δ、E3 Φ）、O2_MOVING_SIGMA_MINUS_BOOK（照印 −1，應 +1）
  DINFH_CLASSES、DINFH_IRREPS、dinfh_chars(γ)（Fig. 5.2.2；Σ^± 的 ± 是鏡面的號，橫向 π 轉的號 ＝ ±×宇稱）＝ dinfh_trace_check（向量與二次式的跡）
§5.3 Euler、Darboux、Hamilton（H502–H503、H507）
  rz、ry、rx、euler_matrix(α, β, γ)（(5.5.76)）、euler_matrix_sym()、euler_from_matrix(R)、euler_alt(α, β, γ, ±1)（(5.3.1)）
    euler_matrix(50, 60, 70) = [[−0.6099, −0.5640, 0.5567], [0.7350, −0.1401, 0.6634], [−0.2962, 0.8138, 0.5000]]（＝ (−130, −60, −110) 的矩陣；
    SU(2) 差一個號：su2_euler(−130, −60, −110) = −su2_euler(50, 60, 70)）
  axis_vec(φ, θ)、darboux_matrix(φ, θ, ω)（四元數）＝ darboux_matrix_cranks（五個曲柄 R(φ00)R(0θ0)R(ω00)R(0−θ0)R(−φ00)）、darboux_from_matrix(R)
  euler_to_darboux(α, β, γ)（(5.5.2) 閉式）：(50, 60, 70) → [80.00, 33.69, 128.68]（(5.3.8)）；darboux_to_euler(φ, θ, ω)（γ − α = 180 − 2φ）
  darboux_sequence(φ, θ, ωs)（Fig. 5.3.5：連續的 SU(2) 分支，β 在 360 < ω < 720 為負）＝ FIG535_BOOK 的 12 張照片（ω = 0, 60, 128.68, 180, 240, 300,
    360, 420, 488.68, 540, 600, 660 → |−10,0,10⟩, |15.7,32.2,35.7⟩, |50,60,70⟩, |80,67.4,100⟩, |114.8,57.4,134.8⟩, |144.3,32.2,164.3⟩, |170,0,190⟩,
    |195.7,−32.2,215.7⟩, |230,−60,250⟩, |260,−67.4,280⟩, |294.8,−57.4,314.8⟩, |324.3,−32.2,344.3⟩；γ − α = 20°）；EX538、EX539
  su2_euler(α, β, γ)（(5.4.30)；(50, 60, 70)：[[0.4330−0.7500i, −0.4924−0.0868i], [0.4924−0.0868i, 0.4330+0.7500i]]）、su2_darboux(φ, θ, ω)（(5.5.1)）、
    quat_euler、PAULI、Q_UNITS（(5.5.3) q = iσ）
  lab_frame_product(steps)、body_frame_product(steps)（steps = [("z", 度), …]）、lab_op(R)、body_op(R)（位置態 S → RS、S R⁻¹；(5.3.6)–(5.3.7)）
  su2_product(t1, t2)（(5.5.10)）、hamilton_example()：先 R_x(90°) 再 R_y(90°) ＝ 繞 (1, 1, −1)/√3 轉 120°；半角弧 45°、45° 接成 60°（dict：omega、axis、hp、arcs）
§5.4 生成元、J、𝒟^j（H504–H506）
  GX、GY、GZ（(5.4.5b)）、comm(A, B)、gen_exp(axis, ω)、gz_projectors()（(5.4.7) 的三個投影 ＝ GZ_PROJ_BOOK）
  small_commutator(ε=1/3)（Fig. 5.4.2）→ dict(R, angle = −6.31°（繞 z）, axis（偏離 z 13.38°）, second = −ε² = −6.37°, path（從 ŷ 出發 5 點）, approx)
  jmat(j, exact=False) → {Jx, Jy, Jz, Jp, Jm, J2}；J2_JX_BOOK（照印 j = 2 的 𝒟²(J_x)，中間四個 √3/2 應 √6/2）
  wigner_d(j, β)（Wigner 公式）＝ D_schwinger(j, β)（兩個振子的多項式）、wigner_D(j, α, β, γ) ＝ D_expm(j, α, β, γ)；D_axis(j, 軸, ω)；d_sym(j)（sympy）
  D1_SYM（(5.4.33)）、D32_SYM（(5.4.34)）、D2_SYM（(5.4.35) 更正）、D2_BOOK（照印；D2_TYPOS 四格）、d2_book_orthogonality(0.7) = 1.3174（應 1）
    wigner_d(1, 60) = [[0.7500, −0.6124, 0.2500], [0.6124, 0.5000, −0.6124], [0.2500, 0.6124, 0.7500]]；𝒟^{1/2}(360°00) = −1
  d_symmetry_checks(j)（習題 5.4.6）、D_from_u(j, u)（習題 5.4.5）、CIRC_T（圓向量 −(x+iy)/√2, z, (x−iy)/√2）
§5.5B Stern–Gerlach（H508）
  sg_probs(j, n, β)、sg_amps(j, n, β)（m = j…−j）、mean_jz（＝ n cos β）、cone_angle(j)（(5.5.17)：(精確, 更正的近似, 照印的近似)；j = 10⁴：1.1459°）、
  cone_theta(j, m)（Fig. 5.4.4；j = 2：35.26, 65.91, 90.00, 114.09, 144.74°）
    j = 2、n = 1：β = 90°：0.2500, 0.2500, 0, 0.2500, 0.2500；β = 60°：0.4219, 0, 0.2813, 0.2500, 0.0469；β = 180°：全到 m = −1
    j = 6、n = 6、β = 60°：振幅 0.1780, 0.3560, 0.4820, 0.5080, 0.4400, 0.3213, 0.2004, 0.1071, 0.0489, 0.0188, 0.0060, 0.0015, 0.0002
    j = 20、n = 20、β = 60°：⟨m⟩ = 10.0000、σ = 2.7386（= √(j/2) sin β）、最大機率 0.1444（m = 10）
§5.5C 量子轉子（H509）
  rotor_ops(J)（|J m n⟩：lab = J ⊗ 1、body = 1 ⊗ J*）、rotor_levels(kind, J, …)（數值）＝ rotor_formula(kind, J, …)；kind = "sphere"、"symmetric"（C）、
    "asymmetric"（A, B, C）、"zeeman"（g）、"stark"（mu）；asym_matrix_formula(J, A, B, C)（(5.5.54)）；ROTOR_CARD（B = 1、C = 2 或 0.5、g = μ = 0.2）
    球形 0, 2, 6, 12（簡併 1, 9, 25, 49）；長形 C = 2、J = 2：6, 7, 10；扁形 C = 0.5：4, 5.5, 6；Zeeman J = 1：1.8, 2.0, 2.2；Stark J = 2：6.0, 6.2, 6.8
§5.5D 球諧（H510）：ylm(l, m, θ, φ) ＝ ylm_from_D、legendre、addition_theorem(l, β, α, θ, φ)（ℓ = 3、(40, 30)、(70, 100)：−0.4457）、rotated_ylm
§5.6 有限對稱下的分裂（H511–H513）
  chi(j, ω) ＝ chi_exact ＝ 𝒟^j 的跡；O_CLASS_ANGLES = [0, 120, 180, 90, 180]；χ²(O 的五類) = 5, −1, 1, −1, 1
  o3_to_o(ℓ)（(A1, A2, E, T1, T2)，特徵標）＝ o3_to_o_matrices(ℓ)（24 個 𝒟^ℓ 的跡）＝ O3_O_BOOK（(5.6.5b) ℓ = 0…20）；o3_label(ℓ)（"E_g\oplus T_{2g}"）、SPEC_LABELS
  invariant_count(k) → 三種算法（k = 0…8：1, 1, 2, 3, 4, 5, 7, 8, 10）
  d 軌域：COSET_EULER、D2_COL_BOOK（(5.6.11)）、d_orbital_states()（陪集代表六項和 → T2_1..3、E1、E2 ＋ "sums"）＝ d_orbital_projector_full ＝ D_ORBITALS_BOOK
    （T2_3 的和 (1, 0, 0, 0, −1)；照印 T2_SUM_BOOK3 = (1, 0, 0, 0, 1)）、D_ORBITAL_POLY（yz、xz、xy、2z²−x²−y²、x²−y²）
  多極：X4（X^4_0、X^4_{±4}）、V4 = x⁴+y⁴+z⁴ − 3r⁴/5（(100)：0.4000、(111)：−0.2667）、v4_identities()、oct_potential(r)、oct_multipole(r)、
    multipole_coeffs() → x^0_0 = 6.0000、x^4_0 = 3.5000、x^4_{±4} = 2.0917（Q = R = 1）
  SF6：sf6_levels(J, t=1, B=0) → [(E, 標記, 簡併)]（由低到高；J 可半整數）：J = 2：T2 −3.6、E 2.4；J = 3：A2 −38.4、T2 −8.4、T1 15.6；
    sf6_clusters(J)（J = 20：[A1, T1, E]、[T2, T1]、[E, T2, A2] … [T1, E, T2]）
  D6：D6_CLASS_ANGLES、D6_IRREPS（A1 A2 E2 B1 B2 E1）、o3_to_d6(ℓ) ＝ o3_to_d6_mrule(ℓ) ＝ O3_D6_BOOK（(5.6.37)）、d180(N)（(5.6.40)）
  助記輪：WHEEL_D6_EVEN／ODD、WHEEL_O（0 號在左、逆時針經底部；A2 在頂）、WHEEL_O_START、wheel_read、wheel_d6(J)、wheel_o(J)、DIMS、_count
    wheel_o(J) 讀出的次序 ＝ sf6_levels(J) 由高到低（J ≤ 25 全部相符）
§5.7 射線表示（H513–H514）
  SpinAlgebra（.mul、.inv、.table(rows, cols)、.ga_mul、.conj_sum(h)、.class_ops()、.trace_char(j)）
  D6S、D6S_AXES（ρ1 120°、ρ2 60°、ρ3 0°、ρ1' 210°、ρ2' 150°、ρ3' 90°）、D6S_COLS、D6S_ROWS、D6S_TABLE_BOOK（(5.7.4)，121 格 ＝ d6s_table_check()）
    ρ3'ρ2 = −h、ρ2ρ3' = h⁵（(5.7.3)）；d6s_class_algebra()：c2 = h²−h⁴、c3 = h−h⁵、c2² = 2−c2、c3² = 2+c2、c2c3 = c3、c3³ = 3c3、冪等元 P
    D6S_CHARS（類 1, h², h：E1 (2, 1, √3)、E2 (2, 1, −√3)、E3 (2, −2, 0)）、ray_freq、half_correlation("D6", j) ＝ half_d6_mrule(j) ＝ wheel_d6_half(j)
    （j = ½…11/2：E1；E1+E3；E1+E2+E3；E1+2E2+E3；E1+2E2+2E3；2(E1+E2+E3)）；D6S_FREQ_BOOK、D6S_TRACE60_BOOK（照印的錯）、WHEEL_D6_HALF
  D2S、D2S_TABLE_BOOK（(5.7.17)）、D2S_DE_BOOK（(5.7.19)）、cn_spin_chars(n)（−1 的 n 次方根）、dn_spin_algebra(n)、dn_ray_irreps(n)
  OS（O 的射線代數）、OS_TABLE_522_BOOK（(5.7.22) 的 25 格 ＝ os_table_check()）、os_class_algebra()（c_r² = 8 + 2c_r、c_rc_R = 4c_R、c_R² = 6 + 3c_r、c_R³ = 18c_R）、
    OS_CHARS（E1 (2, 1, √2)、E2 (2, 1, −√2)、G (4, −1, 0)）、half_correlation("O", j) ＝ wheel_o_half(j)（j = ½…11/2：E1；G；E2+G；E1+E2+G；E1+2G；E1+E2+2G）、
    OS_FREQ_BOOK（照印 5/2 列少 G）、WHEEL_O_HALF、WHEEL_O_HALF_START（讀出的次序 ＝ 半整數 sf6_levels 由高到低）
§5.8 R4 與 U3（H515–H516）
  r4_generators()（(5.8.15) 4×4）、R4_BASIS、r4_regular(name) ＝ R4_REG_LZ_BOOK／R4_REG_KZ_BOOK（(5.8.17b)、(5.8.18)）、r4_ladders()（M±、N±，更正的 (5.8.20a)）
  coulomb_shell(n)（M = J ⊗ 1、N = 1 ⊗ J、L = M + N、K = M − N，n² 維）、coulomb_energy(n) = −(k²m/2)/n²（k = −1、m = 1：−0.5000, −0.1250, −0.0556, −0.0312）、
    coulomb_weights(j_m, j_n)（Fig. 5.8.3）、kepler_check()（古典軌道：ε 守恆、ε² = 1 + 2HL²/k²m = 0.2809）、runge_lenz_poisson()、quantum_xr_p()、
    quantum_eps_denominator()（(5.8.6)、(5.8.8) 的更正以試驗函數核對）
  E3(i, j)、U3_G、U3_K（(5.8.37)–(5.8.39)）、su3_cartan(N=6)（r(1) = (0.2887, 0.5000)、r(2) = (0.2887, −0.5000)、r(3) = (0.5774, 0)）、
    su3_weights(n1, n2)（半標準 Young 表；權 = (a1 − a3, (a1 − 2a2 + a3)/√3)）、su3_dim、osc_shell(N)、osc_fock_levels(3)（1.5, 2.5, 3.5, 4.5；1, 3, 6, 10）
  COULOMB_CLOSURES = [2, 10, 28, 60]、OSC_CLOSURES = [2, 8, 20, 40, 70]、MAGIC_ATOMIC_BOOK、MAGIC_NUCLEAR_BOOK
附錄 E 與習題：appendix_e_check(j)、problem_511、problem_512（φ = 120°：(111) 軸的非正交「跡」3）、problem_513（ω = 0, 60, 90, 120, 180；D5h、D7h、D4d 不行）、
  dn_chars(n)、dinfh_to_dnh(n) ＝ dinfh_to_dnh_mrule(n)（習題 5.2.1）、problem_533（Euler (45, 54.7356, −90)、軸角 [157.5, 126.4751, 69.7356]）、problem_542、
  problem_553（R1²：(180, 180, 0)、[0, 90, 180]；r1：(0, 90, 90)、[45, 54.74, 120]；i1：(0, 90, 180)、[0, 45, 180]；R1：(−90, 90, 90)、[0, 90, 90]）、
  problem_554（(a)(b)(c) 的 J = 1、2 能階）、problem_555（對稱陀螺的 Euler 形式）、problem_561（O3↓D6 到 ℓ = 12）、o3_to_y(ℓ) ＝ o3_to_y_mats、
  o3_to_td(ℓ, ±1) ＝ o3_to_td_chars（習題 5.6.2）
card_values() ＝ CARD_EXPECTED（plan_ch05 §6 卡片上的每個數字）；check()：每個結果兩種以上獨立算法（約 5 秒）；summary()。
直接執行：python utils/sym_rot.py（在 codes/ 裡）印一句成功訊息；--summary 印卡片數字。
"""
import itertools
import math
from fractions import Fraction
from functools import lru_cache

import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp
from scipy.linalg import expm
from scipy.special import lpmv

try:
    from . import sym_mat as SM
    from . import sym_nonabel as SN
    from . import sym_oct as SO
except ImportError:                      # 直接執行 python utils/sym_rot.py
    import sym_mat as SM
    import sym_nonabel as SN
    import sym_oct as SO

R_ = sp.Rational
IU = sp.I
SQ2, SQ3, SQ6 = sp.sqrt(2), sp.sqrt(3), sp.sqrt(6)
BETA, PHI = sp.symbols("beta phi", real=True)


def _close(a, b, tol: float = 1e-9) -> bool:
    return np.allclose(np.asarray(a, dtype=complex), np.asarray(b, dtype=complex), atol=tol, rtol=0)


def _f4(x) -> str:
    return f"{float(x):.4f}"


def _f2(x) -> str:
    return f"{float(x):.2f}"


def _num(M) -> np.ndarray:
    return np.array(sp.Matrix(M).evalf(), dtype=complex)


def _wrap(a: float) -> float:
    """角度（度）收進 (−180, 180]。"""
    a = (a + 180.0) % 360.0 - 180.0
    return 180.0 if abs(a + 180.0) < 1e-9 else a


# ======================================================================== §5.1–5.2 O3、O2、R2、D∞h（H501）
INV3 = -np.eye(3)                                    # 反演 I
SIGMA_XY = np.diag([1.0, 1.0, -1.0])                 # 對 xy 平面的反射
I_Z = SIGMA_XY @ INV3                                 # i_z = σ_xy I = 繞 z 的 180°


def n_params(n: int) -> int:
    """(5.1.12)：O_n 的參數個數 n(n−1)/2 = n² − n(n+1)/2。"""
    return n * (n - 1) // 2


N_PARAMS = [n_params(n) for n in range(1, 7)]        # O1…O6：0, 1, 3, 6, 10, 15


def o2_rot(phi_deg) -> np.ndarray:
    """(5.1.10a)：det = +1 的轉動。"""
    c, s = math.cos(math.radians(phi_deg)), math.sin(math.radians(phi_deg))
    return np.array([[c, -s], [s, c]])


def o2_refl(phi_deg) -> np.ndarray:
    """(5.1.10b)：det = −1；＝ o2_rot(φ)·diag(1, −1) ＝ 對斜角 φ/2 的直線反射。"""
    return o2_rot(phi_deg) @ np.diag([1.0, -1.0])


def is_orthogonal(M) -> bool:
    M = np.asarray(M, float)
    return np.allclose(M.T @ M, np.eye(M.shape[0])) and np.allclose(M @ M.T, np.eye(M.shape[0]))


def o2_standing(m: int) -> sp.Matrix:
    """(5.2.3)：T·diag(e^{imφ}, e^{−imφ})·T† ＝ [[cos mφ, −sin mφ], [sin mφ, cos mφ]]（φ = PHI 符號）。"""
    T = sp.Matrix([[1, 1], [-IU, IU]]) / SQ2
    D = sp.diag(sp.exp(IU * m * PHI), sp.exp(-IU * m * PHI))
    return sp.simplify((T * D * T.H).applyfunc(lambda e: sp.expand(e.rewrite(sp.cos))))


O2_LABELS = [("A_1", r"\Sigma^+", 0), ("A_2", r"\Sigma^-", 0), ("E_1", r"\Pi", 1), ("E_2", r"\Delta", 2), ("E_3", r"\Phi", 3)]
O2_MOVING_SIGMA_MINUS_BOOK = -1                      # (5.2.4) 行進波表照印 𝒟^{Σ−}(φ) = −1（應為 +1，筆誤）

DINFH_CLASSES = ["1", r"R(\gamma00)", r"R(-\gamma\pi\gamma)", "I", r"IR(\gamma00)", r"IR(-\gamma\pi\gamma)"]
DINFH_IRREPS = [r"\Sigma^+_g", r"\Sigma^+_u", r"\Sigma^-_g", r"\Sigma^-_u", r"\Pi_g", r"\Pi_u", r"\Delta_g", r"\Delta_u", r"\Phi_g"]


def dinfh_chars(gamma_deg) -> dict:
    """Fig. 5.2.2 的 D∞h 特徵標（行 = DINFH_CLASSES）。Σ^± 的 ± 是對鏡面 IR(−γπγ) 的號；橫向 π 轉 R(−γπγ) 的號 ＝ (±)×宇稱
    （書上最後一行的表頭漏了 I）。"""
    out = {}
    for lab in DINFH_IRREPS:
        par = 1 if lab.endswith("g") else -1
        if lab.startswith(r"\Sigma"):
            pm = 1 if "+" in lab else -1
            dim, rot, c2 = 1, 1.0, pm * par
        else:
            m = {r"\Pi": 1, r"\Delta": 2, r"\Phi": 3}[lab.split("_")[0]]
            dim, rot, c2 = 2, 2 * math.cos(math.radians(m * gamma_deg)), 0
        out[lab] = [dim, rot, c2, par * dim, par * rot, par * c2]
    return out


def dinfh_trace_check(gamma_deg: float = 40.0) -> bool:
    """第二種算法：向量 (x, y, z) ＝ Σ⁺_u ⊕ Π_u、二次式 (x²+y², z², xz, yz, x²−y², xy) ＝ 2Σ⁺_g ⊕ Π_g ⊕ Δ_g 的跡。"""
    ops = _dinfh_ops(gamma_deg)
    ch = dinfh_chars(gamma_deg)
    ok = all(abs(np.trace(g) - (ch[r"\Sigma^+_u"][k] + ch[r"\Pi_u"][k])) < 1e-9 for k, g in enumerate(ops))
    for k, g in enumerate(ops):
        S = np.kron(g, g)                                       # 二階張量 → 對稱部分的跡 = (tr² + tr(g²))/2
        sym_tr = (np.trace(g) ** 2 + np.trace(g @ g)) / 2
        ok &= abs(sym_tr - (2 * ch[r"\Sigma^+_g"][k] + ch[r"\Pi_g"][k] + ch[r"\Delta_g"][k])) < 1e-9
        assert S.shape == (9, 9)
    return ok


def _dinfh_ops(gamma_deg):
    """D∞h 的六個代表操作（3×3）：1、R_z(γ)、x 軸 180°、I、I R_z(γ)、I·(x 軸 180°)（＝ 對 yz 面… 的鏡面）。"""
    Rz = rz(gamma_deg)
    C2x = np.diag([1.0, -1.0, -1.0])
    return [np.eye(3), Rz, C2x, INV3, INV3 @ Rz, INV3 @ C2x]


# ======================================================================== §5.3 Euler、Darboux、四元數（H502–H503、H507）
def rz(a) -> np.ndarray:
    c, s = math.cos(math.radians(a)), math.sin(math.radians(a))
    return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1.0]])


def ry(b) -> np.ndarray:
    c, s = math.cos(math.radians(b)), math.sin(math.radians(b))
    return np.array([[c, 0, s], [0, 1.0, 0], [-s, 0, c]])


def rx(a) -> np.ndarray:
    c, s = math.cos(math.radians(a)), math.sin(math.radians(a))
    return np.array([[1.0, 0, 0], [0, c, -s], [0, s, c]])


def euler_matrix(a, b, c) -> np.ndarray:
    """(5.3.3)、(5.5.76)：R(αβγ) = R(α00)R(0β0)R(00γ) = R_z(α)R_y(β)R_z(γ)（度；最右邊先做）。"""
    return rz(a) @ ry(b) @ rz(c)


def euler_matrix_sym() -> sp.Matrix:
    """(5.5.76) 的 3×3（符號 α, β, γ）。"""
    a, b, c = sp.symbols("alpha beta gamma", real=True)
    Z = lambda t: sp.Matrix([[sp.cos(t), -sp.sin(t), 0], [sp.sin(t), sp.cos(t), 0], [0, 0, 1]])   # noqa: E731
    Y = sp.Matrix([[sp.cos(b), 0, sp.sin(b)], [0, 1, 0], [-sp.sin(b), 0, sp.cos(b)]])
    return sp.expand(Z(a) * Y * Z(c))


def euler_from_matrix(Rm) -> tuple:
    """3×3 → (α, β, γ)（度），0 ≤ β ≤ 180（(5.3.2)）；β = 0 或 180 時只有 α ± γ 有意義，取 γ = 0。"""
    Rm = np.asarray(Rm, float)
    b = math.degrees(math.acos(max(-1.0, min(1.0, Rm[2, 2]))))
    if math.sin(math.radians(b)) > 1e-9:
        a = math.degrees(math.atan2(Rm[1, 2], Rm[0, 2]))
        c = math.degrees(math.atan2(Rm[2, 1], -Rm[2, 0]))
    elif Rm[2, 2] > 0:
        a, c = math.degrees(math.atan2(Rm[1, 0], Rm[0, 0])), 0.0
    else:
        a, c = math.degrees(math.atan2(-Rm[1, 0], -Rm[0, 0])), 0.0
    return _wrap(a), b, _wrap(c)


def euler_alt(a, b, c, sign: int = 1) -> tuple:
    """(5.3.1)：同一個位置的另一組 Euler 角 (α ± 180, −β, γ ± 180)（sign = ±1 選 ±）。"""
    return a + 180 * sign, -b, c + 180 * sign


def axis_vec(phi, theta) -> np.ndarray:
    """(5.5.5c)：Darboux 軸 (cos φ sin θ, sin φ sin θ, cos θ)。"""
    p, t = math.radians(phi), math.radians(theta)
    return np.array([math.cos(p) * math.sin(t), math.sin(p) * math.sin(t), math.cos(t)])


def darboux_matrix(phi, theta, omega) -> np.ndarray:
    """R[φθω]（軸角，Rodrigues ＝ 四元數 SN.quat）。"""
    return SN.quat_to_rot(SN.quat(axis_vec(phi, theta), omega))


def darboux_matrix_cranks(phi, theta, omega) -> np.ndarray:
    """§5.5A 的五個曲柄：R(φ00)R(0θ0)R(ω00)R(0−θ0)R(−φ00)。"""
    return rz(phi) @ ry(theta) @ rz(omega) @ ry(-theta) @ rz(-phi)


def darboux_from_matrix(Rm) -> tuple:
    """3×3 → [φ, θ, ω]（0 ≤ ω ≤ 180）。"""
    ax, ang = SN.axis_angle(Rm)
    th = math.degrees(math.acos(max(-1.0, min(1.0, ax[2]))))
    ph = math.degrees(math.atan2(ax[1], ax[0])) if math.hypot(ax[0], ax[1]) > 1e-12 else 0.0
    return ph, th, ang


def su2_euler(a, b, c) -> np.ndarray:
    """(5.4.30)：𝒟^{1/2}(αβγ)（度）。"""
    A, B, C = (math.radians(t) / 2 for t in (a, b, c))
    return np.array([[np.exp(-1j * (A + C)) * math.cos(B), -np.exp(-1j * (A - C)) * math.sin(B)],
                     [np.exp(1j * (A - C)) * math.sin(B), np.exp(1j * (A + C)) * math.cos(B)]])


def su2_darboux(phi, theta, omega) -> np.ndarray:
    """(5.5.1)／(5.5.4)：𝒟^{1/2}[φθω] = cos(ω/2)1 − i sin(ω/2) ω̂·σ。"""
    n = axis_vec(phi, theta)
    w = math.radians(omega) / 2
    return math.cos(w) * np.eye(2) - 1j * math.sin(w) * (n[0] * PAULI[0] + n[1] * PAULI[1] + n[2] * PAULI[2])


PAULI = (np.array([[0, 1], [1, 0]], complex), np.array([[0, -1j], [1j, 0]]), np.array([[1, 0], [0, -1]], complex))


def quat_euler(a, b, c) -> np.ndarray:
    """Euler 角的四元數 q(z, α)·q(y, β)·q(z, γ)（SN.qmul）；SN.su2(q) ＝ su2_euler。"""
    return SN.qmul(SN.quat((0, 0, 1), a), SN.qmul(SN.quat((0, 1, 0), b), SN.quat((0, 0, 1), c)))


def euler_to_darboux(a, b, c) -> tuple:
    """(5.5.2) 的閉式：(α,β,γ) → [φ, θ, ω]（0 ≤ ω < 360）。"""
    S, Dd, B = (math.radians(a) + math.radians(c)) / 2, (math.radians(c) - math.radians(a)) / 2, math.radians(b) / 2
    cw = math.cos(S) * math.cos(B)
    v = np.array([math.sin(Dd) * math.sin(B), math.cos(Dd) * math.sin(B), math.sin(S) * math.cos(B)])   # = sin(ω/2) n̂
    sw = float(np.linalg.norm(v))
    om = math.degrees(2 * math.atan2(sw, cw))
    if sw < 1e-12:
        return 0.0, 0.0, om
    n = v / sw
    return math.degrees(math.atan2(n[1], n[0])), math.degrees(math.acos(max(-1.0, min(1.0, n[2])))), om


def darboux_to_euler(phi, theta, omega) -> tuple:
    """(5.5.2) 的反算（取 cos(β/2) ≥ 0、β 可為負：書上 Fig. 5.3.5 的分支）：γ − α = 180 − 2φ。"""
    w = math.radians(omega) / 2
    t = math.radians(theta)
    sb = math.sin(w) * math.sin(t)
    B = math.asin(max(-1.0, min(1.0, sb)))
    S = math.atan2(math.sin(w) * math.cos(t), math.cos(w))
    Dd = math.radians(90 - phi)
    return math.degrees(S - Dd), math.degrees(2 * B), math.degrees(S + Dd)


def darboux_sequence(phi, theta, omegas) -> list:
    """Fig. 5.3.5：R[φθω]|000⟩ 隨 ω 連續變化的 Euler 角（SU(2) 分支：β 在 360° < ω < 720° 為負，720° 才回到原狀）。"""
    out, prev = [], None
    for om in omegas:
        a, b, c = darboux_to_euler(phi, theta, om)
        S = a + c
        if prev is not None:
            S += 360.0 * round((prev - S) / 360.0)
        prev = S
        Dd = c - a
        out.append(((S - Dd) / 2, b, (S + Dd) / 2))
    return out


EX538 = (50, 60, 70)                                  # (5.3.8)
EX539 = (-130, -60, -110)                             # (5.3.9)：同一位置的另一組（ω + 360°）
FIG535_BOOK = {0: (-10, 0, 10), 60: (15.7, 32.20, 35.7), 128.68: (50, 60, 70), 180: (80, 67.4, 100), 240: (114.8, 57.4, 134.8),
               300: (144.3, 32.2, 164.3), 360: (170, 0, 190), 420: (195.7, -32.2, 215.7), 488.68: (230, -60, 250),
               540: (260, -67.4, 280), 600: (294.8, -57.4, 314.8), 660: (324.3, -32.2, 344.3)}


def lab_frame_product(steps) -> np.ndarray:
    """依序繞「實驗室」的軸轉：steps = [(軸, 度), …]，先做第一個 ⇒ R = R_n ⋯ R_1。"""
    R = np.eye(3)
    for ax, d in steps:
        R = _rot_ax(ax, d) @ R
    return R


def body_frame_product(steps) -> np.ndarray:
    """依序繞「本體」（跟著轉的）軸轉 ⇒ R = R_1 ⋯ R_n（順序反過來，(5.3.5)）。"""
    R = np.eye(3)
    for ax, d in steps:
        R = R @ _rot_ax(ax, d)
    return R


_AX = {"x": (1, 0, 0), "y": (0, 1, 0), "z": (0, 0, 1)}


def _rot_ax(ax, d) -> np.ndarray:
    return SN.quat_to_rot(SN.quat(_AX.get(ax, ax), d))


def lab_op(Rm):
    """實驗室算符作用在位置態 |S⟩（S = 3×3）：|S⟩ → |R S⟩。"""
    return lambda S: np.asarray(Rm) @ S


def body_op(Rm):
    """本體算符 R̄：|S⟩ → |S R⁻¹⟩（＝ 把宇宙反轉，(5.3.7)）。"""
    return lambda S: S @ np.asarray(Rm).T


def hamilton_example() -> dict:
    """H507 的例子：先 R_x(90°) 再 R_y(90°) ＝ 繞 (1, 1, −1)/√3 轉 120°（(5.5.10) ＝ 四元數 ＝ 半角弧，SN.hamilton_product）。"""
    hp = SN.hamilton_product((1, 0, 0), 90, (0, 1, 0), 90)
    w1, w2 = np.array([1.0, 0, 0]), np.array([0, 1.0, 0])
    c = math.cos(math.radians(45)) ** 2 - math.sin(math.radians(45)) ** 2 * (w2 @ w1)          # (5.5.10a)
    v = (math.cos(math.radians(45)) * math.sin(math.radians(45)) * (w1 + w2)
         + math.sin(math.radians(45)) ** 2 * np.cross(w2, w1))                                   # (5.5.10b)
    om = math.degrees(2 * math.acos(c))
    arcs = tuple(math.degrees(math.acos(np.clip(u @ w, -1, 1))) for u, w in ((hp["N1"], hp["P"]), (hp["P"], hp["N2p"]), (hp["N1"], hp["N2p"])))
    return dict(omega=om, axis=v / np.linalg.norm(v), hp=hp, arcs=arcs)


def su2_product(t1, t2) -> tuple:
    """(5.5.6)–(5.5.10)：R[t2]R[t1]（t = (φ, θ, ω)）→ (cos(ω''/2), ω̂'' sin(ω''/2))（先 t1）。"""
    w1, w2 = math.radians(t1[2]) / 2, math.radians(t2[2]) / 2
    n1, n2 = axis_vec(t1[0], t1[1]), axis_vec(t2[0], t2[1])
    c = math.cos(w2) * math.cos(w1) - math.sin(w2) * math.sin(w1) * (n2 @ n1)
    v = math.cos(w2) * math.sin(w1) * n1 + math.cos(w1) * math.sin(w2) * n2 + math.sin(w2) * math.sin(w1) * np.cross(n2, n1)
    return c, v


# ======================================================================== §5.4 生成元、J 矩陣、𝒟^j（H504–H506）
GX = sp.Matrix([[0, 0, 0], [0, 0, -1], [0, 1, 0]])    # (5.4.5b)：書上的 𝒢_x, 𝒢_y, 𝒢_z
GY = sp.Matrix([[0, 0, 1], [0, 0, 0], [-1, 0, 0]])
GZ = sp.Matrix([[0, -1, 0], [1, 0, 0], [0, 0, 0]])


def comm(A, B):
    return A * B - B * A if isinstance(A, sp.MatrixBase) else A @ B - B @ A


def gen_exp(axis, w_deg) -> np.ndarray:
    """e^{ω n̂·𝒢}（(5.4.6)；scipy expm）＝ 繞 n̂ 逆時針轉 ω。"""
    n = np.asarray(axis, float) / np.linalg.norm(axis)
    G = n[0] * _num(GX).real + n[1] * _num(GY).real + n[2] * _num(GZ).real
    return expm(math.radians(w_deg) * G)


def gz_projectors() -> dict:
    """(5.4.7)：𝒢_z 的譜分解 {i: 𝒫₊, −i: 𝒫₋, 0: 𝒫₀}（SM.projectors）；e^{ω𝒢_z} = e^{iω}𝒫₊ + e^{−iω}𝒫₋ + 𝒫₀。"""
    return SM.projectors(GZ)


GZ_PROJ_BOOK = {IU: sp.Matrix([[R_(1, 2), IU / 2, 0], [-IU / 2, R_(1, 2), 0], [0, 0, 0]]),
                -IU: sp.Matrix([[R_(1, 2), -IU / 2, 0], [IU / 2, R_(1, 2), 0], [0, 0, 0]]),
                0: sp.Matrix([[0, 0, 0], [0, 0, 0], [0, 0, 1]])}

SMALL_EPS = 1 / 3                                    # Fig. 5.4.2 的 ε ≈ 1/3


def small_commutator(eps: float = SMALL_EPS) -> dict:
    """(5.4.8)–(5.4.10)、Fig. 5.4.2：R = e^{−ε𝒢y}e^{−ε𝒢x}e^{ε𝒢y}e^{ε𝒢x}（最右先做）。
    → dict(R, angle（度，繞 z 的有號轉角）, axis, second（二階預測 −ε² rad 的度數）, path（從 ŷ 出發的 5 個點）)。"""
    Gx, Gy = _num(GX).real, _num(GY).real
    seq = [expm(eps * Gx), expm(eps * Gy), expm(-eps * Gx), expm(-eps * Gy)]
    R = seq[3] @ seq[2] @ seq[1] @ seq[0]
    pts = [np.array([0.0, 1.0, 0.0])]
    for M in seq:
        pts.append(M @ pts[-1])
    ax, ang = SN.axis_angle(R)
    sgn = 1.0 if ax[2] >= 0 else -1.0
    return dict(R=R, angle=sgn * ang, axis=sgn * ax, second=-math.degrees(eps ** 2), path=pts,
                approx=np.eye(3) - eps ** 2 * _num(comm(GX, GY)).real)


def _ms(j) -> list:
    tj = int(round(2 * float(j)))
    return [Fraction(tj - 2 * k, 2) for k in range(tj + 1)]


def jmat(j, exact: bool = False) -> dict:
    """(5.4.21)–(5.4.25)：基底 m = j, j−1, …, −j（列＝行同序），ħ = 1。→ {Jx, Jy, Jz, Jp, Jm, J2}（numpy 或 sympy）。"""
    ms = _ms(j)
    jj = Fraction(int(round(2 * float(j))), 2)
    n = len(ms)
    if exact:
        Jp = sp.zeros(n)
        for k in range(1, n):
            m = ms[k]
            Jp[k - 1, k] = sp.sqrt(sp.Rational((jj - m) * (jj + m + 1)))
        Jz = sp.diag(*[sp.Rational(m) for m in ms])
        Jm = Jp.T
        Jx, Jy = (Jp + Jm) / 2, (Jp - Jm) / (2 * IU)
        return dict(Jx=Jx, Jy=sp.simplify(Jy), Jz=Jz, Jp=Jp, Jm=Jm, J2=sp.simplify(Jx * Jx + Jy * Jy + Jz * Jz))
    Jp = np.zeros((n, n))
    for k in range(1, n):
        m = ms[k]
        Jp[k - 1, k] = math.sqrt(float((jj - m) * (jj + m + 1)))
    Jz = np.diag([float(m) for m in ms])
    Jm = Jp.T.copy()
    Jx, Jy = (Jp + Jm) / 2, (Jp - Jm) / 2j
    return dict(Jx=Jx.astype(complex), Jy=Jy, Jz=Jz.astype(complex), Jp=Jp, Jm=Jm, J2=(Jx @ Jx + Jy @ Jy + Jz @ Jz))


J2_JX_BOOK = sp.Matrix([[0, 1, 0, 0, 0], [1, 0, SQ3 / 2, 0, 0], [0, SQ3 / 2, 0, SQ3 / 2, 0], [0, 0, SQ3 / 2, 0, 1], [0, 0, 0, 1, 0]])   # 照印（中間四個應 √6/2）


@lru_cache(maxsize=None)
def _fact(n: int) -> int:
    return math.factorial(n)


def wigner_d(j, beta_deg) -> np.ndarray:
    """(5.4.45)（α = γ = 0）：d^j_{mn}(β)，列 m、行 n（都從 j 到 −j）。Wigner 的公式（整數運算：j + m 等都是整數）。"""
    tj = int(round(2 * float(j)))
    c, s = math.cos(math.radians(beta_deg) / 2), math.sin(math.radians(beta_deg) / 2)
    n = tj + 1
    jm = [tj - k for k in range(n)]                  # j + m（m = j − k）
    d = np.zeros((n, n))
    for a in range(n):
        jpm, jmm = jm[a], tj - jm[a]                 # j + m、j − m
        for b in range(n):
            jpn, jmn = jm[b], tj - jm[b]
            pre = math.sqrt(_fact(jpm) * _fact(jmm) * _fact(jpn) * _fact(jmn))
            nm = jpn - jpm                           # n − m
            tot = 0.0
            for k in range(max(0, -nm), min(jpm, jmn) + 1):
                den = _fact(jpm - k) * _fact(jmn - k) * _fact(nm + k) * _fact(k)
                tot += (-1) ** k / den * c ** (tj - nm - 2 * k) * s ** (nm + 2 * k)
            d[a, b] = pre * tot
    return d


def wigner_D(j, a, b, c) -> np.ndarray:
    """(5.4.45)：𝒟^j_{mn}(αβγ) = e^{−imα} d^j_{mn}(β) e^{−inγ}（度）。"""
    ms = np.array([float(m) for m in _ms(j)])
    return np.exp(-1j * ms * math.radians(a))[:, None] * wigner_d(j, b) * np.exp(-1j * ms * math.radians(c))[None, :]


def D_expm(j, a, b, c) -> np.ndarray:
    """(5.4.28)：e^{−iαJz}e^{−iβJy}e^{−iγJz}（scipy expm）。"""
    J = jmat(j)
    return expm(-1j * math.radians(a) * J["Jz"]) @ expm(-1j * math.radians(b) * J["Jy"]) @ expm(-1j * math.radians(c) * J["Jz"])


def D_axis(j, axis, w_deg) -> np.ndarray:
    """繞任意軸轉 ω 的 𝒟^j：e^{−iω n̂·J}。"""
    J = jmat(j)
    n = np.asarray(axis, float) / np.linalg.norm(axis)
    return expm(-1j * math.radians(w_deg) * (n[0] * J["Jx"] + n[1] * J["Jy"] + n[2] * J["Jz"]))


def D_schwinger(j, beta_deg) -> np.ndarray:
    """(5.4.36)–(5.4.44)：Schwinger 的兩個振子。把 (c a↑ + s a↓)^{j+n}(−s a↑ + c a↓)^{j−n} 展開成 a↑、a↓ 的多項式
    （numpy 多項式乘法），a↑^{j+m} a↓^{j−m} 的係數乘 √((j+m)!(j−m)!/((j+n)!(j−n)!)) 就是 d^j_{mn}(β)。"""
    ms = _ms(j)
    jj = Fraction(int(round(2 * float(j))), 2)
    c, s = math.cos(math.radians(beta_deg) / 2), math.sin(math.radians(beta_deg) / 2)
    n = len(ms)
    d = np.zeros((n, n))
    for b, nn in enumerate(ms):
        up, dn = int(jj + nn), int(jj - nn)
        poly = np.array([1.0])                       # 以 a↑ 的次方為指標（a↓ 的次方 = 總次數 − 它）
        for _ in range(up):
            poly = np.convolve(poly, [s, c])         # c a↑ + s a↓：a↑¹ 的係數 c、a↑⁰ 的係數 s
        for _ in range(dn):
            poly = np.convolve(poly, [c, -s])        # −s a↑ + c a↓
        for a, m in enumerate(ms):
            k = int(jj + m)
            d[a, b] = poly[k] * math.sqrt(_fact(k) * _fact(int(jj - m)) / (_fact(up) * _fact(dn)))
    return d


def d_sym(j) -> sp.Matrix:
    """d^j(β) 的符號矩陣（Wigner 公式，cos(β/2)、sin(β/2)）。"""
    ms = _ms(j)
    jj = Fraction(int(round(2 * float(j))), 2)
    c, s = sp.cos(BETA / 2), sp.sin(BETA / 2)
    n = len(ms)
    M = sp.zeros(n)
    for a, m in enumerate(ms):
        for b, nn in enumerate(ms):
            pre = sp.sqrt(_fact(int(jj + m)) * _fact(int(jj - m)) * _fact(int(jj + nn)) * _fact(int(jj - nn)))
            tot = 0
            for k in range(0, int(2 * jj) + 1):
                e1, e2, e3 = jj + m - k, jj - nn - k, nn - m + k
                if min(e1, e2, e3) < 0:
                    continue
                tot += sp.Integer(-1) ** k / (_fact(int(e1)) * _fact(int(e2)) * _fact(int(e3)) * _fact(k)) \
                    * c ** int(2 * jj + m - nn - 2 * k) * s ** int(nn - m + 2 * k)
            M[a, b] = pre * tot
    return M


_cb, _sb = sp.cos(BETA), sp.sin(BETA)
_P, _Mn = (1 + _cb) / 2, (1 - _cb) / 2
D1_SYM = sp.Matrix([[_P, -_sb / SQ2, _Mn], [_sb / SQ2, _cb, -_sb / SQ2], [_Mn, _sb / SQ2, _P]])              # (5.4.33) 中間的矩陣
_c2, _s2 = sp.cos(BETA / 2), sp.sin(BETA / 2)
D32_SYM = sp.Matrix([[_c2 ** 3, -SQ3 * _c2 ** 2 * _s2, SQ3 * _c2 * _s2 ** 2, -_s2 ** 3],
                     [SQ3 * _c2 ** 2 * _s2, _c2 ** 3 - 2 * _c2 * _s2 ** 2, _s2 ** 3 - 2 * _c2 ** 2 * _s2, SQ3 * _c2 * _s2 ** 2],
                     [SQ3 * _c2 * _s2 ** 2, -_s2 ** 3 + 2 * _c2 ** 2 * _s2, _c2 ** 3 - 2 * _c2 * _s2 ** 2, -SQ3 * _c2 ** 2 * _s2],
                     [_s2 ** 3, SQ3 * _c2 * _s2 ** 2, SQ3 * _c2 ** 2 * _s2, _c2 ** 3]])                             # (5.4.34)
_r38, _r32 = sp.sqrt(R_(3, 8)), sp.sqrt(R_(3, 2))
D2_SYM = sp.Matrix([                                                                                              # (5.4.35) 更正
    [_P ** 2, -_P * _sb, _r38 * _sb ** 2, -_Mn * _sb, _Mn ** 2],
    [_P * _sb, _P * (2 * _cb - 1), -_r32 * _sb * _cb, _Mn * (2 * _cb + 1), -_Mn * _sb],
    [_r38 * _sb ** 2, _r32 * _sb * _cb, (3 * _cb ** 2 - 1) / 2, -_r32 * _sb * _cb, _r38 * _sb ** 2],
    [_Mn * _sb, _Mn * (2 * _cb + 1), _r32 * _sb * _cb, _P * (2 * _cb - 1), -_P * _sb],
    [_Mn ** 2, _Mn * _sb, _r38 * _sb ** 2, _P * _sb, _P ** 2]])
D2_BOOK = D2_SYM.copy()                                                                                           # 照印（4 個錯）
D2_BOOK[0, 1], D2_BOOK[0, 3], D2_BOOK[2, 3], D2_BOOK[3, 0] = _P * _sb, _P * _sb, _r32 * _sb * _cb, _P * _sb
D2_TYPOS = [(0, 1), (0, 3), (2, 3), (3, 0)]                                                                       # 0 起算的 (列, 行)


def d2_book_orthogonality(beta: float = 0.7) -> float:
    """照印的 (5.4.35)：β = 0.7 rad 時 𝒟𝒟ᵀ 的 (1,1) 元（應為 1）。"""
    M = np.array(D2_BOOK.subs(BETA, beta).evalf(), dtype=float)
    return float((M @ M.T)[0, 0])


def d_symmetry_checks(j, b: float = 37.0) -> dict:
    """習題 5.4.6 的相位對稱（全部數值核對）。"""
    d, dm, dpi = wigner_d(j, b), wigner_d(j, -b), wigner_d(j, 180 - b)
    ms = _ms(j)
    jj = Fraction(int(round(2 * float(j))), 2)
    n = len(ms)
    ok = dict(neg=True, transpose=True, minus_beta=True, pi_minus=True)
    for a in range(n):
        for c in range(n):
            m, nn = ms[a], ms[c]
            sg = (-1) ** int(m - nn)
            ok["neg"] &= abs(d[n - 1 - a, n - 1 - c] - sg * d[a, c]) < 1e-9
            ok["transpose"] &= abs(d[c, a] - sg * d[a, c]) < 1e-9
            ok["minus_beta"] &= abs(dm[a, c] - d[c, a]) < 1e-9
            ok["pi_minus"] &= abs(dpi[a, c] - (-1) ** int(jj + m) * d[a, n - 1 - c]) < 1e-9
    return ok


# ======================================================================== §5.5B Stern–Gerlach（H508）
def sg_probs(j, n, beta_deg) -> np.ndarray:
    """|𝒟^j_{mn}(0β0)|²，m = j…−j（Fig. 5.5.2–5.5.5）。"""
    k = [float(x) for x in _ms(j)].index(float(n))
    return wigner_d(j, beta_deg)[:, k] ** 2


def sg_amps(j, n, beta_deg) -> np.ndarray:
    k = [float(x) for x in _ms(j)].index(float(n))
    return wigner_d(j, beta_deg)[:, k]


def mean_jz(j, n, beta_deg) -> float:
    """(5.5.18)：Σ m |𝒟^j_{mn}|² ＝ n cos β。"""
    ms = np.array([float(m) for m in _ms(j)])
    return float(ms @ sg_probs(j, n, beta_deg))


def cone_angle(j) -> tuple:
    """(5.5.17)：(2 arccos(j/√(j(j+1))), 更正的近似 2 arccos(j/(j+½)), 照印的 arccos(j/(j+½)))（度）。"""
    j = float(j)
    return (math.degrees(2 * math.acos(j / math.sqrt(j * (j + 1)))), math.degrees(2 * math.acos(j / (j + 0.5))),
            math.degrees(math.acos(j / (j + 0.5))))


def cone_theta(j, m) -> float:
    """Fig. 5.4.4 的不確定角：cos θ^j_m = m/√(j(j+1))（度）。"""
    return math.degrees(math.acos(float(m) / math.sqrt(float(j) * (float(j) + 1))))


# ======================================================================== §5.5C 量子轉子（H509）
def rotor_ops(J) -> dict:
    """|J m n⟩ 空間（(2J+1)² 維，kron(m, n)）的實驗室 J_x,J_y,J_z 與本體（Jahn 反轉）J_x̄,J_ȳ,J_z̄ ＝ 1 ⊗ J_i*（(5.5.45)）。"""
    M = jmat(J)
    one = np.eye(M["Jz"].shape[0])
    lab = {k: np.kron(M[k], one) for k in ("Jx", "Jy", "Jz")}
    body = {k: np.kron(one, M[k].conj()) for k in ("Jx", "Jy", "Jz")}
    return dict(lab=lab, body=body, J2=sum(lab[k] @ lab[k] for k in lab))


def _levels(H, tol: float = 1e-7) -> list:
    w = np.sort(np.linalg.eigvalsh((H + H.conj().T) / 2))
    out = []
    for x in w:
        if out and abs(x - out[-1][0]) < tol:
            out[-1][1] += 1
        else:
            out.append([x, 1])
    return [(float(e), d) for e, d in out]


def rotor_levels(kind: str, J, **k) -> list:
    """數值對角化（(2J+1)² 維）→ [(能量, 簡併度)]。kind：
    "sphere"（B J²，(5.5.48)）、"symmetric"（B J² + (C − B)J_z̄²，(5.5.51)）、"asymmetric"（A J_x̄² + B J_ȳ² + C J_z̄²，(5.5.53)）、
    "zeeman"（B J² + g J_z，(5.5.55)）、"stark"（B J² + μ J_z²，(5.5.57)）。"""
    o = rotor_ops(J)
    L, Bd = o["lab"], o["body"]
    if kind == "sphere":
        H = k.get("B", 1) * o["J2"]
    elif kind == "symmetric":
        H = k.get("B", 1) * o["J2"] + (k["C"] - k.get("B", 1)) * Bd["Jz"] @ Bd["Jz"]
    elif kind == "asymmetric":
        H = k["A"] * Bd["Jx"] @ Bd["Jx"] + k["B"] * Bd["Jy"] @ Bd["Jy"] + k["C"] * Bd["Jz"] @ Bd["Jz"]
    elif kind == "zeeman":
        H = k.get("B", 1) * o["J2"] + k["g"] * L["Jz"]
    elif kind == "stark":
        H = k.get("B", 1) * o["J2"] + k["mu"] * L["Jz"] @ L["Jz"]
    else:
        raise ValueError(kind)
    return _levels(H)


def rotor_formula(kind: str, J, **k) -> list:
    """(5.5.49)、(5.5.52)、(5.5.56)、(5.5.58) 的公式（含簡併度）。"""
    B = k.get("B", 1)
    base = B * J * (J + 1)
    d = 2 * J + 1
    if kind == "sphere":
        return [(base, d * d)]
    if kind == "symmetric":
        return sorted([(base + (k["C"] - B) * n * n, d * (1 if n == 0 else 2)) for n in range(J + 1)])
    if kind == "zeeman":
        return sorted([(base + k["g"] * m, d) for m in range(-J, J + 1)])
    if kind == "stark":
        return sorted([(base + k["mu"] * m * m, d * (1 if m == 0 else 2)) for m in range(J + 1)])
    raise ValueError(kind)


def asym_matrix_formula(J, A, B, C) -> np.ndarray:
    """(5.5.54)：H_asym 在 |J n⟩（n = J…−J）的矩陣（三對角跳 2）。"""
    ns = [J - k for k in range(2 * J + 1)]
    n_ = len(ns)
    H = np.zeros((n_, n_))
    for a, n in enumerate(ns):
        H[a, a] = (A + B) / 2 * (J * (J + 1) - n * n) + C * n * n
        if n + 2 <= J:
            H[a - 2, a] = (A - B) / 4 * math.sqrt((J + n + 1) * (J + n + 2) * (J - n) * (J - n - 1))
        if n - 2 >= -J:
            H[a + 2, a] = (A - B) / 4 * math.sqrt((J + n) * (J + n - 1) * (J - n + 1) * (J - n + 2))
    return H


ROTOR_CARD = dict(B=1.0, C_prolate=2.0, C_oblate=0.5, g=0.2, mu=0.2)   # 本系列在 H509 用的常數（書上只給公式）


# ======================================================================== §5.5D 球諧（H510）
def ylm(l, m, theta, phi):
    """(5.5.71)：Y^l_m(θφ)（弧度；含 Condon–Shortley 號，與 (5.6.14) 相同）。"""
    am = abs(m)
    norm = math.sqrt((2 * l + 1) / (4 * math.pi) * _fact(l - am) / _fact(l + am))
    y = norm * lpmv(am, l, np.cos(theta)) * np.exp(1j * am * np.asarray(phi))      # lpmv 已含 (−1)^m
    return y if m >= 0 else (-1) ** am * np.conj(y)


def ylm_from_D(l, m, theta, phi):
    """(5.5.65)：Y^l_m(θφ) = √((2l+1)/4π)·𝒟^{l*}_{m0}(φ, θ, ·)（Euler 角 α = φ、β = θ）。"""
    D = wigner_D(l, math.degrees(phi), math.degrees(theta), 0.0)
    return math.sqrt((2 * l + 1) / (4 * math.pi)) * np.conj(D[[float(x) for x in _ms(l)].index(float(m)), l])


def legendre(l, x) -> float:
    return float(lpmv(0, l, x))


def addition_theorem(l, b_deg, a_deg, th_deg, ph_deg) -> tuple:
    """(5.5.75)：(P_l(cos Θ), (4π/(2l+1)) Σ Y*_m(βα) Y_m(θφ))。"""
    u = axis_vec(a_deg, b_deg)
    v = axis_vec(ph_deg, th_deg)
    lhs = legendre(l, float(u @ v))
    rhs = 4 * math.pi / (2 * l + 1) * sum(np.conj(ylm(l, m, math.radians(b_deg), math.radians(a_deg)))
                                          * ylm(l, m, math.radians(th_deg), math.radians(ph_deg)) for m in range(-l, l + 1))
    return lhs, complex(rhs)


def rotated_ylm(l, m, euler, th, ph) -> tuple:
    """(5.5.70a)：Σ_{m'} 𝒟^l_{m'm}(αβγ) Y_{m'}(θφ) ＝ Y_m 在轉回去的點 R⁻¹r̂ 的值（弧度輸入 θ, φ）。"""
    D = wigner_D(l, *euler)
    ms = [int(x) for x in _ms(l)]
    lhs = sum(D[a, ms.index(m)] * ylm(l, mp, th, ph) for a, mp in enumerate(ms))
    r = euler_matrix(*euler).T @ np.array([math.sin(th) * math.cos(ph), math.sin(th) * math.sin(ph), math.cos(th)])
    th2, ph2 = math.acos(max(-1.0, min(1.0, r[2]))), math.atan2(r[1], r[0])
    return complex(lhs), complex(ylm(l, m, th2, ph2))


# ======================================================================== §5.6 有限對稱下的分裂（H511–H513）
def chi(j, w_deg) -> float:
    """(5.6.4)：χ^j(ω) = sin((j+½)ω)/sin(ω/2)（ω = 0 → 2j+1）。"""
    if abs(math.sin(math.radians(w_deg) / 2)) < 1e-12:
        return (2 * float(j) + 1) * (1 if round(w_deg / 360) % 2 == 0 or float(j) == int(float(j)) else -1)
    return math.sin((float(j) + 0.5) * math.radians(w_deg)) / math.sin(math.radians(w_deg) / 2)


def chi_exact(j, w_deg):
    """Σ_m e^{−imω}（sympy 精確）。"""
    return sp.nsimplify(sp.simplify(sum(sp.cos(sp.Rational(m) * sp.pi * sp.Rational(w_deg) / 180) for m in _ms(j))))


O_CLASS_ANGLES = [0, 120, 180, 90, 180]              # O 的類 1, r, R², R, i 的轉角
O_IRREPS = SO.O_IRREP_NAMES                          # A_1 A_2 E T_1 T_2


def o3_to_o(l) -> tuple:
    """(5.6.1)、(5.6.5b)：𝒟^ℓ↓O 的頻率 (f^{A1}, f^{A2}, f^E, f^{T1}, f^{T2})（特徵標公式，精確）。"""
    tr = [chi_exact(l, w) for w in O_CLASS_ANGLES]
    out = []
    for a in range(5):
        f = sum(SO.O_CHARS[a, k] * SO.O_CLASS_SIZES[k] * tr[k] for k in range(5)) / 24
        out.append(int(sp.nsimplify(f)))
    return tuple(out)


def o3_to_o_matrices(l) -> tuple:
    """第二種算法：24 個 O 元素的 𝒟^ℓ(g) ＝ e^{−iω n̂·J}（軸角 SO.O_TURNS），P 算子的跡。"""
    f = np.zeros(5)
    for g in SO.O_NAMES:
        ax, w = SO.O_TURNS[g]
        tr = np.trace(D_axis(l, ax, w)).real
        k = next(i for i, C in enumerate(SO.O.classes) if g in C)
        f += np.array([float(SO.O_CHARS[a, k]) for a in range(5)]) * tr
    return tuple(int(round(x / 24)) for x in f)


O3_O_BOOK = {0: (1, 0, 0, 0, 0), 1: (0, 0, 0, 1, 0), 2: (0, 0, 1, 0, 1), 3: (0, 1, 0, 1, 1), 4: (1, 0, 1, 1, 1),
             5: (0, 0, 1, 2, 1), 6: (1, 1, 1, 1, 2), 7: (0, 1, 1, 2, 2), 8: (1, 0, 2, 2, 2), 9: (1, 1, 1, 3, 2),
             10: (1, 1, 2, 2, 3), 11: (0, 1, 2, 3, 3), 12: (2, 1, 2, 3, 3), 13: (1, 1, 2, 4, 3), 14: (1, 1, 3, 3, 4),
             15: (1, 2, 2, 4, 4), 16: (2, 1, 3, 4, 4), 17: (1, 1, 3, 5, 4), 18: (2, 2, 3, 4, 5), 19: (1, 2, 3, 5, 5),
             20: (2, 1, 4, 5, 5)}                    # (5.6.5b) 照印（全對）
SPEC_LABELS = "s p d f g h i k l m n o q r t u".split()   # (5.6.5a) 的光譜標記 ℓ = 0…15


def o3_label(l, parity=None) -> str:
    """𝒟^ℓ↓O_h 的 TeX（宇稱 g/u 預設跟 ℓ 的奇偶，如單電子軌域）。"""
    p = parity or ("g" if l % 2 == 0 else "u")
    parts = []
    for f, nm in zip(o3_to_o(l), O_IRREPS):
        if f:
            base = nm.replace("_", "")
            lab = f"{base[0]}_{{{base[1:]}{p}}}" if len(base) > 1 else f"{base}_{p}"
            parts.append((f"{f}" if f > 1 else "") + lab)
    return r"\oplus ".join(parts)


def invariant_count(k: int) -> tuple:
    """(5.6.20)–(5.6.21)：2k 次立方不變多項式的個數，三種算法：
    Σ_α f^α(f^α+1)/2（f = 𝒟^k↓O）、Σ_{ℓ偶 ≤ 2k} f^{A1}(ℓ)（扣除 ℓ 奇的贗不變量）、把 k 分成至多 3 份的分割數。"""
    f = o3_to_o(k)
    a = sum(x * (x + 1) // 2 for x in f)
    b = sum(o3_to_o(l)[0] for l in range(0, 2 * k + 1, 2))
    c = sum(1 for p in itertools.combinations_with_replacement(range(k + 1), 3) if sum(p) == k)
    return a, b, c


# d 軌域（(5.6.7)–(5.6.15)）
COSET_EULER = {"1": (0, 0, 0), "R2^2": (0, 180, 0), "R1": (-90, 90, 90), "R1^3": (-90, -90, 90), "R2": (0, 90, 0), "R2^3": (0, -90, 0)}
D2_COL_BOOK = {"1": [1, 0, 0, 0, 0], "R2^2": [0, 0, 0, 0, 1], "R1": [R_(1, 4), -IU / 2, -sp.sqrt(R_(3, 8)), IU / 2, R_(1, 4)],
               "R1^3": [R_(1, 4), IU / 2, -sp.sqrt(R_(3, 8)), -IU / 2, R_(1, 4)], "R2": [R_(1, 4), R_(1, 2), sp.sqrt(R_(3, 8)), R_(1, 2), R_(1, 4)],
               "R2^3": [R_(1, 4), -R_(1, 2), sp.sqrt(R_(3, 8)), -R_(1, 2), R_(1, 4)]}            # (5.6.11)
T2_SUM_BOOK3 = [1, 0, 0, 0, 1]                       # (5.6.12) 前第三個和照印（最後一個應 −1）


def _t2_book(g) -> np.ndarray:
    """書上 §5.6B 的 T₂：𝒟^{T2}(g) = χ^{A2}(g)𝒟^{T1}(g)（不是第 4 章的 kosher 版）。"""
    k = next(i for i, C in enumerate(SO.O.classes) if g in C)
    return float(SO.O_CHARS[1, k]) * _num(SO.irrep("T_1", g)).real


def d_orbital_states() -> dict:
    """(5.6.9)–(5.6.13)：只對六個陪集代表求和（係數 E：1/3，T₂：1/2），再歸一化。→ {"E1","E2","T2_1","T2_2","T2_3": 5 維向量（m = 2…−2）}
    ＋ "sums"（未歸一化的三個 T₂ 和）。"""
    cols = {g: D_expm(2, *COSET_EULER[g])[:, 0] for g in COSET_EULER}
    out, sums = {}, {}
    for t in range(3):
        v = sum(_t2_book(g)[t, 2] * cols[g] for g in COSET_EULER)
        sums[t + 1] = v
        out[f"T2_{t + 1}"] = v / np.linalg.norm(v)
    for e in range(2):
        v = sum(float(_num(SO.irrep("E", g))[e, 1].real) * cols[g] for g in COSET_EULER)
        out[f"E{e + 1}"] = v / np.linalg.norm(v)
    out["sums"] = sums
    return out


D_ORBITALS_BOOK = {"T2_1": [0, -1, 0, -1, 0], "T2_2": [0, -1j, 0, 1j, 0], "T2_3": [1, 0, 0, 0, -1], "E1": [0, 0, 1, 0, 0], "E2": [1, 0, 0, 0, 1]}  # ×1/√2 等
D_ORBITAL_POLY = {"T2_1": "yz", "T2_2": "xz", "T2_3": "xy", "E1": "2z^2-x^2-y^2", "E2": "x^2-y^2"}


def d_orbital_projector_full(label: str) -> np.ndarray:
    """第二種算法：對 24 個元素求和 P = (ℓ/24)Σ 𝒟^α_{tt}(g)* 𝒟²(g)，作用在 |2 2⟩。"""
    t = {"T2_1": 0, "T2_2": 1, "T2_3": 2, "E1": 0, "E2": 1}[label]
    v = np.zeros(5, complex)
    for g in SO.O_NAMES:
        ax, w = SO.O_TURNS[g]
        D2 = D_axis(2, ax, w)
        if label.startswith("T2"):
            c = _t2_book(g)[t, 2]
            v += 3 / 24 * c * D2[:, 0]
        else:
            c = _num(SO.irrep("E", g))[t, 1].real
            v += 2 / 24 * c * D2[:, 0]
    return v


# 多極展開（(5.6.17)–(5.6.29)）
_x, _y, _z = sp.symbols("x y z", real=True)
_r2 = _x ** 2 + _y ** 2 + _z ** 2
X4 = {0: (35 * _z ** 4 - 30 * _z ** 2 * _r2 + 3 * _r2 ** 2) / 8, 4: sp.sqrt(R_(35, 128)) * (_x + IU * _y) ** 4,
      -4: sp.sqrt(R_(35, 128)) * (_x - IU * _y) ** 4}                     # X^4_m = √(4π/9) r⁴ Y^4_m（附錄 F）
V4 = _x ** 4 + _y ** 4 + _z ** 4 - R_(3, 5) * _r2 ** 2                    # (5.6.16)


def v4_identities() -> tuple:
    """(5.6.28) 方括號的兩種寫法、(5.6.29) 的球座標寫法（sympy 化簡為 0）。"""
    a = sp.expand((2 / sp.sqrt(70)) * (X4[4] + X4[-4]) + R_(2, 5) * X4[0] - V4)
    th, ph, rr = sp.symbols("theta phi r", positive=True)
    sub = {_x: rr * sp.sin(th) * sp.cos(ph), _y: rr * sp.sin(th) * sp.sin(ph), _z: rr * sp.cos(th)}
    b = sp.simplify(sp.expand_trig(V4.subs(sub) - rr ** 4 * (35 * sp.cos(th) ** 4 - 30 * sp.cos(th) ** 2 + 3
                                                             + 5 * sp.sin(th) ** 4 * sp.cos(4 * ph)) / 20))
    return sp.simplify(a), b


def oct_potential(r, Q: float = 1.0, R: float = 1.0) -> float:
    """(5.6.23)：Σ Q/|R r̂_j − r|（六個電荷在 ±x, ±y, ±z）。"""
    r = np.asarray(r, float)
    return float(sum(Q / np.linalg.norm(R * np.array(s, float) - r) for s in SO.SITES))


def oct_multipole(r, Q: float = 1.0, R: float = 1.0) -> float:
    """(5.6.28)：6Q/R + (35Q/4R⁵)(x⁴ + y⁴ + z⁴ − 3r⁴/5)。"""
    x, y, z = r
    rr = x * x + y * y + z * z
    return 6 * Q / R + 35 * Q / (4 * R ** 5) * (x ** 4 + y ** 4 + z ** 4 - 0.6 * rr * rr)


def multipole_coeffs(Q: float = 1.0, R: float = 1.0) -> dict:
    """(5.6.25b)：x^ℓ_m = Σ_j (q_j / r_j^{2ℓ+1}) X^{ℓ*}_m(r_j)（ℓ ≤ 4，用 Y^ℓ_m 數值算）。"""
    out = {}
    for l in range(5):
        for m in range(-l, l + 1):
            s = 0j
            for v in SO.SITES:
                v = np.array(v, float)
                th, ph = math.acos(v[2]), math.atan2(v[1], v[0])
                X = math.sqrt(4 * math.pi / (2 * l + 1)) * R ** l * ylm(l, m, th, ph)
                s += Q / R ** (2 * l + 1) * np.conj(X)
            if abs(s) > 1e-12:
                out[(l, m)] = complex(s)
    return out


# SF6 型離心扭曲（(5.6.30)–(5.6.36)、Fig. 5.6.7）
@lru_cache(maxsize=None)
def _o_mats(twoJ: int) -> dict:
    """24 個 O 元素的 𝒟^J（e^{−iω n̂·J}；半整數 J 用 SO.O_TURNS 的轉角當 SU(2) 代表）。"""
    return {g: D_axis(Fraction(twoJ, 2), *SO.O_TURNS[g]) for g in SO.O_NAMES}


def _o_label_of_space(P, J) -> str:
    """把一個本徵子空間（投影算子 P）依 O 的特徵標分解成標記（"T_2"、"A_1+E"；半整數 J 用射線特徵標 E_1、E_2、G）。"""
    if abs(float(J) - round(float(J))) < 1e-9:
        reps = ["1", "r1", "R1^2", "R1", "i1"]
        mats = _o_mats(int(round(2 * float(J))))
        tr = [np.trace(P @ mats[g]).real for g in reps]
        labs = []
        for a in range(5):
            f = sum(float(SO.O_CHARS[a, k]) * SO.O_CLASS_SIZES[k] * tr[k] for k in range(5)) / 24
            labs += [O_IRREPS[a]] * int(round(f))
        return "+".join(labs)
    full = _ray_chars_full(OS, OS_CHARS, OS_REPS)
    mats = _o_mats(int(round(2 * float(J))))
    trs = {g: np.trace(P @ mats[g]).real for g in OS.names}
    labs = []
    for a in ("E_1", "E_2", "G"):
        f = sum(float(sp.N(full[a][g])) * trs[g] for g in OS.names) / 24
        labs += [a] * int(round(f))
    return "+".join(labs)


def sf6_levels(J, t: float = 1.0, B: float = 0.0) -> list:
    """B J² + t(J_x̄⁴ + J_ȳ⁴ + J_z̄⁴ − (3/5)J⁴) 在 |J n⟩（2J+1 維）的能階 → [(能量, O 標記, 簡併度)]，由低到高
    （J 可為半整數：標記 E_1、E_2、G，每條還有 Kramers 二重）。"""
    M = jmat(J)
    n = M["Jz"].shape[0]
    jj = float(J) * (float(J) + 1)
    H4 = sum(np.linalg.matrix_power(M[k], 4) for k in ("Jx", "Jy", "Jz")) - 0.6 * jj ** 2 * np.eye(n)
    H = B * jj * np.eye(n) + t * H4
    w, V = np.linalg.eigh((H + H.conj().T) / 2)
    out, i = [], 0
    while i < len(w):
        k = i
        while k + 1 < len(w) and abs(w[k + 1] - w[i]) < 1e-6 * max(1.0, abs(w[i])):
            k += 1
        P = V[:, i:k + 1] @ V[:, i:k + 1].conj().T
        out.append((float(w[i]), _o_label_of_space(P, J), k - i + 1))
        i = k + 1
    return out


def sf6_clusters(J: int, gap_ratio: float = 0.01) -> list:
    """大 J 的能階叢（相鄰能階差 < gap_ratio × 總寬度就算同一叢）→ [[標記…], …]（由高到低）。J = 20：
    [A1, T1, E]（6 個態，≈ 0₄↑O）、[T2, T1]（1₄↑O）、[E, T2, A2]（2₄↑O）… 最低 [T1, E, T2]（8 個態，1₃↑O）。"""
    lev = sf6_levels(J)[::-1]
    width = lev[0][0] - lev[-1][0]
    out = [[lev[0]]]
    for a, b in zip(lev, lev[1:]):
        if a[0] - b[0] < gap_ratio * width:
            out[-1].append(b)
        else:
            out.append([b])
    return [[x[1] for x in c] for c in out]


# R3 ⊃ D6（(5.6.37)–(5.6.41)）
D6_CLASS_ANGLES = [0, 120, 180, 180, 60, 180]        # SN.D6_CHARS_BOOK 的行：1, h², ρ, h³, h, ρ'
D6_CLASS_SIZES = [1, 2, 3, 1, 2, 3]
D6_IRREPS = SN.D6_IRREP_NAMES                        # A_1 A_2 E_2 B_1 B_2 E_1（書上 (5.6.37) 的欄序）


def o3_to_d6(l, parity: int = 1) -> tuple:
    """(5.6.37)：𝒟^ℓ↓D6 的頻率（欄 A1 A2 E2 B1 B2 E1）。"""
    T = sp.Matrix(SN.D6_CHARS_BOOK)
    tr = [chi_exact(l, w) for w in D6_CLASS_ANGLES]
    return tuple(int(sp.nsimplify(sum(T[a, k] * D6_CLASS_SIZES[k] * tr[k] for k in range(6)) / 12)) for a in range(6))


def o3_to_d6_mrule(l) -> tuple:
    """第二種算法（(5.6.37) 的註）：m = 0 → A1（ℓ 偶）或 A2（ℓ 奇）；±m 對依 |m| mod 6：±1,±5 → E1、±2,±4 → E2、3 → B1+B2、0 → A1+A2。"""
    f = dict.fromkeys(D6_IRREPS, 0)
    f["A_1" if l % 2 == 0 else "A_2"] += 1
    for m in range(1, l + 1):
        r = m % 6
        if r in (1, 5):
            f["E_1"] += 1
        elif r in (2, 4):
            f["E_2"] += 1
        elif r == 3:
            f["B_1"] += 1
            f["B_2"] += 1
        else:
            f["A_1"] += 1
            f["A_2"] += 1
    return tuple(f[n] for n in D6_IRREPS)


O3_D6_BOOK = {0: (1, 0, 0, 0, 0, 0), 1: (0, 1, 0, 0, 0, 1), 2: (1, 0, 1, 0, 0, 1), 3: (0, 1, 1, 1, 1, 1), 4: (1, 0, 2, 1, 1, 1),
              5: (0, 1, 2, 1, 1, 2), 6: (2, 1, 2, 1, 1, 2), 7: (1, 2, 2, 1, 1, 3), 8: (2, 1, 3, 1, 1, 3), 9: (1, 2, 3, 2, 2, 3),
              10: (2, 1, 4, 2, 2, 3)}               # (5.6.37) 照印（全對）


def d180(N) -> np.ndarray:
    """(5.6.39)–(5.6.40)：𝒟^N(0 180° 0) = δ_{n,−m}(−1)^{N+m}。"""
    return wigner_d(N, 180)


# 助記輪（Fig. 5.6.8、5.6.9、5.7.3；讀法：從起點沿方向讀到維度和 = 2J+1）
WHEEL_D6_EVEN = ["A_1", "E_1", "E_2", "B_1", "B_2", "E_2", "E_1", "A_2"]
WHEEL_D6_ODD = ["A_2", "E_1", "E_2", "B_2", "B_1", "E_2", "E_1", "A_1"]
WHEEL_O = ["A_1", "T_1", "E", "T_2", "T_1", "E", "T_2", "A_2", "T_2", "T_1"]   # Fig. 5.6.9 外圈（A1 在左）
WHEEL_O_START = {0: (0, 1), 2: (5, 1), 1: (4, -1), 3: (9, -1)}               # J mod 4 → (起點, 方向)
DIMS = {"A_1": 1, "A_2": 1, "B_1": 1, "B_2": 1, "E": 2, "E_1": 2, "E_2": 2, "E_3": 2, "T_1": 3, "T_2": 3, "G": 4}


def wheel_read(slots, start: int, step: int, dim: int) -> list:
    """從 slots[start] 起，每次走 step 格，讀到維度和 = dim 為止（不剛好就 AssertionError）。"""
    out, tot, k = [], 0, start
    while tot < dim:
        out.append(slots[k % len(slots)])
        tot += DIMS[slots[k % len(slots)]]
        k += step
    assert tot == dim
    return out


def wheel_d6(J: int) -> list:
    return wheel_read(WHEEL_D6_EVEN if J % 2 == 0 else WHEEL_D6_ODD, 0, 1, 2 * J + 1)


def wheel_o(J: int) -> list:
    s, d = WHEEL_O_START[J % 4]
    return wheel_read(WHEEL_O, s, d, 2 * J + 1)


def _count(labels, names) -> tuple:
    return tuple(labels.count(n) for n in names)


# ======================================================================== §5.7 射線表示（H513–H514）
class SpinAlgebra:
    """有限轉動群的自旋 ½「射線代數」：元素用 SU(2) 矩陣代表（名稱 → 2×2），乘積 ab = ±c（先 b 後 a）。
    .mul(a, b) → (號, 名)、.inv(a) → (號, 名)（代數的逆：a·a⁻¹ = 1）、.table(rows, cols) → TeX 字串表、
    .conj_sum(h) → Σ_g g h g⁻¹（dict 名 → 係數）、.class_ops() → 線性獨立的全對易算子（每個是 dict）、
    .ga_mul(x, y)（群代數乘法，x、y 是 dict）、.chars(Pdict) → 由冪等元讀特徵標。"""

    def __init__(self, mats: dict, tex: dict | None = None, label: str = ""):
        self.names = list(mats)
        self.U = {g: np.asarray(M, complex) for g, M in mats.items()}
        self.tex = tex or {g: g for g in self.names}
        self.label = label
        self.order = len(self.names)
        self._keys = {}
        for g in self.names:
            self._keys[self._key(self.U[g])] = (1, g)
            self._keys[self._key(-self.U[g])] = (-1, g)
        self._mul = {}

    @staticmethod
    def _key(M) -> tuple:
        M = np.asarray(M, complex)
        return tuple(np.round(np.concatenate([M.real.ravel(), M.imag.ravel()]), 6) + 0.0)

    def find(self, M) -> tuple:
        hit = self._keys.get(self._key(M))
        if hit is not None and _close(M, hit[0] * self.U[hit[1]], 1e-8):
            return hit
        for g in self.names:
            if _close(M, self.U[g], 1e-8):
                return 1, g
            if _close(M, -self.U[g], 1e-8):
                return -1, g
        raise ValueError("不在代數裡")

    def mul(self, a, b) -> tuple:
        if (a, b) not in self._mul:
            self._mul[(a, b)] = self.find(self.U[a] @ self.U[b])
        return self._mul[(a, b)]

    def inv(self, a) -> tuple:
        return self.find(np.linalg.inv(self.U[a]))

    def signed_tex(self, sg, g) -> str:
        return ("-" if sg < 0 else "") + self.tex[g]

    def table(self, rows, cols) -> list:
        return [[self.signed_tex(*self.mul(r, c)) for c in cols] for r in rows]

    def ga_mul(self, x: dict, y: dict) -> dict:
        out = {}
        for a, ca in x.items():
            for b, cb in y.items():
                s, g = self.mul(a, b)
                out[g] = out.get(g, 0) + s * ca * cb
        return {g: sp.nsimplify(c) for g, c in out.items() if sp.nsimplify(c) != 0}

    def conj_sum(self, h) -> dict:
        out = {}
        for g in self.names:
            si, gi = self.inv(g)
            s1, x = self.mul(g, h)
            s2, y = self.mul(x, gi)
            out[y] = out.get(y, 0) + si * s1 * s2
        return {g: c for g, c in out.items() if c != 0}

    def class_ops(self) -> list:
        """全對易算子：每個 h 的 Σ g h g⁻¹ 除以 °N_h（＝ 每個不同的共軛元素出現一次），去掉 0 與重複。"""
        ops, seen = [], []
        for h in self.names:
            cs = self.conj_sum(h)
            if not cs:
                continue
            nh = sum(1 for g in self.names if _close(self.U[g] @ self.U[h], self.U[h] @ self.U[g], 1e-8))
            c = {g: sp.Rational(v, nh) for g, v in cs.items()}
            key = frozenset(c)
            if key in seen:
                continue
            seen.append(key)
            ops.append((h, c))
        return ops

    def trace_char(self, j) -> dict:
        """χ^j(g) ＝ 𝒟^j(g) 的跡（g 的 SU(2) 代表 → 轉角 ω ∈ [0, 720)）。"""
        out = {}
        for g in self.names:
            c = max(-1.0, min(1.0, float(np.trace(self.U[g]).real) / 2))
            w = math.degrees(2 * math.acos(c))
            out[g] = chi(j, w)
        return out


def _spin(axis, deg) -> np.ndarray:
    return SN.su2(SN.quat(axis, deg))


def _hor(az) -> tuple:
    return (math.cos(math.radians(az)), math.sin(math.radians(az)), 0.0)


D6S_AXES = {"p1": 120, "p2": 60, "p3": 0, "p1'": 210, "p2'": 150, "p3'": 90}       # 180° 橫軸的方位角（定號）
D6S = SpinAlgebra({**{"1": np.eye(2)}, **{("h" if k == 1 else f"h^{k}"): _spin((0, 0, 1), 60 * k) for k in range(1, 6)},
                   **{p: _spin(_hor(az), 180) for p, az in D6S_AXES.items()}},
                  {"1": r"\mathbf 1", "h": "h", **{f"h^{k}": f"h^{k}" for k in range(2, 6)},
                   "p1": r"\rho_1", "p2": r"\rho_2", "p3": r"\rho_3", "p1'": r"\rho_1'", "p2'": r"\rho_2'", "p3'": r"\rho_3'"}, "D6 spin")
D6S_COLS = ["1", "h^2", "h^4", "p1", "p2", "p3", "h^3", "h", "h^5", "p1'", "p2'", "p3'"]
D6S_ROWS = ["1", "h^4", "h^2", "p1", "p2", "p3", "h^3", "h^5", "h", "p1'", "p2'", "p3'"]
_T574 = """h^4: -1 -h^2 -p2 -p3 p1 -h h^5 -h^3 -p2' -p3' p1'
h^2: h^4 -1 -p3 p1 p2 h^5 h^3 -h -p3' p1' p2'
p1: p2 p3 -1 -h^2 -h^4 -p1' p3' -p2' h^3 h^5 -h
p2: p3 -p1 h^4 -1 -h^2 -p2' -p1' -p3' h h^3 h^5
p3: -p1 -p2 h^2 h^4 -1 -p3' -p2' p1' -h^5 h h^3
h^3: h^5 -h p1' p2' p3' -1 h^4 -h^2 -p1 -p2 -p3
h^5: -h -h^3 -p3' p1' p2' -h^2 -1 -h^4 p3 -p1 -p2
h: h^3 h^5 p2' p3' -p1' h^4 h^2 -1 -p2 -p3 p1
p1': p2' p3' -h^3 -h^5 h p1 -p3 p2 -1 -h^2 -h^4
p2': p3' -p1' -h -h^3 -h^5 p2 p1 p3 h^4 -1 -h^2
p3': -p1' -p2' h^5 -h -h^3 p3 p2 -p1 h^2 h^4 -1"""
D6S_TABLE_BOOK = {ln.split(":")[0]: ln.split(":")[1].split() for ln in _T574.splitlines()}   # (5.7.4)：列 × 行（行 h², h⁴, …, ρ3'）


def d6s_table_check() -> int:
    """(5.7.4) 的 121 個格子逐一和 SU(2) 乘積比對 → 相符的個數。"""
    ok = 0
    for r, entries in D6S_TABLE_BOOK.items():
        for c, e in zip(D6S_COLS[1:], entries):
            s, g = D6S.mul(r, c)
            want = (-1 if e.startswith("-") else 1, e.lstrip("-"))
            ok += (s, g) == want
    return ok


def d6s_class_algebra() -> dict:
    """(5.7.8)–(5.7.13)：c2 = h² − h⁴、c3 = h − h⁵；c2² = 2 − c2、c3² = 2 + c2、c2c3 = c3；c3³ = 3c3；冪等元與特徵標。"""
    one = {"1": 1}
    c2, c3 = {"h^2": 1, "h^4": -1}, {"h": 1, "h^5": -1}
    m = D6S.ga_mul
    P1 = {"1": R_(1, 3), "h^2": R_(1, 6), "h^4": -R_(1, 6), "h": SQ3 / 6, "h^5": -SQ3 / 6}
    P2 = {"1": R_(1, 3), "h^2": R_(1, 6), "h^4": -R_(1, 6), "h": -SQ3 / 6, "h^5": SQ3 / 6}
    P3 = {"1": R_(1, 3), "h^2": -R_(1, 3), "h^4": R_(1, 3)}
    return dict(one=one, c2=c2, c3=c3, c2c2=m(c2, c2), c3c3=m(c3, c3), c2c3=m(c2, c3), c3cube=m(m(c3, c3), c3),
                P=[P1, P2, P3], PP=[[m(a, b) for b in (P1, P2, P3)] for a in (P1, P2, P3)])


D6S_CHARS = {"E_1": (2, 1, SQ3), "E_2": (2, 1, -SQ3), "E_3": (2, -2, 0)}     # (5.7.13)：類 1, h², h
D6S_REPS = ["1", "h^2", "h"]


_FULL_CACHE = {}


def _ray_chars_full(alg: SpinAlgebra, table: dict, reps: list) -> dict:
    """把類代表的特徵標延伸到所有元素（用 ±共軛：χ(g h g⁻¹) 依號；不在任何非零類的元素 χ = 0）。"""
    key = (id(alg), id(table))
    if key in _FULL_CACHE:
        return _FULL_CACHE[key]
    out = {a: dict.fromkeys(alg.names, 0) for a in table}
    _FULL_CACHE[key] = out
    for k, h in enumerate(reps):
        for g in alg.names:
            si, gi = alg.inv(g)
            s1, x = alg.mul(g, h)
            s2, y = alg.mul(x, gi)
            for a in table:
                out[a][y] = si * s1 * s2 * table[a][k]
    return out


def ray_freq(alg: SpinAlgebra, chars_full: dict, j) -> dict:
    """f^α(𝒟^j↓G) = (1/°G)Σ_g χ^α(g)* χ^j(g)（g 的 SU(2) 代表；乘積與號無關）。"""
    tr = alg.trace_char(j)
    return {a: int(round(sum(float(sp.N(chars_full[a][g])) * tr[g] for g in alg.names) / alg.order)) for a in chars_full}


D6S_FREQ_BOOK = {"1/2": (0, 0, 0), "3/2": (1, 0, 1), "5/2": (1, 1, 1), "7/2": (2, 1, 1), "9/2": (2, 1, 2), "11/2": (2, 2, 2)}   # (5.7.14) 照印
D6S_TRACE60_BOOK = {"7/2": SQ3, "9/2": SQ3}           # (5.7.14) 照印（應 −√3）

# D2 的射線代數＝四元數（(5.7.17)–(5.7.20)）
D2S = SpinAlgebra({"1": np.eye(2), "Rx": _spin((1, 0, 0), 180), "Ry": _spin((0, 1, 0), 180), "Rz": _spin((0, 0, 1), 180)},
                  {"1": r"\mathbf 1", "Rx": "R_x", "Ry": "R_y", "Rz": "R_z"}, "D2 spin")
D2S_TABLE_BOOK = [["1", "Rx", "Ry", "Rz"], ["Rx", "-1", "Rz", "-Ry"], ["Ry", "-Rz", "-1", "Rx"], ["Rz", "Ry", "-Rx", "-1"]]   # (5.7.17)
D2S_DE_BOOK = {"Rx": [[0, -1j], [-1j, 0]], "Ry": [[0, -1], [1, 0]], "Rz": [[-1j, 0], [0, 1j]]}                                # (5.7.19)
Q_UNITS = {k: 1j * P for k, P in zip("xyz", PAULI)}                                                                           # (5.5.3) q = iσ


def cn_spin_chars(n: int) -> list:
    """Fig. 5.7.5：C_n 的自旋表示 χ(h) = −1 的 n 次方根 e^{iπ(2k+1)/n}（h = 繞 z 轉 360°/n 的 SU(2) 代表，hⁿ = −1）。"""
    return [complex(np.exp(1j * math.pi * (2 * k + 1) / n)) for k in range(n)]


def dn_spin_algebra(n: int) -> SpinAlgebra:
    mats = {"1": np.eye(2)}
    for k in range(1, n):
        mats[f"h^{k}"] = _spin((0, 0, 1), 360 * k / n)
    for k in range(n):
        mats[f"p{k}"] = _spin(_hor(180 * k / n), 180)
    return SpinAlgebra(mats, label=f"D{n} spin")


def dn_ray_irreps(n: int) -> tuple:
    """D_n 射線不可約表示的維度：由類代數的維數（非零的全對易算子個數）與 Σℓ² = 2n 定出；書上：n 偶 n/2 個二維，
    n 奇 (n−1)/2 個二維 ＋ 一對一維。→ (類算子個數, 維度列表)。"""
    alg = dn_spin_algebra(n)
    k = len(alg.class_ops())
    dims = [2] * (n // 2) if n % 2 == 0 else [2] * ((n - 1) // 2) + [1, 1]
    return k, dims


# 八面體的射線代數（(5.7.21)–(5.7.26)）
OS = SpinAlgebra({g: _spin(*SO.O_TURNS[g]) for g in SO.O_NAMES},
                 {g: SO.O.tex[g] for g in SO.O_NAMES}, "O spin")
OS_TABLE_522_BOOK = {("r1", "-r1^2"): "1", ("r1", "-R2^3"): "R1", ("r2", "-r2^2"): "1", ("r3", "-r3^2"): "1", ("r3", "R2"): "R1",
                     ("r4", "-r4^2"): "1", ("-r1^2", "r1"): "1", ("-r1^2", "-r1^2"): "-r1", ("-r2^2", "r2"): "1",
                     ("-r2^2", "-r3^2"): "r1", ("-r2^2", "R3"): "R1", ("-r3^2", "r3"): "1", ("-r3^2", "-r4^2"): "r1",
                     ("-r4^2", "r4"): "1", ("-r4^2", "-r2^2"): "r1", ("-r4^2", "-R3^3"): "R1", ("R1", "R2"): "r1",
                     ("R1", "-R1^3"): "1", ("R2", "R3"): "r1", ("R2", "-R2^3"): "1", ("R3", "R1"): "r1", ("R3", "-R3^3"): "1",
                     ("-R1^3", "R1"): "1", ("-R2^3", "R2"): "1", ("-R3^3", "R3"): "1"}   # (5.7.22) 的非空格（列 × 行）


def _signed(alg, s: str) -> np.ndarray:
    return (-1 if s.startswith("-") else 1) * alg.U[s.lstrip("-")]


def os_table_check() -> int:
    ok = 0
    for (a, b), c in OS_TABLE_522_BOOK.items():
        ok += _close(_signed(OS, a) @ _signed(OS, b), _signed(OS, c))
    return ok


def os_class_algebra() -> dict:
    """(5.7.21)、(5.7.23)：c_r = Σr_k − Σr_k²、c_R = ΣR_k − ΣR_k³；c_r² = 8 + 2c_r、c_rc_R = 4c_R、c_R² = 6 + 3c_r、c_R³ = 18c_R。"""
    cr = {**{f"r{k}": 1 for k in range(1, 5)}, **{f"r{k}^2": -1 for k in range(1, 5)}}
    cR = {**{f"R{k}": 1 for k in range(1, 4)}, **{f"R{k}^3": -1 for k in range(1, 4)}}
    m = OS.ga_mul
    return dict(cr=cr, cR=cR, crcr=m(cr, cr), crcR=m(cr, cR), cRcR=m(cR, cR), cRcube=m(m(cR, cR), cR))


OS_CHARS = {"E_1": (2, 1, SQ2), "E_2": (2, 1, -SQ2), "G": (4, -1, 0)}        # (5.7.25)：類 1, r1, R1
OS_REPS = ["1", "r1", "R1"]
OS_FREQ_BOOK = {"1/2": (1, 0, 0), "3/2": (0, 0, 1), "5/2": (0, 1, 0), "7/2": (1, 1, 1), "9/2": (1, 0, 2), "11/2": (1, 1, 2)}   # (5.7.26b) 照印


def _jlab(j) -> str:
    f = Fraction(j).limit_denominator(2)
    return f"{f.numerator}/{f.denominator}" if f.denominator != 1 else str(f.numerator)


def half_correlation(group: str, j) -> tuple:
    """𝒟^j↓D6（E1, E2, E3）或 ↓O（E1, E2, G）的頻率（半整數 j，SU(2) 跡）。"""
    if group == "D6":
        fr = ray_freq(D6S, _ray_chars_full(D6S, D6S_CHARS, D6S_REPS), j)
        return tuple(fr[a] for a in ("E_1", "E_2", "E_3"))
    fr = ray_freq(OS, _ray_chars_full(OS, OS_CHARS, OS_REPS), j)
    return tuple(fr[a] for a in ("E_1", "E_2", "G"))


def half_d6_mrule(j) -> tuple:
    """半整數 j 的 D6：±m 對依 m mod 6：±½ → E1、±3/2 → E3、±5/2 → E2（Fig. 5.7.3 的輪；圖上 E2、E3 的名字對調了）。"""
    f = {"E_1": 0, "E_2": 0, "E_3": 0}
    for k in range(int(round(float(j) + 0.5))):
        m = Fraction(2 * k + 1, 2) % 6
        f[{Fraction(1, 2): "E_1", Fraction(11, 2): "E_1", Fraction(3, 2): "E_3", Fraction(9, 2): "E_3",
           Fraction(5, 2): "E_2", Fraction(7, 2): "E_2"}[m]] += 1
    return f["E_1"], f["E_2"], f["E_3"]


WHEEL_D6_HALF = ["E_1", "E_3", "E_2", "E_2", "E_3", "E_1"]   # 1₁₂, 3₁₂, 5₁₂, 7₁₂, 9₁₂, 11₁₂（依 (5.7.13) 的名字；圖上 E2、E3 對調）
# Fig. 5.7.6 的八面體半整數輪（本系列重建：八格、0 號在左、順時針經頂部；E1／E2／G 的位置照圖說——(a) E1 左右、E2 頂上兩格、G 其餘四格；
# (b) E1、E2 對調）。讀法同整數輪：從起點順時針讀到維度和 = 2j+1；讀出的次序 ＝ V⁽⁴⁾ 能階由高到低（j ≤ 21½ 全部核對）。
WHEEL_O_HALF = {"a": ["E_1", "G", "E_2", "E_2", "G", "E_1", "G", "G"], "b": ["E_2", "G", "E_1", "E_1", "G", "E_2", "G", "G"]}
WHEEL_O_HALF_START = {"1/2": ("a", 5), "5/2": ("a", 1), "3/2": ("b", 7), "7/2": ("b", 3)}    # 2j mod 8 → (輪, 起點)


def wheel_d6_half(j) -> list:
    return wheel_read(WHEEL_D6_HALF, 0, 1, int(round(2 * float(j))) + 1)


def wheel_o_half(j) -> list:
    key = {1: "1/2", 5: "5/2", 3: "3/2", 7: "7/2"}[int(round(2 * float(j))) % 8]
    w, s = WHEEL_O_HALF_START[key]
    return wheel_read(WHEEL_O_HALF[w], s, 1, int(round(2 * float(j))) + 1)


# ======================================================================== §5.8 R4 與 U3（H515–H516）
def r4_generators() -> dict:
    """(5.8.15)：4×4 的 L_x, L_y, L_z（前三維）與 K_x, K_y, K_z（第四維）。"""
    def L(a, b):
        M = sp.zeros(4)
        M[a, b], M[b, a] = -IU, IU
        return M
    return dict(Lx=L(1, 2), Ly=L(2, 0), Lz=L(0, 1), Kx=L(0, 3), Ky=L(1, 3), Kz=L(2, 3))


R4_BASIS = ["Lx", "Kx", "Ly", "Ky", "Lz", "Kz"]


def r4_regular(name: str) -> sp.Matrix:
    """(5.8.17b)、(5.8.18)：Lie 代數的正則表示 [X, e_b] = Σ_a 𝓡_ab e_a（基底 L_x, K_x, L_y, K_y, L_z, K_z）。"""
    G = r4_generators()
    A = sp.Matrix([[G[b][i, k] for b in R4_BASIS] for i in range(4) for k in range(4)])
    R = sp.zeros(6)
    for col, b in enumerate(R4_BASIS):
        C = comm(G[name], G[b])
        rhs = sp.Matrix([C[i, k] for i in range(4) for k in range(4)])
        sol = (A.T * A).inv() * A.T * rhs
        R[:, col] = sol
    return sp.simplify(R)


R4_REG_LZ_BOOK = sp.Matrix([[0, 0, -IU, 0, 0, 0], [0, 0, 0, -IU, 0, 0], [IU, 0, 0, 0, 0, 0], [0, IU, 0, 0, 0, 0], [0] * 6, [0] * 6])
R4_REG_KZ_BOOK = sp.Matrix([[0, 0, 0, -IU, 0, 0], [0, 0, -IU, 0, 0, 0], [0, IU, 0, 0, 0, 0], [IU, 0, 0, 0, 0, 0], [0] * 6, [0] * 6])


def r4_ladders() -> dict:
    """(5.8.20a) 更正：M± = L± + K±、N± = L± − K±（L± = L_x ± iL_y）。"""
    G = r4_generators()
    Lp, Lm = G["Lx"] + IU * G["Ly"], G["Lx"] - IU * G["Ly"]
    Kp, Km = G["Kx"] + IU * G["Ky"], G["Kx"] - IU * G["Ky"]
    return dict(Mp=Lp + Kp, Np=Lp - Kp, Nm=Lm - Km, Mm=Lm + Km)


def coulomb_shell(n: int) -> dict:
    """n 殼層（n² 維）：M = J^{(j)} ⊗ 1、N = 1 ⊗ J^{(j)}，j = (n−1)/2；L = M + N、K = M − N（(5.8.21)）。"""
    j = Fraction(n - 1, 2)
    J = jmat(j)
    one = np.eye(n)
    M = {a: np.kron(J[a], one) for a in ("Jx", "Jy", "Jz")}
    N = {a: np.kron(one, J[a]) for a in ("Jx", "Jy", "Jz")}
    L = {a: M[a] + N[a] for a in M}
    K = {a: M[a] - N[a] for a in M}
    M2 = sum(M[a] @ M[a] for a in M)
    N2 = sum(N[a] @ N[a] for a in N)
    return dict(j=j, M=M, N=N, L=L, K=K, M2=M2, N2=N2, dim=n * n)


def coulomb_energy(n: int, k: float = -1.0, m: float = 1.0) -> float:
    """(5.8.33)–(5.8.34)：−(k²m/2)/(2M² + 2N² + 1) ＝ −(k²m/2)/n²（本系列取 k = −1、m = 1：−1/(2n²)）。"""
    return -(k * k * m / 2) / n ** 2


def coulomb_weights(jm, jn) -> list:
    """Fig. 5.8.3：𝒟^{j_m j_n} 的權 (L_z, K_z) = (m + n, m − n)。"""
    return [(float(a + b), float(a - b)) for a in _ms(jm) for b in _ms(jn)]


def kepler_check(r0=(1.0, 0.0, 0.0), v0=(0.0, 1.2, 0.3), k: float = -1.0, m: float = 1.0) -> dict:
    """古典 Kepler 軌道（數值積分）：ε = r/r − L×p/km 守恆、ε² = 1 + (2H/k²m)L²（(5.8.3)–(5.8.5)）。"""
    def f(t, y):
        r = y[:3]
        return np.concatenate([y[3:], k * r / (m * np.linalg.norm(r) ** 3)])          # F = −∇(k/r) = k r/r³
    sol = solve_ivp(f, (0, 20), np.concatenate([r0, v0]), rtol=1e-11, atol=1e-12, t_eval=np.linspace(0, 20, 41))
    eps, e2 = [], []
    for y in sol.y.T:
        r, p = y[:3], m * y[3:]
        L = np.cross(r, p)
        H = p @ p / (2 * m) + k / np.linalg.norm(r)
        ev = r / np.linalg.norm(r) - np.cross(L, p) / (k * m)
        eps.append(ev)
        e2.append((ev @ ev, 1 + 2 * H / (k * k * m) * (L @ L)))
    eps = np.array(eps)
    return dict(spread=float(np.max(np.abs(eps - eps[0]))), e2=e2, eps=eps[0])


def runge_lenz_poisson() -> dict:
    """古典 Poisson 括號（sympy）：{ε₁, ε₂} = −(2H/k²m)L₃ ⇒ 量子 [ε₁, ε₂] = iħ{ , } = −i(2H/k²m)L₃（(5.8.10) 第二行）；
    {L_i, ε_j} = ε_ijk ε_k（(5.8.7)）。"""
    x = sp.symbols("x1:4", real=True)
    p = sp.symbols("p1:4", real=True)
    k, m = sp.symbols("k m", nonzero=True)
    r = sp.sqrt(sum(a * a for a in x))
    L = [x[1] * p[2] - x[2] * p[1], x[2] * p[0] - x[0] * p[2], x[0] * p[1] - x[1] * p[0]]
    Lxp = [L[1] * p[2] - L[2] * p[1], L[2] * p[0] - L[0] * p[2], L[0] * p[1] - L[1] * p[0]]
    eps = [x[i] / r - Lxp[i] / (k * m) for i in range(3)]
    H = sum(a * a for a in p) / (2 * m) + k / r

    def pb(A, B):
        return sum(sp.diff(A, x[i]) * sp.diff(B, p[i]) - sp.diff(A, p[i]) * sp.diff(B, x[i]) for i in range(3))
    d1 = sp.simplify(pb(eps[0], eps[1]) + 2 * H / (k * k * m) * L[2])
    d2 = sp.simplify(pb(L[0], eps[1]) - eps[2])
    d3 = sp.simplify(pb(eps[0], H))
    return dict(e1e2=d1, L1e2=d2, conserved=d3)


def quantum_xr_p(point=(0.3, -0.4, 0.5)) -> tuple:
    """(5.8.8) 更正：[x_i/r, p_j]ψ = i(δ_ij − x_i x_j/r²)ψ/r（p = −i∂；用試驗函數 ψ 在一點數值核對 i = 1、j = 1, 2）。
    → (更正式的最大誤差, 照印 i(δ_ij + x_i x_j)/r 的最大誤差)。"""
    X = sp.symbols("X1:4", real=True)
    r = sp.sqrt(sum(a * a for a in X))
    psi = sp.exp(-(X[0] ** 2 + 2 * X[1] ** 2 + 3 * X[2] ** 2) / 2) * (1 + X[0] + X[1] * X[2])
    sub = dict(zip(X, point))
    err_ok, err_book = 0.0, 0.0
    for i in range(3):
        for j in range(3):
            lhs = (X[i] / r) * (-IU) * sp.diff(psi, X[j]) - (-IU) * sp.diff(X[i] * psi / r, X[j])
            ok = IU * ((1 if i == j else 0) - X[i] * X[j] / r ** 2) * psi / r
            bk = IU * ((1 if i == j else 0) + X[i] * X[j]) * psi / r
            err_ok = max(err_ok, abs(complex(sp.N((lhs - ok).subs(sub)))))
            err_book = max(err_book, abs(complex(sp.N((lhs - bk).subs(sub)))))
    return err_ok, err_book


def quantum_eps_denominator(point=(0.3, -0.4, 0.5)) -> tuple:
    """(5.8.6a)（對稱化）＝ (5.8.6b) 的分母 km 版；照印的 2km 不相等。→ (km 版誤差, 2km 版誤差)。"""
    X = sp.symbols("X1:4", real=True)
    k, m = -1.3, 0.7
    psi = sp.exp(-(X[0] ** 2 + 2 * X[1] ** 2 + 3 * X[2] ** 2) / 2) * (1 + X[0] + X[1] * X[2])
    P = [lambda f, a=a: -IU * sp.diff(f, X[a]) for a in range(3)]
    Lop = [lambda f: X[1] * P[2](f) - X[2] * P[1](f), lambda f: X[2] * P[0](f) - X[0] * P[2](f),
           lambda f: X[0] * P[1](f) - X[1] * P[0](f)]
    sym = Lop[1](P[2](psi)) - Lop[2](P[1](psi)) - P[1](Lop[2](psi)) + P[2](Lop[1](psi))      # (5.8.6a) 的分子
    comp = Lop[1](P[2](psi)) - Lop[2](P[1](psi)) - IU * P[0](psi)
    sub = dict(zip(X, point))
    a = complex(sp.N((sym / (2 * k * m)).subs(sub)))
    b1 = complex(sp.N((comp / (k * m)).subs(sub)))
    b2 = complex(sp.N((comp / (2 * k * m)).subs(sub)))
    return abs(a - b1), abs(a - b2)


# U3／SU3
def E3(i: int, j: int) -> sp.Matrix:
    """(5.8.39)：基本矩陣 𝓔_ij（1 起算）。"""
    M = sp.zeros(3)
    M[i - 1, j - 1] = 1
    return M


U3_G = {"32": E3(3, 2) - E3(2, 3), "13": E3(1, 3) - E3(3, 1), "21": E3(2, 1) - E3(1, 2)}   # G_jk = −iL_jk（(5.8.37)）
U3_K = {"32": E3(3, 2) + E3(2, 3), "13": E3(1, 3) + E3(3, 1), "21": E3(2, 1) + E3(1, 2)}   # (5.8.38)


def su3_cartan(N: int = 6) -> dict:
    """(5.8.41)–(5.8.42)：H1 = (E11 − E33)/√(2N)、H2 = (E11 − 2E22 + E33)/√(6N)、E_{±1} = E12/√6, E21/√6 …（N = 6）。
    根 r(α) = ([H1, E_α], [H2, E_α]) 的係數 → dict(H1, H2, E{α}, roots{α})。"""
    H1 = (E3(1, 1) - E3(3, 3)) / sp.sqrt(2 * N)
    H2 = (E3(1, 1) - 2 * E3(2, 2) + E3(3, 3)) / sp.sqrt(6 * N)
    Eal = {1: E3(1, 2), 2: E3(2, 3), 3: E3(1, 3), -1: E3(2, 1), -2: E3(3, 2), -3: E3(3, 1)}
    Eal = {a: M / sp.sqrt(6) for a, M in Eal.items()}
    roots = {}
    for a, M in Eal.items():
        k = next((i, j) for i in range(3) for j in range(3) if M[i, j] != 0)
        roots[a] = tuple(sp.simplify(comm(H, M)[k] / M[k]) for H in (H1, H2))
    return dict(H1=H1, H2=H2, E=Eal, roots=roots)


def su3_weights(n1: int, n2: int) -> dict:
    """𝒟^{n1 n2} 的權重數：列舉形狀 (n1+n2, n2) 的半標準 Young 表（數字 1–3），權 (a1 − a3, (a1 − 2a2 + a3)/√3)
    （＝ (2√3 H1, 2√3 H2)；H1 軸 ＝ L_z）。→ {(x, y): 重數}（x、y 取 4 位小數當鍵）。"""
    rows = [n1 + n2, n2]
    out = {}

    def fill(r, prev_row, acc):
        if r == len(rows):
            a = [sum(row.count(v) for row in acc) for v in (1, 2, 3)]
            key = (round(a[0] - a[2], 4), round((a[0] - 2 * a[1] + a[2]) / math.sqrt(3), 4))
            out[key] = out.get(key, 0) + 1
            return
        for row in itertools.combinations_with_replacement((1, 2, 3), rows[r]):
            if prev_row is not None and any(row[c] <= prev_row[c] for c in range(len(row))):
                continue
            fill(r + 1, row, acc + [row])
    fill(0, None, [])
    return out


def su3_dim(n1: int, n2: int) -> int:
    return (n1 + 1) * (n2 + 1) * (n1 + n2 + 2) // 2


def osc_shell(N: int) -> dict:
    """三維振子第 N 殼層：(n1, n2, n3) 的個數 ＝ (N+1)(N+2)/2、ℓ = N, N−2, …（Σ(2ℓ+1) 相同）、能量 N + 3/2（ħω = 1）。"""
    states = [s for s in itertools.product(range(N + 1), repeat=3) if sum(s) == N]
    ls = list(range(N, -1, -2))
    return dict(states=states, count=len(states), ls=ls, lsum=sum(2 * l + 1 for l in ls), E=N + 1.5)


def osc_fock_levels(nmax: int = 3) -> list:
    """(5.8.35a)：截斷的 Fock 空間裡 Σa†a + 3/2 的本徵值與重數（ω = 1）。"""
    a = np.diag(np.sqrt(np.arange(1, nmax + 1)), 1)
    one = np.eye(nmax + 1)
    ops = [np.kron(np.kron(a, one), one), np.kron(np.kron(one, a), one), np.kron(np.kron(one, one), a)]
    H = sum(o.T @ o for o in ops) + 1.5 * np.eye((nmax + 1) ** 3)
    return [(e, d) for e, d in _levels(H) if e <= nmax + 1.5 + 1e-9]


COULOMB_SHELLS = [2 * n * n for n in range(1, 5)]                               # 2, 8, 18, 32（含自旋 2n²）
COULOMB_CLOSURES = list(itertools.accumulate(COULOMB_SHELLS))                   # 2, 10, 28, 60
OSC_SHELLS = [(N + 1) * (N + 2) for N in range(5)]                             # 2, 6, 12, 20, 30
OSC_CLOSURES = list(itertools.accumulate(OSC_SHELLS))                          # 2, 8, 20, 40, 70
MAGIC_ATOMIC_BOOK = [2, 10, 18]                                                 # Fig. 5.8.6(a) 圈出的
MAGIC_NUCLEAR_BOOK = [2, 8, 20, 28, 50, 82, 126]                                # 正文（28 以上要自旋–軌道）


# ======================================================================== 附錄 E 與習題
def appendix_e_check(j=Fraction(5, 2)) -> bool:
    """(E.2)–(E.6)、(E.17)：[J3, J±] = ±J±、[J+, J−] = 2J3、J² = J+J− + J3² − J3 = J−J+ + J3² + J3、⟨m+1|J+|m⟩ = √(j(j+1) − m(m+1))。"""
    M = jmat(j)
    Jp, Jm, Jz, J2 = M["Jp"], M["Jm"], M["Jz"], M["J2"]
    ok = _close(comm(Jz, Jp), Jp) and _close(comm(Jz, Jm), -Jm) and _close(comm(Jp, Jm), 2 * Jz)
    ok &= _close(J2, Jp @ Jm + Jz @ Jz - Jz) and _close(J2, Jm @ Jp + Jz @ Jz + Jz)
    jf = float(j)
    for k, m in enumerate(_ms(j)[1:], start=1):
        ok &= abs(Jp[k - 1, k] - math.sqrt(jf * (jf + 1) - float(m) * (float(m) + 1))) < 1e-12
    return bool(ok)


def problem_511(n: int = 4, seed: int = 5) -> dict:
    """習題 5.1.1：正交矩陣的列也正交歸一；φ₁·Ψ₁ = (R²)₁₁、φ₁·Ψ₂ = (R²)₁₂；么正矩陣：⟨φ₁|φ₂⟩ = 0、⟨φ₁|Ψ_k⟩ = (R̄R)₁ₖ。"""
    rng = np.random.default_rng(seed)
    Q, _ = np.linalg.qr(rng.normal(size=(n, n)))
    U, _ = np.linalg.qr(rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n)))
    return dict(rows_orthonormal=bool(np.allclose(Q @ Q.T, np.eye(n))), phi1_psi1=bool(np.isclose(Q[0] @ Q[:, 0], (Q @ Q)[0, 0])),
                phi1_psi2=bool(np.isclose(Q[0] @ Q[:, 1], (Q @ Q)[0, 1])), unitary_rows=bool(np.allclose(U @ U.conj().T, np.eye(n))),
                u_phi1_psi2=bool(np.isclose(np.vdot(U[0], U[:, 1]), (U.conj() @ U)[0, 1])))


def problem_512(phi: float = 120.0) -> dict:
    """習題 5.1.2：跡 = 1 + 2cos φ（任何軸）、非正常 −(1 + 2cos φ)；非正交基底 a = e1+e2、b = e1+e3、c = e2+e3 的
    「跡」Σ⟨a|O|a⟩ = 2 Tr O + Σ_{i≠j} O_ij（不是不變量）。→ {軸: (Tr O, Tr IO, Σ⟨a|O|a⟩, 公式)}。"""
    out = {}
    B = np.array([[1, 1, 0], [1, 0, 1], [0, 1, 1]], float)          # 行 = |a⟩, |b⟩, |c⟩
    for name, ax in (("z", (0, 0, 1)), ("x", (1, 0, 0)), ("111", (1, 1, 1))):
        Om = SN.quat_to_rot(SN.quat(ax, phi))
        out[name] = (round(float(np.trace(Om)), 6), round(float(np.trace(INV3 @ Om)), 6), round(float(np.trace(B.T @ Om @ B)), 6),
                     round(float(np.trace(Om) + Om.sum()), 6))
    return out


def problem_513() -> dict:
    """習題 5.1.3：晶格基底裡的跡是整數 ⇒ 1 + 2cos ω ∈ ℤ ⇒ ω = 0, 60, 90, 120, 180°（非正常 −(1 + 2cos ω) 同理）。
    D2h、D3h、D4h、D6h、D2d、D3d 可以；D5h（72°）、D7h（51.43°）、D4d（S8 的 45°）不行。"""
    ok_ang = [w for w in range(0, 181) if abs(2 * math.cos(math.radians(w)) - round(2 * math.cos(math.radians(w)))) < 1e-9]
    need = {"D2h": [180], "D3h": [120], "D4h": [90], "D5h": [72], "D6h": [60], "D7h": [360 / 7], "D2d": [90], "D3d": [60], "D4d": [45]}
    return dict(angles=ok_ang, groups={g: all(abs(w - round(w)) < 1e-9 and round(w) in ok_ang for w in ws) for g, ws in need.items()})


def dn_chars(n: int) -> tuple:
    """D_n 的特徵標表（類：E、C_n^p（p = 1…⌊n/2⌋）、C2'、[C2'']）→ (類 [(種類, 角度或撇, 大小)], {標記: 特徵標})。
    B1 對 C2'（含 x 軸）為 +1（與第 4 章 D2 的 B1 對 R1² 為 +1 一致）。"""
    cls = [("rot", 0, 1)] + [("rot", 360 * p / n, 1 if 2 * p == n else 2) for p in range(1, n // 2 + 1)]
    cls += [("c2", "'", n // 2), ("c2", "''", n // 2)] if n % 2 == 0 else [("c2", "'", n)]
    ch = {"A_1": [1] * len(cls), "A_2": [1 if c[0] == "rot" else -1 for c in cls]}
    if n % 2 == 0:
        ch["B_1"] = [(-1) ** round(c[1] * n / 360) if c[0] == "rot" else (1 if c[1] == "'" else -1) for c in cls]
        ch["B_2"] = [(-1) ** round(c[1] * n / 360) if c[0] == "rot" else (-1 if c[1] == "'" else 1) for c in cls]
    for k in range(1, (n - 1) // 2 + 1):
        ch[f"E_{k}"] = [2 * math.cos(math.radians(k * c[1])) if c[0] == "rot" else 0 for c in cls]
    return cls, ch


def dinfh_to_dnh(n: int) -> dict:
    """習題 5.2.1：D∞h ↓ D_nh（n 偶：D_n × C_i，標記加 g/u；n 奇：D_n × σ_h，標記加 ′/″）。特徵標算法：D_n 部分用轉角與
    橫向 π 轉的特徵標；第二個陪集用 I（n 偶）或 σ_h = I·C₂(z)（n 奇，轉角加 180°）。→ {D∞h 標記: "A_{1g}+…"}。"""
    cls, ch = dn_chars(n)
    order = 2 * n
    out = {}
    for lab in DINFH_IRREPS + [r"\Phi_u"]:
        par = 1 if lab.endswith("g") else -1
        if lab.startswith(r"\Sigma"):
            pm = 1 if "+" in lab else -1
            mm, c2 = 0, pm * par
        else:
            mm, c2 = {r"\Pi": 1, r"\Delta": 2, r"\Phi": 3}[lab.split("_")[0]], 0

        def rot(w, mm=mm):
            return 1.0 if mm == 0 else 2 * math.cos(math.radians(mm * w))
        part = []
        for a, chi_a in ch.items():
            f1 = sum(c[2] * chi_a[k] * (rot(c[1]) if c[0] == "rot" else c2) for k, c in enumerate(cls)) / order
            shift = 0 if n % 2 == 0 else 180
            f2 = sum(c[2] * chi_a[k] * par * (rot(c[1] + shift) if c[0] == "rot" else c2) for k, c in enumerate(cls)) / order
            tags = ("g", "u") if n % 2 == 0 else ("'", "''")
            for f, t in ((round((f1 + f2) / 2), tags[0]), (round((f1 - f2) / 2), tags[1])):     # °D_nh = 2·°D_n
                if f:
                    part.append((f"{f}" if f > 1 else "") + (f"{a[0]}_{{{a[2:]}{t}}}" if t in ("g", "u") else f"{a}{t}"))
        out[lab] = "+".join(part)
    return out


def dinfh_to_dnh_mrule(n: int) -> dict:
    """第二種算法：Λ = |m|：m ≡ 0 (mod n) → A1+A2（m = 0：Σ⁺ → A1、Σ⁻ → A2，依 C2' 的號 ±·宇稱）；m ≡ n/2 → B1+B2；其他 → E_{min(m, n−m)}；
    n 偶加宇稱 g/u，n 奇加 σ_h 的號 宇稱·(−1)^m（′ 或 ″）。"""
    out = {}
    for lab in DINFH_IRREPS + [r"\Phi_u"]:
        par = 1 if lab.endswith("g") else -1
        if lab.startswith(r"\Sigma"):
            m, pm = 0, (1 if "+" in lab else -1)
            base = ["A_1" if pm * par > 0 else "A_2"]
        else:
            m = {r"\Pi": 1, r"\Delta": 2, r"\Phi": 3}[lab.split("_")[0]]
            r = m % n
            base = ["A_1", "A_2"] if r == 0 else (["B_1", "B_2"] if 2 * r == n else [f"E_{min(r, n - r)}"])
        if n % 2 == 0:
            t = "g" if par > 0 else "u"
            out[lab] = "+".join(f"{a[0]}_{{{a[2:]}{t}}}" for a in base)
        else:
            t = "'" if par * (-1) ** m > 0 else "''"
            out[lab] = "+".join(f"{a}{t}" for a in base)
    return out


def problem_533() -> dict:
    """習題 5.3.3：(b) 給定矩陣的 Euler 角；(c) 軸角（兩種算法：矩陣的軸角、(5.5.2) 的閉式）。"""
    s2, s3, s6 = math.sqrt(2), math.sqrt(3), math.sqrt(6)
    Rm = np.array([[1 / s2, 1 / s6, 1 / s3], [-1 / s2, 1 / s6, 1 / s3], [0, -2 / s6, 1 / s3]])
    e = euler_from_matrix(Rm)
    d = darboux_from_matrix(Rm)
    return dict(R=Rm, euler=e, darboux=d, det=float(np.linalg.det(Rm)), back=bool(np.allclose(euler_matrix(*e), Rm)),
                back2=bool(np.allclose(darboux_matrix(*d), Rm)), via552=euler_to_darboux(*e))


CIRC_T = np.array([[-1, 0, 1], [-1j, 0, -1j], [0, math.sqrt(2), 0]]) / math.sqrt(2)   # 行 = −(x+iy)/√2, z, (x−iy)/√2（(5.4.33)、(5.6.19)）


def problem_542(a=20, b=50, c=70) -> dict:
    """習題 5.4.2：𝒟¹(αβγ) ＝ T†R(αβγ)T（T 的行是圓向量 CIRC_T）；跡相等 1 + 2cos ω。"""
    D1 = wigner_D(1, a, b, c)
    R = euler_matrix(a, b, c)
    return dict(equiv=_close(CIRC_T.conj().T @ R @ CIRC_T, D1), traces=(complex(np.trace(D1)), float(np.trace(R))))


def D_from_u(j, u) -> np.ndarray:
    """習題 5.4.5：由 𝒟^{1/2}(u) = [[u11, u12], [u21, u22]] 得
    𝒟^j_{mn}(u) = Σ_k √((j+m)!(j−m)!(j+n)!(j−n)!) u11^k u12^{j+m−k} u21^{j+n−k} u22^{k−m−n} / (k!(j+m−k)!(j+n−k)!(k−m−n)!)。"""
    ms = _ms(j)
    jj = Fraction(int(round(2 * float(j))), 2)
    u11, u12, u21, u22 = u[0, 0], u[0, 1], u[1, 0], u[1, 1]
    n_ = len(ms)
    out = np.zeros((n_, n_), complex)
    for a, m in enumerate(ms):
        for b, n in enumerate(ms):
            pre = math.sqrt(_fact(int(jj + m)) * _fact(int(jj - m)) * _fact(int(jj + n)) * _fact(int(jj - n)))
            s = 0j
            for k in range(0, int(2 * jj) + 1):
                e = [k, int(jj + m) - k, int(jj + n) - k, k - int(m + n)]
                if min(e) < 0:
                    continue
                s += u11 ** e[0] * u12 ** e[1] * u21 ** e[2] * u22 ** e[3] / (_fact(e[0]) * _fact(e[1]) * _fact(e[2]) * _fact(e[3]))
            out[a, b] = pre * s
    return out


def problem_553() -> dict:
    """習題 5.5.3(b)：R1²、r1、i1、R1 的 Euler 角 (α, β, γ) 與軸角 [φ, θ, ω]（度；β = 180° 時只有 α − γ 有意義，取 γ = 0）。"""
    out = {}
    for g in ("R1^2", "r1", "i1", "R1"):
        ax, w = SO.O_TURNS[g]
        Rm = SN.quat_to_rot(SN.quat(ax, w))
        e = euler_from_matrix(Rm)
        n = np.asarray(ax, float) / np.linalg.norm(ax)
        d = (math.degrees(math.atan2(n[1], n[0])), math.degrees(math.acos(n[2])), float(w))
        out[g] = (tuple(round(x, 2) + 0.0 for x in e), tuple(round(x, 2) + 0.0 for x in d))
    return out


def problem_554() -> dict:
    """習題 5.5.4：A = 1/2I₁、B = 1/2I₂、C = 1/2I₃；J = 1、2 的能階（(5.5.54) 的矩陣；rotor_levels 數值核對）。"""
    cases = {"a": (1, 1, 2), "b": (10, 10, 1), "c": (1, 2, 3)}
    out = {}
    for k, (I1, I2, I3) in cases.items():
        A, B, C = 1 / (2 * I1), 1 / (2 * I2), 1 / (2 * I3)
        out[k] = {J: [round(float(e), 4) for e in np.linalg.eigvalsh(asym_matrix_formula(J, A, B, C))] for J in (1, 2)}
    return out


def problem_555() -> sp.Expr:
    """習題 5.5.5：B(J_x̄² + J_ȳ²) + C J_z̄² 用 (5.5.83b) 換成 B[(J_α − cos β J_γ)²/sin²β + J_β²] + C J_γ²（sympy：差為 0）。"""
    a, b, g, Ja, Jb, Jg, B, C = sp.symbols("alpha beta gamma J_alpha J_beta J_gamma B C", real=True)
    Jx = -sp.cos(g) / sp.sin(b) * Ja + sp.sin(g) * Jb + sp.cot(b) * sp.cos(g) * Jg
    Jy = sp.sin(g) / sp.sin(b) * Ja + sp.cos(g) * Jb - sp.cot(b) * sp.sin(g) * Jg
    H = B * (Jx ** 2 + Jy ** 2) + C * Jg ** 2
    target = B * ((Ja - sp.cos(b) * Jg) ** 2 / sp.sin(b) ** 2 + Jb ** 2) + C * Jg ** 2
    return sp.simplify(H - target)


def problem_561(lmax: int = 12) -> dict:
    """習題 5.6.1：O(3) ↓ D6 到 ℓ = 12（欄 A1 A2 E2 B1 B2 E1）。"""
    return {l: o3_to_d6(l) for l in range(lmax + 1)}


Y_CLASS_ANGLES = [0, 72, 144, 120, 180]


def o3_to_y(l) -> tuple:
    """習題 5.6.2(b)：𝒟^ℓ ↓ I（欄 A, T1, T3, G, H；SO.Y_CHARS）；I_h 的宇稱照 ℓ±。"""
    tr = [chi_exact(l, w) for w in Y_CLASS_ANGLES]
    return tuple(int(sp.nsimplify(sp.simplify(sum(SO.Y_CHARS[a, k] * SO.Y_CLASS_SIZES[k] * tr[k] for k in range(5)) / 60)))
                 for a in range(5))


def o3_to_y_mats(l) -> tuple:
    """第二種算法：60 個轉動的 𝒟^ℓ 跡（數值，軸角公式 χ^ℓ(ω)）配 Y 的類（轉角 → 類；72° 與 144° 分開）。"""
    sizes = SO.Y_CLASS_SIZES
    f = []
    for a in range(5):
        s = sum(complex(sp.N(SO.Y_CHARS[a, k])) * sizes[k] * chi(l, w) for k, w in enumerate(Y_CLASS_ANGLES))
        f.append(int(round((s / 60).real)))
    return tuple(f)


def o3_to_td(l, parity: int) -> tuple:
    """習題 5.6.2(a)：𝒟^{ℓ±} ↓ T_d（欄 A1 A2 E T1 T2）：+ 與 ↓O 相同；− 時 A1↔A2、T1↔T2（T_d 的 IR、Ii 類變號）。"""
    f = o3_to_o(l)
    return f if parity > 0 else (f[1], f[0], f[2], f[4], f[3])


def o3_to_td_chars(l, parity: int) -> tuple:
    """第二種算法：直接用 T_d 的元素（SO.TD：I·R_k、I·i_k）與 χ^{ℓ±}(IR) = ±χ^ℓ(R)，配 O 的特徵標（T_d ≅ O）。"""
    f = np.zeros(5)
    for g in SO.TD.names:
        base = g.split(".")[-1]
        ax, w = SO.O_TURNS[base]
        tr = chi(l, w) * (parity if g.startswith("I.") else 1)
        k = next(i for i, C in enumerate(SO.O.classes) if base in C)
        f += np.array([float(SO.O_CHARS[a, k]) for a in range(5)]) * tr
    return tuple(int(round(x / 24)) for x in f)


# ======================================================================== 卡片上的數字（plan_ch05 §6；check() 逐字比對）
def _row4(v) -> str:
    return " ".join(_f4(x) for x in v)


def card_values() -> dict:
    """plan_ch05.md §6 每張卡上的數字（字串；小數照卡上位數）。"""
    out = {}
    out["H501_params"] = ", ".join(str(x) for x in N_PARAMS)
    out["H501_refl60"] = _row4(o2_refl(60).ravel())
    out["H502_R"] = " / ".join(_row4(r) for r in euler_matrix(*EX538))
    out["H502_su2"] = " / ".join(" ".join(f"{z.real:.4f}{z.imag:+.4f}i" for z in r) for r in su2_euler(*EX538))
    out["H503_darboux"] = " ".join(_f2(x) for x in euler_to_darboux(*EX538))
    out["H503_seq"] = "; ".join(f"{w:g}:" + ",".join(f"{x:.1f}" for x in e)
                                for w, e in zip(FIG535_BOOK, darboux_sequence(80, 33.69, list(FIG535_BOOK))))
    sc = small_commutator()
    out["H504_comm"] = f"{sc['angle']:.2f} {sc['second']:.2f} {math.degrees(math.acos(sc['axis'][2])):.2f}"
    out["H505_cone"] = " ".join(_f2(cone_theta(2, m)) for m in (2, 1, 0, -1, -2)) + f" | {math.sqrt(6):.4f}"
    out["H506_D1"] = " / ".join(_row4(r) for r in wigner_d(1, 60))
    out["H506_book"] = _f4(d2_book_orthogonality(0.7))
    he = hamilton_example()
    out["H507_prod"] = f"{he['omega']:.2f} " + _row4(he["axis"]) + " | " + " ".join(_f2(a) for a in he["arcs"])
    out["H507_cos"] = _f4(math.cos(math.radians(60)) * math.cos(math.radians(30)))
    out["H508_j2b90"] = _row4(sg_probs(2, 1, 90))
    out["H508_j2b60"] = _row4(sg_probs(2, 1, 60))
    out["H508_j6"] = _row4(sg_amps(6, 6, 60))
    p = sg_probs(20, 20, 60)
    ms = np.arange(20, -21, -1)
    out["H508_j20"] = f"{float(ms @ p):.4f} {math.sqrt(float(ms ** 2 @ p) - float(ms @ p) ** 2):.4f} {p.max():.4f}"
    out["H508_cone"] = " ".join(f"{cone_angle(j)[0]:.4f}" for j in (1e4, 1e6)) + f" | {cone_angle(6)[0]:.2f}"
    out["H509_sphere"] = " ".join(f"{e:g}({d})" for J in range(4) for e, d in rotor_formula("sphere", J))
    out["H509_prolate"] = " ".join(f"{e:g}" for e, d in rotor_formula("symmetric", 2, C=2.0))
    out["H509_oblate"] = " ".join(f"{e:g}" for e, d in rotor_formula("symmetric", 2, C=0.5))
    pr = problem_554()["c"]
    out["H509_asym"] = " ".join(_f4(x) for x in pr[1]) + " | " + " ".join(_f4(x) for x in pr[2])
    out["H509_zeeman"] = " ".join(f"{e:.1f}" for e, d in rotor_formula("zeeman", 1, g=0.2))
    out["H509_stark"] = " ".join(f"{e:.1f}" for e, d in rotor_formula("stark", 2, mu=0.2))
    out["H510_add"] = _f4(addition_theorem(3, 40, 30, 70, 100)[0])
    out["H511_chi2"] = " ".join(str(int(round(chi(2, w)))) for w in O_CLASS_ANGLES)
    out["H511_corr"] = "; ".join(o3_label(l) for l in range(5))
    out["H512_inv"] = " ".join(str(invariant_count(k)[0]) for k in range(9))
    mc = multipole_coeffs()
    out["H512_coef"] = f"{mc[(0, 0)].real:.4f} {mc[(4, 0)].real:.4f} {mc[(4, 4)].real:.4f}"
    out["H512_V4"] = f"{float(V4.subs({_x: 1, _y: 0, _z: 0})):.4f} {float(V4.subs({_x: 1 / math.sqrt(3), _y: 1 / math.sqrt(3), _z: 1 / math.sqrt(3)})):.4f}"
    out["H512_N2N3"] = " ".join(f"{e:.1f}{l}" for e, l, d in sf6_levels(2) + sf6_levels(3))
    out["H512_J20"] = " ".join(l.replace("_", "") for e, l, d in sf6_levels(20)[::-1])
    out["H513_d6"] = "; ".join("".join(map(str, o3_to_d6(l))) for l in range(9))
    out["H513_ray"] = "; ".join("".join(map(str, half_correlation("D6", Fraction(2 * k + 1, 2)))) for k in range(6))
    out["H514_ray"] = "; ".join("".join(map(str, half_correlation("O", Fraction(2 * k + 1, 2)))) for k in range(6))
    out["H515_coulomb"] = " ".join(_f4(coulomb_energy(n)) for n in range(1, 5)) + " | " + " ".join(str(n * n) for n in range(1, 5))
    out["H515_kepler"] = _f4(kepler_check()["e2"][0][1])
    rt = su3_cartan()["roots"]
    out["H516_roots"] = " / ".join(f"{float(rt[a][0]):.4f} {float(rt[a][1]):.4f}" for a in (1, 2, 3))
    out["H516_osc"] = " ".join(str(osc_shell(N)["count"]) for N in range(5)) + " | " + ", ".join(map(str, COULOMB_CLOSURES)) \
        + " | " + ", ".join(map(str, OSC_CLOSURES))
    return out


CARD_EXPECTED = {
    "H501_params": "0, 1, 3, 6, 10, 15",
    "H501_refl60": "0.5000 0.8660 0.8660 -0.5000",
    "H502_R": "-0.6099 -0.5640 0.5567 / 0.7350 -0.1401 0.6634 / -0.2962 0.8138 0.5000",
    "H502_su2": "0.4330-0.7500i -0.4924-0.0868i / 0.4924-0.0868i 0.4330+0.7500i",
    "H503_darboux": "80.00 33.69 128.68",
    "H503_seq": "0:-10.0,0.0,10.0; 60:15.7,32.2,35.7; 128.68:50.0,60.0,70.0; 180:80.0,67.4,100.0; 240:114.8,57.4,134.8; "
                "300:144.3,32.2,164.3; 360:170.0,0.0,190.0; 420:195.7,-32.2,215.7; 488.68:230.0,-60.0,250.0; "
                "540:260.0,-67.4,280.0; 600:294.8,-57.4,314.8; 660:324.3,-32.2,344.3",
    "H504_comm": "-6.31 -6.37 13.38",
    "H505_cone": "35.26 65.91 90.00 114.09 144.74 | 2.4495",
    "H506_D1": "0.7500 -0.6124 0.2500 / 0.6124 0.5000 -0.6124 / 0.2500 0.6124 0.7500",
    "H506_book": "1.3174",
    "H507_prod": "120.00 0.5774 0.5774 -0.5774 | 45.00 45.00 60.00",
    "H507_cos": "0.4330",
    "H508_j2b90": "0.2500 0.2500 0.0000 0.2500 0.2500",
    "H508_j2b60": "0.4219 0.0000 0.2813 0.2500 0.0469",
    "H508_j6": "0.1780 0.3560 0.4820 0.5080 0.4400 0.3213 0.2004 0.1071 0.0489 0.0188 0.0060 0.0015 0.0002",
    "H508_j20": "10.0000 2.7386 0.1444",
    "H508_cone": "1.1459 0.1146 | 44.42",
    "H509_sphere": "0(1) 2(9) 6(25) 12(49)",
    "H509_prolate": "6 7 10",
    "H509_oblate": "4 5.5 6",
    "H509_asym": "0.4167 0.6667 0.7500 | 1.2324 1.4167 1.6667 2.4167 2.4343",
    "H509_zeeman": "1.8 2.0 2.2",
    "H509_stark": "6.0 6.2 6.8",
    "H510_add": "-0.4457",
    "H511_chi2": "5 -1 1 -1 1",
    "H511_corr": r"A_{1g}; T_{1u}; E_g\oplus T_{2g}; A_{2u}\oplus T_{1u}\oplus T_{2u}; A_{1g}\oplus E_g\oplus T_{1g}\oplus T_{2g}",
    "H512_inv": "1 1 2 3 4 5 7 8 10",
    "H512_coef": "6.0000 3.5000 2.0917",
    "H512_V4": "0.4000 -0.2667",
    "H512_N2N3": "-3.6T_2 2.4E -38.4A_2 -8.4T_2 15.6T_1",
    "H512_J20": "A1 T1 E T2 T1 E T2 A2 T2 T1 A1 T1 E T2 T1 E T2",
    "H513_d6": "100000; 010001; 101001; 011111; 102111; 012112; 212112; 122113; 213113",
    "H513_ray": "100; 101; 111; 121; 122; 222",
    "H514_ray": "100; 001; 011; 111; 102; 112",
    "H515_coulomb": "-0.5000 -0.1250 -0.0556 -0.0312 | 1 4 9 16",
    "H515_kepler": "0.2809",
    "H516_roots": "0.2887 0.5000 / 0.2887 -0.5000 / 0.5774 0.0000",
    "H516_osc": "1 3 6 10 15 | 2, 10, 28, 60 | 2, 8, 20, 40, 70",
}


# ======================================================================== check()
def _check_51_53():
    assert N_PARAMS == [n * n - n * (n + 1) // 2 for n in range(1, 7)] == [0, 1, 3, 6, 10, 15]
    for M, d in ((SIGMA_XY, -1), (INV3, -1), (I_Z, 1), (rz(37), 1), (INV3 @ rz(37), -1)):
        assert is_orthogonal(M) and round(np.linalg.det(M)) == d
        v = [M @ e for e in np.eye(3)]
        assert abs(v[0] @ np.cross(v[1], v[2]) - d) < 1e-12                       # (5.1.8)
    assert np.allclose(I_Z, rz(180))
    for m in (1, 2, 3):
        assert sp.simplify(o2_standing(m) - sp.Matrix([[sp.cos(m * PHI), -sp.sin(m * PHI)], [sp.sin(m * PHI), sp.cos(m * PHI)]])) == sp.zeros(2)
    Rf = o2_refl(60)
    assert np.allclose(Rf @ np.array([math.cos(math.radians(30)), math.sin(math.radians(30))]),
                       [math.cos(math.radians(30)), math.sin(math.radians(30))]) and round(np.linalg.det(Rf)) == -1
    assert not np.allclose(o2_refl(60) @ o2_rot(30), o2_rot(30) @ o2_refl(60))   # O2 非 Abel
    assert O2_MOVING_SIGMA_MINUS_BOOK != 1 and dinfh_trace_check(40) and dinfh_trace_check(110)
    for n in range(2, 7):
        assert dinfh_to_dnh(n) == dinfh_to_dnh_mrule(n)
    # Euler、Darboux：(5.3.8)、(5.3.9)、Fig. 5.3.5
    ph, th, om = euler_to_darboux(*EX538)
    assert (_f2(ph), _f2(th), _f2(om)) == ("80.00", "33.69", "128.68")
    assert np.allclose(darboux_from_matrix(euler_matrix(*EX538)), (ph, th, om))
    assert np.allclose(darboux_matrix(ph, th, om), euler_matrix(*EX538)) and np.allclose(darboux_matrix_cranks(ph, th, om), euler_matrix(*EX538))
    seq = darboux_sequence(80, 33.69, list(FIG535_BOOK))
    for (w, bk), e in zip(FIG535_BOOK.items(), seq):
        assert all(abs(a - b) < 0.06 for a, b in zip(e, bk)), (w, e, bk)
        assert abs((e[2] - e[0]) - 20) < 1e-9                                       # γ − α = 20°
        assert np.allclose(euler_matrix(*e), darboux_matrix(80, 33.69, w))
    assert np.allclose(euler_matrix(*EX539), euler_matrix(*EX538)) and np.allclose(su2_euler(*EX539), -su2_euler(*EX538))
    assert np.allclose(su2_darboux(80, 33.69, 488.68), -su2_darboux(80, 33.69, 128.68))          # (5.3.9)、(5.3.10)
    assert np.allclose(su2_darboux(ph, th, om), su2_euler(*EX538), atol=1e-12)
    assert np.allclose(euler_matrix(*euler_alt(*EX538, -1)), euler_matrix(*EX538)) and euler_alt(*EX538, -1) == EX539
    rng = np.random.default_rng(3)
    for _ in range(6):
        a, b, c = rng.uniform(-170, 170), rng.uniform(1, 179), rng.uniform(-170, 170)
        assert np.allclose(euler_from_matrix(euler_matrix(a, b, c)), (a, b, c))
        assert np.allclose(SN.su2(quat_euler(a, b, c)), su2_euler(a, b, c))
        d = euler_to_darboux(a, b, c)
        assert np.allclose(darboux_matrix(*d), euler_matrix(a, b, c)) and np.allclose(su2_darboux(*d), su2_euler(a, b, c))
        assert np.allclose(darboux_to_euler(*d), (a, b, c)) or np.allclose(euler_matrix(*darboux_to_euler(*d)), euler_matrix(a, b, c))
        # 實驗室與本體：本體軸依序 z̄(α)、ȳ(β)、z̄(γ) ＝ 實驗室 z(γ)、y(β)、z(α)（(5.3.5)）
        assert np.allclose(body_frame_product([("z", a), ("y", b), ("z", c)]), euler_matrix(a, b, c))
        assert np.allclose(lab_frame_product([("z", c), ("y", b), ("z", a)]), euler_matrix(a, b, c))
        S, R1, R2 = euler_matrix(*rng.uniform(0, 90, 3)), euler_matrix(a, b, c), euler_matrix(c, b, a)
        assert np.allclose(lab_op(R1)(body_op(R2)(S)), body_op(R2)(lab_op(R1)(S)))                 # (5.3.6)
        assert np.allclose(lab_op(R1)(np.eye(3)), np.linalg.inv(body_op(R1)(np.eye(3))))          # (5.3.7)：R̄⁻¹|1⟩ = R|1⟩
    for g, (e, d) in problem_553().items():
        ax, w = SO.O_TURNS[g]
        assert np.allclose(euler_matrix(*e), SN.quat_to_rot(SN.quat(ax, w)), atol=1e-3) and np.allclose(darboux_matrix(*d), euler_matrix(*e), atol=1e-3)
    assert problem_553() == {"R1^2": ((180.0, 180.0, 0.0), (0.0, 90.0, 180.0)), "r1": ((0.0, 90.0, 90.0), (45.0, 54.74, 120.0)),
                             "i1": ((0.0, 90.0, 180.0), (0.0, 45.0, 180.0)), "R1": ((-90.0, 90.0, 90.0), (0.0, 90.0, 90.0))}
    p = problem_533()
    assert p["back"] and p["back2"] and abs(p["det"] - 1) < 1e-12
    assert tuple(_f4(x) for x in p["euler"]) == ("45.0000", "54.7356", "-90.0000")
    assert tuple(_f4(x) for x in p["darboux"]) == ("157.5000", "126.4751", "69.7356") and np.allclose(p["via552"], p["darboux"])


def _check_54():
    assert comm(GX, GY) == GZ and comm(GY, GZ) == GX and comm(GZ, GX) == GY                 # (5.4.10)–(5.4.11)
    P = gz_projectors()
    assert all(sp.simplify(P[k] - GZ_PROJ_BOOK[k]) == sp.zeros(3) for k in GZ_PROJ_BOOK)
    w = sp.Symbol("w", real=True)
    E = SM.mat_func(GZ, lambda e: sp.exp(e * w))
    assert sp.simplify(E.rewrite(sp.cos) - sp.Matrix([[sp.cos(w), -sp.sin(w), 0], [sp.sin(w), sp.cos(w), 0], [0, 0, 1]])) == sp.zeros(3)
    assert sp.factor(GZ.charpoly(sp.Symbol("e")).as_expr()) == sp.factor(sp.Symbol("e") ** 3 + sp.Symbol("e"))
    assert np.allclose(gen_exp((1, 2, 2), 75), darboux_matrix(*darboux_from_matrix(gen_exp((1, 2, 2), 75))))
    assert np.allclose(gen_exp((0, 0, 1), 33), rz(33)) and np.allclose(gen_exp((0, 1, 0), 33), ry(33))
    sc = small_commutator()
    assert abs(sc["angle"] - sc["second"]) / abs(sc["second"]) < 0.02 and np.abs(sc["R"] - sc["approx"]).max() < 0.03
    sc2 = small_commutator(0.01)
    assert abs(sc2["angle"] - sc2["second"]) < 1e-4 and np.allclose(sc["path"][-1], sc["R"] @ sc["path"][0])
    for j in (Fraction(1, 2), 1, Fraction(3, 2), 2, 3, Fraction(7, 2)):
        M = jmat(j)
        jf = float(j)
        assert _close(comm(M["Jx"], M["Jy"]), 1j * M["Jz"]) and _close(comm(M["Jy"], M["Jz"]), 1j * M["Jx"])
        assert _close(M["J2"], jf * (jf + 1) * np.eye(int(2 * jf + 1)))
        assert _close(np.diag(M["Jz"]).real, _ms_num(j)) and _close(M["Jp"], M["Jx"] + 1j * M["Jy"])
    Je = jmat(2, exact=True)
    assert Je["Jx"][1, 2] == SQ6 / 2 and J2_JX_BOOK[1, 2] == SQ3 / 2 and Je["Jx"][0, 1] == 1
    J32, J1 = jmat(Fraction(3, 2), exact=True), jmat(1, exact=True)
    assert [J32["Jp"][k, k + 1] for k in range(3)] == [SQ3, 2, SQ3] and J1["Jx"] == SQ2 / 2 * sp.Matrix([[0, 1, 0], [1, 0, 1], [0, 1, 0]])
    assert J1["Jy"] == SQ2 / 2 * sp.Matrix([[0, -IU, 0], [IU, 0, -IU], [0, IU, 0]])
    assert sum(1 for i in range(5) for k in range(5) if sp.simplify(Je["Jx"][i, k] - J2_JX_BOOK[i, k]) != 0) == 4
    for j in (Fraction(1, 2), 1, Fraction(3, 2), 2, Fraction(5, 2), 3, 4):
        a, b, c = 23.0, 71.0, -38.0
        assert _close(wigner_D(j, a, b, c), D_expm(j, a, b, c), 1e-10) and _close(wigner_d(j, b), D_schwinger(j, b), 1e-12)
        D = wigner_D(j, a, b, c)
        assert _close(D @ D.conj().T, np.eye(D.shape[0]))
        assert _close(D.conj(), wigner_D(j, -c, -b, -a).T)                                  # 更正的 (5.4.46) 後一式
        ms = _ms_num(j)
        assert _close(D.conj(), ((-1.0) ** np.round(ms[:, None] - ms[None, :])) * D[::-1, ::-1])   # 習題 5.4.6 最後一式
        assert d_symmetry_checks(j) == dict(neg=True, transpose=True, minus_beta=True, pi_minus=True)
        u = su2_euler(a, b, c)
        assert _close(D_from_u(j, u), D, 1e-10)                                             # 習題 5.4.5
    assert not _close(wigner_D(1, 23, 71, -38).conj(), wigner_D(1, -23, -71, 38).T)        # 書上的 (−α−β−γ)
    assert _close(wigner_d(20, 60), D_expm(20, 0, 60, 0), 1e-9)
    bb = 0.7
    for S_, j in ((D1_SYM, 1), (D32_SYM, Fraction(3, 2)), (D2_SYM, 2)):
        assert _close(np.array(S_.subs(BETA, bb).evalf(), dtype=float), wigner_d(j, math.degrees(bb)))
    assert sp.simplify(d_sym(1) - D1_SYM) == sp.zeros(3)
    assert _f4(d2_book_orthogonality(0.7)) == "1.3174"
    Dbk = np.array(D2_BOOK.subs(BETA, bb).evalf(), dtype=float)
    assert sorted(zip(*np.nonzero(np.abs(Dbk - wigner_d(2, math.degrees(bb))) > 1e-9))) == D2_TYPOS
    assert _close(wigner_D(0.5, 360, 0, 0), -np.eye(2)) and _close(wigner_D(0.5, 180, 0, 0) @ wigner_D(0.5, 180, 0, 0), -np.eye(2))
    assert _close(wigner_D(1, 360, 0, 0), np.eye(3)) and _close(wigner_D(1.5, 0, 360, 0), -np.eye(4))
    assert _close(wigner_d(0.5, 60), [[math.cos(math.pi / 6), -0.5], [0.5, math.cos(math.pi / 6)]])
    # Schwinger：兩個振子的 Fock 空間（截斷 n ≤ 4）裡 J± = a†a、J_z = (n↑ − n↓)/2（(5.4.37)、(5.4.39)）
    nmax = 4
    a = np.diag(np.sqrt(np.arange(1, nmax + 1)), 1)
    one = np.eye(nmax + 1)
    au, ad = np.kron(a, one), np.kron(one, a)
    Jm, Jp = ad.T @ au, au.T @ ad
    Jz = (au.T @ au - ad.T @ ad) / 2
    keep = [i for i in range((nmax + 1) ** 2) if i // (nmax + 1) + i % (nmax + 1) <= nmax]
    sub = np.ix_(keep, keep)
    assert _close(comm(Jz, Jp)[sub], Jp[sub]) and _close(comm(Jz, Jm)[sub], -Jm[sub])
    for nu, nd in ((2, 1), (1, 2), (3, 0)):
        v = np.zeros((nmax + 1) ** 2)
        v[nu * (nmax + 1) + nd] = 1
        j, m = (nu + nd) / 2, (nu - nd) / 2
        w2 = Jm @ v
        assert abs(np.linalg.norm(w2) - math.sqrt((j + m) * (j - m + 1))) < 1e-12               # (5.4.39) ＝ (5.4.23b)
    assert problem_542()["equiv"] and appendix_e_check() and appendix_e_check(2) and appendix_e_check(Fraction(7, 2))


def _ms_num(j) -> np.ndarray:
    return np.array([float(m) for m in _ms(j)])


def _check_55():
    # (5.5.1)–(5.5.5)
    for t in ((80, 33.69, 128.68), (210, 120, 300), (0, 90, 180)):
        cr = su2_euler(t[0], 0, 0) @ su2_euler(0, t[1], 0) @ su2_euler(t[2], 0, 0) @ su2_euler(0, -t[1], 0) @ su2_euler(-t[0], 0, 0)
        assert _close(cr, su2_darboux(*t)) and _close(SN.su2(SN.quat(axis_vec(t[0], t[1]), t[2])), su2_darboux(*t))
    a, b, c = EX538
    A, B, C = (math.radians(x) / 2 for x in (a, b, c))
    ph, th, om = (math.radians(x) for x in euler_to_darboux(a, b, c))
    W = om / 2
    assert abs(math.cos(A + C) * math.cos(B) - math.cos(W)) < 1e-12 and abs(math.sin(A + C) * math.cos(B) - math.sin(W) * math.cos(th)) < 1e-12
    assert abs(math.cos(C - A) * math.sin(B) - math.sin(W) * math.sin(th) * math.sin(ph)) < 1e-12
    assert abs(math.sin(C - A) * math.sin(B) - math.sin(W) * math.sin(th) * math.cos(ph)) < 1e-12        # (5.5.2a–d)
    J = jmat(0.5)
    for k, nm in enumerate(("Jx", "Jy", "Jz")):
        assert _close(Q_UNITS["xyz"[k]], 2j * J[nm])                                               # (5.5.3)
    rng = np.random.default_rng(7)
    A_, B_ = rng.normal(size=3), rng.normal(size=3)
    sA, sB = sum(A_[i] * PAULI[i] for i in range(3)), sum(B_[i] * PAULI[i] for i in range(3))
    assert _close(sA @ sB, (A_ @ B_) * np.eye(2) + 1j * sum(np.cross(A_, B_)[i] * PAULI[i] for i in range(3)))   # (5.5.7)
    for t1, t2 in (((0, 90, 90), (90, 90, 90)), ((30, 40, 70), (200, 110, 150))):
        cc, v = su2_product(t1, t2)
        U = su2_darboux(*t2) @ su2_darboux(*t1)
        assert _close(U, cc * np.eye(2) - 1j * sum(v[i] * PAULI[i] for i in range(3)))            # (5.5.10)
        hp = SN.hamilton_product(axis_vec(t1[0], t1[1]), t1[2], axis_vec(t2[0], t2[1]), t2[2])
        assert np.allclose(hp["q"], [cc, *v]) and (np.allclose(hp["q_arc"], hp["q"]) or np.allclose(hp["q_arc"], -hp["q"]))
    he = hamilton_example()
    assert (_f2(he["omega"]), _row4(he["axis"])) == ("120.00", "0.5774 0.5774 -0.5774") and np.allclose(he["arcs"], (45, 45, 60))
    assert np.allclose(he["hp"]["axis"], he["axis"]) and abs(he["hp"]["deg"] - 120) < 1e-9
    assert D2S.mul("Ry", "Rx") == (-1, "Rz") and D2S.mul("Rx", "Ry") == (1, "Rz")
    # Stern–Gerlach
    assert _row4(sg_probs(0.5, 0.5, 60)) == "0.7500 0.2500" and _close(sg_probs(0.5, 0.5, 60), [math.cos(math.radians(30)) ** 2, math.sin(math.radians(30)) ** 2])
    for j, n, bt in ((2, 1, 50), (6, 6, 60), (20, 20, 60), (20, 17, 30), (Fraction(7, 2), Fraction(1, 2), 80)):
        assert abs(sg_probs(j, n, bt).sum() - 1) < 1e-10 and abs(mean_jz(j, n, bt) - float(n) * math.cos(math.radians(bt))) < 1e-9
    assert _close(sg_probs(2, 1, 180), [0, 0, 0, 1, 0])
    cb, sb = math.cos(math.radians(30)), math.sin(math.radians(30))
    assert _close(sg_amps(6, 6, 60), [math.sqrt(math.comb(12, 6 + m)) * cb ** (6 + m) * sb ** (6 - m) for m in range(6, -7, -1)])
    p = sg_probs(20, 20, 60)
    assert abs(math.sqrt(float(np.arange(20, -21, -1) ** 2 @ p) - 100) - math.sqrt(10) * math.sin(math.radians(60))) < 1e-9
    ca = cone_angle(1e4)
    assert f"{ca[0]:.2f}" == "1.15" and f"{cone_angle(1e6)[0]:.3f}" == "0.115" and abs(ca[1] - ca[0]) < 1e-3 and abs(2 * ca[2] - ca[1]) < 1e-9
    # 轉子
    for J in range(4):
        assert all(abs(x[0] - y[0]) < 1e-9 and x[1] == y[1] for x, y in zip(rotor_levels("sphere", J), rotor_formula("sphere", J)))
        for kind, kw in (("symmetric", dict(C=2.0)), ("symmetric", dict(C=0.5)), ("zeeman", dict(g=0.2)), ("stark", dict(mu=0.2))):
            num, frm = rotor_levels(kind, J, **kw), rotor_formula(kind, J, **kw)
            assert len(num) == len(frm) and all(abs(x[0] - y[0]) < 1e-9 and x[1] == y[1] for x, y in zip(num, sorted(frm)))
        for (A_, B_, C_) in ((0.5, 0.25, 1 / 6), (0.3, 0.3, 0.9)):
            num = [e for e, d in rotor_levels("asymmetric", J, A=A_, B=B_, C=C_) for _ in range(d // (2 * J + 1))]
            assert np.allclose(num, np.linalg.eigvalsh(asym_matrix_formula(J, A_, B_, C_)))
    o = rotor_ops(2)
    L_, Bd = o["lab"], o["body"]
    assert _close(comm(Bd["Jx"], Bd["Jy"]), -1j * Bd["Jz"]) and _close(comm(L_["Jx"], L_["Jy"]), 1j * L_["Jz"])   # (5.5.46)
    assert all(_close(comm(L_[a], Bd[b]), 0 * L_[a]) for a in L_ for b in Bd) and _close(sum(Bd[k] @ Bd[k] for k in Bd), o["J2"])
    assert problem_554() == {"a": {1: [0.75, 0.75, 1.0], 2: [2.0, 2.0, 2.75, 2.75, 3.0]}, "b": {1: [0.1, 0.55, 0.55], 2: [0.3, 0.75, 0.75, 2.1, 2.1]},
                             "c": {1: [0.4167, 0.6667, 0.75], 2: [1.2324, 1.4167, 1.6667, 2.4167, 2.4343]}}
    assert problem_555() == 0
    # Euler 積分（(5.5.22)–(5.5.28)）：Gauss–Legendre 在 cos β、α、γ 均勻 → ∫d(αβγ) = 1、(2j+1)∫𝒟*𝒟 = δ
    xg, wg = np.polynomial.legendre.leggauss(8)
    angs = np.arange(8) * 45.0
    acc = np.zeros((3, 3, 3, 3), complex)
    for x, wq in zip(xg, wg):
        bdeg = math.degrees(math.acos(x))
        for al in angs:
            for ga in angs:
                D = wigner_D(1, al, bdeg, ga)
                acc += wq / 2 / 64 * np.einsum("ab,cd->abcd", D.conj(), D)
    assert _close(3 * np.einsum("abab->ab", acc), np.ones((3, 3))) and abs(3 * acc[0, 0, 1, 1]) < 1e-12
    # 球諧（(5.5.60)–(5.5.75)）
    for l in range(5):
        for m in range(-l, l + 1):
            for th, phh in ((0.7, 1.1), (2.2, -0.4)):
                assert abs(ylm(l, m, th, phh) - ylm_from_D(l, m, th, phh)) < 1e-12
    th, phh = 0.7, 1.1
    x, y, z = math.sin(th) * math.cos(phh), math.sin(th) * math.sin(phh), math.cos(th)
    n2 = math.sqrt(4 * math.pi / 5)
    want = {2: math.sqrt(3 / 8) * (x + 1j * y) ** 2, 1: -math.sqrt(1.5) * (x + 1j * y) * z, 0: (2 * z * z - x * x - y * y) / 2,
            -1: math.sqrt(1.5) * (x - 1j * y) * z, -2: math.sqrt(3 / 8) * (x - 1j * y) ** 2}
    assert all(abs(n2 * ylm(2, m, th, phh) - v) < 1e-12 for m, v in want.items())                # (5.6.14)
    T15 = {"T2_1": (-(ylm(2, 1, th, phh) + ylm(2, -1, th, phh)) / math.sqrt(2), 1j * math.sqrt(3) * y * z),
           "T2_2": ((-1j * ylm(2, 1, th, phh) + 1j * ylm(2, -1, th, phh)) / math.sqrt(2), 1j * math.sqrt(3) * x * z),
           "T2_3": ((ylm(2, 2, th, phh) - ylm(2, -2, th, phh)) / math.sqrt(2), 1j * math.sqrt(3) * x * y),
           "E2": ((ylm(2, 2, th, phh) + ylm(2, -2, th, phh)) / math.sqrt(2), math.sqrt(3) * (x * x - y * y) / 2)}
    assert all(abs(n2 * a_ - b_) < 1e-12 for a_, b_ in T15.values())                                # (5.6.15)
    for l in (1, 2, 3, 5):
        lhs, rhs = addition_theorem(l, 40, 30, 70, 100)
        assert abs(lhs - rhs) < 1e-12
    for l, m in ((1, 0), (2, 1), (3, -2)):
        a_, b_ = rotated_ylm(l, m, (20, 50, 70), 0.7, 1.1)
        assert abs(a_ - b_) < 1e-12
    D1c = wigner_D(1, 35, 50, 0)[:, 1].conj()
    assert _close(D1c, [-np.exp(1j * math.radians(35)) * math.sin(math.radians(50)) / math.sqrt(2), math.cos(math.radians(50)),
                        np.exp(-1j * math.radians(35)) * math.sin(math.radians(50)) / math.sqrt(2)])                  # (5.5.64a)
    # (5.5.81)–(5.5.83)：角速度矩陣與其逆、動量 = 逆矩陣的轉置
    a_, b_, g_ = (math.radians(t) for t in (25, 50, 80))
    Mlab = np.array([[0, -math.sin(a_), math.cos(a_) * math.sin(b_)], [0, math.cos(a_), math.sin(a_) * math.sin(b_)], [1, 0, math.cos(b_)]])
    Mlab_inv = np.array([[-math.cos(a_) / math.tan(b_), -math.sin(a_) / math.tan(b_), 1], [-math.sin(a_), math.cos(a_), 0],
                         [math.cos(a_) / math.sin(b_), math.sin(a_) / math.sin(b_), 0]])
    Mbody = np.array([[-math.sin(b_) * math.cos(g_), math.sin(g_), 0], [math.sin(b_) * math.sin(g_), math.cos(g_), 0], [math.cos(b_), 0, 1]])
    Mbody_inv = np.array([[-math.cos(g_) / math.sin(b_), math.sin(g_) / math.sin(b_), 0], [math.sin(g_), math.cos(g_), 0],
                          [math.cos(g_) / math.tan(b_), -math.sin(g_) / math.tan(b_), 1]])
    assert np.allclose(Mlab @ Mlab_inv, np.eye(3)) and np.allclose(Mbody @ Mbody_inv, np.eye(3))
    assert np.allclose(euler_matrix(25, 50, 80).T @ Mlab, Mbody)                                    # ω_body = Rᵀ ω_lab


def _check_56():
    for l in range(21):
        assert o3_to_o(l) == O3_O_BOOK[l]
        f = o3_to_o(l)
        assert f[0] + f[1] + 2 * f[2] + 3 * f[3] + 3 * f[4] == 2 * l + 1
    for l in range(9):
        assert o3_to_o_matrices(l) == o3_to_o(l)
    for l in range(9):
        assert tuple(a - b for a, b in zip(o3_to_o(l + 12), o3_to_o(l))) == (1, 1, 2, 3, 3)          # (5.6.6)
    for j in (2, 3, Fraction(5, 2), 7):
        for w in (0, 60, 90, 120, 180, 250):
            assert abs(chi(j, w) - float(chi_exact(j, w))) < 1e-9 and abs(chi(j, w) - np.trace(D_axis(j, (1, 2, 3), w)).real) < 1e-9
    assert tuple(int(round(chi(2, w))) for w in O_CLASS_ANGLES) == (5, -1, 1, -1, 1)
    for k in range(9):
        a, b, c = invariant_count(k)
        assert a == b == c
    assert invariant_count(6) == (7, 7, 7) and o3_to_o(9)[0] == 1 and all(o3_to_o(l)[0] == 0 for l in (1, 3, 5, 7))
    d = d_orbital_states()
    for key, v in D_ORBITALS_BOOK.items():
        w = np.array(v, complex) / np.linalg.norm(v)
        assert _close(d[key], w) and _close(d_orbital_projector_full(key) / np.linalg.norm(d_orbital_projector_full(key)), w)
    assert _close(d["sums"][3], [1, 0, 0, 0, -1]) and not _close(d["sums"][3], T2_SUM_BOOK3)
    for g, e in COSET_EULER.items():
        assert np.allclose(euler_matrix(*e), _num(SO.O.mats[g]).real)
        assert _close(wigner_D(2, *e)[:, 0], np.array([complex(sp.N(x)) for x in D2_COL_BOOK[g]]))   # (5.6.11)
    assert v4_identities() == (0, 0)
    mc = multipole_coeffs()
    assert set(mc) == {(0, 0), (4, 0), (4, 4), (4, -4)} and abs(mc[(4, 0)] - 3.5) < 1e-12 and abs(mc[(4, 4)] - math.sqrt(70) / 4) < 1e-12
    for r in ([0.1, 0.05, -0.07], [0.0, 0.0, 0.12]):
        r = np.array(r)
        e1 = abs(oct_potential(r) - oct_multipole(r))
        e2 = abs(oct_potential(r / 2) - oct_multipole(r / 2))
        assert e1 < 1e-4 and 40 < e1 / e2 < 90                                                     # 下一項 ∝ r⁶
    lv = sf6_levels(2) + sf6_levels(3)
    assert [l for e, l, dd in lv] == ["T_2", "E", "A_2", "T_2", "T_1"] and [round(e, 1) for e, l, dd in lv] == [-3.6, 2.4, -38.4, -8.4, 15.6]
    for J in range(31):
        assert _count(wheel_o(J), O_IRREPS) == o3_to_o(J) and _count(wheel_d6(J), D6_IRREPS) == o3_to_d6(J) == o3_to_d6_mrule(J)
    for J in range(26):
        assert [l for e, l, dd in sf6_levels(J)][::-1] == wheel_o(J)                                 # 輪的次序 ＝ 能階次序（t > 0）
    cl = sf6_clusters(20, 0.01)
    assert cl[0] == ["A_1", "T_1", "E"] and cl[1] == ["T_2", "T_1"] and cl[2] == ["E", "T_2", "A_2"] and cl[-1] == ["T_1", "E", "T_2"]
    for l in range(11):
        assert o3_to_d6(l) == O3_D6_BOOK[l]
    for N in range(7):
        dd = d180(N)
        ms = _ms(N)
        want = np.array([[(-1) ** int(N + m) if m == -n else 0 for n in ms] for m in ms], float)
        assert _close(dd, want)                                                                   # (5.6.40)


def _check_57():
    assert d6s_table_check() == 121
    ca = d6s_class_algebra()
    assert ca["c2c2"] == {"1": 2, "h^2": -1, "h^4": 1} and ca["c3c3"] == {"1": 2, "h^2": 1, "h^4": -1} and ca["c2c3"] == {"h": 1, "h^5": -1}
    assert ca["c3cube"] == {"h": 3, "h^5": -3}                                                   # c3³ = 3c3
    for a in range(3):
        for b in range(3):
            assert {g: sp.simplify(v) for g, v in ca["PP"][a][b].items() if sp.simplify(v) != 0} == \
                ({g: sp.simplify(v) for g, v in ca["P"][a].items()} if a == b else {})
    tot = {}
    for P_ in ca["P"]:
        for g, v in P_.items():
            tot[g] = sp.simplify(tot.get(g, 0) + v)
    assert {g: v for g, v in tot.items() if v != 0} == {"1": 1}
    ops = D6S.class_ops()
    assert [(h, c) for h, c in ops] == [("1", {"1": 1}), ("h", {"h": 1, "h^5": -1}), ("h^2", {"h^2": 1, "h^4": -1})]
    assert D6S.conj_sum("h") == {"h": 6, "h^5": -6} and D6S.conj_sum("p1") == {} and D6S.conj_sum("h^3") == {}
    full = _ray_chars_full(D6S, D6S_CHARS, D6S_REPS)
    for a, Pa in zip(("E_1", "E_2", "E_3"), ca["P"]):
        for g, v in Pa.items():
            assert sp.simplify(v - sp.Rational(2, 12) * full[a][g]) == 0                           # P^α = (ℓ/°G)Σχ* g
    tr = D6S.trace_char(Fraction(1, 2))
    assert all(abs(tr[g] - float(full["E_1"][g])) < 1e-9 for g in D6S.names)
    assert sum(int(v[0]) ** 2 for v in D6S_CHARS.values()) == 12
    for k in range(12):
        j = Fraction(2 * k + 1, 2)
        assert half_correlation("D6", j) == half_d6_mrule(j) == _count(wheel_d6_half(j), ("E_1", "E_2", "E_3"))
        assert half_correlation("O", j) == _count(wheel_o_half(j), ("E_1", "E_2", "G"))
    for k in range(6):
        j = Fraction(2 * k + 1, 2)
        f, fp = half_correlation("D6", j), half_correlation("D6", j + 6)
        assert tuple(b - a for a, b in zip(f, fp)) == (2, 2, 2)                                   # 週期 6（更正的 (5.7.15)）
    bad = [jl for jl, v in D6S_FREQ_BOOK.items() if half_correlation("D6", Fraction(jl)) != v]
    assert bad == ["1/2", "7/2", "9/2"] and abs(chi(3.5, 60) + math.sqrt(3)) < 1e-12 and abs(chi(4.5, 60) + math.sqrt(3)) < 1e-12
    assert half_correlation("D6", Fraction(7, 2)) != tuple(a + b for a, b in zip((1, 1, 1), half_correlation("D6", Fraction(1, 2))))   # (5.7.15) 照印不成立
    assert [[s for s in r] for r in [[e.replace(r"\mathbf 1", "1").replace("R_", "R") for e in row] for row in D2S.table(["1", "Rx", "Ry", "Rz"], ["1", "Rx", "Ry", "Rz"])]] == D2S_TABLE_BOOK
    assert len(D2S.class_ops()) == 1 and all(_close(D2S.U[k], v) for k, v in D2S_DE_BOOK.items())
    q = Q_UNITS
    assert _close(q["x"] @ q["y"], -q["z"]) and _close(q["x"] @ q["x"], -np.eye(2)) and _close(D2S.U["Rz"], -q["z"])   # (5.7.20) 更正
    assert _close(D2S.U["Rx"] @ D2S.U["Ry"], D2S.U["Rz"])
    for n in range(2, 8):
        ch = cn_spin_chars(n)
        assert all(abs(c ** n + 1) < 1e-12 for c in ch) and (any(abs(c + 1) < 1e-12 for c in ch) == (n % 2 == 1))
        h = _spin((0, 0, 1), 360 / n)
        assert _close(np.linalg.matrix_power(h, n), -np.eye(2))
        k, dims = dn_ray_irreps(n)
        assert k == len(dims) and sum(x * x for x in dims) == 2 * n
    assert [dn_ray_irreps(n)[0] for n in range(2, 8)] == [1, 3, 2, 4, 3, 5]
    assert np.allclose(sorted(cn_spin_chars(3), key=lambda z: z.imag), [np.exp(-1j * math.pi / 3), -1, np.exp(1j * math.pi / 3)])
    assert os_table_check() == 25
    oc = os_class_algebra()
    cr, cR = oc["cr"], oc["cR"]
    assert oc["crcr"] == {"1": 8, **{g: 2 * v for g, v in cr.items()}} and oc["crcR"] == {g: 4 * v for g, v in cR.items()}
    assert oc["cRcR"] == {"1": 6, **{g: 3 * v for g, v in cr.items()}} and oc["cRcube"] == {g: 18 * v for g, v in cR.items()}
    assert len(OS.class_ops()) == 3
    fullO = _ray_chars_full(OS, OS_CHARS, OS_REPS)
    PE1 = {"1": sp.Rational(2, 12), **{g: sp.Rational(1, 12) * v for g, v in cr.items()}, **{g: SQ2 / 12 * v for g, v in cR.items()}}
    PG = {"1": sp.Rational(4, 6), **{g: -sp.Rational(1, 6) * v for g, v in cr.items()}}
    for a, Pa in (("E_1", PE1), ("G", PG)):
        lam = 2 if a == "E_1" else 4
        assert all(sp.simplify(v - sp.Rational(lam, 24) * fullO[a][g]) == 0 for g, v in Pa.items())     # (5.7.24)
        PP = OS.ga_mul(Pa, Pa)
        assert all(sp.simplify(PP.get(g, 0) - Pa.get(g, 0)) == 0 for g in set(PP) | set(Pa))
    trO = OS.trace_char(Fraction(1, 2))
    assert all(abs(trO[g] - float(sp.N(fullO["E_1"][g]))) < 1e-9 for g in OS.names)
    assert sum(int(v[0]) ** 2 for v in OS_CHARS.values()) == 24
    bad = [jl for jl, v in OS_FREQ_BOOK.items() if half_correlation("O", Fraction(jl)) != v]
    assert bad == ["5/2"]
    for k in range(22):
        j = Fraction(2 * k + 1, 2)
        assert [l for e, l, dd in sf6_levels(j)][::-1] == wheel_o_half(j)


def _check_58():
    G = r4_generators()
    lab = {"Lx": ("L", 0), "Ly": ("L", 1), "Lz": ("L", 2), "Kx": ("K", 0), "Ky": ("K", 1), "Kz": ("K", 2)}
    for a, (ta, ia) in lab.items():
        for b, (tb, ib) in lab.items():
            k = 3 - ia - ib
            if ia == ib:
                assert comm(G[a], G[b]) == sp.zeros(4)
                continue
            sgn = 1 if (ia, ib) in ((0, 1), (1, 2), (2, 0)) else -1
            tgt = "K" if ta != tb else "L"
            name = {0: "x", 1: "y", 2: "z"}[k]
            assert comm(G[a], G[b]) == sgn * IU * G[tgt + name]                             # (5.8.16)
    RL, RK = r4_regular("Lz"), r4_regular("Kz")
    assert RL == R4_REG_LZ_BOOK and RK == R4_REG_KZ_BOOK
    PJ1, PK1 = RL * (RL + sp.eye(6)) / 2, RK * (RK + sp.eye(6)) / 2
    PJm, PKm = RL * (RL - sp.eye(6)) / 2, RK * (RK - sp.eye(6)) / 2
    assert PJ1[2, 0] == IU / 2 and PJ1[0, 2] == -IU / 2 and PK1[3, 0] == IU / 2 and PJm[2, 0] == -IU / 2
    cols = {(1, 1): [1, 1, IU, IU], (1, -1): [1, -1, IU, -IU], (-1, 1): [1, -1, -IU, IU], (-1, -1): [1, 1, -IU, -IU]}
    for (s1, s2), v in cols.items():
        Pj = PJ1 if s1 > 0 else PJm
        Pk = PK1 if s2 > 0 else PKm
        assert list((Pj * Pk)[:4, 0]) == [sp.Rational(1, 4) * x for x in v]
    lad = r4_ladders()
    want = {"Mp": (1, 1), "Np": (1, -1), "Nm": (-1, 1), "Mm": (-1, -1)}                       # (5.8.20b)、Fig. 5.8.2
    for k, (rl, rk) in want.items():
        assert comm(G["Lz"], lad[k]) == rl * lad[k] and comm(G["Kz"], lad[k]) == rk * lad[k]
    Mz = (G["Lz"] + G["Kz"]) / 2
    assert comm(lad["Mp"] / 2, lad["Mm"] / 2) == 2 * Mz                                       # (5.8.27)ₓ
    for n in range(1, 5):
        S = coulomb_shell(n)
        L, K, M, N = S["L"], S["K"], S["M"], S["N"]
        assert _close(comm(L["Jx"], L["Jy"]), 1j * L["Jz"]) and _close(comm(L["Jx"], K["Jy"]), 1j * K["Jz"]) and _close(comm(K["Jx"], K["Jy"]), 1j * L["Jz"])
        assert all(_close(comm(M[a], N[b]), 0 * M[a]) for a in M for b in N)                     # (5.8.23)
        assert _close(sum(K[a] @ L[a] for a in K), 0 * S["M2"]) and _close(S["M2"], S["N2"])        # (5.8.31)
        E = -0.5 * np.linalg.inv(2 * S["M2"] + 2 * S["N2"] + np.eye(n * n))
        assert _close(E, coulomb_energy(n) * np.eye(n * n)) and S["dim"] == n * n
        L2 = sum(L[a] @ L[a] for a in L)
        w = sorted(np.linalg.eigvalsh(L2))
        assert np.allclose(w, sorted(l * (l + 1) for l in range(n) for _ in range(2 * l + 1)))
    assert sorted(coulomb_weights(0.5, 0.5)) == [(-1.0, 0.0), (0.0, -1.0), (0.0, 1.0), (1.0, 0.0)] and len(coulomb_weights(1, 1)) == 9
    kc = kepler_check()
    assert kc["spread"] < 1e-8 and all(abs(a - b) < 1e-9 for a, b in kc["e2"])
    pb = runge_lenz_poisson()
    assert pb == dict(e1e2=0, L1e2=0, conserved=0)
    e_ok, e_book = quantum_xr_p()
    assert e_ok < 1e-12 and e_book > 0.1
    d_km, d_2km = quantum_eps_denominator()
    assert d_km < 1e-12 and d_2km > 1e-3
    for (i, j_), (k, l_) in itertools.product(itertools.product(range(1, 4), repeat=2), repeat=2):
        assert comm(E3(i, j_), E3(k, l_)) == (E3(i, l_) if j_ == k else sp.zeros(3)) - (E3(k, j_) if l_ == i else sp.zeros(3))   # (5.8.40)
    assert E3(1, 2) == (U3_K["21"] - U3_G["21"]) / 2 and E3(2, 1) == (U3_K["21"] + U3_G["21"]) / 2
    assert E3(1, 3) == (U3_K["13"] + U3_G["13"]) / 2 and E3(3, 2) == (U3_K["32"] + U3_G["32"]) / 2   # (5.8.39)
    cart = su3_cartan()
    rt = cart["roots"]
    assert rt[1] == (SQ3 / 6, R_(1, 2)) and rt[2] == (SQ3 / 6, -R_(1, 2)) and rt[3] == (SQ3 / 3, 0)
    assert all(sp.simplify(rt[-a][0] + rt[a][0]) == 0 and sp.simplify(rt[-a][1] + rt[a][1]) == 0 for a in (1, 2, 3))
    for a in (1, 2, 3):
        Ca = comm(cart["E"][a], cart["E"][-a])
        assert sp.simplify(Ca - (rt[a][0] * cart["H1"] + rt[a][1] * cart["H2"])) == sp.zeros(3)      # (5.8.42)
    rr = lambda a: sp.Matrix(rt[a])                                                           # noqa: E731
    for a, b, c in ((1, 2, 3), (-1, -2, -3), (3, -2, 1), (3, -1, 2)):
        C = comm(cart["E"][a], cart["E"][b])
        jn = 1
        km = jn + 2 * (rr(a).dot(rr(b))) / rr(b).dot(rr(b))
        coef = sp.sqrt(jn * (km + 1) / 2) * sp.sqrt(rr(a).dot(rr(a)))
        assert sp.simplify(C - coef * cart["E"][c]) == sp.zeros(3) or sp.simplify(C + coef * cart["E"][c]) == sp.zeros(3)   # (5.8.43)（差號）
    for n1 in range(4):
        for n2 in range(3):
            w = su3_weights(n1, n2)
            assert sum(w.values()) == su3_dim(n1, n2)
    assert max(su3_weights(1, 1).values()) == 2 and max(su3_weights(2, 1).values()) == 2
    hw = (1 / math.sqrt(12), 1 / math.sqrt(36))                                                  # e1 的權 (H1, H2)
    r1v, r2v = np.array([float(x) for x in rt[1]]), np.array([float(x) for x in rt[2]])
    assert abs(2 * np.dot(hw, r1v) / np.dot(r1v, r1v) - 1) < 1e-12 and abs(2 * np.dot(hw, r2v) / np.dot(r2v, r2v)) < 1e-12   # (5.8.44)：𝒟^{10}
    for N in range(6):
        sh = osc_shell(N)
        assert sh["count"] == sh["lsum"] == (N + 1) * (N + 2) // 2 == su3_dim(N, 0) == len(su3_weights(N, 0)) if N == 0 else \
            sh["count"] == sh["lsum"] == (N + 1) * (N + 2) // 2 == su3_dim(N, 0) == sum(su3_weights(N, 0).values())
    assert osc_fock_levels(3) == [(1.5, 1), (2.5, 3), (3.5, 6), (4.5, 10)]
    assert COULOMB_CLOSURES == [2, 10, 28, 60] and OSC_CLOSURES == [2, 8, 20, 40, 70]


def _check_problems():
    p = problem_511()
    assert all(p.values())
    assert problem_512() == {"z": (0.0, 0.0, 0.0, 0.0), "x": (0.0, 0.0, 0.0, 0.0), "111": (0.0, 0.0, 3.0, 3.0)} or \
        all(abs(a - b) < 1e-9 for k in ("z", "x", "111") for a, b in zip(problem_512()[k], {"z": (0, 0, 0, 0), "x": (0, 0, 0, 0), "111": (0, 0, 3, 3)}[k]))
    q = problem_512(70)
    assert all(abs(v[0] - (1 + 2 * math.cos(math.radians(70)))) < 1e-6 and abs(v[1] + v[0]) < 1e-6 and abs(v[2] - v[3]) < 1e-6 for v in q.values())
    p3 = problem_513()
    assert p3["angles"] == [0, 60, 90, 120, 180]
    assert p3["groups"] == {"D2h": True, "D3h": True, "D4h": True, "D5h": False, "D6h": True, "D7h": False, "D2d": True, "D3d": True, "D4d": False}
    assert dinfh_to_dnh(2)[r"\Pi_u"] == "B_{1u}+B_{2u}" and dinfh_to_dnh(3)[r"\Sigma^+_u"] == "A_2''" and dinfh_to_dnh(6)[r"\Phi_g"] == "B_{1g}+B_{2g}"
    for l in range(13):
        assert problem_561()[l] == o3_to_d6_mrule(l)
        assert o3_to_y(l) == o3_to_y_mats(l)
        f = o3_to_y(l)
        assert f[0] + 3 * f[1] + 3 * f[2] + 4 * f[3] + 5 * f[4] == 2 * l + 1
        for par in (1, -1):
            assert o3_to_td(l, par) == o3_to_td_chars(l, par)
    want_y = [(1, 0, 0, 0, 0), (0, 1, 0, 0, 0), (0, 0, 0, 0, 1), (0, 0, 1, 1, 0), (0, 0, 0, 1, 1), (0, 1, 1, 0, 1), (1, 1, 0, 1, 1)]
    assert [o3_to_y(l) for l in range(7)] == want_y
    assert o3_to_td(1, -1) == (0, 0, 0, 0, 1) and o3_to_td(2, 1) == (0, 0, 1, 0, 1) and o3_to_td(3, -1) == (1, 0, 0, 1, 1)


def _check_cards():
    cv = card_values()
    for k, v in CARD_EXPECTED.items():
        assert cv[k] == v, (k, cv[k], v)


def summary() -> dict:
    return card_values()


def check():
    """每個上卡片／講義的數字用兩種以上獨立算法核對（約 15 秒）。"""
    _check_51_53()
    _check_54()
    _check_55()
    _check_56()
    _check_57()
    _check_58()
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
    print(f"sym_rot（第 5 章）：O3／O2 與 D∞h、Euler ↔ Darboux ↔ 四元數（Fig. 5.3.5 的 720° 序列）、生成元與譜分解、任意 j 的 J 與 𝒟^j"
          f"（Wigner 公式＝矩陣指數＝Schwinger 振子，含 (5.4.35) 的更正）、Stern–Gerlach 機率與古典錐、量子轉子（球形、對稱、不對稱、Zeeman、Stark）、"
          f"球諧與加法定理、O3↓O_h／D6 相關與助記輪、d 軌域、多極展開、SF6 能階叢、D6／D2／O 的射線代數、R4 的 Coulomb 對稱與 U3 振子殼層、"
          f"習題 5.1.1–5.6.2 全部核對通過（{time.time() - t0:.1f} 秒）。")
