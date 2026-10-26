# 由 combkit/tools/sync_from_harter_sym.py 從 harter_sym 專案的 codes/utils/sym_mat.py 複製；請改 harter_sym 那份再同步。
r"""
sym_mat.py —— Harter《Principles of Symmetry, Dynamics, and Spectroscopy》（Wiley 1993）第 1 章
A Review of Matrix Algebra and Quantum Mechanics 的數學（純 Python＋numpy＋sympy＋scipy，不 import manim）。
書上的式號、筆誤與〔注〕見 ref/md/ch01_formulas.md（相對於 C:\mypython\manim\projects\harter_sym）。

約定
  * 精確值用 sympy（Rational、sqrt、I）；浮點用 numpy complex。f-string 不能直接 :.4f 一個 sympy 數或 Fraction，要 float()。
  * 分析器角度 α（度）：分析器的「1」通道方向（在與束垂直的平面上）轉 α。書上 (1.1.3)：
        T(b'←b) = U(θ)，U(θ) = [[cos θ/2, −i sin θ/2], [−i sin θ/2, cos θ/2]]，θ = α_b − α_b'。
    所以 T(β←α) = U(α − β)；U(a)U(b) = U(a+b)，(1.1.6b) 自動成立。通過率 |⟨1'|1⟩|² = cos²(θ/2)。
  * 久期多項式一律取首一的 S(ε) = det(ε1 − M)（書上寫 det|M − ε1|，n 為奇數時差一個負號；見〔注〕）。
  * 書上的 2×2 例子 MB = [[4,1],[3,2]]：ε = 1、5；本系列把書上的 𝒫₂ 寫成 𝒫₅（依本徵值命名）。
    𝔭_j = 未歸一化的 Hamilton–Cayley 因子 (1.2.9)；𝒫_j = 冪等投影算符 (1.2.12)。ħ = 1。

§1.1A 分析器（H101–H102）
    U(θ)、U_num(θ)、T(β, α) = T(β←α)、T_real(θ)（正交的實特例 [[c, −s], [s, c]]）、pass_fraction(θ) = cos²(θ/2)（精確）
        U(60) = [[√3/2, −i/2], [−i/2, √3/2]]；|⟨1'|1⟩|² = 3/4、|⟨2'|1⟩|² = 1/4
        PASS_ANGLES 0、30、60、90、120、150、180° → 1.0000、0.9330、0.7500、0.5000、0.2500、0.0670、0.0000
    analyzer_chain(stages, start=(角度, 通道)) → (p1, p2)（密度矩陣 ρ → TρT†；mode 1/2＝只留該通道、"coherent"＝相干合併、"which"＝機率相加）
    axiom4(θ_in, θ_mid, θ_out, k) → (|Σ振幅|², Σ|振幅|²)
    series_table()（SERIES 0° → 90° → 180°，入射 |1_0°⟩，180° 分析器的 (p1, p2)）：
        direct (0, 1)、coherent (0, 1)（Axiom 4：合併後像沒分過）、which (1/2, 1/2)、stop2（擋掉 90° 的 2 通道）(1/4, 1/4)
        干涉項 = 相干 − 不相干 = −1/2（k = 1）
    problem_111(θ)（習題 1.1.1：90° 分開再合併 → 0° 計數；stopper 擋 2''）＝ problem_111_closed(θ)：
        入射 |1⟩：1、1/4；|2⟩：0、1/4；|1'⟩：cos²(θ/2)、(1 + sin θ)/4 —— θ = 60°：3/4 = 0.7500、(2+√3)/8 = 0.4665
    analyzer_chain_num(stages, start)（浮點版，給動畫）；zeno(n) = cos^{2n}(π/4n) ＝ zeno_chain(n)（習題 1.1.2）：n = 1、2、3、4、10 → 0.5000、0.7286、0.8122、0.8562、0.9401（→ 1）
§1.1B–C bra、ket、基本算子、換基底（H103）
    SX、SY、SZ、E11、E22、basis(k, n)、elementary(i, j, n) = |i⟩⟨j|、ket_primed(θ, j)（(1.1.13)：(√3/2, −i/2) at 60°）
    change_basis(H, T) = THT†、outer(ψ, φ) = |ψ⟩⟨φ|、dag(M)、is_unitary(M)、to_np(M)
        U(60)E11U(60)† = [[3/4, i√3/4], [−i√3/4, 1/4]]（= |1⟩⟨1| 在 60° 基底，也是 H108 的純態 ρ）
        U(60)σzU(60)† = cos θ σz − sin θ σy = [[1/2, i√3/2], [−i√3/2, −1/2]]：本徵值 ±1、跡 0、行列式 −1 不變
        e₁₂|2⟩ = |1⟩、⟨1|e₁₂ = ⟨2|、e₁₂|1⟩ = 0、⟨2|e₁₂ = 0（但 ⟨1|e₁₂ ≠ 0：(1.1.19c) 的條件應是 k ≠ i）
§1.1C 兩態拍頻（H104）：BEAT_AB = (A, B) = (3, 1)、H2(A, B) = [[A, B], [B, A]]：本徵值 4、2，本徵向量 (1, ±1)/√2
    populations(t) = (cos²Bt, sin²Bt) ＝ populations_numeric(ts)（solve_ivp）＝ |psi_t(t)|²（投影算符）＝ scipy expm
        t = 1：0.2919、0.7081；t = π/4：0.5000、0.5000；t = π/2 = 1.5708：0.0000、1.0000；週期 π = 3.1416；ΔE = 2
§1.2 對角化（H105–H107）：MB、NB = [[5, −1], [−3, 7]]、DEG4（反對角 1 的 4×4）、JORDAN = [[2,1,0],[0,2,0],[0,0,3]]
    secular(M)（Poly，首一）、eigen_mult(M)、distinct_eigs(M)、minimal_poly(M)、poly_at(poly, M)
        MB：ε² − 6ε + 5，ε = 1、5；MB² = [[19, 6], [18, 7]]，MB² − 6MB + 5 = 0
    hc_factors(M) = {ε: 𝔭}：𝔭₁ = [[−1, 1], [3, −3]]、𝔭₅ = [[3, 1], [3, 1]]；𝔭₁² = −4𝔭₁、𝔭₅² = 4𝔭₅、𝔭₁𝔭₅ = 0
    projectors(M) = {ε: 𝒫}：𝒫₁ = ¼[[1, −1], [−3, 3]]、𝒫₅ = ¼[[3, 1], [3, 1]]；冪等、互相正交、和 = 1；MB = 𝒫₁ + 5𝒫₅
    RIGHT_VEC {1: (1, −3), 5: (1, 1)}、LEFT_VEC {1: (1, −1), 5: (3, 1)}；eig_book(M) = {ε: (bra, ket)}（(1.2.13)）：
        ⟨e₁| = (1/2, −1/2)、|e₁⟩ = (1/2, −3/2)；⟨e₂| = (√3/2, √3/6)、|e₂⟩ = (√3/2, √3/2)；方向角 −45.00°、−71.57°、18.43°、45.00°；
        兩個 ket 夾 116.57°（不垂直），但 ⟨e₁|⊥|e₂⟩、⟨e₂|⊥|e₁⟩
    angle_of(v)、ellipse_points(M, n)、singular_values(M) = 5.3983、0.9262（乘積 5 = det）
    iterate_angles(MB, (1,0), 4) = 0.00、36.87、43.45、44.69、44.94°（→ 45°）
    mat_func(M, f) = Σ f(ε)𝒫：MB⁵⁰ = 𝒫₁ + 5⁵⁰𝒫₅ = ¼[[1 + 3·5⁵⁰, 5⁵⁰ − 1], [3·5⁵⁰ − 3, 3 + 5⁵⁰]]（(1,1) 元 ≈ 6.6613e+34）；
        √MB = 𝒫₁ + √5𝒫₅ ≈ [[1.9271, 0.3090], [0.9271, 1.3090]]（共 4 個 ±𝒫₁ ± √5𝒫₅）；
        e^{0.1MB} = e^{0.1}𝒫₁ + e^{0.5}𝒫₅ ≈ [[1.5128, 0.1359], [0.4077, 1.2411]]
    lagrange_basis(xs, x)：Σℓ_j = 1、Σx_jℓ_j = x
    DEG4：S = ε⁴ − 2ε² + 1 = (ε−1)²(ε+1)²、最小方程 ε² − 1；𝒫_{±1} = (1 ± DEG4)/2（秩 2、2）
    NB 與 MB 對易，NB = 8𝒫₁ + 4𝒫₅；joint_projectors(A, B) = {(α, β): 𝒫_α(A)𝒫_β(B)}：(1, 8) → 𝒫₁、(5, 4) → 𝒫₅
    adj_vectors(M, ε) = adj(M − ε1)：ε = 1 → [[1, −1], [−3, 3]] = 4𝒫₁；ε = 5 → [[−3, −1], [−3, −1]] = −4𝒫₅（一般 adj(ε_j1 − M) = 𝔭_j）
§1.3 散射與密度矩陣（H108）
    S_EX = (1/5)[[3, 4i], [4i, 3]]：么正、|S₁₁|² = 9/25 = 0.36、|S₂₁|² = 16/25 = 0.64、本徵值 (3 ± 4i)/5（相位 ±53.13°）、本徵向量 (1, ±1)/√2
    flux_in_out(S, A_in)、rho_pure(ψ)、purity(ρ)：純態 tr ρ² = 1；混合 diag(3/4, 1/4)：5/8 = 0.625；1/2·1：1/2
§1.4 耦合振子（Fig. 1.4.1，H108）：OSC = {"equal": (1,1,1,1,1), "beat": (1,1,1,1,1/4), "unequal": (1,2,1,1,1)}（m1, m2, k1, k2, k12）
    force_matrix、mass_matrix、accel_matrix（a = m⁻¹𝓕）、omega2(case)（精確）、omega2_generalized(case) → (ω², E)（m-歸一化）
        equal：ω² = 1、3（ω = 1.0000、1.7321；模態 (1, 1)、(1, −1)）；beat：ω² = 1、3/2（ω = 1.0000、1.2247）；
        unequal：a = [[2, −1], [−1/2, 1]]，ω² = (3 ∓ √3)/2 = 0.6340、2.3660（ω = 0.7962、1.5382）
    osc_analytic(t, case, x0, v0) ＝ osc_numeric(ts, …)；osc_energy；beat_times("beat") = (13.98, 27.96)（能量全到 m₂、拍的週期）
    modes_symbolic() → (a, [ω²₋, ω²₊], 符號)；spring_force_matrix(springs, n)（Fig. 1.4.2 的規則＝位能的 Hessian）
§1.5 連續空間與預解式（H109）
    d1_matrix(N, L)（反對稱）、d2_matrix、fourier_matrix(N)（么正，F†D₁F 對角）、d1_eigen(N, k)：N = 8：0.9003i、N = 32：0.9936i（k = 1，精確 i）
    delta_gauss(x, σ)、DELTA_SIGMAS (1, 0.5, 0.1)、delta_integral(σ) = ∫δ_σ cos = e^{−σ²/2} = 0.6065、0.8825、0.9950
    delta_deriv_integrals(σ) = (∫f ∂_yδ, ∫f ∂_xδ)，f = sin：σ = 0.1：−0.9950、+0.9950（(1.5.10) 第二式的符號）
    adjoint_fd(μ, ν, λ, N) → (M†, 公式 (1.5.12))；resolvent(M, λ)；resolvent_adj(M) = (adj(λ1 − M), S(λ)) = ([[λ−2, 1], [3, λ−4]], λ² − 6λ + 5)
        resolvent(MB, 3) = −¼[[1, 1], [3, −1]] = 𝒫₁/2 − 𝒫₅/2
    contour(M, 圓心, 半徑, f, n=64)：(1/2πi)∮f(λ)(λ − M)⁻¹dλ —— 繞 1 → 𝒫₁、繞 5 → 𝒫₅（誤差 < 1e−12）、兩個都繞 → 1、f = λ → M
    contour_PN(M, λ_a, r)：JORDAN 在 2：𝒫 = diag(1, 1, 0)、𝒩 = e₁₂（𝒩² = 0）；M = Σ(λ_a𝒫_a + 𝒩_a)
附錄 B：inverse2(M, book=False)（書上左下分母：MB 給 −3/2，應為 −3/5）、cramer(M, y)：MB、y = (1, 0) → (2/5, −3/5)
習題：P121 → ε = 4、−1，𝒫₄ = (1/5)[[3, −1], [−6, 2]]、𝒫₋₁ = (1/5)[[2, 1], [6, 3]]，sqrt_all(M) 4 個；P122 → 4、9、36，√ = [[3,1,−1],[1,3,−1],[−1,−1,5]]；
    secular_coeffs(M) 三種算法（P122：−49、504、−1296）；P124（ρ = 1）：ε⁴ − 4ε³、MEq ε² − 4ε；problem_125() = (1/3)[[1,−2,1],[−2,1,1],[1,1,−2]]，
    M¹⁰⁰ = 1 − J/3；P126_M、P126_N 的聯合投影 (1,4)、(3,1)、(3,4)；P127_H、P127_A、P127_U（本徵值 ±1、±i/√2、e^{iπ/12}、e^{−7iπ/12}）；
    shift_matrix(n)（習題 1.2.8）
check()：每個上卡片／講義的數字用兩種以上算法核對（約 2 秒）；summary() → {名稱: 字串}。
直接執行：python utils/sym_mat.py（在 codes/ 裡）印一句成功訊息；--summary 印出全部數字。
"""
import math
from itertools import combinations

import numpy as np
import sympy as sp
from scipy import integrate, linalg

I = sp.I
EPS = sp.Symbol("epsilon")
LAM = sp.Symbol("lambda")
HBAR = 1


# =============================================================================
# 小工具
# =============================================================================
def _f4(x) -> str:
    return f"{float(x):.4f}"


def to_np(M) -> np.ndarray:
    """sympy 矩陣 → numpy complex 陣列。"""
    return np.array([[complex(sp.N(x, 30)) for x in row] for row in sp.Matrix(M).tolist()], dtype=complex)


def dag(M):
    """轉置共軛 M†（(1.1.7)）；sympy 或 numpy 都可以。"""
    return M.H if isinstance(M, sp.MatrixBase) else np.conj(M).T


def _zero(M) -> bool:
    return sp.simplify(M) == sp.zeros(*M.shape)


def is_unitary(M, tol: float = 1e-12) -> bool:
    """(1.1.8c)：M†M = 1 = MM†。"""
    n = M.shape[0]
    if isinstance(M, sp.MatrixBase):
        return _zero(M.H * M - sp.eye(n)) and _zero(M * M.H - sp.eye(n))
    return np.allclose(dag(M) @ M, np.eye(n), atol=tol) and np.allclose(M @ dag(M), np.eye(n), atol=tol)


def _close(a, b, tol: float = 1e-10) -> bool:
    return np.allclose(np.asarray(a, dtype=complex), np.asarray(b, dtype=complex), atol=tol, rtol=0)


# =============================================================================
# §1.1A 分析器與轉換矩陣（H101–H102）
# =============================================================================
def _half(theta):
    """θ/2（θ 以度計）→ 弧度，精確。"""
    return sp.pi * sp.nsimplify(theta) / 360


def U(theta) -> sp.Matrix:
    """(1.1.3)：U(θ) = [[cos θ/2, −i sin θ/2], [−i sin θ/2, cos θ/2]]，θ 以度計（精確）。U(a)U(b) = U(a+b)。"""
    h = _half(theta)
    return sp.Matrix([[sp.cos(h), -I * sp.sin(h)], [-I * sp.sin(h), sp.cos(h)]])


def U_num(theta: float) -> np.ndarray:
    h = math.radians(theta) / 2
    return np.array([[math.cos(h), -1j * math.sin(h)], [-1j * math.sin(h), math.cos(h)]])


def T(beta, alpha) -> sp.Matrix:
    """T(β←α)_{ij} = ⟨i_β|j_α⟩：從角度 α 的分析器基底換到角度 β 的（書上 T(b'←b) = U(θ)，θ = α_b − α_b'）。"""
    return U(sp.nsimplify(alpha) - sp.nsimplify(beta))


def T_real(theta) -> sp.Matrix:
    """實的特例（正交矩陣，(1.1.10)）：[[cos θ/2, −sin θ/2], [sin θ/2, cos θ/2]]（繞另一個軸轉；通過率相同）。"""
    h = _half(theta)
    return sp.Matrix([[sp.cos(h), -sp.sin(h)], [sp.sin(h), sp.cos(h)]])


def pass_fraction(theta):
    """|⟨1'|1⟩|² = cos²(θ/2)（精確）。"""
    return sp.nsimplify(sp.simplify(sp.cos(_half(theta)) ** 2))


PASS_ANGLES = (0, 30, 60, 90, 120, 150, 180)
E11 = sp.Matrix([[1, 0], [0, 0]])
E22 = sp.Matrix([[0, 0], [0, 1]])


def analyzer_chain(stages, start=(0, 1)):
    """
    一串分析器。start = (角度, 通道)：入射束全在該分析器的態 |通道⟩。stages = [(角度, mode), …]，最後一個是計數用的分析器。
    mode：1 或 2 = 只留該通道（另一束被 stopper 擋掉）；"coherent" = 兩束都留、再相干合併（Axiom 4，振幅相加）；
    "which" = 兩束都留但有人量了走哪一條（相位亂掉，機率相加）。用密度矩陣算：ρ → TρT†。回傳最後一個分析器兩個出口的機率 (p1, p2)。
    """
    ang, ch = start
    psi = sp.Matrix([1, 0]) if ch == 1 else sp.Matrix([0, 1])
    rho = psi * psi.H
    for k, (a, mode) in enumerate(stages):
        Tm = T(a, ang)
        rho = Tm * rho * Tm.H
        ang = a
        if k == len(stages) - 1:
            break
        if mode == 1:
            rho = E11 * rho * E11
        elif mode == 2:
            rho = E22 * rho * E22
        elif mode == "which":
            rho = sp.diag(rho[0, 0], rho[1, 1])
        elif mode != "coherent":
            raise ValueError(mode)
    return tuple(sp.nsimplify(sp.simplify(sp.expand(rho[i, i].rewrite(sp.cos)))) for i in (0, 1))


def axiom4(theta_in, theta_mid, theta_out, k: int = 1):
    """Axiom 4：⟨k''|1'⟩ = Σ_j ⟨k''|j⟩⟨j|1'⟩。回傳 (|Σ 振幅|², Σ|振幅|²) = (相干, 不相干)；兩者之差＝干涉項。"""
    A, B = T(theta_mid, theta_in), T(theta_out, theta_mid)
    terms = [B[k - 1, j] * A[j, 0] for j in range(2)]
    s = sum(terms)
    coh = sp.nsimplify(sp.simplify(sp.expand(s * sp.conjugate(s))))
    inc = sp.nsimplify(sp.simplify(sp.expand(sum(t * sp.conjugate(t) for t in terms))))
    return coh, inc


SERIES = (0, 90, 180)       # H101：0° → 90° → 180°


def series_table():
    """0° 的 |1⟩ 進 90° 分析器、再進 180° 分析器：四種情形的 (p1, p2)。"""
    return {
        "direct": analyzer_chain([(180, "count")]),
        "coherent": analyzer_chain([(90, "coherent"), (180, "count")]),
        "which": analyzer_chain([(90, "which"), (180, "count")]),
        "stop2": analyzer_chain([(90, 1), (180, "count")]),
    }


def problem_111(theta=60):
    """習題 1.1.1：先分成 right/left（90°）、再合併、最後 up/dn（0°）計數；stopper 擋 2''。回傳 {入射態: (無 stopper, 有 stopper)}。"""
    out = {}
    for name, start in (("1", (0, 1)), ("2", (0, 2)), ("1'", (theta, 1))):
        out[name] = (analyzer_chain([(90, "coherent"), (0, "count")], start)[0],
                     analyzer_chain([(90, 1), (0, "count")], start)[0])
    return out


def problem_111_closed(theta):
    t = math.radians(theta)
    return {"1": (1.0, 0.25), "2": (0.0, 0.25), "1'": (math.cos(t / 2) ** 2, (1 + math.sin(t)) / 4)}


def zeno(n: int) -> float:
    """習題 1.1.2：n 個分析器依次轉 (k/n)·90°，每個只留通道 1：最後是 right 的比例 cos^{2n}(π/4n)（n → ∞ 趨於 1）。"""
    return math.cos(math.pi / (4 * n)) ** (2 * n)


def analyzer_chain_num(stages, start=(0.0, 1)) -> tuple:
    """analyzer_chain 的浮點版（numpy，角度可為任意浮點；給場景的動畫用）。"""
    ang, ch = start
    psi = np.array([1, 0], dtype=complex) if ch == 1 else np.array([0, 1], dtype=complex)
    rho = np.outer(psi, psi.conj())
    for k, (a, mode) in enumerate(stages):
        Tm = U_num(ang - a)
        rho = Tm @ rho @ dag(Tm)
        ang = a
        if k == len(stages) - 1:
            break
        if mode in (1, 2):
            Pm = np.zeros((2, 2))
            Pm[mode - 1, mode - 1] = 1
            rho = Pm @ rho @ Pm
        elif mode == "which":
            rho = np.diag(np.diag(rho))
        elif mode != "coherent":
            raise ValueError(mode)
    return float(rho[0, 0].real), float(rho[1, 1].real)


def zeno_chain(n: int) -> float:
    stages = [(90 * k / n, 1) for k in range(1, n)] + [(90, "count")]
    return analyzer_chain_num(stages)[0]


# =============================================================================
# §1.1B–C bra、ket、基本算子、換基底（H103）
# =============================================================================
SX = sp.Matrix([[0, 1], [1, 0]])
SY = sp.Matrix([[0, -I], [I, 0]])
SZ = sp.Matrix([[1, 0], [0, -1]])


def basis(k: int, n: int = 2) -> sp.Matrix:
    v = sp.zeros(n, 1)
    v[k - 1] = 1
    return v


def elementary(i: int, j: int, n: int = 2) -> sp.Matrix:
    """(1.1.21) e_ij = |i⟩⟨j|。"""
    return basis(i, n) * basis(j, n).T


def ket_primed(theta, j: int = 1) -> sp.Matrix:
    """(1.1.13)：|j⟩ 在 θ 轉過的（primed）基底下的分量 = U(θ) 的第 j 行；j = 1：(cos θ/2, −i sin θ/2)。"""
    return U(theta)[:, j - 1]


def change_basis(H, Tm):
    """(1.1.23c) H' = T H T†。"""
    return sp.simplify(Tm * H * Tm.H) if isinstance(H, sp.MatrixBase) else Tm @ H @ dag(Tm)


def outer(psi, phi):
    """(1.1.22) |Ψ⟩⟨Φ| 的矩陣 = Ψ ⊗ Φ†。"""
    return psi * phi.H


# =============================================================================
# §1.1C、§1.2 開頭：兩態 Hamilton 與拍頻（H104）
# =============================================================================
BEAT_AB = (3, 1)            # H = [[A, B], [B, A]]，ħ = 1


def H2(A=3, B=1) -> sp.Matrix:
    return sp.Matrix([[A, B], [B, A]])


def populations(t: float, A: float = 3, B: float = 1):
    """|Ψ(0)⟩ = |1⟩：P₁(t) = cos²(Bt/ħ)、P₂ = sin²(Bt/ħ)（閉式）。"""
    return math.cos(B * t / HBAR) ** 2, math.sin(B * t / HBAR) ** 2


def psi_t(t: float, A: float = 3, B: float = 1, psi0=(1, 0)) -> np.ndarray:
    """e^{−iHt/ħ} 用投影算符：e^{−i(A+B)t}P₊ + e^{−i(A−B)t}P₋，P± = (1 ± σx)/2（(1.2.2)：本徵態只轉相位）。"""
    Pp = np.array([[1, 1], [1, 1]]) / 2
    Pm = np.array([[1, -1], [-1, 1]]) / 2
    Ut = np.exp(-1j * (A + B) * t / HBAR) * Pp + np.exp(-1j * (A - B) * t / HBAR) * Pm
    return Ut @ np.array(psi0, dtype=complex)


def populations_numeric(ts, A: float = 3, B: float = 1):
    """數值積分 (1.1.17) iħ dΨ/dt = HΨ（solve_ivp，實部虛部分開）。回傳 (P₁ 陣列, P₂ 陣列)。"""
    Hm = np.array([[A, B], [B, A]], dtype=float)

    def rhs(_t, y):
        psi = y[:2] + 1j * y[2:]
        d = -1j * (Hm @ psi) / HBAR
        return np.concatenate([d.real, d.imag])

    sol = integrate.solve_ivp(rhs, (0, max(ts)), [1, 0, 0, 0], t_eval=ts, rtol=1e-11, atol=1e-12)
    psi = sol.y[:2] + 1j * sol.y[2:]
    return np.abs(psi[0]) ** 2, np.abs(psi[1]) ** 2


# =============================================================================
# §1.2 對角化：久期方程、Hamilton–Cayley、投影算符（H105–H107）
# =============================================================================
MB = sp.Matrix([[4, 1], [3, 2]])          # 書上的例子（(1.2.4b)_x 起）
NB = sp.Matrix([[5, -1], [-3, 7]])        # 與 MB 對易（p. 27）：N = 8P₁ + 4P₅
DEG4 = sp.Matrix([[0, 0, 0, 1], [0, 0, 1, 0], [0, 1, 0, 0], [1, 0, 0, 0]])   # (1.2.16) 的重根例子
JORDAN = sp.Matrix([[2, 1, 0], [0, 2, 0], [0, 0, 3]])                        # 不可對角化（書外，(1.5.20)–(1.5.21) 用）


def secular(M) -> sp.Poly:
    """久期多項式（首一）S(ε) = det(ε1 − M)；MB：ε² − 6ε + 5。"""
    n = M.shape[0]
    return sp.Poly((EPS * sp.eye(n) - M).det(), EPS)


def eigen_mult(M) -> dict:
    """{本徵值: 重數}（由久期多項式的根）。"""
    return sp.roots(secular(M))


def distinct_eigs(M) -> list:
    return sorted(eigen_mult(M), key=lambda z: (sp.re(z), sp.im(z)))


def hc_factors(M) -> dict:
    """(1.2.9) 未歸一化的 𝔭_j = Π_{l≠j}(M − ε_l 1)（相異本徵值）；MB：𝔭₁ = [[−1, 1], [3, −3]]、𝔭₅ = [[3, 1], [3, 1]]。"""
    n = M.shape[0]
    ev = distinct_eigs(M)
    out = {}
    for e in ev:
        p = sp.eye(n)
        for l in ev:
            if l != e:
                p = p * (M - l * sp.eye(n))
        out[e] = sp.simplify(p)
    return out


def projectors(M) -> dict:
    """(1.2.12)/(1.2.15) 𝒫_j = Π_{l≠j}(M − ε_l)/Π_{l≠j}(ε_j − ε_l)（對可對角化的矩陣；重根時就是用最小方程）。"""
    ev = distinct_eigs(M)
    out = {}
    for e, p in hc_factors(M).items():
        d = sp.Integer(1)
        for l in ev:
            if l != e:
                d *= (e - l)
        out[e] = sp.simplify(p / d)
    return out


def minimal_poly(M) -> sp.Poly:
    """最小方程（MEq）：每個本徵值的冪次＝使 nullity((M − ε)^k) 等於代數重數的最小 k。"""
    n = M.shape[0]
    expr = sp.Integer(1)
    for e, m in eigen_mult(M).items():
        A = M - e * sp.eye(n)
        k = 1
        while n - (A ** k).rank() < m:
            k += 1
        expr *= (EPS - e) ** k
    return sp.Poly(sp.expand(expr), EPS)


def poly_at(poly: sp.Poly, M) -> sp.Matrix:
    """把多項式的 ε 換成矩陣（常數項乘 1）——Hamilton–Cayley (1.2.6)。"""
    n = M.shape[0]
    out = sp.zeros(n)
    for (k,), c in poly.terms():
        out += c * M ** k
    return sp.simplify(out)


def mat_func(M, f) -> sp.Matrix:
    """譜分解 f(M) = Σ f(ε_j) 𝒫_j（f 吃 sympy 數）。"""
    out = sp.zeros(M.shape[0])
    for e, P in projectors(M).items():
        out += f(e) * P
    return sp.simplify(out)


def eig_book(M) -> dict:
    """(1.2.13)：從 𝒫_j 取第一個非零對角元 N_j 所在的列與行，bra = 列/√N_j、ket = 行/√N_j。回傳 {ε: (bra 1×n, ket n×1)}。"""
    out = {}
    for e, P in projectors(M).items():
        a = next(i for i in range(P.shape[0]) if P[i, i] != 0)
        Nj = P[a, a]
        out[e] = (sp.simplify(P[a, :] / sp.sqrt(Nj)), sp.simplify(P[:, a] / sp.sqrt(Nj)))
    return out


RIGHT_VEC = {1: (1, -3), 5: (1, 1)}       # MB 的右（ket）本徵向量（整數形）
LEFT_VEC = {1: (1, -1), 5: (3, 1)}        # 左（bra）本徵向量


def adj_vectors(M, e) -> sp.Matrix:
    """(1.2.25)：adj(M − ε_j 1) 的行是 ket、列是 bra 本徵向量；一般 adj(ε_j 1 − M) = 𝔭_j。"""
    return (M - e * sp.eye(M.shape[0])).adjugate()


def joint_projectors(A, B) -> dict:
    """(1.2.24)：對易的 A、B：𝒫^{α,β} = 𝒫_α(A)𝒫_β(B)，只回傳非零的 {(α, β): 矩陣}。"""
    out = {}
    for a, Pa in projectors(A).items():
        for b, Pb in projectors(B).items():
            P = sp.simplify(Pa * Pb)
            if P != sp.zeros(*P.shape):
                out[(a, b)] = P
    return out


def lagrange_basis(xs, x):
    """Lagrange 內插的係數 ℓ_j(x) = Π_{k≠j}(x − x_k)/(x_j − x_k)（(1.2.12) 的「數」版本）。"""
    out = []
    for j, xj in enumerate(xs):
        v = sp.Integer(1)
        for k, xk in enumerate(xs):
            if k != j:
                v *= (x - xk) / sp.nsimplify(xj - xk)
        out.append(v)
    return out


def ellipse_points(M, n: int = 60) -> np.ndarray:
    """單位圓經 M 變成的橢圓（n 點，2×n 實陣列）。"""
    t = np.linspace(0, 2 * np.pi, n, endpoint=False)
    return np.array(to_np(M).real) @ np.vstack([np.cos(t), np.sin(t)])


def singular_values(M):
    return tuple(np.linalg.svd(to_np(M).real, compute_uv=False))


def iterate_angles(M, v=(1, 0), k: int = 4):
    """M^m v 的方向角（度），m = 0..k：MB、v = (1,0)：0、36.87、43.45、44.69、44.94 → 45°（ε = 5 的方向）。"""
    A = to_np(M).real
    w = np.array(v, dtype=float)
    out = []
    for _ in range(k + 1):
        out.append(math.degrees(math.atan2(w[1], w[0])))
        w = A @ w
    return out


def angle_of(v) -> float:
    return math.degrees(math.atan2(float(v[1]), float(v[0])))


# =============================================================================
# §1.3 散射矩陣、密度矩陣（H108 前半）
# =============================================================================
S_EX = sp.Matrix([[3, 4 * I], [4 * I, 3]]) / 5      # 書外的么正 S 矩陣例子：反射 9/25、穿透 16/25


def flux_in_out(S, A_in) -> tuple:
    """(1.3.2)：Σ|A^I|² 與 Σ|A^o|²（A^o = S A^I）。"""
    A_in = np.asarray(A_in, dtype=complex)
    A_out = to_np(S) @ A_in
    return float(np.vdot(A_in, A_in).real), float(np.vdot(A_out, A_out).real)


def rho_pure(psi) -> sp.Matrix:
    """(1.3.6) ρ = |Ψ⟩⟨Ψ|。"""
    return sp.simplify(psi * psi.H)


def purity(rho):
    """tr ρ²：純態 1、混合態 < 1。"""
    return sp.nsimplify(sp.simplify((rho * rho).trace()))


# =============================================================================
# §1.4 耦合振子（Fig. 1.4.1，H108）
# =============================================================================
OSC = {                    # (m1, m2, k1, k2, k12)
    "equal": (1, 1, 1, 1, 1),                       # 等質量、等彈簧：ω² = k/m、3k/m
    "beat": (1, 1, 1, 1, sp.Rational(1, 4)),        # 弱耦合：ω² = 1、3/2，拍頻
    "unequal": (1, 2, 1, 1, 1),                     # m₂ = 2m₁：a = m⁻¹F 不對稱
}


def force_matrix(k1, k2, k12) -> sp.Matrix:
    """(1.4.2) 力矩陣 𝓕（F = −𝓕χ）。"""
    return sp.Matrix([[k1 + k12, -k12], [-k12, k2 + k12]])


def mass_matrix(m1, m2) -> sp.Matrix:
    """(1.4.3) 質量矩陣 m。"""
    return sp.diag(m1, m2)


def accel_matrix(m1, m2, k1, k2, k12) -> sp.Matrix:
    """(1.4.4) 加速度矩陣 a = m⁻¹𝓕（χ̈ = −aχ；書上內文寫 a = −m⁻¹𝓕 是筆誤）。"""
    return sp.simplify(mass_matrix(m1, m2).inv() * force_matrix(k1, k2, k12))


def omega2(case: str = "equal") -> list:
    """a 的本徵值 ω²（由小到大，精確）。"""
    return distinct_eigs(accel_matrix(*OSC[case]))


def omega2_generalized(case: str = "equal"):
    """廣義本徵問題 (1.4.8b)–(1.4.13)：𝓕e = ω² m e，scipy.linalg.eigh(F, m) → (ω², E)，E 的行 m-歸一化（EᵀmE = 1）。"""
    m1, m2, k1, k2, k12 = (float(x) for x in OSC[case])
    F = np.array([[k1 + k12, -k12], [-k12, k2 + k12]])
    m = np.diag([m1, m2])
    w2, E = linalg.eigh(F, m)
    return w2, E


def osc_analytic(t, case: str = "beat", x0=(1.0, 0.0), v0=(0.0, 0.0)) -> np.ndarray:
    """用簡正模解 χ(t)：q_j = ⟨e_j|m|x⟩（(1.4.19b)），每個模態以 ω_j 振盪。t 可為陣列；回傳 2×len(t)。"""
    m1, m2 = (float(x) for x in OSC[case][:2])
    m = np.diag([m1, m2])
    w2, E = omega2_generalized(case)
    w = np.sqrt(w2)
    q0, p0 = E.T @ m @ np.array(x0), E.T @ m @ np.array(v0)
    t = np.atleast_1d(t)
    q = q0[:, None] * np.cos(np.outer(w, t)) + (p0 / w)[:, None] * np.sin(np.outer(w, t))
    return E @ q


def osc_numeric(ts, case: str = "beat", x0=(1.0, 0.0), v0=(0.0, 0.0)) -> np.ndarray:
    """直接數值積分 χ̈ = −aχ（solve_ivp）。"""
    a = to_np(accel_matrix(*OSC[case])).real
    sol = integrate.solve_ivp(lambda _t, y: np.concatenate([y[2:], -a @ y[:2]]), (0, max(ts)),
                              list(x0) + list(v0), t_eval=ts, rtol=1e-11, atol=1e-12)
    return sol.y[:2]


def osc_energy(x, v, case: str = "beat") -> float:
    """H = T + V = ½⟨ẋ|m|ẋ⟩ + ½⟨x|𝓕|x⟩（(1.4.16)、(1.4.19)）。"""
    m1, m2, k1, k2, k12 = (float(x_) for x_ in OSC[case])
    F = np.array([[k1 + k12, -k12], [-k12, k2 + k12]])
    return 0.5 * (m1 * v[0] ** 2 + m2 * v[1] ** 2) + 0.5 * float(np.asarray(x) @ F @ np.asarray(x))


def beat_times(case: str = "beat"):
    """(能量全部到 m₂ 的時間 π/(ω₂ − ω₁), 拍的週期 2π/(ω₂ − ω₁))。"""
    w = [math.sqrt(float(x)) for x in omega2(case)]
    d = w[1] - w[0]
    return math.pi / d, 2 * math.pi / d


def modes_symbolic():
    """書上的一般例子：m₁、m₂、k₁、k₂、k₁₂ 符號。回傳 (a 矩陣, [ω²₋, ω²₊])。"""
    m1, m2, k1, k2, k12 = sp.symbols("m_1 m_2 k_1 k_2 k_12", positive=True)
    a = accel_matrix(m1, m2, k1, k2, k12)
    tr, det = a.trace(), sp.simplify(a.det())
    disc = sp.sqrt(sp.simplify(tr ** 2 - 4 * det))
    return a, [sp.simplify((tr - disc) / 2), sp.simplify((tr + disc) / 2)], (m1, m2, k1, k2, k12)


def spring_force_matrix(springs, n: int) -> np.ndarray:
    """Fig. 1.4.2 的「目視」規則：⟨i|F|j⟩ = Σ_彈簧 k·c_i·c_j（c_i = 彈簧方向在座標軸 i 上的投影，有正負）。
    springs = [(k, {軸 i: c_i}), …]。等於位能 V = ½Σk(Σc_iχ_i)² 的 Hessian（(1.4.18)）。"""
    F = np.zeros((n, n))
    for k, c in springs:
        v = np.zeros(n)
        for i, ci in c.items():
            v[i] = ci
        F += k * np.outer(v, v)
    return F


# =============================================================================
# §1.5 連續向量空間、微分算子、預解式（H109）
# =============================================================================
def d1_matrix(N: int, L: float = 2 * math.pi) -> np.ndarray:
    """週期邊界的中央差分 d/dx（N×N）：反對稱 ⇒ 反 Hermitian；本徵向量是平面波，本徵值 i·sin(kh)/h。"""
    h = L / N
    D = np.zeros((N, N))
    for j in range(N):
        D[j, (j + 1) % N] = 1 / (2 * h)
        D[j, (j - 1) % N] = -1 / (2 * h)
    return D


def d2_matrix(N: int, L: float = 2 * math.pi) -> np.ndarray:
    h = L / N
    D = np.zeros((N, N))
    for j in range(N):
        D[j, j] = -2 / h ** 2
        D[j, (j + 1) % N] = D[j, (j - 1) % N] = 1 / h ** 2
    return D


def fourier_matrix(N: int) -> np.ndarray:
    """⟨x_j|k⟩ = e^{2πijk/N}/√N（(1.5.6) 的離散版）：么正，Fourier 就是換基底。"""
    j = np.arange(N)
    return np.exp(2j * np.pi * np.outer(j, j) / N) / math.sqrt(N)


def d1_eigen(N: int, k: int = 1, L: float = 2 * math.pi) -> complex:
    h = L / N
    return 1j * math.sin(2 * math.pi * k / L * h) / h


def delta_gauss(x, s: float):
    """δ 的高斯序列 δ_σ(x) = e^{−x²/2σ²}/(√(2π)σ)。"""
    return np.exp(-np.asarray(x) ** 2 / (2 * s * s)) / (math.sqrt(2 * math.pi) * s)


DELTA_SIGMAS = (1.0, 0.5, 0.1)


def delta_integral(s: float) -> tuple:
    """∫δ_σ(x) cos x dx：(數值積分, 閉式 e^{−σ²/2}) → σ → 0 時趨於 cos 0 = 1（(1.5.3b)）。"""
    num = integrate.quad(lambda x: delta_gauss(x, s) * math.cos(x), -12 * s, 12 * s, limit=200)[0]
    return num, math.exp(-s * s / 2)


def delta_deriv_integrals(s: float) -> tuple:
    """(1.5.10) 的符號：f = sin、x = 0。回傳 (∫dy f(y) ∂_y δ_σ(x − y), ∫dy f(y) ∂_x δ_σ(x − y))
    = (−e^{−σ²/2}, +e^{−σ²/2}) → (−f'(x), +f'(x))。書上第二式 ∫dy f(y) dδ(x,y)/dy = f'(x) 少了負號（應對 x 微分）。"""
    dy = integrate.quad(lambda y: math.sin(y) * (-y / s ** 2) * delta_gauss(y, s), -12 * s, 12 * s, limit=200)[0]
    dx = integrate.quad(lambda y: math.sin(y) * (y / s ** 2) * delta_gauss(y, s), -12 * s, 12 * s, limit=200)[0]
    return dy, dx


def adjoint_fd(mu, nu, lam, N: int = 64, L: float = 2 * math.pi):
    """(1.5.7)、(1.5.12) 的有限差分版：M = diag(μ)D₂ + diag(ν)D₁ + diag(λ)；
    回傳 (M†, D₂diag(μ*) − D₁diag(ν*) + diag(λ*))——兩者完全相等（D₂ 對稱、D₁ 反對稱）。"""
    x = np.arange(N) * L / N
    D1, D2 = d1_matrix(N, L), d2_matrix(N, L)
    m, n_, l = mu(x), nu(x), lam(x)
    M = np.diag(m) @ D2 + np.diag(n_) @ D1 + np.diag(l)
    Mdag = D2 @ np.diag(np.conj(m)) - D1 @ np.diag(np.conj(n_)) + np.diag(np.conj(l))
    return dag(M), Mdag


def resolvent(M, lam) -> sp.Matrix:
    """(1.5.13) Q_λ(M) = (λ1 − M)⁻¹（精確）。"""
    return sp.simplify((lam * sp.eye(M.shape[0]) - M).inv())


def resolvent_adj(M):
    """(1.5.16) Q_λ = adj(λ1 − M)/S(λ)：回傳 (adj(λ1 − M), S(λ))（λ 是符號 LAM）。"""
    A = LAM * sp.eye(M.shape[0]) - M
    return A.adjugate(), sp.expand(A.det())


def contour(M, center: complex, radius: float, f=None, n: int = 64) -> np.ndarray:
    """(1.5.15)/(1.5.18) 數值圍道積分 (1/2πi)∮ f(λ)(λ1 − M)⁻¹ dλ，圓心 center、半徑 radius、n 點梯形法（對解析函數指數收斂）。"""
    A = to_np(M)
    k = A.shape[0]
    out = np.zeros((k, k), dtype=complex)
    for p in range(n):
        ph = 2 * math.pi * p / n
        lam = center + radius * np.exp(1j * ph)
        dlam = 1j * radius * np.exp(1j * ph) * (2 * math.pi / n)
        w = 1.0 if f is None else f(lam)
        out += w * np.linalg.inv(lam * np.eye(k) - A) * dlam
    return out / (2j * math.pi)


def contour_PN(M, center: float, radius: float = 0.5, n: int = 64):
    """(1.5.18)、(1.5.20)：(P_a, N_a) = ((1/2πi)∮dλ/(λ − M), (1/2πi)∮(λ − λ_a)dλ/(λ − M))。"""
    P = contour(M, center, radius, n=n)
    Nn = contour(M, center, radius, f=lambda z: z - center, n=n)
    return P, Nn


# =============================================================================
# 附錄 B：2×2 反矩陣、伴隨矩陣、Cramer
# =============================================================================
def inverse2(M, book: bool = False) -> sp.Matrix:
    """2×2 反矩陣；book=True 照書上左下分母印成 M₁₁M₂₂ − M₂₂M₂₁ 的寫法。"""
    a, b, c, d = M[0, 0], M[0, 1], M[1, 0], M[1, 1]
    det = a * d - b * c
    den21 = a * d - d * c if book else det
    return sp.Matrix([[d / det, -b / det], [-c / den21, a / det]])


def cramer(M, y) -> sp.Matrix:
    """(B.10)：x_j = det(M 的第 j 行換成 y)/det M。"""
    n = M.shape[0]
    d = M.det()
    out = []
    for j in range(n):
        A = M.copy()
        A[:, j] = sp.Matrix(y)
        out.append(A.det() / d)
    return sp.Matrix(out)


# =============================================================================
# 習題 1.1.3–1.2.8
# =============================================================================
P121 = sp.Matrix([[2, -1], [-6, 1]])
P122 = sp.Matrix([[11, 7, -9], [7, 11, -9], [-9, -9, 27]])
P124 = sp.ones(4, 4)                       # 配對算符 ρ（純量 ρ 取 1）
P125_DATA = (((1, -1, 0), 1), ((1, 1, -2), -1), ((1, 1, 1), 0))
P126_M = sp.Matrix([[2, 1, 0], [1, 2, 0], [0, 0, 3]])
P126_N = sp.Matrix([[3, -1, -1], [-1, 3, -1], [-1, -1, 3]])
P127_H = sp.Matrix([[1, -I], [I, -1]]) / sp.sqrt(2)
P127_A = sp.Matrix([[0, I], [I, 0]]) / sp.sqrt(2)
P127_U = sp.Matrix([[1, 1], [I, -I]]) / sp.sqrt(2)


def sqrt_all(M) -> list:
    """相異本徵值時的全部平方根 Σ ±√ε_j 𝒫_j（2ⁿ 個）。"""
    Ps = list(projectors(M).items())
    out = []
    for signs in range(2 ** len(Ps)):
        R = sp.zeros(M.shape[0])
        for b, (e, P) in enumerate(Ps):
            R += (-1 if (signs >> b) & 1 else 1) * sp.sqrt(e) * P
        out.append(sp.simplify(R))
    return out


def secular_coeffs(M) -> tuple:
    """習題 1.2.3：S(ε) = εⁿ + a₁εⁿ⁻¹ + … + a_n。回傳三種算法：(多項式係數, (−1)^k e_k(本徵值), (−1)^k Σ k 階主子式)。"""
    n = M.shape[0]
    poly = secular(M).all_coeffs()[1:]
    ev = []
    for e, m in eigen_mult(M).items():
        ev += [e] * m
    if not all(e.is_rational for e in ev):                     # 無理根：用數值（避免 sympy 化簡太久）
        ev = [complex(sp.N(e, 30)) for e in ev]
    by_eig = [(-1) ** k * sum(np.prod(c) if isinstance(c[0], complex) else sp.prod(c) for c in combinations(ev, k))
              for k in range(1, n + 1)]
    by_eig = [sp.nsimplify(round(x.real, 9)) if isinstance(x, complex) else x for x in by_eig]
    by_minor = [(-1) ** k * sum(M.extract(list(c), list(c)).det() for c in combinations(range(n), k)) for k in range(1, n + 1)]
    return [sp.simplify(x) for x in poly], [sp.simplify(x) for x in by_eig], [sp.simplify(x) for x in by_minor]


def problem_125() -> sp.Matrix:
    """由本徵解倒推 M = Σ λ |v⟩⟨v|/⟨v|v⟩（Hermitian，向量互相正交）。"""
    M = sp.zeros(3)
    for v, lam in P125_DATA:
        v = sp.Matrix(v)
        M += lam * v * v.T / (v.T * v)[0]
    return M


def shift_matrix(n: int) -> np.ndarray:
    """習題 1.2.8：截斷到 n 維的「升階移位」V = Σ|k+1⟩⟨k|（a† 的方向部分）。無限維時 V†V = 1 但 VV† = 1 − |0⟩⟨0|。"""
    V = np.zeros((n, n))
    for k in range(n - 1):
        V[k + 1, k] = 1
    return V


# =============================================================================
# 摘要與自我測試
# =============================================================================
def summary() -> dict:
    """卡片與講義上的數字（字串）。"""
    s3 = sp.sqrt(3)
    eb = eig_book(MB)
    w_eq = [math.sqrt(float(x)) for x in omega2("equal")]
    w_bt = [math.sqrt(float(x)) for x in omega2("beat")]
    w_un = [math.sqrt(float(x)) for x in omega2("unequal")]
    return {
        "pass_table": [_f4(pass_fraction(a)) for a in PASS_ANGLES],
        "series": {k: tuple(str(x) for x in v) for k, v in series_table().items()},
        "p111_60": [_f4(x) for x in problem_111(60)["1'"]],
        "zeno": [_f4(zeno(n)) for n in (1, 2, 3, 4, 10)],
        "sz_primed_60": str(change_basis(SZ, U(60)).tolist()),
        "beat_t1": [_f4(x) for x in populations(1.0)],
        "beat_period": _f4(math.pi / BEAT_AB[1]),
        "eig_angles": [f"{angle_of(v):.2f}" for v in (eb[1][0], eb[1][1], eb[5][0], eb[5][1])],
        "ket_angle": f"{angle_of(RIGHT_VEC[5]) - angle_of(RIGHT_VEC[1]):.2f}",
        "singular": [_f4(x) for x in singular_values(MB)],
        "iterate": [f"{x:.2f}" for x in iterate_angles(MB)],
        "M50_11": f"{float(mat_func(MB, lambda e: e ** 50)[0, 0]):.4e}",
        "sqrtM": [_f4(x) for x in to_np(mat_func(MB, sp.sqrt)).real.ravel()],
        "expM": [_f4(x) for x in to_np(mat_func(MB, lambda e: sp.exp(e / 10))).real.ravel()],
        "S_probs": [_f4(abs(complex(S_EX[i, 0])) ** 2) for i in (0, 1)],
        "S_phase": f"{math.degrees(math.atan2(4, 3)):.2f}",
        "purity": [str(purity(rho_pure(ket_primed(60)))), str(purity(sp.diag(sp.Rational(3, 4), sp.Rational(1, 4))))],
        "omega": [_f4(x) for x in w_eq + w_bt + w_un],
        "omega2_unequal": [_f4(x) for x in omega2("unequal")],
        "beat_times": [f"{x:.2f}" for x in beat_times("beat")],
        "delta": [_f4(delta_integral(s)[1]) for s in DELTA_SIGMAS],
        "delta_deriv_01": [_f4(x) for x in delta_deriv_integrals(0.1)],
        "d1_eigen": [_f4(d1_eigen(N).imag) for N in (8, 32)],
        "inv_book_21": str(inverse2(MB, book=True)[1, 0]),
        "sqrt3": _f4(s3 / 2),
    }


SUMMARY_EXPECT = {
    "pass_table": ["1.0000", "0.9330", "0.7500", "0.5000", "0.2500", "0.0670", "0.0000"],
    "series": {"direct": ("0", "1"), "coherent": ("0", "1"), "which": ("1/2", "1/2"), "stop2": ("1/4", "1/4")},
    "p111_60": ["0.7500", "0.4665"],
    "zeno": ["0.5000", "0.7286", "0.8122", "0.8562", "0.9401"],
    "sz_primed_60": "[[1/2, sqrt(3)*I/2], [-sqrt(3)*I/2, -1/2]]",
    "beat_t1": ["0.2919", "0.7081"],
    "beat_period": "3.1416",
    "eig_angles": ["-45.00", "-71.57", "18.43", "45.00"],
    "ket_angle": "116.57",
    "singular": ["5.3983", "0.9262"],
    "iterate": ["0.00", "36.87", "43.45", "44.69", "44.94"],
    "M50_11": "6.6613e+34",
    "sqrtM": ["1.9271", "0.3090", "0.9271", "1.3090"],
    "expM": ["1.5128", "0.1359", "0.4077", "1.2411"],
    "S_probs": ["0.3600", "0.6400"],
    "S_phase": "53.13",
    "purity": ["1", "5/8"],
    "omega": ["1.0000", "1.7321", "1.0000", "1.2247", "0.7962", "1.5382"],
    "omega2_unequal": ["0.6340", "2.3660"],
    "beat_times": ["13.98", "27.96"],
    "delta": ["0.6065", "0.8825", "0.9950"],
    "delta_deriv_01": ["-0.9950", "0.9950"],
    "d1_eigen": ["0.9003", "0.9936"],
    "inv_book_21": "-3/2",
    "sqrt3": "0.8660",
}


def check():
    R = sp.Rational
    s2, s3, s5 = sp.sqrt(2), sp.sqrt(3), sp.sqrt(5)
    rng = np.random.default_rng(1993)
    assert summary() == SUMMARY_EXPECT

    # ---- H101：分析器、Axiom 4、習題 1.1.1
    T60 = U(60)
    assert T60 == sp.Matrix([[s3 / 2, -I / 2], [-I / 2, s3 / 2]]) and _close(to_np(T60), U_num(60))
    assert pass_fraction(60) == R(3, 4) and abs(abs(U_num(60)[0, 0]) ** 2 - 0.75) < 1e-14 and abs(abs(U_num(60)[1, 0]) ** 2 - 0.25) < 1e-14
    assert _zero(U(30) * U(60) - U(90)) and _zero(T(180, 90) * T(90, 0) - T(180, 0))            # (1.1.6b)
    for a in PASS_ANGLES:
        assert _close(float(pass_fraction(a)), math.cos(math.radians(a) / 2) ** 2, 1e-14)
    st = series_table()
    assert st == {"direct": (0, 1), "coherent": (0, 1), "which": (R(1, 2), R(1, 2)), "stop2": (R(1, 4), R(1, 4))}
    for k in (1, 2):                                       # 振幅相加（Axiom 4）＝ 相干；機率相加 ＝ which
        coh, inc = axiom4(0, 90, 180, k)
        assert (coh, inc) == (st["coherent"][k - 1], st["which"][k - 1])
    amp = U_num(180 - 90)[0, 0] * U_num(90 - 0)[0, 0]       # ⟨1_180|1_90⟩⟨1_90|1_0⟩
    assert abs(abs(amp) ** 2 - 0.25) < 1e-14
    assert axiom4(0, 90, 180, 1)[0] - axiom4(0, 90, 180, 1)[1] == R(-1, 2)
    for th in (30, 60, 75):
        p = problem_111(th)
        c = problem_111_closed(th)
        assert all(_close([float(x) for x in p[key]], c[key], 1e-12) for key in p)
    p60 = problem_111(60)
    assert p60["1"] == (1, R(1, 4)) and p60["2"] == (0, R(1, 4)) and p60["1'"][0] == R(3, 4)
    assert sp.simplify(p60["1'"][1] - (2 + s3) / 8) == 0
    assert U(90) == sp.Matrix([[s2 / 2, -I * s2 / 2], [-I * s2 / 2, s2 / 2]]) and _f4(s2 / 2) == "0.7071"

    # ---- H102：么正性、正交的實特例、擋板不是么正、習題 1.1.2
    assert is_unitary(T60) and T60.H == sp.Matrix([[s3 / 2, I / 2], [I / 2, s3 / 2]])
    for th in rng.uniform(0, 360, 5):
        Un = U_num(th)
        v = rng.normal(size=2) + 1j * rng.normal(size=2)
        assert is_unitary(Un) and abs(np.linalg.norm(Un @ v) - np.linalg.norm(v)) < 1e-12
    Tr = T_real(60)
    assert Tr == sp.Matrix([[s3 / 2, -R(1, 2)], [R(1, 2), s3 / 2]]) and _zero(Tr.T * Tr - sp.eye(2))
    assert all(sp.simplify(T_real(a)[0, 0] ** 2 - pass_fraction(a)) == 0 for a in (30, 60, 90))
    assert not is_unitary(E11) and E11 * E11 == E11
    for n in (1, 2, 3, 4, 10):
        assert abs(zeno_chain(n) - zeno(n)) < 1e-12
    assert abs(float(analyzer_chain([(45, 1), (90, "count")])[0]) - zeno(2)) < 1e-14
    for key, val in series_table().items():
        st_ = {"direct": [(180, "count")], "coherent": [(90, "coherent"), (180, "count")],
               "which": [(90, "which"), (180, "count")], "stop2": [(90, 1), (180, "count")]}[key]
        assert _close([float(x) for x in val], analyzer_chain_num(st_), 1e-14)
    assert zeno(1000) > 0.999 and abs(zeno(1) - 0.5) < 1e-15

    # ---- H103：ket、基本算子、換基底
    k1 = ket_primed(60)
    assert k1 == sp.Matrix([s3 / 2, -I / 2]) and k1 == T60 * basis(1)                         # (1.1.13)
    assert (k1.H * k1)[0] == 1                                                                # 3/4 + 1/4
    e12 = elementary(1, 2)
    assert e12 * basis(2) == basis(1) and basis(1).T * e12 == basis(2).T                       # (1.1.19a,b)
    assert e12 * basis(1) == sp.zeros(2, 1) and basis(2).T * e12 == sp.zeros(1, 2)             # (1.1.19c)，k ≠ i
    assert basis(1).T * e12 != sp.zeros(1, 2)                                                 # 書上的 k ≠ j 不夠
    P1p = change_basis(E11, T60)
    assert P1p == sp.Matrix([[R(3, 4), s3 * I / 4], [-s3 * I / 4, R(1, 4)]]) and _zero(P1p - outer(k1, k1))
    szp = change_basis(SZ, T60)
    assert _zero(szp - (SZ / 2 - s3 / 2 * SY)) and szp.eigenvals() == {1: 1, -1: 1}
    assert szp.trace() == 0 and sp.simplify(szp.det()) == -1
    Hr = rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3))
    Hr = Hr + dag(Hr)
    Q, _ = np.linalg.qr(rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3)))
    assert _close(np.linalg.eigvalsh(Hr), np.linalg.eigvalsh(change_basis(Hr, Q)), 1e-10)
    Hs = sp.Matrix([[1, 2 + I], [2 - I, -3]])
    assert Hs == sum((Hs[i, j] * elementary(i + 1, j + 1) for i in range(2) for j in range(2)), sp.zeros(2))   # (1.1.18)

    # ---- H104：兩態拍頻
    A, B = BEAT_AB
    assert H2(A, B).eigenvals() == {4: 1, 2: 1}
    assert H2(A, B) * sp.Matrix([1, 1]) == 4 * sp.Matrix([1, 1]) and H2(A, B) * sp.Matrix([1, -1]) == 2 * sp.Matrix([1, -1])
    ts = np.array([0, math.pi / 4, 1.0, math.pi / 2, 2.5, math.pi])
    p1n, p2n = populations_numeric(ts, A, B)
    for i, t in enumerate(ts):
        pc = populations(t, A, B)
        ps = np.abs(psi_t(t, A, B)) ** 2
        pe = np.abs(linalg.expm(-1j * t * np.array([[A, B], [B, A]])) @ np.array([1, 0])) ** 2
        assert _close(pc, ps, 1e-12) and _close(pc, pe, 1e-12) and _close(pc, (p1n[i], p2n[i]), 1e-8)
        ev = abs(np.vdot(np.array([1, 1]) / math.sqrt(2), psi_t(t, A, B))) ** 2
        assert abs(ev - 0.5) < 1e-12                                                          # 本徵態只轉相位
    assert [_f4(x) for x in populations(math.pi / 4)] == ["0.5000", "0.5000"]
    assert [_f4(x) for x in populations(math.pi / 2)] == ["0.0000", "1.0000"] and _f4(math.pi / 2) == "1.5708"
    assert _close(populations(math.pi), (1, 0), 1e-12)

    # ---- H105：久期方程、Hamilton–Cayley、左右本徵向量
    S = secular(MB)
    assert S == sp.Poly(EPS ** 2 - 6 * EPS + 5, EPS) == sp.Poly(EPS ** 2 - MB.trace() * EPS + MB.det(), EPS)
    assert sp.Poly((MB - EPS * sp.eye(2)).det(), EPS) == S and distinct_eigs(MB) == [1, 5]
    assert MB ** 2 == sp.Matrix([[19, 6], [18, 7]]) and 6 * MB == sp.Matrix([[24, 6], [18, 12]])
    assert poly_at(S, MB) == sp.zeros(2) and (MB - sp.eye(2)) * (MB - 5 * sp.eye(2)) == sp.zeros(2)
    for e in (1, 5):
        r, l = sp.Matrix(RIGHT_VEC[e]), sp.Matrix([LEFT_VEC[e]])
        assert MB * r == e * r and l * MB == e * l
    assert sp.Matrix([LEFT_VEC[1]]) * sp.Matrix(RIGHT_VEC[5]) == sp.zeros(1) == sp.Matrix([LEFT_VEC[5]]) * sp.Matrix(RIGHT_VEC[1])
    w, V = np.linalg.eig(to_np(MB).real)
    for lam, v in zip(w, V.T):
        e = int(round(lam))
        assert abs(abs(np.dot(v, RIGHT_VEC[e])) - np.linalg.norm(RIGHT_VEC[e])) < 1e-10
    eb = eig_book(MB)
    assert eb[1] == (sp.Matrix([[R(1, 2), -R(1, 2)]]), sp.Matrix([R(1, 2), -R(3, 2)]))
    assert eb[5] == (sp.Matrix([[s3 / 2, s3 / 6]]), sp.Matrix([s3 / 2, s3 / 2]))
    for a in (1, 5):
        for b in (1, 5):
            assert sp.simplify((eb[a][0] * eb[b][1])[0]) == (1 if a == b else 0)            # 對偶 ⟨e_i|e_j⟩ = δ
        assert _zero(eb[a][1] * eb[a][0] - projectors(MB)[a])                                 # (1.2.14)
    sv = singular_values(MB)
    ev2 = np.sqrt(np.linalg.eigvalsh(to_np(MB.T * MB).real))[::-1]
    assert _close(sv, ev2, 1e-12) and abs(sv[0] * sv[1] - 5) < 1e-12
    exact_angles = [math.degrees(math.atan2(*reversed([float(x) for x in MB ** k * sp.Matrix([1, 0])]))) for k in range(5)]
    assert _close(iterate_angles(MB), exact_angles, 1e-10)
    pts = ellipse_points(MB, 8)
    assert pts.shape == (2, 8) and _close(pts[:, 0], [4, 3], 1e-12)

    # ---- H106：投影算符、譜分解、f(M)
    pf = hc_factors(MB)
    assert pf == {1: sp.Matrix([[-1, 1], [3, -3]]), 5: sp.Matrix([[3, 1], [3, 1]])}             # (1.2.9)_x
    assert pf[1] * pf[1] == -4 * pf[1] and pf[5] * pf[5] == 4 * pf[5] and pf[1] * pf[5] == sp.zeros(2)   # (1.2.11)
    P = projectors(MB)
    assert P[1] == sp.Matrix([[1, -1], [-3, 3]]) / 4 and P[5] == sp.Matrix([[3, 1], [3, 1]]) / 4
    assert P[1] ** 2 == P[1] and P[5] ** 2 == P[5] and P[1] * P[5] == sp.zeros(2) == P[5] * P[1]
    assert P[1] + P[5] == sp.eye(2) and P[1] + 5 * P[5] == MB                                   # (1.2.18)(1.2.21)
    book_1218 = sp.Matrix([[R(1, 4), -R(1, 4)], [-R(3, 4), R(1, 4)]]) + sp.Matrix([[R(3, 4), R(1, 4)], [R(1, 4), R(1, 4)]])
    assert book_1218 == sp.Matrix([[1, 0], [-R(1, 2), R(1, 2)]]) != sp.eye(2)                   # (1.2.18)_x 的印刷錯
    v = sp.Matrix([1, 0])
    assert P[1] * v == sp.Matrix([R(1, 4), -R(3, 4)]) and P[5] * v == sp.Matrix([R(3, 4), R(3, 4)])
    M50 = mat_func(MB, lambda e: e ** 50)
    f = 5 ** 50
    assert M50 == MB ** 50 == sp.Matrix([[1 + 3 * f, f - 1], [3 * f - 3, 3 + f]]) / 4
    assert str(f) == "88817841970012523233890533447265625" and M50[0, 0] == 66613381477509392425417900085449219
    rM = mat_func(MB, sp.sqrt)
    assert rM == sp.Matrix([[1 + 3 * s5, s5 - 1], [3 * s5 - 3, 3 + s5]]) / 4 and sp.simplify(rM ** 2 - MB) == sp.zeros(2)
    assert _close(to_np(rM), linalg.sqrtm(to_np(MB).real), 1e-10)
    roots4 = [sg1 * P[1] + sg5 * s5 * P[5] for sg1 in (1, -1) for sg5 in (1, -1)]
    assert all(sp.simplify(r_ ** 2 - MB) == sp.zeros(2) for r_ in roots4)
    eM = mat_func(MB, lambda e: sp.exp(e / 10))
    assert _close(to_np(eM), linalg.expm(0.1 * to_np(MB).real), 1e-12)
    x = sp.Symbol("x")
    lb = lagrange_basis([1, 5], x)
    assert sp.simplify(sum(lb)) == 1 and sp.simplify(lb[0] * 1 + lb[1] * 5 - x) == 0
    fake = [(MB - 7 * sp.eye(2)) / (2 - 7), (MB - 2 * sp.eye(2)) / (7 - 2)]                    # ε 隨便代：完備仍成立
    assert fake[0] + fake[1] == sp.eye(2) and fake[0] ** 2 != fake[0]

    # ---- H107：重根與最小方程、對易矩陣、伴隨矩陣、附錄 C
    assert secular(DEG4) == sp.Poly((EPS - 1) ** 2 * (EPS + 1) ** 2, EPS) == sp.Poly(EPS ** 4 - 2 * EPS ** 2 + 1, EPS)
    mp = minimal_poly(DEG4)
    assert mp == sp.Poly(EPS ** 2 - 1, EPS) and poly_at(mp, DEG4) == sp.zeros(4)
    assert DEG4 - sp.eye(4) != sp.zeros(4) and DEG4 + sp.eye(4) != sp.zeros(4)
    PD = projectors(DEG4)
    assert PD[1] == (DEG4 + sp.eye(4)) / 2 and PD[-1] == (sp.eye(4) - DEG4) / 2                   # (1.2.16)
    assert PD[1].rank() == PD[-1].rank() == 2 and PD[1].trace() == PD[-1].trace() == 2
    assert PD[1] + PD[-1] == sp.eye(4) and PD[1] * PD[-1] == sp.zeros(4)
    wD = np.linalg.eigvalsh(to_np(DEG4).real)
    assert _close(sorted(wD), [-1, -1, 1, 1], 1e-12)
    rows = [PD[1][0, :] * s2, PD[1][1, :] * s2]
    assert rows[0] == sp.Matrix([[1, 0, 0, 1]]) / s2 and rows[1] == sp.Matrix([[0, 1, 1, 0]]) / s2
    assert minimal_poly(JORDAN) == sp.Poly((EPS - 2) ** 2 * (EPS - 3), EPS)
    Nil = (JORDAN - 2 * sp.eye(3)) * (JORDAN - 3 * sp.eye(3))                                  # (C.2a)
    assert Nil != sp.zeros(3) and Nil * Nil == sp.zeros(3)
    assert MB * NB == NB * MB == sp.Matrix([[17, 3], [9, 11]]) and NB == 8 * P[1] + 4 * P[5]
    assert distinct_eigs(NB) == [4, 8]
    jp = joint_projectors(MB, NB)
    assert jp == {(1, 8): P[1], (5, 4): P[5]}
    jp6 = joint_projectors(P126_M, P126_N)
    assert set(jp6) == {(1, 4), (3, 1), (3, 4)} and sum(jp6.values(), sp.zeros(3)) == sp.eye(3)
    Tj = sp.Matrix([[1, -1, 0], [1, 1, 1], [1, 1, -2]])
    Tj = sp.diag(1 / s2, 1 / s3, 1 / sp.sqrt(6)) * Tj
    assert _zero(Tj * Tj.T - sp.eye(3)) and _zero(Tj * P126_M * Tj.T - sp.diag(1, 3, 3)) and _zero(Tj * P126_N * Tj.T - sp.diag(4, 1, 4))
    assert adj_vectors(MB, 1) == sp.Matrix([[1, -1], [-3, 3]]) == 4 * P[1]
    assert adj_vectors(MB, 5) == sp.Matrix([[-3, -1], [-3, -1]]) == -4 * P[5]
    for Mx in (MB, P122):
        for e, p in hc_factors(Mx).items():
            assert (e * sp.eye(Mx.shape[0]) - Mx).adjugate() == p                               # adj(ε_j1 − M) = 𝔭_j

    # ---- H108：散射、密度矩陣、耦合振子
    assert is_unitary(S_EX) and S_EX.eigenvals() == {R(3, 5) + 4 * I / 5: 1, R(3, 5) - 4 * I / 5: 1}
    for sg in (1, -1):
        vv = sp.Matrix([1, sg])
        assert S_EX * vv == (R(3, 5) + sg * 4 * I / 5) * vv
    assert [abs(S_EX[i, 0]) ** 2 for i in (0, 1)] == [R(9, 25), R(16, 25)]
    for _ in range(3):
        fi, fo = flux_in_out(S_EX, rng.normal(size=2) + 1j * rng.normal(size=2))
        assert abs(fi - fo) < 1e-12                                                           # (1.3.2)
    Hc, Ac = S_EX + S_EX.H, S_EX - S_EX.H                                                     # (C.6)–(C.7)
    assert Hc == R(6, 5) * sp.eye(2) and Ac == R(8, 5) * I * SX and Hc * Ac == Ac * Hc
    rho = rho_pure(k1)
    assert rho == P1p and rho.trace() == 1 and _zero(rho * rho - rho) and purity(rho) == 1
    rho_mix = analyzer_chain([(0, "which"), (0, "count")], start=(60, 1))                      # 60° 的 |1⟩、量了上下
    assert rho_mix == (R(3, 4), R(1, 4)) and purity(sp.diag(*rho_mix)) == R(5, 8) and purity(sp.eye(2) / 2) == R(1, 2)
    for case in OSC:
        m1, m2, k1_, k2_, k12 = OSC[case]
        F = force_matrix(k1_, k2_, k12)
        a = accel_matrix(*OSC[case])
        assert F == F.T and a == mass_matrix(m1, m2).inv() * F
        w2, E = omega2_generalized(case)
        assert _close(w2, [float(x) for x in omega2(case)], 1e-12)
        assert _close(sorted(np.linalg.eigvals(to_np(a).real).real), w2, 1e-12)
        mm = np.diag([float(m1), float(m2)])
        assert _close(E.T @ mm @ E, np.eye(2), 1e-12) and _close(E.T @ to_np(F).real @ E, np.diag(w2), 1e-12)   # (1.4.13)
        for j in range(2):                                                                    # ⟨a_j| ∝ ⟨e_j|m：左本徵向量
            row = E[:, j] @ mm
            assert _close(row @ to_np(a).real, w2[j] * row, 1e-12)
        assert np.all(np.linalg.eigvals(-to_np(a).real).real < 0)                             # −m⁻¹𝓕 不能給 ω²
    assert distinct_eigs(-accel_matrix(*OSC["equal"])) == [-3, -1]
    assert omega2("equal") == [1, 3] and omega2("beat") == [1, R(3, 2)]
    assert omega2("unequal") == [R(3, 2) - s3 / 2, R(3, 2) + s3 / 2]
    assert accel_matrix(*OSC["unequal"]) == sp.Matrix([[2, -1], [-R(1, 2), 1]])
    a_s, w2s, sy = modes_symbolic()
    for case in OSC:
        sub = dict(zip(sy, OSC[case]))
        assert [sp.simplify(z.subs(sub)) for z in w2s] == omega2(case)
    ts = np.linspace(0, 30, 61)
    for case in ("beat", "unequal"):
        xa = osc_analytic(ts, case, (1.0, 0.0), (0.0, 0.0))
        xn = osc_numeric(ts, case, (1.0, 0.0), (0.0, 0.0))
        assert np.max(np.abs(xa - xn)) < 1e-7
    tq = beat_times("beat")[0]
    x_end = osc_analytic(tq, "beat")[:, 0]
    assert abs(x_end[0]) < 0.25 and abs(osc_energy((1, 0), (0, 0)) - 0.625) < 1e-12           # m₁ 在 t = 13.98 幾乎靜止
    h = 1e-6
    vel = (osc_analytic(tq + h, "beat") - osc_analytic(tq - h, "beat"))[:, 0] / (2 * h)
    assert abs(osc_energy(x_end, vel) - 0.625) < 1e-6
    w2, E = omega2_generalized("unequal")                                                     # H = ½Σ(p² + ω²q²)
    xr, vr = rng.normal(size=2), rng.normal(size=2)
    mm = np.diag([1.0, 2.0])
    q, pp = E.T @ mm @ xr, E.T @ mm @ vr
    assert abs(0.5 * np.sum(pp ** 2 + w2 * q ** 2) - osc_energy(xr, vr, "unequal")) < 1e-12     # (1.4.19)
    springs = [(1.3, {0: 0.6, 1: 0.8}), (0.7, {0: -0.28, 1: 0.96}), (2.0, {0: 1.0, 2: -0.5, 3: 0.866})]
    Fs = spring_force_matrix(springs, 4)

    def V(chi):
        return 0.5 * sum(k * sum(c * chi[i] for i, c in cs.items()) ** 2 for k, cs in springs)
    Hess = np.zeros((4, 4))
    for i in range(4):
        for j in range(4):
            ei, ej = np.eye(4)[i] * 1e-4, np.eye(4)[j] * 1e-4
            Hess[i, j] = (V(ei + ej) - V(ei - ej) - V(-ei + ej) + V(-ei - ej)) / (4e-8)
    assert _close(Fs, Hess, 1e-6) and _close(Fs, Fs.T, 0)
    assert abs(Fs[0, 1] - (1.3 * 0.6 * 0.8 + 0.7 * -0.28 * 0.96)) < 1e-12                       # χ₂ 係數 = Σk cosφ₁cosφ₂

    # ---- H109：d/dx、δ 序列、預解式與圍道積分
    for N in (8, 32):
        D1 = d1_matrix(N)
        assert _close(D1.T, -D1, 0) and _close(1j * D1, dag(1j * D1), 0)
        Fm = fourier_matrix(N)
        assert is_unitary(Fm)
        Dd = dag(Fm) @ D1 @ Fm
        assert _close(Dd, np.diag(np.diag(Dd)), 1e-12) and _close(np.diag(Dd)[1], d1_eigen(N), 1e-12)
        assert np.max(np.abs(np.linalg.eigvals(D1).real)) < 1e-10
        assert _close(Fm @ dag(Fm), np.eye(N), 1e-12)                                           # Σ_k⟨x|k⟩⟨k|y⟩ = δ
    for s in DELTA_SIGMAS:
        num, cl = delta_integral(s)
        assert abs(num - cl) < 1e-10
        assert abs(integrate.quad(lambda x_: delta_gauss(x_, s), -12 * s, 12 * s)[0] - 1) < 1e-10
        dy, dx = delta_deriv_integrals(s)
        assert abs(dy + cl) < 1e-10 and abs(dx - cl) < 1e-10
    xg = np.arange(400) * 2 * math.pi / 400                                                    # 有限差分：∂_x 那一種就是 D₁
    assert np.max(np.abs(d1_matrix(400) @ np.sin(xg) - np.cos(xg))) < 1e-3
    Md, Mf = adjoint_fd(lambda x_: 1 + 0.5 * np.sin(x_), lambda x_: 1j * np.cos(x_), lambda x_: 0.3 + 0j * x_)
    assert _close(Md, Mf, 1e-9)
    Q3 = resolvent(MB, 3)
    adjA, Sl = resolvent_adj(MB)
    assert Q3 == -sp.Matrix([[1, 1], [3, -1]]) / 4 == P[1] / 2 - P[5] / 2
    assert adjA == sp.Matrix([[LAM - 2, 1], [3, LAM - 4]]) and Sl == LAM ** 2 - 6 * LAM + 5
    assert Q3 == adjA.subs(LAM, 3) / Sl.subs(LAM, 3)
    lam_ = sp.Symbol("mu")
    assert _zero(resolvent(MB, lam_) - (P[1] / (lam_ - 1) + P[5] / (lam_ - 5)))
    assert _close(contour(MB, 1, 1), to_np(P[1]), 1e-12) and _close(contour(MB, 5, 1), to_np(P[5]), 1e-12)
    assert _close(contour(MB, 3, 4), np.eye(2), 1e-12) and _close(contour(MB, 3, 4, f=lambda z: z), to_np(MB), 1e-11)
    assert _close(contour(MB, 3, 4, f=lambda z: np.exp(0.1 * z), n=128), linalg.expm(0.1 * to_np(MB).real), 1e-11)
    c50 = contour(MB, 3, 3, f=lambda z: z ** 50, n=256)
    assert np.max(np.abs(c50 - to_np(M50))) / float(M50[0, 0]) < 1e-10
    P2, N2 = contour_PN(JORDAN, 2)
    P3, N3 = contour_PN(JORDAN, 3)
    assert _close(P2, np.diag([1, 1, 0]), 1e-12) and _close(N2, to_np(elementary(1, 2, 3)), 1e-12)
    assert _close(N2 @ N2, np.zeros((3, 3)), 1e-12) and _close(N3, np.zeros((3, 3)), 1e-12)
    assert _close(2 * P2 + N2 + 3 * P3 + N3, to_np(JORDAN), 1e-12)                              # (1.5.20) 更正後
    assert not _close(2 * (2 * P2 + N2) + 3 * (3 * P3 + N3), to_np(JORDAN), 1e-3)               # 書上第一行多了 λ_a
    assert _close(P2 + N2 + P3 + N3, np.eye(3) + to_np(elementary(1, 2, 3)), 1e-12)             # 書上末行少了 λ_a：1 + e₁₂ ≠ M
    assert P[1] + 25 * P[5] == MB ** 2                                                          # MB 照書上第一行 → M²

    # ---- 附錄 B
    assert inverse2(MB) == MB.inv() == sp.Matrix([[2, -1], [-3, 4]]) / 5 and inverse2(MB, book=True)[1, 0] == -R(3, 2)
    assert cramer(MB, [1, 0]) == MB.inv() * sp.Matrix([1, 0]) == sp.Matrix([R(2, 5), -R(3, 5)])
    M4 = sp.Matrix(4, 4, lambda i, j: int(rng.integers(-3, 4)) + (5 if i == j else 0))
    y4 = sp.Matrix([1, 2, -1, 3])
    assert M4 * M4.adjugate() == M4.det() * sp.eye(4) and cramer(M4, y4) == M4.inv() * y4      # (B.8)(B.10)
    A_bad = M4.copy()
    A_bad[:, 1] = y4
    A_bad[0, 3] = M4[1, 3]                                                                     # 書上 (B.10) 第二個行列式的 b₄
    assert A_bad.det() / M4.det() != (M4.inv() * y4)[1] or M4[1, 3] == M4[0, 3]

    # ---- 習題
    Ar, Br = (rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3)) for _ in range(2))
    xv, yv = rng.normal(size=3) + 1j * rng.normal(size=3), rng.normal(size=3) + 1j * rng.normal(size=3)
    assert _close(dag(Ar @ Br), dag(Br) @ dag(Ar), 1e-12) and _close(dag(Ar @ xv[:, None]), xv.conj()[None, :] @ dag(Ar), 1e-12)
    assert abs(np.conj(xv.conj() @ Ar @ yv) - yv.conj() @ dag(Ar) @ xv) < 1e-12 and dag(elementary(1, 2)) == elementary(2, 1)
    P1x = projectors(P121)
    assert secular(P121) == sp.Poly(EPS ** 2 - 3 * EPS - 4, EPS)
    assert distinct_eigs(P121) == [-1, 4] and P1x[4] == sp.Matrix([[3, -1], [-6, 2]]) / 5 and P1x[-1] == sp.Matrix([[2, 1], [6, 3]]) / 5
    assert 4 * P1x[4] - P1x[-1] == P121 and len(sqrt_all(P121)) == 4
    assert all(sp.simplify(r_ ** 2 - P121) == sp.zeros(2) for r_ in sqrt_all(P121))
    assert sqrt_all(P121)[0] == sp.Matrix([[6 + 2 * I, -2 + I], [-12 + 6 * I, 4 + 3 * I]]) / 5
    assert distinct_eigs(P122) == [4, 9, 36] and mat_func(P122, sp.sqrt) == sp.Matrix([[3, 1, -1], [1, 3, -1], [-1, -1, 5]])
    assert mat_func(P122, sp.sqrt) ** 2 == P122 and len(sqrt_all(P122)) == 8
    c1, c2, c3 = secular_coeffs(P122)
    assert c1 == c2 == c3 == [-49, 504, -1296]
    Mr = sp.Matrix(3, 3, lambda i, j: int(rng.integers(-4, 5)))
    Tm_ = sp.Matrix([[1, 2, 0], [0, 1, 3], [1, 0, 1]])
    cr = secular_coeffs(Mr)
    assert cr[0] == cr[1] == cr[2] == secular(Tm_ * Mr * Tm_.inv()).all_coeffs()[1:]   # (c) 不變
    assert secular(P124) == sp.Poly(EPS ** 4 - 4 * EPS ** 3, EPS) and minimal_poly(P124) == sp.Poly(EPS ** 2 - 4 * EPS, EPS)
    PJ = projectors(P124)
    assert PJ[4] == P124 / 4 and PJ[0] == sp.eye(4) - P124 / 4
    rt = P124 / 2
    Nn = sp.Matrix([1, -1, 0, 0]) * sp.Matrix([[0, 0, 1, -1]])
    assert rt ** 2 == P124 and (rt + Nn) ** 2 == P124 and (-rt + 7 * Nn) ** 2 == P124            # 無限多個平方根
    M5 = problem_125()
    assert M5 == sp.Matrix([[1, -2, 1], [-2, 1, 1], [1, 1, -2]]) / 3 and M5 == M5.H
    for vv, lam in P125_DATA:
        assert M5 * sp.Matrix(vv) == lam * sp.Matrix(vv)
    assert M5 ** 100 == sp.eye(3) - sp.ones(3, 3) / 3 == mat_func(M5, lambda e: e ** 100)
    for X, ev_ in ((P127_H, {1, -1}), (P127_A, {I / s2, -I / s2})):
        assert set(X.eigenvals()) == ev_ and all(_zero(Pp - Pp.H) for Pp in projectors(X).values())
    assert sp.simplify(P127_U.det()) == -I
    Un = to_np(P127_U)
    assert is_unitary(Un) and is_unitary(P127_U) and _close(sorted(np.degrees(np.angle(np.linalg.eigvals(Un)))), [-105, 15], 1e-9)
    wU, VU = np.linalg.eig(Un)
    assert abs(np.vdot(VU[:, 0], VU[:, 1])) < 1e-12                                            # 么正 ⇒ 本徵向量正交
    assert P[1].H != P[1]                                                                     # 非 Hermitian 的 MB 沒有這性質
    Qr, _ = np.linalg.qr(rng.normal(size=(4, 4)) + 1j * rng.normal(size=(4, 4)))
    assert _close(dag(Qr) @ Qr, np.eye(4), 1e-12) and _close(Qr @ dag(Qr), np.eye(4), 1e-12)     # 有限維 ⇒ 兩邊都成立
    Vs = shift_matrix(6)
    assert _close(dag(Vs) @ Vs, np.diag([1, 1, 1, 1, 1, 0]), 0) and _close(Vs @ dag(Vs), np.diag([0, 1, 1, 1, 1, 1]), 0)


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
    print(f"sym_mat（第 1 章）：分析器與 Axiom 4（密度矩陣＝振幅求和）、么正性、bra/ket 與換基底、兩態拍頻（閉式＝投影算符＝expm＝ODE）、"
          f"久期方程與 Hamilton–Cayley、投影算符（精確）、M⁵⁰／√M／e^Mt（譜分解＝直接）、最小方程、對易矩陣、伴隨矩陣、散射與密度矩陣、"
          f"耦合振子（解析＝數值）、d/dx 與 δ 序列、預解式圍道積分＝投影算符、習題 1.1.1–1.2.8 全部核對通過（{time.time() - t0:.1f} 秒）。")
