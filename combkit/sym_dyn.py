# 由 combkit/tools/sync_from_harter_sym.py 從 harter_sym 專案的 codes/utils/sym_dyn.py 複製；請改 harter_sym 那份再同步。
r"""
sym_dyn.py —— Harter《Principles of Symmetry, Dynamics, and Spectroscopy》（Wiley 1993）第 8 章
Symmetry Analysis for Semiclassical and Quantum Mechanics: Dynamics with High Quanta 的數學
（純 Python＋numpy＋sympy＋scipy，不 import manim；J 矩陣、𝒟^j 用同資料夾的 sym_rot，自旋向量、曲柄向量、轉子能階用 sym_coup，
不重寫）。書上的式號、筆誤與〔注〕見 ref/md/ch08_formulas.md（相對於 C:\mypython\manim\projects\harter_sym）。

約定
  * 諧振子一律 ħ = μ = ω = 1（書上 §8.3 也取 ω = 1）：a = (x + ip)/√2、a† = (x − ip)/√2（(8.2.3)），x = (a + a†)/√2、p = i(a† − a)/√2；
    α = (x0 + ip0)/√2（(8.2.15b)），α_t = α e^{−it}（順時針）。截斷的升降算子 ladder(N)：a[n−1, n] = √n。
  * 平移與 BCH：D(x0) = e^{−ix0 p}、B(p0) = e^{ip0 x}、Q(x0, p0) = e^{i(p0 x − x0 p)} = e^{αa† − α*a} = e^{−|α|²/2}e^{αa†}e^{−α*a}
    （書上 (8.2.15a)、(8.2.21) 把 −α*a 印成 +α*a；照印的 e^{αa† + α*a} 不是么正：⟨0|·|0⟩ = e^{+|α|²/2}）。
  * 自由落體／拋射體 H = p²/2m − f y（(8.1.21)，f 是 +y 方向的力）；FALL = {eps: 128 J, F: −4 N, m: 1 kg}（書上沒給質量，本系列取 1 kg）。
  * 兩能階：照書上 (8.5.4a)、(8.5.5) 的矩陣，|1⟩ 的對角元 +ε/2；曲柄 ω = (2Re H21, 2Im H21, H11 − H22)（(8.5.6)，＝ SC.h_to_omega），
    自旋向量 S = ⟨σ⟩/2 ＝ SC.spin_vector，dS/dt = ω × S。失諧 Δ = ε − Ω；轉動座標 ψ^R = diag(e^{iΩt/2}, e^{−iΩt/2})ψ（(8.5.10)–(8.5.11)）。
    熱平衡 ρ ∝ e^{−βH} ⇒ ⟨S_z⟩ = −½ tanh(βε/2)（|1⟩ 在上；書 p. 752 把能量低的叫 1，得 +½ tanh，只差命名）。
  * Jaynes–Cummings：照 (8.7.6) 與 Fig. 8.7.1，下能階 1 = ↓（能量 0）、上能階 2 = ↑（能量 ε）；基底 |n⟩|s⟩ 的序號 2n + s（s = 0 ↓、1 ↑）；
    H'_JCM = Ω a†a + ε(σz + 1)/2 + (g/2)(a†σ− + aσ+)（實的版本，(8.7.10)）；第 N 塊 {|N−1⟩|↑⟩, |N⟩|↓⟩}。
  * §8.4、§8.6 的物理常數用 SI（CODATA 2018）；書上 (8.6.32b) 是 SI 形式 A = 4ω³q²|r|²/(3ħc³·4πε0)，Gauss 單位把 q²/4πε0 換成 q²。
  * 卡片上的小數一律 f"{float(x):.4f}" 之類（Fraction、sympy 數要 float()）；不要把模組代號 SC、RT、SM 或 D、E、S、T、I、G、J 當變數名。

§8.1A 接觸轉換與 Legendre 轉換（H801）
  legendre_numeric(f, Xs, xs)（Y(X) = max_x [xX − f(x)]，(8.1.2c)）、line_family(f, xjs)（(slope, intercept) = (x_j, −f(x_j))）、
  line_envelope(f, fp, xjs)（切點 (X_j, Y_j) = (f'(x_j), x_j f'(x_j) − f(x_j))，(8.1.2a)）＝ 直線族逐點取極大
    LINE_CARD：y = x²/2 ⇒ Y = X²/2；x_j = 1 → (1, 0.5)、x_j = −1.5 → (−1.5, 1.125)
  circle_envelope(f, fp, R, xs)（一般接觸轉換 S = (X−x)² + (Y−y)² = R² 的包絡＝平行曲線）＝ 到曲線距離恆為 R 的數值核對
  LAGR_CARD = {m: 2, qdot: 1}：L = q̇² ⇒ p = ∂L/∂q̇ = 2、H = p²/4 = 1、切線截距 −H = −1（(8.1.3)、Fig. 8.1.3）
§8.1B–C 作用量函數（H802）
  fall_state(theta, t)（初速 v0 = √(2ε/m) = 16 m/s、仰角 θ）、SP_closed(r0, r1, T)（(8.1.28)）＝ SP_initial(v0vec, y0, T)（(8.1.29)）＝ SP_numeric（∫L dt）
  SH_closed(theta, T)（(8.1.24)）＝ SH_numeric（∫p·dr），fall_checks()：∂S_P/∂T = −H、dS_P/dT = L（(8.1.30)；書 (8.1.30b)ₓ 的 2fy0T 應 2fẏ0T）、∂S_H/∂ε = T、
    (8.1.27) 括號應 (m/2)ẋ0²
  fall_card() → v0 16、頂點 32 m（T = 4 s）、45° 射程 64 m（T = 5.6569 s，S_H = 965.4365、S_P = 241.3591 J·s）、焦散 y = 32 − x²/128、
    垂直上拋到頂點 S_H = 1024/3 = 341.3333 J·s（所以 S_H ≥ 350 的波前碰到焦散成「蝙蝠頭」）
  trajectory(theta, n, tmax)、caustic(x)、wavefront(S, thetas)（S_H = 常數的點：每條軌跡解三次方程）、FIG816_LEVELS（50 … 1000）、
  fig816_curves(levels, n, ymin)（畫圖用 {S_H: [段]}）、
  orthogonality(S, thetas)（波前切線 · p ＝ 0，(8.1.19)）
§8.1D–E 半古典極限與 Huygens（H803）
  semiclassical_check()（sympy：ψ = e^{iS/ħ} 代入 ⇒ (8.1.38)，ħ → 0 得 HJ (8.1.17)）、wkb_ratio(n, x)（ħ|p'|/p²：n = 5、x = 1 → 0.0316；轉折點 √11 = 3.3166 發散）
  huygens_card()（自由空間 S_H = p|r|：最佳路徑 20、經 r'_10（偏 75.52°）的折線 10 + 20 = 30，Fig. 8.1.7）
§8.2 相干態（H804）
  ladder(N)、xp_ops(N)、disp_D(x0, N)、boost_B(p0, N)、Q_op(x0, p0, N)、alpha_of(x0, p0)、bch_check(x0, p0)（(8.2.12)–(8.2.15) 三種分解；照印 + 號不么正）
  coherent_cn(alpha, nmax)（(8.3.1b)）＝ Q|0⟩ 的數值、coherent_stats(alpha) → ⟨n⟩ = |α|²、Δn = |α|（(8.3.2)）
    |α|² = 5：c_0…c_10 = 0.0821 0.1835 0.2902 0.3747 0.4189 0.4189 0.3824 0.3232 0.2555 0.1904 0.1347；⟨n⟩ 5、Δn 2.2361
  ho_eigfun(n, x)（(8.2.7)，遞迴）、packet_closed(x, t, x0, p0)（(8.2.28)）＝ packet_sum（Σ c_n ψ_n e^{−i(n+½)t}）＝ packet_split_operator（FFT 數值演化）
  packet_moments(t, x0, p0) → x_t、p_t、Δx Δp = 0.5（(8.2.22)–(8.2.24)）、energies(alpha) → E_c = 5、E_q = 5.5（(8.2.25)–(8.2.26)，書印 E_0）、
  SH_ho(t, x0, p0)（(8.2.31)）＝ 數值 ∫p dx
§8.3 快照合成本徵函數（H805）
  snapshot_sum(x, x0, m, N, sign)（x0 ＝ 初始位移 √(2|α|²) ＝ √10；(8.3.3)–(8.3.4)；sign = −1 是書上 (8.3.3) 的號）、snapshot_error(x0, m, N, sign) →
    L2 相對誤差 ＝ √(Σ c²_{m+kN})/c_m：|α|² = 5、m = 5、N = 10 → 0.0299（＝ c15/c5，(8.3.9)）、N = 12 → 0.0091（書：< 1%）；照印的號 → 1.4276
§8.3A 自相關、FTAC、EBK（H806）
  ac_ho(t, alpha)（閉式）＝ Σ|c_n|² e^{−i(n+½)t}、ftac(ts, ac, Es)（(8.3.15)；T ＝ ts 的長度）、ho_ftac_card() → 峰在 n + ½、最高兩峰 4.5、5.5（|c|² = 0.1755）
  ho_maslov() → 半古典 AC（MU 波包 × e^{iS_p}）的峰 E_FT = n、E_SC = n + ½：E_FT = E_SC − 2πα/τ（α = ½；書上 (8.3.22) 的號印反、(8.3.20) 少 α）
  ROTOR = {J: 10, A: 0.2, B: 0.4, C: 0.6, beta0: −30}（Fig. 8.3.3）、rotor_H()、rotor_exact() → K = 8 叢 53.1461、53.1489（書 53.146、53.149）、
    初始波包 ⟨H⟩ = 53.5000（J² 約定的古典值 50.000 ＝ 書上 E_classical）、rotor_localized(K) → 0.129 0.922 −0.347 0.106 −0.031 0.009（(8.3.32)；
    書上 K' = 2、0 印 −0.062、0.018，正好兩倍）
  rotor_sc(T, nt)（自旋相干態的 TDVP 半古典傳播＋作用量相位 S = ∫(J cos θ φ̇ − E)dt，(8.3.25)–(8.3.30)，T = 40）→
    FTAC 的 K = 8 峰 E_FT = 53.410、FT 態的 E_SC = 53.135（書 53.134）、⟨K'|E_FT(8)⟩ = 0.119 0.929 −0.341 0.076 −0.011 0.001
    （書 (8.3.31) 0.118 0.929 −0.339 0.075 −0.012 0.001）；精確 AC 的 FTAC 峰 53.147；rotor_ebk(K)（∮J_z dφ = 2πK，|J| = J + ½）→ 53.214
    （書上 Fig. 8.3.3(b) 的橫軸「Frequency」不是能量，7.94200 無法對應，見 ch08_formulas〔注〕）
§8.4A A·p 與 E·r（H807）
  gauge_check()（sympy：ψ^α = e^{iqA·r/ħ}ψ^β 把 (8.4.9) 變成 (8.4.15)；書 (8.4.14) 少 q 就不成立）、hamiltonian_check()（(8.4.5)：p = μv + qA）、
  boost_shift(E0, ω)（(8.4.22)：qE0/μω² = 1.758e−15 m）、ap_er_ratio(ω_NB, ω)（|A·p 矩陣元|/|E·r 矩陣元| = ω_NB/ω）
§8.4B 微擾與譜函數（H808）
  S_amp(Δ, τ)（(8.4.38a)）＝ 數值積分、spectral_card(τ = 4) → 峰 τ² = 16、零點 ±2π/τ = ±1.5708、面積 2πτ = 25.1327、
    第一旁帶 Δ = 2.2467（4.4934·2/τ）高 0.0472τ²（書：約 ±3π/τ = ±2.3562、0.045τ²）
  first_order_check()（弱驅動兩能階的數值解 ＝ (8.4.37) 兩個譜函數；(e^{iΔt} − 1)/Δ = iS）、golden_rule_check()（準連續帶的 P(t)/t ＝ 2π|V|²ρ）、
  ho_strengths(N) → f_{n0} = δ_{n1}（(8.4.48)）、h_strengths() → 1s→2p 0.4162、3p 0.0791、4p 0.0290、束縛和 0.5650、連續 0.4350、總和 1.0000
    ＝ h_box_strengths()（盒中徑向方程的數值對角化）；classical_factor_check()（|4(q/μ)E0|² 應 4|(q/μ)E0|²）
§8.5 兩能階（H809–H811）
  TL_CARD = {eps: 1, r: 0.2}（本系列）、crank(t, r, Ω, ε)（(8.5.6)）、h_lab(t, r, Ω, ε, rwa)（(8.5.5)／(8.5.4a)）、evolve_lab(...)（solve_ivp）、
  spin_precession(...)（dS/dt = ω × S 數值 ＝ SC.spin_vector(ψ)）、rotating_check()（(8.5.11)–(8.5.15)）
  rabi(Δ, r, t)（(8.5.17)）＝ rabi_numeric、rabi_card() → Δ = 0：t = π/r = 15.7080 時 1；Δ = r：最大 0.5000（t = 11.1072）；包絡半高 Δ = ±r（功率展寬 FWHM 2r = 0.4）
  floquet_population(Ω, ε, r)、bloch_siegert(ε, r)（全 Hamilton 量 (8.5.4a) 的 Floquet 共振：兩模時間平均布居 ½）→ ε = 1、r = 0.2：
    位移 0.010025（＝ r²/4ε + r⁴/64ε³；卡上 0.0100 ＝ r²/4ε；書 (8.5.25) 的啟發式 0.0099）、(8.5.27) 分數 r²/4ε²（書例 2.50e−13）
  bloch_rhs、bloch_steady(Δ, r, γ1, γ2, Sz0)（解析 ＝ 長時間數值）→ BLOCH_CARD：Δ = 0、r = 0.2、γ1 = 0.05、γ2 = 0.1：S_z = −0.0556、上能階布居 0.4444（飽和參數 r²/γ1γ2 = 8）
  thermal_sz(βε)（βε = 2 → −0.3808）、dressed(Δ, r)（(8.5.39)–(8.5.44)）、ac_stark(Ω, ε, r)（(8.5.45)：Fig. 8.5.3 的 ε = 1、r = 0.5：共振分裂 r = 0.5（圖上印 2r）、
    AC 位移 δ/2 = 0.25、DC 位移 δ(ε)/2 = 0.0590）、circle_geometry(Δ, r)（Fig. 8.5.4，r = 0.5、Δ = 1.2、0.2、−0.1、−1.2）、
  adiabatic(r, D0, T)（線性掃 Δ：數值反轉 ＝ 1 − e^{−πr²/2v}；r = 0.5、Δ 20 → −20、T = 400（v = 0.1）→ 0.98）
§8.6 量子化的場（H812）
  photon_elements(n)（√n、√(n+1)）、mode_energy_check()（(8.6.3)：2ε0ω²|a|²V，磁場項要 k²）、angular_sum()（∫dΩ Σ_α|e_α·r̂|² = 8π/3）、
  einstein_A(ω, r2)（SI）＝ einstein_A_gauss（cgs）＝ einstein_A_au；h2p_A() → |⟨1s|z|2p0⟩|² = 2¹⁵/3¹⁰ = 0.5549 a0²、A = 6.2683e8 s⁻¹
    （核質量無窮大；乘約化質量因子 → 6.2649e8 ＝ 文獻值）、τ = 1.595 ns；照印 (8.6.31) 少 1/2π → 3.9385e9；golden_A_numeric()（(8.6.28) 對 ω 積分 ＝ A）、number_vs_coherent()（⟨n|a|n⟩ = 0、⟨α|a|α⟩ = α）
§8.7 Jaynes–Cummings（H813）
  jc_matrix(Nmax, Ω, ε, g)、jc_block(N, Ω, ε, g)、jc_levels(...)（數值對角化 ＝ 2×2 塊公式）、jc_splitting(N, g) → g = 0.1：N = 1, 4, 9, 16 → 0.1, 0.2, 0.3, 0.4
  JC_CARD = {Omega: 10, Delta: 0.3, r: 0.4}、mollow(Ω, Δ, r) → Ω_F = 10.5、Ω_L = 10、Ω_3 = 9.5（(8.7.11)、(8.7.13)；遠離共振 ε = 10.3、2Ω − ε = 9.7；
    照印的 −δ/2 給 Ω − Δ = 9.7）、mollow_jc(N, g)（N = 100、g = 0.04 的 JC 能階差：10.4992、9.5008）、jc_ladder(Δ, g, Nmax)（Fig. 8.7.1）
全書回顧：RECAP（八章各一列：主題、關鍵結果、動畫範圍）、recap_check()
card_values() ＝ CARD_EXPECTED（plan_ch08 §6 卡片上與導讀裡的每個數字）；check()：每個結果兩種以上獨立算法（約 20 秒）；summary()。
直接執行：python utils/sym_dyn.py（在 codes/ 裡）印一句成功訊息；--summary 印卡片數字；--freeze 印 CARD_EXPECTED。
"""
import math
import warnings

import numpy as np
import sympy as sp
from scipy.integrate import quad, solve_ivp
from scipy.linalg import eigh_tridiagonal, expm
from scipy.optimize import brentq, minimize_scalar

try:
    from . import sym_coup as SC
    from . import sym_mat as SM          # noqa: F401（給場景用：DY.SM）
    from . import sym_rot as RT
except ImportError:                      # 直接執行 python utils/sym_dyn.py
    import sym_coup as SC
    import sym_mat as SM                 # noqa: F401
    import sym_rot as RT

PI = math.pi


def _f(x, d=4) -> str:
    return f"{float(x):.{d}f}"


def _e(x, d=4) -> str:
    return f"{float(x):.{d}e}"


def _close(a, b, tol=1e-9) -> bool:
    return bool(np.allclose(np.asarray(a, complex), np.asarray(b, complex), atol=tol, rtol=0))


# ======================================================================== §8.1A 接觸轉換、Legendre
def legendre_numeric(f, Xs, xs):
    """(8.1.2c)：Y(X) = max_x [xX − f(x)]（凸的 f；xs 是細格點）。"""
    Xs, xs = np.asarray(Xs, float), np.asarray(xs, float)
    return np.max(np.outer(Xs, xs) - f(xs)[None, :], axis=1)


def line_family(f, xjs) -> list:
    """Fig. 8.1.2：每點 (x_j, y_j) 生成直線 Y = x_j X − y_j → [(斜率, 截距)]。"""
    return [(float(x), float(-f(x))) for x in xjs]


def line_envelope(f, fp, xjs) -> np.ndarray:
    """(8.1.2a,b)：切點 X_j = f'(x_j)、Y_j = x_j X_j − f(x_j)（f(x) = x²/2 時 X = x）。"""
    xjs = np.asarray(xjs, float)
    X = fp(xjs)
    return np.stack([X, xjs * X - f(xjs)], axis=1)


LINE_CARD = dict(f=lambda x: x ** 2 / 2, fp=lambda x: x, xjs=(1.0, -1.5))


def circle_envelope(f, fp, R: float, xs) -> np.ndarray:
    """一般接觸轉換 S(x,y:X,Y) = (X−x)² + (Y−y)² = R²（每點生成半徑 R 的圓）的包絡＝平行曲線（取上側）。"""
    xs = np.asarray(xs, float)
    s = fp(xs)
    nrm = np.sqrt(1 + s ** 2)
    return np.stack([xs - R * s / nrm, f(xs) + R / nrm], axis=1)


def _dist_to_curve(P, f, lo=-6.0, hi=6.0) -> float:
    g = lambda x: (x - P[0]) ** 2 + (f(x) - P[1]) ** 2
    xs = np.linspace(lo, hi, 4001)
    i = int(np.argmin(g(xs)))
    r = minimize_scalar(g, bounds=(xs[max(i - 2, 0)], xs[min(i + 2, len(xs) - 1)]), method="bounded", options=dict(xatol=1e-12))
    return math.sqrt(r.fun)


LAGR_CARD = dict(m=2.0, qdot=1.0)          # 本系列：L = (m/2)q̇²，m = 2


def lagrange_card() -> dict:
    m, qd = LAGR_CARD["m"], LAGR_CARD["qdot"]
    L = lambda v: m * v ** 2 / 2
    p = m * qd                                                     # (8.1.3a)
    H = p * qd - L(qd)                                             # (8.1.3c)
    Hn = float(legendre_numeric(L, [p], np.linspace(-5, 5, 200001))[0])
    return dict(p=p, H=H, H_numeric=Hn, H_formula=p ** 2 / (2 * m), intercept=-H)


# ======================================================================== §8.1B–C 作用量、拋射體
FALL = dict(eps=128.0, F=-4.0, m=1.0)       # 書上 ε = 128 J、F = −4 N；質量本系列取 1 kg
FIG816_LEVELS = (50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 600, 700, 800, 1000)


def fall_v0() -> float:
    return math.sqrt(2 * FALL["eps"] / FALL["m"])


def fall_state(theta_deg: float, t, v0=None):
    """從原點以仰角 θ（度）、速率 v0 射出：(x, y, ẋ, ẏ)（(8.1.26)）。"""
    v0 = fall_v0() if v0 is None else v0
    m, F = FALL["m"], FALL["F"]
    th = math.radians(theta_deg)
    t = np.asarray(t, float)
    vx, vy0 = v0 * math.cos(th), v0 * math.sin(th)
    return vx * t, vy0 * t + F / (2 * m) * t ** 2, vx + 0 * t, vy0 + F / m * t


def SP_closed(r0, r1, T: float) -> float:
    """(8.1.28)：S_P = m(Δx)²/2T + m(Δy)²/2T + (fT/2)(y1 + y0) − f²T³/24m。"""
    m, f = FALL["m"], FALL["F"]
    dx, dy = r1[0] - r0[0], r1[1] - r0[1]
    return m * dx ** 2 / (2 * T) + m * dy ** 2 / (2 * T) + f * T / 2 * (r1[1] + r0[1]) - f ** 2 * T ** 3 / (24 * m)


def SP_initial(v0vec, y0: float, T: float) -> float:
    """(8.1.29)：S_P = mẋ0²T/2 + mẏ0²T/2 + fTy0 + fT²ẏ0 + f²T³/3m。"""
    m, f = FALL["m"], FALL["F"]
    return m * v0vec[0] ** 2 * T / 2 + m * v0vec[1] ** 2 * T / 2 + f * T * y0 + f * T ** 2 * v0vec[1] + f ** 2 * T ** 3 / (3 * m)


def SP_numeric(theta_deg: float, T: float) -> float:
    m, f = FALL["m"], FALL["F"]

    def L(t):
        x, y, vx, vy = fall_state(theta_deg, t)
        return m * (vx ** 2 + vy ** 2) / 2 + f * y
    return quad(L, 0, T, epsabs=1e-12, epsrel=1e-12)[0]


def SH_closed(theta_deg: float, T: float) -> float:
    """(8.1.24)：s_x = mẋ0(x1 − x0)、s_y = (m²/3f)(ẏ1³ − ẏ0³)。"""
    m, f = FALL["m"], FALL["F"]
    x1, y1, vx, vy1 = fall_state(theta_deg, T)
    vy0 = float(fall_state(theta_deg, 0.0)[3])
    return float(m * vx * x1 + m ** 2 / (3 * f) * (vy1 ** 3 - vy0 ** 3))


def SH_numeric(theta_deg: float, T: float) -> float:
    """∫ p·dr = ∫ m v² dt。"""
    m = FALL["m"]
    g = lambda t: m * (fall_state(theta_deg, t)[2] ** 2 + fall_state(theta_deg, t)[3] ** 2)
    return quad(g, 0, T, epsabs=1e-12, epsrel=1e-12)[0]


def _SH_cubic(theta_deg: float, t):
    """S_H(t) = m[v0²t + ẏ0(f/m)t² + (f/m)²t³/3]（沿軌跡，單調遞增）。"""
    m, f = FALL["m"], FALL["F"]
    v0 = fall_v0()
    vy0 = v0 * math.sin(math.radians(theta_deg))
    g = f / m
    return m * (v0 ** 2 * t + vy0 * g * t ** 2 + g ** 2 * t ** 3 / 3)


def trajectory(theta_deg: float, n: int = 60, tmax: float = 8.0) -> np.ndarray:
    ts = np.linspace(0, tmax, n)
    x, y, _, _ = fall_state(theta_deg, ts)
    return np.stack([x, y], axis=1)


def caustic(x):
    """拋物線焦散（安全拋物線）y = ε/|F| − |F|x²/4ε = 32 − x²/128。"""
    eps, F = FALL["eps"], abs(FALL["F"])
    return eps / F - F * np.asarray(x, float) ** 2 / (4 * eps)


def wavefront(Sval: float, thetas) -> np.ndarray:
    """S_H(0:r) = Sval 的點（Fig. 8.1.6(b)）：沿每條仰角 θ 的軌跡解 S_H(t) = Sval。→ (n, 2)。"""
    pts = []
    for th in thetas:
        t = brentq(lambda tt: _SH_cubic(th, tt) - Sval, 0.0, 100.0, xtol=1e-14)
        x, y, _, _ = fall_state(th, t)
        pts.append((float(x), float(y)))
    return np.array(pts)


def orthogonality(Sval: float, thetas, h: float = 1e-4) -> float:
    """波前切線（對仰角的中心差分）與動量的夾角餘弦的最大絕對值（(8.1.19)：p·dr = 0 ⇒ 0）。"""
    worst = 0.0
    for th in thetas:
        a, b = wavefront(Sval, [th - h])[0], wavefront(Sval, [th + h])[0]
        tan_ = b - a
        t = brentq(lambda tt: _SH_cubic(th, tt) - Sval, 0.0, 100.0, xtol=1e-14)
        _, _, vx, vy = fall_state(th, t)
        v = np.array([float(vx), float(vy)])
        worst = max(worst, abs(tan_ @ v) / (np.linalg.norm(tan_) * np.linalg.norm(v)))
    return worst


def fig816_curves(levels=FIG816_LEVELS, n: int = 60, ymin: float = -8.0) -> dict:
    """畫圖用：{S_H: [段…]}，每段是 (k, 2) 的波前折線（仰角 −90°…270° 取 n 個；y < ymin 的部分切掉，所以一個 level 可能有好幾段）。"""
    ths = np.linspace(-89.5, 269.5, n)
    out = {}
    for S_ in levels:
        pts = wavefront(S_, ths)
        segs, cur = [], []
        for p in pts:
            if p[1] >= ymin:
                cur.append(p)
            elif cur:
                segs.append(np.array(cur))
                cur = []
        if cur:
            segs.append(np.array(cur))
        out[S_] = segs
    return out


def fall_checks() -> dict:
    """(8.1.28) ＝ (8.1.29) ＝ ∫L dt；(8.1.24) ＝ ∫p·dr；∂S_P/∂T = −H、dS_P/dT = L、∂S_H/∂ε = T；書上 (8.1.27)、(8.1.30b)ₓ 的錯。"""
    m, f, eps = FALL["m"], FALL["F"], FALL["eps"]
    out = {}
    th, T = 37.0, 3.3
    x1, y1, vx, vy1 = [float(u) for u in fall_state(th, T)]
    vy0 = float(fall_state(th, 0.0)[3])
    sp1 = SP_closed((0, 0), (x1, y1), T)
    sp2 = SP_initial((vx, vy0), 0.0, T)
    sp3 = SP_numeric(th, T)
    out["SP"] = abs(sp1 - sp2) < 1e-9 and abs(sp1 - sp3) < 1e-8
    out["SH"] = abs(SH_closed(th, T) - SH_numeric(th, T)) < 1e-8 and abs(SH_closed(th, T) - (sp1 + eps * T)) < 1e-8   # (8.1.25)
    h = 1e-5
    dSdT = (SP_closed((0, 0), (x1, y1), T + h) - SP_closed((0, 0), (x1, y1), T - h)) / (2 * h)
    out["dSP_dT_partial"] = abs(dSdT + eps) < 1e-5                                       # (8.1.30a)
    tot = (SP_numeric(th, T + h) - SP_numeric(th, T - h)) / (2 * h)
    Lval = m * (vx ** 2 + vy1 ** 2) / 2 + f * y1
    out["dSP_dT_total"] = abs(tot - Lval) < 1e-5                                         # (8.1.30b)
    vx0 = vx
    book_30b = m * vx0 ** 2 / 2 + m * vy0 ** 2 / 2 + f * 0.0 + 2 * f * 0.0 * T + f ** 2 * T ** 2 / m   # 照印 2f y0 T（y0 = 0 時為 0）
    right_30b = m * vx0 ** 2 / 2 + m * vy0 ** 2 / 2 + f * 0.0 + 2 * f * vy0 * T + f ** 2 * T ** 2 / m
    out["30b_right"] = abs(right_30b - Lval) < 1e-9 and abs(book_30b - Lval) > 1.0
    # (8.1.27) 第一行的括號 = ε：照印 (m/2)ẋ0 不是能量
    first = m * vx0 * x1 + m ** 2 / (3 * f) * ((vy0 + f * T / m) ** 3 - vy0 ** 3) - (m / 2 * vx0 ** 2 + m / 2 * vy0 ** 2 - f * 0.0) * T
    first_book = m * vx0 * x1 + m ** 2 / (3 * f) * ((vy0 + f * T / m) ** 3 - vy0 ** 3) - (m / 2 * vx0 + m / 2 * vy0 ** 2 - f * 0.0) * T
    out["27_right"] = abs(first - sp1) < 1e-8 and abs(first_book - sp1) > 1.0
    # ∂S_H/∂ε = T（(8.1.31)）：固定終點、改變能量，S_H(ε) 用兩個能量的閉式
    out["dSH_deps"] = _dSH_deps_check()
    return out


def _SH_endpoint(eps: float, r1, high: bool = False) -> tuple:
    """固定能量 ε、終點 r1：解出經過 r1 的仰角（低拋或高拋），回傳 (S_H, T)。"""
    m, f = FALL["m"], FALL["F"]
    v0 = math.sqrt(2 * eps / m)
    g = -f / m

    def miss(th):
        c = math.cos(th)
        T = r1[0] / (v0 * c)
        return v0 * math.sin(th) * T - g * T ** 2 / 2 - r1[1]
    lo, hi = (math.atan2(r1[1], r1[0]) + 1e-9, PI / 4 + math.atan2(r1[1], r1[0]) / 2) if not high else \
        (PI / 4 + math.atan2(r1[1], r1[0]) / 2, PI / 2 - 1e-9)
    th = brentq(miss, lo, hi, xtol=1e-15)
    T = r1[0] / (v0 * math.cos(th))
    vx, vy0 = v0 * math.cos(th), v0 * math.sin(th)
    vy1 = vy0 - g * T
    return m * vx * r1[0] + m ** 2 / (3 * f) * (vy1 ** 3 - vy0 ** 3), T


def _dSH_deps_check() -> bool:
    r1, eps, h = (20.0, 5.0), 128.0, 1e-4
    a, T = _SH_endpoint(eps, r1)
    b, _ = _SH_endpoint(eps + h, r1)
    c, _ = _SH_endpoint(eps - h, r1)
    return abs((b - c) / (2 * h) - T) < 1e-6


def fall_card() -> dict:
    eps, F, m = FALL["eps"], FALL["F"], FALL["m"]
    v0 = fall_v0()
    g = -F / m
    T45 = 2 * v0 * math.sin(PI / 4) / g
    sh45 = SH_closed(45.0, T45)
    apex_SH = SH_closed(90.0, v0 / g)
    return dict(v0=v0, apex=v0 ** 2 / (2 * g), T_apex=v0 / g, range=v0 ** 2 / g, T45=T45, SH45=sh45, SP45=sh45 - eps * T45,
                SH45_num=SH_numeric(45.0, T45), apex_SH=apex_SH, apex_SH_num=SH_numeric(90.0, v0 / g),
                caustic_c=eps / abs(F), caustic_k=4 * eps / abs(F))


# ======================================================================== §8.1D–E 半古典極限、Huygens
def semiclassical_check() -> dict:
    """ψ = ψ0 e^{iS/ħ} 代入 iħψ_t = −ħ²ψ_xx/2m + Vψ：(1/ψ)(…) = −S_t − S_x²/2m − V + iħS_xx/2m（(8.1.38)）；ħ → 0 得 HJ (8.1.17)。"""
    x, t = sp.symbols("x t", real=True)
    hb, m = sp.symbols("hbar m", positive=True)
    Sf = sp.Function("S")(x, t)
    V = sp.Function("V")(x)
    psi = sp.exp(sp.I * Sf / hb)
    res = sp.simplify((sp.I * hb * sp.diff(psi, t) + hb ** 2 / (2 * m) * sp.diff(psi, x, 2) - V * psi) / psi)
    book = -sp.diff(Sf, t) - sp.diff(Sf, x) ** 2 / (2 * m) - V + sp.I * hb / (2 * m) * sp.diff(Sf, x, 2)
    ok = sp.simplify(res - book) == 0
    hj = sp.limit(res, hb, 0)
    return dict(ok=bool(ok), hj_ok=bool(sp.simplify(hj - (-sp.diff(Sf, t) - sp.diff(Sf, x) ** 2 / (2 * m) - V)) == 0))


def wkb_ratio(n: int, x: float) -> float:
    """ħ|dp/dx|/p²（(8.1.39)）：諧振子 p = √(2(n+½) − x²)（ħ = μ = ω = 1）。"""
    p = math.sqrt(2 * (n + 0.5) - x ** 2)
    return abs(-x / p) / p ** 2


def huygens_card() -> dict:
    """Fig. 8.1.7 的數字版（自由空間 S_H = p|r|、p = 1）：r_10 = (10, 0)、r_20 = (20, 0)；r'_10 在 10 號圓上偏 θ 使 |r'_10 − r_20| = 20。"""
    th = math.acos(0.25)
    rp = np.array([10 * math.cos(th), 10 * math.sin(th)])
    dev = 10 + float(np.linalg.norm(rp - np.array([20.0, 0.0])))
    return dict(theta=math.degrees(th), optimal=20.0, deviant=dev)


# ======================================================================== §8.2 相干態
NBASIS = 90


def ladder(N: int = NBASIS) -> np.ndarray:
    """截斷的 a：a[n−1, n] = √n（(8.2.4b)）。"""
    return np.diag(np.sqrt(np.arange(1, N)), 1).astype(complex)


def xp_ops(N: int = NBASIS):
    a = ladder(N)
    ad = a.conj().T
    return (a + ad) / math.sqrt(2), 1j * (ad - a) / math.sqrt(2)


def disp_D(x0: float, N: int = NBASIS) -> np.ndarray:
    """(8.2.8)：D(x0) = e^{−ix0 p}。"""
    return expm(-1j * x0 * xp_ops(N)[1])


def boost_B(p0: float, N: int = NBASIS) -> np.ndarray:
    """(8.2.10)：B(p0) = e^{ip0 x}。"""
    return expm(1j * p0 * xp_ops(N)[0])


def Q_op(x0: float, p0: float, N: int = NBASIS) -> np.ndarray:
    """(8.2.11)：Q = e^{i(p0 x − x0 p)}。"""
    x, p = xp_ops(N)
    return expm(1j * (p0 * x - x0 * p))


def alpha_of(x0: float, p0: float) -> complex:
    """(8.2.15b)：α = (x0 + ip0)/√2。"""
    return (x0 + 1j * p0) / math.sqrt(2)


def bch_check(x0: float = 1.3, p0: float = -0.7, N: int = NBASIS, k: int = 25) -> dict:
    """(8.2.12)–(8.2.15)：Q ＝ e^{−ix0p0/2}BD ＝ e^{ix0p0/2}DB ＝ e^{αa† − α*a} ＝ e^{−|α|²/2}e^{αa†}e^{−α*a}（比左上 k×k，避開截斷邊緣）；
    照印的 e^{αa† + α*a} 不是么正。"""
    a = ladder(N)
    ad = a.conj().T
    al = alpha_of(x0, p0)
    Qm = Q_op(x0, p0, N)
    Dm, Bm = disp_D(x0, N), boost_B(p0, N)
    f1 = np.exp(-1j * x0 * p0 / 2) * Bm @ Dm
    f2 = np.exp(1j * x0 * p0 / 2) * Dm @ Bm
    f3 = expm(al * ad - np.conj(al) * a)
    f4 = np.exp(-abs(al) ** 2 / 2) * expm(al * ad) @ expm(-np.conj(al) * a)
    book = expm(al * ad + np.conj(al) * a)
    s = slice(0, k)
    ok = all(_close(Qm[s, s], F_[s, s], 1e-8) for F_ in (f1, f2, f3, f4))
    unit_book = np.linalg.norm(book[:k, :k] @ book[:k, :k].conj().T - np.eye(k))
    return dict(ok=ok, book_not_unitary=unit_book > 1.0)


def coherent_cn(alpha: complex, nmax: int) -> np.ndarray:
    """(8.3.1b)：c_n = e^{−|α|²/2} αⁿ/√n!。"""
    n = np.arange(nmax + 1)
    lg = np.array([math.lgamma(k + 1) for k in n])
    return np.exp(-abs(alpha) ** 2 / 2) * np.exp(n * np.log(abs(alpha) + 1e-300) - lg / 2) * np.exp(1j * n * np.angle(alpha)) if alpha != 0 \
        else (n == 0).astype(complex)


def coherent_numeric(x0: float, p0: float, N: int = NBASIS) -> np.ndarray:
    v0 = np.zeros(N, complex)
    v0[0] = 1
    return Q_op(x0, p0, N) @ v0


def coherent_stats(alpha: complex, N: int = NBASIS) -> dict:
    c = coherent_cn(alpha, N - 1)
    nn = np.arange(N)
    pr = abs(c) ** 2
    mean = float(pr @ nn)
    var = float(pr @ nn ** 2) - mean ** 2
    return dict(mean=mean, dn=math.sqrt(var), norm=float(pr.sum()))


def ho_eigfun(n: int, x) -> np.ndarray:
    """(8.2.7)：ψ_n(x)（μ = ω = ħ = 1），用穩定的遞迴 ψ_{k+1} = √(2/(k+1)) x ψ_k − √(k/(k+1)) ψ_{k−1}。"""
    x = np.asarray(x, float)
    p0 = np.pi ** -0.25 * np.exp(-x ** 2 / 2)
    if n == 0:
        return p0
    p1 = math.sqrt(2) * x * p0
    for k in range(1, n):
        p0, p1 = p1, math.sqrt(2 / (k + 1)) * x * p1 - math.sqrt(k / (k + 1)) * p0
    return p1


def packet_moments(t: float, x0: float, p0: float = 0.0) -> dict:
    """(8.2.22)：x_t = x0 cos t + p0 sin t、p_t = p0 cos t − x0 sin t。"""
    return dict(x=x0 * math.cos(t) + p0 * math.sin(t), p=p0 * math.cos(t) - x0 * math.sin(t))


def packet_closed(x, t: float, x0: float, p0: float = 0.0, sign: int = 1) -> np.ndarray:
    """(8.2.28)：χ = e^{ix_tp_t/2 − it/2} N0 e^{−(x−x_t)²/2 + i·sign·p_t(x−x_t)}（sign = −1 是書上 (8.3.3) 的號）。"""
    mo = packet_moments(t, x0, p0)
    xt, pt = mo["x"], mo["p"]
    x = np.asarray(x, float)
    return np.exp(1j * xt * pt / 2 - 1j * t / 2) * np.pi ** -0.25 * np.exp(-(x - xt) ** 2 / 2 + 1j * sign * pt * (x - xt))


def packet_sum(x, t: float, x0: float, p0: float = 0.0, nmax: int = 60) -> np.ndarray:
    """(8.3.1a)：Σ c_n ψ_n(x) e^{−i(n+½)t}。"""
    c = coherent_cn(alpha_of(x0, p0), nmax)
    return sum(c[n] * ho_eigfun(n, x) * np.exp(-1j * (n + 0.5) * t) for n in range(nmax + 1))


def packet_split_operator(t: float, x0: float, p0: float = 0.0, L: float = 24.0, n: int = 1024, steps: int = 400):
    """FFT 分裂算子數值演化（獨立方法）：從 χ(x, 0) 出發 → (x, χ(x, t))。"""
    x = np.linspace(-L / 2, L / 2, n, endpoint=False)
    dx = x[1] - x[0]
    k = 2 * np.pi * np.fft.fftfreq(n, d=dx)
    psi = packet_closed(x, 0.0, x0, p0)
    dt = t / steps
    half = np.exp(-1j * x ** 2 / 2 * dt / 2)
    kin = np.exp(-1j * k ** 2 / 2 * dt)
    for _ in range(steps):
        psi = half * psi
        psi = np.fft.ifft(kin * np.fft.fft(psi))
        psi = half * psi
    return x, psi


def packet_uncertainty(t: float, x0: float, p0: float = 0.0) -> float:
    """數值 ΔxΔp（在 x 格點上直接積分）＝ ½。"""
    x = np.linspace(-14, 14, 4001)
    dx = x[1] - x[0]
    psi = packet_closed(x, t, x0, p0)
    rho = abs(psi) ** 2
    mx = np.sum(x * rho) * dx
    vx = np.sum((x - mx) ** 2 * rho) * dx
    k = 2 * np.pi * np.fft.fftfreq(len(x), d=dx)               # 動量分布用 FFT（譜方法，比差分準）
    pk = abs(np.fft.fft(psi)) ** 2
    pk /= pk.sum()
    mp = float(pk @ k)
    return math.sqrt(vx * (float(pk @ k ** 2) - mp ** 2))


def energies(alpha: complex) -> dict:
    """(8.2.25)–(8.2.26)：E_c = |α|²、E_q = E_c + ½（書印 E_0 + ħω/2）。"""
    st = coherent_stats(alpha)
    return dict(Ec=abs(alpha) ** 2, Eq=st["mean"] + 0.5)


def SH_ho(t: float, x0: float, p0: float = 0.0) -> float:
    """(8.2.31)：∫p_t dx_t = (x_tp_t − x0p0)/2 + E_c t。"""
    mo = packet_moments(t, x0, p0)
    Ec = (x0 ** 2 + p0 ** 2) / 2
    return (mo["x"] * mo["p"] - x0 * p0) / 2 + Ec * t


def SH_ho_numeric(t: float, x0: float, p0: float = 0.0) -> float:
    g = lambda s: packet_moments(s, x0, p0)["p"] * (-x0 * math.sin(s) + p0 * math.cos(s))
    return quad(g, 0, t, epsabs=1e-12)[0]


COH_CARD = dict(alpha2=5.0)        # 書上 Fig. 8.2.2(b)、8.3.2：|α|² = 5（x0 = √10，p0 = 0）


def coh_x0() -> float:
    return math.sqrt(2 * COH_CARD["alpha2"])


# ======================================================================== §8.3 快照合成
def snapshot_sum(x, x0: float, m: int, N: int, sign: int = 1) -> np.ndarray:
    """(8.3.3)–(8.3.4)：(1/N) Σ_p χ(x, t_p) e^{i(m+½)t_p}，t_p = 2πp/N（sign = −1 用書上 (8.3.3) 的 −ip_t 號）。"""
    acc = 0
    for p in range(N):
        tp = 2 * PI * p / N
        acc = acc + packet_closed(x, tp, x0, 0.0, sign) * np.exp(1j * (m + 0.5) * tp)
    return acc / N


def snapshot_error(x0: float, m: int, N: int, sign: int = 1) -> dict:
    """L2 相對誤差 ‖Σ − c_mψ_m‖/‖c_mψ_m‖（格點積分）與 (8.3.8) 的預測 √(Σ_k c²_{m+kN})/c_m、最大誤差／峰值。"""
    x = np.linspace(-12, 12, 6001)
    dx = x[1] - x[0]
    c = coherent_cn(alpha_of(x0, 0.0), 120)
    target = c[m] * ho_eigfun(m, x)
    sm = snapshot_sum(x, x0, m, N, sign)
    l2 = math.sqrt(np.sum(abs(sm - target) ** 2) * dx / (np.sum(abs(target) ** 2) * dx))
    pred = math.sqrt(sum(abs(c[m + k * N]) ** 2 for k in range(1, 12) if m + k * N <= 120)) / abs(c[m])
    return dict(l2=l2, pred=pred, ratio=abs(c[m + N] / c[m]), maxerr=float(np.max(abs(sm - target)) / np.max(abs(target))))


# ======================================================================== §8.3A 自相關、FTAC、EBK
def ac_ho(t, alpha: complex):
    """⟨α|e^{−iHt}|α⟩ = e^{−it/2} exp(|α|²(e^{−it} − 1))。"""
    t = np.asarray(t, float)
    return np.exp(-1j * t / 2) * np.exp(abs(alpha) ** 2 * (np.exp(-1j * t) - 1))


def ac_ho_sum(t, alpha: complex, nmax: int = 60):
    c = coherent_cn(alpha, nmax)
    t = np.asarray(t, float)
    return sum(abs(c[n]) ** 2 * np.exp(-1j * (n + 0.5) * t) for n in range(nmax + 1))


def ftac(ts, ac, Es) -> np.ndarray:
    """(8.3.15)：(1/2T)∫_{−T}^{T} e^{iEt} AC(t) dt，AC(−t) = AC(t)*；ts 從 0 到 T 等間隔（梯形）。"""
    ts = np.asarray(ts, float)
    T = ts[-1]
    w = np.full(len(ts), ts[1] - ts[0])
    w[0] = w[-1] = w[0] / 2
    Es = np.asarray(Es, float)
    ph = np.exp(1j * np.outer(Es, ts))
    return np.real(ph @ (w * ac)) / T


def _peaks(Es, F, lo, hi) -> list:
    Es, F = np.asarray(Es), np.asarray(F)
    out = []
    for i in range(1, len(F) - 1):
        if lo <= Es[i] <= hi and F[i] > F[i - 1] and F[i] >= F[i + 1]:
            out.append((float(Es[i]), float(F[i])))
    return out


def _refine_peak(fun, e0: float, h: float, n: int = 241) -> tuple:
    """在 [e0 − h, e0 + h] 先用細格點找最高點（FTAC 的峰寬約 1/T，旁瓣很多），再在一格內精確化。"""
    es = np.linspace(e0 - h, e0 + h, n)
    vals = [fun(e) for e in es]
    i = int(np.argmax(vals))
    step = es[1] - es[0]
    r = minimize_scalar(lambda e: -fun(e), bounds=(es[i] - step, es[i] + step), method="bounded", options=dict(xatol=1e-10))
    return float(r.x), float(-r.fun)


def ho_ftac_card(T: float = 20 * PI) -> dict:
    al = math.sqrt(COH_CARD["alpha2"])
    ts = np.linspace(0, T, 8001)
    ac = ac_ho(ts, al)
    fun = lambda e: float(ftac(ts, ac, [e])[0])
    peaks = [_refine_peak(fun, n + 0.5, 0.015, 31) for n in range(0, 12)]       # 主峰半寬 1/T ≈ 0.016
    return dict(peaks=peaks, ac_pi=float(abs(ac_ho(PI, al))), ac_2pi=float(abs(ac_ho(2 * PI, al))))


def ho_maslov(T: float = 20 * PI) -> dict:
    """半古典 AC：|MU(t)⟩ = 平移的 Gauss（(8.2.28) 的大括號）× e^{iS_p}、S_p = ∫(p ẋ − E_c)dt ⇒ AC_sc = AC·e^{it/2}，
    FTAC 的峰在 E_FT = n（整數）；FT 態的 ⟨H⟩ = E_SC = n + ½ ⇒ E_FT = E_SC − 2πα/τ（α = ½、τ = 2π）。"""
    x0 = coh_x0()
    al = alpha_of(x0, 0.0)
    ts = np.linspace(0, T, 8001)
    # MU 波包的重疊直接在 x 格點上算（獨立於閉式），作用量用 (8.2.31) − E_c t
    xg = np.linspace(-12, 12, 1201)
    dx = xg[1] - xg[0]
    mu0 = np.pi ** -0.25 * np.exp(-(xg - x0) ** 2 / 2)
    acs = []
    for t in ts[::20]:
        mo = packet_moments(t, x0)
        mut = np.pi ** -0.25 * np.exp(-(xg - mo["x"]) ** 2 / 2 + 1j * mo["p"] * (xg - mo["x"]))
        Sp = SH_ho(t, x0) - abs(al) ** 2 * t
        acs.append(np.sum(np.conj(mu0) * mut) * dx * np.exp(1j * Sp))
    acs = np.array(acs)
    tsub = ts[::20]
    fun = lambda e: float(ftac(tsub, acs, [e])[0])
    eft = _refine_peak(fun, 5.0, 0.3)[0]
    # E_SC：FT 態 ∫ e^{iE_FT t}|ψ_sc(t)⟩ dt 在數態基底（ψ_sc = e^{iS_p}·e^{−ix_tp_t/2}|α_t⟩）
    nmax = 60
    acc = np.zeros(nmax + 1, complex)
    w = np.full(len(ts), ts[1] - ts[0])
    w[0] = w[-1] = w[0] / 2
    for t, wt in zip(np.concatenate([-ts[::-1], ts[1:]]), np.concatenate([w[::-1], w[1:]])):
        mo = packet_moments(t, x0)
        Sp = SH_ho(t, x0) - abs(al) ** 2 * t
        vec = coherent_cn(al * np.exp(-1j * t), nmax) * np.exp(-1j * mo["x"] * mo["p"] / 2) * np.exp(1j * Sp)
        acc += wt * np.exp(1j * eft * t) * vec
    pr = abs(acc) ** 2
    esc = float(pr @ (np.arange(nmax + 1) + 0.5) / pr.sum())
    return dict(E_FT=eft, E_SC=esc, maslov=(esc - eft) * 2 * PI / (2 * PI), ac_ratio_ok=_close(acs, ac_ho(tsub, al) * np.exp(1j * tsub / 2), 1e-6))


ROTOR = dict(J=10, A=0.2, B=0.4, C=0.6, beta0=-30.0)       # Fig. 8.3.3


def _jm(J: int):
    d = RT.jmat(J)
    return [np.array(d[k], dtype=complex) for k in ("Jx", "Jy", "Jz")]


def rotor_H(J: int | None = None) -> np.ndarray:
    J = ROTOR["J"] if J is None else J
    Jx, Jy, Jz = _jm(J)
    return ROTOR["A"] * Jx @ Jx + ROTOR["B"] * Jy @ Jy + ROTOR["C"] * Jz @ Jz


def rotor_packet0() -> np.ndarray:
    """(8.3.24)：R(0, β0, 0)|J J⟩（RT.wigner_D 的第一行，β0 = −30°）。"""
    return np.array(RT.wigner_D(ROTOR["J"], 0, ROTOR["beta0"], 0))[:, 0]


def rotor_exact() -> dict:
    """K = 8 叢的兩個量子能量、初始波包 ⟨H⟩、J² 約定的古典能量（書 E_classical = 50.000）、精確 AC 的 FTAC 峰（叢的加權中心）。"""
    H = rotor_H()
    w, v = np.linalg.eigh(H)
    lv = SC.rigid_levels(ROTOR["J"], ROTOR["A"], ROTOR["B"], ROTOR["C"])
    psi = rotor_packet0()
    wt = abs(v.conj().T @ psi) ** 2
    pair = sorted(float(e) for e in w if abs(e - 53.1475) < 0.01)
    b = math.radians(ROTOR["beta0"])
    n = np.array([math.sin(b), 0.0, math.cos(b)])
    Ecl = ROTOR["J"] ** 2 * (ROTOR["A"] * n[0] ** 2 + ROTOR["B"] * n[1] ** 2 + ROTOR["C"] * n[2] ** 2)
    T = 60.0
    ts = np.linspace(0, T, 12001)
    ac = (wt[None, :] * np.exp(-1j * np.outer(ts, w))).sum(axis=1)
    fun = lambda e: float(ftac(ts, ac, [e])[0])
    eft = _refine_peak(fun, 53.15, 0.3)[0]
    cent = float(np.sum([wt[i] * w[i] for i in range(len(w)) if abs(w[i] - 53.1475) < 0.01]) /
                 np.sum([wt[i] for i in range(len(w)) if abs(w[i] - 53.1475) < 0.01]))
    return dict(pair=pair, pair_sc=sorted(e for e, _ in lv if abs(e - 53.1475) < 0.01), E_packet=float(np.real(psi.conj() @ H @ psi)),
                E_classical=Ecl, levels=[float(e) for e in w], weights=[float(c) for c in wt], E_FT_exact=eft, centroid=cent)


def rotor_localized(K: int = 8) -> list:
    """(8.3.32)：K = 8 叢兩個本徵態中 ⟨J_z⟩ 最大的組合（局域在 +z 的 C 軸）在 K' = 10, 8, …, 0 的分量。"""
    H = rotor_H()
    w, v = np.linalg.eigh(H)
    idx = [i for i in range(len(w)) if abs(w[i] - 53.1475) < 0.01]
    V = v[:, idx]
    Jz = _jm(ROTOR["J"])[2]
    mw, mv = np.linalg.eigh(V.conj().T @ Jz @ V)
    loc = V @ mv[:, -1]
    ms = list(range(ROTOR["J"], -ROTOR["J"] - 1, -1))
    loc = loc * np.exp(-1j * np.angle(loc[ms.index(K)]))
    return [float(np.real(loc[ms.index(k)])) for k in (10, 8, 6, 4, 2, 0)]


def _cs_state(theta: float, phi: float, U, mval, Jz_diag) -> np.ndarray:
    """|n⟩ = e^{−iφJz} e^{−iθJy}|J J⟩（Jy = U diag(m) U†）。"""
    J = int(round(max(mval)))
    e0 = np.zeros(2 * J + 1, complex)
    e0[0] = 1
    v = U @ (np.exp(-1j * theta * mval) * (U.conj().T @ e0))
    return np.exp(-1j * phi * Jz_diag) * v


def rotor_sc(T: float = 40.0, nt: int = 4001) -> dict:
    """自旋相干態的時間相依變分（TDVP）半古典傳播（(8.3.25)–(8.3.30) 的做法）：
    ṅ = (1/J)∇E × n，E(n) = ⟨n|H|n⟩ = J(J−½)Σ A_i n_i² + (J/2)ΣA_i；相位 S = ∫(J cos θ φ̇ − E)dt；|Ψ(t)⟩ = e^{iS}|n(t)⟩。
    → FTAC 的 K = 8 峰 E_FT、FT 態的 E_SC = ⟨H⟩（(8.3.17)）與係數 ⟨K'|E_FT(8)⟩（(8.3.31)）。"""
    J = ROTOR["J"]
    Av = np.array([ROTOR["A"], ROTOR["B"], ROTOR["C"]])
    H = rotor_H()
    Jx, Jy, Jz = _jm(J)
    mval, U = np.linalg.eigh(Jy)
    Jz_diag = np.real(np.diag(Jz))
    c1, c0 = J * (J - 0.5), J / 2 * Av.sum()

    def rhs(t, y):
        n = y[:3]
        grad = 2 * c1 * Av * n
        nd = np.cross(grad, n) / J
        rho2 = n[0] ** 2 + n[1] ** 2
        phid = (n[0] * nd[1] - n[1] * nd[0]) / rho2
        E = c1 * np.sum(Av * n * n) + c0
        return np.concatenate([nd, [J * n[2] * phid - E, phid]])
    b = math.radians(ROTOR["beta0"])
    n0 = np.array([math.sin(b), 0.0, math.cos(b)])
    phi0 = math.atan2(n0[1], n0[0])
    y0 = np.concatenate([n0, [0.0, phi0]])
    ts = np.linspace(0, T, nt)
    fw = solve_ivp(rhs, (0, T), y0, t_eval=ts, rtol=1e-11, atol=1e-12, method="DOP853")
    bw = solve_ivp(rhs, (0, -T), y0, t_eval=-ts, rtol=1e-11, atol=1e-12, method="DOP853")

    def states(sol):
        out = []
        for k in range(sol.y.shape[1]):
            n = sol.y[:3, k]
            th = math.acos(max(-1.0, min(1.0, n[2] / np.linalg.norm(n))))
            out.append(np.exp(1j * sol.y[3, k]) * _cs_state(th, sol.y[4, k], U, mval, Jz_diag))
        return np.array(out)
    psf, psb = states(fw), states(bw)
    psi0 = psf[0]
    ac = psf @ psi0.conj()
    E_cs = float(np.real(psi0.conj() @ H @ psi0))
    fun = lambda e: float(ftac(ts, ac, [e])[0])
    eft = _refine_peak(fun, 53.1, 0.6)[0]
    w = np.full(nt, ts[1] - ts[0])
    w[0] = w[-1] = w[0] / 2
    vec = (w * np.exp(1j * eft * ts)) @ psf + (w * np.exp(-1j * eft * ts)) @ psb - w[0] * psi0
    esc = float(np.real(vec.conj() @ H @ vec) / np.real(vec.conj() @ vec))
    vec = vec / np.linalg.norm(vec)
    ms = list(range(J, -J - 1, -1))
    vec = vec * np.exp(-1j * np.angle(vec[ms.index(8)]))
    coeffs = [float(np.real(vec[ms.index(k)])) for k in (10, 8, 6, 4, 2, 0)]
    nrm = float(np.max(abs(np.linalg.norm(fw.y[:3], axis=0) - 1)))
    Ecl = c1 * np.sum(Av[:, None] * fw.y[:3] ** 2, axis=0) + c0
    return dict(E_cs=E_cs, E_FT=eft, E_SC=esc, coeffs=coeffs, norm_drift=nrm, E_drift=float(np.max(abs(Ecl - E_cs))),
                ts=ts, ac=ac, n_path=fw.y[:3].T)


def rotor_ebk(K: int = 8, Jeff: float | None = None) -> float:
    """EBK（轉動型軌道沒有轉折點，α = 0）：(1/2π)∮ J_z dφ = K，|J| = J + ½；E(n) = |J|²(A n_x² + B n_y² + C n_z²)。"""
    Jm = ROTOR["J"] + 0.5 if Jeff is None else Jeff
    A, B, C = ROTOR["A"], ROTOR["B"], ROTOR["C"]

    def action(E):
        def jz(ph):
            a = A * math.cos(ph) ** 2 + B * math.sin(ph) ** 2
            c2 = (E / Jm ** 2 - a) / (C - a)
            return Jm * math.sqrt(max(c2, 0.0))
        with warnings.catch_warnings():                 # 被積函數在 c2 → 0 有根號尖點：精度仍 ~1e−9
            warnings.simplefilter("ignore")
            return quad(jz, 0, 2 * PI, limit=400)[0] / (2 * PI)
    return brentq(lambda E: action(E) - K, Jm ** 2 * (B + 1e-9) + 1e-6, Jm ** 2 * C - 1e-9, xtol=1e-12)


# ======================================================================== §8.4A A·p 與 E·r
def gauge_check() -> dict:
    """一維偶極近似：ψ^α = e^{iqA(t)x/ħ}ψ^β；若 iħψ^β_t = −ħ²ψ^β_xx/2μ + Vψ^β − qE xψ^β（E = −dA/dt，(8.4.15)），
    則 iħψ^α_t = (p − qA)²ψ^α/2μ + Vψ^α（(8.4.9)）。書上 (8.4.14) 少 q（e^{iAx/ħ}）就不成立。"""
    x, t = sp.symbols("x t", real=True)
    hb, mu, q = sp.symbols("hbar mu q", positive=True)
    A = sp.Function("A")(t)
    V = sp.Function("V")(x)
    pb = sp.Function("psi")(x, t)
    Ef = -sp.diff(A, t)
    pbt = (-hb ** 2 / (2 * mu) * sp.diff(pb, x, 2) + V * pb - q * Ef * x * pb) / (sp.I * hb)

    def residual(qq):
        pa = sp.exp(sp.I * qq * A * x / hb) * pb
        lhs = sp.I * hb * (sp.diff(pa, t).subs(sp.Derivative(pb, t), pbt))
        P = lambda f: -sp.I * hb * sp.diff(f, x) - q * A * f
        rhs = P(P(pa)) / (2 * mu) + V * pa
        return sp.simplify(sp.expand((lhs - rhs) * sp.exp(-sp.I * qq * A * x / hb)))
    return dict(ok=residual(q) == 0, book_fails=residual(1) != 0)


def hamiltonian_check() -> bool:
    """(8.4.4)–(8.4.5)：L = μv²/2 − qΦ + q v·A ⇒ p = μv + qA、H = v·p − L = (p − qA)²/2μ + qΦ ＝ (8.4.5c) 的展開。"""
    v, A_, Phi, mu, q, p = sp.symbols("v A Phi mu q p", real=True)
    L = mu * v ** 2 / 2 - q * Phi + q * v * A_
    pv = sp.diff(L, v)
    H = sp.simplify((v * pv - L).subs(v, sp.solve(sp.Eq(p, pv), v)[0]))
    ok1 = sp.simplify(H - ((p - q * A_) ** 2 / (2 * mu) + q * Phi)) == 0
    ok2 = sp.simplify(H - (p ** 2 / (2 * mu) - q / (2 * mu) * (2 * p * A_) + q ** 2 / (2 * mu) * A_ ** 2 + q * Phi)) == 0
    return bool(ok1 and ok2)


QE = 1.602176634e-19
ME = 9.1093837015e-31
HBAR = 1.054571817e-34
EPS0 = 8.8541878128e-12
CLIGHT = 299792458.0
A0 = 5.29177210903e-11
EH = 4.3597447222071e-18
TAU_AU = HBAR / EH                  # 2.4189e−17 s


def boost_shift(E0: float = 1e6, omega: float = 1e16, q: float = 1.6e-19, mu: float = 9.1e-31) -> float:
    """(8.4.22)：E 座標相對 A 座標的古典平移振幅 qE0/μω²（書上的數：1.758e−15 m，書寫約 2e−15；q 印成 1.6e19）。"""
    return q * E0 / (mu * omega ** 2)


def ap_er_ratio(w_nb: float, w: float) -> float:
    """(8.4.51) 對 (8.4.35a)：天真的 A·p 矩陣元與 E·r 矩陣元的大小比 = ω_NB/ω（只在共振時相等）。"""
    return w_nb / w


# ======================================================================== §8.4B 譜函數、微擾、振子強度
def S_amp(Delta, tau):
    """(8.4.38a)：∫_0^τ e^{iΔt}dt = (2/Δ)e^{iτΔ/2} sin(τΔ/2)（Δ = 0 時 τ）。"""
    Delta = np.asarray(Delta, float)
    with np.errstate(divide="ignore", invalid="ignore"):
        out = np.where(abs(Delta) < 1e-14, tau + 0j, 2 / np.where(Delta == 0, 1, Delta) * np.exp(1j * tau * Delta / 2) * np.sin(tau * Delta / 2))
    return out


def S2(Delta, tau):
    return abs(S_amp(Delta, tau)) ** 2


def spectral_card(tau: float = 4.0) -> dict:
    """|S|² 的峰 τ²（(8.4.41)）、零點 ±2π/τ（(8.4.42)）、面積 2πτ（(8.4.40)）、第一旁帶（u = Δτ/2 解 tan u = u）。"""
    area = 2 * quad(lambda d: float(S2(d, tau)), 0, 200, limit=2000)[0] + 2 * 2 / 200   # 尾巴：兩側各 ∫_{200}^∞ 2/Δ² dΔ（sin² 平均 ½）
    u1 = brentq(lambda u: math.tan(u) - u, PI + 0.1, 1.5 * PI - 1e-6)
    dside = 2 * u1 / tau
    zero = brentq(lambda d: float(np.real(S_amp(d, tau) * np.exp(-1j * tau * d / 2))), 1.0, 2.0)
    num = quad(lambda s: math.cos(0.7 * s), 0, tau)[0], quad(lambda s: math.sin(0.7 * s), 0, tau)[0]
    return dict(peak=float(S2(0.0, tau)), zero=zero, area=area, area_formula=2 * PI * tau, side_D=dside,
                side_h=float(S2(dside, tau)) / tau ** 2, book_side_D=3 * PI / tau, book_side_h=float(S2(3 * PI / tau, tau)) / tau ** 2,
                closed_vs_quad=abs(complex(*num) - complex(S_amp(0.7, tau))) < 1e-12)


def first_order_check() -> dict:
    """弱驅動兩能階（ω_NB = 1、ω = 0.9、E·r 的實餘弦）：數值 |c_N(t)|² ＝ (8.4.37) 兩項（c ≈ (q r E0/2ħ)·i·[e^{iφ}S(Δ⁻) + e^{−iφ}S(Δ⁺)]）；
    (e^{iΔt} − 1)/Δ = iS(Δ, t)（書上 (8.4.37) 少 i）。"""
    wnb, w, V0, phi, tmax = 1.0, 0.9, 0.002, 0.3, 30.0

    def rhs(t, y):
        c = y[:2] + 1j * y[2:]
        Hi = -V0 * math.cos(w * t - phi)                     # −qE0 cos(ωt − φ) r_NB，r_NB = 1
        dc = -1j * np.array([Hi * np.exp(-1j * wnb * t) * c[1], Hi * np.exp(1j * wnb * t) * c[0]])
        return np.concatenate([dc.real, dc.imag])
    sol = solve_ivp(rhs, (0, tmax), [1, 0, 0, 0], rtol=1e-11, atol=1e-13)
    cN = sol.y[1, -1] + 1j * sol.y[3, -1]
    pert = 1j * V0 / 2 * (np.exp(1j * phi) * S_amp(wnb - w, tmax) + np.exp(-1j * phi) * S_amp(wnb + w, tmax))
    d = 0.37
    return dict(ok=abs(cN - pert) / abs(pert) < 5e-3, iS=abs((np.exp(1j * d * 2.5) - 1) / d - 1j * S_amp(d, 2.5)) < 1e-12)


def golden_rule_check(V: float = 1e-3, dD: float = 0.01, t: float = 30.0) -> float:
    """準連續帶（Δ_k = k·dΔ）：P(t) = Σ|V|²|S(Δ_k, t)|² ≈ 2π|V|²t/dΔ（ρ = 1/dΔ）→ 回傳比值。"""
    ks = np.arange(-200000, 200001)
    P = float(np.sum(V ** 2 * S2(ks * dD, t)))
    return P / (2 * PI * V ** 2 * t / dD)


def ho_strengths(N: int = 30) -> list:
    """(8.4.47b)：f_{n0} = 2ω_{n0}|x_{n0}|²（μ = ħ = ω = 1）＝ δ_{n1}（(8.4.48)）。"""
    x, _ = xp_ops(N)
    return [2 * n * abs(x[n, 0]) ** 2 for n in range(1, 8)]


def h_fn(n: int) -> float:
    """氫 1s → np 的振子強度 f_n = 2⁸n⁵(n−1)^{2n−4}/3(n+1)^{2n+4}。"""
    return math.exp(8 * math.log(2) + 5 * math.log(n) + (2 * n - 4) * math.log(n - 1) - math.log(3) - (2 * n + 4) * math.log(n + 1))


def h_strengths(nmax: int = 20000) -> dict:
    bound = sum(h_fn(n) for n in range(2, nmax + 1)) + (2 ** 8 / 3) * math.exp(-4) / (2 * nmax ** 2)
    g = lambda k: (2 ** 8 / 3) * math.exp(-4 / k * math.atan(k)) / ((1 + k * k) ** 4 * (1 - math.exp(-2 * PI / k))) * k
    cont = quad(g, 1e-9, 50, limit=200)[0] + quad(g, 50, np.inf)[0]
    return dict(f2=h_fn(2), f3=h_fn(3), f4=h_fn(4), bound=bound, cont=cont, total=bound + cont)


def h_box_strengths(Rmax: float = 80.0, n: int = 8000) -> dict:
    """徑向方程 −½u'' + [l(l+1)/2r² − 1/r]u = Eu 的有限差分（盒 Rmax）：1s → 所有 l = 1 態，f = (2/3)ΔE|∫u_1s r u_p dr|²。"""
    h = Rmax / (n + 1)
    r = h * np.arange(1, n + 1)
    off = np.full(n - 1, -0.5 / h ** 2)

    def solve(l):
        d = 1 / h ** 2 + l * (l + 1) / (2 * r ** 2) - 1 / r
        return eigh_tridiagonal(d, off)
    e0, u0 = solve(0)
    e1, u1 = solve(1)
    g = u0[:, 0]
    M = u1.T @ (r * g)
    f = 2 / 3 * (e1 - e0[0]) * M ** 2
    return dict(E1s=float(e0[0]), f2=float(f[0]), f3=float(f[1]), total=float(f.sum()))


def classical_factor_check() -> dict:
    """§8.4B(b) 的古典振子：|e^{−ia} − e^{−ib}|² = 4 sin²((a−b)/2) ⇒ |x|² 的係數是 4|(q/μ)E0|²（書印 |4(q/μ)E0|²，大 4 倍）。"""
    a, b = 1.234, 0.456
    lhs = abs(np.exp(-1j * a) - np.exp(-1j * b)) ** 2
    return dict(ok=abs(lhs - 4 * math.sin((a - b) / 2) ** 2) < 1e-14, book_over_right=16 / 4)


# ======================================================================== §8.5 兩能階
TL_CARD = dict(eps=1.0, r=0.2)              # 本系列（Rabi、Bloch–Siegert）
FIG853 = dict(eps=1.0, r=0.5)               # 書上 Fig. 8.5.3 的 ε = 1.0、r = 0.5（Fig. 8.5.4 沒給 r，本系列也取 0.5）
FIG854_DELTAS = (1.2, 0.2, -0.1, -1.2)


def crank(t, r: float, Om: float, eps: float) -> np.ndarray:
    """(8.5.6)：ω = (r cos Ωt, r sin Ωt, ε)。"""
    return np.array([r * math.cos(Om * t), r * math.sin(Om * t), eps])


def h_lab(t: float, r: float, Om: float, eps: float, rwa: bool = True) -> np.ndarray:
    """(8.5.5)（rwa = True）或 (8.5.4a)（實的餘弦曲柄）。"""
    if rwa:
        return np.array([[eps / 2, r / 2 * np.exp(-1j * Om * t)], [r / 2 * np.exp(1j * Om * t), -eps / 2]])
    c = r * math.cos(Om * t)
    return np.array([[eps / 2, c], [c, -eps / 2]], complex)


def evolve_lab(ts, r: float, Om: float, eps: float, psi0=(1, 0), rwa: bool = True) -> np.ndarray:
    def rhs(t, y):
        psi = y[:2] + 1j * y[2:]
        d = -1j * h_lab(t, r, Om, eps, rwa) @ psi
        return np.concatenate([d.real, d.imag])
    p0 = np.asarray(psi0, complex)
    sol = solve_ivp(rhs, (0, ts[-1]), np.concatenate([p0.real, p0.imag]), t_eval=ts, rtol=1e-11, atol=1e-12, method="DOP853")
    return (sol.y[:2] + 1j * sol.y[2:]).T


def spin_precession(ts, r: float, Om: float, eps: float, S0=(0, 0, 0.5)) -> np.ndarray:
    """dS/dt = ω(t) × S（數值）。"""
    sol = solve_ivp(lambda t, S: np.cross(crank(t, r, Om, eps), S), (0, ts[-1]), np.asarray(S0, float), t_eval=ts,
                    rtol=1e-11, atol=1e-12, method="DOP853")
    return sol.y.T


def rotating_check(r: float = 0.2, Om: float = 0.9, eps: float = 1.0, t: float = 7.3) -> dict:
    """(8.5.11)：轉動座標的 H 常數；(8.5.13a) 的曲柄 (r, 0, Δ)；(8.5.15a) ＝ 數值；SC.h_to_omega 給 (8.5.6)。"""
    Rm = np.diag([np.exp(1j * Om * t / 2), np.exp(-1j * Om * t / 2)])
    HR = Rm @ h_lab(t, r, Om, eps) @ Rm.conj().T
    ok1 = _close(HR, 0.5 * np.array([[eps, r], [r, -eps]]))
    Delta = eps - Om
    HRp = HR - Om * np.diag([0.5, -0.5])
    wv = SC.h_to_omega(HRp)["omega_vec"]
    ok2 = _close(wv, [r, 0, Delta])
    om = math.hypot(Delta, r)
    psiR = np.array([math.cos(om * t / 2) - 1j * Delta / om * math.sin(om * t / 2), -1j * r / om * math.sin(om * t / 2)])
    lab = evolve_lab(np.array([0, t]), r, Om, eps)[-1]
    ok3 = _close(Rm @ lab, psiR, 1e-8)
    ok4 = _close(SC.h_to_omega(h_lab(t, r, Om, eps))["omega_vec"], crank(t, r, Om, eps))
    return dict(ok=ok1 and ok2 and ok3 and ok4)


def rabi(Delta, r: float, t):
    """(8.5.17)：|c2|² = r²/(Δ² + r²) sin²(√(Δ² + r²) t/2)。"""
    om = np.hypot(Delta, r)
    return r ** 2 / om ** 2 * np.sin(om * np.asarray(t, float) / 2) ** 2


def rabi_numeric(Delta: float, r: float, ts, eps: float = 1.0) -> np.ndarray:
    psi = evolve_lab(ts, r, eps - Delta, eps)
    return abs(psi[:, 1]) ** 2


def rabi_card() -> dict:
    r = TL_CARD["r"]
    t_pi = PI / r
    t_half = PI / math.hypot(r, r)
    return dict(t_pi=t_pi, p_res=float(rabi(0.0, r, t_pi)), t_half=t_half, p_half=float(rabi(r, r, t_half)), fwhm=2 * r,
                env_half=r ** 2 / (r ** 2 + r ** 2))


def _u_step(eps: float, c: float, dt: float) -> np.ndarray:
    """H = (ε/2)σ_z + cσ_x 的 e^{−iH dt}（2×2 閉式）。"""
    hz, hx = eps / 2, c
    h = math.hypot(hz, hx)
    cs, sn = math.cos(h * dt), math.sin(h * dt) / h
    return np.array([[cs - 1j * sn * hz, -1j * sn * hx], [-1j * sn * hx, cs + 1j * sn * hz]])


def floquet_population(Om: float, eps: float, r: float, nst: int = 400) -> float:
    """全 Hamilton 量 (8.5.4a) 的一週期傳播子（中點、每步閉式指數）→ Floquet 模在一週期內 |1⟩ 布居的時間平均（兩模的較大者）。"""
    Tp = 2 * PI / Om
    dt = Tp / nst
    Us = [np.eye(2, dtype=complex)]
    for k in range(nst):
        Us.append(_u_step(eps, r * math.cos(Om * (k + 0.5) * dt), dt) @ Us[-1])
    _, v = np.linalg.eig(Us[-1])
    U = np.array(Us[:-1])
    pops = [float(np.mean(abs((U @ v[:, j])[:, 0]) ** 2)) for j in range(2)]
    return max(pops)


def bloch_siegert(eps: float = 1.0, r: float = 0.2) -> dict:
    """Floquet 共振（兩個 Floquet 模的時間平均布居都是 ½：|1⟩、|2⟩ 完全混合）→ Ω_res − ε（數值）；
    對照 r²/4ε（(8.5.26a)）與書上 (8.5.25) 的啟發式 (ε + √(4ε² + 3r²))/3 − ε。"""
    g = lambda Om: floquet_population(Om, eps, r) - 0.5
    Os = np.linspace(eps - 0.5 * r ** 2 / eps, eps + 1.5 * r ** 2 / eps, 41)
    vals = [g(O) for O in Os]
    i = int(np.argmin(vals))
    step = Os[1] - Os[0]
    res = minimize_scalar(g, bounds=(Os[i] - step, Os[i] + step), method="bounded", options=dict(xatol=1e-11))
    return dict(shift=float(res.x - eps), formula=r ** 2 / (4 * eps), book_8525=(eps + math.sqrt(4 * eps ** 2 + 3 * r ** 2)) / 3 - eps,
                next_order=r ** 2 / (4 * eps) + r ** 4 / (64 * eps ** 3), pop_min=float(res.fun) + 0.5,
                fraction=r ** 2 / (4 * eps ** 2), book_example=(1e10) ** 2 / (4 * (1e16) ** 2))


def bs_sign_check(r: float = 0.2, Om: float = 0.95, eps: float = 1.0, t: float = 0.4) -> dict:
    """(8.5.20)：轉動座標中全 H 的 H21 = ½(re^{−2iΩt} + r) ⇒ ω_y = 2Im H21 = −r sin 2Ωt（書 (8.5.21) 印 +）。"""
    Rm = np.diag([np.exp(1j * Om * t / 2), np.exp(-1j * Om * t / 2)])
    HR = Rm @ h_lab(t, r, Om, eps, rwa=False) @ Rm.conj().T - Om * np.diag([0.5, -0.5])
    wv = SC.h_to_omega(HR)["omega_vec"]
    return dict(ok=_close(wv, [r + r * math.cos(2 * Om * t), -r * math.sin(2 * Om * t), eps - Om]))


BLOCH_CARD = dict(Delta=0.0, r=0.2, g1=0.05, g2=0.1, Sz0=-0.5)     # 本系列


def bloch_rhs(t, S, Delta, r, g1, g2, Sz0):
    """(8.5.35a)（轉動座標、ω = (r, 0, Δ)）：dS/dt = ω × S − [γ2 S_x, γ2 S_y, γ1(S_z − S_z0)]。"""
    w = np.array([r, 0.0, Delta])
    return np.cross(w, S) - np.array([g2 * S[0], g2 * S[1], g1 * (S[2] - Sz0)])


def bloch_steady(Delta, r, g1, g2, Sz0) -> np.ndarray:
    """解析穩態：S_z = S_z0/(1 + r²γ2/γ1(γ2² + Δ²))、S_y = −rγ2S_z/(γ2² + Δ²)、S_x = −ΔS_y/γ2。"""
    Sz = Sz0 / (1 + r ** 2 * g2 / (g1 * (g2 ** 2 + Delta ** 2)))
    Sy = -r * g2 * Sz / (g2 ** 2 + Delta ** 2)
    return np.array([-Delta * Sy / g2, Sy, Sz])


def bloch_numeric(Delta, r, g1, g2, Sz0, T: float = 600.0) -> np.ndarray:
    sol = solve_ivp(bloch_rhs, (0, T), [0, 0, Sz0], args=(Delta, r, g1, g2, Sz0), rtol=1e-10, atol=1e-12, method="DOP853")
    return sol.y[:, -1]


def thermal_sz(beta_eps: float) -> dict:
    """ρ = e^{−βH}/Z、H = (ε/2)σ_z（|1⟩ 在上）：⟨S_z⟩ = ½Tr(ρσ_z) = −½ tanh(βε/2)（書 (8.5.36) 後的式子把下能階叫 1：+½ tanh）。"""
    Hm = np.diag([0.5, -0.5]) * beta_eps
    rho = expm(-Hm)
    rho /= np.trace(rho)
    num = float(np.real(np.trace(rho @ np.diag([1, -1]))) / 2)
    return dict(numeric=num, formula=-0.5 * math.tanh(beta_eps / 2))


def dressed(Delta: float, r: float) -> dict:
    """(8.5.39)–(8.5.44)：θ = atan2(r, Δ)、R_D 對角化 H^{R'} 成 diag(ω/2, −ω/2)；|1⟩^D = (cos θ/2, sin θ/2) 的本徵值是 +ω/2。"""
    th = math.atan2(r, Delta)
    Rd = np.array([[math.cos(th / 2), -math.sin(th / 2)], [math.sin(th / 2), math.cos(th / 2)]])
    Hp = np.array([[Delta / 2, r / 2], [r / 2, -Delta / 2]])
    om = math.hypot(Delta, r)
    return dict(theta=math.degrees(th), ok=_close(Rd.T @ Hp @ Rd, np.diag([om / 2, -om / 2])), omega=om)


def ac_stark(Om, eps: float, r: float):
    """(8.5.45)：上 ε + δ/2、下 Ω − δ/2，δ = √(Δ² + r²) − Δ、Δ = ε − Ω（＝ [[ε, r/2], [r/2, Ω]] 的本徵值）。"""
    Om = np.asarray(Om, float)
    Dl = eps - Om
    dl = np.hypot(Dl, r) - Dl
    return eps + dl / 2, Om - dl / 2


def ac_stark_card() -> dict:
    eps, r = FIG853["eps"], FIG853["r"]
    up, lo = ac_stark(eps, eps, r)
    w = np.linalg.eigvalsh(np.array([[eps, r / 2], [r / 2, eps]]))
    up0, lo0 = ac_stark(0.0, eps, r)
    wdc = np.linalg.eigvalsh(np.array([[eps, r / 2], [r / 2, 0.0]]))
    return dict(split=float(up - lo), split_num=float(w[1] - w[0]), ac_shift=float(up - eps), dc_shift=float(up0 - eps),
                dc_formula=(math.hypot(eps, r) - eps) / 2, dc_num=float(wdc[1] - eps))


def circle_geometry(Delta: float, r: float) -> dict:
    """Fig. 8.5.4：δ = ω − Δ、δ' = ω + Δ（δδ' = r²）、δ² + r² = 2ωδ（(8.5.46b)）、tan(θ/2) = δ/r、|1⟩^D 振幅 r/√(2δω)、δ/√(2δω)。"""
    om = math.hypot(Delta, r)
    dl, dlp = om - Delta, om + Delta
    th = math.atan2(r, Delta)
    hyp = math.sqrt(2 * dl * om)
    return dict(omega=om, delta=dl, deltap=dlp, theta=math.degrees(th), hyp=hyp, amp1=r / hyp, amp2=dl / hyp,
                ok=abs(dl * dlp - r ** 2) < 1e-12 and abs(dl ** 2 + r ** 2 - 2 * om * dl) < 1e-12
                and abs(math.tan(th / 2) - dl / r) < 1e-12 and abs(math.cos(th / 2) - r / hyp) < 1e-12)


ADIAB = dict(r=0.5, D0=20.0, T=400.0)         # 本系列：Δ 由 +20 線性掃到 −20（v = 0.1）


def adiabatic(r: float = 0.5, D0: float = 20.0, T: float = 400.0) -> dict:
    """轉動座標 H = Δ(t)σ_z/2 + rσ_x/2，Δ 從 +D0 線性掃到 −D0：數值反轉 ＝ 1 − e^{−πr²/2v}（Landau–Zener，v = 2D0/T）。"""
    v = 2 * D0 / T
    th0 = math.atan2(r, D0)
    psi0 = np.array([-math.sin(th0 / 2), math.cos(th0 / 2)], complex)     # 下穿衣態（主要是 |2⟩，S_z = −½ 一側）

    def rhs(t, y):
        Dl = D0 - v * t
        psi = y[:2] + 1j * y[2:]
        Hm = np.array([[Dl / 2, r / 2], [r / 2, -Dl / 2]])
        d = -1j * Hm @ psi
        return np.concatenate([d.real, d.imag])
    sol = solve_ivp(rhs, (0, T), np.concatenate([psi0.real, psi0.imag]), rtol=1e-10, atol=1e-12, method="DOP853")
    psi = sol.y[:2, -1] + 1j * sol.y[2:, -1]
    return dict(inversion=float(abs(psi[0]) ** 2), lz=1 - math.exp(-PI * r ** 2 / (2 * v)))


# ======================================================================== §8.6 量子化的場
def photon_elements(n: int, N: int = 40) -> dict:
    a = ladder(N)
    return dict(absorb=float(abs(a[n - 1, n])), emit=float(abs(a.conj().T[n + 1, n])))


def mode_energy_check() -> dict:
    """(8.6.3)：⟨ε0E²/2 + B²/2μ0⟩ = 2ε0ω²|a|²（E0 = 2|a|ω、B0 = 2|a|k、k² = ω²ε0μ0）；照印的磁場項 4|a|²/2μ0 少 k²（因次不對）。"""
    a, w, e0, m0 = sp.symbols("a omega epsilon_0 mu_0", positive=True)
    k = w * sp.sqrt(e0 * m0)
    U = (e0 / 2 * (2 * a * w) ** 2 + (2 * a * k) ** 2 / (2 * m0)) / 2
    one = sp.solve(sp.Eq(2 * e0 * w ** 2 * a ** 2, sp.Symbol("hbar", positive=True) * w), a)
    return dict(ok=sp.simplify(U - 2 * e0 * w ** 2 * a ** 2) == 0, one_photon=str(one[0]))


def angular_sum() -> float:
    """∫dΩ_k Σ_α |e_α·r̂|² = ∫dΩ (1 − cos²θ) = 8π/3（Σ_α|e_α·r|² = |r|² − |k̂·r|²）。"""
    return 2 * PI * quad(lambda th: math.sin(th) ** 3, 0, PI)[0]


def angular_sum_mc(n: int = 200000, seed: int = 3) -> float:
    """同上，用隨機 k̂ 與兩個正交偏振向量直接相加（獨立方法）。"""
    rng = np.random.default_rng(seed)
    k = rng.normal(size=(n, 3))
    k /= np.linalg.norm(k, axis=1)[:, None]
    tmp = np.where(abs(k[:, :1]) < 0.9, np.array([[1.0, 0, 0]]), np.array([[0, 1.0, 0]]))
    e1 = np.cross(k, tmp)
    e1 /= np.linalg.norm(e1, axis=1)[:, None]
    e2 = np.cross(k, e1)
    rz = np.array([0, 0, 1.0])
    return float(4 * PI * np.mean((e1 @ rz) ** 2 + (e2 @ rz) ** 2))


def einstein_A(omega: float, r2: float, q: float = QE) -> float:
    """(8.6.32b)（SI）：A = 4ω³q²|r|²/(3ħc³·4πε0)。"""
    return 4 * omega ** 3 / (3 * HBAR * CLIGHT ** 3) * q ** 2 / (4 * PI * EPS0) * r2


def einstein_A_gauss(omega: float, r2_cm2: float) -> float:
    """Gauss 單位：A = 4ω³e²|r|²/3ħc³（e = 4.80320e−10 esu、ħ 與 c 用 cgs）。"""
    e = QE * CLIGHT * 10                          # esu
    hb, c = HBAR * 1e7, CLIGHT * 100
    return 4 * omega ** 3 * e ** 2 * r2_cm2 / (3 * hb * c ** 3)


def h2p_A() -> dict:
    """氫 2p → 1s：|⟨1s|z|2p0⟩|² = 2¹⁵/3¹⁰ a0²（閉式 ＝ 徑向數值積分 × 角因子 1/3）、ω = (3/8)E_h/ħ。"""
    r2_exact = 2 ** 15 / 3 ** 10
    rad = quad(lambda r: 2 * math.exp(-r) * r * math.exp(-r / 2) / (2 * math.sqrt(6)) * r ** 3, 0, np.inf)[0]
    r2_num = rad ** 2 / 3
    om = 0.375 * EH / HBAR
    A_si = einstein_A(om, r2_exact * A0 ** 2)
    A_cgs = einstein_A_gauss(om, r2_exact * (A0 * 100) ** 2)
    alpha = 7.2973525693e-3
    A_au = 4 / 3 * 0.375 ** 3 * r2_exact * alpha ** 3 / TAU_AU
    mred = 1 / (1 + ME / 1.67262192369e-27)              # 約化質量：ω ∝ μ、|r|² ∝ 1/μ² ⇒ A ∝ μ
    return dict(r2=r2_exact, r2_num=r2_num, A=A_si, A_cgs=A_cgs, A_au=A_au, A_red=A_si * mred, tau_ns=1e9 / A_si,
                book_printed=A_si * 2 * PI, omega=om)


def golden_A_numeric(t_au: float = 4000.0) -> float:
    """(8.6.28)：P(t) = Σ_α∫dΩ∫dω ρ(ω)|⟨H_I⟩|²|S|²/ħ²（原子單位，偶極沿 z）→ P/t ＝ A（t 大時）。"""
    alpha = 7.2973525693e-3
    c = 1 / alpha
    w0 = 0.375
    r2 = 2 ** 15 / 3 ** 10
    # ρ(ω)/V·(ħω/2ε0)·q²·(8π/3)|r|²·n=1：原子單位 ε0 = 1/4π、ħ = q = 1
    pref = lambda w: w ** 2 / ((2 * PI) ** 3 * c ** 3) * (w / (2 / (4 * PI))) * (8 * PI / 3) * r2
    f = lambda w: pref(w) * float(S2(w - w0, t_au))
    span = 400 / t_au
    P = quad(f, w0 - span, w0 + span, points=[w0], limit=4000)[0] + quad(f, 1e-6, w0 - span, limit=400)[0] + quad(f, w0 + span, 3 * w0, limit=400)[0]
    return P / t_au / TAU_AU


def number_vs_coherent(n: int = 5, N: int = 60) -> dict:
    """§8.6B(c)–(d)：數態 ⟨n|a|n⟩ = 0（「無能」）、⟨n|x²|n⟩ = n + ½；相干態 ⟨α|a|α⟩ = α（(8.6.35)）。"""
    a = ladder(N)
    v = np.zeros(N, complex)
    v[n] = 1
    x = (a + a.conj().T) / math.sqrt(2)
    al = math.sqrt(5)
    c = coherent_cn(al, N - 1)
    return dict(num_a=float(abs(v.conj() @ a @ v)), num_x2=float(np.real(v.conj() @ x @ x @ v)), coh_a=complex(c.conj() @ a @ c))


# ======================================================================== §8.7 Jaynes–Cummings
JC_CARD = dict(Omega=10.0, Delta=0.3, r=0.4)


def jc_matrix(Nmax: int, Om: float, eps: float, g: float, book_phase: bool = False) -> np.ndarray:
    """H'_JCM（(8.7.10) 的實版本；book_phase = True 用 (8.7.5) 的 i(g/2)(a†σ− − aσ+)）；基底序號 2n + s（s = 0 ↓、1 ↑）。"""
    dim = 2 * (Nmax + 1)
    Hm = np.zeros((dim, dim), complex)
    for n in range(Nmax + 1):
        Hm[2 * n, 2 * n] = n * Om
        Hm[2 * n + 1, 2 * n + 1] = n * Om + eps
        if n >= 1:                                   # ⟨n,↓|a†σ−|n−1,↑⟩ = √n
            val = (1j if book_phase else 1.0) * g / 2 * math.sqrt(n)
            Hm[2 * n, 2 * (n - 1) + 1] = val
            Hm[2 * (n - 1) + 1, 2 * n] = np.conj(val)
    return Hm


def jc_block(N: int, Om: float, eps: float, g: float) -> np.ndarray:
    """(8.7.10)：{|N−1⟩|↑⟩, |N⟩|↓⟩}：[[(N−1)Ω + ε, g√N/2], [g√N/2, NΩ]]。"""
    return np.array([[(N - 1) * Om + eps, g * math.sqrt(N) / 2], [g * math.sqrt(N) / 2, N * Om]])


def jc_block_eig(N: int, Om: float, eps: float, g: float) -> tuple:
    c = (N - 1) * Om + (eps + Om) / 2
    w = math.sqrt((eps - Om) ** 2 + g ** 2 * N) / 2
    return c - w, c + w


def jc_levels(Nmax: int, Om: float, eps: float, g: float) -> dict:
    """完整截斷矩陣的數值本徵值 ＝ 各 2×2 塊公式（＋基態 0）；book_phase 的 i 不改本徵值。"""
    w = np.sort(np.linalg.eigvalsh(jc_matrix(Nmax, Om, eps, g)))
    wb = np.sort(np.linalg.eigvalsh(jc_matrix(Nmax, Om, eps, g, book_phase=True)))
    form = [0.0] + [e for N in range(1, Nmax + 1) for e in jc_block_eig(N, Om, eps, g)] + [Nmax * Om + eps]
    return dict(ok=_close(w, np.sort(form), 1e-9) and _close(w, wb, 1e-9), levels=w)


def jc_splitting(N: int, g: float, Om: float = 1.0) -> float:
    """共振（Δ = 0）時第 N 塊的分裂 g√N（數值 2×2 對角化）。"""
    w = np.linalg.eigvalsh(jc_block(N, Om, Om, g))
    return float(w[1] - w[0])


def mollow(Om: float, Delta: float, r: float) -> dict:
    """(8.7.11)、(8.7.13)：Ω_F = Ω + ω、Ω_3 = Ω − ω（ω = √(Δ² + r²)）；遠離共振 Ω_F → ε = Ω + Δ、Ω_3 → 2Ω − ε；照印 CSRS 的 −δ/2 給 Ω − Δ。"""
    om = math.hypot(Delta, r)
    dl = om - Delta
    N = 50
    up = lambda N_: N_ * Om + Delta + dl / 2
    lo = lambda N_: N_ * Om - dl / 2
    fl = up(N) - lo(N - 1)                          # ε^D(2,N,N−1) − ε^D(1,N−1,N−2)
    cs = lo(N) - up(N - 1)                          # 更正的 +δ/2
    cs_book = lo(N) - ((N - 1) * Om + Delta - dl / 2)
    return dict(F=Om + om, L=Om, three=Om - om, F_levels=fl, three_levels=cs, three_book=cs_book, eps=Om + Delta, three_far=Om - Delta,
                res_lines=(Om - r, Om + r))


def mollow_jc(N: int = 100, g: float = 0.04) -> dict:
    """JC 本身（r = g√N 隨 N 變）：螢光 ＝ 第 N 塊上能階 − 第 N−1 塊下能階、CSRS ＝ 第 N 塊下能階 − 第 N−1 塊上能階（用 (8.7.10) 的 NΩ1 + … 形式）。"""
    Om, Delta = JC_CARD["Omega"], JC_CARD["Delta"]
    eps = Om + Delta
    blk = lambda n: np.linalg.eigvalsh(jc_block(n, Om, eps, g))
    bN, bM = blk(N), blk(N - 1)
    return dict(F=float(bN[1] - bM[0]), three=float(bN[0] - bM[1]), r=g * math.sqrt(N))


def jc_ladder(Delta: float, g: float = 0.1, Nmax: int = 3, Om: float = 1.0) -> list:
    """Fig. 8.7.1：[(N, 下穿衣能量, 上穿衣能量, 裸 |N⟩|↓⟩ 能量, 裸 |N−1⟩|↑⟩ 能量)]（N = 1…Nmax；N = 0 是基態 0）。"""
    eps = Om + Delta
    out = [(0, 0.0, 0.0, 0.0, float("nan"))]
    for N in range(1, Nmax + 1):
        lo, hi = jc_block_eig(N, Om, eps, g)
        out.append((N, lo, hi, N * Om, (N - 1) * Om + eps))
    return out


# ======================================================================== 全書回顧
RECAP = [
    (1, "矩陣代數與量子力學", "振幅、投影算符與本徵值：分析器實驗→ Hamilton–Cayley 譜分解", r"M=\begin{pmatrix}4&1\\3&2\end{pmatrix}:\ \varepsilon=1,5;\ \cos^2 30^\circ=3/4", "H101–H109"),
    (2, "Abel 對稱", "C₂、C_n 的投影算符解耦合振子；環上的 Bloch 波", r"\omega_\pm=1,\ \tfrac{11}{10};\ \text{beat }10\pi=31.4159", "H201–H212"),
    (3, "非 Abel 對稱", "類代數、冪等元、特徵標；D₃ 與氨分子", r"\mathrm{NH_3}:\ 3337,\ 908,\ 3414,\ 1628\ \mathrm{cm^{-1}}", "H301–H312"),
    (4, "高階有限對稱與誘導表示", "O、O_h 的子群鏈、Frobenius 互反、XF₆ 振動", r"\mathrm{SF_6}:\ T_{1u}^+\ 955,\ E_g\ 642\ \mathrm{cm^{-1}}", "H401–H413"),
    (5, "連續轉動群", "Euler／Darboux、𝒟^j、轉子與球諧；O 對稱分裂", r"(50,60,70)\to[80.00,\ 33.69,\ 128.68]", "H501–H516"),
    (6, "表示的乘積", "外積群、張量、應力應變、Lorentz 振子與光", r"\text{Lorentz: }\rho=90^\circ,\ R=0.6366\ \mathrm{m}", "H601–H614"),
    (7, "連續群的乘積", "耦合係數、3j／6j、Wigner–Eckart、RE 曲面、分子軌域", r"a=1422.17\ \mathrm{MHz},\ \lambda=21.106\ \mathrm{cm}", "H701–H714"),
    (8, "半古典與量子動力學", "作用量與波前、相干態、Rabi 與穿衣態、量子化的光", r"c_5=0.4189,\ A_{2p}=6.27\times10^{8}\ \mathrm{s^{-1}},\ g\sqrt N", "H801–H813"),
]


def recap_check() -> dict:
    """回顧表裡能便宜重算的數：第 1 章 [[4,1],[3,2]] 的本徵值、第 5 章 (50,60,70) 的 Darboux 角、第 7 章 21 cm 的 a、第 2 章拍的時間、第 8 章 c5 與 A。"""
    ev = sorted(float(np.real(e)) for e in np.linalg.eigvals(np.array([[4.0, 1], [3, 2]])))
    dar = [float(x) for x in RT.euler_to_darboux(50, 60, 70)]
    hn = SC.hyperfine_numbers()
    return dict(ch1=ev == [1.0, 5.0] or _close(ev, [1, 5]), ch5=_close(np.round(dar, 2), [80.00, 33.69, 128.68], 1e-9),
                ch7=f"{hn['a_MHz']:.2f}" == "1422.17" and f"{hn['wavelength_cm']:.3f}" == "21.106",
                ch2=_f(10 * PI) == "31.4159", ch8=_f(abs(coherent_cn(math.sqrt(5), 6)[5])) == "0.4189" and f"{h2p_A()['A']:.2e}" == "6.27e+08")


# ======================================================================== 卡片數字
def card_values() -> dict:
    """plan_ch08 §6 每張卡上的數字與導讀第 8 章的數字（字串；小數照卡上位數）。"""
    out = {}
    # H801
    f, fp = LINE_CARD["f"], LINE_CARD["fp"]
    env = line_envelope(f, fp, LINE_CARD["xjs"])
    out["H801_env"] = " / ".join(f"({_f(a, 3)}, {_f(b, 3)})" for a, b in env)
    lc = lagrange_card()
    out["H801_legendre"] = f"p={_f(lc['p'], 1)} H={_f(lc['H'], 1)} intercept={_f(lc['intercept'], 1)}"
    # H802
    fc = fall_card()
    out["H802_fall"] = f"v0={_f(fc['v0'], 1)} apex={_f(fc['apex'], 1)} T={_f(fc['T_apex'], 1)} range={_f(fc['range'], 1)}"
    out["H802_caustic"] = f"y={_f(fc['caustic_c'], 0)}-x^2/{_f(fc['caustic_k'], 0)}"
    out["H802_45"] = f"T={_f(fc['T45'])} SH={_f(fc['SH45'])} SP={_f(fc['SP45'])}"
    out["H802_apexSH"] = _f(fc["apex_SH"])
    out["H802_levels"] = " ".join(str(v) for v in FIG816_LEVELS)
    # H803
    out["H803_wkb"] = f"{_f(wkb_ratio(5, 1.0))} {_f(math.sqrt(11))}"
    hc = huygens_card()
    out["H803_huygens"] = f"{_f(hc['theta'], 2)} {_f(hc['optimal'], 1)} {_f(hc['deviant'], 1)}"
    # H804
    al = math.sqrt(COH_CARD["alpha2"])
    c = coherent_cn(al, 15)
    out["H804_cn"] = " ".join(_f(abs(c[n])) for n in range(11))
    st = coherent_stats(al)
    out["H804_stats"] = f"{_f(st['mean'])} {_f(st['dn'])}"
    out["H804_x0"] = _f(coh_x0())
    mo = packet_moments(PI / 2, coh_x0())
    out["H804_quarter"] = f"{_f(mo['x'] + 0.0)} {_f(mo['p'])}"
    en = energies(al)
    out["H804_E"] = f"{_f(en['Ec'], 1)} {_f(en['Eq'], 1)}"
    out["H804_bch"] = f"{_f(math.exp(-COH_CARD['alpha2'] / 2))} {_f(math.exp(COH_CARD['alpha2'] / 2))}"
    out["H804_uncert"] = _f(0.5)
    # H805
    se = snapshot_error(coh_x0(), 5, 10)
    se12 = snapshot_error(coh_x0(), 5, 12)
    seb = snapshot_error(coh_x0(), 5, 10, sign=-1)
    out["H805_err"] = f"{_f(se['ratio'])} {_f(se['l2'])} {_f(se12['l2'])}"
    out["H805_book"] = _f(seb["l2"])
    out["H805_p5"] = _f(abs(c[5]) ** 2)
    # H806
    hf = ho_ftac_card()
    out["H806_ho"] = f"{_f(hf['peaks'][4][0], 2)} {_f(hf['peaks'][5][0], 2)} {_f(hf['peaks'][5][1])}"
    out["H806_acpi"] = _e(hf["ac_pi"], 2)
    hm = ho_maslov()
    out["H806_maslov"] = f"{_f(hm['E_FT'], 2)} {_f(hm['E_SC'], 2)}"
    re_ = rotor_exact()
    out["H806_K8"] = " ".join(_f(e) for e in re_["pair"])
    out["H806_packet"] = f"{_f(re_['E_packet'])} {_f(re_['E_classical'], 3)}"
    out["H806_local"] = " ".join(_f(v, 3) for v in rotor_localized())
    rs = rotor_sc()
    out["H806_sc"] = f"{_f(rs['E_FT'], 3)} {_f(rs['E_SC'], 3)}"
    out["H806_coeffs"] = " ".join(_f(v, 3) for v in rs["coeffs"])
    out["H806_ebk"] = _f(rotor_ebk(8), 3)
    out["H806_ftexact"] = _f(re_["E_FT_exact"], 3)
    # H807
    out["H807_shift"] = _e(boost_shift(), 3)
    out["H807_ratio"] = _f(ap_er_ratio(1.0, 0.5), 1)
    # H808
    scd = spectral_card()
    out["H808_S"] = f"{_f(scd['peak'], 1)} {_f(scd['zero'])} {_f(scd['area'])}"
    out["H808_side"] = f"{_f(scd['side_D'])} {_f(scd['side_h'])} {_f(scd['book_side_D'])} {_f(scd['book_side_h'])}"
    hs = h_strengths()
    out["H808_H"] = f"{_f(hs['f2'])} {_f(hs['f3'])} {_f(hs['f4'])} {_f(hs['bound'])} {_f(hs['cont'])} {_f(hs['total'])}"
    out["H808_ho"] = " ".join(_f(v, 1) for v in ho_strengths()[:3])
    # H809
    rk = rotating_check()
    out["H809_crank"] = f"(r, 0, Delta) {'ok' if rk['ok'] else 'NG'}"
    # H810
    rc = rabi_card()
    out["H810_rabi"] = f"{_f(rc['t_pi'])} {_f(rc['p_res'])} {_f(rc['t_half'])} {_f(rc['p_half'])} {_f(rc['fwhm'], 1)}"
    bs = bloch_siegert(TL_CARD["eps"], TL_CARD["r"])
    out["H810_bs"] = f"{_f(bs['shift'])} {_f(bs['formula'])} {_f(bs['book_8525'])}"
    out["H810_bsfrac"] = _e(bs["book_example"], 2)
    # H811
    tz = thermal_sz(2.0)
    out["H811_thermal"] = _f(tz["formula"])
    b = BLOCH_CARD
    ss = bloch_steady(b["Delta"], b["r"], b["g1"], b["g2"], b["Sz0"])
    out["H811_bloch"] = f"{_f(ss[2])} {_f(0.5 + ss[2])} {_f(b['r'] ** 2 / (b['g1'] * b['g2']), 1)}"
    ak = ac_stark_card()
    out["H811_stark"] = f"{_f(ak['split'], 2)} {_f(ak['ac_shift'], 2)} {_f(ak['dc_shift'])}"
    out["H811_circles"] = " / ".join(
        f"{_f(cg['delta'])} {_f(cg['deltap'])} {_f(cg['theta'], 2)}" for cg in (circle_geometry(dd, FIG853["r"]) for dd in FIG854_DELTAS))
    ad = adiabatic(**ADIAB)
    out["H811_adiab"] = f"{_f(ad['inversion'], 2)} {_f(ad['lz'], 2)}"
    # H812
    pe = photon_elements(3)
    out["H812_elements"] = f"{_f(pe['absorb'])} {_f(pe['emit'])}"
    out["H812_ang"] = _f(angular_sum())
    hA = h2p_A()
    out["H812_A"] = f"{_f(hA['r2'])} {_e(hA['A'])} {_e(hA['A_red'])} {_f(hA['tau_ns'], 3)} {_e(hA['book_printed'])}"
    # H813
    out["H813_split"] = " ".join(_f(jc_splitting(N, 0.1), 1) for N in (1, 4, 9, 16))
    ml = mollow(JC_CARD["Omega"], JC_CARD["Delta"], JC_CARD["r"])
    out["H813_mollow"] = f"{_f(ml['F'])} {_f(ml['L'])} {_f(ml['three'])}"
    out["H813_far"] = f"{_f(ml['eps'], 1)} {_f(ml['three_far'], 1)} {_f(ml['three_book'], 1)}"
    out["H813_res"] = f"{_f(ml['res_lines'][0], 1)} {_f(ml['res_lines'][1], 1)}"
    mj = mollow_jc()
    out["H813_jc"] = f"{_f(mj['F'])} {_f(mj['three'])}"
    out["H813_recap"] = " | ".join(f"{row[0]} {row[4]}" for row in RECAP)
    return out


CARD_EXPECTED = {
    'H801_env': '(1.000, 0.500) / (-1.500, 1.125)',
    'H801_legendre': 'p=2.0 H=1.0 intercept=-1.0',
    'H802_fall': 'v0=16.0 apex=32.0 T=4.0 range=64.0',
    'H802_caustic': 'y=32-x^2/128',
    'H802_45': 'T=5.6569 SH=965.4365 SP=241.3591',
    'H802_apexSH': '341.3333',
    'H802_levels': '50 100 150 200 250 300 350 400 450 500 600 700 800 1000',
    'H803_wkb': '0.0316 3.3166',
    'H803_huygens': '75.52 20.0 30.0',
    'H804_cn': '0.0821 0.1835 0.2902 0.3747 0.4189 0.4189 0.3824 0.3232 0.2555 0.1904 0.1347',
    'H804_stats': '5.0000 2.2361',
    'H804_x0': '3.1623',
    'H804_quarter': '0.0000 -3.1623',
    'H804_E': '5.0 5.5',
    'H804_bch': '0.0821 12.1825',
    'H804_uncert': '0.5000',
    'H805_err': '0.0299 0.0299 0.0091',
    'H805_book': '1.4276',
    'H805_p5': '0.1755',
    'H806_ho': '4.50 5.50 0.1755',
    'H806_acpi': '4.54e-05',
    'H806_maslov': '5.00 5.50',
    'H806_K8': '53.1461 53.1489',
    'H806_packet': '53.5000 50.000',
    'H806_local': '0.129 0.922 -0.347 0.106 -0.031 0.009',
    'H806_sc': '53.410 53.135',
    'H806_coeffs': '0.119 0.929 -0.341 0.076 -0.011 0.001',
    'H806_ebk': '53.214',
    'H806_ftexact': '53.147',
    'H807_shift': '1.758e-15',
    'H807_ratio': '2.0',
    'H808_S': '16.0 1.5708 25.1327',
    'H808_side': '2.2467 0.0472 2.3562 0.0450',
    'H808_H': '0.4162 0.0791 0.0290 0.5650 0.4350 1.0000',
    'H808_ho': '1.0 0.0 0.0',
    'H809_crank': '(r, 0, Delta) ok',
    'H810_rabi': '15.7080 1.0000 11.1072 0.5000 0.4',
    'H810_bs': '0.0100 0.0100 0.0099',
    'H810_bsfrac': '2.50e-13',
    'H811_thermal': '-0.3808',
    'H811_bloch': '-0.0556 0.4444 8.0',
    'H811_stark': '0.50 0.25 0.0590',
    'H811_circles': '0.1000 2.5000 22.62 / 0.3385 0.7385 68.20 / 0.6099 0.4099 101.31 / 2.5000 0.1000 157.38',
    'H811_adiab': '0.98 0.98',
    'H812_elements': '1.7321 2.0000',
    'H812_ang': '8.3776',
    'H812_A': '0.5549 6.2683e+08 6.2649e+08 1.595 3.9385e+09',
    'H813_split': '0.1 0.2 0.3 0.4',
    'H813_mollow': '10.5000 10.0000 9.5000',
    'H813_far': '10.3 9.7 9.7',
    'H813_res': '9.6 10.4',
    'H813_jc': '10.4992 9.5008',
    'H813_recap': '1 H101–H109 | 2 H201–H212 | 3 H301–H312 | 4 H401–H413 | 5 H501–H516 | 6 H601–H614 | 7 H701–H714 | 8 H801–H813',
}


# ======================================================================== check()
def _check_81():
    f, fp = LINE_CARD["f"], LINE_CARD["fp"]
    Xs = np.linspace(-2, 2, 41)
    assert _close(legendre_numeric(f, Xs, np.linspace(-6, 6, 120001)), Xs ** 2 / 2, 1e-8)      # y = x²/2 自對偶
    env = line_envelope(f, fp, np.linspace(-2, 2, 9))
    lines = line_family(f, np.linspace(-3, 3, 60001))
    sl = np.array([a for a, _ in lines])
    ic = np.array([b for _, b in lines])
    assert all(abs(np.max(sl * X + ic) - Y) < 1e-8 for X, Y in env)                              # 包絡 ＝ 直線族的上包絡
    ce = circle_envelope(f, fp, 0.7, np.linspace(-2, 2, 9))
    assert all(abs(_dist_to_curve(P, f) - 0.7) < 1e-8 for P in ce)
    lc = lagrange_card()
    assert abs(lc["H"] - lc["H_numeric"]) < 1e-8 and abs(lc["H"] - lc["H_formula"]) < 1e-12
    fcs = fall_checks()
    assert all(fcs.values()), fcs
    fc = fall_card()
    assert abs(fc["SH45"] - fc["SH45_num"]) < 1e-8 and abs(fc["apex_SH"] - 1024 / 3) < 1e-9 and abs(fc["apex_SH_num"] - 1024 / 3) < 1e-7
    # 焦散 ＝ 軌跡族的包絡（每個 x 取所有仰角的最高點）
    xs = np.linspace(-60, 60, 12)                      # 不含 x = 0
    v0, g = fall_v0(), -FALL["F"] / FALL["m"]
    ths = np.linspace(0.5, 179.5, 20001)
    for xx in xs:
        c = np.cos(np.radians(ths))
        T = xx / (v0 * c)
        ok = T > 0
        y = v0 * np.sin(np.radians(ths[ok])) * T[ok] - g * T[ok] ** 2 / 2
        assert abs(np.max(y) - caustic(xx)) < 1e-3
    assert orthogonality(200.0, np.linspace(-150, 150, 13)) < 1e-6
    assert orthogonality(450.0, np.linspace(10, 170, 9)) < 1e-6
    sc = semiclassical_check()
    assert sc["ok"] and sc["hj_ok"]
    assert abs(wkb_ratio(5, 1.0) - 1 / 10 ** 1.5) < 1e-12
    hc = huygens_card()
    assert abs(hc["deviant"] - 30) < 1e-12


def _check_82():
    assert all(bch_check().values())
    assert all(bch_check(-0.4, 1.1).values())
    x0 = coh_x0()
    al = alpha_of(x0, 0.0)
    cnum = coherent_numeric(x0, 0.0)
    assert _close(cnum[:30], coherent_cn(al, 29), 1e-9)
    st = coherent_stats(al)
    vec = cnum
    nn = np.arange(NBASIS)
    mean = float(abs(vec) ** 2 @ nn)
    assert abs(st["mean"] - 5) < 1e-10 and abs(mean - 5) < 1e-8 and abs(st["dn"] - math.sqrt(5)) < 1e-10
    xg = np.linspace(-9, 9, 721)
    for t in (0.0, 0.7, PI / 2, 2.9):
        assert _close(packet_closed(xg, t, x0), packet_sum(xg, t, x0), 1e-9)
    xs_, ps_ = packet_split_operator(2.0, x0)
    assert _close(ps_, packet_closed(xs_, 2.0, x0), 5e-5)                      # dt = 0.005 的二階分裂誤差 ~1.5e−5
    assert abs(packet_uncertainty(1.1, x0) - 0.5) < 1e-4
    # ⟨x⟩、⟨p⟩ ＝ x_t、p_t（數態基底的期望值）
    xo, po = xp_ops(NBASIS)
    for t in (0.4, 2.0):
        v = coherent_cn(al * np.exp(-1j * t), NBASIS - 1)
        mo = packet_moments(t, x0)
        assert abs(np.real(v.conj() @ xo @ v) - mo["x"]) < 1e-9 and abs(np.real(v.conj() @ po @ v) - mo["p"]) < 1e-9
    assert abs(SH_ho(1.3, x0) - SH_ho_numeric(1.3, x0)) < 1e-9
    en = energies(al)
    assert abs(en["Eq"] - en["Ec"] - 0.5) < 1e-10


def _check_83():
    x0 = coh_x0()
    se = snapshot_error(x0, 5, 10)
    assert abs(se["l2"] - se["pred"]) < 1e-6 and abs(se["ratio"] - math.sqrt(5 ** 10 * math.factorial(5) / math.factorial(15))) < 1e-12
    assert _f(se["ratio"]) == "0.0299"
    se12 = snapshot_error(x0, 5, 12)
    assert abs(se12["l2"] - se12["pred"]) < 1e-6 and se12["l2"] < 0.01
    seb = snapshot_error(x0, 5, 10, sign=-1)
    assert seb["l2"] > 1.0
    al = math.sqrt(5)
    ts = np.linspace(0, 9, 50)
    assert _close(ac_ho(ts, al), ac_ho_sum(ts, al), 1e-10)
    hf = ho_ftac_card()
    c = coherent_cn(al, 15)
    ts_ = np.linspace(0, 20 * PI, 8001)
    at_levels = ftac(ts_, ac_ho(ts_, al), np.arange(12) + 0.5)          # 在 E = n + ½ 其他譜線的 sinc 正好為零
    assert _close(at_levels, abs(c[:12]) ** 2, 1e-6)
    for n, (e, h) in enumerate(hf["peaks"]):
        if 3 <= n <= 7:
            assert abs(e - (n + 0.5)) < 3e-3 and abs(h - abs(c[n]) ** 2) < 1e-3
    hm = ho_maslov()
    assert hm["ac_ratio_ok"] and abs(hm["E_FT"] - 5.0) < 1e-3 and abs(hm["E_SC"] - 5.5) < 1e-3
    re_ = rotor_exact()
    assert _close(re_["pair"], re_["pair_sc"], 1e-9) and abs(re_["E_classical"] - 50.0) < 1e-12 and abs(re_["E_packet"] - 53.5) < 1e-9
    assert abs(re_["E_FT_exact"] - re_["centroid"]) < 2e-3
    loc = rotor_localized()
    assert [round(v, 3) for v in loc[:4]] == [0.129, 0.922, -0.347, 0.106]
    rs = rotor_sc()
    assert abs(rs["E_cs"] - 53.5) < 1e-9 and rs["norm_drift"] < 1e-8 and rs["E_drift"] < 1e-7 and abs(rs["E_SC"] - 53.1475) < 0.05
    assert max(abs(a - b) for a, b in zip(rs["coeffs"], (0.118, 0.929, -0.339, 0.075, -0.012, 0.001))) < 0.005   # 書上 (8.3.31)


def _check_84():
    gc = gauge_check()
    assert gc["ok"] and gc["book_fails"]
    assert hamiltonian_check()
    assert f"{boost_shift():.3e}" == "1.758e-15"
    scd = spectral_card()
    assert abs(scd["peak"] - 16) < 1e-12 and abs(scd["zero"] - PI / 2) < 1e-10 and abs(scd["area"] - scd["area_formula"]) < 2e-3
    assert scd["closed_vs_quad"]
    foc = first_order_check()
    assert foc["ok"] and foc["iS"]
    assert abs(golden_rule_check() - 1) < 2e-3
    hs = ho_strengths()
    assert abs(hs[0] - 1) < 1e-12 and all(abs(v) < 1e-12 for v in hs[1:])
    h = h_strengths()
    assert abs(h["total"] - 1) < 1e-6
    hb = h_box_strengths()
    assert abs(hb["E1s"] + 0.5) < 1e-3 and abs(hb["f2"] - h["f2"]) < 2e-3 and abs(hb["f3"] - h["f3"]) < 2e-3 and abs(hb["total"] - 1) < 2e-3
    assert classical_factor_check()["ok"]


def _check_85():
    r, eps = TL_CARD["r"], TL_CARD["eps"]
    ts = np.linspace(0, 30, 301)
    Om = 0.93
    psi = evolve_lab(ts, r, Om, eps, psi0=(0.6, 0.8))
    Sq = np.array([SC.spin_vector(p) for p in psi])
    Sc = spin_precession(ts, r, Om, eps, S0=Sq[0])
    assert _close(Sq, Sc, 1e-7)
    assert rotating_check()["ok"]
    for Dl in (0.0, 0.1, -0.25, 0.6):
        assert _close(rabi_numeric(Dl, r, ts), rabi(Dl, r, ts), 1e-7)
    big = 3.0
    assert abs(float(rabi(big, 0.01, 5.0)) - 0.01 ** 2 / 4 * float(S2(big, 5.0))) < 1e-8      # (8.5.18)
    rc = rabi_card()
    assert abs(rc["p_res"] - 1) < 1e-12 and abs(rc["p_half"] - 0.5) < 1e-12 and abs(rc["env_half"] - 0.5) < 1e-12
    bs = bloch_siegert(eps, r)
    assert abs(bs["shift"] - bs["next_order"]) < 3e-6 and abs(bs["pop_min"] - 0.5) < 1e-6 and _f(bs["shift"]) == "0.0100"
    assert abs(bloch_siegert(eps, 0.1)["shift"] - 0.1 ** 2 / 4) < 2e-6
    assert bs_sign_check()["ok"]
    b = BLOCH_CARD
    for Dl in (0.0, 0.15):
        assert _close(bloch_numeric(Dl, b["r"], b["g1"], b["g2"], b["Sz0"]), bloch_steady(Dl, b["r"], b["g1"], b["g2"], b["Sz0"]), 1e-7)
    tz = thermal_sz(2.0)
    assert abs(tz["numeric"] - tz["formula"]) < 1e-12
    for Dl in (1.2, 0.2, -0.1, -1.2):
        assert dressed(Dl, 0.5)["ok"] and circle_geometry(Dl, 0.5)["ok"]
    ak = ac_stark_card()
    assert abs(ak["split"] - ak["split_num"]) < 1e-12 and abs(ak["split"] - FIG853["r"]) < 1e-12
    assert abs(ak["dc_shift"] - ak["dc_formula"]) < 1e-12 and abs(ak["dc_shift"] - ak["dc_num"]) < 1e-12
    Oms = np.linspace(0, 2, 21)
    up, lo = ac_stark(Oms, 1.0, 0.5)
    for O_, u, l in zip(Oms, up, lo):
        w = np.linalg.eigvalsh(np.array([[1.0, 0.25], [0.25, O_]]))
        assert abs(w[1] - max(u, l)) < 1e-12 and abs(w[0] - min(u, l)) < 1e-12
    ad = adiabatic(**ADIAB)
    assert abs(ad["inversion"] - ad["lz"]) < 0.003


def _check_86():
    pe = photon_elements(3)
    assert abs(pe["absorb"] - math.sqrt(3)) < 1e-12 and abs(pe["emit"] - 2) < 1e-12
    me = mode_energy_check()
    assert me["ok"]
    assert abs(angular_sum() - 8 * PI / 3) < 1e-12 and abs(angular_sum_mc() - 8 * PI / 3) < 0.02
    hA = h2p_A()
    assert abs(hA["r2"] - hA["r2_num"]) < 1e-10
    assert abs(hA["A"] / hA["A_cgs"] - 1) < 1e-9 and abs(hA["A"] / hA["A_au"] - 1) < 1e-8
    assert f"{hA['A']:.4e}" == "6.2683e+08" and f"{hA['A_red']:.4e}" == "6.2649e+08"
    assert abs(golden_A_numeric() / hA["A"] - 1) < 3e-3
    nc = number_vs_coherent()
    assert nc["num_a"] < 1e-15 and abs(nc["num_x2"] - 5.5) < 1e-12 and abs(nc["coh_a"] - math.sqrt(5)) < 1e-9


def _check_87():
    for Dl in (0.3, 0.0, -0.2):
        assert jc_levels(30, 1.0, 1.0 + Dl, 0.1)["ok"]
    for N in (1, 4, 9, 16):
        assert abs(jc_splitting(N, 0.1) - 0.1 * math.sqrt(N)) < 1e-12
    ml = mollow(JC_CARD["Omega"], JC_CARD["Delta"], JC_CARD["r"])
    assert abs(ml["F"] - ml["F_levels"]) < 1e-9 and abs(ml["three"] - ml["three_levels"]) < 1e-9
    assert abs(ml["three_book"] - (JC_CARD["Omega"] - JC_CARD["Delta"])) < 1e-9
    mj = mollow_jc()
    H = jc_matrix(102, JC_CARD["Omega"], JC_CARD["Omega"] + JC_CARD["Delta"], 0.04)
    w = np.sort(np.linalg.eigvalsh(H))
    # 從完整矩陣取第 100、99 塊（每塊兩個能階，能量 ~ NΩ）
    blkN = sorted(e for e in w if abs(e - 100 * JC_CARD["Omega"]) < 2)
    blkM = sorted(e for e in w if abs(e - 99 * JC_CARD["Omega"]) < 2)
    assert abs((blkN[1] - blkM[0]) - mj["F"]) < 1e-8 and abs((blkN[0] - blkM[1]) - mj["three"]) < 1e-8
    rc = recap_check()
    assert all(rc.values()), rc


def _check_cards():
    cv = card_values()
    for k, v in CARD_EXPECTED.items():
        assert cv[k] == v, (k, cv[k], v)
    assert set(cv) == set(CARD_EXPECTED), set(cv) ^ set(CARD_EXPECTED)


def summary() -> dict:
    return card_values()


def check():
    """每個上卡片／講義的數字用兩種以上獨立算法核對（約 20 秒）。"""
    _check_81()
    _check_82()
    _check_83()
    _check_84()
    _check_85()
    _check_86()
    _check_87()
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
    print(f"sym_dyn（第 8 章）：包絡與 Legendre、拋射體的 S_P／S_H 與正交波前、半古典極限、BCH 與相干態（Poisson、波包解析＝數值）、"
          f"快照合成 ⟨x|5⟩（誤差 0.0299）、自相關與 FTAC（諧振子、J = 10 轉子 K = 8 叢）、Maslov 與 EBK、A·p／E·r 規範、譜函數與 TRK、"
          f"Rabi、Bloch–Siegert（Floquet）、Bloch 方程、穿衣態與 AC Stark、Einstein A、Jaynes–Cummings 與 Mollow 全部核對通過"
          f"（{time.time() - t0:.1f} 秒）。")
