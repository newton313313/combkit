# 由 combkit/tools/sync_from_harter_sym.py 從 harter_sym 專案的 codes/utils/sym_prod.py 複製；請改 harter_sym 那份再同步。
r"""
sym_prod.py —— Harter《Principles of Symmetry, Dynamics, and Spectroscopy》（Wiley 1993）第 6 章
Theory and Applications of Symmetry Representation Products 的數學（純 Python＋numpy＋sympy＋scipy，不 import manim；
有限群工具借用同資料夾的 sym_nonabel（FiniteGroup、C3V、QUAT、rot2/refl2/rot3）、sym_oct（O、O_h、T、irreps 與子群鏈、六井、XF6）、
sym_rot（O3 ↓ O、χ^j），不重寫）。書上的式號、筆誤與〔注〕見 ref/md/ch06_formulas.md（相對於 C:\mypython\manim\projects\harter_sym）。

約定
  * 群乘積 ab = 先 b 後 a；ħ = 1；卡片上的小數一律 f"{float(x):.4f}"（Fraction、sympy 數要 float()）。精確的表用 exact(x)／exact_mat(M)
    （浮點 → ±√(p/q) 的 sympy 數）。數值的 irreps：IrrepSet（.elems .names .D[α][g] .dim .chi .stack .decompose(χ) → {α: 頻率}）。
  * irrep 組：oset(basis)（O：A_1 A_2 E T_1 T_2）、ohset(basis)（O_h 十個，"T_{1u}"…）；basis = "tetragonal"（O⊃D4⊃D2，kosher T2 (4.2.16)，
    ＝ sym_oct 的 TETRA）、"nonkosher"（(4.2.15)：𝒟^{T2} = χ^{A2}𝒟^{T1}，基底 yz, xz, xy；(6.1.19c)、(6.3.6) 用這個）、"trigonal"（O⊃D3⊃C2）、
    "trigonal_book"（(6.1.36) 的 T2：E 那一對比 sym_oct 差號）、"c4"（行進波）。E 的基底 (|E1⟩, |E2⟩) ~ ((2z²−x²−y²)/√6, (x²−y²)/√2)。
    d4set()（D4 = {1, R3², R3, R3³, R1², R2², i3, i4}：A_1 B_1 A_2 B_2 E，書上 (6.1.21) 的順序）、d2set()（D2：A_1 A_2 B_1 B_2，(4.2.10)：
    B1 對 R1² 為 +1）、d3set()（D3 = {1, r1, r1², i2, i4, i5}，C2 = {1, i4}）、c3set()（0_3 1_3 2_3，ε^k）、tset()（T：A ^1E ^2E T）、
    qset()（四元數群的二維 irrep "E"）、c3vset()（SN.C3V 的平面幾何：σ₃ 沿 x、σ₂ 60°、σ₁ 120°；A_1 A_2 E）。
  * **耦合係數 cg(S, a, b, gammas) → (C, cols)**：C 的列 = |a i⟩|b j⟩（字典序 i·ℓ^b + j，(6.1.8)），行 cols = [(γ, ω, m)]（0 起算）。
    相位約定＝書上的 (6.1.18b)：每一份 γ 取使 N^{γab}_{nkl} ≠ 0 的第一組 (n,k,l)，C_{kl,n} = √N > 0，其餘行 = P^γ_{mn}|kl⟩/√N。
    用這個約定，書上的 (6.1.19c)、(6.1.23)、(6.1.26)、(6.1.30)、(6.1.35a,b)、(6.1.36) **逐格相同**（cg_book(name)、compare_table）。
  * 六井（Fig. 4.3.1、6.1.1）：SITES |1⟩…|6⟩ = +z, −z, +x, −x, +y, −y；Hamilton 量第一列 [H, −T, −S, −S, −S, −S]（(6.4.9)）⇒
    A1g = H − T − 4S、T1u = H + T、Eg = H − T + 2S；SIX_STATES[(γ, m)]（T1u 是 x, y, z；Eg ~ 2z²−x²−y², x²−y²）。
  * Lorentz：ẍ + 2Γẋ + ω0²x = A cos ω_s t（更正的 (6.5.1)）；e^{−iωt} 約定：複振幅 R e^{iρ}，ρ 是落後。

§6.1A 兩粒子與外積群（H601）
  PAIRS（36 個 (i, j)，(6.1.1)）、site_perm(g)（6×6）、outer_op(a, b) = kron（(6.1.4)）、outer_trace_check()（(6.1.6)、(6.1.9)）、
  outer_frequencies()（2304 個元素 → 九個 (α, β) 各一次：單粒子 A1g ⊕ Eg ⊕ T1u）
  kron_example()（(6.1.7ₓ)：𝒟^{T1g}(r2) = [[0,0,−1],[1,0,0],[0,−1,0]] ⊗ 𝒟^{Eu}(i5) = [[−1/2, −√3/2], [−√3/2, 1/2]]，6×6 ＝ 書上）
  (1, R1)|1,1⟩ = |1,6⟩（本系列的頂點編號；書上 (6.1.12) 寫 |1,3⟩ 是書上的編號）
  six_hamiltonian(H, S, T)、six_levels(H, S, T)、TWO_CARD = {H: 0, SA: 1, SB: 1/4}（本系列）：A 的能階 −4, 0, 2，B 的 −1, 0, 0.5；
  two_particle_levels(SA, SB) → [(E, α, β, 簡併)] ＝ two_particle_numeric（36×36 對角化）：2.5(4) 2(6) 1(2) 0.5(6) 0(9) −1(3) −3.5(2) −4(3) −5(1)
  （Fig. 6.1.2(a) 的上下次序需 SB < SA/3）
§6.1B 相互作用與耦合係數（H602）
  diag_frequencies(a, b)（(6.1.14)）：T1u⊗T1u = A1g + Eg + T1g + T2g；char_product_rows()（χ^{T1u}² = 9 0 1 1 1 9 0 1 1 1）
  cg_norms(S, a, b, γ)（(6.1.16c)）：T1g 的 N_111 = N_112 = N_113 = N_121 = N_122 = 0、N_123 = 1/2；cg_reduces（(6.1.19a,b)）
  CG_619C_BOOK、CG_630_BOOK、CG_623_BOOK、CG_626_BOOK、CG_635A_BOOK、CG_635B_BOOK、CG_636_BOOK（照書的表，numpy）
  cg_book("619c"/"630"/"623"/"626"/"635a"/"635b"/"636")、compare_table(C, cols, book) → {equal, flips, max}
  PAIR_V = {v0: 1, v1: 1/2, v2: 1/4, K: 0.1, W: 0.05}（本系列：直接項 V(0)、V(√2)、V(2)、近鄰交換 K、近鄰成對穿隧 W）、pair_V(**PAIR_V)、
  pair_V_invariant()（與 (a,a) 對易、與 (1, R1) 不對易；⟨1,1|V|1,1⟩ = 1、⟨1,6|V|1,6⟩ = 0.5）、
  t1u_split() ＝ t1u_split_numeric()：A1g 0.7250、Eg 0.5750、T2g 0.6000、T1g 0.4000（只有直接項：0.6250、0.6250、0.5000、0.5000）
§6.1C 子群鏈（H603）
  d2_table() ＝ D2_TABLE_BOOK（(6.1.20)：B1⊗B2 = A2…）；solve_6129()（(6.1.27b) 的 a…k：b = a/√2、d = −√2c、e = −√3c、g = −f、h = √2f、
  j = i、k = √2i；歸一化 a = √6/3、b = c = √3/3、d = −√6/3、e = −1、f = i = j = √2/2、g = −√2/2、h = k = 1 ⇒ ＝ (6.1.30)）；
  solve_6129(True)（加上書上的 d = c/√2 ⇒ 方程矛盾）
§6.2 一般關係（H604）
  fs_indicator(S, α)（(1/°G)Σχ(g²)）＝ irrep_type(S, α)（共軛／不變雙線性型）：O 全 1；C3：1, 0, 0；T：1, 0, 0, 1；Q：−1（第三型）
  sym_squares(S, α) → (對稱, 反對稱) ＝ cg_symmetry(S, α)：[E]² = A1 + E / A2、[T1]² = [T2]² = A1 + E + T2 / T1
  scalar_cg_check（(6.2.5)–(6.2.6)：δ_ij/√ℓ，T1u⊗T1u 的 A1g 行 √3/3）、unitarity_check（(6.2.12)–(6.2.14)）、factorization_lemma（(6.2.15)）
  T 群：T⊗T = A + ¹E + ²E + 2T（(6.2.10)；cg 給兩份正交的 T）；C3：1_3 ⊗ 2_3 = 0_3
§6.3A–C 向量、張量、多項式（H605）
  unit_tensors() → {(γ, m): 3×3}（(6.3.6)）、tensor_components(M)（(6.3.13a)）、tensor_from_components、tensor_transform_check（(6.3.5)）
  POLY_639、poly_from_tensors()、poly_transform_check（(6.3.8)）、o3_labels()：0⁺ → A1g、1⁺ → T1g、2⁺ → Eg + T2g；
  x² − y² 屬 Eg：𝒟^{Eg}(R3) = diag(1, −1)
§6.3D 應力與應變（H606）
  stress_examples()：等向 p·1 → 𝒯^{A1g} = √3 = 1.7321；純張力 ẑẑ → 𝒯^{A1g} = 0.5774、𝒯^{Eg}_1 = −0.8165；𝒯_12 = 𝒯_21 = 1 → 𝒯^{T2g}_3 = 1.4142；
  扭矩 𝒯_23 = −𝒯_32 = 1 → 𝒯^{T1g}_1 = 1.4142；face_forces(T)（(6.3.12)）、strain_of_displacement(s)（(6.3.14)）
§6.3E 彈性常數（H607）
  elastic_counts()（暴力法）＝ elastic_counts_chars()（張量定理）：沒有對稱 81 → 36 → 21；O_h：81 個分量 4 個不變量、36 → 3、21 → 3；
  等向 2；D4h 7 → 6；C3v 8 → 6
  elastic_tensor(kA, kE, kT2)（(6.3.26)）、elastic_from_cartesian（(6.3.28)、更正的 (6.3.28)ₓ）、moduli(kA, kE, kT2)（B、μ、Y、σ、A、k_{11,11}、
  k_{11,22}、k_{12,12}；Fraction 精確）＝ moduli_numeric（對 Cartesian 張量「做實驗」）
  spring_model(j, k, l)（(6.3.41)，單位 1/4a）＝ spring_model_bonds（鈣鈦礦晶格逐根彈簧加總）；MODULI_CARD：
    (1,1,1)：k = 9, 3, 4 → B 3、μ 2、Y 27/7 = 3.8571、σ 2/7 = 0.2857、A −0.5、k_{11,11} 5、k_{11,22} = k_{12,12} = 2（Cauchy）
    (2,1,1)：k = 10, 4, 4（k^E = k^{T2}，等向）→ B 10/3 = 3.3333、μ 2、Y 5、σ 0.25、A 0
  CUBIC_KITTEL、anisotropy(C11, C12, C44)（(6.3.43)）：鑽石 −0.35、Si −0.73、Ge −0.79、NaCl 0.88、Au −1.30、Pb −1.52、W −0.02（＝ Fig. 6.3.9 的標記）
§6.3F 場（H608）
  c3v_field(g, "E"/"B")（3×3；軸向量乘 det）、c3v_field_irreps()：E → A1 + E、B → A2 + E（跡 3 0 1 ／ 3 0 −1）；fig6311()（rE₁ = −½E₁ + (√3/2)E₂、
  σ₃E₂ = −E₂、σ₃B₁ = −B₁、σ₃B₂ = B₂、σ₃B₃ = −B₃）；p3_12("B") = −B̂₂（更正的 (6.3.47)）、p3_12_book() = (0, −3)（照印 3/2）；
  striction_counts() ＝ striction_counts_chars()：C3v 電致 4、磁致 2；O_h 0、0
§6.4 張量算子（H609）
  six_invariants()（(6.4.8)，更正的號）、SIX_INV_ROWS_BOOK（照印）、SIX_CARD = {H: 0, S: 1, T: 0.5} → 能階 −4.5, 0.5, 1.5；six_from_levels（(6.4.10)）
  param_count(freqs)（(6.4.11a)）＝ commutant_count(mats)：六井 3、SF6 振動 7、21 座標全部 11
  six_dipole_elements()、we_check()（三種躍遷的 _we_fit：E→T1u、A1g→T1u、T1u→Eg；T1u→Eg 的 R 更正 1.0000、照印 1/ℓ^α 0.6667）、
  WE_617_BOOK（(6.4.17) 比例 −1/2, −1/2, 1, √3/2, −√3/2, 0）、we_intensities()（極化 1：1:3、極化 3：1:0、(6.4.19) 的比 0.5）、
  we_complex_check()（T 群的 ¹E：共軛在 (6.4.15a) 的 C 上）；數值：⟨T1 3|z|E1⟩ = 0.8165、⟨T1 1|x|E1⟩ = −0.4082、⟨T1 1|x|E2⟩ = 0.7071、⟨T1 k|x_k|A1⟩ = 0.5774
§6.5A Lorentz 振子（H610）
  FIG651（ω0/2π = 5 Hz、Γ = 0.5、f_s = 4, 4.5, 4.9, 5, 5.25, 6 Hz）、LORENTZ_A = 20（本系列的刺激振幅）、lorentz_steady → (R, ρ)：
  ρ = 4.05°, 8.57°, 38.23°, 90.00°, 161.94°, 175.04°（書 4.0 … 175.0）、R = 0.0561, 0.1055, 0.4020, 0.6366, 0.1880, 0.0459 m；
  lorentz_transient（(6.5.7)–(6.5.8)）、lorentz_closed ＝ lorentz_numeric（scipy）、lorentz_trace(f_s, tmax, n) → (ts, 刺激, 響應)（畫圖用）
  response(ω)、smith(Δ, S, Γ)、smith_check()（Γ 圓直徑 S/Γ = 10）、fig652()（Im 最大 10.03 @0.995、Re 5.56 @0.894、−4.55 @1.095、DC 2）
§6.5B 物質對光（H611）
  epsilon(ω, ω0, ωp, Γ, book)、fig654()（ε(0) = 2、ω_c = √2、不透明帶 (1, √2)、ε(1.2) = −1.2727、ε(2) = 0.6667）；ε(1)（Γ = 0.1）= 1 + 5i，照印 1 + 2.5i
  BOOK_CONST、book_examples()（α_e 7.95e−30、α_v 6.6e−30、光頻 −4.2e−33 m³；Γ 4.4e7 s⁻¹、3/2Γ = 3.4e−8 s、2.0e7 次；紅外 3.4e−4 s、2.0e9 次；
  r_e 2.81e−15 m、πr² 2.5e−29 m²）、rayleigh()（450／700 nm：ω⁴ 律 5.86、含共振因子 6.17）、radiation_integral() = 4/3、
  phasor_points()（ρ = 0, 90, 180, −90 → 和波 1.0198 落後 11.31°、0.8、1.0198 超前、1.2）
§6.5C 紅外與 Raman（H612）
  xf6_selection()（振動 = A1g + Eg + T2g + 2T1u + T2u；IR 2T1u（ν3, ν4）、Raman A1g + Eg + T2g（ν1, ν2, ν5）、靜默 T2u（ν6），互斥）＝
  xf6_selection_numeric()（剛性離子偶極、鍵極化率的導數）；XF6_NU_LABELS；C^{T2u Eg T1u}_{323} = −1
§6.6 泛音（H613）
  sym_power_chi（Newton）＝ sym_power_chi_monomial（單項式基底）、overtone(α, n)：T2u：2 → A1g + Eg + T2g、3 → A1u + T1u + 2T2u、
  4 → 2A1g + 2Eg + T1g + 2T2g；Eg：2 → A1g + Eg、3 → A1g + A2g + Eg；pascal(ℓ, n)；combination("E_g", "T_{1u}") = T1u + T2u；henon_heiles_check()
§6.7 Jahn–Teller（H614）
  JT_CARD = {j: 1, k: 1, r: −0.2}、jt_matrix（(6.7.6a)）、jt_sheets（(6.7.8)）、jt_min（ρ0 = 0.5、−0.25）、renner_extrema（jr < 0：極小 60°, 180°, 300°，
  ρ = 0.625、−0.3125；鞍點 ρ = 0.4167、−0.2083；位障 0.1042）、jt_eigvec_track（繞一圈 −1，＝ (6.7.9)）、
  jt_vibronic_levels(j, k, r)（二維振子基底對角化 → [(E, 簡併, "E"/"A_1"/"A_2")]）＝ jt_radial_levels(l)（更正的 (6.7.15)）：
  r = 0：0.9758 (l=½)、2.0929 (3/2)、2.6872 (½)、3.2661 (5/2)…，全是二重；r = −0.2：E 0.9478、A1 1.7963、A2 2.3147、E 2.5684、E 3.1513；
  jt_effective(ρ, l, sheet, book_fig)、jt_effective_min（更正：l = ½ 的下層極小 0.2172 @ρ = 0.874；Fig. 6.7.5 畫的是 2 倍離心：0.5000 @ρ = 1）
習題：problem_61()（四張 E⊗E 表）、problem_62()（D3 的向量、張量、多項式；Hénon–Heiles）、problem_63()（X3）、problem_64()（X4）、
  problem_65()（三方破缺的偶極強度）、problem_66()（D6h：極向量 A2u + E1u、軸向量 A2g + E1g、X6 與苯的振動、IR／Raman 數）；d6h_group、d6h_chi
card_values() ＝ CARD_EXPECTED（plan_ch06 §6 卡片上與導讀裡的每個數字）；check()：每個結果兩種以上獨立算法（約 6 秒）；summary()。
直接執行：python utils/sym_prod.py（在 codes/ 裡）印一句成功訊息；--summary 印卡片數字。
"""
import itertools
import math
from fractions import Fraction
from functools import lru_cache

import numpy as np
import sympy as sp
from scipy.integrate import quad, solve_ivp
from scipy.linalg import eig_banded, eigh, null_space

try:
    from . import sym_mat as SM          # noqa: F401（給場景用：SP.SM）
    from . import sym_nonabel as SN
    from . import sym_oct as SO
    from . import sym_rot as RT
except ImportError:                      # 直接執行 python utils/sym_prod.py
    import sym_mat as SM                 # noqa: F401
    import sym_nonabel as SN
    import sym_oct as SO
    import sym_rot as RT

R_ = sp.Rational
SQ2, SQ3, SQ6 = sp.sqrt(2), sp.sqrt(3), sp.sqrt(6)
EPS3N = complex(np.exp(2j * np.pi / 3))  # ε = e^{2πi/3}（浮點）


def _close(a, b, tol: float = 1e-9) -> bool:
    return np.allclose(np.asarray(a, dtype=complex), np.asarray(b, dtype=complex), atol=tol, rtol=0)


def _f4(x) -> str:
    return f"{float(x):.4f}"


def _f2(x) -> str:
    return f"{float(x):.2f}"


def _num(M) -> np.ndarray:
    return np.array(sp.Matrix(M).evalf(), dtype=complex)


def exact(x, maxden: int = 96):
    """浮點（或複數）→ sympy：每個實、虛部寫成 ±√(p/q)（q ≤ maxden）；表格上的耦合係數都是這種數。"""
    z = complex(x)

    def one(v):
        if abs(v) < 1e-10:
            return sp.Integer(0)
        fr = Fraction(v * v).limit_denominator(maxden)
        if abs(float(fr) - v * v) > 1e-8:
            raise ValueError(f"{v} 不是 ±√(p/q)")
        s = sp.sqrt(sp.Rational(fr.numerator, fr.denominator))
        return s if v > 0 else -s

    return one(z.real) + sp.I * one(z.imag)


def exact_mat(M) -> sp.Matrix:
    M = np.asarray(M)
    return sp.Matrix(M.shape[0], M.shape[1], lambda i, j: exact(M[i, j]))


# ======================================================================== 有限群的一組不可約表示（數值）
class IrrepSet:
    """有限群 G 的一組不可約表示：elems（元素名）、D[α][g]（ℓ×ℓ complex）、names（irrep 的順序）。
    .dim(α) .chi(α)（依 elems 的特徵標陣列）.stack(α)（°G×ℓ×ℓ）.freq(χ) .decompose(χ) → {α: 頻率}"""

    def __init__(self, label: str, elems, D: dict, names=None):
        self.label = label
        self.elems = list(elems)
        self.names = list(names or D)
        self.D = {a: {g: np.asarray(D[a][g], dtype=complex) for g in self.elems} for a in self.names}
        self.order = len(self.elems)
        self._stack = {a: np.array([self.D[a][g] for g in self.elems]) for a in self.names}

    def dim(self, a) -> int:
        return self._stack[a].shape[1]

    def stack(self, a) -> np.ndarray:
        return self._stack[a]

    def chi(self, a) -> np.ndarray:
        return np.einsum("gii->g", self._stack[a])

    def freq(self, chi) -> dict:
        return {a: complex(np.vdot(self.chi(a), chi) / self.order) for a in self.names}

    def decompose(self, chi) -> dict:
        out = {}
        for a, f in self.freq(np.asarray(chi, dtype=complex)).items():
            if abs(f) > 1e-8:
                assert abs(f - round(f.real)) < 1e-8, (a, f)
                out[a] = int(round(f.real))
        return out

    def product_chi(self, a, b) -> np.ndarray:
        return self.chi(a) * self.chi(b)


def is_rep(G, D: dict) -> bool:
    """D[g]（數值）是否為 SN.FiniteGroup G 的表示（D(a)D(b) = D(ab)）。"""
    return all(_close(D[a] @ D[b], D[G.table[a, b]], 1e-9) for a in D for b in D)


def _coef(M, Q) -> float:
    """M Q Mᵀ = c Q 的 c（一維 irrep 由二次式或座標讀出）。"""
    R = M @ Q @ M.T
    return float(np.sum(R * Q) / np.sum(Q * Q))


_S_KOSHER = np.diag([1.0, -1.0, 1.0])
_S_TRIG = np.diag([1.0, -1.0, -1.0])


@lru_cache(maxsize=None)
def _dnum(alpha: str, g: str, basis: str) -> np.ndarray:
    if basis == "nonkosher":
        M = _num(SO.irrep(alpha, g, "tetragonal"))
        if alpha.replace("{", "").replace("}", "").startswith("T_2"):
            M = _S_KOSHER @ M @ _S_KOSHER
        return M
    if basis == "trigonal_book":                 # (6.1.36) 用的 T2 三方基底：E 那一對與 sym_oct 的 TRIG_U_T2 差一個號
        M = _num(SO.irrep(alpha, g, "trigonal"))
        if alpha.replace("{", "").replace("}", "").startswith("T_2"):
            M = _S_TRIG @ M @ _S_TRIG
        return M
    return _num(SO.irrep(alpha, g, basis))


@lru_cache(maxsize=None)
def oset(basis: str = "tetragonal") -> IrrepSet:
    """O 的五個 irrep（A_1 A_2 E T_1 T_2）：basis = "tetragonal"（O⊃D4⊃D2，kosher T2 (4.2.16)＝ sym_oct 的 TETRA）、
    "nonkosher"（(4.2.15)：𝒟^{T2} = χ^{A2}𝒟^{T1}，基底 yz, xz, xy）、"trigonal"（O⊃D3⊃C2）、"c4"（O⊃D4⊃C4 行進波）。"""
    D = {a: {g: _dnum(a, g, basis) for g in SO.O_NAMES} for a in SO.O_IRREP_NAMES}
    return IrrepSet(f"O/{basis}", SO.O_NAMES, D)


@lru_cache(maxsize=None)
def ohset(basis: str = "tetragonal") -> IrrepSet:
    """O_h 的十個 irrep（OH_IRREP_NAMES）；basis 同 oset。"""
    D = {a: {g: _dnum(a, g, basis) for g in SO.OH.names} for a in SO.OH_IRREP_NAMES}
    return IrrepSet(f"Oh/{basis}", SO.OH.names, D)


@lru_cache(maxsize=None)
def d4set() -> IrrepSet:
    """D4 = {1, R3², R3, R3³, R1², R2², i3, i4}（O ⊃ D4 ⊃ D2 的中間一環）：A1、B1（x²−y²）、A2（z）、B2（xy）、E（x, y）（書上 (6.1.21) 的順序）。"""
    D = {"A_1": {}, "B_1": {}, "A_2": {}, "B_2": {}, "E": {}}
    for g in SO.D4Z:
        M = SO.O.num[g].real
        D["A_1"][g] = [[1]]
        D["B_1"][g] = [[_coef(M, np.diag([1.0, -1, 0]))]]
        D["A_2"][g] = [[M[2, 2]]]
        D["B_2"][g] = [[_coef(M, np.array([[0, 1.0, 0], [1, 0, 0], [0, 0, 0]]))]]
        D["E"][g] = M[:2, :2]
    return IrrepSet("D4", SO.D4Z, D)


@lru_cache(maxsize=None)
def d2set() -> IrrepSet:
    """D2 = {1, R1², R2², R3²}：A1、A2（R3² 為 +1，像 z）、B1（R1² 為 +1，像 x）、B2（R2² 為 +1，像 y）（(4.2.10)）。"""
    D = {"A_1": {}, "A_2": {}, "B_1": {}, "B_2": {}}
    for g in SO.D2N:
        M = SO.O.num[g].real
        D["A_1"][g] = [[1]]
        D["A_2"][g] = [[M[2, 2]]]
        D["B_1"][g] = [[M[0, 0]]]
        D["B_2"][g] = [[M[1, 1]]]
    return IrrepSet("D2", SO.D2N, D)


@lru_cache(maxsize=None)
def d3set() -> IrrepSet:
    """D3 = {1, r1, r1², i2, i4, i5}（O ⊃ D3 ⊃ C2，C2 = {1, i4}）：A1、A2、E（三方基底 (v1, v2) 的塊，𝒟^E(r1) = rot(120°)）。"""
    tr = oset("trigonal")
    D = {"A_1": {}, "A_2": {}, "E": {}}
    for g in SO.D3T:
        M = tr.D["T_1"][g].real
        D["A_1"][g] = [[1]]
        D["A_2"][g] = [[M[2, 2]]]
        D["E"][g] = M[:2, :2]
    return IrrepSet("D3", SO.D3T, D)


@lru_cache(maxsize=None)
def c3set() -> IrrepSet:
    """C3 = {1, r1, r1²}：0_3、1_3、2_3，𝒟^{k_3}(r1) = ε^k、ε = e^{2πi/3}（第 2–3 章的約定）。"""
    p = {"1": 0, "r1": 1, "r1^2": 2}
    D = {f"{k}_3": {g: [[EPS3N ** (k * p[g])]] for g in SO.C3T} for k in range(3)}
    return IrrepSet("C3", SO.C3T, D)


@lru_cache(maxsize=None)
def tset() -> IrrepSet:
    """T = {1, r_k, r_k², R_k²}（12 個）：A、¹E、²E（ε、ε* 一對複共軛）、T（三維）。"""
    cls = {g: (0 if g == "1" or g.startswith("R") else (2 if g.endswith("^2") else 1)) for g in SO.T_NAMES}
    D = {"A": {}, "^1E": {}, "^2E": {}, "T": {}}
    for g in SO.T_NAMES:
        D["A"][g] = [[1]]
        D["^1E"][g] = [[EPS3N ** cls[g]]]
        D["^2E"][g] = [[EPS3N ** (2 * cls[g])]]
        D["T"][g] = SO.O.num[g].real
    return IrrepSet("T", SO.T_NAMES, D)


@lru_cache(maxsize=None)
def qset() -> IrrepSet:
    """四元數群 Q = {±1, ±i, ±j, ±k} 的二維 irrep（第三型的例子）。"""
    return IrrepSet("Q", SN.QUAT.names, {"E": {g: SN.QUAT.num[g] for g in SN.QUAT.names}})


@lru_cache(maxsize=None)
def c3vset() -> IrrepSet:
    """C3v（第 3 章的平面幾何 SN.C3V：σ₃ 沿 x 軸、σ₂ 在 60°、σ₁ 在 120°）：A1、A2（行列式）、E（2×2 本身＝ (3.3.7)）。"""
    D = {"A_1": {}, "A_2": {}, "E": {}}
    for g in SN.C3V.names:
        M = SN.C3V.num[g].real
        D["A_1"][g] = [[1]]
        D["A_2"][g] = [[round(np.linalg.det(M))]]
        D["E"][g] = M
    return IrrepSet("C3v", SN.C3V.names, D)


GROUP_OF = {"O": SO.O, "Oh": SO.OH, "D4": SO.O, "D2": SO.O, "D3": SO.O, "C3": SO.O, "T": SO.O, "Q": SN.QUAT, "C3v": SN.C3V}


# ======================================================================== §6.1A 兩粒子態與外積群（H601）
SITES = np.array(SO.SITES, dtype=float)       # |1⟩…|6⟩ = +z, −z, +x, −x, +y, −y
PAIRS = [(i, j) for i in range(1, 7) for j in range(1, 7)]      # (6.1.1)：|1,1⟩, |1,2⟩, …, |6,6⟩（字典序）


@lru_cache(maxsize=None)
def site_perm(g: str) -> np.ndarray:
    """O_h 元素 g 在六個頂點上的置換矩陣（⟨i|g|j⟩ = 1 若 g·site_j = site_i；＝ SO.site_perm 的數值版）。"""
    M = SO.OH.num[g].real
    P = np.zeros((6, 6))
    for s in range(6):
        P[int(np.argmin(np.linalg.norm(SITES - M @ SITES[s], axis=1))), s] = 1
    return P


def outer_op(a: str, b: str) -> np.ndarray:
    """(6.1.4)：(a, b)|k, l⟩ = a|k⟩ b|l⟩ → 36×36 = kron(P(a), P(b))（列、行都依 PAIRS）。"""
    return np.kron(site_perm(a), site_perm(b))


def kron_example() -> dict:
    """(6.1.7ₓ)：𝒟^{T1g}(r2) ⊗ 𝒟^{Eu}(i5)（四方基底；i5 是正常轉動，𝒟^{Eu}(i5) = 𝒟^E(i5)）。"""
    A = _dnum("T_{1g}", "r2", "tetragonal").real
    B = _dnum("E_u", "i5", "tetragonal").real
    return dict(T1=A, E=B, kron=np.kron(A, B))


KRON_617X_BOOK = dict(T1=np.array([[0, 0, -1], [1, 0, 0], [0, -1, 0]], float),
                      E=np.array([[-0.5, -math.sqrt(3) / 2], [-math.sqrt(3) / 2, 0.5]]))


def outer_frequencies() -> dict:
    """36 個兩粒子態在外積群 O_h^A × O_h^B（2304 個元素）下的分解：用 (6.1.9) 的 χ^{αβ}(a,b) = χ^α(a)χ^β(b)
    與 36×36 矩陣的跡逐一相乘求和 → {(α, β): 頻率}（＝ 單粒子 A1g ⊕ Eg ⊕ T1u 兩兩相乘，九個各一次）。"""
    S = ohset()
    tr6 = np.array([np.trace(site_perm(g)) for g in S.elems])
    out = {}
    for a in S.names:
        for b in S.names:
            f = np.einsum("a,b,a,b->", S.chi(a).conj(), S.chi(b).conj(), tr6, tr6) / S.order ** 2
            if abs(f) > 1e-9:
                out[(a, b)] = int(round(f.real))
    return out


def outer_trace_check(n: int = 60, seed: int = 6) -> bool:
    """(6.1.9) 的直接驗證：隨機 n 對 (a, b)，Tr kron(𝒟^α(a), 𝒟^β(b)) = χ^α(a)χ^β(b)，且 (a,b)(a',b') = (aa', bb')（(6.1.6)）。"""
    rng = np.random.default_rng(seed)
    S = ohset()
    names = SO.OH.names
    for _ in range(n):
        a, b, a2, b2 = (names[i] for i in rng.integers(0, 48, 4))
        al, be = (S.names[i] for i in rng.integers(0, 10, 2))
        if not math.isclose(np.trace(np.kron(S.D[al][a], S.D[be][b])).real,
                            (np.trace(S.D[al][a]) * np.trace(S.D[be][b])).real, abs_tol=1e-9):
            return False
        if not _close(outer_op(a, b) @ outer_op(a2, b2), outer_op(SO.OH.table[a, a2], SO.OH.table[b, b2])):
            return False
    return True


def six_hamiltonian(H: float = 0.0, S: float = 1.0, T: float = 0.0) -> np.ndarray:
    """(6.4.9)／Fig. 6.1.2 的六井 Hamilton 量：第一列 [H, −T, −S, −S, −S, −S]（對谷 −T、近鄰 −S）＝ SO.six_H(H, −T, −S)。"""
    return np.array(SO.six_H(H, -T, -S).evalf(), dtype=float)


def six_levels(H: float = 0.0, S: float = 1.0, T: float = 0.0) -> dict:
    """六井的三個能階（Fig. 6.1.2：T = 0 時 A1g = H − 4S、T1u = H、Eg = H + 2S）。"""
    return {"A_{1g}": H - T - 4 * S, "T_{1u}": H + T, "E_g": H - T + 2 * S}


SIX_STATES = {   # 六井的對稱本徵態（(4.3.20)–(4.3.27)；T1u 用實的 x, y, z；Eg 用 (2z²−x²−y², x²−y²) 的型）
    ("A_{1g}", 0): np.ones(6) / math.sqrt(6),
    ("T_{1u}", 0): np.array([0, 0, 1, -1, 0, 0]) / math.sqrt(2),
    ("T_{1u}", 1): np.array([0, 0, 0, 0, 1, -1]) / math.sqrt(2),
    ("T_{1u}", 2): np.array([1, -1, 0, 0, 0, 0]) / math.sqrt(2),
    ("E_g", 0): np.array([2, 2, -1, -1, -1, -1]) / math.sqrt(12),
    ("E_g", 1): np.array([0, 0, 1, 1, -1, -1]) / 2.0}

TWO_CARD = dict(H=0.0, SA=1.0, SB=0.25)          # 本系列取的數（書上只說 B 的間距較小）：照 Fig. 6.1.2(a) 的上下次序需 SB < SA/3


def two_particle_levels(SA: float = 1.0, SB: float = 0.25, H: float = 0.0) -> list:
    """(6.1.10)–(6.1.11)：不相互作用時 ε^{αβ} = ε^α_A + ε^β_B，簡併 ℓ^αℓ^β → [(E, α, β, 簡併)]（由高到低，Fig. 6.1.2(a)）。"""
    la, lb = six_levels(H, SA), six_levels(H, SB)
    dims = {"A_{1g}": 1, "T_{1u}": 3, "E_g": 2}
    out = [(la[a] + lb[b], a, b, dims[a] * dims[b]) for a in la for b in lb]
    return sorted(out, key=lambda t: -t[0])


def two_particle_numeric(SA: float = 1.0, SB: float = 0.25, H: float = 0.0) -> np.ndarray:
    """36×36 的 H_A ⊗ 1 + 1 ⊗ H_B 直接對角化（第二種算法）。"""
    HA, HB = six_hamiltonian(H, SA), six_hamiltonian(H, SB)
    return np.sort(np.linalg.eigvalsh(np.kron(HA, np.eye(6)) + np.kron(np.eye(6), HB)))[::-1]


# ======================================================================== §6.1B 相互作用、對角子群與耦合係數（H602）
def diag_frequencies(a: str, b: str, S: IrrepSet | None = None) -> dict:
    """(6.1.14)：f^γ(𝒟^{αβ}) = (1/°G)Σ χ^{γ*}χ^αχ^β（對角子群 O_h^{AB}）。"""
    S = S or ohset()
    return S.decompose(S.chi(a) * S.chi(b))


def char_product_rows() -> dict:
    """(6.1.14) 下面的特徵標乘積表（O_h 十個類 1 r R² R i I Ir IR² IR Ii）：A1g、Eg、T1g、T2g 與 χ^{T1u}χ^{T1u}。"""
    rows = {}
    for a in ("A_{1g}", "E_g", "T_{1g}", "T_{2g}", "T_{1u}"):
        r = SO.OH_IRREP_NAMES.index(a)
        rows[a] = [int(SO.OH_CHARS[r, c]) for c in range(10)]
    rows["T1u^2"] = [x * x for x in rows.pop("T_{1u}")]
    return rows


def cg(S: IrrepSet, a: str, b: str, gammas=None):
    """耦合係數（(6.1.16)–(6.1.18b) 的投影法）→ (C, cols)：C 是 (ℓ^aℓ^b)×(Σℓ^γ) 的么正矩陣，列 = 乘積基底 |a i⟩|b j⟩
    （字典序 i·ℓ^b + j，(6.1.8)），行 cols = [(γ, ω, m)]（m、ω 從 0 起；ω 是重複 irrep 的份數）。
    **相位約定（＝ 書上的）**：每一份 γ 取使 N^{γab}_{nkl} = ⟨kl|P^γ_{nn}|kl⟩ ≠ 0 的第一組 (n, k, l)（n 最外圈、字典序），
    C_{kl,n} = √N > 0，其餘行 C_{·,m} = P^γ_{mn}|kl⟩/√N。重複的 γ：下一份從下一組 (n, k, l) 出發，先對前幾份做 Gram–Schmidt。"""
    Da, Db = S.stack(a), S.stack(b)
    la, lb = Da.shape[1], Db.shape[1]
    M = np.einsum("gik,gjl->gijkl", Da, Db).reshape(S.order, la * lb, la * lb)
    freqs = S.decompose(np.einsum("gii->g", M))
    blocks, cols = [], []
    for c in (gammas or S.names):
        f = freqs.get(c, 0)
        if not f:
            continue
        Dc = S.stack(c)
        lc = Dc.shape[1]
        P = (lc / S.order) * np.einsum("gmn,gab->mnab", Dc.conj(), M)
        copies = []
        for n in range(lc):
            for kl in range(la * lb):
                if len(copies) == f:
                    break
                v = P[n, n][:, kl].copy()
                for cp in copies:
                    v = v - cp[:, n] * np.vdot(cp[:, n], v)
                nv = np.linalg.norm(v)
                if nv < 1e-7:
                    continue
                vn = v / nv
                copies.append(np.stack([P[m, n] @ vn for m in range(lc)], axis=1))
            if len(copies) == f:
                break
        for w, blk in enumerate(copies):
            blocks.append(blk)
            cols += [(c, w, m) for m in range(lc)]
    return np.hstack(blocks), cols


def cg_norms(S: IrrepSet, a: str, b: str, gamma: str) -> list:
    """(6.1.16c) 的歸一化常數 N^{γab}_{nkl}（依 (n, k, l) 字典序），= (ℓ^γ/°G)Σ 𝒟^{γ*}_{nn}𝒟^a_{kk}𝒟^b_{ll}。"""
    Da, Db, Dc = S.stack(a), S.stack(b), S.stack(gamma)
    lc = Dc.shape[1]
    out = []
    for n in range(lc):
        for k in range(Da.shape[1]):
            for l_ in range(Db.shape[1]):
                out.append(((n + 1, k + 1, l_ + 1),
                            float((lc / S.order) * np.sum(Dc[:, n, n].conj() * Da[:, k, k] * Db[:, l_, l_]).real)))
    return out


def cg_reduces(S: IrrepSet, a: str, b: str, C, cols) -> bool:
    """(6.1.19a)／(6.2.13)：對所有 g，C†(𝒟^a ⊗ 𝒟^b)(g)C = ⊕𝒟^γ(g)（依 cols 的塊）。"""
    for gi, g in enumerate(S.elems):
        L = C.conj().T @ np.kron(S.D[a][g], S.D[b][g]) @ C
        Rm = np.zeros_like(L)
        k = 0
        while k < len(cols):
            c = cols[k][0]
            lc = S.dim(c)
            Rm[k:k + lc, k:k + lc] = S.D[c][g]
            k += lc
        if not _close(L, Rm, 1e-9):
            return False
    return True


def factorization_lemma(S: IrrepSet, a: str, b: str, gamma: str) -> bool:
    """(6.2.15)：(ℓ^γ/°G)Σ_g 𝒟^a_{ik}𝒟^b_{jl}𝒟^{γ*}_{mn} = Σ_ω C^{abγ_ω}_{ijm} C^{abγ_ω *}_{kln}（對所有 i j k l m n）。"""
    Da, Db, Dc = S.stack(a), S.stack(b), S.stack(gamma)
    la, lb, lc = Da.shape[1], Db.shape[1], Dc.shape[1]
    lhs = (lc / S.order) * np.einsum("gik,gjl,gmn->ijmkln", Da, Db, Dc.conj())
    C, cols = cg(S, a, b, [gamma])
    rhs = np.zeros_like(lhs)
    for w in sorted({c[1] for c in cols}):
        idx = [k for k, c in enumerate(cols) if c[1] == w]
        blk = C[:, idx].reshape(la, lb, lc)
        rhs += np.einsum("ijm,kln->ijmkln", blk, blk.conj())
    return _close(lhs, rhs, 1e-9)


s2_, s3_, s6_ = 1 / math.sqrt(2), 1 / math.sqrt(3), 1 / math.sqrt(6)
CG_619C_BOOK = np.array([   # (6.1.19c)：T1u⊗T1u（行 A1g | Eg 1 2 | T1g 1 2 3 | T2g 1 2 3；非 kosher T2；列 11 12 13 21 … 33）
    [s3_, s6_, -s2_, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, s2_, 0, 0, s2_], [0, 0, 0, 0, -s2_, 0, 0, s2_, 0],
    [0, 0, 0, 0, 0, -s2_, 0, 0, s2_], [s3_, s6_, s2_, 0, 0, 0, 0, 0, 0], [0, 0, 0, s2_, 0, 0, s2_, 0, 0],
    [0, 0, 0, 0, s2_, 0, 0, s2_, 0], [0, 0, 0, -s2_, 0, 0, s2_, 0, 0], [s3_, -2 * s6_, 0, 0, 0, 0, 0, 0, 0]])
CG_630_BOOK = CG_619C_BOOK.copy()      # (6.1.30)：kosher 基底（O ⊃ D4 ⊃ D2）——只有 T2 的第 2 行（B2）變號
CG_630_BOOK[:, 7] *= -1
CG_623_BOOK = np.array([[s2_, s2_, 0, 0], [0, 0, s2_, s2_], [0, 0, -s2_, s2_], [s2_, -s2_, 0, 0]])   # (6.1.23)：D4 E⊗E（行 A1 B1 A2 B2）
CG_626_BOOK = {"A_2": np.array([[0, -1.0], [1, 0]]), "B_1": np.array([[1.0, 0], [0, -1]]), "B_2": np.array([[0, 1.0], [1, 0]])}  # (6.1.26)：E⊗A2、E⊗B1、E⊗B2
CG_635A_BOOK = np.array([[s2_, 0, s2_, 0], [0, s2_, 0, -s2_], [0, -s2_, 0, -s2_], [s2_, 0, -s2_, 0]])  # (6.1.35a)：D3 E⊗E（行 A1 A2 E_A E_B）
CG_635B_BOOK = np.array([[0, -1.0], [1, 0]])                                                         # (6.1.35b)：D3 E⊗A2
CG_636_BOOK = np.array([   # (6.1.36) 照印（表頭寫 T2⊗T1；行 A1 | E 1 2 | T2 1 2 3 | T1 1 2 3）
    [s3_, s6_, 0, s6_, s3_, 0, 0, 0, 0], [0, 0, -s6_, 0, 0, -s3_, 0, 0, s2_], [0, 0, -s3_, 0, 0, s6_, 0, -s2_, 0],
    [0, 0, -s6_, 0, 0, -s3_, 0, 0, -s2_], [s3_, -s6_, 0, s6_, -s3_, 0, 0, 0, 0], [0, s3_, 0, 0, -s6_, 0, s2_, 0, 0],
    [0, 0, -s3_, 0, 0, s6_, 0, s2_, 0], [0, s3_, 0, 0, -s6_, 0, -s2_, 0, 0], [s3_, 0, 0, -2 * s6_, 0, 0, 0, 0, 0]])
CG_617_N_BOOK = {(1, 1, 1): 0, (1, 1, 3): 0, (1, 2, 1): 0, (1, 2, 2): 0, (1, 2, 3): 0.5}   # (6.1.16) 之後：T1g 的 N_{nkl}


def cg_book(name: str):
    """依書上的鏈算出的表 → (C, cols)。name：
    "619c"（O_h 四方、非 kosher T2：T1u⊗T1u）、"630"（O 四方 kosher：T1⊗T1）、"623"（D4 E⊗E）、"626"（D4 E⊗A2／B1／B2 → dict）、
    "635a"（D3 E⊗E）、"635b"（D3 E⊗A2）、"636"（O 三方：T1⊗T1，行序照書 A1 E T2 T1）。"""
    if name == "619c":
        return cg(ohset("nonkosher"), "T_{1u}", "T_{1u}", ["A_{1g}", "E_g", "T_{1g}", "T_{2g}"])
    if name == "630":
        return cg(oset("tetragonal"), "T_1", "T_1", ["A_1", "E", "T_1", "T_2"])
    if name == "623":
        return cg(d4set(), "E", "E", ["A_1", "B_1", "A_2", "B_2"])
    if name == "626":
        return {x: cg(d4set(), "E", x, ["E"]) for x in ("A_2", "B_1", "B_2")}
    if name == "635a":
        return cg(d3set(), "E", "E", ["A_1", "A_2", "E"])
    if name == "635b":
        return cg(d3set(), "E", "A_2", ["E"])
    if name == "636":
        return cg(oset("trigonal_book"), "T_1", "T_1", ["A_1", "E", "T_2", "T_1"])
    raise KeyError(name)


def compare_table(C, cols, book) -> dict:
    """把算出的 C 與書上的表逐格比：{"equal": 全同, "flips": 需要整體變號的 γ 份（書上的相位約定不同）, "max": 最大差}。"""
    C = np.asarray(C).real
    book = np.asarray(book, float)
    flips = []
    Cm = C.copy()
    k = 0
    while k < len(cols):
        c, w = cols[k][0], cols[k][1]
        idx = [i for i, x in enumerate(cols) if x[0] == c and x[1] == w]
        if not np.allclose(Cm[:, idx], book[:, idx], atol=1e-9) and np.allclose(-Cm[:, idx], book[:, idx], atol=1e-9):
            Cm[:, idx] *= -1
            flips.append(c)
        k = max(idx) + 1
    return dict(equal=bool(np.allclose(C, book, atol=1e-9)), flips=flips, max=float(np.max(np.abs(Cm - book))))


# ---------------------------------------------------------------- 兩粒子相互作用：T1u⊗T1u 的分裂（Fig. 6.1.2(b)）
PAIR_V = dict(v0=1.0, v1=0.5, v2=0.25, K=0.1, W=0.05)
# 本系列取的相互作用（書上沒有數字，Fig. 6.1.2(b) 只是示意）：直接項 V(d)（同一井 d = 0、近鄰 √2、對谷 2）、
# 近鄰交換 K（|s,t⟩ ↔ |t,s⟩）、近鄰成對穿隧 W（|s,s⟩ ↔ |t,t⟩）。三種都只在對角子群 (a, a) 下不變。


def pair_V(v0: float = 1.0, v1: float = 0.5, v2: float = 0.25, K: float = 0.0, W: float = 0.0) -> np.ndarray:
    """兩粒子相互作用（36×36，列行依 PAIRS）：直接項（對角，依兩井距離 0、√2、2 取 v0、v1、v2）＋ K·近鄰交換 ＋ W·近鄰成對穿隧。"""
    n = len(PAIRS)
    Vm = np.zeros((n, n))
    idx = {p: k for k, p in enumerate(PAIRS)}
    for (i, j), kk in idx.items():
        r = np.linalg.norm(SITES[i - 1] - SITES[j - 1])
        Vm[kk, kk] = v0 if r < 0.1 else (v1 if r < 1.5 else v2)
        if 0.1 < r < 1.5:
            Vm[idx[(j, i)], kk] += K
        if i == j:
            for t in range(1, 7):
                if 0.1 < np.linalg.norm(SITES[i - 1] - SITES[t - 1]) < 1.5:
                    Vm[idx[(t, t)], kk] += W
    return Vm


def t1u_split(v=None) -> dict:
    """V 對 (T1u, T1u) 能階的一階分裂：用 (6.1.19c) 的耦合係數做出 |γ m⟩ = Σ C_{ij,m}|T1u i⟩|T1u j⟩，⟨γ m|V|γ m⟩ → {γ: 位移}。
    PAIR_V：A1g = (v0+v2)/2 + 2W、Eg = (v0+v2)/2 − W、T2g = v1 + K、T1g = v1 − K。"""
    C, cols = cg_book("619c")
    t = [SIX_STATES[("T_{1u}", k)] for k in range(3)]
    Vm = pair_V(**(v or PAIR_V))
    out = {}
    for k, (c, w, m) in enumerate(cols):
        st = sum(C[3 * i + j, k] * np.kron(t[i], t[j]) for i in range(3) for j in range(3))
        out.setdefault(c, []).append(float((st.conj() @ Vm @ st).real))
    for c in out:
        assert max(out[c]) - min(out[c]) < 1e-9, c
        out[c] = out[c][0]
    return out


def t1u_split_numeric(v=None, lam: float = 0.0) -> np.ndarray:
    """第二種算法：V 限制到 9 維 (T1u)⊗(T1u) 子空間後直接對角化（不用耦合係數）；lam > 0 時改為整個 36×36 的
    H_A ⊗ 1 + 1 ⊗ H_B + lam·V（TWO_CARD 的 SA、SB）對角化，取 0 附近的 9 個本徵值再除以 lam（一階）。"""
    Vm = pair_V(**(v or PAIR_V))
    if lam:
        HA, HB = six_hamiltonian(0, TWO_CARD["SA"]), six_hamiltonian(0, TWO_CARD["SB"])
        ev = np.linalg.eigvalsh(np.kron(HA, np.eye(6)) + np.kron(np.eye(6), HB) + lam * Vm)
        return np.sort(ev[np.abs(ev) < 0.25])[:9] / lam
    t = [SIX_STATES[("T_{1u}", k)] for k in range(3)]
    B = np.array([np.kron(t[i], t[j]) for i in range(3) for j in range(3)]).T
    return np.sort(np.linalg.eigvalsh(B.T @ Vm @ B))


def pair_V_invariant() -> tuple:
    """V 與所有 (a, a) 對易（48 個），但不與 (1, R1) 對易（(6.1.12)：⟨1,1|V|1,1⟩ ≠ ⟨1,3|V|1,3⟩）。"""
    Vm = pair_V(**PAIR_V)
    diag_ok = all(_close(outer_op(g, g) @ Vm, Vm @ outer_op(g, g)) for g in SO.OH.names)
    off = not _close(outer_op("1", "R1") @ Vm, Vm @ outer_op("1", "R1"))
    return diag_ok, off


# ======================================================================== §6.1C 子群鏈：D2 乘法表、(6.1.27)–(6.1.29) 的未知數（H603）
def d2_table() -> dict:
    """(6.1.20)：D2 一維 irreps 的乘法表 {(α, β): γ}（由特徵標相乘）。"""
    S = d2set()
    out = {}
    for a in S.names:
        for b in S.names:
            (g, f), = S.decompose(S.chi(a) * S.chi(b)).items()
            out[(a, b)] = g
    return out


D2_TABLE_BOOK = {("A_1", x): x for x in ("A_1", "A_2", "B_1", "B_2")}
D2_TABLE_BOOK.update({("A_2", "A_1"): "A_2", ("A_2", "A_2"): "A_1", ("A_2", "B_1"): "B_2", ("A_2", "B_2"): "B_1",
                      ("B_1", "A_1"): "B_1", ("B_1", "A_2"): "B_2", ("B_1", "B_1"): "A_1", ("B_1", "B_2"): "A_2",
                      ("B_2", "A_1"): "B_2", ("B_2", "A_2"): "B_1", ("B_2", "B_1"): "A_2", ("B_2", "B_2"): "A_1"})


def solve_6129(book_d: bool = False) -> dict:
    """(6.1.27b)–(6.1.29)：把 T1⊗T1 耦合矩陣的結構（未知數 a…k）代入 (𝒟^{T1}⊗𝒟^{T1})(R1)·C = C·(A1⊕E⊕T1⊕T2)(R1)（kosher）
    並解線性方程 → dict（a…k 的關係、歸一化後的值、與 (6.1.30) 是否相同）。book_d=True：再加上書上的 d = c/√2，看方程是否矛盾。"""
    a, c, f, i_ = sp.symbols("a c f i", positive=True)
    b, d, e, g, h, j_, k_ = sp.symbols("b d e g h j k")
    r2 = 1 / SQ2
    C = sp.Matrix([
        [a * r2, c * r2, e * r2, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, h * r2, 0, 0, k_ * r2], [0, 0, 0, 0, -f, 0, 0, -i_, 0],
        [0, 0, 0, 0, 0, -h * r2, 0, 0, k_ * r2], [a * r2, c * r2, -e * r2, 0, 0, 0, 0, 0, 0], [0, 0, 0, f, 0, 0, i_, 0, 0],
        [0, 0, 0, 0, -g, 0, 0, -j_, 0], [0, 0, 0, g, 0, 0, j_, 0, 0], [b, d, 0, 0, 0, 0, 0, 0, 0]])
    T1 = SO.irrep("T_1", "R1")
    blocks = sp.diag(sp.Matrix([[1]]), SO.irrep("E", "R1"), T1, SO.irrep("T_2", "R1"))
    eqs = list(sp.expand(sp.kronecker_product(T1, T1) * C - C * blocks))
    eqs = [q for q in eqs if q != 0]
    if book_d:
        eqs.append(d - c / SQ2)
    sol = sp.solve(eqs, [b, d, e, g, h, j_, k_], dict=True)
    if not sol:
        return dict(consistent=False)
    s = sol[0]
    out = dict(consistent=True, relations={str(x): sp.nsimplify(s[x]) for x in s})
    if book_d:
        out["forced"] = {str(x): s[x] for x in s}
        return out
    Cs = C.subs(s)
    vals = {}
    for col, var in ((0, a), (1, c), (3, f), (6, i_)):
        nrm = sp.sqrt(sp.simplify(sum(x ** 2 for x in Cs[:, col])))
        vals[str(var)] = sp.simplify(1 / (nrm / var))
    full = {a: vals["a"], c: vals["c"], f: vals["f"], i_: vals["i"]}
    Cn = Cs.subs(full)
    Cn = Cn.subs({e: sp.simplify(s[e].subs(full))}) if e in s else Cn
    num = np.array(Cn.evalf(), dtype=float)
    out["values"] = {str(x): sp.nsimplify(sp.simplify(v.subs(full))) for x, v in s.items()}
    out["values"].update({k: sp.nsimplify(v) for k, v in vals.items()})
    out["equals_630"] = bool(np.allclose(num, CG_630_BOOK, atol=1e-9))
    return out


# ======================================================================== §6.2 耦合係數的一般關係（H604）
def fs_indicator(S: IrrepSet, a: str) -> int:
    """Frobenius–Schur 指標 (1/°G)Σχ(g²)（用 𝒟(g)² 的跡）：+1 第一型（實）、0 第二型（複共軛成對）、−1 第三型（必複但與共軛等價）。"""
    St = S.stack(a)
    v = np.einsum("gij,gji->", St, St) / S.order
    return int(round(v.real))


def irrep_type(S: IrrepSet, a: str) -> int:
    """第二種算法（共軛檢查）：χ 不全為實 → 第二型；否則解不變雙線性型 𝒟(g)ᵀJ𝒟(g) = J，J 對稱 → 第一型、反對稱 → 第三型。"""
    if np.max(np.abs(S.chi(a).imag)) > 1e-9:
        return 2
    St = S.stack(a)
    ell = St.shape[1]
    A = np.vstack([np.kron(Dg.T, Dg.T) - np.eye(ell * ell) for Dg in St])
    ns = null_space(A)
    assert ns.shape[1] == 1, ns.shape
    J = ns[:, 0].reshape(ell, ell)
    return 1 if _close(J, J.T, 1e-8) else (3 if _close(J, -J.T, 1e-8) else 0)


def sym_squares(S: IrrepSet, a: str) -> tuple:
    """對稱化與反對稱化平方 (6.2.7)：χ_±(g) = (χ(g)² ± χ(g²))/2 → (對稱的分解, 反對稱的分解)。"""
    St = S.stack(a)
    c1 = np.einsum("gii->g", St)
    c2 = np.einsum("gij,gji->g", St, St)
    return S.decompose((c1 ** 2 + c2) / 2), S.decompose((c1 ** 2 - c2) / 2)


def cg_symmetry(S: IrrepSet, a: str) -> dict:
    """第二種算法：由 C^{aaγ} 的行看 C_{ij,m} = ±C_{ji,m}（(6.2.7a,b)）→ {γ: +1 對稱 / −1 反對稱}。"""
    C, cols = cg(S, a, a)
    ell = S.dim(a)
    perm = [j * ell + i for i in range(ell) for j in range(ell)]
    out = {}
    for k, (c, w, m) in enumerate(cols):
        sgn = 1 if _close(C[perm, k], C[:, k]) else (-1 if _close(C[perm, k], -C[:, k]) else 0)
        out.setdefault((c, w), set()).add(sgn)
    return {key: s.pop() for key, s in out.items() if len(s) == 1}


def scalar_cg_check(S: IrrepSet) -> bool:
    """(6.2.2)、(6.2.5)–(6.2.6)：C^{α α* 0}_{ij·} = δ_ij/√ℓ（第一型 α* = α；第二型取共軛的那個 irrep）。"""
    for a in S.names:
        conj = next(b for b in S.names if _close(S.chi(b), S.chi(a).conj()))
        if S.decompose(S.chi(a) * S.chi(conj)).get(S.names[0], 0) != 1:
            return False
        C, cols = cg(S, a, conj, [S.names[0]])
        ell = S.dim(a)
        if not _close(C[:, 0], np.eye(ell).ravel() / math.sqrt(ell)):
            return False
    return True


def unitarity_check(S: IrrepSet, a: str, b: str) -> bool:
    """(6.2.12a,b)：C†C = 1（正交歸一）、CC† = 1（完備）；(6.2.14) 的反關係 𝒟^a⊗𝒟^b = C(⊕𝒟^γ)C†。"""
    C, cols = cg(S, a, b)
    n = C.shape[0]
    return _close(C.conj().T @ C, np.eye(n)) and _close(C @ C.conj().T, np.eye(n)) and cg_reduces(S, a, b, C, cols)


# ======================================================================== §6.3A–C 單位向量、張量、多項式（H605）
TENSOR_COLS = ["A_{1g}", "E_g", "E_g", "T_{1g}", "T_{1g}", "T_{1g}", "T_{2g}", "T_{2g}", "T_{2g}"]
TENSOR_KEYS = [("A_{1g}", 0), ("E_g", 0), ("E_g", 1), ("T_{1g}", 0), ("T_{1g}", 1), ("T_{1g}", 2),
               ("T_{2g}", 0), ("T_{2g}", 1), ("T_{2g}", 2)]


def unit_tensors() -> dict:
    """(6.3.6)：𝐓^γ_m = Σ C^{T1uT1uγ}_{ij,m} x̂_i x̂_j（(6.1.19c) 的書上相位，非 kosher T2）→ {(γ, m): 3×3}。"""
    C, cols = cg_book("619c")
    return {(c, m): C[:, k].real.reshape(3, 3) for k, (c, w, m) in enumerate(cols)}


def tensor_components(Mt) -> dict:
    """(6.3.13a)：𝒯^γ_m = Σ C_{ij,m}𝒯_{ij}（任何 3×3 張量 → 九個對稱化分量）。"""
    Mt = np.asarray(Mt, float)
    return {key: float(np.sum(Tm * Mt)) for key, Tm in unit_tensors().items()}


def tensor_from_components(comp: dict) -> np.ndarray:
    """(6.3.13b)：𝒯_{ij} = Σ C_{ij,m}𝒯^γ_m。"""
    U = unit_tensors()
    return sum(v * U[k] for k, v in comp.items())


def tensor_transform_check() -> bool:
    """(6.3.5)：𝒱(g)𝐓^γ_m𝒱(g)ᵀ = Σ_{m'}𝒟^γ_{m'm}(g)𝐓^γ_{m'}（48 個 g；非 kosher 的 O_h irreps）。"""
    U = unit_tensors()
    S = ohset("nonkosher")
    for g in SO.OH.names:
        V = SO.OH.num[g].real
        for (c, m), Tm in U.items():
            rhs = sum(S.D[c][g][mp, m].real * U[(c, mp)] for mp in range(S.dim(c)))
            if not _close(V @ Tm @ V.T, rhs):
                return False
    return True


XS = sp.symbols("x y z", real=True)
POLY_639 = {   # (6.3.9)：O_h 的二次對稱化多項式（T1g 為零）
    ("A_{1g}", 0): (XS[0] ** 2 + XS[1] ** 2 + XS[2] ** 2) / SQ3, ("E_g", 0): (XS[0] ** 2 + XS[1] ** 2 - 2 * XS[2] ** 2) / SQ6,
    ("E_g", 1): (-XS[0] ** 2 + XS[1] ** 2) / SQ2, ("T_{2g}", 0): SQ2 * XS[1] * XS[2], ("T_{2g}", 1): SQ2 * XS[0] * XS[2],
    ("T_{2g}", 2): SQ2 * XS[0] * XS[1]}


def poly_from_tensors() -> dict:
    """把 (6.3.6) 的 x̂_i 換成座標 x_i：p^γ_m = Σ C_{ij,m}x_ix_j（(6.3.7)–(6.3.9)；T1g 的反對稱積為零）。"""
    out = {}
    for key, Tm in unit_tensors().items():
        p = sum(exact(Tm[i, j]) * XS[i] * XS[j] for i in range(3) for j in range(3))
        out[key] = sp.expand(p)
    return out


def poly_transform_check() -> bool:
    """(6.3.8)：座標跟基向量一樣變：p_m(𝒱ᵀx) = Σ_{m'}𝒟_{m'm}(g)p_{m'}(x)（在隨機點上比，48 個 g）。"""
    S = ohset("nonkosher")
    rng = np.random.default_rng(3)
    pts = rng.normal(size=(5, 3))
    f = {k: sp.lambdify(XS, v, "numpy") for k, v in POLY_639.items()}
    for g in SO.OH.names:
        V = SO.OH.num[g].real
        for (c, m) in POLY_639:
            for x in pts:
                lhs = f[(c, m)](*(V.T @ x))
                rhs = sum(S.D[c][g][mp, m].real * f[(c, mp)](*x) for mp in range(S.dim(c)))
                if abs(lhs - rhs) > 1e-9:
                    return False
    return True


def o3_labels() -> dict:
    """1 ⊗ 1 = 0 ⊕ 1 ⊕ 2（O3）↓ O_h：0⁺ → A1g、1⁺ → T1g、2⁺ → Eg ⊕ T2g（RT.o3_to_o）；並用 χ¹(ω)² = χ⁰ + χ¹ + χ² 驗證。"""
    out = {}
    for l_, name in ((0, "0^+"), (1, "1^+"), (2, "2^+")):
        f = RT.o3_to_o(l_)
        out[name] = {a: n for a, n in zip(["A_{1g}", "A_{2g}", "E_g", "T_{1g}", "T_{2g}"], f) if n}
    ok = all(abs(RT.chi(1, w) ** 2 - (RT.chi(0, w) + RT.chi(1, w) + RT.chi(2, w))) < 1e-9 for w in (0, 37, 90, 120, 180))
    out["ok"] = ok
    return out


# ======================================================================== §6.3D 應力與應變（H606）
def stress_examples() -> dict:
    """Fig. 6.3.2(b)／(6.3.33) 的例子 → 每個例子的九個對稱化分量（只列非零的）：
    "pressure"（等向張力 p·1）、"uniaxial"（純張力 (F/A)ẑẑ，F/A = 1）、"shear12"（𝒯_12 = 𝒯_21 = 1）、"torque"（𝒯_23 = −𝒯_32 = 1）。"""
    ex = {"pressure": np.eye(3), "uniaxial": np.diag([0.0, 0, 1]), "shear12": np.array([[0, 1.0, 0], [1, 0, 0], [0, 0, 0]]),
          "torque": np.array([[0, 0, 0], [0, 0, 1.0], [0, -1, 0]])}
    return {k: {key: v for key, v in tensor_components(M).items() if abs(v) > 1e-12} for k, M in ex.items()}


def face_forces(Tm) -> list:
    """(6.3.12)：𝐅(j) = 𝒯·x̂_j（立方體第 j 面（+x_j 側）每單位面積受的力；相對的面是 −𝐅(j)）。"""
    Tm = np.asarray(Tm, float)
    return [Tm[:, j].copy() for j in range(3)]


def strain_of_displacement(s_func, h: float = 1e-6) -> np.ndarray:
    """(6.3.14a)：𝒮_{ij} = ∂s_i/∂x_j（數值微分；均勻應變時與點無關）。"""
    x0 = np.array([0.3, -0.2, 0.1])
    Sm = np.zeros((3, 3))
    for j in range(3):
        e = np.zeros(3)
        e[j] = h
        Sm[:, j] = (np.asarray(s_func(x0 + e)) - np.asarray(s_func(x0 - e))) / (2 * h)
    return Sm


# ======================================================================== §6.3E 彈性常數（H607）
def _mandel_basis() -> list:
    """對稱 3×3 矩陣的正交歸一（Frobenius）基底：e11, e22, e33, (e23+e32)/√2, (e13+e31)/√2, (e12+e21)/√2。"""
    out = []
    for i in range(3):
        m = np.zeros((3, 3))
        m[i, i] = 1
        out.append(m)
    for i, j in ((1, 2), (0, 2), (0, 1)):
        m = np.zeros((3, 3))
        m[i, j] = m[j, i] = 1 / math.sqrt(2)
        out.append(m)
    return out


MANDEL = _mandel_basis()


def sym2_rep(Rm) -> np.ndarray:
    """對稱二階張量（應力、應變）上的 6×6 表示：ρ(g)_{ab} = ⟨e_a, 𝒱 e_b 𝒱ᵀ⟩（Mandel 基底）。"""
    Rm = np.asarray(Rm, float)
    return np.array([[np.sum(ea * (Rm @ eb @ Rm.T)) for eb in MANDEL] for ea in MANDEL])


def _invariant_dim(mats, space: str) -> int:
    """在 space 上不變的張量個數（暴力法：對群的生成元解 X = g·X 的線性方程）。
    space："81"（k_{hijk} 不加任何對稱）、"36"（6×6 任意 K：應力、應變對稱）、"21"（6×6 對稱 K：再加互易）。"""
    rows = []
    for Rm in mats:
        if space == "81":
            G4 = np.kron(np.kron(Rm, Rm), np.kron(Rm, Rm))
            rows.append(G4 - np.eye(81))
        else:
            r = sym2_rep(Rm)
            rows.append(np.kron(r, r) - np.eye(36))
    A = np.vstack(rows)
    if space == "21":
        B = []
        for a in range(6):
            for b in range(a, 6):
                v = np.zeros((6, 6))
                v[a, b] = v[b, a] = 1
                B.append(v.ravel())
        B = np.array(B).T
        A = A @ B
    return int(A.shape[1] - np.linalg.matrix_rank(A, tol=1e-9))


def _gens(name: str) -> list:
    num = SO.OH.num
    if name == "Oh":
        return [num["R1"].real, num["r1"].real, num["I"].real]
    if name == "iso":                         # O_h 再加繞 z 的 45°：生成 O3 的稠密子群（等向）
        return _gens("Oh") + [SN_rot3((0, 0, 1), 45)]
    if name == "D4h":
        return [num["R3"].real, num["R1^2"].real, num["I"].real]
    if name == "C3v":                         # 三方軸 (111) 上的 C3v：r1 與鏡面 I·i4
        return [num["r1"].real, num["I.i4"].real]
    if name == "none":
        return [np.eye(3)]
    raise KeyError(name)


def SN_rot3(axis, deg) -> np.ndarray:
    return np.array(SN.rot3(axis, deg).evalf(), dtype=float)


def elastic_counts() -> dict:
    """彈性常數的個數（暴力法）：沒有對稱 81 → 36 → 21；O_h：81 個分量裡 4 個不變量、36 → 3、21 → 3；
    等向（O_h ＋ 45°）：2；D4h：7 → 6（互易）；C3v：8 → 6。"""
    out = {"none": {s: _invariant_dim(_gens("none"), s) for s in ("81", "36", "21")}}
    for name in ("Oh", "iso", "D4h", "C3v"):
        out[name] = {s: _invariant_dim(_gens(name), s) for s in ("81", "36", "21")}
    return out


def elastic_counts_chars() -> dict:
    """第二種算法（張量定理）：Sym²(向量) 的頻率 f^γ → 36 個常數變 Σ(f^γ)²、21 個（互易）變 Σf^γ(f^γ+1)/2。"""
    out = {}
    for name, S, vec in (("Oh", ohset(), "T_{1u}"),):
        sym, _ = sym_squares(S, vec)
        out[name] = (sum(f * f for f in sym.values()), sum(f * (f + 1) // 2 for f in sym.values()), sym)
    # D4h、C3v：用三維矩陣的特徵標（不需要 irreps 的矩陣）
    for name, elems in (("D4h", SO.D4H), ("C3v", SO.C3V)):
        mats = [SO.OH.num[g].real for g in elems]
        chi_s = np.array([(np.trace(M) ** 2 + np.trace(M @ M)) / 2 for M in mats])
        n36 = int(round(np.sum(chi_s * chi_s) / len(mats)))
        chi_sym2 = np.array([np.trace(np.kron(sym2_rep(M), sym2_rep(M)) @ _sym_proj36()) for M in mats])
        n21 = int(round(np.sum(chi_sym2) / len(mats)))
        out[name] = (n36, n21)
    return out


@lru_cache(maxsize=None)
def _sym_proj36() -> np.ndarray:
    P = np.zeros((36, 36))
    for a in range(6):
        for b in range(6):
            P[a * 6 + b, a * 6 + b] += 0.5
            P[a * 6 + b, b * 6 + a] += 0.5
    return P


def elastic_tensor(kA: float, kE: float, kT2: float, kT1: float = 0.0) -> np.ndarray:
    """(6.3.26)：k_{hi,jk} = Σ_γΣ_m C_{hi,m}C_{jk,m}k^γ（(6.1.19c) 的係數）→ 3×3×3×3。"""
    C, cols = cg_book("619c")
    kk = {"A_{1g}": kA, "E_g": kE, "T_{2g}": kT2, "T_{1g}": kT1}
    K = sum(kk[c] * np.outer(C[:, n].real, C[:, n].real) for n, (c, w, m) in enumerate(cols))
    return K.reshape(3, 3, 3, 3)


def elastic_from_cartesian(K4) -> dict:
    """(6.3.28)（更正的 (6.3.28)ₓ）：k^γ = Σ C_{hi,m}C_{jk,m}k_{hi,jk}（取 m = 第一個）→ {A1, E, T2}；並回傳 Cartesian 的三個數。"""
    C, cols = cg_book("619c")
    K = np.asarray(K4).reshape(9, 9)
    first = {}
    for n, (c, w, m) in enumerate(cols):
        if m == 0:
            first[c] = float(C[:, n].real @ K @ C[:, n].real)
    K4 = np.asarray(K4)
    return dict(kA=first["A_{1g}"], kE=first["E_g"], kT2=first["T_{2g}"],
                k1111=float(K4[0, 0, 0, 0]), k1122=float(K4[0, 0, 1, 1]), k1212=float(K4[0, 1, 0, 1]),
                kA_x=float(K4[0, 0, 0, 0] + 2 * K4[0, 0, 1, 1]), kE_x=float(K4[0, 0, 0, 0] - K4[0, 0, 1, 1]), kT2_x=float(2 * K4[0, 1, 0, 1]))


def moduli(kA, kE, kT2) -> dict:
    """(6.3.30)、(6.3.32)、(6.3.36)、(6.3.43)（精確值可給 Fraction）：B = k^{A1}/3、μ = k^{T2}/2、Y = 3k^Ek^{A1}/(2k^{A1}+k^E)、
    σ = (k^{A1}−k^E)/(2k^{A1}+k^E)、A = (k^E − k^{T2})/(k^{T2}/2)；Cartesian k_{11,11} = k^{A1}/3 + 2k^E/3、k_{11,22} = (k^{A1}−k^E)/3、k_{12,12} = k^{T2}/2。"""
    kA, kE, kT2 = (Fraction(x) if isinstance(x, int) else x for x in (kA, kE, kT2))
    return dict(B=kA / 3, mu=kT2 / 2, Y=3 * kE * kA / (2 * kA + kE), sigma=(kA - kE) / (2 * kA + kE),
                A=(kE - kT2) / (kT2 / 2), k1111=kA / 3 + 2 * kE / 3, k1122=(kA - kE) / 3, k1212=kT2 / 2)


def moduli_numeric(kA: float, kE: float, kT2: float) -> dict:
    """第二種算法：用 Cartesian 張量直接「做實驗」——等向應變 ε·1 → B = 𝒯_11/(ΔV/V)；剪應變 𝒮_12 = 𝒮_21 = ε → μ = 𝒯_12/(2ε)；
    單軸應力 σẑẑ → 解 𝒮 = K⁻¹𝒯：Y = σ/𝒮_33、Poisson σ = −𝒮_11/𝒮_33。"""
    K4 = elastic_tensor(kA, kE, kT2)
    Km = np.array([[np.einsum("ij,ijkl,kl->", ea, K4, eb) for eb in MANDEL] for ea in MANDEL])

    def stress(Sm):
        return np.einsum("ijkl,kl->ij", K4, Sm)

    eps = 1e-3
    B = stress(eps * np.eye(3))[0, 0] / (3 * eps)
    mu = stress(eps * np.array([[0, 1.0, 0], [1, 0, 0], [0, 0, 0]]))[0, 1] / (2 * eps)
    t = np.zeros(6)
    t[2] = 1.0
    s = np.linalg.solve(Km, t)
    return dict(B=float(B), mu=float(mu), Y=float(1 / s[2]), sigma=float(-s[0] / s[2]))


MODULI_CARD = {"spring111": (1, 1, 1), "spring211": (2, 1, 1)}      # (j, k, l)，單位 1/4a


def spring_model(j, k, l_) -> tuple:
    """(6.3.41a–c)：SrTiO₃ 的剛性離子彈簧模型（單位 1/4a）：k^{A1} = j + 4k + 4l、k^E = j + k + l、k^{T2} = 2k + 2l。"""
    return (j + 4 * k + 4 * l_, j + k + l_, 2 * k + 2 * l_)


def spring_model_bonds(j: float, k: float, l_: float, a: float = 1.0) -> dict:
    """第二種算法（暴力法）：鈣鈦礦晶格（Ti 在原點、O 在 a·x̂_i、Sr 在 (a,a,a)，晶格常數 2a）每根中心力彈簧的能量 ½κ(d̂·𝒮·d)²
    加起來 → k_{hijk} = (1/V)Σκ d_hd_id_jd_k/|d|²（V = 8a³；每個原子都在反演中心，沒有內部弛豫）→ 乘 4a 後回傳 k^γ 與 Cartesian 常數。"""
    atoms = [("Ti", np.zeros(3))] + [("O", a * np.eye(3)[i]) for i in range(3)] + [("Sr", a * np.ones(3))]
    kap = {frozenset(["Ti", "O"]): (j, a), frozenset(["O"]): (k, a * math.sqrt(2)), frozenset(["Sr", "O"]): (l_, a * math.sqrt(2))}
    K4 = np.zeros((3, 3, 3, 3))
    shifts = [2 * a * np.array(s) for s in itertools.product((-1, 0, 1), repeat=3)]
    for t1, p1 in atoms:
        for t2, p2 in atoms:
            key = frozenset([t1, t2])
            if key not in kap:
                continue
            kk, dist = kap[key]
            for sh in shifts:
                dvec = p2 + sh - p1
                if abs(np.linalg.norm(dvec) - dist) < 1e-9:
                    K4 += 0.5 * kk * np.einsum("h,i,j,k->hijk", dvec, dvec, dvec, dvec) / dist ** 2
    K4 /= 8 * a ** 3
    r = elastic_from_cartesian(K4 * 4 * a)
    return r


CUBIC_KITTEL = {   # 室溫的 C11、C12、C44（10^11 N/m²；Kittel《Introduction to Solid State Physics》表 4.1、4.2；本系列核對 Fig. 6.3.9 用）
    "diamond": (10.76, 1.25, 5.76), "Si": (1.66, 0.64, 0.80), "Ge": (1.29, 0.48, 0.67), "NaCl": (0.487, 0.124, 0.126),
    "Au": (1.923, 1.631, 0.420), "Pb": (0.495, 0.423, 0.149), "Na": (0.073, 0.062, 0.042), "Li": (0.135, 0.114, 0.088),
    "W": (5.233, 2.045, 1.607)}
FIG639_A_BOOK = {"diamond": -0.3, "Si": -0.7, "Ge": -0.8, "NaCl": 0.9, "Au": -1.3, "Pb": -1.5, "Na": -1.7, "Li": -1.8}


def anisotropy(c11, c12, c44) -> float:
    """(6.3.43)：A = (k_{11,11} − k_{11,22} − 2k_{12,12})/k_{12,12}（k_{11,11} = C11、k_{11,22} = C12、k_{12,12} = C44）。"""
    return (c11 - c12 - 2 * c44) / c44


# ======================================================================== §6.3F 對稱化的電場、磁場（H608）
def c3v_field(g: str, kind: str = "E") -> np.ndarray:
    """C3v 操作對場的 3×3 作用（SN.C3V 的平面幾何＋z 軸）：極向量 E → 𝒱(g)；軸向量 B（電流環）→ det(𝒱)·𝒱(g)。"""
    M = np.array(SN.block3(SN.C3V.mats[g]).evalf(), dtype=float)
    return M if kind == "E" else np.linalg.det(M) * M


def c3v_field_irreps() -> dict:
    """E → A1 ⊕ E（(6.3.46)）、B → A2 ⊕ E（(6.3.48)）：用 3×3 的跡分解。"""
    S = c3vset()
    return {k: S.decompose(np.array([np.trace(c3v_field(g, k)) for g in S.elems])) for k in ("E", "B")}


def p3_12(kind: str = "B", vec=(1, 0, 0)) -> np.ndarray:
    """(6.3.47)：P^3_{12} = (2/6)Σ𝒟^E_{12}(g)* g 作用在場上；P^3_{12}B̂₁ = −B̂₂、P^3_{22}B̂₁ = B̂₁（更正的 √3/2）。"""
    S = c3vset()
    v = np.array(vec, float)
    return sum((2 / 6) * S.D["E"][g][0, 1].conjugate().real * (c3v_field(g, kind) @ v) for g in S.elems)


def p3_12_book(s=sp.Rational(3, 2)):
    """照印 (6.3.47)（把每個 √3/2 印成 3/2）逐項算 → B̂₂ 的係數（s = 3/2 → −3；s = √3/2 → −1）與 B̂₁ 的係數。"""
    terms = [(-s, (-R_(1, 2), s)), (s, (-R_(1, 2), -s)), (-s, (R_(1, 2), s)), (s, (R_(1, 2), -s))]
    b1 = sp.Rational(1, 3) * sum(c * t[0] for c, t in terms)
    b2 = sp.Rational(1, 3) * sum(c * t[1] for c, t in terms)
    return sp.simplify(b1), sp.simplify(b2)


def fig6311() -> dict:
    """Fig. 6.3.11 的四個式子（更正後）：rE₁ = −½E₁ + (√3/2)E₂、σ₃E₁ = E₁、σ₃E₂ = −E₂；rB₁ = −½B₁ + (√3/2)B₂、σ₃B₁ = −B₁、σ₃B₂ = B₂ …"""
    out = {}
    for kind in ("E", "B"):
        for g in ("r", "s3"):
            M = c3v_field(g, kind)
            out[(kind, g)] = (M[:, 0].round(12), M[:, 1].round(12), M[:, 2].round(12))
    return out


def striction_counts() -> dict:
    """(6.3.50)–(6.3.52) 的線性電致／磁致伸縮常數個數（暴力法：對稱 ij 的三階張量 ε_{ij,k} 在群下不變；軸向量多乘 det）：
    C3v：電 4、磁 2；O_h：0、0。"""
    out = {}
    for name in ("C3v", "Oh"):
        for kind in ("E", "B"):
            rows = []
            for Rm in _gens(name):
                det = np.linalg.det(Rm)
                Mk = Rm if kind == "E" else det * Rm
                rows.append(np.kron(sym2_rep(Rm), Mk) - np.eye(18))
            A = np.vstack(rows)
            out[(name, kind)] = int(18 - np.linalg.matrix_rank(A, tol=1e-9))
    return out


def striction_counts_chars() -> dict:
    """第二種算法：Hom(場, 應變) 的維數 Σ_γ f^γ(Sym²V) f^γ(場)（C3v：Sym²V = 2A1 + 2E）。"""
    S = c3vset()
    sq = np.array([(np.trace(c3v_field(g, "E")) ** 2 + np.trace(c3v_field(g, "E") @ c3v_field(g, "E"))) / 2 for g in S.elems])
    fS = S.decompose(sq)
    out = {"Sym2": fS}
    for kind in ("E", "B"):
        fF = S.decompose(np.array([np.trace(c3v_field(g, kind)) for g in S.elems]))
        out[kind] = sum(fS.get(a, 0) * fF.get(a, 0) for a in S.names)
    return out


# ======================================================================== §6.4 張量算子、參數定理、Wigner–Eckart（H609）
def six_invariants() -> dict:
    """(6.4.8)：V^{(A1g)_1} = |A1g⟩⟨A1g|、V^{(A1g)_2} = Σ|T1u j⟩⟨T1u j|/√3、V^{(A1g)_3} = Σ|Eg j⟩⟨Eg j|/√2（6×6）。"""
    P = {"A": np.outer(SIX_STATES[("A_{1g}", 0)], SIX_STATES[("A_{1g}", 0)]),
         "T": sum(np.outer(SIX_STATES[("T_{1u}", k)], SIX_STATES[("T_{1u}", k)]) for k in range(3)) / math.sqrt(3),
         "E": sum(np.outer(SIX_STATES[("E_g", k)], SIX_STATES[("E_g", k)]) for k in range(2)) / math.sqrt(2)}
    return P


SIX_INV_ROWS_BOOK = {"A": [1 / 6] * 6, "T": [0.5, 0.5, 0, 0, 0, 0], "E": [1 / 3, 1 / 3, 1 / 6, 1 / 6, 1 / 6, 1 / 6]}   # 照印（少了負號；T、E 還要再除 √3、√2）
SIX_CARD = dict(H=0.0, S=1.0, T=0.5)       # 本系列取的數（H609：由三個能階反推三個參數）


def six_from_levels(eA, eT, eE) -> dict:
    """(6.4.10)：H = ε^{A1g}V₁ + √3ε^{T1u}V₂ + √2ε^{Eg}V₃ 的第一列 ⇒ H = (εA + 3εT + 2εE)/6、T = −(εA − 3εT + 2εE)/6、S = −(εA − εE)/6。"""
    return dict(H=(eA + 3 * eT + 2 * eE) / 6, T=-(eA - 3 * eT + 2 * eE) / 6, S=-(eA - eE) / 6)


def param_count(freqs: dict) -> int:
    """參數定理 (6.4.11a)：n = Σ f^α(f^α + 1)/2（Hermitian）。"""
    return sum(f * (f + 1) // 2 for f in freqs.values())


def commutant_count(mats, hermitian: bool = True) -> int:
    """第二種算法：與所有 Γ(g) 對易的（實對稱）矩陣個數 ＝ 解 Γ(g)XΓ(g)ᵀ = X 的線性方程（暴力法）。"""
    n = mats[0].shape[0]
    rows = [np.kron(M, M) - np.eye(n * n) for M in mats]
    A = np.vstack(rows)
    if hermitian:
        B = []
        for a in range(n):
            for b in range(a, n):
                v = np.zeros((n, n))
                v[a, b] = v[b, a] = 1
                B.append(v.ravel())
        A = A @ np.array(B).T
    return int(A.shape[1] - np.linalg.matrix_rank(A, tol=1e-8))


def xf6_param_counts() -> dict:
    """SF6 型 21 個座標：全部 Σf(f+1)/2 = 11（T1u 三次）；去掉平移、轉動後的 15 維振動空間 = 7（書 p. 496）；兩法。"""
    gens = ["R1", "r1", "I"]
    mats = [SO.xf6_rep(g) for g in gens]
    full_chars = SO.xf6_classify()["all"]
    f = dict(zip(SO.OH_IRREP_NAMES, full_chars))
    f = {a: n for a, n in f.items() if n}
    fv = dict(f)
    fv["T_{1u}"] -= 1
    fv["T_{1g}"] -= 1
    fv = {a: n for a, n in fv.items() if n}
    # 振動子空間：與平移、轉動正交
    tr = [SO.xf6_translation(i) for i in range(3)]
    rot = []
    for ax in range(3):
        v = np.zeros(21)
        for c, (_, at, d) in enumerate(SO.XF6_COORDS_LIST):
            pos = np.zeros(3) if at == 6 else np.array(SO.F_SITES[at], float)
            v[c] = np.cross(np.eye(3)[ax], pos) @ np.array(d, float)
        rot.append(v)
    Q = null_space(np.array(tr + rot))
    vib = [Q.T @ M @ Q for M in mats]
    return dict(full=param_count(f), full_numeric=commutant_count(mats), vib=param_count(fv), vib_numeric=commutant_count(vib),
                freqs=f, vib_freqs=fv)


def six_dipole_elements() -> dict:
    """六井的電偶極矩陣元 ⟨γ k|x_i|β j⟩（位置算子對角，井在單位距離、q = 1）：{(γ, k, i, β, j): 值}，只列非零。"""
    X = [np.diag(SITES[:, i]) for i in range(3)]
    out = {}
    for (g1, k), v1 in SIX_STATES.items():
        for (g2, j), v2 in SIX_STATES.items():
            for i in range(3):
                val = float(v1 @ X[i] @ v2)
                if abs(val) > 1e-12:
                    out[(g1, k, i, g2, j)] = val
    return out


def _we_fit(S: IrrepSet, alpha: str, beta: str, gamma: str, Mel) -> dict:
    """Mel[i, j, k] = ⟨γ k|V^α_i|β j⟩ → 用投影法的 C^{αβγ}_{ij,k} 擬合：R = (1/ℓ^γ)ΣC·⟨⟩（更正的 (6.4.15b)）與照印的 (1/ℓ^α)ΣC·⟨⟩；
    ok = (⟨⟩ = C*·R)、ok_noconj = (⟨⟩ = C·R)、ok_book = (⟨⟩ = C*·R_book)。"""
    C, cols = cg(S, alpha, beta, [gamma])
    la, lb, lg = S.dim(alpha), S.dim(beta), S.dim(gamma)
    Cb = C.reshape(la, lb, lg)
    tot = np.sum(Cb * Mel)
    R, Rb = tot / lg, tot / la
    return dict(R=complex(R), R_book=complex(Rb), ok=_close(Mel, Cb.conj() * R), ok_noconj=_close(Mel, Cb * R),
                ok_book=_close(Mel, Cb.conj() * Rb))


def _six_mel(gamma_states, beta_states, ops) -> np.ndarray:
    return np.array([[[np.vdot(gs, V @ bs) for gs in gamma_states] for bs in beta_states] for V in ops])


def we_check() -> dict:
    """Wigner–Eckart (6.4.15)–(6.4.17)：六井的數值偶極矩陣元 ＝ 耦合係數 × 約化矩陣元（四方基底，投影法的 C）。
    三種躍遷：Eg → T1u、A1g → T1u（ℓ^γ = ℓ^α = 3）與 T1u → Eg（ℓ^γ = 2 ≠ ℓ^α = 3：照印的 1/ℓ^α 給錯的 R）。
    回傳 dict(E_T, A_T, T_E（各為 _we_fit 的結果）, ratios（以 ⟨T1 3|V3|E 1⟩ 為 1 的 (6.4.17) 比例）, ref, RA_abs)。"""
    X = [np.diag(SITES[:, i]) for i in range(3)]
    T = [SIX_STATES[("T_{1u}", k)] for k in range(3)]
    Eg = [SIX_STATES[("E_g", k)] for k in range(2)]
    A = [SIX_STATES[("A_{1g}", 0)]]
    S = ohset()
    out = dict(E_T=_we_fit(S, "T_{1u}", "E_g", "T_{1u}", _six_mel(T, Eg, X)),
               A_T=_we_fit(S, "T_{1u}", "A_{1g}", "T_{1u}", _six_mel(T, A, X)),
               T_E=_we_fit(S, "T_{1u}", "T_{1u}", "E_g", _six_mel(Eg, T, X)))
    el = six_dipole_elements()
    ref = el[("T_{1u}", 2, 2, "E_g", 0)]
    out["ratios"] = {(k + 1, j + 1): el.get(("T_{1u}", k, k, "E_g", j), 0.0) / ref for k in range(3) for j in range(2)}
    out["ref"] = ref
    out["RA_abs"] = el[("T_{1u}", 0, 0, "A_{1g}", 0)]
    out["offdiag_zero"] = all(abs(el.get(("T_{1u}", k, i, "E_g", j), 0.0)) < 1e-12 for k in range(3) for i in range(3) for j in range(2) if i != k)
    return out


WE_617_BOOK = {(1, 1): -0.5, (2, 1): -0.5, (3, 1): 1.0, (1, 2): math.sqrt(3) / 2, (2, 2): -math.sqrt(3) / 2, (3, 2): 0.0}   # (6.4.17)：(k, j) 以 ⟨T1 3|V3|E 1⟩ 為 1


def we_intensities() -> dict:
    """Fig. 6.4.1 的強度比：極化 1（或 2）時 (T1u → Eg 1):(T1u → Eg 2) = 1:3、極化 3 時 1:0；(6.4.19) 的比 (R_A/R_E)²（本模型 = 1/2）。"""
    w = we_check()
    r = w["ratios"]
    return dict(pol1=(r[(1, 1)] ** 2, r[(1, 2)] ** 2), pol3=(r[(3, 1)] ** 2, r[(3, 2)] ** 2), AE=(w["RA_abs"] / w["ref"]) ** 2)


def we_complex_check() -> dict:
    """〔注〕複數耦合係數時共軛在哪裡：T 群（12 個轉動）裡六井的 T1u → T、Eg → ¹E ⊕ ²E（複的一維）；偶極 T ⊗ T → ¹E 的 C 是複數。
    數值 ⟨¹E|V_i|T j⟩ ＝ C^{TT¹E *}_{ij}·R（共軛在 (6.4.15a) 的 C 上；(6.4.15b) 的 ΣC⟨⟩ 不共軛）；不共軛就不對。"""
    S = tset()
    X = [np.diag(SITES[:, i]) for i in range(3)]
    T = [SIX_STATES[("T_{1u}", k)] for k in range(3)]
    P1 = sum(np.conj(S.D["^1E"][g][0, 0]) * site_perm(g) for g in S.elems) / S.order
    e1 = P1 @ SIX_STATES[("E_g", 0)]
    e1 = e1 / np.linalg.norm(e1)
    r = _we_fit(S, "T", "T", "^1E", _six_mel([e1], T, X))
    C, _ = cg(S, "T", "T", ["^1E"])
    r["complex_C"] = bool(np.max(np.abs(C.imag)) > 1e-6)
    return r


# ======================================================================== §6.5A Lorentz 振子（H610）
FIG651 = dict(f0=5.0, gamma=0.5, A=20.0, fs=(4.0, 4.5, 4.9, 5.0, 5.25, 6.0))    # ω0/2π = 5 Hz、Γ = 0.5 s⁻¹、刺激 20 m/s²（＝ LORENTZ_A）
FIG651_RHO_BOOK = (4.0, 8.6, 38.2, 90.0, 161.9, 175.0)


def lorentz_steady(ws: float, w0: float, gamma: float, A: float = 1.0) -> tuple:
    """(6.5.5e,f)：穩態 x_S = R cos(ω_s t − ρ)（刺激 A cos ω_s t）→ (R, ρ 度，0–180)。"""
    R = A / math.hypot(w0 ** 2 - ws ** 2, 2 * gamma * ws)
    rho = math.degrees(math.atan2(2 * gamma * ws, w0 ** 2 - ws ** 2))
    return R, rho


def lorentz_transient(ws, w0, gamma, A=1.0, x0=0.0, v0=0.0) -> dict:
    """(6.5.7)–(6.5.8)：暫態 x_T = a e^{−Γt}cos(ω_Γ t − α)，ω_Γ = √(ω0² − Γ²)。用 a cos α、a sin α 兩式（共振時書上的 tan 式是 0/0）；
    另給書上的形式 (a_book, α_book)（arctan 取主值，a 可為負；ρ ≠ 90° 時兩者給同一條曲線）。"""
    R, rho = lorentz_steady(ws, w0, gamma, A)
    r = math.radians(rho)
    wg = math.sqrt(w0 ** 2 - gamma ** 2)
    c = x0 - R * math.cos(r)
    s = (v0 - R * ws * math.sin(r) + gamma * c) / wg
    a, al = math.hypot(c, s), math.atan2(s, c)
    out = dict(R=R, rho=rho, wg=wg, a=a, alpha=math.degrees(al))
    if abs(math.cos(r)) > 1e-12:
        alb = math.atan((v0 - R * ws * math.sin(r)) / ((x0 - R * math.cos(r)) * wg) + gamma / wg)
        out.update(a_book=(x0 - R * math.cos(r)) / math.cos(alb), alpha_book=math.degrees(alb))
    return out


def lorentz_closed(t, ws, w0, gamma, A=1.0, x0=0.0, v0=0.0):
    """(6.5.6) 的實部：x(t) = R cos(ω_s t − ρ) + a e^{−Γt}cos(ω_Γ t − α)。"""
    tr = lorentz_transient(ws, w0, gamma, A, x0, v0)
    t = np.asarray(t, float)
    return (tr["R"] * np.cos(ws * t - math.radians(tr["rho"]))
            + tr["a"] * np.exp(-gamma * t) * np.cos(tr["wg"] * t - math.radians(tr["alpha"])))


def lorentz_numeric(t, ws, w0, gamma, A=1.0, x0=0.0, v0=0.0):
    """第二種算法：scipy 直接積分 ẍ + 2Γẋ + ω0²x = A cos ω_s t（(6.5.1) 更正的 ẍ）。"""
    t = np.asarray(t, float)
    sol = solve_ivp(lambda tt, y: [y[1], A * math.cos(ws * tt) - 2 * gamma * y[1] - w0 ** 2 * y[0]],
                    (0, float(t[-1])), [x0, v0], t_eval=t, rtol=1e-10, atol=1e-12, method="DOP853")
    return sol.y[0]


def fig651_panels() -> list:
    """Fig. 6.5.1 的六幅：[(f_s, R, ρ 度, 最大 |x|（0–10 s）)]。"""
    p = FIG651
    w0 = 2 * math.pi * p["f0"]
    t = np.linspace(0, 10, 4001)
    out = []
    for fs in p["fs"]:
        ws = 2 * math.pi * fs
        R, rho = lorentz_steady(ws, w0, p["gamma"], p["A"])
        out.append((fs, R, rho, float(np.max(np.abs(lorentz_closed(t, ws, w0, p["gamma"], p["A"]))))))
    return out


def lorentz_trace(fs: float, tmax: float = 10.0, n: int = 801) -> tuple:
    """畫 Fig. 6.5.1 一幅用的取樣：(ts, 刺激 A cos ω_s t, 響應 x(t))（解析式；n 點）。"""
    p = FIG651
    ts = np.linspace(0, tmax, n)
    ws = 2 * math.pi * fs
    return ts, p["A"] * np.cos(ws * ts), lorentz_closed(ts, ws, 2 * math.pi * p["f0"], p["gamma"], p["A"])


def response(w, w0: float = 1.0, gamma: float = 0.1, A: float = 2.0) -> complex:
    """(6.5.5c)：R e^{iρ} = A/(ω0² − ω² − 2iΓω)（Fig. 6.5.2：ω0 = 1、Γ = 0.1、A = 2 ⇒ S = A/2ω0 = 1）。"""
    return A / (w0 ** 2 - w ** 2 - 2j * gamma * w)


def smith(delta: float, S: float = 1.0, gamma: float = 0.1) -> complex:
    """(6.5.9a)：近共振 R e^{iρ} ≈ S/(Δ − iΓ)。"""
    return S / (delta - 1j * gamma)


def smith_check(S: float = 1.0, gamma: float = 0.1) -> dict:
    """Smith 圓 (6.5.10)–(6.5.11)：Δ 變動時 S/(Δ−iΓ) 落在直徑 S/Γ（沿虛軸）的圓上；Γ 變動時落在直徑 S/Δ（沿實軸）的圓上；
    R = (S/Δ)cos ρ = (S/Γ)sin ρ；完整響應 vs 近似在 ω = 0.9…1.1 的最大相對誤差。"""
    ds = np.linspace(-1, 1, 41)
    z = np.array([smith(d, S, gamma) for d in ds])
    on_g = np.max(np.abs(np.abs(z - 1j * S / (2 * gamma)) - S / (2 * gamma)))
    d0 = 0.05
    zg = np.array([smith(d0, S, g) for g in np.linspace(0.01, 1, 30)])
    on_d = np.max(np.abs(np.abs(zg - S / (2 * d0)) - S / (2 * d0)))
    rel = []
    for wv in np.linspace(0.9, 1.1, 21):
        rel.append(abs(response(wv, 1, gamma, 2 * S) - smith(1 - wv, S, gamma)) / abs(response(wv, 1, gamma, 2 * S)))
    rr = smith(0.07, S, gamma)
    rho = np.angle(rr)
    return dict(on_gamma=float(on_g), on_delta=float(on_d), R_cos=float(abs(rr) - S / 0.07 * math.cos(rho)),
                R_sin=float(abs(rr) - S / gamma * math.sin(rho)), maxrel=float(max(rel)), diameter=S / gamma)


def fig652() -> dict:
    """Fig. 6.5.2（ω0 = 1、Γ = 0.1、A = 2）：Im R 的最大（≈ A/2Γω0 = 10）與 Re R 的極大、極小（在 ω0 ∓ Γ 附近，≈ ±5）、DC 響應 A/ω0² = 2。"""
    ws = np.linspace(0.5, 1.5, 200001)
    z = response(ws)
    return dict(im_max=float(np.max(z.imag)), w_im=float(ws[np.argmax(z.imag)]), re_max=float(np.max(z.real)),
                w_re_max=float(ws[np.argmax(z.real)]), re_min=float(np.min(z.real)), w_re_min=float(ws[np.argmin(z.real)]),
                dc=float(response(0.0).real))


# ======================================================================== §6.5B 物質對光（H611）
def epsilon(w, w0: float = 1.0, wp: float = 1.0, gamma: float = 0.1, book: bool = False) -> complex:
    """(6.5.21)：ε(ω) = 1 + ω_p²/(ω0² − ω² − 2iΓω)；book=True 照印（虛部分子少了 2）。"""
    den = (w0 ** 2 - w ** 2) ** 2 + 4 * gamma ** 2 * w ** 2
    im = (gamma if book else 2 * gamma) * wp ** 2 * w / den
    return 1 + wp ** 2 * (w0 ** 2 - w ** 2) / den + 1j * im


def fig654() -> dict:
    """Fig. 6.5.4（ω0 = ω_p = 1）：ε(0) = 2、截止 ω_c = √(ω0² + ω_p²) = √2；Γ = 0 時 Re ε < 0 恰在 (ω0, ω_c)；幾個點的值。"""
    wc = math.sqrt(2)
    grid = np.linspace(0.01, 3, 3000)
    neg = grid[np.array([epsilon(w, gamma=0).real < 0 for w in grid])]
    return dict(eps0=float(epsilon(0.0).real), wc=wc, opaque=(float(neg.min()), float(neg.max())),
                e12=float(epsilon(1.2, gamma=0).real), e2=float(epsilon(2.0, gamma=0).real), n2=math.sqrt(epsilon(2.0, gamma=0).real))


BOOK_CONST = dict(q=1.60e-19, m=9.11e-31, eps0=8.84e-12, c=3.0e8)        # 書上用的值（ε0 照印 8.84；標準 8.854×10⁻¹²）


def book_examples(eps0: float | None = None) -> dict:
    """§6.5B 的數值例（書上的常數）：α_e（氫，ω0 = 2×10¹⁶）、α_v（30 amu、ω0 = 9.4×10¹³＝500 cm⁻¹）與光頻時的 α_v、
    Γ（(6.5.30)，600 THz 綠光與 6 THz 紅外）、衰減 95% 的時間 3/2Γ 與期間的振盪次數、古典電子半徑與截面 πr²、Rayleigh 藍／紅比。"""
    q, m, c = BOOK_CONST["q"], BOOK_CONST["m"], BOOK_CONST["c"]
    e0 = eps0 or BOOK_CONST["eps0"]
    out = {}
    out["alpha_e"] = q ** 2 / (m * e0) / (2e16) ** 2
    out["alpha_v"] = q ** 2 / (5.0e-26 * e0) / (9.4e13) ** 2
    out["alpha_v_opt"] = q ** 2 / (5.0e-26 * e0) / ((9.4e13) ** 2 - (3.7e15) ** 2)
    out["w500"] = 2 * math.pi * c * 100 * 500
    out["w136"] = 13.6 * 1.602e-19 / 1.0546e-34
    for name, f in (("vis", 600e12), ("ir", 6e12)):
        w = 2 * math.pi * f
        G = q ** 2 * w ** 2 / (12 * math.pi * c ** 3 * e0 * m)
        out["G_" + name] = G
        out["t_" + name] = 3 / (2 * G)
        out["n_" + name] = f * 3 / (2 * G)
    out["w_vis"] = 2 * math.pi * 600e12
    re = q ** 2 / (4 * math.pi * e0 * m * c ** 2)
    out["r_e"] = re
    out["sigma"] = math.pi * re ** 2
    return out


RAYLEIGH_CARD = dict(blue=450e-9, red=700e-9, w0=2e16)


def rayleigh(blue=450e-9, red=700e-9, w0=2e16) -> dict:
    """Rayleigh 散射 (6.5.32)–(6.5.33)：藍／紅的功率比 = (λ_r/λ_b)⁴（ω ≪ ω0）；含共振因子 |ω0² − ω²|⁻² 的完整 Lorentz 比；
    第二種算法：直接把 (6.5.32) 的 ⟨P⟩ 相除（Γ = 0）。"""
    c = BOOK_CONST["c"]
    wb, wr = 2 * math.pi * c / blue, 2 * math.pi * c / red

    def P(w):
        return w ** 4 / abs(w0 ** 2 - w ** 2) ** 2

    return dict(law4=(red / blue) ** 4, full=P(wb) / P(wr), wb=wb, wr=wr)


def radiation_integral() -> float:
    """(6.5.28)：∫₀^π sin³θ dθ = 4/3（數值）⇒ P = q²ω⁴x²/(6πc³ε0)·cos²。"""
    return quad(lambda th: math.sin(th) ** 3, 0, math.pi)[0]


def phasor_points(eps: float = 0.2) -> list:
    """Fig. 6.5.8 的響應羅盤（e^{−iωt} 約定：複振幅乘 e^{+iδ} 是落後 δ）：振子落後 ρ、下游響應場再落後 π/2：
    E^R ∝ ε e^{i(ρ+π/2)}，和波 = 1 + E^R → [(ρ 度, 和波的振幅, 和波的相位落後 度)]：ρ = 0 延遲、90 衰減、180 超前、−90 放大。"""
    out = []
    for rho in (0, 90, 180, -90):
        z = 1 + eps * np.exp(1j * math.radians(rho + 90))
        out.append((rho, float(abs(z)), float(math.degrees(np.angle(z)))))
    return out


# ======================================================================== §6.5C 紅外與 Raman（H612）
XF6_NU_LABELS = {"A_{1g}": "nu_1", "E_g": "nu_2", "T_{1u}+": "nu_3", "T_{1u}-": "nu_4", "T_{2g}": "nu_5", "T_{2u}": "nu_6"}


def xf6_selection() -> dict:
    """XF6 的選擇規則（特徵標法）：振動 = 全部 21 − 平移 T1u − 轉動 T1g = A1g + Eg + T2g + T2u + 2T1u；
    紅外 = 振動 ∩ 偶極 {T1u}；Raman = 振動 ∩ 極化率 Sym²(T1u) = {A1g, Eg, T2g}；靜默 = 其餘（T2u）；互斥（g/u）。"""
    allf = {a: n for a, n in zip(SO.OH_IRREP_NAMES, SO.xf6_classify()["all"]) if n}
    vib = dict(allf)
    vib["T_{1u}"] -= 1
    vib["T_{1g}"] -= 1
    vib = {a: n for a, n in vib.items() if n}
    pol, _ = sym_squares(ohset(), "T_{1u}")
    ir = {a: n for a, n in vib.items() if a == "T_{1u}"}
    ram = {a: n for a, n in vib.items() if a in pol}
    silent = {a: n for a, n in vib.items() if a not in ir and a not in ram}
    return dict(all=allf, vib=vib, polarizability=pol, ir=ir, raman=ram, silent=silent, exclusive=not (set(ir) & set(ram)))


def xf6_selection_numeric(qF: float = -1.0, qX: float = 6.0, d_bond: float = 1e-4) -> dict:
    """第二種算法（數值）：每個 O_h 同型子空間（SO.xf6_isotypic）裡，剛性離子偶極 𝐩 = Σq_a𝐮_a 的導數、
    鍵極化率模型 α = Σ_鍵 𝐝𝐝ᵀ（X–F 與 F–F 鍵）的一階導數是否為零 → {α: (偶極導數的範數, 極化率導數的範數)}。"""
    B = np.zeros((21, 21))
    for c, (_, at, v) in enumerate(SO.XF6_COORDS_LIST):
        B[3 * at:3 * at + 3, c] = v
    pos0 = np.vstack([np.array(SO.F_SITES, float), np.zeros(3)])
    q = np.array([qF] * 6 + [qX])
    bonds = [(6, a) for a in range(6)] + [(a, b) for a in range(6) for b in range(a + 1, 6) if abs(pos0[a] @ pos0[b]) < 1e-9]

    def alpha(u21):
        pos = pos0 + (B @ u21).reshape(7, 3)
        return sum(np.outer(pos[b] - pos[a], pos[b] - pos[a]) for a, b in bonds)

    out = {}
    for a in SO.OH_IRREP_NAMES:
        P = SO.xf6_isotypic(a)
        if np.linalg.matrix_rank(P, tol=1e-8) == 0:
            continue
        U = np.linalg.svd(P)[0][:, :int(round(np.trace(P)))]
        dip = max(np.linalg.norm(q @ (B @ U[:, k]).reshape(7, 3)) for k in range(U.shape[1]))
        pol = max(np.linalg.norm(alpha(d_bond * U[:, k]) - alpha(-d_bond * U[:, k])) / (2 * d_bond) for k in range(U.shape[1]))
        out[a] = (float(dip), float(pol))
    return out


# ======================================================================== §6.6 多重激發（泛音）（H613）
def sym_power_chi(S: IrrepSet, a: str, n: int) -> np.ndarray:
    """Sym^n(𝒟^a) 的特徵標：Newton 遞迴 n·h_n = Σ_{k=1}^n p_k h_{n−k}，p_k(g) = Tr 𝒟(g)^k。"""
    St = S.stack(a)
    p = [None] + [np.einsum("gii->g", np.linalg.matrix_power(St, k)) for k in range(1, n + 1)]
    h = [np.ones(S.order, dtype=complex)]
    for m in range(1, n + 1):
        h.append(sum(p[k] * h[m - k] for k in range(1, m + 1)) / m)
    return h[n]


def sym_power_chi_monomial(S: IrrepSet, a: str, n: int, seed: int = 1) -> np.ndarray:
    """第二種算法：在 n 次單項式基底上直接求 𝒟(g) 的作用矩陣（隨機點上最小平方），取跡。"""
    ell = S.dim(a)
    mons = [m for m in itertools.product(range(n + 1), repeat=ell) if sum(m) == n]
    rng = np.random.default_rng(seed)
    pts = rng.normal(size=(3 * len(mons) + 5, ell))

    def ev(Y):
        return np.array([[np.prod(y ** np.array(m)) for m in mons] for y in Y])

    P0 = ev(pts)
    out = []
    for Dg in S.stack(a):
        Pg = ev(pts @ Dg.real)              # p_m(𝒟ᵀy)：座標跟基向量一樣變
        Mg = np.linalg.lstsq(P0, Pg, rcond=None)[0]
        out.append(np.trace(Mg))
    return np.array(out)


def overtone(alpha: str, n: int, S: IrrepSet | None = None) -> dict:
    """α 模的 n 次泛音（Sym^n）的 O_h 分解；u 的 irrep 時宇稱 = u^n（如 T2u：n = 2 → A1g + Eg + T2g、n = 3 → A1u + T1u + 2T2u）。"""
    S = S or ohset()
    return S.decompose(sym_power_chi(S, alpha, n))


def pascal(ell: int, n: int) -> int:
    """(4.4.64)：ℓ 維 irrep 的 n 次泛音簡併度 C(ℓ + n − 1, n)（ℓ = 2：2, 3, 4…；ℓ = 3：3, 6, 10, 15…）。"""
    return math.comb(ell + n - 1, n)


def combination(alpha: str, beta: str) -> dict:
    """組合頻 α ⊗ β（α ≠ β，沒有對易對稱的刪減）：Fig. 6.6.1(b) 的 Eg ⊗ T1u = T1u ⊕ T2u。"""
    S = ohset()
    return S.decompose(S.chi(alpha) * S.chi(beta))


def henon_heiles_check() -> bool:
    """Hénon–Heiles V = x³ − 3xy²（＝ ρ³cos3φ）在 C3v（σ₃ 沿 x）下不變；3x²y − y³（ρ³ sin 3φ）屬 A2（鏡面變號）。"""
    x, y = sp.symbols("x y", real=True)
    V = x ** 3 - 3 * x * y ** 2
    W = 3 * x ** 2 * y - y ** 3
    for g in SN.C3V.names:
        M = SN.C3V.mats[g]
        xp, yp = M.T * sp.Matrix([x, y])
        det = M.det()
        if sp.expand(V.subs({x: xp, y: yp}, simultaneous=True) - V) != 0:
            return False
        if sp.expand(W.subs({x: xp, y: yp}, simultaneous=True) - det * W) != 0:
            return False
    return True


# ======================================================================== §6.7 Jahn–Teller（H614）
JT_CARD = dict(j=1.0, k=1.0, r=-0.2)       # Fig. 6.7.5 的 j = k = 1；Renner 的 r 是本系列取的（jr < 0：與書上 p. 546 的極小位置、Fig. 6.7.4(d) 的 A1、A2 次序相同）


def jt_matrix(x1: float, x2: float, j: float = 1.0, k: float = 1.0, r: float = 0.0) -> np.ndarray:
    """(6.7.6a)：E⊗e 的 2×2 振電位能矩陣（電子基底 |E1⟩、|E2⟩）。"""
    p2 = x1 * x1 + x2 * x2
    return np.array([[j * x1 + k * p2 + r * (x1 * x1 - x2 * x2), -j * x2 + 2 * r * x1 * x2],
                     [-j * x2 + 2 * r * x1 * x2, -j * x1 + k * p2 - r * (x1 * x1 - x2 * x2)]])


def jt_sheets(rho, phi_deg, j: float = 1.0, k: float = 1.0, r: float = 0.0) -> tuple:
    """(6.7.8)：ε± = kρ² ± (j²ρ² + 2jrρ³cos3φ + r²ρ⁴)^{1/2}（r = 0 時 kρ² ± |j|ρ）。"""
    rho = np.asarray(rho, float)
    c3 = np.cos(3 * np.radians(phi_deg))
    root = np.sqrt(np.maximum(j * j * rho ** 2 + 2 * j * r * rho ** 3 * c3 + r * r * rho ** 4, 0))
    return k * rho ** 2 + root, k * rho ** 2 - root


def jt_min(j: float = 1.0, k: float = 1.0) -> tuple:
    """墨西哥帽（r = 0）的谷底：ρ0 = |j|/2k、深度 ε−(ρ0) = −j²/4k。"""
    return abs(j) / (2 * k), -j * j / (4 * k)


def renner_extrema(j: float = 1.0, k: float = 1.0, r: float = 0.2) -> dict:
    """Renner 項（r ≠ 0）：jr > 0 時下層的三個極小在 φ = 0°, 120°, 240°（ρ = j/2(k−r)、ε = −j²/4(k−r)），鞍點在 60°, 180°, 300°
    （ρ = j/2(k+r)、ε = −j²/4(k+r)），位障 j²r/2(k²−r²)；jr < 0 時兩者對調（書上的說法相當於 jr < 0）。第二種算法：格點上直接找。"""
    sg = 1 if j * r > 0 else -1
    mins = [0, 120, 240] if sg > 0 else [60, 180, 300]
    out = dict(min_phi=mins, rho_min=abs(j) / (2 * (k - abs(r))), e_min=-j * j / (4 * (k - abs(r))),
               rho_saddle=abs(j) / (2 * (k + abs(r))), e_saddle=-j * j / (4 * (k + abs(r))))
    out["barrier"] = out["e_saddle"] - out["e_min"]
    rhos = np.linspace(0.01, 1.5, 1500)
    phis = np.arange(0, 360, 1.0)
    grid = np.array([jt_sheets(rhos, p, j, k, r)[1] for p in phis])
    prof = grid.min(axis=1)
    loc = [int(phis[i]) for i in range(360) if prof[i] < prof[i - 1] and prof[i] <= prof[(i + 1) % 360]]
    out["grid_min_phi"] = loc
    out["grid_min"] = float(prof.min())
    out["grid_saddle"] = float(prof.max())
    return out


def jt_eigvec_track(j: float = 1.0, phis=None) -> tuple:
    """(6.7.9)：r = 0 時沿谷底繞一圈追蹤下層電子態（數值本徵向量，連續選號）→ (φ 陣列, 向量陣列, 繞完的號 −1)；
    第二種算法：|ε−⟩ = sin(φ/2)|E1⟩ + cos(φ/2)|E2⟩ 的解析式。"""
    phis = np.linspace(0, 360, 73) if phis is None else np.asarray(phis)
    vecs = []
    prev = None
    for p in phis:
        w, v = np.linalg.eigh(jt_matrix(math.cos(math.radians(p)), math.sin(math.radians(p)), j, 0.0, 0.0))
        u = v[:, 0]
        if prev is None:
            u = u * np.sign(u[1]) if abs(u[1]) > 1e-9 else u
        elif u @ prev < 0:
            u = -u
        vecs.append(u)
        prev = u
    vecs = np.array(vecs)
    sign = float(np.sign(vecs[-1] @ vecs[0]))
    ana = np.array([[math.sin(math.radians(p) / 2), math.cos(math.radians(p) / 2)] for p in phis])
    return phis, vecs, sign, ana


def _ho_ops(nmax: int, w: float):
    a = np.diag(np.sqrt(np.arange(1, nmax + 1)), 1)
    x = (a + a.T) / math.sqrt(2 * w)
    hn = w * (np.arange(nmax + 1) + 0.5)
    return x, hn


def jt_vibronic_levels(j: float = 1.0, k: float = 1.0, r: float = 0.0, nmax: int = 26, count: int = 12) -> list:
    """振電能階（數值對角化）：H = ½(p₁² + p₂²) + kρ² + j(x₁σ_z − x₂σ_x) + r[(x₁² − x₂²)σ_z + 2x₁x₂σ_x]（ħ = 質量 = 1，
    (6.7.6a)＋(6.7.14)），二維振子基底 |n₁n₂⟩（ω = √(2k)，每維 0…nmax）⊗ 電子 |E1⟩, |E2⟩ →
    [(E, 簡併度, C3v 標記 "E"/"A_1"/"A_2")]（鏡面 x₂ → −x₂ ⊗ σ_z 的本徵值分 A1、A2）。"""
    w = math.sqrt(2 * k)
    x, hn = _ho_ops(nmax, w)
    n = nmax + 1
    I1 = np.eye(n)
    X1, X2 = np.kron(x, I1), np.kron(I1, x)
    H0 = np.diag(np.add.outer(hn, hn).ravel())
    sz, sx = np.diag([1.0, -1.0]), np.array([[0, 1.0], [1, 0]])
    H = np.kron(H0, np.eye(2)) + j * (np.kron(X1, sz) - np.kron(X2, sx))
    if r:
        H += r * (np.kron(X1 @ X1 - X2 @ X2, sz) + np.kron(2 * X1 @ X2, sx))
    ev, V = eigh(H, subset_by_index=[0, 2 * count + 1])
    par = np.diag((-1.0) ** np.arange(n))
    sig = np.kron(np.kron(I1, par), sz)
    out = []
    i = 0
    while i < len(ev) and len(out) < count:
        jdeg = 1
        while i + jdeg < len(ev) and abs(ev[i + jdeg] - ev[i]) < 1e-7:
            jdeg += 1
        if jdeg == 2:
            lab = "E"
        else:
            s = float(V[:, i] @ sig @ V[:, i])
            lab = "A_1" if s > 0 else "A_2"
        out.append((float(ev[i]), jdeg, lab))
        i += jdeg
    return out


def jt_radial_levels(l_: float = 0.5, j: float = 1.0, k: float = 1.0, count: int = 4, L: float = 8.0, N: int = 2400,
                     book: bool = False) -> np.ndarray:
    """第二種算法：r = 0 時固定半整數 l 的耦合徑向方程（更正的 (6.7.15)：位能在 ½[⋯] 外；交叉項 ∓il/2ρ²，號不影響能譜）有限差分 →
    最低 count 個能量。把 (R⁺, R⁻) 轉 45°：兩個通道的離心項變成 m = l ± ½ 的二維徑向算子（交錯格點 ρ_i = (i − ½)h，
    對 m = 0 也是二階精度），jρ 變成通道間的耦合。book=True 照印（位能也在 ½ 裡，等於 k、j 減半）。"""
    h = L / N
    rho = h * (np.arange(1, N + 1) - 0.5)
    half = h * np.arange(0, N + 1)                      # ρ_{i−½}
    kk, jj = (k / 2, j / 2) if book else (k, j)
    lap_d = -(half[1:] + half[:-1]) / (h * h * rho)
    lap_o = half[1:-1] / (h * h * np.sqrt(rho[:-1] * rho[1:]))
    diag = np.empty(2 * N)
    for c, m in enumerate((l_ + 0.5, abs(l_ - 0.5))):
        diag[c::2] = -0.5 * lap_d + m * m / (2 * rho ** 2) + kk * rho ** 2
    band = np.zeros((3, 2 * N))
    band[2] = diag
    off1 = np.zeros(2 * N)
    off1[1::2] = jj * rho
    band[1, 1:] = off1[1:]
    off2 = np.zeros(2 * N)
    off2[2:] = -0.5 * np.repeat(lap_o, 2)
    band[0, 2:] = off2[2:]
    return eig_banded(band, lower=False, select="i", select_range=(0, count - 1), eigvals_only=True)


def jt_effective(rho, l_: float, sheet: int = -1, j: float = 1.0, k: float = 1.0, book_fig: bool = False):
    """Born–Oppenheimer 的有效位能：(l² + ¼)/2ρ² + kρ² ± jρ（更正後，ħ = 質量 = 1）；book_fig=True 是 Fig. 6.7.5 畫的
    (l² + ¼)/ρ² + kρ² ± jρ（＝ 照印 (6.7.15) 括號裡的量）。"""
    rho = np.asarray(rho, float)
    cent = (l_ ** 2 + 0.25) / rho ** 2
    return (cent if book_fig else cent / 2) + k * rho ** 2 + sheet * j * rho


def jt_effective_min(l_: float, sheet: int = -1, j: float = 1.0, k: float = 1.0, book_fig: bool = False) -> tuple:
    rhos = np.linspace(0.2, 4, 380001)
    v = jt_effective(rhos, l_, sheet, j, k, book_fig)
    i = int(np.argmin(v))
    return float(rhos[i]), float(v[i])


# ======================================================================== 習題 6.1–6.6
def problem_61() -> dict:
    """習題 6.1：子群鏈標記的 E⊗E 耦合係數表（投影法、書上的相位約定）：
    (a) D3 ⊃ C2（駐波 ＝ (6.1.35a)）、(b) D3 ⊃ C3（行進波 1_3、2_3）、(c) D4 ⊃ C2（C2 = {1, R1²}；＝ (6.1.23)）、(d) D4 ⊃ C4（行進波 1_4、3_4）。
    回傳 {key: (精確的 C, 行標記)}。"""
    out = {}
    C, cols = cg(d3set(), "E", "E", ["A_1", "A_2", "E"])
    out["a"] = (exact_mat(C), cols)
    U3 = np.array([[-1, 1], [-1j, -1j]]) / math.sqrt(2)          # (x, y) → (1_3, 2_3)（仿 (4.2.36)）
    S3 = d3set()
    D = {a_: {g: (U3.conj().T @ S3.D[a_][g] @ U3 if a_ == "E" else S3.D[a_][g]) for g in S3.elems} for a_ in S3.names}
    C, cols = cg(IrrepSet("D3/C3", S3.elems, D), "E", "E", ["A_1", "A_2", "E"])
    out["b"] = (exact_mat(C), cols)
    C, cols = cg(d4set(), "E", "E", ["A_1", "B_1", "A_2", "B_2"])
    out["c"] = (exact_mat(C), cols)
    U4 = np.array(SO.TMOV2.evalf(), dtype=complex)
    S4 = d4set()
    D = {a_: {g: (U4.conj().T @ S4.D[a_][g] @ U4 if a_ == "E" else S4.D[a_][g]) for g in S4.elems} for a_ in S4.names}
    C, cols = cg(IrrepSet("D4/C4", S4.elems, D), "E", "E", ["A_1", "B_1", "A_2", "B_2"])
    out["d"] = (exact_mat(C), cols)
    return out


def _d3_xaxis():
    """習題 6.2 的 D3：C3 繞 z、C2 繞 x（x 沿對稱軸指向頂點）。"""
    c3 = SN.rot3((0, 0, 1), 120)
    c2 = SN.rot3((1, 0, 0), 180)
    return SN.FiniteGroup.from_generators({"c": c3, "p": c2})


def problem_62() -> dict:
    """習題 6.2（x 沿 C2 軸指向頂點、z 垂直平面的 D3）：(a)(b) 一次：(x, y) ∈ E、z ∈ A2；二次：x² + y²、z² ∈ A1，
    (x² − y², −2xy)、(xz, −yz)…∈ E（9 維張量 = 2A1 + A2 + 3E，其中反對稱的 A2 + E 不給多項式）；
    (c) 三次 A1：x³ − 3xy² = ρ³cos3φ（Hénon–Heiles）；(d) 三次 A2：3x²y − y³ = ρ³sin3φ。回傳各多項式的 irrep（逐元素代入核對）。"""
    G = _d3_xaxis()
    x, y, z = sp.symbols("x y z", real=True)

    def label(p):
        vals = []
        for g in G.names:
            M = G.mats[g]
            sub = dict(zip((x, y, z), list(M.T * sp.Matrix([x, y, z]))))
            q = sp.expand(p.subs(sub, simultaneous=True))
            ratio = sp.simplify(q / p) if q != 0 else 0
            vals.append(ratio)
        if all(v == 1 for v in vals):
            return "A_1"
        if all(v in (1, -1) for v in vals):
            return "A_2"
        return "?"

    polys = {"z": z, "x2+y2": x ** 2 + y ** 2, "z2": z ** 2, "HH": x ** 3 - 3 * x * y ** 2, "HH2": 3 * x ** 2 * y - y ** 3,
             "z(x2+y2)": z * (x ** 2 + y ** 2)}
    out = {k: label(p) for k, p in polys.items()}
    vec = np.array([[float(np.trace(np.array(G.mats[g].evalf(), dtype=float))) for g in G.names]])
    chis = {}
    for g in G.names:
        M = np.array(G.mats[g].evalf(), dtype=float)
        chis[g] = (np.trace(M), np.trace(M) ** 2, (np.trace(M) ** 2 + np.trace(M @ M)) / 2)
    out["tensor9"] = "2A_1 + A_2 + 3E"
    out["poly2"] = "2A_1 + 2E"
    out["_chis"] = chis
    out["_vec"] = vec
    return out


def problem_63() -> dict:
    """習題 6.3（第 3 章 X3 模型，C3v，四個能階 A1、A2、E、E，六個態；ω 由 SN.six_dof_omega2）：從基態 A1 出發，
    躍遷算子屬 γ 時末態必須屬 γ ⇒ (a) E1_x（x ∈ E）：兩個 E 能階；(b) M1_z（R_z ∈ A2）：A2；(c) E2_xy（xy ∈ E）：兩個 E。"""
    S = c3vset()
    om = SN.six_dof_omega2()
    ops = {"E1_x": "E", "M1_z": "A_2", "E2_xy": "E"}
    out = {}
    for name, op in ops.items():
        fin = [a for a in S.names if S.decompose(S.chi("A_1") * S.chi(op) * S.chi(a).conj()).get("A_1", 0)]
        out[name] = fin
    out["omega"] = {k: (tuple(math.sqrt(float(v)) for v in val) if isinstance(val, tuple) else math.sqrt(float(val)))
                    for k, val in om.items()}
    return out


def _c4v_2d():
    names, mats = [], []
    for kk in range(4):
        names.append(f"C4^{kk}")
        mats.append(np.array(SN.rot2(90 * kk).evalf(), dtype=float))
    for deg in (0, 45, 90, 135):
        names.append(f"s{deg}")
        mats.append(np.array(SN.refl2(deg).evalf(), dtype=float))
    return names, mats


def problem_64() -> dict:
    """習題 6.4（X4 正方形模型，C4v ≅ D4，質點在 x、y 軸上，σ_v 過質點）：8 個平面座標 = A1 + A2 + B1 + B2 + 2E；
    從基態：(a) E1_x → E；(b) M1_z（R_z ∈ A2）→ A2；(c) E2_xy（xy ∈ B2）→ B2（σ_v 改過鍵中點時 B1、B2 互換）。"""
    names, mats = _c4v_2d()
    sites = [np.array(s, float) for s in ((1, 0), (0, 1), (-1, 0), (0, -1))]

    def chi_ir(kind, M):
        if kind == "A1":
            return 1.0
        if kind == "A2":
            return np.linalg.det(M)
        if kind == "E":
            return np.trace(M)
        Q = np.diag([1, -1.0]) if kind == "B1" else np.array([[0, 1.0], [1, 0]])     # B1 ~ x² − y²、B2 ~ xy
        return _coef(M, Q)

    irs = ["A1", "A2", "B1", "B2", "E"]
    tab = {ir: np.array([chi_ir(ir, M) for M in mats]) for ir in irs}
    disp = np.array([sum(np.trace(M) for s in sites if np.allclose(M @ s, s)) for M in mats])
    dec = {ir: int(round(float(tab[ir] @ disp) / 8)) for ir in irs}
    xy = np.array([_coef(M, np.array([[0, 1.0], [1, 0]])) for M in mats])
    lab_xy = next(ir for ir in irs if np.allclose(tab[ir], xy))
    rz = np.array([np.linalg.det(M) for M in mats])
    lab_rz = next(ir for ir in irs if np.allclose(tab[ir], rz))
    return dict(modes=dec, E1_x="E", M1_z=lab_rz, E2_xy=lab_xy)


def problem_65() -> dict:
    """習題 6.5（Fig. 6.4.1 改成三方 O_h ↓ D3d）：T1u → A2u（沿 v3 = (111)/√3）⊕ Eu（v1, v2），A1g → A1g、Eg → Eg（不分裂）。
    六井模型的偶極強度（單位：q = 1、井距 1）：{(初態, 末態, 極化): |⟨⟩|²}：A1g → T1u(v_a) 用 v_a 極化 1/3；
    T1u(v3) → Eg：v3 極化 0、v1/v2 極化 1/3；T1u(v1 或 v2) → Eg：v3 極化 1/3、同平面極化 1/6（Eg 兩個夥伴相加）。"""
    T43 = np.array(SO.T43.evalf(), dtype=float)
    X = [np.diag(SITES[:, i]) for i in range(3)]
    t = np.array([SIX_STATES[("T_{1u}", k)] for k in range(3)]).T
    e = [SIX_STATES[("E_g", k)] for k in range(2)]
    A = SIX_STATES[("A_{1g}", 0)]
    out = {}
    for a_ in range(3):
        st = t @ T43[:, a_]
        for b in range(3):
            op = sum(T43[i, b] * X[i] for i in range(3))
            out[("A1g", f"v{a_ + 1}", f"v{b + 1}")] = float((st @ op @ A) ** 2)
            out[(f"v{a_ + 1}", "Eg", f"v{b + 1}")] = float(sum((st @ op @ ee) ** 2 for ee in e))
    return {kk: Fraction(v).limit_denominator(48) for kk, v in out.items()}


def d6h_group() -> tuple:
    """D6h 的 24 個元素：(k, s, i) ＝ I^i·R_z(60k)·C2x^s（s = 1 的是繞方位角 30k° 的 C2：k 偶 C2'（過原子）、k 奇 C2''）→ (名, 3×3)。"""
    out = []
    for i in (0, 1):
        for s in (0, 1):
            for kk in range(6):
                M = SN_rot3((0, 0, 1), 60 * kk) @ (SN_rot3((1, 0, 0), 180) if s else np.eye(3))
                if i:
                    M = -M
                out.append(((kk, s, i), M))
    return out


D6H_IRREPS = ["A_{1g}", "A_{2g}", "B_{1g}", "B_{2g}", "E_{1g}", "E_{2g}", "A_{1u}", "A_{2u}", "B_{1u}", "B_{2u}", "E_{1u}", "E_{2u}"]


def d6h_chi(ir: str, key: tuple) -> float:
    """D6h 的特徵標（慣例：B1 對 C2'（過原子）為 +1）：key = (k, s, i)。"""
    kk, s, i = key
    base, par = ir.replace("{", "").replace("}", "")[:-1], ir[-2]
    if s == 0:
        v = {"A_1": 1, "A_2": 1, "B_1": (-1) ** kk, "B_2": (-1) ** kk, "E_1": 2 * math.cos(math.radians(60 * kk)),
             "E_2": 2 * math.cos(math.radians(120 * kk))}[base]
    else:
        through_atom = kk % 2 == 0
        v = {"A_1": 1, "A_2": -1, "B_1": 1 if through_atom else -1, "B_2": -1 if through_atom else 1, "E_1": 0, "E_2": 0}[base]
    if i and par == "u":
        v = -v
    return float(v)


def problem_66() -> dict:
    """習題 6.6（D6h）：(a) 極向量 z ∈ A2u、(x, y) ∈ E1u；(b) 軸向量 R_z ∈ A2g、(R_x, R_y) ∈ E1g；
    (c) 六角 X6：12 個振動 = A1g + B2g + 2E2g + B1u + B2u + E1u + E2u；(d) 苯 C6H6：30 個 = 2A1g + A2g + 2B2g + E1g + 4E2g + A2u + 2B1u + 2B2u + 3E1u + 2E2u，
    古典振動問題最大的方塊是 E2g 的 4×4；(e) 紅外 A2u + E1u（X6：1 個，苯：4 個）、Raman A1g + E1g + E2g（X6：3 個，苯：7 個），互斥。"""
    G = d6h_group()
    n = len(G)

    def dec(chi):
        out = {}
        for ir in D6H_IRREPS:
            f = sum(d6h_chi(ir, key) * c for (key, M), c in zip(G, chi)) / n
            if abs(f) > 1e-9:
                out[ir] = int(round(f))
        return out

    vec = [np.trace(M) for _, M in G]
    ax = [np.linalg.det(M) * np.trace(M) for _, M in G]

    def gamma3n(radii):
        pts = [rr * np.array([math.cos(math.radians(60 * a)), math.sin(math.radians(60 * a)), 0]) for rr in radii for a in range(6)]
        return [sum(1 for p in pts if np.allclose(M @ p, p)) * np.trace(M) for _, M in G]

    def vib(radii):
        f = dec(gamma3n(radii))
        for ir in ("A_{2u}", "E_{1u}", "A_{2g}", "E_{1g}"):
            f[ir] -= 1
        return {k: v for k, v in f.items() if v}

    x6, bz = vib([1.0]), vib([1.39, 2.48])
    ir_set, ram_set = ("A_{2u}", "E_{1u}"), ("A_{1g}", "E_{1g}", "E_{2g}")
    count = lambda f, s: sum(f.get(k, 0) for k in s)   # noqa: E731
    return dict(polar=dec(vec), axial=dec(ax), X6=x6, benzene=bz, X6_ir=count(x6, ir_set), X6_raman=count(x6, ram_set),
                bz_ir=count(bz, ir_set), bz_raman=count(bz, ram_set), bz_block=max(bz.values()),
                n_vib=(sum(v * (2 if k.startswith("E") else 1) for k, v in x6.items()),
                       sum(v * (2 if k.startswith("E") else 1) for k, v in bz.items())))


# ======================================================================== 卡片上的數字
LORENTZ_A = FIG651["A"]   # 本系列取的刺激振幅 20 m/s²：Fig. 6.5.1 的縱軸 ±100，由響應的大小（共振時 10 s 後約 0.63 m 出框）推得約 20


def _ex(x) -> str:
    """精確數 → 卡片字串（sympy 的 str：1/2、sqrt(3)/2、-sqrt(6)/3）。"""
    return str(x) if isinstance(x, (sp.Basic, Fraction, int)) else str(exact(x))


def _tab(M) -> str:
    M = np.asarray(M)
    return " / ".join(" ".join(_ex(x) for x in row) for row in M)


def _dec(d: dict) -> str:
    return " + ".join((f"{n}" if n > 1 else "") + a for a, n in d.items())


def _fg(x) -> str:
    """短小數：2.5、-3.5、0.25；四位以上照 4 位（10/3 → 3.3333）。"""
    x = float(x)
    if abs(x) < 5e-13:
        x = 0.0
    return f"{x:.4f}".rstrip("0").rstrip(".") if abs(x - round(x, 4)) < 1e-12 else f"{x:.4f}"


def _e2(x) -> str:
    return f"{x:.2e}"


def card_values() -> dict:
    """plan_ch06.md §6 每張卡上的數字與導讀第 6 章的數字（字串；小數照卡上位數）。"""
    out = {}
    # H601
    out["H601_count"] = f"{len(PAIRS)} {SO.OH.order ** 2}"
    v = np.zeros(36)
    v[0] = 1
    out["H601_R1"] = str(PAIRS[int(np.argmax(outer_op("1", "R1") @ v))])
    out["H601_kron"] = _tab(kron_example()["kron"])
    out["H601_chi"] = f"{int(round(np.trace(_dnum('T_{1u}', 'r1', 'tetragonal')).real)) ** 2} {len(outer_frequencies())}"
    la, lb = six_levels(0, TWO_CARD["SA"]), six_levels(0, TWO_CARD["SB"])
    out["H601_single"] = " ".join(_fg(la[a]) for a in ("A_{1g}", "T_{1u}", "E_g")) + " | " + \
        " ".join(_fg(lb[a]) for a in ("A_{1g}", "T_{1u}", "E_g"))
    out["H601_levels"] = " ".join(f"{_fg(e)}({d})" for e, a, b, d in two_particle_levels(TWO_CARD["SA"], TWO_CARD["SB"]))
    # H602
    Vm = pair_V(**PAIR_V)
    out["H602_V"] = " ".join(_fg(Vm[PAIRS.index(pq), PAIRS.index(pq)]) for pq in ((1, 1), (1, 6), (1, 2)))
    out["H602_chiprod"] = " ".join(str(x) for x in char_product_rows()["T1u^2"])
    out["H602_freq"] = _dec(diag_frequencies("T_{1u}", "T_{1u}"))
    out["H602_N"] = " ".join(_ex(n) for _, n in cg_norms(ohset("nonkosher"), "T_{1u}", "T_{1u}", "T_{1g}")[:6])
    out["H602_cg"] = _tab(cg_book("619c")[0].real)
    sp1 = t1u_split()
    sp0 = t1u_split(dict(PAIR_V, K=0.0, W=0.0))
    order = ("A_{1g}", "E_g", "T_{2g}", "T_{1g}")
    out["H602_split"] = " ".join(_f4(sp1[a]) for a in order) + " | " + " ".join(_f4(sp0[a]) for a in order)
    # H603
    irs = ("A_1", "A_2", "B_1", "B_2")
    out["H603_d2"] = " / ".join(" ".join(d2_table()[(a, b)] for b in irs) for a in irs)
    out["H603_623"] = _tab(cg_book("623")[0].real)
    out["H603_626"] = " | ".join(_tab(cg_book("626")[x][0].real) for x in ("A_2", "B_1", "B_2"))
    so = solve_6129()
    out["H603_rel"] = " ".join(f"{k}={so['relations'][k]}" for k in ("b", "d", "e", "g", "h", "j", "k"))
    out["H603_val"] = " ".join(f"{k}={so['values'][k]}" for k in ("a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k"))
    out["H603_bookd"] = str(solve_6129(True)["consistent"])
    out["H603_635a"] = _tab(cg_book("635a")[0].real)
    out["H603_635b"] = _tab(cg_book("635b")[0].real)
    out["H603_636"] = _tab(cg_book("636")[0].real)
    # H604
    out["H604_scalar"] = _ex(cg_book("619c")[0][0, 0].real)
    out["H604_fs"] = " | ".join(f"{S_.label}: " + " ".join(str(fs_indicator(S_, a)) for a in S_.names)
                                for S_ in (oset(), c3set(), tset(), qset()))
    out["H604_c3"] = _dec(c3set().decompose(c3set().chi("1_3") * c3set().chi("2_3")))
    O_ = oset()
    out["H604_sq"] = " | ".join(f"{a}: " + _dec(sym_squares(O_, a)[0]) + " / " + _dec(sym_squares(O_, a)[1]) for a in ("E", "T_1", "T_2"))
    Tt = tset()
    out["H604_TxT"] = _dec(Tt.decompose(Tt.chi("T") ** 2))
    # H605
    o3 = o3_labels()
    out["H605_o3"] = " | ".join(f"{k}: {_dec(o3[k])}" for k in ("0^+", "1^+", "2^+"))
    pt = poly_from_tensors()
    out["H605_poly"] = " | ".join(str(pt[k]) for k in POLY_639)
    out["H605_R3"] = _tab(_dnum("E_g", "R3", "tetragonal").real)
    # H606
    st = stress_examples()
    out["H606_stress"] = " | ".join(f"{k}: " + " ".join(_f4(x) for x in st[k].values()) for k in ("pressure", "uniaxial", "shear12", "torque"))
    # H607
    ec = elastic_counts()
    out["H607_counts"] = f"{ec['none']['81']} {ec['none']['36']} {ec['none']['21']} | Oh {ec['Oh']['81']} {ec['Oh']['36']} {ec['Oh']['21']} | " \
                         f"iso {ec['iso']['21']} | D4h {ec['D4h']['36']} {ec['D4h']['21']} | C3v {ec['C3v']['36']} {ec['C3v']['21']}"
    for key, jkl in MODULI_CARD.items():
        kA, kE, kT = spring_model(*jkl)
        m = moduli(kA, kE, kT)
        out["H607_" + key] = f"{kA} {kE} {kT} | " + " ".join(_fg(m[x]) for x in ("B", "mu", "Y", "sigma", "A", "k1111", "k1122", "k1212"))
    out["H607_kittel"] = " ".join(f"{k} {anisotropy(*CUBIC_KITTEL[k]):.2f}" for k in ("diamond", "Si", "Ge", "NaCl", "Au", "Pb", "W"))
    # H608
    fi = c3v_field_irreps()
    out["H608_reps"] = f"E: {_dec(fi['E'])} | B: {_dec(fi['B'])}"
    out["H608_chars"] = " | ".join(" ".join(_fg(np.trace(c3v_field(g, kd))) for g in ("1", "r", "s3")) for kd in ("E", "B"))
    f6 = fig6311()
    out["H608_r"] = " ".join(_f4(x) for x in f6[("E", "r")][0][:2]) + " | " + " ".join(_f4(x) for x in f6[("B", "r")][0][:2])
    out["H608_s3"] = " ".join(_fg(f6[(kd, "s3")][c][c]) for kd in ("E", "B") for c in range(3))
    out["H608_P12"] = " ".join(_fg(x) for x in p3_12("B")) + " | " + " ".join(str(x) for x in p3_12_book()) + " | " + \
        " ".join(str(x) for x in p3_12_book(SQ3 / 2))
    stc = striction_counts()
    out["H608_strict"] = f"C3v {stc[('C3v', 'E')]} {stc[('C3v', 'B')]} | Oh {stc[('Oh', 'E')]} {stc[('Oh', 'B')]}"
    # H609
    inv = six_invariants()
    out["H609_inv"] = " | ".join(_tab(inv[k][0:1] * f) for k, f in (("A", 1), ("T", math.sqrt(3)), ("E", math.sqrt(2))))
    lv = six_levels(**SIX_CARD)
    back = six_from_levels(lv["A_{1g}"], lv["T_{1u}"], lv["E_g"])
    out["H609_levels"] = " ".join(_fg(lv[a]) for a in ("A_{1g}", "T_{1u}", "E_g")) + " | " + " ".join(_fg(back[k]) for k in ("H", "T", "S"))
    xp = xf6_param_counts()
    out["H609_param"] = f"{param_count({'A': 1, 'T': 1, 'E': 1})} {xp['vib']} {xp['full']}"
    S_ = ohset()
    out["H609_multi"] = _dec(S_.decompose(S_.chi("T_{1u}") ** 3)) + " | " + _dec(overtone("T_{1u}", 3))
    w = we_check()
    el = six_dipole_elements()
    keys = (("T_{1u}", 2, 2, "E_g", 0), ("T_{1u}", 0, 0, "E_g", 0), ("T_{1u}", 0, 0, "E_g", 1), ("T_{1u}", 1, 1, "E_g", 1),
            ("T_{1u}", 2, 2, "E_g", 1), ("T_{1u}", 0, 0, "A_{1g}", 0))
    out["H609_mel"] = " ".join(_f4(el.get(kk, 0.0)) for kk in keys)
    out["H609_ratio"] = " ".join(_f4(w["ratios"][(k, j)]) for j in (1, 2) for k in (1, 2, 3))
    out["H609_red"] = f"{_f4(abs(w['E_T']['R']))} {_f4(abs(w['A_T']['R']))} | {_f4(abs(w['T_E']['R']))} {_f4(abs(w['T_E']['R_book']))}"
    wi = we_intensities()
    out["H609_int"] = f"{_f4(wi['pol1'][0])} {_f4(wi['pol1'][1])} | {_f4(wi['pol3'][0])} {_f4(wi['pol3'][1])} | {_f4(wi['AE'])}"
    # H610
    p = FIG651
    w0 = 2 * math.pi * p["f0"]
    out["H610_w"] = f"{_f4(w0)} {_f4(math.sqrt(w0 ** 2 - p['gamma'] ** 2))}"
    out["H610_rho"] = " ".join(_f2(r) for _, _, r, _ in fig651_panels())
    out["H610_R"] = " ".join(_f4(lorentz_steady(2 * math.pi * fs, w0, p["gamma"], LORENTZ_A)[0]) for fs in p["fs"])
    out["H610_res"] = f"{1 / (2 * p['gamma'] * w0):.5f}"
    f2 = fig652()
    out["H610_smith"] = f"{_f4(smith_check()['diameter'])} | {f2['im_max']:.2f} {f2['w_im']:.3f} | {f2['re_max']:.2f} {f2['w_re_max']:.3f} | " \
                        f"{f2['re_min']:.2f} {f2['w_re_min']:.3f} | {_f4(f2['dc'])}"
    # H611
    f4 = fig654()
    e1, e1b = epsilon(1.0), epsilon(1.0, book=True)
    out["H611_eps"] = f"{_f4(f4['eps0'])} {_f4(f4['wc'])} {_f4(f4['e12'])} {_f4(f4['e2'])} {_f4(f4['n2'])} | " \
                      f"{_fg(e1.real)}+{_fg(e1.imag)}i {_fg(e1b.real)}+{_fg(e1b.imag)}i"
    be = book_examples()
    out["H611_alpha"] = f"{_e2(be['alpha_e'])} {be['alpha_v']:.1e} {be['alpha_v_opt']:.1e}"
    out["H611_gamma"] = f"{be['G_vis']:.1e} {be['t_vis']:.1e} {be['n_vis']:.1e} | {be['t_ir']:.1e} {be['n_ir']:.1e}"
    out["H611_cross"] = f"{_e2(be['r_e'])} {be['sigma']:.1e}"
    ry = rayleigh(**RAYLEIGH_CARD)
    out["H611_rayleigh"] = f"{ry['law4']:.2f} {ry['full']:.2f}"
    out["H611_phasor"] = " | ".join(f"{_f4(a)} {_f2(ph)}" for _, a, ph in phasor_points())
    # H612
    sel = xf6_selection()
    out["H612_vib"] = _dec(sel["vib"])
    out["H612_sets"] = f"IR {_dec(sel['ir'])} | Raman {_dec(sel['raman'])} | silent {_dec(sel['silent'])}"
    out["H612_nu"] = " ".join(str(SO.SF6_NU[a]) for a in ("A_{1g}", "E_g", "T_{1u}+", "T_{1u}-", "T_{2g}", "T_{2u}"))
    C, cols = cg(ohset(), "T_{2u}", "E_g", ["T_{1u}"])
    out["H612_cg323"] = _ex(C.real.reshape(3, 2, 3)[2, 1, 2])
    # H613
    out["H613_pascal"] = " ".join(str(pascal(2, n)) for n in (1, 2, 3, 4)) + " | " + " ".join(str(pascal(3, n)) for n in (1, 2, 3, 4))
    out["H613_T2u"] = " | ".join(f"{n}: " + _dec(overtone("T_{2u}", n)) for n in (2, 3, 4))
    out["H613_Eg"] = " | ".join(f"{n}: " + _dec(overtone("E_g", n)) for n in (2, 3))
    out["H613_comb"] = _dec(combination("E_g", "T_{1u}"))
    # H614
    out["H614_sq"] = _dec(sym_squares(S_, "E_g")[0]) + " | " + _dec(sym_squares(S_, "T_{1g}")[0])
    rho0, e0 = jt_min(JT_CARD["j"], JT_CARD["k"])
    out["H614_moat"] = f"{_f4(rho0)} {_f4(e0)}"
    out["H614_sign"] = _fg(jt_eigvec_track(JT_CARD["j"])[2])
    rn = renner_extrema(JT_CARD["j"], JT_CARD["k"], JT_CARD["r"])
    out["H614_renner"] = " ".join(str(x) for x in rn["min_phi"]) + f" | {_f4(rn['rho_min'])} {_f4(rn['e_min'])} | " \
                         f"{_f4(rn['rho_saddle'])} {_f4(rn['e_saddle'])} | {_f4(rn['barrier'])}"
    lv0 = jt_vibronic_levels(JT_CARD["j"], JT_CARD["k"], 0.0, count=8)
    out["H614_levels"] = " ".join(_f4(e) for e, d, lab in lv0[:6])
    out["H614_l"] = " | ".join(" ".join(_f4(x) for x in jt_radial_levels(l_, JT_CARD["j"], JT_CARD["k"], 2, N=4800)) for l_ in (0.5, 1.5, 2.5))
    lvr = jt_vibronic_levels(JT_CARD["j"], JT_CARD["k"], JT_CARD["r"], count=6)
    out["H614_rlevels"] = " ".join(f"{lab} {_f4(e)}" for e, d, lab in lvr[:5])
    rm, vm = jt_effective_min(0.5)
    rb, vb = jt_effective_min(0.5, book_fig=True)
    out["H614_eff"] = f"{rm:.3f} {_f4(vm)} | {rb:.3f} {_f4(vb)}"
    out["H614_book_l"] = " ".join(_f4(x) for x in jt_radial_levels(0.5, JT_CARD["j"], JT_CARD["k"], 1, N=4800, book=True))
    # 習題
    p61 = problem_61()
    out["EX61"] = " || ".join(f"{k}: " + " / ".join(" ".join(str(x) for x in p61[k][0].row(i)) for i in range(4)) for k in "abcd")
    p62 = problem_62()
    out["EX62"] = " ".join(f"{k}:{p62[k]}" for k in ("z", "x2+y2", "z2", "HH", "HH2", "z(x2+y2)")) + f" | {p62['tensor9']} | {p62['poly2']}"
    p63 = problem_63()
    out["EX63"] = " ".join(f"{k}:{'+'.join(p63[k])}" for k in ("E1_x", "M1_z", "E2_xy")) + " | " + \
        f"{_f4(p63['omega']['A1'])} {_f4(p63['omega']['A2'])} {_f4(p63['omega']['E'][0])} {_f4(p63['omega']['E'][1])}"
    p64 = problem_64()
    out["EX64"] = _dec(p64["modes"]) + f" | {p64['E1_x']} {p64['M1_z']} {p64['E2_xy']}"
    p65 = problem_65()
    out["EX65"] = " ".join(f"{a}>{b}@{c}:{v}" for (a, b, c), v in p65.items() if v)
    p66 = problem_66()
    out["EX66"] = f"{_dec(p66['polar'])} | {_dec(p66['axial'])} | X6 {_dec(p66['X6'])} | C6H6 {_dec(p66['benzene'])} | " \
                  f"{p66['X6_ir']} {p66['X6_raman']} {p66['bz_ir']} {p66['bz_raman']} {p66['bz_block']}"
    return out


CARD_EXPECTED = {
    'H601_count': '36 2304',
    'H601_R1': '(1, 6)',
    'H601_kron': '0 0 0 0 1/2 sqrt(3)/2 / 0 0 0 0 sqrt(3)/2 -1/2 / -1/2 -sqrt(3)/2 0 0 0 0 / -sqrt(3)/2 1/2 0 0 0 0 / 0 0 1/2 sqrt(3)/2 0 0 / 0 0 sqrt(3)/2 -1/2 0 0',
    'H601_chi': '0 9',
    'H601_single': '-4 0 2 | -1 0 0.5',
    'H601_levels': '2.5(4) 2(6) 1(2) 0.5(6) 0(9) -1(3) -3.5(2) -4(3) -5(1)',
    'H602_V': '1 0.5 0.25',
    'H602_chiprod': '9 0 1 1 1 9 0 1 1 1',
    'H602_freq': 'A_{1g} + E_g + T_{1g} + T_{2g}',
    'H602_N': '0 0 0 0 0 1/2',
    'H602_cg': 'sqrt(3)/3 sqrt(6)/6 -sqrt(2)/2 0 0 0 0 0 0 / 0 0 0 0 0 sqrt(2)/2 0 0 sqrt(2)/2 / 0 0 0 0 -sqrt(2)/2 0 0 sqrt(2)/2 0 / 0 0 0 0 0 -sqrt(2)/2 0 0 sqrt(2)/2 / sqrt(3)/3 sqrt(6)/6 sqrt(2)/2 0 0 0 0 0 0 / 0 0 0 sqrt(2)/2 0 0 sqrt(2)/2 0 0 / 0 0 0 0 sqrt(2)/2 0 0 sqrt(2)/2 0 / 0 0 0 -sqrt(2)/2 0 0 sqrt(2)/2 0 0 / sqrt(3)/3 -sqrt(6)/3 0 0 0 0 0 0 0',
    'H602_split': '0.7250 0.5750 0.6000 0.4000 | 0.6250 0.6250 0.5000 0.5000',
    'H603_d2': 'A_1 A_2 B_1 B_2 / A_2 A_1 B_2 B_1 / B_1 B_2 A_1 A_2 / B_2 B_1 A_2 A_1',
    'H603_623': 'sqrt(2)/2 sqrt(2)/2 0 0 / 0 0 sqrt(2)/2 sqrt(2)/2 / 0 0 -sqrt(2)/2 sqrt(2)/2 / sqrt(2)/2 -sqrt(2)/2 0 0',
    'H603_626': '0 -1 / 1 0 | 1 0 / 0 -1 | 0 1 / 1 0',
    'H603_rel': 'b=sqrt(2)*a/2 d=-sqrt(2)*c e=-sqrt(3)*c g=-f h=sqrt(2)*f j=i k=sqrt(2)*i',
    'H603_val': 'a=sqrt(6)/3 b=sqrt(3)/3 c=sqrt(3)/3 d=-sqrt(6)/3 e=-1 f=sqrt(2)/2 g=-sqrt(2)/2 h=1 i=sqrt(2)/2 j=sqrt(2)/2 k=1',
    'H603_bookd': 'False',
    'H603_635a': 'sqrt(2)/2 0 sqrt(2)/2 0 / 0 sqrt(2)/2 0 -sqrt(2)/2 / 0 -sqrt(2)/2 0 -sqrt(2)/2 / sqrt(2)/2 0 -sqrt(2)/2 0',
    'H603_635b': '0 -1 / 1 0',
    'H603_636': 'sqrt(3)/3 sqrt(6)/6 0 sqrt(6)/6 sqrt(3)/3 0 0 0 0 / 0 0 -sqrt(6)/6 0 0 -sqrt(3)/3 0 0 sqrt(2)/2 / 0 0 -sqrt(3)/3 0 0 sqrt(6)/6 0 -sqrt(2)/2 0 / 0 0 -sqrt(6)/6 0 0 -sqrt(3)/3 0 0 -sqrt(2)/2 / sqrt(3)/3 -sqrt(6)/6 0 sqrt(6)/6 -sqrt(3)/3 0 0 0 0 / 0 sqrt(3)/3 0 0 -sqrt(6)/6 0 sqrt(2)/2 0 0 / 0 0 -sqrt(3)/3 0 0 sqrt(6)/6 0 sqrt(2)/2 0 / 0 sqrt(3)/3 0 0 -sqrt(6)/6 0 -sqrt(2)/2 0 0 / sqrt(3)/3 0 0 -sqrt(6)/3 0 0 0 0 0',
    'H604_scalar': 'sqrt(3)/3',
    'H604_fs': 'O/tetragonal: 1 1 1 1 1 | C3: 1 0 0 | T: 1 0 0 1 | Q: -1',
    'H604_c3': '0_3',
    'H604_sq': 'E: A_1 + E / A_2 | T_1: A_1 + E + T_2 / T_1 | T_2: A_1 + E + T_2 / T_1',
    'H604_TxT': 'A + ^1E + ^2E + 2T',
    'H605_o3': '0^+: A_{1g} | 1^+: T_{1g} | 2^+: E_g + T_{2g}',
    'H605_poly': 'sqrt(3)*x**2/3 + sqrt(3)*y**2/3 + sqrt(3)*z**2/3 | sqrt(6)*x**2/6 + sqrt(6)*y**2/6 - sqrt(6)*z**2/3 | -sqrt(2)*x**2/2 + sqrt(2)*y**2/2 | sqrt(2)*y*z | sqrt(2)*x*z | sqrt(2)*x*y',
    'H605_R3': '1 0 / 0 -1',
    'H606_stress': 'pressure: 1.7321 | uniaxial: 0.5774 -0.8165 | shear12: 1.4142 | torque: 1.4142',
    'H607_counts': '81 36 21 | Oh 4 3 3 | iso 2 | D4h 7 6 | C3v 8 6',
    'H607_spring111': '9 3 4 | 3 2 3.8571 0.2857 -0.5 5 2 2',
    'H607_spring211': '10 4 4 | 3.3333 2 5 0.25 0 6 2 2',
    'H607_kittel': 'diamond -0.35 Si -0.73 Ge -0.79 NaCl 0.88 Au -1.30 Pb -1.52 W -0.02',
    'H608_reps': 'E: A_1 + E | B: A_2 + E',
    'H608_chars': '3 0 1 | 3 0 -1',
    'H608_r': '-0.5000 0.8660 | -0.5000 0.8660',
    'H608_s3': '1 -1 1 -1 1 -1',
    'H608_P12': '0 -1 0 | 0 -3 | 0 -1',
    'H608_strict': 'C3v 4 2 | Oh 0 0',
    'H609_inv': '1/6 1/6 1/6 1/6 1/6 1/6 | 1/2 -1/2 0 0 0 0 | 1/3 1/3 -1/6 -1/6 -1/6 -1/6',
    'H609_levels': '-4.5 0.5 1.5 | 0 0.5 1',
    'H609_param': '3 7 11',
    'H609_multi': 'A_{1u} + A_{2u} + 2E_u + 4T_{1u} + 3T_{2u} | A_{2u} + 2T_{1u} + T_{2u}',
    'H609_mel': '0.8165 -0.4082 0.7071 -0.7071 0.0000 0.5774',
    'H609_ratio': '-0.5000 -0.5000 1.0000 0.8660 -0.8660 0.0000',
    'H609_red': '0.8165 0.5774 | 1.0000 0.6667',
    'H609_int': '0.2500 0.7500 | 1.0000 0.0000 | 0.5000',
    'H610_w': '31.4159 31.4119',
    'H610_rho': '4.05 8.57 38.23 90.00 161.94 175.04',
    'H610_R': '0.0561 0.1055 0.4020 0.6366 0.1880 0.0459',
    'H610_res': '0.03183',
    'H610_smith': '10.0000 | 10.03 0.995 | 5.56 0.894 | -4.55 1.095 | 2.0000',
    'H611_eps': '2.0000 1.4142 -1.2727 0.6667 0.8165 | 1+5i 1+2.5i',
    'H611_alpha': '7.95e-30 6.6e-30 -4.2e-33',
    'H611_gamma': '4.4e+07 3.4e-08 2.0e+07 | 3.4e-04 2.0e+09',
    'H611_cross': '2.81e-15 2.5e-29',
    'H611_rayleigh': '5.86 6.17',
    'H611_phasor': '1.0198 11.31 | 0.8000 0.00 | 1.0198 -11.31 | 1.2000 0.00',
    'H612_vib': 'A_{1g} + E_g + T_{2g} + 2T_{1u} + T_{2u}',
    'H612_sets': 'IR 2T_{1u} | Raman A_{1g} + E_g + T_{2g} | silent T_{2u}',
    'H612_nu': '775 644 965 617 524 363',
    'H612_cg323': '-1',
    'H613_pascal': '2 3 4 5 | 3 6 10 15',
    'H613_T2u': '2: A_{1g} + E_g + T_{2g} | 3: A_{1u} + T_{1u} + 2T_{2u} | 4: 2A_{1g} + 2E_g + T_{1g} + 2T_{2g}',
    'H613_Eg': '2: A_{1g} + E_g | 3: A_{1g} + A_{2g} + E_g',
    'H613_comb': 'T_{1u} + T_{2u}',
    'H614_sq': 'A_{1g} + E_g | A_{1g} + E_g + T_{2g}',
    'H614_moat': '0.5000 -0.2500',
    'H614_sign': '-1',
    'H614_renner': '60 180 300 | 0.6250 -0.3125 | 0.4167 -0.2083 | 0.1042',
    'H614_levels': '0.9758 2.0929 2.6872 3.2661 3.9258 4.1785',
    'H614_l': '0.9758 2.6872 | 2.0929 4.1785 | 3.2661 5.5522',
    'H614_rlevels': 'E 0.9478 A_1 1.7963 A_2 2.3147 E 2.5684 E 3.1513',
    'H614_eff': '0.874 0.2172 | 1.000 0.5000',
    'H614_book_l': '0.7738',
    'EX61': 'a: sqrt(2)/2 0 sqrt(2)/2 0 / 0 sqrt(2)/2 0 -sqrt(2)/2 / 0 -sqrt(2)/2 0 -sqrt(2)/2 / sqrt(2)/2 0 -sqrt(2)/2 0 || b: 0 0 0 -1 / sqrt(2)/2 sqrt(2)/2 0 0 / sqrt(2)/2 -sqrt(2)/2 0 0 / 0 0 1 0 || c: sqrt(2)/2 sqrt(2)/2 0 0 / 0 0 sqrt(2)/2 sqrt(2)/2 / 0 0 -sqrt(2)/2 sqrt(2)/2 / sqrt(2)/2 -sqrt(2)/2 0 0 || d: 0 sqrt(2)/2 0 sqrt(2)/2 / sqrt(2)/2 0 sqrt(2)/2 0 / sqrt(2)/2 0 -sqrt(2)/2 0 / 0 sqrt(2)/2 0 -sqrt(2)/2',
    'EX62': 'z:A_2 x2+y2:A_1 z2:A_1 HH:A_1 HH2:A_2 z(x2+y2):A_2 | 2A_1 + A_2 + 3E | 2A_1 + 2E',
    'EX63': 'E1_x:E M1_z:A_2 E2_xy:E | 4.8990 2.0000 4.0000 2.4495',
    'EX64': 'A1 + A2 + B1 + B2 + 2E | E A2 B2',
    'EX65': 'A1g>v1@v1:1/3 v1>Eg@v1:1/6 v1>Eg@v2:1/6 v1>Eg@v3:1/3 v2>Eg@v1:1/6 A1g>v2@v2:1/3 v2>Eg@v2:1/6 v2>Eg@v3:1/3 v3>Eg@v1:1/3 v3>Eg@v2:1/3 A1g>v3@v3:1/3',
    'EX66': 'A_{2u} + E_{1u} | A_{2g} + E_{1g} | X6 A_{1g} + B_{2g} + 2E_{2g} + B_{1u} + B_{2u} + E_{1u} + E_{2u} | C6H6 2A_{1g} + A_{2g} + 2B_{2g} + E_{1g} + 4E_{2g} + A_{2u} + 2B_{1u} + 2B_{2u} + 3E_{1u} + 2E_{2u} | 1 3 4 7 4',
}


# ======================================================================== 自我測試（每個結果兩種以上獨立算法）
KRON_617X_FULL_BOOK = np.array([[0, 0, 0, 0, .5, 3 ** .5 / 2], [0, 0, 0, 0, 3 ** .5 / 2, -.5], [-.5, -3 ** .5 / 2, 0, 0, 0, 0],
                                [-3 ** .5 / 2, .5, 0, 0, 0, 0], [0, 0, .5, 3 ** .5 / 2, 0, 0], [0, 0, 3 ** .5 / 2, -.5, 0, 0]])  # (6.1.7ₓ) 右邊照印


def _check_61():
    S = ohset()
    f6 = S.decompose(np.array([np.trace(site_perm(g)) for g in S.elems]))
    assert f6 == {"A_{1g}": 1, "E_g": 1, "T_{1u}": 1}
    assert outer_frequencies() == {(a, b): 1 for a in f6 for b in f6}             # 2304 個 (a,b) 的跡 ＝ 頻率相乘
    assert outer_trace_check()
    ke = kron_example()
    assert _close(ke["T1"], KRON_617X_BOOK["T1"]) and _close(ke["E"], KRON_617X_BOOK["E"])
    assert _close(ke["kron"], KRON_617X_FULL_BOOK)
    for H_, SA, SB in ((0, 1, 0.25), (0.3, 0.7, 0.2)):
        lev = two_particle_levels(SA, SB, H_)
        flat = sorted([e for e, a, b, d in lev for _ in range(d)], reverse=True)
        assert _close(flat, two_particle_numeric(SA, SB, H_), 1e-9)
    H6 = six_hamiltonian(0.2, 1.0, 0.3)
    lv = six_levels(0.2, 1.0, 0.3)
    ev = np.linalg.eigvalsh(H6)
    assert _close(sorted(ev), sorted([lv["A_{1g}"]] + [lv["T_{1u}"]] * 3 + [lv["E_g"]] * 2))
    assert _close(SIX_STATES[("T_{1u}", 0)] @ H6 @ SIX_STATES[("T_{1u}", 0)], lv["T_{1u}"])


def _check_62():
    sizes = [len(c) for c in SO.OH.classes]
    reps = [c[0] for c in SO.OH.classes]
    for a, b in (("T_{1u}", "T_{1u}"), ("E_g", "T_{1u}"), ("E_g", "E_g")):
        ia, ib = SO.OH_IRREP_NAMES.index(a), SO.OH_IRREP_NAMES.index(b)
        cls = {}
        for r, gname in enumerate(SO.OH_IRREP_NAMES):
            f = sum(sizes[c] * float(SO.OH_CHARS[r, c] * SO.OH_CHARS[ia, c] * SO.OH_CHARS[ib, c]) for c in range(10)) / 48
            if abs(f) > 1e-9:
                cls[gname] = int(round(f))
        assert cls == diag_frequencies(a, b), (a, b)
    assert reps[0] == "1"
    assert diag_frequencies("E_g", "E_g") == {"A_{1g}": 1, "A_{2g}": 1, "E_g": 1}
    assert diag_frequencies("E_g", "T_{1u}") == {"T_{1u}": 1, "T_{2u}": 1}
    C, cols = cg_book("619c")
    assert compare_table(C, cols, CG_619C_BOOK)["equal"]
    assert cg_reduces(ohset("nonkosher"), "T_{1u}", "T_{1u}", C, cols)
    Rk = np.kron(_dnum("T_{1u}", "R1", "nonkosher"), _dnum("T_{1u}", "R1", "nonkosher"))
    blk = C.conj().T @ Rk @ C
    assert _close(blk[6:, 6:], np.array([[-1, 0, 0], [0, 0, 1], [0, -1, 0]]))       # (6.1.19b) 右下的非 kosher T2g(R1)
    nn = dict(cg_norms(ohset("nonkosher"), "T_{1u}", "T_{1u}", "T_{1g}"))
    assert all(abs(nn[k] - v) < 1e-12 for k, v in CG_617_N_BOOK.items())
    a1 = t1u_split()
    p = PAIR_V
    form = {"A_{1g}": (p["v0"] + p["v2"]) / 2 + 2 * p["W"], "E_g": (p["v0"] + p["v2"]) / 2 - p["W"],
            "T_{2g}": p["v1"] + p["K"], "T_{1g}": p["v1"] - p["K"]}
    assert all(abs(a1[k] - form[k]) < 1e-12 for k in form)
    expect = sorted([a1["A_{1g}"]] + [a1["E_g"]] * 2 + [a1["T_{1g}"]] * 3 + [a1["T_{2g}"]] * 3)
    assert _close(t1u_split_numeric(), expect, 1e-12)
    assert _close(t1u_split_numeric(lam=1e-6), expect, 1e-4)
    assert pair_V_invariant() == (True, True)


def _check_63():
    for name, book in (("619c", CG_619C_BOOK), ("630", CG_630_BOOK), ("623", CG_623_BOOK), ("635a", CG_635A_BOOK),
                       ("635b", CG_635B_BOOK), ("636", CG_636_BOOK)):
        C, cols = cg_book(name)
        assert compare_table(C, cols, book)["equal"], name
    for x, (C, cols) in cg_book("626").items():
        assert compare_table(C, cols, CG_626_BOOK[x])["equal"], x
    C, cols = cg(oset("trigonal"), "T_1", "T_1", ["A_1", "E", "T_2", "T_1"])     # sym_oct 的 T2 基底：只有 T2 的 E 那兩行差號
    diff = np.abs(np.abs(C.real) - np.abs(CG_636_BOOK)).max()
    bad = [k for k in range(9) if not np.allclose(C[:, k].real, CG_636_BOOK[:, k], atol=1e-9)]
    assert diff < 1e-9 and bad == [4, 5]
    tr = oset("trigonal")
    assert "A_1" not in tr.decompose(tr.chi("T_2") * tr.chi("T_1"))                # T2⊗T1 沒有 A1 ⇒ (6.1.36) 不是 T2⊗T1
    assert tr.decompose(tr.chi("T_1") ** 2) == {"A_1": 1, "E": 1, "T_1": 1, "T_2": 1}
    assert d2_table() == D2_TABLE_BOOK
    S2 = d2set()
    for a in S2.names:                                                             # 第二種算法：一維表示逐元素相乘
        for b in S2.names:
            prod = S2.chi(a) * S2.chi(b)
            assert _close(prod, S2.chi(D2_TABLE_BOOK[(a, b)]))
    so = solve_6129()
    assert so["consistent"] and so["equals_630"] and not solve_6129(True)["consistent"]
    assert sp.simplify(so["relations"]["d"] + SQ2 * sp.Symbol("c", positive=True)) == 0


def _check_64():
    for S_ in (oset(), c3set(), tset(), qset(), c3vset(), d4set(), d3set()):
        for a in S_.names:
            fs = fs_indicator(S_, a)
            assert {1: 1, 0: 2, -1: 3}[fs] == irrep_type(S_, a), (S_.label, a)
    O_ = oset()
    for a in O_.names:
        sym, anti = sym_squares(O_, a)
        cs = cg_symmetry(O_, a)
        assert {k[0] for k, v in cs.items() if v == 1} == set(sym) and {k[0] for k, v in cs.items() if v == -1} == set(anti)
    S = ohset()
    sym, anti = sym_squares(S, "T_{1u}")
    assert sym == {"A_{1g}": 1, "E_g": 1, "T_{2g}": 1} and anti == {"T_{1g}": 1}
    Tt = tset()
    assert Tt.decompose(Tt.chi("T") ** 2) == {"A": 1, "^1E": 1, "^2E": 1, "T": 2}
    C, cols = cg(Tt, "T", "T")
    assert len([c for c in cols if c[0] == "T"]) == 6 and cg_reduces(Tt, "T", "T", C, cols)
    assert _close(C.conj().T @ C, np.eye(9))
    for g in Tt.elems:                                                             # ¹E 是 T 的表示（同態）
        for h in Tt.elems:
            gh = SO.O.table[g, h]
            assert _close(Tt.D["^1E"][g] @ Tt.D["^1E"][h], Tt.D["^1E"][gh])
    assert scalar_cg_check(oset()) and scalar_cg_check(c3set()) and scalar_cg_check(tset()) and scalar_cg_check(ohset())
    assert unitarity_check(O_, "T_1", "T_1") and unitarity_check(O_, "T_2", "E") and unitarity_check(ohset("nonkosher"), "T_{1u}", "T_{2g}")
    assert factorization_lemma(O_, "T_1", "T_1", "E") and factorization_lemma(O_, "T_1", "T_2", "T_1")
    assert factorization_lemma(Tt, "T", "T", "T")                                  # 重複 irrep：Σ_ω 兩份
    C1, _ = cg_book("619c")
    assert abs(C1[0, 0] - 1 / math.sqrt(3)) < 1e-12                                # (6.2.5) 的 1/√ℓ


def _check_65():
    assert tensor_transform_check() and poly_transform_check()
    pt = poly_from_tensors()
    for k, v in POLY_639.items():
        assert sp.simplify(pt[k] - v) == 0, k
    assert all(pt[("T_{1g}", m)] == 0 for m in range(3))
    U = unit_tensors()
    keys = list(U)
    G = np.array([[np.sum(U[a] * U[b]) for b in keys] for a in keys])
    assert _close(G, np.eye(9))
    assert o3_labels()["ok"]
    x, y, z = XS
    for g in ("R3", "R1", "r1", "i4"):                                             # x² − y² 屬 Eg：逐個操作代入
        M = sp.Matrix(SO.O.mats[g])
        sub = dict(zip(XS, list(M.T * sp.Matrix(XS))))
        img = sp.expand((x ** 2 - y ** 2).subs(sub, simultaneous=True))
        Dm = _dnum("E_g", g, "tetragonal").real
        e1, e2 = (2 * z ** 2 - x ** 2 - y ** 2) / SQ3, x ** 2 - y ** 2            # √2·(E 的兩個二次式基底)
        assert sp.simplify(img - (exact(Dm[0, 1]) * e1 + exact(Dm[1, 1]) * e2)) == 0, g
    T43 = np.array(SO.T43.evalf(), dtype=float)                                     # (6.3.3b)：v3 屬 D3d 的 A2u
    for g in SO.D3D:
        M = SO.OH.num[g].real
        img = M @ T43[:, 2]
        sgn = img @ T43[:, 2]
        assert abs(abs(sgn) - 1) < 1e-9
        inv = g.startswith("I")
        c2 = (not inv and g in ("i2", "i4", "i5")) or (inv and g in ("I.r1", "I.r1^2", "I"))
        assert (sgn < 0) == c2, g


def _check_66():
    st = stress_examples()
    assert abs(st["pressure"][("A_{1g}", 0)] - math.sqrt(3)) < 1e-12
    assert abs(st["uniaxial"][("A_{1g}", 0)] - 1 / math.sqrt(3)) < 1e-12 and abs(st["uniaxial"][("E_g", 0)] + 2 / math.sqrt(6)) < 1e-12
    assert set(st["shear12"]) == {("T_{2g}", 2)} and set(st["torque"]) == {("T_{1g}", 0)}
    rng = np.random.default_rng(66)
    M = rng.normal(size=(3, 3))
    assert _close(tensor_from_components(tensor_components(M)), M)
    Ssym = rng.normal(size=(3, 3))
    Ssym = Ssym + Ssym.T
    assert _close(strain_of_displacement(lambda r: Ssym @ r), Ssym, 1e-7)
    assert all(abs(v) < 1e-9 for k, v in tensor_components(Ssym).items() if k[0] == "T_{1g}")
    assert _close(face_forces(M)[1], M[:, 1])


def _check_67():
    ec, ch = elastic_counts(), elastic_counts_chars()
    assert ec["none"] == {"81": 81, "36": 36, "21": 21} and ec["Oh"] == {"81": 4, "36": 3, "21": 3}
    assert (ec["Oh"]["36"], ec["Oh"]["21"]) == ch["Oh"][:2]
    assert (ec["D4h"]["36"], ec["D4h"]["21"]) == ch["D4h"] == (7, 6)
    assert (ec["C3v"]["36"], ec["C3v"]["21"]) == ch["C3v"] == (8, 6)
    assert ec["iso"]["21"] == 2
    S = ohset()
    assert S.decompose(S.chi("T_{1u}") ** 4).get("A_{1g}") == ec["Oh"]["81"]
    for kA, kE, kT in ((9.0, 3.0, 4.0), (7.3, 2.1, 5.9)):
        K4 = elastic_tensor(kA, kE, kT)
        r = elastic_from_cartesian(K4)
        assert abs(r["kA"] - kA) + abs(r["kE"] - kE) + abs(r["kT2"] - kT) < 1e-9
        assert abs(r["kA_x"] - kA) + abs(r["kE_x"] - kE) + abs(r["kT2_x"] - kT) < 1e-9   # 更正的 (6.3.28)ₓ
        assert abs(K4[0, 0, 0, 1]) < 1e-12 and abs((K4[0, 0, 0, 0] + 2 * K4[0, 0, 0, 1]) - kA) > 1                 # 照印 k_{11,12}：不對
        assert abs(K4[0, 0, 0, 0] - (kA / 3 + 2 * kE / 3)) < 1e-12 and abs(K4[0, 0, 1, 1] - (kA - kE) / 3) < 1e-12
        assert abs(K4[0, 1, 0, 1] - kT / 2) < 1e-12 and abs(K4[0, 0, 1, 1] - K4[1, 1, 0, 0]) < 1e-12
        for g in ("R1", "r1", "I.i3"):
            Rm = SO.OH.num[g].real
            assert _close(np.einsum("ai,bj,ck,dl,ijkl->abcd", Rm, Rm, Rm, Rm, K4), K4)
        m, mn = moduli(kA, kE, kT), moduli_numeric(kA, kE, kT)
        assert all(abs(float(m[k]) - mn[k]) < 1e-8 for k in ("B", "mu", "Y", "sigma"))
        sig = float(m["sigma"])
        assert abs(float(m["Y"]) / (1 - 2 * sig) - kA) < 1e-9 and abs(float(m["Y"]) / (1 + sig) - kE) < 1e-9     # (6.3.37)
    for jkl in MODULI_CARD.values():
        f = spring_model(*jkl)
        b = spring_model_bonds(*map(float, jkl))
        assert abs(b["kA"] - f[0]) + abs(b["kE"] - f[1]) + abs(b["kT2"] - f[2]) < 1e-9
        assert abs(b["k1122"] - b["k1212"]) < 1e-12                                # Cauchy 關係 (6.3.42)
    m = moduli(*spring_model(2, 1, 1))
    assert m["A"] == 0 and m["sigma"] == Fraction(1, 4) and m["B"] == m["Y"] / (3 * (1 - 2 * m["sigma"])) and m["mu"] == m["Y"] / (2 * (1 + m["sigma"]))
    K4 = elastic_tensor(10, 4, 4)
    Rr = SN_rot3((1, 2, 3), 37)
    assert _close(np.einsum("ai,bj,ck,dl,ijkl->abcd", Rr, Rr, Rr, Rr, K4), K4)      # k^E = k^{T2}：任意轉動都不變（等向）
    K4b = elastic_tensor(10, 3, 4)
    assert not _close(np.einsum("ai,bj,ck,dl,ijkl->abcd", Rr, Rr, Rr, Rr, K4b), K4b)
    for k, A in FIG639_A_BOOK.items():
        assert abs(anisotropy(*CUBIC_KITTEL[k]) - A) < 0.06, k


def _check_68():
    fi = c3v_field_irreps()
    assert fi == {"E": {"A_1": 1, "E": 1}, "B": {"A_2": 1, "E": 1}}
    assert _close(p3_12("B"), [0, -1, 0]) and _close(p3_12("E"), [0, 0, 0]) and _close(p3_12("E", (0, 1, 0)), [1, 0, 0])
    S = c3vset()
    p22 = sum((2 / 6) * S.D["E"][g][1, 1].real * (c3v_field(g, "B") @ np.array([1.0, 0, 0])) for g in S.elems)
    assert _close(p22, [1, 0, 0])                                                   # P^3_22 B̂1 = B̂1
    basis = [-np.array([0, 1.0, 0]), np.array([1.0, 0, 0])]                         # (6.3.48)：(B̂^E_1, B̂^E_2) = (−B̂2, B̂1)
    for g in S.elems:
        for m in range(2):
            img = c3v_field(g, "B") @ basis[m]
            assert _close(img, sum(S.D["E"][g][mp, m].real * basis[mp] for mp in range(2)))
    assert p3_12_book() == (0, -3) and p3_12_book(SQ3 / 2) == (0, -1)
    f6 = fig6311()
    assert _close(f6[("E", "r")][0], [-0.5, math.sqrt(3) / 2, 0]) and _close(f6[("E", "s3")][1], [0, -1, 0])
    assert _close(f6[("B", "s3")][0], [-1, 0, 0]) and _close(f6[("B", "s3")][1], [0, 1, 0]) and _close(f6[("B", "s3")][2], [0, 0, -1])
    sc, scc = striction_counts(), striction_counts_chars()
    assert (sc[("C3v", "E")], sc[("C3v", "B")]) == (scc["E"], scc["B"]) == (4, 2)
    assert sc[("Oh", "E")] == sc[("Oh", "B")] == 0
    assert scc["Sym2"] == {"A_1": 2, "E": 2}


def _check_69():
    inv = six_invariants()
    for k, f in (("A", 1), ("T", math.sqrt(3)), ("E", math.sqrt(2))):
        assert _close(np.abs(inv[k][0] * f), SIX_INV_ROWS_BOOK[k])                  # 絕對值同書，號不同
    assert inv["T"][0, 1] < 0 and inv["E"][0, 2] < 0
    assert _close(sum(inv[k][0] * f for k, f in (("A", 1), ("T", math.sqrt(3)), ("E", math.sqrt(2)))), [1, 0, 0, 0, 0, 0])
    for H_, S_, T_ in ((0, 1, 0.5), (0.4, -0.7, 0.2)):
        lv = six_levels(H_, S_, T_)
        H6 = six_hamiltonian(H_, S_, T_)
        assert _close(lv["A_{1g}"] * inv["A"] + lv["T_{1u}"] * math.sqrt(3) * inv["T"] + lv["E_g"] * math.sqrt(2) * inv["E"], H6)
        b = six_from_levels(lv["A_{1g}"], lv["T_{1u}"], lv["E_g"])
        assert abs(b["H"] - H_) + abs(b["T"] - T_) + abs(b["S"] - S_) < 1e-12
    assert param_count({"A": 1, "T": 1, "E": 1}) == commutant_count([site_perm(g) for g in ("R1", "r1", "I")]) == 3
    xp = xf6_param_counts()
    assert xp["full"] == xp["full_numeric"] == 11 and xp["vib"] == xp["vib_numeric"] == 7
    assert commutant_count([np.eye(6)]) == 21
    w = we_check()
    assert w["E_T"]["ok"] and w["A_T"]["ok"] and w["T_E"]["ok"] and not w["T_E"]["ok_book"] and w["offdiag_zero"]
    assert all(abs(w["ratios"][k] - v) < 1e-12 for k, v in WE_617_BOOK.items())
    wc = we_complex_check()
    assert wc["complex_C"] and wc["ok"] and not wc["ok_noconj"] and not wc["ok_book"]
    wi = we_intensities()
    assert _close(wi["pol1"], (0.25, 0.75)) and _close(wi["pol3"], (1, 0)) and abs(wi["AE"] - 0.5) < 1e-12
    S = ohset()
    assert S.decompose(S.chi("T_{1u}") * S.chi("A_{1g}") * S.chi("E_g")).get("A_{1g}", 0) == 0   # A1g ↔ Eg 偶極禁戒
    assert S.decompose(S.chi("T_{1g}") * S.chi("E_g")).get("A_{1g}", 0) == 0
    assert S.decompose(S.chi("T_{1u}") * S.chi("T_{1u}") * S.chi("T_{1u}")).get("A_{1g}", 0) == 0


def _check_6105():
    p = FIG651
    w0 = 2 * math.pi * p["f0"]
    ts = np.linspace(0, 10, 1201)
    for fs, rb in zip(p["fs"], FIG651_RHO_BOOK):
        ws = 2 * math.pi * fs
        R, rho = lorentz_steady(ws, w0, p["gamma"], LORENTZ_A)
        assert abs(rho - rb) < 0.06, fs
        xc = lorentz_closed(ts, ws, w0, p["gamma"], LORENTZ_A)
        xn = lorentz_numeric(ts, ws, w0, p["gamma"], LORENTZ_A)
        assert np.max(np.abs(xc - xn)) < 1e-7, fs
        tr = lorentz_transient(ws, w0, p["gamma"], LORENTZ_A)
        if "a_book" in tr:
            xb = R * np.cos(ws * ts - math.radians(rho)) + tr["a_book"] * np.exp(-p["gamma"] * ts) * np.cos(tr["wg"] * ts - math.radians(tr["alpha_book"]))
            assert np.max(np.abs(xb - xc)) < 1e-9
        z = response(ws, w0, p["gamma"], LORENTZ_A)
        assert abs(abs(z) - R) < 1e-12 and abs(math.degrees(np.angle(z)) - rho) < 1e-9
    tr = lorentz_transient(w0, w0, p["gamma"], LORENTZ_A)
    assert abs(tr["alpha"] + 90) < 1e-9 and "a_book" not in tr
    sc = smith_check()
    assert sc["on_gamma"] < 1e-12 and sc["on_delta"] < 1e-12 and abs(sc["R_cos"]) < 1e-12 and abs(sc["R_sin"]) < 1e-12 and sc["maxrel"] < 0.1
    f2 = fig652()
    assert abs(f2["w_re_max"] - 0.9) < 0.01 and abs(f2["w_re_min"] - 1.1) < 0.01 and abs(f2["im_max"] - 10) < 0.05 and f2["dc"] == 2.0
    assert abs(abs(smith(0.1)) ** 2 - 1 / (0.1 ** 2 + 0.1 ** 2)) < 1e-12              # Lorentz 形 (6.5.10b)


def _check_6106():
    for wv in (0.0, 0.5, 1.0, 1.3, 2.0):
        for g in (0.0, 0.1, 0.2):
            if wv == 1.0 and g == 0.0:
                continue
            e = epsilon(wv, 1.0, 1.0, g)
            assert abs(e - (1 + 1 / (1 - wv ** 2 - 2j * g * wv))) < 1e-12
            eb = epsilon(wv, 1.0, 1.0, g, book=True)
            assert abs(eb.imag - e.imag / 2) < 1e-12
    f4 = fig654()
    assert abs(f4["opaque"][0] - 1) < 2e-3 and abs(f4["opaque"][1] - math.sqrt(2)) < 2e-3 and f4["eps0"] == 2.0
    be = book_examples()
    for k, v in (("alpha_e", 7.95e-30), ("alpha_v", 6.6e-30), ("G_vis", 4.4e7), ("t_vis", 3.4e-8), ("t_ir", 3.4e-4), ("sigma", 2.5e-29)):
        assert abs(be[k] / v - 1) < 0.015, k
    assert abs(be["alpha_v_opt"] / -4e-33 - 1) < 0.1 and abs(be["n_vis"] / 2e7 - 1) < 0.02 and abs(be["n_ir"] / 2e9 - 1) < 0.02
    assert abs(be["w_vis"] / 3.77e15 - 1) < 1e-3 and abs(be["w500"] / 9.4e13 - 1) < 3e-3 and abs(be["w136"] / 2e16 - 1) < 0.04
    assert abs(radiation_integral() - 4 / 3) < 1e-12
    q, m, c, e0 = BOOK_CONST["q"], BOOK_CONST["m"], BOOK_CONST["c"], BOOK_CONST["eps0"]
    w, w0, G, E0 = 3e15, 2e16, 1e8, 1.0                                             # (6.5.33) 的三個因子相乘 ＝ (6.5.32)
    P32 = q ** 2 * w ** 4 / (12 * math.pi * c ** 3 * e0) * abs(q * E0 / m / (w0 ** 2 - w ** 2 - 2j * G * w)) ** 2
    rcl = q ** 2 / (4 * math.pi * e0 * m * c ** 2)
    P33 = (0.5 * e0 * c * E0 ** 2) * (math.pi * rcl ** 2) * (8 * w ** 4 / (3 * abs(w0 ** 2 - w ** 2 - 2j * G * w) ** 2))
    assert abs(P32 / P33 - 1) < 1e-12
    ry = rayleigh(**RAYLEIGH_CARD)
    assert abs(ry["law4"] - (ry["wb"] / ry["wr"]) ** 4) < 1e-9 and ry["full"] > ry["law4"]
    ph = phasor_points()
    assert ph[0][2] > 0 and abs(ph[1][1] - 0.8) < 1e-12 and ph[2][2] < 0 and abs(ph[3][1] - 1.2) < 1e-12


def _check_6107():
    sel, num = xf6_selection(), xf6_selection_numeric()
    assert sel["vib"] == {"A_{1g}": 1, "E_g": 1, "T_{2g}": 1, "T_{1u}": 2, "T_{2u}": 1} and sel["exclusive"]
    for a, (dip, pol) in num.items():
        assert (dip > 1e-6) == (a == "T_{1u}"), a
        assert (pol > 1e-6) == (a in sel["polarizability"]), a
    assert set(sel["ir"]) == {"T_{1u}"} and set(sel["raman"]) == {"A_{1g}", "E_g", "T_{2g}"} and set(sel["silent"]) == {"T_{2u}"}
    assert sum(n * (1 if a.startswith("A") else (2 if a.startswith("E") else 3)) for a, n in sel["vib"].items()) == 15


def _check_6108():
    S = ohset()
    for a in ("T_{2u}", "E_g", "T_{1u}"):
        for n in (2, 3, 4):
            c1, c2 = sym_power_chi(S, a, n), sym_power_chi_monomial(S, a, n)
            assert _close(c1, c2, 1e-8), (a, n)
            assert abs(c1[0] - pascal(S.dim(a), n)) < 1e-9
    assert overtone("T_{2u}", 3) == {"A_{1u}": 1, "T_{1u}": 1, "T_{2u}": 2}
    assert combination("E_g", "T_{1u}") == {"T_{1u}": 1, "T_{2u}": 1}
    assert [pascal(3, n) for n in range(1, 5)] == [math.comb(n + 2, 2) for n in range(1, 5)]
    assert henon_heiles_check()


def _check_6109():
    j, k, r = JT_CARD["j"], JT_CARD["k"], JT_CARD["r"]
    rng = np.random.default_rng(67)
    for _ in range(20):
        x1, x2 = rng.normal(size=2)
        rho, ph = math.hypot(x1, x2), math.degrees(math.atan2(x2, x1))
        ep, em = jt_sheets(rho, ph, j, k, r)
        ev = np.linalg.eigvalsh(jt_matrix(x1, x2, j, k, r))
        assert _close(ev, [em, ep], 1e-10)
    rhos = np.linspace(0, 2, 200001)
    em = jt_sheets(rhos, 37.0, j, k, 0.0)[1]
    rho0, e0 = jt_min(j, k)
    assert abs(rhos[np.argmin(em)] - rho0) < 1e-4 and abs(em.min() - e0) < 1e-9
    phis, vecs, sign, ana = jt_eigvec_track(j)
    assert sign == -1 and _close(vecs, ana, 1e-9)
    rn = renner_extrema(j, k, r)
    assert rn["grid_min_phi"] == rn["min_phi"] == [60, 180, 300]
    assert abs(rn["grid_min"] - rn["e_min"]) < 1e-6 and abs(rn["grid_saddle"] - rn["e_saddle"]) < 1e-6
    assert renner_extrema(j, k, -r)["grid_min_phi"] == [0, 120, 240]
    lv = jt_vibronic_levels(j, k, 0.0, count=8)
    assert all(d == 2 for e, d, lab in lv)
    rad = sorted(np.concatenate([jt_radial_levels(l_, j, k, 3, N=4800) for l_ in (0.5, 1.5, 2.5, 3.5)]))
    assert _close([e for e, d, lab in lv[:6]], rad[:6], 2e-5)
    lb = jt_vibronic_levels(j / 2, k / 2, 0.0, count=2)[0][0]                       # 照印 (6.7.15) ＝ k、j 減半
    assert abs(lb - jt_radial_levels(0.5, j, k, 1, N=4800, book=True)[0]) < 2e-5
    lr = jt_vibronic_levels(j, k, r, count=6)
    lr2 = jt_vibronic_levels(j, k, -r, count=6)
    assert lr[0][2] == "E" and _close([e for e, d, lab in lr], [e for e, d, lab in lr2], 1e-9)
    assert [lab for e, d, lab in lr][:5] == ["E", "A_1", "A_2", "E", "E"] and [lab for e, d, lab in lr2][1:3] == ["A_2", "A_1"]
    rb, vb = jt_effective_min(0.5, book_fig=True)
    assert abs(rb - 1) < 1e-4 and abs(vb - 0.5) < 1e-9                              # (l²+¼)/ρ² + ρ² − ρ：ρ = 1 時導數 −1 + 2 − 1 = 0
    rm, vm = jt_effective_min(0.5)
    assert abs(-0.5 / rm ** 3 + 2 * rm - 1) < 1e-3


def _check_problems():
    p61 = problem_61()
    assert np.allclose(np.array(p61["a"][0].evalf(), dtype=float), CG_635A_BOOK) and np.allclose(np.array(p61["c"][0].evalf(), dtype=float), CG_623_BOOK)
    for key in ("b", "d"):
        Cm = np.array(p61[key][0].evalf(), dtype=complex)
        assert _close(Cm.conj().T @ Cm, np.eye(4))
    p62 = problem_62()
    assert (p62["z"], p62["x2+y2"], p62["z2"], p62["HH"], p62["HH2"]) == ("A_2", "A_1", "A_1", "A_1", "A_2")
    p63 = problem_63()
    assert p63["E1_x"] == ["E"] and p63["M1_z"] == ["A_2"] and p63["E2_xy"] == ["E"]
    p64 = problem_64()
    assert sum(n * (2 if a == "E" else 1) for a, n in p64["modes"].items()) == 8 and (p64["M1_z"], p64["E2_xy"]) == ("A2", "B2")
    p65 = problem_65()
    for b in ("v1", "v2", "v3"):
        assert sum(p65[(a, "Eg", b)] for a in ("v1", "v2", "v3")) == Fraction(2, 3)
        assert sum(p65[("A1g", a, b)] for a in ("v1", "v2", "v3")) == Fraction(1, 3)
    assert p65[("v3", "Eg", "v3")] == 0
    p66 = problem_66()
    assert p66["n_vib"] == (12, 30) and (p66["X6_ir"], p66["X6_raman"], p66["bz_ir"], p66["bz_raman"]) == (1, 3, 4, 7)
    G = d6h_group()                                                                 # D6h 特徵標表的正交性（第二種算法）
    for a in D6H_IRREPS:
        for b in D6H_IRREPS:
            v = sum(d6h_chi(a, key) * d6h_chi(b, key) for key, M in G) / 24
            assert abs(v - (1 if a == b else 0)) < 1e-9
    for key, M in G:
        assert abs(d6h_chi("E_{1u}", key) + d6h_chi("A_{2u}", key) - np.trace(M)) < 1e-9


def _check_cards():
    cv = card_values()
    for k, v in CARD_EXPECTED.items():
        assert cv[k] == v, (k, cv[k], v)


def summary() -> dict:
    return card_values()


def check():
    """每個上卡片／講義的數字用兩種以上獨立算法核對（約 15 秒）。"""
    _check_61()
    _check_62()
    _check_63()
    _check_64()
    _check_65()
    _check_66()
    _check_67()
    _check_68()
    _check_69()
    _check_6105()
    _check_6106()
    _check_6107()
    _check_6108()
    _check_6109()
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
    print(f"sym_prod（第 6 章）：八面體上兩粒子的外積群與對角子群、投影法 Clebsch–Gordan 係數（O_h、O⊃D4⊃D2、O⊃D3⊃C2 七張表逐格同書）、"
          f"(6.1.29) 的未知數、Frobenius–Schur 三型與 T 群的重複 irrep、單位張量與多項式、應力應變與彈性常數（81→36→21→3、彈簧模型、"
          f"等向、C3v）、C3v 的 E 與 B、伸縮常數、參數定理與 Wigner–Eckart、Lorentz 振子與 Smith 圓、介電函數與輻射、XF6 紅外／Raman、"
          f"泛音、E⊗e Jahn–Teller 與半整數振電能階、習題 6.1–6.6 全部核對通過（{time.time() - t0:.1f} 秒）。")
