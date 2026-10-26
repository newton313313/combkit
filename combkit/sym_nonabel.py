# 由 combkit/tools/sync_from_harter_sym.py 從 harter_sym 專案的 codes/utils/sym_nonabel.py 複製；請改 harter_sym 那份再同步。
r"""
sym_nonabel.py —— Harter《Principles of Symmetry, Dynamics, and Spectroscopy》（Wiley 1993）第 3 章
Basic Theory and Applications of Symmetry Representations（Non-Abelian Symmetry Groups）的數學
（純 Python＋numpy＋sympy＋scipy，不 import manim；借用同資料夾的 sym_abel（平面波轉子、單擺環、外積）與 sym_mat）。
書上的式號、筆誤與〔注〕見 ref/md/ch03_formulas.md（相對於 C:\mypython\manim\projects\harter_sym）。

約定
  * 群乘積 ab = 先 b 後 a（書 p. 152）；乘法表 table[a, b] = ab（列 a、行 b）。精確值用 sympy；卡片上的小數一律 f"{float(x):.4f}"
    （f-string 不能直接 :.4f 一個 sympy 數或 Fraction）。運動方程 ẍ = −𝒜x（本徵值 = ω²）。ħ = 1。
  * 元素的程式名稱（TeX 在 G.tex[名稱]）：C3V：1, r, r^2, s1, s2, s3（σ_j）；D3：1, r, r^2, p1, p2, p3（ρ_j）；
    D4：1, R^2, R, R^3, R1^2, R2^2, i3, i4（(3.6.1a) 的順序）；C4V：1, R^2, R, R^3, s1..s4；D6：1, h^2, h^4, p1, p2, p3, h^3, h^5, h,
    p1', p2', p3'（(3.6.16)）；D3H：D3 的名稱＋"sh"、"sh.r"…；D3D、D4H、D2H：＋"I"、"I.r"…；D2：1, Rz, Ry, Rx；QUAT：1, -1, i, -i, j, -j, k, -k。
  * 幾何（Fig. 3.1.2，鏡面固定在紙上）：σ₃ 的鏡線沿 x 軸（0°）、σ₂ 在 60°、σ₁ 在 120°（＝ −60°）；r ＝ 逆時針轉 120°；
    三個擺錘 |1⟩、|r⟩、|r²⟩ 在 0°、120°、240°（σ₃|1⟩ = |1⟩、σ₂|1⟩ = |r⟩、σ₁|1⟩ = |r²⟩）。𝒟^E(g) 就是 g 在 xy 平面上的 2×2 矩陣（(3.3.7)）。
    D₃：ρ_j = Iσ_j ＝ 繞水平軸（方向 φ_j + 90°）轉 180°（ρ₃ 繞 y 軸）。
  * C_n 的不可約表示一律 𝒟^{k_n}(r) = ε_n^k = e^{+2πik/n}（第 2 章、(3.3.18)、(3.4.24b) 的約定；P^{k_n} = (1/n)Σ_p ε_n^{−kp} r^p，
    r P^{k_n} = ε_n^k P^{k_n}）。書上 (3.4.24a) 寫 ε = e^{−2πi/3} 是筆誤；書上 §3.6 的 C₄ 用相反的號 𝒟^{1₄}(R) = −i，
    所以書上 (3.6.4b) 的 P^E_{1₄} 是本系列的 P^E_{3₄}。轉子的波 e^{imφ} 照第 2 章 H210 標成 m mod n（(rψ)(φ) = ψ(φ + 2π/n)，
    e^{imφ} 的本徵值 e^{2πim/n}）：m = 1 的行波叫 1₄（與書上 Fig. 3.6.3 同名）。
  * 實驗室與本體座標（§3.1C）：本體的「先 a 後 b」＝ 實驗室的 ab（順序反過來）。

有限群（任何由精確矩陣定義的群；examples/ex31 也用 from_generators）
  FiniteGroup(mats, tex, label)；FiniteGroup.from_generators({名: 矩陣})（元素名＝最短的字，"rs" = r·s，"r^2"）
    .names .order .e .table[a, b] .mul(*gs) .inv(g) .conj(h, g) = hgh⁻¹ .order_of(g) .table_rows(order) .is_abelian()
    .subgroup_of(g)（[1, g, g², …]）.is_subgroup(H) .is_normal(H) .cosets(H, "left"/"right") .double_cosets(H, K) .commutator_subgroup()
    .classes（依 names 中第一次出現）.class_of(g) .class_sizes() .normalizer(g)（°c = °G/°N，(3.2.7)）
    .struct[i, j]（結構常數 γ^k_{ij}）.class_mul(x, y) .class_matrix(i)（類代數的正則表示，習題 3.5.1）.class_alg_rows(symbols)
    .min_poly_class(i)（首一 Poly，變數 c）.class_idempotents_of(i) .class_idempotents（ℙ^α 的類係數向量）.central_idempotents（群代數）
    .char_table_class_algebra() → (ℓ 列表, 特徵標表)（§3.5A：ℓ = (°G p₁)^{1/2}、χ* = p°G/ℓ）
    群代數（dict {名: 係數}）：.ga .el(g) .ga_mul .ga_add(*xs, coeffs) .ga_scale .ga_dag（x† = Σ x* g⁻¹）.ga_eq .class_op(i) .from_class_vec
    .regular(g)（𝓡_{hf}(g) = δ_{g,hf⁻¹}，(3.4.3)）.regular_of(x) .vec_of(x)（|x⟩ = x|1⟩ 的分量）
    .cyclic_projectors(gen)（[P^{0_n}, P^{1_n}, …]）.split(ℙ, projs) .elementary(Ps, ref) → ({(i, j): P_{ij}}, ref)（0 起算；D_{1j}(ref) > 0）
    .irrep_from_elementary(Pij) → {g: Matrix} .irreps(sub_gen, ref)（全部不可約表示，依 class_idempotents 順序）.is_rep(D)
    .char_table_traces(irreps) .frequencies(chars)（(3.5.11)）.correlation(H, H 的特徵標表, 表)（相關表）.iso_to(other, mapping)
  ga_tex(G, x, factor=None, named=None)（群代數元素 → TeX，提出公因數；ℙ³ → \tfrac{1}{3}(2\mathbf 1-r-r^2)；named={EPS3: r"\varepsilon", …}）
  reorder(T, ref_rows)、book_order(G, ref_rows) → (ℓ, 表, perm)（換成書上的列順序）、outer_table(T1, T2) = kron（(3.6.17)）、cn_chars(n)
  幾何：cosd、sind、rot2(deg)、refl2(deg)（鏡線與 x 軸夾 deg）、rot3(axis, deg)、block3(M2)、e2(deg)（浮點單位向量）

§3.1 群（H301–H302）
  C3V（2×2 幾何矩陣）；C3V_TABLE_BOOK（(3.1.1)）＝ C3V.table_rows()：σ₁σ₂ = r、σ₂σ₁ = r²；MIRROR_DEG = {s1: 120, s2: 60, s3: 0}
  D3（3×3）≅ C3v（D3_TO_C3V）、(3.1.2) σ_j = Iρ_j（INV3 = −1）；C3V_TEX、C3V_IRREP_NAMES = [A_1, A_2, E]
  BOB_START（Fig. 3.1.2 的 a、b、c，環半徑 1）、config_after(g, pts) → {名: g·p}、fig312() → {"s2s1": (中間, 最後, r²·起點), "s1s2": (…, r·起點)}
  Fig. 3.1.4：refl2(φ)·refl2(0) = rot2(2φ)（φ = 60°：R(120°) = [[−1/2, −√3/2], [√3/2, −1/2]]）
  Hamilton 轉：quat(axis, deg)、qmul(p, q)（先 q）、quat_to_rot、su2、axis_angle、rotate_vec、reflection3(N)、
    hamilton_arc(axis, deg, start) → (N1, N2)（弧長 deg/2；R = σ_{N2}σ_{N1}）、arc_points(N1, N2, n, radius)、arc_quat(N1, N2)、
    hamilton_product(axis1, deg1, axis2, deg2) → dict(N1, P（交點）, N2p, axis, deg, q, q_arc)（先 axis1）：
      先 R_z(90°) 再 R_x(90°)：兩段 45° 弧接成 60° 弧 ＝ 轉 120°，軸 (1, −1, 1)/√3；順序反過來：軸 (1, 1, 1)/√3
    d3_arcs()（Fig. 3.1.7：r 的 60° 赤道弧、ρ_j 的 90° 弧，單位球）；D3_AXIS_DEG = {p1: 210, p2: 150, p3: 90}；
      ρ₁r = ρ₂（弧相接＝四元數）；SU(2) 裡 ρ₁ 軸取 210° 時 q(ρ₁)q(r) = +q(ρ₂)，ρ₁ 軸反向時 = −q(ρ₂)（Fig. 3.1.7 的 (−)）
    conjugate_rotation(T_axis, T_deg, R_axis, R_deg) → (軸, 角, T·n̂)（TRT⁻¹：同角、換軸；Fig. 3.2.1；ρ₂ρ₁ρ₂⁻¹ = ρ₃）
  body_product(G, first, second)（Fig. 3.1.8：本體 σ̲₂σ̲₁ = r̲，實驗室 σ₂σ₁ = r²）
  D4（3×3，(3.6.1a) 的順序）、D4_TABLE_BOOK、D4_2D（xy 平面的作用 ＝ (3.6.7)）、C4V ≅ D4（D4_TO_C4V）、D2D ≅ D4（(3.6.32)）
§3.2 類與類代數（H303–H304）
  C3V.classes = [[1], [r, r²], [σ₁, σ₂, σ₃]]；N(r) = C₃、N(σ₁) = {1, σ₁}；°c = 6/°N = 1, 2, 3；cosets：σ₁C₃ = {σ₁, σ₂, σ₃}、rC_v = {r, σ₃}、r²C_v = {r², σ₂}
  class_alg_rows ＝ C3V_CLASS_ALG_BOOK（(3.2.1b)）；c₃² = 3c₁ + 3c₂、c₃³ = 9c₃；min_poly_class(2) = c³ − 9c = (c − 3)(c + 3)c
  class_idempotents ＝ C3V_CLASS_IDEMP_BOOK：ℙ¹ = (c₁+c₂+c₃)/6、ℙ² = (c₁+c₂−c₃)/6、ℙ³ = (2c₁−c₂)/3（(3.2.15)）
  c3v_class_spectral()：c_j^{(α)}（列 α、行 j）= [[1, 2, 3], [1, 2, −3], [1, −1, 0]] ⇒ (3.2.17)
§3.3 第二階段（H305）
  c3v_split_cv() = (P³₁, P³₂)（C_v'' = {1, σ₃}）：P³₁ = (2 − r − r² − σ₁ − σ₂ + 2σ₃)/6、P³₂ = (2 − r − r² + σ₁ + σ₂ − 2σ₃)/6
  c3v_elementary("cv")（參考元素 σ₂：P³₁σ₂P³₂ = (−r + r² − σ₁ + σ₂)/4 = (√3/2)P³₁₂）、c3v_elementary("c3")（C₃ 分裂，參考 σ₃）
  c3v_irreps(kind) → [D^{A1}, D^{A2}, D^E]：kind="cv" 的 D^E(g) = C3V.mats[g]（(3.3.7)）；kind="c3" ＝ T_CIRC†·D^E·T_CIRC（(3.4.24)：
    𝒟^{3°}(r) = diag(ε, ε*)、σ₁ → [[0, ε*], [ε, 0]]、σ₂ → [[0, ε], [ε*, 0]]、σ₃ → [[0, 1], [1, 0]]，ε = EPS3 = e^{2πi/3}）
  (3.3.18)：P^{1₃}ℙ³ = (1 + ε*r + εr²)/3、P^{2₃}ℙ³ = (1 + εr + ε*r²)/3；P^{0₃}ℙ³ = 0
§3.4 正則表示與六自由度模型（H306–H307）
  C3V.regular(g)（6×6，跡 6, 0, 0, 0, 0, 0）；c3v_regular_T()：行 |P̂¹⟩、|P̂²⟩、|P̂³₁₁⟩、|P̂³₂₁⟩、|P̂³₁₂⟩、|P̂³₂₂⟩ ＝ GROUP_VECTORS_BOOK（(3.4.8)）；
    T†𝓡(g)T = D¹ ⊕ D² ⊕ D³ ⊕ D³（(3.4.16)）；r|P̂³₁₁⟩ = −½|P̂³₁₁⟩ + (√3/2)|P̂³₂₁⟩（(3.4.12c)）
  SIX = {k1: 4, k2: 8, m: 1, alpha: 15}（α = 15°：k₂ 沿三角形邊的延長線）；SIX_ANGLES（座標 |g⟩ 的方向，度）、SIX_MASS（在哪個質點）
  six_dof_first_row(k1, k2, alpha, m)（(3.4.21b)）：SIX → 12, 1, 1, 2 − √3 = 0.2679, 2 + √3 = 3.7321, 6
  six_dof_matrix(…)（⟨g|a|h⟩ = ⟨1|a|g⁻¹h⟩）＝ six_dof_geometric(…)（彈簧幾何，浮點）
  six_dof_blocks(…) → (a¹, a², a³ 2×2)（(3.4.22)）＝ six_dof_blocks_book（(3.4.23b)）：SIX → 24、4、[[15, 3], [3, 7]]
  six_dof_omega2(…) → {A1: 24, A2: 4, E: (16, 6)}：ω = 4.8990、2.0000、4.0000、2.4495；k₂ = 0（Fig. 3.4.4）：12、0、(6, 0) ⇒ ω = 3.4641、2.4495、0（三重）
  six_dof_mode(…, which) → (|x⟩, |y⟩, ω)（E 模的型 1、型 2 夥伴，正則座標；which 0 = 高頻）；six_dof_displacements(x6) → 三個質點的 2D 位移
  fig344_vectors("fixed" / "book")（Fig. 3.4.4；書上型 2 兩欄錯）
  six_dof_circular(…, which) → dict(rmin, rmax, sense, omega)：k₂ = 0 時三個質點走半徑 1/√3 的圓、都順時針（|x⟩ + i|y⟩ 屬於 2₃）；
    SIX 時是長短軸比 3 的橢圓
§3.5 特徵標（H308）
  C3V.char_table_class_algebra() = ([1, 1, 2], [[1, 1, 1], [1, 1, −1], [2, −1, 0]]) ＝ char_table_traces(c3v_irreps()) ＝ C3V_CHARS_BOOK
  pendulum_rep(g)（三擺徑向擺動的置換表示；跡 3、0、1）、pendulum_frequencies() = [1, 0, 1]（A1 ⊕ E）、
  pendulum_trace_eigs(a, b) = [b, None, 3a + b]（(3.5.13) 用 °c_g）；PEND3 = {a: 1, b: 1} ⇒ [1, None, 4]；切向擺動時跡 3、0、−1 ⇒ A2 ⊕ E
§3.6 D₄、D₆、D_nh（H309–H311）
  D4.class_alg_rows（類名 D4_CLASS_TEX）：c_R² = 2 + 2c₂、c_Rc₁ = 2c₃、c_Rc₃ = 2c₁、c₁² = 2 + 2c₂、c₁c₃ = 2c_R、c₃² = 2 + 2c₂（(3.6.1b)）
  c_R³ = 4c_R、c₁³ = 4c₁（各只給 3 個冪等元）；book_order(D4, D4_CHARS_BOOK)：A1、B1、A2、B2、E（(3.6.3)）；D4.irreps("R1^2", "i3") 的 E ＝ D4_2D
  ℙ^E 用 C₂ = {1, R₁²}：(1 − R² + R₁² − R₂²)/4、(1 − R² − R₁² + R₂²)/4（(3.6.6b)）；用 C₄：(1 − R² − iR + iR³)/4（R → +i，本系列 1₄）、
    (1 − R² + iR − iR³)/4（R → −i，3₄）；圓形 E：𝒟(R) = diag(i, −i)、𝒟(R₁²) = [[0, 1], [1, 0]]
  correlation：D₄↓C₄：A1, A2 → 0₄；B1, B2 → 2₄；E → 1₄ + 3₄；D₄↓D₂（第 2 章標記）：E → B1 + B2；D₄↓{1, R₁²}：A1, B1 → +；A2, B2 → −；E → + 與 −
  D6、D6_CHARS_BOOK（列 A1 A2 E2 B1 B2 E1）＝ outer_table(C2_CHARS, D3_CHARS)；𝒟^{E₁}(h) = rot2(60)、𝒟^{E₂}(h) = 𝒟^{E₂}(h⁴) = rot2(240)
  D6↓C6、D6↓{1, ρ₃}（(3.6.27)）；D3H、D3D（D3H_CHARS_BOOK、D3H_IRREP_NAMES、D3D_IRREP_NAMES）、D2H（D2_CHARS_BOOK、D2_LABELS）、D4H（10×10）
  轉子 V(φ) = −V cos nφ（谷在 φ = 0；ROTOR = {weak: 2, strong4: 16, strong6: 36, zeeman_b: 1/10}）：
    dn_label_free(n, m, "c"/"s")、dn_free_levels(n, count)、rotor_dn(n, V)、dn_label_vector(n, c, ms)、dn_levels(n, V, count) → [(能量, 標記)]、
    dn_level_paths(n, Vs, count) → (array, labels)、dn_gap_2x2(n, V)、dn_density(n, V, label, φs)、dn_zeeman(n, V, b, count) → [(能量, κ)]
    D4 弱（V = 2）：A1 −0.1242、E 0.8358、B1 2.9697、B2 4.9678、E 9.0975、A2 15.9792（2×2：3、5，能隙 2；數值能隙 1.9980）
    D4 強（V = 16）：A1 −6.0558、E −5.8194、B1 −5.5627 ｜ B2 9.5168、E 11.4095、A2 14.6889
    D6 弱：A1 −0.0555、E1 0.9376、E2 3.9006、B1 7.9863、B2 9.9859、E2 16.0708、E1 25.0311、A2 35.9907（2×2：8、10）
    D6 強（V = 36）：A1 −13.6256、E1 −13.3648、E2 −12.8112、B1 −12.5161 ｜ B2 21.4128、E2 23.3741、E1 28.5724、A2 33.0501
    B1 的 |ψ|² 在谷（φ = 0）最大、在山頂（π/n）為 0；B2 相反。Zeeman（n = 4、V = 2、b = 1/10）：A1 −0.1145、(1₄) 0.6552、(3₄) 1.0340（b = 0：E 0.8358）
§3.7 氨（H312）：NH3 = {m: 1, M: 14, D: 20}、NH3_NU（3337, 908, 3414, 1628 cm⁻¹）、NH3_OMEGA_BOOK、NH3_KJ = 0.3、C_LIGHT、M_H_BOOK
  nh3_geometry(D) → (4×3 位置, 12 個 (原子, 單位向量))、nh3_force(k, j, D)（(3.7.13a)）、nh3_mass、nh3_rep(g)、nh3_projector(α, i, j)
  nh3_cvectors()（(3.7.8)–(3.7.11)，c³ 用圖上的歸一化）、nh3_constraints(v)（(3.7.12)：T = A = 0）
  nh3_reduced(k, j, …) ＝ nh3_reduced_formulas(…)（更正後；book=True 照印）、nh3_secular(k, j, …) → dict(S1, P1, S3, P3, w1p, w1m, w3p, w3m)
  nh3_numeric(…)（12×12 的 ω²：6 個 0 ＋ A1 兩個 ＋ E 兩個二重）
  k/j = 0.3（D = 20°，ω 以 √(j/m) 為單位）：S₁ = 1.9251、P₁ = 0.1278、S₃ = 1.5446、P₃ = 0.2939；ω₁₊ = 1.3624、ω₁₋ = 0.2624、ω₃₊ = 1.1499、ω₃₋ = 0.4714
  nh3_omega_obs()：ω（10¹⁴ rad/s）= 6.2857、1.7104、6.4308、3.0666（書上 6.290、1.712、6.435、3.069 是用 c = 3×10¹⁰）
  nh3_fit()：R = 0.2972、D = 16.22°（書上 17°）、j = 847.5 nt/m（書 850）、k = −53.4（更正；書上的式子給 −62.1，書 −62）
  nh3_ratio_fit()：ω₃₋/ω₁₋ = 1.7930 ⇒ k/j = 0.2956；以兩個低頻定尺度：ω₁₊ 預測 8.884（+41.3%）、ω₃₊ 7.513（+16.8%）
  nh3_crossing() = 0.0348（ω₁₊ = ω₃₊ 的偶然簡併；書 0.035）；nh3_scan(k/j)；nh3_modes(k/j) → [(名, ω, 12 維位移)]（"1+"、"1-"、"3+"、"3-"）
習題：problem_311（陪集）、problem_312、problem_322、problem_331（D₄ 類代數與冪等元、D₂／C₄ 分裂）、problem_332（四元數群；Q ≇ D₄）、
  d4_square_matrix(e, s, d, m)、problem_341(e, s, d, m)（P341：(g) e = 0、s = 5、d = 3：ω² = 0, 0, 0, 6, 10, 10, 10, 16；
  (h) e = 1：ω = 0, 0.6892（二重）, 2.6458, 3.1623, 3.2442（二重）, 4.1231）、problem_351、problem_352、problem_355、problem_356
check()：每個上卡片／講義的數字用兩種以上算法核對（約 3 秒）；summary() → {名稱: 字串}。
直接執行：python utils/sym_nonabel.py（在 codes/ 裡）印一句成功訊息；--summary 印出主要數字。
"""
import math
from functools import cached_property

import numpy as np
import sympy as sp
from scipy.optimize import brentq

try:
    from . import sym_abel as SA
    from . import sym_mat as SM          # noqa: F401（給場景用：SN.SM）
except ImportError:                      # 直接執行 python utils/sym_nonabel.py
    import sym_abel as SA
    import sym_mat as SM                 # noqa: F401

R_ = sp.Rational
I = sp.I
SQ3 = sp.sqrt(3)
EPS3 = R_(-1, 2) + SQ3 / 2 * I            # ε = e^{2πi/3}（本系列的約定，見 docstring）


def _close(a, b, tol: float = 1e-9) -> bool:
    return np.allclose(np.asarray(a, dtype=complex), np.asarray(b, dtype=complex), atol=tol, rtol=0)


def _f4(x) -> str:
    return f"{float(x):.4f}"


def _zero(x) -> bool:
    return sp.expand(x) == 0


# ======================================================================== 幾何：轉動、反射（精確）
def cosd(deg):
    return sp.cos(sp.pi * sp.nsimplify(deg) / 180)


def sind(deg):
    return sp.sin(sp.pi * sp.nsimplify(deg) / 180)


def rot2(deg) -> sp.Matrix:
    """xy 平面逆時針轉 deg 度：[[cos, −sin], [sin, cos]]（(3.1.4) 的 ⟨R[2φ]⟩、(3.6.19)）。"""
    c, s = cosd(deg), sind(deg)
    return sp.Matrix([[c, -s], [s, c]])


def refl2(deg) -> sp.Matrix:
    """對「與 x 軸夾 deg 度的直線」（3D：含 z 軸的 φz 面）反射：[[cos 2φ, sin 2φ], [sin 2φ, −cos 2φ]]（Fig. 3.1.4(b)）。"""
    c, s = cosd(2 * sp.nsimplify(deg)), sind(2 * sp.nsimplify(deg))
    return sp.Matrix([[c, s], [s, -c]])


def rot3(axis, deg) -> sp.Matrix:
    """繞 axis（不必歸一）轉 deg 度的 3×3 精確矩陣（Rodrigues）。"""
    n = sp.Matrix(axis)
    n = n / sp.sqrt(n.dot(n))
    c, s = cosd(deg), sind(deg)
    K = sp.Matrix([[0, -n[2], n[1]], [n[2], 0, -n[0]], [-n[1], n[0], 0]])
    return (c * sp.eye(3) + s * K + (1 - c) * n * n.T).applyfunc(lambda x: sp.nsimplify(sp.expand(x)))


def e2(deg) -> np.ndarray:
    """平面單位向量 (cos, sin)（浮點）。"""
    t = math.radians(float(deg))
    return np.array([math.cos(t), math.sin(t)])


def block3(M2) -> sp.Matrix:
    """2×2 → diag(M2, 1)（鉛直軸不動）。"""
    M = sp.eye(3)
    M[:2, :2] = sp.Matrix(M2)
    return M


# ======================================================================== 有限群
def _key(M) -> tuple:
    a = np.asarray(M, dtype=complex).ravel()
    return tuple(np.round(a.real, 7) + 0.0) + tuple(np.round(a.imag, 7) + 0.0)


class FiniteGroup:
    """由矩陣（精確 sympy）定義的有限群。乘積 ab = 先 b 後 a（書上的約定，p. 152）。

    names 的順序就是乘法表、正則表示、類（依第一次出現）的順序。群代數的元素用 dict {名稱: 係數}。"""

    def __init__(self, mats: dict, tex: dict | None = None, label: str = ""):
        self.label = label
        self.names = list(mats)
        self.mats = {g: sp.Matrix(mats[g]) for g in self.names}
        self.tex = {g: g for g in self.names}
        self.tex.update(tex or {})
        self.order = len(self.names)
        self.idx = {g: i for i, g in enumerate(self.names)}
        self.num = {g: np.array(self.mats[g].evalf(), dtype=complex) for g in self.names}
        look = {_key(M): g for g, M in self.num.items()}
        if len(look) != self.order:
            raise ValueError(f"{label}：有重複的元素")
        self.table = {}
        for a in self.names:
            for b in self.names:
                k = _key(self.num[a] @ self.num[b])
                if k not in look:
                    raise ValueError(f"{label} 不封閉：{a}·{b}")
                self.table[a, b] = look[k]
        self.e = look[_key(np.eye(self.num[self.names[0]].shape[0]))]
        self._inv = {a: next(b for b in self.names if self.table[a, b] == self.e) for a in self.names}

    # ---------------------------------------------------------- 群
    @classmethod
    def from_generators(cls, gens: dict, tex: dict | None = None, label: str = "", identity: str = "1"):
        """由生成元（精確矩陣）做封閉：元素名是最短的字（"r s" 寫成 "rs"，連續相同字母寫成冪 "r^2"；ab = 先 b）。"""
        gens = {g: sp.Matrix(M) for g, M in gens.items()}
        n = next(iter(gens.values())).shape[0]
        found = {_key(np.eye(n)): (identity, sp.eye(n), [])}
        frontier = [(sp.eye(n), [])]
        while frontier:
            new = []
            for M, w in frontier:
                for g, G in gens.items():
                    P = (G * M).applyfunc(sp.expand)
                    k = _key(np.array(P.evalf(), dtype=complex))
                    if k not in found:
                        found[k] = (None, P, [g] + w)
                        new.append((P, [g] + w))
            frontier = new
        mats = {}
        for nm, P, w in found.values():
            mats[nm or _word_name(w)] = P
        return cls(mats, tex=tex, label=label)

    def mul(self, *gs) -> str:
        """g1·g2·…·gk（最右邊的最先作用）。"""
        out = self.e
        for g in gs:
            out = self.table[out, g]
        return out

    def inv(self, g: str) -> str:
        return self._inv[g]

    def conj(self, h: str, g: str) -> str:
        """h g h⁻¹（(3.2.4)、Definition p. 159）。"""
        return self.mul(h, g, self.inv(h))

    def order_of(self, g: str) -> int:
        k, x = 1, g
        while x != self.e:
            x, k = self.table[x, g], k + 1
        return k

    def table_rows(self, order=None) -> list:
        """乘法表（列 a、行 b、格子 ab）。"""
        order = order or self.names
        return [[self.table[a, b] for b in order] for a in order]

    def is_abelian(self) -> bool:
        return all(self.table[a, b] == self.table[b, a] for a in self.names for b in self.names)

    def is_subgroup(self, H) -> bool:
        return self.e in H and all(self.table[a, b] in H for a in H for b in H)

    def subgroup_of(self, g: str) -> list:
        """g 生成的循環子群 [1, g, g², …]（依冪次）。"""
        out, x = [self.e], g
        while x != self.e:
            out.append(x)
            x = self.table[x, g]
        return out

    def cosets(self, H, side: str = "left") -> list:
        """左陪集 gH（side="left"）或右陪集 Hg（"right"），依 names 順序找代表元。"""
        seen, out = set(), []
        for g in self.names:
            if g in seen:
                continue
            C = [self.table[g, h] for h in H] if side == "left" else [self.table[h, g] for h in H]
            C = sorted(set(C), key=self.idx.get)
            out.append(C)
            seen |= set(C)
        return out

    def double_cosets(self, H, K) -> list:
        """雙陪集 HgK（習題 3.1.1(e)(f)）。"""
        seen, out = set(), []
        for g in self.names:
            if g in seen:
                continue
            C = sorted({self.mul(h, g, k) for h in H for k in K}, key=self.idx.get)
            out.append(C)
            seen |= set(C)
        return out

    def is_normal(self, H) -> bool:
        return all(self.conj(g, h) in H for g in self.names for h in H)

    def commutator_subgroup(self) -> list:
        """所有 a⁻¹b⁻¹ab 生成的子群（習題 3.5.4）。"""
        S = {self.mul(self.inv(a), self.inv(b), a, b) for a in self.names for b in self.names}
        while True:
            T = S | {self.table[a, b] for a in S for b in S}
            if T == S:
                return sorted(S, key=self.idx.get)
            S = T

    # ---------------------------------------------------------- 類
    @cached_property
    def classes(self) -> list:
        """共軛類（依 names 中第一次出現的順序；類內也依 names 順序）。"""
        seen, out = set(), []
        for g in self.names:
            if g in seen:
                continue
            C = sorted({self.conj(h, g) for h in self.names}, key=self.idx.get)
            out.append(C)
            seen |= set(C)
        return out

    def class_of(self, g: str) -> int:
        return next(k for k, C in enumerate(self.classes) if g in C)

    def normalizer(self, g: str) -> list:
        """N_g = 與 g 對易的元素（p. 160）；Lagrange：°c_g = °G / °N_g（(3.2.7)）。"""
        return [h for h in self.names if self.table[h, g] == self.table[g, h]]

    def class_sizes(self) -> list:
        return [len(C) for C in self.classes]

    @cached_property
    def struct(self) -> dict:
        """結構常數 (3.2.11)：struct[i, j][k] = γ^k_{ij}，c_i c_j = Σ_k γ^k_{ij} c_k。"""
        out = {}
        reps = [C[0] for C in self.classes]
        for i, Ci in enumerate(self.classes):
            for j, Cj in enumerate(self.classes):
                cnt = {}
                for a in Ci:
                    for b in Cj:
                        p = self.table[a, b]
                        cnt[p] = cnt.get(p, 0) + 1
                out[i, j] = tuple(cnt.get(r, 0) for r in reps)
        return out

    def class_mul(self, x, y) -> sp.Matrix:
        """類代數裡兩個元素（係數向量）相乘。"""
        n = len(self.classes)
        out = [0] * n
        for i in range(n):
            if x[i] == 0:
                continue
            for j in range(n):
                if y[j] == 0:
                    continue
                for k, g in enumerate(self.struct[i, j]):
                    if g:
                        out[k] += x[i] * y[j] * g
        return sp.Matrix([sp.expand(v) for v in out])

    def class_matrix(self, i: int) -> sp.Matrix:
        """類代數的正則表示（習題 3.5.1）：(C_i)_{kj} = γ^k_{ij}。"""
        n = len(self.classes)
        return sp.Matrix(n, n, lambda k, j: self.struct[i, j][k])

    def class_alg_rows(self, symbols=None) -> list:
        """類代數表 (3.2.1b)、(3.6.1b)：rows[i][j] = c_i c_j（sympy 式，用 symbols 或 c_1..c_n）。"""
        n = len(self.classes)
        cs = symbols or sp.symbols(f"c_1:{n + 1}")
        return [[sum(g * cs[k] for k, g in enumerate(self.struct[i, j])) for j in range(n)] for i in range(n)]

    def min_poly_class(self, i: int, x=sp.Symbol("c")):
        """類算子 c_i 的最小方程 (3.2.12)–(3.2.13)（首一 Poly）：由 c_i 的冪在類代數中第一次線性相依求得。"""
        n = len(self.classes)
        one = sp.Matrix([1] + [0] * (n - 1))
        ci = sp.Matrix([1 if k == i else 0 for k in range(n)])
        pw = [one]
        while True:
            nxt = self.class_mul(pw[-1], ci)
            A = sp.Matrix.hstack(*pw)
            sol = None
            try:
                sol, params = A.gauss_jordan_solve(nxt)
                if params.shape[0]:
                    sol = sol.subs({p: 0 for p in params})
            except ValueError:
                sol = None
            if sol is not None:
                return sp.Poly(x ** len(pw) - sum(sol[k] * x ** k for k in range(len(pw))), x)
            pw.append(nxt)

    def class_idempotents_of(self, i: int) -> dict:
        """c_i 的最小方程的每個根 λ → 冪等元 Π(c_i − λ')/(λ − λ')（(3.2.14)，類係數向量）。"""
        n = len(self.classes)
        P = self.min_poly_class(i)
        roots = list(sp.roots(P).keys())
        ci = sp.Matrix([1 if k == i else 0 for k in range(n)])
        one = sp.Matrix([1] + [0] * (n - 1))
        out = {}
        for lam in roots:
            v = one
            for mu in roots:
                if mu != lam:
                    v = self.class_mul(v, (ci - mu * one) / (lam - mu))
            out[lam] = v.applyfunc(sp.nsimplify)
        return out

    @cached_property
    def class_idempotents(self) -> list:
        """全對易冪等元 ℙ^α（類係數向量）：逐一用各類算子的最小方程分裂、兩兩相乘（§3.2B；D₄ 要兩個方程，p. 191）。
        排序：先依維度 ℓ，再依特徵標（大的在前）；具體的群用 book_order 換成書上的順序。"""
        n = len(self.classes)
        ids = [sp.Matrix([1] + [0] * (n - 1))]
        for i in range(1, n):
            if len(ids) == n:
                break
            new = []
            for P in ids:
                for Q in self.class_idempotents_of(i).values():
                    R = self.class_mul(P, Q)
                    if any(x != 0 for x in R):
                        new.append(R)
            ids = new
        if len(ids) != n:
            raise ValueError("類代數沒有分裂完全")
        chars = [self._chars_of_idempotent(P) for P in ids]
        order = sorted(range(n), key=lambda a: (chars[a][0], [-float(sp.re(x)) for x in chars[a][1]],
                                                [-float(sp.im(x)) for x in chars[a][1]]))
        return [ids[a] for a in order]

    def _chars_of_idempotent(self, P) -> tuple:
        """(3.5.5)–(3.5.7)：ℓ = (°G p_1)^{1/2}、χ_k* = p_k °G/ℓ。"""
        ell = sp.sqrt(self.order * P[0])
        ell = sp.nsimplify(ell)
        return ell, [sp.nsimplify(sp.conjugate(p * self.order / ell)) for p in P]

    def char_table_class_algebra(self) -> tuple:
        """(維度 ℓ 的列表, 特徵標表 Matrix[α, 類])：只用類代數（§3.5A）。"""
        rows = [self._chars_of_idempotent(P) for P in self.class_idempotents]
        return [r[0] for r in rows], sp.Matrix([r[1] for r in rows])

    # ---------------------------------------------------------- 群代數
    def ga(self, d: dict) -> dict:
        return {g: sp.expand(c) for g, c in d.items() if not _zero(c)}

    def ga_mul(self, x: dict, y: dict) -> dict:
        out = {}
        for a, ca in x.items():
            for b, cb in y.items():
                p = self.table[a, b]
                out[p] = out.get(p, 0) + ca * cb
        return self.ga(out)

    def ga_add(self, *xs, coeffs=None) -> dict:
        out = {}
        for k, x in enumerate(xs):
            c = 1 if coeffs is None else coeffs[k]
            for g, v in x.items():
                out[g] = out.get(g, 0) + c * v
        return self.ga(out)

    def ga_scale(self, x: dict, c) -> dict:
        return self.ga({g: c * v for g, v in x.items()})

    def ga_dag(self, x: dict) -> dict:
        """x† = Σ x_g* g⁻¹（么正表示）。"""
        return self.ga({self.inv(g): sp.conjugate(v) for g, v in x.items()})

    def ga_eq(self, x: dict, y: dict) -> bool:
        return not self.ga_add(x, y, coeffs=[1, -1])

    def el(self, g: str) -> dict:
        return {g: sp.Integer(1)}

    def class_op(self, i: int) -> dict:
        """類算子 c_i = 類中元素之和（(3.2.1a)）。"""
        return {g: sp.Integer(1) for g in self.classes[i]}

    def from_class_vec(self, v) -> dict:
        return self.ga({g: v[i] for i, C in enumerate(self.classes) for g in C})

    @cached_property
    def central_idempotents(self) -> list:
        """ℙ^α 寫成群代數元素（與 class_idempotents 同順序）。"""
        return [self.from_class_vec(P) for P in self.class_idempotents]

    def regular(self, g: str) -> sp.Matrix:
        """正則表示 (3.4.3)：𝓡_{hf}(g) = δ_{g, h f⁻¹}（列 ⟨h|、行 |f⟩，依 names 順序）；g|f⟩ = |gf⟩。"""
        n = self.order
        return sp.Matrix(n, n, lambda i, j: 1 if self.table[g, self.names[j]] == self.names[i] else 0)

    def regular_of(self, x: dict) -> sp.Matrix:
        return sum((c * self.regular(g) for g, c in x.items()), sp.zeros(self.order))

    def vec_of(self, x: dict) -> sp.Matrix:
        """|x⟩ = x|1⟩ 的分量（(3.4.6a)：⟨g|x⟩ = x 裡 g 的係數）。"""
        return sp.Matrix([x.get(g, 0) for g in self.names])

    # ---------------------------------------------------------- 第二階段：子群分裂、基本算子、不可約表示
    def cyclic_projectors(self, gen: str) -> list:
        """循環子群 {1, g, g², …} 的投影算子 P^{k_n} = (1/n) Σ_p ε_n^{−kp} g^p（g P^{k} = ε_n^k P^{k}，ε_n = e^{2πi/n}）。"""
        H = self.subgroup_of(gen)
        n = len(H)
        return [self.ga({H[p]: sp.conjugate(SA.eps(n, k * p)) / n for p in range(n)}) for k in range(n)]

    def split(self, P: dict, projs: list) -> list:
        """ℙ^α 乘上子群的投影算子，留下不為零的（(3.3.2)、(3.3.18)、(3.6.4)、(3.6.6)）。"""
        out = []
        for p in projs:
            Q = self.ga_mul(P, p)
            if Q:
                out.append(Q)
        return out

    def elementary(self, Ps: list, ref: str | None = None) -> tuple:
        """由互相正交的不可約冪等元 [P_1, …, P_ℓ] 做基本算子 P_{ij}（(3.3.11)）：P_{1j} = P_1 g P_j / D_{1j}(g)，
        g = ref（沒給就取第一個使 P_1 g P_2 ≠ 0 的元素），D_{1j}(ref) = √λ > 0（P_1 g P_j (…)† = λ P_1）；
        P_{j1} = P_{1j}†、P_{ij} = P_{i1} P_{1j}。回傳 ({(i, j): 群代數元素}（0 起算）, ref)。"""
        ell = len(Ps)
        if ell == 1:
            return {(0, 0): Ps[0]}, ref
        if ref is None:
            ref = next(g for g in self.names if self.ga_mul(self.ga_mul(Ps[0], self.el(g)), Ps[1]))
        P1j = {0: Ps[0]}
        for j in range(1, ell):
            N = self.ga_mul(self.ga_mul(Ps[0], self.el(ref)), Ps[j])
            NN = self.ga_mul(N, self.ga_dag(N))
            g0 = next(iter(Ps[0]))
            lam = sp.nsimplify(sp.simplify(NN[g0] / Ps[0][g0]))
            P1j[j] = self.ga_scale(N, 1 / sp.sqrt(lam))
        Pi1 = {j: self.ga_dag(P1j[j]) for j in range(ell)}
        Pi1[0] = Ps[0]
        out = {}
        for i in range(ell):
            for j in range(ell):
                out[i, j] = Ps[0] if i == j == 0 else self.ga_mul(Pi1[i], P1j[j])
        return out, ref

    def irrep_from_elementary(self, Pij: dict) -> dict:
        """P_i g P_j = D_{ij}(g) P_{ij}（(3.3.11)）⇒ {g: D(g)}。"""
        ell = int(round(math.sqrt(len(Pij))))
        out = {}
        for g in self.names:
            D = sp.zeros(ell)
            for i in range(ell):
                for j in range(ell):
                    X = self.ga_mul(self.ga_mul(Pij[i, i], self.el(g)), Pij[j, j])
                    if X:
                        h = next(iter(Pij[i, j]))
                        D[i, j] = sp.nsimplify(sp.simplify(X.get(h, 0) / Pij[i, j][h]))
            out[g] = D
        return out

    def irreps(self, sub_gen: str | None = None, ref: str | None = None) -> list:
        """全部不可約表示（依 class_idempotents 的順序）：一維的由 gℙ = D(g)ℙ，多維的用 sub_gen 生成的循環子群分裂。"""
        out = []
        for P in self.central_idempotents:
            g0 = next(iter(P))
            ell = int(sp.sqrt(self.order * P.get(self.e, 0)))
            if ell == 1:
                out.append({g: sp.Matrix([[sp.nsimplify(sp.simplify(self.ga_mul(self.el(g), P)[g0] / P[g0]))]])
                            for g in self.names})
                continue
            Ps = self.split(P, self.cyclic_projectors(sub_gen))
            if len(Ps) != ell:
                raise ValueError(f"用 {sub_gen} 分不開 ℓ = {ell} 的冪等元")
            Pij, _ = self.elementary(Ps, ref)
            out.append(self.irrep_from_elementary(Pij))
        return out

    def is_rep(self, D: dict) -> bool:
        """D(g)D(h) = D(gh) 對所有 g、h（(3.3.13)）。"""
        return all(sp.simplify(D[a] * D[b] - D[self.table[a, b]]) == sp.zeros(*D[a].shape)
                   for a in self.names for b in self.names)

    def char_table_traces(self, irreps: list) -> sp.Matrix:
        """χ^α_k = Tr D^α(類 k 的代表元)（(3.5.1)）。"""
        return sp.Matrix([[sp.nsimplify(sp.simplify(D[C[0]].trace())) for C in self.classes] for D in irreps])

    def frequencies(self, chars, table: sp.Matrix | None = None) -> list:
        """(3.5.11)：f^α = (1/°G) Σ_k °c_k χ^α_k* χ_k（chars：表示在每個類的跡）。"""
        T = table if table is not None else self.char_table_class_algebra()[1]
        sz = self.class_sizes()
        return [sp.nsimplify(sp.simplify(sum(sz[k] * sp.conjugate(T[a, k]) * chars[k] for k in range(len(sz))) / self.order))
                for a in range(T.rows)]

    def correlation(self, sub: list, sub_chars: sp.Matrix, table: sp.Matrix | None = None) -> sp.Matrix:
        """相關（subduction）表 (3.6.15)、(3.6.27)：n_{αβ} = (1/°H) Σ_{h∈H} χ^α(h) χ^β(h)*（sub_chars 的行依 sub 的順序）。"""
        T = table if table is not None else self.char_table_class_algebra()[1]
        return sp.Matrix(T.rows, sub_chars.rows, lambda a, b: sp.nsimplify(sp.simplify(
            sum(T[a, self.class_of(h)] * sp.conjugate(sub_chars[b, p]) for p, h in enumerate(sub)) / len(sub))))

    def iso_to(self, other: "FiniteGroup", mapping: dict) -> bool:
        """mapping（名稱 → 名稱）是同構：保持乘法表。"""
        return all(mapping[self.table[a, b]] == other.table[mapping[a], mapping[b]] for a in self.names for b in self.names)


def ga_tex(G: FiniteGroup, x: dict, factor=None, named: dict | None = None) -> str:
    """群代數元素 → TeX，提出公因數：{1: 1/3, r: −1/6, r²: −1/6} → \\tfrac{1}{6}(2\\mathbf 1-r-r^2)（依 names 順序）。
    factor 不給時取第一個係數的絕對值再調成讓其餘係數都是（Gauss）整數。named：{係數: TeX}，例如
    ga_tex(SN.C3V, P, R_(1, 3), {EPS3: r"\varepsilon", sp.conjugate(EPS3): r"\varepsilon^*"}) → \tfrac{1}{3}(\mathbf 1+\varepsilon^* r+\varepsilon r^2)。"""
    items = [(g, x[g]) for g in G.names if g in x and not _zero(x[g])]
    if not items:
        return "0"
    if factor is None:
        c0 = sp.Abs(items[0][1])
        den = 1
        for _, v in items:
            for part in (sp.re(v / c0), sp.im(v / c0)):
                den = sp.ilcm(den, sp.fraction(sp.nsimplify(part))[1]) if sp.nsimplify(part).is_rational else den
        factor = sp.nsimplify(c0 / den)
    out = ""
    for g, v in items:
        c = sp.nsimplify(sp.expand(v / factor))
        name = G.tex[g]
        hit = next((s for k, s in (named or {}).items() if sp.expand(sp.nsimplify(k) - c) == 0), None)
        if hit is not None:
            out += "+" + hit + " " + name
            continue
        if c == 1:
            term = "+" + name
        elif c == -1:
            term = "-" + name
        elif c.is_real:
            term = ("+" if c > 0 else "-") + sp.latex(abs(c)) + name
        else:
            s = sp.latex(c).replace(" ", "")
            term = ("+(" + s + ")" if sp.re(c) != 0 else ("" if s.startswith("-") else "+") + s) + name
        out += term
    out = out.lstrip("+")
    if factor == 1:
        return out
    return (r"\tfrac{%s}{%s}" % sp.fraction(factor) if factor.is_rational else sp.latex(factor)) + "(" + out + ")"


def _word_name(w: list) -> str:
    out, k = [], 0
    while k < len(w):
        j = k
        while j < len(w) and w[j] == w[k]:
            j += 1
        out.append(w[k] if j - k == 1 else f"{w[k]}^{j - k}")
        k = j
    return "".join(out)


def reorder(T: sp.Matrix, ref_rows) -> list:
    """把特徵標表 T 的列換成 ref_rows 的順序：回傳 perm，T[perm[a], :] = ref_rows[a]。"""
    perm = []
    for row in ref_rows:
        a = next(a for a in range(T.rows) if all(sp.simplify(T[a, k] - row[k]) == 0 for k in range(T.cols)))
        perm.append(a)
    return perm


def outer_table(T1: sp.Matrix, T2: sp.Matrix) -> sp.Matrix:
    """外積的特徵標表 (3.6.17)：kron(T1, T2)（列 (α, β)、行 (類 1, 類 2)，T1 的指標在外）。"""
    return SA.outer_char_table(T1, T2)


def cn_chars(n: int) -> sp.Matrix:
    """C_n：χ^{k_n}(g^p) = ε_n^{kp}（列 k、行 p；本系列的約定 𝒟^{k_n}(r) = e^{2πik/n}）。"""
    return SA.cn_char_table(n)


# ======================================================================== §3.1 具體的群
C3V_NAMES = ["1", "r", "r^2", "s1", "s2", "s3"]
C3V_TEX = {"1": r"\mathbf 1", "r": "r", "r^2": "r^2", "s1": r"\sigma_1", "s2": r"\sigma_2", "s3": r"\sigma_3"}
MIRROR_DEG = {"s1": 120, "s2": 60, "s3": 0}          # σ_j 的鏡線與 x 軸的夾角（Fig. 3.1.2：σ₃ 水平、σ₂ 右上、σ₁ 右下）
C3V = FiniteGroup({"1": sp.eye(2), "r": rot2(120), "r^2": rot2(240),
                   "s1": refl2(120), "s2": refl2(60), "s3": refl2(0)}, C3V_TEX, "C3v")
C3V_TABLE_BOOK = [["1", "r", "r^2", "s1", "s2", "s3"],          # (3.1.1) 照抄（列 × 行）
                  ["r", "r^2", "1", "s3", "s1", "s2"],
                  ["r^2", "1", "r", "s2", "s3", "s1"],
                  ["s1", "s2", "s3", "1", "r", "r^2"],
                  ["s2", "s3", "s1", "r^2", "1", "r"],
                  ["s3", "s1", "s2", "r", "r^2", "1"]]
C3V_CLASS_ALG_BOOK = [["c1", "c2", "c3"], ["c2", "2*c1+c2", "2*c3"], ["c3", "2*c3", "3*c1+3*c2"]]   # (3.2.1b)
C3V_CHARS_BOOK = [[1, 1, 1], [1, 1, -1], [2, -1, 0]]           # (3.5.8)：A1、A2、E 對 1、(r, r²)、(σ₁σ₂σ₃)
C3V_IRREP_NAMES = ["A_1", "A_2", "E"]

D3_AXIS_DEG = {"p1": 210, "p2": 150, "p3": 90}        # ρ_j 的軸 ⊥ σ_j 面：φ_j + 90°（ρ₃ 繞 y 軸，Fig. 3.1.3）
D3 = FiniteGroup({"1": sp.eye(3), "r": rot3((0, 0, 1), 120), "r^2": rot3((0, 0, 1), 240),
                  **{p: rot3((cosd(a), sind(a), 0), 180) for p, a in D3_AXIS_DEG.items()}},
                 {"1": r"\mathbf 1", "r": "r", "r^2": "r^2", "p1": r"\rho_1", "p2": r"\rho_2", "p3": r"\rho_3"}, "D3")
D3_TO_C3V = {"1": "1", "r": "r", "r^2": "r^2", "p1": "s1", "p2": "s2", "p3": "s3"}
INV3 = -sp.eye(3)

D4_NAMES = ["1", "R^2", "R", "R^3", "R1^2", "R2^2", "i3", "i4"]
D4_TEX = {"1": r"\mathbf 1", "R^2": "R^2", "R": "R", "R^3": "R^3", "R1^2": "R_1^2", "R2^2": "R_2^2", "i3": "i_3", "i4": "i_4"}
D4 = FiniteGroup({"1": sp.eye(3), "R^2": rot3((0, 0, 1), 180), "R": rot3((0, 0, 1), 90), "R^3": rot3((0, 0, 1), 270),
                  "R1^2": rot3((1, 0, 0), 180), "R2^2": rot3((0, 1, 0), 180),
                  "i3": rot3((1, 1, 0), 180), "i4": rot3((1, -1, 0), 180)}, D4_TEX, "D4")
D4_TABLE_BOOK = [["1", "R^2", "R", "R^3", "R1^2", "R2^2", "i3", "i4"],     # (3.6.1a) 照抄
                 ["R^2", "1", "R^3", "R", "R2^2", "R1^2", "i4", "i3"],
                 ["R", "R^3", "R^2", "1", "i3", "i4", "R2^2", "R1^2"],
                 ["R^3", "R", "1", "R^2", "i4", "i3", "R1^2", "R2^2"],
                 ["R1^2", "R2^2", "i4", "i3", "1", "R^2", "R^3", "R"],
                 ["R2^2", "R1^2", "i3", "i4", "R^2", "1", "R", "R^3"],
                 ["i3", "i4", "R1^2", "R2^2", "R", "R^3", "1", "R^2"],
                 ["i4", "i3", "R2^2", "R1^2", "R^3", "R", "R^2", "1"]]
D4_CLASS_TEX = [r"\mathbf 1", "c_2", "c_R", "c_1", "c_3"]       # (3.6.1b) 的類名：{1}、{R²}、{R, R³}、{R₁², R₂²}、{i₃, i₄}
D4_CHARS_BOOK = [[1, 1, 1, 1, 1], [1, 1, -1, 1, -1], [1, 1, 1, -1, -1], [1, 1, -1, -1, 1], [2, -2, 0, 0, 0]]   # (3.6.3)
D4_IRREP_NAMES = ["A_1", "B_1", "A_2", "B_2", "E"]
D4_2D = {"1": sp.eye(2), "R^2": rot2(180), "R": rot2(90), "R^3": rot2(270), "R1^2": refl2(0), "R2^2": refl2(90),
         "i3": refl2(45), "i4": refl2(-45)}                    # (3.6.7) 的駐波 E ＝ xy 平面的作用
C4V = FiniteGroup({**{g: D4.mats[g] for g in ["1", "R^2", "R", "R^3"]},
                   "s1": INV3 * D4.mats["R1^2"], "s2": INV3 * D4.mats["R2^2"], "s3": INV3 * D4.mats["i3"],
                   "s4": INV3 * D4.mats["i4"]}, {"1": r"\mathbf 1", "s1": r"\sigma_1", "s2": r"\sigma_2",
                                                 "s3": r"\sigma_3", "s4": r"\sigma_4"}, "C4v")
D4_TO_C4V = {"1": "1", "R^2": "R^2", "R": "R", "R^3": "R^3", "R1^2": "s1", "R2^2": "s2", "i3": "s3", "i4": "s4"}

# D₆ = D₃ × C₂（(3.6.16)，書上的元素順序）
_h = lambda k: rot3((0, 0, 1), 60 * k)                         # noqa: E731
D6 = FiniteGroup({"1": sp.eye(3), "h^2": _h(2), "h^4": _h(4), "p1": D3.mats["p1"], "p2": D3.mats["p2"], "p3": D3.mats["p3"],
                  "h^3": _h(3), "h^5": _h(5), "h": _h(1), "p1'": D3.mats["p1"] * _h(3), "p2'": D3.mats["p2"] * _h(3),
                  "p3'": D3.mats["p3"] * _h(3)},
                 {"1": r"\mathbf 1", "p1": r"\rho_1", "p2": r"\rho_2", "p3": r"\rho_3", "p1'": r"\rho_1'",
                  "p2'": r"\rho_2'", "p3'": r"\rho_3'"}, "D6")
D3_CHARS = sp.Matrix(C3V_CHARS_BOOK)
C2_CHARS = sp.Matrix([[1, 1], [1, -1]])
D6_CHARS_BOOK = [[1, 1, 1, 1, 1, 1], [1, 1, -1, 1, 1, -1], [2, -1, 0, 2, -1, 0],     # (3.6.17)：A1 A2 E2 B1 B2 E1
                 [1, 1, 1, -1, -1, -1], [1, 1, -1, -1, -1, 1], [2, -1, 0, -2, 1, 0]]
D6_IRREP_NAMES = ["A_1", "A_2", "E_2", "B_1", "B_2", "E_1"]


def direct_product(G: FiniteGroup, extra: str, M, tex_extra: str, label: str) -> FiniteGroup:
    """G × {1, X}（X 與 G 全部對易，例如 σ_h、I）：名稱 g、X.g（(3.6.28)–(3.6.30)）。"""
    mats = dict(G.mats)
    tex = dict(G.tex)
    for g in G.names:
        nm = extra if g == G.e else f"{extra}.{g}"
        mats[nm] = sp.Matrix(M) * G.mats[g]
        tex[nm] = tex_extra if g == G.e else f"{tex_extra} {G.tex[g]}"
    return FiniteGroup(mats, tex, label)


SIGMA_H = sp.diag(1, 1, -1)
D3H = direct_product(D3, "sh", SIGMA_H, r"\sigma_h", "D3h")   # (3.6.36)
D3D = direct_product(D3, "I", INV3, "I", "D3d")                # (3.6.34)
D4H = direct_product(D4, "I", INV3, "I", "D4h")                # (3.6.29)
D2 = FiniteGroup({"1": sp.eye(3), "Rz": rot3((0, 0, 1), 180), "Ry": rot3((0, 1, 0), 180), "Rx": rot3((1, 0, 0), 180)},
                 {"1": r"\mathbf 1", "Rz": "R_z^2", "Ry": "R_y^2", "Rx": "R_x^2"}, "D2")
D2H = direct_product(D2, "I", INV3, "I", "D2h")                # (3.6.33)
D3H_CHARS_BOOK = [[1, 1, 1, 1, 1, 1], [1, 1, -1, 1, 1, -1], [2, -1, 0, 2, -1, 0],
                  [1, 1, 1, -1, -1, -1], [1, 1, -1, -1, -1, 1], [2, -1, 0, -2, 1, 0]]
D3H_IRREP_NAMES = ["A_1'", "A_2'", "E'", "A_1''", "A_2''", "E''"]
D3D_IRREP_NAMES = ["A_{1g}", "A_{2g}", "E_g", "A_{1u}", "A_{2u}", "E_u"]
D2_CHARS_BOOK = [[1, 1, 1, 1], [1, -1, 1, -1], [1, 1, -1, -1], [1, -1, -1, 1]]   # p. 204（對 1、R_z²、R_y²、R_x²）
D2_LABELS_BOOK = ["A_1", "B_2", "A_2", "B_2"]                  # 書上第 2 列也印 B₂（筆誤）
D2_LABELS = ["A_1", "B_1", "A_2", "B_2"]
D2H_IRREP_NAMES = ["A_{1g}", "B_{1g}", "A_{2g}", "B_{2g}", "A_{1u}", "B_{1u}", "A_{2u}", "B_{2u}"]
D4H_IRREP_NAMES = [n + s for s in ("_g", "_u") for n in ("A_1", "B_1", "A_2", "B_2", "E")]

# 四元數群（習題 3.3.2、3.4.2）：SU(2) 的 ±1、±iσ_x…；i = R_x(180°) 的自旋 ½ 版本
_sx, _sy, _sz = sp.Matrix([[0, 1], [1, 0]]), sp.Matrix([[0, -I], [I, 0]]), sp.Matrix([[1, 0], [0, -1]])
QUAT = FiniteGroup({"1": sp.eye(2), "-1": -sp.eye(2), "i": -I * _sx, "-i": I * _sx, "j": -I * _sy, "-j": I * _sy,
                    "k": -I * _sz, "-k": I * _sz},
                   {"1": r"\mathbf 1", "-1": r"\underline{\mathbf 1}", "i": r"\mathbf i", "-i": r"\underline{\mathbf i}",
                    "j": r"\mathbf j", "-j": r"\underline{\mathbf j}", "k": r"\mathbf k", "-k": r"\underline{\mathbf k}"}, "Q")


D2D = FiniteGroup({"1": sp.eye(3), "R^2": D4.mats["R^2"], "IR": INV3 * D4.mats["R"], "IR^3": INV3 * D4.mats["R^3"],
                   "IR1^2": INV3 * D4.mats["R1^2"], "IR2^2": INV3 * D4.mats["R2^2"], "i3": D4.mats["i3"], "i4": D4.mats["i4"]},
                  label="D2d")                                  # (3.6.32)
D2D_TO_D4 = {"1": "1", "R^2": "R^2", "IR": "R", "IR^3": "R^3", "IR1^2": "R1^2", "IR2^2": "R2^2", "i3": "i3", "i4": "i4"}

# Fig. 3.1.2：三個擺錘（a 在 0°、b 在 120°、c 在 240°）凍結在某個非平衡位置（單位：環的半徑 1）
BOB_START = {"a": (1.0, 0.28), "b": (-0.72, 0.62), "c": (-0.36, -1.02)}


def config_after(g: str, pts: dict | None = None, G: FiniteGroup | None = None) -> dict:
    """群元素 g（鏡面固定在紙上）作用在擺錘的位置上：{名字: g·p}（浮點 2D）。"""
    G = G or C3V
    pts = pts or BOB_START
    M = np.array(G.mats[g].evalf(), dtype=float)[:2, :2]
    return {k: M @ np.asarray(v, float) for k, v in pts.items()}


def fig312() -> dict:
    """Fig. 3.1.2 的兩列：σ₂σ₁（先 σ₁）、σ₁σ₂（先 σ₂）的中間與最後位置；最後一格＝ r²、r 作用的結果。"""
    mid1, mid2 = config_after("s1"), config_after("s2")
    return {"s2s1": (mid1, config_after("s2", mid1), config_after("r^2")),
            "s1s2": (mid2, config_after("s1", mid2), config_after("r"))}


def book_order(G: FiniteGroup, ref_rows) -> tuple:
    """(ℓ 列表, 特徵標表) 依書上的列順序（ref_rows）。"""
    ells, T = G.char_table_class_algebra()
    perm = reorder(T, ref_rows)
    return [ells[a] for a in perm], sp.Matrix([list(T.row(a)) for a in perm]), perm


# ======================================================================== §3.1 Hamilton 轉、實驗室與本體座標
def reflection3(N) -> np.ndarray:
    """法向量 N 的鏡面反射 1 − 2NNᵀ（浮點）。"""
    N = np.asarray(N, float)
    N = N / np.linalg.norm(N)
    return np.eye(3) - 2 * np.outer(N, N)


def quat(axis, deg) -> np.ndarray:
    """繞 axis 轉 deg 度的單位四元數 (cos θ/2, sin θ/2 n̂)。"""
    n = np.asarray(axis, float)
    n = n / np.linalg.norm(n)
    t = math.radians(deg) / 2
    return np.array([math.cos(t), *(math.sin(t) * n)])


def qmul(p, q) -> np.ndarray:
    """Hamilton 乘積 pq（先 q 後 p 的轉動）。"""
    p0, pv, q0, qv = p[0], np.asarray(p[1:]), q[0], np.asarray(q[1:])
    return np.array([p0 * q0 - pv @ qv, *(p0 * qv + q0 * pv + np.cross(pv, qv))])


def quat_to_rot(q) -> np.ndarray:
    """單位四元數 → 3×3 轉動矩陣。"""
    w, x, y, z = q
    return np.array([[1 - 2 * (y * y + z * z), 2 * (x * y - w * z), 2 * (x * z + w * y)],
                     [2 * (x * y + w * z), 1 - 2 * (x * x + z * z), 2 * (y * z - w * x)],
                     [2 * (x * z - w * y), 2 * (y * z + w * x), 1 - 2 * (x * x + y * y)]])


def su2(q) -> np.ndarray:
    """四元數 → SU(2)：U = q₀1 − i(q₁σ_x + q₂σ_y + q₃σ_z)（第 1 章 (1.1.3) 的半角）。"""
    sx, sy, sz = (np.array(M.tolist(), dtype=complex) for M in (_sx, _sy, _sz))
    return q[0] * np.eye(2) - 1j * (q[1] * sx + q[2] * sy + q[3] * sz)


def axis_angle(Rm) -> tuple:
    """3×3 轉動 → (單位軸, 角度 0–180°)。"""
    Rm = np.asarray(Rm, float)
    ang = math.degrees(math.acos(max(-1.0, min(1.0, (np.trace(Rm) - 1) / 2))))
    if ang < 1e-9:
        return np.array([0.0, 0.0, 1.0]), 0.0
    if abs(ang - 180) < 1e-7:
        w, V = np.linalg.eigh((Rm + Rm.T) / 2)
        return V[:, np.argmax(w)], 180.0
    v = np.array([Rm[2, 1] - Rm[1, 2], Rm[0, 2] - Rm[2, 0], Rm[1, 0] - Rm[0, 1]])
    return v / np.linalg.norm(v), ang


def rotate_vec(v, axis, deg) -> np.ndarray:
    return quat_to_rot(quat(axis, deg)) @ np.asarray(v, float)


def hamilton_arc(axis, deg, start) -> tuple:
    """Hamilton 轉向量（Fig. 3.1.5）：從 N1 = start（⊥ 軸、單位）沿大圓轉 deg/2 到 N2；R[ω] = σ_{N2} σ_{N1}。"""
    N1 = np.asarray(start, float)
    N1 = N1 / np.linalg.norm(N1)
    return N1, rotate_vec(N1, axis, deg / 2)


def arc_points(N1, N2, n: int = 30, radius: float = 1.0) -> np.ndarray:
    """N1 → N2 的大圓弧上 n 點（短弧，球半徑 radius）。"""
    N1, N2 = np.asarray(N1, float), np.asarray(N2, float)
    om = math.acos(max(-1.0, min(1.0, N1 @ N2)))
    if om < 1e-12:
        return np.array([N1 * radius] * n)
    t = np.linspace(0, 1, n)[:, None]
    return radius * (np.sin((1 - t) * om) * N1 + np.sin(t * om) * N2) / math.sin(om)


def arc_quat(N1, N2) -> np.ndarray:
    """弧 N1 → N2 所代表的轉動的四元數 (N1·N2, N1×N2)（轉 2∠(N1, N2)，軸 N1×N2）。"""
    N1, N2 = np.asarray(N1, float), np.asarray(N2, float)
    return np.array([N1 @ N2, *np.cross(N1, N2)])


def hamilton_product(axis1, deg1, axis2, deg2) -> dict:
    """R[ω'']= R[ω']R[ω]（先 R[ω] = (axis1, deg1)，Fig. 3.1.6）：兩段半角弧在兩個大圓的交點 P 頭尾相接，
    第一段 N1 → P、第二段 P → N2'；乘積的弧 N1 → N2'。回傳 dict(N1, P, N2p, axis, deg（由弧）, q（由四元數）, q_arc)。"""
    a1, a2 = np.asarray(axis1, float), np.asarray(axis2, float)
    a1, a2 = a1 / np.linalg.norm(a1), a2 / np.linalg.norm(a2)
    P = np.cross(a1, a2)
    if np.linalg.norm(P) < 1e-12:                               # 同軸：取任一 ⊥ 向量
        P = np.cross(a1, [1, 0, 0] if abs(a1[0]) < 0.9 else [0, 1, 0])
    P = P / np.linalg.norm(P)
    N1 = rotate_vec(P, a1, -deg1 / 2)
    N2p = rotate_vec(P, a2, deg2 / 2)
    q = qmul(quat(a2, deg2), quat(a1, deg1))
    qa = arc_quat(N1, N2p)
    ax, ang = axis_angle(quat_to_rot(qa))
    return dict(N1=N1, P=P, N2p=N2p, axis=ax, deg=ang, q=q, q_arc=qa)


def conjugate_rotation(T_axis, T_deg, R_axis, R_deg) -> tuple:
    """R' = T R T⁻¹（(3.2.4)、Fig. 3.2.1）：同一個轉角，軸換成 T·n̂。回傳 (四元數算的軸, 角, T·n̂)。"""
    qT = quat(T_axis, T_deg)
    q = qmul(qmul(qT, quat(R_axis, R_deg)), qT * np.array([1, -1, -1, -1]))
    ax, ang = axis_angle(quat_to_rot(q))
    return ax, ang, rotate_vec(R_axis, T_axis, T_deg)


def body_product(G: FiniteGroup, first: str, second: str) -> str:
    """本體座標（§3.1C、Fig. 3.1.8）：先做 first，鏡面跟著物體走；再做「跟著走過的」second，即 T·M₂·T⁻¹（T = first）。
    結果的名稱（＝實驗室的 first·second：順序反過來）。"""
    T = G.num[first]
    total = (T @ G.num[second] @ np.linalg.inv(T)) @ T
    return next(g for g in G.names if _close(G.num[g], total))


def d3_arcs() -> dict:
    """Fig. 3.1.7(a)：D₃ 的 Hamilton 弧（單位球）：r：赤道上 60° 弧；ρ_j：過極點、90° 弧（⊥ 轉軸的大圓）。"""
    out = {"r": hamilton_arc((0, 0, 1), 120, (1, 0, 0)), "r^2": hamilton_arc((0, 0, 1), 240, (1, 0, 0))}
    for p, a in D3_AXIS_DEG.items():
        ax = np.array([math.cos(math.radians(a)), math.sin(math.radians(a)), 0.0])
        out[p] = hamilton_arc(ax, 180, (0, 0, 1))
    return out


# ======================================================================== §3.2–3.4 C3v 的冪等元、基本算子、正則表示
C3V_CLASS_IDEMP_BOOK = [[R_(1, 6), R_(1, 6), R_(1, 6)], [R_(1, 6), R_(1, 6), -R_(1, 6)], [R_(2, 3), -R_(1, 3), 0]]  # (3.5.4)ₓ


def c3v_class_spectral() -> sp.Matrix:
    """(3.2.17)：c_j = Σ_α c_j^{(α)} ℙ^α 的係數 c_j^{(α)}（列 α、行 j）＝ °c_j χ^α_j / ℓ^α。"""
    ells, T = C3V.char_table_class_algebra()
    sz = C3V.class_sizes()
    return sp.Matrix(T.rows, T.cols, lambda a, j: sz[j] * T[a, j] / ells[a])


def c3v_split_cv() -> tuple:
    """(3.3.1)–(3.3.3)：用 C_v'' = {1, σ₃} 的 P^± = (1 ± σ₃)/2 分裂 ℙ³ ⇒ (P³₁, P³₂)。"""
    P3 = C3V.central_idempotents[2]
    return tuple(C3V.split(P3, C3V.cyclic_projectors("s3")))


def c3v_elementary(kind: str = "cv") -> dict:
    """C3v 的 E 型基本算子 {(i, j): 群代數}（0 起算）：kind="cv" 用 {1, σ₃}、參考元素 σ₂（(3.3.7)、(3.4.8)）；
    kind="c3" 用 C₃ = {1, r, r²}（(3.3.18)；本系列 P^{k₃} 讓 r P^{k₃} = ε^k P^{k₃}）、參考元素 σ₃（(3.4.24c)：𝒟^{3°}(σ₃) = [[0,1],[1,0]]）。"""
    P3 = C3V.central_idempotents[2]
    if kind == "cv":
        Ps = C3V.split(P3, C3V.cyclic_projectors("s3"))
        return C3V.elementary(Ps, "s2")[0]
    Ps = C3V.split(P3, C3V.cyclic_projectors("r"))
    return C3V.elementary(Ps, "s3")[0]


def c3v_irreps(kind: str = "cv") -> list:
    """[D^{A1}, D^{A2}, D^E]（{g: Matrix}）：E 由基本算子求得（kind 見 c3v_elementary）。"""
    ones = [{g: sp.Matrix([[C3V_CHARS_BOOK[a][C3V.class_of(g)]]]) for g in C3V.names} for a in (0, 1)]
    return ones + [C3V.irrep_from_elementary(c3v_elementary(kind))]


T_CIRC = sp.Matrix([[1, 1], [-I, I]]) / sp.sqrt(2)            # (3.4.24b)


def c3v_regular_T() -> sp.Matrix:
    """(3.4.7)–(3.4.8)、(3.4.16)：行 = |P̂¹⟩、|P̂²⟩、|P̂³₁₁⟩、|P̂³₂₁⟩、|P̂³₁₂⟩、|P̂³₂₂⟩（N^α = ℓ^α/°G，(3.4.20)）。"""
    cols = [C3V.vec_of(C3V.central_idempotents[0]) * sp.sqrt(6), C3V.vec_of(C3V.central_idempotents[1]) * sp.sqrt(6)]
    Pij = c3v_elementary("cv")
    for i, j in [(0, 0), (1, 0), (0, 1), (1, 1)]:
        cols.append(C3V.vec_of(Pij[i, j]) * sp.sqrt(3))
    return sp.Matrix.hstack(*cols).applyfunc(sp.nsimplify)


GROUP_VECTORS_BOOK = {                                          # (3.4.8)：分量依 1, r, r², σ₁, σ₂, σ₃
    "P1": (sp.Matrix([1, 1, 1, 1, 1, 1]) / sp.sqrt(6)), "P2": (sp.Matrix([1, 1, 1, -1, -1, -1]) / sp.sqrt(6)),
    "P11": sp.Matrix([2, -1, -1, -1, -1, 2]) / (2 * SQ3), "P22": sp.Matrix([2, -1, -1, 1, 1, -2]) / (2 * SQ3),
    "P12": sp.Matrix([0, -1, 1, -1, 1, 0]) / 2, "P21": sp.Matrix([0, 1, -1, -1, 1, 0]) / 2}


# ======================================================================== §3.5 C3v 的三擺：頻率公式、跡公式求本徵值
PEND3 = SA.RING3                                               # a = b = 1（第 2 章 H206 的數字）


def pendulum_rep(g: str) -> sp.Matrix:
    """三擺（擺錘在 0°、120°、240°，徑向擺動）的 3×3 置換表示：σ₃|1⟩ = |1⟩、σ₂|1⟩ = |r⟩、σ₁|1⟩ = |r²⟩（p. 189）。"""
    pos = [rot2(120 * p) * sp.Matrix([1, 0]) for p in range(3)]
    M = sp.zeros(3)
    for j, v in enumerate(pos):
        w = C3V.mats[g] * v
        i = next(i for i, u in enumerate(pos) if sp.simplify(u - w) == sp.zeros(2, 1))
        M[i, j] = 1
    return M


def pendulum_frequencies() -> list:
    """(3.5.11)ₓ：三擺座標的 f^{A1}, f^{A2}, f^E（跡 3、0、1）。"""
    tr = [pendulum_rep(C[0]).trace() for C in C3V.classes]
    return C3V.frequencies(tr)


def pendulum_trace_eigs(a=PEND3["a"], b=PEND3["b"]) -> list:
    """(3.5.13)–(3.5.14)：a^α = (1/°G) Σ_類 °c_g χ^α_g* Tr⟨a g⟩（只用第一列）：A1 → b、E → 3a + b（A2 頻率 0，回傳 None）。"""
    A = SA.ring_matrix(3, a, b)
    ells, T = C3V.char_table_class_algebra()
    sz = C3V.class_sizes()
    f = pendulum_frequencies()
    out = []
    for al in range(3):
        if f[al] == 0:
            out.append(None)
            continue
        s = sum(sz[k] * T[al, k] * (A * pendulum_rep(C[0])).trace() for k, C in enumerate(C3V.classes))
        out.append(sp.nsimplify(s / (C3V.order * f[al])))
    return out


# ======================================================================== §3.4 三角形六自由度彈簧質量模型（Fig. 3.4.1）
SIX = dict(k1=4, k2=8, m=1, alpha=15)          # 本系列的數字；α = 15°：k₂ 彈簧沿三角形邊的延長線（Fig. 3.4.1）
SIX_ANGLES = {"1": -45, "r": 75, "r^2": 195, "s1": 285, "s2": 165, "s3": 45}   # 座標 |g⟩ 的方向（度）＝ g·u₁，u₁ = e(−45°)
SIX_MASS = {"1": 0, "s3": 0, "r": 1, "s2": 1, "r^2": 2, "s1": 2}             # 座標 |g⟩ 在哪個質點（0°、120°、240°）


def six_dof_first_row(k1=SIX["k1"], k2=SIX["k2"], alpha=SIX["alpha"], m=SIX["m"]) -> list:
    """(3.4.21b)：⟨1|a|g⟩，g = 1, r, r², σ₁, σ₂, σ₃（精確）。"""
    k1, k2, m = sp.nsimplify(k1), sp.nsimplify(k2), sp.nsimplify(m)
    s = sind(2 * sp.nsimplify(alpha))
    return [(k1 + k2) / m, k1 / (4 * m), k1 / (4 * m), k1 * (2 - SQ3) / (4 * m), k1 * (2 + SQ3) / (4 * m),
            (k1 + 2 * k2 * s) / (2 * m)]


def six_dof_matrix(k1=SIX["k1"], k2=SIX["k2"], alpha=SIX["alpha"], m=SIX["m"]) -> sp.Matrix:
    """6×6 加速度矩陣：⟨g|a|h⟩ = ⟨1|a|g⁻¹h⟩（a 與所有 g 對易，只要第一列）。"""
    row = dict(zip(C3V.names, six_dof_first_row(k1, k2, alpha, m)))
    return sp.Matrix(6, 6, lambda i, j: row[C3V.mul(C3V.inv(C3V.names[i]), C3V.names[j])])


def six_dof_geometric(k1=SIX["k1"], k2=SIX["k2"], alpha=SIX["alpha"], m=SIX["m"]) -> np.ndarray:
    """同一個矩陣，直接由彈簧幾何算（浮點）：k₁ 沿三邊、k₂ 從每個質點往外牆，與徑向夾 45° − α；座標方向 SIX_ANGLES。"""
    u = {g: e2(a) for g, a in SIX_ANGLES.items()}
    pos = [e2(120 * p) for p in range(3)]
    K = np.zeros((6, 6))
    names = C3V.names
    for a in range(3):
        for b in range(a + 1, 3):
            n = (pos[b] - pos[a]) / np.linalg.norm(pos[b] - pos[a])
            for i, g in enumerate(names):
                for j, h in enumerate(names):
                    ma, mb = SIX_MASS[g], SIX_MASS[h]
                    if {ma, mb} == {a, b}:
                        K[i, j] -= k1 * (u[g] @ n) * (u[h] @ n)
                    elif ma == mb and ma in (a, b):
                        K[i, j] += k1 * (u[g] @ n) * (u[h] @ n)
    beta = 45 - float(alpha)
    for i, g in enumerate(names):
        for j, h in enumerate(names):
            if SIX_MASS[g] == SIX_MASS[h]:
                for sg in (1, -1):
                    n = e2(120 * SIX_MASS[g] + sg * beta)
                    K[i, j] += k2 * (u[g] @ n) * (u[h] @ n)
    return K / float(m)


def six_dof_blocks(k1=SIX["k1"], k2=SIX["k2"], alpha=SIX["alpha"], m=SIX["m"]) -> tuple:
    """(3.4.22)：a^(1) = Σ_g ⟨1|a|g⟩、a^(2) = Σ χ^{A2}(g)⟨1|a|g⟩、a^(3)_{jl} = Σ_g D³_{jl}(g)*⟨1|a|g⟩（2×2，兩份相同）。"""
    row = six_dof_first_row(k1, k2, alpha, m)
    irr = c3v_irreps("cv")
    out = []
    for D in irr:
        M = sum((sp.conjugate(D[g]) * row[i] for i, g in enumerate(C3V.names)), sp.zeros(*D["1"].shape))
        out.append(M.applyfunc(lambda x: sp.nsimplify(sp.simplify(x))))
    return out[0][0, 0], out[1][0, 0], out[2]


def six_dof_blocks_book(k1=SIX["k1"], k2=SIX["k2"], alpha=SIX["alpha"], m=SIX["m"]) -> tuple:
    """(3.4.23b) 照抄。"""
    k1, k2, m = sp.nsimplify(k1), sp.nsimplify(k2), sp.nsimplify(m)
    s = sind(2 * sp.nsimplify(alpha))
    return ((3 * k1 + k2 * (1 + s)) / m, k2 * (1 - s) / m,
            sp.Matrix([[3 * k1 / 4 + k2 * (1 + s), 3 * k1 / 4], [3 * k1 / 4, 3 * k1 / 4 + k2 * (1 - s)]]) / m)


def six_dof_omega2(k1=SIX["k1"], k2=SIX["k2"], alpha=SIX["alpha"], m=SIX["m"]) -> dict:
    """{"A1": ω², "A2": ω², "E": (ω²₊, ω²₋)}（精確）。"""
    a1, a2, a3 = six_dof_blocks(k1, k2, alpha, m)
    ev = sorted(a3.eigenvals().keys(), key=lambda x: -float(x))
    return {"A1": a1, "A2": a2, "E": tuple(sp.nsimplify(sp.simplify(x)) for x in ev)}


def six_dof_mode(k1=SIX["k1"], k2=SIX["k2"], alpha=SIX["alpha"], m=SIX["m"], which: int = 0) -> tuple:
    """E 模（which = 0：高頻 ω₊、1：低頻 ω₋）的兩個夥伴 (|x⟩, |y⟩, ω)（浮點 6 維，正則座標）：
    |x⟩ = c₁|P̂³₁₁⟩ + c₂|P̂³₁₂⟩（型 1）、|y⟩ = c₁|P̂³₂₁⟩ + c₂|P̂³₂₂⟩（型 2），(c₁, c₂) 是 a^(3) 的本徵向量（Fig. 3.4.4）。"""
    a3 = np.array(six_dof_blocks(k1, k2, alpha, m)[2].evalf(), dtype=float)
    w, V = np.linalg.eigh(a3)
    idx = np.argsort(-w)[which]
    c = V[:, idx] * (1 if V[0, idx] >= 0 else -1)
    T = np.array(c3v_regular_T().evalf(), dtype=float)
    x = c[0] * T[:, 2] + c[1] * T[:, 4]
    y = c[0] * T[:, 3] + c[1] * T[:, 5]
    return x, y, math.sqrt(max(w[idx], 0.0))


def six_dof_displacements(x6) -> np.ndarray:
    """6 個正則座標 → 三個質點的 2D 位移（列：0°、120°、240° 的質點）。"""
    out = np.zeros((3, 2))
    for i, g in enumerate(C3V.names):
        out[SIX_MASS[g]] += float(np.real(x6[i])) * e2(SIX_ANGLES[g])
    return out


def six_dof_circular(k1=SIX["k1"], k2=SIX["k2"], alpha=SIX["alpha"], m=SIX["m"], which: int = 0, n: int = 24) -> dict:
    """(3.4.24)–(3.4.25)、Fig. 3.4.5：x(t) = Re[e^{−iωt}(|x⟩ + i|y⟩)] = cos ωt|x⟩ + sin ωt|y⟩。
    回傳 dict(rmin, rmax（三個質點各自的最小、最大距離；k₂ = 0 時相等＝圓）, sense=[+1 逆時針／−1 順時針], omega)。
    一般 k₂：每個質點走橢圓（徑向半軸 ∝ c₁、切向 ∝ c₂），三個橢圓一樣大、互相轉 120°。"""
    x, y, om = six_dof_mode(k1, k2, alpha, m, which)
    ts = np.linspace(0, 2 * math.pi, n, endpoint=False)
    D = np.array([six_dof_displacements(math.cos(t) * x + math.sin(t) * y) for t in ts])   # (n, 3, 2)
    r = np.linalg.norm(D, axis=2)
    sense = [int(np.sign(np.cross(D[0, p], D[1, p]))) for p in range(3)]
    return dict(rmin=r.min(axis=0), rmax=r.max(axis=0), sense=sense, omega=om)


def fig344_vectors(which: str = "fixed") -> dict:
    """Fig. 3.4.4（k₂ = 0）的六個駐波向量（分量依 1, r, r², σ₁, σ₂, σ₃）：which="book" 照印、"fixed" 更正（型 2 兩欄）。"""
    s3 = SQ3
    v = lambda *a: sp.Matrix(a) / (2 * sp.sqrt(6))             # noqa: E731
    out = {"P1": GROUP_VECTORS_BOOK["P1"], "P2": GROUP_VECTORS_BOOK["P2"],
           "x+": v(2, -1 - s3, -1 + s3, -1 - s3, -1 + s3, 2), "x-": v(2, -1 + s3, -1 - s3, -1 + s3, -1 - s3, 2)}
    if which == "book":
        out.update({"y+": v(2, -1 + s3, -1 - s3, -1 - s3, -1 + s3, 2), "y-": v(2, -1 + s3, -1 - s3, -1 - s3, -1 - s3, 2)})
    else:
        out.update({"y+": v(2, -1 + s3, -1 - s3, 1 - s3, 1 + s3, -2), "y-": v(-2, 1 + s3, 1 - s3, -1 - s3, -1 + s3, 2)})
    return out


# ======================================================================== §3.6 D₄、D₆ 位能裡的轉子（Bohr 波與 Bloch 波）
ROTOR = dict(weak=2, strong4=16, strong6=36, zeeman_b=R_(1, 10))   # V(φ) = −V cos nφ（E = ħ²/2μr² = 1；谷在 φ = 0）


def dn_label_free(n: int, m: int, kind: str) -> str:
    """自由轉子的駐波 cos mφ（kind="c"）／sin mφ（"s"）的 D_n 標記（由 m mod n）：
    m ≡ 0：cos → A1、sin → A2；n 偶且 m ≡ n/2：cos → B1（波腹在谷）、sin → B2；其他 → E_{min(κ, n−κ)}（n = 4 寫 E）。"""
    kap = m % n
    if kap == 0:
        return "A_1" if kind == "c" else "A_2"
    if n % 2 == 0 and kap == n // 2:
        return "B_1" if kind == "c" else "B_2"
    e = min(kap, n - kap)
    return "E" if n <= 4 else f"E_{e}"


def dn_free_levels(n: int, count: int) -> list:
    """自由轉子由低到高的前 count 個態 [(能量 m², 標記)]（E 二重算兩個）。"""
    out = [(0, "A_1")]
    m = 1
    while len(out) < count:
        a, b = dn_label_free(n, m, "c"), dn_label_free(n, m, "s")
        out += [(m * m, a), (m * m, b)]
        m += 1
    return out[:count]


def rotor_dn(n: int, V, mmax: int = 40) -> tuple:
    """D_n 位能 −V cos nφ 的轉子：(本徵值, 本徵向量(行), ms)（平面波基底，借 sym_abel.rotor_matrix）。"""
    ms, H = SA.rotor_matrix(V, mmax, n)
    w, U = np.linalg.eigh(H)
    return w, U, ms


def dn_label_vector(n: int, c, ms) -> str:
    """由數值本徵向量的對稱定標記：κ = m mod n（R：φ → φ + 2π/n）；κ ≡ 0、n/2 時看 φ → −φ（谷的 ρ）的宇稱。"""
    c = np.asarray(c)
    wts = {}
    for coef, m in zip(c, ms):
        k = min(m % n, (-m) % n)
        wts[k] = wts.get(k, 0) + abs(coef) ** 2
    kap = max(wts, key=wts.get)
    if kap != 0 and not (n % 2 == 0 and kap == n // 2):
        return "E" if n <= 4 else f"E_{kap}"
    idx = {m: i for i, m in enumerate(ms)}
    par = sum(c[i] * c[idx[-m]] for i, m in enumerate(ms)) / sum(abs(c) ** 2)
    if kap == 0:
        return "A_1" if par > 0 else "A_2"
    return "B_1" if par > 0 else "B_2"


def dn_levels(n: int, V, count: int = 8, mmax: int = 40) -> list:
    """最低 count 個能階 [(能量, 由本徵向量定的標記)]。"""
    w, U, ms = rotor_dn(n, V, mmax)
    return [(float(w[i]), dn_label_vector(n, U[:, i], ms)) for i in range(count)]


def dn_level_paths(n: int, Vs, count: int = 8, mmax: int = 40) -> dict:
    """能階隨 V 變化（給「弱 → 強」的連線）：array[len(Vs), count]（第 i 條＝由低到高第 i 個），標記 labels[len(Vs)][count]。
    最低 2n 個態在 0 ≤ V ≤ 強位能之間不交叉，標記沿著每一條不變（check() 核對）。回傳 (levels, labels)。"""
    lev = [dn_levels(n, V, count, mmax) for V in Vs]
    return np.array([[e for e, _ in row] for row in lev]), [[lab for _, lab in row] for row in lev]


def dn_gap_2x2(n: int, V) -> tuple:
    """弱位能：m = ±n/2 的 2×2 [[(n/2)², G], [G, (n/2)²]]，G = −V/2 ⇒ (B1, B2) = ((n/2)² + G, (n/2)² − G)，能隙 2|G| = V。"""
    G = -sp.nsimplify(V) / 2
    e0 = sp.Integer((n // 2) ** 2)
    return e0 + G, e0 - G


def dn_density(n: int, V, label: str, phis) -> np.ndarray:
    """標記為 label 的（最低的那個）態的 |ψ(φ)|²（給「B1 波腹在谷、B2 在山頂」）。"""
    w, U, ms = rotor_dn(n, V)
    i = next(i for i in range(len(w)) if dn_label_vector(n, U[:, i], ms) == label)
    psi = np.array([np.sum(U[:, i] * np.exp(1j * ms * p)) for p in np.atleast_1d(phis)])
    return np.abs(psi) ** 2


def dn_zeeman(n: int, V, b, count: int = 4) -> list:
    """加磁場：H = (m − b)² + 位能（C_n 對稱，ρ 不再是對稱）；每個 C_n 不可約表示 κ = m mod n 各自對角化。
    回傳最低 count 個 [(能量, κ)]（(3.6.14)：E → (1_n) ⊕ ((n−1)_n)）。"""
    out = []
    for kap in range(n):
        msk = np.array([kap + n * t for t in range(-12, 13)])
        H = np.diag((msk - float(b)) ** 2)
        for i in range(len(msk) - 1):
            H[i, i + 1] = H[i + 1, i] = -float(V) / 2
        out += [(float(e), kap) for e in np.linalg.eigvalsh(H)[:3]]
    return sorted(out)[:count]


# ======================================================================== §3.7 氨 NH₃
NH3 = dict(m=1, M=14, D=20)                         # H、N 的質量（amu）與角度 D（N–H 鍵與 H₃ 平面的夾角，度）
NH3_NU = dict(p1=3337, m1=908, p3=3414, m3=1628)    # Herzberg（書 p. 214、Fig. 3.7.3）：ν₁₊、ν₁₋、ν₃₊、ν₃₋（cm⁻¹）
NH3_OMEGA_BOOK = dict(p1=6.290, m1=1.712, p3=6.435, m3=3.069)    # 書上的 ω（10¹⁴ rad/s；用 c = 3×10¹⁰ cm/s）
C_LIGHT = 2.99792458e10                              # cm/s
M_H_BOOK = 1.67e-27                                  # kg（書 p. 216）
NH3_KJ = 0.3                                         # Fig. 3.7.4 的工作點 k/j（D = 20°）


def nh3_geometry(D=NH3["D"]) -> tuple:
    """鍵長 1：H_p = (cos D cos 120p°, cos D sin 120p°, 0)、N = (0, 0, sin D)。
    座標（Fig. 3.7.1）：|H1⟩ = 沿 N→H 鍵 (cos D, 0, −sin D)、|H2⟩ = 切向 ŷ、|H3⟩ ⊥ 鍵向上 (sin D, 0, cos D)，
    r^p|H_i⟩ 轉 120p°；|N1⟩ = x̂、|N2⟩ = ŷ、|N3⟩ = ẑ。回傳 (positions 4×3, [(原子, 單位向量)]×12)（順序 H1,rH1,r²H1,H2,…,N1,N2,N3）。"""
    d = math.radians(float(D))
    c, s = math.cos(d), math.sin(d)
    Rz = lambda p: np.array(quat_to_rot(quat((0, 0, 1), 120 * p)))   # noqa: E731
    pos = np.array([Rz(p) @ np.array([c, 0, 0]) for p in range(3)] + [np.array([0, 0, s])])
    us = [np.array([c, 0, -s]), np.array([0.0, 1, 0]), np.array([s, 0, c])]
    coords = [(p, Rz(p) @ u) for u in us for p in range(3)] + [(3, np.eye(3)[i]) for i in range(3)]
    return pos, coords


def nh3_force(k, j, D=NH3["D"]) -> np.ndarray:
    """12×12 力常數矩陣（H–H 彈簧 k、N–H 彈簧 j，無預應力：Σ k n nᵀ）＝ (3.7.13a)。"""
    pos, coords = nh3_geometry(D)
    springs = [(0, 1, k), (1, 2, k), (0, 2, k), (0, 3, j), (1, 3, j), (2, 3, j)]
    F = np.zeros((12, 12))
    for a, b, kk in springs:
        n = (pos[b] - pos[a]) / np.linalg.norm(pos[b] - pos[a])
        for i, (pa, ua) in enumerate(coords):
            sa = 1 if pa == a else (-1 if pa == b else 0)
            if not sa:
                continue
            for jj, (pb, ub) in enumerate(coords):
                sb = 1 if pb == a else (-1 if pb == b else 0)
                if sb:
                    F[i, jj] += kk * sa * sb * (ua @ n) * (ub @ n)
    return F


def nh3_mass(m=NH3["m"], M=NH3["M"]) -> np.ndarray:
    return np.diag([float(m)] * 9 + [float(M)] * 3)


def nh3_rep(g: str, D=NH3["D"]) -> np.ndarray:
    """C3v 在 12 個座標上的表示（正交 12×12）：g 作用在原子位置與向量上（σ_j 是含 z 軸、鏡線 MIRROR_DEG 的鉛直面）。"""
    pos, coords = nh3_geometry(D)
    G3 = np.array(block3(C3V.mats[g]).evalf(), dtype=float)
    M = np.zeros((12, 12))
    for jj, (pa, ua) in enumerate(coords):
        pb = next(q for q in range(4) if np.allclose(G3 @ pos[pa], pos[q]))
        v = G3 @ ua
        for i, (pc, uc) in enumerate(coords):
            if pc == pb:
                M[i, jj] = uc @ v
    return M


def nh3_projector(alpha: int, i: int, j: int, D=NH3["D"]) -> np.ndarray:
    """P^α_{ij} = (ℓ^α/°G) Σ_g D^α_{ij}(g)* 𝓡(g)（(3.4.19)；α = 0, 1, 2 ↔ A1, A2, E；i, j 0 起算）。"""
    D_ = c3v_irreps("cv")[alpha]
    ell = D_["1"].shape[0]
    return sum(float(sp.re(sp.conjugate(D_[g][i, j]))) * nh3_rep(g, D) for g in C3V.names) * ell / 6


def nh3_cvectors(m=NH3["m"], M=NH3["M"], D=NH3["D"]) -> dict:
    """(3.7.8)–(3.7.11) 的 |c¹⟩、|c¹'⟩、|c³₁⟩（圖上的歸一化：(3.7.11) 再除 √2）、|c³'₁⟩（12 維浮點）。"""
    d = math.radians(float(D))
    e = np.eye(12)
    P1H1, P1H3 = (e[0] + e[1] + e[2]) / math.sqrt(3), (e[6] + e[7] + e[8]) / math.sqrt(3)
    P11H1, P11H3 = (2 * e[0] - e[1] - e[2]) / math.sqrt(6), (2 * e[6] - e[7] - e[8]) / math.sqrt(6)
    P12H2 = (-e[4] + e[5]) / math.sqrt(2)
    c1 = M * P1H1 + math.sqrt(3) * m * math.sin(d) * e[11]
    c1p = M * P1H3 - math.sqrt(3) * m * math.cos(d) * e[11]
    c3 = (math.cos(d) * P11H1 - P12H2 + math.sin(d) * P11H3) / math.sqrt(2)
    c3p = M * P11H1 - 3 * m * math.cos(d) * e[9] / math.sqrt(6)
    return {"c1": c1, "c1'": c1p, "c3": c3, "c3'": c3p}


def nh3_constraints(v, m=NH3["m"], M=NH3["M"], D=NH3["D"]) -> tuple:
    """(3.7.12)：T = Σ m c、A = Σ m r × c（v：12 維位移）。"""
    pos, coords = nh3_geometry(D)
    disp = np.zeros((4, 3))
    for x, (p, u) in zip(v, coords):
        disp[p] += x * u
    ms = np.array([m, m, m, M], float)
    T = (ms[:, None] * disp).sum(axis=0)
    A = sum(ms[p] * np.cross(pos[p], disp[p]) for p in range(4))
    return T, A


def nh3_reduced(k, j, m=NH3["m"], M=NH3["M"], D=NH3["D"]) -> dict:
    """{"A1": (F, m, a), "E": (F, m, a)}：在 {|c¹⟩, |c¹'⟩}、{|c³₁⟩, |c³'₁⟩} 基底的 2×2（數值，由 12×12 直接投影）。"""
    F, Mm = nh3_force(k, j, D), nh3_mass(m, M)
    cv = nh3_cvectors(m, M, D)
    out = {}
    for key, (u, v) in {"A1": (cv["c1"], cv["c1'"]), "E": (cv["c3"], cv["c3'"])}.items():
        B = np.vstack([u, v]).T
        Fr, Mr = B.T @ F @ B, B.T @ Mm @ B
        out[key] = (Fr, Mr, np.linalg.solve(Mr, Fr))
    return out


def nh3_reduced_formulas(k, j, m=NH3["m"], M=NH3["M"], D=NH3["D"], book: bool = False) -> dict:
    """(3.7.20)–(3.7.24)、(3.7.27) 的 2×2：book=True 照印，False 更正後。回傳同 nh3_reduced 的格式。"""
    d = math.radians(float(D))
    c, s = math.cos(d), math.sin(d)
    H, K, J = 1 + 1.5 * m / M * c * c, 1 + 3 * m / M * s * s, 1 + 3 * m / M * c * c
    r2 = math.sqrt(2)
    if book:
        F3 = np.array([[1.5 * k + j / 2 * c * c, M * c * (1.5 * k + j * H)], [M * c * (1.5 * k + j * H), M * M * (0.75 * k * c * c + j * H * H)]])
        F1 = np.array([[M * M * (3 * k * c * c + j * J), 3 * J * M * M * (k - j * m / M) * s * c],
                       [3 * J * M * M * (k - j * m / M) * s * c, M * M * (3 * k + j * (2 * m / M) ** 2 * c * c) * s * s]])
        m1 = np.array([[m * M * M * K, -3 * m / M * s * c], [-3 * m / M * s * c, J]])
        a3 = np.array([[3 * k / m, 3 * k / (m * math.sqrt(8)) * M * c], [j / (m * M * math.sqrt(8)) * c, H * j / m]])
    else:
        F3 = np.array([[1.5 * k + j / 2 * c * c, M / r2 * c * (1.5 * k + j * H)], [M / r2 * c * (1.5 * k + j * H), M * M * (0.75 * k * c * c + j * H * H)]])
        F1 = np.array([[M * M * (3 * k * c * c + j * K * K), 3 * M * M * (k - j * m * K / M) * s * c],
                       [3 * M * M * (k - j * m * K / M) * s * c, M * M * (3 * k + j * (3 * m / M) ** 2 * c * c) * s * s]])
        m1 = m * M * M * np.array([[K, -3 * m / M * s * c], [-3 * m / M * s * c, J]])
        a3 = np.array([[1.5 * k / m, 3 * k / (m * math.sqrt(8)) * M * c], [j / (m * M * r2) * c, H * j / m]])
    m3 = np.array([[m, m * M / r2 * c], [m * M / r2 * c, m * M * M * H]])
    a1 = np.array([[3 * k / m * c * c + K * j / m, 3 * (k / m - j / M) * s * c], [3 * k / m * s * c, 3 * k / m * s * s]])
    return {"A1": (F1, m1, a1), "E": (F3, m3, a3)}


def nh3_secular(k, j, m=NH3["m"], M=NH3["M"], D=NH3["D"]) -> dict:
    """(3.7.25)、(3.7.29)：S₃ = (3k/2 + Hj)/m、P₃ = (3kj/2m²)(H − cos²D/2)、S₁ = (3k + Kj)/m、P₁ = 3kj(1/m² + 3/mM) sin²D；
    根 λ = ω² = (S ± √(S² − 4P))/2。回傳 dict(S1, P1, S3, P3, w1p, w1m, w3p, w3m)（ω，不是 ω²）。"""
    d = math.radians(float(D))
    c2, s2 = math.cos(d) ** 2, math.sin(d) ** 2
    H, K = 1 + 1.5 * m / M * c2, 1 + 3 * m / M * s2
    S3, P3 = (1.5 * k + H * j) / m, 1.5 * k * j / m ** 2 * (H - c2 / 2)
    S1, P1 = (3 * k + K * j) / m, 3 * k * j * (1 / m ** 2 + 3 / (m * M)) * s2
    rt = lambda S, P: (math.sqrt((S + math.sqrt(S * S - 4 * P)) / 2), math.sqrt(max((S - math.sqrt(S * S - 4 * P)) / 2, 0)))  # noqa: E731
    w1p, w1m = rt(S1, P1)
    w3p, w3m = rt(S3, P3)
    return dict(S1=S1, P1=P1, S3=S3, P3=P3, w1p=w1p, w1m=w1m, w3p=w3p, w3m=w3m)


def nh3_numeric(k, j, m=NH3["m"], M=NH3["M"], D=NH3["D"]) -> np.ndarray:
    """12×12 的 ω²（由小到大；前 6 個是平移與轉動的 0）。"""
    return np.sort(np.linalg.eigvals(np.linalg.solve(nh3_mass(m, M), nh3_force(k, j, D))).real)


def nh3_omega_obs(nu=None) -> dict:
    """ω = 2πcν（10¹⁴ rad/s 為單位前的 rad/s）。"""
    nu = nu or NH3_NU
    return {key: 2 * math.pi * C_LIGHT * v for key, v in nu.items()}


def nh3_fit(m=NH3["m"], M=NH3["M"], D_assumed=NH3["D"], m_kg=M_H_BOOK, book_k: bool = False) -> dict:
    """(3.7.30)–(3.7.31) 與 p. 216：R = (ω₁₊ω₁₋/ω₃₊ω₃₋)²、D = arcsin √[R(1+3m/M)/(4(1+3m/M) − R(1−3m/M))]；
    j/m = (2S₃ − S₁)/(1 + 3(m/M)cos 2D)、3k/2m = (S₁H − S₃K)/(1 + 3(m/M)cos 2D)（書上第二式漏了分母：book_k=True）。"""
    om = nh3_omega_obs()
    R = (om["p1"] * om["m1"] / (om["p3"] * om["m3"])) ** 2
    q = 3 * m / M
    s2 = R * (1 + q) / (4 * (1 + q) - R * (1 - q))
    Dfit = math.degrees(math.asin(math.sqrt(s2)))
    d = math.radians(float(D_assumed))
    S1, S3 = om["p1"] ** 2 + om["m1"] ** 2, om["p3"] ** 2 + om["m3"] ** 2
    H, K = 1 + 1.5 * m / M * math.cos(d) ** 2, 1 + 3 * m / M * math.sin(d) ** 2
    den = 1 + q * math.cos(2 * d)
    jm = (2 * S3 - S1) / den
    km = (S1 * H - S3 * K) * 2 / 3 / (1 if book_k else den)
    return dict(R=R, D=Dfit, j=jm * m_kg, k=km * m_kg, kj=km / jm)


def nh3_scan(kj, D=NH3["D"], m=NH3["m"], M=NH3["M"]) -> tuple:
    """Fig. 3.7.4：ω（以 √(j/m) 為單位）對 k/j：(ω₁₊, ω₁₋, ω₃₊, ω₃₋)。"""
    r = nh3_secular(kj, 1.0, m, M, D)
    return r["w1p"], r["w1m"], r["w3p"], r["w3m"]


def nh3_ratio_fit(D=NH3["D"]) -> dict:
    """在 k/j 上找 ω₃₋/ω₁₋ = 實測比（書：「k/j = 0.3 附近」）；用兩個低頻定尺度，預測 ω₁₊、ω₃₊（10¹⁴ rad/s）。"""
    om = nh3_omega_obs()
    obs = om["m3"] / om["m1"]
    x0 = brentq(lambda x: nh3_scan(x, D)[3] / nh3_scan(x, D)[1] - obs, 0.05, 0.6)
    w = nh3_scan(x0, D)
    unit = om["m1"] / w[1]
    return dict(kj=x0, ratio=obs, unit=unit / 1e14, w1p=w[0] * unit / 1e14, w3p=w[2] * unit / 1e14,
                w1m=w[1] * unit / 1e14, w3m=w[3] * unit / 1e14)


def nh3_crossing(D=NH3["D"]) -> float:
    """ω₁₊ = ω₃₊ 的偶然簡併（書：k/j = 0.035）。"""
    return brentq(lambda x: nh3_scan(x, D)[0] - nh3_scan(x, D)[2], 0.005, 0.1)


def nh3_modes(kj=NH3_KJ, D=NH3["D"], m=NH3["m"], M=NH3["M"]) -> list:
    """四種真振動（給動畫）：[(名稱, ω（√(j/m) 單位）, 12 維位移（最大分量 1）)]，名稱 "1+"、"1-"、"3+"、"3-"（E 取型 1 夥伴）。"""
    F, Mm = nh3_force(kj, 1.0, D), nh3_mass(m, M)
    cv = nh3_cvectors(m, M, D)
    out = []
    for key, (u, v) in {"1": (cv["c1"], cv["c1'"]), "3": (cv["c3"], cv["c3'"])}.items():
        B = np.vstack([u, v]).T
        w, X = np.linalg.eig(np.linalg.solve(B.T @ Mm @ B, B.T @ F @ B))
        for sg, idx in (("+", int(np.argmax(w.real))), ("-", int(np.argmin(w.real)))):
            vec = B @ X[:, idx].real
            vec = vec / np.abs(vec).max()
            out.append((key + sg, math.sqrt(w[idx].real), vec))
    return out


# ======================================================================== 習題
def problem_311() -> dict:
    """3.1.1：C_v = {1, σ₃} 的右／左陪集、C₃ 的右／左陪集；{1,σ₁}、{1,σ₂} 的雙陪集；{1,σ₁}、C₃ 的雙陪集。"""
    Cv, C3 = ["1", "s3"], ["1", "r", "r^2"]
    return dict(a=C3V.cosets(Cv, "right"), b=C3V.cosets(Cv, "left"), c=C3V.cosets(C3, "right"), d=C3V.cosets(C3, "left"),
                e=C3V.double_cosets(["1", "s1"], ["1", "s2"]), f=C3V.double_cosets(["1", "s1"], C3))


def problem_312() -> dict:
    """3.1.2：{1, σ₃} 不是正規子群；C₃ 是（左陪集 ＝ 右陪集）。"""
    return dict(a=C3V.is_normal(["1", "s3"]), b=C3V.is_normal(["1", "r", "r^2"]))


def problem_322() -> dict:
    """3.2.2：°N(g)·°C_g = °G（D₂、D₃、D₄ 每個元素）。回傳 {群: [(g, °N, °C)]}。"""
    return {G.label: [(g, len(G.normalizer(g)), len(G.classes[G.class_of(g)])) for g in G.names] for G in (D2, D3, D4)}


def problem_331() -> dict:
    """3.3.1：D₄ 的類代數 (3.6.1b)、冪等元 (3.6.2)、譜分解係數 c_j^{(α)}、用 D₂ 與 C₄ 分裂 ℙ^E。"""
    ells, T, perm = book_order(D4, D4_CHARS_BOOK)
    ids = [D4.class_idempotents[a] for a in perm]
    sz = D4.class_sizes()
    spec = sp.Matrix(5, 5, lambda a, j: sz[j] * T[a, j] / ells[a])
    PE = D4.central_idempotents[perm[4]]
    pD2 = []
    for sx in (1, -1):
        for sy in (1, -1):
            pD2.append(D4.ga({"1": R_(1, 4), "R1^2": R_(sx, 4), "R2^2": R_(sy, 4), "R^2": R_(sx * sy, 4)}))
    return dict(rows=D4.class_alg_rows(), idemp=ids, spectral=spec, split_D2=D4.split(PE, pD2),
                split_C4=D4.split(PE, D4.cyclic_projectors("R")), minpolys=[D4.min_poly_class(2), D4.min_poly_class(3)])


def problem_332() -> dict:
    """3.3.2、3.4.2：四元數群：類、類代數、特徵標表（＝ D₄ 的），E 用 {1, k, −1, −k} 分裂（複數）；Q ≇ D₄（2 階元素 1 個 vs 5 個）。"""
    ells, T = QUAT.char_table_class_algebra()
    irr = QUAT.irreps("k")
    order2 = {G.label: sum(1 for g in G.names if G.order_of(g) == 2) for G in (QUAT, D4)}
    return dict(classes=QUAT.classes, ells=ells, chars=T, irreps=irr, order2=order2,
                split=QUAT.split(QUAT.central_idempotents[-1], QUAT.cyclic_projectors("k")))


P341 = dict(g=dict(e=0, s=5, d=3, m=1), h=dict(e=1, s=5, d=3, m=1))


def d4_square_matrix(e, s, d, m=1) -> sp.Matrix:
    """習題 3.4.1：正方形四角的質點（右上 = e(45°)），邊彈簧 s、對角彈簧 d、外接彈簧 e（沿對角線向外）。
    8 個座標以 D₄ 元素命名：|g⟩ = g|1⟩，|1⟩ = 右上質點的 x̂（|i₃⟩ = ŷ）。精確 8×8 加速度矩陣。"""
    u = {g: D4_2D[g] * sp.Matrix([1, 0]) for g in D4_NAMES}
    p = {g: D4_2D[g] * sp.Matrix([1, 1]) for g in D4_NAMES}            # 質點位置（邊長 2）
    masses = []
    for g in D4_NAMES:
        if not any(p[g] == q for q in masses):
            masses.append(p[g])
    which = {g: next(i for i, q in enumerate(masses) if q == p[g]) for g in D4_NAMES}
    springs = []
    for a in range(4):
        for b in range(a + 1, 4):
            dv = masses[b] - masses[a]
            L2 = dv.dot(dv)
            springs.append((a, b, s if L2 == 4 else d, dv / sp.sqrt(L2)))
    K = sp.zeros(8)
    for i, g in enumerate(D4_NAMES):
        for j, h in enumerate(D4_NAMES):
            a, b = which[g], which[h]
            for (x, y, kk, n) in springs:
                if {a, b} == {x, y}:
                    K[i, j] -= kk * u[g].dot(n) * u[h].dot(n)
                elif a == b and a in (x, y):
                    K[i, j] += kk * u[g].dot(n) * u[h].dot(n)
            if a == b:
                n = masses[a] / sp.sqrt(2)
                K[i, j] += e * u[g].dot(n) * u[h].dot(n)
    return (K / m).applyfunc(sp.nsimplify)


def problem_341(e, s, d, m=1) -> dict:
    """3.4.1：塊對角約化（(3.4.22)：a^α_{jl} = Σ_g D^α_{jl}(g)* ⟨1|a|g⟩，E 用駐波 (3.6.7)）與 8×8 數值對角化。
    回傳 dict(blocks={A1, B1, A2, B2: 數, E: 2×2}, omega2=排序的 ω²（精確，含重數）, numeric=數值)。"""
    A = d4_square_matrix(e, s, d, m)
    row = {g: A[0, j] for j, g in enumerate(D4_NAMES)}
    ells, T, perm = book_order(D4, D4_CHARS_BOOK)
    blocks = {}
    for a, nm in enumerate(D4_IRREP_NAMES[:4]):
        blocks[nm] = sp.nsimplify(sum(T[a, D4.class_of(g)] * row[g] for g in D4_NAMES))
    blocks["E"] = sum((D4_2D[g] * row[g] for g in D4_NAMES), sp.zeros(2)).applyfunc(sp.nsimplify)
    w2 = [blocks[k] for k in D4_IRREP_NAMES[:4]]
    for ev, mult in blocks["E"].eigenvals().items():
        w2 += [sp.nsimplify(ev)] * (2 * mult)
    w2 = sorted(w2, key=float)
    num = np.sort(np.linalg.eigvalsh(np.array(A.evalf(), dtype=float)))
    return dict(blocks=blocks, omega2=w2, numeric=num)


def problem_351() -> dict:
    """3.5.1：C3v 類代數的 3×3 正則表示 R(c_j)；本徵值 c_j^α ⇒ χ^α_j = ℓ^α c_j^α / °c_j。"""
    mats = [C3V.class_matrix(i) for i in range(3)]
    spec = c3v_class_spectral()
    ells, T = C3V.char_table_class_algebra()
    sz = C3V.class_sizes()
    chi = sp.Matrix(3, 3, lambda a, j: ells[a] * spec[a, j] / sz[j])
    return dict(mats=mats, eig=[sorted(M.eigenvals().keys()) for M in mats], chi=chi, ok=chi == T)


def problem_352() -> dict:
    """3.5.2：Σ_g D^α_{ij}(g)* D^β_{kl}(g) = (°G/ℓ^α) δδδ；Σ_α χ^α_j* χ^α_k = (°G/°c_j) δ_{jk}；Σ_j °c_j χ^α_j* χ^β_j = °G δ。"""
    irr = c3v_irreps("cv")
    ok1 = True
    for a, Da in enumerate(irr):
        for b, Db in enumerate(irr):
            la, lb = Da["1"].shape[0], Db["1"].shape[0]
            for i in range(la):
                for j in range(la):
                    for k in range(lb):
                        for l_ in range(lb):
                            v = sp.simplify(sum(sp.conjugate(Da[g][i, j]) * Db[g][k, l_] for g in C3V.names))
                            want = R_(6, la) if (a == b and i == k and j == l_) else 0
                            ok1 &= sp.simplify(v - want) == 0
    ells, T = C3V.char_table_class_algebra()
    sz = C3V.class_sizes()
    comp = sp.Matrix(3, 3, lambda j, k: sum(sp.conjugate(T[a, j]) * T[a, k] for a in range(3)))
    orth = sp.Matrix(3, 3, lambda a, b: sum(sz[j] * sp.conjugate(T[a, j]) * T[b, j] for j in range(3)))
    return dict(great=ok1, completeness=comp, orthogonality=orth)


def problem_355() -> dict:
    """3.5.5：忠實表示的最低維數（群：最小的忠實不可約表示；群代數：Σ ℓ^α）、互相對易的群元素最多幾個、群代數的秩 Σ ℓ^α。"""
    out = {}
    for G in (C3V, C4V):
        ells, T = G.char_table_class_algebra()
        ab = max(len(H) for H in _abelian_subgroups(G))
        out[G.label] = dict(group=min(int(l_) for l_ in ells if l_ > 1), algebra=int(sum(ells)), commuting=ab, rank=int(sum(ells)))
    return out


def _abelian_subgroups(G: FiniteGroup) -> list:
    """所有由一或兩個元素生成的 Abel 子群（小群夠用）。"""
    out = []
    for a in G.names:
        for b in G.names:
            if G.table[a, b] != G.table[b, a]:
                continue
            H = {G.e, a, b}
            while True:
                H2 = H | {G.table[x, y] for x in H for y in H}
                if H2 == H:
                    break
                H = H2
            if all(G.table[x, y] == G.table[y, x] for x in H for y in H):
                out.append(sorted(H, key=G.idx.get))
    return out


def problem_356() -> dict:
    """3.5.6：°G 被子群的階、類的大小、不可約表示的維數整除（C3v、D₄、D₆、Q 逐一核對）。"""
    out = {}
    for G in (C3V, D4, D6, QUAT):
        ells, T = G.char_table_class_algebra()
        out[G.label] = (all(G.order % len(C) == 0 for C in G.classes), all(G.order % int(l_) == 0 for l_ in ells),
                        all(G.order % len(H) == 0 for H in _abelian_subgroups(G)))
    return out


# ======================================================================== 摘要
def summary() -> dict:
    """{名稱: 字串}：場景卡與導讀上的主要數字。"""
    S = {}
    S["C3v 類"] = str(C3V.classes)
    S["C3v 類代數"] = str(C3V.class_alg_rows())
    S["c3 最小方程"] = str(C3V.min_poly_class(2).as_expr())
    S["C3v 特徵標"] = str(C3V.char_table_class_algebra())
    w = six_dof_omega2()
    S["六自由度 SIX ω²"] = str(w)
    S["六自由度 SIX ω"] = str([_f4(sp.sqrt(w["A1"])), _f4(sp.sqrt(w["A2"])), _f4(sp.sqrt(w["E"][0])), _f4(sp.sqrt(w["E"][1]))])
    S["六自由度 第一列"] = str([_f4(x) for x in six_dof_first_row()])
    w0 = six_dof_omega2(SIX["k1"], 0)
    S["Fig 3.4.4 (k2=0) ω"] = str([_f4(sp.sqrt(w0["A1"])), _f4(sp.sqrt(w0["E"][0]))])
    for wh in (0, 1):
        S[f"圓周 k2=0 which={wh}"] = str(six_dof_circular(SIX["k1"], 0, which=wh))
        S[f"圓周 SIX which={wh}"] = str(six_dof_circular(which=wh))
    for n in (4, 6):
        S[f"D{n} weak"] = str([(f"{e:.4f}", l_) for e, l_ in dn_levels(n, ROTOR["weak"], 2 * n)])
        S[f"D{n} strong"] = str([(f"{e:.4f}", l_) for e, l_ in dn_levels(n, ROTOR[f"strong{n}"], 2 * n)])
        S[f"D{n} 2x2"] = str(dn_gap_2x2(n, ROTOR["weak"]))
    S["Zeeman D4"] = str([(f"{e:.4f}", k) for e, k in dn_zeeman(4, ROTOR["weak"], ROTOR["zeeman_b"], 6)])
    S["NH3 ω obs"] = str({k: f"{v / 1e14:.4f}" for k, v in nh3_omega_obs().items()})
    S["NH3 fit"] = str(nh3_fit())
    S["NH3 fit book k"] = str(nh3_fit(book_k=True))
    S["NH3 ratio fit"] = str(nh3_ratio_fit())
    S["NH3 crossing"] = f"{nh3_crossing():.4f}"
    S["NH3 secular k/j=0.3"] = str({k: f"{v:.4f}" for k, v in nh3_secular(NH3_KJ, 1.0).items()})
    S["NH3 scan"] = str({x: [f"{v:.4f}" for v in nh3_scan(x)] for x in (0, 0.035, 0.1, 0.2, 0.3)})
    for key, par in P341.items():
        r = problem_341(**par)
        S[f"習題 3.4.1 {key}"] = str(r["blocks"]) + " ω = " + str([_f4(sp.sqrt(x)) for x in r["omega2"]])
    return S


# ======================================================================== 自我測試
def summary_numbers() -> dict:
    """check() 用：轉子能階的 4 位小數字串。"""
    out = {}
    for n in (4, 6):
        out[f"D{n} weak"] = [f"{e:.4f}" for e, _ in dn_levels(n, ROTOR["weak"], 2 * n)]
        out[f"D{n} strong"] = [f"{e:.4f}" for e, _ in dn_levels(n, ROTOR[f"strong{n}"], 2 * n)]
    return out


def _check_groups():
    G = C3V
    # ---------- §3.1 群、乘法表、同構、Fig. 3.1.2
    assert G.table_rows() == C3V_TABLE_BOOK                                # 矩陣相乘 ＝ (3.1.1) 照抄
    assert G.mul("s1", "s2") == "r" and G.mul("s2", "s1") == "r^2" and not G.is_abelian()
    assert refl2(120) * refl2(60) == rot2(120) and refl2(60) * refl2(120) == rot2(240)
    inv_map = {v: k for k, v in D3_TO_C3V.items()}
    assert D3.iso_to(C3V, D3_TO_C3V) and D3.table_rows() == [[inv_map[x] for x in r] for r in C3V_TABLE_BOOK]
    for p, s in (("p1", "s1"), ("p2", "s2"), ("p3", "s3")):                 # (3.1.2) σ_j = Iρ_j
        assert INV3 * D3.mats[p] == block3(C3V.mats[s])
    f = fig312()
    for key in f:
        mid, end, direct = f[key]
        assert all(_close(end[b], direct[b]) for b in end)
    assert _close(f["s2s1"][1]["a"], np.array(rot2(240).evalf(), dtype=float) @ np.array(BOB_START["a"]))
    phi = sp.Symbol("phi")
    Sf = sp.Matrix([[sp.cos(2 * phi), sp.sin(2 * phi)], [sp.sin(2 * phi), -sp.cos(2 * phi)]])
    Rf = sp.Matrix([[sp.cos(2 * phi), -sp.sin(2 * phi)], [sp.sin(2 * phi), sp.cos(2 * phi)]])
    assert sp.simplify(Sf * sp.diag(1, -1) - Rf) == sp.zeros(2)            # Fig. 3.1.4：σ_{φz}σ_{xz} = R[2φ]
    # ---------- Hamilton 轉：弧 ＝ 四元數 ＝ 3×3 ＝ SU(2) ＝ 兩次反射
    rng = np.random.default_rng(7)
    for _ in range(6):
        a1, a2 = rng.normal(size=3), rng.normal(size=3)
        d1, d2 = rng.uniform(10, 350, 2)
        hp = hamilton_product(a1, d1, a2, d2)
        assert _close(hp["q_arc"], hp["q"])                                 # 頭尾相接的弧 ＝ 四元數乘積（含正負號）
        R2R1 = quat_to_rot(quat(a2, d2)) @ quat_to_rot(quat(a1, d1))
        assert _close(quat_to_rot(hp["q"]), R2R1)
        assert _close(reflection3(hp["N2p"]) @ reflection3(hp["N1"]), R2R1)  # 兩次反射
        assert _close(su2(quat(a2, d2)) @ su2(quat(a1, d1)), su2(hp["q"]))
        half = math.degrees(math.acos(float(np.clip(hp["N1"] @ hp["N2p"], -1, 1))))
        assert abs(2 * half - hp["deg"]) < 1e-7 or abs(360 - 2 * half - hp["deg"]) < 1e-7
    ex = hamilton_product((0, 0, 1), 90, (1, 0, 0), 90)                     # H302：先 R_z(90°) 再 R_x(90°)
    assert abs(ex["deg"] - 120) < 1e-9 and _close(ex["axis"], np.array([1, -1, 1]) / math.sqrt(3))
    ex2 = hamilton_product((1, 0, 0), 90, (0, 0, 1), 90)
    assert abs(ex2["deg"] - 120) < 1e-9 and _close(ex2["axis"], np.array([1, 1, 1]) / math.sqrt(3))
    arcs = d3_arcs()
    assert abs(math.degrees(math.acos(arcs["r"][0] @ arcs["r"][1])) - 60) < 1e-9
    assert all(abs(math.degrees(math.acos(float(np.clip(arcs[p][0] @ arcs[p][1], -1, 1)))) - 90) < 1e-9 for p in D3_AXIS_DEG)
    ax = {p: np.array([math.cos(math.radians(a)), math.sin(math.radians(a)), 0]) for p, a in D3_AXIS_DEG.items()}
    hp = hamilton_product((0, 0, 1), 120, ax["p1"], 180)                    # ρ₁ r = ρ₂（Fig. 3.1.7）
    assert abs(hp["deg"] - 180) < 1e-9 and abs(abs(hp["axis"] @ ax["p2"]) - 1) < 1e-9 and D3.mul("p1", "r") == "p2"
    assert _close(qmul(quat(ax["p1"], 180), quat((0, 0, 1), 120)), quat(ax["p2"], 180))        # 這組軸向：+ρ₂
    assert _close(qmul(quat(-ax["p1"], 180), quat((0, 0, 1), 120)), -quat(ax["p2"], 180))      # ρ₁ 軸反向：(−)ρ₂
    assert body_product(G, "s1", "s2") == "r" == G.mul("s1", "s2") and G.mul("s2", "s1") == "r^2"   # Fig. 3.1.8
    axc, angc, Tn = conjugate_rotation(ax["p2"], 180, ax["p1"], 180)       # ρ₂ρ₁ρ₂⁻¹ = ρ₃：同角、軸 = T·n
    assert abs(angc - 180) < 1e-9 and abs(abs(axc @ Tn) - 1) < 1e-9 and abs(abs(Tn @ ax["p3"]) - 1) < 1e-9
    assert D3.conj("p2", "p1") == "p3" and G.conj("s2", "s1") == "s3" and G.conj("r", "s3") == "s1"
    # ---------- §3.2 類、Lagrange、類代數、冪等元
    assert G.classes == [["1"], ["r", "r^2"], ["s1", "s2", "s3"]]
    assert G.normalizer("r") == ["1", "r", "r^2"] and G.normalizer("s1") == ["1", "s1"]
    assert [G.order // len(G.normalizer(C[0])) for C in G.classes] == G.class_sizes() == [1, 2, 3]
    assert G.cosets(["1", "r", "r^2"]) == [["1", "r", "r^2"], ["s1", "s2", "s3"]]
    assert G.cosets(["1", "s1"]) == [["1", "s1"], ["r", "s3"], ["r^2", "s2"]]
    c1, c2, c3 = sp.symbols("c1 c2 c3")
    book = [[sp.sympify(x, locals={"c1": c1, "c2": c2, "c3": c3}) for x in row] for row in C3V_CLASS_ALG_BOOK]
    assert G.class_alg_rows([c1, c2, c3]) == book                            # (3.2.1b)
    for h in G.names:                                                      # (3.2.2) h⁻¹ c h = c
        for i in range(3):
            assert G.ga_eq(G.ga_mul(G.ga_mul(G.el(G.inv(h)), G.class_op(i)), G.el(h)), G.class_op(i))
    cc3 = sp.Matrix([0, 0, 1])
    sq = G.class_mul(cc3, cc3)
    assert list(sq) == [3, 3, 0] and list(G.class_mul(sq, cc3)) == [0, 0, 9]     # c₃² = 3c₁ + 3c₂、c₃³ = 9c₃
    x = sp.Symbol("c")
    assert G.min_poly_class(2).as_expr() == x ** 3 - 9 * x
    C3m = G.class_matrix(2)
    assert C3m ** 3 - 9 * C3m == sp.zeros(3) and C3m ** 2 - 9 * sp.eye(3) != sp.zeros(3) and C3m ** 2 != 9 * C3m
    ids = G.class_idempotents
    assert [list(v) for v in ids] == C3V_CLASS_IDEMP_BOOK                  # (3.2.15) ＝ (3.5.4)ₓ
    Pl = G.central_idempotents
    for a in range(3):
        for b in range(3):
            assert G.ga_eq(G.ga_mul(Pl[a], Pl[b]), Pl[a] if a == b else {})
        assert G.ga_eq(G.ga_dag(Pl[a]), Pl[a])                             # (3.2.22)
    assert G.ga_eq(G.ga_add(*Pl), G.el("1"))
    spec = c3v_class_spectral()
    assert spec == sp.Matrix([[1, 2, 3], [1, 2, -3], [1, -1, 0]])          # (3.2.17)：c₂ = 2ℙ¹+2ℙ²−ℙ³、c₃ = 3ℙ¹−3ℙ²
    for j in range(3):
        assert G.ga_eq(G.ga_add(*Pl, coeffs=[spec[a, j] for a in range(3)]), G.class_op(j))
    # ---------- §3.3 分裂、基本算子、不可約表示
    P31, P32 = c3v_split_cv()
    assert ga_tex(G, Pl[2]) == r"\tfrac{1}{3}(2\mathbf 1-r-r^2)" and ga_tex(G, P31) == r"\tfrac{1}{6}(2\mathbf 1-r-r^2-\sigma_1-\sigma_2+2\sigma_3)"
    assert G.ga_eq(P31, G.ga_scale({"1": 2, "r": -1, "r^2": -1, "s1": -1, "s2": -1, "s3": 2}, R_(1, 6)))   # (3.3.3a)
    assert G.ga_eq(P32, G.ga_scale({"1": 2, "r": -1, "r^2": -1, "s1": 1, "s2": 1, "s3": -2}, R_(1, 6)))    # (3.3.3b)
    pm = G.cyclic_projectors("s3")
    assert not G.ga_mul(pm[1], Pl[0]) and not G.ga_mul(pm[0], Pl[1])         # (3.3.2c) 的兩個 0
    Pij = c3v_elementary("cv")
    assert G.ga_eq(G.ga_mul(Pij[0, 1], Pij[1, 0]), Pij[0, 0]) and not G.ga_mul(Pij[0, 1], Pij[0, 1])
    assert G.ga_eq(G.ga_mul(G.ga_mul(P31, G.el("s2")), P32),
                   G.ga_scale({"r": -1, "r^2": 1, "s1": -1, "s2": 1}, R_(1, 4)))   # p. 173：P³₁σ₂P³₂
    irr = c3v_irreps("cv")
    assert all(irr[2][g] == G.mats[g] for g in G.names)                    # 基本算子的 𝒟³ ＝ 幾何的 2×2 ＝ (3.3.7)
    assert irr[2]["r"] == sp.Matrix([[R_(-1, 2), -SQ3 / 2], [SQ3 / 2, R_(-1, 2)]])
    assert all(G.is_rep(D) for D in irr)
    for g in G.names:                                                      # (3.3.16)：g = Σ 𝒟^α_{ij}(g) P^α_{ij}
        terms = [G.ga_scale(Pl[0], irr[0][g][0, 0]), G.ga_scale(Pl[1], irr[1][g][0, 0])]
        terms += [G.ga_scale(Pij[i, j], irr[2][g][i, j]) for i in range(2) for j in range(2)]
        assert G.ga_eq(G.ga_add(*terms), G.el(g))
    for i in range(2):                                                     # (3.4.19)：P_{jk} = (ℓ/°G) Σ 𝒟*_{jk}(g) g
        for j in range(2):
            assert G.ga_eq(Pij[i, j], G.ga({g: R_(2, 6) * sp.conjugate(irr[2][g][i, j]) for g in G.names}))
    assert sum(int(D["1"].shape[0]) ** 2 for D in irr) == 6                 # (3.3.17)
    Pc = G.split(Pl[2], G.cyclic_projectors("r"))                           # (3.3.18)
    eps = EPS3
    assert len(Pc) == 2 and not G.ga_mul(Pl[2], G.cyclic_projectors("r")[0])
    assert G.ga_eq(Pc[0], G.ga({"1": R_(1, 3), "r": sp.conjugate(eps) / 3, "r^2": eps / 3}))
    assert G.ga_eq(Pc[1], G.ga({"1": R_(1, 3), "r": eps / 3, "r^2": sp.conjugate(eps) / 3}))
    assert not G.ga_eq(Pc[1], G.ga({"1": R_(1, 3), "r": eps / 3 + sp.conjugate(eps) / 3}))   # 書上照印的 ε*r（應為 ε*r²）
    irc = c3v_irreps("c3")[2]
    for g in G.names:                                                      # (3.4.24)：T†𝒟³T ＝ 𝒟^{3°}（ε = e^{+2πi/3}）
        assert sp.simplify(sp.expand_complex(T_CIRC.H * G.mats[g] * T_CIRC - irc[g])) == sp.zeros(2)
    assert irc["r"] == sp.diag(eps, sp.conjugate(eps)) and irc["s1"] == sp.Matrix([[0, sp.conjugate(eps)], [eps, 0]])
    assert irc["s2"] == sp.Matrix([[0, eps], [sp.conjugate(eps), 0]]) and irc["s3"] == sp.Matrix([[0, 1], [1, 0]])
    assert sp.simplify(sp.expand_complex((T_CIRC.H * G.mats["r"] * T_CIRC)[0, 0] - sp.exp(-2 * sp.pi * I / 3))) != 0   # 書上的 e^{−2πi/3} 不對
    # ---------- §3.4 正則表示與約化
    R_r = G.regular("r")
    assert R_r == sp.Matrix([[0, 0, 1, 0, 0, 0], [1, 0, 0, 0, 0, 0], [0, 1, 0, 0, 0, 0],
                             [0, 0, 0, 0, 1, 0], [0, 0, 0, 0, 0, 1], [0, 0, 0, 1, 0, 0]])     # (3.4.4)
    assert [G.regular(g).trace() for g in G.names] == [6, 0, 0, 0, 0, 0]                    # (3.4.5)
    assert all(G.regular(G.mul(a, b)) == G.regular(a) * G.regular(b) for a in G.names for b in G.names)
    Tm = c3v_regular_T()
    keys = ["P1", "P2", "P11", "P21", "P12", "P22"]
    assert all(sp.simplify(Tm[:, k] - GROUP_VECTORS_BOOK[key]) == sp.zeros(6, 1) for k, key in enumerate(keys))   # (3.4.8)
    assert sp.simplify(Tm.T * Tm) == sp.eye(6)                              # (3.4.10a)
    for g in G.names:                                                      # (3.4.16)
        assert sp.simplify(Tm.T * G.regular(g) * Tm) == sp.diag(irr[0][g], irr[1][g], irr[2][g], irr[2][g])
    assert sp.simplify(R_r * Tm[:, 2] - (R_(-1, 2) * Tm[:, 2] + SQ3 / 2 * Tm[:, 3])) == sp.zeros(6, 1)     # (3.4.12c)
    assert [Pl[a].get("1") for a in range(3)] == [R_(1, 6), R_(1, 6), R_(2, 3)]    # ⟨1|ℙ^α|1⟩ = (ℓ^α)²/°G
    return irr, Pl, ids, R_r


def _check_six_and_chars(irr, Pl, ids, R_r):
    G = C3V
    eps = EPS3
    A6 = six_dof_matrix()
    assert _close(np.array(A6.evalf(), dtype=float), six_dof_geometric())   # 第一列公式 ＝ 彈簧幾何（整個 6×6）
    assert _close(np.array(six_dof_matrix(1, 0.37, 30).evalf(), dtype=float), six_dof_geometric(1, 0.37, 30))
    for g in G.names:
        v = np.array(G.mats[g].evalf(), dtype=float) @ e2(-45)
        assert _close(v, e2(SIX_ANGLES[g])) and _close(np.array(G.mats[g].evalf(), dtype=float) @ e2(0), e2(120 * SIX_MASS[g]))
    k1s, k2s, ms_, al = sp.symbols("k1 k2 m alpha", positive=True)
    gen, bk = six_dof_blocks(k1s, k2s, al, ms_), six_dof_blocks_book(k1s, k2s, al, ms_)
    assert sp.simplify(gen[0] - bk[0]) == 0 and sp.simplify(gen[1] - bk[1]) == 0 and sp.simplify(gen[2] - bk[2]) == sp.zeros(2)  # (3.4.23b)
    assert six_dof_omega2() == {"A1": 24, "A2": 4, "E": (16, 6)}
    assert _close(np.sort(np.linalg.eigvalsh(six_dof_geometric())), [4, 6, 6, 16, 16, 24])
    assert [_f4(sp.sqrt(v)) for v in (24, 4, 16, 6)] == ["4.8990", "2.0000", "4.0000", "2.4495"]
    assert [_f4(x) for x in six_dof_first_row()] == ["12.0000", "1.0000", "1.0000", "0.2679", "3.7321", "6.0000"]
    assert six_dof_omega2(SIX["k1"], 0) == {"A1": 12, "A2": 0, "E": (6, 0)} and _f4(sp.sqrt(12)) == "3.4641"
    A0 = np.array(six_dof_matrix(SIX["k1"], 0).evalf(), dtype=float)
    assert _close(np.sort(np.linalg.eigvalsh(A0)), [0, 0, 0, 6, 6, 12])    # 三個非真振動（平移 x、y 與轉動）
    fx, fb = fig344_vectors("fixed"), fig344_vectors("book")
    lam = {"P1": 12, "P2": 0, "x+": 6, "y+": 6, "x-": 0, "y-": 0}
    A0s = six_dof_matrix(SIX["k1"], 0)
    for key, v in fx.items():
        assert sp.simplify(A0s * v - lam[key] * v) == sp.zeros(6, 1) and sp.simplify(v.dot(v) - 1) == 0
    assert sp.simplify(A0s * fb["y+"] - 6 * fb["y+"]) != sp.zeros(6, 1)     # 書上型 2 的「+」欄不是本徵向量
    assert sp.simplify(fb["y-"].dot(fb["y-"]) - 1) != 0                     # 書上型 2 的「−」欄沒有歸一（(24 + 4√3)/24）
    x6, y6, _ = six_dof_mode(SIX["k1"], 0, which=0)
    assert _close(x6, np.array(fx["x+"].evalf(), dtype=float).ravel()) and _close(y6, np.array(fx["y+"].evalf(), dtype=float).ravel())
    Rr = np.array(R_r.evalf(), dtype=float)
    assert _close(Rr @ (x6 + 1j * y6), complex(sp.conjugate(eps)) * (x6 + 1j * y6))   # |x⟩ + i|y⟩ 屬於 2₃（r → ε*）
    circ = six_dof_circular(SIX["k1"], 0, which=0)
    assert _close(circ["rmin"], [1 / math.sqrt(3)] * 3) and _close(circ["rmax"], [1 / math.sqrt(3)] * 3) and circ["sense"] == [-1, -1, -1]
    circ2 = six_dof_circular(which=0)
    assert _close(circ2["rmax"] / circ2["rmin"], [3, 3, 3]) and circ2["sense"] == [-1, -1, -1]
    # ---------- §3.5 特徵標、頻率、跡公式
    ells, T = G.char_table_class_algebra()
    assert ells == [1, 1, 2] and T == sp.Matrix(C3V_CHARS_BOOK) == G.char_table_traces(irr)   # (3.5.8)：類代數 ＝ 跡
    assert [sp.sqrt(6 * v[0]) for v in ids] == [1, 1, 2]                    # (3.5.7)ₓ
    for a in range(3):                                                     # (3.5.3)
        assert G.ga_eq(Pl[a], G.ga({g: ells[a] * T[a, G.class_of(g)] / 6 for g in G.names}))
    assert [pendulum_rep(C[0]).trace() for C in G.classes] == [3, 0, 1] and pendulum_frequencies() == [1, 0, 1]
    asym, bsym = sp.symbols("a b")
    te = pendulum_trace_eigs(asym, bsym)
    assert sp.simplify(te[0] - bsym) == 0 and te[1] is None and sp.simplify(te[2] - 3 * asym - bsym) == 0
    assert set(SA.ring_matrix(3, asym, bsym).eigenvals().keys()) == {bsym, 3 * asym + bsym}
    assert pendulum_trace_eigs() == [1, None, 4]
    A3 = SA.ring_matrix(3, asym, bsym)
    tr = [sp.expand((A3 * pendulum_rep(C[0])).trace()) for C in G.classes]
    assert tr == [6 * asym + 3 * bsym, -3 * asym, bsym]                    # (3.5.14) 之後
    assert sp.simplify(sum(T[0, k] * tr[k] for k in range(3)) / 6 - bsym) != 0   # 照印 (3.5.13)（沒有類的大小）不對
    tang = [3, 0, -1]                                                      # 切向擺動：Tr 𝓡(σ) = −1 ⇒ A2 ⊕ E
    assert G.frequencies(tang) == [0, 1, 1]


def _check_d4_d6():
    x = sp.Symbol("c")
    assert D4.table_rows() == D4_TABLE_BOOK and D4.iso_to(C4V, D4_TO_C4V) and D2D.iso_to(D4, D2D_TO_D4)
    c = sp.symbols("c1:6")
    rows4 = D4.class_alg_rows(c)
    assert rows4[2][2] == 2 * c[0] + 2 * c[1] and rows4[2][3] == 2 * c[4] and rows4[2][4] == 2 * c[3] and \
        rows4[3][3] == 2 * c[0] + 2 * c[1] and rows4[3][4] == 2 * c[2] and rows4[4][4] == 2 * c[0] + 2 * c[1]   # (3.6.1b)
    assert D4.min_poly_class(2).as_expr() == x ** 3 - 4 * x == D4.min_poly_class(3).as_expr()
    assert len(D4.class_idempotents_of(2)) == 3 < 5                         # 一個方程不夠（p. 191）
    ells4, T4, perm4 = book_order(D4, D4_CHARS_BOOK)
    assert ells4 == [1, 1, 1, 1, 2] and T4 == sp.Matrix(D4_CHARS_BOOK)
    irr4 = D4.irreps("R1^2", "i3")
    assert D4.char_table_traces(irr4) == D4.char_table_class_algebra()[1]
    idem_book = [[R_(1, 8)] * 5, [R_(1, 8), R_(1, 8), -R_(1, 8), R_(1, 8), -R_(1, 8)],
                 [R_(1, 8), R_(1, 8), R_(1, 8), -R_(1, 8), -R_(1, 8)], [R_(1, 8), R_(1, 8), -R_(1, 8), -R_(1, 8), R_(1, 8)],
                 [R_(1, 2), -R_(1, 2), 0, 0, 0]]
    assert [list(D4.class_idempotents[a]) for a in perm4] == idem_book      # (3.6.2)
    DE = irr4[perm4[4]]
    assert all(DE[g] == D4_2D[g] for g in D4_NAMES)                          # (3.6.7) ＝ 平面上的作用
    PE = D4.central_idempotents[perm4[4]]
    sC2 = D4.split(PE, D4.cyclic_projectors("R1^2"))
    assert D4.ga_eq(sC2[0], D4.ga_scale({"1": 1, "R^2": -1, "R1^2": 1, "R2^2": -1}, R_(1, 4)))   # (3.6.6b)
    sC4 = D4.split(PE, D4.cyclic_projectors("R"))
    assert D4.ga_eq(sC4[0], D4.ga_scale({"1": 1, "R^2": -1, "R": -I, "R^3": I}, R_(1, 4)))       # 本系列的 P^E_{1₄}＝書上的 P^E_{3₄}
    assert D4.ga_eq(sC4[1], D4.ga_scale({"1": 1, "R^2": -1, "R": I, "R^3": -I}, R_(1, 4)))
    circE = D4.irrep_from_elementary(D4.elementary(sC4, "R1^2")[0])
    assert circE["R"] == sp.diag(I, -I) and circE["R1^2"] == sp.Matrix([[0, 1], [1, 0]]) and D4.is_rep(circE)
    assert D4.subgroup_of("R") == ["1", "R", "R^2", "R^3"] and D4.is_subgroup(["1", "R", "R^2", "R^3"])
    assert not D4.is_subgroup(["1", "R", "R^2"])                            # 書上 {1, R, R², R⁴} = {1, R, R²}：不是子群
    corr = D4.correlation(D4.subgroup_of("R"), cn_chars(4), T4)
    assert corr == sp.Matrix([[1, 0, 0, 0], [0, 0, 1, 0], [1, 0, 0, 0], [0, 0, 1, 0], [0, 1, 0, 1]])   # (3.6.14)–(3.6.15)
    D2sub = ["1", "R1^2", "R2^2", "R^2"]
    D2tab = sp.Matrix([[SA.D2_CHARS[n][h] for h in ("1", "Rx", "Ry", "Rz")] for n in ("A1", "B2", "B1", "A2")])
    assert D4.correlation(D2sub, D2tab, T4)[4, :] == sp.Matrix([[0, 1, 1, 0]])   # (3.6.8)：E ↓ D₂ = B₁ ⊕ B₂（第 2 章的 D₂ 標記）
    assert all(DE[h].is_diagonal() for h in D2sub)
    assert D4.correlation(["1", "R1^2"], C2_CHARS, T4) == sp.Matrix([[1, 0], [1, 0], [0, 1], [0, 1], [1, 1]])
    # ---------- §3.6 轉子：標記兩種算法、能隙、波腹位置、強位能的簇、Zeeman
    for n in (4, 6):
        lev = dn_levels(n, ROTOR["weak"], 2 * n + 1)                        # 2n + 1 個：到 m = ±n 的殼層剛好完整
        free = dn_free_levels(n, 2 * n + 1)
        shells = {}
        for (e, lab), (e0, lab0) in zip(lev, free):
            shells.setdefault(e0, [[], []])
            shells[e0][0].append(lab)
            shells[e0][1].append(lab0)
        assert all(sorted(a) == sorted(b) for a, b in shells.values())     # 數值本徵向量的標記 ＝ m mod n 的標記
        lo, hi = dn_gap_2x2(n, ROTOR["weak"])
        assert (lo, hi) == ((n // 2) ** 2 - 1, (n // 2) ** 2 + 1)
        eB1 = next(e for e, l_ in lev if l_ == "B_1")
        eB2 = next(e for e, l_ in lev if l_ == "B_2")
        assert abs((eB2 - eB1) - 2) < 3e-3 and abs(eB1 - lo) < 0.04 and abs(eB2 - hi) < 0.04
        dB1 = dn_density(n, ROTOR["weak"], "B_1", [0, math.pi / n])
        dB2 = dn_density(n, ROTOR["weak"], "B_2", [0, math.pi / n])
        assert dB1[0] > 1.9 and dB1[1] < 1e-12 and dB2[0] < 1e-12 and dB2[1] > 1.8   # B₁ 波腹在谷、B₂ 在山頂
        dA1 = dn_density(n, ROTOR["weak"], "A_1", [0, math.pi / n])
        assert dA1[0] > dA1[1]                                              # A₁ 也偏在谷
        Vs = np.linspace(0, ROTOR[f"strong{n}"], 13)
        lv, lb = dn_level_paths(n, Vs[1:], 2 * n)
        assert all(row == lb[0] for row in lb)                              # 弱 → 強，標記沿著每一條不變
    assert [l_ for _, l_ in dn_levels(4, ROTOR["strong4"], 8)] == ["A_1", "E", "E", "B_1", "B_2", "E", "E", "A_2"]
    assert [l_ for _, l_ in dn_levels(6, ROTOR["strong6"], 12)] == ["A_1", "E_1", "E_1", "E_2", "E_2", "B_1",
                                                                    "B_2", "E_2", "E_2", "E_1", "E_1", "A_2"]
    S_ = summary_numbers()
    assert S_["D4 weak"] == ["-0.1242", "0.8358", "0.8358", "2.9697", "4.9678", "9.0975", "9.0975", "15.9792"]
    assert S_["D4 strong"] == ["-6.0558", "-5.8194", "-5.8194", "-5.5627", "9.5168", "11.4095", "11.4095", "14.6889"]
    assert S_["D6 weak"][5:7] == ["7.9863", "9.9859"] and S_["D6 strong"][0] == "-13.6256" and S_["D6 strong"][5] == "-12.5161"
    assert S_["D6 strong"][6] == "21.4128" and S_["D6 strong"][11] == "33.0501"
    zl = dn_zeeman(4, ROTOR["weak"], ROTOR["zeeman_b"], 3)
    assert [k for _, k in zl] == [0, 1, 3] and [f"{e:.4f}" for e, _ in zl] == ["-0.1145", "0.6552", "1.0340"]
    z0 = dn_zeeman(4, ROTOR["weak"], 0, 3)
    assert abs(z0[1][0] - z0[2][0]) < 1e-12 and f"{z0[1][0]:.4f}" == "0.8358"   # b = 0：E 不分裂（＝ dn_levels）
    # ---------- §3.6 D₆ 與外積、D_nh、D_nd、相關表
    ells6, T6, perm6 = book_order(D6, D6_CHARS_BOOK)
    assert T6 == outer_table(C2_CHARS, D3_CHARS) == sp.Matrix(D6_CHARS_BOOK)   # 類代數 ＝ 外積 ＝ (3.6.17)
    dE = {"1": rot2(0), "h^2": rot2(120), "h^4": rot2(240), "p1": refl2(120), "p2": refl2(60), "p3": refl2(0)}
    E1, E2 = {}, {}
    for g in D6.names:                                                     # (3.6.18)、(3.6.20)：𝒟^E(of D₃)·𝒟^{B or A}(of C₂)
        base = {"h^3": "1", "h^5": "h^2", "h": "h^4"}.get(g, g.rstrip("'"))
        sgn = -1 if (g in ("h^3", "h^5", "h") or g.endswith("'")) else 1
        E1[g], E2[g] = sgn * dE[base], dE[base]
    assert D6.is_rep(E1) and D6.is_rep(E2) and E1["h"] == rot2(60) and E1["h^4"] == rot2(240) == E2["h"] == E2["h^4"]
    assert D6.char_table_traces([E1, E2]) == sp.Matrix([D6_CHARS_BOOK[5], D6_CHARS_BOOK[2]])
    w_ = sp.Symbol("w")
    assert sp.Matrix([[sp.cos(2 * w_), sp.sin(2 * w_)], [-sp.sin(2 * w_), sp.cos(2 * w_)]]).subs(w_, sp.pi / 3) == E2["h"]   # (3.6.21)
    assert all(rot2(60 * k) == E1[g] for g, k in (("h", 1), ("h^2", 2), ("h^3", 3), ("h^4", 4), ("h^5", 5)))   # (3.6.19)
    assert D6.correlation(D6.subgroup_of("h"), cn_chars(6), T6) == sp.Matrix(
        [[1, 0, 0, 0, 0, 0], [1, 0, 0, 0, 0, 0], [0, 0, 1, 0, 1, 0], [0, 0, 0, 1, 0, 0], [0, 0, 0, 1, 0, 0], [0, 1, 0, 0, 0, 1]])  # (3.6.27a)
    assert D6.correlation(["1", "p3"], C2_CHARS, T6) == sp.Matrix([[1, 0], [0, 1], [1, 1], [1, 0], [0, 1], [1, 1]])   # (3.6.27b)
    for Gp in (D3H, D3D):                                                  # (3.6.36)、(3.6.34)
        _, Tp, _ = book_order(Gp, D3H_CHARS_BOOK)
        assert Tp == outer_table(C2_CHARS, D3_CHARS) == sp.Matrix(D3H_CHARS_BOOK)
    ref2h = [list(r) for r in outer_table(C2_CHARS, sp.Matrix(D2_CHARS_BOOK)).tolist()]
    _, T2h, _ = book_order(D2H, ref2h)
    assert T2h.shape == (8, 8) and D2H.is_abelian()
    assert list(T2h.row(1)) == [1, -1, 1, -1, 1, -1, 1, -1]                # B_{1g}（(3.6.33)）＝ D₂ 表第 2 列（書上誤標 B₂）
    ell4h, T4h = D4H.char_table_class_algebra()
    assert sorted(ell4h) == [1] * 8 + [2, 2] and T4h.shape == (10, 10)
    assert sorted(map(tuple, T4h.tolist())) == sorted(map(tuple, outer_table(C2_CHARS, T4).tolist()))


def _check_nh3():
    G = C3V
    F = nh3_force(0.3, 1.0)
    d = math.radians(NH3["D"])
    cD, sD = math.cos(d), math.sin(d)
    k_, j_ = 0.3, 1.0
    r3 = math.sqrt(3)
    rowH1 = [1.5 * k_ * cD ** 2 + j_, .75 * k_ * cD ** 2, .75 * k_ * cD ** 2, 0, r3 / 4 * k_ * cD, -r3 / 4 * k_ * cD,
             1.5 * k_ * sD * cD, .75 * k_ * sD * cD, .75 * k_ * sD * cD, -j_ * cD, 0, j_ * sD]
    assert _close(F[0], rowH1)                                             # (3.7.13a) 第一列（「s4」→ 4）
    assert _close(F[3], [0, -r3 / 4 * k_ * cD, r3 / 4 * k_ * cD, k_ / 2, -k_ / 4, -k_ / 4,
                         0, -r3 / 4 * k_ * sD, r3 / 4 * k_ * sD, 0, 0, 0])
    assert abs(F[9, 9] - 1.5 * j_ * cD ** 2) < 1e-12 and abs(F[11, 11] - 3 * j_ * sD ** 2) < 1e-12
    for g in G.names:
        Rg = nh3_rep(g)
        assert _close(Rg @ F, F @ Rg) and _close(Rg.T @ Rg, np.eye(12))
    e12 = np.eye(12)
    s2, s3, s6, h = math.sqrt(2), math.sqrt(3), math.sqrt(6), math.sqrt(1.5)
    assert _close(nh3_projector(0, 0, 0) @ e12[0] * s3, (e12[0] + e12[1] + e12[2]) / s3)                # (3.7.1)
    assert _close(nh3_projector(2, 0, 0) @ e12[0] * h, (2 * e12[0] - e12[1] - e12[2]) / s6)
    assert _close(nh3_projector(2, 1, 0) @ e12[0] * h, (e12[1] - e12[2]) / s2)
    assert _close(nh3_projector(1, 0, 0) @ e12[0], 0)                       # (3.7.2)：P² 消滅 |H1⟩
    assert _close(nh3_projector(1, 0, 0) @ e12[3] * s3, (e12[3] + e12[4] + e12[5]) / s3)                # (3.7.5)
    assert _close(nh3_projector(2, 0, 1) @ e12[3] * h, (-e12[4] + e12[5]) / s2)          # −|rH2⟩（書上印 rH1）
    assert _close(nh3_projector(2, 1, 1) @ e12[3] * h, (2 * e12[3] - e12[4] - e12[5]) / s6)
    assert _close(nh3_projector(2, 0, 0) @ e12[9], e12[9]) and _close(nh3_projector(2, 1, 0) @ e12[9], e12[10])  # (3.7.6)
    assert _close(nh3_projector(0, 0, 0) @ e12[11], e12[11])              # (3.7.7)
    for v in nh3_cvectors().values():                                      # (3.7.12)：零平移、零角動量
        Tq, Aq = nh3_constraints(v)
        assert _close(Tq, 0) and _close(Aq, 0)
    cv = nh3_cvectors()
    assert abs(cv["c3"] @ nh3_mass() @ cv["c3"] - NH3["m"]) < 1e-12          # 圖上的歸一化 ⇒ ⟨c³|m|c³⟩ = m（(3.7.21)）
    rng = np.random.default_rng(3)
    for _ in range(3):
        kk, jj, mm, MM, DD = rng.uniform(0.2, 2), rng.uniform(0.5, 2), rng.uniform(0.5, 2), rng.uniform(5, 20), rng.uniform(10, 40)
        red = nh3_reduced(kk, jj, mm, MM, DD)
        fo, bo = nh3_reduced_formulas(kk, jj, mm, MM, DD), nh3_reduced_formulas(kk, jj, mm, MM, DD, book=True)
        for key in ("A1", "E"):
            assert all(_close(red[key][i], fo[key][i]) for i in range(3))   # 直接投影 ＝ 更正後的公式
        assert not _close(red["E"][0], bo["E"][0]) and _close(red["E"][1], bo["E"][1]) and not _close(red["E"][2], bo["E"][2])
        assert not _close(red["A1"][0], bo["A1"][0]) and not _close(red["A1"][1], bo["A1"][1]) and _close(red["A1"][2], bo["A1"][2])
        sec = nh3_secular(kk, jj, mm, MM, DD)
        num = nh3_numeric(kk, jj, mm, MM, DD)
        assert _close(num[:6], 0, 1e-8)                                     # 6 個零頻（平移、轉動）
        assert _close(num[6:], sorted([sec["w1m"] ** 2, sec["w3m"] ** 2, sec["w3m"] ** 2, sec["w3p"] ** 2,
                                       sec["w3p"] ** 2, sec["w1p"] ** 2]), 1e-8)
        assert abs(np.trace(red["E"][2]) - sec["S3"]) < 1e-9 and abs(np.linalg.det(red["A1"][2]) - sec["P1"]) < 1e-9
    om = nh3_omega_obs()
    assert {k: f"{v / 1e14:.4f}" for k, v in om.items()} == {"p1": "6.2857", "m1": "1.7104", "p3": "6.4308", "m3": "3.0666"}
    assert all(abs(2 * math.pi * 3e10 * NH3_NU[k] / 1e14 - NH3_OMEGA_BOOK[k]) < 6e-4 for k in NH3_NU)   # 書上用 c = 3×10¹⁰
    fit, fitb = nh3_fit(), nh3_fit(book_k=True)
    assert f"{fit['R']:.4f}" == "0.2972" and f"{fit['D']:.2f}" == "16.22" and f"{fit['j']:.1f}" == "847.5"
    assert f"{fit['k']:.1f}" == "-53.4" and f"{fitb['k']:.1f}" == "-62.1"
    q = 3 / 14
    Rof = lambda dd: 2 * (1 + q) * math.sin(dd) ** 2 / (1 + (q / 2 - 0.5) * math.cos(dd) ** 2)   # noqa: E731（(3.7.30)）
    assert abs(Rof(math.radians(fit["D"])) - fit["R"]) < 1e-12              # 反代：D_fit 讓 (3.7.30) 等於實測的 R
    assert abs(Rof(math.radians(17)) - fit["R"]) > 0.02                     # 書上的 17° 給不出實測的 R
    sec20 = nh3_secular(fit["k"] / M_H_BOOK, fit["j"] / M_H_BOOK, 1.0, 14.0, 20)   # 反代：擬合的 k、j 讓 S₁、S₃ 等於實測
    assert abs(sec20["S1"] / (om["p1"] ** 2 + om["m1"] ** 2) - 1) < 1e-12 and abs(sec20["S3"] / (om["p3"] ** 2 + om["m3"] ** 2) - 1) < 1e-12
    rf = nh3_ratio_fit()
    assert f"{rf['kj']:.4f}" == "0.2956" and f"{rf['w1p']:.3f}" == "8.884" and f"{rf['w3p']:.3f}" == "7.513"
    assert abs(rf["w3m"] - om["m3"] / 1e14) < 1e-9 and abs(rf["w1m"] - om["m1"] / 1e14) < 1e-9
    assert f"{100 * (rf['w1p'] / (om['p1'] / 1e14) - 1):.1f}" == "41.3" and f"{100 * (rf['w3p'] / (om['p3'] / 1e14) - 1):.1f}" == "16.8"
    xs = np.linspace(0.005, 0.1, 400)                                      # 偶然簡併：二分法 ＝ 掃描的變號
    dif = np.array([nh3_scan(x)[0] - nh3_scan(x)[2] for x in xs])
    i0 = int(np.flatnonzero(np.diff(np.sign(dif)))[0])
    assert xs[i0] <= nh3_crossing() <= xs[i0 + 1] and f"{nh3_crossing():.4f}" == "0.0348"
    s3_ = nh3_secular(NH3_KJ, 1.0)
    assert [f"{s3_[k]:.4f}" for k in ("S1", "P1", "S3", "P3", "w1p", "w1m", "w3p", "w3m")] == \
        ["1.9251", "0.1278", "1.5446", "0.2939", "1.3624", "0.2624", "1.1499", "0.4714"]
    assert _close(np.sqrt(nh3_numeric(NH3_KJ, 1.0)[6:]), [s3_["w1m"], s3_["w3m"], s3_["w3m"], s3_["w3p"], s3_["w3p"], s3_["w1p"]], 1e-9)
    assert [f"{v:.4f}" for v in nh3_scan(0)] == ["1.0125", "0.0000", "1.0462", "0.0000"]   # k → 0：√K、√H
    assert f"{M_H_BOOK * 6.3e14 ** 2:.0f}" == "663"                       # 書上的粗估 j ≈ 660 nt/m
    modes = nh3_modes()
    assert [nm for nm, _, _ in modes] == ["1+", "1-", "3+", "3-"]
    assert _close([w for _, w, _ in modes], [s3_["w1p"], s3_["w1m"], s3_["w3p"], s3_["w3m"]])
    Fm, Mm_ = nh3_force(NH3_KJ, 1.0), nh3_mass()
    for _, w, v in modes:
        assert _close(Fm @ v, w * w * Mm_ @ v, 1e-9)


def _check_problems():
    p311 = problem_311()
    assert p311["a"] == [["1", "s3"], ["r", "s1"], ["r^2", "s2"]] and p311["b"] == [["1", "s3"], ["r", "s2"], ["r^2", "s1"]]
    assert p311["c"] == p311["d"] == [["1", "r", "r^2"], ["s1", "s2", "s3"]]
    assert p311["e"] == [["1", "r", "s1", "s2"], ["r^2", "s3"]] and p311["f"] == [C3V.names]
    assert problem_312() == {"a": False, "b": True}
    assert all(nn * cc == Gx.order for Gx, lst in zip((D2, D3, D4), problem_322().values()) for _, nn, cc in lst)
    p331 = problem_331()
    assert p331["spectral"] == sp.Matrix([[1, 1, 2, 2, 2], [1, 1, -2, 2, -2], [1, 1, 2, -2, -2], [1, 1, -2, -2, 2], [1, -1, 0, 0, 0]])
    assert len(p331["split_D2"]) == 2 and len(p331["split_C4"]) == 2
    p332 = problem_332()
    assert p332["ells"] == [1, 1, 1, 1, 2] and p332["order2"] == {"Q": 1, "D4": 5}
    assert sorted(map(tuple, p332["chars"].tolist())) == sorted(map(tuple, D4.char_table_class_algebra()[1].tolist()))
    assert QUAT.is_rep(p332["irreps"][-1]) and p332["irreps"][-1]["k"] == sp.diag(I, -I)
    assert _close(su2(quat((1, 0, 0), 180)), QUAT.num["i"]) and _close(su2(quat((0, 0, 1), 180)), QUAT.num["k"])   # 3.3.2(a)
    pg, ph = problem_341(**P341["g"]), problem_341(**P341["h"])
    assert pg["omega2"] == [0, 0, 0, 6, 10, 10, 10, 16] and _close(pg["numeric"], [0, 0, 0, 6, 10, 10, 10, 16])
    assert [_f4(sp.sqrt(v)) for v in ph["omega2"]] == ["0.0000", "0.6892", "0.6892", "2.6458", "3.1623", "3.2442", "3.2442", "4.1231"]
    assert _close(ph["numeric"], [float(v) for v in ph["omega2"]])
    p351 = problem_351()
    assert p351["ok"] and p351["eig"] == [[1], [-1, 2], [-3, 0, 3]]
    p352 = problem_352()
    assert p352["great"] and p352["completeness"] == sp.diag(6, 3, 2) and p352["orthogonality"] == 6 * sp.eye(3)
    assert C3V.commutator_subgroup() == ["1", "r", "r^2"] and D4.commutator_subgroup() == ["1", "R^2"]
    assert problem_355() == {"C3v": dict(group=2, algebra=4, commuting=3, rank=4), "C4v": dict(group=2, algebra=6, commuting=4, rank=6)}
    assert all(all(v) for v in problem_356().values())
    gen = FiniteGroup.from_generators({"r": rot2(120), "s": refl2(0)})   # 由生成元封閉 ＝ C3v
    assert gen.order == 6 and gen.char_table_class_algebra()[1] == sp.Matrix(C3V_CHARS_BOOK)
    gen4 = FiniteGroup.from_generators({"R": rot2(90), "s": refl2(0)})
    assert gen4.order == 8 and _canon_table(gen4) == _canon_table(D4)


def _canon_table(G: FiniteGroup) -> list:
    """特徵標表的標準形（行依 (類大小, 元素的階) 排、列排序）：比較兩個群的表是否相同（差一個重新命名）。"""
    T = G.char_table_class_algebra()[1]
    sig = [(len(C), G.order_of(C[0])) for C in G.classes]
    cols = sorted(range(len(sig)), key=lambda k: sig[k])
    return sorted(tuple(T[a, c] for c in cols) for a in range(T.rows))


def check():
    """每個上卡片／講義的數字用兩種以上獨立算法核對（約 5 秒）。"""
    irr, Pl, ids, R_r = _check_groups()
    _check_six_and_chars(irr, Pl, ids, R_r)
    _check_d4_d6()
    _check_nh3()
    _check_problems()


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
    print(f"sym_nonabel（第 3 章）：C3v／D3／D4／C4v／D6／D3h／D3d／D2h／D4h／Q 的乘法表、類、類代數與冪等元、子群分裂的基本算子與不可約表示"
          f"（＝幾何矩陣）、特徵標表（類代數＝跡＝外積）、正則表示的約化、頻率與跡公式、相關表、Hamilton 轉（弧＝四元數＝SU(2)＝兩次反射）、"
          f"六自由度彈簧質量模型與圓周運動、D4／D6 轉子能帶的對稱標記、NH3 的 12×12 簡正模＝2×2 久期方程與擬合、習題 3.1.1–3.5.6 全部核對通過"
          f"（{time.time() - t0:.1f} 秒）。")
