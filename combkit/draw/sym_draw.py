# 由 combkit/tools/sync_from_harter_sym.py 從 harter_sym 專案的 codes/utils/sym_draw.py 複製；請改 harter_sym 那份再同步。
"""
sym_draw.py —— Harter PSDS 第 1–8 章共用的畫圖元件（ManimGL，`from manimlib import *`）。數學在 sym_mat.py（先算、再畫）。

放置方式：元件都用參數（font_size、unit、width…）決定大小，不要事後 .scale()（線寬不會跟著縮）。
顏色：KET_COLOR 青（右／ket 本徵向量）、BRA_COLOR 橘（左／bra）、P1_COLOR 藍（ε = 1 的投影）、P2_COLOR 金（ε = 5）、
      BEAM_COLOR 淺藍（粒子束）、BOX_COLOR 灰（分析器）、MASS_COLOR、SPRING_COLOR；WRONG_COLOR／RIGHT_COLOR 同 walk_draw。

    entry_tex(x)                                  sympy 數／int／str → TeX 字串（sympy 用 sp.latex）
    mat_tex(rows, colors=None, font_size=34, h_buff=0.45, v_buff=0.28)
                                                  圓括號矩陣 VGroup(格子們, 左括號, 右括號)；rows 是 list of list（或 sympy Matrix）；
                                                  colors：{(i, j): 顏色}、{"row": {i: 顏色}}、{"col": {j: 顏色}} 可混用（0 起算）
    mat_entry(m, i, j)、mat_row(m, i)、mat_col(m, j)   取出格子（0 起算）
    labeled_mat(lhs_tex, rows, **kw)              "lhs =" ＋矩陣 → VGroup(Tex, 矩陣)
    analyzer(angle=0, width=1.6, height=1.0, label=None, recombine=False)
                                                  Stern–Gerlach 分析器盒（盒內小箭頭指出分析器軸的角度，度）；recombine=True 是反轉的合併盒
    port(an, side, k=0)                           side "in"/"out"；k = 0 中央、1 上（通道 1）、2 下（通道 2）。
                                                  一般分析器：入口 k=0、出口 k=1,2；合併盒：入口 k=1,2、出口 k=0
    beam(start, end, intensity=1.0, color=BEAM_COLOR, max_width=14)
                                                  粒子束：線寬與不透明度隨 intensity（0..1）；intensity = 0 時是很淡的細線
    split_beams(an, p1, length=1.4, spread=0.35)  從分析器兩個出口往右的兩束（強度 p1、1 − p1）→ VGroup(上束, 下束)
    vec_plane(unit=1.0, xr=3, yr=3)               淡色座標平面（NumberPlane；plane.c2p 照用）
    vec_arrow(plane, v, color, label_tex=None)    原點 → v 的箭頭（label 放在箭頭尖外側）→ VGroup(箭頭[, 標籤])
    eigen_line(plane, v, color, extent=3)         通過原點、沿 v 的虛線（本徵方向）
    circle_image(plane, M, n=72, alpha=1.0)       單位圓經 (1 − α)1 + αM 變形後的曲線（α = 0 是單位圓）——用 UpdateFromAlphaFunc 做「變形」動畫
    contour_around(axes, center, r, color)        plot_axes 上繞實軸上一點的逆時針圓＋箭頭（H109 的圍道）
    SpringChain(n=2, spacing=2.2, mass_widths=None, unit=1.0, coils=7, amp=0.16)
                                                  Fig. 1.4.1 的牆—彈簧—質量鏈；.set_displacements([x1, x2]) 移動質量並重畫彈簧（用在 updater）；
                                                  .masses、.springs、.walls、.marks（平衡位置虛線）；.eq[i] = 第 i 個質量的平衡位置
    （另外從 combkit.draw.walk_draw 匯入：tex_table、table_cell、erratum_box、plot_axes、curve、pole_marker）

第 2 章新增（BOB_COLOR 擺錘、ROD_COLOR 擺桿、WEAK_COLOR／STRONG_COLOR 軟／硬彈簧、PHASOR_COLOR 相量）：
    PendulumRow(n=2, spacing=2.0, length=2.2, bob_radius=0.18, bob_radii=None, bob_colors=None, unit=1.0,
                spring_height=0.45, coils=6, amp=0.1, spring_colors=None, bar_pad=0.6)
                                                  一排單擺（Fig. 2.2.1、2.5.1）：頂上一根橫桿、n 個支點、擺桿、擺錘；相鄰擺桿在 spring_height
                                                  （擺長的比例）處用彈簧相連。.set_displacements([x1, …]) 把第 i 擺轉 unit·x_i 弧度並重畫彈簧
                                                  （用在 updater）；.bobs、.rods、.springs、.bar、.marks（鉛直虛線）、.pivots（支點座標）
    PendulumRing(n=3, radius=2.0, bob_radius=0.16, unit=0.5, mode="tangential", clockwise=False, start_angle=PI/2,
                 coils=5, amp=0.07, spring_colors=None, spring_widths=None, show_circle=True, labels=False)
                                                  俯視的單擺環（Fig. 2.6.1、2.7.1、2.7.6、2.12.3）：第 p 個錘在角度 start_angle ± 2πp/n；
                                                  .set_displacements(xs) 沿切向（mode="tangential"）或徑向（"radial"）移動 unit·x_p；
                                                  彈簧 p 連 p 與 p+1（spring_colors／spring_widths 可逐根給，對稱破缺用）；
                                                  .bobs、.springs、.marks（平衡位置小刻度）、.circle、.rest[p]（平衡位置）、.labels（labels=True 時 1..n）
    spring(a, b, coils=6, amp=0.1, color=SPRING_COLOR, width=2.5)、spring_points(a, b, coils, amp)
                                                  任意方向的一根彈簧（D₂ 四擺的對角彈簧等）；動畫時 updater 裡 set_points_as_corners(spring_points(…))
    PhasorClock(radius=0.6, color=PHASOR_COLOR, ticks=12, label_tex=None, font_size=24, show_real=True)
                                                  相量鐘面（Fig. 2.3.7、2.6.2）：.set_value(z) 把指針指到 center + radius·(Re z, Im z)（|z| ≤ 1）；
                                                  show_real=True 時在水平直徑上標出實部（＝位置）的小點；.face、.hand、.dot、.real_dot
    band_dots(axes, ks, values, color=YELLOW, radius=0.05)       能帶／色散的點（≤ 60 個；Fig. 2.7.4、2.7.7、2.12.1）→ VGroup
    band_curve(axes, f, k0=-PI, k1=PI, n=60, color=YELLOW, stroke_width=3)
                                                  y = f(k) 的曲線（一條 VMobject、n 點）
    level_lines(axes, values, x0, x1, color=YELLOW, stroke_width=3)
                                                  能階的水平短線（x 從 x0 到 x1，y = 各值；Fig. 2.7.3 的 ω² 譜）
    quad_contours(plane, F, levels=(0.5, 1.0, 2.0), n=60, color=TEAL, stroke_width=2.5)
                                                  ½xᵀ𝓕x = level 的等位橢圓（𝓕 對稱正定，2×2；Fig. 2.3.4）→ VGroup（每個 level 一條閉曲線）

第 3 章新增（MIRROR_COLOR 鏡線、N_COLOR／H_COLOR 氮／氫、LABEL_COLORS {A_1, A_2, B_1, B_2, E, E_1, E_2: 顏色}；
SIX_ORDER、SIX_COORD_DEG、SIX_COORD_MASS ＝ sym_nonabel 的 C3V.names、SIX_ANGLES、SIX_MASS——這裡另抄一份，因為 combkit 的 sym_draw 不 import sym_nonabel）：
    mirror_lines(radius=2.0, angles=(120, 60, 0), labels=(σ₁, σ₂, σ₃), color=MIRROR_COLOR, font_size=26, center=ORIGIN)
                                                  固定在紙上的鏡線（Fig. 3.1.2）→ VGroup(lines, labels)（.lines、.labels）
    BobTriangle(radius=1.6, positions=None, names=("a", "b", "c"), bob_radius=0.2, font_size=26, show_rest=True)
                                                  俯視的三個擺錘；.state = {名: (x, y)}（以 radius 為單位，預設 ＝ SN.BOB_START）；
                                                  .set_state(state)、.anim_to(state) → 動畫列表（新 state 用 SN.config_after(g, tri.state) 算）；
                                                  .bobs（每個是 VGroup(圓, 字母)）、.rest、.ring、.point(p)
    SpringTriangle(radius=1.3, wall_radius=2.3, unit=0.5, mass_radius=0.15, alpha_deg=15, show_k2=True, coils=5, amp=0.07)
                                                  Fig. 3.4.1 的六自由度模型；.set_displacements(x6)（正則座標 1, r, r², σ₁, σ₂, σ₃；給 updater）、
                                                  .coord_arrows(length, font_size)（⟨g|x⟩ 箭頭）、.displacement_arrows(x6, scale)（靜態位移箭頭）、
                                                  .disp(x6)、.masses、.k1、.k2、.wall、.rest
    sym_table(col_labels, row_labels, rows, cell_w=0.62, row_h=0.42, font_size=24, corner="", first_w=None)
                                                  乘法表／類代數表／特徵標表／相關表（格子：TeX 字串、int、sympy）；table_cell(t, i, j, len(col_labels)+1)
    perm_grid(M, cell=0.3, on_color=YELLOW, off_color=GREY_E, row_labels=None, col_labels=None, font_size=18, show_values=False)
                                                  矩陣的非零格子圖（𝓡(g)、T†𝓡T）；.cells[i][j]
    3D 一律投影到 2D 畫（不要用 Sphere／Surface／TrueDot／3D 相機：這台 wgpu 渲染器在它們緊接 VMobject 之後畫時會 GPUValidationError）：
    View3D(theta=−30, phi=65, scale=1.0, center=ORIGIN)   正交投影（phi = 0 俯視、90 平視）；.project(v)、.depth(v)、.matrix()
    SphereView(radius=2.0, view=None, lats=(−60, −30, 0, 30, 60), n_lon=6)
                                                  Hamilton 轉的球（經緯線前亮後淡）；.add_arc(N1, N2, color, width=6, tip=True)（半角弧，
                                                  N1、N2 用 SN.hamilton_arc／hamilton_product 算；ShowCreation(arc[0])）、.add_axis(axis, color)、
                                                  .add_point(v, color)、.project(v)、.set_view(theta, phi)（鏡頭繞轉：updater 裡呼叫）
    NH3Model(D_deg=20, bond=1.6, unit=0.5, n_radius=0.26, h_radius=0.16, view=None, show_hh=True)
                                                  球棍氨分子（投影）；.set_displacements(x12)（座標 ＝ SN.nh3_geometry；給 updater）、
                                                  .set_view(theta, phi)、.positions3(x12)、.atoms（N, H1, H2, H3）、.bonds、.hh；前後用 z_index 排
                                                  ——SphereView、NH3Model 的位置用 View3D(center=…) 給，不要事後 move_to（重畫時會回到 view.center）
    polar_curve(f, radius=1.5, amp=0.3, center, n=60)   環上的波 r(φ) = radius + amp·f(φ)（Fig. 3.6.2、3.6.5）

第 4 章新增（DIAG_COLORS {1..4: 顏色}、POS_COLOR／NEG_COLOR／IPOS_COLOR／INEG_COLOR 振幅的相位、X_COLOR／F_COLOR、OCT_COLORS（LABEL_COLORS 加 T_1、T_2）；
OCT_SITES、OCT_DIAGONALS、XF6_COORDS ＝ sym_oct 的 SITES、DIAGONALS、XF6_COORDS_LIST——另抄一份，因為 combkit 的 sym_draw 不 import sym_oct）：
    irrep_color(label)                           不可約表示標記 → 顏色（"T_{1u}"、"E_g"、"A_1"；g/u 不分）
    rot_matrix(axis, deg)                        繞 axis 逆時針 deg 度的 3×3（浮點）；轉動動畫 R(t) = rot_matrix(axis, t·deg) @ R₀
    CubeView(size=1.3, view=None, show_cube=True, show_octa=False, show_diagonals=True, diag_labels=True, font_size=22)
                                                  投影的立方體＋四條上色的體對角線（端點標 1–4）＋可選內接八面體；.set_rotation(Rm)（物體轉）、
                                                  .set_view(theta, phi)（鏡頭）、.add_axis(axis, color, length, label_tex)（實驗室的轉軸）、
                                                  .add_arrow(v, color, label_tex, body=False)；.edges、.diagonals、.diag_labels、.octa、.R
    OctaSites(radius=1.6, view=None, ball=0.2, labels=True, show_edges=True, font_size=22, label_texts=None)
                                                  八面體六頂點 |1⟩…|6⟩（+z, −z, +x, −x, +y, −y）；.set_amplitudes(v, min_frac, value_texs)
                                                  （球大小 ∝ |v_i|、顏色依相位 +、−、+i、−i）、.set_view、.point(i)、.balls、.edges、.labels、.values
    level_diagram(columns, links=(), col_w=1.0, gap=1.4, unit=0.5, font_size=22, deg_gap=0.07, colors=None, headers=None)
                                                  能階分裂／相關圖：columns ＝ [[(能量, 標記 TeX, 簡併度), …], …]，links ＝ [((欄, i), (欄, j))]；
                                                  .levels[c][i]、.labels[c][i]、.links、.headers
    XF6Model(bond=1.5, unit=0.6, view=None, x_radius=0.26, f_radius=0.17, show_ff=True, coords=None)
                                                  球棍 XF6：.set_displacements(x21)（座標＝ sym_oct.XF6_NAMES；給 updater）、.arrows(x21, scale, color)、
                                                  .set_view、.positions3、.atoms（X, F0…F5）、.bonds、.ff
    ——CubeView、OctaSites、XF6Model 都用 View3D(center=…) 定位置，不要事後 move_to。

第 5 章新增（BODY_COLORS {0, 1, 2: x̄ ȳ z̄ 的顏色}、LAB_COLOR、CRANK_COLOR、AMP_POS_COLOR／AMP_NEG_COLOR、PROB_COLOR）：
    TriadBall(radius=1.4, view=None, show_lab=True, lab_labels=("x","y","z"), body_labels=(r"\bar x", r"\bar y", r"\bar z"), axis_len=1.35, font_size=24)
                                                  投影的球＋跟著轉的赤道／子午線＋本體三腳架＋實驗室虛線軸；.set_rotation(Rm)（可非正常）、.set_view、
                                                  .add_crank(axis, color, label_tex, length)（實驗室固定的曲柄）、.tip(k)、.trail(k, Rs)、.project(v)
    VectorCone(j, m, unit=0.5, view=None, beta=0, color=TEAL, show_axis=True, fill=0.12)
                                                  角動量向量錐（斜高 √(j(j+1))、高 m）；.set_beta(β) 繞 y 傾斜（給 updater）；.rim、.sides、.axis、.foot
    cone_family(j, unit=0.5, view=None, colors=None)  Fig. 5.4.4 的一整族圓錐（.cones[k] ↔ m = j − k）
    ProbBars(values, labels=None, width=0.32, gap=0.1, height=2.0, signed=False, color=PROB_COLOR, font_size=20, label_every=1)
                                                  機率／振幅長條；.set_values(v)（給 updater）、.bars、.labels、.baseline、.top(i)
    ladder_column(levels, width=1.0, unit=0.5, font_size=20)  能階梯：levels = [(E, 左 TeX, 右 TeX, 顏色)]（簡併度寫在右邊，不畫成多條線）
    weight_diagram(weights, unit=0.6, roots=None, axes=True, axis_labels=("H_1","H_2"), extent=2.2, …)
                                                  根圖／權圖：weights = {(x, y): 重數}（重數 > 1 加同心圈）、roots = [(x, y, TeX)] 畫成箭頭
    mnemonic_wheel(labels, radius=1.5, start_deg=180, direction=−1, font_size=24, colors=None, angles=None)
                                                  助記輪：.ring、.slots[k]、.angle(k)、.point(k, r)；wheel_marker(wheel, k, text, inner=True) 箭頭標
    ——TriadBall、VectorCone 也用 View3D(center=…) 定位置；ProbBars、ladder_column、weight_diagram、mnemonic_wheel 是平面元件，可以 shift／next_to。

第 6 章新增（A_COLOR／B_COLOR 兩個粒子、STRESS_COLOR、STRAIN_COLOR、EFIELD_COLOR、BFIELD_COLOR、STIM_COLOR／RESP_COLOR）：
    PairOcta(radius=1.2, view=None, ball=0.13, …)   OctaSites ＋ 粒子 A（大光暈）、B（小點）；.set_pair(i, j)（1 起算，|i,j⟩）、.halo、.dot
    PairGrid(n=6, cell=0.42, font_size=16, labels=True, values=None)
                                                  6×6 的 |i,j⟩ 格子（列 A 的 i、行 B 的 j）；.cells[i][j]、.texts[i][j]、.color_cells(M, palette)
    StrainCube(size=1.0, view=None, color=GREY_A, show_ref=True)
                                                  投影立方體：.set_strain(S)（頂點 (1+S)v；給 updater）、.face_arrows(T, scale, color, faces)
                                                  （面 2k/2k+1 = ±x_k，箭頭 ±T·x̂_k）、.face_center(f)、.vertex(v)
    field_glyph(kind=out/in/+/-, radius=0.12, color)   ⊙、⊗、帶號電荷
    current_loop(center, rx=0.5, ry=0.18, sense=1, color)       斜看的電流環＋箭頭
    trace_curve(axes, ts, ys, color, stroke_width)              plot_axes 上的時間曲線（n ≤ 1000，靜態；資料 sym_prod.lorentz_trace）
    SmithChart(S=1, gamma=0.1, unit=0.25, deltas=())            Γ 圓、Δ 圓、.point(Δ)、.phasor(Δ)（近共振 S/(Δ − iΓ)）
    hat_contours(center, unit=2, j, k, r, levels, n=60, colors) Jahn–Teller 下層 ε− 的等高線（自己用二分法算）
    hat_profile(axes, j, k, r, phi_deg, rho_max, n, colors)     ε± 沿 φ／φ+180° 的截面（Fig. 6.7.3）
    ——PairOcta、StrainCube 用 View3D(center=…) 定位置；其餘是平面元件。

第 7 章新增（P_COLOR／E_COLOR 質子、電子、J1_COLOR／J2_COLOR／J3_COLOR 三個角動量、POL_COLOR 偏振橢圓、STOKES_COLOR、
POT_G_COLOR／POT_U_COLOR 鍵結／反鍵結位能曲線）：
    SpinTop(radius=0.3, color=P_COLOR, arrow_len=0.55, angle=0, label_tex=None, font_size=24, show_orbit=True)
                                                  小陀螺：球＋自旋箭頭（angle：與向上的夾角，度；0 ↑、180 ↓）＋進動虛線橢圓；
                                                  .set_angle(deg)（給 updater）、.ball、.arrow、.orbit、.label
    SpinPair(gap=1.5, radii=(0.34, 0.22), colors=(P_COLOR, E_COLOR), labels=("p", "e"), font_size=24)
                                                  兩個 SpinTop（質子、電子；Fig. 7.1.1）；.set_ms(m1, m2)（±½ → ↑／↓；0 → 水平）、.tops
    vector_triangle(j1, j2, j3, unit=0.5, start=ORIGIN, angle=0, colors=(J1, J2, J3), labels=True, font_size=24)
                                                  向量三角形（Fig. 7.1.5）：邊長 √(j(j+1))·unit，j1 接 j2 的頭尾＝ j3（j3 沿 angle 度）→ VGroup(j1, j2, j3, 標籤)
    CouplingCone(j1, j2, J, M, unit=0.35, view=None, phase=0, colors=(J1, J2, J3), show_range=True)
                                                  Fig. 7.1.6 的錐：J（與 z 夾 Θ，cos Θ = M/|J|）＋ j1 繞 J 的錐（固定夾角）＋ j2（j1 尖 → J 尖）＋ z 軸上
                                                  古典範圍 m1 = |j1|cos(Θ ± θ13) 的兩個刻度；.set_phase(deg)（j1 繞 J 轉；給 updater）、.m1_range
    RESurface(paths, radii=1.6, view=None, colors=None, width=2.0, front=1.0, back=0.22, mesh=None, axes=True,
              axis_labels=(\\bar J_x, \\bar J_y, \\bar J_z), font_size=22)
                                                  RE 曲面的投影（Fig. 7.4.1、7.4.5、7.4.11、7.5.11–7.5.13）：paths ＝ [[n×3 單位向量折線…]（每個能量一組）]
                                                  （用 sym_coup.sphere_contours／rigid_paths／octa_paths／re46_paths 算），radii ＝ 每組的半徑（RE 曲面
                                                  是徑向圖，等能量線半徑固定）；mesh ＝ [n×3 點]（已乘半徑的淡色經緯線，見 RESurface.mesh_lines）；
                                                  .set_view(theta, phi)、.curves[k]、.axes；RESurface.mesh_lines(radius_fn, n_lat, n_lon, n)
    PoincareSphere(radius=1.5, view=None, labels=("A", "B", "C"), font_size=24)
                                                  SphereView ＋ A、B、C 三軸（A→x、B→y、C→z）；.add_stokes(S, color) → VGroup(箭頭, 點)、
                                                  .set_stokes(g, S)（給 updater）、.axis_labels
    polarization_ellipse(a, b, phi=0, center=ORIGIN, unit=1.0, color=POL_COLOR, n=60, arrow=True, width=3)
                                                  (7.5.20) 的橢圓：E = (a cos Φ cos φ − b sin Φ sin φ, a cos Φ sin φ + b sin Φ cos φ)·unit；
                                                  b > 0 逆時針（右旋，R > L）；→ VGroup(曲線[, 方向箭頭])
    ellipse_from_stokes(S, center=ORIGIN, unit=1.0, color=POL_COLOR, **kw)
                                                  Stokes (S_A, S_B, S_C) → ψ = ½asin(S_C/I)、φ = ½atan2(S_B, S_A)、a = √I cos ψ、b = √I sin ψ 的橢圓
    projected_curve(points3, view, color=TEAL, width=2.5, front=1.0, back=0.3)
                                                  3D 折線（振子環面上的 (x1, x2, p1) 軌跡等）投影、前後不同不透明度 → VGroup
    box_frame(X, Y, center=ORIGIN, unit=1.0, angle=0, color=GREY_B)   2X×2Y 的虛線盒（Fig. 7.5.5、7.5.7；angle 度）
    potential_curves(axes, xs, ys_list, colors, labels=None, label_x=None, font_size=22, width=3)
                                                  plot_axes 上的多條位能曲線＋右端 TeX 標籤（H₂⁺ 的 V±、H₂）→ VGroup(curves, labels)
    min_marker(axes, x, y, color=YELLOW, label_tex=None, font_size=20, drop=True)
                                                  極小點：點＋落到橫軸的虛線＋標籤 → VGroup
    ——CouplingCone、RESurface、PoincareSphere 用 View3D(center=…) 定位置，不要事後 move_to；其餘是平面元件。

第 8 章新增（TRAJ_COLOR 軌跡、FRONT_COLOR 波前、ENV_COLOR 包絡／焦散、PACKET_COLOR／PACKET_ENV_COLOR 波包、SPIN_COLOR 自旋向量、
DRESSED_COLOR／BARE_COLOR 穿衣與裸能階）：
    polyline(axes, pts, color, width, dashed=False)   plot_axes 上的一條折線（pts (n, 2)；包絡、焦散、單條軌跡、Rabi 曲線…）
    line_family(axes, lines, x0, x1, color, width, highlight=None)
                                                  直線族 Y = sX + b（lines ＝ sym_dyn.line_family；Fig. 8.1.2；highlight = 序號集合）
    trajectory_fronts(axes, trajs, fronts, …)     Fig. 8.1.6：軌跡族（sym_dyn.trajectory）＋ 等 S_H 波前（sym_dyn.fig816_curves 的 {level: [段]}）
                                                  → VGroup(軌跡們, 波前們)
    PhasePhasor(radius=1.2, center, amp=1, blob=0.22, color, labels, font_size)
                                                  相空間相量（Fig. 8.2.1）：虛線軌道、軸、箭頭、Gauss 圓斑；.set_alpha(z)（給 updater）、.point(z)
    snapshot_stack(xs, rows, unit_x, unit_y, gap, top, …)
                                                  波包快照列（Fig. 8.2.2、8.3.2(a)）：每列 Re χ 實線＋|χ| 虛線＋基線；rows 用 sym_dyn.packet_closed
    CrankSphere(radius=1.5, view, labels=("x","y","z"))
                                                  曲柄 ω 與自旋 S（Fig. 8.5.1）：.add_vector(v, color, label_tex)（|v| = ½ 在球面上，不歸一化）、
                                                  .set_vector(g, v)（給 updater）、.add_trail(points3)、CrankSphere.cone_rim(axis, 半頂角, 長度, n)
    dressed_ladder(blocks, unit, base, …)         JC 穿衣能階梯（Fig. 8.7.1）：blocks ＝ sym_dyn.jc_ladder(Δ, g, Nmax)；.dressed[k]
    ——CrankSphere 用 View3D(center=…) 定位置；其餘是平面元件。譜函數曲線族、Rabi 曲線族用第 7 章的 potential_curves（多條＋標籤）。

注意：Tex 的 color= 參數無效，一律 .set_color()；同一個 play 不要同時動畫父物件與它的子物件；Tex 裡不放中文（用 zh()）。
"""
import itertools

import numpy as np
import sympy as sp
from manimlib import *

from combkit.draw.styles import zh
from combkit.draw.walk_draw import (RIGHT_COLOR, WRONG_COLOR, curve, erratum_box, plot_axes,  # noqa: F401
                                    pole_marker, table_cell, tex_table)

KET_COLOR = TEAL
BRA_COLOR = ORANGE
P1_COLOR = BLUE
P2_COLOR = YELLOW
BEAM_COLOR = "#8FD3FF"
BOX_COLOR = GREY_B
MASS_COLOR = "#E07A5F"
SPRING_COLOR = GREY_A


# ---------------------------------------------------------------- 矩陣
def entry_tex(x) -> str:
    if isinstance(x, str):
        return x
    if isinstance(x, sp.Basic):
        return sp.latex(x)
    return str(x)


def mat_tex(rows, colors=None, font_size: float = 34, h_buff: float = 0.45, v_buff: float = 0.28,
            paren_color=WHITE) -> VGroup:
    if isinstance(rows, sp.MatrixBase):
        rows = rows.tolist()
    nr, nc = len(rows), len(rows[0])
    colors = colors or {}
    cells = [[Tex(entry_tex(x), font_size=font_size) for x in row] for row in rows]
    for i in range(nr):
        for j in range(nc):
            c = colors.get((i, j)) or colors.get("row", {}).get(i) or colors.get("col", {}).get(j)
            if c is not None:
                cells[i][j].set_color(c)
    col_w = [max(cells[i][j].get_width() for i in range(nr)) for j in range(nc)]
    row_h = [max(max(cells[i][j].get_height() for j in range(nc)), 0.3 * font_size / 34) for i in range(nr)]
    xs = np.cumsum([0] + [w + h_buff for w in col_w[:-1]]) + np.array(col_w) / 2
    ys = -(np.cumsum([0] + [h + v_buff for h in row_h[:-1]]) + np.array(row_h) / 2)
    grid = VGroup()
    for i in range(nr):
        for j in range(nc):
            grid.add(cells[i][j].move_to([xs[j], ys[i], 0]))
    grid.move_to(ORIGIN)
    H = grid.get_height() + 0.25
    sag = 0.12 + 0.03 * H
    r = (H * H / 4 + sag * sag) / (2 * sag)
    th = 2 * np.arcsin(H / (2 * r))
    left = Arc(start_angle=PI - th / 2, angle=th, radius=r).set_stroke(paren_color, 2.5)
    right = Arc(start_angle=-th / 2, angle=th, radius=r).set_stroke(paren_color, 2.5)
    left.move_to(grid.get_left() + 0.12 * LEFT, aligned_edge=RIGHT)
    right.move_to(grid.get_right() + 0.12 * RIGHT, aligned_edge=LEFT)
    m = VGroup(grid, left, right)
    m.n_cols = nc
    return m


def mat_entry(m: VGroup, i: int, j: int):
    return m[0][i * m.n_cols + j]


def mat_row(m: VGroup, i: int) -> VGroup:
    return VGroup(*[mat_entry(m, i, j) for j in range(m.n_cols)])


def mat_col(m: VGroup, j: int) -> VGroup:
    return VGroup(*[m[0][k] for k in range(j, len(m[0]), m.n_cols)])


def labeled_mat(lhs_tex: str, rows, font_size: float = 34, **kw) -> VGroup:
    lhs = Tex(lhs_tex + "=", font_size=font_size)
    return VGroup(lhs, mat_tex(rows, font_size=font_size, **kw)).arrange(RIGHT, buff=0.18)


# ---------------------------------------------------------------- Stern–Gerlach 分析器與粒子束
def analyzer(angle: float = 0, width: float = 1.6, height: float = 1.0, label: str | None = None,
             recombine: bool = False, color=BOX_COLOR, font_size: float = 26) -> VGroup:
    box = RoundedRectangle(width=width, height=height, corner_radius=0.12)
    box.set_stroke(color, 2.5).set_fill(GREY_E, 0.85)
    if recombine:
        box.set_stroke(color, 2.5).set_fill("#2A2A3A", 0.85)
    poles = VGroup(Polygon(LEFT * 0.28, RIGHT * 0.28, DOWN * 0.18).shift(UP * (height / 2 - 0.08)),
                   Rectangle(width=0.56, height=0.1).shift(DOWN * (height / 2 - 0.1)))
    poles.set_stroke(width=0).set_fill(GREY_C, 0.9)
    a = np.radians(angle)
    d = 0.28 * height * np.array([np.sin(a), np.cos(a), 0])
    axis = Arrow(-d, d, buff=0, thickness=2.5, max_tip_length_to_length_ratio=0.4).set_color(YELLOW)
    parts = VGroup(box, poles, axis)
    if label is not None:
        parts.add(Tex(label, font_size=font_size).set_color(GREY_A).next_to(box, DOWN, buff=0.12))
    parts.recombine = recombine
    return parts


def port(an: VGroup, side: str, k: int = 0) -> np.ndarray:
    box = an[0]
    dy = {0: 0.0, 1: 0.25, 2: -0.25}[k] * box.get_height()
    x = box.get_left()[0] if side == "in" else box.get_right()[0]
    return np.array([x, box.get_center()[1] + dy, 0.0])


def beam(start, end, intensity: float = 1.0, color=BEAM_COLOR, max_width: float = 14) -> Line:
    p = float(np.clip(intensity, 0, 1))
    return Line(start, end).set_stroke(color, width=max(1.0, max_width * p), opacity=0.15 + 0.85 * p)


def split_beams(an: VGroup, p1: float, length: float = 1.4, spread: float = 0.35, **kw) -> VGroup:
    out = VGroup()
    for k, p, s in ((1, p1, spread), (2, 1 - p1, -spread)):
        a = port(an, "out", k)
        out.add(beam(a, a + length * RIGHT + s * UP, p, **kw))
    return out


# ---------------------------------------------------------------- 平面上的向量與 2×2 矩陣
def vec_plane(unit: float = 1.0, xr: float = 3, yr: float = 3) -> NumberPlane:
    plane = NumberPlane(x_range=(-xr, xr, 1), y_range=(-yr, yr, 1), width=2 * xr * unit, height=2 * yr * unit,
                        background_line_style=dict(stroke_color=GREY_D, stroke_width=1, stroke_opacity=0.6),
                        faded_line_ratio=0)
    plane.axes.set_stroke(GREY_B, 2)
    return plane


def vec_arrow(plane, v, color=KET_COLOR, label_tex: str | None = None, font_size: float = 30) -> VGroup:
    v = [float(x) for x in v]
    tip = plane.c2p(v[0], v[1])
    arr = Arrow(plane.c2p(0, 0), tip, buff=0, thickness=4).set_color(color)
    g = VGroup(arr)
    if label_tex:
        d = tip - plane.c2p(0, 0)
        d = d / (np.linalg.norm(d) + 1e-9)
        g.add(Tex(label_tex, font_size=font_size).set_color(color).move_to(tip + 0.38 * d))
    return g


def eigen_line(plane, v, color=GREY_B, extent: float = 3.0) -> DashedLine:
    v = np.array([float(x) for x in v][:2])
    v = v / np.linalg.norm(v) * extent
    return DashedLine(plane.c2p(*(-v)), plane.c2p(*v), dash_length=0.12).set_stroke(color, 2)


def circle_image(plane, M, n: int = 72, alpha: float = 1.0, color=WHITE, stroke_width: float = 3) -> VMobject:
    A = np.array([[complex(x).real for x in row] for row in (M.tolist() if hasattr(M, "tolist") else M)], dtype=float)
    A = (1 - alpha) * np.eye(2) + alpha * A
    t = np.linspace(0, 2 * np.pi, n + 1)
    pts = A @ np.vstack([np.cos(t), np.sin(t)])
    c = VMobject().set_points_as_corners([plane.c2p(x, y) for x, y in pts.T])
    return c.set_stroke(color, stroke_width)


def contour_around(axes, center: float, r: float, color=TEAL, stroke_width: float = 3) -> VGroup:
    c0 = axes.c2p(center, 0)
    rad = np.linalg.norm(axes.c2p(center + r, 0) - c0)
    circ = Circle(radius=rad).move_to(c0).set_stroke(color, stroke_width)
    tip = ArrowTip(angle=PI / 4 + PI / 2).set_fill(color, 1).set_height(0.16)
    tip.move_to(c0 + rad * np.array([np.cos(PI / 4), np.sin(PI / 4), 0]))
    return VGroup(circ, tip)


# ---------------------------------------------------------------- 牆—彈簧—質量鏈（Fig. 1.4.1）
def _zigzag(a, b, coils: int, amp: float) -> list:
    a, b = np.array(a, dtype=float), np.array(b, dtype=float)
    L = b - a
    pts = [a, a + 0.12 * L]
    for k in range(2 * coils):
        s = 0.12 + 0.76 * (k + 0.5) / (2 * coils)
        pts.append(a + s * L + (amp if k % 2 == 0 else -amp) * UP)
    return pts + [a + 0.88 * L, b]


class SpringChain(VGroup):
    def __init__(self, n: int = 2, spacing: float = 2.2, mass_widths=None, unit: float = 1.0, coils: int = 7,
                 amp: float = 0.16, wall_height: float = 1.4, mass_color=MASS_COLOR, spring_color=SPRING_COLOR, **kw):
        super().__init__(**kw)
        self.n, self.unit, self.coils, self.amp = n, unit, coils, amp
        widths = mass_widths or [0.55] * n
        x_wall = (n + 1) * spacing / 2
        w0, w1 = np.array([-x_wall, 0.0, 0.0]), np.array([x_wall, 0.0, 0.0])
        self.walls = VGroup()
        for s, p in ((-1, w0), (1, w1)):
            wall = Line(p + wall_height / 2 * UP, p + wall_height / 2 * DOWN).set_stroke(GREY_B, 4)
            hatch = VGroup(*[Line(p + y * UP, p + y * UP + 0.18 * s * RIGHT + 0.18 * DOWN).set_stroke(GREY_C, 2)
                             for y in np.linspace(-wall_height / 2 + 0.2, wall_height / 2, 6)])
            self.walls.add(VGroup(wall, hatch))
        eq = self.eq
        self.masses = VGroup(*[Square(side_length=w).set_stroke(WHITE, 2).set_fill(mass_color, 0.9).move_to(eq[i])
                               for i, w in enumerate(widths)])
        self.marks = VGroup(*[DashedLine(e + 0.75 * UP, e + 0.75 * DOWN, dash_length=0.08).set_stroke(GREY_C, 1.5)
                              for e in eq])
        self.springs = VGroup(*[VMobject().set_stroke(spring_color, 2.5) for _ in range(n + 1)])
        self.add(self.marks, self.walls, self.springs, self.masses)
        self.set_displacements([0.0] * n)

    @property
    def wall_pts(self):
        """兩面牆的中點（隨整條鏈移動）。"""
        return [self.walls[0][0].get_center(), self.walls[1][0].get_center()]

    @property
    def eq(self):
        """各質量的平衡位置（隨整條鏈移動）。"""
        a, b = self.wall_pts
        return [a + (i + 1) / (self.n + 1) * (b - a) for i in range(self.n)]

    def set_displacements(self, xs):
        eq = self.eq
        for i, m in enumerate(self.masses):
            m.move_to(eq[i] + self.unit * float(xs[i]) * RIGHT)
        ends = [self.wall_pts[0]] + [m.get_center() for m in self.masses] + [self.wall_pts[1]]
        for k, spr in enumerate(self.springs):
            a = ends[k] + (self.masses[k - 1].get_width() / 2 * RIGHT if k > 0 else 0)
            b = ends[k + 1] - (self.masses[k].get_width() / 2 * RIGHT if k < self.n else 0)
            spr.set_points_as_corners(_zigzag(a, b, self.coils, self.amp))
        return self


# ================================================================ 第 2 章：單擺列、單擺環、相量鐘、能帶、等位線
BOB_COLOR = "#E07A5F"
ROD_COLOR = GREY_B
WEAK_COLOR = "#7FB7E6"
STRONG_COLOR = "#F2C14E"
PHASOR_COLOR = TEAL


def _zigzag_dir(a, b, coils: int, amp: float) -> list:
    """a → b 的彈簧折線，鋸齒垂直於 a→b（任意方向）。"""
    a, b = np.array(a, dtype=float), np.array(b, dtype=float)
    L = b - a
    d = np.linalg.norm(L[:2])
    nrm = np.array([-L[1], L[0], 0.0]) / (d if d > 1e-9 else 1.0)
    pts = [a, a + 0.12 * L]
    for k in range(2 * coils):
        s = 0.12 + 0.76 * (k + 0.5) / (2 * coils)
        pts.append(a + s * L + (amp if k % 2 == 0 else -amp) * nrm)
    return pts + [a + 0.88 * L, b]


def spring(a, b, coils: int = 6, amp: float = 0.1, color=SPRING_COLOR, width: float = 2.5) -> VMobject:
    """a → b 的一根彈簧（任意方向）；要動就在 updater 裡 .set_points_as_corners(spring_points(a, b, …))。"""
    return VMobject().set_points_as_corners(_zigzag_dir(a, b, coils, amp)).set_stroke(color, width)


def spring_points(a, b, coils: int = 6, amp: float = 0.1) -> list:
    return _zigzag_dir(a, b, coils, amp)


class PendulumRow(VGroup):
    """一排單擺（側視）。第 i 擺的偏角 = unit·x_i（弧度）；相鄰擺桿在 spring_height（擺長的比例）處以彈簧相連。"""

    def __init__(self, n: int = 2, spacing: float = 2.0, length: float = 2.2, bob_radius: float = 0.18,
                 bob_radii=None, bob_colors=None, unit: float = 1.0, spring_height: float = 0.45, coils: int = 6,
                 amp: float = 0.1, spring_colors=None, bar_pad: float = 0.6, **kw):
        super().__init__(**kw)
        self.n, self.unit, self.length, self.h, self.coils, self.amp = n, unit, length, spring_height, coils, amp
        x0 = -(n - 1) * spacing / 2
        piv = [np.array([x0 + i * spacing, length / 2, 0.0]) for i in range(n)]
        self.bar = Line(piv[0] + bar_pad * LEFT, piv[-1] + bar_pad * RIGHT).set_stroke(GREY_A, 5)
        self.marks = VGroup(*[DashedLine(p, p + (length + 0.3) * DOWN, dash_length=0.08).set_stroke(GREY_D, 1.5)
                              for p in piv])
        self.pivot_dots = VGroup(*[Dot(p, radius=0.05).set_fill(GREY_A, 1) for p in piv])
        self.rods = VGroup(*[Line(p, p + length * DOWN).set_stroke(ROD_COLOR, 3) for p in piv])
        radii = bob_radii or [bob_radius] * n
        cols = bob_colors or [BOB_COLOR] * n
        self.bobs = VGroup(*[Circle(radius=r).set_stroke(WHITE, 2).set_fill(c, 0.95) for r, c in zip(radii, cols)])
        scol = spring_colors or [SPRING_COLOR] * max(n - 1, 0)
        self.springs = VGroup(*[VMobject().set_stroke(c, 2.5) for c in scol])
        self.add(self.marks, self.bar, self.springs, self.rods, self.pivot_dots, self.bobs)
        self.set_displacements([0.0] * n)

    @property
    def pivots(self):
        """各支點座標（隨整排移動）。"""
        return [d.get_center() for d in self.pivot_dots]

    def bob_position(self, i: int, x: float):
        th = self.unit * float(x)
        return self.pivots[i] + self.length * np.array([np.sin(th), -np.cos(th), 0.0])

    def set_displacements(self, xs):
        pv = self.pivots
        for i in range(self.n):
            tip = self.bob_position(i, xs[i])
            self.rods[i].put_start_and_end_on(pv[i], tip)
            self.bobs[i].move_to(tip)
        for k, spr in enumerate(self.springs):
            a = pv[k] + self.h * (self.rods[k].get_end() - pv[k])
            b = pv[k + 1] + self.h * (self.rods[k + 1].get_end() - pv[k + 1])
            spr.set_points_as_corners(_zigzag_dir(a, b, self.coils, self.amp))
        return self


class PendulumRing(VGroup):
    """俯視的單擺環。第 p 個錘在角度 start_angle + s·2πp/n（clockwise=True 時 s = −1）；位移沿切向或徑向。"""

    def __init__(self, n: int = 3, radius: float = 2.0, bob_radius: float = 0.16, unit: float = 0.5,
                 mode: str = "tangential", clockwise: bool = False, start_angle: float = PI / 2, coils: int = 5,
                 amp: float = 0.07, spring_colors=None, spring_widths=None, show_circle: bool = True,
                 labels: bool = False, bob_color=BOB_COLOR, font_size: float = 22, **kw):
        super().__init__(**kw)
        self.n, self.radius, self.unit, self.mode, self.coils, self.amp = n, radius, unit, mode, coils, amp
        sgn = -1 if clockwise else 1
        self.angles = [start_angle + sgn * 2 * PI * p / n for p in range(n)]
        self.center_dot = Dot(ORIGIN, radius=0.01).set_opacity(0)
        self.circle = DashedVMobject(Circle(radius=radius), num_dashes=max(24, 6 * n)).set_stroke(GREY_D, 1.5)
        if not show_circle:
            self.circle.set_opacity(0)
        rest = self._rest_from(ORIGIN)
        self.marks = VGroup(*[Line(r - 0.12 * self._radial(p), r + 0.12 * self._radial(p)).set_stroke(GREY_C, 2)
                              for p, r in enumerate(rest)])
        self.bobs = VGroup(*[Circle(radius=bob_radius).set_stroke(WHITE, 2).set_fill(bob_color, 0.95).move_to(r)
                             for r in rest])
        cols = spring_colors or [SPRING_COLOR] * n
        wids = spring_widths or [2.5] * n
        self.springs = VGroup(*[VMobject().set_stroke(c, w) for c, w in zip(cols, wids)])
        self.labels = VGroup()
        if labels:
            for p, r in enumerate(rest):
                self.labels.add(Tex(str(p + 1), font_size=font_size).set_color(GREY_A).move_to(r + 0.45 * self._radial(p)))
        self.add(self.center_dot, self.circle, self.marks, self.springs, self.bobs, self.labels)
        self.set_displacements([0.0] * n)

    def _radial(self, p: int):
        a = self.angles[p]
        return np.array([np.cos(a), np.sin(a), 0.0])

    def _tangent(self, p: int):
        a = self.angles[p]
        return np.array([-np.sin(a), np.cos(a), 0.0])

    def _rest_from(self, c):
        return [np.array(c) + self.radius * self._radial(p) for p in range(self.n)]

    @property
    def rest(self):
        """各錘的平衡位置（隨整個環移動）。"""
        return self._rest_from(self.center_dot.get_center())

    def set_displacements(self, xs):
        rest = self.rest
        for p in range(self.n):
            d = self._tangent(p) if self.mode == "tangential" else self._radial(p)
            self.bobs[p].move_to(rest[p] + self.unit * float(xs[p]) * d)
        for p, spr in enumerate(self.springs):
            q = (p + 1) % self.n
            a, b = self.bobs[p].get_center(), self.bobs[q].get_center()
            u = (b - a) / max(np.linalg.norm(b - a), 1e-9)
            r = self.bobs[p].get_width() / 2
            spr.set_points_as_corners(_zigzag_dir(a + r * u, b - r * u, self.coils, self.amp))
        return self


class PhasorClock(VGroup):
    """相量鐘面：指針指到 center + radius·(Re z, Im z)。e^{−iωt} 的相量順時針轉。"""

    def __init__(self, radius: float = 0.6, color=PHASOR_COLOR, ticks: int = 12, label_tex: str | None = None,
                 font_size: float = 24, show_real: bool = True, **kw):
        super().__init__(**kw)
        self.radius = radius
        self.face = Circle(radius=radius).set_stroke(GREY_B, 2).set_fill(GREY_E, 0.6)
        self.ticks = VGroup(*[Line((radius - 0.08) * np.array([np.cos(a), np.sin(a), 0]),
                                   radius * np.array([np.cos(a), np.sin(a), 0])).set_stroke(GREY_C, 1.5)
                              for a in np.linspace(0, 2 * PI, ticks, endpoint=False)])
        self.axis = Line(radius * LEFT, radius * RIGHT).set_stroke(GREY_D, 1)
        self.hand = Line(ORIGIN, radius * RIGHT).set_stroke(color, 4)
        self.dot = Dot(radius * RIGHT, radius=0.05).set_fill(color, 1)
        self.real_dot = Dot(radius * RIGHT, radius=0.045).set_fill(YELLOW, 1)
        if not show_real:
            self.real_dot.set_opacity(0)
        self.label = VGroup()
        if label_tex:
            self.label = Tex(label_tex, font_size=font_size).set_color(GREY_A).next_to(self.face, DOWN, buff=0.12)
        self.add(self.face, self.ticks, self.axis, self.hand, self.dot, self.real_dot, self.label)
        self.set_value(1.0)

    def set_value(self, z):
        z = complex(z)
        c = self.face.get_center()
        tip = c + self.radius * np.array([z.real, z.imag, 0.0])
        if abs(z) > 1e-6:
            self.hand.put_start_and_end_on(c, tip)
            self.hand.set_stroke(opacity=1)
        else:
            self.hand.put_start_and_end_on(c, c + 1e-4 * RIGHT)
            self.hand.set_stroke(opacity=0)
        self.dot.move_to(tip)
        self.real_dot.move_to(c + self.radius * z.real * RIGHT)
        return self


def band_dots(axes, ks, values, color=YELLOW, radius: float = 0.05) -> VGroup:
    """能帶／色散的點 (k, 值)（≤ 60 個）。"""
    return VGroup(*[Dot(axes.c2p(float(k), float(v)), radius=radius).set_fill(color, 1) for k, v in zip(ks, values)])


def band_curve(axes, f, k0: float = -PI, k1: float = PI, n: int = 60, color=YELLOW, stroke_width: float = 3) -> VMobject:
    """y = f(k) 的曲線（n 點，一條 VMobject）。"""
    ks = np.linspace(k0, k1, n)
    return VMobject().set_points_as_corners([axes.c2p(float(k), float(f(k))) for k in ks]).set_stroke(color, stroke_width)


def level_lines(axes, values, x0: float, x1: float, color=YELLOW, stroke_width: float = 3) -> VGroup:
    """能階的水平短線（y = 各值，x 從 x0 到 x1）。"""
    return VGroup(*[Line(axes.c2p(x0, float(v)), axes.c2p(x1, float(v))).set_stroke(color, stroke_width) for v in values])


def quad_contours(plane, F, levels=(0.5, 1.0, 2.0), n: int = 60, color=TEAL, stroke_width: float = 2.5) -> VGroup:
    """½xᵀ𝓕x = level 的等位橢圓（𝓕 對稱正定 2×2，可給 sympy Matrix）。"""
    rows = F.tolist() if hasattr(F, "tolist") else F
    A = np.array([[float(complex(x).real) for x in row] for row in rows], dtype=float)
    w, V = np.linalg.eigh(A)
    t = np.linspace(0, 2 * np.pi, n + 1)
    out = VGroup()
    for lv in levels:
        r = np.sqrt(2 * float(lv) / w)
        pts = V @ np.vstack([r[0] * np.cos(t), r[1] * np.sin(t)])
        out.add(VMobject().set_points_as_corners([plane.c2p(x, y) for x, y in pts.T]).set_stroke(color, stroke_width))
    return out


# ---------------------------------------------------------------- 第 3 章（非 Abel 群）
MIRROR_COLOR = "#9AD1D4"
N_COLOR = "#5B8DEF"
H_COLOR = GREY_A
LABEL_COLORS = {"A_1": BLUE, "A_2": GREEN, "B_1": YELLOW, "B_2": ORANGE, "E": TEAL, "E_1": TEAL, "E_2": PURPLE_B}
SIX_ORDER = ["1", "r", "r^2", "s1", "s2", "s3"]
SIX_COORD_DEG = {"1": -45, "r": 75, "r^2": 195, "s1": 285, "s2": 165, "s3": 45}      # ＝ sym_nonabel.SIX_ANGLES
SIX_COORD_MASS = {"1": 0, "s3": 0, "r": 1, "s2": 1, "r^2": 2, "s1": 2}              # ＝ sym_nonabel.SIX_MASS
SIX_COORD_TEX = {"1": r"\mathbf 1", "r": "r", "r^2": "r^2", "s1": r"\sigma_1", "s2": r"\sigma_2", "s3": r"\sigma_3"}


def _e3(deg) -> np.ndarray:
    t = np.radians(float(deg))
    return np.array([np.cos(t), np.sin(t), 0.0])


def mirror_lines(radius: float = 2.0, angles=(120, 60, 0), labels=(r"\sigma_1", r"\sigma_2", r"\sigma_3"),
                 color=MIRROR_COLOR, font_size: float = 26, center=ORIGIN) -> VGroup:
    """固定在紙上的鏡線（Fig. 3.1.2；角度＝ sym_nonabel.MIRROR_DEG）：過 center 的虛線＋標籤（放在右半邊那一端）。
    → VGroup(lines, labels)，另有 .lines、.labels。"""
    lines, labs = VGroup(), VGroup()
    c = np.array(center, dtype=float)
    for a, lab in zip(angles, labels):
        d = _e3(a)
        lines.add(DashedLine(c - radius * d, c + radius * d, dash_length=0.12).set_stroke(color, 2))
        if lab:
            if d[0] < -1e-9 or (abs(d[0]) < 1e-9 and d[1] < 0):
                d = -d
            labs.add(Tex(lab, font_size=font_size).set_color(color).move_to(c + (radius + 0.3) * d))
    g = VGroup(lines, labs)
    g.lines, g.labels = lines, labs
    return g


class BobTriangle(VGroup):
    """俯視的三個擺錘（Fig. 3.1.2）：名字 a、b、c，平衡位置在 0°、120°、240°（虛線小圓）；
    state = {名字: (x, y)}（以 radius 為單位；預設＝ sym_nonabel.BOB_START）。
    鏡面／轉動固定在紙上：新的 state 由 SN.config_after(g, tri.state) 算，再 self.play(*tri.anim_to(新 state))。"""

    def __init__(self, radius: float = 1.6, positions=None, names=("a", "b", "c"), bob_radius: float = 0.2,
                 font_size: float = 26, bob_color=BOB_COLOR, show_rest: bool = True, **kw):
        super().__init__(**kw)
        self.radius = radius
        self.names = list(names)
        pos = positions or {"a": (1.0, 0.28), "b": (-0.72, 0.62), "c": (-0.36, -1.02)}
        self.center_dot = Dot(ORIGIN, radius=0.01).set_opacity(0)
        self.ring = DashedVMobject(Circle(radius=radius), num_dashes=40).set_stroke(GREY_D, 1.5)
        self.rest = VGroup(*[DashedVMobject(Circle(radius=bob_radius), num_dashes=10).set_stroke(GREY_C, 1.5)
                             .move_to(radius * _e3(120 * p)) for p in range(3)])
        if not show_rest:
            self.rest.set_opacity(0)
        self.bobs = VGroup(*[VGroup(Circle(radius=bob_radius).set_stroke(WHITE, 2).set_fill(bob_color, 0.95),
                                    Tex(n, font_size=font_size).set_color(BLACK)) for n in self.names])
        self.add(self.center_dot, self.ring, self.rest, self.bobs)
        self.set_state(pos)

    def point(self, p) -> np.ndarray:
        """state 座標 → 畫面座標（隨整組移動）。"""
        return self.center_dot.get_center() + self.radius * np.array([float(p[0]), float(p[1]), 0.0])

    def set_state(self, state):
        self.state = {n: np.array(state[n], dtype=float)[:2] for n in self.names}
        for b, n in zip(self.bobs, self.names):
            b.move_to(self.point(self.state[n]))
        return self

    def anim_to(self, state) -> list:
        """移到新的 state 的動畫列表（直線移動；標籤字母不翻轉）。"""
        self.state = {n: np.array(state[n], dtype=float)[:2] for n in self.names}
        return [b.animate.move_to(self.point(self.state[n])) for b, n in zip(self.bobs, self.names)]


class SpringTriangle(VGroup):
    """Fig. 3.4.1 的六自由度彈簧質量模型（俯視）：三個質點在 0°、120°、240°（半徑 radius），三邊彈簧 k₁，
    每個質點兩根 k₂ 接到外牆（與徑向夾 45° − alpha_deg；alpha = 15° 時沿邊的延長線）。
    .set_displacements(x6)：x6 依 1, r, r², σ₁, σ₂, σ₃ 的正則座標（座標 |g⟩ 的方向 SIX_COORD_DEG），位移 = unit·Σx_g u_g（給 updater）。
    .coord_arrows()：六個座標箭頭＋⟨g|x⟩ 標籤；.displacement_arrows(x6, scale)：每個質點的位移箭頭（Fig. 3.4.2、3.4.4 的靜態圖）。
    .masses、.k1、.k2、.wall、.rest（平衡位置）。"""

    def __init__(self, radius: float = 1.3, wall_radius: float = 2.3, unit: float = 0.5, mass_radius: float = 0.15,
                 alpha_deg: float = 15, show_k2: bool = True, coils: int = 5, amp: float = 0.07, **kw):
        super().__init__(**kw)
        self.radius, self.wall_radius, self.unit, self.coils, self.amp = radius, wall_radius, unit, coils, amp
        self.center_dot = Dot(ORIGIN, radius=0.01).set_opacity(0)
        self.wall = Circle(radius=wall_radius).set_stroke(GREY_B, 4)
        beta = 45 - float(alpha_deg)
        self._anchor_rel = []
        for p in range(3):
            P = radius * _e3(120 * p)
            for sg in (1, -1):
                d = _e3(120 * p + sg * beta)
                b = P @ d
                t = -b + np.sqrt(b * b - (P @ P - wall_radius ** 2))
                self._anchor_rel.append((p, P + t * d))
        self.masses = VGroup(*[Circle(radius=mass_radius).set_stroke(WHITE, 2).set_fill(MASS_COLOR, 0.95) for _ in range(3)])
        self.k1 = VGroup(*[VMobject().set_stroke(SPRING_COLOR, 2.5) for _ in range(3)])
        self.k2 = VGroup(*[VMobject().set_stroke(GREY_B, 2) for _ in range(6)])
        if not show_k2:
            self.k2.set_opacity(0)
        self.add(self.center_dot, self.wall, self.k2, self.k1, self.masses)
        self.set_displacements([0.0] * 6)

    @property
    def rest(self) -> list:
        c = self.center_dot.get_center()
        return [c + self.radius * _e3(120 * p) for p in range(3)]

    def disp(self, x6) -> list:
        """6 個正則座標 → 三個質點的位移向量（已乘 unit）。"""
        out = [np.zeros(3) for _ in range(3)]
        for g, x in zip(SIX_ORDER, x6):
            out[SIX_COORD_MASS[g]] += self.unit * float(np.real(x)) * _e3(SIX_COORD_DEG[g])
        return out

    def set_displacements(self, x6):
        rest, d = self.rest, self.disp(x6)
        c = self.center_dot.get_center()
        pos = [rest[p] + d[p] for p in range(3)]
        r = self.masses[0].get_width() / 2
        for p in range(3):
            self.masses[p].move_to(pos[p])
        for p in range(3):
            a, b = pos[p], pos[(p + 1) % 3]
            u = (b - a) / max(np.linalg.norm(b - a), 1e-9)
            self.k1[p].set_points_as_corners(_zigzag_dir(a + r * u, b - r * u, self.coils, self.amp))
        for k, (p, A) in enumerate(self._anchor_rel):
            a, b = pos[p], c + A
            u = (b - a) / max(np.linalg.norm(b - a), 1e-9)
            self.k2[k].set_points_as_corners(_zigzag_dir(a + r * u, b, self.coils, self.amp * 0.8))
        return self

    def coord_arrows(self, length: float = 0.55, font_size: float = 18, color=KET_COLOR) -> VGroup:
        """六個座標箭頭 |g⟩ 與標籤 ⟨g|x⟩（Fig. 3.4.1）。"""
        out = VGroup()
        rest = self.rest
        for g in SIX_ORDER:
            P, d = rest[SIX_COORD_MASS[g]], _e3(SIX_COORD_DEG[g])
            arr = Arrow(P, P + length * d, buff=0, thickness=3, max_tip_length_to_length_ratio=0.3).set_color(color)
            lab = Tex(r"\langle " + SIX_COORD_TEX[g] + r"|x\rangle", font_size=font_size).set_color(color)
            lab.move_to(P + (length + 0.28) * d)
            out.add(VGroup(arr, lab))
        return out

    def displacement_arrows(self, x6, scale: float = 1.0, color=KET_COLOR) -> VGroup:
        """每個質點的位移箭頭（長度 = scale·unit·|位移|；太短的略過）。"""
        out = VGroup()
        for P, d in zip(self.rest, self.disp(x6)):
            if np.linalg.norm(d) * scale > 0.03:
                out.add(Arrow(P, P + scale * d, buff=0, thickness=4, max_tip_length_to_length_ratio=0.35).set_color(color))
        return out


def sym_table(col_labels, row_labels, rows, cell_w: float = 0.62, row_h: float = 0.42, font_size: float = 24,
              corner: str = "", first_w: float | None = None, label_color=YELLOW) -> VGroup:
    """群乘法表、類代數表、特徵標表、相關表：表頭列＝ col_labels、第一行＝ row_labels（TeX），格子是 TeX 字串／int／sympy（entry_tex）。
    取格子用 table_cell(t, i, j, ncols)（i、j 含表頭，0 起算；ncols = len(col_labels) + 1）。左上角在 ORIGIN。"""
    head = [corner if corner else Dot(radius=0.001).set_opacity(0)] + list(col_labels)
    data = [head] + [[rl] + [entry_tex(x) for x in r] for rl, r in zip(row_labels, rows)]
    widths = [first_w or cell_w] + [cell_w] * len(col_labels)
    t = tex_table(data, widths, row_h=row_h, font_size=font_size)
    for i in range(1, len(data)):
        table_cell(t, i, 0, len(widths)).set_color(label_color)
    return t


def perm_grid(M, cell: float = 0.3, on_color=YELLOW, off_color=GREY_E, row_labels=None, col_labels=None,
              font_size: float = 18, show_values: bool = False) -> VGroup:
    """矩陣的格子圖（正則表示 𝓡(g)、塊對角的 T†𝓡T）：非零格填 on_color。→ VGroup(cells, labels, values)，
    .cells[i][j] 是第 i 列第 j 行的方格（0 起算）。row_labels／col_labels 是 TeX（例如 ⟨h|、|f⟩）。"""
    rows = M.tolist() if hasattr(M, "tolist") else M
    cells, values, labels = VGroup(), VGroup(), VGroup()
    grid = []
    for i, row in enumerate(rows):
        line = []
        for j, v in enumerate(row):
            nz = complex(sp.sympify(v)) != 0
            sq = Square(side_length=cell).set_stroke(GREY_C, 1).set_fill(on_color if nz else off_color, 0.9 if nz else 0.25)
            sq.move_to(np.array([j * cell, -i * cell, 0.0]))
            line.append(sq)
            cells.add(sq)
            if show_values and nz:
                values.add(Tex(entry_tex(v), font_size=font_size).set_color(BLACK).move_to(sq))
        grid.append(line)
    if row_labels:
        for i, s in enumerate(row_labels):
            labels.add(Tex(s, font_size=font_size).set_color(GREY_A).next_to(grid[i][0], LEFT, buff=0.1))
    if col_labels:
        for j, s in enumerate(col_labels):
            labels.add(Tex(s, font_size=font_size).set_color(GREY_A).next_to(grid[0][j], UP, buff=0.1))
    g = VGroup(cells, labels, values)
    g.cells = grid
    return g


class View3D:
    """2D 畫面上的 3D 正交投影（不用 3D 相機：這台 wgpu 渲染器的 Surface／DotCloud 緊接在 VMobject 後面畫會出錯）。
    先繞 z 轉 −theta，再繞 x 轉 −phi：phi = 0 是俯視（z 朝外），phi = 90 是平視（z 朝上）；depth 越大越靠近觀眾。
    .project(v) → 畫面座標（center + scale·(x', y')）、.depth(v)。鏡頭繞轉：改 theta、phi 後對各元件呼叫 .set_view(view)。"""

    def __init__(self, theta: float = -30, phi: float = 65, scale: float = 1.0, center=ORIGIN):
        self.theta, self.phi, self.scale = float(theta), float(phi), float(scale)
        self.center = np.array(center, dtype=float)

    def matrix(self) -> np.ndarray:
        t, p = np.radians(self.theta), np.radians(self.phi)
        Rz = np.array([[np.cos(t), np.sin(t), 0], [-np.sin(t), np.cos(t), 0], [0, 0, 1]])
        Rx = np.array([[1, 0, 0], [0, np.cos(p), np.sin(p)], [0, -np.sin(p), np.cos(p)]])
        return Rx @ Rz

    def project(self, v) -> np.ndarray:
        w = self.matrix() @ np.asarray(v, dtype=float)
        return self.center + self.scale * np.array([w[0], w[1], 0.0])

    def depth(self, v) -> float:
        return float((self.matrix() @ np.asarray(v, dtype=float))[2])


def _runs(points3, view, front_op: float, back_op: float, color, width) -> VGroup:
    """3D 折線投影到畫面，前半（depth ≥ 0）與後半分段上不同的不透明度。"""
    out = VGroup()
    cur, side = [], None
    for v in points3:
        s = view.depth(v) >= 0
        if side is None or s == side:
            cur.append(v)
        else:
            cur.append(v)
            out.add(VMobject().set_points_as_corners([view.project(q) for q in cur]).set_stroke(color, width, front_op if side else back_op))
            cur = [v]
        side = s
    if len(cur) > 1:
        out.add(VMobject().set_points_as_corners([view.project(q) for q in cur]).set_stroke(color, width, front_op if side else back_op))
    return out


class SphereView(VGroup):
    """Hamilton 轉的球（Fig. 3.1.5–3.1.7、3.2.1）：經緯線（前面亮、後面淡）＋可加的大圓弧、轉軸、點，全部用 View3D 投影。
    SphereView(radius=2.0, view=None, lats=(−60, −30, 0, 30, 60), n_lon=6, color=GREY_D)；
    .add_arc(N1, N2, color, width=6, tip=True) → 弧（VGroup(線[, 頭])；N1、N2 是單位向量，用 SN.hamilton_arc／hamilton_product 算）、
    .add_axis(axis, color)、.add_point(v, color, radius)、.project(v)（球面上的點 radius·v 的畫面座標）；
    .set_view(theta, phi) 重新投影全部（鏡頭繞轉：ValueTracker ＋ updater 呼叫它）。弧用 ShowCreation(arc[0]) 畫出來。"""

    def __init__(self, radius: float = 2.0, view: View3D | None = None, lats=(-60, -30, 0, 30, 60), n_lon: int = 6,
                 color=GREY_D, width: float = 1.2, **kw):
        super().__init__(**kw)
        self.radius = radius
        self.view = view or View3D()
        self._lats, self._nlon, self._color, self._width = lats, n_lon, color, width
        self.frame_lines = VGroup()
        self.arcs, self.axes, self.points = VGroup(), VGroup(), VGroup()
        self._items = []
        self.add(self.frame_lines, self.axes, self.arcs, self.points)
        self._draw_frame()

    def project(self, v) -> np.ndarray:
        return self.view.project(self.radius * np.asarray(v, dtype=float))

    def _draw_frame(self):
        ts = np.linspace(0, 2 * np.pi, 60)
        curves = VGroup()
        for la in self._lats:
            c, s = np.cos(np.radians(la)), np.sin(np.radians(la))
            curves.add(*_runs([self.radius * np.array([c * np.cos(t), c * np.sin(t), s]) for t in ts], self.view, 0.7, 0.2,
                              self._color, self._width))
        for k in range(self._nlon):
            ph = np.pi * k / self._nlon
            curves.add(*_runs([self.radius * np.array([np.cos(t) * np.cos(ph), np.cos(t) * np.sin(ph), np.sin(t)]) for t in ts],
                              self.view, 0.7, 0.2, self._color, self._width))
        self.frame_lines.set_submobjects(list(curves))

    def _arc_points(self, N1, N2, n: int = 30) -> list:
        N1, N2 = np.asarray(N1, float), np.asarray(N2, float)
        om = np.arccos(np.clip(N1 @ N2, -1, 1))
        if om < 1e-9:
            return [self.project(N1)] * 2
        return [self.project((np.sin((1 - t) * om) * N1 + np.sin(t * om) * N2) / np.sin(om)) for t in np.linspace(0, 1, n)]

    def add_arc(self, N1, N2, color=YELLOW, width: float = 6, tip: bool = True) -> VGroup:
        pts = self._arc_points(N1, N2)
        line = VMobject().set_points_as_corners(pts).set_stroke(color, width)
        g = VGroup(line)
        if tip:
            g.add(Dot(pts[-1], radius=0.07).set_fill(color, 1))
        self.arcs.add(g)
        self._items.append(("arc", g, (np.asarray(N1, float), np.asarray(N2, float))))
        return g

    def add_axis(self, axis, color=GREY_B, length: float = 1.3) -> DashedLine:
        a = np.asarray(axis, float)
        a = a / np.linalg.norm(a)
        ln = DashedLine(self.project(-length * a), self.project(length * a), dash_length=0.12).set_stroke(color, 2)
        self.axes.add(ln)
        self._items.append(("axis", ln, (a, length)))
        return ln

    def add_point(self, v, color=WHITE, radius: float = 0.06) -> Dot:
        d = Dot(self.project(v), radius=radius).set_fill(color, 1)
        self.points.add(d)
        self._items.append(("point", d, np.asarray(v, float)))
        return d

    def set_view(self, theta: float | None = None, phi: float | None = None):
        if theta is not None:
            self.view.theta = float(theta)
        if phi is not None:
            self.view.phi = float(phi)
        self._draw_frame()
        for kind, m, data in self._items:
            if kind == "arc":
                pts = self._arc_points(*data)
                m[0].set_points_as_corners(pts)
                if len(m) > 1:
                    m[1].move_to(pts[-1])
            elif kind == "axis":
                a, length = data
                m.put_start_and_end_on(self.project(-length * a), self.project(length * a))
            else:
                m.move_to(self.project(data))
        return self


def polar_curve(f, radius: float = 1.5, amp: float = 0.3, center=ORIGIN, n: int = 60, color=YELLOW, width: float = 3) -> VMobject:
    """環上的波（Fig. 3.6.2、3.6.5）：r(φ) = radius + amp·f(φ) 的閉曲線（n 點）。"""
    c = np.array(center, dtype=float)
    pts = [c + (radius + amp * float(f(t))) * np.array([np.cos(t), np.sin(t), 0.0]) for t in np.linspace(0, 2 * np.pi, n)]
    return VMobject().set_points_as_corners(pts).set_stroke(color, width)


class NH3Model(VGroup):
    """球棍 NH₃（Fig. 3.7.1），用 View3D 投影到畫面（不用 3D 相機）：N 在上（N_COLOR）、三個 H 在下（H_COLOR），
    D_deg ＝ N–H 鍵與 H₃ 平面的夾角、鍵長 bond。H_p 在 120p°（H1 在 x 軸上，σ₃ 面內）；座標與 sym_nonabel.nh3_geometry 相同
    （H1 沿 N→H 鍵、H2 切向、H3 ⊥ 鍵向上；N1, N2, N3 = x, y, z），以四個原子的幾何中心為原點。
    .set_displacements(x12)：位移 = unit·Σ x_i u_i（給 updater）；.set_view(theta, phi)（鏡頭繞轉）；前後用 z_index 排。
    .atoms（VGroup：N, H1, H2, H3，每個是 VGroup(圓, 亮點)）、.bonds（N–H）、.hh（H–H 虛線，彈簧 k）。"""

    def __init__(self, D_deg: float = 20, bond: float = 1.6, unit: float = 0.5, n_radius: float = 0.26,
                 h_radius: float = 0.16, view: View3D | None = None, show_hh: bool = True, **kw):
        super().__init__(**kw)
        d = np.radians(float(D_deg))
        c, s = np.cos(d), np.sin(d)

        def rz(p, v):
            a = 2 * np.pi * p / 3
            return np.array([[np.cos(a), -np.sin(a), 0], [np.sin(a), np.cos(a), 0], [0, 0, 1]]) @ v

        pos = [bond * rz(p, np.array([c, 0, 0])) for p in range(3)] + [np.array([0, 0, bond * s])]
        ctr = np.mean(pos, axis=0)
        self._rest = [q - ctr for q in pos]
        us = [np.array([c, 0, -s]), np.array([0.0, 1, 0]), np.array([s, 0, c])]
        self._coords = [(p, rz(p, u)) for u in us for p in range(3)] + [(3, np.eye(3)[i]) for i in range(3)]
        self.unit = unit
        self.view = view or View3D(-30, 70)
        self._x = np.zeros(12)

        def ball(r, col):
            return VGroup(Circle(radius=r).set_stroke(WHITE, 1.5).set_fill(col, 1),
                          Circle(radius=0.35 * r).set_stroke(width=0).set_fill(WHITE, 0.45).shift(0.35 * r * (UP + LEFT)))

        self.atoms = VGroup(ball(n_radius, N_COLOR), *[ball(h_radius, H_COLOR) for _ in range(3)])
        self.bonds = VGroup(*[Line(ORIGIN, RIGHT).set_stroke(ROD_COLOR, 6) for _ in range(3)])
        self.hh = VGroup(*[DashedLine(ORIGIN, RIGHT, dash_length=0.1).set_stroke(GREY_C, 2) for _ in range(3)])
        if not show_hh:
            self.hh.set_opacity(0)
        self.add(self.hh, self.bonds, self.atoms)
        self.set_displacements(self._x)

    def positions3(self, x12) -> list:
        """[H1, H2, H3, N] 的 3D 位置（分子座標）。"""
        disp = [np.zeros(3) for _ in range(4)]
        for x, (p, u) in zip(x12, self._coords):
            disp[p] += self.unit * float(np.real(x)) * u
        return [self._rest[p] + disp[p] for p in range(4)]

    def set_displacements(self, x12):
        self._x = np.array([float(np.real(v)) for v in x12])
        P3 = self.positions3(self._x)
        P = [self.view.project(q) for q in P3]
        dep = [self.view.depth(q) for q in P3]
        order = np.argsort(dep)                                   # 遠的先畫
        rank = {int(a): k for k, a in enumerate(order)}
        self.atoms[0].move_to(P[3]).set_z_index(10 + 2 * rank[3])
        for p in range(3):
            self.atoms[p + 1].move_to(P[p]).set_z_index(10 + 2 * rank[p])
            self.bonds[p].put_start_and_end_on(P[3], P[p]).set_z_index(9 + 2 * min(rank[3], rank[p]))
            self.hh[p].put_start_and_end_on(P[p], P[(p + 1) % 3]).set_z_index(1)
        return self

    def set_view(self, theta: float | None = None, phi: float | None = None):
        if theta is not None:
            self.view.theta = float(theta)
        if phi is not None:
            self.view.phi = float(phi)
        return self.set_displacements(self._x)


# ---------------------------------------------------------------- 第 4 章（八面體群、誘導表示、XF6）
DIAG_COLORS = {1: "#FF6B6B", 2: "#F2C14E", 3: "#7FD1AE", 4: "#5B8DEF"}     # 四條體對角線（編號＝ sym_oct.DIAGONALS）
POS_COLOR = "#5B8DEF"                                                        # 振幅 +（實部 > 0）
NEG_COLOR = "#E07A5F"                                                        # 振幅 −
IPOS_COLOR = "#7FD1AE"                                                       # 振幅 +i
INEG_COLOR = "#C792EA"                                                       # 振幅 −i
X_COLOR = "#F2C14E"                                                          # XF6 的中心原子
F_COLOR = "#9AD1D4"                                                          # XF6 的 F
OCT_COLORS = {**LABEL_COLORS, "T_1": "#FF8C69", "T_2": "#C792EA", "A": BLUE, "T": "#FF8C69"}
OCT_SITES = [(0, 0, 1), (0, 0, -1), (1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0)]          # ＝ sym_oct.SITES（|1⟩…|6⟩）
OCT_DIAGONALS = {1: (1, -1, 1), 2: (-1, 1, 1), 3: (1, 1, -1), 4: (1, 1, 1)}                 # ＝ sym_oct.DIAGONALS
XF6_COORDS = [   # ＝ sym_oct.XF6_COORDS_LIST：(名稱, 原子（0–5 = OCT_SITES 上的 F、6 = X）, 位移方向)
    ("A", 0, (0, 0, 1)), ("R1^2A", 1, (0, 0, -1)), ("R1A", 5, (0, -1, 0)), ("R2A", 2, (1, 0, 0)), ("R1^3A", 4, (0, 1, 0)),
    ("R2^3A", 3, (-1, 0, 0)), ("B", 0, (1, 0, 0)), ("r1B", 2, (0, 1, 0)), ("r2B", 3, (0, 1, 0)), ("r1^2B", 4, (0, 0, 1)),
    ("r4^2B", 5, (0, 0, 1)), ("R2^2B", 1, (-1, 0, 0)), ("R1B", 5, (1, 0, 0)), ("R2B", 2, (0, 0, -1)), ("R3B", 0, (0, 1, 0)),
    ("R1^3B", 4, (1, 0, 0)), ("R2^3B", 3, (0, 0, 1)), ("i3B", 1, (0, 1, 0)), ("C1", 6, (1, 0, 0)), ("C2", 6, (0, 1, 0)),
    ("C3", 6, (0, 0, 1))]


def irrep_color(label: str):
    """不可約表示標記 → 顏色（"T_{1u}"、"E_g"、"A_1" 都可以；g/u 不分）。"""
    s = label.replace("{", "").replace("}", "").replace("'", "")
    if len(s) > 1 and s[-1] in "gu" and s[:-1] not in ("", "_"):
        s = s[:-1]
    s = s.rstrip("_")
    if len(s) == 2 and s[1].isdigit():
        s = s[0] + "_" + s[1]
    return OCT_COLORS.get(s, WHITE)


def rot_matrix(axis, deg) -> np.ndarray:
    """繞 axis 逆時針轉 deg 度的 3×3 矩陣（Rodrigues，浮點；動畫用 rot_matrix(axis, t·deg)）。"""
    n = np.asarray(axis, float)
    n = n / np.linalg.norm(n)
    t = np.radians(float(deg))
    K = np.array([[0, -n[2], n[1]], [n[2], 0, -n[0]], [-n[1], n[0], 0]])
    return np.eye(3) + np.sin(t) * K + (1 - np.cos(t)) * K @ K


def _depth_op(view, v, front: float = 1.0, back: float = 0.3) -> float:
    return front if view.depth(v) >= -1e-9 else back


class CubeView(VGroup):
    """投影的立方體（Fig. 4.1.1–4.1.2、4.2.1–4.2.2）：12 條稜、四條上色的體對角線（DIAG_COLORS，端點標 1–4，編號＝ sym_oct.DIAGONALS）、
    可選內接八面體（頂點在面心）。CubeView(size=1.3, view=None, show_cube=True, show_octa=False, show_diagonals=True,
    diag_labels=True, font_size=22)。物體的轉動：.set_rotation(Rm)（3×3 浮點；動畫用 ValueTracker t 與
    rot_matrix(axis, t·deg) @ R_前一個）；鏡頭：.set_view(theta, phi)。.add_axis(axis, color, length, label_tex)＝ 固定在實驗室的虛線轉軸；
    .add_arrow(v, color, label_tex, body=False)＝ 從中心出發的向量（body=True 跟著物體轉）。.edges、.diagonals、.diag_labels、.octa、.R。
    前後用不透明度區分（後面的稜 0.3）。位置用 View3D(center=…)。"""

    def __init__(self, size: float = 1.3, view: View3D | None = None, show_cube: bool = True, show_octa: bool = False,
                 show_diagonals: bool = True, diag_labels: bool = True, font_size: float = 22, **kw):
        super().__init__(**kw)
        self.size = size
        self.view = view or View3D(-25, 70)
        self.R = np.eye(3)
        verts = [np.array(v, float) for v in itertools.product((1, -1), repeat=3)]
        self._edges = [(a, b) for a in verts for b in verts if np.sum(np.abs(a - b)) == 2 and tuple(a) < tuple(b)]
        octv = [np.array(v, float) for v in OCT_SITES]
        self._octa = [(a, b) for i, a in enumerate(octv) for b in octv[i + 1:] if abs(a @ b) < 1e-9]
        self.edges = VGroup(*[Line(ORIGIN, RIGHT).set_stroke(GREY_A, 2.5) for _ in self._edges])
        self.octa = VGroup(*[Line(ORIGIN, RIGHT).set_stroke(TEAL_B, 2) for _ in self._octa])
        self.diagonals = VGroup(*[Line(ORIGIN, RIGHT).set_stroke(DIAG_COLORS[k], 4) for k in range(1, 5)])
        self.diag_labels = VGroup(*[Tex(str(k), font_size=font_size).set_color(DIAG_COLORS[k]) for k in range(1, 5)])
        if not show_cube:
            self.edges.set_opacity(0)
        if not show_octa:
            self.octa.set_opacity(0)
        if not show_diagonals:
            self.diagonals.set_opacity(0)
        if not (diag_labels and show_diagonals):
            self.diag_labels.set_opacity(0)
        self._show = dict(cube=show_cube, octa=show_octa, diag=show_diagonals, lab=diag_labels and show_diagonals)
        self.axes_lab, self.arrows = VGroup(), VGroup()
        self._items = []
        self.add(self.edges, self.octa, self.diagonals, self.diag_labels, self.axes_lab, self.arrows)
        self.redraw()

    def _w(self, v) -> np.ndarray:
        return self.size * (self.R @ np.asarray(v, float))

    def redraw(self):
        pr = self.view.project
        for ln, (a, b) in zip(self.edges, self._edges):
            ln.put_start_and_end_on(pr(self._w(a)), pr(self._w(b)))
            if self._show["cube"]:
                ln.set_stroke(opacity=_depth_op(self.view, self._w((a + b) / 2)))
        for ln, (a, b) in zip(self.octa, self._octa):
            ln.put_start_and_end_on(pr(self._w(a)), pr(self._w(b)))
            if self._show["octa"]:
                ln.set_stroke(opacity=_depth_op(self.view, self._w((a + b) / 2), 0.9, 0.25))
        for k in range(1, 5):
            d = np.array(OCT_DIAGONALS[k], float)
            self.diagonals[k - 1].put_start_and_end_on(pr(self._w(-d)), pr(self._w(d)))
            if self._show["lab"]:
                self.diag_labels[k - 1].move_to(pr(self._w(1.22 * d)))
        for kind, m, data in self._items:
            v, body = data
            w = self._w(v) if body else self.size * np.asarray(v, float)
            if kind == "arrow":
                m[0].put_start_and_end_on(pr(np.zeros(3)), pr(w))
                if len(m) > 1:
                    m[1].move_to(pr(1.18 * w))
        return self

    def set_rotation(self, Rm):
        self.R = np.asarray(Rm, float)
        return self.redraw()

    def set_view(self, theta: float | None = None, phi: float | None = None):
        if theta is not None:
            self.view.theta = float(theta)
        if phi is not None:
            self.view.phi = float(phi)
        for kind, m, data in self._items:
            if kind == "axis":
                a, length = data
                m[0].put_start_and_end_on(self.view.project(-length * self.size * a), self.view.project(length * self.size * a))
                if len(m) > 1:
                    m[1].move_to(self.view.project(1.12 * length * self.size * a))
        return self.redraw()

    def add_axis(self, axis, color=GREY_B, length: float = 1.6, label_tex: str | None = None, font_size: float = 22) -> VGroup:
        a = np.asarray(axis, float)
        a = a / np.linalg.norm(a)
        g = VGroup(DashedLine(self.view.project(-length * self.size * a), self.view.project(length * self.size * a),
                              dash_length=0.1).set_stroke(color, 2.5))
        if label_tex:
            g.add(Tex(label_tex, font_size=font_size).set_color(color).move_to(self.view.project(1.12 * length * self.size * a)))
        self.axes_lab.add(g)
        self._items.append(("axis", g, (a, length)))
        return g

    def add_arrow(self, v, color=KET_COLOR, label_tex: str | None = None, body: bool = False, font_size: float = 22) -> VGroup:
        v = np.asarray(v, float)
        g = VGroup(Arrow(ORIGIN, RIGHT, buff=0, thickness=4).set_color(color))
        if label_tex:
            g.add(Tex(label_tex, font_size=font_size).set_color(color))
        self.arrows.add(g)
        self._items.append(("arrow", g, (v, body)))
        self.redraw()
        return g


def _amp_color(z: complex):
    ph = np.angle(z)
    q = int(np.round(ph / (np.pi / 2))) % 4
    return [POS_COLOR, IPOS_COLOR, NEG_COLOR, INEG_COLOR][q]


class OctaSites(VGroup):
    """八面體的六個頂點（Fig. 4.3.1–4.3.2、4.3.4–4.3.7）：|1⟩…|6⟩ 在 +z, −z, +x, −x, +y, −y（＝ sym_oct.SITES），用 View3D 投影。
    OctaSites(radius=1.6, view=None, ball=0.2, labels=True, show_edges=True, font_size=22, label_texts=None)。
    .set_amplitudes(v, min_frac=0.25, value_texs=None)：本徵向量的六個分量 → 球的大小 ∝ |v_i|、顏色依相位
    （+ POS_COLOR、− NEG_COLOR、+i IPOS_COLOR、−i INEG_COLOR；0 的分量畫成空心小圈）；value_texs 給六個 TeX 字串就標在球旁。
    .set_view(theta, phi)、.point(i)、.balls、.edges、.labels、.values。位置用 View3D(center=…)。"""

    def __init__(self, radius: float = 1.6, view: View3D | None = None, ball: float = 0.2, labels: bool = True,
                 show_edges: bool = True, font_size: float = 22, label_texts=None, **kw):
        super().__init__(**kw)
        self.radius, self.ball, self.font_size = radius, ball, font_size
        self.view = view or View3D(-25, 70)
        self.sites = [radius * np.array(s, float) for s in OCT_SITES]
        self._pairs = [(i, j) for i in range(6) for j in range(i + 1, 6) if abs(self.sites[i] @ self.sites[j]) < 1e-9]
        self.edges = VGroup(*[Line(ORIGIN, RIGHT).set_stroke(GREY_C, 1.5) for _ in self._pairs])
        if not show_edges:
            self.edges.set_opacity(0)
        self._show_edges = show_edges
        self.balls = VGroup(*[Circle(radius=ball).set_stroke(WHITE, 1.5).set_fill(GREY_B, 0.9) for _ in range(6)])
        texts = label_texts or [rf"|{k}\rangle" for k in range(1, 7)]
        self.labels = VGroup(*[Tex(t, font_size=font_size).set_color(GREY_A) for t in texts])
        if not labels:
            self.labels.set_opacity(0)
        self.values = VGroup()
        self._amps = None
        self.add(self.edges, self.balls, self.labels, self.values)
        self.redraw()

    def point(self, i: int) -> np.ndarray:
        return self.view.project(self.sites[i])

    def redraw(self):
        for ln, (i, j) in zip(self.edges, self._pairs):
            ln.put_start_and_end_on(self.point(i), self.point(j))
            if self._show_edges:
                ln.set_stroke(opacity=_depth_op(self.view, (self.sites[i] + self.sites[j]) / 2, 0.6, 0.2))
        dep = [self.view.depth(s) for s in self.sites]
        rank = {int(a): k for k, a in enumerate(np.argsort(dep))}
        for i in range(6):
            self.balls[i].move_to(self.point(i)).set_z_index(5 + rank[i])
            off = self.view.project(1.32 * self.sites[i]) - self.point(i)
            d = off / np.linalg.norm(off) if np.linalg.norm(off) > 0.05 else UP
            self.labels[i].move_to(self.point(i) + d * (1.8 * self.ball + 0.18))
        for i, m in enumerate(self.values):
            m.next_to(self.balls[i], RIGHT, buff=0.06)
        return self

    def set_amplitudes(self, v, min_frac: float = 0.25, value_texs=None):
        a = np.asarray(v, dtype=complex)
        mx = max(np.abs(a).max(), 1e-12)
        for i in range(6):
            b = self.balls[i]
            if abs(a[i]) < 1e-9:
                b.set_width(2 * self.ball * min_frac * 0.8).set_fill(opacity=0).set_stroke(GREY_B, 1.5)
            else:
                r = self.ball * (min_frac + (1 - min_frac) * abs(a[i]) / mx) * 1.6
                b.set_width(2 * r).set_fill(_amp_color(a[i]), 0.95).set_stroke(WHITE, 1.5)
        self._amps = a
        if value_texs is not None:
            self.values.set_submobjects([Tex(t, font_size=self.font_size - 2).set_color(GREY_A) for t in value_texs])
        return self.redraw()

    def set_view(self, theta: float | None = None, phi: float | None = None):
        if theta is not None:
            self.view.theta = float(theta)
        if phi is not None:
            self.view.phi = float(phi)
        return self.redraw()


def level_diagram(columns, links=(), col_w: float = 1.0, gap: float = 1.4, unit: float = 0.5, font_size: float = 22,
                  deg_gap: float = 0.07, colors=None, headers=None, header_size: float = 20, stroke_width: float = 4) -> VGroup:
    """能階分裂／相關圖（Fig. 4.2.3、4.3.2、4.3.5、4.4.3）：columns ＝ 每一欄的 [(能量, 標記 TeX, 簡併度), …]；簡併度 d 畫成 d 條間隔 deg_gap
    的短線（書上的多重線），標記放右邊、顏色用 irrep_color（或 colors={標記: 顏色}）；links ＝ [((欄, 第幾個), (欄, 第幾個)), …] 畫虛線相連；
    headers ＝ 每欄上方的中文字串（zh）。左下角大約在 ORIGIN（第 0 欄、能量 0）。→ VGroup，.levels[c][i]（VGroup 的線）、
    .labels[c][i]、.links、.headers。標記一律 TeX，不放中文。"""
    levels, labels, linkg, heads = [], [], VGroup(), VGroup()
    allg = VGroup()
    for c, col in enumerate(columns):
        x0 = c * (col_w + gap)
        lv_c, lb_c = [], []
        for (E, lab, d) in col:
            y = unit * float(E)
            col_ = (colors or {}).get(lab, irrep_color(lab))
            g = VGroup(*[Line(np.array([x0, y + (k - (d - 1) / 2) * deg_gap, 0]), np.array([x0 + col_w, y + (k - (d - 1) / 2) * deg_gap, 0]))
                         .set_stroke(col_, stroke_width) for k in range(int(d))])
            t = Tex(lab, font_size=font_size).set_color(col_).next_to(g, RIGHT, buff=0.08)
            lv_c.append(g)
            lb_c.append(t)
            allg.add(g, t)
        levels.append(lv_c)
        labels.append(lb_c)
    for (c1, i1), (c2, i2) in links:
        a = levels[c1][i1].get_right() + 0.02 * RIGHT
        if labels[c1][i1].get_width() > 0:
            a = labels[c1][i1].get_right() + 0.05 * RIGHT
        b = levels[c2][i2].get_left() - 0.02 * RIGHT
        linkg.add(DashedLine(a, b, dash_length=0.06).set_stroke(GREY_B, 1.5))
    if headers:
        top = max(m.get_top()[1] for m in allg) + 0.35
        for c, h in enumerate(headers):
            heads.add(zh(h, font_size=header_size).move_to(np.array([c * (col_w + gap) + col_w / 2, top, 0])))
    out = VGroup(allg, linkg, heads)
    out.levels, out.labels, out.links, out.headers = levels, labels, linkg, heads
    return out


class XF6Model(VGroup):
    """球棍 XF6（Fig. 4.4.1–4.4.4、4.4.7），用 View3D 投影：中心 X（X_COLOR）、六個 F（F_COLOR，原子 0–5 在 OCT_SITES·bond）、
    X–F 鍵（彈簧 j）、F–F 虛線（彈簧 k，12 根）。座標＝ sym_oct.XF6_NAMES 的 21 個（coords 預設 XF6_COORDS ＝ sym_oct.XF6_COORDS_LIST）。
    XF6Model(bond=1.5, unit=0.6, view=None, x_radius=0.26, f_radius=0.17, show_ff=True, coords=None)；
    .set_displacements(x21)（位移 ＝ unit·Σ x_c 方向_c；給 updater）、.arrows(x21, scale=1.0, color=KET_COLOR)（靜態位移箭頭 VGroup）、
    .set_view(theta, phi)、.positions3(x21)、.atoms（VGroup：X, F0…F5，每個是 VGroup(圓, 亮點)）、.bonds、.ff；前後用 z_index 排。"""

    def __init__(self, bond: float = 1.5, unit: float = 0.6, view: View3D | None = None, x_radius: float = 0.26,
                 f_radius: float = 0.17, show_ff: bool = True, coords=None, **kw):
        super().__init__(**kw)
        self.coords = coords or XF6_COORDS
        self.unit = unit
        self.view = view or View3D(-25, 70)
        self._rest = [bond * np.array(s, float) for s in OCT_SITES] + [np.zeros(3)]
        self._x = np.zeros(len(self.coords))

        def ball(r, col):
            return VGroup(Circle(radius=r).set_stroke(WHITE, 1.5).set_fill(col, 1),
                          Circle(radius=0.35 * r).set_stroke(width=0).set_fill(WHITE, 0.45).shift(0.35 * r * (UP + LEFT)))

        self.atoms = VGroup(ball(x_radius, X_COLOR), *[ball(f_radius, F_COLOR) for _ in range(6)])
        self.bonds = VGroup(*[Line(ORIGIN, RIGHT).set_stroke(ROD_COLOR, 5) for _ in range(6)])
        self._ffp = [(i, j) for i in range(6) for j in range(i + 1, 6) if abs(self._rest[i] @ self._rest[j]) < 1e-9]
        self.ff = VGroup(*[DashedLine(ORIGIN, RIGHT, dash_length=0.08).set_stroke(GREY_C, 1.5) for _ in self._ffp])
        if not show_ff:
            self.ff.set_opacity(0)
        self.add(self.ff, self.bonds, self.atoms)
        self.set_displacements(self._x)

    def positions3(self, x21) -> list:
        """[F0, …, F5, X] 的 3D 位置。"""
        disp = [np.zeros(3) for _ in range(7)]
        for x, (_, a, v) in zip(x21, self.coords):
            disp[a] += self.unit * float(np.real(x)) * np.asarray(v, float)
        return [self._rest[a] + disp[a] for a in range(7)]

    def set_displacements(self, x21):
        self._x = np.array([float(np.real(v)) for v in x21])
        P3 = self.positions3(self._x)
        P = [self.view.project(q) for q in P3]
        dep = [self.view.depth(q) for q in P3]
        rank = {int(a): k for k, a in enumerate(np.argsort(dep))}
        self.atoms[0].move_to(P[6]).set_z_index(20 + 2 * rank[6])
        for a in range(6):
            self.atoms[a + 1].move_to(P[a]).set_z_index(20 + 2 * rank[a])
            self.bonds[a].put_start_and_end_on(P[6], P[a]).set_z_index(19 + 2 * min(rank[6], rank[a]))
        for ln, (i, j) in zip(self.ff, self._ffp):
            ln.put_start_and_end_on(P[i], P[j]).set_z_index(1)
        return self

    def arrows(self, x21, scale: float = 1.0, color=KET_COLOR, min_len: float = 0.04) -> VGroup:
        out = VGroup()
        P0 = self.positions3(np.zeros(len(self.coords)))
        P1 = self.positions3(scale * np.asarray(x21, float))
        for a in range(7):
            s, e = self.view.project(P0[a]), self.view.project(P1[a])
            if np.linalg.norm(e - s) > min_len:
                out.add(Arrow(s, e, buff=0, thickness=4, max_tip_length_to_length_ratio=0.35).set_color(color).set_z_index(40))
        return out

    def set_view(self, theta: float | None = None, phi: float | None = None):
        if theta is not None:
            self.view.theta = float(theta)
        if phi is not None:
            self.view.phi = float(phi)
        return self.set_displacements(self._x)


# ======================================================================== 第 5 章新增
BODY_COLORS = {0: "#FF6B6B", 1: "#7FD1AE", 2: "#5B8DEF"}      # 本體三腳架 x̄、ȳ、z̄
LAB_COLOR = GREY_B                                               # 實驗室軸
CRANK_COLOR = "#F2C14E"                                          # 曲柄（正在轉的軸）
AMP_POS_COLOR, AMP_NEG_COLOR = POS_COLOR, NEG_COLOR             # 振幅長條的正、負
PROB_COLOR = "#8FD3FF"                                           # 機率長條


def _circle3(normal, radius: float, n: int = 48) -> list:
    """3D 中法向量 normal、半徑 radius、圓心在原點的圓上 n+1 點（閉合）。"""
    nv = np.asarray(normal, float)
    nv = nv / np.linalg.norm(nv)
    a = np.cross(nv, [1.0, 0, 0] if abs(nv[0]) < 0.9 else [0, 1.0, 0])
    a /= np.linalg.norm(a)
    b = np.cross(nv, a)
    return [radius * (np.cos(t) * a + np.sin(t) * b) for t in np.linspace(0, 2 * np.pi, n + 1)]


class TriadBall(VGroup):
    """投影的球＋本體三腳架（Fig. 5.1.1、5.3.1–5.3.5、5.5.10–5.5.11）：外輪廓（正交投影下永遠是圓）、跟著物體轉的赤道與本初子午線
    （前亮後淡）、本體軸 x̄ ȳ z̄ 的箭頭（BODY_COLORS）、固定在實驗室的虛線軸 x y z（LAB_COLOR，show_lab）。
    TriadBall(radius=1.4, view=None, show_lab=True, lab_labels=("x", "y", "z"), body_labels=(r"\bar x", r"\bar y", r"\bar z"),
              axis_len=1.35, font_size=24)；
    .set_rotation(Rm)（3×3，可以是非正常的 −R：三腳架變成左手，赤道與子午線照樣畫）、.set_view(theta, phi)、
    .add_crank(axis, color=CRANK_COLOR, label_tex=None, length=1.6)（實驗室固定的轉軸：粗線＋軸端一小段圓弧＋點；回傳 VGroup，可 FadeOut）、
    .tip(k)（本體第 k 軸箭頭尖的畫面座標）、.trail(k, Rs)（一串轉動下第 k 軸尖的軌跡，VMobject，≤ 60 個 R）、
    .project(v)（球面上 radius·v 的畫面座標）。.outline、.equator、.meridian、.body、.lab、.R。位置用 View3D(center=…)，不要事後 move_to。"""

    def __init__(self, radius: float = 1.4, view: View3D | None = None, show_lab: bool = True, lab_labels=("x", "y", "z"),
                 body_labels=(r"\bar x", r"\bar y", r"\bar z"), axis_len: float = 1.35, font_size: float = 24, **kw):
        super().__init__(**kw)
        self.radius, self.axis_len, self.font_size = radius, axis_len, font_size
        self.view = view or View3D(-30, 70)
        self.R = np.eye(3)
        self.outline = Circle(radius=radius * self.view.scale).set_stroke(GREY_C, 2).move_to(self.view.center)
        self.equator, self.meridian = VGroup(), VGroup()
        self.lab = VGroup()
        if show_lab:
            for k, e in enumerate(np.eye(3)):
                ln = DashedLine(self.view.project(-1.6 * radius * e), self.view.project(1.6 * radius * e), dash_length=0.1)
                ln.set_stroke(LAB_COLOR, 1.8)
                lab = Tex(lab_labels[k], font_size=font_size).set_color(LAB_COLOR).move_to(self.view.project(1.78 * radius * e))
                self.lab.add(VGroup(ln, lab))
        self.body = VGroup()
        for k in range(3):
            g = VGroup(Arrow(ORIGIN, RIGHT, buff=0, thickness=4).set_color(BODY_COLORS[k]))
            if body_labels:
                g.add(Tex(body_labels[k], font_size=font_size).set_color(BODY_COLORS[k]))
            self.body.add(g)
        self.add(self.lab, self.outline, self.equator, self.meridian, self.body)
        self._cranks = []
        self.redraw()

    def project(self, v) -> np.ndarray:
        return self.view.project(self.radius * np.asarray(v, float))

    def redraw(self):
        self.outline.move_to(self.view.center)
        R = self.R
        eq = [R @ p for p in _circle3((0, 0, 1), self.radius)]
        mer = [R @ p for p in _circle3((0, 1, 0), self.radius)]
        self.equator.set_submobjects(list(_runs(eq, self.view, 0.9, 0.25, GREY_A, 2)))
        self.meridian.set_submobjects(list(_runs(mer, self.view, 0.9, 0.25, GREY_A, 2)))
        for k in range(3):
            w = self.axis_len * self.radius * (R @ np.eye(3)[k])
            self.body[k][0].put_start_and_end_on(self.view.project(np.zeros(3)), self.view.project(w))
            self.body[k][0].set_opacity(1.0 if self.view.depth(w) >= -1e-9 else 0.55)
            if len(self.body[k]) > 1:
                self.body[k][1].move_to(self.view.project(1.18 * w))
        return self

    def set_rotation(self, Rm):
        self.R = np.asarray(Rm, float)
        return self.redraw()

    def set_view(self, theta: float | None = None, phi: float | None = None):
        if theta is not None:
            self.view.theta = float(theta)
        if phi is not None:
            self.view.phi = float(phi)
        for k, g in enumerate(self.lab):
            e = np.eye(3)[k]
            g[0].put_start_and_end_on(self.view.project(-1.6 * self.radius * e), self.view.project(1.6 * self.radius * e))
            g[1].move_to(self.view.project(1.78 * self.radius * e))
        for g, (a, length) in self._cranks:
            self._place_crank(g, a, length)
        return self.redraw()

    def _place_crank(self, g, a, length):
        g[0].put_start_and_end_on(self.view.project(-length * self.radius * a), self.view.project(length * self.radius * a))
        ring = [length * 0.92 * self.radius * a + p for p in _circle3(a, 0.18 * self.radius, 24)[:19]]
        g[1].set_points_as_corners([self.view.project(q) for q in ring])
        g[2].move_to(self.view.project(ring[-1]))
        if len(g) > 3:
            g[3].move_to(self.view.project(1.12 * length * self.radius * a))

    def add_crank(self, axis, color=CRANK_COLOR, label_tex: str | None = None, length: float = 1.6) -> VGroup:
        a = np.asarray(axis, float)
        a = a / np.linalg.norm(a)
        g = VGroup(Line(ORIGIN, RIGHT).set_stroke(color, 6), VMobject().set_stroke(color, 3), Dot(radius=0.06).set_fill(color, 1))
        if label_tex:
            g.add(Tex(label_tex, font_size=self.font_size).set_color(color))
        self._place_crank(g, a, length)
        self._cranks.append((g, (a, length)))
        self.add(g)
        return g

    def tip(self, k: int) -> np.ndarray:
        return self.view.project(self.axis_len * self.radius * (self.R @ np.eye(3)[k]))

    def trail(self, k: int, Rs, color=None, width: float = 3) -> VMobject:
        pts = [self.view.project(self.axis_len * self.radius * (np.asarray(R, float) @ np.eye(3)[k])) for R in Rs]
        return VMobject().set_points_as_corners(pts).set_stroke(color or BODY_COLORS[k], width)


class VectorCone(VGroup):
    """角動量向量錐（Fig. 5.4.4、5.5.3、5.5.6）：頂點在原點、軸沿 R_y(β)ẑ、斜高 √(j(j+1))·unit、軸上的高 m·unit（cos θ = m/√(j(j+1))）。
    VectorCone(j, m, unit=0.5, view=None, beta=0, color=TEAL, show_axis=True, fill=0.12)；畫：底圓（投影成橢圓）、兩條輪廓母線、
    中心軸（虛線）、軸上高 m 的點。.set_beta(β)（度；給 updater）、.rim、.sides、.axis、.foot、.view。j、m 可為半整數（float）。
    位置用 View3D(center=…)。"""

    def __init__(self, j, m, unit: float = 0.5, view: View3D | None = None, beta: float = 0.0, color=TEAL, show_axis: bool = True,
                 fill: float = 0.12, **kw):
        super().__init__(**kw)
        self.j, self.m, self.unit, self.color, self.fill = float(j), float(m), unit, color, fill
        self.view = view or View3D(-20, 75)
        self.rim = VMobject().set_stroke(color, 2.5).set_fill(color, fill)
        self.sides = VGroup(Line(ORIGIN, RIGHT).set_stroke(color, 2.5), Line(ORIGIN, RIGHT).set_stroke(color, 2.5))
        self.axis = DashedLine(ORIGIN, RIGHT, dash_length=0.08).set_stroke(GREY_B, 1.8)
        if not show_axis:
            self.axis.set_opacity(0)
        self.foot = Dot(radius=0.045).set_fill(color, 1)
        self.add(self.rim, self.sides, self.axis, self.foot)
        self.set_beta(beta)

    def set_beta(self, beta: float):
        self.beta = float(beta)
        L = np.sqrt(self.j * (self.j + 1)) * self.unit
        h = self.m * self.unit
        r = np.sqrt(max(L * L - h * h, 0.0))
        c, s = np.cos(np.radians(self.beta)), np.sin(np.radians(self.beta))
        Ry = np.array([[c, 0, s], [0, 1, 0], [-s, 0, c]])
        ax = Ry @ np.array([0, 0, 1.0])
        rim3 = [Ry @ (np.array([0, 0, h]) + p) for p in _circle3((0, 0, 1), r, 48)]
        pr = self.view.project
        pts = [pr(q) for q in rim3]
        self.rim.set_points_as_corners(pts)
        o = pr(np.zeros(3))
        foot = pr(Ry @ np.array([0, 0, h]))
        if r > 1e-9:
            u = self.view.matrix() @ ax
            perp = np.array([-u[1], u[0], 0.0])
            perp = perp / np.linalg.norm(perp) if np.linalg.norm(perp) > 1e-9 else np.array([1.0, 0, 0])
            proj = [(np.asarray(p) - foot)[:2] @ perp[:2] for p in pts]
            self.sides[0].put_start_and_end_on(o, pts[int(np.argmax(proj))])
            self.sides[1].put_start_and_end_on(o, pts[int(np.argmin(proj))])
        else:
            self.sides[0].put_start_and_end_on(o, foot)
            self.sides[1].put_start_and_end_on(o, foot)
        sgn = 1.0 if h >= 0 else -1.0
        self.axis.put_start_and_end_on(o, pr(sgn * 1.15 * L * ax))
        self.foot.move_to(foot)
        return self


def cone_family(j, unit: float = 0.5, view: View3D | None = None, colors=None) -> VGroup:
    """Fig. 5.4.4：m = j…−j 的一整族圓錐（同一個頂點）＋ z 軸虛線。colors：{m: 顏色}；預設由暖到冷。→ VGroup，.cones[k] 對應 m = j − k。"""
    view = view or View3D(-20, 75)
    ms = [j - k for k in range(int(round(2 * j)) + 1)]
    pal = color_gradient([RED_B, YELLOW, TEAL, BLUE_B], max(len(ms), 2))
    cones = [VectorCone(j, m, unit=unit, view=view, color=(colors or {}).get(m, pal[k]), show_axis=False, fill=0.06)
             for k, m in enumerate(ms)]
    L = 1.25 * np.sqrt(j * (j + 1)) * unit
    g = VGroup(*cones, DashedLine(view.project((0, 0, -L)), view.project((0, 0, L)), dash_length=0.08).set_stroke(GREY_B, 1.8))
    g.cones = cones
    return g


class ProbBars(VGroup):
    """機率／振幅長條（Fig. 5.5.2–5.5.5）：n 根長條、底線、下面的標記。ProbBars(values, labels=None, width=0.32, gap=0.1,
    height=2.0, signed=False, color=PROB_COLOR, font_size=20, label_every=1)；values 在 [0, 1]（signed=True 時 [−1, 1]，
    正 AMP_POS_COLOR、負 AMP_NEG_COLOR，從中線往上下長）。.set_values(v)（給 updater）、.bars、.labels、.baseline、.top(i)（第 i 根頂端）。
    左下角（signed 時中線左端）在 ORIGIN；整組可以 .shift()（set_values 會照目前的位置重畫）。"""

    def __init__(self, values, labels=None, width: float = 0.32, gap: float = 0.1, height: float = 2.0, signed: bool = False,
                 color=PROB_COLOR, font_size: float = 20, label_every: int = 1, **kw):
        super().__init__(**kw)
        self.n, self.w, self.gap, self.h, self.signed, self.color = len(values), width, gap, height, signed, color
        self.bars = VGroup(*[Rectangle(width=width, height=0.01).set_stroke(WHITE, 0.8).set_fill(color, 0.85) for _ in range(self.n)])
        total = self.n * width + (self.n - 1) * gap
        self.baseline = Line(np.array([-0.1, 0, 0]), np.array([total + 0.1, 0, 0])).set_stroke(GREY_B, 1.5)
        self.labels = VGroup()
        if labels is not None:
            for i, s in enumerate(labels):
                t = Tex(s, font_size=font_size).set_color(GREY_A)
                t.move_to(np.array([self._x(i), -0.25 - (self.h if signed else 0), 0]))
                if i % label_every:
                    t.set_opacity(0)
                self.labels.add(t)
        self.add(self.baseline, self.bars, self.labels)
        self.set_values(values)

    def _x(self, i: int) -> float:
        return i * (self.w + self.gap) + self.w / 2

    def set_values(self, values):
        o = self.baseline.get_left() + np.array([0.1, 0, 0])                      # 目前的原點（跟著 shift）
        for i, (b, v) in enumerate(zip(self.bars, values)):
            v = float(np.real(v))
            hh = max(abs(v) * self.h, 0.01)
            b.stretch_to_fit_height(hh)
            if self.signed:
                b.move_to(o + np.array([self._x(i), (np.sign(v) * hh / 2) if v else 0.0, 0]))
                b.set_fill(AMP_POS_COLOR if v >= 0 else AMP_NEG_COLOR, 0.85)
            else:
                b.move_to(o + np.array([self._x(i), hh / 2, 0]))
        return self

    def top(self, i: int) -> np.ndarray:
        return self.bars[i].get_top()


def ladder_column(levels, width: float = 1.0, unit: float = 0.5, font_size: float = 20, stroke_width: float = 3,
                  left_buff: float = 0.1, right_buff: float = 0.1) -> VGroup:
    """轉子／殼層能階梯（Fig. 5.5.7、5.6.7、5.8.6）：levels = [(能量, 左邊 TeX 或 None, 右邊 TeX 或 None, 顏色), …]；
    每條是寬 width 的水平線（y = unit·能量，x 從 0 到 width），左標記放線左、右標記（例如簡併度）放線右。
    → VGroup，.lines、.left、.right（與 levels 同序；沒有標記的位置放隱形點）。並排幾欄時各自 .shift(RIGHT·x)。"""
    lines, left, right = VGroup(), VGroup(), VGroup()
    for E, lt, rt, col in levels:
        y = unit * float(E)
        ln = Line(np.array([0, y, 0]), np.array([width, y, 0])).set_stroke(col, stroke_width)
        lines.add(ln)
        if lt:
            left.add(Tex(lt, font_size=font_size).set_color(col).next_to(ln, LEFT, buff=left_buff))
        else:
            left.add(Dot(ln.get_left(), radius=0.001).set_opacity(0))
        if rt:
            right.add(Tex(rt, font_size=font_size).set_color(col).next_to(ln, RIGHT, buff=right_buff))
        else:
            right.add(Dot(ln.get_right(), radius=0.001).set_opacity(0))
    g = VGroup(lines, left, right)
    g.lines, g.left, g.right = lines, left, right
    return g


def weight_diagram(weights, unit: float = 0.6, roots=None, axes: bool = True, axis_labels=("H_1", "H_2"), extent: float = 2.2,
                   dot_radius: float = 0.07, color=YELLOW, root_color=TEAL, font_size: float = 20, ring_gap: float = 0.07) -> VGroup:
    """根圖／權圖（Fig. 5.8.2–5.8.5）：weights = {(x, y): 重數} 或 [(x, y), …]（重數 > 1 加同心圓圈）；roots = [(x, y, 標記 TeX 或 None), …]
    畫成從原點出發的箭頭。座標乘 unit；原點在 ORIGIN；軸長 ±extent。→ VGroup，.axes、.dots（與 weights 同序）、.rings、.roots、.root_labels。"""
    items = list(weights.items()) if isinstance(weights, dict) else [(w, 1) for w in weights]
    axg, dots, rings, rts, rlabs = VGroup(), VGroup(), VGroup(), VGroup(), VGroup()
    if axes:
        axg.add(Line(LEFT * extent, RIGHT * extent).set_stroke(GREY_C, 1.5), Line(DOWN * extent, UP * extent).set_stroke(GREY_C, 1.5))
        if axis_labels:
            axg.add(Tex(axis_labels[0], font_size=font_size).set_color(GREY_B).next_to(RIGHT * extent, DOWN, buff=0.1),
                    Tex(axis_labels[1], font_size=font_size).set_color(GREY_B).next_to(UP * extent, LEFT, buff=0.1))
    for (x, y), mlt in items:
        p = np.array([float(x) * unit, float(y) * unit, 0])
        dots.add(Dot(p, radius=dot_radius).set_fill(color, 1))
        for k in range(1, int(mlt)):
            rings.add(Circle(radius=dot_radius + k * ring_gap).set_stroke(color, 1.5).move_to(p))
    for r in roots or []:
        v = np.array([float(r[0]) * unit, float(r[1]) * unit, 0])
        rts.add(Arrow(ORIGIN, v, buff=0, thickness=3).set_color(root_color))
        if len(r) > 2 and r[2]:
            rlabs.add(Tex(r[2], font_size=font_size).set_color(root_color).move_to(v + 0.32 * v / max(np.linalg.norm(v), 1e-9)))
    g = VGroup(axg, rts, rlabs, rings, dots)
    g.axes, g.dots, g.rings, g.roots, g.root_labels = axg, dots, rings, rts, rlabs
    return g


def mnemonic_wheel(labels, radius: float = 1.5, start_deg: float = 180, direction: int = -1, font_size: float = 24, colors=None,
                   ring_color=GREY_B, angles=None) -> VGroup:
    """助記輪（Fig. 5.6.8、5.6.9、5.7.3、5.7.6）：圓環＋沿圓周的 irrep 標記（TeX）。第 k 個標記在角度 start_deg + direction·k·360/n
    （direction = −1 順時針、+1 逆時針；或直接給 angles 列表），顏色 irrep_color（或 colors={標記: 顏色}）。
    → VGroup，.ring、.slots[k]、.angle(k)（度）、.point(k, r)（半徑 r 處的座標，相對中心）。中心在 ORIGIN。
    讀法由 sym_rot 決定（wheel_read、WHEEL_*）；畫面上依序 .slots[k].animate.set_color(…) 點亮。"""
    n = len(labels)
    angs = list(angles) if angles is not None else [start_deg + direction * 360 * k / n for k in range(n)]
    ring = Circle(radius=radius).set_stroke(ring_color, 2.5)
    slots = VGroup()
    for k, lab in enumerate(labels):
        a = np.radians(angs[k])
        col = (colors or {}).get(lab, irrep_color(lab))
        slots.add(Tex(lab, font_size=font_size).set_color(col).move_to(1.22 * radius * np.array([np.cos(a), np.sin(a), 0])))
    g = VGroup(ring, slots)
    g.ring, g.slots, g.radius = ring, slots, radius
    g.angle = lambda k: angs[k % n]
    g.point = lambda k, r=radius: r * np.array([np.cos(np.radians(angs[k % n])), np.sin(np.radians(angs[k % n])), 0])
    return g


def wheel_marker(wheel, k: int, text: str | None = None, inner: bool = True, color=YELLOW, font_size: float = 18) -> VGroup:
    """助記輪上第 k 格的箭頭標（inner=True 在圈內指向該格、False 在圈外）；text 以 $ 開頭結尾的當 TeX，其餘用 zh。
    位置相對於輪心 ORIGIN（輪子 shift 過就對這個 VGroup 做同樣的 shift）。"""
    r = wheel.radius
    a = np.radians(wheel.angle(k))
    u = np.array([np.cos(a), np.sin(a), 0])
    arr = Arrow(0.5 * r * u, 0.88 * r * u, buff=0, thickness=3) if inner else Arrow(1.7 * r * u, 1.42 * r * u, buff=0, thickness=3)
    g = VGroup(arr.set_color(color))
    if text:
        t = Tex(text[1:-1], font_size=font_size) if text.startswith("$") else zh(text, font_size=font_size)
        t.set_color(color).move_to((0.3 if inner else 1.95) * r * u)
        g.add(t)
    return g


# ================================================================ 第 6 章（H601–H614）
A_COLOR = "#FF8C69"            # 粒子 A（大光暈）
B_COLOR = "#7FD1AE"            # 粒子 B（小點）
STRESS_COLOR = "#F2C14E"       # 應力箭頭
STRAIN_COLOR = "#8FD3FF"       # 變形後的立方體
EFIELD_COLOR = "#FFB347"       # 電場
BFIELD_COLOR = "#9AD1D4"       # 磁場／電流環
STIM_COLOR = GREY_B            # Lorentz 刺激曲線
RESP_COLOR = "#5B8DEF"         # Lorentz 響應曲線


class PairOcta(OctaSites):
    """兩個粒子在八面體上（Fig. 6.1.1）：OctaSites 加上粒子 A（大的空心光暈，A_COLOR）與 B（小實心點，B_COLOR）。
    PairOcta(radius=1.2, view=None, ball=0.13, **OctaSites 的參數)；.set_pair(i, j)（1 起算，|i, j⟩ = A 在 i、B 在 j）；.halo、.dot。
    位置用 View3D(center=…)。動畫：self.play(po.animate.set_pair(…)) 不要用，改 FadeOut/FadeIn 或直接 set_pair 後 self.wait。"""

    def __init__(self, radius: float = 1.2, view: View3D | None = None, ball: float = 0.13, **kw):
        super().__init__(radius=radius, view=view, ball=ball, **kw)
        self.halo = Circle(radius=2.2 * ball).set_stroke(A_COLOR, 4).set_fill(A_COLOR, 0.15)
        self.dot = Dot(radius=0.55 * ball).set_fill(B_COLOR, 1).set_stroke(WHITE, 1)
        self.pair = (1, 1)
        self.add(self.halo, self.dot)
        self.set_pair(1, 1)

    def set_pair(self, i: int, j: int):
        self.pair = (i, j)
        self.halo.move_to(self.point(i - 1)).set_z_index(20)
        pj = self.point(j - 1)
        self.dot.move_to(pj + (0.12 * RIGHT if i == j else 0)).set_z_index(21)
        return self

    def redraw(self):
        super().redraw()
        if hasattr(self, "halo"):
            self.set_pair(*self.pair)
        return self


class PairGrid(VGroup):
    """6×6 的兩粒子格子（(6.1.1)、(6.1.8) 的字典序；列 = 粒子 A 的 i、行 = 粒子 B 的 j）：PairGrid(n=6, cell=0.42, font_size=16,
    labels=True, values=None)。.cells[i][j]（Square，0 起算）、.texts[i][j]（values 給 n×n 的 TeX 字串時）、.row_labels、.col_labels；
    .color_cells(M, palette)：M[i][j] → 填色（palette = {值: 顏色}，沒有的值不上色）。"""

    def __init__(self, n: int = 6, cell: float = 0.42, font_size: float = 16, labels: bool = True, values=None, **kw):
        super().__init__(**kw)
        self.n, self.cell = n, cell
        self.cells = [[Square(side_length=cell).set_stroke(GREY_B, 1.2).set_fill(GREY_E, 0.4)
                       .move_to(np.array([(j - (n - 1) / 2) * cell, ((n - 1) / 2 - i) * cell, 0])) for j in range(n)] for i in range(n)]
        self.add(*[c for row in self.cells for c in row])
        self.texts = [[None] * n for _ in range(n)]
        if values is not None:
            for i in range(n):
                for j in range(n):
                    t = Tex(values[i][j], font_size=font_size).move_to(self.cells[i][j])
                    self.texts[i][j] = t
                    self.add(t)
        self.row_labels, self.col_labels = VGroup(), VGroup()
        if labels:
            for k in range(n):
                self.row_labels.add(Tex(str(k + 1), font_size=font_size + 2).set_color(A_COLOR).next_to(self.cells[k][0], LEFT, buff=0.1))
                self.col_labels.add(Tex(str(k + 1), font_size=font_size + 2).set_color(B_COLOR).next_to(self.cells[0][k], UP, buff=0.1))
            self.add(self.row_labels, self.col_labels)

    def color_cells(self, M, palette: dict, opacity: float = 0.75):
        for i in range(self.n):
            for j in range(self.n):
                c = palette.get(M[i][j])
                if c is not None:
                    self.cells[i][j].set_fill(c, opacity)
        return self


class StrainCube(VGroup):
    """投影的立方體（Fig. 6.3.2、6.3.4、6.3.6）：StrainCube(size=1.0, view=None, color=GREY_A, show_ref=True)。
    .set_strain(S)（3×3 浮點；頂點 v → (1 + S)v，show_ref=True 時原來的立方體留成淡虛線；給 updater 用 S(t) = t·S₀）、
    .face_arrows(T, scale=0.5, color=STRESS_COLOR, faces=(0,1,2,3,4,5)) → VGroup（面 2k、2k+1 = +x_k、−x_k 面心；箭頭 = ±scale·T·x̂_k，(6.3.12)）、
    .face_center(f)、.vertex(v)。前後用不透明度。位置用 View3D(center=…)；不要用 get_center 排版（同 CubeView）。"""

    def __init__(self, size: float = 1.0, view: View3D | None = None, color=GREY_A, show_ref: bool = True, **kw):
        super().__init__(**kw)
        self.size, self.view, self.S = size, view or View3D(-25, 70), np.zeros((3, 3))
        verts = [np.array(v, float) for v in itertools.product((1, -1), repeat=3)]
        self._edges = [(a, b) for a in verts for b in verts if np.sum(np.abs(a - b)) == 2 and tuple(a) < tuple(b)]
        self.edges = VGroup(*[Line(ORIGIN, RIGHT).set_stroke(color, 3) for _ in self._edges])
        self.ref = VGroup(*[DashedLine(ORIGIN, RIGHT, dash_length=0.06).set_stroke(GREY_C, 1.5) for _ in self._edges])
        if not show_ref:
            self.ref.set_opacity(0)
        self._show_ref = show_ref
        self.add(self.ref, self.edges)
        self.redraw()

    def vertex(self, v) -> np.ndarray:
        return self.view.project(self.size * (np.eye(3) + self.S) @ np.asarray(v, float))

    def face_center(self, f: int) -> np.ndarray:
        k, sgn = f // 2, (1 if f % 2 == 0 else -1)
        return self.view.project(self.size * (np.eye(3) + self.S) @ (sgn * np.eye(3)[k]))

    def redraw(self):
        pr = self.view.project
        for ln, rl, (a, b) in zip(self.edges, self.ref, self._edges):
            ln.put_start_and_end_on(self.vertex(a), self.vertex(b))
            mid = self.size * (np.eye(3) + self.S) @ ((a + b) / 2)
            ln.set_stroke(opacity=_depth_op(self.view, mid, 1.0, 0.35))
            rl.put_start_and_end_on(pr(self.size * a), pr(self.size * b))
        return self

    def set_strain(self, S):
        self.S = np.asarray(S, float)
        return self.redraw()

    def face_arrows(self, T, scale: float = 0.5, color=STRESS_COLOR, faces=(0, 1, 2, 3, 4, 5), thickness: float = 3) -> VGroup:
        T = np.asarray(T, float)
        g = VGroup()
        for f in faces:
            k, sgn = f // 2, (1 if f % 2 == 0 else -1)
            F = sgn * T[:, k]
            if np.linalg.norm(F) < 1e-9:
                continue
            c3 = self.size * (np.eye(3) + self.S) @ (sgn * np.eye(3)[k])
            a = Arrow(self.view.project(c3), self.view.project(c3 + scale * F), buff=0, thickness=thickness).set_color(color)
            a.set_opacity(_depth_op(self.view, c3, 1.0, 0.45))
            g.add(a)
        return g


def field_glyph(kind: str = "out", radius: float = 0.12, color=BFIELD_COLOR) -> VGroup:
    """場與電荷的小記號（Fig. 6.3.10–6.3.11）：kind = "out"（⊙，指出紙面）、"in"（⊗）、"+"、"-"（帶號的電荷圓）。"""
    c = Circle(radius=radius).set_stroke(color, 2)
    g = VGroup(c)
    if kind == "out":
        g.add(Dot(radius=0.28 * radius).set_fill(color, 1))
    elif kind == "in":
        d = 0.62 * radius
        g.add(Line(np.array([-d, -d, 0]), np.array([d, d, 0])).set_stroke(color, 2),
              Line(np.array([-d, d, 0]), np.array([d, -d, 0])).set_stroke(color, 2))
    else:
        d = 0.6 * radius
        g.add(Line(np.array([-d, 0, 0]), np.array([d, 0, 0])).set_stroke(color, 2))
        if kind == "+":
            g.add(Line(np.array([0, -d, 0]), np.array([0, d, 0])).set_stroke(color, 2))
    return g


def current_loop(center=ORIGIN, rx: float = 0.5, ry: float = 0.18, sense: int = 1, color=BFIELD_COLOR, width: float = 3) -> VGroup:
    """斜看的電流環（橢圓＋前緣一個箭頭；sense = +1 前緣向右＝從上往下看逆時針，B 向上；−1 反向）→ VGroup(環, 箭頭)。"""
    pts = [np.array(center) + np.array([rx * np.cos(t), ry * np.sin(t), 0]) for t in np.linspace(0, 2 * np.pi, 49)]
    ring = VMobject().set_points_as_corners(pts).set_stroke(color, width)
    tip_at = np.array(center) + np.array([0, -ry, 0])
    d = np.array([sense * 0.14, 0, 0])
    arr = Arrow(tip_at - d, tip_at + d, buff=0, thickness=3).set_color(color)
    return VGroup(ring, arr)


def trace_curve(axes, ts, ys, color=RESP_COLOR, stroke_width: float = 2.5) -> VMobject:
    """時間曲線（Fig. 6.5.1）：plot_axes 上的折線（一條 VMobject；n ≤ 1000 點，靜態畫、ShowCreation，不放進 updater）。
    資料用 sym_prod.lorentz_trace(f_s, tmax, n)。"""
    pts = [axes.c2p(float(t), float(y)) for t, y in zip(ts, ys)]
    return VMobject().set_points_as_corners(pts).set_stroke(color, stroke_width)


class SmithChart(VGroup):
    """Smith 圖（Fig. 6.5.3）：SmithChart(S=1, gamma=0.1, unit=0.25, deltas=(), extent=None)：實軸、虛軸、Γ 圓（直徑 S/Γ 沿虛軸）、
    deltas 裡每個 Δ 的 Δ 圓（直徑 S/Δ 沿實軸，虛線）。.point(Δ) → 近共振響應 S/(Δ − iΓ) 的畫面座標；.phasor(Δ, color) → 原點出發的箭頭；
    .gamma_circle、.delta_circles、.axes_lines。中心（原點）在 ORIGIN，可 shift。"""

    def __init__(self, S: float = 1.0, gamma: float = 0.1, unit: float = 0.25, deltas=(), extent=None, **kw):
        super().__init__(**kw)
        self.S, self.gamma, self.unit = S, gamma, unit
        D = S / gamma
        ext = extent or 0.65 * D
        self.axes_lines = VGroup(Line(-ext * unit * RIGHT, ext * unit * RIGHT).set_stroke(GREY_B, 1.5),
                                 Line(-0.1 * D * unit * UP, 1.1 * D * unit * UP).set_stroke(GREY_B, 1.5))
        self.gamma_circle = Circle(radius=D / 2 * unit).move_to(D / 2 * unit * UP).set_stroke(RESP_COLOR, 3)
        self.delta_circles = VGroup()
        for d in deltas:
            r = abs(S / d) / 2 * unit
            self.delta_circles.add(DashedVMobject(Circle(radius=r).move_to(np.sign(d) * r * RIGHT), num_dashes=36)
                                   .set_stroke(STRESS_COLOR, 2))
        self.add(self.axes_lines, self.gamma_circle, self.delta_circles)

    def point(self, delta: float) -> np.ndarray:
        z = self.S / (delta - 1j * self.gamma)
        return self.get_origin() + self.unit * np.array([z.real, z.imag, 0])

    def get_origin(self) -> np.ndarray:
        return self.axes_lines[0].get_center().copy()

    def phasor(self, delta: float, color=YELLOW) -> Arrow:
        return Arrow(self.get_origin(), self.point(delta), buff=0, thickness=3).set_color(color)


def _jt_lower(rho, phi, j, k, r):
    root = np.sqrt(np.maximum(j * j * rho ** 2 + 2 * j * r * rho ** 3 * np.cos(3 * phi) + r * r * rho ** 4, 0))
    return k * rho ** 2 - root


def hat_contours(center=ORIGIN, unit: float = 2.0, j: float = 1.0, k: float = 1.0, r: float = 0.0, levels=(-0.2, -0.1, 0.0, 0.2),
                 n: int = 60, colors=None, stroke_width: float = 2.5) -> VGroup:
    """墨西哥帽（Fig. 6.7.1）下層 ε−(x1, x2) = kρ² − (j²ρ² + 2jrρ³cos3φ + r²ρ⁴)^{1/2} 的等高線（(6.7.8)；自己算，不 import sym_prod）：
    每個 level 每個 φ 用二分法解 ρ（谷底內側、外側各一條）→ VGroup（每條閉曲線 n 點；level 高於谷底中心 0 時只有外側一條）。
    畫面座標 = center + unit·(x1, x2)。colors：{level: 顏色}。r ≠ 0 時 level 要高於鞍點（環狀）；低於鞍點的三個小島畫不出來（自動略過）。"""
    g = VGroup()
    phis = np.linspace(0, 2 * np.pi, n + 1)
    for lev in levels:
        inner, outer = [], []
        for ph in phis:
            rr = np.linspace(1e-4, 4 * abs(j) / max(k - abs(r), 1e-3) + 2, 400)
            vals = _jt_lower(rr, ph, j, k, r)
            imin = int(np.argmin(vals))

            def solve(a, b):
                fa = _jt_lower(a, ph, j, k, r) - lev
                for _ in range(60):
                    m = 0.5 * (a + b)
                    fm = _jt_lower(m, ph, j, k, r) - lev
                    if (fm > 0) == (fa > 0):
                        a, fa = m, fm
                    else:
                        b = m
                return 0.5 * (a + b)

            if vals[imin] < lev:
                outer.append(solve(rr[imin], rr[-1]))
                if vals[0] > lev:
                    inner.append(solve(rr[0], rr[imin]))
        col = (colors or {}).get(lev, TEAL)
        for branch in (inner, outer):
            if len(branch) == len(phis):
                pts = [np.array(center) + unit * np.array([b * np.cos(p), b * np.sin(p), 0]) for b, p in zip(branch, phis)]
                g.add(VMobject().set_points_as_corners(pts).set_stroke(col, stroke_width))
    return g


def hat_profile(axes, j: float = 1.0, k: float = 1.0, r: float = 0.0, phi_deg: float = 0.0, rho_max: float = 1.4, n: int = 60,
                colors=(GREY_A, TEAL), stroke_width: float = 3) -> VGroup:
    """Fig. 6.7.3 的截面：沿 φ（右半，x > 0）與 φ + 180°（左半，x < 0）的 ε+（上層）、ε−（下層）→ VGroup(上層, 下層)（各一條 2n 點的線）。"""
    out = VGroup()
    for sheet, col in ((+1, colors[0]), (-1, colors[1])):
        pts = []
        for x in np.linspace(-rho_max, rho_max, 2 * n):
            ph = np.radians(phi_deg) + (np.pi if x < 0 else 0)
            rho = abs(x)
            root = np.sqrt(max(j * j * rho ** 2 + 2 * j * r * rho ** 3 * np.cos(3 * ph) + r * r * rho ** 4, 0))
            pts.append(axes.c2p(x, k * rho ** 2 + sheet * root))
        out.add(VMobject().set_points_as_corners(pts).set_stroke(col, stroke_width))
    return out


# ======================================================================== 第 7 章
P_COLOR = "#FF8C69"            # 質子（大）
E_COLOR = "#7FD1AE"            # 電子（小）
J1_COLOR = "#FF8C69"           # j1
J2_COLOR = "#7FD1AE"           # j2
J3_COLOR = "#F2C14E"           # j3 = j1 + j2
POL_COLOR = "#8FD3FF"          # 偏振橢圓
STOKES_COLOR = "#F2C14E"       # Stokes 向量
POT_G_COLOR = "#5B8DEF"        # 鍵結 σg
POT_U_COLOR = "#E07A5F"        # 反鍵結 σu


class SpinTop(VGroup):
    """小陀螺（Fig. 7.1.1）：球、自旋箭頭（與向上夾 angle 度）、繞箭頭的進動虛線橢圓、可選的標籤（球下方）。"""

    def __init__(self, radius: float = 0.3, color=P_COLOR, arrow_len: float = 0.55, angle: float = 0.0, label_tex: str | None = None,
                 font_size: float = 24, show_orbit: bool = True, **kw):
        super().__init__(**kw)
        self.radius, self.arrow_len, self.color = radius, arrow_len, color
        self.ball = Circle(radius=radius).set_stroke(color, 2).set_fill(color, 0.35)
        self.arrow = Arrow(ORIGIN, UP, buff=0, stroke_width=4, max_tip_length_to_length_ratio=0.3).set_color(color)
        self.orbit = DashedVMobject(Ellipse(width=0.9 * radius, height=0.3 * radius), num_dashes=12).set_stroke(GREY_B, 1.5)
        if not show_orbit:
            self.orbit.set_opacity(0)
        self.add(self.ball, self.orbit, self.arrow)
        self.label = None
        if label_tex:
            self.label = Tex(label_tex, font_size=font_size).set_color(color).next_to(self.ball, DOWN, buff=0.12)
            self.add(self.label)
        self.set_angle(angle)

    def set_angle(self, angle: float):
        self.angle = float(angle)
        c = self.ball.get_center()
        d = np.array([np.sin(np.radians(self.angle)), np.cos(np.radians(self.angle)), 0.0])
        L = self.radius + self.arrow_len
        self.arrow.put_start_and_end_on(c - 0.35 * self.radius * d, c + L * d)
        self.orbit.move_to(c + 0.82 * L * d)
        return self


class SpinPair(VGroup):
    """兩個 SpinTop（質子大、電子小；Fig. 7.1.1 與 H701 的四個乘積態）。.set_ms(m1, m2)：+½ ↑、−½ ↓、0 → 水平（疊加態的示意）。"""

    def __init__(self, gap: float = 1.5, radii=(0.34, 0.22), colors=(P_COLOR, E_COLOR), labels=("p", "e"), font_size: float = 24, **kw):
        super().__init__(**kw)
        self.tops = [SpinTop(radius=radii[k], color=colors[k], label_tex=labels[k] if labels else None, font_size=font_size)
                     for k in range(2)]
        self.tops[1].shift(gap * RIGHT)
        self.tops[1].set_angle(0)
        self.add(*self.tops)

    def set_ms(self, m1: float, m2: float):
        for top, m in zip(self.tops, (m1, m2)):
            top.set_angle(0 if m > 0 else (180 if m < 0 else 90))
        return self


def _jlen(j) -> float:
    j = float(j)
    return float(np.sqrt(j * (j + 1)))


def vector_triangle(j1, j2, j3, unit: float = 0.5, start=ORIGIN, angle: float = 0.0, colors=(J1_COLOR, J2_COLOR, J3_COLOR),
                    labels: bool = True, font_size: float = 24) -> VGroup:
    """Fig. 7.1.5：|j1| + |j2| ＝ |j3| 的三角形（邊長 √(j(j+1))·unit）；j3 沿 angle 度，j1 在 j3 的左側（逆時針）。
    → VGroup(箭頭 j1, 箭頭 j2, 箭頭 j3[, 標籤們])。"""
    L1, L2, L3 = (_jlen(x) * unit for x in (j1, j2, j3))
    c = np.clip((L1 ** 2 + L3 ** 2 - L2 ** 2) / (2 * L1 * L3), -1, 1)
    a3 = np.radians(angle)
    a1 = a3 + np.arccos(c)
    o = np.array(start, float)
    t1 = o + L1 * np.array([np.cos(a1), np.sin(a1), 0])
    t3 = o + L3 * np.array([np.cos(a3), np.sin(a3), 0])
    arr = [Arrow(o, t1, buff=0, stroke_width=4).set_color(colors[0]), Arrow(t1, t3, buff=0, stroke_width=4).set_color(colors[1]),
           Arrow(o, t3, buff=0, stroke_width=5).set_color(colors[2])]
    g = VGroup(*arr)
    if labels:
        for k, (a, b, txt) in enumerate(((o, t1, "j_1"), (t1, t3, "j_2"), (o, t3, "j_3"))):
            mid = (a + b) / 2
            d = b - a
            nrm = np.array([-d[1], d[0], 0]) / max(np.linalg.norm(d), 1e-9)
            sgn = -1 if k == 2 else 1
            g.add(Tex(txt, font_size=font_size).set_color(colors[k]).move_to(mid + sgn * 0.25 * nrm))
    return g


class CouplingCone(VGroup):
    """Fig. 7.1.6：總角動量 J（長 √(J(J+1))，與 z 軸夾 Θ，cos Θ = M/|J|，在 xz 平面）、j1 繞 J 的錐（夾角 θ13 由餘弦定律）、
    j2 ＝ J − j1；z 軸上的古典範圍 m1 = |j1|cos(Θ ± θ13)（粗線＋兩個刻度）。.set_phase(deg)：j1 繞 J 轉（給 updater）。
    位置用 View3D(center=…)。"""

    def __init__(self, j1, j2, J, M, unit: float = 0.35, view: View3D | None = None, phase: float = 0.0,
                 colors=(J1_COLOR, J2_COLOR, J3_COLOR), show_range: bool = True, **kw):
        super().__init__(**kw)
        self.view = view or View3D(-25, 70)
        self.unit = unit
        L1, L2, LJ = _jlen(j1), _jlen(j2), _jlen(J)
        self.th13 = float(np.arccos(np.clip((LJ ** 2 + L1 ** 2 - L2 ** 2) / (2 * LJ * L1), -1, 1)))
        self.Th = float(np.arccos(np.clip(float(M) / LJ, -1, 1)))
        self.L1, self.LJ = L1, LJ
        self.Jhat = np.array([np.sin(self.Th), 0.0, np.cos(self.Th)])
        e1 = np.array([np.cos(self.Th), 0.0, -np.sin(self.Th)])
        self.e1, self.e2 = e1, np.cross(self.Jhat, e1)
        pr = self.view.project
        zL = 1.2 * max(LJ, L1) * unit
        self.zaxis = DashedLine(pr((0, 0, -zL)), pr((0, 0, zL)), dash_length=0.08).set_stroke(GREY_B, 1.6)
        cen = L1 * np.cos(self.th13) * self.Jhat
        rim = [unit * (cen + L1 * np.sin(self.th13) * (np.cos(t) * e1 + np.sin(t) * self.e2)) for t in np.linspace(0, 2 * np.pi, 49)]
        self.rim = _runs(rim, self.view, 0.9, 0.3, colors[0], 1.8)
        self.Jarrow = Arrow(pr(np.zeros(3)), pr(unit * LJ * self.Jhat), buff=0, stroke_width=5).set_color(colors[2])
        self.j1arrow = Arrow(ORIGIN, RIGHT, buff=0, stroke_width=4).set_color(colors[0])
        self.j2arrow = Arrow(ORIGIN, RIGHT, buff=0, stroke_width=4).set_color(colors[1])
        self.add(self.zaxis, self.rim, self.Jarrow, self.j1arrow, self.j2arrow)
        self.m1_range = (L1 * np.cos(self.Th + self.th13), L1 * np.cos(self.Th - self.th13))
        self.range_marks = VGroup()
        if show_range:
            a, b = (pr((0, 0, unit * m)) for m in self.m1_range)
            self.range_marks.add(Line(a, b).set_stroke(colors[0], 6, 0.7))
            for p in (a, b):
                self.range_marks.add(Line(p + 0.08 * LEFT, p + 0.08 * RIGHT).set_stroke(colors[0], 3))
            self.add(self.range_marks)
        self.set_phase(phase)

    def j1_vector(self, phase: float) -> np.ndarray:
        t = np.radians(phase)
        return self.L1 * (np.cos(self.th13) * self.Jhat + np.sin(self.th13) * (np.cos(t) * self.e1 + np.sin(t) * self.e2))

    def set_phase(self, phase: float):
        self.phase = float(phase)
        pr = self.view.project
        v1 = self.unit * self.j1_vector(phase)
        vJ = self.unit * self.LJ * self.Jhat
        self.j1arrow.put_start_and_end_on(pr(np.zeros(3)), pr(v1))
        self.j2arrow.put_start_and_end_on(pr(v1), pr(vJ))
        return self


class RESurface(VGroup):
    """RE 曲面的投影（不用真 3D）：每個能量一組等高線（單位向量折線 × 該組半徑），前面亮、後面淡；可加淡色經緯網與本體三軸。
    paths：[[n×3, …], …]（sym_coup.sphere_contours 的值、rigid_paths 的第三欄…）；radii：一個數或每組一個；colors：一個或每組一個。"""

    def __init__(self, paths, radii=1.6, view: View3D | None = None, colors=None, width: float = 2.0, front: float = 1.0,
                 back: float = 0.22, mesh=None, axes: bool = True, axis_labels=(r"\bar J_x", r"\bar J_y", r"\bar J_z"),
                 font_size: float = 22, **kw):
        super().__init__(**kw)
        self.view = view or View3D(-30, 65)
        self.paths = [list(ps) for ps in paths]
        n = len(self.paths)
        self.radii = list(radii) if np.iterable(radii) else [float(radii)] * n
        if colors is None:
            colors = color_gradient([BLUE_B, TEAL, YELLOW, RED_B], max(n, 2))
        self.colors = list(colors) if isinstance(colors, (list, tuple)) else [colors] * n
        self._w, self._f, self._b = width, front, back
        self.mesh_data = mesh or []
        self.mesh = VGroup()
        self.curves = VGroup()
        self.axes = VGroup()
        self._axis_labels, self._fs, self._show_axes = axis_labels, font_size, axes
        self.add(self.mesh, self.curves, self.axes)
        self._draw()

    @staticmethod
    def mesh_lines(radius_fn, n_lat: int = 5, n_lon: int = 8, n: int = 48) -> list:
        """淡色經緯網：radius_fn(n×3 單位向量) → 半徑陣列；回傳 [n×3 點]（已乘半徑）。"""
        out = []
        t = np.linspace(0, 2 * np.pi, n)
        for la in np.linspace(-60, 60, n_lat):
            c, s = np.cos(np.radians(la)), np.sin(np.radians(la))
            u = np.stack([c * np.cos(t), c * np.sin(t), np.full(n, s)], 1)
            out.append(u * np.asarray(radius_fn(u))[:, None])
        for k in range(n_lon):
            ph = np.pi * k / n_lon
            u = np.stack([np.cos(t) * np.cos(ph), np.cos(t) * np.sin(ph), np.sin(t)], 1)
            out.append(u * np.asarray(radius_fn(u))[:, None])
        return out

    def _draw(self):
        self.mesh.set_submobjects([m for pts in self.mesh_data for m in _runs(list(pts), self.view, 0.35, 0.12, GREY_C, 1.0)])
        curves = []
        for ps, r, col in zip(self.paths, self.radii, self.colors):
            g = VGroup()
            for p in ps:
                g.add(*_runs([r * np.asarray(q) for q in p], self.view, self._f, self._b, col, self._w))
            curves.append(g)
        self.curves.set_submobjects(curves)
        axs = []
        if self._show_axes:
            R = 1.25 * max(self.radii) if self.radii else 1.5
            for k, e in enumerate(np.eye(3)):
                ln = DashedLine(self.view.project(-R * e), self.view.project(R * e), dash_length=0.1).set_stroke(GREY_B, 1.4)
                lab = Tex(self._axis_labels[k], font_size=self._fs).set_color(GREY_A).move_to(self.view.project(1.12 * R * e))
                axs.append(VGroup(ln, lab))
        self.axes.set_submobjects(axs)

    def set_view(self, theta: float | None = None, phi: float | None = None):
        if theta is not None:
            self.view.theta = float(theta)
        if phi is not None:
            self.view.phi = float(phi)
        self._draw()
        return self


class PoincareSphere(SphereView):
    """Poincaré–Stokes 球（Fig. 7.5.1(b2)）：SphereView ＋ A、B、C 三軸（A→x、B→y、C→z）與標籤；
    .add_stokes(S, color) → VGroup(箭頭, 尖端的點)；.set_stokes(g, S) 更新（給 updater）。北極 +C ＝ 右旋圓偏振、赤道 ＝ 線偏振。"""

    def __init__(self, radius: float = 1.5, view: View3D | None = None, labels=("A", "B", "C"), font_size: float = 24, **kw):
        super().__init__(radius=radius, view=view or View3D(-35, 70), lats=(-60, -30, 0, 30, 60), n_lon=6, **kw)
        self.axis_labels = VGroup()
        for k, e in enumerate(np.eye(3)):
            self.add_axis(e, GREY_B, 1.25)
            self.axis_labels.add(Tex(labels[k], font_size=font_size).set_color(GREY_A).move_to(self.project(1.42 * e)))
        self.add(self.axis_labels)

    def add_stokes(self, S, color=STOKES_COLOR) -> VGroup:
        g = VGroup(Arrow(ORIGIN, RIGHT, buff=0, stroke_width=4).set_color(color), Dot(radius=0.06).set_fill(color, 1))
        self.set_stokes(g, S)
        self.add(g)
        return g

    def set_stokes(self, g: VGroup, S):
        S = np.asarray(S, float)
        u = S / max(np.linalg.norm(S), 1e-12)
        g[0].put_start_and_end_on(self.view.project(np.zeros(3)), self.project(u))
        g[1].move_to(self.project(u))
        return g


def polarization_ellipse(a: float, b: float, phi: float = 0.0, center=ORIGIN, unit: float = 1.0, color=POL_COLOR, n: int = 60,
                         arrow: bool = True, width: float = 3) -> VGroup:
    """(7.5.20a,b)：E(Φ) = (a cos Φ cos φ − b sin Φ sin φ, a cos Φ sin φ + b sin Φ cos φ)·unit；b > 0 時 Φ 增加是逆時針（右旋）。"""
    f = np.radians(phi)
    c = np.array(center, float)
    pts = [c + unit * np.array([a * np.cos(t) * np.cos(f) - b * np.sin(t) * np.sin(f), a * np.cos(t) * np.sin(f) + b * np.sin(t) * np.cos(f), 0])
           for t in np.linspace(0, 2 * np.pi, n)]
    g = VGroup(VMobject().set_points_as_corners(pts).set_stroke(color, width))
    if arrow and abs(b) > 1e-6:
        k = n // 8
        g.add(Arrow(pts[k], pts[k + 2], buff=0, stroke_width=3, max_tip_length_to_length_ratio=0.9).set_color(color))
    return g


def ellipse_from_stokes(S, center=ORIGIN, unit: float = 1.0, color=POL_COLOR, **kw) -> VGroup:
    """Stokes (S_A, S_B, S_C)（(7.5.25) 圓欄）→ ψ = ½asin(S_C/I)、φ = ½atan2(S_B, S_A)、a = √I cos ψ、b = √I sin ψ。"""
    S = np.asarray(S, float)
    I = float(np.linalg.norm(S))
    psi = 0.5 * np.arcsin(np.clip(S[2] / max(I, 1e-12), -1, 1))
    phi = 0.5 * np.degrees(np.arctan2(S[1], S[0]))
    return polarization_ellipse(np.sqrt(I) * np.cos(psi), np.sqrt(I) * np.sin(psi), phi, center, unit, color, **kw)


def projected_curve(points3, view: View3D, color=TEAL, width: float = 2.5, front: float = 1.0, back: float = 0.3) -> VGroup:
    """3D 折線（例如振子相空間 (x1, x2, p1) 的軌跡，Fig. 7.5.9–7.5.10）投影；前（depth ≥ 0）亮、後淡。"""
    return _runs([np.asarray(p, float) for p in points3], view, front, back, color, width)


def box_frame(X: float, Y: float, center=ORIGIN, unit: float = 1.0, angle: float = 0.0, color=GREY_B, width: float = 1.6) -> VMobject:
    """2X×2Y 的虛線盒（Fig. 7.5.5 的線偏振盒；Fig. 7.5.6 的傾斜盒用 angle，度）。"""
    rect = Rectangle(width=2 * X * unit, height=2 * Y * unit).rotate(np.radians(angle)).move_to(center)
    return DashedVMobject(rect, num_dashes=40).set_stroke(color, width)


def potential_curves(axes, xs, ys_list, colors=(POT_G_COLOR, POT_U_COLOR), labels=None, label_x=None, font_size: float = 22,
                     width: float = 3) -> VGroup:
    """plot_axes 上的位能曲線（Fig. 7.6.2）：每條一個 VMobject（建議 ≤ 60 點），標籤放在 x = label_x（預設最右）。→ VGroup(curves, labels)。"""
    curves, labs = VGroup(), VGroup()
    for k, ys in enumerate(ys_list):
        col = colors[k % len(colors)]
        curves.add(VMobject().set_points_as_corners([axes.c2p(float(x), float(y)) for x, y in zip(xs, ys)]).set_stroke(col, width))
        if labels:
            x0 = xs[-1] if label_x is None else label_x
            y0 = float(np.interp(x0, xs, ys))
            labs.add(Tex(labels[k], font_size=font_size).set_color(col).next_to(axes.c2p(x0, y0), RIGHT, buff=0.1))
    return VGroup(curves, labs)


def min_marker(axes, x: float, y: float, color=YELLOW, label_tex: str | None = None, font_size: float = 20, drop: bool = True) -> VGroup:
    """極小點：點＋（可選）落到 y = 0 的虛線＋標籤（點的右下方）。"""
    p = axes.c2p(x, y)
    g = VGroup(Dot(p, radius=0.06).set_fill(color, 1))
    if drop:
        g.add(DashedLine(axes.c2p(x, 0), p, dash_length=0.06).set_stroke(color, 1.5))
    if label_tex:
        g.add(Tex(label_tex, font_size=font_size).set_color(color).next_to(p, DR, buff=0.08))
    return g


# ======================================================================== 第 8 章
TRAJ_COLOR = "#8FD3FF"          # 古典軌跡
FRONT_COLOR = "#F2C14E"         # 等作用量波前
ENV_COLOR = "#E07A5F"           # 包絡／焦散
PACKET_COLOR = "#7FD1AE"        # 波包實部
PACKET_ENV_COLOR = GREY_B       # 波包的 |χ| 包絡
SPIN_COLOR = "#FF8C69"          # 自旋向量 S
DRESSED_COLOR = "#F2C14E"       # 穿衣能階
BARE_COLOR = GREY_B             # 裸能階


def polyline(axes, pts, color=TRAJ_COLOR, width: float = 2.5, dashed: bool = False, num_dashes: int = 40) -> VMobject:
    """plot_axes 上的一條折線（pts：(n, 2) 資料座標；建議 ≤ 60 點、靜態）；dashed=True 畫成虛線。"""
    m = VMobject().set_points_as_corners([axes.c2p(float(x), float(y)) for x, y in pts]).set_stroke(color, width)
    return DashedVMobject(m, num_dashes=num_dashes).set_stroke(color, width) if dashed else m


def line_family(axes, lines, x0: float, x1: float, color=GREY_B, width: float = 1.6, highlight=None, hcolor=YELLOW) -> VGroup:
    """直線族 Y = sX + b（lines = [(s, b)]，sym_dyn.line_family 的輸出；Fig. 8.1.2）在 x ∈ [x0, x1] 的線段；highlight = 要標亮的序號。"""
    g = VGroup()
    for k, (s, b) in enumerate(lines):
        col, w = (hcolor, width + 1.2) if highlight is not None and k in highlight else (color, width)
        g.add(Line(axes.c2p(x0, s * x0 + b), axes.c2p(x1, s * x1 + b)).set_stroke(col, w))
    return g


def trajectory_fronts(axes, trajs, fronts, traj_color=TRAJ_COLOR, front_color=FRONT_COLOR, traj_width: float = 1.6,
                      front_width: float = 2.5) -> VGroup:
    """Fig. 8.1.6：軌跡族（trajs：[(n, 2)]，sym_dyn.trajectory）＋ 等 S_H 波前（fronts：[(n, 2)] 或 {level: [段…]}，sym_dyn.fig816_curves）。
    → VGroup(軌跡們, 波前們)；波前是軌跡的法線（p·dr = 0）。"""
    tg = VGroup(*[polyline(axes, t, traj_color, traj_width) for t in trajs])
    fg = VGroup()
    items = fronts.values() if isinstance(fronts, dict) else [[f] for f in fronts]
    for segs in items:
        for seg in segs:
            if len(seg) > 1:
                fg.add(polyline(axes, seg, front_color, front_width))
    return VGroup(tg, fg)


class PhasePhasor(VGroup):
    """相空間的相量（Fig. 8.2.1）：圓心 center、半徑 radius 的虛線軌道（|α| = amp 對到 radius）、兩軸、從圓心到 α 的箭頭、
    尖端一個半透明 Gauss 圓斑（半徑 blob，代表最小不確定的波包）；.set_alpha(z)（給 updater；α_t = α e^{−it} 是順時針）、
    .arrow、.blob、.orbit、.axes_lines、.labels。座標：右 = Re α（∝ x）、上 = Im α（∝ p）。"""

    def __init__(self, radius: float = 1.2, center=ORIGIN, amp: float = 1.0, blob: float = 0.22, color=PHASOR_COLOR,
                 labels=(r"\mathrm{Re}\,\alpha", r"\mathrm{Im}\,\alpha"), font_size: float = 22):
        super().__init__()
        self.c = np.array(center, float)
        self.radius, self.amp = radius, amp
        ext = radius * 1.35
        self.axes_lines = VGroup(Line(self.c + ext * LEFT, self.c + ext * RIGHT), Line(self.c + ext * DOWN, self.c + ext * UP)).set_stroke(GREY_C, 1.5)
        self.orbit = DashedVMobject(Circle(radius=radius).move_to(self.c), num_dashes=48).set_stroke(GREY_B, 1.5)
        self.labels = VGroup(Tex(labels[0], font_size=font_size).set_color(GREY_A).next_to(self.c + ext * RIGHT, DOWN, buff=0.1),
                             Tex(labels[1], font_size=font_size).set_color(GREY_A).next_to(self.c + ext * UP, RIGHT, buff=0.1))
        self.blob = Circle(radius=blob).set_fill(color, 0.35).set_stroke(color, 1.5)
        self.arrow = Arrow(self.c, self.c + radius * RIGHT, buff=0, stroke_width=4).set_color(color)
        self.add(self.axes_lines, self.orbit, self.labels, self.blob, self.arrow)
        self.set_alpha(amp)

    def point(self, z) -> np.ndarray:
        z = complex(z)
        return self.c + self.radius / self.amp * np.array([z.real, z.imag, 0.0])

    def set_alpha(self, z):
        p = self.point(z)
        self.blob.move_to(p)
        if np.linalg.norm(p - self.c) > 1e-6:
            self.arrow.put_start_and_end_on(self.c, p)
        return self


def snapshot_stack(xs, rows, unit_x: float = 0.45, unit_y: float = 0.5, gap: float = 0.45, top=ORIGIN, re_color=PACKET_COLOR,
                   env_color=PACKET_ENV_COLOR, width: float = 2.0, envelope: bool = True, baseline: bool = True) -> VGroup:
    """波包快照列（Fig. 8.2.2、8.3.2(a)）：rows = [複數陣列（在 xs 上）]，第 k 列畫在 y = top_y − k·gap：實線 Re χ、虛線 |χ|（envelope）、
    淡色基線；x = 0 對到 top_x。每列建議 ≤ 120 點。→ VGroup(列們)，每列 VGroup(Re, [包絡], [基線])。"""
    xs = np.asarray(xs, float)
    t0 = np.array(top, float)
    out = VGroup()
    for k, psi in enumerate(rows):
        y0 = t0[1] - k * gap
        psi = np.asarray(psi, complex)
        row = VGroup()
        row.add(VMobject().set_points_as_corners([np.array([t0[0] + x * unit_x, y0 + v * unit_y, 0]) for x, v in zip(xs, psi.real)])
                .set_stroke(re_color, width))
        if envelope:
            env = VMobject().set_points_as_corners([np.array([t0[0] + x * unit_x, y0 + v * unit_y, 0]) for x, v in zip(xs, abs(psi))])
            row.add(DashedVMobject(env, num_dashes=30).set_stroke(env_color, 1.2))
        if baseline:
            row.add(Line(np.array([t0[0] + xs[0] * unit_x, y0, 0]), np.array([t0[0] + xs[-1] * unit_x, y0, 0])).set_stroke(GREY_D, 1))
        out.add(row)
    return out


class CrankSphere(SphereView):
    """曲柄向量 ω 與自旋向量 S 的 3D 圖（Fig. 8.5.1、7.5.3）：SphereView（球面代表 |v| = ½，即自旋向量的長度）＋ 三軸與標籤；
    .add_vector(v, color, label_tex) → VGroup(箭頭, 點[, 標籤])（v 不歸一化：畫面長 = 2·radius·|v|）、.set_vector(g, v)（給 updater）、
    .add_trail(points3, color)（S 的軌跡，投影、前亮後淡；≤ 200 點、靜態）、
    CrankSphere.cone_rim(axis, half_angle, length, n)（進動錐的口：繞 axis、半頂角 half_angle 度、離原點 length → 3D 點列，交給 add_trail）。
    注意：set_view 只重畫球與軸，第 8 章加的向量要再呼叫 set_vector。"""

    def __init__(self, radius: float = 1.5, view: View3D | None = None, labels=("x", "y", "z"), font_size: float = 22, **kw):
        super().__init__(radius=radius, view=view or View3D(-35, 70), lats=(-60, -30, 0, 30, 60), n_lon=6, **kw)
        self.axis_labels = VGroup()
        for k, e in enumerate(np.eye(3)):
            self.add_axis(e, GREY_B, 1.25)
            self.axis_labels.add(Tex(labels[k], font_size=font_size).set_color(GREY_A).move_to(self.project(1.42 * e)))
        self.add(self.axis_labels)
        self.vectors = VGroup()
        self.add(self.vectors)

    def vec_point(self, v) -> np.ndarray:
        return self.view.project(2 * self.radius * np.asarray(v, float))          # |v| = ½ 落在球面上

    def add_vector(self, v, color=SPIN_COLOR, label_tex: str | None = None, font_size: float = 22) -> VGroup:
        g = VGroup(Arrow(ORIGIN, RIGHT, buff=0, stroke_width=4).set_color(color), Dot(radius=0.05).set_fill(color, 1))
        if label_tex:
            g.add(Tex(label_tex, font_size=font_size).set_color(color))
        self.set_vector(g, v)
        self.vectors.add(g)
        return g

    def set_vector(self, g: VGroup, v):
        v = np.asarray(v, float)
        o, p = self.view.project(np.zeros(3)), self.vec_point(v)
        if np.linalg.norm(p - o) < 1e-4:
            p = o + 1e-3 * RIGHT
        g[0].put_start_and_end_on(o, p)
        g[1].move_to(p)
        if len(g) > 2:
            d = p - o
            g[2].move_to(p + 0.28 * d / max(np.linalg.norm(d), 1e-6))
        return g

    def add_trail(self, points3, color=SPIN_COLOR, width: float = 2.0, front: float = 1.0, back: float = 0.3) -> VGroup:
        tr = _runs([2 * self.radius * np.asarray(q, float) for q in points3], self.view, front, back, color, width)
        self.add(tr)
        return tr

    @staticmethod
    def cone_rim(axis, half_angle: float, length: float = 0.5, n: int = 48) -> list:
        a = np.asarray(axis, float)
        a = a / np.linalg.norm(a)
        t = np.cross(a, [0, 0, 1.0] if abs(a[2]) < 0.9 else [1.0, 0, 0])
        t /= np.linalg.norm(t)
        u = np.cross(a, t)
        h = np.radians(half_angle)
        return [length * (np.cos(h) * a + np.sin(h) * (np.cos(s) * t + np.sin(s) * u)) for s in np.linspace(0, 2 * np.pi, n)]


def dressed_ladder(blocks, unit: float = 1.0, base=ORIGIN, bare_w: float = 0.7, dressed_w: float = 0.8, gap: float = 0.55,
                   font_size: float = 18, labels: bool = True) -> VGroup:
    """Jaynes–Cummings 穿衣能階梯（Fig. 8.7.1）：blocks = [(N, 下穿衣, 上穿衣, 裸 |N⟩|↓⟩, 裸 |N−1⟩|↑⟩)]（sym_dyn.jc_ladder 的輸出；
    N = 0 只有基態）；左欄裸 |N⟩|↓⟩、右欄裸 |N−1⟩|↑⟩、中欄穿衣雙重態，虛線連到最近的穿衣能階。y = base_y + unit·能量。
    → VGroup(左裸, 右裸, 穿衣, 連線, 標籤)；.dressed[k] 是第 k 塊的 VGroup(下[, 上])。"""
    b = np.array(base, float)
    xl, xc, xr = b[0] - dressed_w / 2 - gap - bare_w / 2, b[0], b[0] + dressed_w / 2 + gap + bare_w / 2
    left, right, mid, links, labs = VGroup(), VGroup(), VGroup(), VGroup(), VGroup()

    def hl(x, e, w, col, sw=3):
        return Line(np.array([x - w / 2, b[1] + unit * e, 0]), np.array([x + w / 2, b[1] + unit * e, 0])).set_stroke(col, sw)

    def pt(x, e):
        return np.array([x, b[1] + unit * e, 0])
    for (N, lo, hi, e_dn, e_up) in blocks:
        left.add(hl(xl, e_dn, bare_w, BARE_COLOR))
        if labels:
            labs.add(Tex(r"|%d\rangle|{\downarrow}\rangle" % N, font_size=font_size).set_color(GREY_A).next_to(pt(xl - bare_w / 2, e_dn), LEFT, buff=0.08))
        if N == 0:
            mid.add(VGroup(hl(xc, lo, dressed_w, DRESSED_COLOR)))
            links.add(DashedLine(pt(xl + bare_w / 2, e_dn), pt(xc - dressed_w / 2, lo), dash_length=0.05).set_stroke(GREY_C, 1))
            continue
        right.add(hl(xr, e_up, bare_w, BARE_COLOR))
        if labels:
            labs.add(Tex(r"|%d\rangle|{\uparrow}\rangle" % (N - 1), font_size=font_size).set_color(GREY_A).next_to(pt(xr + bare_w / 2, e_up), RIGHT, buff=0.08))
        mid.add(VGroup(hl(xc, lo, dressed_w, DRESSED_COLOR), hl(xc, hi, dressed_w, DRESSED_COLOR)))
        for e_b, x_b, side in ((e_dn, xl + bare_w / 2, -1), (e_up, xr - bare_w / 2, 1)):
            near = lo if abs(e_b - lo) <= abs(e_b - hi) else hi
            links.add(DashedLine(pt(x_b, e_b), pt(xc + side * dressed_w / 2, near), dash_length=0.05).set_stroke(GREY_C, 1))
    g = VGroup(left, right, mid, links, labs)
    g.dressed = mid
    return g
