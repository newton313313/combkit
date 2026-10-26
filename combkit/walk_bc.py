# 由 combkit/tools/sync_from_walk017.py 從 walk017 專案的 codes/utils/walk_bc.py 複製；請改 walk017 那份再同步。
"""
walk_bc.py —— 《Elements of the Random Walk》第 4 章 Boundary conditions, steady state, and the electrostatic analogy 的數學
（純 Python，不 import manim；只依賴同資料夾的 walk_gf、walk_rec，以便原樣複製進 combkit）

約定（書上式號與筆誤見 ref/md/ch04_formulas.md）：
  一維格子、步長 Δ = a = 1。吸收牆是「踩到就消失」的格點：一面牆在 x = 0（起點 x0 ≥ 1、鏡像點 −x0）；
  兩面牆在 0 與 L（內點 1..L−1，共 L−1 個）——和書上 (4.10) q = nπ/L 一致：格子特徵值 2cos(mπ/L)、m = 1..L−1
  （不是「L 個內點、π/(L+1)」那種約定）。連續極限把 2cos q 換成 2 − q²（(4.14)）。
  §4.3 的補充照書上：格點 n ≥ 0、牆緊鄰在 −1。三維都用簡單立方（sc，每步 6 個方向，z_c = 1/6）。
  穩定態：每步在源放 1 個走路人（s0 = 1）；一維格子的 β = ½、三維 β = D = 1/6（每步）；連續的式子用 D（＝書上的 β）。
  「÷ 2^N」「÷ 6^N」的函數回傳走法數除以無限制的走法數（避免大數）。

§4.1.1 鏡像走路人 (4.1)–(4.8)，Fig. 4.1
    walks(N, x0)                    全部 2^N 條走法（位置 list）；touches_wall(path)：有沒有踩到 ≤ 0
    survivors_brute(N, x0)          {終點 x: N 步都沒踩到 0 的走法數}（暴力列舉）
    survivors_image(N, x0)          同上，(4.1)：C(N; x−x0) − C(N; x+x0)（第二種算法）；survivors_dp(N, x0) 逐步遞迴（第三種）
    reflect_first(path)             第一次踩到 0 之前的那段對 0 鏡射 → 從 −x0 出發的路線
    reflection_pairs(N, x0)         [(碰過牆、終點 > 0 的路線, 鏡像路線)]（一一對應）
    survival_prob(N, x0)            N 步的存活機率（Fraction）= P(1−x0 ≤ S_N ≤ x0)；survival_asym(N, x0) = x0√(2/(πN))
    survival_3d(N, x0)              sc 格子、吸收平面 x = 0：Σ_k C(N,k)(1/3)^k(2/3)^{N−k}·一維存活（x 只在 1/3 的步數裡動）；
                                    survival_3d_dp 用「懶惰的一維漫步」逐步遞迴（第二種）；survival_3d_asym(N, x0) = √(6x0²/(πN))（(4.5) 後）
    wall_reduction_1d(N)            一維版 (4.6)–(4.8)：Σ_{x0≥1}(1 − 存活機率)（Fraction）；wall_reduction_brute(N) = E[max S_k]（暴力列舉）；
                                    wall_reduction_asym(N, d=1) = √(2N/(πd))（d = 3 是書上的 √(2N/3π)）
§4.1.2 兩面牆 (4.9)–(4.15)，Ex. 4.1–4.4
    box_dp(N, L, x0)                {x: 從 x0 走 N 步到 x、沒碰 0 或 L 的走法數}（轉移矩陣＝逐步遞迴）
    box_images(N, L, x0, x)         Ex. 4.3：鏡像點 2kL + x0（正）、2kL − x0（負）；box_image_points(L, x0, K) 列出位置
    box_modes(N, L, x0, x)          (4.15) 的格子精確版 Σ_m (2/L) sin(mπx0/L) sin(mπx/L) (2cos(mπ/L))^N（float）
    box_eigs(L)                     [2cos(mπ/L)]，m = 1..L−1；box_mode(L, m) 正規化的 √(2/L) sin(mπx/L)；box_eig_residual(L)
    box_gf(z, L)                    生成函數矩陣 (I − zT)^{−1}；box_amplitude(z, L, m, x2) 投影到 sin(qx1) 的係數（數值）；
                                    box_amplitude_formula(…) = √(2/L) sin(q x2)/(1 − 2z cos q)（更正後的 (4.12)(4.13)）
§4.1.3 長短時間 (4.16)–(4.27)，Fig. 4.2
    box_total(N, L)                 從任一內點出發、終點不限的存活走法總數（int）；box_survival(N, L) = ÷ (2^N(L−1))（Fraction）
    box_total_modes(N, L)           ÷ 2^N 的格子模態和 Σ_{m 奇} (2/L)cot²(mπ/2L) cos^N(mπ/L)
    box_survival_long(N, L, kind)   對數圖的直線：'lattice' 最低模態、'cont' (8/π²)e^{−π²N/2L²}、'book' (8/π²)e^{−π²N/L²}
    box_total_cont(N, L)            (4.19) ÷ 2^N（連續）；box_long(N, L, book=False) = (8L/π²)e^{−π²N/(2L²)}（book：書上的 e^{−π²N/L²}）
    decay_rate(L, kind)             每步衰減率：'lattice' −ln cos(π/L)、'cont' π²/(2L²)、'book' π²/L²
    odd_inv_squares(M)              Σ 1/(2m+1)²（部分和＋尾巴）→ π²/8；zeta2(M) → π²/6（(4.21)(4.22)）
§4.1.4 球與立方體 (4.28)–(4.37)，Ex. 4.5
    j0(x)、j0_zero(n)               sin x/x 與第 n 個零點（二分法；= nπ）
    sphere_weight(n, R)             (4.34) 8R³/(πn²)；sphere_weight_numeric 數值積分（更正後的 (4.31)）；sphere_weight_book 照書上 (4.31)
    sphere_total_cont(N, R, book)   (4.35) ÷ 6^N（更正：e^{−π²n²N/(6R²)}）；sphere_rate(R, book) = π²/(6R²)（書上 π²/R²）
    ball_sites(R)、ball_totals(R, Nmax)、ball_rate(R)   格子球（x²+y²+z² < R²）的存活總數 ÷ 6^N 與長時間衰減率（核對 1/6）
    cube_total(N, L)                Ex. 4.5：立方體 0 < x, y, z < L 的存活走法總數（int，逐步遞迴）；cube_total_modes(N, L) ÷ 6^N（模態）
    cube_total_cont(N, L)、cube_long(N, L)   連續版與長時間 (512L³/π⁶)e^{−π²N/(2L²)}
    PREFACTOR                       長時間前因子 ÷ 體積：{'interval': 8/π², 'cube': 512/π⁶, 'sphere': 6/π²}
    same_volume_rate_ratio()        同體積：立方體 ÷ 球的衰減率 = 3(3/(4π))^{2/3} = 1.155
§4.2–4.2.1 穩定態、Fick、帳篷 (4.38)–(4.49)，Ex. 4.6
    laplacian_coeff(a)              [(1/6)Σ f(r+aê_i) − f(r)]/(a²∇²f) → 1/6（更正後的 (4.40)；書上 1/3）
    tent_lattice(L, x0)             一維格子穩定態 w(0..L)（Fraction，線性方程精確解）；tent_formula(L, x0, x) = (4.49)（s0 = 1、β = ½）
    absorb_left(L, x0) = (L−x0)/L   Ex. 4.6（賭徒破產）；absorb_left_solve 解調和方程、absorb_left_dp 逐步遞迴
    wall_fluxes(L, x0)              穩定態流進左、右牆的人數 (w(1)/2, w(L−1)/2)
§4.2.2 點源 (4.50)–(4.56)
    heat_kernel(r, t, D)、time_integral(r, D)   ∫_0^∞ (4πDt)^{−3/2}e^{−r²/4Dt} dt（數值）→ 1/(4πDr)
    K_const(a) = √(πa)、K_const_numeric(a)     (4.56) 的 K
    ine_list(nmax, y)               [e^{−y} I_n(y)]（級數／Miller 逆向遞迴）
    sc_green(n1, n2, n3)            sc 的 G(z_c; r) = ∫_0^∞ e^{−t} Π I_{n_i}(t/3) dt＝在原點每步放 1 人時的穩定態 w(r)
    sc_green_asym(r) = 3/(2πr)      更正後的 (4.41)；sc_green_book(r) = 3/(4πr)（照書上 3/a²）；sc_laplace_residual(n1, n2, n3) 核對 (4.39)
§4.2.3 吸收牆旁的點源 (4.57)–(4.60)，Fig. 4.4–4.5
    conc_wall(x, y, z, l)           點源 (l,0,0)＋鏡像 (−l,0,0)，牆 x = 0：(1/4π)(1/|r−r_s| − 1/|r−r_i|)；conc_dipole 遠場 2lx/r³
    plane_flux(X, l)                穿過平面 x = X 的總流量（數值）；plane_flux_exact(X, l)：X > l 為 0、0 < X < l 為 −1
    dipole_flux_integral(x)         (4.60)（數值，= 0）
§4.2.4–4.2.5 吸收球、電容 (4.61)–(4.66)，Fig. 4.6–4.7
    conc_sphere(r, r0, cinf)        c∞(1 − r0/r)；sphere_flux(r0, D, cinf) = 4πDr0c∞；sphere_flux_numeric(r, r0) 任何球面都相同
    kelvin_image(R, r0)             (r0²/R, −r0/R)；hit_prob(R, r0) = r0/R；hit_prob_numeric 鏡像場的球面流量（數值）；sphere_surface_max ≈ 0
    cell_rmax(D, cinf, q)           √(3Dc∞/q)；CELL_EXAMPLE（假設值）→ 3.9×10⁻⁵ m
    capacitance_sphere(r0) = r0；disk_capacitance(a) = 2a/π（書外）；disk_capacitance_numeric、disk_potential_offcenter（圓盤等電位 π²）
§4.2.6 布滿受器的球（Berg & Purcell 1977，書外補充）
    berg_purcell(n, a, r0) = na/(na+πr0)；berg_purcell_series（電容串聯）；bp_half(a, r0) = πr0/a；bp_area_fraction = na²/(4r0²)
    BP_EXAMPLE                      r0 = 1 µm、a = 1 nm（n_half = 3142，占表面 0.079%）
§4.3 補充：吸收邊界 (4:S-1)–(4:S-3)
    half_counts(N, x, y)            格點 n ≥ 0、牆在 −1：從 x 走 N 步到 y 的走法數（逐步遞迴）；half_counts_image = C(N; y−x) − C(N; y+x+2)
    half_gf(z, x, y)                Ex. 4.1 的格子版 (λ^{|x−y|} − λ^{x+y+2})/√(1−4z²)；half_gf_residual(z, x) 核對 (4:S-1)–(4:S-3)
    catalan(k)                      從牆邊 0 出發回到 0 的 2k 步走法數＝Catalan 數 1, 1, 2, 5, 14, 42
    half_visits_critical(x, y)      z → ½：到訪次數 2(min(x, y)+1)（帳篷在 L → ∞ 的極限）

check() 驗證全部（每個結果至少兩種獨立算法、書上的數字、筆誤的更正、plan_ch04.md 場景卡與講義的每個數字）；直接執行本檔會跑 check()（約 1 秒）。
最慢的呼叫：sc_green 第一次建表約 0.15 秒、ball_rate(8) 約 0.15 秒、survival_3d_dp(1000, 3) 約 0.3 秒，其餘都很快。
"""
import math
from fractions import Fraction
from functools import lru_cache

if __package__:
    from . import walk_gf as GF
    from . import walk_rec as R
else:                                   # python utils/walk_bc.py
    import walk_gf as GF
    import walk_rec as R

PI = math.pi
_simpson = R._simpson


# =========================================================================== §4.1.1 鏡像走路人（牆在 0）


def walks(N: int, x0: int) -> list[list[int]]:
    """從 x0 出發的全部 2^N 條一維走法（位置 list，長 N+1）。"""
    return R.walks_1d(N, x0)


def touches_wall(path) -> bool:
    """這條路線有沒有踩到牆（位置 ≤ 0）。"""
    return any(p <= 0 for p in path)


def survivors_brute(N: int, x0: int) -> dict:
    """{終點 x: 從 x0 走 N 步、從沒踩到 0 的走法數}（暴力列舉）。"""
    out = {}
    for p in walks(N, x0):
        if not touches_wall(p):
            out[p[-1]] = out.get(p[-1], 0) + 1
    return dict(sorted(out.items()))


def survivors_image(N: int, x0: int) -> dict:
    """(4.1) 的一維格子版：C(N; x−x0) − C(N; x+x0)（鏡像點 −x0），x ≥ 1。"""
    out = {}
    for x in range(1, x0 + N + 1):
        c = GF.count1(N, x - x0) - GF.count1(N, x + x0)
        if c:
            out[x] = c
    return out


def survivors_dp(N: int, x0: int) -> dict:
    """逐步遞迴（每步把踩到 0 的刪掉；第三種算法）。"""
    cur = {x0: 1}
    for _ in range(N):
        nxt = {}
        for x, c in cur.items():
            for y in (x - 1, x + 1):
                if y > 0:
                    nxt[y] = nxt.get(y, 0) + c
        cur = nxt
    return dict(sorted(cur.items()))


def reflect_first(path: list[int]) -> list[int]:
    """把第一次踩到 0 之前的那一段對 0 鏡射（Fig. 4.1）：碰過牆的路線 ↔ 從 −x0 出發的路線。"""
    k = next(i for i, p in enumerate(path) if p == 0)
    return [-p for p in path[:k]] + list(path[k:])


def reflection_pairs(N: int, x0: int) -> list[tuple]:
    """[(碰過牆、終點 > 0 的路線, 它的鏡像路線)]：鏡像路線從 −x0 出發、終點相同。"""
    return [(p, reflect_first(p)) for p in walks(N, x0) if touches_wall(p) and p[-1] > 0]


def survival_prob(N: int, x0: int) -> Fraction:
    """N 步都沒踩到 0 的機率（精確）= P(1−x0 ≤ S_N ≤ x0)，S_N 是無限制的位移。"""
    return Fraction(sum(GF.count1(N, d) for d in range(1 - x0, x0 + 1)), 2 ** N)


def survival_asym(N: float, x0: float) -> float:
    """一維 (4.5) 的類比：x0·√(2/(πN))（= 2x0 × 高斯密度 (2πN)^{−1/2}）。"""
    return x0 * math.sqrt(2 / (PI * N))


def _log_comb(n: int, k: int) -> float:
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


def _surv_float(N: int, x0: int) -> float:
    tot = 0.0
    for d in range(1 - x0, x0 + 1):
        if abs(d) <= N and (N + d) % 2 == 0:
            tot += math.exp(_log_comb(N, (N + d) // 2) - N * math.log(2))
    return tot


def survival_3d(N: int, x0: int) -> float:
    """sc 格子、吸收平面 x = 0、從 (x0, 0, 0) 出發：N 步都沒踩到平面的機率。
    x 座標只在 1/3 的步數裡移動：Σ_k C(N,k)(1/3)^k(2/3)^{N−k}·(一維 k 步的存活機率)。"""
    tot = 0.0
    for k in range(N + 1):
        w = math.exp(_log_comb(N, k) + k * math.log(1 / 3) + (N - k) * math.log(2 / 3))
        tot += w * _surv_float(k, x0)
    return tot


def survival_3d_dp(N: int, x0: int) -> float:
    """同上，用「懶惰的一維漫步」逐步遞迴（每步 2/3 不動、1/6 向左、1/6 向右；第二種算法）。"""
    width = x0 + N + 2
    p = [0.0] * width
    p[x0] = 1.0
    for _ in range(N):
        q = [0.0] * width
        for x in range(1, width - 1):
            v = p[x]
            if v:
                q[x] += v * 2 / 3
                q[x - 1] += v / 6
                q[x + 1] += v / 6
        q[0] = 0.0
        p = q
    return sum(p)


def survival_3d_asym(N: float, x0: float) -> float:
    """書上 (4.5) 之後：√(6x0²/(πN))。"""
    return math.sqrt(6 * x0 * x0 / (PI * N))


def wall_reduction_1d(N: int) -> Fraction:
    """一維 (4.6)–(4.8)：Σ_{x0≥1} (1 − 存活機率)＝因為牆而少掉的走法數 ÷ 2^N（精確）。"""
    return sum((1 - survival_prob(N, x0) for x0 in range(1, N + 1)), Fraction(0))


def wall_reduction_brute(N: int) -> Fraction:
    """同上的第二種算法：等於 N 步漫步的最大值的期望值 E[max_k S_k]（暴力列舉）。"""
    return Fraction(sum(max(p) for p in R.walks_1d(N, 0)), 2 ** N)


def wall_reduction_asym(N: float, d: int = 1) -> float:
    """(4.8)：一維 √(2N/π)；三維（sc，x 方向的變異數 N/3）√(2N/(3π))。"""
    return math.sqrt(2 * N / (PI * d))


# =========================================================================== §4.1.2 兩面吸收牆（牆在 0 與 L）


def box_dp(N: int, L: int, x0: int) -> dict:
    """{x: 從 x0 走 N 步到 x、從沒碰到 0 或 L 的走法數}（轉移矩陣＝逐步遞迴）。"""
    cur = {x0: 1}
    for _ in range(N):
        nxt = {}
        for x, c in cur.items():
            for y in (x - 1, x + 1):
                if 0 < y < L:
                    nxt[y] = nxt.get(y, 0) + c
        cur = nxt
    return dict(sorted(cur.items()))


def box_images(N: int, L: int, x0: int, x: int) -> int:
    """Ex. 4.3：無窮多個鏡像點 2kL + x0（正）與 2kL − x0（負）。"""
    tot, k = 0, 0
    K = N // (2 * L) + 2
    for k in range(-K, K + 1):
        tot += GF.count1(N, x - x0 - 2 * k * L) - GF.count1(N, x + x0 - 2 * k * L)
    return tot


def box_eigs(L: int) -> list[float]:
    """內點 1..L−1 的轉移矩陣（鄰居的和）的特徵值 2cos(mπ/L)，m = 1..L−1。"""
    return [2 * math.cos(m * PI / L) for m in range(1, L)]


def box_mode(L: int, m: int) -> list[float]:
    """正規化的正弦模態 √(2/L)·sin(mπx/L)，x = 1..L−1（(4.9)(4.10) 的格子版）。"""
    return [math.sqrt(2 / L) * math.sin(m * PI * x / L) for x in range(1, L)]


def box_eig_residual(L: int) -> float:
    """max |T v_m − λ_m v_m|（核對特徵向量；Ex. 4.4 的格子版）。"""
    worst = 0.0
    for m, lam in enumerate(box_eigs(L), start=1):
        v = [0.0] + box_mode(L, m) + [0.0]
        for x in range(1, L):
            worst = max(worst, abs(v[x - 1] + v[x + 1] - lam * v[x]))
    return worst


def box_gf(z: float, L: int) -> list[list[float]]:
    """兩牆之間的生成函數矩陣 G(z; x1, x2) = [(I − zT)^{−1}]，x1, x2 = 1..L−1（高斯消去；索引 0 ↔ x = 1）。"""
    n = L - 1
    A = [[(1.0 if i == j else 0.0) - (z if abs(i - j) == 1 else 0.0) for j in range(n)] + [1.0 if i == k else 0.0
                                                                                          for k in range(n)]
         for i in range(n)]
    for c in range(n):
        piv = A[c][c]
        A[c] = [v / piv for v in A[c]]
        for r in range(n):
            if r != c and A[r][c]:
                f = A[r][c]
                A[r] = [a - f * b for a, b in zip(A[r], A[c])]
    return [row[n:] for row in A]


def box_amplitude(z: float, L: int, m: int, x2: int) -> float:
    """(4.12)(4.13) 的 g(z; q)：把 G(z; x1, x2) 對 x1 投影到 √(2/L) sin(qx1)，q = mπ/L（數值）。"""
    G = box_gf(z, L)
    return sum(math.sqrt(2 / L) * math.sin(m * PI * x1 / L) * G[x1 - 1][x2 - 1] for x1 in range(1, L))


def box_amplitude_formula(z: float, L: int, m: int, x2: int) -> float:
    """更正後的 (4.13)：√(2/L) sin(q x2)/(1 − 2z cos q)（書上印成 sin q x1）。"""
    q = m * PI / L
    return math.sqrt(2 / L) * math.sin(q * x2) / (1 - 2 * z * math.cos(q))


def box_image_points(L: int, x0: int, K: int = 1) -> tuple:
    """Ex. 4.3 的鏡像點（k = −K..K）：(正的 2kL + x0 們, 負的 2kL − x0 們)。"""
    return [2 * k * L + x0 for k in range(-K, K + 1)], [2 * k * L - x0 for k in range(-K, K + 1)]


def box_modes(N: int, L: int, x0: int, x: int) -> float:
    """(4.15) 的格子精確版：Σ_m (2/L) sin(mπx0/L) sin(mπx/L) (2cos(mπ/L))^N。"""
    return sum(2 / L * math.sin(m * PI * x0 / L) * math.sin(m * PI * x / L) * (2 * math.cos(m * PI / L)) ** N
               for m in range(1, L))


def box_total(N: int, L: int) -> int:
    """從 1..L−1 任一點出發、終點不限、N 步都沒碰牆的走法總數（逐步遞迴，精確整數）。"""
    v = [0] + [1] * (L - 1) + [0]
    for _ in range(N):
        v = [0] + [v[x - 1] + v[x + 1] for x in range(1, L)] + [0]
    return sum(v)


def box_total_modes(N: int, L: int) -> float:
    """同上 ÷ 2^N，用格子模態：Σ_{m 奇} (2/L) cot²(mπ/2L) cos^N(mπ/L)（第二種算法）。"""
    return sum(2 / L / math.tan(m * PI / (2 * L)) ** 2 * math.cos(m * PI / L) ** N for m in range(1, L, 2))


def box_total_cont(N: float, L: float, mmax: int = 2000) -> float:
    """(4.19) ÷ 2^N：(8L/π²) Σ_m e^{−N((2m+1)π/L)²/2}/(2m+1)²（連續極限，Δ = 1）。"""
    tot = 0.0
    for m in range(mmax):
        n = 2 * m + 1
        tot += math.exp(-N * (n * PI / L) ** 2 / 2) / n ** 2
    tail = 1 / (4 * mmax) if N == 0 else 0.0               # N = 0 時補上 Σ 的尾巴
    return 8 * L / PI ** 2 * (tot + tail)


def box_long(N: float, L: float, book: bool = False) -> float:
    """(4.24) ÷ 2^N：(8L/π²)e^{−π²N/(2L²)}；book=True 用書上印的 e^{−π²N/L²}。"""
    return 8 * L / PI ** 2 * math.exp(-PI ** 2 * N / (L * L if book else 2 * L * L))


def box_survival(N: int, L: int) -> Fraction:
    """存活比例：box_total ÷ (2^N·(L−1))（起點在 L−1 個內點上平均）。"""
    return Fraction(box_total(N, L), 2 ** N * (L - 1))


def box_survival_long(N: float, L: int, kind: str = "lattice") -> float:
    """存活比例的長時間直線（對數圖用）：'lattice' 格子最低模態 (2/L)cot²(π/2L)cos^N(π/L)/(L−1)；
    'cont' (8/π²)e^{−π²N/(2L²)}（(4.24) ÷ 2^N L）；'book' (8/π²)e^{−π²N/L²}（書上印的指數）。"""
    if kind == "lattice":
        return 2 / L / math.tan(PI / (2 * L)) ** 2 * math.cos(PI / L) ** N / (L - 1)
    return 8 / PI ** 2 * math.exp(-decay_rate(L, kind) * N)


def decay_rate(L: float, kind: str = "lattice") -> float:
    """長時間每步的衰減率（存活比例 ∝ e^{−rate·N}）：'lattice' −ln cos(π/L)、'cont' π²/(2L²)、'book' π²/L²。"""
    if kind == "lattice":
        return -math.log(math.cos(PI / L))
    return PI ** 2 / (2 * L * L) if kind == "cont" else PI ** 2 / (L * L)


def odd_inv_squares(M: int) -> float:
    """Σ_{m<M} 1/(2m+1)² ＋ Euler–Maclaurin 尾巴（→ π²/8，(4.22)）。"""
    s = math.fsum(1 / (2 * m + 1) ** 2 for m in range(M))
    a = 2 * M + 1                                           # 尾巴 Σ_{n≥a, 奇} 1/n² ≈ 1/(2(a−1)) 的精細版
    return s + 1 / (2 * a) + 1 / (2 * a * a) + 1 / (3 * a ** 3)


def zeta2(M: int) -> float:
    """Σ_{n≥1} 1/n²（部分和＋Euler–Maclaurin 尾巴；(4.21)）。"""
    return math.fsum(1 / n ** 2 for n in range(1, M)) + 1 / M + 1 / (2 * M * M) + 1 / (6 * M ** 3)


# =========================================================================== §4.1.4 球與立方體


def j0(x: float) -> float:
    """(4.32) sin x / x。"""
    return 1.0 if x == 0 else math.sin(x) / x


def j0_zero(n: int) -> float:
    """j0 的第 n 個零點（二分法；應等於 nπ，(4.33)）。"""
    a, b = (n - 0.5) * PI, (n + 0.5) * PI
    fa = j0(a)
    for _ in range(200):
        c = (a + b) / 2
        if (j0(c) > 0) == (fa > 0):
            a, fa = c, j0(c)
        else:
            b = c
    return (a + b) / 2


def sphere_weight(n: int, R: float) -> float:
    """(4.34) 的係數 8R³/(πn²)：第 n 個 l = 0 模態對「走法總數」的權重。"""
    return 8 * R ** 3 / (PI * n * n)


def sphere_weight_numeric(n: int, R: float) -> float:
    """同上，數值積分：ψ = c·j0(kr)，c² = 1/∫|j0|²d³r；權重 (∫ψ d³r)²（第二種算法；更正後的 (4.31)）。"""
    k = n * PI / R
    norm = _simpson(lambda r: 4 * PI * r * r * j0(k * r) ** 2, 0.0, R, 2000)
    tot = _simpson(lambda r: 4 * PI * r * r * j0(k * r), 0.0, R, 2000)
    return tot ** 2 / norm


def sphere_weight_book(n: int, R: float) -> float:
    """書上 (4.31) 印的前因子 (8π/R)·(2/R)·(∫r²j0 dr)²（和 (4.34) 對不上）。"""
    k = n * PI / R
    rad = _simpson(lambda r: r * r * j0(k * r), 0.0, R, 2000)
    return 8 * PI / R * 2 / R * rad ** 2


def sphere_total_cont(N: float, R: float, book: bool = False, nmax: int = 4000) -> float:
    """(4.35) ÷ 6^N：(8R³/π) Σ e^{−π²n²N/(6R²)}/n²（更正）；book=True 用書上的 e^{−π²N/R²}（沒有 n²、沒有 1/6）。"""
    tot = 0.0
    for n in range(1, nmax + 1):
        e = PI ** 2 * N / R ** 2 if book else PI ** 2 * n * n * N / (6 * R * R)
        tot += math.exp(-e) / (n * n)
    if N == 0:
        tot += 1 / nmax
    return 8 * R ** 3 / PI * tot


def sphere_rate(R: float, book: bool = False) -> float:
    """球的長時間衰減率（每步）：π²/(6R²)（更正後的 (4.37)）；書上 π²/R²。"""
    return PI ** 2 / (R * R) if book else PI ** 2 / (6 * R * R)


def ball_sites(R: float) -> list[tuple]:
    """sc 格點裡 x²+y²+z² < R² 的點（踩到球外就被吸收）。"""
    n = int(R) + 1
    return [(x, y, z) for x in range(-n, n + 1) for y in range(-n, n + 1) for z in range(-n, n + 1)
            if x * x + y * y + z * z < R * R]


def _region_totals(sites: list[tuple], Nmax: int) -> list[float]:
    idx = {s: i for i, s in enumerate(sites)}
    nb = [[idx[t] for t in ((x + 1, y, z), (x - 1, y, z), (x, y + 1, z), (x, y - 1, z), (x, y, z + 1), (x, y, z - 1))
           if t in idx] for (x, y, z) in sites]
    v = [1.0] * len(sites)
    out = [float(len(sites))]
    for _ in range(Nmax):
        v = [sum(v[j] for j in js) / 6 for js in nb]
        out.append(sum(v))
    return out


def ball_totals(R: float, Nmax: int) -> list[float]:
    """格子球：起點在球內任一點、N 步都沒出球的走法數 ÷ 6^N，N = 0..Nmax。"""
    return _region_totals(ball_sites(R), Nmax)


def ball_rate(R: float, Nmax: int = 120) -> float:
    """格子球的長時間衰減率：−½ ln(T_{N+2}/T_N)（偶數 N，夠大；第二種算法，核對 1/6）。"""
    t = ball_totals(R, Nmax + 2)
    return -0.5 * math.log(t[Nmax + 2] / t[Nmax])


def cube_total(N: int, L: int) -> int:
    """Ex. 4.5：立方體 0 < x, y, z < L（邊長 L、牆在 0 與 L）裡 N 步都沒碰牆的走法總數（逐步遞迴，精確）。"""
    pts = [(x, y, z) for x in range(1, L) for y in range(1, L) for z in range(1, L)]
    v = {p: 1 for p in pts}
    for _ in range(N):
        v = {(x, y, z): sum(v.get(q, 0) for q in ((x + 1, y, z), (x - 1, y, z), (x, y + 1, z), (x, y - 1, z),
                                                     (x, y, z + 1), (x, y, z - 1)))
             for (x, y, z) in pts}
    return sum(v.values())


def cube_total_modes(N: int, L: int) -> float:
    """同上 ÷ 6^N，用格子模態：Σ_{三個奇數 m} Π (2/L)cot²(m_iπ/2L) · ((Σ cos(m_iπ/L))/3)^N。"""
    ms = list(range(1, L, 2))
    w = {m: 2 / L / math.tan(m * PI / (2 * L)) ** 2 for m in ms}
    c = {m: math.cos(m * PI / L) for m in ms}
    return sum(w[a] * w[b] * w[d] * ((c[a] + c[b] + c[d]) / 3) ** N for a in ms for b in ms for d in ms)


def cube_total_cont(N: float, L: float, mmax: int = 60) -> float:
    """Ex. 4.5 的連續極限 ÷ 6^N：(8L/π²)³ Σ_{奇 n1,n2,n3} e^{−π²N(n1²+n2²+n3²)/(6L²)}/(n1n2n3)²。"""
    f = [(math.exp(-PI ** 2 * N * n * n / (6 * L * L)) / (n * n)) for n in range(1, 2 * mmax, 2)]
    s = math.fsum(f) + (1 / (4 * mmax) if N == 0 else 0.0)
    return (8 * L / PI ** 2) ** 3 * s ** 3


def cube_long(N: float, L: float) -> float:
    """立方體長時間 ÷ 6^N：(512L³/π⁶) e^{−π²N/(2L²)}。"""
    return 512 * L ** 3 / PI ** 6 * math.exp(-PI ** 2 * N / (2 * L * L))


PREFACTOR = {"interval": 8 / PI ** 2, "cube": 512 / PI ** 6, "sphere": 6 / PI ** 2}   # 長時間前因子 ÷ 體積


def same_volume_rate_ratio() -> float:
    """同體積的立方體與球：衰減率之比 (π²/2L²)/(π²/6R²) = 3(R/L)² = 3(3/(4π))^{2/3}。"""
    return 3 * (3 / (4 * PI)) ** (2 / 3)


# =========================================================================== §4.2–4.2.1 穩定態、帳篷、賭徒破產


def laplacian_coeff(a: float = 1e-3) -> float:
    """(4.40) 的係數：[(1/6)Σ_i f(r+a ê_i) − f(r)] / (a²∇²f) → 1/6（書上印 1/3）。測試函數 f = e^{−(x²+2y²+3z²)/2}。"""
    def f(x, y, z):
        return math.exp(-(x * x + 2 * y * y + 3 * z * z) / 2)
    x, y, z = 0.3, -0.2, 0.1
    lap = f(x, y, z) * ((x * x - 1) + (4 * y * y - 2) + (9 * z * z - 3))
    nb = sum(f(x + dx, y + dy, z + dz) for dx, dy, dz in
             ((a, 0, 0), (-a, 0, 0), (0, a, 0), (0, -a, 0), (0, 0, a), (0, 0, -a)))
    return (nb / 6 - f(x, y, z)) / (a * a * lap)


def _solve_tridiag(L: int, rhs: list) -> list:
    """解 w(x) − (w(x−1)+w(x+1))/2 = rhs(x)，x = 1..L−1，w(0) = w(L) = 0（Fraction，Thomas 演算法）。"""
    n = L - 1
    a = [Fraction(-1, 2)] * n          # 下對角
    b = [Fraction(1)] * n
    c = [Fraction(-1, 2)] * n          # 上對角
    d = [Fraction(v) for v in rhs]
    for i in range(1, n):
        m = a[i] / b[i - 1]
        b[i] -= m * c[i - 1]
        d[i] -= m * d[i - 1]
    w = [Fraction(0)] * n
    w[-1] = d[-1] / b[-1]
    for i in range(n - 2, -1, -1):
        w[i] = (d[i] - c[i] * w[i + 1]) / b[i]
    return [Fraction(0)] + w + [Fraction(0)]


def tent_lattice(L: int, x0: int) -> list[Fraction]:
    """一維格子的穩定態：每步在 x0 放 1 個走路人，(4.39) 的一維版 w(x) = ½(w(x−1)+w(x+1)) + δ_{x,x0}，牆在 0、L。w(0..L)。"""
    return _solve_tridiag(L, [1 if x == x0 else 0 for x in range(1, L)])


def tent_formula(L: int, x0: int, x, s0=1, beta=Fraction(1, 2)):
    """(4.49)：(s0/β)·x(L−x0)/L（x < x0）、(s0/β)·x0(L−x)/L（x > x0）。格子上 β = a²/(2τ) = ½。"""
    if x <= x0:
        return Fraction(s0) / beta * x * (L - x0) / L
    return Fraction(s0) / beta * x0 * (L - x) / L


def absorb_left(L: int, x0: int) -> Fraction:
    """Ex. 4.6：從 x0 出發、先被左牆 0 吸收的機率 (L − x0)/L（賭徒破產）。"""
    return Fraction(L - x0, L)


def absorb_left_solve(L: int, x0: int) -> Fraction:
    """同上，解調和方程 h(x) = ½(h(x−1)+h(x+1))，h(0) = 1、h(L) = 0（第二種算法；Fraction）。"""
    # h = 1 − x/L 的齊次解；用 _solve_tridiag 解 u = h − (1 − x/L)·0 不方便，改用直接消去
    n = L - 1
    b = [Fraction(1)] * n
    d = [Fraction(0)] * n
    d[0] = Fraction(1, 2)                  # 左邊界 h(0) = 1 移到右邊
    c = Fraction(-1, 2)
    for i in range(1, n):
        m = c / b[i - 1]
        b[i] -= m * c
        d[i] -= m * d[i - 1]
    h = [Fraction(0)] * n
    h[-1] = d[-1] / b[-1]
    for i in range(n - 2, -1, -1):
        h[i] = (d[i] - c * h[i + 1]) / b[i]
    return h[x0 - 1]


def absorb_left_dp(L: int, x0: int, tol: float = 1e-14) -> float:
    """同上，逐步遞迴到幾乎全被吸收（第三種算法；float）。"""
    p = [0.0] * (L + 1)
    p[x0] = 1.0
    left = 0.0
    while sum(p) > tol:
        q = [0.0] * (L + 1)
        for x in range(1, L):
            q[x - 1] += p[x] / 2
            q[x + 1] += p[x] / 2
        left += q[0]
        q[0] = q[L] = 0.0
        p = q
    return left


def wall_fluxes(L: int, x0: int) -> tuple:
    """穩定態流進左、右牆的走路人數（每步）：w(1)/2、w(L−1)/2（Fraction；和 = 1）。"""
    w = tent_lattice(L, x0)
    return w[1] / 2, w[L - 1] / 2


# =========================================================================== §4.2.2 點源：把高斯對時間積分；sc 格子 Green 函數


def heat_kernel(r: float, t: float, D: float = 1.0) -> float:
    """三維高斯 (4πDt)^{−3/2} e^{−r²/(4Dt)}（(4.54) 的正規化版本）。"""
    return (4 * PI * D * t) ** -1.5 * math.exp(-r * r / (4 * D * t))


def time_integral(r: float, D: float = 1.0) -> float:
    """∫_0^∞ heat_kernel dt（數值；t = e^v）→ 1/(4πDr)（(4.52)(4.55)(4.56)）。"""
    return _simpson(lambda v: math.exp(v) * heat_kernel(r, math.exp(v), D), math.log(r * r / D) - 12,
                    math.log(r * r / D) + 40, 6000)


def K_const(a: float) -> float:
    """(4.56) 的 K = ∫_0^∞ T^{−3/2} e^{−1/(aT)} dT = √(πa)。"""
    return math.sqrt(PI * a)


def K_const_numeric(a: float) -> float:
    return _simpson(lambda v: math.exp(-0.5 * v) * math.exp(-1 / (a * math.exp(v))), math.log(1 / a) - 8,
                    math.log(1 / a) + 60, 6000)


def _ine_small(n: int, y: float) -> float:
    s, t, k, q = 0.0, (y / 2) ** n / math.factorial(n), 0, y * y / 4
    if t == 0:
        return 0.0
    while True:
        s += t
        k += 1
        t *= q / (k * (k + n))
        if t < 1e-18 * s:
            return (s + t) * math.exp(-y)


def ine_list(nmax: int, y: float) -> list[float]:
    """[e^{−y} I_n(y)]，n = 0..nmax（y < 1 用級數；否則 Miller 逆向遞迴，用 e^{y} = I_0 + 2Σ I_k 正規化）。"""
    if y < 1:
        return [_ine_small(n, y) for n in range(nmax + 1)]
    M = nmax + int(math.sqrt(80 * y)) + 20
    b = [0.0] * (M + 2)
    b[M] = 1e-30
    for k in range(M, 0, -1):
        b[k - 1] = b[k + 1] + 2 * k / y * b[k]
        if b[k - 1] > 1e250:
            b = [v * 1e-250 for v in b]
    norm = b[0] + 2 * sum(b[1:M + 1])
    return [b[n] / norm for n in range(nmax + 1)]


def _hankel(n: int, K: int) -> list[float]:
    """e^{−y}I_n(y)·√(2πy) ≈ Σ_k c_k y^{−k}（c_k = (−1)^k Π_{j≤k}(4n²−(2j−1)²)/(k!8^k)）。"""
    c, out = 1.0, [1.0]
    for k in range(1, K + 1):
        c *= -(4 * n * n - (2 * k - 1) ** 2) / (k * 8)
        out.append(c)
    return out


_SC_Y, _SC_PTS = 1000.0, 6000


@lru_cache(maxsize=None)
def _sc_table(nmax: int):
    V = math.log1p(_SC_Y)
    h = V / _SC_PTS
    rows = []
    for i in range(_SC_PTS + 1):
        v = i * h
        y = math.expm1(v)
        w = 1 if i in (0, _SC_PTS) else (4 if i % 2 else 2)
        rows.append((w * h / 3 * math.exp(v), ine_list(nmax, y)))
    return rows


@lru_cache(maxsize=None)
def sc_green(n1: int, n2: int, n3: int) -> float:
    """sc 格子的 Green 函數 G(z_c; r)：從原點出發、到訪 r = (n1, n2, n3) 的平均次數
    = ∫_0^∞ e^{−t} I_{n1}(t/3) I_{n2}(t/3) I_{n3}(t/3) dt（數值積分到 t = 3000＋漸近尾巴）。
    也是 (4.39) 在原點每步放 1 個走路人時的穩定態 w(r)。"""
    n1, n2, n3 = sorted((abs(n1), abs(n2), abs(n3)))
    nmax = max(12, n3)
    val = 3 * sum(wt * row[n1] * row[n2] * row[n3] for wt, row in _sc_table(nmax))
    K = 8
    a, b, c = _hankel(n1, K), _hankel(n2, K), _hankel(n3, K)
    prod = R._ser_mul(R._ser_mul(a, b, K), c, K)
    val += 3 * (2 * PI) ** -1.5 * sum(bk * _SC_Y ** (-k - 0.5) / (k + 0.5) for k, bk in enumerate(prod))
    return val


def sc_green_asym(r: float) -> float:
    """遠處 G ≈ 3/(2πr)：更正後的 (4.41) −∇²w = (6/a²)S 的解 6/(4πr)·(a = 1)。"""
    return 3 / (2 * PI * r)


def sc_green_book(r: float) -> float:
    """照書上 (4.41) 的 3/a² 會得到 3/(4πr)（差 2 倍）。"""
    return 3 / (4 * PI * r)


def sc_laplace_residual(n1: int, n2: int, n3: int) -> float:
    """離散 (4.39)：G(r) − (1/6)Σ G(鄰居) − δ_{r,0}（應為 0；核對 sc_green 的數值）。"""
    nb = sum(sc_green(n1 + dx, n2 + dy, n3 + dz) for dx, dy, dz in
             ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)))
    return sc_green(n1, n2, n3) - nb / 6 - (1 if (n1, n2, n3) == (0, 0, 0) else 0)


# =========================================================================== §4.2.3 吸收牆旁的點源：鏡像源、偶極


def conc_wall(x: float, y: float, z: float, l: float) -> float:
    """(4.57)(4.58)：點源在 (l,0,0)、鏡像源在 (−l,0,0)、吸收牆 x = 0；β = s0 = 1：(1/4π)(1/|r−r_s| − 1/|r−r_i|)。"""
    return (1 / math.sqrt((x - l) ** 2 + y * y + z * z) - 1 / math.sqrt((x + l) ** 2 + y * y + z * z)) / (4 * PI)


def conc_dipole(x: float, y: float, z: float, l: float) -> float:
    """(4.58) 末行（r ≫ l）：(1/4π)·2lx/r³。"""
    return 2 * l * x / (x * x + y * y + z * z) ** 1.5 / (4 * PI)


def plane_flux(X: float, l: float) -> float:
    """穿過平面 x = X（往 +x）的總流量 ∫∫ j_x dy dz（j = −∇c，精確的 (4.58) 第一行；數值積分，ρ = X·tan u）。"""
    def jx(rho):
        s = sum(sg * (X - x0) / ((X - x0) ** 2 + rho * rho) ** 1.5 for sg, x0 in ((1, l), (-1, -l)))
        return s / (4 * PI)
    return _simpson(lambda u: jx(abs(X) * math.tan(u)) * 2 * PI * abs(X) * math.tan(u) * abs(X) / math.cos(u) ** 2,
                    0.0, PI / 2 - 1e-9, 4000)


def plane_flux_exact(X: float, l: float) -> float:
    """同上的閉式：每個點源穿過無限大平面的流量是 ±½（立體角 2π）。X > l：½ − ½ = 0；0 < X < l：−½ − ½ = −1（全部流向牆）。"""
    def half(x0):
        return 0.5 if X > x0 else -0.5
    return half(l) - half(-l)


def dipole_flux_integral(x: float) -> float:
    """(4.60)：∫_0^∞ [3x²/(x²+R²)^{5/2} − 1/(x²+R²)^{3/2}] R dR（數值；應為 0）。"""
    return _simpson(lambda u: (3 * x * x / (x * x + (x * math.tan(u)) ** 2) ** 2.5
                               - 1 / (x * x + (x * math.tan(u)) ** 2) ** 1.5) * x * math.tan(u) * x / math.cos(u) ** 2,
                    0.0, PI / 2 - 1e-7, 8000)


# =========================================================================== §4.2.4–4.2.5 吸收球、電容


def conc_sphere(r: float, r0: float, cinf: float = 1.0) -> float:
    """(4.61) c∞(1 − r0/r)。"""
    return cinf * (1 - r0 / r)


def sphere_flux(r0: float, D: float = 1.0, cinf: float = 1.0) -> float:
    """(4.63)/(4.65)(4.66)：吸收球每單位時間吃掉 4πDr0c∞（D 就是書上的 β）。"""
    return 4 * PI * D * r0 * cinf


def sphere_flux_numeric(r: float, r0: float, D: float = 1.0, cinf: float = 1.0) -> float:
    """在半徑 r 的球面上算 4πr²·D·dc/dr（數值微分；任何 r 都相同）。"""
    h = 1e-5 * r
    return 4 * PI * r * r * D * (conc_sphere(r + h, r0, cinf) - conc_sphere(r - h, r0, cinf)) / (2 * h)


def kelvin_image(R_: float, r0: float) -> tuple:
    """吸收球外的點源（距球心 R）：鏡像源在 ρ = r0²/R、強度 −r0/R（Fig. 4.6）。"""
    return r0 * r0 / R_, -r0 / R_


def hit_prob(R_: float, r0: float) -> float:
    """從距球心 R 處出發，終究撞到半徑 r0 吸收球的機率 r0/R（逃脫 (R − r0)/R）。"""
    return r0 / R_


def _source_image_c(r: float, th: float, R_: float, r0: float) -> float:
    rho, q = kelvin_image(R_, r0)
    d1 = math.sqrt(r * r + R_ * R_ - 2 * r * R_ * math.cos(th))
    d2 = math.sqrt(r * r + rho * rho - 2 * r * rho * math.cos(th))
    return (1 / d1 + q / d2) / (4 * PI)


def sphere_surface_max(R_: float, r0: float) -> float:
    """點源＋鏡像源在球面上的濃度最大絕對值（應為 0：球面是吸收邊界）。"""
    return max(abs(_source_image_c(r0, PI * i / 200, R_, r0)) for i in range(201))


def hit_prob_numeric(R_: float, r0: float) -> float:
    """點源（每單位時間 1 個）＋鏡像源：球面上的流入量 ∮ ∂c/∂r dS（數值積分；第二種算法，→ r0/R）。"""
    h = 1e-6 * r0

    def flux(th):
        dcdr = (_source_image_c(r0 + h, th, R_, r0) - _source_image_c(r0 - h, th, R_, r0)) / (2 * h)
        return dcdr * 2 * PI * r0 * r0 * math.sin(th)
    return _simpson(flux, 0.0, PI, 4000)


def cell_rmax(D: float, cinf: float, q: float) -> float:
    """營養流入 4πDr0c∞ ＝ 需求 (4/3)πr0³q → 最大半徑 √(3Dc∞/q)。"""
    return math.sqrt(3 * D * cinf / q)


CELL_EXAMPLE = {"D": 2e-9, "cinf": 0.25, "q": 1.0}     # 假設值（只看數量級）：氧在水中 m²/s、空氣飽和的水 mol/m³、耗氧率 mol/(m³·s)


def capacitance_sphere(r0: float) -> float:
    """(4.66) C = r0（高斯單位）。"""
    return r0


def disk_capacitance(a: float) -> float:
    """（書外）半徑 a 的導體圓盤 C = 2a/π；嵌在不吸收表面上的受器只有一面，等效 a/π，流量 4Dac∞。"""
    return 2 * a / PI


def disk_capacitance_numeric(a: float) -> float:
    """同上，數值：面電荷 σ ∝ 1/√(a²−ρ²)，總電荷 Q 與圓心電位 φ(0) 都用數值積分（ρ = a sin u），C = Q/φ(0)。"""
    Q = _simpson(lambda u: 2 * PI * a * math.sin(u), 0.0, PI / 2, 2000)            # ∫ σ·2πρ dρ，σ = 1/√(a²−ρ²)
    phi0 = _simpson(lambda u: 2 * PI, 0.0, PI / 2, 2000)                            # ∫ σ·2πρ/ρ dρ
    return Q / phi0


def disk_potential_offcenter(a: float, p: float) -> float:
    """同一個電荷分布（σ = 1/√(a²−ρ²)）在圓盤上距圓心 p 處的電位（數值；等於 φ(0) = π²，導體是等電位）。"""
    def ring(rho):                                                                   # ∫_0^{2π} dθ/|p − ρe^{iθ}| = 4K(k)/(ρ+p)
        a_, g_ = 1.0, abs(rho - p) / (rho + p)                                       # 用 k' = |ρ−p|/(ρ+p) 算 AGM（避免 1−k² 相消）
        while abs(a_ - g_) > 1e-15 * a_:
            a_, g_ = (a_ + g_) / 2, math.sqrt(a_ * g_)
        return 4 * (PI / (2 * a_)) / (rho + p)
    f1 = _simpson(lambda s: ring(p * (1 - s * s)) * p * (1 - s * s) / math.sqrt(a * a - (p * (1 - s * s)) ** 2)
                  * 2 * p * s, 1e-6, 1.0, 4000)
    def g2(u):                                                                      # ρ = p + (a−p)sin²u，端點的 cos u 已約掉
        rho = p + (a - p) * math.sin(u) ** 2
        return ring(rho) * rho * 2 * math.sin(u) * math.sqrt(a - p) / math.sqrt(a + rho)
    f2 = _simpson(g2, 1e-6, PI / 2, 4000)
    return f1 + f2


# =========================================================================== §4.2.6 布滿受器的球（Berg & Purcell 1977，書外補充）


def berg_purcell(n: float, a: float, r0: float) -> float:
    """J/J_max = na/(na + πr0)：n 個半徑 a 的圓盤受器（Berg & Purcell 1977）。"""
    return n * a / (n * a + PI * r0)


def berg_purcell_series(n: float, a: float, r0: float) -> float:
    """同上，用「電容串聯」：1/C_eff = 1/r0 + 1/(n·a/π)，J/J_max = C_eff/r0（第二種算法）。"""
    c_eff = 1 / (1 / capacitance_sphere(r0) + 1 / (n * disk_capacitance(a) / 2))
    return c_eff / r0


def bp_half(a: float, r0: float) -> float:
    """半飽和的受器數 n = πr0/a。"""
    return PI * r0 / a


def bp_area_fraction(n: float, a: float, r0: float) -> float:
    """受器占球面的比例 n·πa²/(4πr0²) = na²/(4r0²)。"""
    return n * a * a / (4 * r0 * r0)


BP_EXAMPLE = {"r0": 1e-6, "a": 1e-9}                  # 細菌半徑 1 µm、受器半徑 1 nm（Berg & Purcell 的典型值）


# =========================================================================== §4.3 補充：吸收邊界 G(z; x, −Δ) = 0


def half_counts(N: int, x: int, y: int) -> int:
    """書上 §4.3 的約定：格點 n ≥ 0，牆緊鄰在 −1 的左邊；從 x 走 N 步到 y、從沒踩到 −1 的走法數（逐步遞迴）。"""
    cur = {x: 1}
    for _ in range(N):
        nxt = {}
        for p, c in cur.items():
            for q in (p - 1, p + 1):
                if q >= 0:
                    nxt[q] = nxt.get(q, 0) + c
        cur = nxt
    return cur.get(y, 0)


def half_counts_image(N: int, x: int, y: int) -> int:
    """同上，鏡像點在 −x − 2（對 −1 鏡射）：C(N; y−x) − C(N; y+x+2)。"""
    return GF.count1(N, y - x) - GF.count1(N, y + x + 2)


def half_gf(z: float, x: int, y: int) -> float:
    """Ex. 4.1 的格子版（同一約定）：G'(z; x, y) = (λ^{|x−y|} − λ^{x+y+2})/√(1−4z²)，λ = (1−√(1−4z²))/(2z)。"""
    return GF.G1_exact(z, x - y) - GF.G1_exact(z, x + y + 2)


def half_gf_residual(z: float, x: int, ymax: int = 6) -> float:
    """(4:S-1)–(4:S-3)：在 y ≥ 0 都滿足 G(x,y) = z(G(x,y−1)+G(x,y+1)) + δ_{x,y}，且 G(x,−1) = 0（回傳最大誤差）。"""
    worst = abs(half_gf(z, x, -1))
    for y in range(0, ymax + 1):
        worst = max(worst, abs(half_gf(z, x, y) - z * (half_gf(z, x, y - 1) + half_gf(z, x, y + 1)) - (x == y)))
    return worst


def catalan(k: int) -> int:
    return math.comb(2 * k, k) // (k + 1)


def half_visits_critical(x: int, y: int) -> int:
    """z → z_c = ½：從 x 出發、到訪 y 的平均次數 2(min(x, y) + 1)（帳篷在 L → ∞ 的極限）。"""
    return 2 * (min(x, y) + 1)


# =========================================================================== check()


def _f(x: float, d: int) -> str:
    """場景卡上的寫法（固定小數位）。"""
    return f"{x:.{d}f}"


def check():
    F = Fraction
    # ---- §4.1.1 反射原理：暴力列舉＝鏡像公式 (4.1)＝逐步遞迴；鏡射是一一對應
    for x0 in (1, 2, 3):
        for N in range(0, 11):
            s = survivors_brute(N, x0)
            assert s == survivors_image(N, x0) == survivors_dp(N, x0), (x0, N)
            assert sum(s.values()) == survival_prob(N, x0) * 2 ** N
    for x0 in (1, 3):
        for N in range(1, 11):
            pr = reflection_pairs(N, x0)
            imgs = [tuple(q) for _, q in pr]
            assert all(q[0] == -x0 and q[-1] == p[-1] and touches_wall(q) for p, q in pr)
            assert len(set(imgs)) == len(imgs) == sum(1 for q in walks(N, -x0) if q[-1] > 0)
    W = walks(4, 1)                                                                   # W401 小例子
    assert len(W) == 16 and sum(map(touches_wall, W)) == 10 and survivors_brute(4, 1) == {1: 2, 3: 3, 5: 1}
    assert len(reflection_pairs(4, 1)) == 5 and ([1, 0, 1, 2, 3], [-1, 0, 1, 2, 3]) in reflection_pairs(4, 1)
    assert [GF.count1(4, x - 1) for x in (1, 3, 5)] == [6, 4, 1] and [GF.count1(4, x + 1) for x in (1, 3, 5)] == [4, 1, 0]
    assert sum(1 for p in W if p[-1] > 0) == 11 and 11 - 5 == 6
    assert survival_prob(10, 3) == F(21, 32) == F(672, 1024) and sum(survivors_image(10, 3).values()) == 672
    assert [_f(float(survival_prob(N, 3)), d) for N, d in ((10, 5), (100, 4), (1000, 5))] == ["0.65625", "0.2356", "0.07557"]
    assert [_f(survival_asym(N, 3), d) for N, d in ((10, 4), (100, 4), (1000, 5))] == ["0.7569", "0.2394", "0.07569"]
    assert abs(float(survival_prob(1000, 3)) / survival_asym(1000, 3) - 1) < 2e-3                   # → x0√(2/πN)
    assert abs(survival_3d(100, 3) - survival_3d_dp(100, 3)) < 1e-12 and abs(survival_3d(300, 3) - survival_3d_dp(300, 3)) < 1e-12
    assert abs(survival_3d(1000, 3) - survival_3d_dp(1000, 3)) < 1e-12
    assert [_f(survival_3d(N, 3), 4) for N in (100, 1000)] == ["0.3958", "0.1305"]
    assert [_f(survival_3d_asym(N, 3), 4) for N in (100, 1000)] == ["0.4146", "0.1311"]
    assert all(wall_reduction_1d(N) == wall_reduction_brute(N) for N in range(1, 13))           # (4.6)–(4.8) 一維版
    assert wall_reduction_1d(10) == F(1067, 512)
    assert [_f(float(wall_reduction_1d(N)), 3) for N in (10, 100)] == ["2.084", "7.499"]
    assert [_f(wall_reduction_asym(N), 3) for N in (10, 100)] == ["2.523", "7.979"]
    assert abs(wall_reduction_asym(300, 3) - wall_reduction_asym(100)) < 1e-12                    # √(2N/3π)：x 方向變異數 N/3
    # ---- §4.1.2 兩牆之間：逐步遞迴＝格子模態＝無窮多個鏡像（Ex. 4.3）；(4.12)(4.13) 的更正
    for L in (4, 6, 9):
        assert box_eig_residual(L) < 1e-12
        for x0 in range(1, L):
            for N in range(0, 16):
                d = box_dp(N, L, x0)
                for x in range(1, L):
                    c = d.get(x, 0)
                    assert c == box_images(N, L, x0, x) == round(box_modes(N, L, x0, x)), (L, x0, N, x)
        for N in range(0, 30):
            t = box_total(N, L)
            assert t == sum(sum(box_dp(N, L, x0).values()) for x0 in range(1, L))
            assert abs(t / 2 ** N - box_total_modes(N, L)) < 1e-9 * max(1.0, t / 2 ** N)
    assert [_f(v, 4) for v in box_eigs(6)] == ["1.7321", "1.0000", "0.0000", "-1.0000", "-1.7321"]
    assert [box_dp(N, 6, 3) for N in range(5)] == [{3: 1}, {2: 1, 4: 1}, {1: 1, 3: 2, 5: 1}, {2: 3, 4: 3}, {1: 3, 3: 6, 5: 3}]
    assert box_dp(10, 6, 3) == {1: 81, 3: 162, 5: 81} and 81 + 162 + 81 == 324 and 3 ** 5 == 243
    assert F(2, 6) * (243 + 243) == 162 and F(2, 6) * (F(1, 2) * 243 + F(1, 2) * 243) == 81
    assert box_image_points(6, 3) == ([-9, 3, 15], [-15, -3, 9])
    assert [box_total(N, 6) for N in range(11)] == [5, 8, 14, 24, 42, 72, 126, 216, 378, 648, 1134]
    for z in (0.2, 0.45):                                                             # 筆誤 (4.12)(4.13)：sin q x1 → sin q x2
        for m in (1, 2, 3):
            for x2 in (1, 2, 4):
                assert abs(box_amplitude(z, 6, m, x2) - box_amplitude_formula(z, 6, m, x2)) < 1e-12
    assert abs(box_amplitude(0.3, 6, 1, 2) - math.sqrt(2 / 6) * math.sin(PI * 1 / 6) / (1 - 0.6 * math.cos(PI / 6))) > 0.1
    assert math.sin(0) == 0                                                           # (4.15) 的 n = 0 項本來就是 0
    # ---- §4.1.3 長短時間；Σ_奇 1/m² = π²/8；(4.22)(4.24) 的更正
    assert abs(odd_inv_squares(10 ** 5) - PI ** 2 / 8) < 1e-12 and abs(0.75 * zeta2(10 ** 5) - PI ** 2 / 8) < 1e-12
    assert abs(zeta2(10 ** 5) - PI ** 2 / 6) < 1e-12                                 # (4.21)
    assert abs(math.fsum(1 / (2 * k) ** 2 for k in range(1, 10 ** 5)) + 1 / (4 * 10 ** 5) - PI ** 2 / 24) < 1e-9  # (4.22) 偶數項 Σ_{k≥1}
    assert [_f(math.fsum(1 / (2 * m + 1) ** 2 for m in range(k)), 4) for k in (1, 2, 3, 4)] == ["1.0000", "1.1111", "1.1511", "1.1715"]
    assert _f(PI ** 2 / 8, 4) == "1.2337"
    for L in (10, 20, 40):
        assert abs(box_total_cont(0, L) / L - 1) < 1e-4                               # (4.23)：短時間 = 2^N L
        assert abs(box_total_modes(0, L) - (L - 1)) < 1e-9                            # 格子：L−1 個內點
        assert abs(box_total_cont(3 * L * L, L) / box_long(3 * L * L, L) - 1) < 1e-12  # 長時間只剩第一項
        assert abs(decay_rate(L) / decay_rate(L, "cont") - 1) < 0.02 and decay_rate(L, "book") == 2 * decay_rate(L, "cont")
        lead = sum(2 / L / math.tan(m * PI / (2 * L)) ** 2 * math.cos(m * PI / L) ** (40 * L) for m in (1, L - 1))
        assert abs(box_total_modes(40 * L, L) / lead - 1) < 1e-9                     # m = L−1：二部格子的負特徵值 −2cos(π/L)
    assert [_f(decay_rate(10, k), d) for k, d in (("lattice", 5), ("cont", 5), ("book", 4))] == ["0.05018", "0.04935", "0.0987"]
    assert _f(decay_rate(20), 5) == "0.01239" and _f(decay_rate(10) / decay_rate(20), 2) == "4.05"
    assert [_f(decay_rate(20, "cont"), 5), _f(decay_rate(20, "book"), 4)] == ["0.01234", "0.0247"]
    assert decay_rate(20, "cont") * 4 == decay_rate(10, "cont")
    S10 = {N: box_survival(N, 10) for N in (0, 50, 100, 200)}
    assert S10[0] == 1 and [_f(float(S10[N]), 6) for N in (50, 100)] == ["0.072103", "0.005865"]
    assert f"{float(S10[200]):.3e}" == "3.881e-05"
    r = -math.log(float(box_survival(200, 10) / box_survival(180, 10))) / 20               # 對數圖的斜率＝ln cos(π/10)
    assert abs(r - decay_rate(10)) < 1e-6 and abs(float(S10[200]) / box_survival_long(200, 10) - 1) < 1e-3
    assert [_f(box_survival_long(100, 10, k), 6) for k in ("lattice", "cont")] == ["0.005861", "0.005830"]
    assert f"{box_survival_long(100, 10, 'book'):.2e}" == "4.19e-05" and _f(8 / PI ** 2, 4) == "0.8106"
    assert abs(r - decay_rate(10)) < 1e-6
    assert _f(2 * 100 / PI ** 2, 2) == "20.26"                                                      # e 倍的步數 2L²/π²
    assert _f(PI ** 2 / 2, 3) == "4.935"                                                            # (4.27) 的 ακ
    # ---- §4.1.4 球（j0 零點、權重、(4.31)(4.35)(4.37) 的更正）與立方體（Ex. 4.5）
    assert all(abs(j0_zero(n) - n * PI) < 1e-12 for n in (1, 2, 3, 4))
    for n in (1, 2, 3):
        for Rr in (1.0, 2.5):
            assert abs(sphere_weight_numeric(n, Rr) / sphere_weight(n, Rr) - 1) < 1e-9
    assert [_f(sphere_weight(n, 1), 3) for n in (1, 2)] == ["2.546", "0.637"] and _f(sphere_weight_book(1, 1), 3) == "0.516"
    assert abs(sphere_total_cont(0, 1) - 4 * PI / 3) < 1e-6 and _f(4 * PI / 3, 3) == "4.189"            # (4.36)
    assert abs(sphere_total_cont(30, 1) / (8 / PI * math.exp(-PI ** 2 * 30 / 6)) - 1) < 1e-12       # (4.37) 更正
    rb = ball_rate(8)
    assert abs(rb / sphere_rate(8) - 1) < 0.1 and sphere_rate(8, True) / rb > 6                         # 格子球核對 1/6
    assert abs(ball_rate(6) / sphere_rate(6) - 1) < 0.1 and len(ball_sites(8)) == 2103
    assert [_f(rb, 4), _f(sphere_rate(8), 4), _f(sphere_rate(8, True), 3)] == ["0.0242", "0.0257", "0.154"]
    for N in range(0, 13):
        assert abs(cube_total(N, 6) / 6 ** N - cube_total_modes(N, 6)) < 1e-9 * cube_total_modes(N, 6)
    assert cube_total(0, 6) == 125 and cube_total(2, 6) == 2970 and cube_total(10, 6) > 2 ** 10 * 125   # 三維是 6^N 不是 2^N
    assert abs(cube_total_cont(0, 6) - 216) < 0.01 and abs(cube_total_cont(400, 6) / cube_long(400, 6) - 1) < 1e-12
    assert abs(cube_total_modes(400, 6) / cube_total_modes(398, 6) - math.cos(PI / 6) ** 2) < 1e-12     # 立方體 = 一維的衰減率
    assert [_f(PREFACTOR[k], 3) for k in ("interval", "cube", "sphere")] == ["0.811", "0.533", "0.608"]
    assert _f(same_volume_rate_ratio(), 3) == "1.155"
    Rv = (3 / (4 * PI)) ** (1 / 3)                                                    # 同體積：L = 1 的立方體、半徑 Rv 的球
    assert abs(decay_rate(1, "cont") / sphere_rate(Rv) - same_volume_rate_ratio()) < 1e-12
    # ---- §4.2–4.2.1 穩定態：(4.40)(4.41) 的 1/6；帳篷；賭徒破產（三種算法）
    assert abs(laplacian_coeff() - 1 / 6) < 1e-6 and _f(laplacian_coeff(), 4) == "0.1667"
    for L in (3, 5, 8, 11):
        for x0 in range(1, L):
            w = tent_lattice(L, x0)
            assert all(w[x] == tent_formula(L, x0, x) for x in range(L + 1)), (L, x0)
            assert w[x0] - (w[x0 - 1] + w[x0 + 1]) / 2 == 1
            assert (w[x0 + 1] - w[x0]) - (w[x0] - w[x0 - 1]) == -2                    # (4.48)：斜率跳 −s0/β = −2
            assert absorb_left(L, x0) == absorb_left_solve(L, x0) == wall_fluxes(L, x0)[0]
            assert abs(absorb_left_dp(L, x0) - float(absorb_left(L, x0))) < 1e-12
            assert sum(wall_fluxes(L, x0)) == 1
    assert tent_lattice(5, 2) == [0, F(6, 5), F(12, 5), F(8, 5), F(4, 5), 0]
    assert wall_fluxes(5, 2) == (F(3, 5), F(2, 5)) and absorb_left(5, 2) == F(3, 5)
    assert [float(v) for v in tent_lattice(5, 2)] == [0, 1.2, 2.4, 1.6, 0.8, 0]
    # ---- §4.2.2 點源：時間積分 1/(4πDr)；sc 的 Green 函數 × r → 3/(2π)
    for r_, D_ in ((2.0, 1.0), (0.7, 2.5)):
        assert abs(time_integral(r_, D_) * 4 * PI * D_ * r_ - 1) < 1e-8
    assert _f(time_integral(2.0), 5) == "0.03979" and _f(1 / (8 * PI), 5) == "0.03979"
    for a_ in (0.5, 4.0):
        assert abs(K_const_numeric(a_) / K_const(a_) - 1) < 1e-9
    G = [sc_green(n, 0, 0) for n in range(11)]
    assert abs(G[0] - R.G_C["sc"]) < 1e-10 and abs(G[1] - (G[0] - 1)) < 1e-10                       # 原點的方程：G(0) = 1 + G(1)
    for p in ((0, 0, 0), (1, 0, 0), (2, 0, 0), (3, 0, 0), (1, 1, 0), (2, 1, 0), (1, 1, 1), (5, 0, 0)):
        assert abs(sc_laplace_residual(*p)) < 1e-10, p
    assert [_f(G[n], 4) for n in range(6)] == ["1.5164", "0.5164", "0.2573", "0.1653", "0.1217", "0.0966"]
    assert [_f(n * G[n], 4) for n in (1, 2, 3, 4, 5, 10)] == ["0.5164", "0.5147", "0.4958", "0.4869", "0.4830", "0.4787"]
    assert abs(10 * G[10] / (10 * sc_green_asym(10)) - 1) < 3e-3 and abs(10 * G[10] / (10 * sc_green_book(10)) - 2) < 1e-2
    assert [_f(3 / (2 * PI), 4), _f(3 / (4 * PI), 4)] == ["0.4775", "0.2387"]
    # ---- §4.2.3 吸收牆旁的點源：牆在 x = 0（y–z 平面）；偶極；穿過遠處平面的流量 0
    for y_, z_ in ((0.0, 0.0), (1.3, -0.4), (5.0, 2.0)):
        assert conc_wall(0.0, y_, z_, 1.0) == 0.0
    assert abs(conc_wall(1.0, 0.0, 0.5, 1.0)) > 0.05                                 # x–z 平面（y = 0）上不是 0
    assert abs(conc_wall(30.0, 4.0, 3.0, 1.0) / conc_dipole(30.0, 4.0, 3.0, 1.0) - 1) < 2e-3
    for X in (1.5, 3.0, 10.0):
        assert abs(plane_flux(X, 1.0) - plane_flux_exact(X, 1.0)) < 1e-6 and plane_flux_exact(X, 1.0) == 0
    assert abs(plane_flux(0.5, 1.0) + 1) < 1e-6 and plane_flux_exact(0.5, 1.0) == -1
    assert abs(dipole_flux_integral(2.0)) < 1e-6 and abs(dipole_flux_integral(7.0)) < 1e-6
    # ---- §4.2.4–4.2.5 吸收球、鏡像、命中機率、電容；(4.64) 前的更正
    for r_ in (1.5, 3.0, 10.0):
        assert abs(sphere_flux_numeric(r_, 1.0) - sphere_flux(1.0)) < 1e-6
    assert _f(sphere_flux(1.0), 3) == "12.566"
    for R_, r0 in ((4.0, 1.0), (2.0, 1.0), (10.0, 3.0)):
        assert sphere_surface_max(R_, r0) < 1e-12
        assert abs(hit_prob_numeric(R_, r0) - hit_prob(R_, r0)) < 1e-6
    assert kelvin_image(4.0, 1.0) == (0.25, -0.25) and hit_prob(4.0, 1.0) == 0.25 and 1 - hit_prob(4.0, 1.0) == 0.75
    assert f"{cell_rmax(**CELL_EXAMPLE):.1e}" == "3.9e-05" and cell_rmax(1, 1, 1) == math.sqrt(3)
    q_, rr = 2.0, 1.7                                                                 # 高斯定律：Q = −(1/4π)∮∇φ·dS（∇²φ = −4πρ）
    h = 1e-5
    dphi = (q_ / (rr + h) - q_ / (rr - h)) / (2 * h)
    assert abs(-1 / (4 * PI) * dphi * 4 * PI * rr * rr - q_) < 1e-8
    assert capacitance_sphere(2.0) == 2.0
    # ---- §4.2.6 Berg–Purcell（書外）；圓盤電容 2a/π；p. 92–93 的更正
    assert abs(disk_capacitance_numeric(1.0) - disk_capacitance(1.0)) < 1e-9
    for p in (0.2, 0.5, 0.8):
        assert abs(disk_potential_offcenter(1.0, p) / PI ** 2 - 1) < 1e-7                         # 圓盤是等電位
    a_, r0_ = BP_EXAMPLE["a"], BP_EXAMPLE["r0"]
    for n in (100, 1000, 3142, 10 ** 4, 10 ** 5):
        assert abs(berg_purcell(n, a_, r0_) - berg_purcell_series(n, a_, r0_)) < 1e-12
    assert [_f(berg_purcell(n, a_, r0_), 3) for n in (100, 1000, 3142, 10 ** 4, 10 ** 5)] == \
        ["0.031", "0.241", "0.500", "0.761", "0.970"]
    nh = bp_half(a_, r0_)
    assert round(nh) == 3142 and abs(berg_purcell(nh, a_, r0_) - 0.5) < 1e-12
    assert f"{100 * bp_area_fraction(nh, a_, r0_):.3f}" == "0.079"                                   # 只占 0.079% 的表面
    assert abs(nh * a_ / r0_ - PI) < 1e-9 and nh * a_ * a_ / r0_ ** 2 < 1e-2                        # na ≈ πr0（不是 ≪ r0），na² ≪ r0²
    # ---- §4.3 吸收邊界 G(z; x, −Δ) = 0；Catalan；Ex. 4.1
    for N in range(0, 13):
        for x in range(0, 4):
            for y in range(0, 6):
                assert half_counts(N, x, y) == half_counts_image(N, x, y)
    assert [half_counts(2 * k, 0, 0) for k in range(6)] == [catalan(k) for k in range(6)] == [1, 1, 2, 5, 14, 42]
    for z in (0.2, 0.3, 0.45):
        for x in (0, 1, 3):
            assert half_gf_residual(z, x) < 1e-12
            cur, ser = {x: 1}, 0.0                                                    # 逐步遞迴的級數 Σ C'(N; x, 2) z^N
            for N in range(400):
                ser += cur.get(2, 0) * z ** N
                nxt = {}
                for p_, c_ in cur.items():
                    for q_ in (p_ - 1, p_ + 1):
                        if q_ >= 0:
                            nxt[q_] = nxt.get(q_, 0) + c_
                cur = nxt
            assert abs(ser - half_gf(z, x, 2)) < 1e-9
    z = 0.3                                                                           # 書上的 (4:S-2) 不成立
    assert abs(half_gf(z, 1, 0) - z * half_gf(z, 1, 1) - 1) > 0.1
    assert abs(half_gf(z, 1, 0) - z * half_gf(z, 1, 1) - 0) < 1e-12
    assert all(abs(half_gf(0.5 - 1e-12, x, y) - half_visits_critical(x, y)) < 1e-3 for x in range(3) for y in range(3))
    return True


if __name__ == "__main__":
    import sys
    sys.stdout.reconfigure(encoding="utf-8")
    check()
    print("walk_bc（第 4 章）：反射原理、兩牆之間、球與立方體、穩定態、點源與鏡像、吸收球、Berg–Purcell 全部驗證通過")
