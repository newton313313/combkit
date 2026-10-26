# 由 combkit/tools/sync_from_walk017.py 從 walk017 專案的 codes/utils/walk_shape.py 複製；請改 walk017 那份再同步。
"""walk_shape.py —— 《Elements of the Random Walk》第 6 章 The shape of a random walk 的數學
（Python＋numpy／scipy.integrate，不 import manim；不依賴其他模組，以便原樣複製進 combkit）

約定（書上式號與筆誤見 ref/md/ch06_formulas.md）：
  一條 N 步的漫步走過 N+1 個點 r_0 = 0, r_1, …, r_N（書上 (6.1) 寫「N 個足跡、除以 N」，(6.22) 起改用 N 步、N+1 點；本模組一律 N 步、N+1 點）。
  迴轉半徑張量 T = (1/(N+1)) Σ_j (r_j − r̄)(r_j − r̄)ᵀ（(6.1)(6.2)）；特徵值 R_1² ≥ R_2² ≥ …＝主迴轉半徑的平方，R_i 是長度。
  步伐 η_α = r_α − r_{α−1}；Kramers 矩陣 a_{αβ} = min(α,β)(N+1−max(α,β))/(N+1)²（含對角線，(6.22) 更正），T = Σ a_{αβ} η_α η_βᵀ（(6.21)）。
  a 的特徵值精確為 λ_k = 1/(4(N+1) sin²(kπ/(2N+2)))、特徵向量精確為 sin(kπα/(N+1))（k = 1..N）；大 N 時 λ_k → (N+1)/(π²k²)（(6.19)(6.27)）。
  格子漫步：d 維簡單立方，每步沿 2d 個方向之一走 1（⟨|η|²⟩ = 1）；高斯步伐 (6.28)：每個分量 N(0, 1/d)。
  個別非球度 asph(λ) = Σ_{i<j}(λ_i−λ_j)² / ((d−1)(Σλ)²)（(6.16) 括號內，更正）；平均非球度 A_d = 比值的平均之比（(6.12)）。
  模擬一律固定種子（numpy default_rng）；1/d 與 Fixman 的模擬用高斯步伐、N = 50，數值除以 N+1（書上 Fig. 6.12–6.15 的橫軸）。
  函數多半接受 numpy 陣列；回傳快取（lru_cache）陣列的函數不要就地修改結果。

§6.1 迴轉半徑張量 (6.1)–(6.5)、補充 1（(6:S1-18)）
    sites(steps)                    步伐 (N, d) → 點 (N+1, d)，從原點出發
    center_of_mass(pts)、tensor(pts) (6.1)(6.2) 的定義式（numpy）；tensor_pairs(pts) 兩兩差 (6:S2-2b)（第二種算法）
    tensor_frac(pts)                整數點的精確 T（Fraction 的 list of lists）；com_frac(pts)
    kramers(N)                      a_{αβ}（Fraction，N×N）；kramers_np(N)；kramers_book(N, kind)：'nodiag' 照 (6.22) 只有 α≠β、'noscale' 照 (6.23) 少 1/(N+1)²
    kramers_count(N)                補充 3 的計數：j < m ≤ n ≤ k 的點對 (j, k) 個數＝ m(N+1−n)（暴力列舉，整數矩陣）
    tensor_kramers(steps)、tensor_kramers_frac(steps)   (6.21)：Σ a η ηᵀ（第三種算法）
    principal(T)                    (特徵值由大到小, 特徵向量 columns)；二維時向量的 x ≥ 0
    rotate(pts, theta)              二維點繞重心轉 θ；extent(T, phi) = n·T·n，n = (cos φ, sin φ)（(6:S1-18)；最大＝R_1²）
    SMALL_STEPS                     主例：4 步 東東北西 → 點 (0,0),(1,0),(2,0),(2,1),(1,1)；small_example() → dict
                                    （重心 (6/5, 2/5)、R_g² = 4/5、T = [[14/25, 3/25],[3/25, 6/25]]、特徵值 3/5、1/5、主軸 (3,1)、非球度 1/4）
    lattice_steps(M, N, d, seed)    M 條格子漫步的步伐 (M, N, d)（int8）；batch_tensor(steps) → (M, d, d)；batch_eigs(T) → (M, d) 由大到小
    cloud(M=1000, N=100)            Fig. 6.1：二維 1000 條 100 步（點 (M, N+1, 2)）；cloud_stats() → dict（平均張量 (6.5) 與個別的主半徑）
    fig63(n=6, N=1000)              Fig. 6.3：6 條二維 1000 步（list of 點）；fig63_ratios() = R_1²/R_2²
    fig64(N=1000)                   Fig. 6.4：主例的二維 1000 步（點）；fig64_stats(theta) → dict（重心、T、特徵值、主軸角、轉 θ 之後的 T）
§6.1.4 不變量與非球度 (6.6)–(6.18)，Fig. 6.5，Ex. 6.1
    invariants(T)                   (Tr T, Tr T², Tr T³, det T)；det_from_traces(T) (6.9) 第二行
    asph(vals)                      個別非球度（最後一軸是特徵值；Fraction 也可）；asph_book16(vals, N) 照書上 (6.16)（分母沒平方、加到 N）
    A_d(d) = (4+2d)/(4+5d)          (6.15)（Fraction）；A_from_moments(d, T11², T11T22, T12²) (6.14)；A_S21(d, book) 補充 2 的 (6:S2-1)
    trace_a(N) = ((N+1)²−1)/(6(N+1))、trace_a2(N) = (P²−1)(2P²+7)/(180P²)、diag_a2(N) = (P⁴−1)/(30P³)（P = N+1，Fraction）
    moments_exact(N, d, kind)       (⟨T11²⟩, ⟨T11T22⟩, ⟨T12²⟩)：kind 'gauss'（Wick 配對，精確）或 'lattice'（再加同一步的四次累積量）
    moments_enum(N, d)              暴力列舉全部 (2d)^N 條格子漫步（第二種算法，小 N）
    A_exact(N, d, kind)             有限 N 的精確平均非球度（Fraction）；N → ∞ 為 A_d
    M_p(p)、Aprime(d)               (6.17)(6.18)：A'_d = (d/4)[3 + 4/d − (d/2)M_{d/2}]；Aprime_laplace(d) 第二種算法（1/(TrT)² = ∫t e^{−t TrT}dt，閉式被積函數）
    Aprime_closed(d)                d = 2：(5 − 7ζ(3)/2)/2；d = 4：4 − 3ζ(3)（第三種算法）
    sim3d(M=20000, N=100)           Fig. 6.5：三維格子 2 萬條 100 步的特徵值 (M, 3)（約 0.4 秒）；sim3d_stats() → dict；asph_hist(a, bins=20)
§6.2 d ≫ N 的「樓梯」(6.19)–(6.20)，Ex. 6.2–6.4
    staircase(N)                    第 n 步沿第 n 軸走 1：點 (N+1, N)（整數）；staircase_com(N) = [(N+1−n)/(N+1)]（Ex. 6.3(a)）
    eigs_exact(N)                   a 的精確特徵值（由大到小，numpy）；eigs_asym(N, k) = (N+1)/(π²k²)；sine_mode(N, k) 單位化的 sin(kπα/(N+1))
    staircase_asph(N)               d → ∞ 的非球度 Tr a²/(Tr a)² = (2P²+7)/(5(P²−1))（Fraction）；N → ∞ 為 2/5（(6.20)）
    orth_fraction(N, d)             Ex. 6.2：每步都垂直於之前所有步的比例 Π(1 − k/d)；orth_approx(N, d) = e^{−N²/2d}
    zeta2、zeta4                    Σ1/k² = π²/6、Σ1/k⁴ = π⁴/90（數值部分和＋尾巴）
§6.2.1 Kramers 恆等式、高斯步伐、Wick 配對 (6.21)–(6.37)，Ex. 6.5–6.7
    gauss_moment(m)                 標準高斯的 ⟨x^{2m}⟩（數值積分）＝(2m−1)!!；wick_pairings(m) 暴力列出 2m 個物件的配對；wick_count(m) = (2m−1)!!
    closed_cov(N)                   Ex. 6.5：閉合漫步（Σ η = 0）的步伐共變異（× d）＝ I − J/N（Fraction 矩陣）
    sim_gauss_moments(N, d, M, seed)  高斯步伐的模擬 ⟨T11²⟩、⟨T11T22⟩、⟨T12²⟩（第二種算法）
§6.2.2–6.2.5 1/d 展開 (6.38)–(6.71)，Fig. 6.12–6.14
    resolvent_weight(eps)           ∫ Im[1/(x − iε)] dx = π（(6.41)(6.42)；書上 (6.52) 漏了 π）
    sum_identity(n, M)              (6.60) 的部分和 Σ_{m≠n} n²/(m²−n²)（＋尾巴 n²/M）→ 3/4；shift_sum_exact(i, N) 用精確 λ 的 Σ_{j≠i} λ_j/(λ_i−λ_j)
    shift_1d(d) = 3/(4d)            (6.59)；mean_1d(k, d, N) = (N+1)/(π²k²)·(1+3/4d)（(6.61)）；var_1d(k, d, N) = 2((N+1)/(π²k²))²/d（(6.62)）
    eig_sim(d, N=50, M, seed)       高斯步伐 d 維的模擬（ξ = Uᵀη 的技巧）：(前 3 個特徵值/(N+1) 的 (M,3), Tr T/(N+1) 的 (M,))；
                                    M 預設：d=100、20 為 2000，d=3 為 20000（共約 1 秒）
    xi_identity_error(N, d, seed)   定義式 T 與 ξᵀ diag(λ) ξ 的特徵值差（驗證模擬技巧）
    onedee_table(d)                 dict：k = 1,2,3 的 模擬平均/(6.61)、模擬變異數/(6.62)、模擬平均、預測（都除以 N+1）
§6.2.6 主迴轉半徑的分布、Fixman (6.72)–(6.76)，Fig. 6.15，Ex. 6.8–6.9
    gamma_density(x, alpha, d)      (6.73) 歸一化：x^{d/2−1} e^{−dx/2α}（平均 α、變異數 2α²/d，Ex. 6.8）
    fixman_K(s, d, book=False)      (6.76) 的 K(s)（N → ∞、α_n = 1/(π²n²)，閉式 [√z/sinh√z]^{d/2}，z = 2is/d）；book=True 照書上的號（z = −2is/d）
    fixman_K_product(s, d, n_max)   同上，直接連乘 n ≤ n_max（第二種算法）
    fixman_density(x, d, book=False)  (6.74) 數值反轉（x = R²/(N+1)，向量）；book=True 得 ≈ 0（號錯，密度跑到負的 R²）
    fixman_mean_var(d, N=None)      累積量：平均 Σα_n、變異數 (2/d)Σα_n²（N=None → 1/6、1/(45d)；給 N 用精確 λ/(N+1)）
    fixman_moments_numeric(d)       由反轉的密度積分出 (總和, 平均, 變異數)（第二種算法）
§6.2.7 自迴避 (6.77)：A_saw(d) = (2d+4)/(5d+4) + 0.008(4−d)（只引用）
summary()                         場景卡與講義用的全部數字（dict，字串已照卡上的小數位格式化）

check() 驗證全部（每個結果至少兩種獨立算法、書上的數字、筆誤的更正、plan_ch06.md 場景卡與講義的每個數字）；直接執行本檔會跑 check()（約 2 秒）。
`python walk_shape.py --summary` 印出 summary() 的全部數字。
最慢的呼叫：sim3d 約 0.4 秒、eig_sim 三個 d 合計約 1 秒、fixman_density 120 點約 0.08 秒（前兩者有快取），其餘都很快。
"""
from __future__ import annotations

import math
from fractions import Fraction
from functools import lru_cache
from itertools import product

import numpy as np
from scipy import integrate

# 固定種子（場景與講義都用這些）
SEED_CLOUD = 601
SEED_FIG63 = 603
SEED_FIG64 = 643              # 比值約 8、橫向 53 × 縱向 31（適合 16:9 畫面）
SEED_FIG65 = 605
SEED_GAUSS = 606
SEED_D = {100: 6100, 20: 6020, 3: 6003}
N_D = 50                       # 1/d 與 Fixman 模擬的步數
M_D = {100: 2000, 20: 2000, 3: 20000}

ZETA3 = 1.2020569031595942


# ---------------------------------------------------------------- §6.1 張量
def sites(steps):
    steps = np.asarray(steps)
    out = np.zeros((steps.shape[0] + 1, steps.shape[1]), dtype=float if steps.dtype.kind == "f" else np.int64)
    out[1:] = np.cumsum(steps, axis=0)
    return out


def center_of_mass(pts):
    return np.asarray(pts, dtype=float).mean(axis=0)


def tensor(pts):
    X = np.asarray(pts, dtype=float)
    X = X - X.mean(axis=0)
    return X.T @ X / len(X)


def tensor_pairs(pts):
    """(6:S2-2b)：T = (1/(2n²)) Σ_i Σ_j (r_i − r_j)(r_i − r_j)ᵀ，n = 點數。"""
    X = np.asarray(pts, dtype=float)
    n = len(X)
    D = X[:, None, :] - X[None, :, :]
    return np.einsum("ijk,ijl->kl", D, D) / (2 * n * n)


def com_frac(pts):
    n = len(pts)
    d = len(pts[0])
    return [Fraction(sum(int(p[i]) for p in pts), n) for i in range(d)]


def tensor_frac(pts):
    n = len(pts)
    d = len(pts[0])
    c = com_frac(pts)
    return [[sum((int(p[i]) - c[i]) * (int(p[j]) - c[j]) for p in pts) / n for j in range(d)] for i in range(d)]


def kramers(N):
    P = N + 1
    return [[Fraction(min(a, b) * (P - max(a, b)), P * P) for b in range(1, N + 1)] for a in range(1, N + 1)]


def kramers_np(N):
    P = N + 1
    al = np.arange(1, N + 1)
    return np.minimum.outer(al, al) * (P - np.maximum.outer(al, al)) / P**2


def kramers_book(N, kind="nodiag"):
    """照書上：'nodiag' (6.22) 只寫 α<β 與 α>β（對角線 0）；'noscale' (6.23) x<(N+1−x>)（少 1/(N+1)²）。"""
    a = kramers_np(N)
    if kind == "nodiag":
        return a - np.diag(np.diag(a))
    if kind == "noscale":
        return a * (N + 1) ** 2
    raise ValueError(kind)


def kramers_count(N):
    """補充 3：(r_k − r_j)² 展開後含 η_m η_n（m ≤ n）的點對 j < k 個數（點編號 0..N）；回傳 N×N 整數矩陣。"""
    c = [[0] * N for _ in range(N)]
    for j in range(N + 1):
        for k in range(j + 1, N + 1):
            # r_k − r_j = η_{j+1} + … + η_k
            for m in range(j + 1, k + 1):
                for n in range(j + 1, k + 1):
                    c[m - 1][n - 1] += 1
    return c


def tensor_kramers(steps):
    eta = np.asarray(steps, dtype=float)
    a = kramers_np(len(eta))
    return eta.T @ a @ eta


def tensor_kramers_frac(steps):
    a = kramers(len(steps))
    d = len(steps[0])
    N = len(steps)
    return [[sum(a[x][y] * int(steps[x][i]) * int(steps[y][j]) for x in range(N) for y in range(N))
             for j in range(d)] for i in range(d)]


def principal(T):
    vals, vecs = np.linalg.eigh(np.asarray(T, dtype=float))
    order = np.argsort(vals)[::-1]
    vals, vecs = vals[order], vecs[:, order]
    for c in range(vecs.shape[1]):
        v = vecs[:, c]
        k = int(np.argmax(np.abs(v) > 1e-12))
        if v[k] < 0:
            vecs[:, c] = -v
    return vals, vecs


def rotate(pts, theta):
    X = np.asarray(pts, dtype=float)
    c = X.mean(axis=0)
    R = np.array([[math.cos(theta), -math.sin(theta)], [math.sin(theta), math.cos(theta)]])
    return (X - c) @ R.T + c


def extent(T, phi):
    phi = np.asarray(phi, dtype=float)
    T = np.asarray(T, dtype=float)
    c, s = np.cos(phi), np.sin(phi)
    return T[0, 0] * c * c + 2 * T[0, 1] * s * c + T[1, 1] * s * s


SMALL_STEPS = [(1, 0), (1, 0), (0, 1), (-1, 0)]


def small_example():
    pts = [tuple(int(v) for v in p) for p in sites(SMALL_STEPS)]
    T = tensor_frac(pts)
    tr = T[0][0] + T[1][1]
    det = T[0][0] * T[1][1] - T[0][1] ** 2
    disc = tr * tr - 4 * det                      # = 4/25 → √ = 2/5
    sq = Fraction(math.isqrt(disc.numerator), math.isqrt(disc.denominator))
    assert sq * sq == disc
    l1, l2 = (tr + sq) / 2, (tr - sq) / 2
    d2 = [sum((Fraction(p[i]) - com_frac(pts)[i]) ** 2 for i in range(2)) for p in pts]
    return dict(pts=pts, com=com_frac(pts), d2=d2, Rg2=sum(d2) / len(pts), T=T, eig=(l1, l2),
                vec=(3, 1), asph=((l1 - l2) / (l1 + l2)) ** 2)


def lattice_steps(M, N, d, seed):
    rng = np.random.default_rng(seed)
    s = rng.integers(0, 2 * d, size=(M, N))
    steps = np.zeros((M, N, d), dtype=np.int8)
    ax = s % d
    sign = np.where(s < d, 1, -1).astype(np.int8)
    np.put_along_axis(steps, ax[..., None], sign[..., None], axis=2)
    return steps


def batch_tensor(steps):
    st = np.asarray(steps, dtype=float)
    M, N, d = st.shape
    X = np.zeros((M, N + 1, d))
    X[:, 1:] = np.cumsum(st, axis=1)
    X -= X.mean(axis=1, keepdims=True)
    return np.einsum("mki,mkj->mij", X, X) / (N + 1)


def batch_eigs(T):
    return np.linalg.eigvalsh(T)[:, ::-1]


@lru_cache(maxsize=None)
def cloud(M=1000, N=100, seed=SEED_CLOUD):
    st = lattice_steps(M, N, 2, seed).astype(np.int64)
    out = np.zeros((M, N + 1, 2), dtype=np.int64)
    out[:, 1:] = np.cumsum(st, axis=1)
    return out


@lru_cache(maxsize=None)
def cloud_stats():
    pts = cloud()
    M, n, _ = pts.shape
    st = np.diff(pts, axis=1)
    T = batch_tensor(st)
    ev = batch_eigs(T)
    Tm = T.mean(axis=0)
    evm = principal(Tm)[0]
    return dict(Tmean=Tm, eig_of_mean=evm, ratio_of_mean=evm[0] / evm[1], R1=ev[:, 0].mean(), R2=ev[:, 1].mean(),
                ratio=ev[:, 0].mean() / ev[:, 1].mean(), exact_T11=trace_a(n - 1) / 2,
                asph_mean=float(np.mean(asph(ev))), A_ratio=float(_A_ratio(ev)),
                end_r2=float(np.mean(np.sum(pts[:, -1].astype(float) ** 2, axis=1))))


@lru_cache(maxsize=None)
def fig63(n=6, N=1000, seed=SEED_FIG63):
    st = lattice_steps(n, N, 2, seed).astype(np.int64)
    return tuple(sites(s) for s in st)


def fig63_ratios():
    out = []
    for p in fig63():
        v = principal(tensor(p))[0]
        out.append(v[0] / v[1])
    return out


@lru_cache(maxsize=None)
def fig64(N=1000, seed=SEED_FIG64):
    return sites(lattice_steps(1, N, 2, seed)[0].astype(np.int64))


def fig64_stats(theta=math.pi / 6):
    p = fig64()
    T = tensor(p)
    vals, vecs = principal(T)
    pr = rotate(p, theta)
    Tr = tensor(pr)
    return dict(com=center_of_mass(p), T=T, eig=vals, R=np.sqrt(vals), angle=math.degrees(math.atan2(vecs[1, 0], vecs[0, 0])),
                ratio=vals[0] / vals[1], Trot=Tr, eig_rot=principal(Tr)[0], asph=float(asph(vals)),
                Tkr=tensor_kramers(np.diff(p, axis=0)))


# ---------------------------------------------------------------- §6.1.4 非球度
def invariants(T):
    T = np.asarray(T, dtype=float)
    return np.trace(T), np.trace(T @ T), np.trace(T @ T @ T), np.linalg.det(T)


def det_from_traces(T):
    t1, t2, t3, _ = invariants(T)
    return (t1**3 + 2 * t3 - 3 * t1 * t2) / 6


def asph(vals):
    if isinstance(vals, (list, tuple)) and vals and isinstance(vals[0], Fraction):
        d = len(vals)
        num = sum((vals[i] - vals[j]) ** 2 for i in range(d) for j in range(i + 1, d))
        return num / ((d - 1) * sum(vals) ** 2)
    v = np.asarray(vals, dtype=float)
    d = v.shape[-1]
    num = d * np.sum(v * v, axis=-1) - np.sum(v, axis=-1) ** 2
    return num / ((d - 1) * np.sum(v, axis=-1) ** 2)


def asph_book16(vals, N):
    """書上 (6.16) 括號內照字面：分母 (d−1) Σ_{k=1}^{N} R_k²（沒平方；k > d 的 R_k 當 0）。"""
    v = list(vals)
    d = len(v)
    num = sum((v[i] - v[j]) ** 2 for i in range(d) for j in range(i + 1, d))
    return num / ((d - 1) * sum(v[:N]))


def _A_ratio(ev):
    d = ev.shape[1]
    num = np.mean(d * np.sum(ev * ev, axis=1) - np.sum(ev, axis=1) ** 2)
    return num / ((d - 1) * np.mean(np.sum(ev, axis=1) ** 2))


def A_d(d):
    return Fraction(4 + 2 * d, 4 + 5 * d)


def A_from_moments(d, t11sq, t1122, t12sq):
    return (t11sq - t1122 + d * t12sq) / (t11sq + (d - 1) * t1122)


def A_S21(d, book=False):
    """(6:S2-1) 用補充 2 的比值 ⟨T11²⟩/⟨T11T22⟩ = 9/5、⟨T11²⟩/⟨T12²⟩ = 9/2。book=True 照書上（分母少乘 ⟨T11T22⟩/⟨T12²⟩）。"""
    r1, r2 = Fraction(9, 5), Fraction(9, 2)          # T11²/T11T22, T11²/T12²
    t11t22_over_t12 = r2 / r1                         # = 5/2
    num = (r2 - t11t22_over_t12) + d
    den = r1 + d - 1
    return num / den if book else num / (den * t11t22_over_t12)


def trace_a(N):
    P = N + 1
    return Fraction(P * P - 1, 6 * P)


def trace_a2(N):
    P = N + 1
    return Fraction((P * P - 1) * (2 * P * P + 7), 180 * P * P)


def diag_a2(N):
    P = N + 1
    return Fraction(P**4 - 1, 30 * P**3)


def moments_exact(N, d, kind="gauss"):
    ta, ta2 = trace_a(N), trace_a2(N)
    t11sq = (ta * ta + 2 * ta2) / d**2
    t1122 = ta * ta / d**2
    t12sq = ta2 / d**2
    if kind == "lattice":
        D = diag_a2(N)
        t11sq += D * (Fraction(1, d) - Fraction(3, d**2))
        t1122 -= D / d**2
        t12sq -= D / d**2
    return t11sq, t1122, t12sq


def moments_enum(N, d):
    """暴力列舉全部 (2d)^N 條格子漫步（精確 Fraction）：(⟨T11²⟩, ⟨T11T22⟩, ⟨T12²⟩)。"""
    dirs = []
    for i in range(d):
        for s in (1, -1):
            e = [0] * d
            e[i] = s
            dirs.append(tuple(e))
    s11 = s1122 = s12 = Fraction(0)
    cnt = 0
    for w in product(dirs, repeat=N):
        T = tensor_frac([tuple(v) for v in sites(w)])
        s11 += T[0][0] ** 2
        s1122 += T[0][0] * T[1][1]
        s12 += T[0][1] ** 2
        cnt += 1
    return s11 / cnt, s1122 / cnt, s12 / cnt


def A_exact(N, d, kind="gauss"):
    return A_from_moments(d, *moments_exact(N, d, kind))


def M_p(p):
    def f(x):
        if x < 1e-8:
            return x
        return math.exp((p + 1) * math.log(x) - p * (x + math.log1p(-math.exp(-2 * x)) - math.log(2)))
    return integrate.quad(f, 0, np.inf, limit=200, epsabs=1e-13, epsrel=1e-12)[0]


def Aprime(d):
    return d / 4 * (3 + 4 / d - d / 2 * M_p(d / 2))


def _S1S2(y):
    """S1 = Σ 1/(π²k²+y²)，S2 = Σ 1/(π²k²+y²)²。"""
    if y < 0.05:
        y2 = y * y
        return 1 / 6 - y2 / 90 + y2 * y2 / 945, 1 / 90 - 2 * y2 / 945 + 3 * y2 * y2 / 9450
    f = y / math.tanh(y)
    fp = 1 / math.tanh(y) - y / math.sinh(y) ** 2 if y < 300 else 1.0
    S1 = (f - 1) / (2 * y * y)
    S2 = -fp / (4 * y**3) + (f - 1) / (2 * y**4)
    return S1, S2


def Aprime_laplace(d):
    """⟨TrT²/(TrT)²⟩ = ∫_0^∞ t⟨TrT² e^{−t TrT}⟩dt；高斯：傾斜後每個模態的變異數變小，Z = (y/sinh y)^{d/2}，y = √(2t/d)。"""
    def f(y):
        if y == 0:
            return 0.0
        S1, S2 = _S1S2(y)
        lsh = y + math.log1p(-math.exp(-2 * y)) - math.log(2)
        Z = math.exp(d / 2 * (math.log(y) - lsh))
        return y**3 * Z * ((d + 1) * S2 + S1 * S1)
    E = d / 2 * integrate.quad(f, 0, np.inf, limit=400, epsabs=1e-14, epsrel=1e-12)[0]
    return (d * E - 1) / (d - 1)


def Aprime_closed(d):
    if d == 2:
        return (5 - 3.5 * ZETA3) / 2
    if d == 4:
        return 4 - 3 * ZETA3
    raise ValueError(d)


@lru_cache(maxsize=None)
def sim3d(M=20000, N=100, seed=SEED_FIG65):
    out = np.empty((M, 3))
    B = 5000
    rng_seeds = np.random.SeedSequence(seed).spawn((M + B - 1) // B)
    for i, s in enumerate(range(0, M, B)):
        m = min(B, M - s)
        st = lattice_steps(m, N, 3, rng_seeds[i])
        out[s:s + m] = batch_eigs(batch_tensor(st))
    return out


@lru_cache(maxsize=None)
def sim3d_stats():
    ev = sim3d()
    a = asph(ev)
    m = ev.mean(axis=0)
    return dict(A_ratio=float(_A_ratio(ev)), A_mean=float(a.mean()), hist=asph_hist(a), frac_lt_half=float(np.mean(a < 0.5)),
                means=m, ratios=(m[0] / m[2], m[1] / m[2]), mean_trace=float(ev.sum(axis=1).mean()),
                ratio_of_ratios=(float(np.mean(ev[:, 0] / ev[:, 2])), float(np.mean(ev[:, 1] / ev[:, 2]))))


def asph_hist(a, bins=20):
    h, _ = np.histogram(a, bins=bins, range=(0.0, 1.0))
    return h


# ---------------------------------------------------------------- §6.2 樓梯
def staircase(N):
    pts = np.zeros((N + 1, N), dtype=np.int64)
    for j in range(1, N + 1):
        pts[j] = pts[j - 1]
        pts[j, j - 1] = 1
    return pts


def staircase_com(N):
    return [Fraction(N + 1 - n, N + 1) for n in range(1, N + 1)]


def eigs_exact(N):
    P = N + 1
    k = np.arange(1, N + 1)
    return 1.0 / (4 * P * np.sin(k * np.pi / (2 * P)) ** 2)


def eigs_asym(N, k):
    return (N + 1) / (math.pi**2 * k * k)


def sine_mode(N, k):
    al = np.arange(1, N + 1)
    v = np.sin(k * np.pi * al / (N + 1))
    return v / np.linalg.norm(v)


def staircase_asph(N):
    P = N + 1
    return Fraction(2 * P * P + 7, 5 * (P * P - 1))


def orth_fraction(N, d):
    f = Fraction(1)
    for k in range(N):
        f *= Fraction(d - k, d)
    return f


def orth_approx(N, d):
    return math.exp(-N * N / (2 * d))


def _zeta(s, K=100000):
    k = np.arange(1, K + 1, dtype=float)
    return float(np.sum(k ** (-s))) + K ** (1 - s) / (s - 1) - 0.5 * K ** (-s)


zeta2 = _zeta(2)
zeta4 = _zeta(4)


# ---------------------------------------------------------------- §6.2.1 高斯、Wick
def gauss_moment(m):
    f = lambda x: x ** (2 * m) * math.exp(-x * x / 2) / math.sqrt(2 * math.pi)
    return 2 * integrate.quad(f, 0, np.inf)[0]


def wick_pairings(m):
    def rec(items):
        if not items:
            return [[]]
        a, rest = items[0], items[1:]
        out = []
        for i, b in enumerate(rest):
            for sub in rec(rest[:i] + rest[i + 1:]):
                out.append([(a, b)] + sub)
        return out
    return rec(list(range(2 * m)))


def wick_count(m):
    r = 1
    for k in range(1, 2 * m, 2):
        r *= k
    return r


def closed_cov(N):
    return [[Fraction(1 if i == j else 0) - Fraction(1, N) for j in range(N)] for i in range(N)]


def sim_gauss_moments(N=100, d=3, M=20000, seed=SEED_GAUSS):
    rng = np.random.default_rng(seed)
    eta = rng.standard_normal((M, N, d)) / math.sqrt(d)
    T = batch_tensor(eta)
    return float(np.mean(T[:, 0, 0] ** 2)), float(np.mean(T[:, 0, 0] * T[:, 1, 1])), float(np.mean(T[:, 0, 1] ** 2)), T


# ---------------------------------------------------------------- §6.2.2–6.2.5 1/d
def resolvent_weight(eps=1e-3):
    f = lambda x: eps / (x * x + eps * eps)
    return 2 * integrate.quad(f, 0, np.inf, limit=200)[0]


def sum_identity(n, M=200000):
    m = np.arange(1, M + 1, dtype=float)
    m = m[m != n]
    return float(np.sum(n * n / (m * m - n * n))) + n * n / M


def shift_sum_exact(i, N=4000):
    lam = eigs_exact(N)
    li = lam[i - 1]
    lj = np.delete(lam, i - 1)
    return float(np.sum(lj / (li - lj)))


def shift_1d(d):
    return 3 / (4 * d)


def mean_1d(k, d, N=N_D):
    return eigs_asym(N, k) * (1 + shift_1d(d))


def var_1d(k, d, N=N_D):
    return 2 * eigs_asym(N, k) ** 2 / d


@lru_cache(maxsize=None)
def eig_sim(d, N=N_D, M=None, seed=None):
    M = M if M is not None else M_D.get(d, 2000)
    seed = seed if seed is not None else SEED_D.get(d, 6000 + d)
    rng = np.random.default_rng(seed)
    lam = eigs_exact(N)
    P = N + 1
    top = np.empty((M, 3))
    tr = np.empty(M)
    B = 500 if d > 10 else 5000
    for s in range(0, M, B):
        m = min(B, M - s)
        xi = rng.standard_normal((m, N, d)) / math.sqrt(d)
        if d <= N:
            T = np.matmul((xi * lam[None, :, None]).transpose(0, 2, 1), xi)
        else:
            Y = xi * np.sqrt(lam)[None, :, None]
            T = np.matmul(Y, Y.transpose(0, 2, 1))
        ev = np.linalg.eigvalsh(T)[:, ::-1]
        top[s:s + m] = ev[:, :3] / P
        tr[s:s + m] = ev.sum(axis=1) / P
    return top, tr


def xi_identity_error(N=12, d=5, seed=1):
    rng = np.random.default_rng(seed)
    eta = rng.standard_normal((N, d)) / math.sqrt(d)
    T1 = tensor(sites(eta))
    U = np.column_stack([sine_mode(N, k) for k in range(1, N + 1)])
    xi = U.T @ eta
    T2 = xi.T @ np.diag(eigs_exact(N)) @ xi
    return float(np.max(np.abs(np.linalg.eigvalsh(T1) - np.linalg.eigvalsh(T2))))


@lru_cache(maxsize=None)
def onedee_table(d):
    top, _ = eig_sim(d)
    P = N_D + 1
    mean = top.mean(axis=0)
    var = top.var(axis=0)
    pred_mean = np.array([mean_1d(k, d) / P for k in (1, 2, 3)])
    pred_var = np.array([var_1d(k, d) / P**2 for k in (1, 2, 3)])
    lead = np.array([eigs_asym(N_D, k) / P for k in (1, 2, 3)])
    return dict(mean=mean, var=var, pred_mean=pred_mean, pred_var=pred_var, lead=lead,
                mean_ratio=mean / pred_mean, var_ratio=var / pred_var, std=np.sqrt(var))


# ---------------------------------------------------------------- §6.2.6 Fixman
def gamma_density(x, alpha, d):
    x = np.asarray(x, dtype=float)
    k = d / 2
    th = 2 * alpha / d
    with np.errstate(divide="ignore", invalid="ignore"):
        out = np.where(x > 0, np.exp((k - 1) * np.log(np.where(x > 0, x, 1)) - x / th - math.lgamma(k) - k * math.log(th)), 0.0)
    return out


def fixman_K(s, d, book=False):
    s = np.asarray(s, dtype=float)
    z = (-2j if book else 2j) * s / d
    w = np.sqrt(z.astype(complex))
    small = np.abs(w) < 1e-6
    wsafe = np.where(small, 1.0, w)
    lsh = wsafe + np.log1p(-np.exp(-2 * wsafe)) - math.log(2)       # log sinh w（Re w > 0）
    out = np.exp(d / 2 * (np.log(wsafe) - lsh))
    return np.where(small, 1.0 + 0j, out)


def fixman_K_product(s, d, n_max=200000, book=False):
    n = np.arange(1, n_max + 1, dtype=float)
    al = 1 / (np.pi**2 * n * n)
    sg = -1 if book else 1
    out = []
    for si in np.atleast_1d(s):
        out.append(np.exp(-d / 2 * np.sum(np.log1p(sg * 2j * si * al / d))))
    return np.array(out)


def fixman_density(x, d, book=False, h=0.002):
    """(6.74)：P(x) = (1/2π)∫e^{isx}K(s)ds = (1/π)Re∫_0^∞ e^{isx}K(s)ds，s = u²。"""
    x = np.atleast_1d(np.asarray(x, dtype=float))
    U = 80 / math.sqrt(d) + 5
    u = np.arange(0, U + h, h)
    s = u * u
    g = fixman_K(s, d, book) * 2 * u
    w = np.full(len(u), h)
    w[0] = w[-1] = h / 2
    gw = g * w
    out = np.empty(len(x))
    for i0 in range(0, len(x), 40):
        xs = x[i0:i0 + 40]
        out[i0:i0 + 40] = (np.exp(1j * np.outer(xs, s)) @ gw).real / math.pi
    return out


def fixman_mean_var(d, N=None):
    if N is None:
        return 1 / 6, 2 / (90 * d)
    lam = eigs_exact(N) / (N + 1)
    return float(lam.sum()), float(2 / d * np.sum(lam * lam))


@lru_cache(maxsize=None)
def fixman_moments_numeric(d, xmax=1.6, n=801):
    x = np.linspace(0, xmax, n)
    p = fixman_density(x, d)
    tot = integrate.simpson(p, x=x)
    mu = integrate.simpson(x * p, x=x)
    var = integrate.simpson((x - mu) ** 2 * p, x=x)
    return float(tot), float(mu), float(var)


def A_saw(d):
    return (2 * d + 4) / (5 * d + 4) + 0.008 * (4 - d)


# ---------------------------------------------------------------- 場景卡與講義的數字
def f4(x):
    return f"{x:.4f}"


@lru_cache(maxsize=None)
def summary():
    S = {}
    se = small_example()
    S["small_com"] = se["com"]
    S["small_T"] = se["T"]
    S["small_eig"] = se["eig"]
    cs = cloud_stats()
    S["cloud_Tmean"] = (f"{cs['Tmean'][0, 0]:.2f}", f"{cs['Tmean'][1, 1]:.2f}", f"{cs['Tmean'][0, 1]:.2f}")
    S["cloud_exact_T11"] = f"{float(cs['exact_T11']):.2f}"
    S["cloud_eig_of_mean"] = tuple(f"{v:.2f}" for v in cs["eig_of_mean"])
    S["cloud_R"] = (f"{cs['R1']:.2f}", f"{cs['R2']:.2f}", f"{cs['ratio']:.2f}")
    S["fig63_ratios"] = tuple(f"{r:.1f}" for r in fig63_ratios())
    fs = fig64_stats()
    S["fig64_eig"] = tuple(f"{v:.1f}" for v in fs["eig"])
    S["fig64_R"] = tuple(f"{v:.2f}" for v in fs["R"])
    S["fig64_angle"] = f"{fs['angle']:.1f}"
    S["fig64_ratio"] = f"{fs['ratio']:.2f}"
    S["fig64_T"] = tuple(f"{v:.1f}" for v in (fs["T"][0, 0], fs["T"][1, 1], fs["T"][0, 1]))
    S["fig64_Trot"] = tuple(f"{v:.1f}" for v in (fs["Trot"][0, 0], fs["Trot"][1, 1], fs["Trot"][0, 1]))
    S["fig64_com"] = tuple(f"{v:.1f}" for v in fs["com"])
    S["fig64_asph"] = f4(fs["asph"])
    S["A_table"] = {d: f4(float(A_d(d))) for d in (2, 3, 4, 10)}
    S["Aprime"] = {d: f"{Aprime(d):.6f}" for d in (2, 3, 4)}
    st = sim3d_stats()
    S["sim3d_A_ratio"] = f4(st["A_ratio"])
    S["sim3d_A_mean"] = f4(st["A_mean"])
    S["A_exact_lattice100"] = f4(float(A_exact(100, 3, "lattice")))
    S["A_exact_gauss100"] = f4(float(A_exact(100, 3, "gauss")))
    S["sim3d_hist"] = tuple(int(v) for v in st["hist"])
    S["sim3d_frac_lt_half"] = f"{st['frac_lt_half']:.3f}"
    S["sim3d_ratios"] = tuple(f"{v:.2f}" for v in st["ratios"])
    S["stair3_eig"] = tuple(f4(v) for v in eigs_exact(3))
    S["stair3_asym"] = tuple(f4(eigs_asym(3, k)) for k in (1, 2, 3))
    S["stair_asph"] = {N: f4(float(staircase_asph(N))) for N in (3, 10, 100)}
    S["orth"] = (f4(float(orth_fraction(10, 100))), f4(orth_approx(10, 100)))
    S["kr4_eig"] = tuple(f4(v) for v in eigs_exact(4))
    S["kr4_asym"] = tuple(f4(eigs_asym(4, k)) for k in (1, 2, 3, 4))
    S["A_exact_gauss_N4_d3"] = f4(float(A_exact(4, 3, "gauss")))
    g = sim_gauss_moments()
    S["sim_gauss_moments"] = tuple(f"{v:.2f}" for v in g[:3])
    S["exact_gauss_moments100"] = tuple(f"{float(v):.2f}" for v in moments_exact(100, 3, "gauss"))
    S["sim_gauss_A"] = f4(A_from_moments(3, *g[:3]))
    for d in (100, 20, 3):
        t = onedee_table(d)
        S[f"od{d}_mean"] = tuple(f4(v) for v in t["mean"])
        S[f"od{d}_pred"] = tuple(f4(v) for v in t["pred_mean"])
        S[f"od{d}_std"] = tuple(f4(v) for v in t["std"])
        S[f"od{d}_mean_ratio"] = tuple(f"{v:.2f}" for v in t["mean_ratio"])
        S[f"od{d}_var_ratio"] = tuple(f"{v:.2f}" for v in t["var_ratio"])
    S["lead"] = tuple(f4(eigs_asym(N_D, k) / (N_D + 1)) for k in (1, 2, 3))
    _, tr3 = eig_sim(3)
    S["fix_sim"] = (f4(float(tr3.mean())), f4(float(tr3.var())))
    S["fix_exact_N"] = tuple(f4(v) for v in fixman_mean_var(3, N_D))
    S["fix_num"] = tuple(f4(v) for v in fixman_moments_numeric(3))
    top3, _ = eig_sim(3)
    xg = np.linspace(0.001, 0.8, 800)
    pg = fixman_density(xg, 3)
    S["fix_peak"] = (f"{xg[np.argmax(pg)]:.3f}", f"{pg.max():.2f}")
    S["fix_at"] = tuple(f"{v:.3f}" for v in fixman_density(np.array([0.1, 0.2, 0.3]), 3))
    S["fig615_mean"] = f4(float(top3[:, 0].mean()))
    S["fig615_mode"] = f4(1 / (3 * math.pi**2))
    S["fig615_peak"] = f"{float(gamma_density(1 / (3 * math.pi**2), 1 / math.pi**2, 3)):.2f}"
    S["saw3"] = f4(A_saw(3))
    return S


# ---------------------------------------------------------------- check
def check():
    pi2 = math.pi**2
    # ---- 小例子（W601–W602）：定義式 = 兩兩差 = Kramers（精確）
    se = small_example()
    assert se["pts"] == [(0, 0), (1, 0), (2, 0), (2, 1), (1, 1)]
    assert se["com"] == [Fraction(6, 5), Fraction(2, 5)]
    assert se["d2"] == [Fraction(8, 5), Fraction(1, 5), Fraction(4, 5), Fraction(1), Fraction(2, 5)]
    assert se["Rg2"] == Fraction(4, 5)
    assert se["T"] == [[Fraction(14, 25), Fraction(3, 25)], [Fraction(3, 25), Fraction(6, 25)]]
    assert tensor_kramers_frac(SMALL_STEPS) == se["T"]
    assert np.allclose(tensor_pairs(se["pts"]), np.array(se["T"], dtype=float))
    assert se["eig"] == (Fraction(3, 5), Fraction(1, 5))
    vals, vecs = principal(np.array(se["T"], dtype=float))
    assert np.allclose(vals, [0.6, 0.2]) and np.allclose(vecs[:, 0], np.array([3, 1]) / math.sqrt(10))
    assert se["asph"] == Fraction(1, 4) and asph(list(se["eig"])) == Fraction(1, 4)
    assert abs(math.degrees(math.atan2(1, 3)) - 18.43) < 0.005
    # 補充 1：n·T·n 的最大、最小＝特徵值
    phi = np.linspace(0, math.pi, 100001)
    e = extent(np.array(se["T"], dtype=float), phi)
    assert abs(e.max() - 0.6) < 1e-9 and abs(e.min() - 0.2) < 1e-9
    assert abs(math.degrees(phi[np.argmax(e)]) - 18.43) < 0.01

    # ---- Kramers 矩陣：計數 m(N+1−n)、(6.21) = 定義式、精確的特徵值與 sin 模態
    for N in (1, 2, 3, 4, 6):
        c = kramers_count(N)
        assert all(c[m - 1][n - 1] == min(m, n) * (N + 1 - max(m, n)) for m in range(1, N + 1) for n in range(1, N + 1))
        a = kramers(N)
        assert all(a[i][j] == Fraction(c[i][j], (N + 1) ** 2) for i in range(N) for j in range(N))
        assert sum(a[i][i] for i in range(N)) == trace_a(N)
        assert sum(a[i][j] ** 2 for i in range(N) for j in range(N)) == trace_a2(N)
        assert sum(a[i][i] ** 2 for i in range(N)) == diag_a2(N)
        A = kramers_np(N)
        assert np.allclose(np.sort(np.linalg.eigvalsh(A))[::-1], eigs_exact(N))
        for k in range(1, N + 1):
            v = sine_mode(N, k)
            assert np.allclose(A @ v, eigs_exact(N)[k - 1] * v)
    assert kramers(4) == [[Fraction(v, 25) for v in row] for row in [[4, 3, 2, 1], [3, 6, 4, 2], [2, 4, 6, 3], [1, 2, 3, 4]]]
    assert kramers(3) == [[Fraction(v, 16) for v in row] for row in [[3, 2, 1], [2, 4, 2], [1, 2, 3]]]
    rng = np.random.default_rng(7)
    for d in (2, 3, 5):
        st = rng.standard_normal((30, d))
        p = sites(st)
        assert np.allclose(tensor(p), tensor_kramers(st)) and np.allclose(tensor(p), tensor_pairs(p))
    # 筆誤 (6.22)：沒有對角線 → 小例子的 T11 = T22 = 0（應為 14/25、6/25）
    eta = np.array(SMALL_STEPS, dtype=float)
    Tb = eta.T @ kramers_book(4, "nodiag") @ eta
    assert abs(Tb[0, 0]) < 1e-12 and abs(Tb[1, 1]) < 1e-12
    # 筆誤 (6.23)：少 1/(N+1)² → 特徵值大 (N+1)² 倍
    assert np.allclose(np.sort(np.linalg.eigvalsh(kramers_book(4, "noscale")))[::-1], 25 * eigs_exact(4))
    # (6.26)：大 N 時 λ_k → (N+1)/(π²k²)
    assert abs(eigs_exact(2000)[0] / eigs_asym(2000, 1) - 1) < 1e-6

    # ---- 不變量 (6.6)–(6.10)
    Tm = np.array([[3.0, 1.0, 0.5], [1.0, 2.0, 0.2], [0.5, 0.2, 1.0]])
    t1, t2, t3, det = invariants(Tm)
    ev = np.linalg.eigvalsh(Tm)
    assert abs(t1 - ev.sum()) < 1e-12 and abs(t2 - np.sum(ev**2)) < 1e-12 and abs(t3 - np.sum(ev**3)) < 1e-12
    assert abs(det - det_from_traces(Tm)) < 1e-12 and abs(det - np.prod(ev)) < 1e-12
    lhs = (ev[0] - ev[1]) ** 2 + (ev[0] - ev[2]) ** 2 + (ev[1] - ev[2]) ** 2
    assert abs(lhs - (3 * t2 - t1**2)) < 1e-12                      # (6.10)
    # 形狀的例子：圓 0、線 1、(4,1,1) 與扁盤 (1,1,0) 都是 1/4
    F = Fraction
    assert asph([F(1), F(1), F(1)]) == 0 and asph([F(1), F(0), F(0)]) == 1
    assert asph([F(4), F(1), F(1)]) == F(1, 4) and asph([F(1), F(1), F(0)]) == F(1, 4) and asph([F(3), F(1)]) == F(1, 4)
    # 筆誤 (6.16)：分母沒平方 → 不是無因次（整條漫步放大 2 倍，值變 4 倍）
    v = [0.6, 0.2]
    assert abs(asph_book16(v, 4) - 0.2) < 1e-12 and abs(asph_book16([4 * x for x in v], 4) - 0.8) < 1e-12
    assert abs(asph(np.array(v)) - 0.25) < 1e-12 and abs(asph(4 * np.array(v)) - 0.25) < 1e-12

    # ---- 平均非球度 (6.12)–(6.15)：公式 = 有限 N 精確（N → ∞）= 補充 2 = 暴力列舉 = 模擬
    assert A_d(3) == F(10, 19) and A_d(2) == F(4, 7) and A_d(4) == F(1, 2) and A_d(10) == F(4, 9)
    assert abs(float(A_d(10**9)) - 0.4) < 1e-8
    for d in (2, 3, 4, 7):
        assert A_S21(d) == A_d(d) and A_S21(d, book=True) == A_d(d) * F(5, 2)
        big = A_exact(10**6, d, "gauss")
        assert abs(float(big) - float(A_d(d))) < 1e-5
        t11, t1122, t12 = moments_exact(10**6, d)
        assert abs(float(t11 / t1122) - 1.8) < 1e-5 and abs(float(t11 / t12) - 4.5) < 1e-5       # 9/5、9/2（(6:S2-26)(6:S2-27)）
    assert A_S21(3, book=True) == F(25, 19)
    # 補充 2：三張圖 2 + 16 + 18 = 36（書上寫「加 (6:S2-23)」應為 (6:S2-24)）；[z^N](1−z)^{−7} ~ N⁶/720 → 36/720 = 1/20 = (6.35c)
    assert F(2 + 16 + 18, 720) == F(1, 20)
    # 有限 N 的大 N 係數：⟨T11²⟩ → N²/20d²、⟨T11T22⟩ → N²/36d²、⟨T12²⟩ → N²/90d²（(6.35)–(6.37)）
    t11, t1122, t12 = moments_exact(10**6, 1)
    P2 = (10**6 + 1) ** 2
    assert abs(float(t11) / P2 - 1 / 20) < 1e-5 and abs(float(t1122) / P2 - 1 / 36) < 1e-5 and abs(float(t12) / P2 - 1 / 90) < 1e-5
    # 筆誤 (6.35b)：少了 2 → 7/180 而不是 1/20；(6.35a) 多一個 d²、少平方
    assert abs((zeta2**2 + 2 * zeta4) / math.pi**4 - 1 / 20) < 1e-9
    assert abs((1 / 36 + 1 / 90) - 7 / 180) < 1e-15 and abs(7 / 180 - 1 / 20) > 0.01
    assert abs((zeta2 + 2 * zeta4) / math.pi**4 - 1 / 20) > 0.003
    # 格子漫步的精確動差 = 暴力列舉
    for N, d in ((3, 2), (4, 2), (3, 3)):
        assert moments_enum(N, d) == moments_exact(N, d, "lattice")
    # 個別非球度的平均 (6.17) = Laplace 積分 = ζ(3) 閉式
    assert f"{Aprime(3):.6f}" == "0.394274" and f"{Aprime(3):.4f}" == "0.3943"      # 書上 0.394274…
    for d in (2, 3, 4, 6):
        assert abs(Aprime(d) - Aprime_laplace(d)) < 2e-7
    for d in (2, 4):
        assert abs(Aprime(d) - Aprime_closed(d)) < 1e-9
    assert abs(Aprime_laplace(200) - 0.4) < 0.005

    # ---- d ≫ N 的樓梯 (6.19)(6.20)、Ex. 6.2–6.4
    for N in (1, 2, 3, 5, 8):
        pts = staircase(N)
        assert com_frac([tuple(p) for p in pts]) == staircase_com(N)
        assert tensor_frac([tuple(p) for p in pts]) == kramers(N)          # 樓梯的 T 就是 a
        Tn = np.array(kramers(N), dtype=float)
        assert abs(float(staircase_asph(N)) - float(np.trace(Tn @ Tn) / np.trace(Tn) ** 2)) < 1e-12
        # Ex. 6.3(c)：Σ_j (x_{j,n}−⟨x_n⟩)(x_{j,m}−⟨x_m⟩) = (N+1)·n/(N+1)·(1−m/(N+1))，n ≤ m
        X = pts - pts.mean(axis=0)
        for n_ in range(1, N + 1):
            for m_ in range(n_, N + 1):
                assert abs(X[:, n_ - 1] @ X[:, m_ - 1] - (N + 1) * n_ / (N + 1) * (1 - m_ / (N + 1))) < 1e-12
    # Ex. 6.3(b) 的編號：點 j = 1..N+1 時第 n 軸座標是 1 ⇔ j > n（書上寫 j ≥ n，平均就變成 (N+2−n)/(N+1)）
    N = 3
    assert [Fraction(sum(1 for j in range(1, N + 2) if j >= n), N + 1) for n in (1, 2, 3)] == [F(1), F(3, 4), F(1, 2)]
    assert staircase_com(3) == [F(3, 4), F(1, 2), F(1, 4)]
    assert staircase_asph(3) == F(13, 25) and abs(float(staircase_asph(10**6)) - 0.4) < 1e-9
    e3 = eigs_exact(3)
    assert np.allclose(e3, [(2 + math.sqrt(2)) / 8, 1 / 8, (2 - math.sqrt(2)) / 8])
    assert abs(zeta2 - pi2 / 6) < 1e-12 and abs(zeta4 - math.pi**4 / 90) < 1e-12 and abs(zeta4 / zeta2**2 - 0.4) < 1e-12
    # Ex. 6.2：精確比例 vs e^{−N²/2d}
    assert orth_fraction(2, 3) == F(2, 3)
    assert abs(float(orth_fraction(10, 10**6)) - orth_approx(10, 10**6)) < 1e-4
    # Ex. 6.4：x<(1−x>) 是 −d²/dx² 的反：離散版 (N+1)·a 與三對角 (2,−1) 互逆
    for N in (3, 7):
        Lap = 2 * np.eye(N) - np.eye(N, k=1) - np.eye(N, k=-1)
        assert np.allclose(Lap @ ((N + 1) * kramers_np(N)), np.eye(N))

    # ---- 高斯步伐 (6.28)–(6.37)、Wick（Ex. 6.7）、閉合漫步（Ex. 6.5）
    for m in range(1, 6):
        assert len(wick_pairings(m)) == wick_count(m) == math.factorial(2 * m) // (2**m * math.factorial(m))
        assert abs(gauss_moment(m) - wick_count(m)) < 1e-8
    assert [wick_count(m) for m in range(1, 6)] == [1, 3, 15, 105, 945]
    C = closed_cov(4)
    assert all(sum(row) == 0 for row in C) and C[0][1] == F(-1, 4)
    s11, s1122, s12, _ = sim_gauss_moments()
    ex = [float(v) for v in moments_exact(100, 3, "gauss")]
    assert abs(s11 / ex[0] - 1) < 0.05 and abs(s1122 / ex[1] - 1) < 0.05 and abs(s12 / ex[2] - 1) < 0.05
    assert abs(A_from_moments(3, s11, s1122, s12) - float(A_exact(100, 3))) < 0.02

    # ---- 1/d 展開 (6.41)–(6.62)
    assert abs(resolvent_weight() - math.pi) < 1e-6                     # (6.42) 有 π，(6.52) 漏了
    for n in (1, 2, 3):
        assert abs(sum_identity(n) - 0.75) < 1e-5
        assert abs(shift_sum_exact(n) - 0.75) < 2e-3
    assert xi_identity_error() < 1e-12
    for d in (100, 20, 3):
        t = onedee_table(d)
        tol_m = 0.03 if d > 3 else 0.3
        assert np.all(np.abs(t["mean_ratio"] - 1) < tol_m)
        assert abs(t["var_ratio"][0] - 1) < 0.05                          # 最大的特徵值：寬度 (6.62) 準
        if d == 100:
            assert np.all(np.abs(t["var_ratio"] - 1) < 0.1)
    _, tr3 = eig_sim(3)
    mu_ex, var_ex = fixman_mean_var(3, N_D)
    assert abs(mu_ex - float(trace_a(N_D)) / (N_D + 1)) < 1e-12
    assert abs(var_ex - 2 / 3 * float(trace_a2(N_D)) / (N_D + 1) ** 2) < 1e-12
    assert abs(tr3.mean() / mu_ex - 1) < 0.01 and abs(tr3.var() / var_ex - 1) < 0.05

    # ---- Fixman (6.74)–(6.76)：閉式 K = 連乘；反轉的密度 歸一、平均 1/6、變異數 1/135；書上的號給 0
    s = np.array([0.5, 3.0, 20.0])
    assert np.allclose(fixman_K(s, 3), fixman_K_product(s, 3), rtol=1e-4)
    assert np.allclose(fixman_K(s, 3, book=True), fixman_K_product(s, 3, book=True), rtol=1e-4)
    tot, mu, var = fixman_moments_numeric(3)
    assert abs(tot - 1) < 1e-5 and abs(mu - 1 / 6) < 1e-5 and abs(var - 1 / 135) < 1e-5
    assert np.max(np.abs(fixman_density(np.array([0.05, 0.15, 0.3]), 3, book=True))) < 1e-5
    assert fixman_density(np.array([0.15]), 3)[0] > 3
    # (6.73)：歸一、平均 α、變異數 2α²/d（Ex. 6.8）；峰在 α(1−2/d)
    al = 1 / pi2
    x = np.linspace(1e-9, 3, 300001)
    gd = gamma_density(x, al, 3)
    assert abs(integrate.simpson(gd, x=x) - 1) < 1e-4 and abs(integrate.simpson(x * gd, x=x) - al) < 1e-5
    assert abs(integrate.simpson((x - al) ** 2 * gd, x=x) - 2 * al**2 / 3) < 1e-5
    assert abs(x[np.argmax(gd)] - al / 3) < 1e-4
    assert f"{A_saw(3):.4f}" == "0.5343"

    # ---- 場景卡與講義的數字（字串比對）
    S = summary()
    assert S["cloud_Tmean"] == ("8.35", "8.59", "0.03") and S["cloud_exact_T11"] == "8.42"
    assert S["cloud_eig_of_mean"] == ("8.60", "8.34")
    assert S["cloud_R"] == ("14.14", "2.80", "5.05")
    assert S["fig63_ratios"] == ("2.1", "3.0", "6.7", "5.5", "14.5", "4.6")
    assert S["fig64_eig"] == ("258.8", "31.6") and S["fig64_R"] == ("16.09", "5.62")
    assert S["fig64_angle"] == "28.4" and S["fig64_ratio"] == "8.18"
    assert S["fig64_T"] == ("207.5", "82.9", "95.0") and S["fig64_Trot"] == ("94.1", "196.3", "101.4")
    assert S["fig64_com"] == ("-17.9", "-8.9") and S["fig64_asph"] == "0.6117"
    fs = fig64_stats()
    assert np.allclose(fs["eig_rot"], fs["eig"]) and np.allclose(fs["Tkr"], fs["T"])          # 旋轉不變、(6.21) = 定義式
    assert np.allclose(tensor_pairs(fig64()[::10]), tensor(fig64()[::10]))
    assert S["A_table"] == {2: "0.5714", 3: "0.5263", 4: "0.5000", 10: "0.4444"}
    assert S["Aprime"] == {2: "0.396400", 3: "0.394274", 4: "0.393829"}
    assert S["sim3d_A_ratio"] == "0.5231" and S["sim3d_A_mean"] == "0.3939"
    assert S["A_exact_lattice100"] == "0.5235" and S["A_exact_gauss100"] == "0.5265"
    assert abs(float(S["sim3d_A_ratio"]) - float(S["A_exact_lattice100"])) < 0.005
    assert abs(float(S["sim3d_A_mean"]) - Aprime(3)) < 0.01
    assert S["sim3d_hist"] == (119, 616, 1183, 1631, 1890, 1795, 1808, 1621, 1669, 1529, 1505, 1305, 1065, 1000, 655, 389, 162, 56, 2, 0)
    assert S["sim3d_frac_lt_half"] == "0.693"
    st = sim3d_stats()
    assert abs(st["mean_trace"] / float(trace_a(100)) - 1) < 0.005                             # ⟨Tr T⟩ = Tr a（(6:S2-13)）
    assert S["sim3d_ratios"] == ("11.95", "2.69")                                              # ⟨R1²⟩:⟨R2²⟩:⟨R3²⟩
    assert tuple(f"{v:.2f}" for v in st["ratio_of_ratios"]) == ("13.69", "2.98")               # ⟨R1²/R3²⟩:⟨R2²/R3²⟩ ≈ 書上 (6.54)
    assert abs(1.01 * 9 / 0.75 - 12.12) < 1e-9 and abs(0.91 * 2.25 / 0.75 - 2.73) < 1e-9     # 書上 (6.69)–(6.71) 換算
    assert S["stair3_eig"] == ("0.4268", "0.1250", "0.0732") and S["stair3_asym"] == ("0.4053", "0.1013", "0.0450")
    assert S["stair_asph"] == {3: "0.5200", 10: "0.4150", 100: "0.4002"}
    assert S["orth"] == ("0.6282", "0.6065")
    assert S["kr4_eig"] == ("0.5236", "0.1447", "0.0764", "0.0553") and S["kr4_asym"] == ("0.5066", "0.1267", "0.0563", "0.0317")
    assert S["A_exact_gauss_N4_d3"] == "0.6013"
    assert S["sim_gauss_moments"] == ("56.09", "31.31", "12.56") and S["exact_gauss_moments100"] == ("56.67", "31.48", "12.60")
    assert S["sim_gauss_A"] == "0.5261"
    assert S["lead"] == ("0.1013", "0.0253", "0.0113")
    assert S["od100_mean"] == ("0.1027", "0.0254", "0.0114") and S["od100_pred"] == ("0.1021", "0.0255", "0.0113")
    assert S["od100_std"] == ("0.0142", "0.0035", "0.0015")
    assert S["od100_mean_ratio"] == ("1.01", "1.00", "1.00") and S["od100_var_ratio"] == ("0.98", "0.97", "0.94")
    assert S["od20_mean"] == ("0.1059", "0.0264", "0.0117") and S["od20_pred"] == ("0.1051", "0.0263", "0.0117")
    assert S["od20_std"] == ("0.0324", "0.0075", "0.0029")
    assert S["od20_mean_ratio"] == ("1.01", "1.01", "1.00") and S["od20_var_ratio"] == ("1.02", "0.87", "0.68")
    assert S["od3_mean"] == ("0.1280", "0.0286", "0.0105") and S["od3_pred"] == ("0.1267", "0.0317", "0.0141")
    assert S["od3_std"] == ("0.0830", "0.0142", "0.0045")
    assert S["od3_mean_ratio"] == ("1.01", "0.90", "0.75") and S["od3_var_ratio"] == ("1.01", "0.47", "0.23")
    # 書上 (6.63)–(6.71) 的比值 vs 本模組的模擬（差 ≤ 0.05）
    for got, book in ((onedee_table(100)["var_ratio"], (0.97, 0.93, 0.95)), (onedee_table(3)["var_ratio"], (0.98, 0.48, 0.23)),
                      (onedee_table(3)["mean_ratio"], (1.01, 0.91, 0.75))):
        assert np.all(np.abs(np.array(got) - np.array(book)) <= 0.05)
    assert S["fix_sim"] == ("0.1671", "0.0076") and S["fix_exact_N"] == ("0.1666", "0.0074") and S["fix_num"] == ("1.0000", "0.1667", "0.0074")
    assert S["fix_peak"] == ("0.111", "6.20") and S["fix_at"] == ("6.078", "3.281", "1.011")
    hd, _ = np.histogram(tr3, bins=30, range=(0, 0.6))
    hd = hd / (len(tr3) * 0.02)                                          # 模擬的直方圖密度 vs 反轉的曲線
    assert np.max(np.abs(hd[4:15] - fixman_density(np.arange(4, 15) * 0.02 + 0.01, 3))) < 0.4
    assert S["fig615_mean"] == "0.1280" and S["fig615_mode"] == "0.0338" and S["fig615_peak"] == "7.16"
    assert S["saw3"] == "0.5343"
    # 其餘卡上的數字
    Ts = se["T"]
    assert Ts[0][0] * Ts[1][1] - Ts[0][1] ** 2 == F(3, 25) and Ts[0][0] + Ts[1][1] == F(4, 5)
    assert Ts[0][0] ** 2 + 2 * Ts[0][1] ** 2 + Ts[1][1] ** 2 == F(2, 5)                     # Tr T² = (3/5)² + (1/5)²
    e = extent(np.array(Ts, dtype=float), phi)
    v2 = principal(np.array(Ts, dtype=float))[1][:, 1]
    ang2 = math.degrees(math.atan2(v2[1], v2[0])) % 180
    assert f"{ang2:.2f}" == "108.43" and abs(math.degrees(phi[np.argmin(e)]) - ang2) < 0.01
    vr = principal(fs["Trot"])[1][:, 0]
    assert f"{math.degrees(math.atan2(vr[1], vr[0])):.1f}" == "58.4"
    assert trace_a(4) == F(4, 5) and trace_a(3) == F(5, 8)
    for k, v in ((1, [1, math.sqrt(2), 1]), (2, [1, 0, -1]), (3, [1, -math.sqrt(2), 1])):
        assert abs(abs(sine_mode(3, k) @ (np.array(v) / np.linalg.norm(v))) - 1) < 1e-12
    assert (f"{zeta2:.4f}", f"{zeta2**2:.4f}", f"{zeta4:.4f}") == ("1.6449", "2.7058", "1.0823")
    assert f"{resolvent_weight():.4f}" == "3.1416"
    assert tuple(f"{shift_1d(d):.4f}" for d in (100, 20, 3)) == ("0.0075", "0.0375", "0.2500")
    assert tuple(f"{math.sqrt(2 / d):.4f}" for d in (100, 20, 3)) == ("0.1414", "0.3162", "0.8165")
    assert (f"{1 / 6:.4f}", f"{1 / 135:.4f}", f"{math.sqrt(1 / 135):.4f}") == ("0.1667", "0.0074", "0.0861")
    return True


if __name__ == "__main__":
    import sys
    import time
    sys.stdout.reconfigure(encoding="utf-8")
    if "--summary" in sys.argv:
        for k, v in summary().items():
            print(k, v)
        sys.exit(0)
    t0 = time.time()
    check()
    print(f"walk_shape（第 6 章）：迴轉半徑張量、Kramers 恆等式、非球度、無限多維的樓梯、Wick 配對、1/d 展開、Fixman 分布全部驗證通過（{time.time() - t0:.1f} 秒）")
