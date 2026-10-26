# 由 combkit/tools/sync_from_symosc.py 從 symosc 專案的 codes/utils/symosc.py 複製；請改 symosc 那份再同步。
"""
symosc.py —— 對稱群與耦合振盪子：D3、破缺到 C3／C2、四元群 Q8 的數學（numpy；不 import manim）

畫面上的每個數字都從這裡算，場景的 check_math() 再 assert。直接執行本檔會跑 check()。

記號：振盪子 m ẍ = −K x（m = 1）；對稱操作 D(g) 是置換矩陣，(D(g) x)_{g·j} = x_j；[D(g), K] = 0。
ω = e^{2πi/3}（D3／C3）；Q8 = {±1, ±i, ±j, ±k}，i² = j² = k² = ijk = −1。

D3（三個相同的擺放在正三角形頂點，兩兩用彈簧 k 相連，自己的回復常數 k0）
    D3_NAMES                 ['e', 'r', 'r2', 's1', 's2', 's3']；r：0→1→2→0；s1 固定擺 1（index 0）交換 2、3；s2 固定 2；s3 固定 3
    d3_matrix(name)          3×3 置換矩陣；D3_CLASSES = {'e': [...], 'C3': [...], 'σ': [...]}
    K_d3(k0, k)              剛性矩陣；d3_chars_of(D) 每個類的跡；D3_CHARS 特徵標表（A1、A2、E）
    decompose(chars, table, sizes)   約化：每個不可約表示出現幾次
    projector_d3(irrep)      P^α = (d/6) Σ χ_α(g) D(g)（實特徵標）
    class_sum_values(table, sizes)   類和在每個不可約表示上的值 |類|·χ/d
    d3_freqs(k0, k)          ω_A = √k0、ω_E = √(k0+3k)
    beat(t, k0, k)           x(0) = (1,0,0)、ẋ(0) = 0 的解 (x1, x2, x3)
破缺到 C3（帶相位的耦合 b e^{iφ}；或陀螺項 G = g(P − Pᵀ)）與 C2（2–3 之間的彈簧變成 k+δ）
    P_CYC                    輪換矩陣 P_{j,j+1} = 1；P v_m = ω^m v_m，v_m = (1, ω^m, ω^{2m})/√3
    H_flux(a, b, phi)        a I + b e^{iφ} P + b e^{−iφ} Pᵀ；flux_levels(a, b, phi) = a + 2b cos(φ + 2πm/3)，m = 0, 1, −1
    gyro_freqs(k0, k, g)     m ẍ + G ẋ + K x = 0 的正頻率（m = 0, +1, −1）；裂距 √3 g
    precession(t, delta)     駐波 v_{+1} + v_{−1} 在頻率差 delta = ω_{−1} − ω_{+1} 下的振幅分布 2cos(2πj/3 − δt/2)
                             （圖樣往 j 增加的方向轉：1→2→3；check() 用 ODE 積分核對方向）
    K_c2(k0, k, delta)       換掉一根彈簧；c2_levels → k0、k0+3k（對稱 (2,−1,−1)）、k0+3k+2δ（反對稱 (0,1,−1)）
Q8（8 個擺放在 8 個群元素上，h 和 h·s 用彈簧相連：s = ±i、±j、±k 各用 ki、kj、kk）
    Q8_NAMES                 ['1', 'i', 'j', 'k', '-1', '-i', '-j', '-k']；q8_mul(a, b)；q8_inv(a)
    Q8_REP2                  二維不可約表示 ρ(i) = diag(i, −i)、ρ(j) = [[0,−1],[1,0]]、ρ(k) = ρ(i)ρ(j)
    L(g)、R(g)               左乘 (L(g)x)_h = x_{g⁻¹h}、右乘 (R(g)x)_h = x_{h g}（8×8；兩者交換）
    K_q8(k0, ki, kj, kk)     (k0 + 2ki + 2kj + 2kk) I − Σ_s k_s R(s)
    Q8_CHARS、Q8_CLASSES     特徵標表（A1..A4、E）
    q8_levels(k0, ki, kj, kk)   {A1: k0, A2: k0+4kj+4kk, A3: k0+4ki+4kk, A4: k0+4ki+4kj, E: k0+2(ki+kj+kk)}
    q8_E_basis()             E 的 4 個基底 v1..v4 與標籤 (L(i) 的本徵值 m, R(i) 的本徵值 n)
    frobenius_schur(chars, sizes)   (1/|G|) Σ χ(g²)：E 是 −1（四元數型）
"""
from __future__ import annotations

import cmath
import math

import numpy as np

OMEGA = cmath.exp(2j * math.pi / 3)

# =========================================================================== D3

D3_NAMES = ["e", "r", "r2", "s1", "s2", "s3"]
_D3_PERM = {"e": (0, 1, 2), "r": (1, 2, 0), "r2": (2, 0, 1), "s1": (0, 2, 1), "s2": (2, 1, 0), "s3": (1, 0, 2)}
D3_CLASSES = {"e": ["e"], "C3": ["r", "r2"], "σ": ["s1", "s2", "s3"]}
D3_CHARS = {"A1": (1, 1, 1), "A2": (1, 1, -1), "E": (2, -1, 0)}
D3_SIZES = (1, 2, 3)


def perm_matrix(p) -> np.ndarray:
    """p[j] = 擺 j 被送到的位置：D[p[j], j] = 1。"""
    n = len(p)
    D = np.zeros((n, n))
    for j, pj in enumerate(p):
        D[pj, j] = 1
    return D


def d3_matrix(name: str) -> np.ndarray:
    return perm_matrix(_D3_PERM[name])


def K_d3(k0: float, k: float) -> np.ndarray:
    return (k0 + 2 * k) * np.eye(3) - k * (np.ones((3, 3)) - np.eye(3))


def d3_chars_of(D) -> tuple:
    return tuple(int(round(np.trace(D(D3_CLASSES[c][0])).real)) for c in D3_CLASSES)


def decompose(chars, table: dict, sizes) -> dict:
    order = sum(sizes)
    return {a: round(sum(s * c * x for s, c, x in zip(sizes, chars, ch)) / order) for a, ch in table.items()}


def class_of_d3(name: str) -> int:
    return [i for i, c in enumerate(D3_CLASSES.values()) if name in c][0]


def projector_d3(irrep: str) -> np.ndarray:
    ch = D3_CHARS[irrep]
    return ch[0] / 6 * sum(ch[class_of_d3(g)] * d3_matrix(g) for g in D3_NAMES)


def class_sum_values(table: dict, sizes) -> dict:
    return {a: tuple(s * x / ch[0] for s, x in zip(sizes, ch)) for a, ch in table.items()}


def d3_freqs(k0: float, k: float) -> tuple:
    return math.sqrt(k0), math.sqrt(k0 + 3 * k)


def beat(t: float, k0: float, k: float) -> tuple:
    wa, we = d3_freqs(k0, k)
    a, e = math.cos(wa * t) / 3, math.cos(we * t) / 3
    return a + 2 * e, a - e, a - e


# =========================================================================== 破缺：C3 與 C2

P_CYC = np.roll(np.eye(3), 1, axis=1)          # P[j, j+1] = 1


def v_m(m: int) -> np.ndarray:
    return np.array([OMEGA ** (m * j) for j in range(3)]) / math.sqrt(3)


def H_flux(a: float, b: float, phi: float) -> np.ndarray:
    return a * np.eye(3) + b * cmath.exp(1j * phi) * P_CYC + b * cmath.exp(-1j * phi) * P_CYC.T


def flux_levels(a: float, b: float, phi: float) -> dict:
    return {m: a + 2 * b * math.cos(phi + 2 * math.pi * m / 3) for m in (0, 1, -1)}


def gyro_freqs(k0: float, k: float, g: float) -> dict:
    """m = 1：ω² + 2g sinθ_m ω − κ_m = 0 的正根（x ∝ v_m e^{iωt}）。"""
    out = {}
    for m in (0, 1, -1):
        kap = k0 if m == 0 else k0 + 3 * k
        s = g * math.sin(2 * math.pi * m / 3)
        out[m] = -s + math.sqrt(s * s + kap)
    return out


def gyro_freqs_numeric(k0: float, k: float, g: float) -> list:
    """一階化 [x, v]' = [[0, I], [−K, −G]] [x, v] 的本徵值 iω → 排序後的正 ω。"""
    G = g * (P_CYC - P_CYC.T)
    A = np.block([[np.zeros((3, 3)), np.eye(3)], [-K_d3(k0, k), -G]])
    w = np.linalg.eigvals(A)
    return sorted(round(x.imag, 12) for x in w if x.imag > 0)


def precession(t: float, delta: float) -> list:
    """cos(θ_j + ω₊t) + cos(θ_j − ω₋t) = 2cos(ω̄t)·cos(θ_j − δt/2)：m = +1 是 v₁e^{iω₊t}、m = −1 是 v₋₁e^{iω₋t}。"""
    return [2 * math.cos(2 * math.pi * j / 3 - delta * t / 2) for j in range(3)]


def K_c2(k0: float, k: float, delta: float) -> np.ndarray:
    K = K_d3(k0, k)
    K[1, 1] += delta
    K[2, 2] += delta
    K[1, 2] -= delta
    K[2, 1] -= delta
    return K


def c2_levels(k0: float, k: float, delta: float) -> dict:
    return {"A1": k0, "sym": k0 + 3 * k, "anti": k0 + 3 * k + 2 * delta}


# =========================================================================== Q8

Q8_NAMES = ["1", "i", "j", "k", "-1", "-i", "-j", "-k"]
_U = {"1": np.eye(2, dtype=complex), "i": np.diag([1j, -1j]), "j": np.array([[0, -1], [1, 0]], dtype=complex)}
_U["k"] = _U["i"] @ _U["j"]
Q8_REP2 = {n: (-1 if n.startswith("-") else 1) * _U[n.lstrip("-")] for n in Q8_NAMES}
Q8_CLASSES = {"1": ["1"], "-1": ["-1"], "±i": ["i", "-i"], "±j": ["j", "-j"], "±k": ["k", "-k"]}
Q8_SIZES = (1, 1, 2, 2, 2)
Q8_CHARS = {"A1": (1, 1, 1, 1, 1), "A2": (1, 1, 1, -1, -1), "A3": (1, 1, -1, 1, -1),
            "A4": (1, 1, -1, -1, 1), "E": (2, -2, 0, 0, 0)}


def q8_mul(a: str, b: str) -> str:
    M = Q8_REP2[a] @ Q8_REP2[b]
    return next(n for n in Q8_NAMES if np.allclose(Q8_REP2[n], M))


def q8_inv(a: str) -> str:
    return next(n for n in Q8_NAMES if q8_mul(a, n) == "1")


def _idx(n: str) -> int:
    return Q8_NAMES.index(n)


def L(g: str) -> np.ndarray:
    """(L(g)x)_h = x_{g⁻¹h}：把擺 h 送到 g·h。"""
    return perm_matrix([_idx(q8_mul(g, h)) for h in Q8_NAMES])


def R(g: str) -> np.ndarray:
    """(R(g)x)_h = x_{h g}。"""
    M = np.zeros((8, 8))
    for h in Q8_NAMES:
        M[_idx(h), _idx(q8_mul(h, g))] = 1
    return M


def K_q8(k0: float, ki: float, kj: float, kk: float) -> np.ndarray:
    ks = {"i": ki, "j": kj, "k": kk}
    K = (k0 + 2 * (ki + kj + kk)) * np.eye(8)
    for u, kv in ks.items():
        K -= kv * (R(u) + R("-" + u))
    return K


def q8_class(n: str) -> int:
    return [i for i, c in enumerate(Q8_CLASSES.values()) if n in c][0]


def q8_levels(k0: float, ki: float, kj: float, kk: float) -> dict:
    return {"A1": k0, "A2": k0 + 4 * kj + 4 * kk, "A3": k0 + 4 * ki + 4 * kk, "A4": k0 + 4 * ki + 4 * kj,
            "E": k0 + 2 * (ki + kj + kk)}


def q8_irrep_vector(irrep: str) -> np.ndarray:
    """一維表示：x_h = χ_α(h)。"""
    ch = Q8_CHARS[irrep]
    return np.array([ch[q8_class(h)] for h in Q8_NAMES], dtype=float)


def q8_E_basis() -> list:
    """[(名稱, 向量, m = L(i) 的本徵值, n = R(i) 的本徵值)]；v3 = L(j)v1、v4 = L(j)v2。"""
    v1 = np.zeros(8, dtype=complex)
    v2 = np.zeros(8, dtype=complex)
    for h, a, b in (("1", 1, 1), ("i", -1j, 1j), ("-1", -1, -1), ("-i", 1j, -1j)):
        v1[_idx(h)], v2[_idx(h)] = a, b
    vs = [("v1", v1), ("v2", v2), ("v3", L("j") @ v1), ("v4", L("j") @ v2)]
    out = []
    for name, v in vs:
        m = _eig_of(L("i"), v)
        n = _eig_of(R("i"), v)
        out.append((name, v, m, n))
    return out


def _eig_of(A: np.ndarray, v: np.ndarray) -> complex:
    w = A @ v
    j = int(np.argmax(abs(v)))
    lam = w[j] / v[j]
    assert np.allclose(w, lam * v)
    return complex(round(lam.real, 12), round(lam.imag, 12))


def frobenius_schur(chars, classes: dict, sizes, square) -> float:
    """(1/|G|) Σ_g χ(g²)。"""
    names = [n for c in classes.values() for n in c]
    cls = {n: i for i, c in enumerate(classes.values()) for n in c}
    return sum(chars[cls[square(g)]] for g in names) / len(names)


# =========================================================================== check


def check():
    # ---- D3：對稱、約化、投影、類和、頻率
    K = K_d3(1.0, 0.5)
    for g in D3_NAMES:
        assert np.allclose(d3_matrix(g) @ K, K @ d3_matrix(g))
    chi = d3_chars_of(d3_matrix)
    assert chi == (3, 0, 1)
    assert decompose(chi, D3_CHARS, D3_SIZES) == {"A1": 1, "A2": 0, "E": 1}
    assert np.allclose(projector_d3("A1") @ [1, 0, 0], np.array([1, 1, 1]) / 3)
    assert np.allclose(projector_d3("E") @ [1, 0, 0], np.array([2, -1, -1]) / 3)
    assert np.allclose(projector_d3("A2"), 0)
    assert class_sum_values(D3_CHARS, D3_SIZES) == {"A1": (1, 2, 3), "A2": (1, 2, -3), "E": (1, -1, 0)}
    w = sorted(np.linalg.eigvalsh(K))
    assert np.allclose(w, [1.0, 2.5, 2.5])
    assert np.allclose(K @ [0, 1, -1], 2.5 * np.array([0, 1, -1]))
    wa, we = d3_freqs(1.0, 0.5)
    assert round(wa, 4) == 1.0 and round(we, 4) == 1.5811 and round(2 * math.pi / (we - wa), 2) == 10.81
    # 拍：數值積分 vs 公式；最小振幅 1/3、能量 1/9
    x, v, dt = np.array([1.0, 0, 0]), np.zeros(3), 1e-3
    for _ in range(5000):                           # 速度 Verlet 到 t = 5
        v += -0.5 * dt * (K @ x)
        x += dt * v
        v += -0.5 * dt * (K @ x)
    assert np.allclose(x, beat(5.0, 1.0, 0.5), atol=1e-5)
    assert abs(1 / 3 - 2 / 3) == 1 / 3 and round((1 / 3) ** 2, 4) == 0.1111
    # ---- C3：實對稱的 C3 不變矩陣自動 D3 不變
    for b, c in ((0.3, 0.3), (0.3, 0.7)):
        M = 2 * np.eye(3) + b * P_CYC + c * P_CYC @ P_CYC
        assert np.allclose(M @ P_CYC, P_CYC @ M)
        assert np.allclose(M, M.T) == (b == c) == np.allclose(M @ d3_matrix("s1"), d3_matrix("s1") @ M)
    for m in (0, 1, -1):
        assert np.allclose(P_CYC @ v_m(m), OMEGA ** m * v_m(m))
    a, b = 2.0, -0.5                                 # = K_d3(1, 0.5)
    for phi in (0.0, 0.2, 0.5):
        H = H_flux(a, b, phi)
        lev = flux_levels(a, b, phi)
        for m in (0, 1, -1):
            assert np.allclose(H @ v_m(m), lev[m] * v_m(m))
        assert abs(lev[1] - lev[-1] + 2 * math.sqrt(3) * b * math.sin(phi)) < 1e-12
    assert flux_levels(a, b, 0.0) == {0: 1.0, 1: 2.5, -1: 2.5} or np.allclose(list(flux_levels(a, b, 0.0).values()), [1, 2.5, 2.5])
    gf = gyro_freqs(1.0, 0.5, 0.2)
    assert np.allclose(sorted(gf.values()), gyro_freqs_numeric(1.0, 0.5, 0.2))
    assert abs(gf[-1] - gf[1] - math.sqrt(3) * 0.2) < 1e-12 and round(gf[0], 6) == 1.0
    assert round(gf[1], 4) == 1.4174 and round(gf[-1], 4) == 1.7638
    assert np.allclose(precession(0, 1.0), [2, -1, -1])
    # 進動方向：一階化 [x, v]' = A [x, v] 用 RK4 積分 x(0) = (2,−1,−1)、ẋ(0) = 0，
    # 每個時刻前後 2（≈ 一個載波週期）內振幅最大的擺 = precession() 的 |峰|（1→3→2 每半格、1→2→3 每格）
    d = gf[-1] - gf[1]
    hop = 4 * math.pi / (3 * d)
    G = 0.2 * (P_CYC - P_CYC.T)
    A = np.block([[np.zeros((3, 3)), np.eye(3)], [-K_d3(1.0, 0.5), -G]])
    y, h, xs = np.array([2.0, -1, -1, 0, 0, 0]), 2e-3, []
    for _ in range(int((2 * hop + 2) / h) + 1):
        xs.append(y[:3].copy())
        k1 = A @ y
        k2 = A @ (y + h / 2 * k1)
        k3 = A @ (y + h / 2 * k2)
        k4 = A @ (y + h * k3)
        y = y + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
    xs = np.array(xs)
    for T, peak in ((hop / 2, 2), (hop, 1), (2 * hop, 2)):
        win = xs[int((T - 2) / h):int((T + 2) / h)]
        assert int(np.argmax(np.abs(win).max(axis=0))) == peak == int(np.argmax(np.abs(precession(T, d))))
    # ---- C2
    Kc = K_c2(1.0, 0.5, 0.25)
    assert np.allclose(Kc @ [2, -1, -1], 2.5 * np.array([2, -1, -1]))
    assert np.allclose(Kc @ [0, 1, -1], 3.0 * np.array([0, 1, -1]))
    assert np.allclose(Kc @ d3_matrix("s1"), d3_matrix("s1") @ Kc) and not np.allclose(Kc @ P_CYC, P_CYC @ Kc)
    assert c2_levels(1.0, 0.5, 0.25) == {"A1": 1.0, "sym": 2.5, "anti": 3.0}
    # ---- Q8
    assert q8_mul("i", "j") == "k" and q8_mul("j", "i") == "-k" and q8_mul("i", "i") == "-1"
    assert q8_mul(q8_mul("i", "j"), "k") == "-1"
    for g in Q8_NAMES:
        for h in Q8_NAMES:
            assert np.allclose(L(g) @ R(h), R(h) @ L(g))
    Kq = K_q8(1.0, 0.2, 0.3, 0.45)
    assert np.allclose(Kq, Kq.T)
    for g in Q8_NAMES:
        assert np.allclose(L(g) @ Kq, Kq @ L(g))
    chiL = tuple(int(round(np.trace(L(c[0])))) for c in Q8_CLASSES.values())
    assert chiL == (8, 0, 0, 0, 0)
    assert decompose(chiL, Q8_CHARS, Q8_SIZES) == {"A1": 1, "A2": 1, "A3": 1, "A4": 1, "E": 2}
    for a_, ch in Q8_CHARS.items():                  # 特徵標表 = ρ 的跡（E）與正交
        if a_ == "E":
            assert all(np.isclose(np.trace(Q8_REP2[c[0]]), x) for c, x in zip(Q8_CLASSES.values(), ch))
    lev = q8_levels(1.0, 0.2, 0.3, 0.45)
    assert lev == {"A1": 1.0, "A2": 4.0, "A3": 3.6, "A4": 3.0, "E": 2.9} or \
        np.allclose(list(lev.values()), [1.0, 4.0, 3.6, 3.0, 2.9])
    for irrep in ("A1", "A2", "A3", "A4"):
        v = q8_irrep_vector(irrep)
        assert np.allclose(Kq @ v, lev[irrep] * v)
        assert np.allclose(L("-1") @ v, v)
    assert np.allclose(sorted(np.linalg.eigvalsh(Kq)), sorted([1.0, 4.0, 3.6, 3.0] + [2.9] * 4))
    basis = q8_E_basis()
    assert [(n, m, nn) for n, _, m, nn in basis] == [("v1", 1j, -1j), ("v2", -1j, 1j), ("v3", -1j, -1j), ("v4", 1j, 1j)]
    for _, v, _, _ in basis:
        assert np.allclose(Kq @ v, 2.9 * v) and np.allclose(L("-1") @ v, -v)
    Bm = np.array([v for _, v, _, _ in basis])
    assert np.linalg.matrix_rank(Bm) == 4 and np.allclose(Bm.conj() @ Bm.T, 4 * np.eye(4))
    assert q8_mul(q8_mul("j", "i"), q8_inv("j")) == "-i"
    assert np.allclose(Q8_REP2["j"] @ Q8_REP2["j"], -np.eye(2))
    fs = frobenius_schur(Q8_CHARS["E"], Q8_CLASSES, Q8_SIZES, lambda g: q8_mul(g, g))
    assert fs == -1
    fs_d3 = frobenius_schur(D3_CHARS["E"], D3_CLASSES, D3_SIZES,
                            lambda g: {"e": "e", "r": "r2", "r2": "r", "s1": "e", "s2": "e", "s3": "e"}[g])
    assert fs_d3 == 1
    # 實對稱、和 L 交換的 K = Σ c_g R(g)（c_g = c_{g⁻¹}）在 E 上是常數 (c_1 − c_{−1})
    rng = np.random.default_rng(8)
    c = {g: 0.0 for g in Q8_NAMES}
    for g in Q8_NAMES:
        c[g] = c[q8_inv(g)] = c[g] if c[g] else rng.normal()
    Kr = sum(c[g] * R(g) for g in Q8_NAMES)
    assert np.allclose(Kr, Kr.T)
    for _, v, _, _ in basis:
        assert np.allclose(Kr @ v, (c["1"] - c["-1"]) * v)
    return True


def main():
    """`python utils/symosc.py`、`python -m symosc` 與安裝後的 `symosc-check` 都跑這裡。"""
    import sys
    sys.stdout.reconfigure(encoding="utf-8")
    check()
    print("symosc：D3 約化與拍、破缺到 C3／C2、Q8 的左右標籤與四元數型全部驗證通過")


if __name__ == "__main__":
    main()
