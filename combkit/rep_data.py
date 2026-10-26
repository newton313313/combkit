"""
rep_data.py —— Fulton 第 7 章（場景 F28–F41）的固定算例與小工具（純 Python）

    load(key)                 讀 codes/data/fulton7_math.json 的一段（複製自 fulton7_agent_pack 的 MATH_DATA.json，已驗證）
    act(p, rows)              置換作用在編號上：σ·T 把 i 換成 σ(i)（置換用一行記法 tuple，1 起算）
    compose(p, q)             p∘q（先 q 後 p）；sign(p)、cycles(p)、cycle_type(p)
    tabloid_key(rows)         列 tabloid 的標準形：各列排序後的 tuple（例 ((1,2),(3,))）；tabloid_str → "12/3"
    polytabloid(rows)         v_T = b_T·{T} = {tabloid_key: 係數}
    specht21(vec)             λ=(2,1) 的向量（{tabloid: 係數}）寫成 (v_a, v_b) 座標
    rep_matrix21(p)           S^(2,1) 在 (v_a, v_b) 下的矩陣（第 j 欄 = 第 j 個基底的像），list of rows
    gr_mul(a, b)              群環乘法（元素 = {perm tuple: 係數}）；symmetrizers(rows) → (a_T, b_T, c_T)

直接執行：python -m combkit.rep_data（在 codes/ 裡）會用 tableaux.py 獨立重算並核對資料。
資料慣例同 Fulton：列 = row、行 = column；矩陣以 rows 存，第 j 欄是第 j 個基底向量的像。
"""
from __future__ import annotations

import json
from collections import Counter
from itertools import permutations, product
from pathlib import Path

_HERE = Path(__file__).resolve().parent
_DATA = next(p for p in (_HERE / "data" / "fulton7_math.json", _HERE.parent / "data" / "fulton7_math.json")
             if p.exists())                      # combkit 裡放在套件內的 data/；catalan930 裡在 codes/data/
_CACHE = None


def load(key: str):
    global _CACHE
    if _CACHE is None:
        _CACHE = json.loads(_DATA.read_text(encoding="utf-8"))
    return _CACHE[key]


# ---------------------------------------------------------------- 置換
def compose(p, q):
    """(p∘q)(i) = p(q(i))：先 q 後 p。"""
    return tuple(p[q[i] - 1] for i in range(len(p)))


def sign(p) -> int:
    return -1 if sum(p[i] > p[j] for i in range(len(p)) for j in range(i + 1, len(p))) % 2 else 1


def cycles(p, keep_fixed=True):
    seen, out = set(), []
    for i in range(1, len(p) + 1):
        if i not in seen:
            c, j = [], i
            while j not in seen:
                seen.add(j)
                c.append(j)
                j = p[j - 1]
            if keep_fixed or len(c) > 1:
                out.append(tuple(c))
    return out


def cycle_type(p):
    return tuple(sorted((len(c) for c in cycles(p)), reverse=True))


def transposition(a, b, n):
    p = list(range(1, n + 1))
    p[a - 1], p[b - 1] = b, a
    return tuple(p)


# ---------------------------------------------------------------- 編號與 tabloid
def act(p, rows):
    return tuple(tuple(p[x - 1] for x in r) for r in rows)


def tabloid_key(rows):
    return tuple(tuple(sorted(r)) for r in rows)


def tabloid_str(key) -> str:
    return "/".join("".join(map(str, r)) for r in key)


def columns(rows):
    return [[r[j] for r in rows if len(r) > j] for j in range(len(rows[0]))]


def _stabilizer(groups, n):
    out = []
    for orders in product(*(list(permutations(g)) for g in groups)):
        p = list(range(1, n + 1))
        for g, o in zip(groups, orders):
            for a, b in zip(g, o):
                p[a - 1] = b
        out.append(tuple(p))
    return out


def row_group(rows):
    return _stabilizer([list(r) for r in rows], sum(map(len, rows)))


def col_group(rows):
    return _stabilizer(columns(rows), sum(map(len, rows)))


def polytabloid(rows):
    """v_T = Σ_{q∈C(T)} sgn(q) {qT}，回傳 {tabloid_key: 係數}（略去 0）。"""
    out = Counter()
    for q in col_group(rows):
        out[tabloid_key(act(q, rows))] += sign(q)
    return {k: v for k, v in out.items() if v}


def act_vec(p, vec):
    out = Counter()
    for key, c in vec.items():
        out[tabloid_key(act(p, key))] += c
    return {k: v for k, v in out.items() if v}


# λ = (2,1)：t_i = 第二列是 i 的 tabloid；v_a = t3 - t1、v_b = t2 - t1
T_A, T_B = ((1, 2), (3,)), ((1, 3), (2,))


def t21(i):
    rest = tuple(sorted({1, 2, 3} - {i}))
    return (rest, (i,))


def specht21(vec):
    """零和向量 x1 t1 + x2 t2 + x3 t3 = a v_a + b v_b：a = x3、b = x2（且 x1 = -a - b）。"""
    x = [vec.get(t21(i), 0) for i in (1, 2, 3)]
    assert sum(x) == 0, "不在 S^(2,1) 裡"
    return (x[2], x[1])


def rep_matrix21(p):
    cols = [specht21(act_vec(p, polytabloid(T))) for T in (T_A, T_B)]
    return [[cols[0][0], cols[1][0]], [cols[0][1], cols[1][1]]]


# ---------------------------------------------------------------- 群環
def gr_mul(a, b):
    out = Counter()
    for p, x in a.items():
        for q, y in b.items():
            out[compose(p, q)] += x * y
    return {k: v for k, v in out.items() if v}


def symmetrizers(rows):
    a = {p: 1 for p in row_group(rows)}
    b = {q: sign(q) for q in col_group(rows)}
    return a, b, gr_mul(b, a)


# ---------------------------------------------------------------- 自我檢查
def check():
    from math import factorial

    from combkit.counting import hook_product
    from combkit.tableaux import (character_table, conjugate, kostka, lr_coefficient, standard_tableaux,
                                z_lambda)
    s21 = load("shape21")
    for item in s21["all_numberings"]:
        rows = tuple(map(tuple, item["rows"]))
        v = polytabloid(rows)
        assert {tabloid_str(k).replace("", ""): c for k, c in v.items()} == \
            {k.replace(",", ""): c for k, c in item["v_terms"].items()}
        assert list(specht21(v)) == item["specht_coordinates"]
    m = load("matrices21")
    assert rep_matrix21((2, 1, 3)) == m["rho12"] and rep_matrix21((1, 3, 2)) == m["rho23"]
    cs = load("character_S3")
    _, X = character_table(3)                  # 列 (3),(2,1),(1^3)；行 μ 同序
    assert [list(reversed(r)) for r in X] == [cs["specht"]["3"], cs["specht"]["21"], cs["specht"]["111"]]
    for mu, size, z in zip(cs["types"], cs["sizes"], cs["z"]):
        assert z_lambda(tuple(mu)) == z and factorial(3) // z == size
    tr = {(1, 1, 1): (1, 2, 3), (2, 1): (2, 1, 3), (3,): (2, 3, 1)}
    assert [sum(rep_matrix21(tr[tuple(mu)])[i][i] for i in range(2)) for mu in cs["types"]] == cs["specht"]["21"]
    s32 = load("shape32")
    assert [list(map(list, T)) for T in standard_tableaux((3, 2))] == s32["standard"]
    assert tuple(s32["transpose"]) == conjugate((3, 2)) and hook_product((3, 2)) == s32["hook_product"] == 24
    T32 = ((1, 2, 3), (4, 5))
    assert len(row_group(T32)) == s32["row_group_size"] and len(col_group(T32)) == s32["column_group_size"]
    for comp in load("young221")["components"]:
        assert kostka(tuple(comp["lambda"]), (2, 2, 1)) == comp["K"] == len(comp["semistandard"])
        assert len(standard_tableaux(tuple(comp["lambda"]))) == comp["f"]
    assert sum(c["K"] * c["f"] for c in load("young221")["components"]) == load("young221")["dim_M"] == 30
    pc = load("permutation_char32")
    fixed = {tabloid_key(r) for r in (((1, 2, 5), (3, 4)), ((3, 4, 5), (1, 2)))}
    assert {tabloid_key(tuple(map(tuple, r))) for r in pc["fixed_tabloids"]} == fixed and pc["character"] == 2
    ys = load("young_symmetrizer21")
    a, b, c = symmetrizers(T_A)
    assert {tuple(t["perm"]): t["coef"] for t in ys["c_terms"]} == c
    assert gr_mul(c, c) == {p: 3 * v for p, v in c.items()} and ys["n_lambda"] == 3
    lr = load("LR_example")
    assert lr_coefficient(tuple(lr["lambda"]), tuple(lr["mu"]), tuple(lr["nu"])) == lr["multiplicity"] == 2
    br = load("branching")
    assert [len(standard_tableaux(tuple(l))) for l in br["add21"]] == br["f_add"] and sum(br["f_add"]) == 8
    return True


if __name__ == "__main__":
    check()
    print("rep_data：固定算例與 tableaux.py 的獨立重算全部一致")
