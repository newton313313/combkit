"""
run_selftests.py —— 檢查 combkit 是否正常

    python tests/run_selftests.py           純 Python 部分：counting、tableaux 的自我測試，tikz 能產生圖
    python tests/run_selftests.py --draw    另外載入 combkit.draw 的每個模組（需要 ManimGL）

會自動找出 combkit/ 裡所有有 `if __name__ == "__main__":` 自我測試的模組，以及 combkit/draw/ 的所有模組，
所以新加入的模組不用改這個檔案；只有想檢查輸出裡某句話時，才在 EXPECT 加一行。
"""
import importlib
import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
env = dict(os.environ, PYTHONIOENCODING="utf-8", PYTHONPATH=str(ROOT))
EXPECT = {
    "combkit.counting": ["全部一致", "Raney 公式證明的每一步都通過", "卡塔蘭基礎模組：所有物件的數目與一一對應都通過",
                         "MIT 補充模組：所有計數與一一對應都通過"],
    "combkit.raney": ["Raney (1960)：字、循環引理、L(n;a)、L*、Lagrange 反演、第 6–7 節全部驗證通過",
                      "Raney 引理（整數序列）、加強版、出棧、摸牌全部驗證通過"],
    "combkit.raney_examples": ["Raney 例題產生器：均勻性、隨機不交叉分法、運算式、例題集全部通過"],
    "combkit.walk1d": ["walk1d（第 1 章）：計數、動差、Stirling、高斯近似、遞迴、傅立葉、生成函數、最陡下降全部驗證通過"],
    "combkit.walk_gf": ["walk_gf（第 2 章）：Legendre、統計力學、格子計數、傅立葉、高斯極限、普適性、Cauchy、擴散、補充 1–3 全部驗證通過"],
    "combkit.walk_rec": ["walk_rec（第 3 章）：到訪次數、一維／二維／三維的回家、Watson 積分、足跡全部驗證通過"],
    "combkit.walk_bc": ["walk_bc（第 4 章）：反射原理、兩牆之間、球與立方體、穩定態、點源與鏡像、吸收球、Berg–Purcell 全部驗證通過"],
    "combkit.walk_var": ["walk_var（第 5 章）：偏向漫步、有慣性的漫步、電報方程、二三維的跑與轉、連續時間漫步全部驗證通過"],
    "combkit.walk_shape": ["walk_shape（第 6 章）：迴轉半徑張量、Kramers 恆等式、非球度、無限多維的樓梯、Wick 配對、1/d 展開、Fixman 分布全部驗證通過"],
    "combkit.walk_saw": ["walk_saw（第 7 章）：路徑和的切段、高斯摺積、自迴避漫步的精確列舉、1−f 展開（排容）、弧線圖與 Dyson 方程、冪次數算、平均場、Flory 指數全部驗證通過"],
    "combkit.walk_scale": ["walk_scale（第 8 章）：長尾巴步伐的 χ(q) 與寬度 N^(1/y)、尾巴 r^(−5+δ) 與換變數、生成函數的標度塌縮與實空間 ξ、最陡下降與 Fisher 指數、γ = (2−η)ν 與比值法、重新標度與上臨界維度 4、Flory 指數全部驗證通過"],
    "combkit.walk_crit": ["walk_crit（第 9 章）：平均場與 Onsager 的 β、小 Ising 的精確列舉與 Metropolis 翻轉時間、一維 Ising 的 G、ξ、χ、標度律與 Widom 形式、Landau 雙井與 (9.34)、(z_c−z)^(−s) 的係數與 ⟨R²⟩ ∝ N^(2ν)、輔助場與高斯積分、Wick、Trotter、相交權重與 n → 0 的自迴避全部驗證通過"],
    "combkit.walk_cond": ["walk_cond（第 10 章）：環面上的自迴避列舉（三種算法、與無限平面並排）、平均場積分與 n → 0、註 3 的對照、自旋波矩陣的特徵值與逆矩陣、牆邊的 tanh、相關洞的 Coulomb 形與 Ewald 常數、erfc 形的步數、Redner–Reynolds 的 ρ 與指數全部驗證通過"],
    "combkit.walk_rg": ["walk_rg（第 11 章）：Koch 曲線與 Cantor 集的質量－半徑＝盒子計數＝複本公式、普通漫步與自迴避的碎形維度、Landau 與 Gaussian 的標度塌縮、一維 Gaussian 鏈的抽格（完成平方＝數值積分、環的行列式）、一維 Ising 的抽格（暴力＝轉移矩陣＝連乘）、x 的遞迴與 ξ、動量殼與連續流（解析＝數值）、(11.87) 的驗證、線性化與指數公式（Fraction＝Widom）全部驗證通過"],
    "combkit.walk_eps": ["walk_eps（第 12 章）：兩變數玩具（閉式＝數值積分＝擬合）、殼層體積（指數＝線性＝格點計數）、重新標度、配對列舉＝組合公式（(12:S-4)(12:S-6)、f₃ = 6）、連結叢集（四種算法、三變數 det＝求積＝環的和）、O(n) 縮併 2n+4 與 4(n+8)、不動點（Fraction＝fsolve）、線性化（閉式＝numpy＝差分＝流的時間）、u(ℓ) 與 t(ℓ)（解析＝數值）、ν 的一階公式與 Flory 展開、(12.86)＝(12.39)、冪次數算（Fraction＝數值積分）、全書回顧的數字全部驗證通過"],
    "combkit.sym_mat": ["sym_mat（第 1 章）：分析器與 Axiom 4（密度矩陣＝振幅求和）、么正性、bra/ket 與換基底、兩態拍頻（閉式＝投影算符＝expm＝ODE）、久期方程與 Hamilton–Cayley、投影算符（精確）、M⁵⁰／√M／e^Mt（譜分解＝直接）、最小方程、對易矩陣、伴隨矩陣、散射與密度矩陣、耦合振子（解析＝數值）、d/dx 與 δ 序列、預解式圍道積分＝投影算符、習題 1.1.1–1.2.8 全部核對通過"],
    "combkit.sym_abel": ["sym_abel（第 2 章）：D₂／C_n 群表與公理（含「不是群」的表）", "習題 2.3.3–2.12.3 全部核對通過"],
    "combkit.sym_nonabel": ["sym_nonabel（第 3 章）：C3v／D3／D4／C4v／D6／D3h／D3d／D2h／D4h／Q 的乘法表", "習題 3.1.1–3.5.6 全部核對通過"],
    "combkit.sym_oct": ["sym_oct（第 4 章）：O／Oh／Td／T／Th 由轉動矩陣產生", "習題 4.1.1–4.4.2 全部核對通過"],
    "combkit.symosc": ["symosc：D3 約化與拍、破缺到 C3／C2、Q8 的左右標籤與四元數型全部驗證通過"],
    "combkit.sym_rot": ["sym_rot（第 5 章）：O3／O2 與 D∞h、Euler ↔ Darboux ↔ 四元數", "習題 5.1.1–5.6.2 全部核對通過"],
    "combkit.sym_prod": ["sym_prod（第 6 章）：八面體上兩粒子的外積群與對角子群", "習題 6.1–6.6 全部核對通過"],
    "combkit.sym_coup": ["sym_coup（第 7 章）：R3 耦合係數", "混成與 Hückel 全部核對通過"],
    "combkit.sym_dyn": ["sym_dyn（第 8 章）：包絡與 Legendre", "Jaynes–Cummings 與 Mollow 全部核對通過"],
    "combkit.young_extras": ["young_extras（youngart 桌布）：patience、滑動路徑、擠壓路徑、影子線、快速 R-S-K、pipe 水管、雙曲面直線、Bézout 交點全部驗證通過"],
    "combkit.tableaux": ["Fulton 第 1 章：插入、乘積、jeu de taquin 全部驗證通過",
                         "Fulton 第 2–4 章：Knuth 等價、L(w,k)、R-S-K、矩陣球全部驗證通過",
                         "Fulton 第 5–6 章：Littlewood–Richardson 規則、對稱函數全部驗證通過"],
}
PURE = sorted(f"combkit.{f.stem}" for f in (ROOT / "combkit").glob("*.py")
              if '__name__ == "__main__"' in f.read_text(encoding="utf-8"))
bad = 0
for mod in PURE:
    needles = EXPECT.get(mod, [])
    r = subprocess.run([sys.executable, "-m", mod], cwd=ROOT, env=env, capture_output=True, text=True,
                       encoding="utf-8", errors="replace")
    ok = r.returncode == 0 and all(n in r.stdout for n in needles)
    print(f"[{'OK' if ok else 'FAIL'}] python -m {mod}")
    bad += not ok

sys.path.insert(0, str(ROOT))
from combkit.counting import binary_trees, tree_to_triangulation  # noqa: E402
from combkit.tikz import picture, place, tri_prims  # noqa: E402

tex = picture(place(tri_prims(5, tree_to_triangulation(binary_trees(3)[0], 0, 4), 0.6), 0, 0))
ok = tex.startswith(r"\begin{tikzpicture}") and tex.rstrip().endswith(r"\end{tikzpicture}")
print(f"[{'OK' if ok else 'FAIL'}] combkit.tikz")
bad += not ok

if "--draw" in sys.argv:
    sys.argv = sys.argv[:1]                       # manimlib 在 import 時會解析命令列
    for name in sorted(f.stem for f in (ROOT / "combkit" / "draw").glob("*.py") if f.stem != "__init__"):
        try:
            importlib.import_module(f"combkit.draw.{name}")
            print(f"[OK] import combkit.draw.{name}")
        except Exception as err:                  # noqa: BLE001
            print(f"[FAIL] import combkit.draw.{name}: {err}")
            bad += 1

print("全部通過" if not bad else f"{bad} 項失敗")
sys.exit(1 if bad else 0)
