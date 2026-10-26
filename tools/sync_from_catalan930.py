"""
sync_from_catalan930.py —— 由 catalan930 重新產生 combkit/combkit/（catalan930 本身不會被修改）

catalan930/codes/utils 或 notes/make_catalan_figs.py 改過之後執行：python tools/sync_from_catalan930.py
只覆寫「來自 catalan930 的那些檔案」（請不要在 combkit 裡直接改它們，改 catalan930 那份再同步）；
其他計畫加進 combkit 的模組、以及 __init__.py（版本號）都不會被動到。最後重新產生 API.md。
catalan930/codes/utils 裡新增的模組也會自動帶進來：有 import manimlib 的放 combkit/draw/，其餘放 combkit/。
"""
import ast
import re
from pathlib import Path

ROOT = Path(r"C:\mypython\manim\projects")
UTILS = ROOT / "catalan930" / "codes" / "utils"
FIGS = ROOT / "catalan930" / "notes" / "make_catalan_figs.py"
DST = Path(__file__).resolve().parents[1]
PKG = DST / "combkit"
DRAW = PKG / "draw"

MODS = sorted(f.stem for f in UTILS.glob("*.py") if f.stem != "__init__")
DRAWMODS = [m for m in MODS if "manimlib" in (UTILS / f"{m}.py").read_text(encoding="utf-8")]
PURE = [m for m in MODS if m not in DRAWMODS]


def rewrite(text: str) -> str:
    text = re.sub(rf"\bfrom utils\.({'|'.join(PURE)}) import", r"from combkit.\1 import", text)
    text = re.sub(r"\bfrom utils\.(\w+) import", r"from combkit.draw.\1 import", text)
    # 說明文字裡的 utils.X（例如 "python -m utils.rep_data"）
    text = re.sub(r"\butils\.(\w+)", lambda m: ("combkit." if m.group(1) in PURE else "combkit.draw.") + m.group(1), text)
    assert "utils." not in text.replace("utils.py", ""), "leftover utils reference"
    return text


DRAW.mkdir(parents=True, exist_ok=True)

for name in PURE:
    (PKG / f"{name}.py").write_text(rewrite((UTILS / f"{name}.py").read_text(encoding="utf-8")), encoding="utf-8")
for name in DRAWMODS:
    (DRAW / f"{name}.py").write_text(rewrite((UTILS / f"{name}.py").read_text(encoding="utf-8")), encoding="utf-8")

# tikz.py：make_catalan_figs.py 裡通用的「幾何」與「物件」兩段，加上 picture()
src = FIGS.read_text(encoding="utf-8").splitlines(keepends=True)
start = next(i for i, l in enumerate(src) if l.startswith("# 幾何"))  - 1
stop = next(i for i, l in enumerate(src) if l.startswith("# 圖")) - 1
pic0 = next(i for i, l in enumerate(src) if l.startswith("def picture"))
pic1 = next(i for i in range(pic0 + 1, len(src)) if src[i].startswith("def "))
imports = re.search(r"from utils\.counting import \((.*?)\)", "".join(src), re.S).group(1)
header = '''r"""
tikz.py —— 用 Python 產生 TikZ 程式碼（純 Python，不需要 manim）

做法：先在原點附近建立「圖元」（填色多邊形、線段、點、方塊、圓、文字節點），
算出外框後整體平移到指定位置，再輸出成 TikZ 指令。
卡塔蘭物件（三角剖分、二元樹、Dyck／格子路徑、圓上配對、括號）都由 counting.py 的樹產生，
畫法與 manim 端的 combkit.draw.catalan_draw 相同。

顏色名稱 cA（一般）、cTop（根）、cLeft（左）、cRight（右）、cBad（不合法）要在 LaTeX 導言區定義，例如：
    \\colorlet{cA}{blue!75!black}  \\colorlet{cTop}{orange!90!black}
    \\colorlet{cLeft}{teal!80!black}  \\colorlet{cRight}{purple!80!black}  \\colorlet{cBad}{red!85!black}

來源：catalan930/notes/make_catalan_figs.py（講義第 1 節的圖就是用它產生的）。
"""
from __future__ import annotations

import math

from combkit.counting import (''' + imports + ''')

COLOR = {"root": "cTop", "left": "cLeft", "right": "cRight", None: "cA"}
FILL = {"root": "cTop!38", "left": "cLeft!22", "right": "cRight!20"}


'''
body = "".join(src[start:stop]) + "".join(src[pic0:pic1]).rstrip() + "\n"
(PKG / "tikz.py").write_text(header + body, encoding="utf-8")

if not (PKG / "__init__.py").exists():
  (PKG / "__init__.py").write_text('''"""
combkit —— 組合數學 × ManimGL 動畫工具包（由 catalan930 專案整理出來）

    combkit.counting   卡塔蘭數、廣義卡塔蘭數、Raney、各種卡塔蘭物件、MIT 18.212 補充（純 Python）
    combkit.tableaux   Young 表：列插入、乘積、jeu de taquin、Schur 多項式（純 Python）
    combkit.tikz       用 Python 產生 TikZ 圖（純 Python）
    combkit.draw.*     ManimGL 的畫圖工具與場景基底（需要 ManimGL）

匯入 combkit 本身不會載入 manim；要畫圖時才 import combkit.draw 底下的模組。
"""
__version__ = "1.0.0"
''', encoding="utf-8")
if not (DRAW / "__init__.py").exists():
  (DRAW / "__init__.py").write_text('''"""
combkit.draw —— ManimGL 的畫圖工具（需要已安裝的 ManimGL，本工具包用 3b1b/manim commit fafa083 測過）。
各模組請個別匯入，例如 from combkit.draw.catalan_draw import dyck_path_mob。
"""
''', encoding="utf-8")

DATA_SRC = ROOT / "catalan930" / "codes" / "data"           # rep_data.py 用的固定算例（JSON）
if DATA_SRC.exists():
    (PKG / "data").mkdir(exist_ok=True)
    for f in DATA_SRC.glob("*.json"):
        (PKG / "data" / f.name).write_bytes(f.read_bytes())

import subprocess, sys  # noqa: E401
subprocess.run([sys.executable, str(Path(__file__).with_name("make_api.py"))], check=True)
print(f"synced into {PKG}: {PURE} + draw/{DRAWMODS} + tikz")
