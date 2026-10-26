"""
sync_from_harter_sym.py —— 把 harter_sym 專案（Harter《Principles of Symmetry, Dynamics, and Spectroscopy》導讀與動畫）的共用模組複製進 combkit

    python tools/sync_from_harter_sym.py        （在 combkit 資料夾執行）

harter_sym 的 codes/utils/ 是正本（場景 import 的是 utils.X）；這裡的複本給其他計畫用（import combkit.X）。
sym_mat 只依賴 numpy／scipy／sympy；sym_draw 只依賴 manimlib、combkit.draw.styles 與 combkit.draw.walk_draw（絕對 import，兩邊都能用），
所以內容原樣複製，只在檔頭加一行來源註記。
以後 harter_sym 加新章節的模組，把名稱加進 FILES 再執行一次（第 1 章：sym_mat.py、sym_draw.py；第 2 章：sym_abel.py；
第 3 章：sym_nonabel.py；第 4 章：sym_oct.py；第 5 章：sym_rot.py；第 6 章：sym_prod.py；第 7 章：sym_coup.py；第 8 章：sym_dyn.py（相對 import sym_coup、sym_mat、sym_rot）；sym_draw.py 一直是同一份擴充）。sym_abel、sym_nonabel、sym_oct、sym_rot、sym_prod、sym_coup 用相對 import
取 sym_mat、sym_abel、sym_nonabel、sym_oct、sym_rot、sym_prod（`from . import sym_nonabel`，直接執行時退回 `import sym_nonabel`），兩邊都能用；
sym_draw 不 import sym_nonabel、sym_oct、sym_rot、sym_prod、sym_coup（需要的常數另抄一份）。
最後重新產生 API.md。
"""
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = Path(r"C:\mypython\manim\projects\harter_sym\codes\utils")
FILES = [                                   # (harter_sym 的檔名, combkit 裡的位置)
    ("sym_mat.py", ROOT / "combkit" / "sym_mat.py"),
    ("sym_abel.py", ROOT / "combkit" / "sym_abel.py"),
    ("sym_nonabel.py", ROOT / "combkit" / "sym_nonabel.py"),
    ("sym_oct.py", ROOT / "combkit" / "sym_oct.py"),
    ("sym_rot.py", ROOT / "combkit" / "sym_rot.py"),
    ("sym_prod.py", ROOT / "combkit" / "sym_prod.py"),
    ("sym_coup.py", ROOT / "combkit" / "sym_coup.py"),
    ("sym_dyn.py", ROOT / "combkit" / "sym_dyn.py"),
    ("sym_draw.py", ROOT / "combkit" / "draw" / "sym_draw.py"),
]
NOTE = "# 由 combkit/tools/sync_from_harter_sym.py 從 harter_sym 專案的 codes/utils/{name} 複製；請改 harter_sym 那份再同步。\n"


def main():
    for name, dst in FILES:
        text = (SRC / name).read_text(encoding="utf-8")
        assert "from utils" not in text and "import utils" not in text, f"{name} 依賴 harter_sym 的 utils，不能直接複製"
        dst.write_text(NOTE.format(name=name) + text, encoding="utf-8", newline="\n")
        print(f"{name} → {dst.relative_to(ROOT)}")
    subprocess.run([sys.executable, str(ROOT / "tools" / "make_api.py")], check=True)


if __name__ == "__main__":
    main()
