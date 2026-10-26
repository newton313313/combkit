"""
sync_from_walk017.py —— 把 walk017 專案（《Elements of the Random Walk》導讀與動畫）的共用模組複製進 combkit

    python tools/sync_from_walk017.py        （在 combkit 資料夾執行）

walk017 的 codes/utils/ 是正本（場景 import 的是 utils.X）；這裡的複本給其他計畫用（import combkit.X）。
它們只依賴 numpy／scipy／manimlib、combkit.draw.styles 與同套件的其他 walk 模組（相對 import），所以內容原樣複製，只在檔頭加一行來源註記。
以後 walk017 加新章節的模組，把名稱加進 FILES 再執行一次（第 2 章加了 walk_gf.py、lattice_draw.py；第 3 章加了 walk_rec.py；第 4 章加了 walk_bc.py；第 5 章加了 walk_var.py；第 6 章加了 walk_shape.py；第 7 章加了 walk_saw.py；第 8 章加了 walk_scale.py，它用相對 import 引用同套件的 walk_gf、walk_saw；第 9 章加了 walk_crit.py，用相對 import 引用 walk_saw、walk_scale；第 10 章加了 walk_cond.py，用相對 import 引用 walk_saw、walk_crit、walk_rec、walk_gf；第 11 章加了 walk_rg.py，用相對 import 引用 walk_saw、walk_crit；第 12 章加了 walk_eps.py，用相對 import 引用 walk_saw、walk_crit、walk_rg）。
最後重新產生 API.md。
"""
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = Path(r"C:\mypython\manim\projects\walk017\codes\utils")
FILES = [                                   # (walk017 的檔名, combkit 裡的位置)
    ("walk1d.py", ROOT / "combkit" / "walk1d.py"),
    ("walk_draw.py", ROOT / "combkit" / "draw" / "walk_draw.py"),
    ("walk_gf.py", ROOT / "combkit" / "walk_gf.py"),
    ("lattice_draw.py", ROOT / "combkit" / "draw" / "lattice_draw.py"),
    ("walk_rec.py", ROOT / "combkit" / "walk_rec.py"),
    ("walk_bc.py", ROOT / "combkit" / "walk_bc.py"),
    ("walk_var.py", ROOT / "combkit" / "walk_var.py"),
    ("walk_shape.py", ROOT / "combkit" / "walk_shape.py"),
    ("walk_saw.py", ROOT / "combkit" / "walk_saw.py"),
    ("walk_scale.py", ROOT / "combkit" / "walk_scale.py"),
    ("walk_crit.py", ROOT / "combkit" / "walk_crit.py"),
    ("walk_cond.py", ROOT / "combkit" / "walk_cond.py"),
    ("walk_rg.py", ROOT / "combkit" / "walk_rg.py"),
    ("walk_eps.py", ROOT / "combkit" / "walk_eps.py"),
]
NOTE = "# 由 combkit/tools/sync_from_walk017.py 從 walk017 專案的 codes/utils/{name} 複製；請改 walk017 那份再同步。\n"


def main():
    for name, dst in FILES:
        text = (SRC / name).read_text(encoding="utf-8")
        assert "from utils" not in text and "import utils" not in text, f"{name} 依賴 walk017 的 utils，不能直接複製"
        dst.write_text(NOTE.format(name=name) + text, encoding="utf-8", newline="\n")
        print(f"{name} → {dst.relative_to(ROOT)}")
    subprocess.run([sys.executable, str(ROOT / "tools" / "make_api.py")], check=True)


if __name__ == "__main__":
    main()
