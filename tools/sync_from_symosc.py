"""
sync_from_symosc.py —— 把 symosc 專案（D3／C3／C2／Q8 與耦合振盪子的小演示）的共用模組複製進 combkit

    python tools/sync_from_symosc.py                    （在 combkit 資料夾執行）
    python tools/sync_from_symosc.py <symosc 的 codes/utils>   （symosc 不在預設位置時）

symosc 的 codes/utils/ 是正本（場景 import 的是 utils.X）；這裡的複本給其他計畫用（import combkit.X）。
symosc.py 只依賴 numpy；osc_draw.py 只依賴 manimlib 與 combkit.draw.styles（絕對 import，兩邊都能用），
所以內容原樣複製，只在檔頭加一行來源註記。最後重新產生 API.md。
"""
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(r"C:\mypython\manim\projects\symosc\codes\utils")
FILES = [                                   # (symosc 的檔名, combkit 裡的位置)
    ("symosc.py", ROOT / "combkit" / "symosc.py"),
    ("osc_draw.py", ROOT / "combkit" / "draw" / "osc_draw.py"),
]
NOTE = "# 由 combkit/tools/sync_from_symosc.py 從 symosc 專案的 codes/utils/{name} 複製；請改 symosc 那份再同步。\n"


def main():
    for name, dst in FILES:
        text = (SRC / name).read_text(encoding="utf-8")
        assert "from utils" not in text and "import utils" not in text, f"{name} 依賴 symosc 的 utils，不能直接複製"
        dst.write_text(NOTE.format(name=name) + text, encoding="utf-8", newline="\n")
        print(f"{name} → {dst.relative_to(ROOT)}")
    subprocess.run([sys.executable, str(ROOT / "tools" / "make_api.py")], check=True)


if __name__ == "__main__":
    main()
