"""
make_zip.py —— 把整個 combkit 打包成 zip，放在 combkit 資料夾的上一層

    python tools/make_zip.py            產生 ..\\combkit_v<版本>.zip（版本取自 combkit/__init__.py）
    python tools/make_zip.py --force    同版本的 zip 已存在時也覆寫（預設不覆寫，舊版 zip 會保留）

zip 裡的最上層是 combkit/ 資料夾；不包含試渲染的 runs/、__pycache__、build/、*.egg-info，
以及範例產生的 .tex／.pdf 等檔案。打包前會先跑 tests/run_selftests.py（純 Python 部分），失敗就不打包。
"""
import re
import subprocess
import sys
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SKIP_DIRS = {"runs", "__pycache__", "build", "latex_cache", "renders", "media", "cache"}
SKIP_SUFFIX = {".pyc", ".zip"}
GENERATED = {".tex", ".pdf", ".aux", ".log", ".mp4", ".png"}      # 只在 examples/ 裡略過（範例產生的檔案）

version = re.search(r'__version__ = "([^"]+)"', (ROOT / "combkit" / "__init__.py").read_text(encoding="utf-8")).group(1)
out = ROOT.parent / f"combkit_v{version}.zip"
if out.exists() and "--force" not in sys.argv:
    sys.exit(f"{out} 已存在：請先把 combkit/__init__.py 的 __version__ 加一版，或加 --force 覆寫")

if subprocess.run([sys.executable, str(ROOT / "tests" / "run_selftests.py")]).returncode != 0:
    sys.exit("自我測試失敗，沒有打包")


def keep(f: Path) -> bool:
    rel = f.relative_to(ROOT)
    if any(p in SKIP_DIRS or p.endswith(".egg-info") for p in rel.parts[:-1]) or f.suffix in SKIP_SUFFIX:
        return False
    return not (rel.parts[0] == "examples" and f.suffix in GENERATED)


files = sorted(f for f in ROOT.rglob("*") if f.is_file() and keep(f))
with zipfile.ZipFile(out, "w", zipfile.ZIP_DEFLATED) as z:
    for f in files:
        z.write(f, Path("combkit") / f.relative_to(ROOT))
print(f"{out}  {len(files)} 個檔案  {out.stat().st_size / 1024:.0f} KB")
