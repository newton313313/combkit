"""
make_bundle.py —— 把 combkit 連同目前使用的 ManimGL（3b1b/manim，commit fafa083，1.7.2）打包成一個可以直接安裝的 zip

    python tools/make_bundle.py            產生 ..\\combkit_v<版本>_with_manimgl.zip（已存在就不覆寫）
    python tools/make_bundle.py --force    覆寫同名的 zip

zip 裡最上層是 combkit_with_manim/：
    manimgl/                  C:\\mypython\\manim 裡 git 追蹤的檔案（manimlib、setup.py…；不含 docs、logo、.github、tests）
    combkit/                  與 make_zip.py 相同的檔案（不含試渲染、快取、範例產生的檔）
    docs/                     環境與交接、manim 使用說明、新計畫起步指南（來自 C:\\mypython\\manim）
    先讀我_README.md、install.bat、requirements_pinned.txt
打包前先跑 tests/run_selftests.py（純 Python 部分），失敗就不打包；也確認 ManimGL 的工作樹沒有未提交的修改。
"""
import re
import subprocess
import sys
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MANIM = Path(r"C:\mypython\manim")
TOP = Path("combkit_with_manim")
SKIP_DIRS = {"runs", "__pycache__", "build", "latex_cache", "renders", "media", "cache"}
GENERATED = {".tex", ".pdf", ".aux", ".log", ".mp4", ".png"}
DOCS = ["環境與交接_ENVIRONMENT.md", "manim使用說明.md", "新計畫起步指南_NEW_PROJECT_GUIDE.md"]
PINNED = """addict==2.4.0
appdirs==1.4.4
av==17.1.0
colour==0.1.5
contourpy==1.0.5
diskcache==5.6.3
fonttools==4.25.0
glfw==2.10.2
ipython==8.39.0
isosurfaces==0.1.2
ManimPango==0.6.1
matplotlib==3.7.0
numpy==1.23.5
Pillow==9.4.0
pydub==0.25.1
Pygments==2.17.2
pyperclip==1.11.0
PyWavefront==1.3.3
PyYAML==6.0
rendercanvas==2.7.2
rich==13.7.1
scipy==1.10.0
setuptools==65.6.3
skia-pathops==0.9.2
svgelements==1.9.6
sympy==1.11.1
tqdm==4.64.1
trimesh==5.1.0
typing-extensions
validators==0.35.0
wgpu==0.31.0
"""
INSTALL = r"""@echo off
chcp 65001 >nul
rem 在 combkit_with_manim\ 裡執行；建議先開一個乾淨的 Python 3.10 環境（例：conda create -n manimgl python=3.10）
python -m pip install -r requirements_pinned.txt || goto :fail
python -m pip install --no-deps -e manimgl || goto :fail
python -m pip install --no-deps -e combkit || goto :fail
python combkit\tests\run_selftests.py --draw || goto :fail
echo.
echo 安裝完成，自我測試全部通過。
goto :eof
:fail
echo.
echo 安裝或測試失敗，請看上面的訊息與 先讀我_README.md 的「常見問題」。
exit /b 1
"""


def readme(version, commit):
    return f"""# combkit {version} ＋ ManimGL 1.7.2 安裝包

（本文的相對路徑都以解壓後的 `combkit_with_manim\\` 資料夾為基準。）

## 1. 內容

| 資料夾／檔案 | 內容 |
|---|---|
| `manimgl\\` | ManimGL 1.7.2 原始碼（3Blue1Brown 版，git commit `{commit}`；MIT 授權，見 `manimgl\\LICENSE.md`）。**不是** Manim Community。 |
| `combkit\\` | combkit {version}：組合數學 × ManimGL 的工具包（數學模組、`combkit.draw` 畫圖工具、`examples\\`、`tests\\`、`tools\\`；說明在 `combkit\\README.md`，函式清單在 `combkit\\API.md`） |
| `docs\\` | 原機器的環境與交接說明、ManimGL 使用說明、新計畫起步指南（裡面提到的 `C:\\mypython\\...` 路徑是原機器的位置） |
| `requirements_pinned.txt` | 原機器上能正常運作的套件版本（Python 3.10.9） |
| `install.bat` | 一鍵安裝＋自我測試（Windows） |

不包含：各計畫的資料夾與成品（影片、桌布）、參考書的 PDF 與轉檔（只供原使用者個人學習，請勿散布）。

## 2. 需要先裝好的東西

- **Python 3.10**（原機器是 Anaconda 的 3.10.9；3.11 以上沒測過）。
- **LaTeX**：MiKTeX 或 TeX Live，要有 `latex`／`xelatex`（ManimGL 的 `Tex` 用；預設模板含 amsmath、amssymb、tipa、mathrsfs）。MiKTeX 請開啟「自動安裝缺少的套件」。
- **ffmpeg**（原機器 4.3.1）放在 PATH 上（輸出影片、動畫抽格都要用）。
- 字型：藝術桌布用 Garamond、Corbel（Windows 內建），中文字幕用 Microsoft JhengHei。

## 3. 安裝

在 `combkit_with_manim\\` 裡（建議在新的 conda 環境）執行：

```bat
install.bat
```

它依序做：裝 `requirements_pinned.txt` 的固定版本 → `pip install --no-deps -e manimgl` → `pip install --no-deps -e combkit`
→ `python combkit\\tests\\run_selftests.py --draw`（純 Python 自我測試＋每個 `combkit.draw` 模組都能載入）。
`--no-deps` 是故意的：不讓 pip 把 ManimGL 換成別的版本。之後不要對 manimgl 執行升級或重裝。

## 4. 試跑

```bat
cd combkit\\examples
python ex28_young_wallpaper.py
```

會產生 `ex28_*.png`（兩張 960×540 的 Young 表桌布）。其他範例與執行方式見 `combkit\\README.md` 第 5 節；
做新計畫的步驟見 `docs\\新計畫起步指南_NEW_PROJECT_GUIDE.md`。

## 5. 常見問題（細節見 `docs\\`）

- 每個場景資料夾放一個 `custom_config.yml`，內容必須是純 ASCII（中文 Windows 會用 cp950 讀）；從該資料夾執行 `manimgl`。
- 這版 `manimgl -s`（只存最後一格）不能用；靜態圖用 `combkit.draw.art_kit.render_png()`。
- `Tex` 裡不要用 `\\;`、`\\!`、`\\:`（tipa 讓它們悄悄失效）；`Tex(color=)` 無效，要 `.set_color()`／`.set_fill()`。
- 公式快取壞掉時刪掉場景資料夾的 `latex_cache\\`。
"""


def main():
    version = re.search(r'__version__ = "([^"]+)"', (ROOT / "combkit" / "__init__.py").read_text(encoding="utf-8")).group(1)
    out = ROOT.parent / f"combkit_v{version}_with_manimgl.zip"
    if out.exists() and "--force" not in sys.argv:
        sys.exit(f"{out} 已存在（要覆寫加 --force）")
    if subprocess.run([sys.executable, str(ROOT / "tests" / "run_selftests.py")]).returncode != 0:
        sys.exit("自我測試失敗，沒有打包")
    git = lambda *a: subprocess.run(["git", "-C", str(MANIM), *a], capture_output=True, text=True, check=True).stdout
    commit = git("log", "-1", "--format=%h").strip()
    dirty = [l for l in git("status", "--porcelain", "--untracked-files=no").splitlines() if l.strip()]
    if dirty:
        sys.exit(f"ManimGL 工作樹有未提交的修改，先確認：{dirty[:5]}")
    mfiles = [MANIM / p for p in git("ls-files").splitlines()
              if p.split("/")[0] not in {".github", "docs", "logo", "tests"}]

    def keep(f):
        rel = f.relative_to(ROOT)
        if any(p in SKIP_DIRS or p.endswith(".egg-info") for p in rel.parts[:-1]) or f.suffix in {".pyc", ".zip"}:
            return False
        return not (rel.parts[0] == "examples" and f.suffix in GENERATED)

    cfiles = sorted(f for f in ROOT.rglob("*") if f.is_file() and keep(f))
    with zipfile.ZipFile(out, "w", zipfile.ZIP_DEFLATED) as z:
        for f in mfiles:
            z.write(f, TOP / "manimgl" / f.relative_to(MANIM))
        for f in cfiles:
            z.write(f, TOP / "combkit" / f.relative_to(ROOT))
        for d in DOCS:
            z.write(MANIM / d, TOP / "docs" / d)
        z.writestr(str(TOP / "先讀我_README.md"), readme(version, commit))
        z.writestr(str(TOP / "install.bat"), INSTALL.replace("\n", "\r\n").encode("utf-8"))
        z.writestr(str(TOP / "requirements_pinned.txt"), PINNED)
    n = len(mfiles) + len(cfiles) + len(DOCS) + 3
    print(f"{out}  {n} 個檔案（ManimGL {len(mfiles)}、combkit {len(cfiles)}）  {out.stat().st_size / 1024 / 1024:.1f} MB")


if __name__ == "__main__":
    main()
