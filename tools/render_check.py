"""
render_check.py —— 低畫質試渲染一個場景，並做一張總覽圖（contact sheet）方便檢查排版

    python tools/render_check.py examples/ex01_dyck_paths.py DyckPathsExample   （--tile 480 縮小總覽圖）
    選項：--sheet 4x3（總覽圖的格數）  --frames 3.5,8（另外抽這幾秒的原尺寸畫面）
          --config 路徑（custom_config.yml；預設用場景檔同一資料夾的那一份）  --hd（改輸出 1080p）

每次在 runs/<類別名稱>/ 這個獨立資料夾裡執行 manimgl（有自己的 latex_cache、公式快取 cache/ 與 renders/），
所以可以同時跑好幾個——ManimGL 預設的公式快取是全機共用的 SQLite，多個程式同時寫入曾經把它弄壞；
並用 strict_render.py 的嚴格模式：公式編譯失敗會直接報錯，不會沿用上一個公式。
需要 ManimGL 與 ffmpeg／ffprobe。印出 STATUS、長度、影片與總覽圖的路徑；失敗時印出記錄的最後幾行。
"""
import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path

TOOLS = Path(__file__).resolve().parent
RUNS = TOOLS.parent / "runs"

ap = argparse.ArgumentParser()
ap.add_argument("scene_file")
ap.add_argument("cls")
ap.add_argument("--sheet", default="4x3")
ap.add_argument("--frames", default="")
ap.add_argument("--config")
ap.add_argument("--hd", action="store_true")
ap.add_argument("--tile", type=int, default=640, help="總覽圖每一格的寬度（像素）；480 可讓總覽圖小約三分之一")
a = ap.parse_args()

scene_file = Path(a.scene_file).resolve()
config = Path(a.config).resolve() if a.config else scene_file.with_name("custom_config.yml")
if not config.exists():
    sys.exit(f"STATUS FAIL: no custom_config.yml next to the scene; give --config")
run = RUNS / a.cls
if (run / "renders").exists():
    shutil.rmtree(run / "renders")
run.mkdir(parents=True, exist_ok=True)
cfg = config.read_text(encoding="ascii")
own_cache = f'  cache: "{(run / "cache").as_posix()}"'
if "\ndirectories:" in "\n" + cfg:
    cfg = cfg.replace("directories:\n", "directories:\n" + own_cache + "\n", 1)
else:
    cfg += "\ndirectories:\n" + own_cache + "\n"
(run / "custom_config.yml").write_text(cfg, encoding="ascii")

log = run / "render.log"
tmp = run / "tmp"                                  # 每次試渲染自己的暫存資料夾：Text 的 .svg 暫存檔放在 %TEMP%，
tmp.mkdir(exist_ok=True)                           # 兩個 render_check 同時跑時才不會互相覆蓋
env = dict(os.environ, PYTHONIOENCODING="utf-8", TEMP=str(tmp), TMP=str(tmp))
with open(log, "w", encoding="utf-8") as fh:
    r = subprocess.run([sys.executable, str(TOOLS / "strict_render.py"), str(scene_file), a.cls,
                        "--hd" if a.hd else "-l"], cwd=run, stdout=fh, stderr=subprocess.STDOUT, env=env)
text = log.read_text(encoding="utf-8", errors="replace").replace("\r", "\n")
videos = sorted((run / "renders").rglob(f"{a.cls}.mp4")) if (run / "renders").exists() else []
if r.returncode != 0 or "Traceback" in text or not videos:
    lines = [ln for ln in text.splitlines() if ln.strip() and "it/s]" not in ln and "it [" not in ln]
    print("STATUS FAIL")
    print("\n".join(lines[-40:]))
    sys.exit(1)
video = videos[0]
dur = float(subprocess.run(["ffprobe", "-v", "error", "-show_entries", "format=duration", "-of", "csv=p=0",
                            str(video)], capture_output=True, text=True).stdout.strip())
cols, rows = map(int, a.sheet.lower().split("x"))
sheet = run / "sheet.png"
subprocess.run(["ffmpeg", "-v", "error", "-y", "-i", str(video), "-vf",
                f"fps={cols * rows / max(dur - 0.15, 0.1)},scale={a.tile}:-1,tile={cols}x{rows}:padding=6:color=white",
                "-frames:v", "1", str(sheet)], check=True)
print(f"STATUS OK  class={a.cls}  duration={dur:.1f}s")
print(f"video={video}")
print(f"sheet={sheet}")
for t in [x.strip() for x in a.frames.split(",") if x.strip()]:
    out = run / f"frame_{t}.png"
    subprocess.run(["ffmpeg", "-v", "error", "-y", "-ss", t, "-i", str(video), "-frames:v", "1", str(out)], check=True)
    print(f"frame={out}")
