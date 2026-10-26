# ManimGL 工作環境與檔案結構（給繼任者）

整理日期：2026-10-09，最後更新 2026-10-13。本文件說明這台電腦上 manim 動畫工作的**環境、依賴、檔案結構與注意事項**。
各專案的內容細節請看各自的說明（見第 4 節）。
本文件的正本在 `C:\mypython\manim`；`projects/catalan_portable_105`、`projects/catalan_portable_104`、`projects/catalan_portable_103` 與 `projects/combkit` 裡各有一份副本，
副本裡的路徑仍然指開發機上的位置。

**路徑的寫法：** 本文件中沒有寫磁碟機的相對路徑，一律相對於 `C:\mypython\manim`（本文件所在的資料夾）。
例如 `projects/combkit` 就是 `C:\mypython\manim\projects\combkit`。

---

## 1. 一句話總覽

- 用的是 **ManimGL**（3Blue1Brown 版，`manimgl` 指令、`from manimlib import *`）**v1.7.2**，原始碼就在 `C:\mypython\manim`，以**可編輯模式**安裝。
- **不是** Manim Community（`manim` 指令、`from manim import *`）。網路上大多數教學是社群版，語法不通用。
- 作業系統：Windows 11 Home，**繁體中文**（系統編碼 cp950）。這一點會影響設定檔與中文字的處理，見第 6 節。

---

## 2. 軟體與版本（2026-10-09 實測）

| 項目 | 版本／位置 | 備註 |
|---|---|---|
| Python | 3.10.9，Anaconda base：`C:\Users\User\anaconda3\python.exe` | conda 23.3.1、pip 22.3.1 |
| ManimGL | 1.7.2，git commit **fafa083**（2026-09-09，"Add VideoMobject and Sprite (#2521)"），remote `https://github.com/3b1b/manim.git` | `pip install -e .` 可編輯安裝；`import manimlib` 會載入 `C:\mypython\manim\manimlib` |
| 繪圖後端 | `wgpu` 0.31.0、`rendercanvas` 2.7.2、`glfw` 2.10.2 | 這個版本的 ManimGL 已不用 moderngl |
| 其他主要套件 | numpy 1.23.5、scipy 1.10.0、sympy 1.11.1、Pillow 9.4.0、matplotlib 3.7.0、ManimPango 0.6.1、skia-pathops 0.9.2、svgelements 1.9.6、diskcache 5.6.3、PyYAML 6.0、`av` 17.1.0 | `av` 不在 setup.cfg 裡，要另外由 `requirements.txt` 安裝 |
| ffmpeg／ffprobe | 4.3.1（conda-forge），`C:\Users\User\anaconda3\Library\bin\` | |
| LaTeX | MiKTeX 25.12（使用者安裝）：`C:\Users\User\AppData\Local\Programs\MiKTeX\miktex\bin\x64\` | `latex`（pdfTeX 4.23）、`xelatex`（XeTeX 4.16）、`dvisvgm` 3.4.3 |
| MiKTeX 額外套件 | **`cm-super`**（必要，否則公式出現 "SVG has no content"）；講義另需 `xeCJK`、`tcolorbox`、`tikz` 等（MiKTeX 第一次用到時會自動安裝） | MiKTeX 會警告「尚未檢查更新」，可忽略 |
| 中文字型 | 微軟正黑體 `Microsoft JhengHei`（`C:\Windows\Fonts\msjh*.ttc`） | 動畫字幕與講義都用它 |
| git | 2.44.0.windows.1 | |
| PyMuPDF | 1.28.2 | 只用在 PDF 轉 Markdown（`projects/ref/*/convert_to_md.py`）與講義預覽 |
| 講義用的 LaTeX 套件 | `xeCJK`、`tcolorbox`、`tikz`、`ytableau`、`booktabs`、`hyperref`（MiKTeX 自動安裝） | 這版 `ytableau` 沒有 `enableskew` 選項；空格子要寫 `{}` |

> **不要升級或重裝 ManimGL**（也不要 `pip install manimgl` 或 `git pull`）。所有場景都是用 commit fafa083 測試的。
> 因此 combkit 的 `pyproject.toml` 刻意沒有把 manimgl 列為依賴。

### 2.1 在新電腦上重建環境

1. 安裝 Anaconda（Python 3.10）、git、MiKTeX，再執行 `miktex packages install cm-super`。
2. 取得同一版 ManimGL，二擇一：
   - 用可攜包裡附的原始碼：`projects/catalan_portable_105/third_party/manim`，版本記在 `MANIM_VERSION.txt`；
   - 或者從網路取得：
     ```bash
     git clone https://github.com/3b1b/manim.git
     ```
     ```bash
     git -C manim checkout fafa083
     ```
3. 在 manim 原始碼資料夾裡安裝：
   ```bash
   pip install -e .
   ```
   ```bash
   pip install -r requirements.txt
   ```
4. 安裝 ffmpeg：
   ```bash
   conda install -c conda-forge ffmpeg
   ```
5. 檢查環境。在 `projects/catalan930/codes` 裡執行：
   ```bash
   python check_env.py --render
   ```
   它會檢查 Python、ManimGL、ffmpeg、LaTeX、cm-super、字型、設定檔，跑數學自我測試，並試渲染一段。

可攜包 `projects/catalan_portable_105` 裡有一鍵安裝腳本 `install_windows.bat`，`codes/render_all.bat` 可以一鍵輸出全部影片（104 起不附影片）。

原始安裝過程記錄在 `install.log`，ManimGL 的基本用法見 `manim使用說明.md`。可攜包的 `setup_records/` 也各有一份。

---

## 3. 檔案結構

```text
C:\mypython\manim\                      ← ManimGL 原始碼（git：3b1b/manim，master，commit fafa083）
├── manimlib\                           ManimGL 本體（可編輯安裝；不要修改）
├── example_scenes.py, docs\, ...       官方原始內容
├── install.log                         安裝記錄（2026-09-30）
├── manim使用說明.md                     ManimGL 簡短使用說明（中文）
├── 環境與交接_ENVIRONMENT.md            ← 本文件
├── videos\                             官方範例的試渲染影片
├── latex_cache\                        在這一層執行 manimgl 時產生的公式暫存（可刪）
└── projects\                           ← 我們自己的工作都在這裡（外層 manim 的 git 不追蹤）
    ├── catalan930\                     主專案：卡塔蘭數、Chu 1987、MIT 18.212、Fulton《Young Tableaux》的動畫與講義（2026-10-13 起是獨立的 git repo）
    │   ├── 交接說明_HANDOFF.md          ★ 這個專案的權威說明（內容、指令、場景與影片清單、踩過的坑）
    │   ├── codes\                      動畫程式（必須在這裡執行 manimgl）
    │   │   ├── custom_config.yml       ManimGL 設定（只能用英文）
    │   │   ├── scenes\                 各場景：00a–00m、01–09、08r、M01–M14、F01–F65（Fulton 第 1–10 章與附錄）
    │   │   ├── data\                   fulton7_math.json（第 7 章的固定算例，由 utils/rep_data.py 讀取）
    │   │   ├── utils\                  共用模組：數學（counting、tableaux、rep_data）＋畫圖（相依圖見 codes/README.md）
    │   │   ├── main.py                 每個場景的同名子類別，以及十種連續播放的 demo（含 FultonDemo～FultonDemo5）
    │   │   ├── render_all.py           批次輸出影片
    │   │   ├── check_env.py            環境檢查
    │   │   └── renders\                輸出的影片（共 92 支；不要刪除或覆寫）
    │   ├── notes\                      講義：chu1987_導讀（Chu＋卡塔蘭＋MIT）、young_tableaux_導讀（Fulton 第 1–10 章與附錄 A、B，內文 fulton/），
    │   │                               出圖兼數字檢查程式、各章計畫表 fulton_plan_*.md
    │   └── ref\                        參考文獻 PDF 與抄錄（僅供個人研讀，不可散布）
    ├── catalan_portable_101.zip, catalan_portable_102\, catalan_portable_103\, catalan_portable_104\、105\（及 .zip）
    │                                   可攜版快照（105 最新：全部 87 個場景含 Fulton F01–F45、Raney R01–R06，ManimGL 原始碼與安裝腳本，不含影片）
    ├── ref\                            Fulton 各章 PDF 與 md 轉檔（ref/fulton/md）、Sagan、整合參考與第 7–8 章分鏡大綱（僅供個人研讀）
    ├── fulton7_agent_pack\             使用者提供的第 7 章工作包（固定算例 MATH_DATA.json、場景卡）
    ├── combkit\                        可重複使用的工具包（數學＋畫圖＋TikZ），說明見它的 README.md
    └── combkit_v1.0.0.zip … combkit_v1.8.0.zip   工具包的打包檔（每次加入新演示升一版；舊版保留）
```

本機另有一個 ManimGL 的公式快取 `%LOCALAPPDATA%\manim\manim\Cache`，是全機共用的 SQLite 資料庫。
它的旁邊有一個 `Cache_corrupt_20261006`，是損毀後改名留下的舊快取，可以刪掉。

---

## 4. 該先讀哪些文件

| 想做的事 | 讀這個 |
|---|---|
| 了解 catalan930 的全部內容、指令、場景清單 | `projects/catalan930/交接說明_HANDOFF.md`（最權威）、`projects/catalan930/codes/README.md` |
| 在別的計畫重複使用程式 | `projects/combkit/README.md`（三層架構「先算再畫」、模組說明、範例、§10 加入新演示的步驟）、`projects/combkit/API.md` |
| ManimGL 基本用法（最小範例、常用指令、常見問題） | `manim使用說明.md` |
| 安裝時遇過的問題 | `install.log` |

---

## 5. 日常操作

- **輸出影片：** 在 `projects/catalan930/codes` 裡執行，例如：
  ```bash
  python render_all.py --fulton
  ```
  也可以直接執行 `manimgl scenes/xxx.py 類別名 -w`。影片會寫到該資料夾的 `renders/`。
- **預覽（互動視窗）：** `manimgl 檔案.py 類別名`。這一版的按鍵和官方教學不同：
  - Ctrl+Q 結束；
  - 滾輪縮放；
  - 按住 `f` 移動畫面；
  - `r` 重設視角。
- **公式出錯時清快取：** `manimgl --clear-cache`。
- **試渲染並檢查排版：** 在 `projects/combkit` 裡執行：
  ```bash
  python tools/render_check.py <場景檔> <類別>
  ```
  它產生低畫質影片與總覽圖，每次用獨立資料夾與獨立快取，所以可以同時跑好幾個。
- **編譯講義：** 在 `projects/catalan930/notes` 裡執行兩次 `xelatex chu1987_導讀.tex`；Young 表講義先執行 `python fulton/make_figs.py`（檢查數字並產生圖），再執行兩次 `xelatex young_tableaux_導讀.tex`。
- **測試：** 執行 `python -m utils.counting` 與 `python -m utils.tableaux`（在 `codes` 裡），或在 `projects/combkit` 裡執行：
  ```bash
  python tests/run_selftests.py --draw
  ```
- **工具包同步與打包：** catalan930 的 `utils` 改過之後，在 `projects/combkit` 裡執行：
  ```bash
  python tools/sync_from_catalan930.py
  ```
  然後打包：
  ```bash
  python tools/make_zip.py
  ```

---

## 6. 這個環境特有的注意事項

1. **`custom_config.yml` 只能寫英文（ASCII）。** ManimGL 用系統編碼 cp950 讀它，有中文就會壞。它必須放在執行 `manimgl` 的資料夾。
2. **中文不要放進 `Tex`。** 中文用 `Text`，或用工具包的 `zh()`、`mixed()`。字型是 `Microsoft JhengHei`。
3. **舊公式快取的錯誤。** ManimGL 編譯公式失敗時可能不報錯，直接沿用上一個公式，產生錯的畫面。
   - 常見原因：Tex 開頭或結尾有 `\;`、`\,`；或者用了 `t2c` 上色卻同時寫 `\left( … \right)`。
   - 專案用「嚴格模式」防止這種情況：每次編譯前刪掉 `latex_cache/working.dvi`。
   - 嚴格模式的程式在 `projects/combkit/tools/strict_render.py`。
4. **平行渲染要分開資料夾與快取：**
   - 同一資料夾同時跑兩個 `manimgl`，`latex_cache/working.tex` 會互相覆蓋；
   - 全機共用的 SQLite 快取曾因此損毀，錯誤訊息是 `database disk image is malformed`；
   - 每個程式要用自己的 `directories: cache:` 設定，或者直接用 `render_check.py`；
   - 快取壞掉時把那個資料夾改名或刪掉，ManimGL 會自動重建；
   - 正式輸出（`render_all.py`）一次只跑一個：`Text` 的暫存 .svg 放在共用的 `%TEMP%`，兩個同時跑會互相覆蓋。
5. **用 Bash 執行指令時：**
   - Git Bash 會吃掉反斜線，含 `\` 的檔案請用編輯器寫，不要用 sed 或 echo；
   - `.bat` 要寫成 `.\x.bat` 執行；
   - 呼叫 Windows 指令時加 `MSYS_NO_PATHCONV=1`。
6. **輸出中文到終端機時**先設定 `PYTHONIOENCODING=utf-8`，避免 cp950 編碼錯誤。
7. 其他 ManimGL 的陷阱見 `projects/catalan930/交接說明_HANDOFF.md` 第 7 節，共 18 點，例如：
   - `Tex(...)[i]` 取到的是第 i 個字形，不是第 i 段；
   - 自訂屬性不要和 ManimGL 的屬性撞名；
   - 用 `play_sequence` 連續播放時，輔助函式要寫在模組層級。

---

## 7. 工作慣例（使用者的要求）

- 文件、字幕、README 一律用**繁體中文**。文件裡的相對路徑要寫明是相對於哪個資料夾。
- **已完成的影片不要刪除或覆寫**。新影片另外命名或放進新資料夾。
- 只用 ManimGL；**不要升級或重裝**。
- 參考文獻 PDF（Chu、Selimi、Fulton）與它們的抄錄僅供個人研讀，**不可散布**（不放進公開的壓縮檔或網路）。
- 盡量節省 token 與步驟。
- **今後的新演示也要加進 combkit**：照 `projects/combkit/README.md` §10 的步驟，然後版本號加一版，並重新執行 `tools/make_zip.py`。
