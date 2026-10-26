"""
combkit —— 組合數學＋ManimGL 的工具包（由 catalan930 專案整理出來，之後各計畫的共用程式都收進這裡）

    combkit.counting   卡塔蘭家族的列舉與對應（Dyck 路徑、二元樹、三角剖分、Raney、MIT 18.212 補充）（純 Python）
    combkit.tableaux   Young 表：插入、乘積、jeu de taquin、Schur 多項式、Knuth 等價、L(w,k)、R-S-K、矩陣球、LR 規則、對稱函數（純 Python）
    combkit.rep_data   對稱群表示的資料與計算（tabloid、polytabloid、表示矩陣、特徵標）（純 Python）
    combkit.raney      Raney (1960)：循環引理、字與樹的計數 L(n;a)、Lagrange 反演、整數序列的 Raney 引理（純 Python）
    combkit.raney_examples 例題產生器：均勻隨機的字／樹／不交叉分法、Markdown 例題集（純 Python）
    combkit.rmt_moments 隨機矩陣的動差：Wick 配對、GOE／GUE 精確值、1/N 展開、樹形閉路徑（純 Python）
    combkit.schubert   Fulton 第 9–10 章：Bruhat 序、Schubert 胞腔、Schubert 多項式、pipe dream（純 Python）
    combkit.fulton_appendix Fulton 附錄 A–B：evacuation、鑰匙、charge、−1 現象（純 Python）
    combkit.walk1d     一維隨機漫步（Rudnick & Gaspari 第 1 章）：走法數、二項分布、Stirling、高斯近似（純 Python）
    combkit.walk_gf    格子上的漫步與生成函數（第 2 章）：走法數、結構函數、布里淵區積分、Cauchy 步伐（純 Python）
    combkit.walk_rec   回家與維度（第 3 章）：回原點的生成函數、Watson 積分、逃脫機率（純 Python）
    combkit.walk_bc    邊界條件與穩定態（第 4 章）：鏡像走路人、吸收牆、存活機率、1/r 與 Berg–Purcell（純 Python）
    combkit.walk_var   漫步的變形（第 5 章）：持續漫步、電報方程、Laplace 變換（Python＋numpy／scipy.special）
    combkit.walk_shape 漫步的形狀（第 6 章）：迴轉半徑張量、Kramers 恆等式、非球度、1/d 展開、Fixman（Python＋numpy／scipy.integrate）
    combkit.walk_saw   路徑積分與自迴避（第 7 章）：二維自迴避漫步的精確列舉、1−f 展開、Dyson 方程、Flory（Python＋numpy／scipy）
    combkit.walk_scale 標度入門（第 8 章）：長尾巴步伐、標度塌縮、Fisher 指數、比值法、上臨界維度（Python＋numpy／scipy）
    combkit.walk_crit  臨界現象（第 9 章）：平均場與 Onsager、小 Ising 的列舉、臨界指數、Widom、Landau、n → 0（Python＋numpy／scipy）
    combkit.walk_cond  O(n) 模型的平均場與自旋波（第 10 章）：環面上的自迴避漫步、自旋波矩陣、牆邊的 tanh（Python＋numpy／scipy）
    combkit.walk_rg    標度、碎形與重整化（第 11 章）：Koch／Cantor 的維度、Gaussian 與 Ising 鏈的抽格、動量殼的流（Python＋numpy／scipy）
    combkit.walk_eps   重整化群的 ε 展開（第 12 章）：積掉自由度、連結叢集、O(n) 縮併、Wilson–Fisher 不動點、ν（Python＋numpy／scipy）
    combkit.young_extras youngart 桌布的計算：patience、滑動路徑、擠壓路徑、影子線、快速 R-S-K、pipe 水管（Python＋numpy）
    combkit.sym_mat    Harter《Principles of Symmetry, Dynamics, and Spectroscopy》第 1 章：分析器與 Axiom 4、么正性、換基底、
                       Hamilton–Cayley 投影算符與譜分解、最小方程、對易矩陣、散射與密度矩陣、耦合振子、預解式的圍道積分（Python＋numpy／scipy／sympy）
    combkit.sym_abel   Harter PSDS 第 2 章：Abel 對稱群、投影算符（最小方程＝特徵標）、耦合擺與單擺環、能帶與對稱破缺、奔馳波、
                       點群、C₁₂ 轉子與 Bloch 波、Kronig–Penney、兩能階、方井 S 矩陣與穿越矩陣（Python＋numpy／scipy／sympy）
    combkit.sym_nonabel Harter PSDS 第 3 章：有限群工具（類代數、冪等元、基本算子、不可約表示、特徵標表、正則表示、相關表）、
                       C3v／D4／D6／D_nh、Hamilton 轉（四元數）、彈簧質量的簡正模、D_n 轉子能帶、氨的振動（Python＋numpy／scipy／sympy）
    combkit.sym_oct    Harter PSDS 第 4 章：八面體群 O／Oh／Td／T／Th、四條子群鏈的不可約表示、相關表、誘導表示與 Frobenius 互反、
                       六點穿隧模型與外場、XF6 的簡正模與同位素位移、升降算子與簡併度（Python＋numpy／scipy／sympy）
    combkit.sym_rot    Harter PSDS 第 5 章：O3／O2 與 D∞h、Euler ↔ Darboux ↔ 四元數、生成元、任意 j 的 J 與 Wigner 𝒟^j（三種算法）、
                       Stern–Gerlach、量子轉子、球諧、O3↓Oh／D6 相關與助記輪、SF6 能階叢、射線表示、R4 與 U3（Python＋numpy／scipy／sympy）
    combkit.sym_prod   Harter PSDS 第 6 章：外積群與對角子群、投影法 Clebsch–Gordan 係數（書上的相位、各子群鏈）、irrep 三型、張量與彈性常數、
                       C3v 的場、參數定理與 Wigner–Eckart、Lorentz 振子與介電函數、XF6 紅外／Raman、泛音、E⊗e Jahn–Teller（Python＋numpy／scipy／sympy）
    combkit.sym_coup   Harter PSDS 第 7 章：R3 耦合係數（降算子、Racah、Schwinger 玻色子）、3j／6j、超精細與 (2p)²、張量算子與單位張量、R3 Wigner–Eckart 與晶體場、
                       剛體與八面體轉子的 RE 曲面與能階叢、旋量與 Stokes 向量、Lissajous 與局域模、H2+／H2、混成與 Hückel（Python＋numpy／scipy／sympy）
    combkit.sym_dyn    Harter PSDS 第 8 章：接觸轉換與作用量、拋射體的波前、相干態與 BCH、快照合成本徵函數、自相關與 EBK、A·p／E·r、譜函數與 TRK、
                       Rabi、Bloch–Siegert（Floquet）、Bloch 方程與穿衣態、Einstein A、Jaynes–Cummings 與 Mollow（Python＋numpy／scipy／sympy）
    combkit.symosc     對稱群與耦合振盪子（symosc 專案）：D3 約化與拍、破缺到 C3（帶相位耦合、陀螺項、進動）／C2、Q8 的左右標籤與四元數型（Python＋numpy）
    combkit.tikz       用 Python 產生 TikZ 圖（純 Python）
    combkit.tikz_trees 平面樹的 TikZ 畫法（純 Python）
    combkit.draw.*     ManimGL 的畫圖元件（需要 ManimGL；sym_draw 是 Harter PSDS 的矩陣、分析器、彈簧鏈、單擺環、相量鐘、群表、Hamilton 球、氨分子、立方體、八面體振幅、能階圖、XF6、Euler 球、向量錐、機率長條、權圖、助記輪、兩粒子八面體、應變立方體、Smith 圖、墨西哥帽、自旋陀螺、耦合錐、RE 曲面、Poincaré 球、偏振橢圓、位能曲線、軌跡與波前、相空間相量、波包快照、曲柄球、穿衣能階梯）

匯入 combkit 本身不會載入 manim；只有 import combkit.draw 底下的模組才需要。
（1.27.0：這段說明之前存檔時被轉碼弄壞，已重寫。）
"""
__version__ = "1.35.2"
